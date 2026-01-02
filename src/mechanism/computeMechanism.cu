#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <iomanip>
#include <map>

struct Committors {
  double qf; // forward committor mean
  double qb; // backward committor mean
  double pA; // population of A state
  double pAB; // population of AB state
  double pB; // population of B state

  double Eqf; // forward committor mean
  double Eqb; // backward committor mean
  double EpA; // population of A state
  double EpAB; // population of AB state
  double EpB; // population of B state
};

struct KramersRate {
  double mean;   // geometric mean
  double se;     // standard error
  double bamean; // block-averaged weighted mean
  double base;   // block-averaged weighted SE
};

struct TPTinfo {
  size_t fpt; // first passage time to the boundary
  size_t let; // last exit time from the boundary
  int loc; // 0 if particle is inside of the boundary (small r), 1 if it is out
};

struct Stats {
  double mean;
  double error;
};

struct TwoHistogram{
  std::vector<double> x;
  std::vector<double> pA;
  std::vector<double> pB;
};


struct Sample {
  double r;
  double cos;
};

struct infoB{ // center id, counter id, time when nndist touched B from reactive trajectory
  int centerID;
  int counterID;
  int time; // after removing eqtime frames
};
struct Box {
  float x;
  float y;
  float z;
};
  

using Trajectory = std::vector<Sample>;
using Trajectories = std::vector<std::vector<Sample>>;


// ----------------------- Stat ---------------------------------
double mean(const std::vector<double>& data) {
  double sum = 0.0;
  for (double x : data) sum += x;
  return sum/data.size();
}

double variance ( const std::vector<double>& data) {
  double mu = mean(data);
  double var = 0.0;
  for ( double x : data) var += (x-mu) * (x-mu) ;
  return var/( data.size() -1);
}

double blockAverage( const std::vector<double>& data, int nB) {
  std::vector<double> blockData = data;
  std::vector<double> semList;
  std::vector<int> blockSizes;

  int blockSize = 1;
  int level=0;

  while ( level <= nB) {
    int nBlocks = blockData.size() / 2;
    // new block
    std::vector<double> newBlockData(nBlocks);
    for ( int i =0 ; i<nBlocks; i++ ) {
      newBlockData[i] = 0.5 * ( blockData[2*i] + blockData[2*i+1] );
    }
    double var = variance(newBlockData);
    double sem = std::sqrt(var / (nBlocks)); 

    semList.push_back(sem);
    blockSizes.push_back(blockSize);
  //  std::cout << blockSize << "\t\t" << nBlocks << "\t\t" << std::setprecision(6) << sem << "\n";

    blockData = newBlockData;
    blockSize *= 2;
    ++level;
  }
  return semList.back();
}

Stats weightedMean ( std::vector<double>& means, std::vector<double>& errors){
  double mean, error;
  for( int i =0 ; i< means.size(); i++ ) {
    double weight = 1/errors[i]/errors[i];
    error += weight;
    mean += weight * means[i];
  }
  mean = mean / error;
  error = std::sqrt( 1/ error);
  return Stats{mean, error};
}


Stats unweightedMean ( std::vector<double>& means){
  double average, error;
  average = mean(means);
  error = std::sqrt( variance(means) / static_cast<double> (means.size())) ;
  return Stats{average, error};
}

// ----------------------------



std::vector<std::vector<TPTinfo>> findPassageTime(double rcut, size_t numatoms, size_t numsnap, std::vector<std::vector<double>>& dist) {
  std::vector<std::vector<TPTinfo>> out(numatoms);
  for (auto& row : out ) row.reserve(numsnap);
  for ( size_t traj = 0; traj < numatoms; traj++ ) {
    std::vector<size_t> fpt(numsnap, 0), let(numsnap, 0);
    std::vector<int> loc(numsnap, 0); 
    loc[0] = ( dist[0][traj] > rcut  ? 1 : 0);
    fpt[0] = 0;
    let[0] = 0;
    size_t lasttime=0;
    for( size_t time = 1; time < numsnap; time++ ) {
      const bool above_prev = dist[time-1][traj] > rcut;
      const bool above_curr = dist[time][traj] > rcut;
      const bool crossed = (above_prev != above_curr);
      if (crossed) {
        loc[time] = (loc[time-1] + 1) % 2 ;
        // update first passage and last exit times
        for ( size_t k = lasttime; k<time; k++ ) fpt[k] = time;
        for (size_t k = lasttime+1; k<time; k++ ) let[k] = let[lasttime];
        let[time] = time;
        lasttime = time;
      } else{ // no trajsition
        loc[time] = loc[time-1];
      }
    }
    for ( size_t kk = lasttime; kk < numsnap; kk++) { // treating last end of trajectory
      fpt[kk] = 0;
      let[kk] = let[lasttime];
    }

    auto& row = out[traj];
    for ( size_t time=0; time< numsnap; time++ ) {
      row.emplace_back( TPTinfo{fpt[time], let[time], loc[time]});
    }
  }
  return out;
}


Stats findKramersRate(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary, size_t numatoms, size_t numsnap, int numBlock){
  std::vector<double> meanKr; meanKr.reserve(numatoms*2);
  std::vector<double> baseKr; baseKr.reserve(numatoms*2);
  for ( size_t traj = 0; traj < numatoms ; traj++ ) {
    std::vector<double> nuAB; nuAB.reserve(numsnap);
    std::vector<double> nuBA; nuBA.reserve(numsnap);
    for(size_t time=1; time < numsnap; time++ ) {
      const bool welldefinedAB = (outerBoundary[traj][time-1].let != 0); // assume analysis starts at the first cross
      const bool crossB = (innerBoundary[traj][time].loc == 0 ) && ( innerBoundary[traj][time-1].loc == 1);
      const bool cameFromA = innerBoundary[traj][time-1].let  < outerBoundary[traj][time-1].let; 
      if (welldefinedAB) {
        if (crossB && cameFromA ) nuAB.push_back(1);
        else nuAB.push_back(0);
      }

      const bool welldefinedBA = (innerBoundary[traj][time-1].let != 0);
      const bool crossA = (outerBoundary[traj][time].loc == 1 ) && ( outerBoundary[traj][time-1].loc == 0);
      const bool cameFromB = innerBoundary[traj][time-1].let  > outerBoundary[traj][time-1].let; 
      if (welldefinedBA) {
        if (crossA && cameFromB) nuBA.push_back(1);
        else nuBA.push_back(0);
      }
    }
    /*
    double meanAB = mean(nuAB);
    if (!std::isnan(meanAB)) meanKr.push_back(meanAB);
    else meanKr.push_back(0);

    double meanBA = mean(nuBA);
    if (!std::isnan(meanBA)) meanKr.push_back(meanBA);
    else meanKr.push_back(0);

    double baseAB = blockAverage(nuAB, numBlock);
    if (!std::isnan(baseAB)) baseKr.push_back(baseAB);

    double baseBA = blockAverage(nuBA, numBlock);
    if (!std::isnan(baseBA)) baseKr.push_back(baseBA);
    */
    meanKr.push_back(mean(nuAB));
    meanKr.push_back(mean(nuBA));
    baseKr.push_back(blockAverage(nuAB, numBlock));
    baseKr.push_back(blockAverage(nuBA, numBlock));
  }
  // erase NaN elements
  size_t idx=0;
  baseKr.erase(std::remove_if(baseKr.begin(), baseKr.end(), [&](double x) {return std::isnan(meanKr[idx++]);} ), baseKr.end());
  meanKr.erase(std::remove_if(meanKr.begin(), meanKr.end(), [](double x) {return std::isnan(x);} ), meanKr.end());

  Stats out;
  if( std::none_of(meanKr.begin(), meanKr.end(), [](double x) { return x==0;} ) ){
    out = weightedMean(meanKr, baseKr);
  } else {
    out = unweightedMean(meanKr);
  }
  return out;
}


Committors findCommittors(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary ) {
  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();
  std::vector<double> countqf; countqf.reserve(numatoms * numsnap);
  std::vector<double> countqb; countqb.reserve(numatoms * numsnap);
  std::vector<double> countpA; countpA.reserve(numatoms * numsnap);
  std::vector<double> countpAB;countpAB.reserve(numatoms * numsnap);
  std::vector<double> countpB; countpB.reserve(numatoms * numsnap);

  for( size_t atom =0; atom < numatoms; atom++) {
    for ( size_t time=0; time < numsnap; time++) {
      const bool domainA = outerBoundary[atom][time].loc == 1;
      const bool domainB = innerBoundary[atom][time].loc == 0;
      const bool domainAB = (!domainA) && (!domainB);
      const bool gotoB = innerBoundary[atom][time].fpt < outerBoundary[atom][time].fpt;
      const bool cameFromA = innerBoundary[atom][time].let < outerBoundary[atom][time].let;
      const bool fptdefined = innerBoundary[atom][time].fpt != 0 && outerBoundary[atom][time].fpt != 0;
      const bool letdefined = innerBoundary[atom][time].let != 0 && outerBoundary[atom][time].let != 0;

      if (fptdefined) {
        if (gotoB) countqf.push_back(1);
        else countqf.push_back(0);
      }

      if (letdefined) {
        if (cameFromA) countqb.push_back(1);
        else countqb.push_back(0);
      }
      if ( domainA) countpA.push_back(1);
      else countpA.push_back(0);
      if ( domainB) countpB.push_back(1);
      else countpB.push_back(0);
      if ( domainAB) countpAB.push_back(1);
      else countpAB.push_back(0);
    }
  }
  Committors out;
  out.qf = mean(countqf);
  out.qb = mean(countqb);
  out.pA = mean(countpA);
  out.pAB= mean(countpAB);
  out.pB = mean(countpB);

  out.Eqf = std::sqrt(variance(countqf) / static_cast<double>(countqf.size()) );
  out.Eqb = std::sqrt(variance(countqb) / static_cast<double>(countqb.size()) );
  out.EpA = std::sqrt(variance(countpA) / static_cast<double>(countpA.size()) );
  out.EpAB= std::sqrt(variance(countpAB)/ static_cast<double>(countpAB.size()));
  out.EpB = std::sqrt(variance(countpB) / static_cast<double>(countpB.size()) );

  return out;
}


Stats findqpqm(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary ) {
  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();
  std::vector<double> count; count.reserve(numatoms * numsnap);

  for( size_t atom =0; atom < numatoms; atom++) {
    for ( size_t time=0; time < numsnap; time++) {
      const bool domainA = outerBoundary[atom][time].loc == 1;
      const bool domainB = innerBoundary[atom][time].loc == 0;
      const bool domainAB = (!domainA) && (!domainB);
      const bool gotoB = innerBoundary[atom][time].fpt < outerBoundary[atom][time].fpt;
      const bool cameFromA = innerBoundary[atom][time].let < outerBoundary[atom][time].let;
      const bool fptdefined = innerBoundary[atom][time].fpt != 0 && outerBoundary[atom][time].fpt != 0;
      const bool letdefined = innerBoundary[atom][time].let != 0 && outerBoundary[atom][time].let != 0;

      if (fptdefined && letdefined) {
        if (gotoB && cameFromA) count.push_back(1);
        else count.push_back(0);
      }
    }
  }
  Stats out;
  out.mean = mean(count);
  out.error = std::sqrt(variance(count) / static_cast<double>(count.size()) );

  return out;
}


Committors findCommittorsBA(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary, int numblock) {
  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();
  std::vector<double> meanqf; meanqf.reserve( numatoms);
  std::vector<double> baseqf; baseqf.reserve( numatoms);
  std::vector<double> meanqb; meanqb.reserve( numatoms);
  std::vector<double> baseqb; baseqb.reserve( numatoms);
  
  /*
  std::vector<double> meanpA; meanpA.reserve( numatoms);
  std::vector<double> basepA; basepA.reserve( numatoms);
  std::vector<double> meanpAB; meanpAB.reserve( numatoms);
  std::vector<double> basepAB; basepAB.reserve( numatoms);
  std::vector<double> meanpB; meanpB.reserve( numatoms);
  std::vector<double> basepB; basepB.reserve( numatoms);
  */

  for( size_t atom =0; atom < numatoms; atom++) {
    std::vector<double> countqf; countqf.reserve( numsnap);
    std::vector<double> countqb; countqb.reserve( numsnap);
    /*
    std::vector<double> countpA; countpA.reserve( numsnap);
    std::vector<double> countpB; countpB.reserve( numsnap);
    std::vector<double> countpAB; countpAB.reserve( numsnap);
    */

    for ( size_t time=0; time < numsnap; time++) {
      const bool domainA = outerBoundary[atom][time].loc == 1;
      const bool domainB = innerBoundary[atom][time].loc == 0;
      const bool domainAB = (!domainA) && (!domainB);
      const bool gotoB = innerBoundary[atom][time].fpt < outerBoundary[atom][time].fpt;
      const bool cameFromA = innerBoundary[atom][time].let < outerBoundary[atom][time].let;
      const bool fptdefined = innerBoundary[atom][time].fpt != 0 && outerBoundary[atom][time].fpt != 0;
      const bool letdefined = innerBoundary[atom][time].let != 0 && outerBoundary[atom][time].let != 0;

      if (fptdefined) {
        if (gotoB) countqf.push_back(1);
        else countqf.push_back(0);
      }

      if (letdefined) {
        if (cameFromA) countqb.push_back(1);
        else countqb.push_back(0);
      }
      /*
      if(domainA) countpA.push_back(1);
      else countpA.push_back(0);
      if(domainB) countpB.push_back(1);
      else countpB.push_back(0);
      if(domainAB) countpAB.push_back(1);
      else countpAB.push_back(0);
      */
    }
    meanqf.push_back( mean(countqf));
    baseqf.push_back(blockAverage(countqf, numblock));
    meanqb.push_back( mean(countqb));
    baseqb.push_back(blockAverage(countqb, numblock));
    /*
    meanpA.push_back( mean(countpA));
    basepA.push_back(blockAverage(countpA, numBlock));
    meanpB.push_back( mean(countpB));
    basepB.push_back(blockAverage(countpB, numBlock));
    meanpAB.push_back( mean(countpAB));
    basepAB.push_back(blockAverage(countpAB, numBlock));
    */
  }

  size_t idx=0;
  baseqf.erase(std::remove_if(baseqf.begin(), baseqf.end(), [&](double x) {return std::isnan(meanqf[idx++]);} ), baseqf.end());
  meanqf.erase(std::remove_if(meanqf.begin(), meanqf.end(), [](double x) {return std::isnan(x);} ), meanqf.end());

  idx=0;
  baseqb.erase(std::remove_if(baseqb.begin(), baseqb.end(), [&](double x) {return std::isnan(meanqb[idx++]);} ), baseqb.end());
  meanqb.erase(std::remove_if(meanqb.begin(), meanqb.end(), [](double x) {return std::isnan(x);} ), meanqb.end());


  Stats qf, qb;
  /*
  if( std::none_of(meanqf.begin(), meanqf.end(), [](double x) { return x==0 || x==1;} ) ){
    qf = weightedMean(meanqf, baseqf);
  } else {
    qf = unweightedMean(meanqf);
  }
  if( std::none_of(meanqb.begin(), meanqb.end(), [](double x) { return x==0 || x==1;} ) ){
    qb = weightedMean(meanqb, baseqb);
  } else {
    qb = unweightedMean(meanqb);
  }
  */
  qf = unweightedMean(meanqf);
  qb = unweightedMean(meanqb);
  /*
  auto pA = weightedMean(meanpA, basepA);
  auto pB = weightedMean(meanpB, basepB);
  auto pAB = weightedMean(meanpAB, basepAB);
  */

  Committors out;
  out.qf = qf.mean;
  out.Eqf = qf.error;
  out.qb = qb.mean;
  out.Eqb = qb.error;
  /*
  out.pA = pA.mean;
  out.EpA = pA.error;
  out.pB = pB.mean;
  out.EpB = pB.error;
  out.pAB = pAB.mean;
  out.EpAB = pAB.error;
  */

  return out;
}

Stats findMFPTAB(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary){
  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();
  std::vector<double> pathTime;
  for (int atom=0; atom < numatoms; atom++) {
    size_t currentTime=outerBoundary[atom][0].fpt;
    while ( currentTime != 0) {
      const bool cameFromB = innerBoundary[atom][currentTime-1].let > outerBoundary[atom][currentTime-1].let && outerBoundary[atom][currentTime-1].let != 0; // before crossing A, came from B?
      size_t nexttime = innerBoundary[atom][currentTime].fpt;
      const bool gotoB = nexttime != 0;
      if( cameFromB && gotoB) {
        size_t time = nexttime - currentTime;
        pathTime.push_back(static_cast<double>(time));
      }
      currentTime = outerBoundary[atom][currentTime].fpt;
    }
  }
  auto out = unweightedMean(pathTime); 
  return out;
}

Stats findMTPTAB(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary){
  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();
  std::vector<double> pathTime;
  for (int atom=0; atom < numatoms; atom++) {
    size_t currentTime=outerBoundary[atom][0].fpt;
    while ( currentTime != 0) {
      //const bool cameFromB = innerBoundary[atom][currentTime-1].let > outerBoundary[atom][currentTime-1].let && outerBoundary[atom][currentTime-1].let != 0; // before crossing A, came from B?
      const bool cameFromA = innerBoundary[atom][currentTime-1].let < outerBoundary[atom][currentTime-1].let && innerBoundary[atom][currentTime-1].let != 0; // before crossing B, came from A?
      const bool gotoB = innerBoundary[atom][currentTime].fpt < outerBoundary[atom][currentTime].fpt && innerBoundary[atom][currentTime].fpt != 0; // hitting inner boundary first than outer boundary
      size_t nexttime = innerBoundary[atom][currentTime].fpt;
      if( cameFromA && gotoB) {
        size_t time = nexttime - currentTime;
        pathTime.push_back(static_cast<double>(time));
      }
      currentTime = outerBoundary[atom][currentTime].fpt;
    }
  }
  auto out = unweightedMean(pathTime); 
  return out;
}

TwoHistogram findHistAB(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary,
    std::vector<std::vector<double>>& angl){
  // note, angl [ time ] [traj]  vs innerBoundary [ traj ] [ time]
  // angle : degree
  const double pi = 3.14159265358979323846; 

  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();

  // make histogram
  const double minVal = -1.0;
  const double maxVal = 1.0;
  int nbins=50;
  const double width = (maxVal - minVal)/static_cast<double>(nbins);
  std::vector<double> countsA(nbins, 0.0);
  std::vector<double> countsB(nbins, 0.0);

  double NA=0.0;
  double NB=0.0;

  for (int atom=0; atom < numatoms; atom+=2) { // only for cation 
    size_t currentTime=outerBoundary[atom][0].fpt;
    while ( currentTime != 0) {
      const bool cameFromA = innerBoundary[atom][currentTime-1].let < outerBoundary[atom][currentTime-1].let && innerBoundary[atom][currentTime-1].let != 0; // before crossing B, came from A?
      const bool gotoB = innerBoundary[atom][currentTime].fpt < outerBoundary[atom][currentTime].fpt && innerBoundary[atom][currentTime].fpt != 0; // hitting inner boundary first than outer boundary
      size_t nextTime = innerBoundary[atom][currentTime].fpt;
      if( cameFromA && gotoB) {
        double angleA = std::cos( angl[currentTime][atom] / 180.0 * pi);
        double angleB = std::cos( angl[nextTime][atom] / 180.0 * pi);

        if (angleA <= 1.0 && angleA >= -1.0 ){ 
          size_t idxA = static_cast<std::size_t>((angleA-minVal)/ width);
          if (idxA >= nbins) idxA = nbins-1;
          countsA[idxA] += 1.0;
          NA+=1.0;
        }

        if (angleB <= 1.0 && angleB >= -1.0 ){ 
          size_t idxB = static_cast<std::size_t>((angleB-minVal)/ width);
          if (idxB >= nbins) idxB = nbins-1;
          countsB[idxB] += 1.0;
          NB+=1.0;
        }
      }
      currentTime = outerBoundary[atom][currentTime].fpt;
    }
  }

  TwoHistogram out;
  out.x.resize(nbins);
  out.pA.resize(nbins);
  out.pB.resize(nbins);


  for( size_t i =0; i< nbins; i++) {
    double center = minVal + (i+0.5) * width;
    out.x[i] = center;
    out.pA[i] = countsA[i] /(NA*width);
    out.pB[i] = countsB[i] /(NB*width);
  }
  return out;
}



TwoHistogram findHistBA(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary,
    std::vector<std::vector<double>>& angl){
  // note, angl [ time ] [traj]  vs innerBoundary [ traj ] [ time]
  // angle : degree
  const double pi = 3.14159265358979323846; 

  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();

  // make histogram
  const double minVal = -1.0;
  const double maxVal = 1.0;
  int nbins=50;
  const double width = (maxVal - minVal)/static_cast<double>(nbins);
  std::vector<double> countsA(nbins, 0.0);
  std::vector<double> countsB(nbins, 0.0);

  double NA=0.0;
  double NB=0.0;

  for (int atom=0; atom < numatoms; atom+=2) { // only for cation 
    size_t currentTime=innerBoundary[atom][0].fpt;
    while ( currentTime != 0) {
      const bool cameFromB = innerBoundary[atom][currentTime-1].let > outerBoundary[atom][currentTime-1].let && outerBoundary[atom][currentTime-1].let != 0; // before crossing A, came from B?
      const bool gotoA = innerBoundary[atom][currentTime].fpt > outerBoundary[atom][currentTime].fpt && outerBoundary[atom][currentTime].fpt != 0;
      size_t nextTime = outerBoundary[atom][currentTime].fpt;
      if( cameFromB && gotoA) {
        double angleA = std::cos( angl[nextTime][atom] / 180.0 * pi);
        double angleB = std::cos( angl[currentTime][atom] / 180.0 * pi);

        if (angleA <= 1.0 && angleA >= -1.0 ){ 
          size_t idxA = static_cast<std::size_t>((angleA-minVal)/ width);
          if (idxA >= nbins) idxA = nbins-1;
          countsA[idxA] += 1.0;
          NA+=1.0;
        }

        if (angleB <= 1.0 && angleB >= -1.0 ){ 
          size_t idxB = static_cast<std::size_t>((angleB-minVal)/ width);
          if (idxB >= nbins) idxB = nbins-1;
          countsB[idxB] += 1.0;
          NB+=1.0;
        }
      }
      currentTime = innerBoundary[atom][currentTime].fpt;
    }
  }

  TwoHistogram out;
  out.x.resize(nbins);
  out.pA.resize(nbins);
  out.pB.resize(nbins);


  for( size_t i =0; i< nbins; i++) {
    double center = minVal + (i+0.5) * width;
    out.x[i] = center;
    out.pA[i] = countsA[i] /(NA*width);
    out.pB[i] = countsB[i] /(NB*width);
  }
  return out;
}


std::vector<infoB> getTrajectoryAB(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary,
    std::vector<std::vector<int>>& cA, std::vector<std::vector<int>>& cB){
  // note, angl [ time ] [traj]  vs innerBoundary [ traj ] [ time]
  // angle : degree
  const double pi = 3.14159265358979323846; 

  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();
  
  std::vector<infoB> out;
  out.reserve(20000);

  for (int atom=0; atom < numatoms; atom+=2) { // only for cation 
    size_t currentTime=outerBoundary[atom][0].fpt;
    while ( currentTime != 0) {
      const bool cameFromA = innerBoundary[atom][currentTime-1].let < outerBoundary[atom][currentTime-1].let && innerBoundary[atom][currentTime-1].let != 0; // before crossing B, came from A?
      const bool gotoB = innerBoundary[atom][currentTime].fpt < outerBoundary[atom][currentTime].fpt && innerBoundary[atom][currentTime].fpt != 0; // hitting inner boundary first than outer boundary
      size_t nextTime = innerBoundary[atom][currentTime].fpt;
      if( cameFromA && gotoB) {
        infoB one;
        one.centerID = cA[nextTime][atom];
        one.counterID = cB[nextTime][atom];
        one.time=nextTime;
        out.push_back(one);
      }
      currentTime = outerBoundary[atom][currentTime].fpt;
    }
  }
  return out;
}

std::vector<infoB> getTrajectoryBA(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary,
    std::vector<std::vector<int>>& cA, std::vector<std::vector<int>>& cB){
  // note, angl [ time ] [traj]  vs innerBoundary [ traj ] [ time]
  // angle : degree
  const double pi = 3.14159265358979323846; 

  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();
  
  std::vector<infoB> out;
  out.reserve(20000);

  for (int atom=0; atom < numatoms; atom+=2) { // only for cation 
    size_t currentTime=innerBoundary[atom][0].fpt;
    while ( currentTime != 0) {
      const bool cameFromB = innerBoundary[atom][currentTime-1].let > outerBoundary[atom][currentTime-1].let && outerBoundary[atom][currentTime-1].let != 0; // before crossing A, came from B?
      const bool gotoA = innerBoundary[atom][currentTime].fpt > outerBoundary[atom][currentTime].fpt && outerBoundary[atom][currentTime].fpt != 0;
      size_t nextTime = outerBoundary[atom][currentTime].fpt;
      if( cameFromB && gotoA) {
        infoB one;
        one.centerID = cA[nextTime][atom];
        one.counterID = cB[nextTime][atom];
        one.time=nextTime;
        out.push_back(one);
      }
      currentTime = innerBoundary[atom][currentTime].fpt;
    }
  }
  return out;
}

Stats findMFPTBA(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary){
  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();
  std::vector<double> pathTime;
  for (int atom=0; atom < numatoms; atom++) {
    size_t currentTime=innerBoundary[atom][1].fpt;
    while ( currentTime != 0) {
      const bool cameFromA = innerBoundary[atom][currentTime-1].let < outerBoundary[atom][currentTime-1].let && innerBoundary[atom][currentTime-1].let != 0; // before crossing A, came from B?
      size_t nexttime = outerBoundary[atom][currentTime].fpt;
      const bool gotoA = nexttime != 0;
      if( cameFromA && gotoA) {
        size_t time = nexttime - currentTime;
        pathTime.push_back(static_cast<double>(time));
      }
      currentTime = innerBoundary[atom][currentTime].fpt;
    }
  }
  auto out = unweightedMean(pathTime); 
  return out;
}


struct ReadError : std::runtime_error { 
  using std::runtime_error::runtime_error;
};

std::vector<std::vector<double>> readBTrajFromFile(const std::string& filename) {
  std::ifstream in(filename, std::ios::binary);
  if (!in) throw ReadError("Cannot open file: " + filename);
  std::int32_t rows=0, cols=0;
  in.read(reinterpret_cast<char*>(&rows), sizeof(rows));
  in.read(reinterpret_cast<char*>(&cols), sizeof(cols));
  if(!in) throw ReadError("Failed to read header (rows/cols)");
  if( rows < 0 || cols < 0) throw ReadError("Wrong rows/cols in header");
  if (rows == 0 || cols == 0) return {};

  std::vector<std::vector<double>> matrix(static_cast<size_t>(rows), std::vector<double>(static_cast<size_t>(cols)));
  
  std::vector<float> rowf(static_cast<size_t>(cols));

  for( std::int32_t i=0; i<rows; i++){
    in.read(reinterpret_cast<char*>(rowf.data()), static_cast<std::streamsize>(rowf.size() * sizeof(float)));
    if (!in) throw ReadError("Short read at row ");

    auto& rowd = matrix[static_cast<size_t>(i)];
    std::transform(rowf.begin(), rowf.end(), rowd.begin(),
        [](float x) {return static_cast<double>(x); });
   }
  in.close();
//  std::cout << "file reading done" << std::endl;
  return matrix;
}

std::vector<std::vector<int>> readBTrajFromFileInt(const std::string& filename) {
  std::ifstream in(filename, std::ios::binary);
  if (!in) throw ReadError("Cannot open file: " + filename);
  std::int32_t rows=0, cols=0;
  in.read(reinterpret_cast<char*>(&rows), sizeof(rows));
  in.read(reinterpret_cast<char*>(&cols), sizeof(cols));
  if(!in) throw ReadError("Failed to read header (rows/cols)");
  if( rows < 0 || cols < 0) throw ReadError("Wrong rows/cols in header");
  if (rows == 0 || cols == 0) return {};

  std::vector<std::vector<int>> matrix(static_cast<size_t>(rows), std::vector<int>(static_cast<size_t>(cols)));
  
  std::vector<int> rowf(static_cast<size_t>(cols));

  for( std::int32_t i=0; i<rows; i++){
    in.read(reinterpret_cast<char*>(rowf.data()), static_cast<std::streamsize>(rowf.size() * sizeof(int)));
    if (!in) throw ReadError("Short read at row ");

    auto& rowd = matrix[static_cast<size_t>(i)];
    std::transform(rowf.begin(), rowf.end(), rowd.begin(),
        [](int x) {return static_cast<int>(x); });
   }
  in.close();
//  std::cout << "file reading done" << std::endl;
  return matrix;
}

std::vector<double> linspace(double start, double stop, int num, bool endpoint = true) {
  std::vector<double> result;
  if (num <= 0) {
    return result;
  }

  double step = (stop-start) / (endpoint ? num-1 : num);

  for (int i=0; i< num ; i++){
    double exponent = start + i *step;
    result.push_back(exponent);
  }
  return result;
}



void generateArgumentVector(std::vector<double>& argVector, int numBins, double minVal, double maxVal) {
  double binWidth = (maxVal - minVal) / numBins;
  for (int i = 0; i<numBins; i++){
    argVector.push_back(minVal + binWidth * ( i + 0.5f));
  }
}


void printFileInt(std::string& file, std::vector<int>& x, std::vector<int>& y){
  std::cout << "Outpul File generated: " << file << std::endl;
  std::cout << "Size: " << x.size() <<  std::endl;

  std::ofstream out ( file);
  for (size_t i = 0; i< x.size(); i++) {
    out << std::fixed << std::setprecision(7) << x[i] << " " << y[i]  << std::endl;
  }
  out.close();
}
void printFile(std::string& file, std::vector<double>& x, std::vector<double>& y){
  std::cout << "Outpul File generated: " << file << std::endl;
  std::cout << "Size: " << x.size() <<  std::endl;

  std::ofstream out ( file);
  for (size_t i = 0; i< x.size(); i++) {
    out << std::fixed << std::setprecision(12) << x[i] << " " << y[i]  << std::endl;
  }
  out.close();
}


void readTrajFromFile(std::string& file,
                      std::vector<std::vector<double>>& traj) {
  std::ifstream in(file);

  if (!in.is_open()) {
    std::cerr << "Error: Cannot find file" << std::endl;
  }

  std::string line;
  while (std::getline(in, line)) {
    std::istringstream lineStream(line);
    std::vector<double> row;

    double value;

    while (lineStream >> value) {
      row.push_back(value);
    }
    
    if (!row.empty()){
      traj.push_back(row);
    }
  }
  in.close();
  std::cout << "Read File from " << file << " is done" << std::endl;
  std::cout << "Rows : " << traj.size() << " Cols : " << traj[0].size() << std::endl;
}

std::vector<double> logspace(double start, double stop, int num, bool endpoint = true) {
  std::vector<double> result;
  if (num <= 0) {
    return result;
  }

  double step = (stop-start) / (endpoint ? num-1 : num);

  for (int i=0; i< num ; i++){
    double exponent = start + i *step;
    result.push_back(pow(10, exponent));
  }
  return result;
}

template <typename T>
std::vector<T> removeDuplicates(const std::vector<T>& input) {
  std::set<T> uniqueSet(input.begin(), input.end());
  return std::vector<T>(uniqueSet.begin(), uniqueSet.end());
}


double fieldConvert(const std::string& s) {
  if (s.find_first_not_of('0') == std::string::npos) {
    return 0.0f;
  }
  size_t leadingZeros = s.find_first_not_of('0');

  if ( leadingZeros >= 2) {
    std::string digits = s.substr(leadingZeros);
    double value = std::stof(digits);
    return value / std::pow(10, leadingZeros -1);
  }
  return std::stof(s);
}

double divisionError(double x, double ex, double y, double ey ) { // return error of x/y 
  if( y == 0) return 0;
  double ratio = x/y;
  double relEx = ex/x;
  double relEy = ey/y;
  return ratio * std::sqrt( relEx * relEx + relEy * relEy );
}

float applyPBC(float x, float box){
  float hbox = box/2.0;
  float wrapped = fmod(x + hbox, box);
  if (wrapped < 0) wrapped += box;
  return wrapped - hbox;
}

float distance(const Box& a, const Box& b, const Box& box) {
  float dx, dy, dz, rsq;
  dx=applyPBC( a.x - b.x, box.x);
  dy=applyPBC( a.y - b.y, box.y);
  dz=applyPBC( a.z - b.z, box.z);
  rsq = dx*dx +dy*dy +dz*dz;
  return std::sqrt(rsq);
}

std::vector<std::vector<float>> readBinary(const std::string& filename) {
  std::ifstream in(filename, std::ios::binary);
  if (!in.is_open()) {
    throw std::runtime_error("Failed to open file: " + filename);
  }
  int rows, cols;
  in.read(reinterpret_cast<char*>(&rows), sizeof(int));
  in.read(reinterpret_cast<char*>(&cols), sizeof(int));


  std::vector<std::vector<float>> matrix(rows, std::vector<float>(cols));
  for( int i=0; i<rows; i++){
    in.read(reinterpret_cast<char*>(matrix[i].data()), sizeof(float)*cols);
   }
  in.close();
  return matrix;
}

float computeAngle(const Box& cation, const Box& anion, const Box& box) {
  float dx, dy, dz;
  dx=applyPBC( anion.x - cation.x, box.x);
  dy=applyPBC( anion.y - cation.y, box.y);
  dz=applyPBC( anion.z - cation.z, box.z);
  float norm = std::sqrt( dx*dx + dy*dy +dz*dz);
  if (norm == 0.0){
    std::cout<< "Error! Zero-length vector \n";
  }

  float cosTheta = dz/norm;
  if (cosTheta > 1.0) cosTheta=1.0;
  if (cosTheta <-1.0) cosTheta=-1.0;

  float angleRad = std::acos(cosTheta);
  float angleDeg = angleRad * ( 180.0/M_PI);

  return cosTheta;
}

int main(int argc, char* argv[]) {
  if (argc != 11 ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << "dir_name density field cutoff_in(A) cutoff_out(A) boxX(A) boxY(A) boxZ(A) timestep(ps)  eqtime(ns) \n";
    return 1;
  }
  std::string dirName=argv[1]; // name of directory to output results, in results directory
  std::string rho = argv[2]; 
  //float density = std::stof(rho)* 0.1; // density in M unit
  int n = rho.size();
  long long val = std::stoll(rho);
  float density = static_cast<float>(val) / std::pow(10.0, n-1);
  float boxX = std::stof(argv[6]); // box size in A unit
  float boxY = std::stof(argv[7]);
  float boxZ = std::stof(argv[8]); 
  Box box{boxX, boxY, boxZ};
  float avogadro = 6.0221408;
  int numatoms= static_cast<int>(std::round(2.0*density * boxX * boxY * boxZ *avogadro /10000) ); // cation N + anion N
  int numpion = numatoms/2;
  std::string fieldname = argv[3];
  float field = fieldConvert(fieldname)* 0.023; //kcal /molA,
  float fconversion = 25.7/0.592 ; // kcal/molA to mV/A
  float CUTOFFin = std::stof(argv[4]);
  float CUTOFFout = std::stof(argv[5]);
  float timestep = std::stof(argv[9]); // ps unit
  float eqtime = std::stof(argv[10]); // ns unit


  std::cout << "\nAnalysis Setups\n";
  std::cout << "dump Directory : ../data/dumps" << dirName << "/\n"; 
  std::cout << "output Directory : ../results/" << dirName << "/\n"; 
  std::cout << "numIons: " << numatoms << " numPairs: " << numpion << "\n";
  std::cout << "boxSizes x: " << boxX << " A y: " << boxY << " A z: " << boxZ << " A \n";
  std::cout << "density of ions : " << density << "M \n";
  std::cout << "fieldStrength: " << field << " kcal/molA = " << field * fconversion << " mV/A" << "\n"; 
  std::cout << "domainCutoffs: in " << CUTOFFin << " A\t\tout "<< CUTOFFout << " A\n";
  std::cout << "timeStep: " << timestep << " ps\n";
  std::cout << "eqTime: " << eqtime << " ns\n";
  std::cout << "\n\n";



  // get dist data, flatten it to feed it to paralle job
  std::vector<std::vector<double>> dist;// note numsnap / numatom order
  std::string distFile;
  distFile = std::string("../data/cnnDist/") + dirName + "D" + rho + "E" + fieldname + ".binary";
  dist=readBTrajFromFile(distFile);
  size_t numsnap = dist.size();
  size_t ratoms = dist[0].size();

  // get angle file
  std::vector<std::vector<double>> angl;
  std::string anglFile;
  anglFile = std::string("../data/cnnAngle/") + dirName + "D" + rho  + "E" + fieldname + ".binary";
  angl=readBTrajFromFile(anglFile);
  size_t asnap = angl.size();
  size_t aatoms = angl[0].size();

  // get id
  std::vector<std::vector<int>> cid;
  std::string idFile = std::string("../data/cnnId/") + dirName + "D" + rho + "E" + fieldname + ".binary";
  cid = readBTrajFromFileInt(idFile);
  size_t idsnap = cid.size();
  size_t idatoms = cid[0].size();

  // get cid
  std::vector<std::vector<int>> ccid;
  std::string cidFile = std::string("../data/cnnCid/") + dirName + "D" + rho + "E" + fieldname + ".binary";
  ccid = readBTrajFromFileInt(cidFile);
  size_t cidsnap = ccid.size();
  size_t cidatoms = ccid[0].size();


  if (asnap != numsnap) std::cout << "Error, angle size not match with dist size" << "\n";
  if (idsnap != numsnap) std::cout << "Error, angle size not match with dist size" << "\n";
  if (cidsnap != numsnap) std::cout << "Error, angle size not match with dist size" << "\n";
  if (ratoms != numatoms) std::cout << "Error, angle size not match with dist size" << "\n";
  if (aatoms != numatoms) std::cout << "Error, angle size not match with dist size" << "\n";
  if (idatoms != numatoms) std::cout << "Error, angle size not match with dist size" << "\n";
  if (cidatoms != numatoms) std::cout << "Error, angle size not match with dist size" << "\n";

  std::cout << "TPT-based Hitting probability distribution measurement starting...\n";
  std::cout << "Read Trajectory of length " << numsnap*timestep/1000 << "ns\n\n";

  // evaluate TPT variables for inner boundary
  // Note, each structure has numatoms in row, numsnap in column, opposite to dist
  double rcut = CUTOFFin;
  std::vector<std::vector<TPTinfo>> innerBoundary;
  innerBoundary = findPassageTime(rcut, numatoms, numsnap, dist);

  std::cout << "finding first passage time at " << rcut << " done" << std::endl;

  rcut = CUTOFFout;
  std::vector<std::vector<TPTinfo>> outerBoundary;
  outerBoundary = findPassageTime(rcut, numatoms, numsnap, dist);
  std::cout << "finding first passage time at " << rcut << " done" << std::endl;


  // read dump files (original trajectory files)

  std::cout << "Now read Dumps again \n";
  std::vector<std::vector<float>> traj;
  std::string inputFileName = std::string("../data/dumps") + dirName + "/dumpD" + rho + "E" + fieldname + ".binary";
  std::cout << "READ dump file at " << inputFileName << "\n";
  traj = readBinary(inputFileName);
  if (traj.size() % numatoms != 0) {
    std::cout << "\n\nError in trajectory size or number of atoms\n\n";
  } else {
    numsnap = traj.size()/numatoms;
  }
  std::cout << "Read from dump files " << inputFileName  << " of length " << numsnap * timestep/1000 << " ns\n";
  std::cout << "Removing initial " << eqtime << " ns for relaxation\n";
  int removesnap=static_cast<int>(1000*eqtime/timestep*numatoms);
  traj.erase(traj.begin(), traj.begin()+removesnap);
  numsnap = traj.size()/numatoms;
  std::cout << "Analysis will be on trajectory of length " << numsnap*timestep/1000 << " ns\n\n";

  // trajectory binary accepts only two type of ion ordering.
  //check the order of types
  bool alter;
  if (traj[0][0] != traj[1][0] ) {
    std::cout << "Ions are ordered alternatively (+-+-) \n";
    alter=true;
  } else {
    std::cout << "Ions are ordered sequentially (++--) \n";
    alter=false;
  }

  // work on trajectory AB
  auto partsAB =getTrajectoryAB(innerBoundary, outerBoundary, cid, ccid);
  Trajectories outAB;
  outAB.reserve(partsAB.size());
  for( int i=0; i<partsAB.size(); i++ ) {
    const auto& aset = partsAB[i];
    int time = aset.time;
    int id = aset.centerID;
    int cid = aset.counterID;
    int cA, cB;

    if (alter){
      cA = id;
      cB = cid;
    }
    else {
      cA = id%2==0 ? id/2 : id/2 + numpion;
      cB = cid%2==0 ? cid/2 : cid/2 + numpion;
    }
    bool pass=true;
    int back=1;
    Trajectory trj;
    trj.reserve(30000);
    while (pass) {
      int loc = (time-back)*numatoms;
      Box atomA{traj[loc+cA][1], traj[loc+cA][2], traj[loc+cA][3]};
      Box atomB{traj[loc+cB][1], traj[loc+cB][2], traj[loc+cB][3]};
      float dist = distance(atomA, atomB, box);
      float angl = computeAngle(atomA, atomB, box);
      Sample s;
      s.r = dist;
      s.cos = angl;
      trj.push_back(s);
      if (dist < CUTOFFin) { 
        break;
      }
      if (dist > CUTOFFout) {
        pass=false;
        outAB.push_back(trj);
      }
      back+=1;
    }
  }

  // work on trajectory BA
  auto partsBA =getTrajectoryBA(innerBoundary, outerBoundary, cid, ccid);
  Trajectories outBA;
  outBA.reserve(partsBA.size());
  for( int i=0; i<partsBA.size(); i++ ) {
    const auto& aset = partsBA[i];
    int time = aset.time;
    int id = aset.centerID;
    int cid = aset.counterID;
    int cA, cB;

    if (alter){
      cA = id;
      cB = cid;
    }
    else {
      cA = id%2==0 ? id/2 : id/2 + numpion;
      cB = cid%2==0 ? cid/2 : cid/2 + numpion;
    }
    bool pass=true;
    int back=1;
    Trajectory trj;
    trj.reserve(30000);
    while (pass) {
      int loc = (time-back)*numatoms;
      Box atomA{traj[loc+cA][1], traj[loc+cA][2], traj[loc+cA][3]};
      Box atomB{traj[loc+cB][1], traj[loc+cB][2], traj[loc+cB][3]};
      float dist = distance(atomA, atomB, box);
      float angl = computeAngle(atomA, atomB, box);
      Sample s;
      s.r = dist;
      s.cos = angl;
      trj.push_back(s);
      if (dist > CUTOFFout) { 
        break;
      }
      if (dist < CUTOFFin) {
        pass=false;
        outBA.push_back(trj);
      }
      back+=1;
    }
  }


  /*
  float numreactivetraj=static_cast<float>(out.size());
  float MTPT=0;
  for( const auto& one : out) {
    MTPT+= static_cast<float>(one.size());
  }
  MTPT/=numreactivetraj;
  std::cout << "MTPT : " << MTPT * timestep << " ps  and " << 1/(MTPT*timestep/1000) << " /ns" << "\n";
  std::cout << "Total trajectories : " << ctest << " and survivals " << count << "\n"; 
  std::cout << out.size() << "\t\t" << out[0].size() << "\t\t" << out[1].size() << "\n";
  */

  std::string trajFilerAB = std::string("../results/traj/") + dirName + "D" + rho + "E" + fieldname + "AB.r";
  std::string trajFileaAB = std::string("../results/traj/") + dirName + "D" + rho + "E" + fieldname + "AB.a";
  std::ofstream outrAB(trajFilerAB);
  std::ofstream outaAB(trajFileaAB);

  outrAB << std::fixed << std::setprecision(6);
  outaAB << std::fixed << std::setprecision(6);
  for( size_t i=0; i< outAB.size(); i++) {
    size_t length = outAB[i].size(); 
    for( size_t j=0; j< length ; j++ ) {
      outrAB << outAB[i][j].r << " ";
      outaAB << outAB[i][j].cos << " ";
    }
    outrAB << "\n";
    outaAB << "\n";
  }

  outrAB.close();
  outaAB.close();


  std::string trajFilerBA = std::string("../results/traj/") + dirName + "D" + rho + "E" + fieldname + "BA.r";
  std::string trajFileaBA = std::string("../results/traj/") + dirName + "D" + rho + "E" + fieldname + "BA.a";
  std::ofstream outrBA(trajFilerBA);
  std::ofstream outaBA(trajFileaBA);

  outrBA << std::fixed << std::setprecision(6);
  outaBA << std::fixed << std::setprecision(6);
  for( size_t i=0; i< outBA.size(); i++) {
    size_t length = outBA[i].size(); 
    for( size_t j=0; j< length ; j++ ) {
      outrBA << outBA[i][j].r << " ";
      outaBA << outBA[i][j].cos << " ";
    }
    outrBA << "\n";
    outaBA << "\n";
  }

  outrBA.close();
  outaBA.close();


  /*
  std::string trajFiler = std::string("../results/traj/") + dirName + "D" + rho + "E" + fieldname + ".r";
  std::string trajFilea = std::string("../results/traj/") + dirName + "D" + rho + "E" + fieldname + ".a";
  std::ofstream outr(trajFiler);
  std::ofstream outa(trajFilea);

  outr << std::fixed << std::setprecision(6);
  outa << std::fixed << std::setprecision(6);
  for( size_t i=0; i< 100; i++) {
    size_t length = parts[i].size(); 
    for( size_t j=0; j< length ; j++ ) {
      outr << parts[i][j].r << " ";
      outa << parts[i][j].cos << " ";
    }
    outr << "\n";
    outa << "\n";
  }

  outr.close();
  outa.close();
  */



  return 0;
}

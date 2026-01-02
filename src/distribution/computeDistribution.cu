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


Trajectories getTrajectoryAB(std::vector<std::vector<TPTinfo>>& innerBoundary, std::vector<std::vector<TPTinfo>>& outerBoundary,
    std::vector<std::vector<double>>& dist, std::vector<std::vector<double>>& angl){
  // note, angl [ time ] [traj]  vs innerBoundary [ traj ] [ time]
  // angle : degree
  const double pi = 3.14159265358979323846; 

  size_t numatoms = innerBoundary.size();
  size_t numsnap = innerBoundary[0].size();
  
  Trajectories out;
  out.reserve(20000);

  for (int atom=0; atom < numatoms; atom+=2) { // only for cation 
    size_t currentTime=outerBoundary[atom][0].fpt;
    while ( currentTime != 0) {
      const bool cameFromA = innerBoundary[atom][currentTime-1].let < outerBoundary[atom][currentTime-1].let && innerBoundary[atom][currentTime-1].let != 0; // before crossing B, came from A?
      const bool gotoB = innerBoundary[atom][currentTime].fpt < outerBoundary[atom][currentTime].fpt && innerBoundary[atom][currentTime].fpt != 0; // hitting inner boundary first than outer boundary
      size_t nextTime = innerBoundary[atom][currentTime].fpt;
      if( cameFromA && gotoB) {
        size_t start = currentTime;
        size_t end = nextTime;
        size_t length = end-start+1;
        Trajectory trj;
        trj.reserve(length);
        for(size_t i=0; i < length; i++) {
          Sample s;
          s.r = dist[start+i][atom];
          s.cos = std::cos( angl[start+i][atom] / 180.0 * pi);
          trj.push_back(s);
        }
        out.push_back(trj);
      }
      currentTime = outerBoundary[atom][currentTime].fpt;
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




int main(int argc, char* argv[]) {
  if (argc != 8 ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << " dir_name density field cutoff_in(A) cutoff_out(A) timestep(ps) thermoflag\n";
    return 1;
  }
  std::string dirName=argv[1]; // name of directory to output results, in results directory
  std::string rho = argv[2]; 
  //double density = std::stof(rho)* 0.1; // density in M unit
  int n = rho.size();
  long long val = std::stoll(rho);
  float density = static_cast<float>(val) / std::pow(10.0, n-1);
  std::string fieldname = argv[3];
  double field = fieldConvert(fieldname)* 0.023; //kcal /molA,
  double fconversion = 25.7/0.592 ; // kcal/molA to mV/A
  double CUTOFFin = std::stof(argv[4]);
  double CUTOFFout = std::stof(argv[5]);
  double timestep = std::stof(argv[6]); // ps unit
  int thermoflag=static_cast<int>(std::stoi(argv[7]));// if 1: use partial thermostat dataset, starting with Tdump..


  // get dist data, flatten it to feed it to paralle job
  std::vector<std::vector<double>> dist;// note numsnap / numatom order
  std::string distFile;
  distFile = std::string("../data/cnnDist/") + dirName + "D" + rho + "E" + fieldname + ".binary";
  dist=readBTrajFromFile(distFile);
  size_t numsnap = dist.size();
  size_t numatoms = dist[0].size();

  // get angle file
  std::vector<std::vector<double>> angl;
  std::string anglFile;
  anglFile = std::string("../data/cnnAngle/") + dirName + "D" + rho  + "E" + fieldname + ".binary";
  angl=readBTrajFromFile(anglFile);
  size_t asnap = angl.size();
  size_t aatoms = angl[0].size();
  if (asnap != numsnap) std::cout << "Error, angle size not match with dist size" << "\n";
  if (aatoms != numatoms) std::cout << "Error, angle size not match with dist size" << "\n";

  std::cout << "\nAnalysis Setups\n";
  std::cout << "data location: " << distFile << std::endl ; 
  if(thermoflag == 1) { std::cout << "This trajectory was generated from partial Thermostat \n";}
  else if(thermoflag == 2) { std::cout << "This trajectory was generated from solvent Thermostat \n";}
  std::cout << "numIons: " << numatoms << " numPairs: " << numatoms/2 << "\n";
  std::cout << "density of ions : " << density << "M \n";
  std::cout << "fieldStrength: " << field << " kcal/molA = " << field * fconversion << " mV/A" << "\n"; 
  std::cout << "domainCutoffs: in " << CUTOFFin << " A\t\tout "<< CUTOFFout << " A\n";
  std::cout << "timeStep: " << timestep << " ps\n";
  std::cout << "\n\n";


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


  auto histAB = findHistAB(innerBoundary, outerBoundary, angl);
  auto histBA = findHistBA(innerBoundary, outerBoundary, angl);


  std::string histFile;
  histFile=std::string("../results/distribution/") + dirName + "D" + rho + "E" + fieldname + ".dat";
  std::ofstream out(histFile );
  out << std::fixed << std::setprecision(8);
  for( size_t i=0; i< histAB.x.size(); i++) {
    out << histAB.x[i] << "\t\t" << histAB.pA[i] << "\t\t" << histAB.pB[i] << "\t\t" << histBA.pA[i] << "\t\t" << histBA.pB[i] << "\n";
  }

  out.close();


  auto parts =getTrajectoryAB(innerBoundary, outerBoundary, dist, angl);
  std::cout << parts.size() << "\t\t" << parts[0].size() << "\t\t" << parts[1].size() << "\n"; 
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



  return 0;
}

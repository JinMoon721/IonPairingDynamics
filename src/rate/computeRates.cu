#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <iomanip>
#include <map>

__global__ void computeCommittorsKernel(const int* data1, const int* time1, const int* btime1, const int* data3, const int* time3, const int* btime3, int numSteps, int numTraj, int* counts){
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int den,  currentIndex;
  int n1, n2, n3, n4, n5;
  if (idx < numSteps) {
    den=0;
    n1=n2=n3=n4=n5=0;
    for (int traj = 0 ; traj < numTraj; traj++) {
      currentIndex = idx*numTraj + traj;
      if ( time1[currentIndex] != -1 && btime3[currentIndex] != -1 && time3[currentIndex] != -1 && btime1[currentIndex] != -1 ){ 
        den += 1;
        if ( (time3[currentIndex] > time1[currentIndex] )  ) { // forwardward committor
          n1+=1;
        }
        if ( (btime3[currentIndex] > btime1[currentIndex] )  ) { // backward committor
          n2+=1;
        }

        if (data3[currentIndex] % 2 == 1 ) { // A state
          n3+=1;
        } else if (data1[currentIndex] % 2 == 0) { // B state
          n5+=1;
        } else{
          n4+=1;
        }
      }
    }
    atomicAdd(&counts[0], den);
    atomicAdd(&counts[1], n1);
    atomicAdd(&counts[2], n2);
    atomicAdd(&counts[3], n3);
    atomicAdd(&counts[4], n4);
    atomicAdd(&counts[5], n5);
  }
}


void computeCommittors(const std::vector<int>& h_data1, const std::vector<int>& h_time1, const std::vector<int>& h_btime1, 
                       const std::vector<int>& h_data3, const std::vector<int>& h_time3, const std::vector<int>& h_btime3,  
                       int numSteps, int numTraj, std::vector<float>& result) {
  int *d_btime1, *d_btime3;
  int *d_time1, *d_time3;
  int *d_data1, *d_data3;
  int *d_result;
  cudaMalloc(&d_data1, h_data1.size() * sizeof(int));
  cudaMalloc(&d_data3, h_data3.size() * sizeof(int));
  cudaMalloc(&d_btime1, h_btime1.size() * sizeof(int));
  cudaMalloc(&d_btime3, h_btime3.size() * sizeof(int));
  cudaMalloc(&d_time1, h_time1.size() * sizeof(int));
  cudaMalloc(&d_time3, h_time3.size() * sizeof(int));
  cudaMalloc(&d_result, 6 * sizeof(int));

  cudaMemcpy(d_data1, h_data1.data(), h_data1.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_data3, h_data3.data(), h_data3.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_btime1, h_btime1.data(), h_btime1.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_btime3, h_btime3.data(), h_btime3.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_time1, h_time1.data(), h_time1.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_time3, h_time3.data(), h_time3.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemset(d_result, 0, 6 * sizeof(int) );

  dim3 threadsPerBlock(256);
  dim3 numBlocks((numSteps + threadsPerBlock.x -1) / threadsPerBlock.x);

  computeCommittorsKernel<<<numBlocks, threadsPerBlock>>>( d_data1, d_time1, d_btime1, d_data3, d_time3, d_btime3, numSteps, numTraj, d_result );

  std::vector<int> counts(6, 0);
  cudaMemcpy( counts.data(), d_result, 6*sizeof(int), cudaMemcpyDeviceToHost);
  for ( int i=0; i<5; i++) {
    result[i] = (float) counts[i+1] / counts[0];
  }

  cudaFree(d_data1);
  cudaFree(d_data3);
  cudaFree(d_time1);
  cudaFree(d_time3);
  cudaFree(d_btime1);
  cudaFree(d_btime3);
  cudaFree(d_result);
}

__global__ void passageTimeKernel(float* dist,  int* data, int* time, int* btime, float rcut, size_t numTraj, size_t numSteps) {
  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < numTraj) {
    size_t index, lastindex;
    if ((dist[idx] > rcut)  ) {
      data[idx] = 1;
    } else {
      data[idx]=0;
    }
    time[idx]=0;
    btime[idx]=-1;
    lastindex=0;
    for (size_t i=1; i<numSteps; i++){
      index = i*numTraj+idx;
      if ( ( (dist[index - numTraj] > rcut)  && (  dist[index] <= rcut  ) ) ||
           ( (dist[index] > rcut ) && (  dist[index-numTraj] <= rcut   ))){ 
        data[index] = data[index-numTraj] + 1;
        for (size_t k=lastindex; k< i; k++) {
          time[ k*numTraj + idx] = i;
        }
        for (size_t k=lastindex+1; k< i; k++) {
          btime[ k*numTraj + idx] = btime[ lastindex*numTraj + idx];
        }
        btime[i*numTraj + idx]=i;
        lastindex=i;
      } else{
        data[index] = data[index-numTraj];
      }
    }
    for( size_t k=lastindex; k< numSteps; k++) {
      time[ k*numTraj + idx] = -1;
      btime[ k*numTraj + idx] = btime[lastindex*numTraj + idx];
    }
  }
}


void findPassageTime(float rcut, size_t numTraj, size_t numSteps, std::vector<float>& fdist, std::vector<int>& data, std::vector<int>& time, std::vector<int>& btime){
  size_t dataSize = data.size();
  size_t dataBytes = dataSize * sizeof(int);
  size_t distBytes = dataSize * sizeof(float);

  int *d_data, *d_time, *d_btime;
  float *d_dist;
  cudaMalloc(&d_data, dataBytes);
  cudaMalloc(&d_time, dataBytes);
  cudaMalloc(&d_btime, dataBytes);
  cudaMalloc(&d_dist, distBytes);

  cudaMemset(d_data, 0, dataBytes);
  cudaMemset(d_time, 0, dataBytes);
  cudaMemset(d_btime, 0, dataBytes);
  cudaMemcpy(d_dist, fdist.data(), distBytes, cudaMemcpyHostToDevice);
  
  int threadsPerBlock = 8;
  int blocksPerGrid = (numTraj + threadsPerBlock-1) / threadsPerBlock;
  passageTimeKernel<<<blocksPerGrid, threadsPerBlock>>>(d_dist, d_data, d_time, d_btime, rcut, numTraj, numSteps);

  cudaMemcpy(data.data(), d_data, dataBytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(time.data(), d_time, dataBytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(btime.data(), d_btime, dataBytes, cudaMemcpyDeviceToHost);
  cudaFree(d_data);
  cudaFree(d_time);
  cudaFree(d_btime);
  cudaFree(d_dist);
}


std::vector<std::vector<float>> readBTrajFromFile(const std::string& filename) {
  std::ifstream in(filename, std::ios::binary);
  int rows, cols;
  in.read(reinterpret_cast<char*>(&rows), sizeof(int));
  in.read(reinterpret_cast<char*>(&cols), sizeof(int));

  std::vector<std::vector<float>> matrix(rows, std::vector<float>(cols));
  for( int i=0; i<rows; i++){
    in.read(reinterpret_cast<char*>(matrix[i].data()), sizeof(float)*cols);
   }
  in.close();
//  std::cout << "file reading done" << std::endl;
  return matrix;
}

std::vector<float> linspace(float start, float stop, int num, bool endpoint = true) {
  std::vector<float> result;
  if (num <= 0) {
    return result;
  }

  float step = (stop-start) / (endpoint ? num-1 : num);

  for (int i=0; i< num ; i++){
    float exponent = start + i *step;
    result.push_back(exponent);
  }
  return result;
}



void generateArgumentVector(std::vector<float>& argVector, int numBins, float minVal, float maxVal) {
  float binWidth = (maxVal - minVal) / numBins;
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
void printFile(std::string& file, std::vector<float>& x, std::vector<float>& y){
  std::cout << "Outpul File generated: " << file << std::endl;
  std::cout << "Size: " << x.size() <<  std::endl;

  std::ofstream out ( file);
  for (size_t i = 0; i< x.size(); i++) {
    out << std::fixed << std::setprecision(12) << x[i] << " " << y[i]  << std::endl;
  }
  out.close();
}


void readTrajFromFile(std::string& file,
                      std::vector<std::vector<float>>& traj) {
  std::ifstream in(file);

  if (!in.is_open()) {
    std::cerr << "Error: Cannot find file" << std::endl;
  }

  std::string line;
  while (std::getline(in, line)) {
    std::istringstream lineStream(line);
    std::vector<float> row;

    float value;

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

std::vector<float> logspace(float start, float stop, int num, bool endpoint = true) {
  std::vector<float> result;
  if (num <= 0) {
    return result;
  }

  float step = (stop-start) / (endpoint ? num-1 : num);

  for (int i=0; i< num ; i++){
    float exponent = start + i *step;
    result.push_back(pow(10, exponent));
  }
  return result;
}

template <typename T>
std::vector<T> removeDuplicates(const std::vector<T>& input) {
  std::set<T> uniqueSet(input.begin(), input.end());
  return std::vector<T>(uniqueSet.begin(), uniqueSet.end());
}

float mean(const std::vector<float>& data) {
  float sum = 0.0;
  for (float x : data) sum += x;
  return sum/data.size();
}

float variance ( const std::vector<float>& data) {
  float mu = mean(data);
  float var = 0.0;
  for ( float x : data) var += (x-mu) * (x-mu) ;
  return var/( data.size() -1);
}

float blockAverage( const std::vector<float>& data) {
  std::vector<float> blockData = data;
  std::vector<float> semList;
  std::vector<int> blockSizes;

  int blockSize = 1;
  int level=0;

  while ( blockData.size() >= 2) {
    int nBlocks = blockData.size() / 2;
    // new block
    std::vector<float> newBlockData(nBlocks);
    for ( int i =0 ; i<nBlocks; i++ ) {
      newBlockData[i] = 0.5 * ( blockData[2*i] + blockData[2*i+1] );
    }
    float var = variance(newBlockData);
    float sem = std::sqrt(var / (nBlocks)); 

    semList.push_back(sem);
    blockSizes.push_back(blockSize);
  //  std::cout << blockSize << "\t\t" << nBlocks << "\t\t" << std::setprecision(6) << sem << "\n";

    blockData = newBlockData;
    blockSize *= 2;
    ++level;
  }
  int stabilityWindow=3;
  float tolerance=0.05;
  float finalSEM=0;
  for ( size_t i =0; i+stabilityWindow <= semList.size(); ++i){
    bool stable = true;
    float ref = semList[i];
    for ( size_t j=0; j<stabilityWindow; ++j){
      float diff = std::fabs(semList[i+j] - ref);
      if ( diff/ref > tolerance) {
        stable = false;
        break;
      }
    }
    if (stable){
//      std::cout << "\nSEM plateau detected starting at block step " <<i << " (stable for " << stabilityWindow << " sizes)\n";
      finalSEM=semList[i];
      break;
    }
  }
  if (finalSEM==0) {
    //std::cout << semList[0] << "\t\t" << semList[1];
    std::cout << "\nERROR. SEM plateua not detected. Use the largest error\n";
    finalSEM=*std::max_element(semList.begin(), semList.end());
    if (finalSEM > 1e6 ) {
      finalSEM = mean(data);
      std::cout << "\n\nZero Err occur, replace by mean\n\n";
    }
    //finalSEM = semList.back();
    //std::cout << "Estimated mean: " << finalSEM << "\t\tSEM: " << "\n";
//    std::cout << "Estimated mean: " << mean(data) << "\t\tmaxSEM: " << finalSEM << "\t\tlastSEM: " << semList.back() <<"\n";
  }
//  std::cout << "Estimated mean: " << mean(data) << "\t\tSEM: " << finalSEM << "\n";
  return finalSEM;
}


int main(int argc, char* argv[]) {
  if (argc != 6 ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << " density field cutoff_in(A) cutoff_out(A) timestep(ps) \n";
    return 1;
  }
  std::string density=argv[1];
  std::string field=argv[2];
  float CUTOFFin = std::stof(argv[3]);
  float CUTOFFout = std::stof(argv[4]);
  float timestep=std::stof(argv[5]); // 50 fs  = 0.05 ps 
  
  std::vector<std::vector<float>> dist;
  std::string distFile = "../data/cnnDist/distD" + density + "E" + field + ".binary";
  dist=readBTrajFromFile(distFile);
  std::vector<float> fdist;
  for (size_t i =0; i<dist.size(); i++) {
    fdist.insert(fdist.end(), dist[i].begin(), dist[i].end());
  }
  size_t numTraj = dist[0].size();
  for (auto& row : dist){
    row.clear();
    row.shrink_to_fit();
  }
  dist.clear();
  dist.shrink_to_fit();


  // tranform distance data into indicator data
  size_t numSample = fdist.size();
  size_t numSteps = numSample / numTraj;

  std::cout << "TPT-based Rate Measurement Starts...\n";
  std::cout << "Read Trajectory of length " << numSteps*timestep/1000 << "ns\n\n";


  // first passage time for inner cutoff domain
  float rcut;
  rcut = CUTOFFin;
  std::vector<int> h_data1(numSample, 0);  // indicator function  odd if r> rcut
  std::vector<int> h_time1(numSample, 0); // first entrance time
  std::vector<int> h_btime1(numSample, 0); // last exit time
  findPassageTime(rcut, numTraj, numSteps, fdist, h_data1, h_time1, h_btime1);

  std::cout << "finding first passage time at " << rcut << " done" << std::endl;


  rcut = CUTOFFout;
  std::vector<int> h_data3(numSample, 0); // indicator function  odd if r>rcut
  std::vector<int> h_time3(numSample, 0); // first entrance time 
  std::vector<int> h_btime3(numSample, 0); // last exit time from
  findPassageTime(rcut, numTraj, numSteps, fdist, h_data3, h_time3, h_btime3);
  std::cout << "finding first passage time at " << rcut << " done" << std::endl;

  std::vector<int> h_state(numSample, 0); // 1 if B(CIP), 2 if intermediate, 3 if A (FREE)
  for ( size_t element=0; element < numSample; element++) {
    if( h_data1[element] % 2 == 0) {
      h_state[element] = 1;
    }else if (h_data3[element] % 2 == 1) {
      h_state[element] = 3;
    } else{
      h_state[element] = 2;
    }
  }


  int prevIndex,currentIndex, count;
  count=0;
  std::vector<float> meanKr; // mean Kramers Rate over AB and BA
  std::vector<float> semsKr; // standard error of mean from Block Average
  for ( size_t traj=0; traj < numTraj; traj++) {
    std::vector<float> nuAB;  //indicator function, 1 only if at t-dt~t, traj have "reactive trajectory" crossing B boundary
    std::vector<float> nuBA;  //indicator function, 1 only if at t-dt~t, traj have "reactive trajectory" crossing B boundary
    for ( size_t time=1; time < numSteps; time++ ) {
      currentIndex = time * numTraj + traj;
      //construction of nuAB: only when all btimes and times are defined
      if( h_time1[currentIndex] != -1 && h_time3[currentIndex] != -1 
          && h_btime1[currentIndex] != -1 && h_btime3[currentIndex] != -1 ) {
        prevIndex = (time-1)*numTraj + traj;
        // find transitions A->B
        if( h_state[prevIndex] == 2 && h_state[currentIndex] == 1 
            && h_btime3[prevIndex] > h_btime1[prevIndex] ) {
          nuAB.push_back(1);
          count+=1;
        }else {
          nuAB.push_back(0);
        }
        // find transitions B->A
        if( h_state[prevIndex] == 2 && h_state[currentIndex] == 3 
            && h_btime3[prevIndex] < h_btime1[prevIndex] ) {
          nuBA.push_back(1);
          count+=1;
        }else {
          nuBA.push_back(0);
        }
      }
    }
    meanKr.push_back(mean(nuAB));
    meanKr.push_back(mean(nuBA));
    semsKr.push_back(blockAverage(nuAB));
    semsKr.push_back(blockAverage(nuBA));
//    std::cout << "Traj " << traj << " AB " << mean(nuAB)*nuAB.size() << " BA " << mean(nuBA)* nuAB.size()<< "\n";
  }

  float KramersRate=0;
  float KramersErr=0;
  for( int i =0; i<meanKr.size(); i++ ) {
    KramersRate += meanKr[i]/semsKr[i]/semsKr[i];
    KramersErr += 1/semsKr[i]/semsKr[i];
  }
  KramersRate /= KramersErr*timestep/1000;
  KramersErr = std::sqrt(1/KramersErr) / timestep*1000;
  std::cout << "Kramers Rate : " << KramersRate << " /ns\t\t Err: " << KramersErr  << " count: " << count << "\n";


  // check if A<->B jump occurs 
  count=0;
  for (int traj=0; traj < numTraj; traj++) {
    for ( size_t time =1; time < numSteps-1 ; time++) {
      currentIndex  = time * numTraj + traj;
      prevIndex = currentIndex-numTraj;
      if (std::abs(h_state[currentIndex]-h_state[prevIndex]) == 2) {
        count+=1;
        std::cout << "Error! jump occured,  count : " << count << std::endl;
      }
    }
  }

  std::vector<std::vector<float>> committor(5); // q+, q-, pA, pAB, pB
  std::vector<std::vector<float>> committorErr(5); // q+, q-, pA, pAB, pB

  for ( size_t traj=0; traj < numTraj; traj++) {
    std::vector<std::vector<float>> commit(5);  //indicator function, 1 
    for ( size_t time=0; time < numSteps; time++ ) {
      currentIndex = time * numTraj + traj;
      //consideer only part of trajectory who has full history
      if( h_time1[currentIndex] != -1 && h_time3[currentIndex] != -1 
          && h_btime1[currentIndex] != -1 && h_btime3[currentIndex] != -1 ) {
        // find q+: go to B before A
        if( h_time1[currentIndex] < h_time3[currentIndex] ) { 
          commit[0].push_back(1);
        }else {
          commit[0].push_back(0);
        }
        // find q-: come from A than B
        if( h_btime1[currentIndex] < h_btime3[currentIndex] ) { 
          commit[1].push_back(1);
        }else {
          commit[1].push_back(0);
        }
        // population A
        if ( h_state[currentIndex] == 3) {
          commit[2].push_back(1);
        } else {
          commit[2].push_back(0);
        }
        // population AB
        if ( h_state[currentIndex] == 2) {
          commit[3].push_back(1);
        } else {
          commit[3].push_back(0);
        }
        // population B
        if ( h_state[currentIndex] == 1) {
          commit[4].push_back(1);
        } else {
          commit[4].push_back(0);
        }

      }
    }
    for( int j=0; j<committor.size(); j++ ) {
      committor[j].push_back( mean(commit[j]));
      committorErr[j].push_back( blockAverage(commit[j]));
    }
//    std::cout << "Traj " << traj << " AB " << mean(nuAB)*nuAB.size() << " BA " << mean(nuBA)* nuAB.size()<< "\n";
  }
  std::vector<float> fcommittor(5, 0.0);
  std::vector<float> fcommittorErr(5, 0.0);
  for( size_t i =0; i<committor.size(); i++ ) {
    for( size_t j=0; j<committor[i].size(); j++ ) {
      fcommittor[i] += committor[i][j]/committorErr[i][j]/committorErr[i][j];
      //fcommittor[i] += committor[i][j];
      fcommittorErr[i] += 1.0/committorErr[i][j]/committorErr[i][j];
    }
    fcommittor[i] /= fcommittorErr[i];
    //fcommittor[i] /= committor[i].size();
    fcommittorErr[i] =std::sqrt(1/fcommittorErr[i]) ;
  }

  std::cout << "\n\nCommittors and Population Analysis\n";
  float sumq, sumqerr, sump, sumperr;
  sumq = fcommittor[0] + fcommittor[1];
  sumqerr = std::sqrt( fcommittorErr[0]*fcommittorErr[0] + fcommittorErr[1]*fcommittorErr[1])/2.;
  sump = fcommittor[2] + fcommittor[3] + fcommittor[4];
  sumperr = std::sqrt( fcommittorErr[2]*fcommittorErr[2] + fcommittorErr[3]*fcommittorErr[3] + fcommittorErr[4]*fcommittorErr[4])/3.;

  std::cout << "means:\t\t" << std::fixed << std::setprecision(7) << fcommittor[0] << "\t\t" << fcommittor[1] << "\t\t" << fcommittor[2] << "\t\t"
    << fcommittor[3] << "\t\t" << fcommittor[4] << "\t\t" << sumq << "\t\t" << sump  <<  "\n";

  std::cout << "SEMS:\t\t" << std::fixed << std::setprecision(7) << fcommittorErr[0] << "\t\t" << fcommittorErr[1] << "\t\t" << fcommittorErr[2] << "\t\t"
    << fcommittorErr[3] << "\t\t" << fcommittorErr[4] << "\t\t" <<  sumqerr << "\t\t" << sumperr <<  "\n";


  /*
  // compute committors and populations
  std::vector<float> committors(5, 0); // forward, backward, A, B, C
  if (1) {
    computeCommittors(h_data1, h_time1, h_btime1, h_data3, h_time3, h_btime3,  numSteps, numTraj, committors);
    std::cout << "forward : " << committors[0] << " backward " << committors[1] << " sum " << committors[0]+committors[1] << " pA "  << committors[2] << " pAB " << committors[3] << " pB " << committors[4] << std::endl;
  }
  */

  // compute recombination MFPT with conditioning, few sample > no block average
  int currentTime;
  count=0;
  std::vector<float> pathTime;
  for (int traj=0; traj < numTraj; traj++) {
    currentTime=0;
    currentIndex = currentTime * numTraj + traj;
    for ( int j =0; j<100000; j++) {
      if (currentIndex < numSteps*numTraj ) {
        currentTime = h_time3[currentIndex];
      } else {
        currentTime=0;
        break;
      }
      if (currentTime == -1){
        currentTime=0;
        break;
      }
      currentIndex = (currentTime) * numTraj + traj; // when a traj first hits A boundary  2->3 or 3->2
      prevIndex = (currentTime-1) * numTraj + traj; // right before it hits A boundary
      if ( h_btime1[prevIndex] > h_btime3[prevIndex] &&  h_time1[currentIndex] != -1 && prevIndex >0 ) { // if the trajectory comes from B, not A, and currently located at partial A
        pathTime.push_back(h_time1[currentIndex] - currentTime);
        count +=1;
      }
    }
  }
  float rateAB = 1/(mean(pathTime)*timestep/1000);
  float rateABErr = std::sqrt(variance(pathTime)/pathTime.size()) /(timestep/1000) / mean(pathTime) / mean(pathTime);

  float KramersABErr = KramersRate/fcommittor[1]*std::sqrt( (KramersErr/KramersRate)*(KramersErr/KramersRate) + (fcommittorErr[1]/fcommittor[1])*(fcommittorErr[1]/fcommittor[1]) );
  std::cout << "\nRecombination MFPT rate: " << rateAB << " /ns\tErr: " << rateABErr <<  "\tcount " << count << "\tcompare: " << KramersRate/fcommittor[1] << "\tErr\t" << KramersABErr << std::endl;


  // compute dissociation MFPT with conditioning
  count=0;
  std::vector<float> bpathTime;
  for (int traj=0; traj < numTraj; traj++) {
    currentTime=0;
    currentIndex = currentTime * numTraj + traj;
    for ( int j =0; j<100000; j++) {
      if (currentIndex < numSteps*numTraj ) {
        currentTime = h_time1[currentIndex];
      } else {
        currentTime=0;
        break;
      }
      if (currentTime == -1){
        currentTime=0;
        break;
      }
      currentIndex = (currentTime) * numTraj + traj; // when a traj first hits A boundary  2->3 or 3->2
      prevIndex = (currentTime-1) * numTraj + traj; // right before it hits A boundary
      if ( h_btime1[prevIndex] < h_btime3[prevIndex] &&  h_time3[currentIndex] != -1 ) { // if the trajectory comes from B, not A, and currently located at partial A
        bpathTime.push_back(h_time3[currentIndex] - currentTime);
        count +=1;
      }
    }
  }
  float rateBA = 1/(mean(bpathTime)*timestep/1000);
  float rateBAErr = std::sqrt(variance(bpathTime)/bpathTime.size()) /(timestep/1000) / mean(bpathTime) / mean(bpathTime);

  float KramersBAErr = KramersRate/fcommittor[0]*std::sqrt( (KramersErr/KramersRate)*(KramersErr/KramersRate) + (fcommittorErr[0]/fcommittor[0])*(fcommittorErr[0]/fcommittor[0]) );
  std::cout << "\nDissociation MFPT rate: " << rateBA << " /ns\tErr: " << rateBAErr <<  "\tcount " << count << "\tcompare: " << KramersRate/fcommittor[0] << "\tErr\t" << KramersBAErr << std::endl;

  
  float rho;
  if (density == "05") {
    rho = 0.5;
  } else if (density == "025") {
    rho = 0.25;
  } else {
    rho = 0;
  }
  float bimErr = rateAB/fcommittor[1]/rho * std::sqrt( (rateABErr/rateAB)*(rateABErr/rateAB) + (fcommittorErr[1]/fcommittor[1])* (fcommittorErr[1]/fcommittor[1]) );
  float KassErr = fcommittor[0]/fcommittor[1]/fcommittor[1]/rho * std::sqrt( (fcommittorErr[0]/fcommittor[0] ) *(fcommittorErr[0]/fcommittor[0]) + (fcommittorErr[1]/fcommittor[1])*(fcommittorErr[1]/fcommittor[1]) *2 );
  // generate output result
  std::string resfile = "../results/rate/rateD" + density + "E" + field + ".dat";
  std::ofstream out (resfile );
  out << std::fixed << std::setprecision(5) << density << "\t\t" << field << "\t\t"  
    << CUTOFFin << "\t\t" << CUTOFFout <<"\t\t" << numSteps*timestep/1000  << "\t\t"
    << fcommittor[0] << "\t\t" << fcommittorErr[0] << "\t\t" << fcommittor[1] << "\t\t" << fcommittorErr[1] << "\t\t" 
    << KramersRate << "\t\t" << KramersErr << "\t\t" 
    << rateAB << "\t\t" << rateABErr << "\t\t" << rateBA << "\t\t" << rateBAErr << "\t\t" 
    << rateAB/fcommittor[1]/rho << "\t\t" << bimErr  << "\t\t"
    << fcommittor[0]/fcommittor[1]/fcommittor[1]/rho << "\t\t" << KassErr  << "\t\t"
    << "\n";
  out.close();
  return 0;
}

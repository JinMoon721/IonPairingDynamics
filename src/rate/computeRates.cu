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
      if ( time1[currentIndex] != -1 && btime3[currentIndex] != -1 ){ 
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

__global__ void passageTimeKernelA(float* dist,  int* data, int* time, int* btime, float rcut, size_t numTraj, size_t numSteps) {
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

__global__ void passageTimeKernelB(float* dist, int* data, int* time, int* btime, float rcut, size_t numTraj, size_t numSteps) {
  size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < numTraj) {
    size_t index, lastindex;
    if ((dist[idx] <= rcut) ) {
      data[idx] = 0;
    } else {
      data[idx]=1;
    }
    time[idx]=0;
    btime[idx]=-1;
    lastindex=0;
    for (size_t i=1; i<numSteps; i++){
      index = i*numTraj+idx;
      if ( ( (dist[index - numTraj] > rcut)  && (  dist[index] <= rcut ) ) ||
           ( (dist[index] > rcut ) && (  dist[index-numTraj] <= rcut  ))){ 
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
  if ( rcut > 20) {
    passageTimeKernelA<<<blocksPerGrid, threadsPerBlock>>>(d_dist, d_data, d_time, d_btime, rcut, numTraj, numSteps);
  } else{
    passageTimeKernelB<<<blocksPerGrid, threadsPerBlock>>>(d_dist, d_data, d_time, d_btime, rcut, numTraj, numSteps);
  }

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
  std::cout << "file reading done" << std::endl;
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

__global__ void computeAbsCorrKernel(const int* data1, const int* time1, const int* btime1,
                                       const int* data3, const int* time3, const int* btime3,
                                       const float* dist, const float* pdist, int* nresult, int* dresult, size_t numSteps, size_t numTraj, size_t numLag, const int* dt){
  size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < numSteps) {
    size_t lag,  initIndex, lastIndex, num, den;
    int check, currentTime, nextTime;
    size_t currentIndex;
    for (size_t lagind = 0 ; lagind < numLag; lagind++) {
      lag = dt[lagind];
      if (idx < numSteps - lag) {
        num = 0;
        den = 0;
        for (size_t traj = 0 ; traj < numTraj; traj++) {
          currentTime=idx;
          currentIndex = currentTime*numTraj + traj;
          //if (( btime3[currentIndex] < btime1[currentIndex] ) && (btime3[currentIndex] != -1 ) && (dist[currentIndex] >= 65.0) && (data3[currentIndex] % 2 == 0)) { // conditioned on 1-q- near the boundary
          if ( ( btime3[currentIndex] < btime1[currentIndex] ) && (btime3[currentIndex] != -1 ) && 
              ( ((dist[currentIndex] >= 65.0) && (dist[currentIndex] <= 70.0)) || ( pdist[currentIndex] >= 25.0 && pdist[currentIndex] <= 35.0 && dist[currentIndex] > 70.0)  ) )  { // conditioned on 1-q- near the boundary
            den +=1;
            if (time1[currentIndex] < lag+idx) {
              num +=1;
            }
          }
        }
        atomicAdd(&dresult[lagind], den);
        atomicAdd(&nresult[lagind], num);
      }
    }
  }
}


void computeAbsCorr(const std::vector<int>& dt, const std::vector<int>& h_data1, const std::vector<int>& h_time1,  const std::vector<int>& h_btime1,
                                                  const std::vector<int>& h_data3, const std::vector<int>& h_time3,  const std::vector<int>& h_btime3,
                                                  const std::vector<float>& h_dist, const std::vector<float>& h_pdist, size_t numLag, size_t numSteps, size_t numTraj, 
                                                  std::vector<float>& corr) {
  // allocate GPU memory and initialize
  int *d_dresult, *d_nresult;
  int *d_data1, *d_data3, *d_dt, *d_time1, *d_btime1, *d_time3, *d_btime3;
  float *d_dist, *d_pdist;
  cudaMalloc(&d_data1, h_data1.size() * sizeof(int));
  cudaMalloc(&d_time1, h_time1.size() * sizeof(int));
  cudaMalloc(&d_btime1, h_btime1.size() * sizeof(int));
  cudaMalloc(&d_data3, h_data3.size() * sizeof(int));
  cudaMalloc(&d_time3, h_time3.size() * sizeof(int));
  cudaMalloc(&d_btime3, h_btime3.size() * sizeof(int));
  cudaMalloc(&d_dist, h_dist.size() * sizeof(float));
  cudaMalloc(&d_pdist, h_pdist.size() * sizeof(float));
  cudaMalloc(&d_dt, dt.size() * sizeof(int));
  cudaMalloc(&d_dresult, numLag * sizeof(int));
  cudaMalloc(&d_nresult, numLag * sizeof(int));

  cudaMemset(d_dresult, 0, numLag * sizeof(int));
  cudaMemset(d_nresult, 0, numLag * sizeof(int));
  cudaMemcpy(d_data1, h_data1.data(), h_data1.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_time1, h_time1.data(), h_time1.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_btime1, h_btime1.data(), h_btime1.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_data3, h_data3.data(), h_data3.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_time3, h_time3.data(), h_time3.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_btime3, h_btime3.data(), h_btime3.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_dist, h_dist.data(), h_dist.size() * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_pdist, h_pdist.data(), h_pdist.size() * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_dt, dt.data(), dt.size() * sizeof(int), cudaMemcpyHostToDevice);

  dim3 threadsPerBlock(256);
  dim3 numBlocks((numSteps + threadsPerBlock.x -1) / threadsPerBlock.x);

  computeAbsCorrKernel<<<numBlocks, threadsPerBlock>>>( d_data1, d_time1, d_btime1, d_data3, d_time3, d_btime3, d_dist, d_pdist, d_nresult, d_dresult, numSteps, numTraj, numLag, d_dt);

  std::vector<int> numerator(numLag, 0);
  std::vector<int> denominator(numLag, 0);

  cudaMemcpy( numerator.data(), d_nresult, numLag*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy( denominator.data(), d_dresult, numLag*sizeof(int), cudaMemcpyDeviceToHost);

  std::transform(numerator.begin(), numerator.end(), denominator.begin(), corr.begin(), 
      [](float a, float b) {
        return b != 0.0 ? a/b : 0.0;
        });

  cudaFree(d_data1);
  cudaFree(d_time1);
  cudaFree(d_btime1);
  cudaFree(d_data3);
  cudaFree(d_time3);
  cudaFree(d_btime3);
  cudaFree(d_dist);
  cudaFree(d_pdist);
  cudaFree(d_dresult);
  cudaFree(d_nresult);
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



int main(int argc, char* argv[]) {
  if (argc != 7 ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << " density field cutoff_in(A) cutoff_out(A) box(A) timestep(ps) \n";
    return 1;
  }
  int numtraj=1;
  int numsnap;
  std::string density=argv[1];
  std::string field=argv[2];
  std::map<std::string, float> str2float = {
    {"00", 0.0},
    {"10", 0.23},
    {"50", 1.15},
    {"100", 2.3} 
  };
  float efield = str2float[field]; //kcal/molA

  float CUTOFFin = std::stof(argv[3]);
  float CUTOFFout = std::stof(argv[4]);
  float box=std::stof(argv[5]);
  float timestep=std::stof(argv[6]); // 50 fs  = 0.05 ps 
  
  int numatoms;
  if ( density == "05" ) {
    numatoms = 16;
  } else {
    numatoms=32;
  }

  std::vector<std::vector<float>> dist;
  std::string distFile = "../cnnTraj/trajD" + density + "E" + field + ".binary";
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

  std::cout << "Traj: " << numTraj << " Sample: " << numSample << " #snaps: " << numSteps << " Time : " << numSteps*timestep/1000 << "ns"  << std::endl;


  // first passage time for A: nearest neighbor anion is bare and > 70. B: dipole is also bare
  float rcut;
  rcut = CUTOFFin;
  std::vector<int> h_data1(numSample, 0);  // indicator function  odd in A state, even in A^C
  std::vector<int> h_time1(numSample, 0); // first entrance time
  std::vector<int> h_btime1(numSample, 0); // last exit time
  findPassageTime(rcut, numTraj, numSteps, fdist, h_data1, h_time1, h_btime1);

  std::cout << "finding first passage time at " << rcut << " done" << std::endl;


  rcut = CUTOFFout;
  std::vector<int> h_data3(numSample, 0); // indicator function  even in B state, odd in B^C
  std::vector<int> h_time3(numSample, 0); // first entrance time 
  std::vector<int> h_btime3(numSample, 0); // last exit time from
  findPassageTime(rcut, numTraj, numSteps, fdist, h_data3, h_time3, h_btime3);
  std::cout << "finding first passage time at " << rcut << " done" << std::endl;

  // compute Kramers Rate
  int currentTime, currentIndex, ipassageTime, count;
  float KramersRateAB, KramersRateBA, KramersRate;
  if (1){
    count=0;
    std::vector<int> passageTime;
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
        if ( (h_time1[currentIndex] < h_time3[currentIndex] )  && (h_time1[currentIndex] != -1) ) { // reactive trajectory not including A->B transition
          ipassageTime = h_time1[currentIndex] - currentTime;
          passageTime.push_back(ipassageTime);
          count+=1;
        }
      }
    }
    KramersRateAB = (float) count / numSteps /timestep/  numTraj * 1000;
    std::cout << "Kramers Rate AB : " << KramersRateAB << " /ns " << "count " << count<< std::endl;
  }

  if (1){ //  Kramers Rate BA
    count=0;
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

        currentIndex = (currentTime) * numTraj + traj; 
        if ( (h_time1[currentIndex] > h_time3[currentIndex] )  && (h_time3[currentIndex] != -1)) { 
          ipassageTime = h_time3[currentIndex] - currentTime;
          count+=1;
        }
      }
    }
    KramersRateBA = (float) count / numSteps/ timestep / numTraj * 1000;
    std::cout << "Kramers Rate BA : " << KramersRateBA << " /ns " << "count " << count<< std::endl;
  }
  KramersRate = (KramersRateAB + KramersRateBA) / 2.0;

  // find any trajectory jump from A to B
  if (1){
    count=1;
    //for (int traj=0; traj < numTraj; traj++) {
    for (int traj=0; traj < 1; traj++) {
      for ( size_t time =0; time < numSteps-1 ; time++) {
        currentIndex  = time * numTraj + traj;
        if ( h_data3[currentIndex] % 2 == 1 && h_data1[currentIndex + numTraj] % 2 == 0 ) { // transition from 3 to 1
          std::cout << "Error! jump from A to B occured,  count : " << count << std::endl;
          count+=1;
        }
        //if (fdist[currentIndex] > 300) {
          //std::cout << "Error! dist,  count : " << fdist[currentIndex] << std::endl;
        //}
      }
    }
  }

  // compute committors and populations
  std::vector<float> committors(5, 0); // forward, backward, A, B, C
  if (1) {
    computeCommittors(h_data1, h_time1, h_btime1, h_data3, h_time3, h_btime3,  numSteps, numTraj, committors);
    std::cout << "forward : " << committors[0] << " backward " << committors[1] << " sum " << committors[0]+committors[1] << " pA "  << committors[2] << " pAB " << committors[3] << " pB " << committors[4] << std::endl;
  }

  // compute MFPT with conditioning
  if (1){
    float MFPT;
    int pathTime=0;
    int prevIndex=0;
    count=0;
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
        if ( h_btime1[prevIndex] > h_btime3[prevIndex] &&  h_time1[currentIndex] != -1 ) { // if the trajectory comes from B, not A, and currently located at partial A
          pathTime += h_time1[currentIndex] - currentTime;
          count +=1;
        }
      }
    }
    MFPT = (float) pathTime * timestep / count / 1000;
    std::cout << "MFPT is : " << MFPT << " ns " << count << std::endl;
    std::cout << "compare it with <q->/Kr  : " << committors[1] / KramersRate  << "  " << std::endl;
  }

  // compute MFPT for dissociation
  if (1){
    float MFPT;
    int pathTime=0;
    int prevIndex=0;
    count=0;
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
        currentIndex = (currentTime) * numTraj + traj; // when a traj first hits B boundary  1->2 or 2->1
        prevIndex = (currentTime-1) * numTraj + traj; // right before it hits A boundary
        if ( h_btime1[prevIndex] < h_btime3[prevIndex] && h_time3[currentIndex] != -1  ) { // if the trajectory comes from A, not B
          pathTime += h_time3[currentIndex] - currentTime;
          count +=1;
        }
      }
    }
    MFPT = (float) pathTime * timestep / count / 1000;
    std::cout << "Dissociation MFPT is : " << MFPT << " ns " << count << std::endl;
    std::cout << "compare it with <q+>/Kr  : " << committors[0] / KramersRate  << "  " << std::endl;
  }

  // compute probability of finding vacant ion
  if (1){
    int num=0;
    int den=0;
    for (int traj=0; traj < numTraj; traj++) {
      for (size_t time=0; time < numSteps; time++) {
        den+=1;
        currentIndex = time*numTraj + traj;
        if( fdist[currentIndex] > 5.0 ){
          num +=1;
        }
      }
    }
    std::cout << "vacant ion probability : " << (float) num/den << std::endl;
  }

  /// compute radial histograms

  return 0;
}

  




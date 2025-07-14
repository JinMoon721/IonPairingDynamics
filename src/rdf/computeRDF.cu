#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <iomanip>
#include <map>

void writeHistogram(std::string& file, std::vector<float>& arg, std::vector<int>& hist){
  std::cout << "Histogram File generated: " << file << std::endl;
  std::cout << "Size: " << arg.size() <<  std::endl;
  std::cout << "Arg | log10(g(r)) | log10(r^2 g(r)) " << std::endl;

  std::ofstream out ( file);
  float sum1=0;
  float sum2=0;
  float r=0;
  float dr=0;
  for (size_t i = 0; i< arg.size(); i++) {
    r = arg[i]; // assume uniform bin
    dr = arg[1]-arg[0];
    sum1 += static_cast<float>(hist[i]) * dr / r / r;
    sum2 += static_cast<float>(hist[i]) * dr ;
  }

  float val1, val2;
  for (size_t i = 0; i< arg.size(); i++) {
    r = arg[i] + (arg[1] - arg[0])/2.0 ; // assume uniform bin
    val1 = static_cast<float>(hist[i])/sum1/r/r;
    val2 = static_cast<float>(hist[i])/sum2;
    //out << std::fixed << std::setprecision(5) << arg[i] << "\t" << log10(val1)  << "\t" << log10(val2) << std::endl;
    out << std::fixed << std::setprecision(5) << arg[i] << "\t\t" << -std::log(val1)  << "\t\t" << -std::log(val2) << std::endl;
  }
  out.close();
}

__global__ void computeHistogramKernel(float* data, int* histogram, int numBins, float minVal, float maxVal, int dataSize) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < dataSize) {
    float value = data[idx];
    if (value >= minVal && value < maxVal) {
      int bin = (int) ((value-minVal) / (maxVal - minVal) * numBins);
      atomicAdd(&histogram[bin], 1);
    }
  }
}



void computeHistogramCUDA(const std::vector<float>& data, int numBins, float minVal, float maxVal, std::vector<int>& histogram) {
  int dataSize = data.size();
  size_t dataBytes = dataSize * sizeof(float);
  size_t histBytes = numBins * sizeof(float);

  float* d_data;
  int* d_histogram;
  cudaMalloc(&d_data, dataBytes);
  cudaMalloc(&d_histogram, histBytes);

  cudaMemset(d_histogram, 0, histBytes);
  cudaMemcpy(d_data, data.data(), dataBytes, cudaMemcpyHostToDevice);

  int threadsPerBlock = 256;
  int blocksPerGrid = (dataSize + threadsPerBlock-1) / threadsPerBlock;
  computeHistogramKernel<<<blocksPerGrid, threadsPerBlock>>>(d_data, d_histogram, numBins, minVal, maxVal, dataSize);

  histogram.resize(numBins);
  cudaMemcpy(histogram.data(), d_histogram, histBytes, cudaMemcpyDeviceToHost);
  cudaFree(d_data);
  cudaFree(d_histogram);
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

  std::string density=argv[1];
  std::string field=argv[2];

  std::vector<std::vector<float>> dist;
  std::string distFile = "../data/cnnDist/distD" + density + "E" + field + ".binary"; ; 
  dist=readBTrajFromFile(distFile);

  std::vector<float> fdist;
  for (size_t i =0; i<dist.size(); i++) {
    fdist.insert(fdist.end(), dist[i].begin(), dist[i].end());
  }
  dist.clear();
  dist.shrink_to_fit();


  // tranform distance data into indicator data
  size_t numSample = fdist.size();
  size_t numSteps = numSample ;

  std::cout << "Traj: " <<  " Sample: " << numSample << " Times: " << numSteps << std::endl;


  /// compute radial histograms
  int numBins=1000;
  float minVal = *std::min_element(fdist.begin(), fdist.end());
  float maxVal = *std::max_element(fdist.begin(), fdist.end());
  std::vector<int> histogram;
  std::vector<float> argVec;
  
  computeHistogramCUDA(fdist, numBins, minVal, maxVal, histogram);
  generateArgumentVector(argVec, numBins, minVal, maxVal);

  std::string rhist = "../results/hist/rdfD" + density + "E" +  field + ".dat";
  writeHistogram(rhist, argVec, histogram);

  return 0;
}

  




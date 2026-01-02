#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <iomanip>
#include <map>

#include <stats/stats.hpp>
using namespace stats;


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


std::vector<std::vector<double>> readBTrajFromFile(const std::string& filename) {
  std::ifstream in(filename, std::ios::binary);
  int32_t rows32=0, cols32=0;
  in.read(reinterpret_cast<char*>(&rows32), sizeof(rows32));
  in.read(reinterpret_cast<char*>(&cols32), sizeof(cols32));

  size_t rows = static_cast<size_t>(rows32);
  size_t cols = static_cast<size_t>(cols32);

  std::vector<std::vector<double>> out(rows, std::vector<double>(cols));
  std::vector<float> rowf(cols);
  for( size_t i=0; i<rows; i++){
    in.read(reinterpret_cast<char*>(rowf.data()), static_cast<std::streamsize>( sizeof(float)*cols));
    for (size_t j=0; j<cols; j++) {
      out[i][j] = static_cast<double>(rowf[j]);
    }
   }
  in.close();
  std::cout << "file reading done" << std::endl;
  return out;
}







float fieldConvert(const std::string& s) {
  if (s.find_first_not_of('0') == std::string::npos) {
    return 0.0f;
  }
  size_t leadingZeros = s.find_first_not_of('0');

  if ( leadingZeros >= 2) {
    std::string digits = s.substr(leadingZeros);
    float value = std::stof(digits);
    return value / std::pow(10, leadingZeros -1);
  }
  return std::stof(s);
}

int main(int argc, char* argv[]) {
  if (argc != 8 ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << " dir_name density field cutoff_in(A) cutoff_out(A) timestep(ps) thermoflag\n";
    return 1;
  }
  std::string dirName=argv[1]; // name of directory to output results, in results directory
  std::string rho = argv[2]; 
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
  std::string distFile;
  if (thermoflag==1){ distFile = std::string("../data/cnnDist/") + dirName + "TD" + rho + "E" + fieldname + ".binary";}
  else if (thermoflag==2){ distFile = std::string("../data/cnnDist/") + dirName + "TTD" + rho + "E" + fieldname + ".binary";}
  else { distFile = std::string("../data/cnnDist/") + dirName + "D" + rho + "E" + fieldname + ".binary";}
  auto dist=readBTrajFromFile(distFile);
  size_t numsnap = dist.size();
  size_t numatoms = dist[0].size();

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


  std::cout << "Radial Distribution Function computation Starts...\n";
  std::cout << "Read Trajectory of length " << numsnap*timestep/1000 << "ns\n\n";

  //flatten 2d double std vector into 1d vector
  std::vector<double> fdist;
  for (auto& row : dist) {
    fdist.insert(fdist.end(), row.begin(), row.end());
  }

  /// compute radial histograms
  size_t numBins=1000;
  auto rdf = makeHistogramAuto(fdist, numBins);

  std::string outfile = std::string("../results/rdf/") + dirName + "D" + rho + "E" +  fieldname + ".dat";
  std::ofstream out(outfile);
  writeHistogram(rdf, out);

  return 0;
}

  




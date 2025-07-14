#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include<algorithm>
#include <iomanip>
#include <cmath>
#include<cuda_runtime.h>
#include <utility>
#include <map>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <cstdlib> // for exit()

void writeHistogram(std::string& file, std::vector<float>& arg, std::vector<float>& hist, std::vector<int>& count){
  std::cout << "Histogram File generated: " << file << std::endl;
  std::cout << "Size: " << arg.size() <<  std::endl;
  std::cout << "Arg | log10(g(r)) | log10(r^2 g(r)) " << std::endl;
  float sum=0;
  for ( size_t i =0 ; i < count.size(); i++ ) {
    sum += count[i];
  }


  std::ofstream out ( file);
  for (size_t i = 0; i< arg.size(); i++) {
    //out << std::fixed << std::setprecision(5) << arg[i] << "\t" << log10(val1)  << "\t" << log10(val2) << std::endl;
    out << std::fixed << std::setprecision(5) << arg[i] << "\t" << hist[i]  << "\t" << count[i]/sum  << "\t" << hist[i] * count[i]/sum << std::endl;
  }
  out.close();
}


//using namespace uammd;
using CationTagMap = std::unordered_map<int, int>;
using InverseCationTagMap = std::unordered_map<int, int>;

void readMatrixFromFile(const std::string& file,
                        std::vector<std::vector<float>>& matrix,
                        std::vector<float>& time){
  std::ifstream in(file);

  if (!in.is_open()) {
    std::cerr << "Error: Cannot find file" << std::endl;
    //System::log<System::CRITICAL>("Cannot open file %s", file.c_str());
  }

  std::string line;
  while (std::getline(in, line)) {
    std::istringstream lineStream(line);
    std::vector<float> row;

    float value;

    if (lineStream >> value) {
      time.push_back(value);
    }

    while (lineStream >> value) {
      row.push_back(value);
    }
    
    if (!row.empty()){
      matrix.push_back(row);
    }
  }
  in.close();
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
}

void printMatrix(std::string& file,  std::vector<std::vector<float>>& matrix){
  std::cout << "Output File generated: " << file << std::endl;
  std::cout << "Rows: " << matrix.size() << ", Columns: " << (matrix.empty() ? 0 : matrix[0].size()) << std::endl;

  std::ofstream out ( file);

  for (size_t i = 0; i< matrix.size(); i++) {
    for (size_t j = 0; j< matrix[0].size(); j++) {
      out << std::fixed << std::setprecision(5) << matrix[i][j] << " ";
     }
    out << "\n";
  }
  out.close();
}

void printBinary(std::string& file,  std::vector<std::vector<float>>& matrix){
  std::cout << "Binary Output File generated: " << file << std::endl;
  std::cout << "Rows: " << matrix.size() << ", Columns: " << (matrix.empty() ? 0 : matrix[0].size()) << std::endl;

  std::ofstream out ( file, std::ios::binary);
  int rows = matrix.size();
  int cols = matrix[0].size();
  out.write(reinterpret_cast<const char*>(&rows), sizeof(int));
  out.write(reinterpret_cast<const char*>(&cols), sizeof(int));

  for (const auto& row : matrix) {
    out.write(reinterpret_cast<const char*>(row.data()), sizeof(float) * cols);
  }

  out.close();
}

std::vector<std::vector<float>> readBinary(const std::string& filename) {
  std::ifstream in(filename, std::ios::binary);
  int rows, cols;
  in.read(reinterpret_cast<char*>(&rows), sizeof(int));
  in.read(reinterpret_cast<char*>(&cols), sizeof(int));

  int eqsteps=0;
  std::vector<float> temp(cols);
  for( int i=0; i<eqsteps; i++ ) {
    in.read(reinterpret_cast<char*>(temp.data()), sizeof(float));
  }

  std::vector<std::vector<float>> matrix(rows-eqsteps, std::vector<float>(cols));
  for( int i=0; i<rows-eqsteps; i++){
    in.read(reinterpret_cast<char*>(matrix[i].data()), sizeof(float)*cols);
   }
  in.close();
  return matrix;
}

std::map<int, std::pair<int, int>> indexTopair(int size) {
  std::map<int, std::pair<int, int>> index2pair;
  int k=0;
  for ( int i =0; i<size-1; i++ ) {
    for ( int j =i+1; j<size; j++ ) {
      index2pair[k++] = {i, j};
    }
  }
  return index2pair;
}

float applyPBC(float x, float box){
  float hbox = box/2.0;
  float wrapped = fmod(x + hbox, box);
  if (wrapped < 0) wrapped += box;
  return wrapped - hbox;
}

float distance(const std::vector<float> a, const std::vector<float> b, float box) {
  float dx, dy, dz, rsq;
  dx=applyPBC( a[0] - b[0], box);
  dy=applyPBC( a[1] - b[1], box);
  dz=applyPBC( a[2] - b[2], box*2);
  rsq = dx*dx +dy*dy +dz*dz;
  return std::sqrt(rsq);
}

void findClusters(const std::vector<std::vector<int>>& graph, std::vector<std::vector<int>>& clusters, int N) {
  std::vector<bool> visited(N, false);
  
  for (int i = 0; i < N; ++i) {
    if (!visited[i]) {
      std::vector<int> cluster;
      std::queue<int> q;
      q.push(i);
      visited[i] = true;
      
      while (!q.empty()) {
        int node = q.front(); q.pop();
        cluster.push_back(node);
        
        for (int neighbor : graph[node]) {
          if (!visited[neighbor]) {
            visited[neighbor] = true;
            q.push(neighbor);
          }
        }
      }
      
      clusters.push_back(cluster);
    }
  }
}

float netCharge(const std::vector<int>& cluster) {
  float net=0;
  for( int i =0; i<cluster.size(); i++) {
    if ( cluster[i] % 2 == 0 ) {
      net +=1.0;
    } else {
      net -= 1.0;
    }
  }
  return net;
}

float findNearestAnionOutside (const std::vector<std::vector<float>>& ions, const std::vector<std::vector<int>>& clusters, std::vector<int>& atom2cluster, const std::vector<std::vector<float>>& distanceMatrix, int cation, int& nnAnion) {
  float minDist = 1e9;
  nnAnion=-1;
  
  for ( int anion=1; anion < ions.size(); anion+=2) {
    if ( netCharge(clusters[atom2cluster[anion]]) < 0) { //any clusters with total negative charge
      float d = distanceMatrix[cation][anion];
      if (d<minDist) {
        minDist = d;
        nnAnion=anion;
      }
    }
  }
  if (nnAnion == -1 ) {
    std::cout << "Error, finding nearest anion outside didn't work properply " << "\n";
    nnAnion=1;
    minDist=2.5;
  }
  return minDist;
}

float findNearestAnionInside (const std::vector<std::vector<float>>& ions, const std::vector<int>& cluster, std::vector<int>& atom2cluster, const std::vector<std::vector<float>>& distanceMatrix, int cation, int& nnAnion) {
  float minDist = 1e9;
  nnAnion=-1;
  for (int idx : cluster ) {
    if ( idx %2 == 1) {
      float d = distanceMatrix[cation][idx];
      if (d<minDist) {
        minDist = d;
        nnAnion=idx;
      }
    }
  }
  return minDist;
  if (nnAnion == -1 ) {
    std::cout << "Error, finding nearest anion inside didn't work properply " << "\n";
    nnAnion=1;
    minDist=2;
  }
}

float computeAngle(const std::vector<std::vector<float>>& ions, int cation, int anion, float box) {
  float dx, dy, dz;
  dx=applyPBC( ions[anion][0] - ions[cation][0], box);
  dy=applyPBC( ions[anion][1] - ions[cation][1], box);
  dz=applyPBC( ions[anion][2] - ions[cation][2], box*2);
  float norm = std::sqrt( dx*dx + dy*dy +dz*dz);
  if (norm == 0.0){
    std::cout<< "Zero-length vector \n";
  }

  float cosTheta = dz/norm;
  if (cosTheta > 1.0) cosTheta=1.0;
  if (cosTheta <-1.0) cosTheta=-1.0;

  float angleRad = std::acos(cosTheta);
  float angleDeg = angleRad * ( 180.0/M_PI);

  return angleDeg;
}


std::vector<float> findNN(
    const std::vector<std::vector<float>>& ions,
    const std::vector<std::vector<int>>& clusters,
    const std::vector<std::vector<float>> distanceMatrix,
    int numatoms,
    std::vector<int>& role,
    std::vector<int>& swapcation,
    CationTagMap& tagMap,
    std::vector<int>& anionIdx,
    std::vector<float>& nnangl,
    float box
    ) {

  std::vector<float> nndist(ions.size(), 0.0); // return vector

  // atom index -> cluster index
  std::vector<int> atom2cluster(ions.size(), -1);
  for (size_t c = 0; c < clusters.size(); ++c) {
    for( int atom : clusters[c]) {
      atom2cluster[atom] = c;
    }
  }

  for( int cation=0; cation < ions.size(); cation +=2) {
    int clusterIdx = atom2cluster[cation];
    const auto& cluster = clusters[clusterIdx];

    // cluster analysis with size 1
    if (cluster.size() == 1) {
      nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster,distanceMatrix, cation, anionIdx[cation]); 
      nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);

      if (role[cation] %2 == 0 && role[cation] > 2 ) {
        std::swap(tagMap[cation], tagMap[swapcation[cation]]);
      }
      role[cation] = 1; // update role as singleton
    }
    // cluster analysis with size 2
    else if (cluster.size() == 2) {
      int other = (cluster[0] == cation) ? cluster[1] : cluster[0];
      if (other % 2 == 1) { // + and -
        nndist[cation] = distanceMatrix[cation][other];
        role[cation] = 2; // update role as dipole
        anionIdx[cation] = other;
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
      }
      else { // ++
        nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster,distanceMatrix, cation, anionIdx[cation]); 
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation] = 1; // update role as singleton
      }
    }
    // cluster analysis with size 3
    else if (cluster.size() == 3) { 
      if (cluster[0]%2==0 && cluster[1]%2==0 && cluster[2]%2==0) { // +++
        nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster,distanceMatrix, cation, anionIdx[cation]); 
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=1; // update role as singleton
      } 
      else if ( cluster[0]%2 + cluster[1]%2 + cluster[2]%2 == 2 ) { // +-- always dipole
        for ( int j =0; j<3; j++ ) {
          if (cluster[j] == cation) {
            int a= (j+1)%3;
            int b= (j+2)%3;
            float d1 = distanceMatrix[cation][cluster[a]];
            float d2 = distanceMatrix[cation][cluster[b]];
            if (d1 >= d2 ) {
              nndist[cation] = d2;
              anionIdx[cation] = cluster[b];
              nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
            } else {
              nndist[cation] = d1;
              anionIdx[cation] = cluster[a];
              nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
            }
            role[cation]=2; 
          }
        }
      } else { // ++-  possibility of cation exchange exist
        // find two cations and check if their role is both 4 > make one 3
        int othercation=-1;
        for ( int j =0; j<3; j++ ) {
          if (cluster[j] == cation) {
            int a= (j+1)%3;
            int b= (j+2)%3;
            if (cluster[a]%2 == 0) {
              othercation=cluster[a];
            } else {
              othercation=cluster[b];
            }
          }
        }
        if (cation < othercation && role[cation] %2 == 0 && role[othercation] %2 == 0) { //  if both cation are binary, since it made from ++--
          role[cation]=1;
        } else if (cation < othercation && role[cation] %2 == 1 && role[othercation] %2 == 1) {
          role[cation]=2;
        }

        if( role[cation] %2 == 1) { // if it was singleton or 3-cluster singleton
          nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster,distanceMatrix, cation, anionIdx[cation]); 
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation] = 3; //update role as 3-body singleton
          for ( int j =0; j<3; j++ ) {
            if (cluster[j] == cation) {
              int a= (j+1)%3;
              int b= (j+2)%3;
              if (cluster[a]%2 == 0 ) {
                swapcation[cation] = cluster[a];
              } else {
                swapcation[cation] = cluster[b];
              }
            }
          }
        } 
        else if (role[cation] %2 == 0) { // dipole charactor
          for ( int j =0; j<3; j++ ) {
            if (cluster[j] == cation) {
              int a= (j+1)%3;
              int b= (j+2)%3;
              if (cluster[a]%2 == 0 ) {
                swapcation[cation] = cluster[a];
                nndist[cation] = distanceMatrix[cation][cluster[b]];
                anionIdx[cation] = cluster[b];
                nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
              } else {
                swapcation[cation] = cluster[b];
                nndist[cation] = distanceMatrix[cation][cluster[a]];
                anionIdx[cation] = cluster[a];
                nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
              }
            }
          }
          role[cation] =4;
        }
      }
    }
    //cluster analysis with size 4
    else if (cluster.size() == 4) { 
      if ( cluster[0]%2 + cluster[1]%2 + cluster[2]%2 + cluster[3]%2 == 0 ) { //++++
        nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster,distanceMatrix, cation, anionIdx[cation]); 
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=1;
      }
      else if ( cluster[0]%2 + cluster[1]%2 + cluster[2]%2 + cluster[3]%2 == 3 ) { //+---
        nndist[cation] = findNearestAnionInside(ions, cluster, atom2cluster,distanceMatrix, cation, anionIdx[cation]); 
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=2;
      }
      else if ( cluster[0]%2 + cluster[1]%2 + cluster[2]%2 + cluster[3]%2 == 2 ) { //++--
        nndist[cation] = findNearestAnionInside(ions, cluster, atom2cluster,distanceMatrix, cation, anionIdx[cation]); 
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=2;
      }
      else if ( cluster[0]%2 + cluster[1]%2 + cluster[2]%2 + cluster[3]%2 == 1 ) { //+++-
        if ( role[cation] %2 == 1){
          nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster,distanceMatrix, cation, anionIdx[cation]); 
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation]=5;
        }
        else if ( role[cation] %2 == 0 ) {
          float minDist = 1e9;
          int nnAnion=-1;
          for ( int idx : cluster ) {
            if (idx %2 == 1) { // find anion
              float d = distanceMatrix[cation][idx];
              if (d<minDist) {
                minDist = d;
                nnAnion=idx;
              }
            } else if (idx != cation && idx %2 == 0 && role[idx] %2 == 1) { //find another singleton cation
              swapcation[cation]=idx;
            }
          }
          nndist[cation]=minDist;
          anionIdx[cation] = nnAnion;
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation]=6;
        }
      }
    }
    // cluster analysis with size 5
    else if (cluster.size() == 5) { //
      int charge=0; // number of anion in the cluster
      for ( int j=0; j<5; j++) {
        charge += cluster[j]%2 ;
      }

      if (charge == 0) { // +++++
        nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=1;
      }
      else if (charge == 4 || charge == 3 ) { // +---- or ++---
        nndist[cation] = findNearestAnionInside(ions, cluster, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=2;
      }
      else if ( charge ==2  || charge == 1) { // +++-- or ++++-
        if ( role[cation] %2 == 1){
          nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation]=5;
        }
        else if ( role[cation] %2 == 0 ) {
          float minDist = 1e9;
          int nnAnion=-1;
          for ( int idx : cluster ) {
            if (idx %2 == 1) { // find anion
              float d = distanceMatrix[cation][idx];
              if (d<minDist) {
                minDist = d;
                nnAnion=idx;
              }
            } else if (idx != cation && idx %2 == 0 && role[idx] %2 == 1) { //find another singleton cation
              swapcation[cation]=idx;
            }
          }
          nndist[cation]=minDist;
          anionIdx[cation] = nnAnion;
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation]=6;
        }
      }
    }
    // cluster analysis with size 6
    else if (cluster.size() == 6) { //
      int charge=0; // number of anion in the cluster
      for ( int j=0; j<6; j++) {
        charge += cluster[j]%2 ;
      }

      if (charge == 0) { // ++++++
        nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=1;
      }
      else if (charge == 5 || charge == 4 || charge == 3 ) { // +----- or ++---- or +++---
        nndist[cation] = findNearestAnionInside(ions, cluster, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=2;
      }
      else if ( charge ==2  || charge == 1) { // +++ +-- or +++ ++-
        if ( role[cation] %2 == 1){
          nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation]=5;
        }
        else if ( role[cation] %2 == 0 ) {
          float minDist = 1e9;
          int nnAnion=-1;
          for ( int idx : cluster ) {
            if (idx %2 == 1) { // find anion
              float d = distanceMatrix[cation][idx];
              if (d<minDist) {
                minDist = d;
                nnAnion=idx;
              }
            } else if (idx != cation && idx %2 == 0 && role[idx] %2 == 1) { //find another singleton cation
              swapcation[cation]=idx;
            }
          }
          nndist[cation]=minDist;
          anionIdx[cation]=nnAnion;
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation]=6;
        }
      }
    }
    else if ( cluster.size() % 2 == 0) {
      int charge=0; // number of anion in the cluster
      for ( int j=0; j<cluster.size(); j++) {
        charge += cluster[j]%2 ;
      }

      if (charge == 0) { // ++++++
        nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=1;
      }
      else if (charge >cluster.size()/2 ) { // +-------- or ++------- or +++------ or ++++----- 
        nndist[cation] = findNearestAnionInside(ions, cluster, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=2;
      }
      else { // ++++ --- or ++++ +-- or ++++ ++-
        if ( role[cation] %2 == 1){
          nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation]=5;
        }
        else if ( role[cation] %2 == 0 ) {
          float minDist = 1e9;
          int nnAnion=-1;
          for ( int idx : cluster ) {
            if (idx %2 == 1) { // find anion
              float d = distanceMatrix[cation][idx];
              if (d<minDist) {
                minDist = d;
                nnAnion=idx;
              }
            } else if (idx != cation && idx %2 == 0 && role[idx] %2 == 1) { //find another singleton cation
              swapcation[cation]=idx;
            }
          }
          nndist[cation]=minDist;
          anionIdx[cation] = nnAnion;
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation]=6;
        }
      }
    }
    else  {
      int charge=0; // number of anion in the cluster
      for ( int j=0; j<cluster.size(); j++) {
        charge += cluster[j]%2 ;
      }

      if (charge == 0) { // ++++++
        nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation],box);
        role[cation]=1;
      }
      else if (charge >=cluster.size()/2 ) { // +-------- or ++------- or +++------ or ++++----- 
        nndist[cation] = findNearestAnionInside(ions, cluster, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
        nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
        role[cation]=2;
      }
      else { // ++++ --- or ++++ +-- or ++++ ++-
        if ( role[cation] %2 == 1){
          nndist[cation] = findNearestAnionOutside(ions, clusters, atom2cluster, distanceMatrix, cation, anionIdx[cation]);
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation]=5;
        }
        else if ( role[cation] %2 == 0 ) {
          float minDist = 1e9;
          int nnAnion=-1;
          for ( int idx : cluster ) {
            if (idx %2 == 1) { // find anion
              float d = distanceMatrix[cation][idx];
              if (d<minDist) {
                minDist = d;
                nnAnion=idx;
              }
            } else if (idx != cation && idx %2 == 0 && role[idx] %2 == 1) { //find another singleton cation
              swapcation[cation]=idx;
            }
          }
          nndist[cation]=minDist;
          anionIdx[cation]=nnAnion;
          nnangl[cation] = computeAngle(ions, cation, anionIdx[cation], box);
          role[cation]=6;
        }
      }
    }
  }
  return nndist;
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

  while ( blockData.size() >= 4) {
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
//    std::cout << blockSize << "\t\t" << nBlocks << "\t\t" << std::setprecision(6) << sem << "\n";

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
//    std::cout << "\nERROR. SEM plateua not detected. Use the largest error\n";
    finalSEM=*std::max_element(semList.begin(), semList.end());
  }
//  std::cout << "Estimated SEM: " << finalSEM << "\n";
  return finalSEM;
}

std::vector<float> getColumn(const std::vector<std::vector<float>>& mat, size_t colIndex, float conversion) {
  std::vector<float> column;
  for (const auto& row : mat) {
    if (colIndex < row.size()) {
      column.push_back(row[colIndex]*conversion);
    } else {
      // Handle error or skip if row too short
      column.push_back(0);  // or throw, or continue
    }
  }
  return column;
}


int main(int argc, char* argv[]) {
  if (argc != 9 ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << " density field cutoff_in(A) cutoff_out(A) box(A) timestep(ps) cat(1)/an(2)  eqtime(ns)\n";
    return 1;
  }

  int numtraj=1;
  int numsnap;
  std::string density=argv[1];
  std::string field=argv[2];
  std::map<std::string, float> str2float = {
    {"00", 0.023*0},
    {"02", 0.023*2},
    {"04", 0.023*4},
    {"06", 0.023*6},
    {"08", 0.023*8},
    {"10", 0.023*10},
    {"12", 0.023*12},
    {"14", 0.023*14},
    {"16", 0.023*16},
    {"18", 0.023*18},
    {"20", 0.023*20},
    {"22", 0.023*22},
  };
  float efield = str2float[field]; //kcal/molA

  float CUTOFFin = std::stof(argv[3]);
  float CUTOFFout = std::stof(argv[4]);
  float box=std::stof(argv[5]);
  float timestep=std::stof(argv[6]); // 50 fs  = 0.05 ps 
  float center = std::stof(argv[7]); // 1 if cation is center, 2 else
  float eqtime = std::stof(argv[8]); // equilibration time in ns unit
  int rerun=0;
  
  int numatoms;
  if ( density == "05" ) {
    numatoms = 32;
  } else if (density == "025") {
    numatoms=16;
  } else {
    std::cout << "Error: not valid density\n";
  }

  std::cout << "\nAnalysis Setups\n";
  std::cout << "numIons: " << numatoms << "\n";
  std::cout << "simulationBox: " << box << " A\tx\t" << box << " A\tx\t" << box*2 << " A\n";
  if ( density =="05") {
    std::cout << "densityIon: " << "0.5 M\n";
  } else if (density == "025" ) {
    std::cout << "densityIon: " << "0.25 M\n";
  } else {
    std::cout << "Error: not valid density\n";
  }
  std::cout << "fieldStrength: " << efield << " meV/A\n";
  std::cout << "domainCutoffs: in " << CUTOFFin << " A\t\tout "<< CUTOFFout << " A\n";
  std::cout << "timeStep: " << timestep << " ps\n";
  if (center < 1.5 ) {
    std::cout << "centerIon : Cation\n";
  } else {
    std::cout << "centerIon : Anion\n";
  }
  std::cout << "\n";

  int numpion = numatoms/2;
  std::vector<std::vector<float>> traj;
  //reading coordinate file
  if (rerun == 0) {
    std::string inputFileName1 = "../data/dumps/dumpD" + density + "E" + field + ".binary";
    traj = readBinary(inputFileName1);
    numsnap = traj.size()/numatoms;
    std::cout << "Read from dump files " << inputFileName1 << " of length " << numsnap*timestep << " ps\n";
  } else {
  //reading coordinate file from two dumps
    std::string inputFileName1 = "../data/dumps/dumpD" + density + "E" + field + ".binary";
    std::vector<std::vector<float>> traj1;
    traj1 = readBinary(inputFileName1);

    std::string inputFileName2 = "../data/dumps/dumpD" + density + "E" + field + ".rbinary";
    std::vector<std::vector<float>> traj2;
    traj2 = readBinary(inputFileName2);

    traj = std::move(traj1);
    traj.insert(traj.end(), std::make_move_iterator(traj2.begin()), std::make_move_iterator(traj2.end()));
    traj1.clear();
    traj1.shrink_to_fit();
    traj2.clear();
    traj2.shrink_to_fit();
    numsnap = traj.size()/numatoms;

    std::cout << "Read from dump files " << inputFileName1 << " and " << inputFileName2 << " of length " << numsnap * timestep << " ps\n";
  }


  std::cout << "Steady state relaxation time: " << eqtime << " ns\n";
  traj.erase(traj.begin(), traj.begin()+ static_cast<int>(1000*eqtime/timestep*numatoms));
  numsnap = traj.size()/numatoms;
  std::cout << "Analysis done on trajectory length " << numsnap*timestep/1000 << " ns\n\n";


  std::vector<std::vector<float>> cdist(numsnap, std::vector<float>( numtraj*numpion, 0.0) ); // conditioned nearest neighbor distance, A unit
  std::vector<std::vector<float>> cangl(numsnap, std::vector<float>( numtraj*numpion, 0.0) ); // conditioned nearest neighbor angle, from positive z axis, degree unit
  std::vector<std::vector<int>> anionIndex(numsnap, std::vector<int>( numtraj*numpion, -1) ); // conditioned nearest neighbor index
  
  CationTagMap tagMap;
  for ( int i = 0; i<numatoms; i+=2 ) {
    tagMap[i]=i;
  }

  std::vector<std::vector<int>> clusterStat(10, std::vector<int>(10, 0));


  std::vector<int> role(numatoms, 1);
  std::vector<int> swapcation(numatoms, 0);
  std::vector<float> prevz(numatoms, 0.0);
  std::vector<int> prevcbin(numatoms, -1);
  int every=10; // frequency of measureing displacement vector for computing current, timestep unit

  // histogram of conductivity
  float cmin=0;
  float cmax=30.0;
  int numcurrent=300;
  float bincurrent = (cmax-cmin)/numcurrent; // 0.1 A interval
  std::vector<float> currentvec(numcurrent, 0.0);
  std::vector<float> current(numatoms, 0.0); //  ion conductivity from each snapshot
  std::vector<std::vector<float>> allCurrents; // time-series of conductivity, later used for error analysis
  std::vector<float> catCurrentvec(numcurrent, 0.0);
  std::vector<int> currentcount(numcurrent, 0);
  std::vector<float> currentargument(numcurrent, 0.0);

  for ( int i =0; i<numcurrent; i++) {
    currentargument[i] = (0.5 + i ) * bincurrent;
  }


  for ( int time = 0; time<numsnap; time++){
    //get one snapshot >> ions
    int init = time*numatoms;
    int init2 = 1; // first element : type
    std::vector<std::vector<float>> ions;
    for ( int i=init; i<init+numatoms;i++) {
      std::vector<float> atoms;
      for ( int j=init2; j<init2+3; j++) {
        atoms.push_back(traj[i][j]);
      }
      ions.push_back(atoms);
    }

    // for anion-centered frame, we simply convert role of cation and anion
    if ( center > 1.5) {
      for ( int cation = 0; cation < numatoms; cation += 2) {
        std::swap(ions[cation], ions[cation+1] );
      }
    }

    // construct a graph and distance matrix
    std::vector<std::vector<int>> graph(numatoms);
    std::vector<std::vector<float>> distanceMatrix(numatoms, std::vector<float>(numatoms, 0.0));
    for ( int i=0; i<numatoms; i++ ){
      for ( int j=i+1; j<numatoms; j++ ){
        float d = distance(ions[i], ions[j], box);
        distanceMatrix[i][j] = d;
        distanceMatrix[j][i] = d;
        if (d < CUTOFFin) {
          graph[i].push_back(j);
          graph[j].push_back(i);
        }
      }
    }

    // construct clusters, BFS algorithm
    std::vector<std::vector<int>> clusters;
    findClusters(graph, clusters, numatoms);


    std::vector<int> anionIdx(numatoms, -1);
    std::vector<float> nnangl(numatoms, 0);
    std::vector<float> nndist=findNN( ions, clusters, distanceMatrix, numatoms, role, swapcation, tagMap, anionIdx, nnangl, box);


    std::vector<int> atom2cluster(ions.size(), -1);
    for (size_t c = 0; c < clusters.size(); ++c) {
      for( int atom : clusters[c]) {
        atom2cluster[atom] = c;
      }
    }

    // cluster charge analysis
    /*
    for ( int c=0 ; c<clusters.size() ; c++) {
      int charge=0;
      for ( int e=0; e<clusters[c].size(); e++) {
        if (clusters[c][e] %2 == 0) { 
          charge += 1;
        } else{
          charge -= 1;
        }
      }
      if (clusters[c].size() > 1) {
        int csizeidx = clusters[c].size()-2;
        int chargeidx = (clusters[c].size() + charge)/2;
        clusterStat[csizeidx][chargeidx] += 1;
      }
    }
    */

    InverseCationTagMap inverseTagMap;
    for ( const std::pair< const int, int>& pair : tagMap)  {
      int atom = pair.first;
      int tag = pair.second;
      inverseTagMap[tag] = atom;
    }

    for ( int tag = 0; tag < numatoms; tag+=2) {
      int cation = inverseTagMap[tag];
      float cd = nndist[cation];
      cdist[time][tag/2]=cd;
      anionIndex[time][tag/2] = anionIdx[cation];
      cangl[time][tag/2] = nnangl[cation];
    }

    // computing current vector, using displacement vector accumulated for "every" timestep
    if ( time >= every && time%every == 0) {
      std::vector<float> dispVec(numatoms, 0.0);
      int aidx=0;
      for ( int atom=0; atom<numatoms; atom++) {
        float disp = (ions[atom][2] - prevz[atom] );
        if (std::abs(disp) > box) {
          if (disp > 0) {
            disp-=box*2;
          } else{
            disp += box*2;
          }
        }
        dispVec[atom] = disp;
        float charge = -2*(atom%2)+1;
        current[atom] = charge*disp/timestep/every; // A/ps unit, not normalized by ion numbers
      }
      allCurrents.push_back(current);
      // compute pair current and conductivity, conditioned over its current state, which eventually goes to overall current if statistics are enough
      for (int cation=0; cation<numatoms; cation+=2) {
        aidx = anionIdx[cation];
        if (aidx <=0 ){
          std::cout << "Err, aidx " << aidx << "\n";
        }
        int cbin = prevcbin[cation]; // initialized to -1
        if ( cbin > 0) {
          currentvec[cbin] += (dispVec[cation] - dispVec[aidx] ) / timestep/2/every;
          catCurrentvec[cbin] += (dispVec[cation]) / timestep/every;
          currentcount[cbin] += 1;
        }
        // update cbin for conditioning ## note that t ~ t+dt averaged displacement is conditioned on state at t
        cbin =(int) ( (nndist[cation] -cmin) / (cmax-cmin) * numcurrent);
        if (cbin <0 or cbin >= numcurrent ) { 
          std::cout << "error occurs in current histogram\n";
        } else {
          prevcbin[cation] = cbin;
        }
      }
      for ( int atom=0; atom<numatoms; atom++ ){
        prevz[atom]=ions[atom][2];
      }
    }
  }
  std::string outputc = "../data/cnnDist/distD" + density + "E" + field + ".binary"; ; // nearest neighbor distance
  std::cout << "Output generated : " << outputc << " Rows: " << cdist.size() << " Cols: " << cdist[0].size() << "\n";
  printBinary(outputc, cdist);

  std::string outputa = "../data/cnnAngle/anglD" + density + "E" + field + ".binary"; ; // nearest neighbor angle
  std::cout << "Output generated : " << outputc << " Rows: " << cangl.size() << " Cols: " << cangl[0].size() << "\n";
  printBinary(outputa, cangl);

  std::cout << "\n\n";

  /*
  std::cout << "Cluster charge statistics" << "\n" ;
  for ( int hx=0; hx<clusterStat.size(); hx++) {
    std::cout << "ClusterSize" << hx+2 << "\t";
    for ( int hy=0; hy<hx+3; hy++) {
      std::cout << clusterStat[hx][hy] << "\t";
    }
    std::cout << "\n";
  }
  */

  float conversion = 0; // unit conversion to S cm^2/mol
  if ( center > 1.5 ) {
    conversion = -6.022 * 6.022 * 1.6 * 1.6 /4184 * 10000/efield;
  } else {
    conversion = 6.022 * 6.022 * 1.6 * 1.6 /4184 * 10000/efield;
  }


  std::vector<float> cmeans;
  std::vector<float> csems;
  std::vector<float> ameans;
  std::vector<float> asems;

  for ( int cation=0; cation < numatoms ; cation+=2 ) {
    std::vector<float> onecurrent = getColumn(allCurrents, cation, conversion);
    csems.push_back(blockAverage(onecurrent));
    cmeans.push_back(mean(onecurrent));
//    std::cout << "Ion " << cation << " mean: " << mean(onecurrent) << " sem : " << csems.back() << "\n";
  }
  for ( int anion=1; anion < numatoms ; anion+=2 ) {
    std::vector<float> onecurrent = getColumn(allCurrents, anion, conversion);
    asems.push_back(blockAverage(onecurrent));
    ameans.push_back(mean(onecurrent));
//    std::cout << "Ion " << anion << " mean: " << mean(onecurrent) << " sem : " << asems.back() << "\n";
  }

  float cferr=0;
  float cfmean=0;
  for( int i =0; i<cmeans.size(); i++ ) {
    cfmean += cmeans[i]/csems[i]/csems[i];
    cferr += 1/csems[i]/csems[i];
  }
  cfmean /= cferr;
  cferr = std::sqrt(1/cferr);
  std::cout << "Cation error-weighted mean: " << cfmean << "\tsems:\t" << cferr << "\n";

  float aferr=0;
  float afmean=0;
  for( int i =0; i<ameans.size(); i++ ) {
    afmean += ameans[i]/asems[i]/asems[i];
    aferr += 1/asems[i]/asems[i];
  }
  afmean /= aferr;
  aferr = std::sqrt(1/aferr);
  std::cout << "Anion error-weighted mean: " << afmean << "\tsems:\t" << aferr << "\n";
  std::cout << "Total error-weighted mean: " << (cfmean+afmean)/2 << "\tsems:\t" << std::sqrt((cferr*cferr) + (aferr*aferr))/2 << "\n";


  float CIP=0;
  float SSIP=0;
  float free=0;

  float CIPc=0;
  float SSIPc=0;
  float freec=0;
  for ( int i=0; i< numcurrent; i++ ) {
    if (currentargument[i] < CUTOFFin) {
      CIPc += currentcount[i];
      CIP += currentvec[i] *conversion;
    } else if ( currentargument[i] > CUTOFFout ) {
      freec += currentcount[i];
      free += currentvec[i] * conversion;
    } else {
      SSIPc +=currentcount[i] ;
      SSIP += currentvec[i] * conversion;
    }
  }
  std::cout << "Conductivity CIP " << CIP/CIPc << "\t\tSSIP " << SSIP/SSIPc << "\t\tFree " << free/freec << "\n";
  float all = CIPc + SSIPc + freec;
  std::cout << "Populations  CIP " << CIPc/all << "\t\tSSIP " << SSIPc/all << "\t\tFree " << freec/all << "\n";
  std::cout << "Contribution CIP " << CIP/all << "\t\tSSIP " << SSIP/all << "\t\tFree " << free/all << "\n";

  /*
  for ( int i =0 ; i< numcurrent; i++ ) {
    if ( currentcount[i] != 0 ){
      currentvec[i] = catCurrentvec[i] / currentcount[i] * conversion;
    } else {
      currentvec[i] = 0;
    }
  }
  */
  // generate output result
  std::string resfile = "../results/conductivity/condD" + density + "E" + field;
  std::ofstream out (resfile );
  out << std::fixed << std::setprecision(5) << density << "\t\t" << field << "\t\t"  
    << CUTOFFin << "\t\t" << CUTOFFout <<"\t\t" << numsnap*timestep/1000  << "\t\t"
    << cfmean << "\t\t" << cferr << "\t\t" << afmean << "\t\t" << aferr << "\t\t" 
    << (cfmean + afmean)/2 << "\t\t" << std::sqrt( cferr*cferr+aferr*aferr)/2. << "\t\t" 
    << CIP/CIPc << "\t\t" << SSIP/SSIPc << "\t\t" << free/freec << "\t\t" 
    << CIPc/all << "\t\t" << SSIPc/all << "\t\t" << free/all << "\n";
  out.close();
  //std::string outputcurrent= "./results/condCatCondD" + density + "E" + field +".dat";
  //writeHistogram(outputcurrent, currentargument, currentvec, currentcount);

  return 0;
}

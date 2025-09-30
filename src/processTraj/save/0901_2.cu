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
#include <random>
#include <cassert>
#include <numeric>

////////// Name space ////////////
using TagMap = std::unordered_map<int, int>;
using InverseTagMap = std::unordered_map<int, int>;
using Rows = std::vector<std::vector<int>>;
////////// End Name Space ///////////////

///////// Structures ////////////////////

struct Atom { 
  int id;
  int charge;
  int role; // role as separated ion or bound ion
  int swap; // if 1, need swap
  int tag; // printing order
  float nndist; // conditioned nn distance
  float nnangl; // conditioned nn angle
  int nncounter; // index for pairing counter ion
  float x, y, z;

  Atom(int id_, int charge_, float x_, float y_, float z_) 
    : id(id_), charge(charge_), role(-1), swap(-1), tag(-1), nndist(-1), nnangl(-1), nncounter(-1), x(x_), y(y_), z(z_) {}
};

struct Box {
  float x;
  float y;
  float z;
};

struct ChildClassification {
  int kind;
  std::vector<std::pair<int, int>> parents;
};

////////// End structures //////////////////

////////////    IO
void writeLammpsDump(std::ostream& out, long timestep, const std::vector<Atom>& atoms, const Box& box) {
  std::vector<Atom> sorted = atoms;
  std::sort(sorted.begin(), sorted.end(), [](const Atom& a, const Atom& b){ return a.id < b.id; });
  out << "ITEM: TIMESTEP\n" << timestep << "\n";
  out << "ITEM: NUMBER OF ATOMS\n" << sorted.size() << "\n";
  out << "ITEM: BOX BOUNDS pp pp pp\n";
  out << std::setprecision(10) << std::fixed;
  float bl=0;
  out << bl << " " << box.x << "\n";
  out << bl << " " << box.y << "\n";
  out << bl << " " << box.z << "\n";

  out << "ITEM: ATOMS id mol type charge v_user x y z\n";
  out << std::setprecision(6) << std::fixed;
  for (const auto& a : sorted) {
    int type = (a.charge > 0 ? 1 : 2 ); 
    out << a.id+1 << " " << a.role << " " << type << " " <<  a.charge << " " << a.role << " "  << a.x << " " << a.y << " " << a.z << "\n";
  }
}
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
//////////////// END IO
///////////////// Basic trajectory analysis tools
float applyPBC(float x, float box){
  float hbox = box/2.0;
  float wrapped = fmod(x + hbox, box);
  if (wrapped < 0) wrapped += box;
  return wrapped - hbox;
}

float distance(const Atom& a, const Atom& b, const Box& box) {
  float dx, dy, dz, rsq;
  dx=applyPBC( a.x - b.x, box.x);
  dy=applyPBC( a.y - b.y, box.y);
  dz=applyPBC( a.z - b.z, box.z);
  rsq = dx*dx +dy*dy +dz*dz;
  return std::sqrt(rsq);
}
float computeAngle(const std::vector<Atom>& atoms, int cation, int anion, const Box& box) {
  float dx, dy, dz;
  dx=applyPBC( atoms[anion].x - atoms[cation].x, box.x);
  dy=applyPBC( atoms[anion].y - atoms[cation].y, box.y);
  dz=applyPBC( atoms[anion].z - atoms[cation].z, box.z);
  float norm = std::sqrt( dx*dx + dy*dy +dz*dz);
  if (norm == 0.0){
    std::cout<< "Error! Zero-length vector \n";
  }

  float cosTheta = dz/norm;
  if (cosTheta > 1.0) cosTheta=1.0;
  if (cosTheta <-1.0) cosTheta=-1.0;

  float angleRad = std::acos(cosTheta);
  float angleDeg = angleRad * ( 180.0/M_PI);

  return angleDeg;
}
///////////////// End basic trajectory analysis

//////////////// Construct Cluster 
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
int netCharge(const std::vector<int>& cluster) {
  int net=0;
  for( int i =0; i<cluster.size(); i++) {
    if ( cluster[i] % 2 == 0 ) {
      net +=1;
    } else {
      net -= 1;
    }
  }
  return net;
}
/////////////////// End Cluster

/////////////////// Cluster History Analysis
// Overlap map : prevIdx -> [(currIdx, overlap)]
// currIdx -> [(prevIdx, overlap)]
using P2C = std::unordered_map<int, std::vector<std::pair<int,int>>>;
using C2P = std::unordered_map<int, std::vector<std::pair<int,int>>>;

static void build_overlap_maps(const Rows& prevC, const Rows& currC, P2C& p2c, C2P& c2p) {
  p2c.clear(); c2p.clear();

  // atom -> prev cluster index
  std::unordered_map<int, int> atom2prev;
  atom2prev.reserve(100);
  for ( int i =0; i < prevC.size(); i++) {
    for ( int a : prevC[i]) atom2prev[a] = i;
  }

  for (int j=0; j< currC.size(); j++) {
    const auto& child = currC[j];
    std::unordered_map<int, int> contrib; // prev -> overlap count
    contrib.reserve(8);
    for (int a: child) {
      auto it = atom2prev.find(a);
      if (it != atom2prev.end()) ++contrib[it->second];
      // contrib : frequency histogram of previous cluster index
    }
    for (auto& kv : contrib){
      p2c[kv.first].push_back({j, kv.second});
      c2p[j].push_back({kv.first, kv.second});
    }
  }
}

// BFS the small "overlap component" around a child
static void component_for_child(int childIdx, const P2C& p2c, const C2P& c2p, 
                                std::unordered_set<int>& prevs,
                                std::unordered_set<int>& currs)
{
  // currs : all child index that are connected to childIdx through any parent set
  // prevs : all previous cluster indexes composing current cluster
  prevs.clear(); currs.clear();
  std::queue<std::pair<bool, int>> q; // (isCurr, index)
  currs.insert(childIdx);
  q.push({true, childIdx});
  while(!q.empty()){
    auto p = q.front(); q.pop();
    bool isCurr = p.first;
    int idx = p.second;
    if (isCurr) {
      auto it = c2p.find(idx);
      if ( it==c2p.end()) continue;
      for ( auto& pr : it->second){
        int p= pr.first;
        if (prevs.insert(p).second) {
          auto jt = p2c.find(p);
          if (jt != p2c.end()){
            for (auto& ch : jt->second){
              int c = ch.first;
              if (currs.insert(c).second) q.push({true, c});
            }
          }
        }
      }
    } else {
      auto jt = p2c.find(idx);
      if (jt == p2c.end()) continue;
      for ( auto& ch : jt->second){
        int c=ch.first;
        if (currs.insert(c).second) q.push({true, c});
      }
    }
  }
}

// exact set-equality check
static bool same_set(const std::vector<int>& a, const std::vector<int>& b ){
  return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
}
// classify, 1: merge, 2: split, 3: exchange, 4:continuation
// classify one current cluster indexed by j
static ChildClassification classify_child(int j, const Rows& prevC, const Rows& currC, const P2C& p2c, const C2P& c2p) {
  ChildClassification res;
  res.parents.clear();

  auto it = c2p.find(j);
  int indeg = (it==c2p.end() ? 0 : (int)it->second.size());
  if (indeg == 0) {
    res.kind=0;
    return res;
  }

  res.parents = it->second;

  int parent = (indeg==1 ? it->second[0].first : -1 );
  int parentOutdeg = (indeg==1 ? (int)p2c.at(parent).size() : 0);

  std::unordered_set<int> compPrev, compCurr;
  component_for_child(j, p2c, c2p, compPrev, compCurr);
  bool hasSplit = false, hasMerge = false;
  for (int p : compPrev) {
    auto jt = p2c.find(p);
    if (jt != p2c.end() && (int)jt->second.size() >=2) {hasSplit = true; break;}
  }
  for (int c : compCurr) {
    auto kt = c2p.find(c);
    if (kt != c2p.end() && (int)kt->second.size() >=2) {hasMerge = true; break;}
  }
  if ( hasSplit && hasMerge && compPrev.size() >=2 && compCurr.size() >=2 ) { 
    res.kind = 3;
    return res;
  }
 
  if (indeg >=2 ) {
    res.kind = 1; return res;
  }
  if (parentOutdeg >= 2) {
    res.kind=2; return res;
  }
  if (same_set(prevC[parent], currC[j])) {
    res.kind=4; return res;
  }
  res.kind=1;
  return res;
}

static std::vector<ChildClassification> classify_all(const Rows& prev, const Rows& curr) {
  P2C p2c; C2P c2p; build_overlap_maps(prev, curr, p2c, c2p);

  std::vector<ChildClassification> out;
  out.reserve(curr.size());
  for (int j=0; j< curr.size(); j++) 
    out.push_back( classify_child(j, prev, curr, p2c, c2p));
  return out;
}

static void print_set(const std::vector<int>& v){
  std::cout << "{";
  for (size_t i=0; i<v.size(); i++ ){ if (i) std::cout<<","; std::cout<<v[i];}
  std::cout << "}";
}
////////////// End Cluster history analysis

std::vector<int> findLeastConnectedIons(int n, int sign, const std::vector<int>& cluster, const Rows& graph, const std::vector<Atom>& atoms ){
  if (n<=0) return {};
  
  struct Item { int cnt; int idx;};
  std::vector<Item> items;
  items.reserve(cluster.size());

  for( int u : cluster) {
    if (atoms[u].charge != sign) continue; // only ions who have same charge as cluster
    int cnt = 0; // connected counter ions
    for (int v : graph[u]) {
      if (atoms[v].charge != sign) cnt++;
    }
    items.push_back({cnt, u});
  }
  if (items.empty()) return {};

  auto cmp = [](const Item& a, const Item& b) {
    if (a.cnt != b.cnt) return a.cnt < b.cnt;
    return a.idx < b.idx;
  };

  if( n>= static_cast<int>(items.size())){
    std::sort(items.begin(), items.end(), cmp);
  } else{
    std::nth_element(items.begin(), items.begin() + n, items.end(), cmp);
    std::sort(items.begin(), items.begin() + n, cmp);
    items.resize(n);
  }
  std::vector<int> result;
  result.reserve(std::min<int>(n, items.size()));
  for(const auto& it : items) result.push_back(it.idx);
  return result;
}




void findCounterOutside (int index, std::vector<Atom>& atoms, const std::vector<std::vector<int>>& clusters, const std::vector<std::vector<float>>& distanceMatrix, const Box& box) {
  assert(index >=0 && static_cast<size_t>(index) < atoms.size());
  const int myCharge = atoms[index].charge;
  const int needSign = (myCharge > 0 ? -1 : +1);

  int myCluster = -1;
  for(int cl=0; cl<clusters.size() && myCluster <0; cl++ ){
    for (int i: clusters[cl]){
      if (i == index) { myCluster = cl; break;}
    }
  }
  if (myCluster <0 ) std::cout << "Error, current cluster not found\n";

  int bestCluster = -1;
  float bestClusterDist = std::numeric_limits<float>::infinity();

  for ( int cl = 0; cl < clusters.size(); cl++ ) {
    if ( cl== myCluster) continue;

    int net = netCharge(clusters[cl]);
    if ((needSign > 0 && net <= 0) || (needSign < 0 && net >= 0)) {
      continue;
    }
    float minD = std::numeric_limits<float>::infinity();
    for (int j : clusters[cl]){
      float d  = distanceMatrix[index][j];
      if ( d< minD) minD = d;
    }
    if( minD < bestClusterDist) {
      bestClusterDist = minD;
      bestCluster = cl;
    }
  }
  if (bestCluster < 0 ) {
    atoms[index].nncounter = -1;
    atoms[index].nndist = -1.0;
    std::cout << "Error, counter ion not found\n";
  }

  int bestAtom = -1;
  float bestD = std::numeric_limits<float>::infinity();

  for( int j : clusters[bestCluster]) {
    //if (atoms[j].charge == needSign && atoms[j].role == 1 ) { 
    if (atoms[j].charge == needSign ) { // role 2 can be a candidate 
      float d = distanceMatrix[index][j];
      if (d<bestD) { bestD = d; bestAtom = j;}
    }
  }
  if (bestAtom < 0) {
    std::cout << "Out Error, cannot find best ion for " << index << "\n\n";
  }
  atoms[index].nncounter = bestAtom;
  atoms[index].nndist = bestD;
  atoms[index].nnangl = computeAngle(atoms, index, bestAtom, box);
}

float FfindCounterOutside (int index, std::vector<Atom>& atoms, const std::vector<std::vector<int>>& clusters, const std::vector<std::vector<float>>& distanceMatrix, const Box& box) {
  assert(index >=0 && static_cast<size_t>(index) < atoms.size());
  const int myCharge = atoms[index].charge;
  const int needSign = (myCharge > 0 ? -1 : +1);

  int myCluster = -1;
  for(int cl=0; cl<clusters.size() && myCluster <0; cl++ ){
    for (int i: clusters[cl]){
      if (i == index) { myCluster = cl; break;}
    }
  }
  if (myCluster <0 ) std::cout << "Error, current cluster not found\n";

  int bestCluster = -1;
  float bestClusterDist = std::numeric_limits<float>::infinity();

  for ( int cl = 0; cl < clusters.size(); cl++ ) {
    if ( cl== myCluster) continue;

    int net = netCharge(clusters[cl]);
    if ((needSign > 0 && net <= 0) || (needSign < 0 && net >= 0)) {
      continue;
    }
    float minD = std::numeric_limits<float>::infinity();
    for (int j : clusters[cl]){
      float d  = distanceMatrix[index][j];
      if ( d< minD) minD = d;
    }
    if( minD < bestClusterDist) {
      bestClusterDist = minD;
      bestCluster = cl;
    }
  }
  if (bestCluster < 0 ) {
    atoms[index].nncounter = -1;
    atoms[index].nndist = -1.0;
    std::cout << "Error, counter ion not found\n";
  }

  int bestAtom = -1;
  float bestD = std::numeric_limits<float>::infinity();

  for( int j : clusters[bestCluster]) {
    //if (atoms[j].charge == needSign && atoms[j].role == 1 ) { 
    if (atoms[j].charge == needSign ) { // role 2 can be a candidate 
      float d = distanceMatrix[index][j];
      if (d<bestD) { bestD = d; bestAtom = j;}
    }
  }
  if (bestAtom < 0) {
    std::cout << "Out Error, cannot find best ion for " << index << "\n\n";
  }
  return bestD;
}


void findCounterInside (int index, std::vector<Atom>& atoms, const std::vector<std::vector<int>>& clusters, const std::vector<std::vector<float>>& distanceMatrix, const Box& box) {
  assert(index >=0 && static_cast<size_t>(index) < atoms.size());
  const int myCharge = atoms[index].charge;
  const int needSign = (myCharge > 0 ? -1 : +1);

  int myCluster = -1;
  for(int cl=0; cl<clusters.size() && myCluster <0; cl++ ){
    for (int i: clusters[cl]){
      if (i == index) { myCluster = cl; break;}
    }
  }
  if (myCluster <0 ) std::cout << "Error, current cluster not found\n";

  int bestAtom = -1;
  float bestD = std::numeric_limits<float>::infinity();

  for( int j : clusters[myCluster]) {
    if ( j== index) continue;
    if (atoms[j].charge == needSign && atoms[j].role == 2 ) { 
      float d = distanceMatrix[index][j];
      if (d<bestD) { bestD = d; bestAtom = j;}
    }
  }
  if (bestAtom < 0) {
    std::cout << "In Error, cannot find best ion for " << index << "\n\n";
  }
  atoms[index].nncounter = bestAtom;
  atoms[index].nndist = bestD;
  atoms[index].nnangl = computeAngle(atoms, index, bestAtom, box);
}

void findCounterIons (std::vector<Atom>& atoms, const std::vector<std::vector<int>>& clusters, const std::vector<std::vector<float>>& distanceMatrix, const Box& box) {
  for ( int atom = 0; atom < atoms.size(); atom++ ) {
    if ( atoms[atom].role == 1 ) {
      findCounterOutside(atom, atoms, clusters, distanceMatrix, box); // find counter ion, measure distance and angle
    } else if ( atoms[atom].role == 2 ) {
      findCounterInside(atom, atoms, clusters, distanceMatrix, box);
    } else {
      std::cout << "Role is not defined for " << atom << "\n";
    }
  }
}



void roleAssignment(const std::vector<std::vector<int>>& clusters, std::vector<Atom>& atoms, std::vector<std::vector<int>>& graph, const std::vector<std::vector<float>>& distanceMatrix, const Box& box ){
  for ( const auto& cl : clusters) {
    for ( int member : cl) {
      atoms[member].role = 2;
    }
    int net = netCharge(cl);
    if (net == 0) continue;

    const int q = std::abs(net);
    const int sign = (net > 0 ? +1 : -1);

    std::vector<int> candidates;  //ions having same charge as cluster
    std::vector<int> candidatesNN;  //nearest neighbor opposite charged cluster
    for( int member : cl) {
      if (atoms[member].charge == sign ) candidates.push_back(member);
    }
    // need to assign "q" ions as role 1 from candidates
    if( q==candidates.size()){
      for(int atom : candidates) {
        atoms[atom].role = 1;
      }
    } else {
      std::vector<std::pair<int, float>> id2dist;
      for(int atom : candidates) {
        findCounterOutside(atom, atoms, clusters, distanceMatrix, box);
        id2dist.emplace_back(atom, atoms[atom].nndist);
      }
      std::sort(id2dist.begin(), id2dist.end(), 
          [](const auto& a , const auto& b) {
          return a.second < b.second;//assending order
          });
      for( int n=0; n<q; n++) atoms[id2dist[n].first].role = 1;
    }
  }
}

void updateRoleMerge(const std::vector<std::vector<int>>& clusters, std::vector<Atom>& atoms, std::vector<int>& proles, int time) {
  for ( const auto& cl : clusters) {
    int net = netCharge(cl);
    if (net == 0) continue;

    const int q = std::abs(net);
    const int sign = (net > 0 ? +1 : -1);

    std::vector<int> candidates;
    candidates.reserve(cl.size());
    for( int member : cl) {
      if (atoms[member].charge == sign ) candidates.push_back(member);
    }

    const int take = std::min<int>(q, static_cast<int>(candidates.size()));
    if (take < candidates.size() ) { // when ambiguity in determining role
      std::vector<int> deltaRole;
      for ( int k : candidates) deltaRole.push_back( atoms[k].role - proles[k]);
      int nonzeros = std::count_if(deltaRole.begin(), deltaRole.end(), [](int x) { return x != 0; } );
      int sum = std::accumulate(deltaRole.begin(), deltaRole.end(), 0);
      if ( nonzeros >= 2 ) {
        int first = deltaRole.size();
        int second = deltaRole.size();
        for ( int i = 0; i<deltaRole.size(); i++ ) {
          if (deltaRole[i] != 0) {
            if ( first == deltaRole.size()){ first = i;}
            else { second = i; break; }
          }
        }
        /*
        std::cout << "Time " << time << "\n";
        std::cout << "cluster size " << cl.size() << " candidates " << candidates.size() << " take " << take << " nonzeros " << nonzeros << " sum " << sum << std::endl;
        std::cout << "first : " << first << " second : " << second << "\n";
        std::cout << "first : " << candidates[first] << " second : " << candidates[second] << "\n";
        */
        std::swap(atoms[candidates[first]].role, atoms[candidates[second]].role );
      }
    }
  }
}
void updateTagSplit(const std::vector<std::vector<int>>& clusters, std::vector<Atom>& atoms, std::vector<int>& proles, int time) {
  for ( const auto& cl : clusters) {
    int net = netCharge(cl);
    if (net == 0) continue;

    const int q = std::abs(net);
    const int sign = (net > 0 ? +1 : -1);

    std::vector<int> candidates;
    candidates.reserve(cl.size());
    for( int member : cl) {
      if (atoms[member].charge == sign ) candidates.push_back(member);
    }

    const int take = std::min<int>(q, static_cast<int>(candidates.size()));
    if (take < candidates.size() ) { // when ambiguity in determining role
      std::vector<int> deltaRole;
      for ( int k : candidates) deltaRole.push_back( atoms[k].role - proles[k]);
      int nonzeros = std::count_if(deltaRole.begin(), deltaRole.end(), [](int x) { return x != 0; } );
      int sum = std::accumulate(deltaRole.begin(), deltaRole.end(), 0);
      if ( nonzeros >= 2 ) {
        int first = deltaRole.size();
        int second = deltaRole.size();
        for ( int i = 0; i<deltaRole.size(); i++ ) {
          if (deltaRole[i] != 0) {
            if ( first == deltaRole.size()){ first = i;}
            else { second = i; break; }
          }
        }
        /*
        std::cout << "Time " << time << "\n";
        std::cout << "cluster size " << cl.size() << " candidates " << candidates.size() << " take " << take << " nonzeros " << nonzeros << " sum " << sum << std::endl;
        std::cout << "first : " << first << " second : " << second << "\n";
        std::cout << "first : " << candidates[first] << " second : " << candidates[second] << "\n";
        */
        std::swap(atoms[candidates[first]].role, atoms[candidates[second]].role );
      }
    }
  }
}

std::vector<int> detectJump( const std::vector<Atom>& atoms, const std::vector<float>& pdist ,const std::vector<int> ptags, float cutoff) {
  std::vector<int> idx;
  for ( int i =0; i<pdist.size(); i++){ // here i is tag 
    // find previous particle with tag i
    float prevDist = 0;
    auto it = std::find(ptags.begin(), ptags.end(), i);
    if ( it == ptags.end() ) {
      std::cout << "Previous, cannot find properly tagged particle\n\n";
    } else {
      int index = std::distance(ptags.begin(), it);
      prevDist = pdist[index];
    }

    float currDist = 0;
    auto cit = std::find_if( atoms.begin(), atoms.end(), [i](const Atom& a ){ return a.tag == i; });
    if (cit == atoms.end()) { 
      std::cout << "Current  cannot find properly tagged particle\n\n";
    } else {
      currDist = cit->nndist;
    }
      
    if (std::fabs(currDist - prevDist) > cutoff ) {
      idx.push_back(i);
    }
  }
  return idx;
}



////////// Statistics
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

///////////End statistics

int main(int argc, char* argv[]) {
  if (argc != 11 ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << "dir_name density field cutoff_in(A) cutoff_out(A) boxX(A) boxY(A) boxZ(A) timestep(ps)  eqtime(ns) \n";
    //std::cerr << "numberOfEnsemble is non-unity if there are parallel trajectories. Corresponding trajectories are named with final integer if it is not 1.\n";

    return 1;
  }

  std::string dirName=argv[1]; // name of directory to output results, in results directory
  std::string rho = argv[2]; 
  float density = std::stof(rho)* 0.1; // density in M unit
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
  //int dumpOrder = std::stoi(argv[11]); // 1 if cation and anion alternates, 0 otherwise
  //int numtraj = std::stoi(argv[12]); // if more than 1, nump file names have integer at the end.


  std::cout << "\nAnalysis Setups\n";
  std::cout << "dump Directory : ../data/dumps" << dirName << "/\n"; 
  std::cout << "output Directory : ../results/" << dirName << "/\n"; 
  //std::cout << "Units : A, ns, kcal/mol, mol/L\n";
  std::cout << "numIons: " << numatoms << " numPairs: " << numpion << "\n";
  std::cout << "boxSizes x: " << boxX << " A y: " << boxY << " A z: " << boxZ << " A \n";
  std::cout << "density of ions : " << density << "M \n";
  std::cout << "fieldStrength: " << field << " kcal/molA = " << field * fconversion << " mV/A" << "\n"; 

  std::cout << "domainCutoffs: in " << CUTOFFin << " A\t\tout "<< CUTOFFout << " A\n";
  std::cout << "timeStep: " << timestep << " ps\n";
  std::cout << "eqTime: " << eqtime << " ns\n";
  /*
  if (dumpOrder == 1) {
    std::cout << "dump file is ordered alternatively\n";
  } else if (dumpOrder == 0){
    std::cout << "dump file is ordered sequentially\n";
  }
  */
  //std::cout << "Number of Ensembles : " << numtraj << "\n";

  std::cout << "\n\n";


  // Get dump files for numberOfEnsemble = 1 case
  std::cout << "Dump Reading Starts\n";
  int numsnap;
  std::vector<std::vector<float>> traj;
  //reading coordinate file
  std::string inputFileName = "../data/dumps" + dirName + "/dumpD" + rho + "E" + fieldname + ".binary";
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


  //check the order of types
  bool alter;
  if (traj[0][0] != traj[1][0] ) {
    std::cout << "Ions are ordered alternatively (+-+-) \n";
    alter=true;
  } else {
    std::cout << "Ions are ordered sequentially (++--) \n";
    alter=false;
  }

  /*
  for ( int time = 0; time < 2; time++ ) {
    for ( int index = 0; index < numatoms; index++ ) {
      int init = time*numatoms;
      std::cout << "Atom" << index+1 << "\t" << traj[init+index][0] << "\t" 
        << traj[init+index][1] << "\t" << traj[init+index][2] << "\t" << traj[init+index][3] << "\n";
    }
    std::cout << "\n";
  }
  */

  std::vector<std::vector<float>>  cdist(numsnap, std::vector<float>( numpion, 0.0) ); // conditioned nearest neighbor distance, A unit
  std::vector<std::vector<float>>  cangl(numsnap, std::vector<float>( numpion, 0.0) ); // conditioned nearest neighbor angle, from positive z axis, degree unit
  std::vector<std::vector<int>> anionIndex(numsnap, std::vector<int>( numpion, -1) ); // conditioned nearest neighbor index
  std::vector<std::vector<float>> allCurrent(numsnap-1, std::vector<float>( numatoms, 0.0)); 

  
  std::vector<int> proles(numatoms, -1);
  std::vector<int> ptags(numatoms);
  std::iota(ptags.begin(), ptags.end(), 0);
  std::vector<float> pdist(numatoms, -1);
  std::vector<std::vector<int>> pclusters;
  std::vector<float> prevz(numatoms, 0.0);


  // save visualization for each time step
  std::ofstream ofs("video.lammpstrj");
  if (!ofs) {
    std::cerr << "Error. Cannot open the file\n";
    return 1;
  }
  
  //numsnap = 200000;
  //for ( int time = 0; time<numsnap; time++){
  int st=40000;
  for (int time= st; time<st+10000; time++){
    std::vector<Atom> frame; // id charge role swap tag nndist nnangl nnanion x y z
    //get one snapshot >> ions
    int init = time*numatoms;
    if (alter ) {
      for ( int i=init; i<init+numatoms;i++) {
        int charge = ( (i-init) % 2 == 0 ? +1 : -1);
        frame.push_back({i-init, charge , traj[i][1], traj[i][2], traj[i][3]  });
      }
    } else { // +++++ ----
      for ( int i=init; i<init+numpion;i++) {
        frame.push_back({2*(i-init),+1, traj[i][1], traj[i][2], traj[i][3]  });
        frame.push_back({2*(i-init)+1, -1, traj[i+numpion][1], traj[i+numpion][2], traj[i+numpion][3]  });
      }
    }

    // construct a graph and distance matrix
    std::vector<std::vector<int>> graph(numatoms);
    std::vector<std::vector<float>> distanceMatrix(numatoms, std::vector<float>(numatoms, 0.0));
    for ( int i=0; i<numatoms; i++ ){
      for ( int j=i+1; j<numatoms; j++ ){
        float d = distance(frame[i], frame[j], box);
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

    // Analysis
    if ( time == st ) {
      roleAssignment(clusters, frame, graph, distanceMatrix, box); // frame[i].role determines whether it finds counter ion inside or outside of cluster
      findCounterIons(frame, clusters, distanceMatrix, box); // based on role, find conditioned counter ion (index, dist, angl)
      for(auto& atom : frame) atom.tag = atom.id; // initialize tag

      // save role, tag and cluster info 
      std::transform(frame.begin(), frame.end(), proles.begin(),[](const Atom& a){ return a.role;}); // save roles
      std::transform(frame.begin(), frame.end(), ptags.begin(),[](const Atom& a){ return a.tag;}); // save tags
      std::transform(frame.begin(), frame.end(), pdist.begin(),[](const Atom& a){ return a.nndist;}); // save prev nn dist
      pclusters = clusters; // save cluster info
    } else {
      roleAssignment(clusters, frame, graph, distanceMatrix, box);
      updateRoleMerge(clusters, frame,  proles, time); // when merge occur, swap role to reduce discontinuity

      /*
      auto cls = classify_all(pclusters, clusters);
      std::vector<int> candidates;
      for (int j=0; j< clusters.size(); j++ ) {
        if (cls[j].kind == 2 || cls[j].kind == 3) { // split or exchange
          for( int atom : clusters[j]) {
            candidates.push_back(atom) ;
          }
          std::vector<float> shortest;
          std::vector<int> shortestIdx;
          int net = netCharge(clusters[j]);
          int q = std::abs(net);
          int sign = (net > 0 ? +1 : -1);

          if ( q != 0 && q < clusters[j].size() ) { // ambiguity
            for ( int atom : clusters[j]) {
              if( frame[atom].charge == sign ) {
                shortest.push_back( FfindCounterOutside(atom, frame, clusters, distanceMatrix, box ) );
                shortestIdx.push_back(atom);
              }
            }
            //std::cout << "Time " << time << " Cluster " << j << " charge " << net << " size " << clusters[j].size() << "\n"; 
            //for ( int ii =0; ii < shortest.size(); ii++) {
            //  std::cout << "index " << clusters[j][ii] << " role " << frame[clusters[j][ii]].role << " dist " << shortest[ii] << " ";
            //}
            //std::cout << std::endl;
            if ( shortest[0] > shortest[1] && frame[shortestIdx[0]].role > frame[shortestIdx[1]].role ) { // 
              std::swap(frame[shortestIdx[0]].role, frame[shortestIdx[1]].role);
            }
            else if ( shortest[0] < shortest[1] && frame[shortestIdx[0]].role < frame[shortestIdx[1]].role ) { // 
              std::swap(frame[shortestIdx[0]].role, frame[shortestIdx[1]].role);
            }
          }
        }
      }
      */

      //updateRoleSplit(// current cluster, if it is formed from split, cls[j].kind==2, 
      findCounterIons(frame, clusters, distanceMatrix, box);
      for(auto& atom : frame) atom.tag = ptags[atom.id]; // initialize tag
//      updateTagSplit(...); // swap tag if there is discontinuity occur during split or exchange 
      //// update tag
    


      std::vector<int> deltaRole;
      for ( int cand : candidates) {
        int target = frame[cand].tag;
        auto it = std::find(ptags.begin(), ptags.end(), target);
        int pidx = std::distance( ptags.begin(), it);
        deltaRole.push_back( frame[cand].role - proles[pidx]);
      }
      int nonzeros = std::count_if(deltaRole.begin(), deltaRole.end(), [](int x) { return x != 0; } );
      int sum = std::accumulate(deltaRole.begin(), deltaRole.end(), 0) ;

      int first = deltaRole.size();
      int second = deltaRole.size();

      if (nonzeros >= 2) {
        for ( int i =0; i<deltaRole.size(); i++ ){
          if (deltaRole[i] != 0) {
            if (first == deltaRole.size()) { first = i;}
            else {second = i; break;}
          }
        }
      }

      if (!candidates.empty()) { 
//        std::cout << " nonempty at time " << time << "\t";
        for (int j = 0; j< candidates.size(); j++ ) {
//          std::cout << candidates[j] << " ";
        }
        if ( nonzeros >= 2) {
//          std::cout << "/// nonzeros " << candidates[first] << " and " << candidates[second];
          std::swap(frame[candidates[first]].tag, frame[candidates[second]].tag);
        }
//        std::cout << std::endl;
      }

      /// end update tag

      /*
      if ( time==40949 ) {
        std::swap(frame[30].tag, frame[54].tag);
      }
      */

      std::vector<int> jumps = detectJump(frame, pdist, ptags, 3.0);
      if (!jumps.empty()) {
        for(int mem : jumps){
          auto pit = std::find(ptags.begin(), ptags.end(), mem);
          int pindex = std::distance(ptags.begin(), pit);
          auto cit = std::find_if(frame.begin(), frame.end(), [mem](const Atom& a ) { return a.tag == mem;});
          if ( cit->nndist < CUTOFFin || pdist[pindex] < CUTOFFin ) {
            std::cout <<"Time " << time << " atom index " << cit->id << " : " << pdist[pindex] << "(" << pindex << ") to "  << cit->nndist << " \n";
          }
        }
      }

      /*
      if ( time == 40949 || time == 40948 ) {
        auto cls = classify_all(pclusters, clusters);
        std::cout << "Time : " << time << "\n";
        for (int j=0; j<clusters.size(); j++) {
          if (cls[j].kind != 4) {
            std::cout << "child ";
            print_set(clusters[j]);
            std::cout << " -> " << cls[j].kind;

            if (!cls[j].parents.empty()) {
              std::cout << " parents : ";
              for (const auto& pr : cls[j].parents) {
                int pidx = pr.first;
                int ov = pr.second;
                std::cout << "[prev#" << pidx << " ";
                print_set(pclusters[pidx]);
                std::cout << " ov=" << ov << "] ";
              }
            }
            std::cout << "\n";
            for (int atom : clusters[j]) { 
              std::cout << " atom " << atom << " nn " << pdist[atom] << " > " << frame[atom].nndist << "\t";
            }
            std::cout << "\n";


          }
        }
      }
      */
      std::transform(frame.begin(), frame.end(), proles.begin(),[](const Atom& a){ return a.role;}); // save roles
      std::transform(frame.begin(), frame.end(), ptags.begin(),[](const Atom& a){ return a.tag;}); // save tags
      std::transform(frame.begin(), frame.end(), pdist.begin(),[](const Atom& a){ return a.nndist;}); // save prev nn dist
      pclusters = clusters; // save cluster info
      // save role, tag and cluster info 

    }
    writeLammpsDump(ofs, time, frame, box);

    /*
    if (time > 0)  {
      auto cls = classify_all(pclusters, clusters);
      int count=0;
      for (int jj=0; jj<clusters.size(); jj++){
        if (cls[jj].kind != 4 ){
          for( const auto& pr : cls[jj].parents) {
            int pidx = pr.first;
            if (pclusters[pidx].size() %2 == 1) count++;
          }
        }
      }
      if (1){
        std::cout << "Time : " << time << "\n";
        for (int j=0; j<clusters.size(); j++) {
          if (cls[j].kind != 4) {
            std::cout << "child ";
            print_set(clusters[j]);
            std::cout << " -> " << cls[j].kind;

            if (!cls[j].parents.empty()) {
              std::cout << " parents : ";
              for (const auto& pr : cls[j].parents) {
                int pidx = pr.first;
                int ov = pr.second;
                std::cout << "[prev#" << pidx << " ";
                print_set(pclusters[pidx]);
                std::cout << " ov=" << ov << "] ";
              }
            }
            std::cout << "\n";
          }
        }
      }
    } 
    */

    
  
    //std::vector<int> anionIdx(numatoms, -1);
    //std::vector<float> nnangl(numatoms, 0);
    /*
    updateNN( frame, clusters, distanceMatrix, prole, pswap, tagMap, box);
  
  
    std::vector<int> atom2cluster(ions.size(), -1);
    for (size_t c = 0; c < clusters.size(); ++c) {
      for( int atom : clusters[c]) {
        atom2cluster[atom] = c;
      }
    }

  
    // construct inverse tag map after updating tag
    InverseTagMap inverseTagMap;
    for ( const std::pair< const int, int>& pair : tagMap)  {
      int atom = pair.first;
      int tag = pair.second;
      inverseTagMap[tag] = atom;
    }

    if ( time == numsnap-1) {
      for (size_t c = 0; c < clusters.size(); ++c) {
        std::cout << "Cluster\t" << c+1 << "\tAtoms\t";
        for( int atom : clusters[c]) {
          std::cout << inverseTagMap[atom]+1 << " ";
        }
        std::cout << "\n";
      }
    }
  
    for ( int tag = 0; tag < numatoms; tag+=2) {
      int cation = inverseTagMap[tag];
      float cd = nndist[cation];

      if ( alter) {
        cdist[time][tag/2]=cd;
        anionIndex[time][tag/2] = anionIdx[cation];
        cangl[time][tag/2] = nnangl[cation];
      } else {
        cdist[time][numpion+tag/2]=cd;
        anionIndex[time][numpion+tag/2] = anionIdx[cation];
        cangl[time][numpion+tag/2] = nnangl[cation];
      }
    }
    std::vector<Atom> frame;
    for( int atom =0; atom<numatoms; atom++ ) {
      //frame.push_back({ atom+1, inverseTagMap[atom]+1, atom%2+1, ions[atom][0], ions[atom][1], ions[atom][2]});
      // id type role swap tag nndist nnangl nnanion
      frame.push_back({ atom+1, atom%2+1, role[atom], swapcation[atom], tagMap[atom], nndist[atom], nnangl[atom], anionIdx[atom],ions[atom][0], ions[atom][1], ions[atom][2]});
    }
    writeLammpsDump(ofs, time, frame, box);
    */
  }
  ofs.close();
  std::cout << "Wrote LAMMPS dump to traj.dump\n";




  std::string outputc = "../data/cnnDist/" + dirName + "D" + rho + "E" + fieldname + ".binary"; ; // nearest neighbor distance
  std::cout << "Output generated : " << outputc << " Rows: " << cdist.size() << " Cols: " << cdist[0].size() << "\n";
  printBinary(outputc, cdist);

  std::string outputa = "../data/cnnAngle/" + dirName + "D" + rho + "E" + fieldname + ".binary"; ; // nearest neighbor angle
  std::cout << "Output generated : " << outputa << " Rows: " << cangl.size() << " Cols: " << cangl[0].size() << "\n";
  printBinary(outputa, cangl);


  std::cout << "\n\n";

  return 0;

  /*
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

  for ( int cation=0; cation < numatoms*numtraj ; cation+=2 ) {
    std::vector<float> onecurrent = getColumn(allCurrent, cation, conversion);
    csems.push_back(blockAverage(onecurrent));
    cmeans.push_back(mean(onecurrent));
//    std::cout << "Ion " << cation << " mean: " << mean(onecurrent) << " sem : " << csems.back() << "\n";
  }
  for ( int anion=1; anion < numatoms*numtraj ; anion+=2 ) {
    std::vector<float> onecurrent = getColumn(allCurrent, anion, conversion);
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


  // generate output result
  std::string resfile = "../results/conductivity/condD" + density + "E" + field;
  std::ofstream out (resfile );
  out << std::fixed << std::setprecision(5) << density << "\t\t" << field << "\t\t"  
    << CUTOFFin << "\t\t" << CUTOFFout <<"\t\t" << numsnap*timestep/1000  << "\t\t"
    << cfmean << "\t\t" << cferr << "\t\t" << afmean << "\t\t" << aferr << "\t\t" 
    << (cfmean + afmean)/2 << "\t\t" << std::sqrt( cferr*cferr+aferr*aferr)/2. << "\n";
  out.close();
  //std::string outputcurrent= "./results/condCatCondD" + density + "E" + field +".dat";
  //writeHistogram(outputcurrent, currentargument, currentvec, currentcount);

  return 0;
  */
}

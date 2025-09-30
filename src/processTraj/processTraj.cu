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
//------------------ progress bar ----------------//
#include <chrono>

class ProgressBar {
public:
  ProgressBar(std::size_t total, std::size_t bar_width = 40,
              std::string prefix = "Progress")
      : total_(total), bar_width_(bar_width), prefix_(std::move(prefix)),
        start_(Clock::now()), last_print_(start_) {}

  // Call with 'done' in [0, total], as your loop advances.
  void update(std::size_t done) {
    done = std::min(done, total_);
    const auto now = Clock::now();

    // Throttle printing to ~10 Hz to avoid spamming stdout.
    if (done < total_ &&
        std::chrono::duration_cast<std::chrono::milliseconds>(now - last_print_).count() < 100)
      return;
    last_print_ = now;

    const double frac = total_ ? static_cast<double>(done) / total_ : 1.0;
    const std::size_t filled = static_cast<std::size_t>(std::round(frac * bar_width_));

    // Rate & ETA
    const auto elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(now - start_).count();
    const double rate = elapsed > 0 ? static_cast<double>(done) / elapsed : 0.0; // items/s
    const long long eta_s = (rate > 0.0 && done > 0)
                                ? static_cast<long long>(std::llround((total_ - done) / rate))
                                : -1;

    // Build bar
    std::string bar;
    bar.reserve(bar_width_);
    for (std::size_t i = 0; i < bar_width_; ++i) bar += (i < filled ? '#' : ' '); // or '#'

    std::cout << '\r' << prefix_ << " ["
              << bar << "] "
              << std::setw(3) << static_cast<int>(std::round(frac * 100)) << "%  "
              << done << '/' << total_
              << "  " << std::fixed << std::setprecision(1) << rate << " it/s"
              << "  ETA: " << format_hms(eta_s)
              << std::flush;

    if (done >= total_) std::cout << '\n';
  }

private:
  using Clock = std::chrono::steady_clock;

  std::size_t total_;
  std::size_t bar_width_;
  std::string prefix_;
  Clock::time_point start_, last_print_;

  static std::string format_hms(long long s) {
    if (s < 0) return "--:--:--";
    long long h = s / 3600; s %= 3600;
    long long m = s / 60;   s %= 60;
    std::ostringstream oss;
    oss << std::setw(2) << std::setfill('0') << h << ':'
        << std::setw(2) << std::setfill('0') << m << ':'
        << std::setw(2) << std::setfill('0') << s;
    return oss.str();
  }
};



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
      for( int n=0; n<q; n++) {
        atoms[id2dist[n].first].role = 1;
        if( graph[id2dist[n].first].size() > 1 ) {
          std::swap(atoms[id2dist[n].first].role, atoms[id2dist[q].first].role);
        }
      }
    }
  }
}

void updateTagMerge(const std::vector<std::vector<int>>& clusters, const Rows& pclusters, std::vector<Atom>& atoms, std::vector<int>& proles, int time) {
  auto cls = classify_all(pclusters, clusters);
  for( int icl=0; icl< clusters.size(); icl++ ) {
    if( cls[icl].kind == 4 || cls[icl].kind==1) {
      auto cl = clusters[icl];
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
          //std::swap(atoms[candidates[first]].role, atoms[candidates[second]].role );
          std::swap(atoms[candidates[first]].tag, atoms[candidates[second]].tag );
          //std::cout << "MERGE/CONT " << time << " " << candidates[first] << " " << candidates[second] << "\n";
        }else if (nonzeros > 2) {
          std::cout << "merge, nonzeros > 2\n";
        }
      }
    }
  }
}

void updateTag(const Rows& clusters, const Rows& pclusters, std::vector<Atom>& atoms, std::vector<int>& proles, int time) {
  auto cls = classify_all(pclusters, clusters);
  std::vector<int> parentsIndex; // all parents who dissociate
  for( int cl=0; cl < clusters.size(); cl++ ) {
    if( cls[cl].kind == 2 ) { // if cluster cl experienced dissociation
      for(const auto& pr : cls[cl].parents) {
        parentsIndex.push_back(pr.first);
      }
    }
  }
  if(!parentsIndex.empty()) {
    std::sort(parentsIndex.begin(), parentsIndex.end());
    auto last = std::unique(parentsIndex.begin(), parentsIndex.end());
    parentsIndex.erase(last, parentsIndex.end());
    for( int cl : parentsIndex) {
      std::vector<int> deltaRole;
      for( int candidate : pclusters[cl]) {
        deltaRole.push_back( atoms[candidate].role - proles[candidate]);
      }
      int nonzeros = std::count_if(deltaRole.begin(), deltaRole.end(), [](int x) { return x != 0; } );
      //std::cout << "Time " << time <<" dissociation with # nonzeros " << nonzeros << "\n";
      if (nonzeros >= 2) {
        std::vector<int> neg, pos;
        for ( int i = 0 ; i < deltaRole.size(); i++ ){
          if(deltaRole[i] == -1) neg.push_back(pclusters[cl][i]);
          else if (deltaRole[i] == 1) pos.push_back(pclusters[cl][i]);
        }

        int n = std::min(neg.size(), pos.size());
        for(int k=0; k<n; k++){
          std::swap( atoms[neg[k]].tag, atoms[pos[k]].tag);
          //std::cout << "DISS " << time << " " << neg[k] << " " << pos[k] << "\n";
        }
      } else if ( nonzeros >=3 ) {
        std::cout << "parents has more nonzeros\n";
      }
    }
  }

  std::vector<int> ingredients;
  for( int cl=0; cl < clusters.size(); cl++ ) {
    if( cls[cl].kind == 3 ) { // if cluster cl experienced exchange
      for( int atom : clusters[cl]) {
        ingredients.push_back(atom);
      }
    }
  }
  if(!ingredients.empty()) {
//    std::cout << "EXCH " << time << " ";
//    for( int elem : ingredients) { std::cout << elem << " ";}
//    std::cout << "\n";
    std::vector<int> deltaRole;
    for( int candidate : ingredients) {
      deltaRole.push_back( atoms[candidate].role - proles[candidate]);
    }
    int nonzeros = std::count_if(deltaRole.begin(), deltaRole.end(), [](int x) { return x != 0; } );
    //std::cout << "Time " << time << " exchange with # nonzeros " << nonzeros << "\n";
    if (nonzeros >= 2) {
      std::vector<int> neg, pos;
      for ( int i = 0 ; i < deltaRole.size(); i++ ){
        if(deltaRole[i] == -1) neg.push_back(ingredients[i]);
        else if (deltaRole[i] == 1) pos.push_back(ingredients[i]);
      }

      int n = std::min(neg.size(), pos.size());
      for(int k=0; k<n; k++){
        std::swap( atoms[neg[k]].tag, atoms[pos[k]].tag);
  //      std::cout << "EXCH " << time << " " << neg[k] << " " << pos[k] << "\n";
      }
    } else if ( nonzeros >=3 ) {
      std::cout << "ingredients has more nonzeros\n";
    }
  }
}

struct indexTag {
  int tag;
  int pid;
  int cid;
};

void detectJump( const std::vector<Atom>& atoms, const std::vector<float>& pdist ,const std::vector<int> ptags, float CUTOFFin, float CUTOFFout, int time) {
  std::vector<indexTag> idx;
  for ( int tag =0; tag<pdist.size(); tag++){  
    float prevDist = 0;
    int index;
    auto it = std::find(ptags.begin(), ptags.end(), tag);
    if ( it == ptags.end() ) {
      std::cout << "Previous, cannot find properly tagged particle\n\n";
    } else {
      index = std::distance(ptags.begin(), it);
      prevDist = pdist[index];
    }

    float currDist = 0;
    auto cit = std::find_if( atoms.begin(), atoms.end(), [tag](const Atom& a ){ return a.tag == tag; });
    if (cit == atoms.end()) { 
      std::cout << "Current  cannot find properly tagged particle\n\n";
    } else {
      currDist = cit->nndist;
    }
      
    if (std::fabs(currDist - prevDist) > (CUTOFFout - CUTOFFin) ) {
      idx.emplace_back( indexTag{tag, index, cit->id});
    }
  }
  // idx contains tags of jumped particles
  if(!idx.empty()){ 
    for (auto& mem : idx) {
      if ( pdist[mem.pid] > CUTOFFout && atoms[mem.cid].nndist < CUTOFFin ) {
        std::cout << "Jump from out to in | time " << time << " | pid " << mem.pid << " pdist " << pdist[mem.pid] << " cid " << mem.cid << " cdist " << atoms[mem.cid].nndist << "\n";
      }
      if ( pdist[mem.pid] < CUTOFFin && atoms[mem.cid].nndist > CUTOFFout ) {
        std::cout << "Jump from in to out | time " << time << "  pid " << mem.pid << " pdist " << pdist[mem.pid] << " cid " << mem.cid << " cdist " << atoms[mem.cid].nndist << "\n";
      }
    }
  }
}


////////// Statistics
struct Histogram {
  std::vector<float> edges;
  std::vector<float> counts;
  float binWidth(std::size_t i=0) const {
    return edges.size() > 1? (edges[i+1]-edges[i]) : 0.0;
  }
};

Histogram makeHistogram(const std::vector<float>& data, size_t nbins, float minEdge, float maxEdge){
  Histogram h;
  h.edges.resize(nbins + 1);
  h.counts.assign(nbins, 0.0);

  if (nbins == 0) return h;
  if (maxEdge == minEdge) maxEdge = minEdge + 1;

  const float width = (maxEdge - minEdge)/ static_cast<float>(nbins);
  for( size_t i =0; i<= nbins; i++ ) h.edges[i] = minEdge + i* width;

  for(size_t k=0; k<data.size(); k++) {
    float x = data[k];
    long idx = static_cast<long>(std::floor( (x-minEdge) / width));
    if (idx < 0 ) idx =0;
    if (idx > static_cast<long>(nbins)) idx = static_cast<long>(nbins)-1;
    h.counts[static_cast<size_t>(idx)] += 1.0;
  }
  return h;
}

Histogram makeHistogramAuto(const std::vector<float>& data, size_t nbins) {
  if ( data.empty() || nbins == 0) return {};
  auto mm = std::minmax_element(data.begin(), data.end());
  float mn = *mm.first;
  float mx = *mm.second;
  if ( mx == mn ) {mn -= 0.5; mx += 0.5;}

  return makeHistogram(data, nbins, mn, mx);
}
std::vector<float> normalizePDF(const Histogram& h) {
  float N = std::accumulate(h.counts.begin(), h.counts.end(), 0.0);
  std::vector<float> pdf(h.counts.size(), 0.0);
  if ( N > 0 && h.edges.size() >=2) {
    for (size_t i =0; i< h.counts.size(); i++) {
      float bw = h.edges[i+1] - h.edges[i];
      if (bw > 0 ) pdf[i] = h.counts[i] / (N*bw);
    }
  }
  return pdf;
}

void printHistogram( const Histogram& h, std::ostream& os = std::cout) {
  std::vector<float> pdf = normalizePDF(h);
  os.setf(std::ios::fixed);
  os << std::setprecision(10);
  for( size_t i =0; i< h.edges.size(); i++ ) {
    float x = 0.5 * (h.edges[i] + h.edges[i+1]);
    os << x << "\t" << pdf[i] << "\n";
  }
}

std::vector<float> histPDFatEdges(const Histogram& h) {
  size_t nb = h.counts.size();
  std::vector<float> y; 
  y.assign(nb+1, 0.0);

  float N = std::accumulate(h.counts.begin(), h.counts.end(), 0.0);
  if ( N<=0.0 || h.edges.size() != nb+1 ) return y;

  for (size_t i =0; i< nb ; i++){
    float bw = h.edges[i+1] - h.edges[i];
    y[i] = (bw > 0.0) ? (h.counts[i] / N / bw) : 0.0;
  }
  y[nb] = y[nb - (nb > 0 ? 1: 0)];
  return y;
}

float linInterp(float x0, float y0, float x1, float y1, float x) {
  float dx = x1-x0;
  if (dx==0.0) return 0.5*(y0+y1);
  float t= (x-x0)/dx;
  return y0 + t* (y1-y0);
}

float integratePDF(const Histogram& h, float xi, float xf) {
  if (h.edges.size() < 2 || h.counts.size() +1 != h.edges.size()) return 0.0;
  if (xi > xf ) std::swap(xi, xf);

  const std::vector<float>& x = h.edges;
  std::vector<float> y = histPDFatEdges(h);

  float a = x.front();
  float b = x.back();
  if ( xf <= a || xi >= b) return 0.0;
  xi = std::max(xi, a);
  xf = std::min(xf, b);

  std::vector<float>::const_iterator itL = std::upper_bound(x.begin(), x.end(), xi);
  std::vector<float>::const_iterator itR = std::upper_bound(x.begin(), x.end(), xf);
  size_t j = (itL == x.begin() ? 0 : static_cast<size_t>((itL - x.begin()) -1));
  size_t k = (itR == x.begin() ? 0 : static_cast<size_t>((itR - x.begin()) -1));
  if( j >= x.size() -1) j = x.size() -2;
  if( k >= x.size() -1) k= x.size() -2;
  float y_xi = linInterp(x[j], y[j], x[j+1], y[j+1], xi);
  float y_xf = linInterp(x[k], y[k], x[k+1], y[k+1], xf);
  if ( j== k) { 
    return 0.5 * (y_xi + y_xf) * (xf - xi);
  }

  float area = 0.0;
  area += 0.5 * (y_xi + y[j+1]) * (x[j+1] - xi);

  for ( size_t m = j+1; m< k; m++ ) {
    area += 0.5 * (y[m]+y[m+1]) * (x[m+1] - x[m]);
  }
  area += 0.5 * (y[k] + y_xf) * (xf - x[k]);
  return area;
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
  bool videoflag=false;
  int thermoflag=0;
  if (argc < 11 ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << "dir_name density field cutoff_in(A) cutoff_out(A) boxX(A) boxY(A) boxZ(A) timestep(ps)  eqtime(ns) \n";
    return 1;
  } else if (argc>11) {
    videoflag=static_cast<bool>(std::stoi(argv[11])); // if 1: write video upto first 20000 steps
    thermoflag=static_cast<int>(std::stoi(argv[12]));// if 1: use partial thermostat dataset, starting with Tdump..
  } else {
    std::cerr << "Error: Too much  arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << "dir_name density field cutoff_in(A) cutoff_out(A) boxX(A) boxY(A) boxZ(A) timestep(ps)  eqtime(ns) \n";
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


  std::cout << "Dump Reading Starts\n";
  int numsnap;
  std::vector<std::vector<float>> traj;
  std::string inputFileName;
  if (thermoflag == 1) { inputFileName = std::string("../data/dumps") + dirName + "/TdumpD" + rho + "E" + fieldname + ".binary";}
  else if (thermoflag == 2) { inputFileName = std::string("../data/dumps") + dirName + "/TTdumpD" + rho + "E" + fieldname + ".binary";}
  else {inputFileName = std::string("../data/dumps") + dirName + "/dumpD" + rho + "E" + fieldname + ".binary";}
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

  //collect all time series of reduced trajectory
  std::vector<std::vector<float>>  cdist(numsnap, std::vector<float>( numatoms, 0.0) ); // conditioned nearest neighbor distance, A unit
  std::vector<std::vector<float>>  cangl(numsnap, std::vector<float>( numatoms, 0.0) ); // conditioned nearest neighbor angle, from positive z axis, degree unit

  // placeholder for previous frame Atom properties
  std::vector<int> proles(numatoms, -1);
  std::vector<int> ptags(numatoms);
  std::iota(ptags.begin(), ptags.end(), 0);
  std::vector<float> pdist(numatoms, -1);
  std::vector<std::vector<int>> pclusters;

  std::vector<float> deltaDist;

  // save visualization for each time step
  std::ofstream ofs;
  if (videoflag)  {
    std::cout << "\n\n-------- generate initial 20000 step video ---------------------\n\n";
    ofs.open("video.lammpstrj");
    if (!ofs) {
      std::cerr << "Error. Cannot open the file\n";
      return 1;
    }
  }
  
  int st=0;
  //int ed=143600;
  int ed=numsnap;
  if (videoflag) ed=20000;
  ProgressBar pb(numsnap, /*bar_width=*/50, "Analyzing");
  int every=numsnap/100;
  
  //for (int time= st; time<numsnap; time++){
  for (int time= st; time<ed; time++){
    std::vector<Atom> frame; // id charge role swap tag nndist nnangl nnanion x y z
    //get one frame
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
      for(auto& atom : frame) atom.tag = atom.id; // initialize tag
      roleAssignment(clusters, frame, graph, distanceMatrix, box); // frame[i].role determines whether it finds counter ion inside or outside of cluster
      findCounterIons(frame, clusters, distanceMatrix, box); // based on role, find conditioned counter ion (index, dist, angl)

    } else {
      for(auto& atom : frame) atom.tag = ptags[atom.id]; // update tag
      roleAssignment(clusters, frame, graph, distanceMatrix, box);
      findCounterIons(frame, clusters, distanceMatrix, box);
      updateTagMerge(clusters,pclusters,  frame,  proles, time); // merge or continue
      updateTag(clusters, pclusters, frame, proles, time);  // dissociate or exchange

      detectJump(frame, pdist, ptags, CUTOFFin, CUTOFFout, time); // list of tags

      for(int tagg=0; tagg<numatoms; tagg++) {
        auto it = std::find_if(frame.begin(), frame.end(), [tagg](const Atom& a){return a.tag==tagg;});
        cdist[time][tagg] = it->nndist;
        cangl[time][tagg] = it->nnangl;
      }

      for( int tag =0; tag < numatoms; tag++ ) {
        int pid, cid;
        auto it = std::find( ptags.begin(), ptags.end(), tag);
        if (it != ptags.end()) {
          pid = std::distance(ptags.begin(), it);
        }
        auto cit = std::find_if(frame.begin(), frame.end(), [tag](const Atom& a) {return a.tag == tag;});
        if (cit != frame.end()){
          cid = cit->id;
        }
        deltaDist.push_back(cit->nndist - pdist[pid] );
      }
    }
    // save role, tag and cluster info 
    std::transform(frame.begin(), frame.end(), proles.begin(),[](const Atom& a){ return a.role;}); // save roles
    std::transform(frame.begin(), frame.end(), ptags.begin(),[](const Atom& a){ return a.tag;}); // save tags
    std::transform(frame.begin(), frame.end(), pdist.begin(),[](const Atom& a){ return a.nndist;}); // save prev nn dist
    pclusters = clusters; // save cluster info
    if(videoflag) writeLammpsDump(ofs, time, frame, box);
    if( time%every == 0 ) pb.update(time+1);
  }
  if(videoflag) ofs.close();
  //std::cout << "Wrote LAMMPS dump to traj.dump\n";

  /*
  auto h = makeHistogramAuto( deltaDist, 500);
  std::ofstream out("displacementHist.dat");
  printHistogram(h, out);
  */

  //float cutone = integratePDF(h, -1, 1);
  //std::cout << "-1 to 1 : " << cutone << std::endl;

  std::string outputc;
  if (thermoflag==1) { outputc =  "../data/cnnDist/" + dirName + "TD" + rho + "E" + fieldname + ".binary";}
  else if (thermoflag==2) { outputc =  "../data/cnnDist/" + dirName + "TTD" + rho + "E" + fieldname + ".binary";}
  else { outputc = "../data/cnnDist/" + dirName + "D" + rho + "E" + fieldname + ".binary";}
  std::cout << "Output generated : " << outputc << " Rows: " << cdist.size() << " Cols: " << cdist[0].size() << "\n";
  printBinary(outputc, cdist);
  /*
  std::string text = "./cnnDist.dat";
  printMatrix(text, cdist);
  */

  std::string outputa;
  if( thermoflag==1) { outputa = "../data/cnnAngle/" + dirName + "TD" + rho + "E" + fieldname + ".binary"; }
  else if( thermoflag==2) { outputa = "../data/cnnAngle/" + dirName + "TTD" + rho + "E" + fieldname + ".binary"; }
  else { outputa = "../data/cnnAngle/" + dirName + "D" + rho + "E" + fieldname + ".binary"; }
  std::cout << "Output generated : " << outputa << " Rows: " << cangl.size() << " Cols: " << cangl[0].size() << "\n";
  printBinary(outputa, cangl);

  std::cout << "\n\n";

  return 0;
}

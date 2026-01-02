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

//linear regression
#include <cstddef>
#include <limits>
#include <stdexcept>

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

float wrapPBC(float x, float bx){
  return x - bx*std::floor(x/bx);
}

void checkWrapping( float x, float bx) {
  if ( x < 0 || x > bx ) {
    std::cout << "wrapping pbc image doesn't work properly\n\n";
  }
}

std::vector<float> constructImageIndex( std::vector<float>& z, int numsnap, int numatoms, float bz) {
  std::vector<float> iz(numsnap * numatoms, 0.0);
  for(size_t time=1; time < numsnap; time++ ) {
    for( size_t atom=0; atom < numatoms; atom++) {
      int index =  time * numatoms + atom;
      if( (z[index] - z[index-numatoms]) >  bz/2.0 ) {
        iz[index] = iz[index-numatoms] -1;
      } else if ( (z[index] -z[index-numatoms]) < -bz/2.0) {
        iz[index] = iz[index-numatoms] +1;
      } else {
        iz[index] = iz[index-numatoms];
      }
    }
  }
  return iz;
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

float trapz(const std::vector<float>& x, const std::vector<float>& y) {
  size_t n = x.size();
  float sum = 0;
  for( size_t i=0; i<n-1; i++ ) {
    float dx = x[i+1] - x[i];
    sum += (y[i] + y[i+1]) * dx * 0.5;
  }
  return sum;
}


struct LinRegResult {
  float slope;      // b
  float intercept;  // a
  float r2;         // coefficient of determination
  float sse;        // sum of squared residuals
  float sst;        // total sum of squares (about mean of y)
  std::size_t n;    // number of points
};

inline LinRegResult linear_fit(const std::vector<float>& x,
                               const std::vector<float>& y)
{
  if (x.size() != y.size())
    throw std::invalid_argument("x and y must have the same length");
  const std::size_t n = x.size();
  if (n < 2)
    throw std::invalid_argument("need at least 2 points");

  // Accumulate sums
  long double Sx = 0.0L, Sy = 0.0L, Sxx = 0.0L, Sxy = 0.0L, Syy = 0.0L;
  for (std::size_t i = 0; i < n; ++i) {
    const long double xi = x[i];
    const long double yi = y[i];
    Sx  += xi;
    Sy  += yi;
    Sxx += xi * xi;
    Sxy += xi * yi;
    Syy += yi * yi;
  }

  const long double nL = static_cast<long double>(n);
  const long double denom_x = nL * Sxx - Sx * Sx;
  const long double denom_y = nL * Syy - Sy * Sy;

  if (std::abs(denom_x) == 0.0L) {
    // All x are identical -> vertical line; slope undefined
    return {std::numeric_limits<float>::quiet_NaN(),
            std::numeric_limits<float>::quiet_NaN(),
            std::numeric_limits<float>::quiet_NaN(),
            std::numeric_limits<float>::quiet_NaN(),
            std::numeric_limits<float>::quiet_NaN(),
            n};
  }

  const long double b = (nL * Sxy - Sx * Sy) / denom_x;   // slope
  const long double a = (Sy - b * Sx) / nL;               // intercept

  // Compute SSE and SST to report R^2 via 1 - SSE/SST
  const long double ybar = Sy / nL;
  long double SSE = 0.0L; // residual sum of squares
  long double SST = 0.0L; // total sum of squares
  for (std::size_t i = 0; i < n; ++i) {
    const long double yi = y[i];
    const long double fi = a + b * x[i];
    const long double ri = yi - fi;
    SSE += ri * ri;
    const long double dy = yi - ybar;
    SST += dy * dy;
  }

  float r2;
  if (SST > 0.0L) {
    r2 = static_cast<float>(1.0L - SSE / SST);
  } else {
    // All y identical: define R^2 = 1 if SSE == 0 (perfect fit), else 0
    r2 = (SSE == 0.0L) ? 1.0f : 0.0f;
  }

  return {static_cast<float>(b),
          static_cast<float>(a),
          r2,
          static_cast<float>(SSE),
          static_cast<float>(SST),
          n};
}



float GreenKubo( std::vector<float>& unwrapped, std::vector<float>& charges, float volume, int dt, size_t numsnap, int numatoms, float timestep) {
  //dt : number of timestep to be used for velocity estimate
  std::vector<float> Jz; Jz.reserve(numsnap);
  for( size_t time=dt; time<numsnap-dt; time++) {
    float current = 0;
    for ( int atom=0; atom< numatoms; atom++) {
      int index = numatoms * time + atom;
      float disp = unwrapped[index + dt*numatoms] - unwrapped[index - dt*numatoms]; 
      current += charges[index] * disp/2.0/dt/timestep ; // C A/ps
    }
    Jz.push_back(current);
  }

  int maxTime=2000; // 100 timesteps = 5ps
  std::vector<float> interval;
  std::vector<float> ACF;
  for( int intv = 0; intv < maxTime; intv++){
    interval.push_back(intv*timestep);
    float mean;
    for( size_t time=0; time < Jz.size()-maxTime-1; time++ ) {
      float val = Jz[time+intv] * Jz[time];
      mean += (val - mean) / (time+1);
    }
    ACF.push_back(mean);
  }
  //for(int i =0; i< ACF.size(); i++) std::cout << "t " << interval[i] << " ACF " << ACF[i] << "\n";
  
  float kBT = 0.592*4184; // J/mol
  float numpair = numatoms/2;
  float sigma =  trapz(interval, ACF) / kBT / numpair;
  sigma = sigma *96485 * 96485 / 10000; // unit conversion e^2 A^2/ps

  return sigma;
}

Box EinsteinHelfand( std::vector<float>& unwrapped, std::vector<float>& charges, float volume, int dt, size_t numsnap, int numatoms, float timestep) {
  //dt : number of timestep to be used for velocity estimate

  std::vector<int> interval;
  for(int intv=400; intv < 700; intv+=2 ) {
    interval.push_back(intv);
  }

  std::vector<float> meansqG;
  for( int intv : interval) {
    float Gsq=0; // time integrated current, 
    for ( size_t time=intv; time < numsnap - intv; time++ ) {
      float frameG = 0;
      for( int atom=0; atom < numatoms; atom++ ) {
        int index = numatoms*time + atom;
        float disp = unwrapped[index +intv*numatoms] - unwrapped[index];
        frameG += charges[index] * disp;
      }
      Gsq+= ( frameG*frameG - Gsq) / (time-intv+1);
    }
    meansqG.push_back(Gsq);
  }
  std::vector<float> time;
  for(int i=0; i<meansqG.size(); i++) {
    //std::cout << " x " << interval[i] << " y " << meansqG[i] << "\n";
    time.push_back( interval[i] * timestep);
  }

  LinRegResult res = linear_fit(time, meansqG);

  /*
  std::cout << "\nLinear Regression Summary\n";
  std::cout << "slope b      = " << res.slope << "\n";
  std::cout << "intercept a  = " << res.intercept << "\n";
  std::cout << "R^2          = " << res.r2 << "\n";
  std::cout << "SSE          = " << res.sse << "\n";
  std::cout << "n            = " << res.n << "\n";
  */

  float kBT = 0.592*4184; // J/mol
  float numpair=numatoms/2;
  float sigma = res.slope /2/kBT/numpair;
  sigma = sigma *96485 * 96485 / 10000;

  // computing standard error of the slope
  float SE = res.slope * std::sqrt( (1-res.r2)/(res.n-2) / res.r2);
  SE = SE/2/kBT/numpair* 96485 * 96485 / 10000;
  Box out;
  out.x = sigma;
  out.y = SE;

  return out;
}

Box partialEinsteinHelfand( std::vector<float>& unwrapped, std::vector<float>& charges, float volume, int dt, size_t numsnap, int numatoms, float timestep, int ensemble) {
  //dt : number of timestep to be used for velocity estimate

  std::vector<int> interval;
  for(int intv=400; intv < 700; intv+=2 ) {
    interval.push_back(intv);
  }
  int chunck = numsnap/3;

  std::cout << "\n chunck " << chunck << " ensemble " << ensemble << "\n";

  std::vector<float> meansqG;
  for( int intv : interval) {
    float Gsq=0; // time integrated current, 
    
    for ( size_t time=ensemble*chunck + intv ; time < (ensemble+1)*chunck - intv; time++ ) {
      float frameG = 0;
      for( int atom=0; atom < numatoms; atom++ ) {
        int index = numatoms*time + atom;
        float disp = unwrapped[index +intv*numatoms] - unwrapped[index];
        frameG += charges[index] * disp;
      }
      Gsq+= ( frameG*frameG - Gsq) / (time-(intv+ensemble*chunck)+1);
    }
    meansqG.push_back(Gsq);
  }
  std::vector<float> time;
  for(int i=0; i<meansqG.size(); i++) {
  //  std::cout << " x " << interval[i] << " y " << meansqG[i] << "\n";
    time.push_back( interval[i] * timestep);
  }

  LinRegResult res = linear_fit(time, meansqG);

  /*
  std::cout << "\nLinear Regression Summary\n";
  std::cout << "slope b      = " << res.slope << "\n";
  std::cout << "intercept a  = " << res.intercept << "\n";
  std::cout << "R^2          = " << res.r2 << "\n";
  std::cout << "SSE          = " << res.sse << "\n";
  std::cout << "n            = " << res.n << "\n";
  */

  float kBT = 0.592*4184; // J/mol
  float numpair=numatoms/2;
  float sigma = res.slope /2/kBT/numpair;
  sigma = sigma *96485 * 96485 / 10000;

  // computing standard error of the slope
  float SE = res.slope * std::sqrt( (1-res.r2)/(res.n-2) / res.r2);
  SE = SE/2/kBT/numpair* 96485 * 96485 / 10000;
  Box out;
  out.x = sigma;
  out.y = SE;

  return out;
}

float NernstEinstein( std::vector<float>& unwrapped, std::vector<float>& charges, float volume, size_t numsnap, int numatoms, float timestep){
  std::vector<int> interval;
  for(int intv=40; intv < 400; intv++ ) {
    interval.push_back(intv);
  }

  std::vector<float> meansqDc;
  std::vector<float> meansqDa;
  for( int intv : interval) {
    float Dsqc=0, Dsqa=0; // displacement square for cation and anion
    float countc=0, counta=0;
    for ( size_t time=intv; time < numsnap - intv; time++ ) {
      for( int atom=0; atom < numatoms; atom++ ) {
        int index = numatoms*time + atom;
        float disp = unwrapped[index +intv*numatoms] - unwrapped[index];
        if ( charges[index] == 1 ) {
          Dsqc += (disp * disp);
          countc += 1;
        } else {
          Dsqa += (disp*disp);
          counta += 1;
        }
      }
    }
    meansqDc.push_back(Dsqc/countc);
    meansqDa.push_back(Dsqa/counta);
  }
  std::vector<float> time;
  for(int i=0; i<meansqDc.size(); i++) {
    time.push_back( interval[i] * timestep);
  }

  LinRegResult resc = linear_fit(time, meansqDc);
  LinRegResult resa = linear_fit(time, meansqDa);

  std::cout << "\nLinear Regression Summary\n";
  std::cout << "Cation\n";
  std::cout << "slope b      = " << resc.slope << "\n";
  std::cout << "R^2          = " << resc.r2 << "\n";
  std::cout << "Cation\n";
  std::cout << "slope b      = " << resa.slope << "\n";
  std::cout << "R^2          = " << resa.r2 << "\n";

  float Dc = resc.slope ; // A^2/ps
  float Da = resa.slope ; // A^2/ps

  Dc*= 10;// nm^2/ns
  Da*= 10;// nm^2/ns

  std::cout << "Diffusion Coefficient nm^2/ns: Dc " << Dc  << " Da " << Da  << "\n";

  float RT=8.31446*300; // 
  float lambda = (Dc+Da)* (96485.0*96485.0) / RT /100000 ;
  return lambda;

}
  
Box Response( std::vector<float>& unwrapped, std::vector<float>& charges, float volume, size_t numsnap, int numatoms, float timestep, float fieldStrength) {
  //dt : number of timestep to be used for velocity estimate
  int dt=1;
  std::vector<float> Jz; Jz.reserve(numsnap);
  for( size_t time=dt; time<numsnap-dt; time++) {
    float current = 0;
    for ( int atom=0; atom< numatoms; atom++) {
      int index = numatoms * time + atom;
      float disp = unwrapped[index + dt*numatoms] - unwrapped[index - dt*numatoms]; 
      current += charges[index] * disp/2.0/dt/timestep ; // C A/ps
    }
    Jz.push_back(current);
  }

  float meanFlux = mean(Jz);
  float stdFlux = std::sqrt(variance(Jz))/std::sqrt(Jz.size()); 

  float kBT = 0.592*4184; // J/mol
  float numpair = numatoms/2;
  meanFlux = meanFlux * 96485 * 96485 /kBT/ (1.0/25.0) / 10000 / fieldStrength / numpair;
  stdFlux = stdFlux * 96485 * 96485 /kBT/ (1.0/25.0) / 10000 / fieldStrength / numpair;
  Box out;
  out.x = meanFlux;
  out.y = stdFlux;

  return out;
}






  

///////////End statistics

int main(int argc, char* argv[]) {
  int thermoflag=0;
  if (argc < 11 ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << "dir_name density field cutoff_in(A) cutoff_out(A) boxX(A) boxY(A) boxZ(A) timestep(ps)  eqtime(ns) \n";
    return 1;
  } else if (argc>11) {
    thermoflag=static_cast<int>(std::stoi(argv[11]));// if 1: use partial thermostat dataset, starting with Tdump..
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
  

  std::cout << "\nConductivity Measruement\n";
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

  std::vector<float> deltaDist;

  
  int st=0;
  int ed=numsnap;
  ProgressBar pb(numsnap, /*bar_width=*/50, "Analyzing");
  int every=numsnap/100;
  
  std::vector<float> unwrappedX(numsnap * numatoms, 0.0);
  std::vector<float> wrappedX(numsnap * numatoms, 0.0);
  std::vector<float> imageX;

  std::vector<float> unwrappedY(numsnap * numatoms, 0.0);
  std::vector<float> wrappedY(numsnap * numatoms, 0.0);
  std::vector<float> imageY;

  std::vector<float> unwrappedZ(numsnap * numatoms, 0.0);
  std::vector<float> wrappedZ(numsnap * numatoms, 0.0);
  std::vector<float> imageZ;
  std::vector<float> charges(numsnap* numatoms, 0.0);

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

    // wrapPBC first before assigning image index
    
    for( int atom=0; atom < numatoms; atom++) {
      int index = time*numatoms + atom;
      wrappedX[index] = wrapPBC(frame[atom].x, box.x);
      wrappedY[index] = wrapPBC(frame[atom].y, box.y);
      wrappedZ[index] = wrapPBC(frame[atom].z, box.z);
      //checkWrapping(wrappedZ[index], box.z);
      charges[index]=frame[atom].charge;
    }

    if( time%every == 0 ) pb.update(time+1);
  }

  imageX = constructImageIndex(wrappedX,numsnap, numatoms, box.x);
  imageY = constructImageIndex(wrappedY,numsnap, numatoms, box.y);
  imageZ = constructImageIndex(wrappedZ,numsnap, numatoms, box.z);
  for(size_t element=0; element < imageZ.size(); element++ ) {
    unwrappedX[element] = wrappedX[element] + imageX[element]* box.x;
    unwrappedY[element] = wrappedY[element] + imageY[element]* box.y;
    unwrappedZ[element] = wrappedZ[element] + imageZ[element]* box.z;
  }

  //Box conductivityGK; // Green-Kubo
  //Box conductivityEH; // Einstein-Helfand
  //Box conductivityNE; // Nernst-Einstein

  Box condEHx, condEHy, condEHz;
  float volume = box.x * box.y * box.z; // A^3

  std::cout << std::fixed << std::setprecision(6);

  float molarConductivity=0;
  float molarConductivityE=0;

  if ( field==0) {
  //if (1){
    /*
    Box conductivityGK;
    conductivityGK.x = GreenKubo(unwrappedX, charges, volume, 1,  numsnap, numatoms, timestep);
    conductivityGK.y = GreenKubo(unwrappedY, charges, volume, 1,  numsnap, numatoms, timestep);
    conductivityGK.z = GreenKubo(unwrappedZ, charges, volume, 1,  numsnap, numatoms, timestep);
    std::cout << "\nGreen-Kubo conductivity\n";
    //std::cout << "conductivity : " << conductivityGK << "S / m\n";
    std::cout << "molar conductivity x " << conductivityGK.x   << " S cm^2/mol\n";
    std::cout << "molar conductivity y " << conductivityGK.y   << " S cm^2/mol\n";
    std::cout << "molar conductivity z " << conductivityGK.z   << " S cm^2/mol\n";

    

    condEHx = EinsteinHelfand(unwrappedX, charges, volume, 1, numsnap, numatoms, timestep);
    condEHy = EinsteinHelfand(unwrappedY, charges, volume, 1, numsnap, numatoms, timestep);
    condEHz = EinsteinHelfand(unwrappedZ, charges, volume, 1, numsnap, numatoms, timestep);
    std::cout << "\nEinstein-Helfand conductivity\n";
    std::cout << "molar conductivity x  " << condEHx.x   << " +- " << condEHx.y << " " << "S cm^2/mol\n";
    std::cout << "molar conductivity y  " << condEHy.x   << " +- " << condEHy.y << " " << "S cm^2/mol\n";
    std::cout << "molar conductivity z  " << condEHz.x   << " +- " << condEHz.y << " " << "S cm^2/mol\n";

    //float wx = condEHx.y * condEHx.y;
    //float wy = condEHy.y * condEHy.y;
    //float wz = condEHz.y * condEHz.y;

    //molarConductivity = (condEHx.x / wx + condEHy.x / wy + condEHz.x / wz) / (1/wx + 1/wy + 1/wz);
    //molarConductivityE=std::sqrt( 1/ ( 1/wx + 1/wy + 1/wz));
    molarConductivity = (condEHx.x + condEHy.x + condEHz.x )/3.0;
    float sampleSD = std::sqrt( std::pow(condEHx.x-molarConductivity, 2) +  std::pow(condEHy.x-molarConductivity, 2) +  std::pow(condEHz.x-molarConductivity, 2))/std::sqrt(2);
    molarConductivityE = sampleSD/std::sqrt(3);

    std::cout << "Total mean : " << molarConductivity << " +- " << molarConductivityE << "\n";
    */


    std::vector<float> values;
    for(int ensemble=0; ensemble<3; ensemble++) {
      condEHx = partialEinsteinHelfand(unwrappedX, charges, volume, 1, numsnap, numatoms, timestep, ensemble);
      condEHy = partialEinsteinHelfand(unwrappedY, charges, volume, 1, numsnap, numatoms, timestep, ensemble);
      condEHz = partialEinsteinHelfand(unwrappedZ, charges, volume, 1, numsnap, numatoms, timestep, ensemble);
      values.push_back(condEHx.x);
      values.push_back(condEHy.x);
      values.push_back(condEHz.x);
      std::cout << "values " << condEHx.x << " " << condEHy.x << " " << condEHz.x << "\n";
    }


    std::cout << "split3 Total mean : " << mean(values) << " +- " << std::sqrt( variance(values)/ static_cast<float>(values.size()) ) << "\n";




    /*
    Box conductivityNE;
    conductivityNE.x = NernstEinstein(unwrappedX, charges, volume, numsnap, numatoms, timestep);
    conductivityNE.y = NernstEinstein(unwrappedY, charges, volume, numsnap, numatoms, timestep);
    conductivityNE.z = NernstEinstein(unwrappedZ, charges, volume, numsnap, numatoms, timestep);
    std::cout << "\nNernst-Einstein conductivity\n";
    std::cout << "molar conductivity x " << conductivityNE.x   << " S cm^2/mol\n";
    std::cout << "molar conductivity y " << conductivityNE.y   << " S cm^2/mol\n";
    std::cout << "molar conductivity z " << conductivityNE.z   << " S cm^2/mol\n";
    */
  }
  else {
    Box response; // <J>/E, not differential conductivity
    Box std;
    response.x = Response(unwrappedX, charges, volume, numsnap, numatoms, timestep, field*fconversion).x;
    response.y = Response(unwrappedY, charges, volume, numsnap, numatoms, timestep, field*fconversion).x;
    response.z = Response(unwrappedZ, charges, volume, numsnap, numatoms, timestep, field*fconversion).x;
    std.x = Response(unwrappedX, charges, volume, numsnap, numatoms, timestep, field*fconversion).y;
    std.y = Response(unwrappedY, charges, volume, numsnap, numatoms, timestep, field*fconversion).y;
    std.z = Response(unwrappedZ, charges, volume, numsnap, numatoms, timestep, field*fconversion).y;
    std::cout << "\nResponse x " << response.x   << " +- " << std.x << "  S cm^2/mol\n";
    std::cout << "Response y " << response.y   <<   " +- " << std.y << " S cm^2/mol\n";
    std::cout << "Response z " << response.z   <<   " +- " << std.z << " S cm^2/mol\n";

    molarConductivity = response.z;
    molarConductivityE = std.z;
  }


  /*
  for(size_t time=0; time<1000; time++) {
    int index = time*numatoms +1;
    std::cout << "raw " << unwrappedZ[index] << " wrapped " << wrappedZ[index] << " imageID " << imageZ[index] << "\n";
  }
  */

  /*
  float conductivityEH=0; // Einstein-Helfand
  float flux=0; // current flux for differential conductivity

  if ( field==0) {
    conductivityGK = GreenKubo();
    conductivityEH = EinsteinHelfand();
  }

  flux = measureFlux();

  if (field==0) {
    std::cout << "\nZero Field Conductivity Results\n" << "Green-Kubo : " << conductivityGK << "mS/A" ...


    */


  std::cout << "\n\n";

  // generate output result
  std::string rateFile;
  rateFile=std::string("../results/conductivity/") + dirName + "D" + rho + "E" + fieldname + ".dat";
  std::ofstream out(rateFile );
  out << std::fixed << std::setprecision(8) << density << "\t\t" << field << "\t\t"  
    << CUTOFFin << "\t\t" << CUTOFFout <<"\t\t"  
    << molarConductivity << "\t\t" << molarConductivityE 
    << "\n";
  out.close();
  return 0;
}

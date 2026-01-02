#include <ioTraj/atom.hpp>
#include <ioTraj/dcdreader.hpp>
#include <ioTraj/lammpsreader.hpp>

using namespace ioTraj;

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include<algorithm>
#include <iomanip>
#include <cmath>
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
using Rows = std::vector<std::vector<int>>;

struct Dipole {
  float x, y, z, mu; // mu=amplitude
};
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


////////////    IO
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
void writeLammpsDumpLiPF6ACN(std::ostream& out, long timestep, const std::vector<Atom>& atoms, const std::vector<int> atomTypes, const Box& box) {
  out << "ITEM: TIMESTEP\n" << timestep << "\n";
  out << "ITEM: NUMBER OF ATOMS\n" << atoms.size() << "\n";
  out << "ITEM: BOX BOUNDS pp pp pp\n";
  out << std::setprecision(10) << std::fixed;
  double bl=0;
  out << bl << " " << box.x << "\n";
  out << bl << " " << box.y << "\n";
  out << bl << " " << box.z << "\n";

  std::map<int, std::string> names = {
    {1, "C"},
    {2, "C"},
    {3, "N"},
    {4, "H"},
    {5, "Li"},
    {6, "P"},
    {7, "F"},
  };

  out << "ITEM: ATOMS id type element x y z\n";
  out << std::setprecision(10) << std::fixed;
  for( int atom=0; atom< atoms.size(); atom++ ) {
    out << atom+1 << " " << atomTypes[atom] << " " << names[atomTypes[atom]] << " " << atoms[atom].x << " " <<  atoms[atom].y << " " << atoms[atom].z << "\n"; 
  }
}

void writeLammpsDumpLiPF6H2O(std::ostream& out, long timestep, const std::vector<Atom>& atoms, const std::vector<int> atomTypes, const Box& box) {
  out << "ITEM: TIMESTEP\n" << timestep << "\n";
  out << "ITEM: NUMBER OF ATOMS\n" << atoms.size() << "\n";
  out << "ITEM: BOX BOUNDS pp pp pp\n";
  out << std::setprecision(10) << std::fixed;
  double bl=0;
  out << bl << " " << box.x << "\n";
  out << bl << " " << box.y << "\n";
  out << bl << " " << box.z << "\n";

  std::map<int, std::string> names = {
    {1, "O"},
    {2, "H"},
    {3, "Li"},
    {4, "P"},
    {5, "F"},
  };

  out << "ITEM: ATOMS id type element x y z\n";
  out << std::setprecision(5) << std::fixed;
  for( int atom=0; atom< atoms.size(); atom++ ) {
    out << atom+1 << " " << atomTypes[atom] << " " << names[atomTypes[atom]] << " " << atoms[atom].x << " " <<  atoms[atom].y << " " << atoms[atom].z << "\n"; 
  }
}


///////////////// Basic trajectory analysis tools
float applyPBC(float x, float box){
  float hbox = box/2.0;
  float wrapped = fmod(x + hbox, box);
  if (wrapped < 0) wrapped += box;
  return wrapped - hbox;
}

double wrapPBC(double x, double box){
  return x - box*std::floor(x/box);
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

float sError ( const std::vector<float>& data) {
  float mu = mean(data);
  float var = variance(data);
  return std::sqrt( var /static_cast<float>(data.size()));
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

int safeDivision(int a, int b) {
  if (b==0) {
    throw std::runtime_error("Division by zero");
  }
  if ( a%b != 0) {
    throw std::runtime_error("Non-integer division result");
  }
  return a/b;
}


///////////End statistics

int main(int argc, char* argv[]) {
  if (argc != 8  ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << "dir_name density field eqtime(ns) natomsInSolvent natomsInCation natomsInAnion \n";
    return 1;
  } 

  std::string dirName=argv[1]; // name of directory to output results, in results directory
  std::string rho = argv[2]; 
  int oneSolvent = std::stoi(argv[5]); 
  int oneCation = std::stoi(argv[6]); 
  int oneAnion = std::stoi(argv[7]); 
  std::cout << oneSolvent << " " << oneCation << " " << oneAnion << "\n";

  int n = rho.size();
  long long val = std::stoll(rho);
  float density = static_cast<float>(val) / std::pow(10.0, n-1);
  std::string fieldname = argv[3];
  float field = fieldConvert(fieldname)* 0.023; //kcal /molA,
  float fconversion = 25.7/0.592 ; // kcal/molA to mV/A

  std::string topology = std::string("../data/topology/") + dirName + "D" + rho + ".top";
  std::cout << "-------Reading Lammps Input from " << topology << " ------------\n";
  LammpsData top = readLammpsData(topology);

  std::cout << "Atoms style : " 
            << (top.atomsStyleHint.empty() ? "<none>" : top.atomsStyleHint) << "\n";
  std::cout << "Read atoms: " << top.atoms.size() << "\n";
  int numatoms= static_cast<int>(top.atoms.size());
  Box box=top.box;
  float avogadro = 6.0221408;
  int numions= static_cast<int>(std::round(2.0*density * box.x * box.y * box.z *avogadro /10000) ); // cation N + anion N
  int numpion = numions/2;
  float eqtime = std::stof(argv[4]); // ns unit
  

  std::cout << "\nAnalysis Setups\n";
  std::cout << "dump Directory : ../data/dielectric" << dirName << "/\n"; 
  std::cout << "output Directory : ../results/dielectric/" << "/\n"; 
  std::cout << "numAtoms: " << numatoms <<  " numIons: " << numions << " numPairs: " << numpion << "\n";
  std::cout << "boxSizes x: " << box.x << " A y: " << box.y << " A z: " << box.z << " A \n";
  std::cout << "density of ions : " << density << "M \n";
  std::cout << "fieldStrength: " << field << " kcal/molA = " << field * fconversion << " mV/A" << "\n"; 
  std::cout << "eqTime: " << eqtime << " ns\n";
  std::cout << "-----------End Reading Lammps Input File --------------";
  std::cout << "\n\n";

  //sorting topology.atom
  auto atoms = top.atoms;
  std::sort(atoms.begin(), atoms.end(), [](const AtomRow& a, const AtomRow&b) { return a.index < b.index;});

  std::cout << "first 20 atom \n" ;
  for(int i=0; i<20; i++) {
    std::cout << atoms[i].index <<  " " << atoms[i].mol << " " << atoms[i].type << " " << atoms[i].charge << "\n";
  }
  // set atom info
  std::vector<int> types; types.reserve(atoms.size());
  std::vector<double> charges; charges.reserve(atoms.size());
  for(int i=0; i<atoms.size(); i++) {
    types.push_back(atoms[i].type);
    charges.push_back(atoms[i].charge);
  }

  std::string inputFileName = std::string("../data/videos/videoD" + rho )+ "E" + fieldname + ".dcd";
  std::cout << "----------- Start Reading dcd File from " << inputFileName << "------------------------------------\n";
  DCDReader reader(inputFileName);
  int numatomsfromreader = reader.natoms();
  if ( numatomsfromreader != numatoms) std::cout << "Number of atoms in topology doesn't match with dcd file\n";
  int numsnap = reader.nframes();
  auto frames = reader.read_all();
  std::cout << "Read " << numsnap << " frames from dcd file\n";
  std::cout << "----------- End Reading dcd File ------------------\n";
  // frames : std::vector<Frame>
  // Frame : std::vector<Atom> atoms (x, y, z), Box box


  int numsolvents = safeDivision(numatoms - numpion*oneCation - numpion * oneAnion,oneSolvent);
  std::cout << "One solvent contains " << oneSolvent << " atoms, and total number of moleucles : " << numsolvents << "\n";

  std::string outFile;
  outFile=std::string("../results/videos/") + dirName + "D" + rho + "E" + fieldname + ".lammpstrj";
  std::ofstream out(outFile );

  for( size_t time = numsnap/2; time < numsnap; time++ ) {
    auto frame= frames[time].atoms;
    writeLammpsDumpLiPF6H2O(out, time, frame, types, box);
  }
  out.close();

  

  return 0;
}

# LiHalides
Analysis code for Ion association/dissociation rate measurement, TPT analysis
Electric field strength ranged from 0 to 22 meV/A.

## Structure
./results : all results, including conductivity, rate, 
./data/cnnDist ./data/cnnAngl < place for saving conditioend nn info
./data/dumps<name>/ < add binary trajectory dump files, format : type, x, y, z

## external library
use namespace of ioTraj and atom structures, used for organic solvent project

## note
processTraj : Required to be run firstly, since it generates all conditioned nearest neighbor list.

###
processTraj: 

construct graph-based clusters of ion ( cutoff = CUTOFFin). 
Goal is to generate an order parameter for each ion that will be used as state space definition.
Standard error of means are computed from block-averaging of Peterson (Franchel and Smith). If SEMs find convergence, (consequtive 3 , tolerance of less than 5% ), we find error.

outputs:
1. The conditioned nearest neighbor distance will be saved in data/cnnDist directory.
2. The angle from the positive z axis ( the electric field direction) of the conditioned nearest neighbor counterion, saved in data/cnnAngle
3. Conditioned conductivity on the conditioned nearest neighbor distance will be saved in results/conductivity

File Format
density field cutoffin(A) cutoffout(A) trajLength(ns) catCond catCondErr anCond anCondErr allCond allCondErr CIP SSIP FREE popCIP popSSIP popFREE
note that catCond and anCond are total conductivity of ions present in box, while CIP... are computed using unique pairs of cation

Conductivity is computed in units of S cm^2/mol, separately by cation and anion


use script in scripts/conductivity to control parameters


###
rate:

based on the conditioned nearest neighbor distance saved in data/cnnDist, compute TPT-based rates and committor averages.
CUDA-necessary

outputs :
1. results/rate/rateD00E00.dat
File Format
density | field | cutoffin(A) | cutoffout(A) | trajLength(ns) | 
q+ | q+err | q- | q-err | Kr(/ns) | Krerr(/ns) | krec(/ns) | krecErr(/ns) | kdis(/ns) | kdisErr(/ns) | kbim(/ns/M) | kbimerr(/ns/M) | Kass(/M) | KassErr(/M)

check if Error message comes out, which tells you how the block-averaging works. If it has error, the SEM may not be accurate.

Jul 14 2025/ modify not to use block averaging for MFPT measurement, since in each trajectory there are only few trajectories.
given that each crossing event are very uncorrelated, we measured standard error of mean for all colelctions of reactive trajectory

### 
rdf:

based on the conditioned nearest neighbor distance saved in data/cnnDist, generate 1-dimensional histogram to get radial distribution function

outputs :
1. results/hist/rdfD00E00.dat
File Format
cnn distance (A) | -ln( g(r)  )|  -ln( g(r) r^2)


### datasets details
data/dumpsLiIinH2O 

data/dumpsLiPF6inACN

data/dumpsNaIinACN
mostly aggregate even with strong field ~50

data/dumpsLiBF4inACN
Tom's ACN model with LiBF4

dumps file with partial thermostat, ( apply thermostat only on x and y direction, and only for solvents), are named as Tdumps...

### conductivity
measure conductivity. 
if field = 0
conductivity from Green-Kubo relation is included
conductivity from Einstein-Helfand equation

otherwise
conductivity from differential conductivity (linear response around Neq steady state)



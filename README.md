# LiHalides
Analysis code for Ion association/dissociation rate measurement, TPT analysis
Currently, all data come from LiCl of density 0.25M and 0.5M, with electric field strength ranged from 0 to 22 meV/A.

###
conductivity:

construct graph-based clusters of ion and measure a conditioned nearest neighbor counter ion distance.
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


rate:
based on the conditioned nearest neighbor distance saved in data/cnnTraj, compute TPT-based rates and committor averages.



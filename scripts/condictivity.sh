#!/bin/bash

target="LiIinH2O"
#target="Test"
density="10" ## ion density, /10 M
cutoffin="4.0" ## inner cutoff for domain, A
cutoffout="6" ## outer cutoff for domain, A
boxX="29.84" ## box size A
boxY="29.84" ## box size A
boxZ="59.68" ## box size A
timestep="0.05" ## time step of snapshots, ps
eqtime="0" ## equilibration time, ns
#dumpOrder="1" ## 1 if alternates cation and anion, 0 for sequencial cation and next anion
#numtraj="1" ## if nonunity, dump files are named with integer at the end


fields=("00" "02" "04" "06" "08" "10" "12" "14" "16" "18" "20" "22" ) ## external field, meV/A
fields=("09")

for field in "${fields[@]}";
do
  mkdir -p ../results/conductivity
  ../bin/conductivity $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime 
  echo "$field is done"
done
#./computeConductivity $density $field $cutoffin $cutoffout $box $timestep 2 $eqtime

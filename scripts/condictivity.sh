#!/bin/bash

density="025" ## ion density, /10 M
cutoffin="3." ## inner cutoff for domain, A
cutoffout="6" ## outer cutoff for domain, A
box="29.84" ## box size A, (note z axis *2)
timestep="0.05" ## time step of snapshots, ps
center="1" ## 1 for cation, 2 for anion
eqtime="5" ## equilibration time, ns

fields=("00" "02" "04" "06" "08" "10" "12" "14" "16" "18" "20" "22" ) ## external field, meV/A

for field in "${fields[@]}";
do
  ../bin/conductivity $density $field $cutoffin $cutoffout $box $timestep $center $eqtime
  echo "$field is done"
done
#./computeConductivity $density $field $cutoffin $cutoffout $box $timestep 2 $eqtime

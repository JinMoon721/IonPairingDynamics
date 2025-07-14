#!/bin/bash

density="025" ## ion density, /10 M
cutoffin="3." ## inner cutoff for domain, A
cutoffout="6" ## outer cutoff for domain, A
timestep="0.05" ## time step of snapshots, ps

fields=("00" "02" "04" "06" "08" "10" "12" "14" "16" "18" "20" "22" ) ## external field, meV/A
for field in "${fields[@]}";
do
  ../bin/rate $density $field $cutoffin $cutoffout $timestep
done


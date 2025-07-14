#!/bin/bash

density="05" ## ion density, /10 M
field="10" ## external field, meV/A
cutoffin="3." ## inner cutoff for domain, A
cutoffout="6" ## outer cutoff for domain, A
box="29.84" ## box size A, (note z axis *2)
timestep="0.05" ## time step of snapshots, ps
center="1" ## 1 for cation, 2 for anion
eqtime="5" ## equilibration time, ns


../bin/conductivity $density $field $cutoffin $cutoffout $box $timestep $center $eqtime
#./computeConductivity $density $field $cutoffin $cutoffout $box $timestep 2 $eqtime

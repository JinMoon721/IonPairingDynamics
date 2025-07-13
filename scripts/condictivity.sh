#!/bin/bash

density="025"
field="10"
cutoffin="3."
cutoffout="6"
box="29.84"
timestep="0.05"
center="1"


./computeConductivity $density $field $cutoffin $cutoffout $box $timestep 1
./computeConductivity $density $field $cutoffin $cutoffout $box $timestep 2

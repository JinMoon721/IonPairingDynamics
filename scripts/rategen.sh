#!/bin/bash

cutoffin=("4.3" "4.5" "4.7" "4.9" "5.1" "5.3" "5.5" "5.7" "5.9")
cutoffout=("7.0" "7.3" "7.6" "7.9" "8.2" "8.5" "8.8" "9.1" "9.4" "9.7" "10.0" "10.3" "10.6" "10.9" "11.2")


len1=${#cutoffin[@]}
len2=${#cutoffout[@]}

#rm rate.lst

count=0
target="KPF6inACN"
density="05"
field="00"
timestep="0.05"
for(( i=0; i<len1; i++ ));
do
  for(( j=0; j<len2; j++));
  do
    in=${cutoffin[i]}
    out=${cutoffout[j]}
    count=$((count+1))
    echo "../bin/rate $target $density $field $in $out $timestep $count" >> rate.lst
  done
done


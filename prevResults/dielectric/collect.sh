#!/bin/bash


FIELD=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09" "11" "13" "15" "17" "19" "22" "25" "30" "35" "40" "45" "50"  "60" "70")
len=${#FIELD[@]}

name="LiPF6inACN"
#name="LiIinH2O"
#name="Stock"
density="025"
rm dielectric${name}D$density
for(( i=0; i<len; i++));
do
  field=${FIELD[i]}
  cat ${name}D${density}E${field}.dat >> dielectric${name}D$density
done

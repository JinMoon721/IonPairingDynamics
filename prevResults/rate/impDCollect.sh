#!/bin/bash


DENSITY=("010")
len1=${#DENSITY[@]}

FIELD=("20" "30" "40" "50" "60" )
len2=${#FIELD[@]}

name="ImpD010"

rm rate${name}
for(( i=0; i<len1; i++));
do
  for(( j=0; j<len2; j++));
  do
    density=${DENSITY[i]}
    field=${FIELD[j]}
  cat ${name}D${density}E${field}.dat >> rate${name}
  done
done

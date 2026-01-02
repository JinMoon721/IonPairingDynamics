#!/bin/bash


DENSITY=("001" "002" "005" "010" "020" "050")
len1=${#DENSITY[@]}

FIELD=("00" "10" "20" "30")
len2=${#FIELD[@]}

name="ImpX20"

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

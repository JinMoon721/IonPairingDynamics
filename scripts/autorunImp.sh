#!/bin/bash
set -Eeuo pipefail

declare -A ImpX20=(
  [cutoffin]="4.3" ## inner cutoff for domain, A
  [cutoffout]="8.0" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="5.00" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="6 1 7"
)

measures=(processTraj rate conductivity)

#echo "Choose a measurement : "
#select measure in "${measures[@]}"; do
#  [[ -n "${measure:-}" ]] && break
#  echo "Invalid choice. Try again."
#done

measure=$1


targets=(ImpX20)

load_params() {
  local tgt=$1
  declare -n P="${tgt}"
  if [[ ${#P[@]} -eq 0 ]]; then
    echo "ERROR: unknown target ' ${tgt}. Known: ${targets[*]}" >&2
    return 1
  fi

  for k in "${!P[@]}"; do
    printf -v "$k" "%s" "${P[$k]}"
  done
}


#echo "Choose a target :"
#select target in "${targets[@]}"; do
#  [[ -n "${targets:-}" ]] && break
##    [[ -n "$target" ]] && break
#  echo "Invalid choice. Try again."
#done
target=$2

echo "Selected target trajectory : $target"
load_params "$target"


#DENSITY=("001" "002" "005" "010" "050" "100")
DENSITY=("001" "002" "005" "010" "020"  )
BOX=("132.62" "105.26" "77.56" "61.56" "48.86"  )
FIELD=( "00" "5" "10" "15" "20" "25" "30" "40" "50" "60" "70" "80" "90" "100" )

len1=${#DENSITY[@]}
len2=${#FIELD[@]}

if [[ $measure == "processTraj" ]]; then
  dump="../data/dumps$target/"
  echo "Analysis will be done on trajectory file : $dump"

  mkdir -p ./data/cnnDist
  for(( i=0; i<len1; i++));
  do
    for(( j=0; j<len2; j++));
    do
      field=${FIELD[j]}
      density=${DENSITY[i]}
      box=${BOX[i]}
      ../bin/processTraj $target $density $field $cutoffin $cutoffout $box $box $box $timestep $eqtime $videoflag $thermoflag 
    done
  done
elif [[ $measure == "rate" ]]; then
  dump="../data/cnnDist/"
  echo "Analysis will be done on trajectory file : $dump"
  mkdir -p ../results/rate
  for(( i=0; i<len1; i++));
  do
    for(( j=0; j<len2; j++));
    do
      field=${FIELD[j]}
      density=${DENSITY[i]}
      ../bin/rate $target $density $field $cutoffin $cutoffout $timestep $thermoflag
    done
  done
elif [[ $measure == "conductivity" ]]; then
  dump="../data/dumps$target/"
  echo "Analysis will be done on trajectory file : $dump"
  mkdir -p ../results/conductivity
  for(( i=0; i<len1; i++));
  do
    for(( j=0; j<len2; j++));
    do
      field=${FIELD[j]}
      density=${DENSITY[i]}
      box=${BOX[i]}
      ../bin/conductivity $target $density $field $cutoffin $cutoffout $box $box $box $timestep $eqtime $thermoflag 
    done
  done
fi

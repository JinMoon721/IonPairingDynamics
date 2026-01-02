#!/bin/bash
set -Eeuo pipefail

declare -A LiIinH2O=(
  [cutoffin]="3.1"
  [cutoffout]="6.0" ## outer cutoff for domain, A
  [boxX]="29.84" ## box size A
  [boxY]="29.84" ## box size A
  [boxZ]="29.84" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="5" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="3 1 1" ## number of atoms in solvent, cation, anion
)
declare -A LiPF6inACN=(
  [cutoffin]="4.3" ## inner cutoff for domain, A
  [cutoffout]="8.0" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="5" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="6 1 7"
)
declare -A KPF6inACN=(
  [cutoffin]="4.3" ## inner cutoff for domain, A
  [cutoffout]="8.0" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="5" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="6 1 7"
)
declare -A LiPF6inH2O=(
  [cutoffin]="4.3" ## inner cutoff for domain, A
  [cutoffout]="8.0" ## outer cutoff for domain, A
  [boxX]="28.57" ## box size A
  [boxY]="28.57" ## box size A
  [boxZ]="28.57" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="5" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="6 1 7"
)
declare -A TLiPF6inACN=(
  [cutoffin]="4.3" ## inner cutoff for domain, A
  [cutoffout]="8.0" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="5" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="6 1 7"
)

declare -A Stock=(
  [cutoffin]="3.1"
  [cutoffout]="6.0" ## outer cutoff for domain, A
  [boxX]="29.84" ## box size A
  [boxY]="29.84" ## box size A
  [boxZ]="29.84" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)
declare -A StockD=(
  [cutoffin]="3.1"
  [cutoffout]="6.0" ## outer cutoff for domain, A
  [boxX]="29.84" ## box size A
  [boxY]="29.84" ## box size A
  [boxZ]="29.84" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)

measures=(processTraj rate dielectric conductivity)
#echo "Choose a measurement : "
#select measure in "${measures[@]}"; do
#  [[ -n "${measure:-}" ]] && break
#  echo "Invalid choice. Try again."
#done

measure=$1


targets=(LiIinH2O LiPF6inACN Stock StockD TLiPF6inACN LiPF6inH2O)

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
#  echo "Invalid choice. Try again."
#done
target=$2

echo "Selected target trajectory : $target"
load_params "$target"

densities=("05" "025")
#echo "Choose a density :"
#select density in "${densities[@]}"; do
#  [[ -n "${densities:-}" ]] && break
#  echo "Invalid choice. Try again."
#done
density="05"

echo "Selected target density trajectory : $density"


#FIELD=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09" "11" "13" "15" "17" "19" "22" "25" "30" "35" "40" "45" )
FIELD=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09" "11" "13" "15" "17" "19" "22" "25" "30" "35" "40" "45" )
#FIELD=( "16" "17" "18" "19" "20")
len=${#FIELD[@]}

if [[ $measure == "processTraj" ]]; then
  dump="../data/dumps$target/"
  echo "Analysis will be done on trajectory file : $dump"

  mkdir -p ./data/cnnDist
  for(( i=0; i<len; i++));
  do
    field=${FIELD[i]}
    ../bin/processTraj $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime $videoflag $thermoflag 
  done
elif [[ $measure == "rate" ]]; then
  dump="../data/cnnDist/"
  echo "Analysis will be done on trajectory file : $dump"
  mkdir -p ../results/rate
  for(( i=0; i<len; i++));
  do
    field=${FIELD[i]}
    ../bin/rate $target $density $field $cutoffin $cutoffout $timestep $thermoflag
  done
elif [[ $measure == "dielectric" ]]; then
  dump="../data/dielectric$target/"
  echo "Analysis will be done on trajectory file : $dump"
  mkdir -p ../results/dielectric
  for(( i=0; i<len; i++));
  do
    field=${FIELD[i]}
    ../bin/dielectric $target $density $field $eqtime $natoms
  done
elif [[ $measure == "conductivity" ]]; then
  dump="../data/dumps$target/"
  echo "Analysis will be done on trajectory file : $dump"
  mkdir -p ../results/conductivity
  for(( i=0; i<len; i++));
  do
    field=${FIELD[i]}
    ../bin/conductivity $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime $thermoflag 
  done
fi

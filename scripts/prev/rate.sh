#!/bin/bash

set -Eeuo pipefail

declare -A LiIinH2O=(
  [density]="10"
  [cutoffin]="3.1"
  [cutoffout]="6.0" ## outer cutoff for domain, A
  [boxX]="29.84" ## box size A
  [boxY]="29.84" ## box size A
  [boxZ]="59.68" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)
declare -A LiPF6inACN=(
  [density]="05" ## ion density, /10 M
  [cutoffin]="4.3" ## inner cutoff for domain, A
  [cutoffout]="8.0" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)
declare -A NaIinACN=(
  [density]="05" ## ion density, /10 M
  [cutoffin]="4.0" ## inner cutoff for domain, A
  [cutoffout]="8.0" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)
declare -A LiBF4inACN=(
  [density]="05" ## ion density, /10 M
  [cutoffin]="4." ## inner cutoff for domain, A
  [cutoffout]="8.0" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)
declare -A LiIinACN=(
  [density]="05" ## ion density, /10 M
  [cutoffin]="3.1" ## inner cutoff for domain, A
  [cutoffout]="6.0" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)
declare -A LiPF6inACNionField=(
  [density]="05" ## ion density, /10 M
  [cutoffin]="4.3" ## inner cutoff for domain, A
  [cutoffout]="8.0" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)
declare -A LiPF6inACNacnField=(
  [density]="05" ## ion density, /10 M
  [cutoffin]="4.3" ## inner cutoff for domain, A
  [cutoffout]="8.0" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)
declare -A LiIinH2OionField=(
  [density]="10" ## ion density, /10 M
  [cutoffin]="3.1" ## inner cutoff for domain, A
  [cutoffout]="6.0" ## outer cutoff for domain, A
  [boxX]="29.84" ## box size A
  [boxY]="29.84" ## box size A
  [boxZ]="59.68" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)
declare -A LiIinH2Oh2oField=(
  [density]="10" ## ion density, /10 M
  [cutoffin]="3.1" ## inner cutoff for domain, A
  [cutoffout]="6.0" ## outer cutoff for domain, A
  [boxX]="29.84" ## box size A
  [boxY]="29.84" ## box size A
  [boxZ]="59.68" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)

targets=(LiIinH2O LiPF6inACN NaIinACN LiBF4inACN LiIinACN LiPF6inACNionField LiPF6inACNacnField LiIinH2OionField LiIinH2Oh2oField)

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

target="${1-}"
if [[ -z "${target}" && -n "${TARGET-}" ]]; then target="$TARGET"; fi
if [[ -z "${target}" ]]; then
  echo "Choose target:"
  select t in "${targets[@]}"; do
    [[ -n "${t:-}" ]] && target="$t" && break
  done
fi

load_params "$target"
thermoflag=0
fields=("50")


for field in "${fields[@]}";
do
  ../bin/rate $target $density $field $cutoffin $cutoffout $timestep $thermoflag
done


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

# choose_option VAR "Prompt" option1 option2 ...
choose_option() {
  local __var="$1"; shift
  local prompt="$1"; shift
  local -a opts=("$@")

  # CLI arg or env override (FIELD) if present
  local value="${1-}"; [[ -n "${value}" ]] && { printf -v "$__var" %s "$value"; return; }
  [[ -n "${FIELD-}" ]] && { printf -v "$__var" %s "$FIELD"; return; }

  # menu
  PS3="${prompt} "
  local pick
  select pick in "${opts[@]}"; do
    [[ -n "${pick:-}" ]] && { printf -v "$__var" %s "$pick"; break; }
    echo "Invalid choice."
  done
  PS3=
}


# usage
choose_option field "Choose field (number):" "${fields[@]}"
echo "field=$field"



load_params "$target"
videoflag=0
thermoflag=0

fields=("20")
#density="05"


for field in "${fields[@]}";
do
  mkdir -p ../results/conductivity
  #../bin/conductivity $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime $videoflag $thermoflag 
  ../bin/conductivity $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime $videoflag $thermoflag   
  echo "$field is done"
done

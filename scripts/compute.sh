#!/bin/bash
set -Eeuo pipefail

declare -A LiIinH2O=(
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
declare -A LiTFSIinACN=(
  [cutoffin]="3.5" ## inner cutoff for domain, A
  [cutoffout]="6.6" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="2.8652" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
)


measures=(processTraj rate dielectric conductivity)
echo "Choose a measurement : "
select measure in "${measures[@]}"; do
  [[ -n "${measure:-}" ]] && break
  echo "Invalid choice. Try again."
done


targets=(LiIinH2O LiPF6inACN NaIinACN LiBF4inACN LiIinACN LiPF6inACNionField LiPF6inACNacnField LiIinH2OionField LiIinH2Oh2oField LiTFSIinACN)

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

echo "Choose a target :"
select target in "${targets[@]}"; do
  [[ -n "${targets:-}" ]] && break
  echo "Invalid choice. Try again."
done

echo "Selected target trajectory : $target"
load_params "$target"

if [[ $measure == "processTraj" ]]; then
  dumpdir="../data/dumps$target/"
  pattern="dumpD"

  ## find all trajectories
  shopt -s nullglob
  dynamic=()
  for f in "$dumpdir"/*"$pattern"*; do
    [[ -f $f ]] && dynamic+=( "$(basename -- "$f")")
  done
  shopt -u nullglob

  dumps=()
  dumps+=("${dynamic[@]}")

  if ((${#dumps[@]} == 0)); then
    echo "Dump trajectories are not available (pattern='${pattern}' in dir='${dumpdir}')."
    exit 1
  fi

  echo "Choose a dump file :"
  select dump in "${dumps[@]}"; do
    [[ -n "${dump:-}" ]] && break
    echo "Invalid choice. Try again."
  done

  echo "Analysis will be done on trajectory file : $dump"

  ## get density and field
  [[ $dump =~ D([0-9]+)E([0-9]+) ]] && density="${BASH_REMATCH[1]}" field="${BASH_REMATCH[2]}"

  if [[ $dump == *Tdump* ]]; then
    thermoflag=1
    if [[ $target == *LiIinH2O* ]]; then
      boxZ="29.84"
    fi
  fi

  if [[ $dump == *TTdump* ]]; then
    thermoflag=2
    if [[ $target == *LiIinH2O* ]]; then
      boxZ="29.84"
    fi
  fi

  #mkdir -p ../results/conductivity
  ../bin/processTraj $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime $videoflag $thermoflag 
elif [[ $measure == "rate" ]]; then
  dumpdir="../data/cnnDist/"
  pattern="$target"

  ## find all trajectories
  shopt -s nullglob
  dynamic=()
  for f in "$dumpdir"/*"$pattern"*; do
    [[ -f $f ]] && dynamic+=( "$(basename -- "$f")")
  done
  shopt -u nullglob

  dumps=()
  dumps+=("${dynamic[@]}")

  if ((${#dumps[@]} == 0)); then
    echo "CNN trajectories are not available (pattern='${pattern}' in dir='${dumpdir}')."
    exit 1
  fi

  echo "Choose a data file :"
  select dump in "${dumps[@]}"; do
    [[ -n "${dump:-}" ]] && break
    echo "Invalid choice. Try again."
  done

  echo "Analysis will be done on trajectory file : $dump"

  ## get density and field
  [[ $dump =~ D([0-9]+)E([0-9]+) ]] && density="${BASH_REMATCH[1]}" field="${BASH_REMATCH[2]}"

  if [[ $dump == *TD* ]]; then
    thermoflag=1
    if [[ $target == *LiIinH2O* ]]; then
      boxZ="29.84"
    fi
  fi

  if [[ $dump == *TTD* ]]; then
    thermoflag=2
    if [[ $target == *LiIinH2O* ]]; then
      boxZ="29.84"
    fi
  fi

  mkdir -p ../results/rate
  ../bin/rate $target $density $field $cutoffin $cutoffout $timestep $thermoflag
elif [[ $measure == "dielectric" ]]; then
  dumpdir="../data/dielectric$target/"
  pattern="video"

  ## find all trajectories
  shopt -s nullglob
  dynamic=()
  for f in "$dumpdir"/*"$pattern"*; do
    [[ -f $f ]] && dynamic+=( "$(basename -- "$f")")
  done
  shopt -u nullglob

  dumps=()
  dumps+=("${dynamic[@]}")

  if ((${#dumps[@]} == 0)); then
    echo "Full trajectories are not available (pattern='${pattern}' in dir='${dumpdir}')."
    exit 1
  fi

  echo "Choose a trajectory file :"
  select dump in "${dumps[@]}"; do
    [[ -n "${dump:-}" ]] && break
    echo "Invalid choice. Try again."
  done

  echo "Analysis will be done on trajectory file : $dump"

  ## get density and field
  [[ $dump =~ D([0-9]+)E([0-9]+) ]] && density="${BASH_REMATCH[1]}" field="${BASH_REMATCH[2]}"

  if [[ $dump == *TD* ]]; then
    thermoflag=1
    if [[ $target == *LiIinH2O* ]]; then
      boxZ="29.84"
    fi
  fi

  if [[ $dump == *TTD* ]]; then
    thermoflag=2
    if [[ $target == *LiIinH2O* ]]; then
      boxZ="29.84"
    fi
  fi

  mkdir -p ../results/dielectric
  ../bin/dielectric $target $density $field $eqtime
elif [[ $measure == "conductivity" ]]; then
  dumpdir="../data/dumps$target/"
  pattern="dumpD"

  ## find all trajectories
  shopt -s nullglob
  dynamic=()
  for f in "$dumpdir"/*"$pattern"*; do
    [[ -f $f ]] && dynamic+=( "$(basename -- "$f")")
  done
  shopt -u nullglob

  dumps=()
  dumps+=("${dynamic[@]}")

  if ((${#dumps[@]} == 0)); then
    echo "Dump trajectories are not available (pattern='${pattern}' in dir='${dumpdir}')."
    exit 1
  fi

  echo "Choose a dump file :"
  select dump in "${dumps[@]}"; do
    [[ -n "${dump:-}" ]] && break
    echo "Invalid choice. Try again."
  done

  echo "Analysis will be done on trajectory file : $dump"

  ## get density and field
  [[ $dump =~ D([0-9]+)E([0-9]+) ]] && density="${BASH_REMATCH[1]}" field="${BASH_REMATCH[2]}"

  if [[ $dump == *Tdump* ]]; then
    thermoflag=1
    if [[ $target == *LiIinH2O* ]]; then
      boxZ="29.84"
    fi
  fi

  if [[ $dump == *TTdump* ]]; then
    thermoflag=2
    if [[ $target == *LiIinH2O* ]]; then
      boxZ="29.84"
    fi
  fi

  mkdir -p ../results/conductivity
  ../bin/conductivity $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime $thermoflag 
fi

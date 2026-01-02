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
  [cutoffin]="5.3" ## inner cutoff for domain, A
  [cutoffout]="7.2" ## outer cutoff for domain, A
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
  [cutoffin]="5.3" ## inner cutoff for domain, A
  [cutoffout]="7.2" ## outer cutoff for domain, A
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
  [cutoffin]="5.3" ## inner cutoff for domain, A
  [cutoffout]="7.2" ## outer cutoff for domain, A
  [boxX]="28.57" ## box size A
  [boxY]="28.57" ## box size A
  [boxZ]="28.57" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="5" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="3 1 7"
)

declare -A ImpX20=(
  [cutoffin]="3.5" ## inner cutoff for domain, A
  [cutoffout]="8.5" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="5.00" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="6 1 7"
)
declare -A ImpD010=(
  [cutoffin]="3.0" ## inner cutoff for domain, A
  [cutoffout]="5.5" ## outer cutoff for domain, A
  [boxX]="36.00" ## box size A
  [boxY]="36.00" ## box size A
  [boxZ]="36.00" ## box size A
  [timestep]="0.05" ## time step of snapshots, ps
  [eqtime]="5.00" ## equilibration time, ns
  [videoflag]="0"
  [thermoflag]="0"
  [natoms]="6 1 7"
)

measures=(processTraj rate dielectric conductivity rdf mtpt distribution mechanism singleboundary)
echo "Choose a measurement : "
select measure in "${measures[@]}"; do
  [[ -n "${measure:-}" ]] && break
  echo "Invalid choice. Try again."
done


targets=(LiIinH2O LiPF6inACN LiPF6inH2O ImpD010 ImpX20 KPF6inACN)

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
elif [[ $measure == "singleboundary" ]]; then
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

  mkdir -p ../results/singleboundary
  ../bin/singleboundary $target $density $field $cutoffin $cutoffout $timestep $thermoflag
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
  ../bin/dielectric $target $density $field $eqtime $natoms
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
elif [[ $measure == "rdf" ]]; then
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

  mkdir -p ../results/rdf
  echo $target $density $field
  ../bin/rdf $target $density $field $cutoffin $cutoffout $timestep $thermoflag
elif [[ $measure == "mtpt" ]]; then
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

  mkdir -p ../results/mtpt
  ../bin/mtpt $target $density $field $cutoffin $cutoffout $timestep $thermoflag
elif [[ $measure == "distribution" ]]; then
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

  mkdir -p ../results/distribution
  ../bin/distribution $target $density $field $cutoffin $cutoffout $timestep $thermoflag
elif [[ $measure == "mechanism" ]]; then
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
  mkdir -p ../results/mechanism
  ../bin/mechanism $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime  
fi

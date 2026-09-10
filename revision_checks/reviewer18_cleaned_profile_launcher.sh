#!/bin/sh
set -eu

root=${1:-/home/amizuno/phylo_spatial_tutorial_revision}
workers=${2:-16}
script="$root/revision_checks/reviewer18_cleaned_profile_point.R"
out="$root/revision_checks/reviewer18_cleaned_primary_outputs"
profiles="$out/profiles"
grid="$out/profile_grid.csv"
prepared="$root/revision_checks/reviewer18_influential_effects_outputs/published_cleaned_prepared.rds"
models="$root/revision_checks/reviewer18_influential_effects_outputs"

test -f "$script"
test -f "$grid"
test -f "$prepared"
test -f "$models/spatial_only.rds"
test -f "$models/combined.rds"
mkdir -p "$profiles"

existing=$(find "$profiles" -maxdepth 1 -type f -name '*.csv' | wc -l | tr -d ' ')
if [ "$existing" -ne 0 ]; then
  echo "Refusing to relaunch over $existing existing profile files" >&2
  exit 2
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

awk -F, 'NR > 1 {gsub(/"/, "", $1); gsub(/"/, "", $2); print $1, $2, $3}' "$grid" |
  xargs -P "$workers" -n 3 sh -c '
    script=$1
    prepared=$2
    models=$3
    profiles=$4
    model=$5
    component=$6
    value=$7
    Rscript "$script" "$prepared" "$models/$model.rds" "$profiles" "$model" "$component" "$value"
  ' sh "$script" "$prepared" "$models" "$profiles"

completed=$(find "$profiles" -maxdepth 1 -type f -name '*.csv' | wc -l | tr -d ' ')
echo "PROFILE_GRID_COMPLETE points=$completed workers=$workers"
test "$completed" -eq 43

#!/bin/sh
set -eu

root=${1:-/home/amizuno/phylo_spatial_tutorial_revision}
workers=${2:-9}
script="$root/revision_checks/reviewer18_cleaned_targeted_multistart.R"
prepared="$root/revision_checks/reviewer18_influential_effects_outputs/published_cleaned_prepared.rds"
primary="$root/revision_checks/reviewer18_influential_effects_outputs/spatial_only.rds"
out="$root/revision_checks/reviewer18_cleaned_primary_outputs/targeted_multistart"
mkdir -p "$out"

existing=$(find "$out" -maxdepth 1 -type f -name '*.csv' | wc -l | tr -d ' ')
if [ "$existing" -ne 0 ]; then
  echo "Refusing to relaunch over $existing existing multistart files" >&2
  exit 2
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

for tau2 in 0.8 1.0 1.3; do
  for rho_start in 0.02 0.14 0.3; do
    printf '%s %s\n' "$tau2" "$rho_start"
  done
done |
  xargs -P "$workers" -n 2 sh -c '
    Rscript "$1" "$2" "$3" "$4" "$5" "$6"
  ' sh "$script" "$prepared" "$primary" "$out"

completed=$(find "$out" -maxdepth 1 -type f -name '*.csv' | wc -l | tr -d ' ')
echo "TARGETED_MULTISTART_COMPLETE points=$completed workers=$workers"
test "$completed" -eq 9

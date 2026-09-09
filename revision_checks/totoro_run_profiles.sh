#!/usr/bin/env bash
# Run fixed-parameter REML profile points in bounded parallel batches on Totoro.
# Individual R processes write their own result CSV immediately on completion.
set -euo pipefail

workers=4
root_dir="${1:-$PWD}"
cd "$root_dir"
mkdir -p outputs/profiles

# Values are deliberately concentrated around each primary estimate, with
# enough logarithmic coverage to detect boundaries or a displaced maximum.
tasks=(
  "spatial_only tau2 0"
  "spatial_only tau2 0.25"
  "spatial_only tau2 0.5"
  "spatial_only tau2 0.8"
  "spatial_only tau2 1.0"
  "spatial_only tau2 1.5"
  "spatial_only tau2 2"
  "spatial_only tau2 3"
  "spatial_only tau2 5"
  "spatial_only rho 0.005"
  "spatial_only rho 0.01"
  "spatial_only rho 0.02"
  "spatial_only rho 0.03"
  "spatial_only rho 0.05"
  "spatial_only rho 0.08"
  "spatial_only rho 0.15"
  "spatial_only rho 0.3"
  "spatial_only rho 1"
  "spatial_only rho 3"
  "spatial_only rho 10"
  "spatial_only rho 100"
  "combined tau2 0"
  "combined tau2 0.005"
  "combined tau2 0.01"
  "combined tau2 0.02"
  "combined tau2 0.04"
  "combined tau2 0.06"
  "combined tau2 0.1"
  "combined tau2 0.2"
  "combined tau2 0.4"
  "combined tau2 0.8"
  "combined tau2 1.6"
  "combined rho 0.03"
  "combined rho 0.1"
  "combined rho 0.3"
  "combined rho 1"
  "combined rho 3"
  "combined rho 10"
  "combined rho 30"
  "combined rho 100"
  "combined rho 200"
  "combined rho 385"
  "combined rho 700"
  "combined rho 1500"
  "combined rho 3000"
  "combined rho 6000"
  "combined rho 12000"
)

run_point() {
  local model="$1"
  local component="$2"
  local value="$3"
  local tag="${model}_${component}_${value}"
  Rscript totoro_profile_point.R "outputs/${model}.rds" Roger_etal_2024.csv \
    outputs/profiles "$model" "$component" "$value" \
    > "outputs/profiles/${tag}.log" 2>&1
}

for task in "${tasks[@]}"; do
  read -r model component value <<< "$task"
  while (( $(jobs -rp | wc -l) >= workers )); do
    wait -n
  done
  run_point "$model" "$component" "$value" &
done
wait

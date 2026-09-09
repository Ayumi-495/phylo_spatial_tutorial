#!/usr/bin/env bash
# Bounded profile-likelihood grid for the corrected Scholer metafor audit.
# Each profile point is an independent R process and writes its own RDS/CSV
# before exiting. Eight workers are used, each limited to one BLAS/OpenMP thread.
set -euo pipefail

base_dir="/home/amizuno/phylo_spatial_tutorial_revision_checks/scholer_spatial_audit_2026-09-08"
workers=8

jobs=(
  "spatial_only tau2 0"
  "spatial_only tau2 0.02"
  "spatial_only tau2 0.05"
  "spatial_only tau2 0.1"
  "spatial_only tau2 0.2"
  "spatial_only tau2 0.3226583722"
  "spatial_only tau2 0.5"
  "spatial_only tau2 0.8"
  "spatial_only tau2 1.2"
  "combined tau2 0"
  "combined tau2 0.001"
  "combined tau2 0.003"
  "combined tau2 0.005"
  "combined tau2 0.01"
  "combined tau2 0.02"
  "combined tau2 0.04"
  "combined tau2 0.08"
  "combined tau2 0.16"
  "spatial_only rho 1"
  "spatial_only rho 5"
  "spatial_only rho 20"
  "spatial_only rho 50"
  "spatial_only rho 100"
  "spatial_only rho 168.7093726"
  "spatial_only rho 300"
  "spatial_only rho 600"
  "spatial_only rho 1200"
  "spatial_only rho 3000"
  "combined rho 10"
  "combined rho 30"
  "combined rho 100"
  "combined rho 200"
  "combined rho 400"
  "combined rho 535.5823133"
  "combined rho 800"
  "combined rho 1500"
  "combined rho 3000"
  "combined rho 6000"
  "combined rho 12000"
)

printf '%s\n' "${jobs[@]}" |
  xargs -n 3 -P "${workers}" bash -c '
    set -euo pipefail
    model="$1"; parameter="$2"; value="$3"
    cd "'"${base_dir}"'"
    PROFILE_APPEND=false OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      Rscript scholer_spatial_models.R profile_point scholer_prepared.rds output \
      "${model}" "${parameter}" "${value}"
  ' _

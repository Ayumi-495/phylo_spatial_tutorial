#!/usr/bin/env bash
# Run figure-only variance-component profile intervals on Totoro.
# At most eight serial metafor profile jobs run at once (eight CPU cores).
set -euo pipefail

root="$HOME/phylo_spatial_tutorial_revision"
script="$root/revision_checks/profile_variance_component_ci.R"
workers=8

jobs=(
  "reviewer18_influential_effects_outputs unstructured_only sigma2 1"
  "reviewer18_influential_effects_outputs unstructured_only sigma2 2"
  "reviewer18_influential_effects_outputs spatial_only sigma2 1"
  "reviewer18_influential_effects_outputs spatial_only tau2 1"
  "reviewer18_influential_effects_outputs combined sigma2 1"
  "reviewer18_influential_effects_outputs combined sigma2 2"
  "reviewer18_influential_effects_outputs combined tau2 1"
  "scholer_spatial_audit_outputs unstructured_only sigma2 1"
  "scholer_spatial_audit_outputs unstructured_only sigma2 2"
  "scholer_spatial_audit_outputs spatial_only sigma2 1"
  "scholer_spatial_audit_outputs spatial_only tau2 1"
  "scholer_spatial_audit_outputs combined sigma2 1"
  "scholer_spatial_audit_outputs combined sigma2 2"
  "scholer_spatial_audit_outputs combined tau2 1"
)

running=0
for spec in "${jobs[@]}"; do
  read -r audit model type index <<< "$spec"
  log_dir="$root/revision_checks/$audit/variance_profile_ci/logs"
  mkdir -p "$log_dir"
  (
    Rscript "$script" "$audit" "$model" "$type" "$index"
  ) >"$log_dir/${model}_${type}_${index}.log" 2>&1 &
  ((running += 1))
  if (( running >= workers )); then
    wait -n
    ((running -= 1))
  fi
done
wait

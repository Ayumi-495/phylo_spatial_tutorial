# Archived legacy spatial material

The historical spatial material was removed from `tutorial_v2.qmd` and remains
recoverable in Git history. The replacement spatial tutorial uses only the
audited workflows and numerical outputs in this directory.

The archived block contained:

- global coordinate transformations using EPSG:3857;
- old global `brms` and `glmmTMB` examples that were not matched to the
  verified great-circle `metafor` analysis;
- stale squared-exponential numerical output;
- the Scholer spatial term indexed by `effect_id` while a 454-location matrix
  was supplied;
- the undefined object `ma_sp2_exp`; and
- visualisations based on unaudited fit objects.

Replacement evidence:

- `totoro_spatial_audit_outputs/` for the full global Grau-Andrés comparison;
- `regional_cross_package_audit_outputs/` for the Spain implementation;
- `scholer_structure_audit_outputs/` and `scholer_spatial_audit_outputs/` for
  the corrected Scholer analysis.

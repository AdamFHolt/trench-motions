# Config Matrix

This folder defines reproducible misfit-sweep runs for core run comparisons.

Matrix dimensions currently included:
- Reference frame (`vt_ref`): `hs3`, `nnr`, `sa`
- Formulation suite:
  - `linear` (`formulation: 1`)
  - `powerlaw` (`formulation: 4`, `trans_strain_rate: 1e-13`)

Common settings:
- `include_DP: 1`
- `DP_ref: 2.35e7` (Pa)
- `PSP_slab_pull_factor: 0.25`
- `include_ridge_push: 1`
- `skip_map: true` (use quick plots; avoids GMT dependency)

Outputs are routed via `out_prefix` into `results/matrix/...`.

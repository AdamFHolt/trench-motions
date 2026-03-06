# Session Notes (Trench Motions)

## Current Status
- Repository is cleaned and pushed to `origin/main`.
- Core workflow is simplified around:
  - `make run-matrix-with-summary`
  - `make run-matrix-maps-with-summary` (with `DATASETS_DIR` set)
- `README.md` and `Makefile` have been pared down to fundamentals.

## Key Refactors Completed
- Removed legacy power-law mode from active workflow and core implementation.
- Removed `PSP_slab_pull_factor` from:
  - CLI args
  - configs
  - output naming
  - core function signature usage
- Consolidated matrix configs:
  - `configs/run_params.yaml`
- Added shared prep utilities:
  - `workflow_common.py`
  - used by both `compute_rates_misfit.py` and `compute_rates_single.py`
- Reorganized outputs:
  - `plots/<hs3|nnr|sa>/param-sweep/`
  - `plots/<hs3|nnr|sa>/maps/`
  - `plots/summary/param-sweep/`
  - `plots/summary/maps/`
- Removed results `tmp` outputs (temporary files now use OS temp dirs and are cleaned automatically).

## Minimal Daily Workflow
1. Activate environment:
   - `source env/bin/activate`
2. Run canonical sweep + summary:
   - `make run-matrix-with-summary`
3. (Optional) Run map workflow:
   - `export DATASETS_DIR=/path/to/datasets`
   - `make run-matrix-maps-with-summary`

## Important Paths
- Inputs:
  - `data/`
  - `data/vt/`
- Sweep results:
  - `plots/<hs3|nnr|sa>/param-sweep/`
- Map results:
  - `plots/<hs3|nnr|sa>/maps/`
- Summaries:
  - `plots/summary/param-sweep/`
  - `plots/summary/maps/`

## Suggested Next Steps (Science)
1. Define a small set of parameter experiments to run first (per reference frame: `hs3`, `nnr`, `sa`).
2. Compare summary outputs:
   - `plots/summary/param-sweep/matrix_summary.csv`
   - `rmse_by_run.png`, `sign_match_by_run.png`
3. Select best-fit cases and generate map products for interpretation.
4. Add short interpretation notes (per reference frame) in a new markdown file for your student handoff.

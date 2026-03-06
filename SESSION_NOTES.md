# Session Notes (Trench Motions)

## Current Status
- Repository is cleaned and pushed to `origin/main`.
- Core workflow is simplified around:
  - `make run-matrix-with-summary`
  - `make run-matrix-maps-with-summary`
- `README.md` and `Makefile` have been pared down to fundamentals.

## Key Refactors Completed
- Removed legacy power-law mode from active workflow and core implementation.
- Retooled active formulation numbering to a contiguous set:
  - `1..6` (old `7/8/9` remapped to new `4/5/6`)
- Removed `PSP_slab_pull_factor` from:
  - CLI args
  - configs
  - output naming
  - core function signature usage
- Consolidated run config to a single file:
  - `configs/run_params.yaml`
- Added shared prep utilities:
  - `workflow_common.py`
  - used by both `compute_rates_misfit.py` and `compute_rates_single.py`
- Replaced GMT map plotting with Python plotting:
  - `plot_trench_motions.py`
- Reorganized outputs:
  - `plots/<hs3|nnr|sa>/param-sweep/`
  - `plots/<hs3|nnr|sa>/maps/`
  - `plots/summary/`
- Removed results `tmp` outputs (temporary files now use OS temp dirs and are cleaned automatically).

## Minimal Daily Workflow
1. Activate environment:
   - `source env/bin/activate`
2. Run canonical sweep + summary:
   - `make run-matrix-with-summary`
3. (Optional) Run map workflow:
   - `make run-matrix-maps-with-summary`

## Important Paths
- Inputs:
  - `data/`
  - `data/vt/`
- Sweep results:
  - `plots/<hs3|nnr|sa>/param-sweep/`
- Map results:
  - `plots/<hs3|nnr|sa>/maps/`
- Combined summary outputs:
  - `plots/summary/`

## Suggested Next Steps (Science)
1. Define a small set of parameter experiments to run first (per reference frame: `hs3`, `nnr`, `sa`).
2. Compare summary outputs:
   - `plots/summary/matrix_summary.csv`
   - `rmse_by_run.png`, `sign_match_by_run.png`
3. Select best-fit cases and generate map products for interpretation.
4. Add short interpretation notes (per reference frame) in a new markdown file for your student handoff.

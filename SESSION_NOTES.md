# Session Notes (Trench Motions)

## Current Status
- Repository workflow is Python-only (no active GMT or shell-script execution).
- Core workflow is simplified around:
  - `make run-matrix-with-summary`
  - `make run-matrix-maps-with-summary`
- `README.md` and `Makefile` have been pared down to fundamentals.

## Key Refactors Completed
- Removed legacy power-law mode from active workflow and core implementation.
- Retooled active formulation numbering to a contiguous active set:
  - active set is now `1..5` (`3/4/5` correspond to former `4/5/6`)
- Removed slab-pull-prefactor case from active workflow and core implementation.
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
- Consolidated plotting into:
  - `plotting_functions.py`
  - includes misfit heatmaps, quick diagnostics, and full trench-motion map rendering
- Updated sweep behavior:
  - quick diagnostic plot is generated for best-fit results in every sweep run
  - `--skip-map` now skips only full map rendering
- Removed legacy script files:
  - `plot_trench_motions.py`
  - `plot_trench_motions.sh`
  - `quick_plot.py` (wrapper no longer needed)
  - legacy `archive/legacy/*.sh` run scripts
- Removed generated GMT session artifacts from repo root:
  - `gmt.conf`
  - `gmt.history`
- Reorganized outputs:
  - `plots/<hs3|nnr|sa>/<viscous|plastic|viscous_LspShear|viscous_VspShear|viscous_ShearOP>/param-sweep/`
  - `plots/<hs3|nnr|sa>/<viscous|plastic|viscous_LspShear|viscous_VspShear|viscous_ShearOP>/maps/`
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
  - `plots/<hs3|nnr|sa>/<model>/param-sweep/`
- Map results:
  - `plots/<hs3|nnr|sa>/<model>/maps/`
- Combined summary outputs:
  - `plots/summary/`

## Latest Run Results (Viscous, Formulation 1)
- Executed full viscous runs (with maps) for all reference frames using `configs/run_params.yaml` (`DP_ref=2.35e7`, `include_ridge_push=1`):
  - `hs3`: minimum RMS `4.01 cm/yr`, signs `70/120`
  - `nnr`: minimum RMS `2.27 cm/yr`, signs `94/120`
  - `sa`: minimum RMS `3.08 cm/yr`, signs `84/120`
- Current generated structure is consistent across all refs:
  - `plots/<hs3|nnr|sa>/viscous/param-sweep/`
  - `plots/<hs3|nnr|sa>/viscous/maps/`

## Suggested Next Steps (Science)
1. Define a small set of parameter experiments to run first (per reference frame: `hs3`, `nnr`, `sa`).
2. Compare summary outputs:
   - `plots/summary/matrix_summary.csv`
   - `rmse_by_run.png`, `sign_match_by_run.png`
3. Select best-fit cases and generate map products for interpretation.
4. Add short interpretation notes (per reference frame) in a new markdown file for your student handoff.

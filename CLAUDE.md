# CLAUDE.md — trench-motions

## What this project is
Analytical subduction/trench-motion modeling. Computes predicted vertical trench velocity for ~120 global subduction segments using physical formulations (slab pull, ridge push, dynamic pressure). Compares predictions to observed trench motion data in three reference frames (hs3, nnr, sa). Outputs: RMSE metrics, scatter plots, misfit heatmaps, global maps.

Main dataset: `Lallemand_et_al-2005_G3_dataset_WithAdditionalParams.txt`
Observed velocities: `data/vt/tnew.<hs3|nnr|sa>.dat`

## Development commands

```bash
# Activate environment
source env/bin/activate

# Canonical sweep + summary (no maps, fast)
make run-matrix-with-summary

# Full run with maps
make run-matrix-maps-with-summary

# Single smoke test
python compute_rates_single.py --config configs/single_smoke.yaml --skip-map

# Sweep smoke test
python compute_rates_misfit.py --config configs/misfit_smoke.yaml --vt-ref hs3 --skip-map
```

## Active scripts

| File | Role |
|------|------|
| `compute_rates_misfit.py` | Parameter sweep driver. Accepts `--config` YAML or positional args. |
| `compute_rates_single.py` | Single deterministic run with fixed viscosity values. |
| `functions.py` | Physics engine — `compute_vsp_withDP()` implements all 5 formulations. |
| `workflow_common.py` | Shared data prep: VT column mapping, segment arrays, depth/dip interpolation. |
| `plotting_functions.py` | `save_quick_plot()`, `save_misfit_heatmap()`, `save_trench_motion_map()`. |
| `matrix_summary.py` | Aggregates sweep results across all frames/models into CSV + bar plots. |

## Formulations (active set 1–4)

1. Viscous bending + dynamic pressure (canonical)
2. Plastic bending + DP
3. Viscous, layer thickness ∝ slab length (Lsp)
4. Viscous, layer thickness ∝ predicted velocity (Vsp) — closed-form quadratic

Canonical config (`configs/run_params.yaml`): formulation=1, DP_ref=2.35e7 Pa, include_ridge_push=1.

## Output layout

```
plots/<hs3|nnr|sa>/<model>/param-sweep/   # sweep plots and best-fit diagnostics
plots/<hs3|nnr|sa>/<model>/maps/          # global maps and prediction tables
plots/summary/param-sweep/               # matrix_summary.csv, RMSE/sign-match bar plots
plots/summary/maps/                      # summary across map runs
```

## Latest results (Formulation 1, viscous, DP_ref=2.35e7, with ridge push)

| Frame | Min RMSE | Sign match |
|-------|----------|------------|
| hs3   | 4.01 cm/yr | 70/120   |
| nnr   | 2.27 cm/yr | 94/120   |
| sa    | 3.08 cm/yr | 84/120   |

## Suggested next steps (science)

1. Run parameter experiments across all formulations (1–5) for all three reference frames.
2. Compare summary outputs: `plots/summary/param-sweep/matrix_summary.csv`, `rmse_by_run.png`, `sign_match_by_run.png`.
3. Select best-fit cases per frame and generate map products for interpretation.
4. Add per-reference-frame interpretation notes (e.g. for student handoff).

## Map data

Map plotting auto-detects:
- `data/age.3.6.NaN.grd` — oceanic plate age raster
- `data/PB2002_tdiddy.gmt` — Bird plate boundaries

Or set `DATASETS_DIR` to a root containing `age/age.3.6.NaN.grd` and `plate_boundaries/bird_PB2002/PB2002_tdiddy.gmt`. Maps run without these files, just without the background overlays.

# CLAUDE.md — trench-motions

## What this project is
Analytical subduction/trench-motion modeling. Computes predicted trench velocity for ~120 global subduction segments using physical force-balance formulations (slab pull, ridge push, dynamic pressure). Compares predictions to observed trench motion data in three reference frames (hs3, nnr, sa). Outputs: RMSE metrics, scatter plots, misfit heatmaps, global maps.

Main dataset: `Lallemand_et_al-2005_G3_dataset_WithAdditionalParams.txt`
Observed velocities: `data/vt/tnew.<hs3|nnr|sa>.dat`
Force-balance derivations: `docs/force_balance.md`
Sketch/figures: `plots/sketch/`

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

# Run a single formulation/frame only
make run-matrix FORMULATIONS=1 REF_FRAMES=nnr
```

## Active scripts

| File | Role |
|------|------|
| `compute_rates_misfit.py` | Parameter sweep driver. Accepts `--config` YAML or positional args. Key flags: `--formulation`, `--vt-ref`, `--skip-map`, `--smoke`. |
| `compute_rates_single.py` | Single deterministic run with fixed viscosity values. |
| `functions.py` | Physics engine — `compute_vsp_withDP()` implements formulations 1 and 2. |
| `workflow_common.py` | Shared data prep: VT column mapping, segment arrays, depth/dip interpolation. |
| `plotting_functions.py` | `save_quick_plot()`, `save_misfit_heatmap()`, `save_trench_motion_map()`, `save_vt_param_plot()`. |
| `matrix_summary.py` | Aggregates sweep results across all frames/models into CSV + bar plots. |

## Formulations (active set 1–2, each run with/without ΔP)

| # | Slug (no ΔP / with ΔP) | Bending | Channel thickness h |
|---|------------------------|---------|---------------------|
| 1 | `viscous` / `viscous_wDynP` | Viscous, $\frac{2}{3}(H^3/R^3)\eta_L v_c$ | Fixed, 200 km |
| 2 | `plastic` / `plastic_wDynP` | Plastic, $\frac{1}{6}(H^2/R)\sigma_Y$ | Fixed, 200 km |

Canonical config (`configs/run_params.yaml`): formulation=1, DP_ref=2.35e7 Pa, include_ridge_push=1.
Use `--dp-ref 0` (or `DP_REF_VALUES=0` in make) for the no-ΔP variants.

## Output layout

```
plots/<hs3|nnr|sa>/<model>/param-sweep/   # misfit heatmaps only (withRP only)
plots/<hs3|nnr|sa>/<model>/best-fit/      # bestfit.*.txt, correlation.*.png, vt_param.*.png, map.*.png
plots/summary/                            # matrix_summary.csv, rmse_by_formulation.png, sign_match_by_formulation.png
plots/sketch/                             # hand-drawn / reference figures (tracked in git)
```

`matrix_summary.py` reads `bestfit.*.txt` from `best-fit/` directories (suite `best-fit`).

## Latest results (viscous_wDynP / F1+ΔP, DP_ref=2.35e7, with ridge push)

| Frame | Min RMSE | Sign match |
|-------|----------|------------|
| hs3   | 4.01 cm/yr | 70/120   |
| nnr   | 2.27 cm/yr | 94/120   |
| sa    | 3.08 cm/yr | 84/120   |

Results stale — rerun needed with corrected drag lengths (Lsp=col26, slabL=col8) and 4-model matrix.
Note: active segment count is 120 (not 98 as previously stated).

## Map data

Map plotting auto-detects:
- `data/age.3.6.NaN.grd` — oceanic plate age raster
- `data/PB2002_tdiddy.gmt` — Bird plate boundaries

Or set `DATASETS_DIR` to a root containing `age/age.3.6.NaN.grd` and `plate_boundaries/bird_PB2002/PB2002_tdiddy.gmt`. Maps run without these files, just without the background overlays.

Age grid and plate boundaries are cached in memory across map calls (module-level dicts in `plotting_functions.py`), so they are only read from disk once per process.

---

## Next steps

### ~~i) Matrix summary — cover all formulations and reference frames~~ ✓ done
### ~~i-b) Force budget map physics~~ ✓ resolved
### ~~iii) Vt vs. parameter plots~~ ✓ done

### A) Full rerun — regenerate all outputs with current code
Run `make run-matrix-with-summary` (12 runs: 2 formulations × 2 DP cases × 3 frames).
This will populate `plots/` with corrected drag lengths and produce the 4-bar summary figures.

### B) Dataset audit — which segments are included/excluded ✓ done
159 rows total → **120 included**, 39 excluded. Run `python dataset_audit.py --csv plots/audit/segment_audit.csv`.
- Primary exclusion reason: missing dip (col 6) — 35 segments (Philippines, Nankai, Mexico, S.Chile, Antilles)
- Also: missing slabL (13), slabD (5), age/Lsp (1 each)
- 41 segments had Rmin filled with global mean (243 km) before the filter
- All 120 included segments have observed vt in all 3 reference frames

### C) Collaborator handoff — final checks before passing to student
- Confirm end-to-end run from clean clone (`make run-matrix-with-summary`)
- Check `docs/force_balance.md` matches current code exactly
- Remove any remaining debug prints or stale comments
- Tag a clean release commit

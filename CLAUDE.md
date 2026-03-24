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
| `functions.py` | Physics engine — `compute_vsp_withDP()` implements all 4 formulations. |
| `workflow_common.py` | Shared data prep: VT column mapping, segment arrays, depth/dip interpolation. |
| `plotting_functions.py` | `save_quick_plot()`, `save_misfit_heatmap()`, `save_trench_motion_map()`, `save_vt_param_plot()`. |
| `matrix_summary.py` | Aggregates sweep results across all frames/models into CSV + bar plots. |

## Formulations (active set 1–4)

| # | Name (slug) | Bending | Channel thickness h |
|---|-------------|---------|---------------------|
| 1 | `viscous` | Viscous, $\frac{2}{3}(H^3/R^3)\eta_L v_c$ | Fixed, 200 km |
| 2 | `plastic` | Plastic, $\frac{1}{6}(H^2/R)\sigma_Y$ | Fixed, 200 km |
| 3 | `viscous_LspShear` | Viscous | $h \propto L_{sp}$ (100–250 km) |
| 4 | `viscous_VspShear` | Viscous | $h \propto v_{sp}$ (150–200 km), closed-form quadratic |

Canonical config (`configs/run_params.yaml`): formulation=1, DP_ref=2.35e7 Pa, include_ridge_push=1.

All formulations include dynamic pressure (DP) back-pressure parameterised from a reference Stokes solution. See `docs/force_balance.md` for full derivations.

## Output layout

```
plots/<hs3|nnr|sa>/<model>/param-sweep/   # misfit heatmaps only (withRP only)
plots/<hs3|nnr|sa>/<model>/best-fit/      # bestfit.*.txt, correlation.*.png, vt_param.*.png, map.*.png
plots/summary/                            # matrix_summary.csv, rmse_by_formulation.png, sign_match_by_formulation.png
plots/sketch/                             # hand-drawn / reference figures (tracked in git)
```

`matrix_summary.py` reads `bestfit.*.txt` from `best-fit/` directories (suite `best-fit`).

## Latest results (Formulation 1, viscous, DP_ref=2.35e7, with ridge push)

| Frame | Min RMSE | Sign match |
|-------|----------|------------|
| hs3   | 4.01 cm/yr | 70/120   |
| nnr   | 2.27 cm/yr | 94/120   |
| sa    | 3.08 cm/yr | 84/120   |

## Map data

Map plotting auto-detects:
- `data/age.3.6.NaN.grd` — oceanic plate age raster
- `data/PB2002_tdiddy.gmt` — Bird plate boundaries

Or set `DATASETS_DIR` to a root containing `age/age.3.6.NaN.grd` and `plate_boundaries/bird_PB2002/PB2002_tdiddy.gmt`. Maps run without these files, just without the background overlays.

Age grid and plate boundaries are cached in memory across map calls (module-level dicts in `plotting_functions.py`), so they are only read from disk once per process.

---

## Next steps

### ~~i) Matrix summary — cover all formulations and reference frames~~ ✓ done
Grouped bar plots (`rmse_by_formulation.png`, `sign_match_by_formulation.png`) in `plots/summary/`.

### i-b) Force budget map — finalise (physics resolved, rerun needed)

`save_force_budget_map()` in `plotting_functions.py`, called from `compute_rates_misfit.py` at best-fit params.
Output: `plots/<vt_ref>/<model>/best-fit/force_budget.*.png`

**Physics resolved.** Key facts from `docs/force_balance.md`:
- Physical balance: `F_sp + F_rp = F_bend + F_pdrag + F_sdrag + F_DP`
- `F_DP = η_A · C_DP · v_t` where `v_t = v_sp − v_c` (positive = trench retreat)
- `F_pdrag = 2η_A · (L/h) · v_sp` (plate drag, scales with **absolute plate velocity** v_sp, not v_c)
- `F_sdrag = η_A · (L_sp/h) · v_sp` (slab channel drag, scales with v_sp)
- For **retreating** zones (v_t > 0): F_DP resists (resisting half of pie)
- For **advancing** zones (v_t < 0): F_DP drives (driving half of pie)
- The `+η_A·C_DP·v_c` that appears in the numerator when solving for `v_sp` is an **algebraic artifact** of substituting `v_t = v_sp − v_c`, not a physical driving force

**TODO:** Rerun the sweep to regenerate force budget plots with the corrected plate drag (∝ v_sp).

### ii) Dataset audit — which trenches are included/excluded
Understand which of the ~160 Lallemand segments make it into the active set (~120) and why the others are dropped (missing dip, Rmin, age, vt, etc.). Useful outputs:
- A table or map of included vs. excluded segments
- Distribution plots of key parameters (age, dip, Rmin, Lsp, slabD) for included segments
- Flag segments where observed vt is missing in one or more reference frames

### ~~iii) Vt vs. parameter plots (model curves + data)~~ ✓ done
`save_vt_param_plot()` in `plotting_functions.py`. Outputs to `plots/<vt_ref>/<model>/best-fit/vt_param_*.png`.

### iv) Collaborator handoff — clean up and document
- Confirm all scripts run end-to-end from a clean clone (`make run-matrix-with-summary`)
- Add a brief methods note or caption-ready description for each plot type
- Check that `docs/force_balance.md` is self-contained and matches the current code exactly
- Remove any remaining debug prints or stale comments
- Tag a clean release commit once the above are done

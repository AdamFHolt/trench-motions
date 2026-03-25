# trench-motions 

This repository contains analytical subduction/trench-motion modeling scripts and supporting datasets.

## Start Here 
Run from repository root.

1. Setup environment:
```bash
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```

2. Run the canonical sweep + summary:
```bash
make run-matrix-with-summary
```

3. Look at results:
- Sweep products (heatmaps): `plots/<hs3|nnr|sa>/<model>/param-sweep/`
- Best-fit outputs (vt-param plots, correlation, force-budget map, prediction table): `plots/<hs3|nnr|sa>/<model>/best-fit/`
- Summary table and bar plots: `plots/summary/`

4. Optional sweep only (no summary step):
```bash
make run-matrix
```

5. Optional map-producing run:
```bash
make run-matrix-maps-with-summary
```

## What Is Active
- `compute_rates_misfit.py`: parameter sweep and best-fit selection.
- `compute_rates_single.py`: single-parameter model run.
- `plotting_functions.py`: shared plotting utilities for:
  - misfit heatmaps
  - quick observed-vs-predicted diagnostics
  - full trench-motion maps

## Project Layout
- Inputs: `data/` and `data/vt/`.
- Active outputs:
  - `plots/<hs3|nnr|sa>/<model>/param-sweep/` — misfit heatmaps.
  - `plots/<hs3|nnr|sa>/<model>/best-fit/` — prediction tables, vt-param plots, correlation, force-budget map, global map.
  - `plots/summary/` — summary CSV and grouped bar plots.

## Runtime notes
- Active scripts are Python 3 compatible.
- Python dependencies: `numpy`, `scipy`, `matplotlib`.
- Map plotting uses the same Python dependencies listed in `requirements.txt`.
- In sweep runs, a quick diagnostic plot is generated for the best-fit case in every reference frame.

## Core Commands
From repository root.

Main workflow:
```bash
make run-matrix-with-summary
```

Map workflow:
```bash
make run-matrix-maps-with-summary
```

Single direct run:
```bash
python3 compute_rates_single.py sa 1 23.5e6 1e21 1e22
```

For full CLI options and arguments:
```bash
python3 compute_rates_misfit.py --help
python3 compute_rates_single.py --help
```


## Datasets for map plotting
Map plotting uses the same dataset layout controls via `DATASETS_DIR`.

Set `DATASETS_DIR` to a dataset root containing:
- `age/age.3.6.NaN.grd`
- `plate_boundaries/bird_PB2002/PB2002_tdiddy.gmt`

Example:
```bash
export DATASETS_DIR=/path/to/datasets
```

If `DATASETS_DIR` is not set, map plotting auto-detects a flat layout in `data/`:
- `data/age.3.6.NaN.grd`
- `data/PB2002_tdiddy.gmt`
If these datasets are missing, map plotting still runs, but without the age raster background and PB2002 boundary overlay.

`PB2002_tdiddy.gmt` is tracked in this repository. `age.3.6.NaN.grd` is not (25 MB); obtain it from the EarthByte/Müller et al. oceanic age grid distribution and place it at `data/age.3.6.NaN.grd`.

## Future improvements

- **Update subduction zone parameters from newer sources** (e.g. slab dip, depth, length). The Lallemand et al. (2005) dataset excludes 39 of 159 segments due to missing parameters — a more complete modern compilation would increase coverage and reduce this exclusion rate.
- **Refine plate and trench geometry estimates** — in particular the ridge-to-trench distance ($L_{sp}$, col 26) and the along-strike trench width ($w$, col 25), which were estimated manually. These could be recalculated systematically from up-to-date plate boundary and spreading ridge datasets.

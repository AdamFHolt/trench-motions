# trench-motions (subduction rates)

This repository contains analytical subduction/trench-motion modeling scripts and supporting datasets.

## Start Here (Student Workflow)
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
- Sweep products: `plots/<hs3|nnr|sa>/param-sweep/...`
- Sweep summary table: `plots/summary/param-sweep/matrix_summary.csv`
- Sweep summary plots: `plots/summary/param-sweep/*.png`

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
- `quick_plot.py`: fast diagnostic observed-vs-predicted check.
- `plot_trench_motions.py`: map plotting (used when `--skip-map` is not set).

## Project Layout
- Inputs: `data/` and `data/vt/`.
- Active outputs:
  - `plots/hs3/param-sweep/`, `plots/nnr/param-sweep/`, `plots/sa/param-sweep/` for parameter sweeps.
  - `plots/hs3/maps/`, `plots/nnr/maps/`, `plots/sa/maps/` for map-producing runs.
  - `plots/summary/param-sweep/` and `plots/summary/maps/` for summary products.
  - direct script runs without `--out-prefix` write under `plots/<vt_ref>/maps/`.
- Archived generated outputs: `archive/generated/`.
- Archived reference files: `archive/reference/`.
- Archived legacy scripts: `archive/legacy/`.

## Runtime notes
- Active scripts are Python 3 compatible.
- Python dependencies: `numpy`, `scipy`, `matplotlib`.
- Map plotting dependencies: Python packages already in this repo (`numpy`, `scipy`, `matplotlib`).

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
python3 compute_rates_single.py sa 1 1 23.5e6 1e21 1e22
```

For full CLI options and arguments:
```bash
python3 compute_rates_misfit.py --help
python3 compute_rates_single.py --help
```

## Datasets for map plotting
`plot_trench_motions.py` accepts the same dataset layout controls via `DATASETS_DIR`.

Set `DATASETS_DIR` to a dataset root containing:
- `age/age.3.6.NaN.grd`
- `plate_boundaries/bird_PB2002/PB2002_tdiddy.gmt`

Example:
```bash
export DATASETS_DIR=/path/to/datasets
```

If `DATASETS_DIR` is not set, the script looks in `data/gmt_datasets`.
It also auto-detects a flat layout in `data/`:
- `data/age.3.6.NaN.grd`
- `data/PB2002_tdiddy.gmt`
If these datasets are missing, map plotting still runs, but without the age raster background and PB2002 boundary overlay.

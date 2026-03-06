# trench-motions (subduction rates)

This repository contains analytical subduction/trench-motion modeling scripts and supporting datasets.

## Start Here (Student Workflow)
Run from repository root.

1. Setup environment:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

2. Run the canonical sweep + summary:
```bash
make run-matrix-with-summary
```

3. Look at results:
- Sweep products: `results/param-sweep/...`
- Summary table: `results/param-sweep/summary/matrix_summary.csv`
- Summary plots: `results/param-sweep/summary/*.png`

4. Optional sweep only (no summary step):
```bash
make run-matrix
```

5. Optional map-producing run (requires GMT datasets):
```bash
export DATASETS_DIR=/path/to/datasets
make run-matrix-maps-with-summary
```

## What Is Active
- `compute_rates_misfit.py`: parameter sweep and best-fit selection.
- `compute_rates_single.py`: single-parameter model run.
- `quick_plot.py`: fast non-GMT diagnostic map.
- `create_trench_motion_table.py`: rebuilds `data/vts_hs3-nnr-sa.txt` from `data/vt/*.dat`.
- `plot_trench_motions.sh`: GMT map plotting (used only when `--skip-map` is not set).

## Project Layout
- Inputs: `data/` and `data/vt/`.
- Active outputs:
  - `results/param-sweep/` for parameter sweeps + summaries.
  - `results/maps/` for map-producing matrix runs.
  - `results/one-off/` for direct one-off script runs without `--out-prefix`.
- Archived generated outputs: `archive/generated/`.
- Archived reference files: `archive/reference/`.
- Archived legacy scripts: `archive/legacy/`.

## Runtime notes
- Active scripts are Python 3 compatible.
- Python dependencies: `numpy`, `scipy`, `matplotlib`.
- Map plotting dependencies: GMT tools (`gmtset`, `makecpt`, `psxy`, `grdview`, etc.), `eps2eps`, and ImageMagick `convert`.

## Additional Commands
From repository root.

Show CLI help:
```bash
python3 compute_rates_misfit.py --help
python3 compute_rates_single.py --help
python3 create_trench_motion_table.py --help
```

Environment setup via Makefile:
```bash
make venv
make install
```

Example misfit sweep:
```bash
python3 compute_rates_misfit.py sa 1 1 23.5e6 1
```

Smoke run (fast grid, no GMT map plotting):
```bash
python3 compute_rates_misfit.py sa 1 1 23.5e6 0 --smoke --skip-map
```

Isolated outputs under a custom directory:
```bash
python3 compute_rates_misfit.py sa 1 1 23.5e6 0 --smoke --skip-map --out-prefix results/exp1
```

Equivalent via Makefile:
```bash
make smoke
```

Config-driven smoke:
```bash
python3 compute_rates_misfit.py --config configs/misfit_smoke.yaml
```

Equivalent via Makefile:
```bash
make smoke-config
```

Example single run:
```bash
python3 compute_rates_single.py sa 1 1 23.5e6 1e21 1e22
```

Single-run smoke (skip GMT map plotting):
```bash
python3 compute_rates_single.py sa 1 1 23.5e6 1e21 1e22 --skip-map
```

Isolated single-run outputs:
```bash
python3 compute_rates_single.py sa 1 1 23.5e6 1e21 1e22 --skip-map --out-prefix results/exp1_single
```

Equivalent via Makefile:
```bash
make single-smoke
```

Config-driven single smoke:
```bash
python3 compute_rates_single.py --config configs/single_smoke.yaml
```

Equivalent via Makefile:
```bash
make single-smoke-config
```

Matrix runs (config-driven):
```bash
make run-matrix
```

Map-producing matrix runs (requires GMT datasets):
```bash
export DATASETS_DIR=/path/to/datasets
make run-matrix-maps
```

Shared matrix config files:
- `configs/matrix.yaml` (skip-map true)
- `configs/matrix_maps.yaml` (skip-map false)

Batch matrix summary tables/plots from `results/param-sweep/*`:
```bash
make matrix-summary
```

One-shot run + summary:
```bash
make run-matrix-with-summary
```

One-shot map run + summary:
```bash
export DATASETS_DIR=/path/to/datasets
make run-matrix-maps-with-summary
```

One-off quick plot (no GMT) from existing files:
```bash
make quick-plot
```

Or explicitly:
```bash
python3 quick_plot.py --predicted <predicted_txt> --observed data/vt/tnew.sa.dat --output <output_png>
```

## Argument Reference
- `compute_rates_misfit.py` positional args:
  `vt_ref formulation include_DP DP_ref include_ridge_push`
- `compute_rates_misfit.py` optional flags:
  `--smoke --skip-map --out-prefix <dir> --config <path.yaml> --vt-ref <hs3|nnr|sa>`
- `compute_rates_single.py` positional args:
  `vt_ref formulation include_DP DP_ref asthen_visc lith_visc`
- `compute_rates_single.py` optional flags:
  `--skip-map --out-prefix <dir> --config <path.yaml>`

## GMT datasets for map plotting
`plot_trench_motions.sh` no longer uses a hardcoded absolute path.

Set `DATASETS_DIR` to a dataset root containing:
- `age/age.3.6.NaN.grd`
- `plate_boundaries/bird_PB2002/PB2002_tdiddy.gmt`

Example:
```bash
export DATASETS_DIR=/path/to/datasets
```

If `DATASETS_DIR` is not set, the script looks in `data/gmt_datasets`.

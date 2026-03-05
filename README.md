# trench-motions (subduction rates)

This repository contains analytical subduction/trench-motion modeling scripts and supporting datasets.

## Main scripts
- `compute_rates_misfit.py`: parameter sweep over viscosity/prefactor/yield settings and selection of best-fit models.
- `compute_rates_single.py`: single-parameter model run.
- `functions.py`: core analytical model functions.
- `create_trench_motion_table.py`: rebuilds `data/vts_hs3-nnr-sa.txt` from `data/vt/*.dat`.
- `plot_trench_motions.sh`: GMT-based map plotting for observed vs predicted trench motions.

## Active vs legacy
- Active entrypoints:
  - `compute_rates_misfit.py`
  - `compute_rates_single.py`
  - `create_trench_motion_table.py`
  - `plot_trench_motions.sh`
- Legacy/deprecated:
  - `compute_rate_plots.py`
  - `compute_rate_plots_vtvp.py`
  - `old/` (historical versions)

## Data and outputs
- Inputs: `data/` and `data/vt/`.
- Model outputs: `predictions/`.
- Figure outputs: `plots/`.
- Temporary text outputs: `tmp/` (created automatically by active scripts).

## Runtime notes
- Active scripts are Python 3 compatible.
- Python dependencies: `numpy`, `scipy`, `matplotlib`.
- Map plotting dependencies: GMT tools (`gmtset`, `makecpt`, `psxy`, `grdview`, etc.), `eps2eps`, and ImageMagick `convert`.

## Environment Setup
Quick setup:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Equivalent via Makefile:
```bash
make venv
make install
```

## Running
From the repository root.

Show CLI help:
```bash
python3 compute_rates_misfit.py --help
python3 compute_rates_single.py --help
python3 create_trench_motion_table.py --help
```

Example misfit sweep:
```bash
python3 compute_rates_misfit.py sa 4 1 23.5e6 1e-13 0.25 1
```

Smoke run (fast grid, no GMT map plotting):
```bash
python3 compute_rates_misfit.py sa 1 1 23.5e6 1e-13 0.25 0 --smoke --skip-map
```

Equivalent via Makefile:
```bash
make smoke
```

Example single run:
```bash
python3 compute_rates_single.py sa 4 1 23.5e6 1e-13 0.25 1e21 1e22
```

Single-run smoke (skip GMT map plotting):
```bash
python3 compute_rates_single.py sa 1 1 23.5e6 1e-13 0.25 1e21 1e22 --skip-map
```

Equivalent via Makefile:
```bash
make single-smoke
```

Manual quick plot (no GMT) from existing files:
```bash
make quick-plot
```

Or explicitly:
```bash
python3 quick_plot.py --predicted <predicted_txt> --observed data/vt/tnew.sa.dat --output <output_png>
```

## Argument Reference
- `compute_rates_misfit.py` positional args:
  `vt_ref formulation include_DP DP_ref trans_strain_rate PSP_slab_pull_factor include_ridge_push`
- `compute_rates_misfit.py` optional flags:
  `--smoke --skip-map`
- `compute_rates_single.py` positional args:
  `vt_ref formulation include_DP DP_ref trans_strain_rate PSP_slab_pull_factor asthen_visc lith_visc`
- `compute_rates_single.py` optional flags:
  `--skip-map`

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

# trench-motions (subduction rates)

This repository contains analytical subduction/trench-motion modeling scripts and supporting datasets.

## Main scripts
- `compute_rates_misfit.py`: parameter sweep over viscosity/prefactor/yield settings and selection of best-fit models.
- `compute_rates_single.py`: single-parameter model run.
- `functions.py`: core analytical model functions.
- `create_trench_motion_table.py`: rebuilds `data/vts_hs3-nnr-sa.txt` from `data/vt/*.dat`.
- `plot_trench_motions.sh`: GMT-based map plotting for observed vs predicted trench motions.

## Data and outputs
- Inputs: `data/` and `data/vt/`.
- Model outputs: `predictions/`.
- Figure outputs: `plots/`.
- Temporary text outputs: `tmp/` (created automatically by active scripts).

## Runtime notes
- Code uses Python 2 syntax (`print ...`) and was developed in a legacy environment.
- Python dependencies: `numpy`, `scipy`, `matplotlib`.
- Map plotting dependencies: GMT tools (`gmtset`, `makecpt`, `psxy`, `grdview`, etc.), `eps2eps`, and ImageMagick `convert`.

## Running
From the repository root.

Example misfit sweep:
```bash
python compute_rates_misfit.py sa 4 1 23.5e6 1e-13 0.25 1
```

Example single run:
```bash
python compute_rates_single.py sa 4 1 23.5e6 1e-13 0.25 1e21 1e22
```

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

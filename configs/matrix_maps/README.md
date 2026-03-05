# Matrix Map Configs

This config set mirrors `configs/matrix/` but enables GMT map generation.

Differences:
- `skip_map: false`
- `out_prefix` points to `runs/matrix_maps/...`

Prerequisite for running:
- set `DATASETS_DIR` for `plot_trench_motions.sh`
  (must contain `age/age.3.6.NaN.grd` and `plate_boundaries/bird_PB2002/PB2002_tdiddy.gmt`)

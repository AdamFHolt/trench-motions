# Data directory — trench-motions

## Files

| File | Description |
|------|-------------|
| `Lallemand_et_al-2005_G3_dataset_WithAdditionalParams.txt` | Main model input — 159 subduction segments, whitespace-delimited, no header |
| `Lallemand_et_al-2005_G3_dataset_WithAdditionalParams.xlsx` | Spreadsheet version with column headers; authoritative reference for column meanings |
| `vt/tnew.hs3.dat` | Observed trench velocities, HS3 reference frame (174 rows) |
| `vt/tnew.nnr.dat` | Observed trench velocities, NNR reference frame |
| `vt/tnew.sa.dat` | Observed trench velocities, South America reference frame |
| `vt/details.txt` | Notes on the observed velocity files |
| `age.3.6.NaN.grd` | Oceanic plate age raster (Müller et al.); used as map background only — **not tracked in repo (25 MB); download separately** |
| `PB2002_tdiddy.gmt` | Bird (2002) plate boundaries in GMT 2-point segment format; used as map background only |

---

## Main dataset: `WithAdditionalParams.txt`

Derived from Lallemand et al. (2005) G-Cubed Table 1, with two columns appended.
Read in Python as `np.genfromtxt(path)` — the segment-name string in column 0
becomes NaN; all other columns are floats. **Column indices below are 0-based**
(as used throughout the code).

### Columns from Lallemand / Heuret & Lallemand (2005, 2007)

| Col (0-based) | Symbol | Units | Description |
|---|---|---|---|
| 0 | — | — | Segment name (string; read separately with `dtype=str, usecols=0`) |
| 1 | φ | °N | Latitude of segment midpoint |
| 2 | λ | °E | Longitude of segment midpoint |
| 3 | — | — | Original Lallemand column — see xlsx for definition |
| 4 | Az | ° | Azimuth of convergence |
| 5 | δ₁ | ° | Slab dip at shallow depth (< ~125 km) |
| 6 | δ₂ | ° | Slab dip at deeper level (> ~125 km); **used as `dip` in model** |
| 7 | D_max | km | Maximum depth of seismicity (slab depth); **used as `slabD`** |
| 8 | L | km | Measured slab length along dip; **used as `slabL` (slab drag)** |
| 9 | v_c | mm/yr | Trench-normal convergence velocity; code reads as `data[:,9]/10` → cm/yr |
| 10–11 | v_T | mm/yr | Trench velocity components in Lallemand reference frames (not used; model uses `tnew.*.dat`) |
| 12 | — | — | Original Lallemand column — see xlsx |
| 13 | R_min | km | Minimum radius of curvature of slab; **used as `Rmin`**. NaN for 41 segments — filled with global mean (243 km) before filtering |
| 14–19 | v_T | mm/yr | Additional trench/plate velocity components in various frames (not used in model) |
| 20 | age | Ma | Age of subducting oceanic plate; **used for thermal structure (H, buoyancy, ridge push)** |
| 21–22 | — | — | Original Lallemand columns — see xlsx |
| 23 | — | — | Continental (c) / oceanic (o) flag (string → NaN when read as float; not used) |
| 24 | — | — | Original Lallemand flag — see xlsx |

### Columns added/estimated by Holt

| Col (0-based) | Symbol | Units | Description |
|---|---|---|---|
| 25 | w | km | Along-strike width of the subducting plate (trench length); **used in ΔP scaling** |
| 26 | Lsp | km | Ridge-to-trench distance (length of subducting plate); **used as `Lsp` (plate drag)** |

Both are assigned per subduction system (all segments in the same system share the
same value). See the xlsx for notes on how each value was estimated.


---

## Inclusion filter (applied in `workflow_common.py`)

A segment is included in the model if **all** of the following columns are non-NaN
**after** preprocessing:

| Col | Variable | Notes |
|-----|----------|-------|
| 6 | dip (δ₂) | Not interpolated from δ₁ (option disabled) |
| 7 | slabD | — |
| 8 | slabL | — |
| 13 | Rmin | NaN values filled with mean (243 km) **before** this check |
| 20 | age | — |
| 26 | Lsp | — |

**Result: 120 of 159 segments included.** The 39 excluded segments fail primarily
on missing dip (35 segments: Philippines, Nankai, Mexico, S. Chile, Antilles),
followed by missing slabL (13) and slabD (5). All 120 included segments have
observed v_T in all three reference frames.

---

## Observed velocity files: `vt/tnew.*.dat`

Nine-column whitespace-delimited files, one row per segment:

| Col (0-based) | Description |
|---|---|
| 0 | Longitude (°E) |
| 1 | Latitude (°N) |
| 2 | Tangential azimuth (°) |
| 3 | vE_trench — eastward trench velocity (mm/yr) |
| 4 | vN_trench — northward trench velocity (mm/yr) |
| 5 | vnormal_trench — trench-normal velocity v_T (mm/yr); **used by model**; positive = trench retreat |
| 6 | vE_plate — eastward plate velocity (mm/yr) |
| 7 | vN_plate — northward plate velocity (mm/yr) |
| 8 | vnnormal_plate — plate-normal velocity (mm/yr) |

Segments are matched to the main dataset by rounded (lat, lon) coordinates
(tolerance 0.1°). All 120 model-included segments are present in all three files.

---

## Reference frames

| Label | Description |
|-------|-------------|
| `hs3` | Hot-spot reference frame (HS3-NUVEL1A) |
| `nnr` | No-net-rotation reference frame |
| `sa`  | South America fixed reference frame |

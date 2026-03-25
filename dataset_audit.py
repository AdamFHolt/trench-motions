#!/usr/bin/env python3
"""
Dataset audit: which Lallemand segments are included/excluded and why.

Usage:
    python dataset_audit.py [--csv <path>]

Mirrors the preprocessing in compute_rates_misfit.py exactly:
    use_avg_Rmin = 1, interpolate_shallow_dip = 0

Outputs:
    - Console: summary + per-segment table
    - CSV (optional): one row per segment with status and flags
"""

import os, sys
import numpy as np

DATA_PATH = 'data/Lallemand_et_al-2005_G3_dataset_WithAdditionalParams.txt'
VT_FRAMES = ['hs3', 'nnr', 'sa']
VT_COLS   = {'hs3': 2, 'nnr': 3, 'sa': 4}   # col in data_vt where vt is stored

FILTER_COLS = {
    'dip':   6,
    'slabD': 7,
    'slabL': 8,
    'Rmin':  13,
    'age':   20,
    'Lsp':   26,
}

COL_LABELS = {
    6:  'dip (col 6)',
    7:  'slabD (col 7)',
    8:  'slabL (col 8)',
    13: 'Rmin (col 13)',
    20: 'age (col 20)',
    26: 'Lsp (col 26)',
}


def load_vt_table(data, vt_ref, tol=0.1):
    """Return array[n_rows] of observed vt (NaN where missing)."""
    path = 'data/vt/tnew.{}.dat'.format(vt_ref)
    tnew = np.genfromtxt(path)
    if tnew.ndim == 1:
        tnew = tnew.reshape(1, -1)
    decimals = max(0, int(round(-np.log10(tol)))) if tol < 1 else 0
    vt_by_coord = {}
    for i in range(tnew.shape[0]):
        lat = round(float(tnew[i, 1]), decimals)
        lon = round(float(tnew[i, 0]), decimals)
        vt_by_coord[(lat, lon)] = float(tnew[i, 5])
    out = np.full(data.shape[0], np.nan)
    for i in range(data.shape[0]):
        lat = round(float(data[i, 1]), decimals)
        lon = round(float(data[i, 2]), decimals)
        if (lat, lon) in vt_by_coord:
            out[i] = vt_by_coord[(lat, lon)]
    return out


def main():
    csv_path = None
    args = sys.argv[1:]
    if '--csv' in args:
        idx = args.index('--csv')
        csv_path = args[idx + 1]

    # --- load dataset (mirrors compute_rates_misfit.py) ---
    data = np.genfromtxt(DATA_PATH)
    names = np.genfromtxt(DATA_PATH, dtype=str, usecols=0)
    n_total = data.shape[0]

    # preprocessing: fill NaN Rmin with global mean (use_avg_Rmin=1)
    rmin_mean = float('%.2f' % np.nanmean(data[:, 13]))
    nan_rmin = np.where(np.isnan(data[:, 13]))[0]
    data_proc = data.copy()
    data_proc[nan_rmin, 13] = rmin_mean

    # load observed vt for all frames
    vt = {ref: load_vt_table(data, ref) for ref in VT_FRAMES}

    # --- classify each segment ---
    rows = []
    for i in range(n_total):
        missing = []
        for key, col in FILTER_COLS.items():
            if np.isnan(data_proc[i, col]):
                missing.append(key)
        included = len(missing) == 0

        rmin_was_filled = (i in nan_rmin)

        vt_available = {ref: not np.isnan(vt[ref][i]) for ref in VT_FRAMES}

        rows.append({
            'name':           names[i],
            'lat':            data[i, 1],
            'lon':            data[i, 2],
            'included':       included,
            'missing':        missing,
            'rmin_filled':    rmin_was_filled,
            'vt_hs3':         vt_available['hs3'],
            'vt_nnr':         vt_available['nnr'],
            'vt_sa':          vt_available['sa'],
        })

    # --- summary ---
    n_included = sum(r['included'] for r in rows)
    n_excluded = n_total - n_included

    print("=" * 60)
    print("DATASET AUDIT: {}".format(DATA_PATH))
    print("=" * 60)
    print("Total segments:    {:3d}".format(n_total))
    print("Included:          {:3d}".format(n_included))
    print("Excluded:          {:3d}".format(n_excluded))
    print()
    print("Rmin filled with mean ({:.0f} km) for {:d} segments before filter".format(
        rmin_mean, len(nan_rmin)))
    print()

    # exclusion breakdown by missing column
    from collections import Counter
    excl_counter = Counter()
    for r in rows:
        if not r['included']:
            for m in r['missing']:
                excl_counter[m] += 1
    print("Exclusion reasons (segments can fail multiple criteria):")
    for key in FILTER_COLS:
        if excl_counter[key]:
            print("  {:8s}  {:3d} segments".format(key, excl_counter[key]))
    print()

    # vt coverage for included segments
    print("Observed vt coverage (included segments only):")
    for ref in VT_FRAMES:
        n_have = sum(r['vt_{}'.format(ref)] for r in rows if r['included'])
        n_miss = n_included - n_have
        print("  {:4s}  {:3d} have vt   {:3d} missing".format(ref, n_have, n_miss))
    print()

    # --- per-segment table ---
    print("{:<10s}  {:>6s}  {:>7s}  {:^8s}  {:^5s}{:^5s}{:^5s}  {}".format(
        "Segment", "Lat", "Lon", "Status", "hs3", "nnr", "sa", "Missing"))
    print("-" * 72)
    for r in rows:
        status = "OK" if r['included'] else "EXCL"
        vt_str = "".join(
            (" Y  " if r['vt_{}'.format(ref)] else " -  ")
            for ref in VT_FRAMES
        )
        miss_str = ",".join(r['missing']) if r['missing'] else ""
        if r['rmin_filled'] and r['included']:
            miss_str = "(Rmin filled)"
        print("{:<10s}  {:>6.1f}  {:>7.1f}  {:^8s}  {}  {}".format(
            r['name'], r['lat'], r['lon'], status, vt_str, miss_str))

    # --- CSV ---
    if csv_path:
        os.makedirs(os.path.dirname(csv_path) or '.', exist_ok=True)
        with open(csv_path, 'w') as f:
            f.write("name,lat,lon,included,missing_cols,rmin_filled,vt_hs3,vt_nnr,vt_sa\n")
            for r in rows:
                f.write("{},{:.1f},{:.1f},{},{},{},{},{},{}\n".format(
                    r['name'], r['lat'], r['lon'],
                    int(r['included']),
                    "|".join(r['missing']),
                    int(r['rmin_filled']),
                    int(r['vt_hs3']), int(r['vt_nnr']), int(r['vt_sa']),
                ))
        print("\nCSV written to {}".format(csv_path))


if __name__ == '__main__':
    main()

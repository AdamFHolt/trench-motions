#!/usr/bin/env python3

import argparse
import csv
import glob
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def load_2d(path):
    arr = np.genfromtxt(path)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr


def detect_vt_ref_from_name(name):
    lower = name.lower()
    if 'hs3model' in lower:
        return 'hs3'
    if 'nnrmodel' in lower:
        return 'nnr'
    if 'samodel' in lower:
        return 'sa'
    return None


def build_observed_map(obs_path, decimals):
    obs = load_2d(obs_path)
    if obs.shape[1] < 6:
        raise ValueError('Observed file has unexpected format: {}'.format(obs_path))
    obs_lon = obs[:, 0]
    obs_lat = obs[:, 1]
    obs_vt = obs[:, 5] / 10.0  # mm/yr -> cm/yr

    out = {}
    for i in range(len(obs_vt)):
        key = (round(float(obs_lat[i]), decimals), round(float(obs_lon[i]), decimals))
        out[key] = float(obs_vt[i])
    return out


def compute_metrics(pred_path, obs_map, decimals, neutral):
    pred = load_2d(pred_path)
    if pred.shape[1] < 3:
        raise ValueError('Predicted file has unexpected format: {}'.format(pred_path))

    pred_lat = pred[:, 0]
    pred_lon = pred[:, 1]
    pred_vt = pred[:, 2]

    obs_vals = []
    pred_vals = []
    for i in range(len(pred_vt)):
        key = (round(float(pred_lat[i]), decimals), round(float(pred_lon[i]), decimals))
        if key in obs_map:
            obs_vals.append(obs_map[key])
            pred_vals.append(float(pred_vt[i]))

    if not obs_vals:
        raise ValueError('No coordinate matches found for {}'.format(pred_path))

    obs_vals = np.array(obs_vals)
    pred_vals = np.array(pred_vals)
    residual = pred_vals - obs_vals

    rmse = float(np.sqrt(np.mean(residual ** 2)))
    mae = float(np.mean(np.abs(residual)))
    bias = float(np.mean(residual))

    sign_match = (
        np.sign(pred_vals) == np.sign(obs_vals)
    ) | (
        (np.abs(pred_vals) <= neutral) & (np.abs(obs_vals) <= neutral)
    )
    sign_count = int(np.sum(sign_match))
    sign_pct = 100.0 * float(sign_count) / float(len(sign_match))

    return {
        'n': int(len(obs_vals)),
        'rmse_cm_yr': rmse,
        'mae_cm_yr': mae,
        'bias_cm_yr': bias,
        'sign_match_count': sign_count,
        'sign_match_pct': sign_pct,
    }


def plot_metric_bar(rows, metric_key, ylabel, output_path, title):
    suites = [r['suite'] for r in rows]
    vt_refs = [r['vt_ref'] for r in rows]
    if len(set(vt_refs)) == 1:
        labels = suites
    elif len(set(suites)) == 1:
        labels = vt_refs
    else:
        labels = ['{}|{}'.format(s, v) for s, v in zip(suites, vt_refs)]
    vals = [float(r[metric_key]) for r in rows]
    x = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(max(8, 1.25 * len(rows)), 4.8), constrained_layout=True)
    bars = ax.bar(x, vals, color='#4C78A8', edgecolor='#1F2A44', linewidth=0.8, alpha=0.92)
    ax.set_axisbelow(True)
    ax.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha='right')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    ax.margins(y=0.1)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Summarize matrix run outputs into tables and plots.')
    parser.add_argument('--runs-dir', default='plots', help='Root directory containing per-reference outputs.')
    parser.add_argument('--vt-ref', default='', help='Reference frame override (hs3|nnr|sa).')
    parser.add_argument('--suites', default='param-sweep', help='Comma-separated suite names to summarize (e.g., param-sweep,maps).')
    parser.add_argument('--output-dir', default='plots/summary', help='Directory for summary CSV/plots.')
    parser.add_argument('--tol', type=float, default=0.1, help='Coordinate rounding tolerance in degrees.')
    parser.add_argument('--neutral', type=float, default=0.3, help='Neutral sign threshold in cm/yr.')
    args = parser.parse_args()

    ensure_dir(args.output_dir)
    decimals = max(0, int(round(-np.log10(args.tol)))) if args.tol < 1 else 0

    suites = [s.strip() for s in args.suites.split(',') if s.strip()]
    if not suites:
        raise SystemExit('At least one suite is required (--suites).')

    pred_files = []
    file_suite = {}
    for suite in suites:
        suite_files = sorted(glob.glob(os.path.join(args.runs_dir, '*', suite, 'rms_*.txt')))
        for p in suite_files:
            pred_files.append(p)
            file_suite[p] = suite
    if not pred_files:
        raise SystemExit('No matrix prediction files found under {} (suites={})'.format(args.runs_dir, args.suites))

    rows = []
    for pred_path in pred_files:
        rel = os.path.relpath(pred_path, args.runs_dir)
        parts = rel.split(os.sep)
        suite = file_suite.get(pred_path, parts[1] if len(parts) > 1 else 'unknown')
        vt_ref = args.vt_ref if args.vt_ref else (parts[0] if len(parts) > 0 else detect_vt_ref_from_name(os.path.basename(pred_path)))
        if vt_ref not in ('hs3', 'nnr', 'sa'):
            vt_ref = detect_vt_ref_from_name(os.path.basename(pred_path))
        if vt_ref not in ('hs3', 'nnr', 'sa'):
            raise ValueError('Could not determine vt_ref for {}'.format(pred_path))

        obs_path = os.path.join('data', 'vt', 'tnew.{}.dat'.format(vt_ref))
        obs_map = build_observed_map(obs_path, decimals)
        metrics = compute_metrics(pred_path, obs_map, decimals, args.neutral)

        row = {
            'suite': suite,
            'vt_ref': vt_ref,
            'predicted_file': pred_path,
        }
        row.update(metrics)
        rows.append(row)

    rows.sort(key=lambda r: (r['suite'], r['vt_ref']))

    csv_path = os.path.join(args.output_dir, 'matrix_summary.csv')
    fields = [
        'suite',
        'vt_ref',
        'n',
        'rmse_cm_yr',
        'mae_cm_yr',
        'bias_cm_yr',
        'sign_match_count',
        'sign_match_pct',
        'predicted_file',
    ]
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    plot_metric_bar(
        rows,
        metric_key='rmse_cm_yr',
        ylabel='RMSE (cm/yr)',
        title='RMSE by Run',
        output_path=os.path.join(args.output_dir, 'rmse_by_run.png'),
    )
    plot_metric_bar(
        rows,
        metric_key='sign_match_count',
        ylabel='Correct Sign Locations (#)',
        title='Correct Sign Locations by Run',
        output_path=os.path.join(args.output_dir, 'sign_match_by_run.png'),
    )

    print('Wrote {}'.format(csv_path))
    print('Wrote {}'.format(os.path.join(args.output_dir, 'rmse_by_run.png')))
    print('Wrote {}'.format(os.path.join(args.output_dir, 'sign_match_by_run.png')))


if __name__ == '__main__':
    main()

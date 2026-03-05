#!/usr/bin/env python3

import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent and not os.path.isdir(parent):
        os.makedirs(parent)


def load_2d(path):
    arr = np.genfromtxt(path)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr


def main():
    parser = argparse.ArgumentParser(description='Quick observed-vs-predicted trench motion plot (no GMT).')
    parser.add_argument('--predicted', required=True, help='Predicted file path (lat lon vt_cm_yr azim).')
    parser.add_argument('--observed', required=True, help='Observed file path (tnew.*.dat).')
    parser.add_argument('--output', required=True, help='Output PNG path.')
    parser.add_argument('--title', default='Quick Trench Motion Check', help='Figure title.')
    parser.add_argument('--tol', type=float, default=0.1, help='Coordinate rounding tolerance in degrees.')
    parser.add_argument('--neutral', type=float, default=0.3, help='Neutral sign threshold in cm/yr.')
    args = parser.parse_args()

    pred = load_2d(args.predicted)
    obs = load_2d(args.observed)

    if pred.shape[1] < 3:
        raise ValueError('Predicted file must have at least 3 columns: lat lon vt.')
    if obs.shape[1] < 6:
        raise ValueError('Observed file must have at least 6 columns (expected tnew.*.dat layout).')

    pred_lat = pred[:, 0]
    pred_lon = pred[:, 1]
    pred_vt = pred[:, 2]

    # tnew.*.dat uses lon, lat and vT in column 6 (mm/yr); convert to cm/yr.
    obs_lon = obs[:, 0]
    obs_lat = obs[:, 1]
    obs_vt = obs[:, 5] / 10.0

    decimals = max(0, int(round(-np.log10(args.tol)))) if args.tol < 1 else 0
    obs_map = {}
    for i in range(len(obs_vt)):
        key = (round(float(obs_lat[i]), decimals), round(float(obs_lon[i]), decimals))
        obs_map[key] = float(obs_vt[i])

    x_obs = []
    y_pred = []
    for i in range(len(pred_vt)):
        key = (round(float(pred_lat[i]), decimals), round(float(pred_lon[i]), decimals))
        if key in obs_map:
            x_obs.append(obs_map[key])
            y_pred.append(float(pred_vt[i]))

    if not x_obs:
        raise ValueError('No coordinate matches found between observed and predicted files.')

    x_obs = np.array(x_obs)
    y_pred = np.array(y_pred)
    residual = y_pred - x_obs
    rmse = float(np.sqrt(np.mean(residual ** 2)))

    neutral = args.neutral
    sign_match = (
        np.sign(y_pred) == np.sign(x_obs)
    ) | (
        (np.abs(y_pred) <= neutral) & (np.abs(x_obs) <= neutral)
    )
    sign_pct = 100.0 * float(np.sum(sign_match)) / float(len(sign_match))

    lo = min(np.min(x_obs), np.min(y_pred))
    hi = max(np.max(x_obs), np.max(y_pred))
    pad = 0.1 * (hi - lo + 1e-6)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

    axes[0].scatter(x_obs, y_pred, s=18, alpha=0.8)
    axes[0].plot([lo - pad, hi + pad], [lo - pad, hi + pad], 'k--', lw=1)
    axes[0].set_xlabel('Observed vT (cm/yr)')
    axes[0].set_ylabel('Predicted vT (cm/yr)')
    axes[0].set_title('Observed vs Predicted')
    axes[0].text(
        0.03,
        0.97,
        'N={}\nRMSE={:.2f} cm/yr\nSign match={:.1f}%'.format(len(x_obs), rmse, sign_pct),
        transform=axes[0].transAxes,
        va='top',
    )

    axes[1].hist(residual, bins=20, edgecolor='black', alpha=0.8)
    axes[1].axvline(0.0, color='k', ls='--', lw=1)
    axes[1].set_xlabel('Predicted - Observed (cm/yr)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Residual Distribution')

    fig.suptitle(args.title)
    ensure_parent_dir(args.output)
    fig.savefig(args.output, dpi=200)


if __name__ == '__main__':
    main()

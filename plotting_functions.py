#!/usr/bin/env python3

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
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


def save_quick_plot(predicted_path, observed_path, output_path, title, tol=0.1, neutral=0.3):
    pred = load_2d(predicted_path)
    obs = load_2d(observed_path)

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

    decimals = max(0, int(round(-np.log10(tol)))) if tol < 1 else 0
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

    fig.suptitle(title)
    ensure_parent_dir(output_path)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_misfit_heatmap(
    sign,
    rms,
    n,
    formulation,
    asth_visc_min,
    asth_visc_max,
    lith_visc_min,
    lith_visc_max,
    yield_min,
    yield_max,
    pre_min,
    pre_max,
    signs_asthen_visc,
    signs_y_param,
    rms_asthen_visc,
    rms_y_param,
    num_sign_matches_max,
    rms_min,
    output_path,
):
    fig = plt.figure(figsize=(6.8, 9.2), constrained_layout=True)

    ax = fig.add_subplot(211)
    if formulation == 2:
        im1 = ax.imshow(
            100.0 * (sign / n),
            cmap=cm.RdYlGn,
            extent=[asth_visc_min, asth_visc_max, yield_min / 1e6, yield_max / 1e6],
            vmax=100.0,
            vmin=60.0,
            aspect=(asth_visc_max - asth_visc_min) / ((yield_max / 1e6) - (yield_min / 1e6)),
        )
        ax.set_ylabel("$\mathregular{\sigma_{Y}}$  [MPa]")
    elif formulation in (1, 7, 8, 9):
        im1 = ax.imshow(
            100.0 * (sign / n),
            cmap=cm.RdYlGn,
            origin='lower',
            extent=[asth_visc_min, asth_visc_max, lith_visc_min, lith_visc_max],
            vmax=100.0,
            vmin=60.0,
            aspect=(asth_visc_max - asth_visc_min) / (lith_visc_max - lith_visc_min),
        )
        ax.set_ylabel("$\mathregular{\eta_{L}}$  [Pa.s]")
    elif formulation == 3:
        im1 = ax.imshow(
            100.0 * (sign / n),
            cmap=cm.RdYlGn,
            origin='lower',
            extent=[asth_visc_min, asth_visc_max, pre_min, pre_max],
            vmax=100.0,
            vmin=60.0,
            aspect=(asth_visc_max - asth_visc_min) / (pre_max - pre_min),
        )
        ax.set_ylabel("K")
    else:
        raise ValueError("Unsupported formulation for plotting: {}".format(formulation))

    ax.scatter(
        [signs_asthen_visc],
        [signs_y_param],
        marker='*',
        s=180,
        c='#2f4b7c',
        edgecolors='black',
        linewidths=0.8,
        zorder=5,
    )
    ax.text(0.03, 0.96, "%s with correct sign" % ('%'), size=19, color="black", transform=ax.transAxes, va='top')
    ax.grid(color='white', linewidth=0.4, alpha=0.25)
    cbar1 = plt.colorbar(im1, ax=ax)
    cbar1.set_label('Correct sign (%)')
    annot_string = ''.join(['best: ', str(num_sign_matches_max), '/', str(n), ' signs, RMS = ', str("%.2f" % rms_min), ' cm/yr'])
    plt.annotate(
        annot_string,
        xy=(0.035, 1.05),
        xycoords='axes fraction',
        verticalalignment='center',
        horizontalalignment='left',
        fontsize=13,
    )

    ax = fig.add_subplot(212)
    if formulation == 2:
        im2 = ax.imshow(
            rms,
            cmap=cm.RdYlGn,
            extent=[asth_visc_min, asth_visc_max, yield_min / 1e6, yield_max / 1e6],
            vmax=10.0,
            vmin=0.0,
            aspect=(asth_visc_max - asth_visc_min) / ((yield_max / 1e6) - (yield_min / 1e6)),
        )
        ax.set_ylabel("$\mathregular{\sigma_{Y}}$  [MPa]")
    elif formulation in (1, 7, 8, 9):
        im2 = ax.imshow(
            rms,
            cmap=cm.RdYlGn,
            origin='lower',
            extent=[asth_visc_min, asth_visc_max, lith_visc_min, lith_visc_max],
            vmax=10.0,
            vmin=0.0,
            aspect=(asth_visc_max - asth_visc_min) / (lith_visc_max - lith_visc_min),
        )
        ax.set_ylabel("$\mathregular{\eta_{L}}$  [Pa.s]")
    elif formulation == 3:
        im2 = ax.imshow(
            rms,
            cmap=cm.RdYlGn,
            origin='lower',
            extent=[asth_visc_min, asth_visc_max, pre_min, pre_max],
            vmax=10.0,
            vmin=0.0,
            aspect=(asth_visc_max - asth_visc_min) / (pre_max - pre_min),
        )
        ax.set_ylabel("K")
    else:
        raise ValueError("Unsupported formulation for plotting: {}".format(formulation))

    ax.set_xlabel("$\mathregular{\eta_{A}}$  [Pa.s]")
    ax.scatter(
        [rms_asthen_visc],
        [rms_y_param],
        marker='*',
        s=180,
        c='#2f4b7c',
        edgecolors='black',
        linewidths=0.8,
        zorder=5,
    )
    ax.text(0.03, 0.96, 'RMS', size=19, color="black", transform=ax.transAxes, va='top')
    ax.grid(color='white', linewidth=0.4, alpha=0.25)
    cbar2 = plt.colorbar(im2, ax=ax)
    cbar2.set_label('RMSE (cm/yr)')

    ensure_parent_dir(output_path)
    fig.savefig(output_path, bbox_inches='tight', dpi=400)
    plt.close(fig)

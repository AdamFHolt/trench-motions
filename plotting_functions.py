#!/usr/bin/env python3

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import numpy as np
from scipy.io import netcdf_file


VECTOR_LENGTH_SCALE_MM_YR = 19.2
PROJ_LON0 = 210.0
LON_MIN = 60.0
LON_MAX = 360.0
LAT_MIN = -70.0
LAT_MAX = 70.0


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
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def resolve_dataset_files(datasets_dir=''):
    if datasets_dir:
        age_grd = os.path.join(datasets_dir, 'age', 'age.3.6.NaN.grd')
        pb_file = os.path.join(datasets_dir, 'plate_boundaries', 'bird_PB2002', 'PB2002_tdiddy.gmt')
    else:
        age_grd = os.path.join('data', 'gmt_datasets', 'age', 'age.3.6.NaN.grd')
        pb_file = os.path.join('data', 'gmt_datasets', 'plate_boundaries', 'bird_PB2002', 'PB2002_tdiddy.gmt')

    # Support flat layout in data/
    if not os.path.isfile(age_grd):
        flat_age = os.path.join('data', 'age.3.6.NaN.grd')
        if os.path.isfile(flat_age):
            age_grd = flat_age
    if not os.path.isfile(pb_file):
        flat_pb = os.path.join('data', 'PB2002_tdiddy.gmt')
        if os.path.isfile(flat_pb):
            pb_file = flat_pb

    return age_grd if os.path.isfile(age_grd) else None, pb_file if os.path.isfile(pb_file) else None


def load_age_grid(age_grd):
    with netcdf_file(age_grd, 'r', mmap=False) as f:
        lon = np.array(f.variables['lon'][:], dtype=float)
        lat = np.array(f.variables['lat'][:], dtype=float)
        z = np.array(f.variables['z'][:], dtype=float)
    return lon, lat, z


def parse_multisegment_gmt(path):
    seg_lons = []
    seg_lats = []
    cur_lon = []
    cur_lat = []
    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if cur_lon:
                    seg_lons.append(np.array(cur_lon, dtype=float))
                    seg_lats.append(np.array(cur_lat, dtype=float))
                    cur_lon = []
                    cur_lat = []
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            lon = float(parts[0])
            lat = float(parts[1])
            if lon < 0:
                lon = lon + 360.0
            cur_lon.append(lon)
            cur_lat.append(lat)
    if cur_lon:
        seg_lons.append(np.array(cur_lon, dtype=float))
        seg_lats.append(np.array(cur_lat, dtype=float))

    # Split segments where longitudes jump across the 0/360 seam to avoid
    # drawing long wrap-around artifacts across the map.
    out_lons = []
    out_lats = []
    for lon_seg, lat_seg in zip(seg_lons, seg_lats):
        if len(lon_seg) < 2:
            out_lons.append(lon_seg)
            out_lats.append(lat_seg)
            continue
        split_inds = np.where(np.abs(np.diff(lon_seg)) > 180.0)[0]
        start = 0
        for idx in split_inds:
            end = idx + 1
            if end - start >= 2:
                out_lons.append(lon_seg[start:end])
                out_lats.append(lat_seg[start:end])
            start = end
        if len(lon_seg) - start >= 2:
            out_lons.append(lon_seg[start:])
            out_lats.append(lat_seg[start:])
    return out_lons, out_lats


def arrow_components(speed_mm_yr, azim_deg, scale=VECTOR_LENGTH_SCALE_MM_YR):
    theta = np.deg2rad(azim_deg)
    mag = -speed_mm_yr / scale
    dx = mag * np.sin(theta)
    dy = mag * np.cos(theta)
    return dx, dy


def project_robinson(lon_deg, lat_deg, lon0_deg=PROJ_LON0):
    lon = np.asarray(lon_deg, dtype=float).copy()
    lat = np.asarray(lat_deg, dtype=float).copy()
    lon[lon > 360.0] -= 360.0
    lon[lon < 0.0] += 360.0

    dlon = ((lon - lon0_deg + 180.0) % 360.0) - 180.0
    lam = np.deg2rad(dlon)
    abslat = np.abs(lat)

    # Robinson table at 5-degree intervals (Snyder).
    x_tab = np.array([
        1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600, 0.9427, 0.9216,
        0.8962, 0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322
    ])
    y_tab = np.array([
        0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958,
        0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000
    ])
    lat_clipped = np.clip(abslat, 0.0, 90.0)
    idx = np.floor(lat_clipped / 5.0).astype(int)
    idx = np.clip(idx, 0, len(x_tab) - 2)
    frac = (lat_clipped - 5.0 * idx) / 5.0
    xr = x_tab[idx] + frac * (x_tab[idx + 1] - x_tab[idx])
    yr = y_tab[idx] + frac * (y_tab[idx + 1] - y_tab[idx])

    x = 0.8487 * lam * xr
    y = 1.3523 * np.sign(lat) * yr
    return x, y


def projected_vector_components(lon, lat, speed, azim):
    dlon, dlat = arrow_components(speed, azim)
    x0, y0 = project_robinson(lon, lat)
    x1, y1 = project_robinson(lon + dlon, lat + dlat)
    return x0, y0, x1 - x0, y1 - y0


def make_map_boundary_patch(ax):
    lon_top = np.linspace(LON_MIN, LON_MAX, 480)
    lat_top = np.full_like(lon_top, LAT_MAX)
    lon_right = np.full(240, LON_MAX)
    lat_right = np.linspace(LAT_MAX, LAT_MIN, 240)
    lon_bottom = np.linspace(LON_MAX, LON_MIN, 480)
    lat_bottom = np.full_like(lon_bottom, LAT_MIN)
    lon_left = np.full(240, LON_MIN)
    lat_left = np.linspace(LAT_MIN, LAT_MAX, 240)

    lon = np.concatenate([lon_top, lon_right, lon_bottom, lon_left])
    lat = np.concatenate([lat_top, lat_right, lat_bottom, lat_left])
    x, y = project_robinson(lon, lat)

    verts = np.column_stack([x, y])
    codes = np.full(len(verts), Path.LINETO, dtype=np.uint8)
    codes[0] = Path.MOVETO
    path = Path(verts, codes, closed=True)
    patch = PathPatch(path, facecolor='none', edgecolor='black', lw=0.8, zorder=10)
    ax.add_patch(patch)
    return patch, x, y


def plot_trench_panel(ax, lon, lat, speed, azim, title, age_grid, pb_segments, speed_norm, add_key=False):
    boundary_patch, bx, by = make_map_boundary_patch(ax)

    if age_grid is not None:
        grid_lon, grid_lat, grid_z = age_grid
        lon_mask = grid_lon >= LON_MIN
        glon = grid_lon[lon_mask]
        gz = grid_z[:, lon_mask]
        # Subsample for speed and memory while preserving broad structure.
        slat = grid_lat[::6]
        slon = glon[::6]
        sz = gz[::6, ::6]
        lons2, lats2 = np.meshgrid(slon, slat)
        gx, gy = project_robinson(lons2, lats2)
        mesh = ax.pcolormesh(
            gx,
            gy,
            sz,
            cmap='Greys',
            vmin=0.0,
            vmax=200.0,
            shading='auto',
            alpha=0.42,
            zorder=0,
        )
        mesh.set_clip_path(boundary_patch)

    if pb_segments is not None:
        seg_lons, seg_lats = pb_segments
        for x, y in zip(seg_lons, seg_lats):
            mask = x >= LON_MIN
            if np.sum(mask) < 2:
                continue
            px, py = project_robinson(x[mask], y[mask])
            line, = ax.plot(px, py, color='darkmagenta', lw=0.5, zorder=1)
            line.set_clip_path(boundary_patch)

    # Boost vector visibility on global maps by using a shorter quiver scale
    # (longer arrows) and a white underlay for contrast.
    px, py, u, v = projected_vector_components(lon, lat, speed, azim)
    q_under = ax.quiver(
        px,
        py,
        u,
        v,
        color='white',
        angles='xy',
        scale_units='xy',
        scale=0.30,
        width=0.0052,
        headwidth=4.2,
        headlength=4.8,
        headaxislength=4.2,
        alpha=0.95,
        zorder=2,
    )
    q_under.set_clip_path(boundary_patch)
    q = ax.quiver(
        px,
        py,
        u,
        v,
        speed,
        cmap='RdYlGn',
        norm=speed_norm,
        angles='xy',
        scale_units='xy',
        scale=0.30,
        width=0.0030,
        headwidth=3.6,
        headlength=4.2,
        headaxislength=3.8,
        zorder=3,
    )
    q.set_clip_path(boundary_patch)
    if add_key:
        _, _, key_u, key_v = projected_vector_components(
            np.array([120.0]), np.array([0.0]), np.array([50.0]), np.array([90.0])
        )
        qk = ax.quiverkey(
            q,
            X=0.02,
            Y=1.03,
            U=float(np.sqrt(key_u[0] ** 2 + key_v[0] ** 2)),
            label='50 mm/yr',
            labelpos='E',
            coordinates='axes',
            color='black',
            fontproperties={'size': 10},
        )
        qk.text.set_clip_on(False)
        qk.text.set_bbox({'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.95, 'pad': 0.25})
    ax.set_xlim(float(np.min(bx)), float(np.max(bx)))
    ax.set_ylim(float(np.min(by)), float(np.max(by)))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    ax.grid(alpha=0.12, lw=0.35)
    ax.set_title(title)
    ax.set_ylabel('')
    return q, boundary_patch


def save_trench_motion_map(predicted_base, observed_file, matches_file, mode, datasets_dir=''):
    predicted_txt = predicted_base + '.txt'
    output_png = predicted_base + '.png'
    if not os.path.isfile(predicted_txt):
        raise ValueError('Missing predicted file: {}'.format(predicted_txt))
    if not os.path.isfile(observed_file):
        raise ValueError('Missing observed file: {}'.format(observed_file))
    if not os.path.isfile(matches_file):
        raise ValueError('Missing matches file: {}'.format(matches_file))

    pred = np.genfromtxt(predicted_txt)
    obs = np.genfromtxt(observed_file)
    matches = np.genfromtxt(matches_file)

    if pred.ndim == 1:
        pred = pred.reshape(1, -1)
    if obs.ndim == 1:
        obs = obs.reshape(1, -1)
    if matches.ndim == 1:
        matches = matches.reshape(1, -1)

    age_grd, pb_file = resolve_dataset_files(datasets_dir)
    age_grid = load_age_grid(age_grd) if age_grd else None
    pb_segments = parse_multisegment_gmt(pb_file) if pb_file else None

    # Observed: lon, lat, azimuth, vt(mm/yr in col 6)
    obs_lon = obs[:, 0].copy()
    obs_lon[obs_lon < 0] += 360.0
    obs_lat = obs[:, 1]
    obs_azim = obs[:, 2]
    obs_vt = obs[:, 5]

    # Predicted: lat, lon, vt(cm/yr), azimuth
    pred_lat = pred[:, 0]
    pred_lon = pred[:, 1].copy()
    pred_lon[pred_lon < 0] += 360.0
    pred_vt = pred[:, 2] * 10.0  # cm/yr -> mm/yr
    pred_azim = pred[:, 3]

    speed_norm = matplotlib.colors.Normalize(vmin=-60.0, vmax=60.0)
    fig, axes = plt.subplots(2, 1, figsize=(12, 8.5), constrained_layout=True)

    q1, _ = plot_trench_panel(axes[0], obs_lon, obs_lat, obs_vt, obs_azim, 'Observed', age_grid, pb_segments, speed_norm, add_key=True)
    _, b2 = plot_trench_panel(axes[1], pred_lon, pred_lat, pred_vt, pred_azim, 'Predicted', age_grid, pb_segments, speed_norm, add_key=False)
    axes[0].set_xlabel('')
    axes[1].set_xlabel('')

    cbar_v = fig.colorbar(q1, ax=axes[0], fraction=0.03, pad=0.03)
    cbar_v.set_label('vT (mm/yr)')

    # Project match points for lower-panel overlays.
    m_lon = matches[:, 1].copy()
    m_lon[m_lon < 0] += 360.0
    m_lat = matches[:, 0]
    mx, my = project_robinson(m_lon, m_lat)

    if mode == 'rms':
        # Matches for RMS are stored in cm/yr; convert to mm/yr for display consistency.
        mvals = 10.0 * matches[:, 2]
        m = axes[1].scatter(mx, my, c=mvals, cmap='seismic_r', vmin=0.0, vmax=60.0, s=12, lw=0.0, zorder=4)
        m.set_clip_path(b2)
        cbar_m = fig.colorbar(m, ax=axes[1], fraction=0.03, pad=0.03)
        cbar_m.set_label('misfit (mm/yr)')
    elif mode == 'signs':
        mvals = matches[:, 2]
        cmap_sign = ListedColormap(['#d73027', '#2b59c3'])
        norm_sign = BoundaryNorm([-0.5, 0.5, 1.5], cmap_sign.N)
        m = axes[1].scatter(mx, my, c=mvals, cmap=cmap_sign, norm=norm_sign, s=12, lw=0.0, zorder=4)
        m.set_clip_path(b2)
        cbar_m = fig.colorbar(m, ax=axes[1], fraction=0.03, pad=0.03)
        cbar_m.set_ticks([0, 1])
        cbar_m.set_ticklabels(['no', 'yes'])
        cbar_m.set_label('signs match?')
    else:
        raise ValueError("Unsupported mode '{}'; expected 'signs' or 'rms'.".format(mode))

    ensure_parent_dir(output_png)
    fig.savefig(output_png, dpi=150, bbox_inches='tight', pad_inches=0.02)
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
    elif formulation in (1, 3, 4):
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
    elif formulation in (1, 3, 4):
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
    fig.savefig(output_path, bbox_inches='tight', dpi=150)
    plt.close(fig)

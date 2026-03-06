import numpy as np

from functions import compute_plate_buoyancy, compute_plate_isotherm


def get_vt_col(vt_ref):
    if vt_ref == 'hs3':
        return 2
    if vt_ref == 'nnr':
        return 3
    return 4


def build_vt_table_from_tnew(data, vt_ref, tol=0.1):
    tnew_path = 'data/vt/tnew.{}.dat'.format(vt_ref)
    tnew = np.genfromtxt(tnew_path)
    if tnew.ndim == 1:
        tnew = tnew.reshape(1, -1)
    if tnew.shape[1] < 6:
        raise ValueError('Observed file has unexpected format: {}'.format(tnew_path))

    decimals = max(0, int(round(-np.log10(tol)))) if tol < 1 else 0
    vt_by_coord = {}
    for i in range(0, tnew.shape[0]):
        lat = round(float(tnew[i, 1]), decimals)
        lon = round(float(tnew[i, 0]), decimals)
        vt_by_coord[(lat, lon)] = float(tnew[i, 5])  # mm/yr

    data_vt = np.full((data.shape[0], 5), np.nan)
    data_vt[:, 0] = data[:, 1]  # lat
    data_vt[:, 1] = data[:, 2]  # lon
    vt_col = get_vt_col(vt_ref)

    n_match = 0
    for i in range(0, data.shape[0]):
        lat = round(float(data[i, 1]), decimals)
        lon = round(float(data[i, 2]), decimals)
        key = (lat, lon)
        if key in vt_by_coord:
            data_vt[i, vt_col] = vt_by_coord[key]
            n_match = n_match + 1

    if n_match == 0:
        raise ValueError('No coordinate matches found between model table and {}'.format(tnew_path))
    return data_vt


def preprocess_data_table(data, limit_max_depth, const_slab_depth, use_avg_Rmin, interpolate_shallow_dip):
    for i in range(0, data.shape[0]):
        if limit_max_depth == 1 and data[i, 7] >= 660:
            data[i, 7] = 660.0
        if const_slab_depth == 1:
            data[i, 7] = 660.0

    rmin_mean = float("%.2f" % np.nanmean(data[:, 13]))
    nan_inds = np.where(np.isnan(data[:, 13]))
    if use_avg_Rmin == 1:
        data[nan_inds, 13] = rmin_mean

    if interpolate_shallow_dip == 1:
        sum_diff = 0.0
        n_diff = 0
        for i in range(0, data.shape[0]):
            if np.isnan(data[i, 5]) == 0 and np.isnan(data[i, 6]) == 0:
                sum_diff = sum_diff + (data[i, 6] - data[i, 5])
                n_diff = n_diff + 1
        avg_dip_diff = sum_diff / n_diff
        for i in range(0, data.shape[0]):
            if np.isnan(data[i, 5]) == 0 and np.isnan(data[i, 6]) == 1:
                data[i, 6] = data[i, 5] + avg_dip_diff


def build_segment_arrays(data, data_vt, vt_col, vel_converter, max_age, calc_slabL_using_dip):
    num = 0
    for i in range(0, data.shape[0]):
        if np.isnan(data[i, 26]) == False and np.isnan(data[i, 6]) == False and np.isnan(data[i, 8]) == False \
            and np.isnan(data[i, 13]) == False and np.isnan(data[i, 20]) == False and np.isnan(data[i, 7]) == False:
            num = num + 1

    Lsp = np.zeros((num, 1))
    dip = np.zeros((num, 1))
    slabD = np.zeros((num, 1))
    vc = np.zeros((num, 1))
    Rmin = np.zeros((num, 1))
    age = np.zeros((num, 1))
    vt_actual = np.zeros((num, 1))
    slabL = np.zeros((num, 1))
    latlon = np.zeros((num, 2))
    azims = np.zeros((num, 1))
    w = np.zeros((num, 1))
    slabL_buoy = np.zeros((num, 1))
    Lop = np.zeros((num, 1))
    external_force_factor = np.zeros((num, 1))

    n = 0
    for i in range(0, data.shape[0]):
        if np.isnan(data[i, 26]) == False and np.isnan(data[i, 6]) == False and np.isnan(data[i, 8]) == False \
            and np.isnan(data[i, 13]) == False and np.isnan(data[i, 20]) == False and np.isnan(data[i, 7]) == False:
            latlon[n, 0] = data_vt[i, 0]
            latlon[n, 1] = data_vt[i, 1]
            azims[n, 0] = data[i, 4]
            Lsp[n, 0] = data[i, 26] * 1e3
            Lop[n, 0] = data[i, 27] * 1e3
            dip[n, 0] = data[i, 6]
            vc[n, 0] = (data[i, 9] / 10.0) * vel_converter
            Rmin[n, 0] = data[i, 13] * 1e3
            if data[i, 20] > max_age:
                age[n, 0] = max_age
            else:
                age[n, 0] = data[i, 20]
            w[n, 0] = data[i, 25] * 1e3
            vt_actual[n, 0] = data_vt[i, vt_col] / 10.0
            external_force_factor[n, 0] = data[i, 27]
            slabD[n, 0] = data[i, 7] * 1e3

            if calc_slabL_using_dip == 1:
                slabL[n, 0] = slabD[n, 0] / np.tan(np.deg2rad(dip[n, 0]))
                slabL_buoy[n, 0] = slabD[n, 0] / np.tan(np.deg2rad(dip[n, 0]))
            else:
                slabL[n, 0] = data[i, 8] * 1e3
                slabL_buoy[n, 0] = data[i, 8] * 1e3

            n = n + 1

    return {
        'n': n,
        'num': num,
        'Lsp': Lsp,
        'dip': dip,
        'slabD': slabD,
        'vc': vc,
        'Rmin': Rmin,
        'age': age,
        'vt_actual': vt_actual,
        'slabL': slabL,
        'latlon': latlon,
        'azims': azims,
        'w': w,
        'slabL_buoy': slabL_buoy,
        'Lop': Lop,
        'external_force_factor': external_force_factor,
    }


def build_thermal_terms(age, include_ridge_push, dT, g, rho0, rhoW, alpha, kappa, ma_to_s):
    n = age.shape[0]
    H = np.zeros((n, 1))
    oceanic_buoy = np.zeros((n, 1))
    ridge_push = np.zeros((n, 1))
    for i in range(0, n):
        H[i] = compute_plate_isotherm(age[i], 1300, 1e-6, 1200) * 1e3
        oceanic_buoy[i] = compute_plate_buoyancy(age[i], 1300, 1e-6, 3e-5, 3300)
        if include_ridge_push == 1:
            ridge_push[i] = g * rho0 * alpha * dT * (1.0 + ((2.0 * rho0 * alpha * dT) / (np.pi * (rho0 - rhoW)))) * kappa * (age[i] * ma_to_s)
    return H, oceanic_buoy, ridge_push

#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import os, numpy as np
import sys, tempfile, shutil
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml
from functions import compute_vsp_withDP
from plotting_functions import save_quick_plot, save_trench_motion_map
from workflow_common import get_vt_col, build_vt_table_from_tnew, preprocess_data_table, build_segment_arrays, build_thermal_terms
plt.ioff()


def ensure_parent_dir(path):
	parent = os.path.dirname(path)
	if parent and not os.path.isdir(parent):
		os.makedirs(parent)


USAGE = """Usage:
  python3 compute_rates_single.py <vt_ref> <formulation> <include_DP> <DP_ref> <asthen_visc> <lith_visc> [--skip-map] [--out-prefix <dir>]
  python3 compute_rates_single.py --config <path.yaml> [--skip-map] [--out-prefix <dir>]

Arguments:
  vt_ref                hs3 | nnr | sa
  formulation           integer model id
  include_DP            0 | 1
  DP_ref                float (Pa)
  asthen_visc           float (Pa.s)
  lith_visc             float (Pa.s)

Optional flags:
  --skip-map  Skip map plotting call.
  --out-prefix <dir>  Write generated outputs under a custom base directory.
  --config <path.yaml>  Load arguments from YAML config.
"""


raw_args = sys.argv[1:]
if '--help' in raw_args or '-h' in raw_args:
	print(USAGE)
	sys.exit(0)

cfg = {}
args_wo_config = list(raw_args)
if '--config' in args_wo_config:
	config_ind = args_wo_config.index('--config')
	if config_ind + 1 >= len(args_wo_config):
		print("Missing value for --config")
		print(USAGE)
		sys.exit(2)
	config_path = args_wo_config[config_ind + 1]
	del args_wo_config[config_ind:config_ind + 2]
	with open(config_path, 'r') as f:
		cfg = yaml.safe_load(f) or {}
	if not isinstance(cfg, dict):
		print("Config must deserialize to a mapping/dictionary")
		sys.exit(2)

if cfg:
	missing_keys = [
		k for k in [
			'vt_ref',
			'formulation',
			'include_DP',
			'DP_ref',
			'asthen_visc',
			'lith_visc',
		]
		if k not in cfg
	]
	if missing_keys:
		print("Missing config keys: %s" % ", ".join(missing_keys))
		sys.exit(2)
	vt_ref = str(cfg['vt_ref'])
	formulation = int(cfg['formulation'])
	include_DP = int(cfg['include_DP'])
	DP_ref = float(cfg['DP_ref'])
	asthen_visc = float(cfg['asthen_visc'])
	lith_visc = float(cfg['lith_visc'])
	extra_args = args_wo_config
	skip_map = bool(cfg.get('skip_map', False))
	out_prefix = str(cfg.get('out_prefix', ''))
else:
	if len(args_wo_config) < 6:
		print(USAGE)
		sys.exit(2)
	main_args = args_wo_config[:6]
	extra_args = args_wo_config[6:]
	vt_ref=str(main_args[0])    				# hs3, nnr, sa
	formulation=int(main_args[1])			# 1 = regular, 2 = plastic bending, 3 = regular, hSP \propto LSP, 4 = regular, hSP \propto VSP, 5 = regular with OP
	include_DP=int(main_args[2])				# 1 = include DP force, 0 = do not
	DP_ref=float(main_args[3]) 				# DP values from analytical computations: free slip base: avg DP_0 = 18.3, max DP_0 = 23.5, no slip: avg DP_0 = 73.1, max DP_0 = 93.9 MPa
	asthen_visc = float(main_args[4])
	lith_visc = float(main_args[5])
	skip_map = False
	out_prefix = ''

i = 0
while i < len(extra_args):
	arg = extra_args[i]
	if arg == '--skip-map':
		skip_map = True
		i += 1
	elif arg == '--out-prefix':
		if i + 1 >= len(extra_args):
			print("Missing value for --out-prefix")
			print(USAGE)
			sys.exit(2)
		out_prefix = extra_args[i + 1]
		i += 2
	else:
		print("Unknown argument: %s" % arg)
		print(USAGE)
		sys.exit(2)


def out_path(relative_path):
	base_ref_dir = out_prefix if out_prefix else os.path.join('plots', vt_ref)
	return os.path.join(base_ref_dir, formulation_slug(formulation), 'maps', relative_path)


def misfit_plot_out_path(filename):
	base_ref_dir = out_prefix if out_prefix else os.path.join('plots', vt_ref)
	return os.path.join(base_ref_dir, formulation_slug(formulation), 'param-sweep', filename)


def formulation_slug(formulation_id):
	if formulation_id == 1:
		return 'viscous'
	if formulation_id == 2:
		return 'plastic'
	if formulation_id == 3:
		return 'viscous_LspShear'
	if formulation_id == 4:
		return 'viscous_VspShear'
	raise ValueError("unsupported formulation: {}".format(formulation_id))

if vt_ref not in ['hs3', 'nnr', 'sa']:
	print("Invalid vt_ref: %s" % vt_ref)
	print(USAGE)
	sys.exit(2)
if include_DP not in [0, 1]:
	print("include_DP must be 0 or 1")
	print(USAGE)
	sys.exit(2)
if formulation not in [1, 2, 3, 4]:
	print("unsupported formulation: %s" % formulation)
	print(USAGE)
	sys.exit(2)
if include_DP == 0:
	print("trying to set include_DP=0. Can't do that - just set DP_ref = 0...")
	sys.exit(1)


# calculation parameters
calc_slabL_using_dip = 1 		# 0 = take down-dip slab length from table, 1 = calculate it from dip and slabD
const_slab_depth = 0  			# 0 = use lallemand depths, 1 = all slabs go to 660 km
limit_max_depth = 0  			# 0 = no slab depth limit,  1 = limit depth to 660
use_avg_Rmin = 1      			# 1 = for segments without an Rmin, use the global average.
include_ridge_push = 1 			# 0 = no ridge push, 1 = approximation for ridge push
# get_DP_from_dip = 1 			# 0 = no DP from dip (get DP from analytical slab wall), 1 = DP from dip, -1 = DP from dip (in opposite direction)

# less important calculation parameters:
interpolate_shallow_dip = 0 # doesn't improve things
limit_max_age = 0     		# has a negligible effect

# reference parameters
vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s  
h = 200e3; # m
n = 3.5;
g = 9.81;
ma_to_s = 1e6 * 365 * 24 * 60 * 60;
kappa = 1e-6; alpha = 3.e-5
dT = 1300.
rho0 = 3300.; rhoW = 1000.

# search parameters
lith_visc_min = 20; lith_visc_max = 25
asth_visc_min = 17; asth_visc_max = 22
yield_min = 100e6; 	yield_max = 10000e6;


# DP stuff (SI units), used in all active formulations
visc_asthen_ref =  3e20
w_ref = 2224e3 * 2 # [m], factor 2 to convert from half width to full width,  (2224km = 20 degrees)
trenchv_ref = 5.0 * vel_converter # [cm/yr] -> [m/s], positive (as here) = trench retreat


data_name = 'data/Lallemand_et_al-2005_G3_dataset_WithAdditionalParams.txt'
data = np.genfromtxt(data_name) # see spreadsheet for column details
vt_col = get_vt_col(vt_ref)
data_vt = build_vt_table_from_tnew(data, vt_ref)
preprocess_data_table(data, limit_max_depth, const_slab_depth, use_avg_Rmin, interpolate_shallow_dip)

max_age = 1000
if limit_max_age == 1:
	max_age = 90.

segments = build_segment_arrays(data, data_vt, vt_col, vel_converter, max_age, calc_slabL_using_dip)
n = segments['n']
num = segments['num']
Lsp = segments['Lsp']
dip = segments['dip']
slabD = segments['slabD']
vc = segments['vc']
Rmin = segments['Rmin']
age = segments['age']
vt_actual = segments['vt_actual']
slabL = segments['slabL']
latlon = segments['latlon']
azims = segments['azims']
w = segments['w']
slabL_buoy = segments['slabL_buoy']

H, oceanic_buoy, ride_push = build_thermal_terms(age, include_ridge_push, dT, g, rho0, rhoW, alpha, kappa, ma_to_s)

rms_min = 10.0
num_sign_matches_max = 0
stored_rms_vals = np.zeros((num,1));

visc_asthen = asthen_visc
visc_lith = lith_visc
yield_sigma = 1
pre = 1

vsp_estimate = compute_vsp_withDP(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,slabL_buoy,dip,oceanic_buoy,DP_ref,visc_asthen_ref,\
	w_ref,trenchv_ref,w,slabD,yield_sigma,pre,ride_push)

vt_estimate = (vc/vel_converter) - vsp_estimate
rms_sum = 0; num_sign_matches = 0;
for i in range(0,n):
	diff = float(vt_estimate[i,0] - vt_actual[i,0])
	rms_sum = rms_sum + diff**2;
	stored_rms_vals[i] = np.abs(diff)
	if np.sign(vt_estimate[i]) == np.sign(vt_actual[i]):
		num_sign_matches = num_sign_matches + 1
	elif (vt_estimate[i] > -0.3 and vt_estimate[i] <= 0.3) and (vt_actual[i] > -0.3 and vt_actual[i] <= 0.3):
		num_sign_matches = num_sign_matches + 1
rms_act = float(np.sqrt(rms_sum / float(n)))

if rms_act < rms_min:
	rms_min = rms_act;
	rms_lith_visc = np.log10(visc_lith);
	rms_asthen_visc = np.log10(visc_asthen);
	rms_yield_stress = yield_sigma/1e6;
	rms_predicted_vts = np.concatenate((latlon,vt_estimate,azims), axis=1)
	rms_separated = np.concatenate((latlon,stored_rms_vals), axis=1)
	vsps = np.concatenate((latlon,vsp_estimate), axis=1)

if num_sign_matches > num_sign_matches_max:
	num_sign_matches_max = num_sign_matches;

print("------------")
print("Minimum RMS = %.2f cm/yr, Signs %.0f/%.0f" % (rms_min,num_sign_matches_max,n))


if include_DP == 0:
	DP_string=''
else:
	DP_string = ''.join(['.DP',str("%.3g" % (DP_ref/1e6)),'MPa'])

if include_ridge_push:
	RP_string = '.withRP'
else:
	RP_string = ''

if formulation == 1:  # viscous bending
	plot_name=''.join(['misfits_',str(vt_ref),'model',DP_string,RP_string,'.viscous_bending.png'])
	rms_name=''.join(['rms_',str(vt_ref),'model',DP_string,RP_string,'.l',str(rms_lith_visc),'_a10e',str(rms_asthen_visc),'.viscous_bending'])
elif formulation == 2: # plastic bending
	plot_name=''.join(['misfits_',str(vt_ref),'model',DP_string,RP_string,'.plastic_bending.png'])
	rms_name=''.join(['rms_',str(vt_ref),'model',DP_string,RP_string,'.y',str(rms_yield_stress),'_a',str(rms_asthen_visc),'.plastic_bending'])
elif formulation == 3:  # regular, hSP \propto LSP
	plot_name=''.join(['misfits_',str(vt_ref),'model',DP_string,RP_string,'.viscous_bending_hSPproptoLSP.png'])
	rms_name=''.join(['rms_',str(vt_ref),'model',DP_string,RP_string,'.l',str(rms_lith_visc),'_a10e',str(rms_asthen_visc),'.viscous_bending_hSPproptoLSP'])
elif formulation == 4:  # regular, hSP \propto VSP
	plot_name=''.join(['misfits_',str(vt_ref),'model',DP_string,RP_string,'.viscous_bending_hSPproptoVSP.png'])
	rms_name=''.join(['rms_',str(vt_ref),'model',DP_string,RP_string,'.l',str(rms_lith_visc),'_a10e',str(rms_asthen_visc),'.viscous_bending_hSPproptoVSP'])
plot_name = misfit_plot_out_path(plot_name)
rms_name = out_path(rms_name)

# plot map
rms_txt_name=''.join([rms_name,'.txt'])
ensure_parent_dir(rms_txt_name)
np.savetxt(rms_txt_name, rms_predicted_vts, fmt='%.4f')

vt_observed=''.join(['tnew.',str(vt_ref),'.dat'])
if skip_map:
	print("skipping map plotting (--skip-map)")
	try:
		save_quick_plot(
			predicted_path=rms_txt_name,
			observed_path=os.path.join('data', 'vt', vt_observed),
			output_path=plot_name,
			title='Single-run quick check ({})'.format(vt_ref),
		)
	except Exception as exc:
		print("warning: quick plot generation failed ({})".format(exc))
else:
	tmp_dir = tempfile.mkdtemp(prefix='trench-motions-')
	rms_sep_base = os.path.join(tmp_dir, 'rms_separated')
	np.savetxt(''.join([rms_sep_base, '.txt']), vsps, fmt='%.4f')
	try:
		save_trench_motion_map(
			predicted_base=rms_name,
			observed_file=os.path.join('data', 'vt', vt_observed),
			matches_file=''.join([rms_sep_base, '.txt']),
			mode='rms',
			datasets_dir=os.environ.get('DATASETS_DIR', ''),
		)
	finally:
		shutil.rmtree(tmp_dir, ignore_errors=True)
print("output: %s" % plot_name)

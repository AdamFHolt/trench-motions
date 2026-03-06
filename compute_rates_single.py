#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import os, numpy as np
import subprocess, sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml
from functions import compute_plate_buoyancy, compute_plate_isotherm
from functions import compute_vsp_withDP
plt.ioff()


def ensure_parent_dir(path):
	parent = os.path.dirname(path)
	if parent and not os.path.isdir(parent):
		os.makedirs(parent)


USAGE = """Usage:
  python3 compute_rates_single.py <vt_ref> <formulation> <include_DP> <DP_ref> <PSP_slab_pull_factor> <asthen_visc> <lith_visc> [--skip-map] [--out-prefix <dir>]
  python3 compute_rates_single.py --config <path.yaml> [--skip-map] [--out-prefix <dir>]

Arguments:
  vt_ref                hs3 | nnr | sa
  formulation           integer model id
  include_DP            0 | 1
  DP_ref                float (Pa)
  PSP_slab_pull_factor  float
  asthen_visc           float (Pa.s)
  lith_visc             float (Pa.s)

Optional flags:
  --skip-map  Skip GMT map plotting subprocess call.
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
			'PSP_slab_pull_factor',
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
	PSP_slab_pull_factor = float(cfg['PSP_slab_pull_factor'])
	asthen_visc = float(cfg['asthen_visc'])
	lith_visc = float(cfg['lith_visc'])
	extra_args = args_wo_config
	skip_map = bool(cfg.get('skip_map', False))
	out_prefix = str(cfg.get('out_prefix', ''))
else:
	if len(args_wo_config) < 7:
		print(USAGE)
		sys.exit(2)
	main_args = args_wo_config[:7]
	extra_args = args_wo_config[7:]
	vt_ref=str(main_args[0])    				# hs3, nnr, sa
	formulation=int(main_args[1])			# 1 = regular, 2 = plastic bending, 3 = slab pull pre-factor
	include_DP=int(main_args[2])				# 1 = include DP force, 0 = do not
	DP_ref=float(main_args[3]) 				# DP values from analytical computations: free slip base: avg DP_0 = 18.3, max DP_0 = 23.5, no slip: avg DP_0 = 73.1, max DP_0 = 93.9 MPa
	PSP_slab_pull_factor=float(main_args[4]) 	
	asthen_visc = float(main_args[5])
	lith_visc = float(main_args[6])
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
	if out_prefix:
		return os.path.join(out_prefix, relative_path)
	return os.path.join('results', relative_path)

if vt_ref not in ['hs3', 'nnr', 'sa']:
	print("Invalid vt_ref: %s" % vt_ref)
	print(USAGE)
	sys.exit(2)
if include_DP not in [0, 1]:
	print("include_DP must be 0 or 1")
	print(USAGE)
	sys.exit(2)
if formulation in [4, 5, 6]:
	print("formulations 4-6 (power-law) are no longer supported")
	print(USAGE)
	sys.exit(2)
if formulation not in [1, 2, 3]:
	print("unsupported formulation: %s" % formulation)
	print(USAGE)
	sys.exit(2)
composite = 0
trans_strain_rate = 1.0e-13


# calculation parameters
include_PSP_interactions = 1    # include simple parameterization of PSP force transmission?
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
ma_to_s = 1e6 * 356 * 24 * 60 * 60;
kappa = 1e-6; alpha = 3.e-5
dT = 1500.
rho0 = 3300.; rhoW = 1000.

# PSP force transmission
PSP_age = 45.  		# Ma
PSP_slab_pull = compute_plate_buoyancy(PSP_age,dT,kappa,alpha,rho0) * 600.e3 * 9.81
PSP_force_transmitted = PSP_slab_pull_factor * PSP_slab_pull
	
# search parameters
lith_visc_min = 20; lith_visc_max = 25
asth_visc_min = 17; asth_visc_max = 22
yield_min = 100e6; 	yield_max = 10000e6;
pre_min = 0.0; 		pre_max = 1.0;


# DP stuff (SI units), for formulation == 8
visc_asthen_ref =  3e20
w_ref = 2224e3 * 2 # [m], factor 2 to convert from half width to full width,  (2224km = 20 degrees)
trenchv_ref = 5.0 * vel_converter # [cm/yr] -> [m/s], positive (as here) = trench retreat


data_name = 'data/Lallemand_et_al-2005_G3_dataset_WithTrenchWidths_withInteractions.txt'
data = np.genfromtxt(data_name) # see spreadsheet for column details
data_vt =np.genfromtxt('data/vts_hs3-nnr-sa.txt') # lat, lon, vtn: hs3, vtn: nnr, vtn: sa
if vt_ref == 'hs3':
   vt_col = 2;
elif vt_ref == 'nnr':
   vt_col = 3;
else:
   vt_col = 4; 

for i in range(0,data.shape[0]):
	if limit_max_depth == 1 and data[i,7] >= 660:
		data[i,7] = 660.
	if const_slab_depth == 1:
		data[i,7] = 660.

# get average radius of curvature
Rmin_mean = "%.2f" % np.nanmean(data[:,13])
nan_inds = np.where(np.isnan(data[:,13]))
if use_avg_Rmin == 1:
	data[nan_inds,13] = Rmin_mean


if interpolate_shallow_dip == 1:
	# get average difference between deep and shallow dip
	sum_diff = 0; n_diff = 0;
	for i in range(0,data.shape[0]):
		if np.isnan(data[i,5]) == 0 and np.isnan(data[i,6]) == 0:
			sum_diff = sum_diff + (data[i,6] - data[i,5])
			n_diff	 = n_diff + 1
	avg_dip_diff = sum_diff / n_diff
	# replace in places where deep dip is missing
	for i in range(0,data.shape[0]):
		if np.isnan(data[i,5]) == 0 and np.isnan(data[i,6]) == 1:
			data[i,6] = data[i,5] + avg_dip_diff

# plate model?
max_age = 1000
if limit_max_age == 1:
	max_age = 90.

num = 0;
for i in range(0,data.shape[0]):
	if np.isnan(data[i,26]) == False and np.isnan(data[i,6]) == False and np.isnan(data[i,8]) == False \
		and np.isnan(data[i,13]) == False and np.isnan(data[i,20]) == False and np.isnan(data[i,7]) == False:
		num = num + 1

Lsp=np.zeros((num,1));			dip=np.zeros((num,1))
slabD=np.zeros((num,1));		vc=np.zeros((num,1))
Rmin=np.zeros((num,1));			age=np.zeros((num,1))
vt_actual=np.zeros((num,1));	slabL=np.zeros((num,1))
latlon=np.zeros((num,2));		azims=np.zeros((num,1))
w = np.zeros((num,1));			slabL_buoy=np.zeros((num,1))
Lop=np.zeros((num,1));
external_force_factor = np.zeros((num,1));	ride_push = np.zeros((num,1))

n = 0
for i in range(0,data.shape[0]):
	if np.isnan(data[i,26]) == False and np.isnan(data[i,6]) == False and np.isnan(data[i,8]) == False \
		and np.isnan(data[i,13]) == False and np.isnan(data[i,20]) == False and np.isnan(data[i,7]) == False:
		
		latlon[n,0]=data_vt[i,0]; latlon[n,1]=data_vt[i,1];  
		azims[n,0] = data[i,4]
		Lsp[n,0] = data[i,26] * 1e3
		Lop[n,0] = data[i,27] * 1e3
		dip[n,0] = data[i,6]
		vc[n,0] = (data[i,9] / 10.0) * vel_converter
		Rmin[n,0] = data[i,13] * 1e3
		if data[i,20] > max_age:
			age[n,0] = max_age
		else:
			age[n,0] = data[i,20]
		w[n,0] = data[i,25] * 1e3
		vt_actual[n,0] = data_vt[i,vt_col] / 10.0
		external_force_factor[n,0] = data[i,27]
		slabD[n,0] = data[i,7] * 1e3

		if calc_slabL_using_dip == 1:
			slabL[n,0] = slabD[n,0]/np.tan(np.deg2rad(dip[n,0]))
			slabL_buoy[n,0] = slabD[n,0]/np.tan(np.deg2rad(dip[n,0]))
		else:
			slabL[n,0] = data[i,8] * 1e3
			slabL_buoy[n,0] = data[i,8] * 1e3

		n = n + 1

if include_PSP_interactions == 0:
	external_force_factor = 0. * external_force_factor

H = np.zeros((num,1))
oceanic_buoy = np.zeros((num,1))
for i in range(0,n):
	H[i] = compute_plate_isotherm(age[i],1300,1e-6,1200) * 1e3; # [m]
	oceanic_buoy[i] = compute_plate_buoyancy(age[i],1300,1e-6,3e-5,3300); # [kg/m2]
	if include_ridge_push == 1:
		ride_push[i] = g * rho0 * alpha * dT * ( 1.0 + ((2.0 * rho0 * alpha * dT)/(np.pi * (rho0 - rhoW))) ) * kappa * (age[i] * ma_to_s); # [N/m]

rms_min = 10.0
num_sign_matches_max = 0
for k in range(0,1):

	for j in range(0,1):

		stored_rms_vals = np.zeros((num,1));
		stored_sign_vals = np.zeros((num,1));

		visc_asthen = asthen_visc
		visc_lith = lith_visc
		yield_sigma = 1
		pre = 1

		if include_DP == 1:
			vsp_estimate = compute_vsp_withDP(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,slabL_buoy,dip,oceanic_buoy,DP_ref,visc_asthen_ref,\
				w_ref,trenchv_ref,w,slabD,yield_sigma,n,pre,trans_strain_rate,composite,external_force_factor,PSP_force_transmitted,ride_push,Lop)
		else:
			sys.exit(1)
			# vsp_estimate = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,dip,oceanic_buoy,yield_sigma,trans_strain_rate,n,pre)

		vt_estimate = (vc/vel_converter) - vsp_estimate
		rms_sum = 0; num_sign_matches = 0;
		for i in range(0,n):
			diff = float(vt_estimate[i,0] - vt_actual[i,0])
			rms_sum = rms_sum + diff**2;
			stored_rms_vals[i] = np.abs(diff)
			if np.sign(vt_estimate[i]) == np.sign(vt_actual[i]):
				num_sign_matches = num_sign_matches + 1
				stored_sign_vals[i] = 1
			elif (vt_estimate[i] > -0.5 and vt_estimate[i] <= 0.5) and (vt_actual[i] > -0.5 and vt_actual[i] <= 0.5):
				num_sign_matches = num_sign_matches + 1
				stored_sign_vals[i] = 1
		rms_act = float(np.sqrt(rms_sum / float(n)))


		if rms_act < rms_min:
			rms_min = rms_act;
			rms_lith_visc = np.log10(visc_lith);
			rms_asthen_visc = np.log10(visc_asthen);
			rms_yield_stress = yield_sigma/1e6;
			rms_pre = pre;
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
	DP_string = ''.join(['.DP',str("%.3g" % DP_ref),'MPa'])

if include_PSP_interactions == 1:
	PSP_string = ''.join(['.PSPfact',str(PSP_slab_pull_factor)])
else:
	PSP_string = ''

if include_ridge_push:
	RP_string = '.withRP'
else:
	RP_string = ''

if formulation == 1:  # viscous bending
	plot_name=''.join(['plots/misfits_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.viscous_bending.png'])
	rms_name=''.join(['predictions/rms_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.l',str(rms_lith_visc),'_a10e',str(rms_asthen_visc),'.viscous_bending'])
elif formulation == 2: # plastic bending
	plot_name=''.join(['plots/misfits_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.plastic_bending.png'])
	rms_name=''.join(['predictions/rms_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.y',str(rms_yield_stress),'_a',str(rms_asthen_visc),'.plastic_bending'])
elif formulation == 3: # just slab pull (with prefactor)
	plot_name=''.join(['plots/misfits_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.just-slab-pull.png'])
	rms_name=''.join(['predictions/rms_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.pre',str(rms_pre),'_a',str(rms_asthen_visc),'.just-slab-pull'])
plot_name = out_path(plot_name)
rms_name = out_path(rms_name)
ensure_parent_dir(plot_name)
plt.savefig(plot_name, bbox_inches='tight',dpi=400)

# plot map
rms_txt_name=''.join([rms_name,'.txt'])
ensure_parent_dir(rms_txt_name)
tmp_dir = out_path('tmp')
if not os.path.isdir(tmp_dir):
	os.makedirs(tmp_dir)
rms_sep_base = os.path.join(tmp_dir, 'rms_separated')
np.savetxt(rms_txt_name, rms_predicted_vts, fmt='%.4f')  
np.savetxt(''.join([rms_sep_base, '.txt']), vsps, fmt='%.4f')  

vt_observed=''.join(['tnew.',str(vt_ref),'.dat'])  
if skip_map:
	print("skipping map plotting (--skip-map)")
	quick_plot_name = ''.join([rms_name, '.quick.png'])
	try:
		subprocess.check_call([
			sys.executable,
			'quick_plot.py',
			'--predicted', rms_txt_name,
			'--observed', os.path.join('data', 'vt', vt_observed),
			'--output', quick_plot_name,
			'--title', 'Single-run quick check ({})'.format(vt_ref)
		])
	except subprocess.CalledProcessError as exc:
		print("warning: quick plot generation failed (exit code {})".format(exc.returncode))
else:
	subprocess.check_call(['./plot_trench_motions.sh',rms_name,vt_observed,rms_sep_base,'1'])
print("output: %s" % plot_name)

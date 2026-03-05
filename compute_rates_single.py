#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import scipy, math, os, numpy as np
import subprocess, sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from functions import compute_plate_buoyancy, compute_plate_isotherm
from functions import compute_vsp_withDP
plt.ioff()

# vt data to compare
vt_ref=str(sys.argv[1])    				# hs3, nnr, sa
formulation=int(sys.argv[2])			# 1 = regular, 2 = plastic bending, 3 = slab pull pre-factor, 4 = regular, power-law, 5 = plastic, power-law, 6 = pre-factor, power-law
include_DP=int(sys.argv[3])				# 1 = include DP force, 0 = do not
DP_ref=float(sys.argv[4]) 				# DP values from analytical computations: free slip base: avg DP_0 = 18.3, max DP_0 = 23.5, no slip: avg DP_0 = 73.1, max DP_0 = 93.9 MPa
trans_strain_rate=float(sys.argv[5]) 	# typical value: 1e-13 s-1
PSP_slab_pull_factor=float(sys.argv[6]) 	
asthen_visc = float(sys.argv[7])
lith_visc = float(sys.argv[8])
composite = 0


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
external_force_factor = np.zeros((num,1));	ride_push = np.zeros((num,1))

n = 0
for i in range(0,data.shape[0]):
	if np.isnan(data[i,26]) == False and np.isnan(data[i,6]) == False and np.isnan(data[i,8]) == False \
		and np.isnan(data[i,13]) == False and np.isnan(data[i,20]) == False and np.isnan(data[i,7]) == False:
		
		latlon[n,0]=data_vt[i,0]; latlon[n,1]=data_vt[i,1];  
		azims[n,0] = data[i,4]
		Lsp[n,0] = data[i,26] * 1e3
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
		print ride_push[i]

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

		# power-law
		if formulation >= 4:
			vsp_estimate=np.zeros((num,1));
			for i in range(0,len(vt_actual)):
				if include_DP == 1:
					vsp_estimate[i] = compute_vsp_withDP(formulation,vc[i],h,visc_asthen,visc_lith,H[i],Lsp[i],Rmin[i],slabL[i],slabL_buoy[i],dip[i],oceanic_buoy[i],DP_ref,visc_asthen_ref,\
						w_ref,trenchv_ref,w[i],slabD[i],yield_sigma,n,pre,trans_strain_rate,composite,external_force_factor[i],PSP_force_transmitted,ride_push[i])
				else:
					exit()
					# vsp_estimate[i] = compute_vsp(formulation,vc[i],h,visc_asthen,visc_lith,H[i],Lsp[i],Rmin[i],slabL[i],dip[i],oceanic_buoy[i],yield_sigma,trans_strain_rate,n,pre)

		# non power-law			
		else:
			if include_DP == 1:
				vsp_estimate = compute_vsp_withDP(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,slabL_buoy,dip,oceanic_buoy,DP_ref,visc_asthen_ref,\
					w_ref,trenchv_ref,w,slabD,yield_sigma,n,pre,trans_strain_rate,composite,external_force_factor,PSP_force_transmitted,ride_push)
			else:
				exit()
				# vsp_estimate = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,dip,oceanic_buoy,yield_sigma,trans_strain_rate,n,pre)

		vt_estimate = (vc/vel_converter) - vsp_estimate
		rms_sum = 0; num_sign_matches = 0;
		for i in range(0,n):
			rms_sum = rms_sum + (vt_estimate[i] - vt_actual[i])**2;
			stored_rms_vals[i] = np.abs(vt_estimate[i] - vt_actual[i])
			if np.sign(vt_estimate[i]) == np.sign(vt_actual[i]):
				num_sign_matches = num_sign_matches + 1
				stored_sign_vals[i] = 1
			elif (vt_estimate[i] > -0.5 and vt_estimate[i] <= 0.5) and (vt_actual[i] > -0.5 and vt_actual[i] <= 0.5):
				num_sign_matches = num_sign_matches + 1
				stored_sign_vals[i] = 1
		rms_act = np.sqrt(rms_sum / float(n))


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
			signs_lith_visc = np.log10(visc_lith);
			signs_asthen_visc = np.log10(visc_asthen);
			signs_yield_stress = yield_sigma/1e6;
			signs_pre = pre;
			signs_predicted_vts = np.concatenate((latlon,vt_estimate,azims), axis=1)
			signs_separated = np.concatenate((latlon,stored_sign_vals), axis=1)

print "------------"
print "Minimum RMS = %.2f cm/yr, Signs %.0f/%.0f" % (rms_min,num_sign_matches_max,n)


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
	plot_name=''.join(['plots/new/linear/misfits_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.viscous_bending.png'])
	signs_name=''.join(['predictions/new/linear/signs_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.l',str(signs_lith_visc),'_a10e',str(signs_asthen_visc),'.viscous_bending'])
	rms_name=''.join(['predictions/new/linear/rms_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.l',str(rms_lith_visc),'_a10e',str(rms_asthen_visc),'.viscous_bending'])
elif formulation == 2: # plastic bending
	plot_name=''.join(['plots/new/linear/misfits_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.plastic_bending.png'])
	signs_name=''.join(['predictions/new/linear/signs_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.y',str(signs_yield_stress),'_a',str(signs_asthen_visc),'.plastic_bending'])
	rms_name=''.join(['predictions/new/linear/rms_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.y',str(rms_yield_stress),'_a',str(rms_asthen_visc),'.plastic_bending'])
elif formulation == 3: # just slab pull (with prefactor)
	plot_name=''.join(['plots/new/linear/misfits_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.just-slab-pull.png'])
	signs_name=''.join(['predictions/new/linear/signs_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.pre',str(signs_pre),'_a',str(signs_asthen_visc),'.just-slab-pull'])
	rms_name=''.join(['predictions/new/linear/rms_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.pre',str(rms_pre),'_a',str(rms_asthen_visc),'.just-slab-pull'])
elif formulation == 4: # power law, viscous bending
	plot_name=''.join(['plots/new/powerlaw/misfits_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.power-law',str(trans_strain_rate),'_viscous_bending.png'])
	signs_name=''.join(['predictions/new/powerlaw/signs_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.l',str(signs_lith_visc),'_a10e',str(signs_asthen_visc),'.power-law',str(trans_strain_rate),'_viscous_bending'])
	rms_name=''.join(['predictions/new/powerlaw/rms_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.l',str(rms_lith_visc),'_a10e',str(rms_asthen_visc),'.power-law',str(trans_strain_rate),'_viscous_bending'])
elif formulation == 5: # power law, plastic bending
	plot_name=''.join(['plots/new/powerlaw/misfits_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.power-law',str(trans_strain_rate),'_plastic-bending.png'])   
	signs_name=''.join(['predictions/new/powerlaw/signs_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.l',str(signs_lith_visc),'_a10e',str(signs_asthen_visc),'.power-law',str(trans_strain_rate),'_plastic_bending'])
	rms_name=''.join(['predictions/new/powerlaw/rms_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.l',str(rms_lith_visc),'_a10e',str(rms_asthen_visc),'.power-law',str(trans_strain_rate),'_plastic_bending'])
elif formulation == 6: # power law, just slab pull (with prefactor)
	plot_name=''.join(['plots/new/powerlaw/misfits_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.power-law',str(trans_strain_rate),'_just-slab-pull.png'])
	signs_name=''.join(['predictions/new/powerlaw/signs_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.pre',str(signs_pre),'_a10e',str(signs_asthen_visc),'.power-law',str(trans_strain_rate),'_just-slab-pull'])
	rms_name=''.join(['predictions/new/powerlaw/rms_',str(vt_ref),'model',DP_string,PSP_string,RP_string,'.pre',str(rms_pre),'_a10e',str(rms_asthen_visc),'.power-law',str(trans_strain_rate),'_just-slab-pull'])	

plt.savefig(plot_name, bbox_inches='tight',dpi=400)

# plot map
rms_txt_name=''.join([rms_name,'.txt'])
np.savetxt(rms_txt_name, rms_predicted_vts, fmt='%.4f')  
np.savetxt('tmp/rms_separated.txt', vsps, fmt='%.4f')  

vt_observed=''.join(['tnew.',str(vt_ref),'.dat'])  
FNULL = open(os.devnull, 'w')
subprocess.check_call(['./plot_trench_motions.sh',rms_name,vt_observed,'tmp/rms_separated','1'])#,stdout=FNULL, stderr=subprocess.STDOUT)
print "output: %s" % plot_name

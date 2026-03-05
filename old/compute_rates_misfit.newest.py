#!/usr/bin/python

import scipy, math, os, numpy as np
import subprocess, sys
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from functions import compute_plate_buoyancy, compute_plate_isotherm
from functions import compute_vsp, compute_vsp_withDP

# vt data to compare
vt_ref=str(sys.argv[1])    				# hs3, nnr, sa
formulation=int(sys.argv[2])			# 1 = regular, 2 = plastic bending, 3 = slab pull pre-factor, 4 = regular, power-law, 5 = plastic, power-law, 6 = pre-factor, power-law
include_DP=int(sys.argv[3])				# 1 = include DP force, 0 = do not
DP_ref=float(sys.argv[4]) 				# DP values from analytical computations: free slip base: avg DP_0 = 18.3, max DP_0 = 23.5, no slip: avg DP_0 = 73.1, max DP_0 = 93.9 MPa
trans_strain_rate=float(sys.argv[5]) 	# typical value: 1e-13 s-1

# calculation parameters
calc_slabL_using_dip = 0
depth_cutoff=0		# 0 = all slabs, 1 = just shallow slabs
const_slab_depth=0  # 0 = use lallemand depths, 1 = all slabs go to 660 km
side_cutoff=0 		# 0 = all, 1 = just segments not at slab sides
use_avg_Rmin=0      # 1 = for segments without an Rmin, use the global average.
limit_max_age=0     # limit maximum subducting plate age/buoyancy? (think sea-floor flattening...)
composite=1         # for power-law: 0 = effective visc is power-law visc, 1 = effective visc is composite visc
					# (currently, just for include_DP=1)
					
# search parameters
lith_visc_min = 20; lith_visc_max = 25
asth_visc_min = 17; asth_visc_max = 22
yield_min = 100e6; 	yield_max = 10000e6;
pre_min = 0.0; 		pre_max = 1.0;

# reference parameters
vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s  
h = 200e3; # m
n = 3.5;

# DP stuff (SI units), for formulation == 8
visc_asthen_ref =  3e20
w_ref = 2224e3 * 2 # [m], factor 2 to convert from half width to full width,  (2224km = 20 degrees)
trenchv_ref = 5.0 * vel_converter # [cm/yr] -> [m/s], positive (as here) = trench retreat

# parameter space and error matrices
asthen_viscs = np.linspace(asth_visc_min,asth_visc_max,51)
lith_viscs = np.linspace(lith_visc_min,lith_visc_max,51)
yield_stresses = np.linspace(yield_min,yield_max,51)
prefactors = np.linspace(pre_min,pre_max,51)
asthen_visc = np.zeros((len(lith_viscs),len(asthen_viscs)))
lith_visc = np.zeros((len(lith_viscs),len(asthen_viscs)))
yield_stress = np.zeros((len(yield_stresses),len(asthen_viscs)))
pre_values = np.zeros((len(prefactors),len(asthen_viscs)))
rms = np.zeros((len(lith_viscs),len(asthen_viscs)))
sign = np.zeros((len(lith_viscs),len(asthen_viscs)))
for i in range(0,len(lith_viscs)):
	asthen_visc[i,:] = asthen_viscs[:];
	lith_visc[:,i] = lith_viscs[:];
	yield_stress[:,i] = yield_stresses[:];
	pre_values[:,i] = prefactors[:]

data = np.genfromtxt('data/Lallemand_et_al-2005_G3_dataset_WithTrenchWidths.txt') # see spreadsheet for column details
data_vt =np.genfromtxt('data/vts_hs3-nnr-sa.txt') # lat, lon, vtn: hs3, vtn: nnr, vtn: sa
if vt_ref == 'hs3':
   vt_col = 2;
elif vt_ref == 'nnr':
   vt_col = 3;
else:
   vt_col = 4; 

if depth_cutoff == 0:
	depth_cutoff_km = 100000.0;
elif depth_cutoff == 1:
	depth_cutoff_km = 670.0

if side_cutoff == 0:
	side_cutoff_val = 10;
elif side_cutoff == 1:
	side_cutoff_val = 0.5;

Rmin_mean = "%.2f" % np.nanmean(data[:,13])
nan_inds = np.where(np.isnan(data[:,13]))
if use_avg_Rmin == 1:
	data[nan_inds,13] = Rmin_mean

if limit_max_age == 1:
	max_age = 90.
else: 
	max_age = 500.

num = 0; Lsp_tot = 0;
for i in range(0,data.shape[0]):
	if np.isnan(data[i,26]) == False and np.isnan(data[i,6]) == False and data[i,7] <= depth_cutoff_km and data[i,24] < side_cutoff_val and np.isnan(data[i,8]) == False \
		and np.isnan(data[i,13]) == False and np.isnan(data[i,20]) == False:
		num = num + 1
		Lsp_tot = Lsp_tot + (data[i,26] * 1e3)
Lsp_avg = Lsp_tot/num

Lsp=np.zeros((num,1));			dip=np.zeros((num,1))
slabD=np.zeros((num,1));		vc=np.zeros((num,1))
Rmin=np.zeros((num,1));			age=np.zeros((num,1))
vt_actual=np.zeros((num,1));	slabL=np.zeros((num,1))
latlon=np.zeros((num,2));		azims=np.zeros((num,1))
w = np.zeros((num,1));			slabL_buoy=np.zeros((num,1))

n = 0
for i in range(0,data.shape[0]):
	if np.isnan(data[i,26]) == False and np.isnan(data[i,6]) == False and data[i,7] <= depth_cutoff_km and data[i,24] < side_cutoff_val and np.isnan(data[i,8]) == False \
		and np.isnan(data[i,13]) == False and np.isnan(data[i,20]) == False:
		
		azims[n,0] = data[i,4]
		Lsp[n,0] = data[i,26] * 1e3

		dip[n,0] = data[i,6]
		if const_slab_depth == 1:
			slabD[n,0] = 660e3
		else:
			slabD[n,0] = data[i,7] * 1e3

		if calc_slabL_using_dip == 1:
			slabL[n,0] = slabD[n,0]/np.tan(np.deg2rad(dip[n,0]))
		else:
			slabL[n,0] = data[i,8] * 1e3
		slabL_buoy[n,0] = data[i,8] * 1e3

		vc[n,0] = (data[i,9] / 10.0) * vel_converter
		Rmin[n,0] = data[i,13] * 1e3
		if data[i,20] > max_age:
			age[n,0] = max_age
		else:
			age[n,0] = data[i,20]
		w[n,0] = data[i,25] * 1e3
		vt_actual[n,0] = data_vt[i,vt_col] / 10.0
		latlon[n,0]=data_vt[i,0]; latlon[n,1]=data_vt[i,1];  
		n = n + 1

H = np.zeros((num,1))
oceanic_buoy = np.zeros((num,1))
for i in range(0,n):
	H[i] = compute_plate_isotherm(age[i],1300,1e-6,1100) * 1e3; # [m]
	oceanic_buoy[i] = compute_plate_buoyancy(age[i,0],1300,1e-6,3e-5,3300); # [kg/m2]

rms_min = 10.0
num_sign_matches_max = 0

for k in range(0,len(lith_viscs)):

	for j in range(0,len(asthen_viscs)):


		stored_rms_vals = np.zeros((num,1));
		stored_sign_vals = np.zeros((num,1));

		visc_asthen = 10**asthen_visc[k,j];
		visc_lith = 10**lith_visc[k,j];
		yield_sigma = yield_stress[k,j]; 
		pre = pre_values[k,j]

		if formulation == 1 or formulation == 4:
			print "lith, asthen visc: %e, %e" % (10**lith_visc[k,j],10**asthen_visc[k,j])

		# power-law
		if formulation >= 4:
			vsp_estimate=np.zeros((num,1));
			for i in range(0,len(vt_actual)):
				if include_DP == 1:
					vsp_estimate[i] = compute_vsp_withDP(formulation,vc[i],h,visc_asthen,visc_lith,H[i],Lsp[i],Rmin[i],slabL[i],slabL_buoy[i],dip[i],oceanic_buoy[i],DP_ref,visc_asthen_ref,\
						w_ref,trenchv_ref,w[i],slabD[i],yield_sigma,n,pre,trans_strain_rate,composite)
				else:
					vsp_estimate[i] = compute_vsp(formulation,vc[i],h,visc_asthen,visc_lith,H[i],Lsp[i],Rmin[i],slabL[i],dip[i],oceanic_buoy[i],yield_sigma,trans_strain_rate,n,pre)
		# non power-law			
		else:
			if include_DP == 1:
				vsp_estimate = compute_vsp_withDP(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,slabL_buoy[i],dip,oceanic_buoy,DP_ref,visc_asthen_ref,\
					w_ref,trenchv_ref,w,slabD,yield_sigma,n,pre,trans_strain_rate,composite)
			else:
				vsp_estimate = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,dip,oceanic_buoy,yield_sigma,trans_strain_rate,n,pre)

		vt_estimate = (vc/vel_converter) - vsp_estimate
		rms_sum = 0; num_sign_matches = 0;
		for i in range(0,n):
			rms_sum = rms_sum + (vt_estimate[i] - vt_actual[i])**2;
			stored_rms_vals[i] = np.abs(vt_estimate[i] - vt_actual[i])
			if np.sign(vt_estimate[i]) == np.sign(vt_actual[i]):
				num_sign_matches = num_sign_matches + 1
				stored_sign_vals[i] = 1
		rms_act = np.sqrt(rms_sum / float(n))
		rms[k,j] = rms_act;
		sign[k,j] = num_sign_matches;

		if rms_act < rms_min:
			rms_min = rms_act;
			rms_lith_visc = np.log10(visc_lith);
			rms_asthen_visc = np.log10(visc_asthen);
			rms_yield_stress = yield_sigma/1e6;
			rms_pre = pre;
			rms_predicted_vts = np.concatenate((latlon,vt_estimate,azims), axis=1)
			rms_separated = np.concatenate((latlon,stored_rms_vals), axis=1)

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

fig = plt.figure()

ax = fig.add_subplot(211)
if formulation == 2 or formulation == 5:
	im1 = ax.imshow(100.0 * (sign/n),  cmap=cm.RdYlGn, extent=[asth_visc_min, asth_visc_max, yield_min/1e6, yield_max/1e6], vmax=100.0, vmin=50.0,aspect=(asth_visc_max-asth_visc_min)/((yield_max/1e6)-(yield_min/1e6)))
	ax.set_ylabel("$\mathregular{\sigma_{Y}}$  [MPa]")
elif formulation == 1 or formulation == 4:
	im1 = ax.imshow(100.0 * (sign/n),  cmap=cm.RdYlGn, origin='lower', extent=[asth_visc_min, asth_visc_max, lith_visc_min, lith_visc_max], vmax=100.0, vmin=50.0,aspect=(asth_visc_max-asth_visc_min)/(lith_visc_max-lith_visc_min))
	ax.set_ylabel("$\mathregular{\eta_{L}}$  [Pa.s]")
elif formulation == 3 or formulation == 6:
	im1 = ax.imshow(100.0 * (sign/n),  cmap=cm.RdYlGn, origin='lower', extent=[asth_visc_min, asth_visc_max, pre_min, pre_max], vmax=100.0, vmin=50.0,aspect=(asth_visc_max-asth_visc_min)/(pre_max-pre_min))
	ax.set_ylabel("K")
ax.text(0.07,0.88,"%s with correct sign" % ('%'),size=12.5, color="black",transform = ax.transAxes)
plt.colorbar(im1)
annot_string = ''.join(['best: ',str(num_sign_matches_max),'/',str(n),' signs, RMS = ',str("%.2f" % rms_min),' cm/yr']);
plt.annotate(annot_string, xy=(0.035, 1.05), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7)

ax = fig.add_subplot(212)
if formulation == 2 or formulation == 5:
	im2 = ax.imshow(rms,  cmap=cm.RdYlGn,  extent=[asth_visc_min, asth_visc_max, yield_min/1e6, yield_max/1e6], vmax=10.0, vmin=0.0,aspect=(asth_visc_max-asth_visc_min)/((yield_max/1e6)-(yield_min/1e6)))
	ax.set_ylabel("$\mathregular{\sigma_{Y}}$  [MPa]")
elif formulation == 1 or formulation == 4:
	im2 = ax.imshow(rms,  cmap=cm.RdYlGn, origin='lower', extent=[asth_visc_min, asth_visc_max, lith_visc_min, lith_visc_max], vmax=10.0, vmin=0.0,aspect=(asth_visc_max-asth_visc_min)/(lith_visc_max-lith_visc_min))
	ax.set_ylabel("$\mathregular{\eta_{L}}$  [Pa.s]")
elif formulation == 3 or formulation == 6:
	im2 = ax.imshow(rms,  cmap=cm.RdYlGn, origin='lower', extent=[asth_visc_min, asth_visc_max, pre_min, pre_max], vmax=10.0, vmin=0.0,aspect=(asth_visc_max-asth_visc_min)/(pre_max-pre_min))
	ax.set_ylabel("K")
ax.set_xlabel("$\mathregular{\eta_{A}}$  [Pa.s]")
ax.text(0.07,0.88,'RMS',size=12.5, color="black",transform = ax.transAxes)
plt.colorbar(im2)

if depth_cutoff == 0:
	depths_string=''
else:
	depths_string='_no-deep-slabs'

if side_cutoff == 0:
	side_string=''
else:
	side_string='_non-edge-slabs'

if include_DP == 0:
	DP_string=''
else:
	DP_string = ''.join(['.DP',str("%.3g" % DP_ref),'MPa'])

if composite == 1:
	comp_string=''
else:
	comp_string='NoComp'

if formulation == 1:  # viscous bending
	plot_name=''.join(['plots/new/misfits_',str(vt_ref),'model',DP_string,depths_string, side_string,'.viscous_bending.png'])
	signs_name=''.join(['predictions/new/signs_',str(vt_ref),'model',DP_string,depths_string, side_string,'.l',str(signs_lith_visc),'_a10e',str(signs_asthen_visc),'.viscous_bending'])
	rms_name=''.join(['predictions/new/rms_',str(vt_ref),'model',DP_string,depths_string, side_string,'.l',str(rms_lith_visc),'_a10e',str(rms_asthen_visc),'.viscous_bending'])
elif formulation == 2: # plastic bending
	plot_name=''.join(['plots/new/misfits_',str(vt_ref),'model',DP_string,depths_string, side_string,'.plastic_bending.png'])
	signs_name=''.join(['predictions/new/signs_',str(vt_ref),'model',DP_string,depths_string, side_string,'.y',str(signs_yield_stress),'_a',str(signs_asthen_visc),'.plastic_bending'])
	rms_name=''.join(['predictions/new/rms_',str(vt_ref),'model',DP_string,depths_string, side_string,'.y',str(rms_yield_stress),'_a',str(rms_asthen_visc),'.plastic_bending'])
elif formulation == 3: # just slab pull (with prefactor)
	plot_name=''.join(['plots/new/misfits_',str(vt_ref),'model',DP_string,depths_string, side_string,'.just-slab-pull.png'])
	signs_name=''.join(['predictions/new/signs_',str(vt_ref),'model',DP_string,depths_string, side_string,'.pre',str(signs_pre),'_a',str(signs_asthen_visc),'.just-slab-pull'])
	rms_name=''.join(['predictions/new/rms_',str(vt_ref),'model',DP_string,depths_string, side_string,'.pre',str(rms_pre),'_a',str(rms_asthen_visc),'.just-slab-pull'])
elif formulation == 4: # power law, viscous bending
	plot_name=''.join(['plots/new/misfits_',str(vt_ref),'model',DP_string,depths_string, side_string,'.power-law',str(trans_strain_rate),str(comp_string),'_viscous_bending.png'])
	signs_name=''.join(['predictions/new/signs_',str(vt_ref),'model',DP_string,depths_string, side_string,'.l',str(signs_lith_visc),'_a10e',str(signs_asthen_visc),'.power-law',str(trans_strain_rate),str(comp_string),'_viscous_bending'])
	rms_name=''.join(['predictions/new/rms_',str(vt_ref),'model',DP_string,depths_string, side_string,'.l',str(rms_lith_visc),'_a10e',str(rms_asthen_visc),'.power-law',str(trans_strain_rate),str(comp_string),'_viscous_bending'])
elif formulation == 5: # power law, plastic bending
	plot_name=''.join(['plots/new/misfits_',str(vt_ref),'model',DP_string,depths_string, side_string,'.power-law',str(trans_strain_rate),str(comp_string),'_plastic-bending.png'])   
	signs_name=''.join(['predictions/new/signs_',str(vt_ref),'model',DP_string,depths_string, side_string,'.l',str(signs_lith_visc),'_a10e',str(signs_asthen_visc),'.power-law',str(trans_strain_rate),str(comp_string),'_plastic_bending'])
	rms_name=''.join(['predictions/new/rms_',str(vt_ref),'model',DP_string,depths_string, side_string,'.l',str(rms_lith_visc),'_a10e',str(rms_asthen_visc),'.power-law',str(trans_strain_rate),str(comp_string),'_plastic_bending'])
elif formulation == 6: # power law, just slab pull (with prefactor)
	plot_name=''.join(['plots/new/misfits_',str(vt_ref),'model',DP_string,depths_string, side_string,'.power-law',str(trans_strain_rate),str(comp_string),'_just-slab-pull.png'])
	signs_name=''.join(['predictions/new/signs_',str(vt_ref),'model',DP_string,depths_string, side_string,'.pre',str(signs_pre),'_a10e',str(signs_asthen_visc),'.power-law',str(trans_strain_rate),str(comp_string),'_just-slab-pull'])
	rms_name=''.join(['predictions/new/rms_',str(vt_ref),'model',DP_string,depths_string, side_string,'.pre',str(rms_pre),'_a10e',str(rms_asthen_visc),'.power-law',str(trans_strain_rate),str(comp_string),'_just-slab-pull'])	

plt.savefig(plot_name, bbox_inches='tight',dpi=400)

# plot map
signs_txt_name=''.join([signs_name,'.txt'])
rms_txt_name=''.join([rms_name,'.txt'])
np.savetxt(signs_txt_name, signs_predicted_vts, fmt='%.4f')  
np.savetxt(rms_txt_name, rms_predicted_vts, fmt='%.4f')  
np.savetxt('tmp/signs_separated.txt', signs_separated, fmt='%.4f')  
np.savetxt('tmp/rms_separated.txt', rms_separated, fmt='%.4f')  

vt_observed=''.join(['tnew.',str(vt_ref),'.dat'])  
FNULL = open(os.devnull, 'w')
subprocess.check_call(['./plot_trench_motions.sh',signs_name,vt_observed,'tmp/signs_separated','0'],stdout=FNULL, stderr=subprocess.STDOUT)
subprocess.check_call(['./plot_trench_motions.sh',rms_name,vt_observed,'tmp/rms_separated','1'],stdout=FNULL, stderr=subprocess.STDOUT)
print "output: %s" % plot_name

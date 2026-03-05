#!/usr/bin/python

import scipy, math, os, numpy as np
import subprocess, pygplates
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from functions import compute_plate_buoyancy, compute_plate_isotherm
from functions import compute_vsp

# options
vt_ref='sa'    		# hs3, nnr, sa
formulation=7		# 4 = power-law, 5 = plastic + power-law, 7 = power-law version 2 
depth_cutoff=0		# 0 = all, 1 = just shallow
side_cutoff=0 		# 0 = all, 1 = just segments not at slab sides
use_avg_Lsp=0		# 0 = regular, 1 = use average Lsp (gives better match...)

lith_visc_min = 20; lith_visc_max = 25
asth_visc_min = 17; asth_visc_max = 22
yield_min = 100e6; yield_max = 10000e6;
trans_strain_rate_min = -14; trans_strain_rate_max = -11;
pre_min = 0.0; pre_max = 1.0;
v_trans_min = 5.0; v_trans_max = 6.0;

# parameter space and error matrices
trans_strain_rates = np.linspace(trans_strain_rate_min,trans_strain_rate_max,11)
trans_vels=np.linspace(v_trans_min,v_trans_max,1)

asthen_viscs = np.linspace(asth_visc_min,asth_visc_max,51)
lith_viscs = np.linspace(lith_visc_min,lith_visc_max,51)
yield_stresses = np.linspace(yield_min,yield_max,51)
prefactors = np.linspace(pre_min,pre_max,51)

asthen_visc = np.zeros((len(lith_viscs),len(asthen_viscs)))
lith_visc = np.zeros((len(lith_viscs),len(asthen_viscs)))
yield_stress = np.zeros((len(yield_stresses),len(asthen_viscs)))
pre_values = np.zeros((len(prefactors),len(asthen_viscs)))

for i in range(0,len(lith_viscs)):
	asthen_visc[i,:] = asthen_viscs[:];
	lith_visc[:,i] = lith_viscs[:];
	yield_stress[:,i] = yield_stresses[:];
	pre_values[:,i] = prefactors[:]

# reference parameters, SI units
vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s  
h = 200e3; exp = 3.5;

data = np.genfromtxt('data/Lallemand_et_al-2005_G3_dataset.txt') # see spreadsheet for column details
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

num = 0; Lsp_tot = 0;
for i in range(0,data.shape[0]):
	if np.isnan(data[i,3]) == False and np.isnan(data[i,6]) == False and data[i,7] <= depth_cutoff_km and data[i,24] < side_cutoff_val and np.isnan(data[i,8]) == False and np.isnan(data[i,9]) == False and np.isnan(data[i,13]) == False and np.isnan(data[i,20]) == False and np.isnan(data[i,15]) == False:
		num = num + 1
		Lsp_tot = Lsp_tot + (data[i,3] * 1e3)
Lsp_avg = Lsp_tot/num

# data=np.mat(data)
Lsp=np.zeros((num,1));		 dip=np.zeros((num,1))
D=np.zeros((num,1));         vc=np.zeros((num,1))
Rmin=np.zeros((num,1));      age=np.zeros((num,1))
vt_actual=np.zeros((num,1)); slabL=np.zeros((num,1))
latlon=np.zeros((num,2));	 azims=np.zeros((num,1))

n = 0
for i in range(0,data.shape[0]):
	if np.isnan(data[i,3]) == False and np.isnan(data[i,6]) == False and data[i,7] <= depth_cutoff_km and data[i,24] < side_cutoff_val and  np.isnan(data[i,8]) == False and np.isnan(data[i,9]) == False and np.isnan(data[i,13]) == False and np.isnan(data[i,20]) == False and np.isnan(data[i,15]) == False:
		azims[n,0] = data[i,4]
		if use_avg_Lsp == 1:
			Lsp[n,0] = Lsp_avg
		else:
			Lsp[n,0] = data[i,3] * 1e3
		dip[n,0] = data[i,6]
		slabL[n,0] = data[i,8] * 1e3
		vc[n,0] = (data[i,9] / 10.0) * vel_converter
		Rmin[n,0] = data[i,13] * 1e3
		age[n,0] = data[i,20]
		vt_actual[n,0] = data_vt[i,vt_col] / 10.0
		latlon[n,0]=data_vt[i,0]; latlon[n,1]=data_vt[i,1];  
		n = n + 1

H = np.zeros((num,1))
oceanic_buoy = np.zeros((num,1))
for i in range(0,n):
	H[i] = compute_plate_isotherm(age[i],1300,1e-6,1100) * 1e3; # [m]
	oceanic_buoy[i] = compute_plate_buoyancy(age[i,0],1300,1e-6,3e-5,3300); # [kg/m2]

if formulation == 4 or formulation == 5:
	num_iterations=len(trans_strain_rates)
elif formulation == 7:
	num_iterations=len(trans_vels)

for e in range(0,num_iterations):

	if formulation == 4 or formulation == 5:
		trans_strain_rate = 10 ** trans_strain_rates[e]
		v_trans = 0.0
		print "Parameter space search for transition strain rate = %.4e" % trans_strain_rate
	elif formulation == 7:
		trans_strain_rate = 0.0
		v_trans = trans_vels[e]
		print "Parameter space search for transition velocity = %.2f cm/yr" % v_trans

	
	rms = np.zeros((len(lith_viscs),len(asthen_viscs)))
	sign = np.zeros((len(lith_viscs),len(asthen_viscs)))
	for k in range(0,len(lith_viscs)):
		for j in range(0,len(asthen_viscs)):
			visc_asthen = 10**asthen_visc[k,j];
			visc_lith = 10**lith_visc[k,j];
			yield_sigma = yield_stress[k,j]; 
			pre = pre_values[k,j]

			vsp_estimate = 0.0*vt_actual
			for i in range(0,len(vt_actual)):
				vsp_estimate[i] = compute_vsp(formulation,vc[i],h,visc_asthen,visc_lith,H[i],Lsp[i],Rmin[i],slabL[i],dip[i],oceanic_buoy[i],yield_sigma,trans_strain_rate,exp,pre,v_trans)
			vt_estimate = (vc/vel_converter) - vsp_estimate;

			rms_sum = 0; num_sign_matches = 0;
			for i in range(0,n):
				# if vt_estimate[i] > 20:
				# 	vt_estimate[i] == 20;
				# if vt_estimate[i] < -20:
				# 	vt_estimate[i] = -20;
				rms_sum = rms_sum + (vt_estimate[i] - vt_actual[i])**2;
				if np.sign(vt_estimate[i]) == np.sign(vt_actual[i]):
					num_sign_matches = num_sign_matches + 1
			rms_act = np.sqrt(rms_sum / float(n))
			rms[k,j] = rms_act;
			sign[k,j] = num_sign_matches;

	plt.close("all")

	fig = plt.figure()

	ax = fig.add_subplot(211)
	if formulation == 5:
		im1 = ax.imshow(100.0 * (sign/n),  cmap=cm.RdYlGn, extent=[asth_visc_min, asth_visc_max, yield_min/1e6, yield_max/1e6], vmax=100.0, vmin=0,aspect=(asth_visc_max-asth_visc_min)/((yield_max/1e6)-(yield_min/1e6)))
		ax.set_ylabel("$\mathregular{\sigma_{Y}}$  [MPa]")
	elif formulation == 4:
		im1 = ax.imshow(100.0 * (sign/n),  cmap=cm.RdYlGn, origin='lower', extent=[asth_visc_min, asth_visc_max, lith_visc_min, lith_visc_max], vmax=100.0, vmin=0)
		ax.set_ylabel("$\mathregular{\eta_{L}}$  [Pa.s]")
	elif formulation == 7:
		im1 = ax.imshow(100.0 * (sign/n),  cmap=cm.RdYlGn, origin='lower', extent=[asth_visc_min, asth_visc_max, v_trans_min, v_trans_max], vmax=100.0, vmin=0,aspect=(asth_visc_max-asth_visc_min)/(v_trans_max-v_trans_min))
		ax.set_ylabel("$\mathregular{v_{T}}$  [cm/yr]")
	ax.text(0.07,0.88,"%s with correct sign" % ('%'),size=12.5, color="black",transform = ax.transAxes)
	plt.colorbar(im1)

	ax = fig.add_subplot(212)
	if formulation == 3:
		im2 = ax.imshow(rms,  cmap=cm.RdYlGn,  extent=[asth_visc_min, asth_visc_max, yield_min/1e6, yield_max/1e6], vmax=40.0, vmin=0.0,aspect=(asth_visc_max-asth_visc_min)/((yield_max/1e6)-(yield_min/1e6)))
		ax.set_ylabel("$\mathregular{\sigma_{Y}}$  [MPa]")
	elif formulation == 4:
		im2 = ax.imshow(rms,  cmap=cm.RdYlGn, origin='lower', extent=[asth_visc_min, asth_visc_max, lith_visc_min, lith_visc_max], vmax=40.0, vmin=0.0)
		ax.set_ylabel("$\mathregular{\eta_{L}}$  [Pa.s]")
	elif formulation == 7:
		im2 = ax.imshow(rms,  cmap=cm.RdYlGn, origin='lower', extent=[asth_visc_min, asth_visc_max, v_trans_min, v_trans_max], vmax=40.0, vmin=0.0, aspect=(asth_visc_max-asth_visc_min)/(v_trans_max-v_trans_min))
		ax.set_ylabel("$\mathregular{v_{T}}$  [cm/yr]")
	ax.set_xlabel("$\mathregular{\eta_{A}}$  [Pa.s]")
	ax.text(0.07,0.88,'RMS',size=12.5, color="black",transform = ax.transAxes)
	plt.colorbar(im2)

	if depth_cutoff == 0:
		depths_string=''
	else:
		depths_string='_deep-slabs'

	if side_cutoff == 0:
		side_string=''
	else:
		side_string='_non-edge-slabs'

	if formulation == 4:
		plot_name=''.join(['plots/power-law_misfits/misfits_',str(vt_ref),'model',depths_string, side_string,'_trans10e',str(trans_strain_rates[e]),'s-1.power-law.png'])
	elif formulation == 5:
		plot_name=''.join(['plots/power-law_misfits/misfits_',str(vt_ref),'model',depths_string, side_string,'_trans10e',str(trans_strain_rates[e]),'s-1.plastic-bending_and_power-law.png'])
	elif formulation == 7:
		 plot_name=''.join(['plots/power-law_misfits/misfits_',str(vt_ref),'model',depths_string, side_string,'_vtrans',str(v_trans),'cmyr.power-law-v2.png'])  

	plt.savefig(plot_name, bbox_inches='tight',dpi=400)




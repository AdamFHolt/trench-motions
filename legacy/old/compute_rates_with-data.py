#!/usr/bin/python

import scipy, math, os, numpy as np
import subprocess, pygplates
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FuncFormatter
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from functions import compute_plate_buoyancy, compute_plate_isotherm

data = np.genfromtxt('Lallemand_et_al-2005_G3_dataset.txt') # see spreadsheet for column details

# reference parameters, SI units
vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s 
D = 660e3; g = 9.81; dip = 60.0; Lsp = 5000e3; age = 60.0; 
Rmin = 250e3; 
visc_lith = 500.0e20; visc_asthen = 1.0e21; h = 200e3;
oceanic_buoy = compute_plate_buoyancy(age,1300,1e-6,3e-5,3300); 
H = compute_plate_isotherm(age,1300,1e-6,1300) * 1e3;
vc = 1.0 * vel_converter;

# slab_buoy = (D * oceanic_buoy * g)/np.sin(np.deg2rad(dip));
# bending = (2./3.) * (H**3/Rmin**3) * visc_lith * vc;
# slab_shear = 2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (visc_asthen/h);
# pressure_force = (D * oceanic_buoy * g) * np.cos(np.deg2rad(dip)); # NOT INCLUDED

### effect of plate length (Lsp)
plate_length_effect = np.zeros((101,4))
plate_length_effect[:,0] = np.linspace(400e3,12000e3,101)
vc = 1.0 * vel_converter; 
plate_length_effect[:,1] = (h/(visc_asthen*plate_length_effect[:,0])) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (visc_asthen/h)) ) * (1.0/vel_converter)
vc = 5.0 * vel_converter;
plate_length_effect[:,2] = (h/(visc_asthen*plate_length_effect[:,0])) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (visc_asthen/h))  ) * (1.0/vel_converter)
vc = 10.0 * vel_converter;
plate_length_effect[:,3] = (h/(visc_asthen*plate_length_effect[:,0])) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (visc_asthen/h))  ) * (1.0/vel_converter)


### effect of slab dip
slab_dip_effect = np.zeros((101,4))
slab_dip_effect[:,0] = np.linspace(30,100,101)
vc = 1.0 * vel_converter; 
slab_dip_effect[:,1] = (h/(visc_asthen*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(slab_dip_effect[:,0]))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(slab_dip_effect[:,0]))) * (visc_asthen/h)) ) * (1.0/vel_converter)
vc = 5.0 * vel_converter; 
slab_dip_effect[:,2] = (h/(visc_asthen*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(slab_dip_effect[:,0]))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(slab_dip_effect[:,0]))) * (visc_asthen/h)) ) * (1.0/vel_converter)
vc = 10.0 * vel_converter; 
slab_dip_effect[:,3] = (h/(visc_asthen*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(slab_dip_effect[:,0]))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(slab_dip_effect[:,0]))) * (visc_asthen/h)) ) * (1.0/vel_converter)

### effect of slab age
slab_age_effect = np.zeros((101,4))
slab_age_effect[:,0] = np.linspace(0,150,101)
H = np.zeros((101)); oceanic_buoy = np.zeros((101)); 
for i in range(0,101):
	H[i] = compute_plate_isotherm(slab_age_effect[i,0],1300,1e-6,1300) * 1e3; # [m]
	oceanic_buoy[i] = compute_plate_buoyancy(slab_age_effect[i,0],1300,1e-6,3e-5,3300); # [kg/m2]
vc = 1.0 * vel_converter; 
slab_age_effect[:,1] = (h/(visc_asthen*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (visc_asthen/h)) ) * (1.0/vel_converter)
vc = 5.0 * vel_converter; 
slab_age_effect[:,2] = (h/(visc_asthen*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (visc_asthen/h)) ) * (1.0/vel_converter)
vc = 10.0 * vel_converter; 
slab_age_effect[:,3] = (h/(visc_asthen*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (visc_asthen/h)) ) * (1.0/vel_converter)

oceanic_buoy = compute_plate_buoyancy(60.0,1300,1e-6,3e-5,3300); 
H = compute_plate_isotherm(60.0,1300,1e-6,1300) * 1e3;
### effect of asthen viscosity / thickness
asthen_effect = np.zeros((101,4))
asthen_effect[:,0] = np.linspace(2e20,20e20,101)
vc = 1.0 * vel_converter; 
asthen_effect[:,1] = (h/(asthen_effect[:,0]*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (asthen_effect[:,0]/h)) ) * (1.0/vel_converter)
vc = 5.0 * vel_converter; 
asthen_effect[:,2] = (h/(asthen_effect[:,0]*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (asthen_effect[:,0]/h)) ) * (1.0/vel_converter)
vc = 10.0 * vel_converter; 
asthen_effect[:,3] = (h/(asthen_effect[:,0]*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (asthen_effect[:,0]/h)) ) * (1.0/vel_converter)

### effect of slab depth
depth_effect = np.zeros((101,4))
depth_effect[:,0] = np.linspace(100e3,1200e3,101)
vc = 1.0 * vel_converter; 
depth_effect[:,1] = (h/(visc_asthen*Lsp)) * ( ((depth_effect[:,0] * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (depth_effect[:,0]/np.sin(np.deg2rad(dip))) * (visc_asthen/h)) ) * (1.0/vel_converter)
vc = 5.0 * vel_converter; 
depth_effect[:,2] = (h/(visc_asthen*Lsp)) * ( ((depth_effect[:,0] * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (depth_effect[:,0]/np.sin(np.deg2rad(dip))) * (visc_asthen/h)) ) * (1.0/vel_converter)
vc = 10.0 * vel_converter; 
depth_effect[:,3] = (h/(visc_asthen*Lsp)) * ( ((depth_effect[:,0] * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
	 (2.0 * vc * (depth_effect[:,0]/np.sin(np.deg2rad(dip))) * (visc_asthen/h)) ) * (1.0/vel_converter)

### effect of radius of curvature
radius_effect = np.zeros((101,4))
radius_effect[:,0] = np.linspace(50e3,700e3,101)
vc = 1.0 * vel_converter; 
radius_effect[:,1] = (h/(visc_asthen*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/radius_effect[:,0]**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (visc_asthen/h)) ) * (1.0/vel_converter)
vc = 5.0 * vel_converter; 
radius_effect[:,2] = (h/(visc_asthen*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/radius_effect[:,0]**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (visc_asthen/h)) ) * (1.0/vel_converter)
vc = 10.0 * vel_converter; 
radius_effect[:,3] = (h/(visc_asthen*Lsp)) * ( ((D * oceanic_buoy * g)/np.sin(np.deg2rad(dip))) - ((2./3.) * (H**3/radius_effect[:,0]**3) * visc_lith * vc) - \
	 (2.0 * vc * (D/np.sin(np.deg2rad(dip))) * (visc_asthen/h)) ) * (1.0/vel_converter)

## organize data for plotting

# data=np.mat(data)
data_plate_length=np.zeros((data.shape[0],3))
for i in range(0,data.shape[0]):
	if np.isnan(data[i,15]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,3]) == False:
		data_plate_length[i,0] = data[i,3];
		data_plate_length[i,1] = data[i,15];
		data_plate_length[i,2] = data[i,9];
data_plate_length = data_plate_length[~np.all(data_plate_length == 0, axis=1)]

data_slab_dip=np.zeros((data.shape[0],3))
for i in range(0,data.shape[0]):
	if np.isnan(data[i,15]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,6]) == False:
		data_slab_dip[i,0] = data[i,6];
		data_slab_dip[i,1] = data[i,15];
		data_slab_dip[i,2] = data[i,9];
data_slab_dip = data_slab_dip[~np.all(data_slab_dip == 0, axis=1)]

data_slab_age=np.zeros((data.shape[0],3))
for i in range(0,data.shape[0]):
	if np.isnan(data[i,15]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,20]) == False:
		data_slab_age[i,0] = data[i,20];
		data_slab_age[i,1] = data[i,15];
		data_slab_age[i,2] = data[i,9];
data_slab_age = data_slab_age[~np.all(data_slab_age == 0, axis=1)]

data_slab_depth=np.zeros((data.shape[0],3))
for i in range(0,data.shape[0]):
	if np.isnan(data[i,15]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,7]) == False:
		data_slab_depth[i,0] = data[i,7];
		data_slab_depth[i,1] = data[i,15];
		data_slab_depth[i,2] = data[i,9];
data_slab_depth = data_slab_depth[~np.all(data_slab_depth == 0, axis=1)]

data_Rmin=np.zeros((data.shape[0],3))
for i in range(0,data.shape[0]):
	if np.isnan(data[i,16]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,13]) == False:
		data_Rmin[i,0] = data[i,13];
		data_Rmin[i,1] = data[i,16]; # using heuret's vt to be consistent with heuret's Rmin
		data_Rmin[i,2] = data[i,9];
data_Rmin = data_Rmin[~np.all(data_Rmin == 0, axis=1)]

fig, axs = plt.subplots(3, 2, tight_layout=True)

axs[0,0].set_xlabel("$\mathregular{L_{SP}}$  [km]")
axs[0,0].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[0,0].set_ylim(-20,20); axs[0,0].set_xlim(400,12000)
axs[0,0].tick_params(axis='both', which='major', labelsize=10)
axs[0,0].tick_params(axis='both', which='minor', labelsize=8)
axs[0,0].plot(plate_length_effect[:,0]/1e3,0.0*plate_length_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[0,0].plot(plate_length_effect[:,0]/1e3,(1.0-plate_length_effect[:,1]), linewidth=2, color='orange')
axs[0,0].plot(plate_length_effect[:,0]/1e3,(5.0-plate_length_effect[:,2]), linewidth=2, color='red')
axs[0,0].plot(plate_length_effect[:,0]/1e3,(10.0-plate_length_effect[:,3]), linewidth=2, color='brown')
axs[0,0].text(0.0,1.07,'$\mathregular{V_{C}}$ = 1 cm/yr',size=9, color="orange",transform = axs[0,0].transAxes)
axs[0,0].text(0.325,1.07,'$\mathregular{V_{C}}$ = 5 cm/yr',size=9, color="red",transform = axs[0,0].transAxes)
axs[0,0].text(0.65,1.07,'$\mathregular{V_{C}}$ = 10 cm/yr',size=9, color="brown",transform = axs[0,0].transAxes)
colors=data_plate_length[:,2]/10.
cset1=axs[0,0].scatter(data_plate_length[:,0],data_plate_length[:,1]/10.,vmin=0,vmax=10,s=5,c=colors,cmap=plt.cm.OrRd,edgecolors='none')

axs[0,1].set_xlabel("$\mathregular{dip}$  [deg]")
axs[0,1].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[0,1].set_ylim(-20,20); axs[0,1].set_xlim(30,100)
axs[0,1].tick_params(axis='both', which='major', labelsize=10)
axs[0,1].tick_params(axis='both', which='minor', labelsize=8)
axs[0,1].plot(slab_dip_effect[:,0],0.0*slab_dip_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[0,1].plot(slab_dip_effect[:,0],(1.0-slab_dip_effect[:,1]), linewidth=2, color='orange')
axs[0,1].plot(slab_dip_effect[:,0],(5.0-slab_dip_effect[:,2]), linewidth=2, color='red')
axs[0,1].plot(slab_dip_effect[:,0],(10.0-slab_dip_effect[:,3]), linewidth=2, color='brown')
colors=data_slab_dip[:,2]/10.
cset1=axs[0,1].scatter(data_slab_dip[:,0],data_slab_dip[:,1]/10.,vmin=0,vmax=10,s=5,c=colors,cmap=plt.cm.OrRd,edgecolors='none')

axs[1,0].set_xlabel("$\mathregular{age_{SP}}}$  [Ma]")
axs[1,0].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[1,0].set_ylim(-20,20); axs[1,0].set_xlim(0,150); 
axs[1,0].tick_params(axis='both', which='major', labelsize=10)
axs[1,0].tick_params(axis='both', which='minor', labelsize=8)
axs[1,0].plot(slab_age_effect[:,0],0.0*slab_age_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[1,0].plot(slab_age_effect[:,0],(1.0-slab_age_effect[:,1]), linewidth=2, color='orange')
axs[1,0].plot(slab_age_effect[:,0],(5.0-slab_age_effect[:,2]), linewidth=2, color='red')
axs[1,0].plot(slab_age_effect[:,0],(10.0-slab_age_effect[:,3]), linewidth=2, color='brown')
colors=data_slab_age[:,2]/10.
cset1=axs[1,0].scatter(data_slab_age[:,0],data_slab_age[:,1]/10.,vmin=0,vmax=10,s=5,c=colors,cmap=plt.cm.OrRd,edgecolors='none')

#axs[1,1].set_xlabel("$\mathregular{visc_{A}}$  [Pa.s]")
axs[1,1].set_xlabel("$\mathregular{\eta_{A}}$  [Pa.s]")
axs[1,1].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[1,1].set_ylim(-20,20)
axs[1,1].tick_params(axis='both', which='major', labelsize=10)
axs[1,1].tick_params(axis='both', which='minor', labelsize=8)
axs[1,1].plot(asthen_effect[:,0],0.0*asthen_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[1,1].plot(asthen_effect[:,0],(1.0-asthen_effect[:,1]), linewidth=2, color='orange')
axs[1,1].plot(asthen_effect[:,0],(5.0-asthen_effect[:,2]), linewidth=2, color='red')
axs[1,1].plot(asthen_effect[:,0],(10.0-asthen_effect[:,3]), linewidth=2, color='brown')

axs[2,0].set_xlabel("$\mathregular{D_{SP}}$  [km]")
axs[2,0].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[2,0].set_ylim(-20,20); axs[2,0].set_xlim(100,1200); 
axs[2,0].tick_params(axis='both', which='major', labelsize=10)
axs[2,0].tick_params(axis='both', which='minor', labelsize=8)
axs[2,0].plot(depth_effect[:,0]/1e3,0.0*depth_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[2,0].plot(depth_effect[:,0]/1e3,(1.0-depth_effect[:,1]), linewidth=2, color='orange')
axs[2,0].plot(depth_effect[:,0]/1e3,(5.0-depth_effect[:,2]), linewidth=2, color='red')
axs[2,0].plot(depth_effect[:,0]/1e3,(10.0-depth_effect[:,3]), linewidth=2, color='brown')
colors=data_slab_depth[:,2]/10.
cset1=axs[2,0].scatter(data_slab_depth[:,0],data_slab_depth[:,1]/10.,vmin=0,vmax=10,s=5,c=colors,cmap=plt.cm.OrRd,edgecolors='none')

axs[2,1].set_xlabel("$\mathregular{R_{min}}$  [km]")
axs[2,1].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[2,1].set_ylim(-20,20)
axs[2,1].set_xlim(50,550); 
axs[2,1].tick_params(axis='both', which='major', labelsize=10)
axs[2,1].tick_params(axis='both', which='minor', labelsize=8)
axs[2,1].plot(radius_effect[:,0],0.0*radius_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[2,1].plot(radius_effect[:,0]/1e3,(1.0-radius_effect[:,1]), linewidth=2, color='orange')
axs[2,1].plot(radius_effect[:,0]/1e3,(5.0-radius_effect[:,2]), linewidth=2, color='red')
axs[2,1].plot(radius_effect[:,0]/1e3,(10.0-radius_effect[:,3]), linewidth=2, color='brown')
colors=data_Rmin[:,2]/10.
cset1=axs[2,1].scatter(data_Rmin[:,0],data_Rmin[:,1]/10.,vmin=0,vmax=10,s=5,c=colors,cmap=plt.cm.OrRd,edgecolors='none')

plt.savefig('test_with-data.png', bbox_inches='tight')



#!/usr/bin/python

import scipy, math, os, numpy as np
import subprocess, pygplates
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FuncFormatter
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from functions import compute_plate_buoyancy, compute_plate_isotherm
from functions import compute_vt

# reference parameters, SI units
vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s 
D = 660e3; g = 9.81; dip = 60.0; Lsp = 5000e3; age = 60.0; 
Rmin = 250e3; 
visc_lith = 500.0e20; visc_asthen = 1.0e21; h = 200e3;
oceanic_buoy = compute_plate_buoyancy(age,1300,1e-6,3e-5,3300); 
H = compute_plate_isotherm(age,1300,1e-6,1300) * 1e3;
vc = 1.0 * vel_converter;


### effect of plate length (Lsp)
plate_length_effect = np.zeros((101,4))
plate_length_effect[:,0] = np.linspace(2000e3,12000e3,101)
vc = 1.0 * vel_converter; 
plate_length_effect[:,1] = compute_vt(vc,h,visc_asthen,visc_lith,H,plate_length_effect[:,0],Rmin,D,dip,oceanic_buoy)
vc = 5.0 * vel_converter;
plate_length_effect[:,2] = compute_vt(vc,h,visc_asthen,visc_lith,H,plate_length_effect[:,0],Rmin,D,dip,oceanic_buoy)
vc = 10.0 * vel_converter;
plate_length_effect[:,3] = compute_vt(vc,h,visc_asthen,visc_lith,H,plate_length_effect[:,0],Rmin,D,dip,oceanic_buoy)


### effect of slab dip
slab_dip_effect = np.zeros((101,4))
slab_dip_effect[:,0] = np.linspace(30,100,101)
vc = 1.0 * vel_converter; 
slab_dip_effect[:,1] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,D,slab_dip_effect[:,0],oceanic_buoy)
vc = 5.0 * vel_converter; 
slab_dip_effect[:,2] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,D,slab_dip_effect[:,0],oceanic_buoy)
vc = 10.0 * vel_converter; 
slab_dip_effect[:,3] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,D,slab_dip_effect[:,0],oceanic_buoy)

### effect of slab dip
slab_age_effect = np.zeros((101,4))
slab_age_effect[:,0] = np.linspace(0,120,101)
H = np.zeros((101)); oceanic_buoy = np.zeros((101)); 
for i in range(0,101):
	H[i] = compute_plate_isotherm(slab_age_effect[i,0],1300,1e-6,1300) * 1e3; # [m]
	oceanic_buoy[i] = compute_plate_buoyancy(slab_age_effect[i,0],1300,1e-6,3e-5,3300); # [kg/m2]
vc = 1.0 * vel_converter; 
slab_age_effect[:,1] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,D,dip,oceanic_buoy)
vc = 5.0 * vel_converter; 
slab_age_effect[:,2] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,D,dip,oceanic_buoy)
vc = 10.0 * vel_converter; 
slab_age_effect[:,3] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,D,dip,oceanic_buoy)

oceanic_buoy = compute_plate_buoyancy(60.0,1300,1e-6,3e-5,3300); 
H = compute_plate_isotherm(60.0,1300,1e-6,1300) * 1e3;
### effect of asthen viscosity / thickness
asthen_effect = np.zeros((101,4))
asthen_effect[:,0] = np.linspace(2e20,20e20,101)
vc = 1.0 * vel_converter; 
asthen_effect[:,1] = compute_vt(vc,h,asthen_effect[:,0],visc_lith,H,Lsp,Rmin,D,dip,oceanic_buoy)
vc = 5.0 * vel_converter; 
asthen_effect[:,2] = compute_vt(vc,h,asthen_effect[:,0],visc_lith,H,Lsp,Rmin,D,dip,oceanic_buoy)
vc = 10.0 * vel_converter; 
asthen_effect[:,3] = compute_vt(vc,h,asthen_effect[:,0],visc_lith,H,Lsp,Rmin,D,dip,oceanic_buoy)

### effect of slab depth
depth_effect = np.zeros((101,4))
depth_effect[:,0] = np.linspace(100e3,1000e3,101)
vc = 1.0 * vel_converter; 
depth_effect[:,1] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,depth_effect[:,0],dip,oceanic_buoy)
vc = 5.0 * vel_converter; 
depth_effect[:,2] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,depth_effect[:,0],dip,oceanic_buoy)
vc = 10.0 * vel_converter; 
depth_effect[:,3] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,depth_effect[:,0],dip,oceanic_buoy)

### effect of radius of curvature
radius_effect = np.zeros((101,4))
radius_effect[:,0] = np.linspace(50e3,700e3,101)
vc = 1.0 * vel_converter; 
radius_effect[:,1] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,radius_effect[:,0],D,dip,oceanic_buoy)
vc = 5.0 * vel_converter; 
radius_effect[:,2] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,radius_effect[:,0],D,dip,oceanic_buoy)
vc = 10.0 * vel_converter; 
radius_effect[:,3] = compute_vt(vc,h,visc_asthen,visc_lith,H,Lsp,radius_effect[:,0],D,dip,oceanic_buoy)

# #
fig, axs = plt.subplots(3, 2, tight_layout=True)

axs[0,0].set_xlabel("$\mathregular{L_{SP}}$  [km]")
axs[0,0].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[0,0].set_ylim(-20,20)
axs[0,0].tick_params(axis='both', which='major', labelsize=10)
axs[0,0].tick_params(axis='both', which='minor', labelsize=8)
axs[0,0].plot(plate_length_effect[:,0]/1e3,0.0*plate_length_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[0,0].plot(plate_length_effect[:,0]/1e3,(1.0-plate_length_effect[:,1]), linewidth=2.5, color='orange')
axs[0,0].plot(plate_length_effect[:,0]/1e3,(5.0-plate_length_effect[:,2]), linewidth=2.5, color='red')
axs[0,0].plot(plate_length_effect[:,0]/1e3,(10.0-plate_length_effect[:,3]), linewidth=2.5, color='brown')
axs[0,0].text(0.6,0.125,'$\mathregular{V_{C}}$ = 1 cm/yr',size=9.5, color="orange",transform = axs[0,0].transAxes)
axs[0,0].text(0.6,0.25,'$\mathregular{V_{C}}$ = 5 cm/yr',size=9.5, color="red",transform = axs[0,0].transAxes)
axs[0,0].text(0.6,0.375,'$\mathregular{V_{C}}$ = 10 cm/yr',size=9.5, color="brown",transform = axs[0,0].transAxes)

axs[0,1].set_xlabel("$\mathregular{dip}$  [deg]")
axs[0,1].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[0,1].set_ylim(-20,20)
axs[0,1].tick_params(axis='both', which='major', labelsize=10)
axs[0,1].tick_params(axis='both', which='minor', labelsize=8)
axs[0,1].plot(slab_dip_effect[:,0],0.0*slab_dip_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[0,1].plot(slab_dip_effect[:,0],(1.0-slab_dip_effect[:,1]), linewidth=2.5, color='orange')
axs[0,1].plot(slab_dip_effect[:,0],(5.0-slab_dip_effect[:,2]), linewidth=2.5, color='red')
axs[0,1].plot(slab_dip_effect[:,0],(10.0-slab_dip_effect[:,3]), linewidth=2.5, color='brown')

axs[1,0].set_xlabel("$\mathregular{age_{SP}}}$  [Ma]")
axs[1,0].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[1,0].set_ylim(-20,20)
axs[1,0].tick_params(axis='both', which='major', labelsize=10)
axs[1,0].tick_params(axis='both', which='minor', labelsize=8)
axs[1,0].plot(slab_age_effect[:,0],0.0*slab_age_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[1,0].plot(slab_age_effect[:,0],(1.0-slab_age_effect[:,1]), linewidth=2.5, color='orange')
axs[1,0].plot(slab_age_effect[:,0],(5.0-slab_age_effect[:,2]), linewidth=2.5, color='red')
axs[1,0].plot(slab_age_effect[:,0],(10.0-slab_age_effect[:,3]), linewidth=2.5, color='brown')

#axs[1,1].set_xlabel("$\mathregular{visc_{A}}$  [Pa.s]")
axs[1,1].set_xlabel("$\mathregular{\eta_{A}}$  [Pa.s]")
axs[1,1].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[1,1].set_ylim(-20,20)
axs[1,1].tick_params(axis='both', which='major', labelsize=10)
axs[1,1].tick_params(axis='both', which='minor', labelsize=8)
axs[1,1].plot(asthen_effect[:,0],0.0*asthen_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[1,1].plot(asthen_effect[:,0],(1.0-asthen_effect[:,1]), linewidth=2.5, color='orange')
axs[1,1].plot(asthen_effect[:,0],(5.0-asthen_effect[:,2]), linewidth=2.5, color='red')
axs[1,1].plot(asthen_effect[:,0],(10.0-asthen_effect[:,3]), linewidth=2.5, color='brown')

axs[2,0].set_xlabel("$\mathregular{D_{SP}}$  [km]")
axs[2,0].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[2,0].set_ylim(-20,20)
axs[2,0].tick_params(axis='both', which='major', labelsize=10)
axs[2,0].tick_params(axis='both', which='minor', labelsize=8)
axs[2,0].plot(depth_effect[:,0]/1e3,0.0*depth_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[2,0].plot(depth_effect[:,0]/1e3,(1.0-depth_effect[:,1]), linewidth=2.5, color='orange')
axs[2,0].plot(depth_effect[:,0]/1e3,(5.0-depth_effect[:,2]), linewidth=2.5, color='red')
axs[2,0].plot(depth_effect[:,0]/1e3,(10.0-depth_effect[:,3]), linewidth=2.5, color='brown')

axs[2,1].set_xlabel("$\mathregular{R_{min}}$  [km]")
axs[2,1].set_ylabel("$\mathregular{V_{T}}$  [cm/yr]")
axs[2,1].set_ylim(-20,20)
axs[2,1].set_xlim(50,700)
axs[2,1].tick_params(axis='both', which='major', labelsize=10)
axs[2,1].tick_params(axis='both', which='minor', labelsize=8)
axs[2,1].plot(radius_effect[:,0],0.0*radius_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[2,1].plot(radius_effect[:,0]/1e3,(1.0-radius_effect[:,1]), linewidth=2.5, color='orange')
axs[2,1].plot(radius_effect[:,0]/1e3,(5.0-radius_effect[:,2]), linewidth=2.5, color='red')
axs[2,1].plot(radius_effect[:,0]/1e3,(10.0-radius_effect[:,3]), linewidth=2.5, color='brown')

plt.savefig('test_using-function.png', bbox_inches='tight')



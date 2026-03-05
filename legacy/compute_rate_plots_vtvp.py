#!/usr/bin/python
import sys

sys.stderr.write(
	"compute_rate_plots_vtvp.py is legacy and currently disabled.\n"
	"Reason: it depends on compute_vsp(), which is not present in functions.py.\n"
	"Use compute_rates_misfit.py / compute_rates_single.py for active workflows.\n"
)
sys.exit(1)

import scipy, math, os, numpy as np
import subprocess, pygplates
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FuncFormatter
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from functions import compute_plate_buoyancy, compute_plate_isotherm
from functions import compute_vsp

vt_ref='sa'    # hs3, nnr, sa
depth_cutoff=0 # 0 = none, 1 = separate at 670
formulation=1 # 1 = regular, 2 = no-bending, 3 = plastic bending, 4 = power-law, 5 = plastic + power-law

vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s 
g = 9.81;

# rheology
visc_lith = 5.0e22; visc_asthen = 1.0e21; h = 200e3;
yield_stress = 4000e6; # [MPa], used for formulation = 3
n = 3.5; trans_strain_rate = 1e-13  # [s^-1], used for formulation = 4

data = np.genfromtxt('data/Lallemand_et_al-2005_G3_dataset.txt') # see spreadsheet for column details
data_vt =np.genfromtxt('data/vts_hs3-nnr-sa.txt') # lat, lon, vtn: hs3, vtn: nnr, vtn: sa
if vt_ref == 'hs3':
   vt_col = 2;
elif vt_ref == 'nnr':
   vt_col = 3;
else:
   vt_col = 4; 

## compute average properties for the curves
Rmin_total = 0; age_total = 0; Lsp_total = 0;
slabL_total = 0; dip_total = 0;
Rmin_n = 0; age_n = 0; Lsp_n = 0;
slabL_n = 0; dip_n = 0;
for i in range(0,data.shape[0]):
    if np.isnan(data[i,13]) == False:
        Rmin_total = Rmin_total + data[i,13];
        Rmin_n = Rmin_n + 1
    if np.isnan(data[i,20]) == False:
        age_total = age_total + data[i,20];
        age_n = age_n + 1
    if np.isnan(data[i,3]) == False:
        Lsp_total = Lsp_total + data[i,3];
        Lsp_n = Lsp_n + 1
    if np.isnan(data[i,8]) == False:
        slabL_total = slabL_total + data[i,8];
        slabL_n = slabL_n + 1
    if np.isnan(data[i,6]) == False:
        dip_total = dip_total + data[i,6];   
        dip_n = dip_n + 1
Rmin_avg = (Rmin_total/Rmin_n) * 1e3;
age_avg = (age_total/age_n);
Lsp_avg = (Lsp_total/Lsp_n) * 1e3;
slabL_avg = (slabL_total/slabL_n) * 1e3;
dip_avg = (dip_total/dip_n);

dip = dip_avg; Lsp = Lsp_avg; age = age_avg; 
Rmin = Rmin_avg; slabL = slabL_avg;
oceanic_buoy = compute_plate_buoyancy(age,1300,1e-6,3e-5,3300); 
H = compute_plate_isotherm(age,1300,1e-6,1300) * 1e3;

## data for plotting
data_plate_length=np.zeros((data.shape[0],3))
data_plate_length_deep=np.zeros((data.shape[0],3))
for i in range(0,data.shape[0]):
    if np.isnan(data_vt[i,vt_col]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,3]) == False and data[i,7] <= 680.0:
        data_plate_length[i,0] = data[i,3];
        data_plate_length[i,1] = data_vt[i,vt_col];
        data_plate_length[i,2] = data[i,9];
    elif np.isnan(data_vt[i,vt_col]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,3]) == False and data[i,7] > 680.0:
        data_plate_length_deep[i,0] = data[i,3];
        data_plate_length_deep[i,1] = data_vt[i,vt_col];
        data_plate_length_deep[i,2] = data[i,9];
data_plate_length = data_plate_length[~np.all(data_plate_length == 0, axis=1)]
data_plate_length_deep = data_plate_length_deep[~np.all(data_plate_length_deep == 0, axis=1)]


data_slab_age=np.zeros((data.shape[0],3))
data_slab_age_deep=np.zeros((data.shape[0],3))
for i in range(0,data.shape[0]):
    if np.isnan(data_vt[i,vt_col]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,20]) == False and data[i,7] <= 680.0:
        data_slab_age[i,0] = data[i,20];
        data_slab_age[i,1] = data_vt[i,vt_col];
        data_slab_age[i,2] = data[i,9];
    elif np.isnan(data_vt[i,vt_col]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,20]) == False and data[i,7] > 680.0:
        data_slab_age_deep[i,0] = data[i,20];
        data_slab_age_deep[i,1] = data_vt[i,vt_col];
        data_slab_age_deep[i,2] = data[i,9];
data_slab_age = data_slab_age[~np.all(data_slab_age == 0, axis=1)]
data_slab_age_deep = data_slab_age_deep[~np.all(data_slab_age_deep == 0, axis=1)]


data_slab_length=np.zeros((data.shape[0],3))
data_slab_length_deep=np.zeros((data.shape[0],3))
for i in range(0,data.shape[0]):
    if np.isnan(data_vt[i,vt_col]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,8]) == False and data[i,7] <= 680.0:
        data_slab_length[i,0] = data[i,8];
        data_slab_length[i,1] = data_vt[i,vt_col];
        data_slab_length[i,2] = data[i,9];
    elif np.isnan(data_vt[i,vt_col]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,8]) == False and data[i,7] > 680.0:
        data_slab_length_deep[i,0] = data[i,8];
        data_slab_length_deep[i,1] = data_vt[i,vt_col];
        data_slab_length_deep[i,2] = data[i,9];
data_slab_length = data_slab_length[~np.all(data_slab_length == 0, axis=1)]
data_slab_length_deep = data_slab_length_deep[~np.all(data_slab_length_deep == 0, axis=1)]

data_Rmin=np.zeros((data.shape[0],3))
data_Rmin_deep=np.zeros((data.shape[0],3))
for i in range(0,data.shape[0]):
    if np.isnan(data_vt[i,vt_col]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,13]) == False and data[i,7] <= 680.0:
        data_Rmin[i,0] = data[i,13];
        data_Rmin[i,1] = data_vt[i,vt_col]; 
        data_Rmin[i,2] = data[i,9];
    if np.isnan(data_vt[i,vt_col]) == False and  np.isnan(data[i,9]) == False and np.isnan(data[i,13]) == False and data[i,7] > 680.0:
        data_Rmin_deep[i,0] = data[i,13];
        data_Rmin_deep[i,1] = data_vt[i,vt_col]; 
        data_Rmin_deep[i,2] = data[i,9];
data_Rmin = data_Rmin[~np.all(data_Rmin == 0, axis=1)]
data_Rmin_deep = data_Rmin_deep[~np.all(data_Rmin_deep == 0, axis=1)]

### effect of plate length (Lsp)
plate_length_effect = np.zeros((101,4))
plate_length_effect[:,0] = np.linspace(300e3,12000e3,101)
vc = 1.0 * vel_converter; 
if formulation == 4 or formulation == 5:
    for i in range(0,len(plate_length_effect)):
        plate_length_effect[i,1] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,plate_length_effect[i,0],Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
else:
    plate_length_effect[:,1] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,plate_length_effect[:,0],Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
vc = 5.0 * vel_converter;
if formulation == 4 or formulation == 5:
    for i in range(0,len(plate_length_effect)):
        plate_length_effect[i,2] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,plate_length_effect[i,0],Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
else:
    plate_length_effect[:,2] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,plate_length_effect[:,0],Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
vc = 10.0 * vel_converter;
if formulation == 4 or formulation == 5:
    for i in range(0,len(plate_length_effect)):
        plate_length_effect[i,3] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,plate_length_effect[i,0],Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
else:
    plate_length_effect[:,3] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,plate_length_effect[:,0],Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)

### effect of slab length 
slab_length_effect = np.zeros((101,4))
slab_length_effect[:,0] = np.linspace(200e3,1500e3,101)
vc = 1.0 * vel_converter; 
if formulation == 4 or formulation == 5:
    for i in range(0,len(slab_length_effect)):
        slab_length_effect[i,1] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slab_length_effect[i,0],dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
else:
    slab_length_effect[:,1] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slab_length_effect[:,0],dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
vc = 5.0 * vel_converter; 
if formulation == 4 or formulation == 5:
    for i in range(0,len(slab_length_effect)):
        slab_length_effect[i,2] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slab_length_effect[i,0],dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
else:
    slab_length_effect[:,2] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slab_length_effect[:,0],dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
vc = 10.0 * vel_converter; 
if formulation == 4 or formulation == 5:
    for i in range(0,len(slab_length_effect)):
        slab_length_effect[i,3] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slab_length_effect[i,0],dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
else:
    slab_length_effect[:,3] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slab_length_effect[:,0],dip,oceanic_buoy,yield_stress,trans_strain_rate,n)

### effect of slab age 
slab_age_effect = np.zeros((101,4))
slab_age_effect[:,0] = np.linspace(5,150,101)
H = np.zeros((101)); oceanic_buoy = np.zeros((101)); 
for i in range(0,101):
	H[i] = compute_plate_isotherm(slab_age_effect[i,0],1300,1e-6,1300) * 1e3; # [m]
	oceanic_buoy[i] = compute_plate_buoyancy(slab_age_effect[i,0],1300,1e-6,3e-5,3300); # [kg/m2]
vc = 1.0 * vel_converter; 
if formulation == 4 or formulation == 5:
    for i in range(0,len(slab_age_effect)):
        slab_age_effect[i,1] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H[i],Lsp,Rmin,slabL,dip,oceanic_buoy[i],yield_stress,trans_strain_rate,n)
else:
    slab_age_effect[:,1] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
vc = 5.0 * vel_converter; 
if formulation == 4 or formulation == 5:
    for i in range(0,len(slab_age_effect)):
        slab_age_effect[i,2] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H[i],Lsp,Rmin,slabL,dip,oceanic_buoy[i],yield_stress,trans_strain_rate,n)
else:
    slab_age_effect[:,2] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
vc = 10.0 * vel_converter; 
if formulation == 4 or formulation == 5:
    for i in range(0,len(slab_age_effect)):
        slab_age_effect[i,3] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H[i],Lsp,Rmin,slabL,dip,oceanic_buoy[i],yield_stress,trans_strain_rate,n)
else:
    slab_age_effect[:,3] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)


oceanic_buoy = compute_plate_buoyancy(60.0,1300,1e-6,3e-5,3300); 
H = compute_plate_isotherm(60.0,1300,1e-6,1300) * 1e3;
### effect of radius of curvature
radius_effect = np.zeros((101,4))
radius_effect[:,0] = np.linspace(50e3,700e3,101)
vc = 1.0 * vel_converter; 
if formulation == 4 or formulation == 5:
    for i in range(0,len(radius_effect)):
        radius_effect[i,1] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,radius_effect[i,0],slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
else:
    radius_effect[:,1] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,radius_effect[:,0],slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
vc = 5.0 * vel_converter; 
if formulation == 4 or formulation == 5:
    for i in range(0,len(radius_effect)):
        radius_effect[i,2] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,radius_effect[i,0],slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
else:
    radius_effect[:,2] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,radius_effect[:,0],slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
vc = 10.0 * vel_converter; 
if formulation == 4 or formulation == 5:
    for i in range(0,len(radius_effect)):
        radius_effect[i,3] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,radius_effect[i,0],slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
else:
    radius_effect[:,3] = compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,radius_effect[:,0],slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n)
# #
fig, axs = plt.subplots(2, 2, tight_layout=True)


ymin=-2
ymax=8
color_bound1=3;
color_bound2=7;

axs[0,0].set_xlabel("$\mathregular{L_{SP}}$  [km]")
axs[0,0].set_ylabel("$\mathregular{V_{T}}$ \ $\mathregular{V_{P}}$")
axs[0,0].set_ylim(ymin,ymax)
axs[0,0].set_xlim(200,12000)
axs[0,0].tick_params(axis='both', which='major', labelsize=10)
axs[0,0].tick_params(axis='both', which='minor', labelsize=8)
axs[0,0].plot(plate_length_effect[:,0]/1e3,0.0*plate_length_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[0,0].plot(plate_length_effect[:,0]/1e3,(1.0-plate_length_effect[:,1])/plate_length_effect[:,1], linewidth=1.5, color='orange')
axs[0,0].plot(plate_length_effect[:,0]/1e3,(5.0-plate_length_effect[:,2])/plate_length_effect[:,2], linewidth=1.5, color='red')
axs[0,0].plot(plate_length_effect[:,0]/1e3,(10.0-plate_length_effect[:,3])/plate_length_effect[:,3], linewidth=1.5, color='brown')
axs[0,0].text(0.0,1.05,'$\mathregular{V_{C}}$ = 1 cm/yr',size=9, color="orange",transform = axs[0,0].transAxes)
axs[0,0].text(0.325,1.05,'$\mathregular{V_{C}}$ = 5 cm/yr',size=9, color="red",transform = axs[0,0].transAxes)
axs[0,0].text(0.65,1.05,'$\mathregular{V_{C}}$ = 10 cm/yr',size=9, color="brown",transform = axs[0,0].transAxes)
if depth_cutoff == 0:
    color_col=data_plate_length[:,2]/10.
    colors=np.where(color_col<color_bound1,'orange',np.where(color_col<color_bound2,'red','brown'))
    cset1=axs[0,0].scatter(data_plate_length[:,0],data_plate_length[:,1]/(data_plate_length[:,2]-data_plate_length[:,1]),vmin=0,vmax=10,s=8,c=colors,edgecolors='none')
    color_col=data_plate_length_deep[:,2]/10.
    colors=np.where(color_col<color_bound1,'orange',np.where(color_col<color_bound2,'red','brown'))
    cset2=axs[0,0].scatter(data_plate_length_deep[:,0],data_plate_length_deep[:,1]/(data_plate_length_deep[:,2]-data_plate_length_deep[:,1]),vmin=0,vmax=10,s=8,c=colors,edgecolors='none')
else:
    axs[0,0].scatter(data_plate_length_deep[:,0],data_plate_length_deep[:,1]/10.,s=6,c='gray',edgecolors='gray')
    axs[0,0].scatter(data_plate_length[:,0],data_plate_length[:,1]/10.,s=6,c='black',edgecolors='black')
if vt_ref == 'hs3':
    axs[0,0].text(0.07,0.875,'HS3',size=14, color="black",transform = axs[0,0].transAxes)
elif vt_ref == 'nnr':
    axs[0,0].text(0.07,0.875,'NNR',size=14, color="black",transform = axs[0,0].transAxes)
else:
    axs[0,0].text(0.07,0.875,'S-A',size=14, color="black",backgroundcolor="lightgrey",transform = axs[0,0].transAxes)
# axs11 = inset_axes(axs[0,0], width="40%", height="4%")
# cbar=plt.colorbar(cset1,cax=axs11,orientation="horizontal",ticks=[0,5,10])
# cbar.ax.tick_params(labelsize=8)

axs[0,1].set_xlabel("$\mathregular{L_{slab}}$  [deg]")
axs[0,1].set_ylabel("$\mathregular{V_{T}}$ \ $\mathregular{V_{P}}$")
axs[0,1].set_ylim(ymin,ymax)
axs[0,1].set_xlim(200,1500)
axs[0,1].tick_params(axis='both', which='major', labelsize=10)
axs[0,1].tick_params(axis='both', which='minor', labelsize=8)
axs[0,1].plot(slab_length_effect[:,0]/1e3,0.0*slab_length_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[0,1].plot(slab_length_effect[:,0]/1e3,(1.0-slab_length_effect[:,1])/slab_length_effect[:,1], linewidth=1.5, color='orange')
axs[0,1].plot(slab_length_effect[:,0]/1e3,(5.0-slab_length_effect[:,2])/slab_length_effect[:,2], linewidth=1.5, color='red')
axs[0,1].plot(slab_length_effect[:,0]/1e3,(10.0-slab_length_effect[:,3])/slab_length_effect[:,3], linewidth=1.5, color='brown')
if depth_cutoff == 0:
    color_col=data_slab_length[:,2]/10.
    colors=np.where(color_col<color_bound1,'orange',np.where(color_col<color_bound2,'red','brown'))
    cset1=axs[0,1].scatter(data_slab_length[:,0],data_slab_length[:,1]/(data_slab_length[:,2]-data_slab_length[:,1]),vmin=0,vmax=10,s=8,c=colors,cmap=plt.cm.OrRd,edgecolors='none')
    color_col=data_slab_length_deep[:,2]/10.
    colors=np.where(color_col<color_bound1,'orange',np.where(color_col<color_bound2,'red','brown'))
    cset2=axs[0,1].scatter(data_slab_length_deep[:,0],data_slab_length_deep[:,1]/(10.*(data_slab_length_deep[:,2]-data_slab_length_deep[:,1])),vmin=0,vmax=10,s=8,c=colors,edgecolors='none')
else:
    axs[0,1].scatter(data_slab_length_deep[:,0],data_slab_length_deep[:,1]/10.,s=6,c='gray',edgecolors='gray')
    axs[0,1].scatter(data_slab_length[:,0],data_slab_length[:,1]/10.,s=6,c='black',edgecolors='black')


axs[1,0].set_xlabel("$\mathregular{age_{SP}}}$  [Ma]")
axs[1,0].set_ylabel("$\mathregular{V_{T}}$ \ $\mathregular{V_{P}}$")
axs[1,0].set_ylim(ymin,ymax)
axs[1,0].set_xlim(5,150)
axs[1,0].tick_params(axis='both', which='major', labelsize=10)
axs[1,0].tick_params(axis='both', which='minor', labelsize=8)
axs[1,0].plot(slab_age_effect[:,0],0.0*slab_age_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[1,0].plot(slab_age_effect[:,0],(1.0-slab_age_effect[:,1])/slab_age_effect[:,1], linewidth=1.5, color='orange')
axs[1,0].plot(slab_age_effect[:,0],(5.0-slab_age_effect[:,2])/slab_age_effect[:,2], linewidth=1.5, color='red')
axs[1,0].plot(slab_age_effect[:,0],(10.0-slab_age_effect[:,3])/slab_age_effect[:,3], linewidth=1.5, color='brown')
if depth_cutoff == 0:
    color_col=data_slab_age[:,2]/10.
    colors=np.where(color_col<color_bound1,'orange',np.where(color_col<color_bound2,'red','brown'))
    cset1=axs[1,0].scatter(data_slab_age[:,0],data_slab_age[:,1]/(data_slab_age[:,2]-data_slab_age[:,1]),vmin=0,vmax=10,s=8,c=colors,edgecolors='none')
    color_col=data_slab_age_deep[:,2]/10.
    colors=np.where(color_col<color_bound1,'orange',np.where(color_col<color_bound2,'red','brown'))
    cset2=axs[1,0].scatter(data_slab_age_deep[:,0],data_slab_age_deep[:,1]/(data_slab_age_deep[:,2]-data_slab_age_deep[:,1]),vmin=0,vmax=10,s=8,c=colors,edgecolors='none')
else:
    axs[1,0].scatter(data_slab_age_deep[:,0],data_slab_age_deep[:,1]/10.,s=6,c='gray',edgecolors='gray')
    axs[1,0].scatter(data_slab_age[:,0],data_slab_age[:,1]/10.,s=6,c='black',edgecolors='black')

axs[1,1].set_xlabel("$\mathregular{R_{min}}$  [km]")
axs[1,1].set_ylabel("$\mathregular{V_{T}}$ \ $\mathregular{V_{P}}$")
axs[1,1].set_ylim(ymin,ymax)
axs[1,1].set_xlim(50,600)
axs[1,1].tick_params(axis='both', which='major', labelsize=10)
axs[1,1].tick_params(axis='both', which='minor', labelsize=8)
axs[1,1].plot(radius_effect[:,0]/1e3,0.0*radius_effect[:,0], linewidth=1, color='grey',linestyle='--')
axs[1,1].plot(radius_effect[:,0]/1e3,(1.0-radius_effect[:,1])/radius_effect[:,1], linewidth=1.5, color='orange')
axs[1,1].plot(radius_effect[:,0]/1e3,(5.0-radius_effect[:,2])/radius_effect[:,2], linewidth=1.5, color='red')
axs[1,1].plot(radius_effect[:,0]/1e3,(10.0-radius_effect[:,3])/radius_effect[:,3], linewidth=1.5, color='brown')
if depth_cutoff == 0:
    color_col=data_Rmin[:,2]/10.
    colors=np.where(color_col<color_bound1,'orange',np.where(color_col<color_bound2,'red','brown'))
    cset1=axs[1,1].scatter(data_Rmin[:,0],data_Rmin[:,1]/(data_Rmin[:,2]-data_Rmin[:,1]),vmin=0,vmax=10,s=8,c=colors,edgecolors='none')
    color_col=data_Rmin_deep[:,2]/10.
    colors=np.where(color_col<color_bound1,'orange',np.where(color_col<color_bound2,'red','brown'))
    cset2=axs[1,1].scatter(data_Rmin_deep[:,0],data_Rmin_deep[:,1]/(data_Rmin_deep[:,2]-data_Rmin_deep[:,1]),vmin=0,vmax=10,s=8,c=colors,edgecolors='none')
else:
    axs[1,1].scatter(data_Rmin_deep[:,0],data_Rmin_deep[:,1]/10.,s=6,c='gray',edgecolors='gray')
    axs[1,1].scatter(data_Rmin[:,0],data_Rmin[:,1]/10.,s=6,c='black',edgecolors='black')

if depth_cutoff == 0:
    depths_string='no_depth_cutoff'
else:
    depths_string='depths_separated'

if formulation == 1:
    plot_name=''.join(['plots/ratioVtVp_',str(vt_ref),'model_',depths_string,'_LithVisc',str(visc_lith),'_AsthVisc',str(visc_asthen),'.png'])
elif formulation == 2:
    plot_name=''.join(['plots/ratioVtVp_',str(vt_ref),'model_',depths_string,'_LithVisc',str(visc_lith),'_AsthVisc',str(visc_asthen),'.no_bending.png'])
elif formulation == 3:
    plot_name=''.join(['plots/ratioVtVp_',str(vt_ref),'model_',depths_string,'_LithVisc',str(visc_lith),'_AsthVisc',str(visc_asthen),'.plastic_bending.png'])
elif formulation == 4:
    plot_name=''.join(['plots/ratioVtVp_',str(vt_ref),'model_',depths_string,'_LithVisc',str(visc_lith),'_AsthVisc',str(visc_asthen),'_Trans',str(trans_strain_rate),'.power-law.png'])
elif formulation == 5:
    plot_name=''.join(['plots/ratioVtVp_',str(vt_ref),'model_',depths_string,'_YieldStress',str(yield_stress/1e9),'GPa_AsthVisc',str(visc_asthen),'_Trans',str(trans_strain_rate),'.plastic-bending_and_power-law.png'])

plt.savefig(plot_name, bbox_inches='tight',dpi=400)


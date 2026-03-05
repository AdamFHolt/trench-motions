#!/bin/python
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import scipy
import scipy.integrate

h_dim = 80e3
a_dim = (660e3 - h_dim);
mu = 5.0e20;
w_dim = 0.5 * h_dim * (1.0/np.tan(np.deg2rad(60.0)))
C_fs = mu * ((3.0 * h_dim)/a_dim**3) 
C_ns = mu * ((12.0 * h_dim)/a_dim**3)
vp = 10.0 * (0.01/(365*24*60*60));

width_dim=2000e3
w = w_dim / (0.5*width_dim)
width_dim2=4000e3
w2 = w_dim / (0.5*width_dim2)

Lsp_range=np.arange(500,13000,100)
forces_fs=np.zeros((len(Lsp_range)))
forces_ns=np.zeros((len(Lsp_range)))
forces_fs2=np.zeros((len(Lsp_range)))
forces_ns2=np.zeros((len(Lsp_range)))

for i in range(len(Lsp_range)):

	Lsp_dim=Lsp_range[i] * 1.e3;
	y_dim = np.linspace(0,Lsp_dim,201)

	Lsp=Lsp_dim/(0.5*width_dim)
	y=y_dim/(0.5*width_dim)

	Lsp2=Lsp_dim/(0.5*width_dim2)
	y2=y_dim/(0.5*width_dim2)

	# width 1
	P_fs = (2.0/np.pi)*(vp * C_fs * (0.5*width_dim)) * ( (y + 2*Lsp)*np.arctan(1.0/(y+2*Lsp)) - y*np.arctan(1.0/y) + 0.5*np.log((1.0+(y + 2*Lsp)**2)/(1.0+y**2)) - 0.5*((np.arctan(2*Lsp))/(1.0-w))*(np.sqrt(1+y**2)-y)  )
	P_ns = (2.0/np.pi)*(vp * C_ns * (0.5*width_dim)) * ( (y + 2*Lsp)*np.arctan(1.0/(y+2*Lsp)) - y*np.arctan(1.0/y) + 0.5*np.log((1.0+(y + 2*Lsp)**2)/(1.0+y**2)) - 0.5*((np.arctan(2*Lsp))/(1.0-w))*(np.sqrt(1+y**2)-y)  )
	# width 2
	P_fs2 = (2.0/np.pi)*(vp * C_fs * (0.5*width_dim2)) * ( (y2 + 2*Lsp2)*np.arctan(1.0/(y2+2*Lsp2)) - y2*np.arctan(1.0/y2) + 0.5*np.log((1.0+(y2 + 2*Lsp2)**2)/(1.0+y2**2)) - 0.5*((np.arctan(2*Lsp2))/(1.0-w2))*(np.sqrt(1+y2**2)-y2)  )
	P_ns2 = (2.0/np.pi)*(vp * C_ns * (0.5*width_dim2)) * ( (y2 + 2*Lsp2)*np.arctan(1.0/(y2+2*Lsp2)) - y2*np.arctan(1.0/y2) + 0.5*np.log((1.0+(y2 + 2*Lsp2)**2)/(1.0+y2**2)) - 0.5*((np.arctan(2*Lsp2))/(1.0-w2))*(np.sqrt(1+y2**2)-y2)  )

	dPdy_fs = np.zeros((len(P_fs)));
	dPdy_ns = np.zeros((len(P_ns)));
	dPdy_fs2 = np.zeros((len(P_fs2)));
	dPdy_ns2 = np.zeros((len(P_ns2)));

	for j in range(0,len(P_fs)):
		if j == 0:
			dPdy_fs[j] = (P_fs[j+1] - P_fs[j])/(y_dim[j+1] - y_dim[j]);  # [Pa/m]
			dPdy_ns[j] = (P_ns[j+1] - P_ns[j])/(y_dim[j+1] - y_dim[j]); 
			dPdy_fs2[j] = (P_fs2[j+1] - P_fs2[j])/(y_dim[j+1] - y_dim[j]);  
			dPdy_ns2[j] = (P_ns2[j+1] - P_ns2[j])/(y_dim[j+1] - y_dim[j]); 
	  	elif j == len(P_fs)-1:
			dPdy_fs[j] = (P_fs[j]   - P_fs[j-1])/(y_dim[j] - y_dim[j-1]);
			dPdy_ns[j] = (P_ns[j]   - P_ns[j-1])/(y_dim[j] - y_dim[j-1]);
			dPdy_fs2[j] = (P_fs2[j]   - P_fs2[j-1])/(y_dim[j] - y_dim[j-1]);
			dPdy_ns2[j] = (P_ns2[j]   - P_ns2[j-1])/(y_dim[j] - y_dim[j-1]);
		else:
			dPdy_fs[j] = (P_fs[j+1] - P_fs[j-1])/(y_dim[j+1] - y_dim[j-1]);
			dPdy_ns[j] = (P_ns[j+1] - P_ns[j-1])/(y_dim[j+1] - y_dim[j-1]);\
			dPdy_fs2[j] = (P_fs2[j+1] - P_fs2[j-1])/(y_dim[j+1] - y_dim[j-1]);
			dPdy_ns2[j] = (P_ns2[j+1] - P_ns2[j-1])/(y_dim[j+1] - y_dim[j-1]);

	shear_stress_fs = dPdy_fs * (h_dim - a_dim) # [Pa]
	shear_stress_ns = (dPdy_ns * (h_dim - 0.5*a_dim)) + mu*(vp/a_dim) 
	shear_stress_fs2 = dPdy_fs2 * (h_dim - a_dim) # [Pa]
	shear_stress_ns2 = (dPdy_ns2 * (h_dim - 0.5*a_dim)) + mu*(vp/a_dim) 
	
	forces_fs[i] = scipy.integrate.simps(y=shear_stress_fs, x=y_dim, even='avg') # [N/m]
	forces_ns[i] = scipy.integrate.simps(y=shear_stress_ns, x=y_dim, even='avg')
	forces_fs2[i] = scipy.integrate.simps(y=shear_stress_fs2, x=y_dim, even='avg')
	forces_ns2[i] = scipy.integrate.simps(y=shear_stress_ns2, x=y_dim, even='avg')
 
# plot
fig, axs = plt.subplots(2, 1, tight_layout=False)

axs[0].plot(Lsp_range, forces_fs, linewidth=1.5, color='blue',linestyle='-')
axs[0].plot(Lsp_range, forces_fs2, linewidth=1.5, color='blue',linestyle='--')
axs[0].set_ylabel("$F_{SP}$   [N/m]")


axs[1].plot(Lsp_range, forces_ns, linewidth=1.5, color='red',linestyle='-')
axs[1].plot(Lsp_range, forces_ns2, linewidth=1.5, color='red',linestyle='--')
axs[1].set_ylabel("$F_{SP}$   [N/m]")
axs[1].set_xlabel("$L_{SP}$   [km]")


plt.show()






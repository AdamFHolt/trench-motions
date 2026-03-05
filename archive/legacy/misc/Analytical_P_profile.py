#!/bin/python
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm

BC=1  # 1 = free slip, 0 = no slip

h_dim = 80e3
a_dim = (660e3 - h_dim);
mu = 5.0e20;
w_dim = 0.5 * h_dim * (1.0/np.tan(np.deg2rad(60.0)))
if BC == 1:
	C = mu * ((3.0 * h_dim)/a_dim**3) 
else:
	C = mu * ((12.0 * h_dim)/a_dim**3)

# 10 cm/yr, 10000 km long, 5100 km wide
Lsp_dim=10000e3
width_dim=5000e3
y_dim = np.linspace(0,11000e3,201)
w = w_dim / (0.5*width_dim)
Lsp=Lsp_dim/(0.5*width_dim)
y=y_dim/(0.5*width_dim)

vp = 10.0 * (0.01/(365*24*60*60));
P2 = (2.0/np.pi)*(vp * C * (0.5*width_dim)) * ( (y + 2*Lsp)*np.arctan(1.0/(y+2*Lsp)) - y*np.arctan(1.0/y) + 0.5*np.log((1.0+(y + 2*Lsp)**2)/(1.0+y**2)) - 0.5*((np.arctan(2*Lsp))/(1.0-w))*(np.sqrt(1+y**2)-y)  )
dP2dy = np.zeros((len(P2)));
for i in range(0,len(P2)):
	if i == 0:
		dP2dy[i] = (P2[i+1] - P2[i])/(y_dim[i+1] - y_dim[i]);  # [Pa/m]
  	elif i == len(P2)-1:
		dP2dy[i] = (P2[i] - P2[i-1])/(y_dim[i] - y_dim[i-1]);
	else:
		dP2dy[i] = (P2[i+1] - P2[i-1])/(y_dim[i+1] - y_dim[i-1]);
if BC == 1:
	shear_stress2 = dP2dy * (h_dim - a_dim)
else:
	shear_stress2 = (dP2dy * (h_dim - 0.5*a_dim)) + mu*(vp/a_dim) 

# 10 cm/yr, 5800 km long, 4300 km wide
Lsp_dim=4500e3
width_dim = 5000e3
y_dim2 = np.linspace(0,4500e3,201)
Lsp = Lsp_dim/(0.5*width_dim)
w = w_dim / (0.5*width_dim)
y2=y_dim2/(0.5*width_dim)
vp = 10.0 * (0.01/(365*24*60*60));
P4 = (2.0/np.pi)*(vp * C * (0.5*width_dim)) * ( (y2 + 2*Lsp)*np.arctan(1.0/(y2+2*Lsp)) - y2*np.arctan(1.0/y2) + 0.5*np.log((1.0+(y2 + 2*Lsp)**2)/(1.0+y2**2)) - 0.5*((np.arctan(2*Lsp))/(1.0-w))*(np.sqrt(1+y2**2)-y2)  )
dP4dy = np.zeros((len(P4)));
for i in range(0,len(P4)):
	if i == 0:
		dP4dy[i] = (P4[i+1] - P4[i])/(y_dim2[i+1] - y_dim[i]);  
  	elif i == len(P4)-1:
		dP4dy[i] = (P4[i] - P4[i-1])/(y_dim2[i] - y_dim2[i-1]);
	else:
		dP4dy[i] = (P4[i+1] - P4[i-1])/(y_dim2[i+1] - y_dim2[i-1]);
if BC == 1:
	shear_stress4 = dP4dy * (h_dim - a_dim)
else:
	shear_stress4 = (dP4dy * (h_dim - 0.5*a_dim)) + mu*(vp/a_dim)

# plot
fig, axs = plt.subplots(2, 1, tight_layout=False)

axs[0].plot(y_dim/1e3, P2/1e6, linewidth=1.5, color='blue',linestyle='-')
axs[0].plot(y_dim2/1e3, P4/1e6, linewidth=1.5, color='red',linestyle='-')
axs[0].set_ylabel("P  [MPa]")

axs[1].plot(y_dim/1e3, shear_stress2/1e6, linewidth=1.5, color='blue',linestyle='-')
axs[1].plot(y_dim2/1e3, shear_stress4/1e6, linewidth=1.5, color='red',linestyle='-')
axs[1].set_ylabel("shear stress  [MPa]")
axs[1].set_xlabel("y  [km]")

plt.show()






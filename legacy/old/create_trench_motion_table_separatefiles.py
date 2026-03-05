#!/usr/bin/python

## ABANDONED...

import scipy, math, os, numpy as np

coords = np.genfromtxt('data/Lallemand_et_al-2005_G3_dataset.txt')[:,1:3] 
azims =  np.genfromtxt('data/Lallemand_et_al-2005_G3_dataset.txt')[:,4] # azimuths

data_hs3 = np.concatenate(( coords[:,1], coords[:,0] , azims, np.zeros((data.shape[0],4))), axis=1)
data_nnr = np.concatenate(( coords[:,1], coords[:,0] , azims, np.zeros((data.shape[0],4))), axis=1)
data_sa = np.concatenate((  coords[:,1], coords[:,0] , azims, np.zeros((data.shape[0],4))), axis=1)

hs3 = np.genfromtxt('data/vt/tnew.hs3.dat')
nnr = np.genfromtxt('data/vt/tnew.nnr.dat')
sa = np.genfromtxt('data/vt/tnew.sa.dat')

for i in range(0,len(hs3)):
	for n in range(0,len(data)):
		if abs(data[n,1] - hs3[i,1]) <= 0.1 and abs(data[n,0] - hs3[i,0]) <= 0.1:
			data_hs3[n,3] = hs3[i,5];
			data_hs3[n,4] = 
			data_hs3[n,5] =

np.savetxt('data/vts_hs3-nnr-sa.txt',data,fmt='%.4f')




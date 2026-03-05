#!/usr/bin/env python3

import scipy, math, os, numpy as np

data = np.genfromtxt('data/Lallemand_et_al-2005_G3_dataset.txt')[:,1:3] # see spreadsheet for column details
data = np.concatenate(( data , np.zeros((data.shape[0],3))), axis=1)

hs3 = np.genfromtxt('data/vt/tnew.hs3.dat')
nnr = np.genfromtxt('data/vt/tnew.nnr.dat')
sa = np.genfromtxt('data/vt/tnew.sa.dat')

for i in range(0,len(hs3)):
	for n in range(0,len(data)):
		if abs(data[n,0] - hs3[i,1]) <= 0.1 and abs(data[n,1] - hs3[i,0]) <= 0.1:
			data[n,2] = hs3[i,5];
			data[n,3] = nnr[i,5];
			data[n,4] = sa[i,5];


np.savetxt('data/vts_hs3-nnr-sa.txt',data,fmt='%.4f')




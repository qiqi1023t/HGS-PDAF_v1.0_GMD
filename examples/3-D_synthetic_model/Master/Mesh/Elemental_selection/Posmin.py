#!/usr/bin/python
import string,sys,math,struct,array,csv
import numpy
import time
import pickle
import numpy as np

print('Make The Link Between ALLUVSIN GRID AND HGS GRID...')
# ========= READ ALLUVSIM MODEL ==============

filename = './AlluvsimGrid.bin'
f = open(filename, 'rb')

Ya, Xa, Za, Lithoa = pickle.load(f)

f.close()

Ya_corr = np.array([(ya - 600)*-1 for ya in Ya])
Xa_corr = np.array([(xa - 100) for xa in Xa])

Ya = Ya_corr
Xa = Xa_corr

# ========= READ THE 3D HGS MODEL ==============

filename = './hgs_3dmodel.bin'
f = open(filename, 'rb')

Xr, Yr, Zr, Zr_top, E = pickle.load(f)

f.close()

triangle_ele_filename = "./mesh_HGS.ele" 

# Assign Geology

Posmin = {}

for i in range(len(E)):
    HGSelem = E[i]
    dist2 = numpy.square(Xr[i]-Xa) + numpy.square(Yr[i]-Ya) + numpy.square(Zr[i]-Za)
    posmin= numpy.argmin(dist2)
    Posmin[HGSelem] = posmin

Posmin_bin = (Posmin)

#  save output
filename = './Posmin.bin'
f = open(filename, 'wb')
pickle.dump( Posmin_bin, f, 2)
f.close()











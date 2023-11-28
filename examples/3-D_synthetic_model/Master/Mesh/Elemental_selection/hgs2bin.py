#!/usr/bin/python
import string,sys,math,struct,array,csv
import numpy
import time
import pickle

# !!!!!!!!!!
# THIS PYTHON SCRIPT MUST BE EXECUTE ON A 64bit PYTHON CONSOLE IN ORDER TO PROPERLY READ THE MODEL_VSL.txt file

baryelem = open('./baryelem.txt','r')

# ------------ > Read the barycentre elements 

print('Reading and store information about HydroGeosphere 3D model...')

elements = []
xr = []
yr = []
zr = []
zt = []

while 1:
    line = baryelem.readline()
    if line=='':
        break
    ielem,x,y,z,ztop=line.split(';')
    elements.append(float(ielem))
    xr.append(float(x))
    yr.append(float(y))
    zr.append(float(z))
    zt.append(float(ztop))

hgs_3dmodel = ( xr, yr, zr, zt, elements)

#  save output
filename = './hgs_3dmodel.bin'
f = open(filename, 'wb')
pickle.dump( hgs_3dmodel, f, 2)
f.close()

import string,sys,math
import numpy
import time
import os
import numpy as np
from scipy.stats import gmean

plot = False

if plot == True:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

######################
# GENERIC PLOT OPTIONS
######################

linewidth = 2.0
size_scatter = 1
colobs = 'black'
colsim = 'black'

size = 3 # controls default text sizes
labelsize=3  # fontsize of the x and y labels
titlesize=2  

Tstep_times = list(np.arange(365,460,1))

x_coordobs_plain = [15,35,55,75,95,115,135,155,175,195,215,235,255,275]
y_coordobs_plain = [5,50,100,150,200,250,300,350,400,450,495]

z_coordobs = [85]

XX_obs,YY_obs = np.meshgrid(list(x_coordobs_plain),list(y_coordobs_plain))

Grid_obs = list(zip(XX_obs.ravel(), YY_obs.ravel()))

obs_idx = [16,26,46,57,82,86,109,128]
reduced_Grid_obs = [Grid_obs[idx] for idx in obs_idx]
Grid_obs = reduced_Grid_obs

nobs = len(Grid_obs)

for i in range(nobs):
    Heads = []
    outputfile = open('./SIM/H/obs' + str(i+1) + '.dat','r')
    while 1:
        lineread = outputfile.readline()
        if lineread == '':
            break
        Heads.append(float(lineread.split()[0]))
    filename = './HGS/Flowo.observation_well_flow.obs_' + str(i+1) + '.dat'
    hgsfile = open(filename, 'r')
    lineread = hgsfile.readline()
    lineread = hgsfile.readline()
    lineread = hgsfile.readline()
    T = []
    H = []
    while 1:
        lineread = hgsfile.readline()
        if lineread == '':
            break
        t, h, p, s, q, ho, po, qo, x, y, z, zsurf, node = lineread.split()
        T.append(float(t))
        H.append(float(h))
    sigma = 0.05 # m
    outputfile = open('./SIM/H/obs' + str(i+1) + '.dat','w')
    Heads.append(H[-1])
    for val in Heads:
        outputfile.write(str(val) + '\n')



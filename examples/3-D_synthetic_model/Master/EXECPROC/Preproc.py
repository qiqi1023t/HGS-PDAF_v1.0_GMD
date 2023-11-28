import string,sys,math,os, time
from operator import itemgetter, attrgetter
import datetime
from math import *
import numpy as np
import pickle

from utils import *

plot = False

sep = 1

if plot == True:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

# Read depths

Depths = []
Depths_m =  []

currentfile = open('./Mesh/Depths.dat','r')

while 1:
    line = currentfile.readline()
    if line == '':
        break
    depth = float(line)*101.97 # g/cm2 
    depth_m = float(line) #  m 
    Depths.append(depth)
    Depths_m.append(depth_m)

n_layer = len(Depths)

# ========= READ RIVERBED ELEMENTS ==============

Riverbed_elements = []

currentfile = open('./Mesh/Riverbed_elements.dat','r')

while 1:
    line = currentfile.readline()
    if line == '':
        break
    elem = float(line) 
    Riverbed_elements.append(elem)

currentfile.close()

# ========= READ THE 3D HGS MODEL ==============

filename = './Mesh/Elemental_selection/hgs_3dmodel.bin'
f = open(filename, 'rb')

Xr, Yr, Zr, Zt, E = pickle.load(f)

f.close()

# ========= READ ALLUVSIM INDEXES ==============

filename = './Mesh/Elemental_selection/Posmin.bin'
f = open(filename, 'rb')

Posmin = pickle.load(f)

f.close()

# ========= READ RIVERBED ELEMENTS ==============

Riverbed_elements = []

currentfile = open('./Mesh/Riverbed_elements.dat','r')

while 1:
    line = currentfile.readline()
    if line == '':
        break
    elem = float(line) 
    Riverbed_elements.append(elem)

f.close()

# ========= READ BELOW RIVER ELEMENTS ==============

Belowriver_elements = []

currentfile = open('./Mesh/Belowriver_elements.dat','r')

while 1:
    line = currentfile.readline()
    if line == '':
        break
    elem = float(line) 
    Belowriver_elements.append(elem)

f.close()


####################### READ PARAMETERS ##############################
######################################################################
######################################################################

####################
# Aquifer properties
####################

#..K Aquifer
inputfile = './K.dat'

K_list = []

currentfile = open(inputfile,'r')
while 1:
    line = currentfile.readline()[:-1]
    if line == '':
        break
    val = float(line)
    K_list.append(val)

#..S Aquifer
inputfile = './S.dat'

S_list = []

currentfile = open(inputfile,'r')
while 1:
    line = currentfile.readline()[:-1]
    if line == '':
        break
    val = float(line)/100.0
    S_list.append(val)


#..S background
inputfile = './Parameters/Calibrated/param_S.dat'

S_background = []

currentfile = open(inputfile,'r')
while 1:
    line = currentfile.readline()
    if line == '':
        break
    S_background.append(float(line))


#..Alpha Van Genuchten
inputfile = './Parameters/Calibrated/param_Alpha_VG.dat'

alpha = []

currentfile = open(inputfile,'r')
while 1:
    line = currentfile.readline()
    if line == '':
        break
    alpha.append(float(line))

#..Beta Van Genuchten
inputfile = './Parameters/Calibrated/param_Beta_VG.dat'

beta = []

currentfile = open(inputfile,'r')
while 1:
    line = currentfile.readline()
    if line == '':
        break
    beta.append(float(line))

#..Recharge rate
inputfile = './Parameters/Calibrated/param_R.dat'

Rch = []

currentfile = open(inputfile,'r')
while 1:
    line = currentfile.readline()
    if line == '':
        break
    Rch.append(float(line))

rch = Rch[0] / (1000*365.0) # m/d

###############################
# Other non adjustable parameters for now
###############################

Ss = [1.00E-04]

Sat_res = 0.05

u_ghb = [99.47] # loosing river
d_ghb = [93.22] # loosing river

Rb_anisotropy = [4.0]
Aq_anisotropy = [4.0]

VG_param = {'alpha' : alpha[0], 'beta': beta[0], 'Swr': Sat_res, 'lp': 0.5}

####################### WRITE MPROPS FILE ###########################
#####################################################################
#####################################################################

Stop_exp = True

minP = 0.01
maxP = 20
n_layer_vg = 150.0

if Stop_exp == True:
    max_dP = 0.5
    stop_exp = 1

P = P_func(minP,maxP,n_layer_vg)

if Stop_exp == True:
    # Do a postproc to stop exponential at stop_exp and continue depth calculation following the regular space defined by max_dz
    mask = P <= stop_exp
    P = list(P[mask])
    linear_spacing = list(np.arange(P[-1], maxP, max_dP))[1:]
    P.extend(linear_spacing)

P = list(np.array(P) *-1)

P.insert(0,0)
P.reverse()

Sw = []
Kr = []
Se = []

for p in P:
    c_aquifer_sw, c_aquifer_kr, c_aquifer_se = vg_func(p,VG_param['alpha'],VG_param['beta'],VG_param['Swr'],VG_param['lp'])
    Sw.append(c_aquifer_sw)
    Kr.append(c_aquifer_kr)
    Se.append(c_aquifer_se)

mprops = {'S':S_background[0],'alpha':alpha[0],'beta':beta[0],'theta_r':Sat_res,'Ss':Ss[0]}

write_mprop(mprops,P,Sw,Kr)

if plot == True:
    fig = plt.figure(figsize=(20,15))
    plt.rc('font', family='serif', size=20)
    x = P
    # PLOT SATURATION
    ax1 = fig.add_subplot(111)
    y = Sw
    lns1 = ax1.plot(x, y, linewidth=6, color='red',label = 'Saturation')
    ax1.set_xlim(0,-20)
    ax1.set_ylim(0,1)
    ax1.grid()
    # PLOT RELATIVE PERMEABILITY
    ax2 = ax1.twinx()
    y = Kr
    lns2 = ax2.plot(x, y, linewidth=6, color='blue',label = 'Permeability')
    ax2.set_ylim(0,1)
    plt.tight_layout()
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns,labs,fontsize=30)
    fig.subplots_adjust(hspace=0)
    fig.savefig('./figs/VG_curves.png', dpi=600)

################
# Write elemental values
################

Kfile = open('./HGS/Parameters/K.txt','w')

n = 0

for i in range(len(E)):
    HGSelem = E[i]
    #if HGSelem in Riverbed_elements:
    #    k = Krb[0]
    #    Kfile.write(str(HGSelem) + ' ' + str(k) + ' ' + str(k) + ' ' + str(k/Aq_anisotropy[0]) +'\n')
    #    continue
    k = K_list[n]
    Kfile.write(str(HGSelem) + ' ' + str(k) + ' ' + str(k) + ' ' + str(k/Aq_anisotropy[0]) +'\n')
    n = n+1

Kfile.close()

Sfile = open('./HGS/Parameters/S.txt','w')

n = 0

for i in range(len(E)):
    HGSelem = E[i]
    #if HGSelem in Riverbed_elements:
    #    s = Srb[0]
    #    Sfile.write(str(HGSelem) + ' ' + str(s) + ' ' + str(s) + ' ' + str(s) +'\n')
    #    continue
    s = S_list[n]
    Sfile.write(str(HGSelem) + ' ' + str(s) + ' ' + str(s) + ' ' + str(s) +'\n')
    n = n+1

Kfile.close()

####################### WRITE BC FILE ##############################
#####################################################################
#####################################################################

write_GHB(u_ghb[0],d_ghb[0])

rch_inc(rch)

######################################
# WRITE OBSERVATION POINT FILE
######################################

z_coordobs = [85]

x_coordobs_plain = [15,35,55,75,95,115,135,155,175,195,215,235,255,275]
y_coordobs_plain = [5,50,100,150,200,250,300,350,400,450,495]

XX_obs,YY_obs = np.meshgrid(list(x_coordobs_plain),list(y_coordobs_plain))

Grid_obs = list(zip(XX_obs.ravel(), YY_obs.ravel())) 

obs_idx = [16,26,46,57,82,86,109,128]
reduced_Grid_obs = [Grid_obs[idx] for idx in obs_idx]
Grid_obs = reduced_Grid_obs

if plot == True:
    fig = plt.figure()
    ax = plt.axes()
    plt.axis('scaled')
    #ax.set_xlim(min(Y_grid), max(Y_grid))
    #ax.set_ylim(min(X_grid), max(X_grid))
    ax.scatter(XX_obs+100,YY_obs+100, s=10, marker='+',color='black', linewidths = 0.7)
    ax.invert_yaxis()
    ax.set_xlim(0,500)
    ax.set_ylim(700,0,-10)
    ax.add_patch(Rectangle((100,100),300,500,color='black',alpha=0.2))
    ax.xaxis.tick_top()
    fig.subplots_adjust(hspace=0)
    figure_name = './figs/obs_loc.pdf'
    fig.savefig(figure_name, dpi=600)

outputfile= open('./SIM/ObservationPoints.dat', 'w')
outputfile.write('X;' + 'Y' + '\n')

outputfile_obs = './HGS/Grokfiles/observation_pts.inc'
currentfile = open(outputfile_obs,'w')

n_obs = 1

for z in z_coordobs:
    for i in range(len(Grid_obs)):
        point = Grid_obs[i]
        x = point[0]
        y = point[1]
        outputfile.write(str(x) + ';' + str(y) + '\n')
        currentfile.write('' + '\n')
        currentfile.write('' + '\n')
        currentfile.write('Make observation point' + '\n')
        currentfile.write('obs_' + str(n_obs) + '\n')
        currentfile.write(str(x) + ' ' + str(y) + ' ' + str(z) + '\n')
        currentfile.write('' + '\n')
        currentfile.write('clear chosen zones' + '\n')
        currentfile.write('clear chosen nodes' + '\n')
        currentfile.write('clear chosen elements' + '\n')
        currentfile.write('clear chosen faces' + '\n')
        currentfile.write('clear chosen segments' + '\n')
        n_obs = n_obs+1

outputfile.close()

outputfile= open('./HGS/Grokfiles/ObservationPoints_selection.dat', 'w')

for z in z_coordobs:
    for i in range(len(Grid_obs)):
        point = Grid_obs[i]
        x = point[0]
        y = point[1]
        outputfile.write(str(x) + ';' + str(y) + ';' + str(85) + '\n')

outputfile.close()

sys.exit()

######################################
# WRITE PARTICLE TRACING POINT FILE
######################################

z_coordobs = list(np.arange(75,90,1)) 

x_coordobs_plain = list(np.arange(100,290,1)) 
y_coordobs_plain = [500]

XX_obs,YY_obs,ZZ_obs = np.meshgrid(list(x_coordobs_plain),list(y_coordobs_plain),list(z_coordobs))

Grid_obs = list(zip(XX_obs.ravel(), YY_obs.ravel(), ZZ_obs.ravel())) 

outputfile= open('./HGS/Grokfiles/particles.dat', 'w')

for i in range(len(Grid_obs)):
    current_coords = Grid_obs[i]
    outputfile.write(str(current_coords[0]) + ' ' + str(current_coords[1]) + ' ' + str(current_coords[2]) + ' ' + '1' + ' ' + '0.0' + '\n')

outputfile.close()





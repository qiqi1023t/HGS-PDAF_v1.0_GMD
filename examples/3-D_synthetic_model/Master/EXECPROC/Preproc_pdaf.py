import string,sys,math,os, time
from operator import itemgetter, attrgetter
import datetime
from math import *
import numpy as np
import pickle

sep = 1

Tstep_times = list(np.arange(364,460,1))

# WRITE SEQUENTIAL WELL FLOW FILE

# ========= READ BATCH FILE ==============

batchfile_elements = []

currentfile = open('./HGS/Grokfiles/batchfiles/wells_flow.inc','r')

while 1:
    line = currentfile.readline()[:-1]
    if line == '!endbatch':
        break
    batchfile_elements.append(line)

pumping_list = []

currentfile = open('./HGS/Grokfiles/pumping.dat','r')

while 1:
    line = currentfile.readline()[:-1]
    if line == '':
        break
    pumping = line.split()[1]
    pumping_list.append(pumping)

# ========= WRITE SEQUENTIAL GROKFILES ==============

for t in range(len(Tstep_times)):
    print(t)
    tstep = t+1
    sequentialfile = open('./HGS/Grokfiles/wells/wells_flow_' + str(tstep) + '.inc','w')
    flag = False
    for i in range(len(batchfile_elements)):
        line = batchfile_elements[i]
        if flag == True:
            sequentialfile.write('     0.0' + ' ' + pumping_list[t] + '\n')
            flag = False
            continue
        if line == '    !adaptive':
            flag = True
            continue
        sequentialfile.write(line + '\n')
    sequentialfile.close()


# ========= INITIALIZED OUTPUTFILES ==============

nobs = 8

for i in range(nobs):
    outputfile = open('./SIM/H/obs' + str(i+1) + '.dat','w')
    outputfile.close()


sys.exit()


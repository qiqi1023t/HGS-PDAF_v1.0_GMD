import string,sys,math,os,shutil,time
from operator import itemgetter, attrgetter
import datetime
from math import *
import numpy as np
import pickle


Ensemble_size = 100

def chunks(l, n):
    for i in range(0, len(l),n):
        yield l[i:i + n]


for real in list(range(1,Ensemble_size+1)):
    reals_list = []
    currentfile = open('./Prior_Ensemble.csv','r')
    line = currentfile.readline()
    while 1:
        line = currentfile.readline()
        if line == '':
            break
        n = float(line.split(',')[0])
        if n == real:
            K, S = list(chunks(line.split(',')[1:],204000))
            break
    currentfile.close()
    shutil.copytree('./Master','./' + str(real).zfill(3))
    currentfile = open('./' + str(real).zfill(3) + '/K.dat','w')
    for k in K:
        currentfile.write(str(k) + '\n')
    currentfile.close()
    currentfile = open('./' + str(real).zfill(3) + '/S.dat','w')
    for s in S:
        currentfile.write(str(s) + '\n')
    currentfile.close()

sys.exit()


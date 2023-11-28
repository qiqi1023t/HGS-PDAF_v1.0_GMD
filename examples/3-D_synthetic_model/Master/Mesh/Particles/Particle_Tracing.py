import string,sys,math,os,numpy, time
import numpy as np

nodes_2d = 7015

# ========= READ 2D RIVER NODES ID ==============

Nodes_list = []

currentfile = open('./streamnodes.txt','r')

while 1:
    line = currentfile.readline()[:-1]
    if line == '':
        break
    node = line.split()
    Nodes_list.append(float(node[0]))

currentfile.close()

# ========= Transform 2D IDENTS 3D IDENTS ==============

nsheet = 16

Nodes_3D_list = np.array(Nodes_list) + (nodes_2d*(nsheet-1))

outputfile = './Stream_node_set.dat'
currentfile = open(outputfile,'w')

for node in Nodes_3D_list:
    currentfile.write(str(node) +'\n')

currentfile.close()

# ========= Transform Z coordinates from SuperGrid ==============

inputfile = '../SuperGrid.dat'

currentfile = open(inputfile,'r')

while 1:
    line = currentfile.readline()[:-1]
    if line == '# x':
        break

X = []

while 1:
    line = currentfile.readline()[:-1]
    if line == '# y':
        break
    temp = line.split()
    for i in range(len(temp)):
        X.append(float(temp[i]))

Y = []

while 1:
    line = currentfile.readline()[:-1]
    if line == '# z':
        break
    temp = line.split()
    for i in range(len(temp)):
        Y.append(float(temp[i]))

Z = []

while 1:
    line = currentfile.readline()[:-1]
    if line == '# propnum (cell-centered)':
        break
    temp = line.split()
    for i in range(len(temp)):
        Z.append(float(temp[i]))


# ========= WRITE PARTICLE TRACING INPUT FILE ==============

outputfile = './particles.dat'
currentfile = open(outputfile,'w')

for node in Nodes_3D_list:
    x = X[int(node)-1]
    y = Y[int(node)-1]
    z = Z[int(node)-1]
    currentfile.write(str(x) + ' ' + str(y) + ' ' + str(z) + ' ' + '1 0.0' + '\n')

currentfile.close()



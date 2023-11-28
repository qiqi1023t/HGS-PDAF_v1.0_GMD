import numpy as np

ah2file = open('../SuperMesh.ah2','r')

line = ah2file.readline()
n_nodes = float(line)

for n in range(int(n_nodes)):
    line = ah2file.readline()

n_elem = int(ah2file.readline())

elemfile = open('./mesh_HGS.ele','w')

elemfile.write( str(n_elem) + ' ' + '3' + ' 0' + '\n')

for i in range(n_elem):
    line = ah2file.readline()[:-1]
    elemfile.write( str(i+1) + ' ' + str(line) + '\n')

elemfile.close()






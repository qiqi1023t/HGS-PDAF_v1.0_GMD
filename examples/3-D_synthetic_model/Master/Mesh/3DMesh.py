import numpy as np
import string,sys,math,os, time

sep = 1

# Read depths

Depths = []

currentfile = open('./Depths.dat','r')

while 1:
    line = currentfile.readline()[:-sep]
    if line == '':
        break
    depth = float(line)
    Depths.append(depth)


n_layer = len(Depths)-1


def grid_inc(n_layer):
    outputfile = './3DGrid.inc'
    currentfile = open(outputfile,'w')
    currentfile.write('Generate layers interactive ' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('    base elevation' + '\n' )
    currentfile.write('      elevation from raster file' + '\n' )
    currentfile.write('        ../Mesh/elevation.asc' + '\n' )
    currentfile.write('    offset base' + '\n' )
    currentfile.write('     -30' + '\n' )
    currentfile.write('    end base elevation' + '\n' )
    currentfile.write('' + '\n' )
    currentfile.write('' + '\n' )
    for i in range(n_layer):
        currentfile.write('new layer' + '\n' )
        currentfile.write('    layer name' + '\n' )
        currentfile.write('      Layer' + str(i)  + '\n' )
        currentfile.write('    Minimum layer thickness' + '\n' )
        currentfile.write('      0.0001' + '\n' )
        currentfile.write('    elevation from raster file' + '\n' )
        currentfile.write('      ../Mesh/Layers/Layer' + str(i) + '.asc' + '\n' )
        currentfile.write('end' + '\n' )
        currentfile.write('' + '\n' )
    currentfile.write('new layer' + '\n' )
    currentfile.write('    layer name' + '\n' )
    currentfile.write('      Layer' + str(n_layer)  + '\n' )
    currentfile.write('    Minimum layer thickness' + '\n' )
    currentfile.write('      0.0001' + '\n' )
    currentfile.write('      elevation from raster file' + '\n' )
    currentfile.write('        ../Mesh/elevation.asc' + '\n' )
    currentfile.write('      end' + '\n' )
    currentfile.write('end' + '\n' )
    currentfile.write('end grid generation' + '\n' )
    currentfile.close()

# Write grid.inc
grid_inc(n_layer)


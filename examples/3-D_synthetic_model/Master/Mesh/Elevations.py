import numpy as np
from scipy.interpolate import interp1d
from math import exp
from math import log

Stop_exp = True

# 2D Elevation Mesh Generation

ncols = 400
nrows = 600
xllcorner = -50
yllcorner = -50
cellsize = 1
NoData = -9999

slope = 3/1000.0
channel_depth = 3.0 # m
channel_width = 15 # m

origin_z = 100.0 #m

x = np.linspace(100, 100-(500.0*slope), num= nrows, endpoint=True)

outputfile= open('./elevation.asc', 'w')

# Headers
outputfile.write('ncols' + ' ' + str(ncols) + '\n')
outputfile.write('nrows' + ' ' + str(nrows) + '\n')
outputfile.write('xllcorner' + ' ' + str(xllcorner) + '\n')
outputfile.write('yllcorner' + ' ' + str(yllcorner) + '\n')
outputfile.write('cellsize' + ' ' + str(cellsize) + '\n')
outputfile.write('NODATA_value' + ' ' + str(NoData) + '\n')

for line in range(nrows):
    plaine = np.repeat(x[line],ncols-abs(xllcorner)-channel_width)
    channel = np.repeat(x[line]-channel_depth,channel_width+abs(yllcorner))
    to_write = np.append(plaine,channel)
    for val in range(len(to_write)):
        current_val = str(to_write[val])
        outputfile.write(current_val + ' ')
    outputfile.write('\n')

outputfile.close()

# 3D Mesh Generation

depths = [0.50,1.0,1.5,2.0,3.0,4.0,5.0,8.0,11.0,14.0,17.0,20.0,23.0,26.0]

outputfile= open('Depths.dat', 'w')

for d in depths:
    outputfile.write(str(d) + '\n')

outputfile.write(str(30.0) + '\n')

outputfile.close()

for i in range(len(depths)):
    current_depth = depths[::-1][i]
    outputfile= open('./Layers/Layer' + str(i) + '.asc', 'w')
    # Headers
    outputfile.write('ncols' + ' ' + str(ncols) + '\n')
    outputfile.write('nrows' + ' ' + str(nrows) + '\n')
    outputfile.write('xllcorner' + ' ' + str(xllcorner) + '\n')
    outputfile.write('yllcorner' + ' ' + str(yllcorner) + '\n')
    outputfile.write('cellsize' + ' ' + str(cellsize) + '\n')
    outputfile.write('NODATA_value' + ' ' + str(NoData) + '\n')
    for line in range(nrows):
        plaine = np.repeat(x[line],ncols-abs(xllcorner)-channel_width)
        channel = np.repeat(x[line]-channel_depth,channel_width+abs(yllcorner))
        to_write = np.append(plaine,channel) - current_depth
        for val in range(len(to_write)):
            current_val = str(to_write[val])
            outputfile.write(current_val + ' ')
        outputfile.write('\n')
    outputfile.close()

    

    










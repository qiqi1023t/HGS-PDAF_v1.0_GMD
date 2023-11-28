import numpy as np

def readhgsoutput(prefix, variable, time='001'):
    """Read output from a HydroGeoSphere simulation, return a numpy array
    
    Args:
        prefix: String containing the prefix of the output files.
        variable: String specifying which variable to read.
            Can be one of the following:
                'coordinates_pm'
                'elements_pm'
                'head_pm'
                'sat_pm'
                'v_pm'
                'q_pm'
                'ElemK_pm'
                'coordinates_olf'
                'indices_olf'
                'elements_olf'
                'head_olf'
                'v_olf'
                'ExchFlux_olf'
                'element_centroids_pm'
                'element_centroids_olf'
                'num_layers'
        time: The output time to read"""

    if (variable == 'coordinates_pm'):
        with open (prefix + 'o.' + variable, 'rb') as f:
            f.seek(4)
            numpoints = np.fromfile(file=f, dtype='<i4', count=1)
            f.seek(8, 1)
            data = np.fromfile(file=f, dtype='<f8', count=numpoints*3).reshape(-1,3)

    elif (variable == 'indices_olf'):
        with open (prefix + 'o.coordinates_olf', 'rb') as f:
            f.seek(4)
            # numpoints = np.fromfile(file=f, dtype='<i4', count=1)
            nodesPerElement= np.fromfile(file=f, dtype='<i4', count=1)  #changed Shanghua
            f.seek(8, 1)
            numberOfElements = np.fromfile(file=f, dtype='<i4', count=1)#changed Shanghua
            f.seek(8, 1)
            # changed Shanghua
            data = np.fromfile(file=f, dtype='<i4', count=numberOfElements*nodesPerElement).reshape(-1,1)
            # data = np.fromfile(file=f, dtype='<i4', count=numpoints).reshape(-1,1)
            
    elif (variable == 'coordinates_olf'):
        coordinates_pm = readhgsoutput(prefix, 'coordinates_pm')
        with open (prefix + 'o.' + variable, 'rb') as f:
            f.seek(4)
            numpoints = np.fromfile(file=f, dtype='<i4', count=1)
            f.seek(8, 1)
            ids_olf = np.fromfile(file=f, dtype='<i4', count=numpoints).reshape(-1,1)
            data = coordinates_pm[ids_olf-1].reshape(-1,3);
            
    elif (variable == 'elements_pm'):
        with open (prefix + 'o.' + variable, 'rb') as f:
            f.seek(4)
            nodesPerElement = np.fromfile(file=f, dtype='<i4', count=1)
            f.seek(8, 1)
            numberOfElements = np.fromfile(file=f, dtype='<i4', count=1)
            f.seek(8, 1)
            data = np.fromfile(file=f, dtype='<i4', count=numberOfElements*nodesPerElement).reshape(-1,nodesPerElement)


    elif (variable == 'element_centroids_pm'):
        elements = readhgsoutput(prefix, 'elements_pm')
        coords = readhgsoutput(prefix, 'coordinates_pm')
        data = np.mean(coords[elements-1],1)

    elif (variable == 'element_centroids_olf'):
        elements = readhgsoutput(prefix, 'elements_olf')
        coords = readhgsoutput(prefix, 'coordinates_pm')
        data = np.mean(coords[elements-1],1)

    elif (variable == 'element_coordinates_olf'):
        elements = readhgsoutput(prefix, 'elements_olf')
        coords = readhgsoutput(prefix, 'coordinates_pm')
        data = coords[elements-1]
            
    elif (variable == 'elements_olf'):      
        with open (prefix + 'o.' + variable, 'rb') as f:
            f.seek(4)
            nodesPerElement = np.fromfile(file=f, dtype='<i4', count=1)
            f.seek(8, 1)
            numberOfElements = np.fromfile(file=f, dtype='<i4', count=1)
            f.seek(8, 1)
            data = np.fromfile(file=f, dtype='<i4', count=numberOfElements*nodesPerElement).reshape(-1,nodesPerElement)
            if data[0,3] == 0:
                data = data[:,0:3]
    elif (variable == 'num_layers'):
        nodes_olf = np.shape(readhgsoutput(prefix, 'coordinates_olf'))[0]
        nodes_pm = np.shape(readhgsoutput(prefix, 'coordinates_pm'))[0]
        data = nodes_pm/nodes_olf
    
    # 64bit floats
    elif (variable in {'head_pm' , 'head_olf'}):
        with open (prefix + 'o.' + variable + '.' + time, 'rb') as f:
            headerLength = np.fromfile(file=f, dtype='<i4', count=1)
            f.seek(headerLength+8)
            
            numpoints = np.fromfile(file=f, dtype='<i4', count=1) / 8
            
            # TODO should this be like the following line?
            # numpoints = int(np.fromfile(file=f, dtype='<i4', count=1) / 8)
            data = np.fromfile(file=f, dtype='<f8', count=numpoints).reshape(-1,1)

    # 32bit floats            
    elif (variable in {'sat_pm', 'ExchFlux_olf'}):
        with open (prefix + 'o.' + variable + '.' + time, 'rb') as f:
            headerLength = np.fromfile(file=f, dtype='<i4', count=1)
            f.seek(headerLength+8)
            numpoints = np.fromfile(file=f, dtype='<i4', count=1) / 4
            data = np.fromfile(file=f, dtype='<f4', count=numpoints).reshape(-1,1)      

    # 32bit 3-component vectors 
    elif (variable in {'v_pm', 'q_pm', 'v_olf'}):
        with open (prefix + 'o.' + variable + '.' + time, 'rb') as f:
            count = np.fromfile(file=f, dtype='<i4', count=1)
            f.seek(count+8)
            data = np.fromfile(file=f, dtype='<f4').reshape(-1,5)[:,1:4]

    # 64bit 3-component vectors 
    elif (variable in {'ElemK_pm'}):
        with open (prefix + 'o.' + variable + '.' + time, 'rb') as f:
            count = np.fromfile(file=f, dtype='<i4', count=1)
            f.seek(count+8)
            data = np.fromfile(file=f, dtype=('i4,3<f8,i4'))
            data = np.array( [entry[1] for entry in data])
    
    else:
        raise ValueError('Not a recognised variable: ' + variable)
    return data

def createLinearInterpolator(prefix, variable, time='001'):
    from scipy.interpolate import LinearNDInterpolator
    
    input_coords = readhgsoutput(prefix, 'coordinates_pm')
    input_var = readhgsoutput(prefix, variable, time)
    interpolator = LinearNDInterpolator(input_coords, input_var)
    return interpolator

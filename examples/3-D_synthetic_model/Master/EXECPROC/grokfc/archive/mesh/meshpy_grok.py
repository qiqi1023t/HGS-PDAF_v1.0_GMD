#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
this script attempts to

- generate a mesh using meshpy
- generating the necessary grok files
    - .coordinates
    - .elements
    - o.elemental_k


meshpy
- doku: https://documen.tician.de/meshpy/index.html
- at Andreas Klöckner's webpage: https://mathema.tician.de/software/meshpy/
- on github: https://github.com/inducer/meshpy
    - particularly the function `DoTriMesh` from https://github.com/inducer/meshpy/blob/master/examples/jw_meshtools.py
"""
import sys
import os
import datetime
import numpy as np
import numpy.linalg as la
import matplotlib
import matplotlib.pyplot as plt

import jw_meshtools as jwm
import meshpy.triangle as triangle

import shapely.geometry as shpgeom

### General comment by Philipp = BEWARE: Python counts from 0, Fortran and HGS from 1 ###

__author1__ = "Philipp Selzer (Philipp.Selzer@gmx.net)"
__author2__ = "Claus Haslauer (mail@planetwater.org)"
__version__ = "$Revision: 2.0 $"
__date__ = datetime.date(2017,4,27)
__copyright__ = "Copyright (c) 2017 Claus Haslauer, Philipp Selzer"
__license__ = "Python"


### SCOPE:
### THE FOLLOWING PROGRAM WRITES HYDROGEOSPHERE INPUT FILES (.coordinates, .elements,
### o.elemental_k) FOR MODELS WITH TRIANGULAR PRISMATIC ELEMENTS (PARALLELOGRAMS OF
### TRIANGULATED PLANES) INCLUDING N ARBITRARILIY SHAPED ZONES OF ALTERED K IN N LAYERS.
### THE EXAMPLE PROGRAM DOES NOT ACOOUNT FOR SLOPED INTERFACES. HOWEVER, THE PROGRAMS IS WRITTEN
### IN A WAY THAT SLOPED INTERFACES CAN BE INCLUDED EASILY. (ALL YOU HAVE TO DO IS TO DEFINE YOUR
### SLOPES AND ALTER THE Z-VALUES ACCORDINGLY. THEN THE PROGRAM SHOULD HANDLE THE DEFINED SLOPES
### AUTOMATICALLY WITHOUT PROBLEMS OR FURTHER CHANGES.)


def main():
    #matplotlib.use('Agg')

    plt.ioff()
    
    # --------------------------------------------------------------------------------
    #                                                                       build mesh
    # --------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------
    # EXAMPLE WITH TWO RECTANGLES defining mesh boundary
    #  with inclusion 
    #length = 2.0
    #domain_ll, domain_ur = [0,0],[10,10]
    #inclusion_ll, inclusion_ur = [3.5, 3.5],[6.5, 6.5]
    #p1,v1=jwm.RectangleSegments(domain_ll, domain_ur, num_points=100)
    #p2,v2=jwm.RectangleSegments(inclusion_ll, inclusion_ur, num_points=50)
    #p,v=jwm.AddCurves(p1,v1,p2,v2)
    #mesh = jwm.DoTriMesh(p,v,edge_length=length, min_angle=20) # change by Philipp
    #mesh_points = np.array(mesh.points)
    #mesh_elements = np.array(mesh.elements)

    reactangular_inclusion = True 
    other_inclusions = True
    Two_D = False

    ## DEFINE POLYGONS WITH ALTERED K-VALUE ##

    ld_corner = (3.5,3.5)
    lu_corner = (3.5,6.5)
    ru_corner = (6.5,6.5)
    rd_corner = (3.5,6.5)

    first_polygon = shpgeom.Polygon([ld_corner,lu_corner,ru_corner,rd_corner])

    ld_corner2 = (7.5,7.5)
    lu_corner2 = (7.5,10.0)
    ru_corner2 = (10.0,10.0)
    rd_corner2 = (10.0,7.5)

    second_polygon = shpgeom.Polygon([ld_corner2,lu_corner2,ru_corner2,rd_corner2])

    polygons = [first_polygon, second_polygon]

    # background k for all layers, can also be seen as starting value to get the array filled..
    background_k_xx = 1.0e-04
    background_k_yy = 1.0e-04
    background_k_zz = 1.0e-04

    ## ALTERED K-VALUES IN THE DEFINITION ORDER OF POLYGONS FOR EVERY LAYER##
    ##-> Every layer gets an own row in the matrix. If k should not be altered
    ##   just insert the background k.
    # note: here 3 nodal layers define 2 FE layers
    #                         FE layer 1        FE layer 2
    polygon_ks_xx = np.array([[1.0e-06,1.0e-8],[1.0e-06,1.0e-8]])
    polygon_ks_yy = np.array([[1.0e-06,1.0e-8],[1.0e-06,1.0e-8]])
    polygon_ks_zz = np.array([[1.0e-06,1.0e-8],[1.0e-06,1.0e-8]])


    # plt.subplot(111, aspect='equal')
    # plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_elements)

    # --------------------------------------------------------------------------------
    # EXAMPLE WITH ONE RECTANGLE defining mesh boundary
    length = 3.0
    p,v=jwm.RectangleSegments([0.0,0.0],[10.,10.],edge_length=length)
    mesh  = jwm.DoTriMesh(p,
                          v,
                          edge_length=length,
                          min_angle=20.)
                          # ,
                          # 
                          # 
                      
    mesh_points = np.array(mesh.points)
    mesh_elements = np.array(mesh.elements)
    
    

    print (type(mesh))
    print("===============================")
    plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_elements,) 
    print (mesh_points.shape)
    n_nodes_per_layer = mesh_points.shape[0]
    print("===============================")
    print (mesh_elements.shape)
    n_elem_per_layer = mesh_elements.shape[0]

    # plotting
    fig, axs = plt.subplots(nrows=1,
                            ncols=1,
                            dpi=150)
    axs.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_elements,)
    plt.xlim(0.0, 10.)
    plt.ylim(0.0, 10.)
    plt.savefig('../mesh.png')
    
    n_fe_layers = 2
    n_node_layers = n_fe_layers + 1

    n_nodes_tot = n_nodes_per_layer * n_node_layers 

    ### Define z-coordinates ###
    z_layers_def = np.array([0.0, 1.0, 2.0])
    # z-coordinate as matrix manual:
    z_layer1 = np.ones(n_nodes_per_layer)*z_layers_def[0]
    z_layer2 = np.ones(n_nodes_per_layer)*z_layers_def[1]
    z_layer3 = np.ones(n_nodes_per_layer)*z_layers_def[2]
    z_layers_v = np.array([z_layer1,z_layer2,z_layer3])
    z_layers = z_layers_v.T

    # z-coordinate as matrix quick:     
    z_layers = np.tile(z_layers_def,(n_nodes_per_layer,1))
    
 
    n_elem_tot = n_elem_per_layer * n_fe_layers 

 
    # i/o names
    prefix = 'case_neu' 
    coord_fname = prefix + '.coordinates'
    elem_fname = prefix + '.elements'

    k_fname = prefix + 'o.elemental_k' 

    
    n_lines_max_x = 0
    n_lines_max_y = 0
    n_lines_max_z = 2
    n_species = 0
    n_nodes_per_element = 6
    
    
    # write coordinates file
    coord_fobj = os.path.join('./', coord_fname)
    with open(coord_fobj, 'w') as f:
        f.write("{:d}\n".format(n_nodes_tot)) 
        for cur_fe_layer in range(n_node_layers):
            for cur_node in range(n_nodes_per_layer):
                f.write("{:f} {:f} {:f}\n".format(mesh_points[cur_node, 0], mesh_points[cur_node, 1], z_layers[cur_node][cur_fe_layer]))
        cur_line = "{:d} {:d} {:d} {:d}\n".format(n_lines_max_x, n_lines_max_y, n_lines_max_z, n_species)
        f.write(cur_line)
        f.write("{:d}\n".format(n_elem_tot))
        false_line = ".false.\n"
        f.write(false_line)
        f.write(false_line)
        f.write(false_line)
        last_line = "9\n"
        f.write(last_line)

    mesh_elements_write = mesh_elements + 1
            
            
    # write elements file
    elems_fobj = os.path.join('./', elem_fname)
    with open(elems_fobj, 'w') as g:
        g.write("{:d}\n".format(n_nodes_per_element))
        g.write("{:d}\n".format(n_elem_tot))
        for cur_fe_layer in range(n_fe_layers): 
            for cur_elem in range(n_elem_per_layer):    
                    g.write("{} {} {} {} {} {}\n".format(mesh_elements_write[cur_elem][0] + (cur_fe_layer*n_nodes_per_layer),           
                                                    mesh_elements_write[cur_elem][1] + (cur_fe_layer*n_nodes_per_layer),
                                                    mesh_elements_write[cur_elem][2] + (cur_fe_layer*n_nodes_per_layer),
                                                    mesh_elements_write[cur_elem][0] + (cur_fe_layer*n_nodes_per_layer) + n_nodes_per_layer,
                                                    mesh_elements_write[cur_elem][1] + (cur_fe_layer*n_nodes_per_layer) + n_nodes_per_layer,
                                                    mesh_elements_write[cur_elem][2] + (cur_fe_layer*n_nodes_per_layer) + n_nodes_per_layer))
            
        for cur_id in range(n_elem_tot):
                g.write("{}\n".format(1))  


    # write o.elemental_k file (by Philipp):


    # background k            
    k_xx = background_k_xx*np.ones(n_elem_tot)
    k_yy = background_k_yy*np.ones(n_elem_tot)
    k_zz = background_k_zz*np.ones(n_elem_tot)



    # inclusions k
    shape_coordinates = np.shape(mesh_points)
    len_coordinates = shape_coordinates[0]

    shape_elements = np.shape(mesh_elements)
    len_elements = shape_elements[0]

    for ith_layer in range(n_fe_layers):

        for ith_poly in range(len(polygons)):
    
            indices_in_rectangular_incl_coord = []
    
            for find_p_in_p in range(len_coordinates):
                current_np_point = mesh_points[find_p_in_p,]
                current_shp_point = shpgeom.Point(current_np_point)
                if current_shp_point.within(polygons[ith_poly]):
                    indices_in_rectangular_incl_coord.append(find_p_in_p)
                    #print(mesh_points[find_p_in_p,])

            indices_of_elements_in_polygon = []
        
            for find_p_in_e in range(len_elements):
                first_node = mesh_elements[find_p_in_e,0]
                second_node = mesh_elements[find_p_in_e,1]
                third_node = mesh_elements[find_p_in_e,2]

                first_where = np.where(indices_in_rectangular_incl_coord==first_node)[0]
                second_where = np.where(indices_in_rectangular_incl_coord==second_node)[0]
                third_where = np.where(indices_in_rectangular_incl_coord==third_node)[0]

                first_len = len(first_where)
                second_len = len(second_where)
                third_len = len(third_where)

                if (first_len==1 or second_len==1 or third_len==1):
                    indices_of_elements_in_polygon.append(find_p_in_e)
                    #print(find_p_in_e)
        
            k_indices = np.array([indices_of_elements_in_polygon]) + (ith_layer*n_elem_per_layer)
            k_xx[[k_indices]]=polygon_ks_xx[ith_layer][ith_poly]
            k_yy[[k_indices]]=polygon_ks_yy[ith_layer][ith_poly] 
            k_zz[[k_indices]]=polygon_ks_zz[ith_layer][ith_poly] 


    k_fobj = os.path.join('./', k_fname)
    with open(k_fobj, 'w') as h:
        for cur_k in range(n_elem_tot):
            h.write("{} {} {} {}\n".format(cur_k + 1,
                                           k_xx[cur_k],
                                           k_yy[cur_k],
                                           k_zz[cur_k]))

    print("Done! Yay!")
if __name__ == '__main__':
    main()



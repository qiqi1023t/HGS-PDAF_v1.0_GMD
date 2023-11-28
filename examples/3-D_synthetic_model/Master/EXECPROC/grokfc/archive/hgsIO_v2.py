# -------------------------------------------------------------------------------
# Name:        hydrogeosphere I/O
# Author:      guthke
# Created:     08.03.2012
# Copyright:   (c) guthke 2012
# -------------------------------------------------------------------------------
# !/usr/bin/env python
import os
import numpy as np
import pylab as plt
import scipy.stats as st
import datetime
import shutil

"""
grok file making has to be checked... timestepping...
"""


class hgsIO(object):
    def __init__(
            self,
            folder=None,  # current working directory or specified one
            prefix='cube',  # hgs prefix
            talk_to_me=True,  # print things
            hgs_type='hgs',  # string that stands for the OS command that calls hydrogeosphere
            grok_type='grok',  # string that stands for the OS command that calls grok (preprocessor for hydrogeosphere)
            plot_type='gsplot'  # string that stands for the OS command that calls postprocessor for hydrogeosphere
    ):
        self.prefix = prefix
        self.hgs_type = hgs_type
        self.grok_type = grok_type
        self.plot_type = plot_type
        self.talk_to_me = talk_to_me
        if folder == None:
            self.folder = os.path.curdir
        else:
            self.folder = folder
            if os.path.isdir(folder) == False:
                os.mkdir(folder)
        self.write_batch()

    ##-------------------------------------------------------------------------#
    def mgridstyle_read_elemental_h(self):
        # load domainproperties and mprops
        """
        this has to be checked...
        """
        self.get_domainproperties()
        self.get_mprops()

        # load element node list
        if hasattr(self, 'element_node_list') == False:
            self.get_element_node_list()

        # GET EVERYTHING IN LISTSTYLE (if necessary interpolate to elementcenter)
        # coordinates
        self.get_nodal_xyz()
        x = (self.helper_node_to_elementcentroid(self.nodal_xyz[:, 0]))
        y = (self.helper_node_to_elementcentroid(self.nodal_xyz[:, 1]))
        z = (self.helper_node_to_elementcentroid(self.nodal_xyz[:, 2]))
        self.elemental_xyz = np.array((x, y, z)).T

        # pressure heads
        self.get_nodal_heads()
        self.elemental_head = self.helper_node_to_elementcentroid(self.nodal_head)

        # TRANSFORM TO MGRIDSTYLE
        xyz = np.copy(self.elemental_xyz)
        x = xyz[:, 0]
        y = xyz[:, 1]
        z = xyz[:, 2]
        h = np.copy(self.elemental_head)

        for ii in range(3)[::-1]:
            ind = np.argsort(xyz[:, ii], kind='mergesort')
            xyz = xyz[ind]
            h = h[ind]

        h = h.reshape(self.elementnumbers)

        return h

    ##-------------------------------------------------------------------------#
    def mgridstyle_read_elemental_xyzKhqqq(self):
        # load domainproperties and mprops
        """
        this has to be checked...
        """
        self.get_domainproperties()
        self.get_mprops()

        # load element node list
        if hasattr(self, 'element_node_list') == False:
            self.get_element_node_list()

        # GET EVERYTHING IN LISTSTYLE (if necessary interpolate to elementcenter)
        # coordinates
        self.get_nodal_xyz()
        x = (self.helper_node_to_elementcentroid(self.nodal_xyz[:, 0]))
        y = (self.helper_node_to_elementcentroid(self.nodal_xyz[:, 1]))
        z = (self.helper_node_to_elementcentroid(self.nodal_xyz[:, 2]))
        self.elemental_xyz = np.array((x, y, z)).T

        # K
        self.get_elemental_K()

        # pressure heads
        self.get_nodal_heads()
        self.elemental_head = self.helper_node_to_elementcentroid(self.nodal_head)
        # velocities
        self.get_elemental_velocity()

        # TRANSFORM TO MGRIDSTYLE
        xyz = np.copy(self.elemental_xyz)
        x = xyz[:, 0]
        y = xyz[:, 1]
        z = xyz[:, 2]
        K = np.copy(self.elemental_K)
        h = np.copy(self.elemental_head)
        qx = self.elemental_velocity[:, 0]
        qy = self.elemental_velocity[:, 1]
        qz = self.elemental_velocity[:, 2]
        for ii in range(3)[::-1]:
            ind = np.argsort(xyz[:, ii], kind='mergesort')
            xyz = xyz[ind]
            x = x[ind]
            y = y[ind]
            z = z[ind]
            K = K[ind]
            h = h[ind]
            qx = qx[ind]
            qy = qy[ind]
            qz = qz[ind]

        x = x.reshape(self.elementnumbers)
        y = y.reshape(self.elementnumbers)
        z = z.reshape(self.elementnumbers)
        K = K.reshape(self.elementnumbers)
        h = h.reshape(self.elementnumbers)
        qx = qx.reshape(self.elementnumbers)
        qy = qy.reshape(self.elementnumbers)
        qz = qz.reshape(self.elementnumbers)

        return x, y, z, K, h, qx, qy, qz

    ##-------------------------------------------------------------------------#
    def mgridstyle_read_elemental_tc(self, solutename='conservative tracer'):
        # load domainproperties and mprops
        """
        this has to be checked...
        """
        self.get_domainproperties()
        self.get_mprops()

        # load element node list
        if hasattr(self, 'element_node_list') == False:
            self.get_element_node_list()

        # coordinates
        self.get_nodal_xyz()
        x = (self.helper_node_to_elementcentroid(self.nodal_xyz[:, 0]))
        y = (self.helper_node_to_elementcentroid(self.nodal_xyz[:, 1]))
        z = (self.helper_node_to_elementcentroid(self.nodal_xyz[:, 2]))
        self.elemental_xyz = np.array((x, y, z)).T

        # timestep files
        filelist = self.ret_timestep_files(solutename)

        # concentration and timesteps
        tstps_all = []
        c_all = []
        for fname in filelist:
            c_nodal, t = self.ret_nodal_concentrations(fname)  # read c
            c_elemental = self.helper_node_to_elementcentroid(c_nodal)

            # TRANSFORM TO MGRIDSTYLE
            xyz = np.copy(self.elemental_xyz)
            c = np.copy(c_elemental)
            for ii in range(3)[::-1]:
                ind = np.argsort(xyz[:, ii], kind='mergesort')
                xyz = xyz[ind]
                c = c[ind]
            c = c.reshape(self.elementnumbers)

            tstps_all.append(t)
            c_all.append(c)

        return np.array(tstps_all), np.array(c_all)

    ##-------------------------------------------------------------------------#
    def run(self, gsplot=False):
        cd = os.path.abspath(os.path.curdir)
        os.chdir(self.folder)
        os.system(self.grok_type)
        os.system(self.hgs_type)
        if gsplot == True: os.system(self.plot_type)
        os.chdir(cd)

    ##-------------------------------------------------------------------------#
    def del_all_files_in_folder(self):
        for fname in os.listdir(self.folder):
            fpath = os.path.join(self.folder, fname)
            try:
                os.unlink(fpath)
            except Exception as e:
                print(e)

    ##-------------------------------------------------------------------------#
    def del_runfolder(self):
        shutil.rmtree(self.folder)

    ##-------------------------------------------------------------------------#
    def get_element_node_list(self):
        '''
        gets a list with the id of every node that belongs to an element id
        [numberofelements x 8]
        '''
        if self.talk_to_me == True: print('read element node list:')
        ii = []

        f = open(os.path.join(self.folder, self.prefix + "o.elements"), 'rb')
        i, ii = self.helper_reader(f, ii)

        nodes_per_element = f.read(i)
        nodes_per_element = np.fromstring(nodes_per_element, dtype='int32')[0]
        if self.talk_to_me == True: print('\tnodes_per_element: ', nodes_per_element)
        i, ii = self.helper_reader(f, ii)
        i, ii = self.helper_reader(f, ii)

        anz = i / 4
        liste = np.fromfile(f, dtype='int32', count=anz)
        liste = liste - 1  # hgs --> python indexing
        self.element_node_list = liste.reshape(np.shape(liste)[0] / nodes_per_element, nodes_per_element)
        i, ii = self.helper_reader(f, ii)

        if ii[0] == ii[1]:
            if ii[2] != ii[3]:
                print('\tReading somehow failed!')

    ##-------------------------------------------------------------------------#
    def get_domainproperties(self):
        '''
        reads from scrach_grok
             domainsize      (for each direction),
             elementnumbers  (for each direction),
             elementsize     (for each direction),
        '''
        if self.talk_to_me == True: print('read domainproperties:')
        DomainSize = []
        NumberElements = []
        read_domain_properties = False
        fobj = open(os.path.join(self.folder, 'scratch_grok'), 'r')
        for line in fobj:
            if read_domain_properties == True:
                horst = line.split()
                DomainSize.append(float(horst[0]))
                NumberElements.append(int(horst[1]))
                if len(DomainSize) == 3:
                    self.domainsize = np.array(DomainSize)
                    self.elementnumbers = np.array(NumberElements)
                    self.elementsize = self.domainsize / self.elementnumbers
                    if self.talk_to_me == True:
                        print('\tdomainsize:      ', self.domainsize)
                        print('\telementnumbers:  ', self.elementnumbers)
                        print('\telementsize:     ', self.elementsize)
                    read_domain_properties = False
                    break
            if line[0:23] == 'generate uniform blocks':
                read_domain_properties = True
        fobj.close()

    ##-------------------------------------------------------------------------#
    def get_mprops(self):
        '''
        reads from scach_mprops
        '''
        if self.talk_to_me == True: print('read material properties:')
        read_porosity = False
        read_alphax = False
        read_alphay = False
        read_alphaz = False
        fobj = open(os.path.join(self.folder, 'scratch_mprops'), 'r')
        for line in fobj:
            # POROSITY
            if read_porosity == True:
                self.porosity = float(line.split()[0])
                read_porosity = False
            elif line[0:8] == 'porosity':
                read_porosity = True
            # ALPHA_X
            if read_alphax == True:
                self.alphax = float(line.split()[0])
                read_alphax = False
            elif line[0:25] == 'longitudinal dispersivity':
                read_alphax = True
            # ALPHA_Y
            if read_alphay == True:
                self.alphay = float(line.split()[0])
                read_alphay = False
            elif line[0:23] == 'transverse dispersivity':
                read_alphay = True
            # ALPHA_Z
            if read_alphaz == True:
                self.alphaz = float(line.split()[0])
                read_alphaz = False
            elif line[0:32] == 'vertical transverse dispersivity':
                read_alphaz = True

        if self.talk_to_me == True:
            print('\tporosity:                         ', self.porosity)
            print('\tlongitudinal dispersivity:        ', self.alphax)
            print('\ttransverse dispersivity:          ', self.alphay)
            print('\tvertical transverse dispersivity: ', self.alphaz)
        fobj.close()

    ##-------------------------------------------------------------------------#
    def get_elemental_K(self):
        if self.talk_to_me == True: print('read K elemental...')
        ii = []

        f = open(os.path.join(self.folder, self.prefix + 'o.elemental_k'), 'rb')

        # header
        i, ii = self.helper_reader(f, ii)
        header = f.read(i)
        i, ii = self.helper_reader(f, ii)

        # Kxx Kyy Kzz
        KxxKyyKzz = []
        for ee in range(np.prod(self.elementnumbers)):  # loop over all ids
            i, ii = self.helper_reader(f, ii)
            string = f.read(i)
            KxxKyyKzz.append(np.fromstring(string, dtype='float64'))
            i, ii = self.helper_reader(f, ii)
        self.elemental_K = np.array(KxxKyyKzz)[:, 0]

    ##-------------------------------------------------------------------------#
    def get_nn_ne_per_layer(self):
        '''
        gets number of nodes (nn) and number of elements (ne) per layer
        designed for lerma catchment (carolin, 20150717)
        '''
        import re

        if self.talk_to_me == True: print('read eco file...')
        ii = []

        f = open(os.path.join(self.folder, self.prefix + "o.eco"), 'r')

        text = f.readlines()

        search_nn = '\s*Number of nodes\s*(\d+)'
        search_ne = '\s*Number of elements\s*(\d+)'

        # find nn
        for cur_line in text:
            match_nn = re.match(search_nn, cur_line)
            match_ne = re.match(search_ne, cur_line)
            if match_nn:
                # print(match_nn.group(1))
                nn = int(match_nn.group(1))
            if match_ne:
                ne = int(match_ne.group(1))

        return nn, ne

    ##-------------------------------------------------------------------------#
    def get_nodal_xyz(self):
        '''
        gets coordinates of every node id [n_n x 3]
        '''
        if self.talk_to_me == True: print('read coordinates...')
        ii = []

        f = open(os.path.join(self.folder, self.prefix + "o.coordinates_pm"), 'rb')
        i, ii = self.helper_reader(f, ii)

        n_nodes = f.read(i)
        n_nodes = np.fromstring(n_nodes, dtype='int32')[0]
        i, ii = self.helper_reader(f, ii)
        i, ii = self.helper_reader(f, ii)

        anz = i / 8
        xyz = np.fromfile(f, dtype='float64', count=anz)
        self.nodal_xyz = xyz.reshape(np.shape(xyz)[0] / 3, 3)

        i, ii = self.helper_reader(f, ii)

        if ii[0] == ii[1]:
            if ii[2] != ii[3]:
                print('\tReading somehow failed!')

    ##-------------------------------------------------------------------------#
    def get_nodal_heads(self):
        '''
        gets a one d list with a head for every node id
        '''
        if self.talk_to_me == True: print('read heads...')
        ii = []

        f = open(os.path.join(self.folder, self.prefix + "o.head.001"), 'rb')
        i, ii = self.helper_reader(f, ii)

        horst = f.read(i)

        i, ii = self.helper_reader(f, ii)
        i, ii = self.helper_reader(f, ii)

        anz = i / 8
        self.nodal_head = np.fromfile(f, dtype='float64', count=anz)

        i, ii = self.helper_reader(f, ii)

        if ii[0] == ii[1]:
            if ii[2] != ii[3]:
                print('\tReading somehow failed!')

    ##-------------------------------------------------------------------------#
    def get_elemental_velocity(self):
        '''
        gets a one d list with a velocity vor every element id
        '''
        if self.talk_to_me == True: print('read elemental velocities...')
        vel = np.fromfile(os.path.join(self.folder, self.prefix + "o.vel"), dtype='float32')
        self.elemental_velocity = vel[1:-1].reshape(np.shape(vel)[0] / 3, 3)

    ##-------------------------------------------------------------------------#
    def ret_timestep_files(self, solutename):
        '''
        gets filelist with path+filename (of concentrations at timesteps)
        '''
        if self.talk_to_me == True: print('read concentration files...')
        filelist = os.listdir(self.folder)

        clist = []
        for item in filelist:
            if item[:-3] == self.prefix + 'o.concentration.' + solutename + '.':
                clist.append(os.path.join(self.folder, item))
        return clist

    ##-------------------------------------------------------------------------#
    def ret_nodal_concentrations(self, filepath):
        '''
        returns:    * nodal concentrations (a 1d list with a concetration for every node id)
                    * timestep
        '''
        if self.talk_to_me == True: print('read nodal concentrations', end=' ')
        ii = []

        f = open(filepath, 'rb')

        i, ii = self.helper_reader(f, ii)

        # read header (timesep)
        t = f.read(i)
        t = float(t)
        if self.talk_to_me == True: print('at timestep=%f' % t)

        i, ii = self.helper_reader(f, ii)
        i, ii = self.helper_reader(f, ii)

        anz = i / 8
        c = np.fromfile(f, dtype='float64', count=anz)

        i, ii = self.helper_reader(f, ii)

        f.close()

        return c, t

    ##-------------------------------------------------------------------------#
    def helper_reader(self, fobj, ii=[]):
        i = fobj.read(4)
        i = np.fromstring(i, dtype='int32')[0]
        ii.append(i)
        return i, ii

    ##-------------------------------------------------------------------------#
    def helper_node_to_elementcentroid(self, parameter):
        centroidvalues = []
        for line in self.element_node_list:
            Avg = np.mean(parameter[line])
            centroidvalues.append(Avg)
        centroidvalues = np.array(centroidvalues)
        return centroidvalues

    ##-------------------------------------------------------------------------#
    def write_idKKK_from_Pmgrid(self, KXXgrid, KYYgrid=None, KZZgrid=None):
        '''
        if Kyy or Kzz = None, Kxx is taken instead
        '''
        xyzKxx = self.conv_Pmgrid_to_XYZPlist(KXXgrid)
        if KYYgrid != None:
            xyzKyy = self.conv_Pmgrid_to_XYZPlist(KYYgrid)
        else:
            xyzKyy = None
        if KZZgrid != None:
            xyzKzz = self.conv_Pmgrid_to_XYZPlist(KZZgrid)
        else:
            xyzKzz = None
        self.write_idKKK_from_xyzklist(xyzKxx, xyzKyy, xyzKzz)

    ##-------------------------------------------------------------------------#
    def write_idKKK_from_xyzklist(self, xyzKxx, xyzKyy=None, xyzKzz=None):
        '''
        if Kyy or Kzz = None, Kxx is taken instead
        '''
        if self.talk_to_me == True: print('write idkkk file')
        xyzKxx = self.sort_XYZPlist_HGSstyle(xyzKxx)
        if xyzKyy != None:
            xyzKyy = self.sort_XYZPlist_HGSstyle(xyzKyy)
        else:
            xyzKyy = xyzKxx
        if xyzKzz != None:
            xyzKzz = self.sort_XYZPlist_HGSstyle(xyzKzz)
        else:
            xyzKzz = xyzKxx
        idkkk = np.concatenate((
            np.arange(xyzKxx.shape[0])[:, np.newaxis] + 1,
            xyzKxx[:, -1][:, np.newaxis],
            xyzKyy[:, -1][:, np.newaxis],
            xyzKzz[:, -1][:, np.newaxis]), axis=1)
        np.savetxt(os.path.join(self.folder, 'idkkk.dat'), idkkk)

    ##-------------------------------------------------------------------------#
    def conv_Pmgrid_to_XYZPlist(self, grid, dx=1, dy=1, dz=1):
        """
        self.elementsize statt dx, dy, dz ... mit abfrage obs da ist???
        """
        grid = np.atleast_3d(grid)
        gridslice = [slice(0, grid.shape[i], 1) for i in range(grid.ndim)]
        xyz = np.array(np.mgrid[gridslice])
        xyz = ((xyz.flatten()).reshape(grid.ndim, -1)).T
        xyz[:, 0] = xyz[:, 0] * dx
        xyz[:, 1] = xyz[:, 1] * dy
        xyz[:, 2] = xyz[:, 2] * dz
        xyzp = np.concatenate((xyz, grid.flatten()[:, np.newaxis]), axis=1)
        return xyzp

    ##-------------------------------------------------------------------------#
    def sort_XYZPlist_HGSstyle(self, xyzp):
        for ii in range(min(xyzp.shape[1] - 1, 3)):
            ind = np.argsort(xyzp[:, ii], kind='mergesort')
            xyzp = xyzp[ind]
        return xyzp

    ##-------------------------------------------------------------------------#
    def write_grok(self,
                   domainlength=(100, 50, 1),
                   elementnumber=(100, 50, 1),

                   units='Kilogram-metre-day',  # 'Kilogram-centimetre-hour'
                   # 'Kilogram-metre-day'
                   headgradient=0.01,  # in x-direction

                   xyzttc=np.array([
                       [5.0, 25.0, 0.0, 0.0, 1.0, 1.0],
                       [5.0, 25.0, 1.0, 0.0, 1.0, 1.0]
                   ]),
                   flux='Specified concentration',  # 'Specified concentration' --> dirichlet
                   # 'Specified mass flux'     --> neumann
                   # 'initial concentration'   --> anfangsbedingung (tt of xyzttc useless!!)
                   diffusion=1.0e-10,  # free-solution diffusion coefficient

                   c_ctrl=None,  # Concentration Contol

                   tstps=[10000],  # output timesteps
                   maxtstp=None,  # maximum timestep (if None --> 1% of max(tstps))
                   mintstp=10e-8,  # minimum timestep

                   force_mass=False,  # force a concentration in upper right corner!
                   make_wells=False,
                   ):

        grokstring = ''
        # ###############################
        # TITLE
        # ###############################
        grokstring += '''!--------------------------  Problem description
%s
automatically generated with hgsIO
end title''' % str(datetime.datetime.now())

        # ###############################
        # GRID DEFINITION
        # ###############################
        grokstring += '''
!--------------------------  Grid generation
generate uniform blocks
%f %i 	     ! Domain length and number of blocks in X
%f %i		 ! Domain length and number of blocks in Y
%f %i		 ! Domain length and number of blocks in Z
end grid definition''' % (domainlength[0], elementnumber[0],
                          domainlength[1], elementnumber[1],
                          domainlength[2], elementnumber[2],
                          )

        # ###############################
        # GENERAL SIMULATION CTRL
        # ###############################
        grokstring += '''
!--------------------------  General simulation parameters
Finite difference mode
!finite differenzen methode wird benutzt(anstelle finite elemente methode - default)
Control volume
units: %s
do transport
!echo to output''' % units

        # ###############################
        # POROUS MEDIA PROPERTIES
        # ###############################
        grokstring += '''
!--------------------------  Porous media properties
use domain type
porous media

properties file
%s.mprops

clear chosen zones
choose zone number
1
read properties
porous medium

clear chosen elements
choose elements all
Read elemental k from file
idkkk.dat
clear chosen elements''' % self.prefix

        # ###############################
        # FLOW CTRL
        # ###############################
        # calculate heads at boundaries
        h0 = domainlength[0] * float(headgradient)
        h1 = 0

        x0 = 0
        x1 = domainlength[0]
        grokstring += '''\n\n
!=======================================
!==============  F L O W  ==============
!=======================================
!----------- INITIAL FLOW CONDITIONS -----------------------------------------------------------
!clear chosen nodes
!choose nodes all
!initial head
!30.0

!----------- BC --------------------------------------------------------------------------------
!----------- BC am Rand des Gebietes (1.TYP / DIRICHLET)----------------------------------------
clear chosen nodes
choose nodes x plane
%f             		! xk-Koordinate der Ebene
1.0e-5           		! Distanz-Kriterium (in diesem Abstand werden Punkte ebenfalls beruecksichtigt)
specified head
1							! Anzahl der Ranbedingungen im Zeitverlauf
0.0 %f	 			! Startzeitpunkt, spez.Druckhoehe

clear chosen nodes
choose nodes x plane
%f          	 		! x-Koordinate der Ebene
1.e-5           		! Distanz-Kriterium (in diesem Abstand werden Punkte ebenfalls beruecksichtigt)
specified head
1				! Anzahl der Ranbedingungen im Zeitverlauf
0.0 %f				! Startzeitpunkt, spez.Druckhoehe''' % (x0, h0,  # left boundary condition
                                                              x1, h1)  # right boundary condition

        # ###############################
        # TRANSPORT CONTROL
        # ###############################
        # BUILT CONCENTRATION BOUNDARY STRING OUT OF xyzttc:
        # compute xyzttc, so that it is compatible with grok:
        # - cumulate coordinates...
        # - no overlapping of times etc...
        # - correct shape
        # make a numpy array
        xyzttc = np.array(xyzttc)

        # if only one point: make new axis for same dimensionality
        if len(xyzttc.shape) == 1:
            xyzttc = xyzttc[np.newaxis]

        # sort ing
        for i in range(xyzttc.shape[1])[::-1]:
            xyzttc = xyzttc[np.argsort(xyzttc[:, i], axis=0, kind='mergesort')]

        # check same space positions and
        # write to extra array
        dim = 3
        allx = []
        onex = []
        horst = xyzttc[0, :3]
        for i in range(xyzttc.shape[0]):
            if (xyzttc[i, :dim] == horst[:dim]).all():
                onex.append(xyzttc[i])
            else:
                allx.append(np.array(onex))
                onex = []
                onex.append(xyzttc[i])
            if i == xyzttc.shape[0] - 1:
                allx.append(np.array(onex))
                onex = []
            horst = xyzttc[i, :3]
        # allx is for now the new shit!

        # cumulate same space time fractions of concentration intensities
        finalarray = []  # this will be the array without overlapping
        for item in allx:
            x = item[0, 0]
            y = item[0, 1]
            z = item[0, 2]
            t0t1c = []
            ts = (item[:, 3:5]).flatten()
            ts = np.sort(list(set(np.round(ts, 10))))  # jeder eintrag sollte nur einmal vorkommen!
            for i in range(ts.shape[0] - 1):
                ind = np.where((item[:, 3] < ts[i + 1]) & (item[:, 4] > ts[i]))[0]
                if ind.shape[0] != 0:
                    t0 = ts[i]
                    t1 = ts[i + 1]
                    c = item[ind, -1].sum()
                    t0t1c.append(np.array([t0, t1, c]))
            t0t1c = np.array(t0t1c)
            # now add to final array
            finalarray.append([x, y, z, t0t1c])

        # -----------------------------
        # MAKE CONCENTRATION BOUNDARIES:
        cboundstr = '\n'

        for coord in finalarray:
            # initial concentration:
            if flux == 'initial concentration':
                cboundstr += 'clear chosen nodes\nchoose node\n%f %f %f\n%s\n%f\n\n' % (coord[0],
                                                                                        coord[1],
                                                                                        coord[2],
                                                                                        flux,
                                                                                        coord[3][0][-1]
                                                                                        )


                # specified mass flux & specified concentration:
            else:
                cboundstr += 'clear chosen nodes\nchoose node\n%f %f %f\n%s\n%i\n' % (coord[0],
                                                                                      coord[1],
                                                                                      coord[2],
                                                                                      flux,
                                                                                      coord[3].shape[0]
                                                                                      )
                for step in coord[3]:
                    cboundstr += '%f, %f, %1.9f\n' % (step[0], step[1], step[2])
                cboundstr += '\n'

        # --------------------------------------------------------------------------
        if force_mass == True:
            # FORCE SYSTEM TO HAVE MASS IN IT THAT IT NOT TERMINATES BEFORE FIRST OR
            # LAST TIMESTEP. Source is in upper right corner!
            # check out the minimum starting time
            mint0 = xyzttc[:, 3].min()
            # and if > 0 make some source for having mass in system
            if mint0 > 0.:
                cboundstr += "\n! This is for havin' mass in system!\n"
                cboundstr += 'clear chosen nodes\nchoose node\n%f %f %f\n%s\n1\n0.0, %f, %1.9f' % (domainlength[0],
                                                                                                   domainlength[1],
                                                                                                   domainlength[2],
                                                                                                   flux,
                                                                                                   # this is a string
                                                                                                   mint0,
                                                                                                   # endtime is first starting time
                                                                                                   xyzttc[:,
                                                                                                   5].min() / 1000.
                                                                                                   # conc is .1% of lowest conc
                                                                                                   )
                cboundstr += "\n! This is for havin' mass in system!\n"
            # check out the maximum ending time
            maxt1 = xyzttc[:, 4].max()
            # and if > 0 make some source for having mass in system
            if maxt1 < max(tstps):
                cboundstr += "\n! This is for havin' mass in system!\n"
                cboundstr += 'clear chosen nodes\nchoose node\n%f %f %f\n%s\n1\n%f, %f, %1.9f' % (domainlength[0],
                                                                                                  domainlength[1],
                                                                                                  domainlength[2],
                                                                                                  flux,
                                                                                                  # this is a string
                                                                                                  maxt1,
                                                                                                  # last endingtime time is startingtime
                                                                                                  max(tstps),
                                                                                                  # endtime is last timestep
                                                                                                  xyzttc[:,
                                                                                                  5].min() / 1000.
                                                                                                  # conc is .1% of lowest conc
                                                                                                  )
                cboundstr += "\n! This is for havin' mass in system!\n\n"
        # --------------------------------------------------------------------------

        # SOLUTE DEFINITION
        grokstring += '''\n\n
!===============================================
!============== T R A N S P O R T ==============
!===============================================
Solute
name
conservative tracer
free-solution diffusion coefficient
%e
End Solute\n\n''' % diffusion

        # CONCENTRATION INPUT
        grokstring += cboundstr

        # TRANSPORT CONVERGENCE PARAMETER
        grokstring += '''
Upstream weighting of velocities
0 0 0
!0 = no upstream; 1 = fully upstream; x y z direction

Transport time weighting
0.5
!0 = explicit; 0.5 = Crank Nicholson; 1 = fully implicit'''

        # ------------
        # TIMESTEPPING
        # ------------
        if maxtstp == None:
            maxtstp = np.array(tstps).max() / 100.0
        tstpstr = ''
        for t in tstps:
            tstpstr += '%f\n' % t

        grokstring += '''\n\n
!=================\n!= TIME STEPPING =\n!=================
Maximum timestep
%s
Minimum timestep
%s

output times
%send
    ''' % (maxtstp, mintstp, tstpstr)

        if c_ctrl != None:
            grokstring += '''
Concentration control
%f''' % c_ctrl

            # ###############################
            # SOLVER
            # ###############################
        grokstring += '''\n\n
!=====================\n!= SOLVER PARAMETERS =\n!=====================
!flow solver convergence criteria
!1.0e-12
!transport solver convergence criteria
!1.0d-12
!Flux limiter for transport'''

        # ###############################
        # OUTPUT CTRL
        # ###############################
        grokstring += '''\n\n
!==================\n!= OUTPUT CONTROL =\n!==================
!echo incidences
!echo coordinates
'''
        if make_wells == True:
            # WELL STRING
            wells = '! OUTPUT NODES:'
            for i, xwellcoord in enumerate((0.25 * domainlength[0], 0.5 * domainlength[0], 0.75 * domainlength[0])):
                for j, ywellcoord in enumerate((0.25 * domainlength[1], 0.5 * domainlength[1], 0.75 * domainlength[1])):
                    wells += '''
clear chosen nodes
choose node
%f, %f, %f
choose node
%f, %f, %f
Slice flux output nodes from chosen
horst_node%i%i
Slice flux contributing nodes from chosen
    ''' % (xwellcoord, ywellcoord, 0,
           xwellcoord, ywellcoord, domainlength[2],
           i, j)

            grokstring += '''
    %s
    ''' % wells

            # --------------------------------------------------------------------------
            # finally write grok file
            # --------------------------------------------------------------------------
        fobj = open(os.path.join(self.folder, self.prefix + '.grok'), 'w')
        fobj.write(grokstring)
        fobj.close()

        if self.talk_to_me == True: print('%s.grok file saved' % self.prefix)

    ##-------------------------------------------------------------------------#
    def write_mprops(self,
                     porosity=0.35,
                     alpha_x=2.0,
                     alpha_y=0.2,
                     alpha_z=0.2,
                     ):
        txt = '''Porous medium
!k isotropic
!7.20576
porosity
%f
specific storage
1.0e-6
longitudinal dispersivity
%f
transverse dispersivity
%f
vertical transverse dispersivity
%f
tortuosity
0.8
bulk density
2650.0
end material
    ''' % (porosity, alpha_x, alpha_y, alpha_z)
        mpropsobj = open(os.path.join(self.folder, self.prefix + '.mprops'), 'w')
        mpropsobj.write(txt)
        mpropsobj.close()
        if self.talk_to_me == True: print('%s.mprops saved' % self.prefix)

    def write_batch(self):
        batchobj = open(os.path.join(self.folder, 'batch.pfx'), 'w')
        batchobj.write(self.prefix)
        batchobj.close()
        if self.talk_to_me == True: print('batch.pfx saved')

    def write_array_sizes(self):
        fobj = open(os.path.join(self.folder, 'array_sizes.default'), 'w')
        fobj.write(
            '''
    channel flow: 1d elements
           50000
    channel flow: material zones
              20
    channel flow bc: zero-depth gradient segments
            5000
    dual flow bc: flux faces
           10000
    dual flow bc: flux function panels
              10
    dual flow bc: flux nodes
           10000
    dual flow bc: flux zones
              10
    dual flow bc: head function panels
             100
    dual flow bc: head nodes
           10000
    dual: material zones
              20
    flow: material zones
              20
    flow bc: drain-type flux nodes
               2
    flow bc: evaporation faces
           10000
    flow bc: evaporation nodes
           10000
    flow bc: evaporation zones
              10
    flow bc: evaporation function panels
              10
    flow bc: flux nodes
           10000
    flow bc: flux faces
           10000
    flow bc: flux zones
              10
    flow bc: flux function panels
              10
    flow bc: free drainage nodes
            1000
    flow bc: head nodes
           10000
    flow bc: head function panels
             100
    flow bc: river-type flux nodes
               2
    flow bc: specified nodal flowrate
             501
    flow bc: specified nodal flowrate function panels
             100
    flow bc: hydrostatic node columns
             100
    heat transfer permafrost: thawing table
              50
    heat transfer permafrost: freezing table
              50
    heat transfer permafrost: thawing-freezing table
              50
    heat transfer permafrost: temperature function panels
             300
    fractures: 2d elements
          100000
    fractures: zones
             300
    general: list
             300
    mesh: node connections
             100
    mesh: node sheets in z for layered grids
              50
    mesh: x grid lines (rectangular)
            1000
    mesh: y grid lines (rectangular)
            1000
    mesh: z grid lines (rectangular)
            1000
    mesh: number of layers
             100
    mesh: number of sublayers per layer
             100
    observation wells: nodes
             501
    output: flux volume nodes
            1000
    output: flux volumes
              10
    output: nodes
             100
    output: times
            1000
    permafrost : elements
           10000
    permafrost : function panels
             100
    seepage face: 3d elements intersecting
            1000
    seepage face: nodes
            1000
    solution: target times
            50000
    surface flow: 2d elements
           50000
    surface flow: boundary segments
            5000
    surface flow: hydrographs
              20
    surface flow: hydrograph nodes
             100
    surface flow: material zones
              20
    surface flow bc: critical depth segments
            5000
    surface flow bc: zero-depth gradient segments
            5000
    tile drains: 1d elements
           10000
    tile drains: 3d elements intersecting
           10000
    tile drains: concentration function panels
             100
    tile drains: nodes
            1000
    transport: species
               5
    transport: species kinetic reactions
               2
    transport bc: concentration nodes
           10000
    transport bc: concentration function panels
             100
    transport bc: flux nodes
          100000
    transport bc: flux function panels
            1000
    transport bc: immiscible phase dissolution nodes
            1000
    transport bc: third-type concentration faces
           10000
    transport bc: third-type concentration function panels
             100
    transport bc: zero-order source function panels
             100
    transport bc: first-order source function panels
             100
    wells: 1d elements
            1000
    wells: 2d fracture elements intersecting
            1000
    wells: 3d elements intersecting
            1000
    wells: flux function panels
              10
    tiles: flux function panels
              20
    wells: injection concentration function panels
             100
    wells: nodes
             100
    stress : stressed nodes
           10000
    stress : stress function panels
             100
    end
    ''')
        fobj.close()
        if self.talk_to_me == True: print('array_sizes.default saved')

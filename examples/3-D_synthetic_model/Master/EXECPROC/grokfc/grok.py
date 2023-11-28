"""
Functions to deal with GROK and HGS
Created on Wed Mar 25 20:39:31 2015
@author: HydroPy Tubingen
"""
import sys
import os
import numpy as np
import datetime
import re
import platform
import grokfc.diradmin as dm


# import shutil
# import subprocess as sp
# import multiprocessing as mp
# import time
# import scipy
# from itertools import repeat
# from statistics import mean
# import scipy.sparse as ssp
# import fieldgen as fg
# import dirmanip as dm
# import gridmanip as gm
# import mystats as mysst
# import kalman as klm
# import pdb


class Grokfile(object):  # todo: test when passing a normal string rather than an os.path or r'' path
    """ Manipulate strings within an existent grok file
    Arguments:
        model_dir:  str, full path to the directory containing the grok file
        grokname:   str, grok file identifier
    """

    def __init__(self, model_dir, grokname='', mkbatch=False, mkletmerun=False):
        self.model_dir = model_dir
        self.grokname = grokname
        self.mkbatch = mkbatch
        self.mkletmerun = mkletmerun
        if self.mkbatch is True:
            self.createprefix()
            self.createletmerun()

    def creategrok(self, description='No description'):
        """ Create a modelname.grok file, including main headers for a basic hgs model
        Returns:

        """

        grokfile, tempgrokfile = self.findgrok()
        assert grokfile == 'NotaFile', 'Grok file < %s.grok > already exists.' % self.grokname
        print('Creating it...')
        mainheaders = \
            """start title
        !--- end title

        !--- grid generation
        !--- end grid generation

        !--- general simulation parameters
        !--- end general simulation parameters

        !--- porous media properties
        !--- end porous media properties

        !--- overland properties
        !--- end overland properties

        !--- initial conditions flow
        !--- end initial conditions flow

        !--- solute properties
        !--- end solute properties

        !--- initial conditions transport
        !--- end initial conditions transport

        !--- well properties
        !--- end well properties

        !--- boundary conditions flow
        !--- end boundary conditions flow

        !--- boundary conditions transport
        !--- end boundary conditions transport

        !--- boundary conditions well
        !--- end boundary conditions well

        !--- tracer injection
        !--- end tracer injection

        !--- output times
        !--- end output times

        !--- output locations
        !--- end output locations
        """  # todo: pass this to a support file...
        # Create the grok file in disk:
        with open(os.path.join(self.model_dir, '%s.grok' % self.grokname), 'w') as newgrok:
            newgrok.write(mainheaders)
        # Add Problem description:
        self.modelTitle(description=description)

        print('working process')

    def findgrok(self):
        """

        Return:
             String for a temporal grok file, based on an existent grok file
        """
        dirObj = dm.Manipdirs(self.model_dir)
        grokfile, dummy = dirObj.getdirs(mystr='%s.grok' % self.grokname, fullpath=True)

        if type(grokfile) == np.ndarray:
            if len(grokfile) > 0:
                grokfile = grokfile[0]
            elif len(grokfile) == 0:
                grokfile = 'NotaFile'
                print('Grok file < %s.grok >  not found' % self.grokname)
        # if not isinstance(grokfile, str):
        #     sys.exit('Error in finding the grok file, check directory')

        return grokfile, grokfile + 'temp'

    def createprefix(self):  # todo: write a function to check if the file exists
        """ Create a prefix file to run GROK, HGS and/ HSPLOT
        Returns:
            The batch file batch.pfx with model prefix
        """
        prefFile = 'batch.pfx'
        if prefFile in os.listdir(self.model_dir):
            os.remove(os.path.join(self.model_dir, prefFile))

        prefix = open(os.path.join(self.model_dir, prefFile), 'w')
        prefix.write(self.grokname)
        prefix.close()

        print('Prefix file < %s > created' % prefFile)

    def createletmerun(self, hsplot=False):
        """ Create a letmerun batch file for running GROK, HGS and HSPLOT
        Arguments:
            hsplot:     bool, if True: include instruction to run hsplot
        Returns:
            Batch file letmerun.bat in the given directory
        """
        batchfile = 'letmerun.bat'
        if batchfile in os.listdir(self.model_dir):
            os.remove(os.path.join(self.model_dir, batchfile))

        letmerun = open(os.path.join(self.model_dir, batchfile), 'w')
        letmerun.write('@echo off\n\ndel *' + self.grokname + 'o* > nul\n\n')
        if platform.system() == 'Linux':
            letmerun.write('grok.x\nhgs.x\n')
        else:
            letmerun.write('grok\nphgs\n')

        if hsplot is True:
            if platform.system() == 'Linux':
                letmerun.write('hsplot.x\n')
            else:
                letmerun.write('hsplot\n')

        letmerun.close()

        print('Batch file < %s > created' % batchfile)

    def addtime2grok(self, outputtimes):
        """ Add/change output times to an existent GROK file
        Arguments:
            outputtimes:    np.array or list of the output time values
        Returns:
            GROK file with updated time step
        """
        if type(outputtimes) == list:
            outputtimes = np.asarray(outputtimes)

        grokfile, tempgrokfile = self.findgrok()

        writing = True
        with open(grokfile) as f:

            with open(tempgrokfile, 'w') as out:
                for line in f:
                    if writing:
                        if "output times" in line:
                            writing = False
                            out.write("output times")
                            for ii in range(0, len(outputtimes)):
                                out.write('\n%s' % (outputtimes[ii]))
                            buffer = [line]
                            out.write("\nend\n")
                        else:
                            out.write(line)
                    elif "end" in line:
                        writing = True
                    else:
                        buffer.append(line)
                else:
                    if not writing:
                        # There wasn't a closing "event", so write buffer contents
                        out.writelines(buffer)
                out.close()
                f.close()
        os.remove(grokfile)
        os.rename(tempgrokfile, grokfile)
        print('Output times in grok file < %s > updated...' % str(grokfile))

    # @staticmethod
    def add2grok(self, pattern_looking_for, string_to_add):
        """
        example:
        add another well to existing well, just before `pattern_looking_for`
            pattern_looking_for = "!--- End Well properties"
            string_to_add = "choose segments polyline\n2\n0.000000, 10.000000, 1.000000\n0.000000, 10.000000, 0.000000\nread properties\nTheis\n"

        """
        # open current grok file and read as one string
        fobj, tempgrokfile = self.findgrok()
        with open(fobj, 'r') as f:
            text = f.read()

        text2replace = string_to_add + '\n' + pattern_looking_for

        result = re.sub(pattern_looking_for,
                        text2replace,
                        text)
        # write text with replacement in place back to file
        with open(fobj, 'w') as f:
            f.writelines(result)

    def addstr2grok(self, startstr='', endstr='', str2add=[]):
        """ Add a string (or list of strings) between to given strings to an existent grok file.
        Arguments:
            startstr:
            endstr:
            str2add:
        :return:
        """

        grokfile, tempgrokfile = self.findgrok()

        fileObj = dm.Manipfiles(grokfile)
        oldstr_idx = fileObj.find_str(startstr, index=True)
        # oldstr_idx = str_infile(grokfile, startstr, index=True)  # todo: reference properly this function

        if oldstr_idx is None:
            print('String < %s > not found in <%s>' % (startstr, os.path.split(grokfile)[-1]))
        elif oldstr_idx is not None:
            writing = True
            with open(grokfile) as f:
                with open(tempgrokfile, 'w') as out:
                    for line in f:
                        if writing:
                            if startstr in line:
                                writing = False
                                out.write(startstr)
                                for ii in range(0, len(str2add)):
                                    out.write('\n%s' % (str2add[ii]))
                                buffer = [line]
                                out.write("\n%s\n" % endstr)
                            else:
                                out.write(line)
                        elif endstr in line:
                            writing = True
                        else:
                            buffer.append(line)
                    else:
                        if not writing:
                            # There wasn't a closing "event", so write buffer contents
                            out.writelines(buffer)
                    out.close()
                    f.close()
            os.remove(grokfile)
            os.rename(tempgrokfile, grokfile)
            # print('Given strings added to < %s > ' % grokfile)

    def replace_strgrok(self, oldstr='', newstr=''):
        """ Replace a specific string from an existent grok file
        Arguments:
            oldstr:
            newstr:
        :return:
        """
        grokfile, tempgrokfile = self.findgrok()

        fileObj = dm.Manipfiles(grokfile)
        oldstr_idx = fileObj.find_str(grokfile, oldstr, index=True)
        # oldstr_idx = str_infile(grokfile, oldstr, index=True)
        # replace = {'Initial head': 'Initial head from file'}
        if oldstr_idx is None:
            print('String < %s > not found in <%s>' % (oldstr, os.path.split(grokfile)[-1]))
        elif oldstr_idx is not None:

            InFileObj = open(grokfile)
            InFile = InFileObj.read()
            InFile_replaced = InFile.replace(oldstr, newstr)

            OutFile = open(tempgrokfile, 'w')
            OutFile.write(InFile_replaced)

            OutFile.close()
            InFileObj.close()

            os.remove(grokfile)
            os.rename(tempgrokfile, grokfile)

    def modelTitle(self, description=''):
        # Write in description of the model:
        if not isinstance(description, str):
            sys.exit('description must be a valid string')
        else:
            # fmt_string = [description, 'end title']
            fmt_string = [description]
            self.addstr2grok(startstr='start title', endstr='end title',
                             str2add=fmt_string)
        print('< Problem description section written ok >')

    def gridgen(self, gridmethod='uniform',
                length=None,
                no_elements=None,
                coordinates=None,
                grade=None,
                mesh_file=None,
                base=None,
                top=None,
                layers=None,
                totecplot=False):
        """ Fill in Grid generation section with the geometry provided. Methods allowed....
        Argument:
            gridmethod:     str, {uniform, grade, gb_uniform}. Default 'uniform'
        generate uniform blocks, Domain length and number of blocks in X, Y, Z
        generate blocks interactive

        Returns:

        """
        if gridmethod == 'uniform':
            fmt_string = ['generate uniform blocks']
            [fmt_string.append('%.3f %d' % (length[ii], no_elements[ii])) for ii in range(0, len(length))]
            fmt_string.append('end grid generation')

        elif gridmethod == 'grade':
            xtra_strings = ['grade x', 'grade y', 'grade z']
            fmt_string = ['generate blocks interactive']

            [fmt_string.append('%s\n%.2f %.2f %.2f %.2f %.2f' % (xtra_strings[ii],
                                                                 coordinates[ii, 0],
                                                                 coordinates[ii, 1],
                                                                 grade[ii, 0],
                                                                 grade[ii, 1],
                                                                 grade[ii, 2])) for ii in
             range(0, coordinates.shape[0])]
            fmt_string.extend(['end generate blocks interactive', 'end'])

        elif gridmethod == 'gb_uniform':
            fmt_string = ['read gb 2d grid', mesh_file, 'generate layers interactive', '\tzone by layer', '']
            fmt_string.extend(['\tBase elevation', '\t\tElevation constant', '\t\t%f' % base, '\tend'])
            fmt_string.extend(
                ['\tnew layer', '\t\tlayer name\n\t\ttop layer', '\t\tElevation constant', '\t\t%f' % top])
            fmt_string.extend(['\t\tuniform sublayering', '\t\t%d' % layers, '\tend', 'end'])
            
        elif gridmethod == 'meshpy':
            fmt_string = ['read 3d grid', 'ascii', 'end grid definition']

        elif gridmethod == 'algomesh':
            fmt_string = ['read algomesh 2d grid', mesh_file, 'end grid definition']


        if totecplot is True:
            fmt_string.extend(['\nmesh to tecplot', 'meshtecplot.dat'])

        self.addstr2grok(startstr='!--- grid generation',
                         endstr='!--- end grid generation',
                         str2add=fmt_string)

        print('< Grid generation section written ok >')

    def simulation(self, units='kms', transient=False, transport=False, cvol=True, fdif=False):
        """

        Args:
            units:  str, units of mass(M), length(L) and time(T). k:kilogram, m:meter, c:centimeter, h:hour, d:day,
                    {kms, kmh, kmd, kmy, kcs, kcmin, kcmin, kch, kcd, kcy} Default 'kms': kilogram-metre-second
            fd:
            cv:
            transport:
            transient:
        Returns:

        """
        unit_dict = {'kms': 'kilogram-metre-second', 'kmh': 'kilogram-metre-hour',
                     'kmd': 'kilogram-metre-day', 'kmy': 'kilogram-metre-year', 'kcs': 'kilogram-centimetre-second',
                     'kcmin': 'kilogram-centimetre-minute', 'kch': 'kilogram-centimetre-hour',
                     'kcd': 'kilogram-centimetre-day', 'kcy': 'kilogram-centimetre-year'}
        var_dict = {'Control volume': cvol, 'Finite difference mode': fdif,
                    'Transient flow': transient, 'do transport': transport}

        if units in unit_dict:
            fmt_string = ['units: %s' % unit_dict[units]]
        else:
            sys.exit('wrong units')

        for key, value in var_dict.items():
            if value:
                fmt_string.append(key)

        self.addstr2grok(startstr='!--- general simulation parameters',
                         endstr='!--- end general simulation parameters',
                         str2add=fmt_string)

        print('< Simulation parameters added >')

    def mediaProps(self, mediaType='porous',
                   propsFile=None,
                   mediaName=None,
                   component=None,
                   typesel=None,
                   coords=None,
                   npts=None,
                   ptol=None,
                   kkkfile=None):
        """

        Returns:

        """
        fmt_string = ['use domain type']

        if mediaType == 'porous':
            startstr, endstr = '!--- porous media properties', '!--- end porous media properties'
            fmt_string.append('porous media')
        elif mediaType == 'overland':
            startstr, endstr = '!--- overland properties', '!--- end overland properties'
            fmt_string.append('surface')
        elif mediaType == 'well':
            startstr, endstr = '!--- well properties', '!--- end well properties'
            fmt_string.append('well')
        else:
            sys.exit('Wrong value for -mediaType-: choose porous, overland, or well')

        # Clear whatever is selected:
        fmt_string.extend(['properties file', propsFile, '\n'])

        # Select the given mesh components:
        aux_string = self.selComponents(component=component, selection=typesel, npts=npts, coords=coords, ptol=ptol)
        for ii in aux_string:
            fmt_string.append(ii)

        # Read properties (with proper name) from properties file:
        fmt_string.extend(['read properties', mediaName])

        # Add elemental file if provided:
        if kkkfile:
            fmt_string.extend(['Read elemental k from file', kkkfile])

        if (typesel != 'element') and (typesel == 'node') and (typesel == 'zone') and \
                (typesel == 'face') and (typesel == 'segment'):
            sys.exit('Wrong value for -typesel-: choose element, node, face, segment or zone')

        self.addstr2grok(startstr=startstr, endstr=endstr, str2add=fmt_string)

        print('< %s properties added >' % mediaType)

    def obs_pts(self, pt_name=None, pt_coord=None, well_name=None, well_coord=None):
        """ Write a list of observation wells/points to output information
            Args:
        pt_name:        array_like, string array with point names. Optional
        pt_coord:       array like, x,y,z coordinates of each obs point in pt_name. Optional
        well_name:      array_like, string array with well names. Optional
        well_coord:     array like, x1,y1,z1,x2,y2,z2 coordinates of each obs well in well_name. Optional
            Returns:
        Updated grok file with well/points included the proper section
        """
        fmt_string = []
        startstr, endstr = '!--- output locations', '!--- end output locations'

        if well_coord.ndim == 1:
            well_coord = np.expand_dims(well_coord, 0)
        if pt_name is not None:
            for ii in range(0, len(pt_name)):
                fmt_string.append('make observation point')
                fmt_string.append('%s\n%f  %f  %f' % (pt_name[ii], pt_coord[ii, 0], pt_coord[ii, 1], pt_coord[ii, 2]))

        if well_name is not None:
            for ii in range(0, len(well_name)):
                fmt_string.append('make observation well')
                fmt_string.append('%s\n%f  %f  %f\n%f  %f  %f'
                                  % (well_name[ii], well_coord[ii, 0], well_coord[ii, 1], well_coord[ii, 2],
                                     well_coord[ii, 3], well_coord[ii, 4], well_coord[ii, 5]))

        self.addstr2grok(startstr=startstr, endstr=endstr, str2add=fmt_string)

    def output_times(self, times=None):
        """ Write output time list to a grok file:
            Args:
        times: np array_like      array with the output times
            Returns:
        Updated grok file with output times in the proper section
        """
        fmt_string = []
        if times is not None:
            startstr, endstr = '!--- output times', '!--- end output times'
            fmt_string.append('output times')
            for ii in range(0, len(times)):
                fmt_string.append(times[ii])
            fmt_string.append('end')
        self.addstr2grok(startstr=startstr, endstr=endstr, str2add=fmt_string)

    @staticmethod
    def selComponents(component=None, selection=None, npts=None, coords=None, ptol=None, deselect=False):

        aux_string = []

        if deselect is not False:
            assert type(deselect) == tuple, 'Unexpected value for deselect argument'
            if deselect[1] == 'all':
                deselect_elements = ['element', 'node', 'face', 'zone', 'segment']
            else:
                deselect_elements = deselect[1:]

            for ii in deselect_elements:
                aux_string.append('clear chosen %ss' % ii)

        if (component == 'node') and (selection is None):
            aux_string = ['choose %s' % component, '%f, %f, %f' % (coords[0], coords[1], coords[2])]
        elif (selection is not None) and (component is not None):
            aux_string = ['choose %ss %s' % (component, selection)]

        if selection is not None:

            if (selection == 'polyline') and (component == 'segment'):
                aux_string.append(npts)
                for ii in range(0, npts):
                    aux_string.append('%f, %f, %f' % (coords[ii, 0], coords[ii, 1], coords[ii, 2]))

            elif selection.endswith('plane'):
                aux_string.extend(['%f' % coords, ptol])

        return aux_string

    def ini_cond(self, type_ic='head', val_ini=None, mode='flow', domaintype=None):
        """

        Args:
            type_ic:
            val_ini:
            mode:
            domaintype:

        Returns:

        """
        fmt_string = ['use domain type', domaintype]
        if mode == 'flow':
            startstr, endstr = '!--- initial conditions flow', '!--- end initial conditions flow'
        elif mode == 'transport':
            startstr, endstr = '!--- initial conditions transport', '!--- end initial conditions transport'

        aux_str = self.selComponents(deselect=(True, 'all'))
        for ii in aux_str:
            fmt_string.append(ii)
        fmt_string.extend(self.selComponents(component='node', selection='all'))

        if type_ic == 'head':
            fmt_string.extend(['Initial head', val_ini])
        elif type_ic == 'conc':
            fmt_string.extend(['Initial concentration', val_ini])

        self.addstr2grok(startstr=startstr, endstr=endstr, str2add=fmt_string)

        print('< %s > initial conditioned added to grok' % type_ic)

    def solute(self, name=None, diffcoef=None):
        """ Add solute properties
            Args:
        name:       string, solute or tracer name
        diffcoef:   float, diffusion coefficient of the solute
            Returns:
        Updated grok file with a new solute definition
        """
        startstr, endstr = '!--- solute properties', '!--- end solute properties'
        fmt_string = ['solute', ' name', ' %s' % name, ' free-solution diffusion coefficient', ' %e' % diffcoef, 'end solute']
        self.addstr2grok(startstr=startstr, endstr=endstr, str2add=fmt_string)

        print('Solute < %s > added to grok ' % name)

    def bound_cond(self, domaintype=None, mode='flow',
                   type_bc='head',
                   component='node',
                   selection=None,
                   coords=None,
                   setname=None,
                   bc_time=None,
                   bc_value=None,
                   npts=None,
                   ptol=None):
        """        !--- Boundary conditions well
        !--- End Boundary conditions well

        Args:
            mode:       flow, transport or well
            type_bc:
            component:
            selection:
            setname:
            bc_time:
            bc_value:

        Returns:
        """
        fmt_string = []
        # First deselect all components
        if mode == 'flow':
            aux_string = '!--- end boundary conditions flow'
        elif mode == 'well':
            aux_string = '!--- end boundary conditions well'
        # self.add2grok(aux_string, string_to_add)
        fmt_string.extend(self.selComponents(deselect=(True, 'all')))
        if domaintype is not None:
            fmt_string.extend(['use domain type', domaintype])

        fmt_string.extend(self.selComponents(component=component, selection=selection, npts=npts, coords=coords, ptol=ptol))

        fmt_string.extend(['create %s set' % component,
                           setname,
                           'boundary condition',
                           '\ttype',
                           '\t%s' % type_bc])

        fmt_string.extend(['\tname', '\t%s' % setname, '\t%s set' % component, '\t%s' % setname])

        if (bc_time is not None) and (bc_value is not None):
            fmt_string.append('\ttime value table')
            for ii in range(0, len(bc_time)):
                fmt_string.append('\t%f %f' % (bc_time[ii], bc_value[ii]))
            fmt_string.extend(['\tend', 'end'])

        fmt_string = [str(i) for i in fmt_string]
        self.add2grok(aux_string, '\n'.join(fmt_string))

        print('Boundary condition < %s > added to grok < %s >' % (type_bc, self.grokname))


class Inputfile(object):

    def __init__(self, model_dir, filename=None, filetype=None):
        self.model_dir = model_dir
        self.filename = filename
        self.filetype = filetype

    def addMaterial(self, name=None, mode='append', **kwargs):
        """

        Args:
            name:
            mode:
            *kwargs:

        Returns:

        """

        myfile = os.path.join(self.model_dir, '%s.%s' % (self.filename, self.filetype))
        fmt_string = ['!--- material', name]

        if mode == 'append':
            openmode = 'a'
        if mode == 'new':
            openmode = 'w'
        # mainheaders = """!--- Material'
        #         !--- End Material
        #         """
        for cur_key, cur_val in kwargs.items():
            fmt_string.extend(['%s' % cur_key, cur_val])

        fmt_string.extend(['end', '!--- end material\n'])

        with open(myfile, openmode) as fileobj:
            for ii in fmt_string:
                fileobj.write(str(ii)+'\n')

        print('Material < %s > added to %s' % (name, myfile))


def writeinGrok(myInDir, myOutDir, myoutfile, probdescFile, grid_coordFile, gridmethod,
                Transient, Transport, Echo2Output, FD, CV, porousfile, kfile, heads_init, transport_init, heads_bc,
                TimeStepCtrl, well_propsFile, pumpwellsFile, tracer_injFile, outputimesFile, obswellsFile):
    """ Type in the main parts of a grok file
    Arguments:
        myInDir:        Str, directory where input files are located
        myOutDir:       Str, directory where GROK file is/will be located
        myoutfile:      Str, name assigned to the GROK file
        probdescFile:   Str, filename containing description of the project
        grid_coordFile: Str, filename containing grid coordinates
        gridmethod:     Str, method to generate grid-'HGS_uniform' 'HGS_grade' 'GB_file'
        Transient:      Boolean, run transient mode or not
        Transport:      Boolean, run transport mode or not
        Echo2Output:    Boolean, create echo output file or not
        FD:             Boolean, finite difference mode or not
        CV:             Boolean, control volume mode or not
        porousfile:     Str, porous media properties file
        kfile:          Str, KKK filename for homogeneous fields, if homogeneous type ''
        heads_init:     Tuple [mode, arg]:
                        mode:'output',arg: string file name
                        mode: 'onevalue', arg: number
                        mode: 'raster', arg: string file name
                        mode:'txt', arg: string file name
        transport_init: Tuple [mode, arg]
                        mode:'output',arg: string file name
                        mode: 'onevalue', arg: number.
                        If transport is False, just type ''
        heads_bc:       tuple [X1, X2], X corrdinate of the X plane to be selected to assign constant head equals to initial value.
                        tuple ['all', 'head equals initial']: all heads stay the same
                        (In the future should support Y and Z planes and possibility
                        to assign boundary conditions from GB files)
        TimeStepCtrl:   Boolean, et the default time step controlers or not
        well_propsFile: Str, filename containing well construction details
        pumpwellsFile:  Str, filename containing well coordinates, rates, etc.
        tracer_injFile: Str, filename containing tracr injection data, if Transport is False, type ''
        outputimesFile: Str, filename containing defined output times
        obswellsFile:   Str, filename containing observation well information

    Returns:
        A new grok file
    """
    # %% PROBLEM DESCRIPTION:
    # Read problem description content from file:
    with open(os.path.join(myInDir, probdescFile), "r") as myfile:
        banner = myfile.read()  # .replace('\n', '')
    # Generate the grok file:
    Ofile = open(os.path.join(myOutDir, myoutfile), 'w')  # Generate the file
    Ofile.write(banner % (datetime.datetime.now()))  # %(time.strftime("%d/%m/%Y"))
    Ofile.flush()
    del banner
    print('<Problem description> written ok...')

    # %% GRID GENERATION: gb or directly in HGS
    banner = '\n\n!---  Grid generation\n'
    Ofile.write('%s' % banner)
    del banner
    if gridmethod == 'HGS_grade':
        Ofile.write('generate blocks interactive\n')
        dim = ['x', 'y', 'z']

        try:
            coordinates = np.loadtxt(os.path.join(myInDir, grid_coordFile))
        except IOError('Unable to load instructions'):
            raise

        if len(dim) != len(coordinates):
            raise IndexError('number of lines in coordinate file must be equal to length of dim vector')

        for cur_dim in range(0, len(dim), 1):
            Ofile.write('grade %s\n' % dim[cur_dim])

            for cur_val in range(0, coordinates.shape[1], 1):
                Ofile.write('%.8e\t' % coordinates[cur_dim, cur_val])
            Ofile.write('\n')

        Ofile.write('end generate blocks interactive\nend\n')
        Ofile.write('Mesh to tecplot\nmeshtecplot.dat\n')

    Ofile.flush()
    print('<Grid generation> written ok...')

    # %% SIMULATION PARAMETERS GROK:
    units = 'units: kilogram-metre-second'
    banner = '\n\n!---  General simulation parameters\n'

    Ofile.write('%s' % banner)
    Ofile.write('%s\n' % units)

    if FD == True: Ofile.write('Finite difference mode\n')
    if CV == True: Ofile.write('Control volume\n')
    if Transient == True: Ofile.write('Transient flow\n')
    if Transport == True: Ofile.write('do Transport\n')
    if Echo2Output == True: Ofile.write('Echo to output\n')

    # %% POROUS MEDIA PROPERTIES:
    Ofile.write('\n!--- Porous media properties\n')
    Ofile.write('use domain type\nporous media\nproperties file\n')
    Ofile.write('./%s\n\nclear chosen zones\nchoose zones all\nread properties\nporous medium\n' % (porousfile))

    if len(kfile) > 0:
        Ofile.write('\nclear chosen elements\nchoose elements all\nRead elemental k from file\n./%s\n' % (kfile))

    if Transport == True:
        name = 'fluorescein'
        difCoef = 1.000e-10
        Ofile.write('\n!--- Solute properties\n')
        Ofile.write('solute\n name\n %s\n free-solution diffusion coefficient\n %.3e\nend solute\n' % (name, difCoef))

    Ofile.flush()
    print('<Simulation param and porous media written ok...')

    # %% INITIAL CONDITIONS:
    del banner
    banner = 'clear chosen nodes\nchoose nodes all\n'
    Ofile.write('\n!--- IC\n!--- Flow:\n%s' % banner)

    if heads_init[0] == 'output':
        try:
            Ofile.write('initial head from output file\n')
            Ofile.write('%s\n' % (heads_init[1]))
        except TypeError('Innapropiate argument type in heads'):
            raise
    elif heads_init[0] == 'onevalue':
        Ofile.write('Initial head\n')
        try:
            Ofile.write('%.3e\n' % (heads_init[1]))
        except TypeError('Innapropiate argument type in heads'):
            raise
    elif heads_init[0] == 'raster':
        try:
            Ofile.write('Map initial head from raster\n')
            Ofile.write('%s\n' % (heads_init[1]))
        except TypeError('Innapropiate argument type in heads'):
            raise
    elif heads_init[0] == 'txt':
        try:
            Ofile.write('Initial head from file\n')
            Ofile.write('%s\n' % (heads_init[1]))
        except TypeError('Innapropiate argument type in heads'):
            raise

    if Transport == True:
        if transport_init[0] == 'onevalue':
            Ofile.write('\n!--- Transport:\n%s' % banner)
            Ofile.write('Initial concentration\n%.3e' % transport_init[1])

        elif transport_init[0] == 'output':
            try:
                Ofile.write('Initial concentration from output file\n')
                Ofile.write('%s\n' % (transport_init[1]))
            except TypeError('Innapropiate argument type in transport_init'):
                raise

    Ofile.flush()
    print('<Initial conditions written ok...')

    # %% BOUNDARY CONDITIONS:
    del banner
    banner = 'clear chosen nodes\nuse domain type\nporous media\n'
    Ofile.write('\n\n!--- BC\n!--- Flow:\n%s' % banner)

    if heads_bc[0] == 'all':
        Ofile.write('\nchoose nodes all\n')
    else:
        Ofile.write('\nchoose nodes x plane\n%.8e\n1.e-5\n' % (heads_bc[0]))
        Ofile.write('\nchoose nodes x plane\n%.8e\n1.e-5\n' % (heads_bc[1]))

    Ofile.write('\ncreate node set\nconstant_head_boundary\n')
    Ofile.write('\nBoundary condition\n type\n head equals initial\n node set\n constant_head_boundary\nEnd\n')

    Ofile.flush()
    print('<Boundary conditions written ok...')

    # %% TIME STEP CONTROLS:
    if TimeStepCtrl == True:
        flowsolver_conv_ = 1.0e-6
        flowsolver_max_iter = 500
        # head_control = 0.01
        initial_timestep = 1.0e-5
        min_timestep = 1.0e-9
        max_timestep = 100.0
        max_timestep_mult = 2.0
        min_timestep_mult = 1.0e-2
        underrelax_factor = 0.5

        Ofile.write('\n!--- Timestep controls\n')
        Ofile.write('flow solver convergence criteria\n%.2e\n' % (flowsolver_conv_))
        Ofile.write('flow solver maximum iterations\n%d\n' % (flowsolver_max_iter))
        # Ofile.write('head control\n%.2e\n'%(head_control))
        Ofile.write('initial timestep\n%.2e\n' % (initial_timestep))
        Ofile.write('minimum timestep\n%.2e\n' % (min_timestep))
        Ofile.write('maximum timestep\n%.2e\n' % (max_timestep))
        Ofile.write('maximum timestep multiplier\n%.2e\n' % (max_timestep_mult))
        Ofile.write('minimum timestep multiplier\n%.2e\n' % (min_timestep_mult))
        Ofile.write('underrelaxation factor\n%.2f\n' % (underrelax_factor))
        Ofile.write('No runtime debug\n\n')
        Ofile.flush()
        print('<Default time step controllers written ok...')

    # %% WELL DEFINITION:

    if pumpwellsFile:
        Ofile.write('\n!--- Pumping wells\n\n')
        try:
            wells_data = np.genfromtxt(os.path.join(myInDir, pumpwellsFile), comments="#", delimiter='',
                                       dtype='|S12, f8, f8, f4, f4, f4, f8, f8, f8, f8')
            len(wells_data)
        except (IOError, TypeError):
            print('There was an error with the <wells> file')
            print('File must contain at least two lines, if only one entrance, just duplicate it')
            Ofile.flush()
            Ofile.close()
            raise

        if wells_data[0] == wells_data[1]:
            length_wellsdata = 1
        else:
            length_wellsdata = len(wells_data)

        for ii in range(0, len(wells_data)):
            wellID = wells_data[ii][0]
            X = wells_data[ii][1]
            Y = wells_data[ii][2]
            top = wells_data[ii][3]
            bot = wells_data[ii][4]
            ext_point = wells_data[ii][5]
            pumprate_ini = wells_data[ii][6]
            pumprate_end = wells_data[ii][7]
            time_ini = wells_data[ii][8]
            time_end = wells_data[ii][9]

            Ofile.write('!-Well: %s\nuse domain type\nwell\nproperties file\n%s\n\n' % (wellID, well_propsFile))
            Ofile.write('clear chosen zones\nclear chosen segments\nchoose segments polyline\n2\n')
            Ofile.write('%s, %s, %s ! xyz of top of well\n' % (X, Y, top))
            Ofile.write('%s, %s, %s ! xyz of bottom of well\n' % (X, Y, bot))
            Ofile.write('new zone\n%d\n\n' % (ii + 1))
            Ofile.write('clear chosen zones\nchoose zone number\n%d\n' % (ii + 1))
            Ofile.write(
                'read properties\ntheis well\nclear chosen nodes\nchoose node\n%s, %s, %s\n' % (X, Y, ext_point))
            Ofile.write(
                'create node set\n %s\n\nboundary condition\n  type\n  flux nodal\n\n  name\n  %s\n  node set\n  %s\n' % (
                    wellID, wellID, wellID))
            Ofile.write('\n  time value table\n  %f   %f\n  %f   %f\n end\nend\n\n' % (
                time_ini, pumprate_ini, time_end, pumprate_end))

        Ofile.flush()
        print('<Pumping wells written ok...')
    if not pumpwellsFile:
        print('<Pumping wells not included in grok>')

    # %% TRACER INJECTION:
    if Transport == True:
        Ofile.write('\n!--- tracer injection\n\n')

        try:
            tracer_data = np.genfromtxt(os.path.join(myInDir, tracer_injFile), comments="#", delimiter='',
                                        dtype='f8, f8, f4, f4, f4, f4')
            len(tracer_data)
        except:
            print('File must contain at least two lines, if only one entrance, just duplicate it')
            print('There was an error with the Tracer Injection File')
            Ofile.flush()
            Ofile.close()
            raise
        if tracer_data[0] == tracer_data[1]:
            length_injdata = 1
        else:
            length_injdata = len(tracer_data)

        for ii in range(0, length_injdata):
            X = tracer_data[ii][0]
            Y = tracer_data[ii][1]
            Z = tracer_data[ii][2]
            time_ini = tracer_data[ii][3]
            time_end = tracer_data[ii][4]
            massflux = tracer_data[ii][5]

            Ofile.write(
                'clear chosen zones\nclear chosen segments\nclear chosen nodes\n\nuse domain type\nporous media\n')
            Ofile.write('\nchoose node\n')
            Ofile.write('%s, %s, %s\n' % (X, Y, Z))
            Ofile.write('\nspecified mass flux\n1\n%s,%s,%s\n\n' % (time_ini, time_end, massflux))

        Ofile.flush()
        print('<Tracer injection written ok...')

    # %% OUTPUT TIMES:
    if outputimesFile:
        try:
            times = np.loadtxt(os.path.join(myInDir, outputimesFile))
        except IOError:
            print('Check the directories of the output-times file')
            raise
        Ofile.write('!--- Outputs\noutput times\n')
        for ii in range(0, len(times)):
            Ofile.write('%s\n' % times[ii])

        Ofile.write('end\n' % times[ii])
        Ofile.flush()
        print('<Output times written ok...')
    if not outputimesFile:
        print('<No output times written in grok...>')

        # %% OBSERVATION WELLS:
    if obswellsFile:
        try:
            obs_wellsData = np.genfromtxt(os.path.join(myInDir, obswellsFile), comments="#", delimiter='',
                                          dtype='|S12, f8, f8, f4')
            len(obs_wellsData)
        except (IOError, TypeError, ValueError, OSError):
            print('Check the directories of the grok and output-times files')
            print('observation wells file must contain at least two lines, if only one entrance, just duplicate it')
            raise

        Ofile.write('\n!--- Observation wells\nuse domain type\nporous media\n\n')
        if obs_wellsData[0] == obs_wellsData[1]:
            length_obsData = 1
        else:
            length_obsData = len(obs_wellsData)

        for ii in range(0, length_obsData):
            ID = obs_wellsData[ii][0]
            X = obs_wellsData[ii][1]
            Y = obs_wellsData[ii][2]
            Z = obs_wellsData[ii][3]

            Ofile.write('make observation point\n')
            Ofile.write('%s\n%f  %f  %f\n' % (ID, X, Y, Z))
        Ofile.flush()
        print('<Observation wells written ok...')

    if not obswellsFile:
        print('<No observation wells written in grok...>')

    Ofile.close()
    print('GROK file finished!')



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
import grokfc.grok as grok
import subprocess as sp
# import shutil

import multiprocessing as mp
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


def fwdHGS(model_dir, modelID, hsplot=False):
    """ Run GROK and HGS and if True, run also hsplot
    Arguments:
        model_dir:      str, full path directory where GROK and HGS should be run
        modelID:        str, model ID
        hsplot:         bool, if false, hsplot is not executed. Default is False
    Returns:
    """
    grokObj = grok.Grokfile(model_dir, grokname=modelID)
    grokName, dummy = grokObj.findgrok()
    del dummy
    print('Model <%s>: Running GROK... ' % grokName[:-5])
    if platform.system() == 'Windows':
        qq = sp.Popen('grok', shell=True, cwd=model_dir, stdout=sp.PIPE, stderr=sp.PIPE)
    else:
        qq = sp.Popen('./grok', shell=True, cwd=model_dir, stdout=sp.PIPE, stderr=sp.PIPE)
    output_qq, error_qq = qq.communicate()
    if 'Normal exit' not in str(output_qq):
        myfile = open(os.path.join(model_dir, '_errorgrok.txt'), 'wb')
        myfile.write(output_qq)
        myfile.close()
        sys.exit('Something went wrong with GROK (Model: %s). Check _errorgrok.txt' % grokName[:-5])

    print('Model <%s>: Running HGS... ' % grokName[:-5])
    if platform.system() == 'Windows':
        qq_2 = sp.Popen('hgs', cwd=model_dir, stdout=sp.PIPE, stderr=sp.PIPE)
    else:
        qq_2 = sp.Popen('./hgs', cwd=model_dir, stdout=sp.PIPE, stderr=sp.PIPE)
    output_qq_2, error_qq_2 = qq_2.communicate()
    if 'NORMAL EXIT' not in str(output_qq_2):
        myfile = open(os.path.join(model_dir, '_errorhgs.txt'), 'wb')
        myfile.write(output_qq_2)
        myfile.close()
        if ('No more mass' in str(output_qq_2)) or ('Warning -1' in str(output_qq_2)):
            print('No more mass stored in the system: Model << %s >> ... ' % grokName[:-5])
            print('Check _errorhgs.txt')
        else:
            sys.exit('Something went wrong with HGS (Model: %s).Check _errorhgs.txt file' % grokName[:-5])

    if hsplot is True:
        print('Model <%s>: Running HSPLOT... ' % grokName[:-5])
        if platform.system() == 'Windows':
            qq_3 = sp.Popen('hsplot', cwd=model_dir, stdout=sp.PIPE, stderr=sp.PIPE)
        else:
            qq_3 = sp.Popen('./hsplot', cwd=model_dir, stdout=sp.PIPE, stderr=sp.PIPE)
        output_qq_3, error_qq_3 = qq_3.communicate()
        if 'Normal exit' not in str(output_qq_3):
            myfile = open(os.path.join(model_dir, '_errorhsplot.txt'), 'wb')
            myfile.write(output_qq_3)
            myfile.close()
            sys.exit('Something went wrong with HSPLOT (Model: %s). Check _errorgrok.txt' % grokName[:-5])
        if platform.system() == 'Windows':
            qq_3.kill()
    if platform.system() == 'Windows':
        qq.kill()
        qq_2.kill()
    print('Successful run. Model: %s' % grokName[:-5])

def fwdHGS_batch(procfolder, mymode, mydir, updateGrok, newtimes, batchfile, remove, files2keep):
    """
    Function designed to run HGS from a batchfile (e.g. "letmerun"). It is designed to be able to run in parallel and
    the directory structure for the EnKF...
    Args:
        procfolder:     str, folder containing all processes of the ensemble
        mymode:         str, 'fl_' for flow model or 'tr_' for transport
        mydir:          str, main directory (one level above procfolder)
        updateGrok:     bool, if it is wanted to update output times in each grok file
        newtimes:       np.array, with the output times
        batchfile:      str, name of the batch file to run the HGS model
        files2keep:
        remove:
    Returns:

    """

    modelfolder = os.listdir(os.path.join(mydir, procfolder))

    for idy, cur_modelfolder in enumerate(modelfolder):

        if mymode in cur_modelfolder:
            cur_dir = os.path.join(mydir, procfolder, cur_modelfolder)

            if updateGrok is True:
                mydict = {'fld_mode': cur_dir}
                addt2grok(mydict, None, newtimes, mode=mymode, mytype='full')
                print('Working with directory << %s >>' % os.path.split(mydir)[0])
            print('Running model %s. Executing batchfile %s' % (cur_modelfolder, batchfile))

            qq = sp.Popen(batchfile, shell=True, cwd=cur_dir, stdout=sp.PIPE, stderr=sp.PIPE)
            output_qq, error_qq = qq.communicate()
            if 'Normal exit' not in str(output_qq):
                myfile = open(os.path.join(cur_dir, '%s.txt' % batchfile), 'wb')
                myfile.write(output_qq)
                myfile.close()
                sys.exit('Something went wrong when executing %s in %s... Check %s.txt file' % (
                    batchfile, cur_modelfolder, batchfile))

        if remove is True:
            dummypaths = dm.rmfiles(files2keep, cur_modelfolder, os.path.join(mydir, procfolder))
        # dummypaths = dm.rmfiles(files2keep, '%s%.5d' % (mymode, zz + 1), os.path.join(allprocDir, curprocess))
            print('Removing files from folder < %s >' % cur_modelfolder)
            del dummypaths

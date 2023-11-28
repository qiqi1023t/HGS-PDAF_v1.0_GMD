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


class Manipdirs(object):
    """ Manipulate directories
    Arguments:
        mydir:          str, dir of interest
    """

    def __init__(self, mydir):
        self.my_dir = mydir

    def getdirs(self, mystr=None, fullpath=False, onlylast=False):
        """ List folders/ files in a defined directory, all wih a common string
            Arguments:
        mystr:          str, common string of files. Default is None
        fullpath:       bool, whether to attach full directory to file names or not
        onlylast:       bool, whether to return only the last value of the list or not
            Returns:
        mylist:         np.array type string with all folders or files of interest
        mydir:          list, directory explored
        """
        mydirs = self.my_dir
        if os.path.isdir(mydirs):
            mylist = os.listdir(mydirs)
        else:
            mydirs = str(input('Dir doesnt exist. Type a valid one: '))
            mylist = os.listdir(mydirs)  # Here I can use the function to constrain user input
        new_list = []
        mylist.sort()

        # %% If filter names according to a string:
        if mystr:
            for ii in mylist:
                if mystr in ii:
                    if fullpath is True:
                        new_list.append(os.path.join(mydirs, ii))
                    else:
                        new_list.append(ii)
            mylist = np.copy(np.array([new_list]))

        elif not mystr:
            if fullpath is True:
                for ii in mylist:
                    new_list.append(os.path.join(mydirs, ii))

                mylist = np.copy(np.array([new_list]))
        mylist.sort()
        new_list.sort()

        mylist = mylist.flatten()
        if onlylast is True:
            return mylist[-1], my_dirs
        elif not onlylast:
            return mylist, mydirs

    def rmfiles(self, keepfiles, mystr=False):
        """ Remove unnecessary files. Useful to save disk space
        Arguments:
            keepfiles:        np.array, strings to identify all files to be KEPT
            mystr:            str, common string in files/dirs. Default to False
            Returns:
            --------
        mypaths:        str, full path of file(s) of the process in execution, not deleted ????
        """
        # Read all folders within proper process number:
        mylist, dummy = self.getdirs(mystr=mystr, fullpath=True, onlylast=False)
        idx = []

        full_indices = np.arange(0, len(mylist), 1)

        for curname in keepfiles:

            mycurfile, dummy = self.getdirs(mystr=curname, fullpath=True, onlylast=False)

            if len(mycurfile) == 1:
                idx.append(np.where(mycurfile == mylist)[0])

            elif len(mycurfile) > 1:
                for ii in mycurfile:
                    idx.append(np.where(np.array(ii) == mylist)[0])

        toremove_indices = np.delete(full_indices, idx)

        try:
            for zz in toremove_indices:
                os.remove(mylist[zz])
        except (IOError, PermissionError, AssertionError):
            print('Not all unnecessary files deleted, remove manually...')
            pass


class Manipfiles(object):
    """ Manipulate content of single (and generic) ascii files.
    Arguments:
        myfile:         str, full path of the file of interest
    """

    def __init__(self, myfile):
        self.myfile = myfile

    def find_str(self, mystr, index=False):
        """ Find mystr + 1 character in the file. The string of interest should
        be able to be converted to "FLOAT" number(s). File format is any.
            Arguments:
        mystr:          str, string BEFORE the one of interest
        index:          bool, if True return the line index where the string was found
            Returns:
            --------
        String after the one given in the arguments. If index is True returns the line index
        """
        next_str = []
        ids = []
        with open(self.myfile, 'r') as fp:
            # counter = 0
            for idx, line in enumerate(fp):
                match = re.search(mystr + '([^,]+)', line)
                if match:
                    ids.append(idx)
                    next_str.append(match.group(1))
        if len(next_str) == 1:
            next_str = next_str[0]
            ids = ids[0]
        else:
            try:
                next_str = np.asarray(next_str, dtype='float')
                ids = np.asarray(ids, dtype='int')
            except:
                next_str = np.asarray(next_str, dtype='str')
                ids = np.asarray(ids, dtype='int')

        if len(next_str) > 0:
            if index is False:
                return next_str
            elif index is True:
                return ids, next_str
        else:
            print('Given string not found in file')


# def add_head2grok(model_dir, grokname='', hhfile=''):
#     # def add_head2grok(grokfc, headsFile1):
#     """
#
#     Args:
#         model_dir: full path to the grok file of interest
#
#     Returns:
#
#     """
#     grokfile, dummy = getdirs(model_dir, mystr='%s.grok' % grokname, fullpath=True)
#     tempgrokfile = os.path.join('%stemp' % grokfile)
#
#     h_str = str_infile(grokfile, 'Initial head from file', index=True)
#
#     if h_str is None:
#         replace = {'Initial head': 'Initial head from file'}
#         InFileObj = open(grokfile)
#         InFile = InFileObj.read()
#         for i in list(replace.keys()):
#             InFile = InFile.replace(i, replace[i])
#         OutFile = open(tempgrokfile, 'w')
#         OutFile.write(InFile)
#         OutFile.close()
#         InFileObj.close()
#         os.remove(grokfile)
#         OutFile = open(tempgrokfile)
#         strfound = str_infile(tempgrokfile, list([replace.values()][0])[0], index=True)
#         fp = open(grokfile, 'w')
#         for i, line in enumerate(OutFile):
#             if (i - 1) == strfound[0]:
#                 fp.write(hhfile)
#             else:
#                 fp.write(line)
#         fp.close()
#         OutFile.close()
#         os.remove(tempgrokfile)
#     else:
#         print('Initial head command was already updated......')  # todo: improve this print statement


# def kkk_bin2ascii(infile, outfile, orig_shape, biggrid_elid='', smallgrid_elid='', get_mean=True, exptransf=True,
#                   logtransf=False, expand=True):
#     """ Read binary files generated during inversion and creates a parameter file (kkk) that HGS can read
#     Args:
#         infile:
#         outfile:
#         orig_shape:
#         biggrid_elid
#         smallgrid_elid
#         get_mean:
#         exptransf:
#         logtransf:
#         expand:
#     Returns:
#     """
#     if type(infile) == str:
#         myparam = np.load(infile)
#     else:
#         myparam = infile
#
#     if get_mean is True:
#         myparam_mn = myparam.mean(axis=1)
#     elif get_mean is False:
#         myparam_mn = myparam
#
#     if exptransf is True:
#         myparam_mn = np.exp(myparam_mn)
#     if logtransf is True:
#         myparam_mn = np.log(myparam_mn)
#
#     if expand is True:
#         myparam_mn = gm.expshr_param(myparam_mn, biggrid_elid, smallgrid_elid, what='expand')
#
#     assert myparam_mn.shape[0] == orig_shape, 'shape of original and modified parameter arrays do NOT match!'
#     np.savetxt(outfile, np.transpose((np.arange(1, len(myparam_mn) + 1, 1), myparam_mn, myparam_mn, myparam_mn)),
#                fmt='%d %.18e %.18e %.18e')
#
#     print('File <%s> successfully stored' % outfile)


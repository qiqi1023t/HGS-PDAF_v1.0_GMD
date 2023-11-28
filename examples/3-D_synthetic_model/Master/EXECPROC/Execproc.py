import string,sys,math,os, time
from operator import itemgetter, attrgetter
import datetime
from math import *
import numpy as np

from grokfc import grok as gr
from grokfc import diradmin as da
from grokfc import runmodel as rm

start = time.time()

do_flow = True

################################
# LAUNCH ALLUVIAL FLOW MODEL...
################################

mod_dir = './HGS/'

if do_flow == True:
    print('LAUNCH ALLUVIAL MODEL FOR FLOW ONLY...')
    cls = lambda: os.system('cls')
    # Create the model folder if does not exists yet
    mygrok = 'Flow'
    batchfile = mod_dir + 'batch.pfx'
    batch = open(batchfile,'w')
    batch.write(mygrok)
    batch.close()
    time.sleep(2)
    # Run Grok and HGS:
    rm.fwdHGS(mod_dir, mygrok, hsplot=True)
    os.system('./flow_initials.sh')
    time.sleep(2)

elapsed = (time.time() - start)
print('ALLUVIAL MODEL DONE IN : ' + str(round(elapsed)/60.0) + ' minutes')
time.sleep(2)


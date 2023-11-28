# HGS-PDAF_v1.0_GMD
November 28,2023
This is the public repository for HGS-PDAF. The codes are developed based on the development made by Dr. Wolfgang Kurtz (DWD, Germany) and Prof. Dr. Oliver S. Schilling (Uni Basel, Switzerland). It is currently being developed and maintained by Dr. Qi Tang at the University of Neuchâtel, Switzerland.
## Introduction
HGS-PDAF is a modular DA framework for integrated hydrological modelling based on PDAF and HGS. It allows updating of states and parameters in integrated flow and transport simulations, assimilate different observation data, alone or jointly.
## Getting started with HGS-PDAF
### Download the source code of HGS-PDAF
- Use _git clone_ to checkout the source code from GitHub repository
- Use _git status_ to check if there is any update in the repository. This allows the user to get the latest version of the code. 
### Download the source code of PDAF
The PDAF source code can be downloaded from:
- The GitHub repository https://github.com/PDAF/PDAF
- The official PDAF website https://pdaf.awi.de/trac/wiki
After the PDAF source code is downloaded to your local directory, make a copy inside the HGS-PDAF directory named as _pdaf_. 
### Compiling HGS-PDAF
- Go to the _pdaf_ directory and follow the instructions on GitHub or the PDAF website to compile PDAF as a library. This is machine and system dependent.
- Go to the _hgsiolib_ directory to compile the HGS I/O subroutines which are to write files that are compatible with HGS. This can be simply by typing:
  _make_.
- Go to the _hgs-pdaf/offline_omi_ directory to compile _hgs-pdaf_ by typing _make_
- When everything is successfully compiled, you will find the executable _hgs-pdaf_ in the _hgs-pdaf/offline_omi_ directory.
## Test HGS-PDAF with the example test case
With this simplified, three dimensional synthetic model, the user can easily test HGS-PDAF. All the model configuration files are included in the _examples_ directory. The following steps help the user to prepare running this test case on a Linux machine:
- Go to the _examples/3-D_synthetic_model_ directory. Use _tar -xzvf prior_ensemble.tar.gz_ to extract the ensemble of 100 realisations of hydraulic conductivity fields from the archive file. The following two files will be extracted:
  - Prior_Ensemble.csv
  - Prior_Ensemble_homogeneous.csv
Prepare the initial ensemble of hydraulic conductivity fields based on these two .csv files. This can be done by simply modifying the python script _initialize.py_. For example, if the user wants to run a test simulation with ensemble size 2, this can be specified in line 9 as _Ensemble_size = 2_. One example script _init_100.py_ is geven for which an ensemble size of 100 is tested. After running this python script, a number of _n_ folders will be created, where _n_ is equal to the ensemble size.
- Go back to the _run_scripts_ directory and modify the run script _runscript_n2.job_ to adapt to the actual ensemble size. In the runs cript, the user needs to define the ensemble size, paths to the source code, model input, and model output. This is an example run script for the _slurm_ system. After the run script is prepared, the user can run this example case. Results will be written into the ouput directory. 

A detailed description on different subroutines and modules of HGS-PDAF as well as the preparation for the example test case can be found in the _docs_ directory. 

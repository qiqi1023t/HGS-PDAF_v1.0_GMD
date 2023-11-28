# HGS-PDAF_v1.0_GMD
November 28,2023
This is the public reporsitory for HGS-PDAF. The codes are developed based on the development made by Dr. Wolfgang Kurtz (DWD, Germany) and Prof. Dr. Oliver S. Schilling (Uni Basel, Switzerland). It is currently being developed and maintained by Dr. Qi Tang at the University of Neuchâtel, Switzerland.
## Introduction
HGS-PDAF is a modular DA framework for integrated hydrological modelling based on PDAF and HGS. It allows updating of states and parameters in integrated flow and transport simulations, assimilate different observation data, alone or jointly.
## Getting started with HGS-PDAF
### Downland the source code of HGS-PDAF
- Use _git clone_ to checkout the source code from GitHub reporsitory
- Use _git status_ to check if there is any update in the reporsitory. This allows the user to get the latest version of the code. 
### Downland the source code of PDAF
The PDAF source code can be downloaded from:
- The GitHub reporsitory https://github.com/PDAF/PDAF
- The official PDAF website https://pdaf.awi.de/trac/wiki
After the PDAF source code is downloaded to your local directory, make a copy inside the HGS-PDAF directory named as _pdaf_. 
### Compiling HGS-PDAF
- Go to the _pdaf_ directory and follow the instructions on GitHub or the PDAF website to compile PDAF as a library. This is maschine and system dependent.
- Go to the _hgsiolib_ directory to compile the HGS I/O subroutines which are to write files that are compatible with HGS. This can be simply by typing:
  _make_.
- Go to the _hgs-pdaf/offline_omi_ directory to compile _hgs-pdaf_ by typing _make_
- When everything is successfully compiled, you will find the executable _hgs-pdaf_ in the _hgs-pdaf/offline_omi_ directory.
## Test HGS-PDAF with the example test case
 


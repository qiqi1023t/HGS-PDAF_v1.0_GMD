# $Id: Makefile for HGS-PDAF 2021-09-14 OS & WK $

#######################################################
# Makefile to build PDAF-HGS based on offline model   #
# Architecture is updated macos_gfortran.h            #
#######################################################

######################################################

# User specifications
# 1. Set BASEDIR, the directory where the PDAF package resides
# 2. Set PDAF_ARCH to include compile definitions
# (See directory BASEDIR/pdaflib/make.arch for files. PDAF_ARCH is filename without .h)

# Root directory of PDAF package

BASEDIR = ../..

# architecture to use

#PDAF_ARCH = macos_gfortran


# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.


# include $(BASEDIR)/pdaflib/make.arch/macos_gfortran.h
include $(BASEDIR)/pdaf/make.arch/$(PDAF_ARCH).h

# Name of executable
EXE = hgs-pdaf

# End of user specifications
######################################################

.SUFFIXES: .F90 .o .f90

######################################################

# Modules used for the model part
MODULES = mod_parallel_pdaf.o \
		timer.o \
		mod_memcount.o \
		parser_mpi.o \

# Module required for assimilation
MOD_ASSIM = mod_assim_pdaf.o \
            mod_assim_hgs_pdaf.o

# Modules required for HGS IO lib
# --> done in separate makefile in /hgsiolib/

# Model routines used with PDAF
OBJ_MODEL = main_offline.o \
		initialize.o \

# Routines of observation handling (PDAF-OMI)
OBJ_USER_PDAFOMI = obs_HEAD_pdafomi.o \
		   obs_SAT_pdafomi.o \
		   callback_obs_pdafomi.o \


# Interface to PDAF - model sided
OBJ_PDAF_INT = init_parallel_pdaf.o \
	       read_config_pdaf.o \
		init_pdaf_offline.o \
		init_pdaf_info.o \
		assimilation_pdaf_offline.o \
		finalize_pdaf.o \
        
# Generic user-supplied routines
OBJ_USER_GEN = init_ens_offline.o \
		collect_state_pdaf_offline.o \
		output_netcdf_pdaf.o \
		prepoststep_ens_offline.o \


## User-supplied routines for SEIK
#OBJ_USER_SEEKSEIK = prodrinva_pdaf.o \
#		init_obsvar_pdaf.o \
#
## User-supplied routines for EnKF
#OBJ_USER_ENKF = add_obs_error_pdaf.o \
#		init_obscovar_pdaf.o
#
## User-supplied routines for localized analysis (LSEIK)
#OBJ_USER_LOCAL = init_n_domains_pdaf.o \
#		init_dim_l_pdaf.o \
#		g2l_state_pdaf.o \
#		l2g_state_pdaf.o \
#		prodrinva_l_pdaf.o \
#		init_obs_l_pdaf.o \
#		init_dim_obs_l_pdaf.o \
#		g2l_obs_pdaf.o \
#		init_obsvar_l_pdaf.o \
#		init_dim_obs_f_pdaf.o \
#		init_obs_f_pdaf.o \
#		obs_op_f_pdaf.o \
#
## User-supplied routines for LEnKF
#OBJ_USER_LENKF = localize_covar_pdaf.o
#
## User-supplied routines for NETF
#OBJ_USER_NETF = likelihood_pdaf.o
#
## User-supplied routines for localized NETF
#OBJ_USER_LNETF = likelihood_l_pdaf.o

# User-supplied routines for localized analysis
OBJ_USER_LOCAL = init_n_domains_pdaf.o \
		init_dim_l_pdaf.o \
		g2l_state_pdaf.o \
		l2g_state_pdaf.o \

# Full list of user-supplied routines for online modes

#OBJ_PDAF_USER = $(OBJ_USER_PDAFOMI) $(OBJ_USER_GEN) $(OBJ_USER_SEEK) $(OBJ_USER_ENS) \
#	$(OBJ_USER_SEEKSEIK) $(OBJ_USER_ENKF) $(OBJ_USER_LOCAL) \
#	$(OBJ_USER_LENKF) $(OBJ_USER_LNETF) $(OBJ_USER_NETF)

OBJ_PDAF_USER = $(OBJ_USER_PDAFOMI) $(OBJ_USER_GEN) $(OBJ_USER_LOCAL)

######################################################

$(EXE) : libpdaf-d.a libhgs.a \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_USER) $(OBJ_PDAF_INT) $(OBJ_MPI)
	$(LD)  $(OPT_LNK) -L$(BASEDIR)/hgsiolib/  -o $@  \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_USER) $(OBJ_PDAF_INT) $(OBJ_MPI) \
	-L$(BASEDIR)/pdaf/lib -L$(BASEDIR)/hgsiolib/ -lpdaf-d -lhgs \
	$(NC_LIB) $(LINK_LIBS)
	@echo "++++++ Done ++++++"

######################################################

libpdaf-d.a:
	@echo "++++++ Generate Filter library ++++++"
	@cd $(BASEDIR)/pdaf/src; make;

######################################################
libhgs.a:
	@echo "++++++ Generate HGS library ++++++"
	@cd $(BASEDIR)/hgsiolib; make;

######################################################
# now compile pdaflib and libhgs together

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS) $(MPI_INC) -c -I$(BASEDIR)/pdaf/include -I$(BASEDIR)/hgsiolib $*.F90
.f90.o :
	$(FC) $(OPT) $(CPP_DEFS) $(MPI_INC) -c -I$(BASEDIR)/pdaf/include -I$(BASEDIR)/hgsiolib $*.f90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

######################################################
# Cleans

cleanall : cleanpdaf clean cleanhgslib

clean :
	rm -f *.o *.mod *.a $(EXE)

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(BASEDIR)/pdaf/src; make clean

cleanhgs: cleanhgslib clean2

clean2 :
	rm -f *.o *.mod *.a $(EXE)
	
cleanhgslib:
	@echo "+++ Clean up HGS directory"
	cd $(BASEDIR)/hgsiolib; make clean

# cleandata:
#	rm -f ens*ana.txt state_ana.txt

######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 $(BASEDIR)/pdaf/make.arch | cut -d"." -f1

#######################################################
# Makefile to build HGS I/O lib for HGS-PDAF          #
# Based on updated macos_gfortran.h architecture      #
#######################################################

#FC      =gfortran

FC = mpif90
FCFLAGS = -g -Wall -O3 -r8 #-fdefault-real-8 #-fallow-argument-mismatch
LIBNAME = libhgs.a

main: hgs_mod.o hgs_fun.o
	ar rc $(LIBNAME) hgs_fun.o hgs_mod.o

hgs_mod.o:
	$(FC) -c $(FCFLAGS) hgs_mod.F90

hgs_fun.o:
	$(FC) -c $(FCFLAGS) hgs_fun.f90

clean:
	rm -f *mod *.o $(LIBNAME)


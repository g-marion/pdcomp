#---------------------#
# Makefile for pdcomp #
#---------------------#
#  Config variables   #

# Fortran compiler
F90=gfortran

# Locations of netcdf include files and libraries
NCINC=-I/usr/local/netcdf-4.8.1/include
NCLIB=-L/usr/local/netcdf-4.8.1/lib

# Locations of FFTW include files and libraries
FFTWINC=-I/usr/include # location of include files for fftw
FFTWLIB=-L/usr/local/lib64 # location of libraries for fftw

#################################################################################
# YOU PROBABLY DON'T WANT TO CHANGE ANYTHING BELOW THIS                         #
#################################################################################

OBFILES = getpp.o bc.o singleton.o pdcomp.o
OUTPUTOPT = -DNETCDF
LINKOPTS  = -lnetcdff

#################################################################################

pdcomp :        $(OBFILES)
		$(F90) -O3 -o ../run/pdcomp.exe $(OBFILES) -fcheck=bounds $(FFTWINC) $(FFTWLIB) -lfftw3 $(NCLIB) -lm $(LINKOPTS)

getpp.o :       getpp.f90
		$(F90) -c getpp.f90 -fcheck=bounds $(OPTS) $(NCINC) $(NCLIB)

bc.o :          bc.f90
		$(F90) -c bc.f90 -fcheck=bounds $(OPTS) $(NCINC) $(NCLIB)

singleton.o :   singleton.f90
		$(F90) -c singleton.f90 -fcheck=bounds $(OPTS) $(NCINC) $(NCLIB)

pdcomp.o : 	pdcomp.f90
		$(F90) -c pdcomp.f90 -fcheck=bounds $(OPTS) $(NCINC) $(NCLIB)

clean:
		rm -rf *o *mod

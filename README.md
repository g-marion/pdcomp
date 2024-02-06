#---------------------------#
# README for pdcomp program #
#---------------------------#

Contents:
1. GENERAL DESCRIPTION
2. IMPORTANT NOTES
3. BEFORE YOU RUN
4. SPECIFIC TO BIGBANG4 USERS

1. GENERAL DESCRIPTION
This program was originally created by George Bryan (NCAR) and modified by
Jeff Trapp (UIUC) and Geoff Marion (CIWRO). Similar code is implemented 
in CM1, but it is not set up for runs using distributed memory 
(the majority of CM1 simulations). This program is intended to bridge 
until that code is better-integrated into CM1, though could be used 
for post-processing if preferred.

Some modifications have been made to the original program for ease
of use (largely inspired by or taken from a similar program by Ryan 
Hastings), including some reconfiguring of the program to use an input
namelist, adding some output variables for pressure perturbations
and horizontal forcings, and adding a Makefile to make compiling the
program a bit easier.

A simple bash script has been added to the run directory to run the 
program for sequential CM1 output files, intended for running on the
head node of whatever machine you're using. Modify as necessary to
for batch job submission script, if the computer used supports that
use case.


2. IMPORTANT NOTES
- Program is currently set up to read/write netcdf output. Original
  program only read/write grads.
- DO NOT RUN PDCOMP OVER A SUBDOMAIN IN CM1. pdcomp uses information
  taken from horizontally averaging some variables at the top of the 
  analysis domain. If you run pdcomp for two subdomains of the same
  model run, they WILL NOT be comparable. It is especially bad to do this
  with the pdcomp domain not extending the full depth of the model domain.
  Just don't do it (unless you're prepared to make very extensive
  modifications of the pdcomp calculations, CM1 output reading subroutines,
  etc.)!
- If you want to run pdcomp for a model run that uses an unsupported
  (i.e., not Morrison or NSSL double-moment) microphysics scheme, it's 
  relatively easy to add that functionality. Simply change the if-statements
  that check if the input ptype is supported and add an additional
  if-statement where qtot is calculated (for calculating the buoyancy
  pressure) to use whatever mixing ratios your scheme uses.
- If you are using pdcomp on a model run that includes boundary layer
  turbulence from random potential temperature perturbations, recomputing 
  the base state variables (e.g., th0, p0) may be necessary in order to
  account for modifications to the base state by said turbulence. There is a 
  loop in getpp.f90 (immediately following 'Checkpoint 1') that is commented 
  out by default that does this by averaging the base state values over some 
  subdomain (hard set variables avgstart,avglen). It's recommended to do the 
  averaging over as large a subdomain as possible to ensure that it is 
  representative of this new base state.


3. BEFORE YOU RUN:
- Modify range of desired CM1 output times to perform analysis in run_pdcomp.bash
- Modify the filepaths in the def.pdcomp.input file and run_pdcomp.bash script
  to reflect your pdcomp and CM1 output directories.
- If not already done, make sure to set the LD_LIBRARY_PATH environment
  variable within your .bashrc (or whatever is relevant for the shell
  your machine is using).
- In the run directory, create a symbolic link to the netcdf.mod module file
  on your machine.


4. SPECIFIC TO BIGBANG4 USERS:
Use the following for LD_LIBRARY_PATH in .bashrc:
export LD_LIBRARY_PATH={$LD_LIBRARY_PATH}:/usr/local/netcdf-4.8.1/lib

Also, add the filepath to the directory for netcdf:
export NETCDF={$NETCDF}:/usr/local/netcdf-4.8.1

To create a symbolic link to the netcdf.mod file:
ln -s /usr/local/netcdf-4.8.1/include 

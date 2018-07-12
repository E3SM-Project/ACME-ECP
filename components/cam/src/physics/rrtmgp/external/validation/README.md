# Validation

This directory contains tools to validate RRTMGP against reference results. 'LBLRTM'/'CHARTS' provide the reference for clear-sky
fluxes with uniform boundary conditions (solar zenith angle, surface properties). 'SHDOMPP' provides the reference for more
general problems, including variable boundary conditions and clouds.

## Synopsis
### Compile  Code (From Parent RRTMGP Directory)
```
cd build; make; cd ../test; make; cd ../validation; make
```
### Run Clear Sky (From Validation Directory)
This will create PDF files comparing RRTMGP vs SHDOM, RRTMGP vs LBLRTM/CHARTS, and SHDOM vs  LBLRTM/CHARTS.
```
source run-val-clear.sh
```
Note: If the code has been compiled previously, consider executing
```
make clean
```
before the make command to ensure a pristine build.

### Run Cloudy Sky (From Validation Directory)
This will create PDF files comparing RRTMGP vs SHDOM for 12 cloudy cases in a common clear-sky profile.
```
source run-val-cloud.sh
```

## Script Details
'SHDOMPP' and 'RRTMGP' are run via the Python script '^/test/util/scripts/run\_tests.py'.
Script 'validation.py' creates PDF files comparing fluxes (bias, RMS) between two sets of files as specified in '.val' files.

Builds in this directory inherit compiler information from '^/build/Makefile.conf'.
Ensure that the RRMTMGP library and all the tools in '^/test' and '^/validation ' have been built.
run_all.sh will generate files with new boundary conditions and cloud optics, then run both SHDOMPP and
RRTMGP on those files.

Files with atmospheric conditions are in '^/test/data'.
Files in which gas optical properties have been computed from these conditions (using RRTMGP) are in '^/validation/data'.
Fluxes are also in '^/validation/data/'. SHDOM fluxes are computed using 8 polar angles.

Clear-sky reference files come from line-by-line calculations using uniform boundary conditions.
Longwave calculations come from 'LBLRTM'. Shortwave calculations from from 'LBLRTM + CHARTS'.
There are two sets of shortwave calculations, one with a black surface and one with uniform
surface albedo = 0.2.

Results for 'SHDOM' and 'RRTMGP' can be computed using the 'run-val-clear.sh' or 'run-val-cloud.sh'
scripts. These compute fluxes for the Garand atmosphere files in ../test/data. Differences between 'SHDOM' and 'LBLRTM' mostly
demonstrate the accuracy of the k-distribution.

The 'run-val-sw-black-sfc.sh' script is analogous to 'run-val-clear.sh' but uses a black surface in the shortwave.

The 'run-cmp-clear.sh' and 'run-cmp-cloud.sh' scripts compare RRTMGP solvers with ECRAD, the radiation code
developed by Robin Hogan at ECMWF. ECRAD requires its own license. Routines in '^/validation/ecrad/build/'
replace RRTMGP kernels with ECRAD kernels. An ECRAD-specific version of 'test_flux_compute_from_optics' can be
build in the same directory with 'make'. One of Robin's routines needs to be modified to accept variable
mu0.

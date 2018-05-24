#!/bin/bash
#
# Generate and compute fluxes for files with randomly-varying boundary conditions. 

# Ensure large enough stack (16 Mb) for reading files
# Ensure compatible Python environment (netCDF4, ...) 

# Create files with generalized boundary conditions starting from Garand atmosphere with uniform boundary conditions
python ../test/util/scripts/run_tests.py --root .. --no_diff --test make-general-clear-problems.ini

#
# Fluxes with SHDOM 
#
# Longwave clear-sky fluxes for uniform and variable BCs (this takes a few minutes)
python ../test/util/scripts/run_tests.py --root .. --no_diff --test run-shdompp-lw-gen-clear.ini
# Shortwave clear-sky fluxes for uniform and variable BCs (about 15 minutes)
python ../test/util/scripts/run_tests.py --root .. --no_diff --test run-shdompp-sw-gen-clear.ini

#
# Fluxes with RRTMGP 
#
# One-angle and three-angle clear-sky fluxes using flux_compute
python ../test/util/scripts/change-n-quad-angles.py 
python ../test/util/scripts/run_tests.py --root .. --no_diff --test run-rrtmgp-gen-clear.ini




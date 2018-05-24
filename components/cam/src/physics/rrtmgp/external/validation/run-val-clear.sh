#!/bin/bash
#
# Compare flux profiles for black surface, varying atmosphere, solar zenith angle of 30. degrees 

# Ensure large enough stack (16 Mb) for reading files
# Ensure compatible Python environment (netCDF4, ...) 

#
# Change input file to match SW calculations to match CHARTS 
#
python change-values.py ../test/data/garand-atmospheres-sw.nc sfc_alb_direct      0.2 
python change-values.py ../test/data/garand-atmospheres-sw.nc sfc_alb_diffuse     0.2 
python change-values.py ../test/data/garand-atmospheres-sw.nc solar_zenith_angle 30. 

#
# Fluxes with SHDOM 
#
# Longwave clear-sky fluxes (this takes a few minutes)
python ../test/util/scripts/run_tests.py --root .. --no_diff --test shdompp-lw-clear.ini
# Shortwave clear-sky fluxes (about 8 minutes)
python ../test/util/scripts/run_tests.py --root .. --no_diff --test shdompp-sw-clear.ini

#
# Fluxes with RRTMGP -- these start from atmospheric profiles, not the output of gas optics
#
# One-angle and three-angle clear-sky fluxes using flux_compute
python ../test/util/scripts/change-n-quad-angles.py 
python ../test/util/scripts/run_tests.py --root .. --no_diff --test rrtmgp-clear.ini

#
# Comparisons -- SHDOM/LBL, RRTMGP/LBL, RRTMGP/SHDOM, for both LW and SW 
#
cd data
ln -sf charts-sw-inputs-outputs-clear_alb_0.2_sza30.nc charts-sw-inputs-outputs-clear.nc
cd ..
python validation.py --val shdom-lblrtm-clear.val rrtmgp-lblrtm-clear.val rrtmgp-shdom-clear.val
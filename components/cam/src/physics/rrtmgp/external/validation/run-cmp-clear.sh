#!/bin/bash
#
# Compare flux profiles for grey surface, varying atmosphere, solar zenith angle of 30. degrees

# Ensure large enough stack (16 Mb) for reading files
# Ensure compatible Python environment (netCDF4, ...)

#
# Default files in RRTMGP tree have first index at surface; ECRAD requires opposite vertical ordering.
#

cp ../test/data/garand-atmospheres-lw.nc data/garand-atmospheres-lw-vertrev.nc
python reverse-vertical.py data/garand-atmospheres-lw-vertrev.nc
cp ../test/data/garand-atmospheres-sw.nc data/garand-atmospheres-sw-vertrev.nc
python reverse-vertical.py data/garand-atmospheres-sw-vertrev.nc

python change-values.py data/garand-atmospheres-sw-vertrev.nc sfc_alb_direct      0.2
python change-values.py data/garand-atmospheres-sw-vertrev.nc sfc_alb_diffuse     0.2
python change-values.py data/garand-atmospheres-sw-vertrev.nc solar_zenith_angle 30.


#
# Run RRTMGP and ECRAD on these reversed files
#
python ../test/util/scripts/run_tests.py --root .. --no_diff --test rrtmgp-ecrad-clear.ini

#
# Compare the files
#
python validation.py --val rrtmgp-ecrad-clear.val

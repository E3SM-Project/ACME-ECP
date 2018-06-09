#!/bin/bash

# Ensure large enough stack (16 Mb) for reading files
# Ensure compatible Python environment (netCDF4, ...)

#
# Frank's add_cloud_optics program has 12 sets of clouds. Test those with a single atmosphere. Last atmos won't have clouds in.
#

# A single Garand atmosphere, by default the first (Tropical?), repeated 13 times
#  (Column 7, counting from 0, has the lowest RMSE in the LW for flux up, dn)
#
python repeat-cols.py ../test/data/garand-atmospheres-lw.nc data/one-clear-atmos-x13-lw.nc --ncol 13
python repeat-cols.py ../test/data/garand-atmospheres-sw.nc data/one-clear-atmos-x13-sw.nc --ncol 13

#
# Add cloud optics - do this before reversing direction because I'm not sure how Frank's code
#   treats direction
#
if [ ! -e data/solver-lw-inputs-cloud.nc ]; then python ../test/run_tests.py --root .. --no_diff --test make-lw-cloud-problems.ini; fi
if [ ! -e data/solver-sw-inputs-cloud.nc ]; then python ../test/run_tests.py --root .. --no_diff --test make-sw-cloud-problems.ini; fi

# Longwave cloudy-sky fluxes (about 8 minutes)
python ../test/run_tests.py --root .. --no_diff --test shdompp-lw-cloud.ini
# Shortwave cloudy-sky  fluxes (about 15 minutes)
python ../test/run_tests.py --root .. --no_diff --test shdompp-sw-cloud.ini

#
# Fluxes with RRTMGP, relying on the cloud files built for SHDOMPP
#
# Cloudy-sky fluxes
python ../test/run_tests.py --root .. --no_diff --test rrtmgp-cloud.ini

#
# Validation plots
#
echo "Making longwave cloud plots"
python compare-clouds-by-profile.py \
  data/shdompp-lw-inputs-outputs-cloud.nc \
  data/rrtmgp-lw-inputs-outputs-cloud.nc \
  data/solver-one-clear-atmos-x13-lw.nc \
  rrtmgp-shdompp-lw-clouds-by-profile.pdf

echo "Making shortwave cloud plots"
python compare-clouds-by-profile.py \
  data/shdompp-sw-inputs-outputs-cloud.nc \
  data/rrtmgp-sw-inputs-outputs-cloud.nc \
  data/solver-one-clear-atmos-x13-sw.nc \
  rrtmgp-shdompp-sw-clouds-by-profile.pdf

#!/bin/bash

# Ensure large enough stack (16 Mb) for reading files
# Ensure compatible Python environment (netCDF4, ...)

#
# Frank's add_cloud_optics program has 12 sets of clouds. Test those with a single atmosphere. Last atmos won't have clouds in.
#

# A single Garand atmosphere, by default the first (Tropical?), repeated 12 times
#  (Column 7, counting from 0, has the lowest RMSE in the LW for flux up, dn)
#
python repeat-cols.py ../test/data/garand-atmospheres-lw.nc data/one-clear-atmos-x13-lw.nc --ncol 13
python repeat-cols.py ../test/data/garand-atmospheres-sw.nc data/one-clear-atmos-x13-sw.nc --ncol 13

#
# Add cloud optics - do this before reversing direction because I'm not sure how Frank's code
#   treats direction
#
if [ ! -e data/solver-lw-inputs-cloud.nc ]; then python ../test/util/scripts/run_tests.py --root .. --no_diff --test make-lw-cloud-problems.ini; fi
if [ ! -e data/solver-sw-inputs-cloud.nc ]; then python ../test/util/scripts/run_tests.py --root .. --no_diff --test make-sw-cloud-problems.ini; fi
#
# Default files in RRTMGP tree have first index at surface; ECRAD requires opposite vertical ordering.
#
cp data/solver-sw-inputs-cloud.nc data/solver-sw-inputs-cloud-vertrev.nc
cp data/solver-lw-inputs-cloud.nc data/solver-lw-inputs-cloud-vertrev.nc
python reverse-vertical.py data/solver-sw-inputs-cloud-vertrev.nc
python reverse-vertical.py data/solver-lw-inputs-cloud-vertrev.nc

#
# Need to have matching clear-sky results to compute cloud radiative effects
#
cp data/solver-one-clear-atmos-x13-lw.nc data/solver-one-clear-atmos-x13-lw-vertrev.nc
cp data/solver-one-clear-atmos-x13-sw.nc data/solver-one-clear-atmos-x13-sw-vertrev.nc
python reverse-vertical.py data/solver-one-clear-atmos-x13-lw-vertrev.nc
python reverse-vertical.py data/solver-one-clear-atmos-x13-sw-vertrev.nc

#
# Run RRTMGP and ECRAD on these reversed files
#
python ../test/util/scripts/run_tests.py --root .. --no_diff --test rrtmgp-ecrad-cloud.ini

#
# Compare the files
#
# python validation.py --val rrtmgp-ecrad-flux-cloud.val
echo "Making longwave cloud plots"
python compare-clouds-by-profile.py \
  data/ecrad-lw-inputs-outputs-cloud-vertrev.nc \
  data/rrtmgp-lw-inputs-outputs-cloud-vertrev.nc \
  data/ecrad-lw-inputs-outputs-clear-vertrev.nc \
  rrtmgp-ecrad-lw-clouds-by-profile.pdf

echo "Making shortwave cloud plots"
python compare-clouds-by-profile.py \
  data/ecrad-sw-inputs-outputs-cloud-vertrev.nc \
  data/rrtmgp-sw-inputs-outputs-cloud-vertrev.nc \
  data/ecrad-sw-inputs-outputs-clear-vertrev.nc \
  rrtmgp-ecrad-sw-clouds-by-profile.pdf

#!/bin/bash
# Compare RRTMGP and ECRAD bulk properties (R, T, source function) for clear-sky and scattering
#   atmosphere for LW and SW

############################
# Clear atmospheres - this part repeats run-cmp-clear.sh
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
############################
# Cloudy atmospheres- this part repeats run-cmp-clour.sh
#   Frank's add_cloud_optics program has 12 sets of clouds. Test those with a single atmosphere. Last atmos won't have clouds in.
#   A single Garand atmosphere, by default the first (Tropical?), repeated 12 times
#  (Column 7, counting from 0, has the lowest RMSE in the LW for flux up, dn)
#
python repeat-cols.py ../test/data/garand-atmospheres-lw.nc data/one-clear-atmos-x13-lw.nc --ncol 13
python repeat-cols.py ../test/data/garand-atmospheres-sw.nc data/one-clear-atmos-x13-sw.nc --ncol 13

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

python change-values.py data/solver-sw-inputs-cloud-vertrev.nc sfc_alb_direct      0.2
python change-values.py data/solver-sw-inputs-cloud-vertrev.nc sfc_alb_diffuse     0.2
python change-values.py data/solver-sw-inputs-cloud-vertrev.nc solar_zenith_angle 30.
############################
#
# Compute source functions and, for scattering atmospheres, the two-stream reflectance and
#   transmittance with ECRAD and RTE
#
../test/util/scripts/run_tests.py --root .. --test rrtmgp-ecrad-two-stream.ini --no_diff

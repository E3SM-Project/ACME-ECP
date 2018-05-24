#!/bin/bash
#
# Compare flux profiles for grey surface, varying atmosphere, solar zenith angle of 30. degrees

# Ensure large enough stack (16 Mb) for reading files
# Ensure compatible Python environment (netCDF4, ...)

root=`pwd`
cd ${root}/validation

# A single Garand atmosphere, by default the first (Tropical?), repeated 13 times
#  (Column 7, counting from 0, has the lowest RMSE in the LW for flux up, dn)
#
python repeat-cols.py ../test/data/garand-atmospheres-lw.nc data/one-clear-atmos-x13-lw.nc --ncol 13
python repeat-cols.py ../test/data/garand-atmospheres-sw.nc data/one-clear-atmos-x13-sw.nc --ncol 13

#
# Add cloud optics
#
if [ ! -e data/solver-lw-inputs-cloud.nc ]; then python ../test/util/scripts/run_tests.py --root .. --no_diff --test make-lw-cloud-problems.ini; fi
if [ ! -e data/solver-sw-inputs-cloud.nc ]; then python ../test/util/scripts/run_tests.py --root .. --no_diff --test make-sw-cloud-problems.ini; fi

for w in lw sw; do
  #
  # Compute fluxes
  #
  cd ${root}/test/flux_compute
  cp ${root}/validation/data/solver-${w}-inputs-cloud.nc ./rrtmgp-inputs-outputs.nc
  ./test_flux_compute_from_optics
  cp rrtmgp-inputs-outputs.nc rrtmgp-${w}-inputs-outputs-cloud.nc

  #
  # Now step-by-step
  #
  cd ${root}/test/two_stream
  cp ${root}/validation/data/solver-${w}-inputs-cloud.nc ./rrtmgp-inputs-outputs.nc
  ./test_two_stream
  cp rrtmgp-inputs-outputs.nc rrtmgp-${w}-inputs-outputs-cloud.nc

  cd ${root}/test/sources
  cp ${root}/test/two_stream/rrtmgp-${w}-inputs-outputs-cloud.nc ./rrtmgp-inputs-outputs.nc
  ./test_sources
  cp rrtmgp-inputs-outputs.nc rrtmgp-${w}-inputs-outputs-cloud.nc

  cd ${root}/test/adding
  cp ${root}/test/sources/rrtmgp-${w}-inputs-outputs-cloud.nc ./rrtmgp-inputs-outputs.nc
  ./test_adding
  cp rrtmgp-inputs-outputs.nc rrtmgp-${w}-inputs-outputs-cloud.nc

  cd ${root}/test/flux_reduce
  cp ${root}/test/adding/rrtmgp-${w}-inputs-outputs-cloud.nc ./rrtmgp-inputs-outputs.nc
  ./test_flux_reduce
  cp rrtmgp-inputs-outputs.nc rrtmgp-${w}-inputs-outputs-cloud.nc
done

cd ${root}

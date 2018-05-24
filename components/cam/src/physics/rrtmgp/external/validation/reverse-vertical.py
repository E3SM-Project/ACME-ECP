# make sure print behaves the same in 2.7 and 3.x
from __future__ import print_function

import os, sys

# package netCDF4 (https://github.com/Unidata/netcdf4-python)
import netCDF4 as nc
import numpy as np

# for reversing the vertical direction (global attribute in netCDF)
revTopAtt = 'top_at_1'


def path_check(path):
  """
  Quick check if path exists.  Use before reading a file.
  """
  if not os.path.exists(path):
    sys.exit('Could not find %s, returning' % path)
# end path_check()

def reverseVertical(inFile):
  """
  Reverse vertical dimension for all applicable variables in given
    netCDF file

  Input
    inFile -- string, path to netCDF to be modified

  Output
    nothing. inFile is overwritten

  Keywords
  """

  # open netCDF4 object and loop over variables in it
  ncObj = nc.Dataset(inFile, 'r+')
  ncVars = list(ncObj.variables)

  for ncVar in ncVars:
    inVar = ncObj.variables[ncVar]

    # these are just layer and level indices and should not be
    # reversed
    ll = ['lev', 'lay']
    if ncVar in ll: continue

    dims = inVar.dimensions; nDims = inVar.ndim

    # determine which axis to invert (either lev or lay nc dimension)
    for l in ll:
      if l in dims: axis = dims.index(l)
    # end l loop

    # is there a vertical dimension (i.e., has axis been assigned)?
    # if not, proceed to next array
    if not 'axis' in locals(): continue

    # not optimized for arrays with more than 3 dimensions
    if axis == 0:
      if nDims == 1:
        outVar = inVar[::-1]
      elif nDims == 2:
        outVar = inVar[::-1, :]
      elif nDims == 3:
        outVar = inVar[::-1, :, :]
      # endif nDims

    elif axis == 1:
      if nDims == 2:
        outVar = inVar[:, ::-1]
      elif nDims == 3:
        outVar = inVar[:, ::-1, :]
      # endif nDims

    elif axis == 2:
      outVar = inVar[:, :, ::-1]
    # end axis conditional

    ncObj.variables[ncVar][:] = outVar

    # so we don't carry axis to the next variable
    del(axis)
  # end loop over variables

  # These variables are referenced to the vertical ordering:
  #   "_inc" refers to increasing index along the vertical axis. So they
  #   need to be swapped.
  if('lev_src_inc' in ncVars and 'lev_src_dec' in ncVars):
      ncObj.renameVariable('lev_src_inc', 'temp')
      ncObj.renameVariable('lev_src_dec', 'lev_src_inc')
      ncObj.renameVariable('temp',        'lev_src_dec')

  # Change the global attribute
  if revTopAtt in ncObj.ncattrs():
    if ncObj.getncattr(revTopAtt) == 1: ncObj.setncattr(revTopAtt, 0)
    else: ncObj.setncattr(revTopAtt, 1)

  ncObj.close()

  return True
# end reverseVertical()

if __name__ == '__main__':

  import argparse

  parser = argparse.ArgumentParser(\
    description='Reverse file in lay/lev dimension')
  parser.add_argument('file', type=str, \
    help='Name of file.')

  args = parser.parse_args()
  path_check(args.file)
  reverseVertical(args.file)

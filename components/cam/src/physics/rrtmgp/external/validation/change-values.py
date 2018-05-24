from __future__ import print_function
# package netCDF4 (https://github.com/Unidata/netcdf4-python)
from netCDF4 import Dataset 
import os 
# ---------------------------------------------------------------------------------

def path_check(path):
  """
  Quick check if path exists.  Use before reading a file.
  """
  if not os.path.exists(path):
    sys.exit('Could not find %s, returning' % path)
# ---------------------------------------------------------------------------------


if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    description='Change values in a file.')
  parser.add_argument('file', type=str, \
    help='Name of file.')
  parser.add_argument('variable', type=str, \
    help='Name of variable.')
  parser.add_argument('value', type=float, default=0.0,\
    help='Value to set the variable to. Default 0.')

  args = parser.parse_args() 
  path_check(args.file) 
  f = Dataset(args.file, 'a')   
  
  f.variables[args.variable][:] = args.value
  f.close() 
  
  
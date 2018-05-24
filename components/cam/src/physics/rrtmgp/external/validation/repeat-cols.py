from __future__ import print_function
# package netCDF4 (https://github.com/Unidata/netcdf4-python)
import netCDF4 as nc
import numpy as np
import os 
import sys 
# ---------------------------------------------------------------------------------

def path_check(path):
  """
  Quick check if path exists.  Use before reading a file.
  """
  if not os.path.exists(path):
    sys.exit('Could not find %s, returning' % path)

# ---------------------------------------------------------------------------------
# Copy a variable definition including attributes from one file to another 
#
def copyVarDef(nc_in, nc_out, name, newname=None) :   
	if newname is None : 
		newname = name
	nc_out.createVariable(newname, nc_in.variables[name].dtype, nc_in.variables[name].dimensions) 
	nc_out.variables[newname].setncatts(nc_in.variables[name].__dict__)
	
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    description='Create testing files with constant atmosphere and a single variable change systematically between limits')
  parser.add_argument('input_file', type=str, \
    help='Name of input file.')
  parser.add_argument('output_file', type=str, \
    help='Name of output file.')
  parser.add_argument('--col_to_rep', type=int, default=0,\
    help='Which column (0-based) to repeat? Default 0.')
  parser.add_argument('--ncol', type=int, default=0,\
    help='Number of times to repeat the columns')
  parser.add_argument('--variables', type=str, nargs='+', \
    help='Name of any variable(s) whose values will change. All variables vary across the same range' + \
         'and are constant across a column')
  parser.add_argument('--minval', type=float, default=0.0,\
    help='Minimum value for variable. Default 0.')
  parser.add_argument('--maxval', type=float, default=1.0,\
    help='Minimum value for variable. Default 1.')

  ##############
  args = parser.parse_args() 
  # Open prototype input file 
  path_check(args.input_file) 
  ifile = nc.Dataset(args.input_file) 
  
  #
  # Error checking 
  #
  if args.col_to_rep > (ifile.dimensions['col'].size-1): 
      sys.exit('Asking to repeat column ' + str(args.col_to_rep) + \
               ' but max index is ' + str(ifile.dimensions['col'].size -1))
  
  #
  # Create new file with same dimensions and variables but different number of columns 
  #
  ofile = nc.Dataset(args.output_file, mode='w', FORMAT='NETCDF4_CLASSIC')  
  ncol = args.ncol
  if (ncol == 0): ncol = ifile.dimensions['col'].size
  d = ofile.createDimension('col',  ncol )
  for d in ifile.dimensions: 
    if not d == 'col': ofile.createDimension(d, ifile.dimensions[d].size) 

  #
  # For all variables: copy the definition to the new file with changed size 
  #
  for v in ifile.variables: copyVarDef(ifile, ofile, v) 
  #
  # ... and copy the specific column values
  #
  for v in ifile.variables: 
    if not 'col' in ifile.variables[v].dimensions: 
      ofile.variables[v][:] = ifile.variables[v][:]
    else: 
      # col is either the first or the last dimension 
      # Maybe there's an analog to Fortran SPREAD() but I don't fine one easily 
      if 'col' == ifile.variables[v].dimensions[0]:
        for i in range(0, ncol): ofile.variables[v][i, ...] =  ifile.variables[v][args.col_to_rep, ...]
      else: 
        for i in range(0, ncol): ofile.variables[v][..., i] =  ifile.variables[v][...,args.col_to_rep ]

  #
  # For the variables we want to change: vary the value across the 'col' dimension
  #
  if args.variables is not None: 
    for v in args.variables: 
      # Is the variable present in the file? 
      if not v in ifile.variables: sys.exit('Variable ' + v + ' not present in file ' + args.input_file)
      val = np.linspace(args.minval, args.maxval, ncol)
      # col is either the first or the last dimension 
      if 'col' == ifile.variables[v].dimensions[0]:
        for i in range(0, ncol): ofile.variables[v][i, ...] = val[i]
      else: 
        for i in range(0, ncol): ofile.variables[v][..., i] = val[i]

  ofile.close() 
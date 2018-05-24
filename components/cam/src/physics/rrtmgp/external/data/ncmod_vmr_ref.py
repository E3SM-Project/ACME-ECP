#!/usr/bin/env python

import os, sys, argparse
import numpy as np
import netCDF4 as nc

parser = argparse.ArgumentParser(\
  description='Edit some fields in a netCDF and redirect them ' + \
  'to a new netCDF.')
parser.add_argument('infile', type=str, help='netCDF to modify')
parser.add_argument('outfile', type=str, help='modified netCDF')
args = parser.parse_args()

inNC = args.infile; outNC = args.outfile
if not os.path.exists(inNC): sys.exit('Could not find %s' % inNC)

# open open the data and a file pointer for the output netCDF
inObj = nc.Dataset(inNC, 'r')
outObj = nc.Dataset(outNC, 'w')

# out has same dims as in
for dim in inObj.dimensions:
  if dim == 'key_absorber': continue
  dimObj = inObj.dimensions[dim]
  outObj.createDimension(dim, dimObj.size)
# end dim loop

# Copy global attributes
for name in inObj.ncattrs():
    outObj.setncattr(name, inObj.getncattr(name))

# except for one new dimension for outfile
outObj.createDimension(u'absorber_ext', \
  inObj.dimensions[u'absorber'].size+1)

# change dimensions of vmr_ref from
# (temperature x key_absorber x atmos_layer) to
# (temperature x absorber x atmos_layer)
for var in inObj.variables:
  # grab variable attributes and use them for new file
  outVar = inObj.variables[var]
  outType = outVar.dtype

  dat = np.array(outVar)

  if var == 'vmr_ref':
    dim1 = outVar.dimensions[0]
    dim3 = outVar.dimensions[2]
    outDim = (dim1, u'absorber_ext', dim3)

    # add 0.781 (per lblatm, N2) array to absorber_ext dimension
    absExt = np.zeros((outObj.dimensions[dim1].size, 1, \
      outObj.dimensions[dim3].size)) + 0.781
    dat = np.append(dat, absExt, axis=1)
  else:
    outDim = outVar.dimensions
  # endif var

  ncVar = outObj.createVariable(var, outType, outDim)
  outObj.variables[var].setncatts({k: inObj.variables[var].getncattr(k) for k in inObj.variables[var].ncattrs()})
  ncVar[:] = dat
# end var loop

# add in new value into vmr_ref

inObj.close()
outObj.close()

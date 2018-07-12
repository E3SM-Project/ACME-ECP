#
# Dirt-simple Python to see which fields differ between two files
#   Saving this not because we'll need it but because I've learned a tiny 
#   bit of Python that would be painful to reconstruct  
#
import numpy as np
from netCDF4 import Dataset  

def delta(ds1, ds2, varName) :
  dx = ds1.variables[varName][:] - ds2.variables[varName][:] 
  return (np.min(dx), np.max(dx))


def cmp2files(f1, f2):
  ds1 = Dataset(f1, mode = 'r') 
  ds2 = Dataset(f2, mode = 'r') 
  varnames = [v for v in ds1.variables]
  for v in varnames:
    mn,mx = delta(ds1, ds2, v)
    if(np.any(np.not_equal([mn,mx], 0.))) : 
      print(v)

cmp2files('../test/gas_optics/ref/rrtmgp-lw-inputs-outputs.nc', 'rrtmgp-lw-inputs-outputs-gen-clear.nc')
cmp2files('../test/gas_optics/ref/rrtmgp-lw-inputs-outputs.nc', 'rrtmgp-lw-inputs-outputs-cloud.nc')
cmp2files('../test/gas_optics/ref/rrtmgp-lw-inputs-outputs.nc', 'rrtmgp-lw-inputs-outputs-gen-cloud.nc')

cmp2files('../test/gas_optics/ref/rrtmgp-sw-inputs-outputs.nc', 'rrtmgp-sw-inputs-outputs-gen-clear.nc')
cmp2files('../test/gas_optics/ref/rrtmgp-sw-inputs-outputs.nc', 'rrtmgp-sw-inputs-outputs-cloud.nc')
cmp2files('../test/gas_optics/ref/rrtmgp-sw-inputs-outputs.nc', 'rrtmgp-sw-inputs-outputs-gen-cloud.nc')


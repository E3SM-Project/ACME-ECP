import netCDF4
import sys
import numpy as np

nc1 = netCDF4.Dataset(sys.argv[1])
nc2 = netCDF4.Dataset(sys.argv[2])

print("Var Name".ljust(20)+":  "+"Min 2-norm".ljust(20)+"  ,  "+"Avg 2-norm".ljust(20)+"  ,  "+"Max 2-norm")
for v in nc1.variables.keys() :
  if (nc2.variables[v].dtype == np.float64 or nc2.variables[v].dtype == np.float32) :
    diff2 = ( nc2.variables[v][:] - nc1.variables[v][:] )**2 / (nc1.variables[v][:]**2 + 1.e-14)
    norm2_avg = np.sum(diff2)/np.prod(diff2.shape)
    norm2_max = np.amax(diff2)
    norm2_min = np.amin(diff2)
    print(v.ljust(20)+":  %20.10e  ,  %20.10e  ,  %20.10e"%(norm2_min,norm2_avg,norm2_max))


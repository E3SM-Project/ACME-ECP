import netCDF4
import sys
import numpy as np

nc1 = netCDF4.Dataset(sys.argv[1])
nc2 = netCDF4.Dataset(sys.argv[2])

print("Var Name".ljust(20)+":  "+"rel 2-norm".ljust(20)+"  ,  "+"rel inf-norm".ljust(20)+"  ,  "+"max abs")
for v in nc1.variables.keys() :
  if (nc2.variables[v].dtype == np.float64 or nc2.variables[v].dtype == np.float32) :
    a1 = nc1.variables[v][:]
    a2 = nc2.variables[v][:]
    adiff = abs(a2-a1)

    norm2 = np.sum( adiff**2 )
    norm2_denom = np.sum( a1**2 )
    if (norm2_denom != 0) :
      norm2 = norm2 / norm2_denom

    normi = np.amax( adiff )
    normi_denom = np.amax(a1) - np.amin(a1)
    if (normi_denom != 0) :
      normi = normi / normi_denom

    maxabs = np.amax( adiff )

    print(v.ljust(20)+":  %20.10e  ,  %20.10e  ,  %20.10e"%(norm2,normi,maxabs))


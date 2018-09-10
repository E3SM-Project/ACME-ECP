import netCDF4
import sys
import numpy as np

nc1 = netCDF4.Dataset(sys.argv[1])
nc2 = netCDF4.Dataset(sys.argv[2])
nc3 = netCDF4.Dataset(sys.argv[3])

print("Var Name".ljust(20)+":  "+"|1-2|".ljust(20)+"  ,  "+"|2-3|".ljust(20)+"  ,  "+"|2-3|/|1-2|")
for v in nc1.variables.keys() :
  if (nc2.variables[v].dtype == np.float64 or nc2.variables[v].dtype == np.float32) :
    a1 = nc1.variables[v][:]
    a2 = nc2.variables[v][:]
    a3 = nc3.variables[v][:]
    adiff12 = abs(a2-a1)
    adiff23 = abs(a3-a2)

    norm12 = np.sum( adiff12**2 )
    norm23 = np.sum( adiff23**2 )
    norm_denom = np.sum( a1**2 )
    if (norm_denom != 0) :
      norm12 = norm12 / norm_denom
      norm23 = norm23 / norm_denom

    normRatio = norm23
    if (norm12 != 0) :
      normRatio = norm23 / norm12
    else :
      if (norm23 == 0) :
        normRatio = 1
      else :
        normRatio = 1e50

    if (normRatio > 2) :
      print(v.ljust(20)+":  %20.10e  ,  %20.10e  ,  %20.10e"%(norm12,norm23,norm23/norm12))


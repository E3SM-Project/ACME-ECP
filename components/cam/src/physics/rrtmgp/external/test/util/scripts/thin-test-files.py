from netCDF4 import Dataset
import os 

def copyVar(nc_in, nc_out, name) :   
	nc_out.createVariable(name, nc_in.variables[name].dtype, nc_in.variables[name].dimensions) 
	nc_out.variables[name].setncatts(nc_in.variables[name].__dict__)
	nc_out.variables[name][:] = nc_in.variables[name][:] 

oldFName = '../validation/data/lblrtm-lw-flux-inputs-outputs-garand-all-r674.nc'
newFName = 'garand-atmospheres-lw.nc'
oldF = Dataset(oldFName,  'r')
newF = Dataset(newFName, 'w') 
varsToDel = ['flux_up', 'flux_dn', 'flux_net', 'heating_rate', 'flux_dif_dn', 'flux_dir_dn', 'band_flux_up', 'band_flux_dn', 'band_flux_net', 'band_heating_rate', 'band_flux_dir_dn', 'band_flux_dif_dn', 'angle_secant', 'angle_weight']

#
# Copy dimensions
#
for d in oldF.dimensions: 
  newF.createDimension(oldF.dimensions[d].name, oldF.dimensions[d].size) 

for v in oldF.variables: 
  if not v in varsToDel: copyVar(oldF, newF, v)

newF.close()
oldF.close()
# os.rename(newFName, oldFName) 

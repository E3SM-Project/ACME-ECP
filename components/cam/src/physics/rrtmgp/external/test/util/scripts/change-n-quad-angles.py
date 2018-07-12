from netCDF4 import Dataset  
# ---------------------------------------------------------------------------------
# Copy a variable and all its attributes from one netCDF file to another 
#
def copyVar(nc_in, nc_out, name, newname=None) :   
	if newname is None : 
		newname = name
	nc_out.createVariable(newname, nc_in.variables[name].dtype, nc_in.variables[name].dimensions) 
	nc_out.variables[newname].setncatts(nc_in.variables[name].__dict__)
	nc_out.variables[newname][:] = nc_in.variables[name][:] 
# ---------------------------------------------------------------------------------


nangles = 3
# Script will get executed in ^/test or ^/validation
inFile = '../test/data/garand-atmospheres-lw.nc' 
outFile = 'data/garand-atmospheres-lw-ang-' + str(nangles) + '.nc'

inDS  = Dataset(inFile) 
outDS = Dataset(outFile, 'w') 

for attname in inDS.ncattrs():
  setattr(outDS,attname,getattr(inDS,attname))

# Copy dimensions 
for dimname,dim in inDS.dimensions.iteritems():
  if not dimname == 'angle' : 
    outDS.createDimension(dimname,len(dim))

# Copy variables
for varname,ncvar in inDS.variables.iteritems():
  if not varname == 'angle': copyVar(inDS, outDS, varname) 

# Number of angles 
outDS.createDimension('angle',nangles)
v = outDS.createVariable('angle', 'i4', ['angle']) 
v[:] = range(1, nangles + 1) 

inDS.close()
outDS.close()

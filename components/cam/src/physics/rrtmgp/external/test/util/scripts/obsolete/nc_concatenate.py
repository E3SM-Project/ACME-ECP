# make sure print behaves the same in 2.7 and 3.x
from __future__ import print_function

# package netCDF4
# https://github.com/Unidata/netcdf4-python
import netCDF4

import sys
import numpy as np
import argparse

# setup: CSCS environment
# module switch PrgEnv-cray PrgEnv-gnu
# module load Python/2.7.10-CrayGNU-5.2.40
# module load cray-netcdf
# virtualenv-2.7 --system-site-packages --python=$( which python2.7 ) venv
# source venv/bin/activate
# easy_install pip==1.2.1
# pip install netCDF4
# python -c 'from netCDF4 import *'

var_concatenate = {}


def main(outfile, files, dimension_name):
    from netCDF4 import Dataset
    import numpy as np

    print('concartinating over dimension:', dimension_name)
    print('creating output file:', outfile)
    dsout = Dataset(outfile, "w", format="NETCDF4", zlib=True)

    # for each input file
    for filename_in in files:
        # open input file
        print('reading input file:', filename_in)
        dsin = netCDF4.Dataset(filename_in, 'r')

        # Copy dimensions
        for dname, the_dim in dsin.dimensions.iteritems():
            if dname not in dsout.dimensions:
                dim_val = None
                if dname == dimension_name:
                    dim_val = len(files)
                elif the_dim.isunlimited():
                    dim_val = None
                else:
                    dim_val = len(the_dim)
                print('creating dimension:', dname, dim_val)
                dsout.createDimension(dname, dim_val)

        # Copy or concartenate variables
        for vname, varin in dsin.variables.iteritems():

            if dimension_name not in varin.dimensions:  # this variable should be copied
                # Copy variable
                if vname not in dsout.variables:  # only if not already there
                    print('copying variable:', vname)
                    varin = dsin.variables[vname]
                    varout = dsout.createVariable(
                        vname, varin.datatype, varin.dimensions)
                    varout.setncatts({k: varin.getncattr(k)
                                      for k in varin.ncattrs()})
                    varout[:] = varin[:]

            else:  # this variable should be concartenated
                # create variable
                if vname not in dsout.variables:  # only if not already there
                    print('creating variable:', vname)
                    varin = dsin.variables[vname]
                    varout = dsout.createVariable(
                        vname, varin.datatype, varin.dimensions)
                    varout.setncatts({k: varin.getncattr(k)
                                      for k in varin.ncattrs()})

                # concartenate variable
                varout = dsout.variables[vname]
                if vname not in var_concatenate:
                    var_concatenate[vname] = []
                var_concatenate[vname].append(varin)
                if len(var_concatenate[vname]) == len(files):
                    index = varout.dimensions.index(dimension_name)
                    if varin[:].shape[index] != 1:
                        raise ValueError('cannot concartinate: dimension ' +
                                         dimension_name + ' is not 1 for variable ' + vname)
                    print('concartinating variable', vname, 'on index', index)
                    varout[:] = np.concatenate(
                        var_concatenate[vname], axis=index)
        dsin.close()
    dsout.close()

if __name__ == '__main__':
    # command line arguments
    parser = argparse.ArgumentParser(
        description='Merge NetCDF input profiles.')
    parser.add_argument('-o', '--output', type=str, nargs=1, required=True,
                        help='output file')
    parser.add_argument('-d', '--dimension', type=str, nargs=1, default=[u'time'],
                        help='dimension to concartinate; default is "time"')
    parser.add_argument('file', nargs='+',
                        help='netCDF files')
    args = parser.parse_args()
    outfile = args.output[0]
    dimension_name = args.dimension[0]
    files = args.file
    # call main
    main(outfile, files, dimension_name)

# make sure print behaves the same in 2.7 and 3.x
from __future__ import print_function

import sys, os

# setup: CSCS environment
# module switch PrgEnv-cray PrgEnv-gnu
# module load Python/2.7.10-CrayGNU-5.2.40
# module load cray-netcdf
# virtualenv-2.7 --system-site-packages --python=$( which python2.7 ) venv
# source venv/bin/activate
# easy_install pip==1.2.1
# pip install netCDF4
# python -c 'from netCDF4 import *'


def diff2files(filename1, filename2, varname, digits, factor):
    # package netCDF4
    # https://github.com/Unidata/netcdf4-python
    import netCDF4
    import numpy as np

    # Reading file 1 (output data)
    file1 = netCDF4.Dataset(filename1)
    var1 = file1.variables[varname]
    shape1 = var1.shape
    #print(filename1, varname + str(shape1))
    data1 = var1[:]
    file1.close()

    # Reading file 2 (expected data)
    file2 = netCDF4.Dataset(filename2)
    var2 = file2.variables[varname]
    shape2 = var2.shape
    #print(filename2, varname + str(shape2))
    data2 = var2[:]
    file2.close()

    # shape the same?
    if shape1 != shape2:
        #print('different shapes')
        return 0

    # inf values present?
    if np.isinf(data1).any():
        """
        print('inf values found in ' + filename1 + ':' + varname)
        for r in np.argwhere(np.isinf(data1)):
            print(r)
        """
        return 0
    if np.isinf(data2).any():
        """
        print('inf values found in ' + filename2 + ':' + varname)
        for r in np.argwhere(np.isinf(data2)):
            print(r)
        """
        return 0

    # NaN values present?
    if np.isnan(data1).any():
        """
        print('NaN values found in ' + filename1 + ':' + varname)
        for r in np.argwhere(np.isnan(data1)):
            print(r)
        """
        return 0
    if np.isnan(data2).any():
        """
        print('NaN values found in ' + filename2 + ':' + varname)
        for r in np.argwhere(np.isnan(data2)):
            print(r)
        """
        return 0

    # compare by number of digits
    if digits:
        # compare element by element
        diff_count = 0
        itr1 = np.nditer(data1, flags=['multi_index'])
        itr2 = np.nditer(data2, flags=['multi_index'])
        while not itr1.finished:
            v1 = itr1[0]
            v2 = itr2[0]
            v12 = (v1 + v2) / 2.
            #numformat = "%."+str(digits)+"f"
            numformat = "%." + str(digits) + "e"
            s1 = numformat % abs(v1)
            s2 = numformat % abs(v2)
            s12 = numformat % abs(v12)
            r1 = s1.replace(".", "")[0:digits] + '|' + s1[len(s1) - 4:len(s1)]
            r2 = s2.replace(".", "")[0:digits] + '|' + s2[len(s2) - 4:len(s2)]
            r12 = s12.replace(".", "")[0:digits] + \
                '|' + s12[len(s12) - 4:len(s12)]
            if (np.isnan(v1) or np.isnan(v2) or ((r1 != r12) and (r2 != r12))):
                pos1 = pos2 = 0
                for x in range(digits + 1):
                    if (r1[x] != r12[x]):
                        pos1 = x
                        break
                for x in range(digits + 1):
                    if (r2[x] != r12[x]):
                        pos2 = x
                        break
                pos = max(pos1, pos2) + 1
                numformat2 = "%." + str(digits) + "e"
                #print('digit-diff:', pos, varname + str(itr1.multi_index),
                #      'values:', numformat2 % v1, numformat2 % v2)
                #print('digit-diff:',pos, varname+str(itr1.multi_index), 'values:', r1, r12, r2)
                diff_count += 1

            itr1.iternext()
            itr2.iternext()

        """
        if (diff_count > 0):
            print('total differences: ' + str(diff_count))
        else:
            print('identical ' + str(digits) + ' digits')
        """
        return diff_count

    # compare by relative difference
    if factor:
        # Compare output data to reference data.
        def myDiff(a, b):
            # both values are 0; they are identical
            if a == 0 and b == 0:
                return 0
            # a is 0, but b is non-zero
            if a == 0:
                return np.inf
            # b is 0, but a is non-zero
            if b == 0:
                return np.inf
            # compute the absolute relative error
            return abs(a / b - 1.0)
        vMyDiff = np.vectorize(myDiff)
        diff = vMyDiff(data1, data2)

        # RLP: 10-Oct-2016
        err = np.argwhere(diff >= factor)
        diff_count = len(err)
        nErr = err.size; nDiff = diff.size
        if diff_count > 0:
          nPoints = diff.size
          flatDiff = diff.flatten()
          diffIdx = np.where(flatDiff >= factor)[0]
          relDiffPer = diff_count/float(nPoints) * 100

          outStr1 = '%s: %.2f%% of values differ by more than %f; ' % \
            (varname, relDiffPer, factor)
          outStr2 = 'Relative Differences Range: [%e, %e]' % \
            (flatDiff[diffIdx].min()*100, flatDiff[diffIdx].max()*100)
          print(outStr1 + outStr2)
        # end nErr

        """
        for r in err:
            print(', '.join(['index:' + str(r), 'left value: ' + str(data1[tuple(r)]),
                             'right value: ' + str(data2[tuple(r)]), 'relative difference: ' + str(diff[tuple(r)])]))
             print(str(r)+str(data1[tuple(r)]))

        print('factor:', factor)
        if (diff_count > 0):
            print('total differences: ' + str(diff_count))
        else:
            print('identical')
        """
        return diff_count


# MAIN -------------------------------------------------------------------
if __name__ == "__main__":
    # command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Compare NetCDF values.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--digits', type=int, nargs=1,
                       help='allowed error digits')
    group.add_argument('-f', '--factor', type=float, nargs=1,
                       help='allowed relative error factor')
    parser.add_argument('-v', '--var', type=str, nargs=1, required=True,
                        help='dataset name')
    parser.add_argument('file', type=file, nargs=2,
                        help='netCDF file')
    parser.add_argument('-s', '--suppress', action='store_true', \
      help='Suppress standard output.')
    args = parser.parse_args()
    filename1 = args.file[0].name
    filename2 = args.file[1].name
    varname = args.var[0]
    digits = args.digits[0] if args.digits else None
    factor = args.factor[0] if args.factor else None

    if args.suppress: sys.stdout = os.devnull

    # calling main function
    diff_count = diff2files(filename1, filename2, varname, digits, factor)
    if (diff_count): sys.exit(1)

    if args.suppress: sys.stdout = sys.__stdout__


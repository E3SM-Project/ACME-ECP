#! /usr/bin/env python
#
# author: Andre Wehe, 2016
#

from __future__ import print_function
import sys, os

def findNcFiles(folder):
    # finds *.nc recursivley in a folder
    import fnmatch
    matches = []
    for root, dirnames, filenames in os.walk(folder):
        for filename in fnmatch.filter(filenames, '*.nc'):
            fullname = os.path.join(root, filename)[len(folder) + 1:]
            matches.append(fullname)
    return matches


def getVariables(filename):
    # package netCDF4
    # https://github.com/Unidata/netcdf4-python
    import netCDF4
    # Reading file 1 (output data)
    file = netCDF4.Dataset(filename, 'r')
    vars = [v for v in file.variables]
    file.close()
    return vars


def diff2folders(folder1, folder2, digits, factor):
    import nc_diff
    files1 = set(findNcFiles(folder1))
    files2 = set(findNcFiles(folder2))

    # warn for missing files
    for f in files2 - files1:
        print('WARNING: file ' + f + ' does not exist in folder ' + folder1)
    for f in files1 - files2:
        print('WARNING: file ' + f + ' does not exist in folder ' + folder2)

    # process common files
    diff_count = 0
    diff_vars = []
    intersection = files1.intersection(files2)
    for f in intersection:
        print('comparing file ' + f)
        fullname1 = os.path.join(folder1, f)
        fullname2 = os.path.join(folder2, f)
        vars1 = set(getVariables(os.path.join(folder1, f)))
        vars2 = set(getVariables(os.path.join(folder2, f)))

        # warn for missing variables
        """
        for v in vars2 - vars1:
            print('WARNING: variable ' + v +
                  ' does not exist in file ' + fullname1)
        for v in vars1 - vars2:
            print('WARNING: variable ' + v +
                  ' does not exist in file ' + fullname2)
        """

        # process common files
        var_intersection = vars1.intersection(vars2)
        for v in var_intersection:
            #print('comparing variable ' + v)
            dc = nc_diff.diff2files(fullname1, fullname2, v, digits, factor)
            diff_count += dc
            if (dc):
                diff_count += 1
                diff_vars.append(f + ':' + v)

    if diff_count == 0:
        print('all files identical')
    return diff_count


# MAIN -------------------------------------------------------------------
if __name__ == "__main__":
    # command line arguments
    import argparse
    parser = argparse.ArgumentParser(
        description='Recursively compare 2 folders of NetCDF files.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--digits', type=int, nargs=1,
                       help='allowed error digits')
    group.add_argument('-f', '--factor', type=float, nargs=1,
                       help='allowed relative error factor')
    parser.add_argument('folder', type=str, nargs=2,
                        help='folder with netCDF files')
    parser.add_argument('-s', '--suppress', action='store_true', \
      help='Suppress standard output.')
    args = parser.parse_args()
    folder1 = args.folder[0]
    folder2 = args.folder[1]
    digits = args.digits[0] if args.digits else None
    factor = args.factor[0] if args.factor else None

    if args.suppress: sys.stdout = os.devnull

    # calling main function
    diff_count = diff2folders(folder1, folder2, digits, factor)
    if (diff_count): sys.exit(1)

    if args.suppress: sys.stdout = sys.__stdout__


#!/usr/bin/env python

# make sure print behaves the same in 2.7 and 3.x
from __future__ import print_function

import os, sys, shutil
import subprocess as sub
if sys.version_info[0] < 3:
    import ConfigParser
else:
    import configparser as ConfigParser

# package netCDF4 (https://github.com/Unidata/netcdf4-python)
import netCDF4 as nc
import numpy as np

# for reversing the vertical direction (global attribute in netCDF)
revTopAtt = 'top_at_1'

def path_check(path):
  """
  Quick check if path exists.  Use before reading a file.
  """
  if not os.path.exists(path):
    sys.exit('Could not find %s, returning' % path)
# end path_check()

def spawn(cmd, noSplit=False, errStop=True):
  """
  Simplifies the call to a shell command in a Python session

  Call:
    results = spawn(cmd)

  Input:
    cmd -- a simple string that would be used at the Unix command line

  Keywords:
    noSplit -- boolean, if True no string split is performed on the
      standard output
    errStop -- boolean, if True the function will exit upon
      encountering any standard error that results from spawning cmd

  Returns:
    stOut, stErr -- lists of standard output and standard error
  """

  call = sub.Popen(cmd, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
  callout, callerr = call.communicate()
  # Python 3 returns byte strings that need to be decoded.
  callout, callerr = callout.decode('utf-8'), callerr.decode('utf-8')
  rCode = call.returncode

  # http://stackoverflow.com/questions/3630389/python-error-codes
  # http://stackoverflow.com/questions/18731791/determining-if-a-python-subprocess-segmentation-faults
  # https://linux.die.net/man/7/signal
  if rCode < 0:
    sys.exit('Fatal error in running %s. Err code %d' % (cmd, rCode) )

  if errStop and len(callerr) > 0: sys.exit(callerr)

  if noSplit:
    return callout, callerr
  else:
    stOut = callout.split()
    return stOut, callerr
  # end noSplit
# end spawn()

##################### From Andre's nc_diff.py ######################

def ncVarDiff(filename1, filename2, varname, digits, factor,
  validation=False):

  """
  Run the regression test (do files differ? what variables differ?
  by how much?) for given test step (unit)

  Call:
    diff_count =
      ncVarDiff(filename1, filename2, varname, digits, factor)

  Input:
    filename1 -- string, path to reference netCDF
    filename2 -- string, path to test netCDF
    varname -- string, name of netCDF variable to compare
    digits -- int, allowed error digits
    factor -- float, allowed relative error factor

  Keywords:
    validation -- boolean, run the validation differencing script,
      which is a more detailed regression test

  Returns:
    diff_count -- int, the number of differences in varname between
      filename1 and filename2
  """

  # Reading file 1 (reference data)
  file1 = nc.Dataset(filename1)
  var1 = file1.variables[varname]
  shape1 = var1.shape
  if validation: print(filename1, varname + str(shape1))
  data1 = var1[:]
  file1.close()

  # Reading file 2 (test data)
  file2 = nc.Dataset(filename2)
  var2 = file2.variables[varname]
  shape2 = var2.shape
  if validation: print(filename2, varname + str(shape2))
  data2 = var2[:]
  file2.close()

  # shape the same?
  if shape1 != shape2:
    #print('different shapes')
    return 0
  # endif shape1/2

  # inf values present?
  if np.isinf(data1).any():
    if validation:
      print('inf values found in ' + filename1 + ':' + varname)
      for r in np.argwhere(np.isinf(data1)): print(r)
    # endif validation
    return 0
  # endif data1 inf

  if np.isinf(data2).any():
    if validation:
      print('inf values found in ' + filename2 + ':' + varname)
      for r in np.argwhere(np.isinf(data2)): print(r)
    # endif validation
    return 0
  # endif data2 inf

  # NaN values present?
  if np.isnan(data1).any():
    if validation:
      print('NaN values found in ' + filename1 + ':' + varname)
      for r in np.argwhere(np.isnan(data1)): print(r)
    # endif validation
    return 0
  # endif data1 NaN

  if np.isnan(data2).any():
    if validation:
      print('NaN values found in ' + filename2 + ':' + varname)
      for r in np.argwhere(np.isnan(data2)): print(r)
    # endif validation
    return 0
  # endif data2 NaN

  # compare by number of digits
  if digits:
    # not a whole lot of testing with this keyword -- eliminate?
    # compare element by element
    diff_count = 0
    itr1 = np.nditer(data1, flags=['multi_index'])
    itr2 = np.nditer(data2, flags=['multi_index'])

    while not itr1.finished:
      v1 = itr1[0]
      v2 = itr2[0]
      v12 = (v1 + v2) / 2.
      numformat = "%." + str(digits) + "e"
      s1 = numformat % abs(v1)
      s2 = numformat % abs(v2)
      s12 = numformat % abs(v12)
      r1 = '%s|%s' % \
        (s1.replace(".", "")[0:digits], s1[len(s1) - 4:len(s1)])
      r2 = '%s|%s' % \
        (s2.replace(".", "")[0:digits], s2[len(s2) - 4:len(s2)])
      r12 = '%s|%s' % \
        (s12.replace(".", "")[0:digits], s12[len(s12) - 4:len(s12)])

      if (np.isnan(v1) or np.isnan(v2) or ((r1 != r12) and (r2 != r12))):
        pos1 = pos2 = 0
        for x in range(digits + 1):
          if (r1[x] != r12[x]):
            pos1 = x
            break
          # endif r1 != r12
        # end x loop

        for x in range(digits + 1):
          if (r2[x] != r12[x]):
            pos2 = x
            break
          # endif r2 != r12
        # end x loop

        pos = max(pos1, pos2) + 1
        if validation:
          # this has not been tested
          numformat2 = "%." + str(digits) + "e"
          print('digit-diff:', pos, varname + str(itr1.multi_index),
            'values:', numformat2 % v1, numformat2 % v2)
          print('digit-diff:',pos, varname+str(itr1.multi_index),
            'values:', r1, r12, r2)
        # end validation

        diff_count += 1
      # endif NaN v1 or v2

      itr1.iternext()
      itr2.iternext()

    # endwhile

    if validation:
      if (diff_count > 0):
        print('total differences: %d' % diff_count)
      else:
        print('identical %d digits' % digits)
    # end validation
  # endif digits

  # compare by relative difference
  if factor:

    # Compare output data to reference data.
    def myDiff(a, b):
      # both values are 0; they are identical
      if a == 0 and b == 0: return 0.0

      # a is 0, but b is non-zero
      if a == 0: return np.nan

      # b is 0, but a is non-zero
      if b == 0: return np.nan

      # compute the absolute relative error
      return abs(b / a - 1.0)
    # end myDiff()

    vMyDiff = np.vectorize(myDiff)
    diff = vMyDiff(data1, data2)
    err = np.argwhere(diff >= factor)
    diff_count = len(err)

    # if we in validation mode, print every set of values that are
    # different by more than "factor" variable
    if validation:
      for r in err:
        print(', '.join(\
          ['index:' + str(r), 'ref value: ' + str(data1[tuple(r)]), \
           'test value: ' + str(data2[tuple(r)]), \
           'relative difference: ' + str(diff[tuple(r)])]))

      print('factor:', factor)
      if (diff_count > 0):
        print('total differences: ' + str(diff_count))
      else:
        print('identical')
    else:
      # if we're not in validation mode, calculate and print some
      # summary statistics
      if diff_count > 0:
        nDiffFactor = diff.size
        flatDiff = diff.flatten()
        diffIdx = np.where(flatDiff >= factor)[0]
        relDiffPer = diff_count/float(nDiffFactor) * 100

        outStr1 = '%s: %.2f%% of values differ by more than %f; ' % \
          (varname, relDiffPer, factor)
        outStr2 = 'Percentage Differences Range: [%e, %e]' % \
          (flatDiff[diffIdx].min()*100, flatDiff[diffIdx].max()*100)
        print(outStr1)
        print(outStr2)
      # end diff_count
    # endif validation
  # endif factor

  return diff_count

# end ncVarDiff()

##################### End Andre's nc_diff.py #######################

################## From Andre's nc_diff_folders.py ###################

def getVariables(filename):
  """
  Read in all variables from netCDF (filename)
  """

  ncObj = nc.Dataset(filename, 'r')
  varsNC = [v for v in ncObj.variables]
  ncObj.close()
  return varsNC
# end getVariables

def ncDiff(testDir, ref, test, relDiffCut=0.0001, validation=False):
  """
  Run the regression test (do files differ? what variables differ?
  by how much?) for given test step (unit)

  Call:
    status = ncDiff(testDir, ref, test)

  Input:
    testDir -- string, directory in which unit test shell
      script exists
    ref -- string, path to reference netCDF file
    test -- string, path to test netCDF file

  Keywords:
    relDiffCut -- float, percentage difference by which any reference-
      test is considered significantly different (ie, anything above
      this cutoff)
    validation -- boolean, run the validation differencing script,
      which is a more detailed regression test

  Returns:
    1 if the files have any differences, 0 if not
  """

  print('TEST: %s' % os.path.basename(testDir))
  curDir = os.getcwd()

  # test/util contains the nc_diff library
  sys.path.append('../util')

  path_check(testDir)
  os.chdir(testDir)

  diffCount = 0
  print('Comparing %s and %s ' % (ref, test) )
  varsRef = set(getVariables(ref))
  varsTest = set(getVariables(test))

  # warn for missing variables
  for v in varsTest - varsRef:
    print('WARNING: variable %s does not exist in file %s' % (v, ref))
  for v in varsRef - varsTest:
    print('WARNING: variable %s does not exist in file %s' % (v, test))

  # process common variables
  varIntersection = varsRef.intersection(varsTest)
  for varNC in varIntersection:
    if validation: print('Comparing variable %s' % varNC)

    # none of the software i've seen specifies the "digits"
    # argument, so i'm setting it to None always
    dc = ncVarDiff(ref, test, varNC, None, relDiffCut, \
      validation=validation)

    if dc > 0:
      #print('%s has %d indices that are different' % (varNC, dc ) )
      diffCount += 1
    # endif dc
  # end varNC loop

  os.chdir(curDir)

  return 1 if diffCount > 0 else 0
# end ncDiff()

################## End Andre's nc_diff_folders.py ###################

def reverseVertical(inFile):
  """
  Reverse vertical dimension for all applicable variables in given
    netCDF file

  Input
    inFile -- string, path to netCDF to be modified

  Output
    nothing. inFile is overwritten

  Keywords
  """

  # open netCDF4 object and loop over variables in it
  ncObj = nc.Dataset(inFile, 'r+')
  ncVars = list(ncObj.variables)

  for ncVar in ncVars:
    inVar = ncObj.variables[ncVar]

    # these are just layer and level indices and should not be
    # reversed
    ll = ['lev', 'lay']
    if ncVar in ll: continue

    dims = inVar.dimensions; nDims = inVar.ndim

    # determine which axis to invert (either lev or lay nc dimension)
    for l in ll:
      if l in dims: axis = dims.index(l)
    # end l loop

    # is there a vertical dimension (i.e., has axis been assigned)?
    # if not, proceed to next array
    if not 'axis' in locals(): continue

    # not optimized for arrays with more than 3 dimensions
    if axis == 0:
      if nDims == 1:
        outVar = inVar[::-1]
      elif nDims == 2:
        outVar = inVar[::-1, :]
      elif nDims == 3:
        outVar = inVar[::-1, :, :]
      # endif nDims

    elif axis == 1:
      if nDims == 2:
        outVar = inVar[:, ::-1]
      elif nDims == 3:
        # stupid level source arrays...
        """
        if ncVar == 'lev_src_dec':
          # move the "zero vector" to the end after inversion
          outVar = inVar[:, ::-1, :]
          goodOut = outVar[:, 1:, :]
          zeroOut = outVar[:, 0, :]
          outVar = np.zeros_like(outVar)
          outVar[:,:,:] = np.nan
          outVar[:, -1, :] = np.array(zeroOut)
          outVar[:, :-1, :] = np.array(goodOut)
        elif ncVar == 'lev_src_inc':
          # move the "zero vector" to the beginning inversion
          outVar = inVar[:, ::-1, :]
          goodOut = outVar[:, :-1, :]
          zeroOut = outVar[:, -1, :]
          outVar = np.zeros_like(outVar)
          outVar[:,:,:] = np.nan
          outVar[:, 0, :] = np.array(zeroOut)
          outVar[:, 1:, :] = np.array(goodOut)
        else:
          outVar = inVar[:, ::-1, :]
        """
        outVar = inVar[:, ::-1, :]
      # endif nDims

    elif axis == 2:
      outVar = inVar[:, :, ::-1]
    # end axis conditional

    ncObj.variables[ncVar][:] = outVar

    # so we don't carry axis to the next variable
    del(axis)
  # end loop over variables

  # These variables are referenced to the vertical ordering:
  #   "_inc" refers to increasing index along the vertical axis. So they
  #   need to be swapped.
  if('lev_src_inc' in ncVars and 'lev_src_dec' in ncVars):
    ncObj.renameVariable('lev_src_inc', 'temp')
    ncObj.renameVariable('lev_src_dec', 'lev_src_inc')
    ncObj.renameVariable('temp',        'lev_src_dec')

  ncObj.close()

  return True
# end reverseVertical()

def unitTest(inDict, verbose=False, reverse=False, skipDiff=False, \
  workCopy='inverted.nc'):
  """

  Call:
    unitTest(inDict)

  Input:
    inDict -- dictionary with the following key/value pairs
      directory: string, directory where unit tests are performed
      top_directory: string, top level directory that contains
        build/ and test/ subdirectories
      refNC: string, netCDF with reference model output
      testNC: string, netCDF with test model output
      coefficients: string, netCDF with LW or SW opt prop coefficients

  Keywords:
    verbose -- boolean; print executable output to standard output
    reverse -- boolean; if the vertical direction is inverted, it is
      assumed that this function was called for the uninverted file,
      and we are re-running the function such that the file staging
      does not need to be repeated. so if this keyword is True, no
      file staging is performed
    skipDiff -- boolean; skip the reference-test netCDF comparison
    workCopy -- string; filename of inverted netCDF

  Returns:
    diffCount -- 1 if netCDF files are at all different, 0 otherwise
  """

  inDir = inDict['directory']
  uTest = inDict['test_name']
  print('%s' % uTest)

  topDir = inDict['top_directory']
  testExe = inDict['executable']; path_check(testExe)
  exeOptions = inDict['options']
  if exeOptions:
    exe = '%s/%s %s' % (inDir, testExe, exeOptions)
  else:
    exe = '%s/%s' % (inDir, testExe)

  # just in case exeOptions is empty; we don't want any spaces at the
  # end of exe
  exe = exe.strip()

  # stage (copy, link, mv, etc.) some files and run the unit test
  coeffNC = inDict['coefficients']
  if coeffNC:
    path_check(coeffNC)

    # is the coeff NC already in the working directory?
    tempCoeffNC = 'coefficients.nc'
    if os.path.exists(tempCoeffNC): os.remove(tempCoeffNC)
    os.symlink(coeffNC, tempCoeffNC)
  # endif

  inNC = inDict['refNC'] if inDict['chainNC'] is None else \
    inDict['chainNC']
  outNC = inDict['workNC']
  testNC = inDict['testNC']

  shutil.copyfile(inNC, outNC)

  if reverse:
    # run the test exe on the working NC, copy the working NC to
    # an inverted file, revert working NC back to un-inverted
    # pressure grid, move to test, compare ref and test

    print('REVERSING VERTICAL DIMENSION')
    status = reverseVertical(outNC)
    if status:
      # top_at_1 should be set so as not confuse the RRTMGP
      # executables about where the TOA is
      ncObj = nc.Dataset(outNC, 'r+')
      if revTopAtt in ncObj.ncattrs(): ncObj.delncattr('top_at_1')
      ncObj.setncattr('top_at_1', 1)
      ncObj.close()

      # keep a copy of the inverted vertical file (workNC is
      # inverted back and overwritten before the tests, but we
      # may want to examine the inverted file as well)
      shutil.copyfile(outNC, workCopy)
    # endif status
  # endif reverse

  sOut, sErr = spawn(exe, noSplit=True)
  if verbose: print(sOut)

  # invert the "working" netCDF back to original pressure grid
  # before comparing to reference NC
  if reverse:
    print('REVERTING VERTICAL DIRECTION TO INITIAL STATE')
    status = reverseVertical(outNC)

    # resetting top_at_1
    ncObj = nc.Dataset(outNC, 'r+')
    if revTopAtt in ncObj.ncattrs(): ncObj.delncattr('top_at_1')
    ncObj.close()
  # endif reverse

  # in the original unit tests, we moved these guys to the test
  # directory, but for now let's copy so we can also copy outNC to
  # the next unit test (this is done in configSetup() w/
  # "replace" set)
  shutil.copyfile(outNC, testNC)

  # now do an NC-diff test
  if skipDiff:
    diffCount = np.nan
  else:
    # just in case we chained results -- we still want to compare
    # refNC to testNC, not chainNC to testNC
    inNC = inDict['refNC']
    diffCount = ncDiff(inDir, inNC, testNC, \
      relDiffCut=inDict['relative_diff_cut'], \
      validation=inDict['validation_switch'])
  # endif skipDiff

  return diffCount
# end unitTest()

def configSetup(configFile, chain=False, replace=False, \
  relDiffCut=0.0001, validInfo=False, revLayers=False, \
  build=False, rootDir='../', failQuit=False, **kwargs):
  """
  Run unit tests as specified by the user in a configuration file

  Call:
    configSetup(configFile)

  Input:
    configFile -- string, path to configuration (.ini) file that
      contains assignments for:
        executable -- string, name of test executable
        directory -- string, path in which executable exists
        results_src -- string, netCDF with reference results
        results_dst -- string, netCDF with test results
        results -- string, name of output netCDF file after a test
          model run
        coefficients -- string or list of strings, path to coefficients
          netCDF (input). this is only required for gas_optics

          multiple values can be assigned to this field (separated by
          a comma), but this option has not yet been extensively tested
        options -- string, arguments to be used with executable (just
          as they would be entered in the command line interface,
          without the executable name). this is optional

        The paths that are assigned to the variables are
        expected to be relative to the current working directory

  Keywords:
    chain -- boolean, if multiple values are assigned to the variables
      in configFile, then setting this keyword places the output from
      unit test n into the working directory of unit test n+1 so that
      the output from n is input into n+1
    replace -- boolean, replace the current reference netCDF with the
      test results that are produced by a run of a given test
      executable. this will be done for all tests that are performed
      in the chain defined by the config file.
    relDiffCut -- float, relDiff = |(a/b) - 1| with a being a value
      from the reference file and b being a value from the test file.
      this cutoff is the threshold by which the test and reference
      are considered significantly different
      (i.e., relDiff > relDiffCut)
    validInfo -- boolean; by default, the difference tests print out
      the name of the test and either "all files identical" or
      statistics for whatever variables differ significantly (see
      relDiffCut keyword). by setting this keyword, *every* test value
      that is significantly different from the reference value is
      printed to standard output, so this is much more extensive
      diagnostic
    revLayers -- boolean; reverse the vertical dimension and
      corresponding arrays
    build -- boolean, build the RRTMGP library and any executables
      that are to be used in the regression tests
    rootDir -- string, absolute path to RRTMGP root (which includes
      build, data, extensions, src, and test subdirs). this is
      necessary for the script to know where the executables are
    failQuit -- boolean, quit program as soon as any netCDF
      differences are found

    Overloaded Keywords (**kwargs, passed to unitTest)
      verbose -- boolean; print executable output to standard output
      skipDiff -- boolean; bypass the reference-test netCDF comparison

  Returns:
    Nothing
  """

  cParse = ConfigParser.ConfigParser()
  cParse.read(configFile)
  cpSections = cParse.sections()

  if build:
    # first build the RRTMGP library
    os.chdir('%s/build' % rootDir)

    buildOut, buildErr = spawn('make clean', errStop=False, \
      noSplit=True)
    print(buildOut)
    buildOut, buildErr = spawn('make', errStop=False, \
      noSplit=True)
    print(buildOut)

    os.chdir(rootDir)
  # end RRTMGP library build

  # loop over each section, which represents a unit test
  fileDiff = 0
  for iCPS, cps in enumerate(cpSections):

    # read in each assigned variable
    # these guys can be split with .split(',') if we add more
    exe = cParse.get(cps, 'executable')
    exeDir = cParse.get(cps, 'directory')
    refNC = cParse.get(cps, 'results_src')
    workNC = cParse.get(cps, 'results')
    testNC = cParse.get(cps, 'results_dst')

    exeOpt = cParse.get(cps, 'options') if \
      cParse.has_option(cps, 'options') else None

    coeffs = '%s/%s' % (rootDir, cParse.get(cps, 'coefficients')) if \
      cParse.has_option(cps, 'coefficients') else None

    eDirFull = '%s/%s' % (rootDir, exeDir)
    fullRefNC = '%s/%s' % (eDirFull, refNC)
    fullWorkNC = '%s/%s' % (eDirFull, workNC)
    fullTestNC = '%s/%s' % (eDirFull, testNC)

    chainNC = str(prevTest) if chain and iCPS > 0 else None

    unitDict = {'executable': exe, \
      'directory': eDirFull, 'top_directory': rootDir, \
      'refNC': fullRefNC, 'testNC': fullTestNC, 'workNC': fullWorkNC, \
      'chainNC': chainNC, 'coefficients': coeffs, \
      'options': exeOpt, 'test_name': cps, \
      'relative_diff_cut': relDiffCut, 'validation_switch': validInfo}

    os.chdir(eDirFull)

    if build:
      # now build the executable
      os.chdir('build')

      buildOut, buildErr = spawn('make', errStop=False, noSplit=True)
      print(buildOut)

      os.chdir('..')
    # end build

    diffCount = unitTest(unitDict, verbose=kwargs['verbose'],
      skipDiff=kwargs['skipNCD'])

    # redo the test with inverted vertical dimension
    if revLayers:
      cpWorkNC = '%s/inverted_%s' % (eDirFull, workNC)
      diffCount = unitTest(unitDict, verbose=kwargs['verbose'],
        skipDiff=kwargs['skipNCD'], reverse=True, workCopy=cpWorkNC)
    # endif revLayers

    if diffCount == 0: print('No differences in %s' % cps)
    if diffCount > 0 and failQuit:
      sys.exit('Differences found, returning')

    fileDiff += diffCount

    if replace: shutil.copyfile(os.path.basename(workNC), fullRefNC)
    os.chdir(rootDir)

    prevTest = str(fullTestNC)

  # end loop over sections

  if fileDiff == 0: print('all files identical')

# end configSetup

def configSetupSHDOMPP(configFile, rootDir='../', \
  **kwargs):
  """
  Run executables as specified by an input configuration file
  for SHDOMPP validation

  Call:
    configSetup(configFile)

  Input:
    configFile -- string, path to configuration (.ini) file that
      contains assignments for:
        executable -- string, name of test executable
        directory -- string, path in which executable exists,
          relative to VALIDATION subdirectory in rootDir
        results_src -- string, netCDF with reference results,
          relative to TEST subdirectory in rootDir
        results_dst -- string, netCDF with test results
          relative to TEST subdirectory in rootDir
        results -- string, name of output netCDF file after a test
          model run, relative to "directory" input
        options -- string, arguments to be used with executable (just
          as they would be entered in the command line interface,
          without the executable name). any paths should be
          relative to rootDir. this is optional

          multiple values can be assigned to this field (separated by
          a comma), but this option has not yet been extensively
          tested

  Output:

  Keywords:
    rootDir -- string, absolute path to RRTMGP root
      This is the directory in which build, data, extensions, src,
      test, and validation reside). this is necessary for the script
      to know where the executables are

    Overloaded Keywords (**kwargs, passed to unitTest)
      verbose -- boolean; print executable output to standard output
      skipDiff -- boolean; bypass the reference-test netCDF comparison
  """

  validDir = '%s/validation' % rootDir; path_check(validDir)

  cParse = ConfigParser.ConfigParser()
  cParse.read(configFile)
  cpSections = cParse.sections()

  # loop over each section, which represents a unit test
  fileDiff = 0
  for iCPS, cps in enumerate(cpSections):

    # read in each assigned variable
    # these guys can be split with .split(',') if we add more
    exe = cParse.get(cps, 'executable')
    exeDir = cParse.get(cps, 'directory')
    refNC = cParse.get(cps, 'results_src')
    workNC = cParse.get(cps, 'results')
    testNC = cParse.get(cps, 'results_dst')

    # if the data/ subdir is part of the options, we need to prepend
    # the rootDir because it is assumed that the paths in options
    # are relative to root
    exeOpt = cParse.get(cps, 'options') if \
      cParse.has_option(cps, 'options') else None
    if 'data' in exeOpt:
      exeOpt = exeOpt.replace('data', '%s/data' % rootDir)

    eDirFull = '%s/%s' % (validDir, exeDir)
    fullRefNC = '%s/%s' % (rootDir, refNC)
    fullWorkNC = '%s/%s/%s' % (validDir, exeDir, workNC)
    fullTestNC = '%s/%s' % (rootDir, testNC)

    unitDict = {'executable': exe, \
      'directory': eDirFull, 'top_directory': rootDir, \
      'refNC': fullRefNC, 'testNC': fullTestNC, 'workNC': fullWorkNC, \
      'coefficients': None, 'options': exeOpt, 'test_name': cps, \
      'relative_diff_cut': None, 'validation_switch': None}

    os.chdir(eDirFull)
    diffCount = unitTest(unitDict, verbose=kwargs['verbose'],
      skipDiff=True)
  # end loop over sections
# end configSetupSHDOMPP()

def runOpticalProps(inDir, replace=False, verbose=False, build=False):
  """
  Call:
    runOpticalProps(inDir)


  Input:
    inDir -- string, path to optical_props/ test

  Keywords:
    replace -- boolean, move output rrtmgp-inputs-outputs.nc to
      ref/ dir instead of default test/ dir (have not yet implemented)
      STILL NEED TO IMPLEMENT THIS
    verbose -- boolean; print executable output to standard output
    build -- boolean, build the optical properties unit test

  Returns:
    Nothing
  """

  # we're assuming the RRTMGP library was already built
  if args.build:
    curDir = os.getcwd()

    # now build the executable
    os.chdir('%s/build' % inDir)

    buildOut, buildErr = spawn('make', errStop=False, noSplit=True)
    print(buildOut)

    os.chdir(curDir)
  # end build

  path_check(inDir)
  curDir = os.getcwd()
  os.chdir(inDir)
  exe = '%s/test_optical_props' % inDir
  print('Optical Properties')
  sOut, sErr = spawn(exe, noSplit=True)
  if verbose: print(sOut)
  sOut, sErr = spawn('mv *.nc test/')
  print()
  os.chdir(curDir)

  return
# end runOpticalProps()

def cleanUp(inFile, rootDir='../', removeTest=True):
  """
  Remove all intermediate (staging) files from unit test execution

  Call:
    cleanUp(inFile)

  Input:
    inFile -- string, path to configuration file used in
      configSetup()

  Keywords:
    rootDir -- string, path to RRTMGP root (which includes build,
      data, extensions, src, and test subdirs). this is necessary for
      the script to know where the executables are
    removeTest -- boolean, removes test netCDF files as well as the
      coefficients netCDF

  Returns:
    Nothing
  """

  print('Cleaning up intermediate files in %s' % inFile)

  cParse = ConfigParser.ConfigParser()
  cParse.read(inFile)
  cpSections = cParse.sections()

  # strings that point to staging/intermdiate files
  # relative paths first
  relPaths = []

  # loop over each section, which represents a unit test
  for iCPS, cps in enumerate(cpSections):
    print('  Removing files in %s' % cps)
    cpsDir = cParse.get(cps, 'directory')

    relPaths.append('%s/%s' % (cpsDir, cParse.get(cps, 'results')))

    if removeTest: \
      relPaths.append('%s/%s' % \
        (cpsDir, cParse.get(cps, 'results_dst')))

    if cParse.has_option(cps, 'coefficients'):
      relPaths.append('%s/coefficients.nc' % cpsDir)

  # end loop over sections

  absPaths = ['%s/%s' % (rootDir, rPath) for rPath in relPaths]

  for path in absPaths:
    if os.path.exists(path): os.remove(path)

  return True
# end cleanUp()

if __name__ == '__main__':
  import socket
  if socket.gethostname() == 'rotor':
    sys.path.append('/home/rpernak/python_lib/')

  import argparse

  parser = argparse.ArgumentParser(\
    description='Run selected unit tests for RRTMGP builds.')
  parser.add_argument('--environment', type=str, \
    help='If set to a string corresponding to one of the ' + \
    'available environments ("conda info --envs" at CLI), ' + \
    'a shell command will be executed to set the environment ' + \
    'to this string.')
  parser.add_argument('--ref_replace', action='store_true', \
    help='Move rrtmgp-inputs-outputs.nc that results from ' + \
    'test executable runs into the ref/ subdirectory instead ' + \
    'of the default test/ subdir.')
  parser.add_argument('--test_config_file', type=str, nargs='+', \
    help='Name of config file(s) that the user can use to setup ' + \
    'the regression test schemes. The default is to process ' + \
    'all of the .ini files in the working directory except ' + \
    'for the user-defined chain (user_define_workflow_config.ini).')
  parser.add_argument('--unit_chain', action='store_true', \
    help='Used in conjuction with --test_config_file. If set, ' + \
    'it is assumed that the config file specifies that many ' + \
    'unit tests are to be performed but independent of each ' + \
    'other (the default action). This forces the output of ' + \
    'each unit test to be used as input into the next test.')
  parser.add_argument('--optical_props', action='store_true', \
    help='Run the optical properties test.')
  parser.add_argument('--cleanup', action='store_true', \
    help='Remove intermediate files specified in the input ' + \
    'configuration files, then exit from program.')
  parser.add_argument('--rel_diff_cut', type=float, default=0.0001,\
    help='Percentage difference over which any reference-' + \
    'test difference is considered significantly different ' + \
    '(i.e., anything above this cutoff)')
  parser.add_argument('--verbose', action='store_true', \
    help='If set, prints the standard output from each executable.')
  parser.add_argument('--very_verbose', action='store_true', \
    help='Instead of only returning statistics on the ' + \
    'differences, return the filename, variable name, and array ' + \
    'indices of every difference that exists in the regression tests.')
  parser.add_argument('--reverse', action='store_true', \
    help='Reverse the vertical dimension and the corresponding ' + \
    'arrays.')
  parser.add_argument('--build', action='store_true', \
    help='Build the RRTMGP library and executables before ' + \
    'running any tests.')
  parser.add_argument('--root_dir', type=str, default='../', \
    help='This script runs with a number of assumptions on ' + \
    'where directories and executables exist. This keyword ' + \
    'specifies what the RRTMGP root directory is, and then ' + \
    'all of the paths in the configuration files (assumed to be ' + \
    'in root_dir) will be relative to the test root.')
  parser.add_argument('--no_diff', action='store_true', \
    help='If set, runs the executables in the configuration file ' + \
    'but does not perform the subsequent reference-test netCDF ' + \
    'comparison.')
  parser.add_argument('--quit_on_fail', action='store_true', \
    help='Quit the script as soon as a significant difference ' + \
    'is found between a reference netCDF and its corresponding ' + \
    'test netCDF.')
  parser.add_argument('--validation', action='store_true', \
    help='Grab .ini files from validation/ subdirectory ' + \
    'maintained by Frank Evans and Robert Pincus rather than the ' + \
    'default test/ subdir.')
  args = parser.parse_args()

  baseDir = args.root_dir

  # default root_dir should be three levels above the scripts
  # directory in which run_tests.py resides
  rootRel = '../../..' if baseDir is None else str(baseDir)
  cwd = os.getcwd()

  # get the absolute path of the root and replace baseDir with it
  os.chdir(rootRel)
  baseDir = os.getcwd()
  os.chdir(cwd)
  # endif baseDir

  # build, test, and validation directories must exist in the current
  # working directory (CWD)
  path_check(baseDir)
  testDir = '%s/test' % baseDir; path_check(testDir)
  if args.build: path_check('%s/build' % baseDir)

  iniSub = 'validation' if args.validation else 'test'
  iniDir = '%s/%s' % (baseDir, iniSub); path_check(iniDir)

  newRef = args.ref_replace

  # set Python environment
  pEnv = args.environment
  if pEnv: sOut, sErr = spawn('source activate %s' % pEnv)

  # user-defined test procedure
  cFiles = args.test_config_file

  # default is to loop over all config files in test dir
  if cFiles is None: cFiles, sErr = spawn('ls %s/*.ini' % iniDir)

  if args.cleanup:
    for cFile in cFiles:
      base = os.path.basename(cFile)
      cleanUp(cFile, rootDir=baseDir)
    # end loop over config files

    sys.exit('Intermediate files have been removed')
  # end cleanUp

  for cFile in cFiles:
    base = os.path.basename(cFile)

    print('Working on %s' % cFile)
    path_check(cFile)
    if args.validation:
      configSetupSHDOMPP(cFile, rootDir=baseDir, verbose=args.verbose)
    else:
      configSetup(\
        cFile, replace=newRef, chain=args.unit_chain, \
        relDiffCut=args.rel_diff_cut, validInfo=args.very_verbose, \
        verbose=args.verbose, revLayers=args.reverse, \
        build=args.build, rootDir=baseDir, skipNCD=args.no_diff, \
        failQuit=args.quit_on_fail)
    print("\n")
    # endif validation
  # end loop over config files

  # are we doing anything with this unit test anymore?
  """
  if args.optical_props:
    # we eventually might wanna change this so that optical props
    # is run in configSetup() like everything else
    runOpticalProps('%s/optical_props' % testDir, replace=newRef, \
      verbose=args.verbose, build=args.build)
  # end optical_props
  """
# end main()

#
# What we want -- plots of {bias, RMSE, max or 95% err?...} of n RRTMGP results vs. reference LBLRTM result
#  or instead of n, maybe test and ref
#  For simplicity could do only broadband fluxes up, down, direct down if present, *maybe* broadband heating rates
#    this implies that we're only concerned with flux_compute
#  How to specify files? Directories within validation/ in parallel to test/ would be simplest. Maybe then we rename files
#    and store model/source as an attribute
# No - files provide two-three directories, two files, list of variables
#   Plots are RMS differences across columns vs height (and maybe across heights vs column)
import os, ConfigParser
import numpy   as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages

def path_check(path):
  """
  Quick check if path exists.  Use before reading a file.
  """
  if not os.path.exists(path):
    sys.exit('Could not find %s, returning' % path)
# end path_check()

def compareTwoFiles(refFile, testFile, figFile, vars):
  """
  Plots the RMSE, measured as test-ref accumulated across the column dimension, for
  a set of variables.

  Call:
    compareTwoFiles(refFile, testFile, figFile, vars)

  Input:
    refFile -- full path to netCDF file containing reference results
    testFile -- full path to netCDF file containing test results
    figFile -- full path and name of PDF file to be produced
    vars -- set of variable names to plot over
  """
  plt.style.use('seaborn-dark-palette')

  ref = nc.Dataset(refFile)
  tst = nc.Dataset(testFile)
  pp = PdfPages(figFile)

  for v in vars:
    # ref.variables['flux_dn'].dims() = ('levl', 'col')
    # ref.variables['flux_dn'].shape = (43, 42)
    # np.mean(ref.variables['flux_dn'], axis = 1) is a vector of length 43
    a = ref.variables[v].dimensions.index('col')
    # Also want to average across band if present
    a = a if not 'band' in ref.variables[v].dimensions else (a, ref.variables[v].dimensions.index('band'))
    delta = tst.variables[v][:] - ref.variables[v][:]
    if 'heating_rate' in v: delta = delta * (24. * 60. * 60.) # Convert to K/day
    bias = np.mean(delta,                 axis = a)
    rmse = np.sqrt(np.mean(delta * delta, axis = a))
    y = np.arange(0, rmse.shape[0]) # arrange plot in the vertical
    # Pressures are uniform across columns; use the first column, label with hPa
    y = ref.variables['p_lev'][:][:,0] if 'lev' in ref.variables[v].dimensions else ref.variables['p_lay'][:][:,0]
    y = y * .01 
    #
    # Mean profile
    #
    fig, ax = plt.subplots()
    if 'heating_rate' in v:
      ax.plot(np.mean(ref.variables[v][:] * (24. * 60. * 60.), axis = a), y)
    else:
      ax.plot(np.mean(ref.variables[v][:],                     axis = a), y)
    ax.set_title("Mean profile: " + v)
    fig.gca().invert_yaxis()
    pp.savefig()
    plt.close(fig)
    #
    # Bias and RMSE
    #
    fig, ax = plt.subplots()
    ax.plot(np.abs(bias), y, linestyle='--')
    ax.plot(       rmse,  y)
    ax.set_title("Bias (dashed) and RMSE in " + v)
    fig.gca().invert_yaxis()
    pp.savefig()
    plt.close(fig)

  pp.close()
  ref.close()
  tst.close()
 #  end compareTwoFiles()

if __name__ == '__main__':
  import argparse, glob

  parser = argparse.ArgumentParser(\
    description='Plot RMS differences between two file.')
  parser.add_argument('--validation_config_file', type=str, nargs='+', \
    help='Name of config file(s) that the user can use to describe ' + \
    'the desired validation plots. The default is to process ' + \
    'all of the .val files in the working directory.')
  args = parser.parse_args()

  rootDir = '..'
  cwd = os.getcwd()

  # user-defined test procedure
  cFiles = args.validation_config_file
  cParse = ConfigParser.ConfigParser()
  # default is to loop over all config files in test dir
  if cFiles is None: cFiles = glob.glob('*.val')

  for cFile in cFiles:
    base = os.path.basename(cFile)
    print('Making plots in %s' % cFile)
    path_check(cFile)
    cParse.read(cFile)
    cpSections = cParse.sections()
    for iCPS, cps in enumerate(cpSections):
      # read in each assigned variable
      # these guys can be split with .split(',') if we add more
      refDir = cParse.get(cps, 'ref_dir')
      refFile = cParse.get(cps, 'ref_file')
      testDir = cParse.get(cps, 'test_dir')
      testFile = cParse.get(cps, 'test_file')
      figDir = cParse.get(cps, 'fig_dir')
      figFile = cParse.get(cps, 'fig_file')
      varList = cParse.get(cps, 'variables').split(',')
      # Strip any leading spaces from variable names
      vars = []
      for v in varList: vars.append(v.lstrip())

      fullRefFile  = '%s/%s/%s' % (rootDir, refDir,  refFile)
      fullTestFile = '%s/%s/%s' % (rootDir, testDir, testFile)
      path_check(fullRefFile)
      path_check(fullTestFile)
      fullFigFile  = '%s/%s/%s' % (rootDir, figDir,  figFile)

      compareTwoFiles(fullRefFile, fullTestFile, fullFigFile, vars)

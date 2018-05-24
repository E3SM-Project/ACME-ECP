#
# Validation of cloudy-sky cases, meaning comparison of RRTGMP vs 2-stream SHDOMPP
#
import numpy as np
import os, sys
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ---------------------------------------------------------------------------------
def path_check(path):
  """
  Quick check if path exists.  Use before reading a file.
  """
  if not os.path.exists(path):
    sys.exit('Could not find %s, returning' % path)
# ---------------------------------------------------------------------------------

def compFluxesByProfile(ref_file, tst_file, clr_file, pdf_file):
  """
  Compare two sets of clouds profile by profile.
  """
  ref   = Dataset(ref_file) #SHDOM results
  tst   = Dataset(tst_file) #RRTMGP results

  nlev  = ref.variables['p_lev'][:].shape[0]
  ncol  = ref.variables['p_lev'][:].shape[1]
  # Surface is at first layer in these files; we could check this programmatically
  sfc, toa = (0, nlev-1)
  # By inspection, this is the band that contains 0.67 microns in the SW
  gpt_idx = 161
  # Pressure axis for plots
  p_lev   = -tst.variables['p_lev'][:][...,0]/100.
  # Midpoint of the pressures
  mid_vps = p_lev.min()/2.
  mid_idx = (np.abs(p_lev-mid_vps)).argmin()

  #
  # Cloud optical thickness
  #
  clr = Dataset(clr_file) # Clear-sky file
  cld_tau = ref.variables['tau'][:][gpt_idx, ...] - clr.variables['tau'][:][gpt_idx, ...]
  clr.close()

  #
  # Cloud radiative effect
  #
  rCREdn = tst.variables['flux_dn'][:]
  sCREdn = ref.variables['flux_dn'][:]
  rCREup = tst.variables['flux_up'][:]
  sCREup = ref.variables['flux_up'][:]
  for i in range(0, ncol):
    rCREdn[...,i]  = rCREdn[...,i] - tst.variables['flux_dn'][...,ncol-1]
    sCREdn[...,i]  = sCREdn[...,i] - ref.variables['flux_dn'][...,ncol-1]
    rCREup[...,i]  = rCREup[...,i] - tst.variables['flux_up'][...,ncol-1]
    sCREup[...,i]  = sCREup[...,i] - ref.variables['flux_up'][...,ncol-1]

  plt.style.use('seaborn-dark-palette')
  pp = PdfPages(pdf_file)
  # Loop over cloudy profiles
  for i in range(0, ncol):
    fig, axes = plt.subplots(2,2, sharey='row')
    #
    # Plot optical depth profile
    #
    axes[0,0].plot(np.append(0, cld_tau[...,i]), p_lev)
    axes[0,0].text(2., -100, 'Optical depth', color = 'black')
    #
    # Flux profiles
    #
    axes[0,1].plot(ref.variables['flux_dn'][...,i], p_lev)
    axes[0,1].plot(ref.variables['flux_up'][...,i], p_lev)
    axes[0,1].plot(tst.variables['flux_dn'][...,i], p_lev, 'r--', color = "C0")
    axes[0,1].plot(tst.variables['flux_up'][...,i], p_lev, 'r--', color = "C1")
    axes[0,1].text(tst.variables['flux_dn'][mid_idx,i], mid_vps, "Flux dn", color = "C0", ha = 'right')
    axes[0,1].text(tst.variables['flux_up'][mid_idx,i], mid_vps, "Flux up", color = "C1", ha = 'left')
    #
    # Flux error
    #
    axes[1,0].plot(tst.variables['flux_dn'][...,i] - ref.variables['flux_dn'][...,i], p_lev)
    axes[1,0].plot(tst.variables['flux_up'][...,i] - ref.variables['flux_up'][...,i], p_lev)
    axes[1,0].plot((tst.variables['flux_dn'][...,i] - tst.variables['flux_up'][...,i]) - \
                   (ref.variables['flux_dn'][...,i] - ref.variables['flux_up'][...,i]), \
                                                                                         p_lev)
    axes[1,0].text(0., -100, "Error")
    #
    # CRE error assuming clear sky is in the last columns
    #
    axes[1,1].plot( rCREdn[...,i] - sCREdn[...,i],  p_lev)
    axes[1,1].plot( rCREup[...,i] - sCREup[...,i],  p_lev)
    axes[1,1].plot((rCREdn[...,i] - rCREup[...,i]) - \
                   (sCREdn[...,i] - sCREup[...,i]), p_lev)
    axes[1,1].text(0., -100, "CRE error")
    pp.savefig()
    plt.tight_layout()
    plt.close()

  #
  # Plot of error in CRE at TOA and surface vs. optical thickness
  #
  plt.scatter(cld_tau.sum(axis=0), rCREdn[sfc,...] - sCREdn[sfc,...], )
  plt.scatter(cld_tau.sum(axis=0), rCREup[toa,...] - sCREup[toa,...])
  plt.text(cld_tau.sum(axis=0).max()/3.,      1., "dn@sfc", color = "C0")
  plt.text(cld_tau.sum(axis=0).max()/3. * 2, -1., "up@toa", color = "C1")
  plt.title('CRE error vs optical thickness')
  pp.savefig()
  plt.close()

  pp.close()
  ref.close()
  tst.close()

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(\
    description='Makes profile-by-profile plots of fluxes and ' +
    'cloud radiative effects. Needs a file with clear-sky results to compute CRE.' )
  parser.add_argument('ref_file', type=str, \
    help='Name of reference file.')
  parser.add_argument('tst_file', type=str, \
    help='Name of test file.')
  parser.add_argument('clr_file', type=str, \
    help='Name of file with clear-sjy results.')
  parser.add_argument('PDF_file', type=str, \
    help='File where the plots go.')

  args = parser.parse_args()
  path_check(args.ref_file)
  path_check(args.tst_file)
  path_check(args.clr_file)
  compFluxesByProfile(args.ref_file, args.tst_file, args.clr_file, args.PDF_file)

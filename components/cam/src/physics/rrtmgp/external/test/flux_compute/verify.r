setwd('Codes/RRTMGP/branches/test-revision/test/flux_compute') 

f1 = nc_open('../../data/rrtmgp-lw-flux-inputs-outputs-garand-all-r674.nc') 
f2 = nc_open('rrtmgp-lw-flux-inputs-outputs.nc') 

compare = function(f1, f2, varName, thresh = 0) { 
  ref = ncvar_get(f1, varName)
  tst = ncvar_get(f2, varName)
  ref[abs(ref) < thresh * max(abs(ref))] = NA
  tst[abs(ref) < thresh * max(abs(tst))] = NA
  return(range((tst-ref)/ref, na.rm = T))
} 
  
compare(f1, f2, "flux_dn")   
compare(f1, f2, "flux_up")   
compare(f1, f2, "band_flux_dn")   
compare(f1, f2, "band_flux_up")   
compare(f1, f2, "heating_rate")   
compare(f1, f2, "band_heating_rate")   

compare(f1, f2, "flux_dn", thresh = .01)   
compare(f1, f2, "flux_up", thresh = .01)   
compare(f1, f2, "band_flux_dn", thresh = .01)   
compare(f1, f2, "band_flux_up", thresh = .01)   
compare(f1, f2, "heating_rate", thresh = .01)   
compare(f1, f2, "band_heating_rate", thresh = .01)   


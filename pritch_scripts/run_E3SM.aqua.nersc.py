#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
import os
import subprocess as sp
newcase,config,build,clean,submit,continue_run = False,False,False,False,False,False

# clean        = True
# newcase      = True
# config       = True
build        = True
submit       = True
# continue_run = True
queue = 'regular'  # regular or debug on Cori (30 min time limit for debug queue)

ne,npg = 30,2        # use ne = 4 or 30, leav2e npg=0 for now
compset = 'F-EAMv1-AQP1'   # FC5AV1C-L / F-EAMv1-AQP1 / FSP1V1 / FSP2V1
# compset = 'F-EAMv1-AQP8'

cld = 'SP1'    # ZM / SP1 / SP2
crm_nx,crm_ny,crm_dx = 64,1,1000
# crm_nx,crm_ny,crm_dx = 64,16,1000

if queue=='debug' : 
   stop_opt,stop_n,resub,walltime = 'ndays',1,1,'0:30:00'
   # stop_opt,stop_n,resub,walltime = 'nsteps',2,1,'0:10:00'
else:
   stop_opt,stop_n,resub,walltime = 'ndays',5,0,'4:00:00'
   # if continue_run : resub=25


res = 'ne'+str(ne) if npg==0 else  'ne'+str(ne)+'pg'+str(npg)
cldc = cld+'_'+str(crm_nx)+'x'+str(crm_ny)+'_'+str(crm_dx)+'m' if 'SP' in cld else cld


# case = 'E3SM_'+cldc+'_'+res+'_'+compset+'_00'  # runs for pg2 paper
# case = 'E3SM_'+cldc+'_'+res+'_'+compset+'_01'  # runs for pg2 paper - alt setting to enhance imprinting

# case = 'E3SM_TEST-HVx1_'+cldc+'_'+res+'_'+compset+'_01'  # control for hyperviscosity test
# case = 'E3SM_TEST-HVx2_'+cldc+'_'+res+'_'+compset+'_01'  # np4 test w/ enhanced hyperviscosity
# case = 'E3SM_TEST-RESTART_'+cldc+'_'+res+'_'+compset+'_00' 
# case = 'E3SM_TEST-NO-THRDS_'+cldc+'_'+res+'_'+compset+'_00' 
# case = 'E3SM_TEST-NO-PERTURB_'+cldc+'_'+res+'_'+compset+'_00' 

case = 'E3SM_TEST-DAMP_'+cldc+'_'+res+'_'+compset+'_00' # 30 min damping

# case = case+'_debug-on'
# case = case+'_checks-on'

top_dir  = '/global/homes/w/whannah/E3SM/'
case_dir = top_dir+'Cases/'

# src_dir  = top_dir+'E3SM_SRC1/'
src_dir  = top_dir+'E3SM_SRC2/'
# src_dir  = top_dir+'E3SM_SRC3/'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

num_dyn = ne*ne*6
dtime = 20*60           # use 20 min for SP (default is 30 min for E3SM @ ne30)

if case == 'E3SM_'+cldc+'_'+res+'_'+compset+'_01' : dtime = 10*60

if ne==120: dtime = 5*60
ncpl  = 86400 / dtime

# Enforce max node limit on debug queue
if queue=='debug' and num_dyn>(64*32) : num_dyn = 64*32

#---------------------------------------------------------------------------------------------------
# Create new case
#---------------------------------------------------------------------------------------------------
if newcase :
   grid = res+'_'+res
   cmd = src_dir+'cime/scripts/create_newcase -case '+case_dir+case
   cmd = cmd + ' -compset '+compset+' -res '+grid#+' --pecount '+str(num_dyn)+'x1'
   if '_gnu_' in case : cmd = cmd + ' --compiler gnu '
   os.system(cmd)
#---------------------------------------------------------------------------------------------------
# Configure
#---------------------------------------------------------------------------------------------------
os.chdir(case_dir+case+'/')
if config : 

   # if compset not in ['FSP1V1','FSP2V1'] :
   #    if cld=='ZM' :
   #       cam_opt  = ' -phys cam5  -rad rrtmg -nlev 72 '  \
   #                 +' -clubb_sgs -microphys mg2 '        \
   #                 +' -chem linoz_mam4_resus_mom_soag '  \
   #                 +' -rain_evap_to_coarse_aero '        \
   #                 +' -bc_dep_to_snow_updates '

   if 'SP' in cld and compset not in ['FSP1V1','FSP2V1'] :
      # set options common to all SP setups
      nlev_crm = 58
      crm_dt   = 5
      cam_opt = ' -phys cam5 -use_SPCAM  -rad rrtmg -nlev 72 -microphys mg2 ' \
               +' -crm_nz '+str(nlev_crm) +' -crm_adv MPDATA '                \
               +' -crm_nx '+str(crm_nx)   +' -crm_ny '+str(crm_ny)            \
               +' -crm_dx '+str(crm_dx)   +' -crm_dt '+str(crm_dt)            \
               +' -crm_nx_rad 4 -crm_ny_rad 1 '
      # 1-moment microphysics
      if cld=='SP1': cam_opt = cam_opt + ' -SPCAM_microp_scheme sam1mom -chem none '
      # 2-moment microphysics
      if cld=='SP2': cam_opt = cam_opt + ' -SPCAM_microp_scheme m2005  '      \
                                       +' -chem linoz_mam4_resus_mom_soag '   \
                                       +' -rain_evap_to_coarse_aero '         \
                                       +' -bc_dep_to_snow_updates '

      cpp_opt = ' -DSP_DIR_NS -DSP_MCICA_RAD '

      if 'TEST-NO-PERTURB' in case : cpp_opt = cpp_opt+' -DSP_DISABLE_PERTURB '
      if 'TEST-DAMP_'      in case : cpp_opt = cpp_opt+' -DSP_CRM_DAMPING '

      cam_opt = cam_opt+' -cppdefs \''+cpp_opt+'\' '
      if 'AQP' in compset : cam_opt = cam_opt+' -aquaplanet '

      os.system('./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS  -val  \"'+cam_opt+'\"' )

   # if compset in ['FSP1V1','FSP2V1'] : 
   #    os.system('./xmlchange --append -file env_build.xml -id CAM_CONFIG_OPTS  -val  \" -crm_nx_rad 4 \" ' )

   # if 'ESMT' in case :
   #    cpp_opt = ' -DSP_DIR_NS -DSP_MCICA_RAD '
   #    cpp_opt = cpp_opt+' -DSP_ESMT -DSP_USE_ESMT '
   #    os.system('./xmlchange --append -file env_build.xml -id CAM_CONFIG_OPTS  -val  \" -cppdefs \''+cpp_opt+'\'  \" ' )

   #-------------------------------------------------------
   # Set tasks and threads - disable threading for SP
   #-------------------------------------------------------
   if ne<=60 : 
      os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val '+str(num_dyn)+' ')
   else:
      os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 5400 ')

   # if compset in ['FSP1V1','FSP2V1'] : os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1 ')
   # if 'SP' in cld : os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1 ')
   os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1 ')
   
   if clean : os.system('./case.setup --clean')
   os.system('./case.setup --reset')
#---------------------------------------------------------------------------------------------------
# Build
#---------------------------------------------------------------------------------------------------
if build : 
   if 'debug-on' in case : os.system('./xmlchange -file env_build.xml -id DEBUG -val TRUE ')
   if clean : os.system('./case.build --clean')
   os.system('./case.build')
#---------------------------------------------------------------------------------------------------
# Write the namelist options and submit the run
#---------------------------------------------------------------------------------------------------
if submit : 

   if ne==30 and npg==1 :
      domain_file_lnd = '/global/cscratch1/sd/whannah/acme_scratch/init_files/domain.lnd.ne30pg1_gx1v6.191206.nc'
      domain_file_ocn = '/global/cscratch1/sd/whannah/acme_scratch/init_files/domain.ocn.ne30pg1_gx1v6.191206.nc'
   if ne==30 and npg==3 :
      domain_file_lnd = '/global/cscratch1/sd/whannah/acme_scratch/init_files/domain.lnd.ne30pg3_gx1v6.190918.nc'
      domain_file_ocn = '/global/cscratch1/sd/whannah/acme_scratch/init_files/domain.ocn.ne30pg3_gx1v6.190918.nc'
   if ne==30 and npg==4 :
      domain_file_lnd = '/global/cscratch1/sd/whannah/acme_scratch/init_files/domain.lnd.ne30pg4_gx1v6.190919.nc'
      domain_file_ocn = '/global/cscratch1/sd/whannah/acme_scratch/init_files/domain.ocn.ne30pg4_gx1v6.190919.nc'
   if ne==120 and npg==1 :
      domain_file_lnd = '/global/cscratch1/sd/whannah/acme_scratch/init_files/domain.lnd.ne120pg1_gx1v6.191206.nc'
      domain_file_ocn = '/global/cscratch1/sd/whannah/acme_scratch/init_files/domain.ocn.ne120pg1_gx1v6.191206.nc'
   if 'domain_file_ocn' in vars() :
      os.system("./xmlchange -file env_run.xml -id ATM_DOMAIN_FILE  -val "+domain_file_lnd)
      os.system("./xmlchange -file env_run.xml -id LND_DOMAIN_FILE  -val "+domain_file_lnd)
      os.system("./xmlchange -file env_run.xml -id OCN_DOMAIN_FILE  -val "+domain_file_ocn)
      os.system("./xmlchange -file env_run.xml -id ICE_DOMAIN_FILE  -val "+domain_file_ocn)
   

   if ne==120 and npg==1 : os.system('./xmlchange -file env_run.xml EPS_AGRID=1e-11' )
   if ne==120 and npg==2 : os.system('./xmlchange -file env_run.xml EPS_AGRID=1e-11' )
   if             npg==4 : os.system('./xmlchange -file env_run.xml EPS_AGRID=1e-11' )
   #-------------------------------------------------------
   # First query some stuff about the case
   #-------------------------------------------------------
   (din_loc_root   , err) = sp.Popen('./xmlquery DIN_LOC_ROOT    -value', \
                                     stdout=sp.PIPE, shell=True, \
                                     universal_newlines=True).communicate()
   (cam_config_opts, err) = sp.Popen('./xmlquery CAM_CONFIG_OPTS -value', \
                                     stdout=sp.PIPE, shell=True, \
                                     universal_newlines=True).communicate()
   cam_config_opts = ' '.join(cam_config_opts.split())   # remove extra spaces to simplify string query
   #-------------------------------------------------------
   # Namelist options
   #-------------------------------------------------------
   nfile = 'user_nl_cam'
   file = open(nfile,'w') 
   #------------------------------
   # Specify history output frequency and variables
   #------------------------------
   if stop_opt=='nsteps' :
      file.write(' nhtfrq    = 0,1,1 \n')
      file.write(' mfilt     = 1,10,10 \n')
   else:
      # file.write(' nhtfrq    = 0,-3,-3 \n')
      # file.write(' mfilt     = 1,8,8 \n')
      file.write(' nhtfrq    = 0,-3 \n')
      file.write(' mfilt     = 1,8 \n')
   file.write(" fincl1 = 'DYN_T','DYN_Q','DYN_U','DYN_OMEGA','DYN_PS' \n")
   file.write(" fincl2    = 'PS','TS'")
   file.write(             ",'PRECT','TMQ'")
   file.write(             ",'LHFLX','SHFLX'")             # surface fluxes
   file.write(             ",'TGCLDLWP','TGCLDIWP'")       # liq and ice water path
   file.write(             ",'FSNT','FLNT'")               # Net TOM heating rates
   file.write(             ",'FLNS','FSNS'")               # Surface rad for total column heating
   # file.write(             ",'FSNTC','FLNTC'")             # clear sky heating rates for CRE
   # file.write(             ",'LWCF','SWCF'")               # cloud radiative foricng
   # file.write(             ",'TAUX','TAUY'")               # surface stress
   # file.write(             ",'CLDLOW','CLDMED','CLDHGH','CLDTOT' ")
   # file.write(             ",'CLOUD','CLDLIQ','CLDICE'")   # 3D cloud fields
   # file.write(             ",'T','Q','Z3' ")               # 3D thermodynamic budget components
   # file.write(             ",'U','V','OMEGA'")             # 3D velocity components
   # file.write(             ",'QRL','QRS'")                 # 3D radiative heating profiles
   file.write(             ",'OMEGA','DYN_OMEGA','DIV'")
   # if 'SP' in cld :
      # file.write(         ",'SPDT','SPDQ'")               # CRM heating/moistening tendencies
      # file.write(         ",'SPTLS','SPQTLS' ")           # CRM large-scale forcing
      # file.write(         ",'SPQPEVP','SPMC'")            # CRM rain evap and total mass flux
      # file.write(         ",'SPMCUP','SPMCDN'")           # CRM saturated mass fluxes
      # file.write(         ",'SPMCUUP','SPMCUDN'")         # CRM unsaturated mass fluxes
      # file.write(         ",'SPTKE'")
      # if any(x in cam_config_opts for x in ["SP_ESMT","SP_USE_ESMT","SPMOMTRANS"]) : 
      #    file.write(",'ZMMTU','ZMMTV','uten_Cu','vten_Cu' ")
      # if "SP_USE_ESMT" in cam_config_opts : file.write(",'U_ESMT','V_ESMT'")
      # if "SPMOMTRANS"  in cam_config_opts : 
   # if 'ESMT' in case:
   #    file.write(         ",'SPDT','SPDQ'")
   #    file.write(",'ZMMTU','ZMMTV','uten_Cu','vten_Cu' ")
   #    file.write(",'U_ESMT','V_ESMT'")
   #    file.write(",'UCONVMOM','VCONVMOM'")
   file.write('\n')

   # file.write(" fincl3    = 'OMEGA:I','DYN_OMEGA:I' ")

   #------------------------------
   # Prescribed aerosol settings
   #------------------------------
   if 'chem none' in cam_config_opts and compset not in ['FSP1V1','FSP2V1']  :
      prescribed_aero_path = '/atm/cam/chem/trop_mam/aero'
      prescribed_aero_file = 'mam4_0.9x1.2_L72_2000clim_c170323.nc'
      file.write(' use_hetfrz_classnuc = .false. \n')
      file.write(' aerodep_flx_type = \'CYCLICAL\' \n')
      file.write(' aerodep_flx_datapath = \''+din_loc_root+prescribed_aero_path+'\' \n')
      file.write(' aerodep_flx_file = \''+prescribed_aero_file+'\' \n')
      file.write(' aerodep_flx_cycle_yr = 01 \n')
      file.write(' prescribed_aero_type = \'CYCLICAL\' \n')
      file.write(' prescribed_aero_datapath=\''+din_loc_root+prescribed_aero_path+'\' \n')
      file.write(' prescribed_aero_file = \''+prescribed_aero_file+'\' \n')
      file.write(' prescribed_aero_cycle_yr = 01 \n')
   #------------------------------
   # Other namelist stuff
   #------------------------------
   # if npg==0: 
   file.write(' srf_flux_avg = 1 \n')              # Sfc flux smoothing (for SP stability)
   if ne<=60 : file.write(' dyn_npes = '+str(num_dyn)+' \n')   # limit dynamics tasks
   
   file.write(' se_fv_phys_remap_alg = 1 \n')

   # if '_TEST-LO_' in case : 
   #    file.write(' se_fv_phys_remap_alg = 0 \n')
   # else:
   #    file.write(' se_fv_phys_remap_alg = 1 \n')

   

   if 'checks-on' in case : 
      file.write(' state_debug_checks = .true. \n')
   else:
      file.write(' state_debug_checks = .false. \n')
   
   if '-HVx2_' in case and ne==30 : 
      file.write(" nu      = 2.0e15 \n")
      file.write(" nu_p    = 2.0e15 \n")
      file.write(" nu_div  = 5.0e15 \n")
      file.write(" hypervis_subcycle = 6 \n")
      file.write(' se_nsplit = 6 \n')

   if dtime == (20*60) : file.write(' rsplit = 2 \n')


   if ne==120 and npg==1 : file.write(' ncdata = \'/global/cscratch1/sd/whannah/acme_scratch/init_files/cami_aquaplanet_ne120_L72.nc\' \n')
   
   file.close()
      
   #-------------------------------------------------------
   # Set some run-time stuff
   #-------------------------------------------------------
   os.system('./xmlchange -file env_run.xml      ATM_NCPL='   +str(ncpl)   )
   os.system('./xmlchange -file env_run.xml      STOP_OPTION='+stop_opt    )
   os.system('./xmlchange -file env_run.xml      STOP_N='     +str(stop_n) )
   os.system('./xmlchange -file env_run.xml      RESUBMIT='   +str(resub)  )

   def xml_check_and_set(file_name,var_name,value):
      if var_name in open(file_name).read(): 
         os.system('./xmlchange -file '+file_name+' '+var_name+'='+str(value) )
   
   xml_check_and_set('env_workflow.xml','JOB_QUEUE',           queue)
   xml_check_and_set('env_batch.xml',   'JOB_QUEUE',           queue)
   xml_check_and_set('env_batch.xml',   'USER_REQUESTED_QUEUE',queue)
   xml_check_and_set('env_workflow.xml','JOB_WALLCLOCK_TIME',     walltime)
   xml_check_and_set('env_batch.xml',   'JOB_WALLCLOCK_TIME',     walltime)
   xml_check_and_set('env_batch.xml',   'USER_REQUESTED_WALLTIME',walltime)
   
   # os.system('./xmlchange -file env_workflow.xml JOB_QUEUE='           +queue )
   # os.system('./xmlchange -file env_batch.xml    JOB_QUEUE='           +queue )
   # os.system('./xmlchange -file env_batch.xml    USER_REQUESTED_QUEUE='+queue )
   # os.system('./xmlchange -file env_workflow.xml JOB_WALLCLOCK_TIME='     +walltime )
   # os.system('./xmlchange -file env_batch.xml    JOB_WALLCLOCK_TIME='     +walltime )
   # os.system('./xmlchange -file env_batch.xml    USER_REQUESTED_WALLTIME='+walltime )

   if continue_run :
      os.system('./xmlchange -file env_run.xml CONTINUE_RUN=TRUE ')   
   else:
      os.system('./xmlchange -file env_run.xml CONTINUE_RUN=FALSE ')
   #-------------------------------------------------------
   # Submit the run
   #-------------------------------------------------------
   os.system('./case.submit')

#---------------------------------------------------------------------------------------------------
# Done!
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n') # Print the case name again
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

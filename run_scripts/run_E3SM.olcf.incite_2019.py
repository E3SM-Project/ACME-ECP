#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
import os
import subprocess as sp
import datetime
newcase,config,build,clean,submit,continue_run = False,False,False,False,False,False

# directory info
case_dir = os.getenv("HOME")+'/E3SM/Cases/'
src_dir  = os.getenv("HOME")+'/E3SM/E3SM_SRC1/'

# clean        = True
newcase      = True
config       = True
build        = True
submit       = True
# continue_run = True

ne,npg         = 120,2
num_nodes      = 922             # ne4=>1, ne30=>15, ne120=>??
tasks_per_node = 18              # only for GPU cases
pcols          = 32              # should be slightly larger than #CRM/(total mpi tasks)
arch           = 'GPU'           # GPU / CPU
compset        = 'FC5AV1C-H01A' 

stop_opt,stop_n,resub = 'ndays',5,0

cld = 'SP1'    # ZM / SP1 / SP2
crm_nx,crm_ny,crm_dx = 64,1,1000

res = 'ne'+str(ne) if npg==0 else  'ne'+str(ne)+'pg'+str(npg)
cldc = '_'+cld+'_'+str(crm_nx)+'x'+str(crm_ny)+'_'+str(crm_dx)+'m' if 'SP' in cld else '_'+cld
if compset not in ['FSP1V1','FSP2V1'] : cldc = ''

# timestamp = '20191021'  # added more output
# timestamp = '20191022'  # reduced irad to 30 min (forgot to set set!), add land output, and use 64bit_data
# timestamp = '20191023' # reduce pcol => 32 and nx_rad => 4, and increase crm_accel_factor => 4
timestamp = '20191024' # set irad to 6 and change NTASK for non-ATM components

case = 'E3SM_TEST-INCITE_'+arch+'_'+res+'_'+compset+'_'+timestamp 

# Impose wall limits for Summit
if num_nodes>=  1: walltime =  '2:00'
if num_nodes>= 46: walltime =  '6:00'
if num_nodes>= 92: walltime = '12:00'
if num_nodes>=922: walltime = '24:00'

if 'TEST-INCITE' in case : walltime =  '2:00'
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

nlev_crm = 58     # limit CRM levels 
crm_dt   = 5      # CRM time step
dtime    = 5*60   # GCM physics time step
ncpl     = 86400 / dtime

if compset=='FSP1V1' : cld='SP1'
if compset=='FSP2V1' : cld='SP2'
#---------------------------------------------------------------------------------------------------
# Create new case
#---------------------------------------------------------------------------------------------------
if newcase :
   grid = res+'_'+res
   cmd = src_dir+'cime/scripts/create_newcase -case '+case_dir+case
   cmd = cmd + ' -compset '+compset+' -res '+grid
   if arch=='CPU' : cmd = cmd + ' -mach summit-cpu -compiler pgi    -pecount '+str(num_nodes*64)+'x1 '
   if arch=='GPU' : cmd = cmd + ' -mach summit     -compiler pgigpu -pecount '+str(num_nodes*tasks_per_node)+'x1 '
   os.system(cmd)

   # Change run directory to be next to bld directory
   os.chdir(case_dir+case+'/')
   os.system('./xmlchange -file env_run.xml   RUNDIR=\'$CIME_OUTPUT_ROOT/$CASE/run\' ' )
   
   # Change the max tasks per node
   os.system('./xmlchange -file env_mach_pes.xml      MAX_TASKS_PER_NODE='+str(tasks_per_node) )
   os.system('./xmlchange -file env_mach_pes.xml   MAX_MPITASKS_PER_NODE='+str(tasks_per_node) )
   
#---------------------------------------------------------------------------------------------------
# Configure
#---------------------------------------------------------------------------------------------------
os.chdir(case_dir+case+'/')
if config : 
   #----------------------------------------------------------------------------
   # Special CAM_CONFIG_OPTS for SP runs not using SP compsets
   #----------------------------------------------------------------------------
   if 'SP' in cld and compset not in ['FSP1V1','FSP2V1'] :
      # set options common to all SP setups
      
      cam_opt = ' -phys cam5 -use_SPCAM  -rad rrtmg -nlev 72 -microphys mg2 ' \
               +' -crm_nz '+str(nlev_crm) +' -crm_adv MPDATA '                \
               +' -crm_nx '+str(crm_nx)   +' -crm_ny '+str(crm_ny)            \
               +' -crm_dx '+str(crm_dx)   +' -crm_dt '+str(crm_dt)            \
               +' -crm_nx_rad 4 -crm_ny_rad 1 -bc_dep_to_snow_updates '
      # 1-moment microphysics
      if cld=='SP1': cam_opt = cam_opt + ' -SPCAM_microp_scheme sam1mom -chem none ' 
      # 2-moment microphysics
      if cld=='SP2': cam_opt = cam_opt + ' -SPCAM_microp_scheme m2005  '      \
                                       +' -chem linoz_mam4_resus_mom_soag '   \
                                       +' -rain_evap_to_coarse_aero '         \
                                       +' -bc_dep_to_snow_updates '
      cam_opt = cam_opt+' -cppdefs \' -DSP_DIR_NS -DSP_MCICA_RAD \' '

      os.system('./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS  -val  \"'+cam_opt+'\"' )
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------

   # Change pcols for GPU runs
   if arch=='GPU' : os.system('./xmlchange --append -file env_build.xml -id CAM_CONFIG_OPTS  -val  \" -pcols '+str(pcols)+' \" ' )

   # disable threading for SP
   if 'SP' in cld : os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1 ')
   
   # reduce task count for the non-atmos components
   ntask_atm = num_nodes*tasks_per_node
   os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val '+str(ntask_atm/4)+' ')
   os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val '+str(ntask_atm/4)+' ')
   os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val '+str(ntask_atm/4)+' ')
   os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 1 ')
   os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val 1 ')
   os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ESP -val 1 ')
   os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val 1 ')
   os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_IAC -val 1 ')

   os.system('./xmlchange ATM_PIO_NETCDF_FORMAT=\"64bit_data\" ')

   if clean : os.system('./case.setup --clean')
   os.system('./case.setup --reset')

#---------------------------------------------------------------------------------------------------
# Build
#---------------------------------------------------------------------------------------------------
if build : 
   # os.system('./xmlchange -file env_build.xml -id DEBUG -val TRUE ')   # enable debug mode
   if clean : os.system('./case.build --clean')
   os.system('./case.build')
#---------------------------------------------------------------------------------------------------
# Write the namelist options and submit the run
#---------------------------------------------------------------------------------------------------
if submit : 
   # Change inputdata from default due to permissions issue
   os.system('./xmlchange -file env_run.xml  DIN_LOC_ROOT=/gpfs/alpine/scratch/hannah6/cli115/inputdata ')
   
   #-------------------------------------------------------
   # Query some stuff about the case
   #-------------------------------------------------------
   (din_loc_root   , err) = sp.Popen('./xmlquery DIN_LOC_ROOT    -value', \
                                     stdout=sp.PIPE,shell=True,universal_newlines=True).communicate()
   (cam_config_opts, err) = sp.Popen('./xmlquery CAM_CONFIG_OPTS -value', \
                                     stdout=sp.PIPE,shell=True,universal_newlines=True).communicate()
   cam_config_opts = ' '.join(cam_config_opts.split())   # remove extra spaces to simplify string query
   #-------------------------------------------------------
   # Namelist options
   #-------------------------------------------------------
   nfile = 'user_nl_cam'
   file = open(nfile,'w') 
   #------------------------------
   # history output frequency and variables
   #------------------------------
   file.write(' nhtfrq = 0,-1,-3,-6 \n') 
   file.write(' mfilt = 1,120,40,20 \n')
   file.write(" avgflag_pertape = 'A','A','A','I' \n")     

   # Add dynamics grid variables to h0 files
   if npg>0 : file.write(" fincl1 = 'DYN_T','DYN_Q','DYN_U','DYN_OMEGA','DYN_PS' \n")
   # hourly 2D fields
   file.write(" fincl2 = 'PRECT','TMQ','LHFLX','SHFLX','TS','PS'")
   file.write(         ",'TGCLDLWP','TGCLDIWP','SWCF','LWCF'")    
   file.write(         ",'FLNT','FSNT','FSNS','FLNS'")            # TOA rad fields
   file.write(         ",'FSDS','FSDSC','FLDS','FLDSC'")          # SFC rad for CRF @ sfc
   file.write(         ",'TREFHT','QREFHT','U10'")                # for diurnal cycle analysis
   file.write(         ",'UBOT','VBOT','Z500','Z300','PSL'")      # for TC tracking 
   file.write("\n")
   # 3-hourly 3D fields for budget terms 
   file.write(" fincl3 = 'PS','T','Q','Z3','U','V','OMEGA','CLDLIQ','CLDICE'")
   file.write(         ",'SPTVFLUX','SPQVFLUX','SPTKE','SPTKES' ")
   file.write("\n")
   # CRM-fields near SGP site (halo of 9 crms for ne120np4)
   file.write(" fincl4 = 'CRM_U:I','CRM_W:I','CRM_QV:I'")
   file.write(         ",'CRM_QC:I','CRM_QI:I','CRM_QPC:I','CRM_QPI:I'")
   file.write("\n")
   file.write(" fincl4lonlat = '97.8w:97.2w_36.27n:36.75n' \n")

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
   # Other atm namelist stuff
   #------------------------------
   # num_dyn = ne*ne*6
   # file.write(' dyn_npes = '+str(num_dyn)+' \n')   # limit dynamics tasks 
   # file.write(' srf_flux_avg = 1 \n')              # Sfc flux smoothing (for SP stability)

   # Write initialization files at the end of each submission
   file.write(' inithist = \'ENDOFRUN\' \n')

   # hig-order physgrid mapping algorithm
   file.write(' se_fv_phys_remap_alg = 1 \n')

   # mean-state acceleration
   file.write(' use_crm_accel    = .true. \n')
   file.write(' crm_accel_uv     = .true. \n')
   file.write(' crm_accel_factor = 4 \n')

   # radiation every 30 minutes
   file.write(' iradlw = 6 \n')
   file.write(' iradsw = 6 \n')

   # file.write(' state_debug_checks = .true. \n')
   file.close()
   #------------------------------
   # new land initial condition file (not in inputdata yet) 
   #------------------------------
   nfile = 'user_nl_clm'
   file = open(nfile,'w') 
   # # file.write(' finidat = \'/gpfs/alpine/scratch/hannah6/cli115/init_files/clmi.ICLM45BC.ne30_ne30.d0241119c.clm2.r.nc\' \n')
   # file.write(' finidat = \'/gpfs/alpine/scratch/hannah6/cli115/init_files/???????\' \n')

   file.write(' hist_nhtfrq = 0, -1, -24 \n')
   file.write(" hist_fincl2 = 'EFLX_SOIL_GRND', 'FCEV', 'FCTR', 'FGEV' ")
   file.write(              ",'FSH', 'FSH_G', 'FSH_V' ")
   file.write(              ",'QINTR', 'QDRIP', 'QH2OSFC', 'QTOPSOIL', 'QINFL' ")
   file.write(              ",'QDRAI', 'QOVER', 'QRUNOFF' ")
   file.write(              ", 'SOILWATER_10CM', 'TSOI_10CM' ")
   file.write(              ",'SWdown', 'SWup', 'LWdown', 'LWup'")
   file.write('\n')
   file.write(" hist_fincl3 = 'H2OCAN', 'H2OSNO', 'H2OSOI', 'WA', 'TWS', 'ZWT', 'TSOI' ")
   file.write('\n')
   file.close()
   
   #-------------------------------------------------------
   # Set some run-time stuff
   #-------------------------------------------------------
   os.system('./xmlchange -file env_run.xml      ATM_NCPL='          +str(ncpl)   )
   os.system('./xmlchange -file env_run.xml      STOP_OPTION='       +stop_opt    )
   os.system('./xmlchange -file env_run.xml      STOP_N='            +str(stop_n) )
   os.system('./xmlchange -file env_run.xml      RESUBMIT='          +str(resub)  )
   os.system('./xmlchange -file env_workflow.xml JOB_WALLCLOCK_TIME='+walltime    )

   if ne==120 and npg==2 : os.system('./xmlchange -file env_run.xml      EPS_AGRID=1e-11' )

   if continue_run :
      os.system('./xmlchange -file env_run.xml CONTINUE_RUN=TRUE ')   
   else:
      os.system('./xmlchange -file env_run.xml CONTINUE_RUN=FALSE ')
   #-------------------------------------------------------
   # Submit the run
   #-------------------------------------------------------
   os.system('./case.submit')

#---------------------------------------------------------------------------------------------------
# Print the case name again
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n') 
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

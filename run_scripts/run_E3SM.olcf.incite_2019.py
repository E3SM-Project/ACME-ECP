#!/usr/bin/env python
# script for running SP-E3SM simulations using the 2019 INICTE allocation (CLI115)
# Branch for this campaign: https://github.com/E3SM-Project/ACME-ECP/tree/whannah/incite_2019
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
# newcase      = True
# config       = True
# build        = True
submit       = True
continue_run = True

stop_opt,stop_n,resub,walltime = 'nmonths',1,4,'5:00'
# stop_opt,stop_n,resub,walltime = 'ndays',183,0,'24:00'

ne,npg         = 120,2
compset        = 'FC5AV1C-H01A' 
crm_nx         = 64
crm_dx         = 1000
arch           = 'GPU'           # GPU / CPU
num_nodes      = 922             # ne4=>1, ne30=>15, ne120=>??
tasks_per_node = 18              # GPU=>18 / CPU=>64
pcols          = 32              # should be slightly larger than #CRM/(total mpi tasks)

phys = f'SP1_{crm_nx}x1_{crm_dx}m'
res  = f'ne{ne}' if npg==0 else  f'ne{ne}pg{npg}'

timestamp = '20191026'

case = '.'.join(['INCITE2019',arch,res,compset,phys,timestamp])

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

crm_accel_fac = 4  # CRM mean-state acceleration factor

dtime    = 5*60   # GCM physics time step
ncpl     = 86400 / dtime
#---------------------------------------------------------------------------------------------------
# Create new case
#---------------------------------------------------------------------------------------------------
if newcase :
   grid = res+'_'+res
   cmd = src_dir+'cime/scripts/create_newcase --case '+case_dir+case
   cmd = cmd + ' --compset '+compset+' --res '+grid
   cmd = cmd + ' --pecount '+str(num_nodes*tasks_per_node)+'x1 '
   if arch=='CPU' : cmd = cmd + ' --machine summit-cpu --compiler pgi    '
   if arch=='GPU' : cmd = cmd + ' --machine summit     --compiler pgigpu '
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
   cam_opt = ' -phys cam5 -use_SPCAM  -rad rrtmg -nlev 72 -microphys mg2 ' \
            +' -crm_nz 58 -crm_adv MPDATA -crm_dt 5 '                      \
            +' -crm_nx '+str(crm_nx)+' -crm_ny 1 -crm_dx '+str(crm_dx)     \
            +' -crm_nx_rad 4 -crm_ny_rad 1 -bc_dep_to_snow_updates '       \
            + ' -SPCAM_microp_scheme sam1mom -chem none '                  \
            +' -cppdefs \' -DSP_DIR_NS -DSP_MCICA_RAD \' '
   os.system('./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS  -val  \"'+cam_opt+'\"' )
   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   # Change pcols for GPU runs
   if arch=='GPU' : os.system('./xmlchange --append -file env_build.xml -id CAM_CONFIG_OPTS  -val  \" -pcols '+str(pcols)+' \" ' )

   # disable threading for SP
   os.system('./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1 ')
   
   # reduce task count for the non-atmos components
   ntask_atm = num_nodes*tasks_per_node
   os.system('./xmlchange -file env_mach_pes.xml NTASKS_OCN='+str(ntask_atm/4)+' ')
   os.system('./xmlchange -file env_mach_pes.xml NTASKS_ICE='+str(ntask_atm/4)+' ')
   os.system('./xmlchange -file env_mach_pes.xml NTASKS_LND='+str(ntask_atm/4)+' ')
   os.system('./xmlchange -file env_mach_pes.xml NTASKS_GLC=1,NTASKS_WAV=1,NTASKS_ESP=1,NTASKS_ROF=1,NTASKS_IAC=1 ')

   # 64_data format is needed for ne120 output
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
   (din_loc_root, err) = sp.Popen('./xmlquery DIN_LOC_ROOT -value',  \
                                    stdout=sp.PIPE,shell=True,       \
                                    universal_newlines=True).communicate()
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
   # Write initialization files at the end of each submission
   file.write(' inithist = \'ENDOFRUN\' \n')

   # hig-order physgrid mapping algorithm
   file.write(' se_fv_phys_remap_alg = 1 \n')

   # mean-state acceleration
   file.write(' use_crm_accel    = .true. \n')
   file.write(' crm_accel_uv     = .true. \n')
   file.write(' crm_accel_factor = '+str(crm_accel_fac)+' \n')

   # radiation every 30 minutes
   file.write(' iradlw = 6 \n')
   file.write(' iradsw = 6 \n')

   # Adjust dycore time step to keep 1.25 min dynamics time step [i.e. 15/(6*2) = 5/(2*2) = 1.25 ]
   file.write(" rsplit    = 2 \n")  # default is 2 
   file.write(" se_nsplit = 2 \n")  # default is 6

   # file.write(' state_debug_checks = .true. \n')
   file.close()
   #------------------------------
   # new land initial condition file (not in inputdata yet) 
   #------------------------------
   nfile = 'user_nl_clm'
   file = open(nfile,'w') 
   # file.write(' finidat = \'/gpfs/alpine/scratch/hannah6/cli115/init_files/clmi.ICLM45BC.ne30_ne30.d0241119c.clm2.r.nc\' \n')
   file.write(' finidat = \'/gpfs/alpine/scratch/hannah6/cli115/init_files/E3SM_PG-LAND-SPINUP_ne120pg2_FC5AV1C-H01A_00.clm2.r.0004-02-25-00000.nc\' \n')
   # file.write(' finidat = \'/ccs/proj/cli115/hannah6/init_files/E3SM_PG-LAND-SPINUP_ne120pg2_FC5AV1C-H01A_00.clm2.r.0004-02-25-00000.nc\' \n')
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
   os.system('./xmlchange -file env_workflow.xml PROJECT=cli115')

   # Restart Frequency
   os.system('./xmlchange -file env_run.xml      REST_OPTION='       +stop_opt    )
   os.system('./xmlchange -file env_run.xml      REST_N='            +str(stop_n) )
   # os.system('./xmlchange -file env_run.xml      REST_OPTION=nmonths')
   # os.system('./xmlchange -file env_run.xml      REST_N=2')

   # An alternate grid checking threshold is needed for ne120pg2 (still not sure why...)
   if ne==120 and npg==2 : os.system('./xmlchange -file env_run.xml  EPS_AGRID=1e-11' )

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

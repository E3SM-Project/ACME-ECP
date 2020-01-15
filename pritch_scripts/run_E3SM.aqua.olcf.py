#!/usr/bin/env python
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
import os
import subprocess as sp
newcase,config,build,clean,submit,continue_run = False,False,False,False,False,False

project = 'cli115'

# directory info
case_dir = os.getenv("HOME")+'/E3SM/Cases/'
src_dir  = os.getenv("HOME")+'/E3SM/E3SM_SRC1/'

# clean        = True
newcase      = True
config       = True
build        = True
submit       = True
# continue_run = True

stop_opt,stop_n,resub = 'ndays',1,0

ne,npg         = 30,2
compset        = 'F-EAMv1-AQP1'
crm_nx         = 64               # <<< change this one!
crm_ny         = 64
crm_dx         = 1000
arch           = 'GPU'           # GPU / CPU
num_nodes      = 15              # ne4=>1, ne30=>15, ne120=>922+
tasks_per_node = 18              # GPU=>18 / CPU=>64
pcols          = 32              # should be slightly larger than #CRM/(total mpi tasks)

if crm_nx==64 and crm_ny==64 : num_nodes,pcols   = 25,2
if crm_nx==16 and crm_ny==64 : num_nodes,pcols   = 25,20

if crm_nx==1 : 
   crm_nx = crm_ny
   crm_ny = 1
   phys = f'SP1_{crm_nx}x{crm_ny}R_{crm_dx}m'
else:
   phys = f'SP1_{crm_nx}x{crm_ny}_{crm_dx}m'
res  = f'ne{ne}' if npg==0 else  f'ne{ne}pg{npg}'

case = '.'.join(['3D_CONV_TEST',arch,res,compset,phys])
# case = '.'.join(['E3SM.3D_TEST',arch,res,compset,phys])


# Impose wall limits for Summit
if num_nodes>=  1: walltime =  '2:00'
if num_nodes>= 46: walltime =  '6:00'
if num_nodes>= 92: walltime = '12:00'
if num_nodes>=922: walltime = '24:00'

# case = case+'.debug-on'

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

crm_accel_fac = 3  # CRM mean-state acceleration factor

dtime    = 20*60   # GCM physics time step
ncpl     = 86400 / dtime
#---------------------------------------------------------------------------------------------------
# Create new case
#---------------------------------------------------------------------------------------------------
if newcase :
   grid = res+'_'+res
   cmd = src_dir+'cime/scripts/create_newcase -case '+case_dir+case
   cmd = cmd + ' -compset '+compset+' -res '+grid
   if arch=='CPU' : cmd = cmd + ' -mach summit-cpu -compiler pgi    -pecount '+str(num_nodes*84)+'x1 '
   if arch=='GPU' : cmd = cmd + ' -mach summit     -compiler pgigpu -pecount '+str(num_nodes*36)+'x1 '
   os.system(cmd)

   #-------------------------------------------------------
   # Change run directory to be next to bld directory
   #-------------------------------------------------------
   os.chdir(case_dir+case+'/')
   os.system('./xmlchange -file env_run.xml   RUNDIR=\'$CIME_OUTPUT_ROOT/$CASE/run\' ' )
#---------------------------------------------------------------------------------------------------
# Configure
#---------------------------------------------------------------------------------------------------
os.chdir(case_dir+case+'/')
if config : 
   #----------------------------------------------------------------------------
   # Special CAM_CONFIG_OPTS for SP runs not using SP compsets
   #----------------------------------------------------------------------------   
   if f'SP1_{crm_nx}x{crm_ny}R' in case:
      cpp_defs = cpp_defs + ' -DSP_DIR_NS '
   else:
      cpp_defs = ' -DSP_MCICA_RAD '
   
   cam_opt = ' -phys cam5 -use_SPCAM  -rad rrtmg -nlev 72 -microphys mg2 ' \
            +' -crm_nz 58 -crm_adv MPDATA -crm_dt 5 '                      \
            +f' -crm_nx {crm_nx} -crm_ny {crm_ny} -crm_dx {crm_dx} '       \
            +' -crm_nx_rad 1 -crm_ny_rad 1 -bc_dep_to_snow_updates '       \
            + ' -SPCAM_microp_scheme sam1mom -chem none '                  \
            +f' -cppdefs \' {cpp_defs} \' '   
   if 'AQP' in compset : cam_opt = cam_opt+' -aquaplanet '
   os.system(f'./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS  -val  \"{cam_opt}\" ' )
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
   if 'debug-on' in case : os.system('./xmlchange -file env_build.xml -id DEBUG -val TRUE ')
   if clean : os.system('./case.build --clean')
   os.system('./case.build')
#---------------------------------------------------------------------------------------------------
# Write the namelist options and submit the run
#---------------------------------------------------------------------------------------------------
if submit : 
   # Change inputdata from default due to permissions issue
   os.system('./xmlchange -file env_run.xml  DIN_LOC_ROOT=/gpfs/alpine/scratch/hannah6/cli115/inputdata ')
   
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
   file.write(' nhtfrq    = 0,-3 \n') 
   file.write(' mfilt     = 1,8 \n')     
   # if npg>0 : file.write(" fincl1    = 'DYN_T','DYN_Q','DYN_U','DYN_OMEGA'")
   file.write(" fincl2    = 'PS','TS'")
   file.write(             ",'PRECT','TMQ'")
   file.write(             ",'LHFLX','SHFLX'")             # surface fluxes
   file.write(             ",'TGCLDLWP','TGCLDIWP'")       # liq and ice water path
   file.write(             ",'FSNT','FLNT'")               # Net TOM heating rates
   file.write(             ",'FLNS','FSNS'")               # Surface rad for total column heating
   file.write(             ",'FSNTC','FLNTC'")             # clear sky heating rates for CRE
   # file.write(             ",'LWCF','SWCF'")               # cloud radiative foricng
   # file.write(             ",'TAUX','TAUY'")               # surface stress
   # file.write(             ",'CLDLOW','CLDMED','CLDHGH','CLDTOT' ")
   # file.write(             ",'CLOUD','CLDLIQ','CLDICE'")   # 3D cloud fields
   # file.write(             ",'T','Q','Z3' ")               # 3D thermodynamic budget components
   # file.write(             ",'U','V','OMEGA'")             # 3D velocity components
   # file.write(             ",'QRL','QRS'")                 # 3D radiative heating profiles
   # file.write(             ",'DYN_T','DYN_Q','DYN_OMEGA'")
   # if 'SP' in cld :
   #    file.write(         ",'SPDT','SPDQ'")               # CRM heating/moistening tendencies
   #    file.write(         ",'SPTLS','SPQTLS' ")           # CRM large-scale forcing
      # file.write(         ",'SPQPEVP','SPMC'")            # CRM rain evap and total mass flux
      # file.write(         ",'SPMCUP','SPMCDN'")           # CRM saturated mass fluxes
      # file.write(         ",'SPMCUUP','SPMCUDN'")         # CRM unsaturated mass fluxes
      # file.write(         ",'SPTKE'")
      # if any(x in cam_config_opts for x in ["SP_ESMT","SP_USE_ESMT","SPMOMTRANS"]) : 
      #    file.write(",'ZMMTU','ZMMTV','uten_Cu','vten_Cu' ")
      # if "SP_USE_ESMT" in cam_config_opts : file.write(",'U_ESMT','V_ESMT'")
      # if "SPMOMTRANS"  in cam_config_opts : 
   # if 'ESMT' in case:
   #    file.write(",'ZMMTU','ZMMTV','uten_Cu','vten_Cu' ")
   #    file.write(",'U_ESMT','V_ESMT'")
   #    file.write(",'UCONVMOM','VCONVMOM'")
   file.write('\n')

   # file.write(" fincl3    = 'OMEGA:I','DYN_OMEGA:I' \n ")

   #------------------------------
   # Prescribed aerosol settings
   #------------------------------
   if 'chem none' in cam_config_opts and compset not in ['FSP1V1','FSP2V1'] :
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
   # file.write(' srf_flux_avg = 1 \n')              # Sfc flux smoothing (for SP stability)

   file.write(' se_fv_phys_remap_alg = 1 \n')

   # mean-state acceleration
   file.write(' use_crm_accel    = .true. \n')
   file.write(' crm_accel_uv     = .true. \n')
   file.write(' crm_accel_factor = '+str(crm_accel_fac)+' \n')

   # adjust dycore time step
   file.write(" rsplit    = 2 \n")  # default is 2 
   file.write(" se_nsplit = 2 \n")  # default is 6
   
   file.close()

   #------------------------------
   # new land initial condition file (not in inputdata yet) 
   #------------------------------
   # nfile = 'user_nl_clm'
   # file = open(nfile,'w') 
   # if res=='ne30' : file.write(' finidat = \'/gpfs/alpine/scratch/hannah6/cli115/init_files/clmi.ICLM45BC.ne30_ne30.d0241119c.clm2.r.nc\' \n')
   # file.close()
   
   #-------------------------------------------------------
   # Set some run-time stuff
   #-------------------------------------------------------
   

   os.system('./xmlchange -file env_run.xml      ATM_NCPL='          +str(ncpl)   )
   os.system('./xmlchange -file env_run.xml      STOP_OPTION='       +stop_opt    )
   os.system('./xmlchange -file env_run.xml      STOP_N='            +str(stop_n) )
   os.system('./xmlchange -file env_run.xml      RESUBMIT='          +str(resub)  )
   os.system('./xmlchange -file env_workflow.xml JOB_WALLCLOCK_TIME='+walltime    )
   os.system('./xmlchange -file env_workflow.xml PROJECT='           +project)

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

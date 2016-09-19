module setparm_mod

contains

subroutine setparm
	
!       initialize parameters:

use vars
!use micro_params
use params
!use isccp, only : isccp_zero
!use isccpTables, only : isccp_tables_init
use microphysics, only: micro_setparm
use grid, only: doclubb, doclubbnoninter

implicit none
	
integer icondavg, ierr

!NAMELIST /PARAMETERS/ dodamping, doupperbound, docloud, doprecip, &
!                dolongwave, doshortwave, dosgs, &
!                docoriolis, dosurface, dolargescale, doradforcing, &
!		nadams,fluxt0,fluxq0,tau0,tabs_s,z0,tauls,nelapse, &
!		dt, dx, dy, fcor, ug, vg, nstop, caseid, &
!		nstat, nstatfrq, nprint, nrestart, doradsimple, &
!		nsave3D, nsave3Dstart, nsave3Dend, dosfcforcing, &
!		donudging_uv, donudging_tq, dosmagor, doscalar,  &
!		timelargescale, longitude0, latitude0, day0, nrad, &
!		CEM,LES,OCEAN,LAND,SFC_FLX_FXD,SFC_TAU_FXD, soil_wetness, &
!                doensemble, nensemble, doxy, dowallx, dowally, &
!                nsave2D, nsave2Dstart, nsave2Dend, qnsave3D, & 
!                docolumn, save2Dbin, save2Davg, save3Dbin, &
!                save2Dsep, save3Dsep, dogzip2D, dogzip3D, restart_sep, &
!	        doseasons, doperpetual, doradhomo, dosfchomo, doisccp, &
!	        dodynamicocean, ocean_type, &
!		dosolarconstant, solar_constant, zenith_angle, rundatadir, &
!                dotracers, output_sep, perturb_type, &
!                doSAMconditionals, dosatupdnconditionals, &
!                doscamiopdata, iopfile, dozero_out_day0, &
!                nstatmom, nstatmomstart, nstatmomend, savemomsep, savemombin, &
!                nmovie, nmoviestart, nmovieend, nrestart_skip, &
!                bubble_x0,bubble_y0,bubble_z0,bubble_radius_hor, &
!               bubble_radius_ver,bubble_dtemp,bubble_dq, dosmoke, &
!               doclubb, doclubbnoninter, doclubb_sfc_fluxes, & ! added by dschanen UWM
!               docam_sfc_fluxes           ! added by mhwang 
	
!----------------------------
!  Set defaults:

	dodamping 	= .false.
	doupperbound   	= .false.
	docloud   	= .false.
	doprecip        = .false.
	dolongwave	= .false.
	doshortwave	= .false.
	doradsimple 	= .false.
	dosgs		= .false.
	dosmagor	= .false.
	doscalar	= .false.
	dosubsidence	= .false.
	docoriolis	= .false.
	dosurface	= .false.
	dolargescale    = .false.
	doradforcing    = .false.
	dosfcforcing    = .false.
	donudging_uv	= .false.
	donudging_tq	= .false.
	doensemble	= .false.
	doxy    	= .false.
	dowallx    	= .false.
	dowally    	= .false.
	docup		= .false.
	docolumn	= .false.
	doseasons	= .false.
	doperpetual	= .false.
	doradhomo	= .false.
	dosfchomo	= .false.
	doisccp		= .false.
	dodynamicocean 	= .false.
        dosolarconstant = .false.
        dotracers       = .false.
        dosmoke         = .false.
        ! Begin code addition for CLUBB
        ! This code is a part of a larger set of code that declares the
        ! doclubb flag so that if a user includes
        ! doclubb in a namelist but compiles without the -DCLUBB flag,
        ! we can return an intelligible error message    
        doclubb         = .false.
        doclubbnoninter = .false.
        doclubb_sfc_fluxes = .false.
        docam_sfc_fluxes = .false.
        nclubb          = 1 
        ! End of code addition for CLUBB
	CEM		= .false.
	LES		= .false.
	OCEAN		= .false.
	LAND		= .false.
	SFC_FLX_FXD	= .false.
	SFC_TAU_FXD	= .false.
		
	nadams		= 3
	dt		= 0
	dx		= 0
	dy		= 0
	longitude0	= 0.
	latitude0	= 0.
	fcor	        = -999.
	day0		= 0.
	nrad		= 1
	ug		= 0.
	vg		= 0.
	fluxt0		= 0.
	fluxq0		= 0.
	tau0		= 0.
	z0		= 0.035
        soil_wetness    = 1.
	timelargescale  = 0.
	tauls		= 7200.
	tabs_s 		= 0.
	nstop 		= 0
	nelapse		= 999999999
	caseid		= 'les00000'
	nstat		= 1000
	nstatfrq	= 50
	nprint		= 1000
	nrestart 	= 0
	restart_sep 	= .false.
        nrestart_skip   = 0
	output_sep 	= .false.
	save3Dbin	= .false.
	save2Dsep	= .false.
	save3Dsep	= .false.
	nsave3D		= 1
	nsave3Dstart	= 99999999
	nsave3Dend	= 999999999
	dogzip2D	= .false.
	dogzip3D	= .false.
	save2Dbin	= .false.
	save2Davg	= .false.
	nsave2D		= 1
	nsave2Dstart	= 99999999
	nsave2Dend	= 999999999
        savemombin      = .false.
        savemomsep      = .false.
        nstatmom        = 1
        nstatmomstart    = 99999999
        nstatmomend      = 999999999
        nmovie           = 1
        nmoviestart      = 99999999
        nmovieend        = 999999999
	nensemble	= 0
 	qnsave3D	= 0.
	ocean_type	= 0 
        rundatadir      = './RUNDATA'
        perturb_type    = 0
        bubble_x0       = 0.
        bubble_y0       = 0.
        bubble_z0       = 0.
        bubble_radius_hor=0. 
        bubble_radius_ver=0. 
        bubble_dtemp     =0. 
        bubble_dq        =0.


        ! Specify solar constant and zenith angle for perpetual insolation.
        ! Note that if doperpetual=.true. and dosolarconstant=.false.
        ! the insolation will be set to the daily-averaged value on day0.

        solar_constant = 685. ! Values from Tompkins & Craig, J. Climate (1998)
        zenith_angle = 51.7

        !bloss: add option for core updraft, core downdraft conditional statistics
        doSAMconditionals = .true.

        !bloss: add option for additional conditional averages:
        !        cloudy updrafts, cloudy downdrafts and cloud-free.
        dosatupdnconditionals = .false.
        ! Allow sounding, forcing and surface data to be read in
        ! from a SCAM IOP input file in netcdf format.
        doscamiopdata = .false.
        iopfile = trim(case) // '.nc' ! default name: CASENAME.nc
        dozero_out_day0 = .false.

!----------------------------------
!  Read namelist variables from the standard input:
!------------

!open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
!read (55,PARAMETERS)
!close(55)

	doprecip        = .true.
	dosgs		= .true.
	dosmagor	= .true.
	dosurface	= .true.
	dodamping 	= .true.
	dt		= CRM_DT
	dx		= CRM_DX
	dy		= CRM_DY
#ifndef CLUBB_CRM
        doclubb         = .false.   ! then docloud must be .true.
        docloud         = .true.
#else
!        doclubb         = .false.   ! then docloud must be .true.
!        docloud         = .true.
        doclubb         = .true.    ! then docloud must be .false.
        docloud         = .false.
        doclubbnoninter = .false.
        doclubb_sfc_fluxes = .false.
        docam_sfc_fluxes = .true.   ! update variables in cam, neither in sam nor in clubb +++mhwang
        nclubb          = 3 

#ifdef sam1mom
! for sam1mom, nclubb needs to be 1. 
! see comments in ./MICRO_SAM1MOM/microphysics.F90
        nclubb          = 1 
#endif

#endif
        rank            = 0   ! in MMF model, rank = 0
!------------------------------------
!  Set parameters 


        ! Allow only special cases for separate output:

        output_sep = output_sep.and.RUN3D
        if(output_sep)  save2Dsep = .true.

	if(RUN2D) dy=dx

	if(RUN2D.and.YES3D.eq.1) then
	  print*,'Error: 2D run and YES3D is set to 1. Exitting...'
	  call task_abort()
	endif
	if(RUN3D.and.YES3D.eq.0) then
	  print*,'Error: 3D run and YES3D is set to 0. Exitting...'
	  call task_abort()
	endif
#ifdef CLUBB_CRM
        if ( dx >= 1000. .and. LES ) then
          print*,'Error: Horizonatal grid spacing is >= 1000. meters'
          print*,'but LES is true.  Use CEM mode for coarse resolutions.'
          call task_abort()
        end if
#endif

	pi = acos(-1.)
	if(fcor.eq.-999.) fcor= 4*pi/86400.*sin(latitude0*pi/180.)
	fcorz = sqrt(4.*(2*pi/(3600.*24.))**2-fcor**2)	  
	coszrs = 0.637 ! default mean solar zenith angle
	
	if(ny.eq.1) dy=dx

	na = 1
	nb = 2
	nc = 3
	nstep = 0
        time = 0.
	dtn = dt

	notopened2D = .true.
	notopened3D = .true.

!        call isccp_tables_init()   ! initialize isccp tables
!        call isccp_zero()
        call micro_setparm() ! read in microphysical options from prm file.

        if(dosmoke) then
           epsv=0.
        else    
           epsv=0.61
        endif   

        if(navgmom_x.lt.0.or.navgmom_y.lt.0) then  
            nstatmom        = 1
            nstatmomstart    = 99999999
            nstatmomend      = 999999999
        end if

        masterproc = rank.eq.0
          
end subroutine setparm
end module setparm_mod

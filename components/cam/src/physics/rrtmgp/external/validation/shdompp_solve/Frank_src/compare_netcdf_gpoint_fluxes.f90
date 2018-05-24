 ! compare_netcdf_gpoint_fluxes.f90
 !   Reads the band and g-point downward and upward flux arrays from two
 ! RRTMGP netcdf files and prints out the down and up flux profiles 
 ! for the two files and the differences for a specified range of columns
 ! and specified band and g-point number.
 !
 ! gfortran -O2 -o compare_netcdf_gpoint_fluxes compare_netcdf_gpoint_fluxes.f90
 !   -lnetcdf -lnetcdff -I/usr/local/include/ -L/usr/local/lib/


SUBROUTINE Read_netcdf_size (infile, Nlev, Ncol, Nband, Ngpt)
 ! Returns the size of the dimensions
  USE netcdf
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: infile
  INTEGER, INTENT(OUT) :: Nlev, Ncol, Nband, Ngpt
  INTEGER :: status(2), ncid, dimid

   ! Get ID of netcdf file   
  status(1) = NF90_OPEN (infile,0,ncid)
  if (status(1) /= 0) CALL ErrorNC('error opening file: ',infile)

  status(1) = NF90_INQ_DIMID (ncid, 'lev', dimid)
  status(2) = NF90_INQUIRE_DIMENSION (ncid, dimid, len=Nlev)
  if (ANY(status(:) /= 0)) CALL ErrorNC('error getting length of lev dimension: ',infile)

  status(1) = NF90_INQ_DIMID (ncid, 'col', dimid)
  status(2) = NF90_INQUIRE_DIMENSION (ncid, dimid, len=Ncol)
  if (ANY(status(:) /= 0)) CALL ErrorNC('error getting length of col dimension: ',infile)

  status(1) = NF90_INQ_DIMID (ncid, 'band', dimid)
  status(2) = NF90_INQUIRE_DIMENSION (ncid, dimid, len=Nband)
  if (ANY(status(:) /= 0)) CALL ErrorNC('error getting length of band dimension: ',infile)

  status(1) = NF90_INQ_DIMID (ncid, 'gpt', dimid)
  status(2) = NF90_INQUIRE_DIMENSION (ncid, dimid, len=Ngpt)
  if (ANY(status(:) /= 0)) CALL ErrorNC('error getting length of gpt dimension: ',infile)

  status(1) = NF90_CLOSE (ncid)
  if (status(1) /= 0) CALL ErrorNC('error closing file: ',infile)
END SUBROUTINE Read_netcdf_size


SUBROUTINE Read_netcdf_fluxes (infile, Ncol, Nlev, Nband, Ngpt, &
                               preslev, band_waveno, &
                               band_flux_dn, band_flux_up, &
                               gpt_flux_dn, gpt_flux_up)
 ! Opens the netcdf file and reads the variables output from this routine:
 !  pressure at layer edges (p_lev, preslev, mb), 
 !  wavenumber ranges for bands (band_lims_wvn, band_waveno, cm^-1-).,
 !  flux down by band (band_flux_dn, W/m^2), flux up by band (band_flux_up),
 !  g-point flux down (gpt_flux_dn), and g-point flux up (gpt_flux_up).
  USE netcdf
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: infile
  INTEGER, INTENT(IN) :: Ncol, Nlev, Nband, Ngpt
  REAL,    INTENT(OUT) :: preslev(Ncol,Nlev), band_waveno(Nband,2)
  REAL,    INTENT(OUT) :: band_flux_dn(Nband,Ncol,Nlev), band_flux_up(Nband,Ncol,Nlev)
  REAL,    INTENT(OUT) :: gpt_flux_dn(Ncol,Nlev,Ngpt), gpt_flux_up(Ncol,Nlev,Ngpt)
  INTEGER :: status(3), ncid, arrayid

   ! Get ID of netcdf file   
  status(1) = NF90_OPEN (infile,0,ncid)
  if (status(1) /= 0) CALL ErrorNC('error opening file: ',infile)

  status(1) = NF90_INQ_VARID (ncid, 'p_lev', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, preslev)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting p_lev: ',infile)
  preslev(:,:) = 0.01*preslev(:,:)

  status(1) = NF90_INQ_VARID (ncid, 'band_lims_wvn', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, band_waveno)
  if (ANY(status(1:2) /= 0)) then
    band_waveno(:,:) = 0.0
  endif

  status(1) = NF90_INQ_VARID (ncid, 'band_flux_dn', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, band_flux_dn)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting band_flux_dn: ',infile)

  status(1) = NF90_INQ_VARID (ncid, 'band_flux_up', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, band_flux_up)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting band_flux_up: ',infile)

  status(1) = NF90_INQ_VARID (ncid, 'gpt_flux_dn', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, gpt_flux_dn)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting gpt_flux_dn: ',infile)

  status(1) = NF90_INQ_VARID (ncid, 'gpt_flux_up', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, gpt_flux_up)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting gpt_flux_up: ',infile)

  status(1) = NF90_CLOSE (ncid)
  if (status(1) /= 0) CALL ErrorNC('error closing file: ',infile)
END SUBROUTINE Read_netcdf_fluxes


SUBROUTINE ErrorNC (errstring, filename)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: errstring, filename
  PRINT *, 'NETCDF ',TRIM(errstring), ' ', TRIM(filename)
  STOP
END SUBROUTINE ErrorNC



PROGRAM compare_netcdf_gpoint_fluxes
  IMPLICIT NONE
  CHARACTER(LEN=96) :: fluxfile1, fluxfile2, outfile
  INTEGER :: Nlev, Ncol, Nband, Ngpt, Nlev2, Ncol2, Nband2, Ngpt2
  INTEGER :: icol1, icol2, ibnd, igpt, icol, ilev
  REAL, ALLOCATABLE :: preslev1(:,:), preslev2(:,:)
  REAL, ALLOCATABLE :: band_waveno1(:,:), band_waveno2(:,:)
  REAL, ALLOCATABLE :: band_flux_dn1(:,:,:), band_flux_up1(:,:,:)
  REAL, ALLOCATABLE :: band_flux_dn2(:,:,:), band_flux_up2(:,:,:)
  REAL, ALLOCATABLE :: gpt_flux_dn1(:,:,:), gpt_flux_up1(:,:,:)
  REAL, ALLOCATABLE :: gpt_flux_dn2(:,:,:), gpt_flux_up2(:,:,:)
  REAL, ALLOCATABLE :: flux_dn_bias(:)
  REAL, ALLOCATABLE :: flux_up_bias(:)


  PRINT *, 'First flux netcdf input file name'
  READ (*,'(A)') fluxfile1
  WRITE (*,*) TRIM(fluxfile1)

  PRINT *, 'Second flux netcdf input file name'
  READ (*,'(A)') fluxfile2
  WRITE (*,*) TRIM(fluxfile2)

  PRINT *, 'Flux profile comparison output file name'
  READ (*,'(A)') outfile
  WRITE (*,*) TRIM(outfile)

  PRINT *, 'Range of column numbers to output flux comparisons'
  READ (*,*) icol1, icol2
  WRITE (*,*) icol1, icol2

  PRINT *, 'Band number to compare (0 for none)'
  READ (*,*) ibnd
  WRITE (*,*) ibnd

  PRINT *, 'g-point number to compare (0 for none)'
  READ (*,*) igpt
  WRITE (*,*) igpt


  CALL Read_netcdf_size (fluxfile1, Nlev, Ncol, Nband, Ngpt)
  CALL Read_netcdf_size (fluxfile2, Nlev2, Ncol2, Nband2, Ngpt2)
  IF (Nlev2 /= Nlev .OR. Ncol2 /= Ncol .OR. &
      Nband2 /= Nband .OR. Ngpt2 /= Ngpt) THEN
    PRINT *, 'Input netcdf files have different array sizes'
    STOP
  ENDIF

  ALLOCATE (preslev1(Ncol,Nlev), preslev2(Ncol,Nlev))
  ALLOCATE (band_waveno1(Nband,2), band_waveno2(Nband,2))
  ALLOCATE (band_flux_dn1(Nband,Ncol,Nlev), band_flux_up1(Nband,Ncol,Nlev))
  ALLOCATE (band_flux_dn2(Nband,Ncol,Nlev), band_flux_up2(Nband,Ncol,Nlev))
  ALLOCATE (gpt_flux_dn1(Ncol,Nlev,Ngpt), gpt_flux_up1(Ncol,Nlev,Ngpt))
  ALLOCATE (gpt_flux_dn2(Ncol,Nlev,Ngpt), gpt_flux_up2(Ncol,Nlev,Ngpt))

  CALL Read_netcdf_fluxes (fluxfile1, Ncol, Nlev, Nband, Ngpt, &
                       preslev1, band_waveno1,&
                       band_flux_dn1, band_flux_up1, gpt_flux_dn1, gpt_flux_up1)
  CALL Read_netcdf_fluxes (fluxfile2, Ncol, Nlev, Nband, Ngpt, &
                       preslev2, band_waveno2,&
                       band_flux_dn2, band_flux_up2, gpt_flux_dn2, gpt_flux_up2)
  IF (ANY(preslev1(:,:) /= preslev2(:,:))) THEN
    stop 'Pressure levels in two input files disagree'
  ENDIF
  IF (ANY(band_waveno1(:,:) /= band_waveno2(:,:)) .AND. &
       SUM(band_waveno1(:,:))>1.0 .AND. SUM(band_waveno2(:,:))>1.0) THEN
    stop 'Band wavenumber ranges in two input files disagree'
  ENDIF


  OPEN (UNIT=1, FILE=outfile, STATUS='UNKNOWN')
  WRITE (1,'(2(A,I2))') '! Flux comparisons for columns ',icol1,' to ',icol2
  WRITE (1,'(A,I2,A,I4,1X,I4,A,I3)') '! Band compared: ',ibnd, &
        ' (', NINT(band_waveno1(ibnd,:)),')    g-point compared: ',igpt
  WRITE (1,'(A,A)') '!  Input 1: ', TRIM(fluxfile1)
  WRITE (1,'(A,A)') '!  Input 2: ', TRIM(fluxfile2)
  WRITE (1,'(A)') '!                  Downwelling Flux                   Upwelling Flux'
  WRITE (1,'(A)') '! Pres         1          2        2 - 1         1          2        2 - 1'
  DO icol = icol1, icol2
    IF (ibnd > 0) THEN
      WRITE (1,'(A)') '!'
      WRITE (1,'(A,I3)') '! Band flux comparison for column: ',icol
      DO ilev = 1, Nlev
        WRITE (1,'(1X,F7.2,1X,3(1X,F10.6),1X,3(1X,F10.6))') preslev1(icol,ilev), &
           band_flux_dn1(ibnd,icol,ilev), band_flux_dn2(ibnd,icol,ilev), &
           band_flux_dn2(ibnd,icol,ilev)-band_flux_dn1(ibnd,icol,ilev), &
           band_flux_up1(ibnd,icol,ilev), band_flux_up2(ibnd,icol,ilev), &
           band_flux_up2(ibnd,icol,ilev)-band_flux_up1(ibnd,icol,ilev)
      ENDDO
    ENDIF
    IF (igpt > 0) THEN
      WRITE (1,'(A)') '!'
      WRITE (1,'(A,I3)') '! g-point flux comparison for column: ',icol
      DO ilev = 1, Nlev
        WRITE (1,'(1X,F7.2,1X,3(1X,F10.6),1X,3(1X,F10.6))') preslev1(icol,ilev), &
           gpt_flux_dn1(icol,ilev,igpt), gpt_flux_dn2(icol,ilev,igpt), &
           gpt_flux_dn2(icol,ilev,igpt)-gpt_flux_dn1(icol,ilev,igpt), &
           gpt_flux_up1(icol,ilev,igpt), gpt_flux_up2(icol,ilev,igpt), &
           gpt_flux_up2(icol,ilev,igpt)-gpt_flux_up1(icol,ilev,igpt)
      ENDDO
    ENDIF
  ENDDO
  CLOSE (1)

  DEALLOCATE (preslev1, preslev2)
  DEALLOCATE (gpt_flux_dn2, gpt_flux_up2, gpt_flux_dn1, gpt_flux_up1)
  DEALLOCATE (band_flux_dn2, band_flux_up2, band_flux_dn1, band_flux_up1)
END PROGRAM compare_netcdf_gpoint_fluxes

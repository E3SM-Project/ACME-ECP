 ! compare_netcdf_fluxes.f90
 !
 ! gfortran -O2 -o compare_netcdf_fluxes compare_netcdf_fluxes.f90
 !   -lnetcdf -lnetcdff -I/usr/local/include/ -L/usr/local/lib/


SUBROUTINE Read_netcdf_size (infile, Nlev, Ncol, Nband)
 ! Returns the size of the dimensions
  USE netcdf
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: infile
  INTEGER, INTENT(OUT) :: Nlev, Ncol, Nband
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
  if (ANY(status(:) /= 0)) CALL ErrorNC('error getting length of band  dimension: ',infile)

  status(1) = NF90_CLOSE (ncid)
  if (status(1) /= 0) CALL ErrorNC('error closing file: ',infile)
END SUBROUTINE Read_netcdf_size


SUBROUTINE Read_netcdf_fluxes (infile, Ncol, Nlev, Nband, preslev, band_waveno,&
                               flux_dn, flux_up, flux_net, &
                               band_flux_dn, band_flux_up, band_flux_net)
 ! Opens the netcdf file and reads the variables output from this routine:
 !  Pressure at layer edges (p_lev, preslev, mb), 
 !  Wavenumber ranges for bands (band_lims_wvn, band_waveno, cm^-1-).,
 !  Broadband flux down (flux_dn, W/m^2), Broadband flux up (flux_up, W/m^2),
 !  Broadband flux net (down-up) (flux_net, W/m^2), Flux down by band
 !  (band_flux_dn, W/m^2), Flux up by band (band_flux_up, W/m^2), and
 !  Flux net by band (down-up) (band_flux_net, W/m^2).
  USE netcdf
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: infile
  INTEGER, INTENT(IN) :: Ncol, Nlev, Nband
  REAL,    INTENT(OUT) :: preslev(Ncol,Nlev), band_waveno(Nband,2)
  REAL,    INTENT(OUT) :: flux_dn(Ncol,Nlev), flux_up(Ncol,Nlev), flux_net(Ncol,Nlev)
  REAL,    INTENT(OUT) :: band_flux_dn(Nband,Ncol,Nlev), band_flux_up(Nband,Ncol,Nlev)
  REAL,    INTENT(OUT) :: band_flux_net(Nband,Ncol,Nlev)
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

  status(1) = NF90_INQ_VARID (ncid, 'flux_dn', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, flux_dn)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting flux_dn: ',infile)

  status(1) = NF90_INQ_VARID (ncid, 'flux_up', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, flux_up)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting flux_up: ',infile)

  status(1) = NF90_INQ_VARID (ncid, 'flux_net', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, flux_net)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting flux_net: ',infile)

  status(1) = NF90_INQ_VARID (ncid, 'band_flux_dn', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, band_flux_dn)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting band_flux_dn: ',infile)

  status(1) = NF90_INQ_VARID (ncid, 'band_flux_up', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, band_flux_up)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting band_flux_up: ',infile)

  status(1) = NF90_INQ_VARID (ncid, 'band_flux_net', arrayid)
  status(2) = NF90_GET_VAR (ncid, arrayid, band_flux_net)
  if (ANY(status(1:2) /= 0)) CALL ErrorNC('error getting band_flux_net: ',infile)

  status(1) = NF90_CLOSE (ncid)
  if (status(1) /= 0) CALL ErrorNC('error closing file: ',infile)
END SUBROUTINE Read_netcdf_fluxes


SUBROUTINE ErrorNC (errstring, filename)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: errstring, filename
  PRINT *, 'NETCDF ',TRIM(errstring), ' ', TRIM(filename)
  STOP
END SUBROUTINE ErrorNC



PROGRAM compare_netcdf_fluxes
!  USE compare_netcdf_fluxes_sub
  IMPLICIT NONE
  CHARACTER(LEN=96) :: fluxfile1, fluxfile2, profoutfile, bandoutfile
  INTEGER :: Nlev, Ncol, Nband, Nlev2, Ncol2, Nband2
  INTEGER :: ilev, ibnd, compband
  REAL, ALLOCATABLE :: preslev1(:,:), preslev2(:,:)
  REAL, ALLOCATABLE :: band_waveno1(:,:), band_waveno2(:,:)
  REAL, ALLOCATABLE :: flux_dn1(:,:), flux_up1(:,:), flux_net1(:,:)
  REAL, ALLOCATABLE :: flux_dn2(:,:), flux_up2(:,:), flux_net2(:,:)
  REAL, ALLOCATABLE :: band_flux_dn1(:,:,:), band_flux_up1(:,:,:)
  REAL, ALLOCATABLE :: band_flux_dn2(:,:,:), band_flux_up2(:,:,:)
  REAL, ALLOCATABLE :: band_flux_net1(:,:,:), band_flux_net2(:,:,:)
  REAL, ALLOCATABLE :: preslevavg(:)
  REAL, ALLOCATABLE :: flux_dn_rms(:), flux_dn_bias(:)
  REAL, ALLOCATABLE :: flux_up_rms(:), flux_up_bias(:)
  REAL, ALLOCATABLE :: flux_net_rms(:), flux_net_bias(:)
  REAL, ALLOCATABLE :: band_flux_dn_rms(:), band_flux_dn_bias(:)
  REAL, ALLOCATABLE :: band_flux_up_rms(:), band_flux_up_bias(:)
  REAL, ALLOCATABLE :: band_flux_net_rms(:), band_flux_net_bias(:)


  PRINT *, 'First flux netcdf input file name'
  READ (*,'(A)') fluxfile1
  WRITE (*,*) TRIM(fluxfile1)

  PRINT *, 'Second flux netcdf input file name'
  READ (*,'(A)') fluxfile2
  WRITE (*,*) TRIM(fluxfile2)

  PRINT *, 'Broadband flux profile comparison output file name'
  READ (*,'(A)') profoutfile
  WRITE (*,*) TRIM(profoutfile)

  PRINT *, 'Band flux comparison output file name'
  READ (*,'(A)') bandoutfile
  WRITE (*,*) TRIM(bandoutfile)

  PRINT *, 'Band number for band flux comparison (or 0 for rms over levels for each bands)'
  READ (*,*) compband
  WRITE (*,*) compband


  CALL Read_netcdf_size (fluxfile1, Nlev, Ncol, Nband)
  CALL Read_netcdf_size (fluxfile2, Nlev2, Ncol2, Nband2)
  IF (Nlev2 /= Nlev .OR. Ncol2 /= Ncol .OR. Nband2 /= Nband) THEN
    PRINT *, 'Input netcdf files have different array sizes'
    STOP
  ENDIF

  ALLOCATE (preslev1(Ncol,Nlev), preslev2(Ncol,Nlev))
  ALLOCATE (band_waveno1(Nband,2), band_waveno2(Nband,2))
  ALLOCATE (flux_dn1(Ncol,Nlev), flux_up1(Ncol,Nlev), flux_net1(Ncol,Nlev))
  ALLOCATE (flux_dn2(Ncol,Nlev), flux_up2(Ncol,Nlev), flux_net2(Ncol,Nlev))
  ALLOCATE (band_flux_dn1(Nband,Ncol,Nlev), band_flux_up1(Nband,Ncol,Nlev), &
            band_flux_net1(Nband,Ncol,Nlev))
  ALLOCATE (band_flux_dn2(Nband,Ncol,Nlev), band_flux_up2(Nband,Ncol,Nlev), &
            band_flux_net2(Nband,Ncol,Nlev))

  CALL Read_netcdf_fluxes (fluxfile1, Ncol, Nlev, Nband, preslev1, band_waveno1, &
                           flux_dn1, flux_up1, flux_net1, &
                           band_flux_dn1, band_flux_up1, band_flux_net1)
  CALL Read_netcdf_fluxes (fluxfile2, Ncol, Nlev, Nband, preslev2, band_waveno2, &
                           flux_dn2, flux_up2, flux_net2, &
                           band_flux_dn2, band_flux_up2, band_flux_net2)
  IF (ANY(preslev1(:,:) /= preslev2(:,:))) THEN
    stop 'Pressure levels in two input files disagree'
  ENDIF
  IF (ANY(band_waveno1(:,:) /= band_waveno2(:,:)) .AND. &
       SUM(band_waveno1(:,:))>1.0 .AND. SUM(band_waveno2(:,:))>1.0) THEN
    stop 'Band wavenumber ranges in two input files disagree'
  ENDIF

  ALLOCATE (preslevavg(Nlev))
  ALLOCATE (flux_dn_rms(Nlev), flux_dn_bias(Nlev))
  ALLOCATE (flux_up_rms(Nlev), flux_up_bias(Nlev))
  ALLOCATE (flux_net_rms(Nlev), flux_net_bias(Nlev))
  preslevavg(:) = SUM(0.5*(preslev1(:,:)+preslev2(:,:)),DIM=1)/Ncol

   ! Do the broadband comparisons
  DO ilev = 1, Nlev
    flux_dn_rms(ilev) = SQRT(SUM((flux_dn2(:,ilev)-flux_dn1(:,ilev))**2)/Ncol)
    flux_dn_bias(ilev) = SUM(flux_dn2(:,ilev)-flux_dn1(:,ilev))/Ncol
    flux_up_rms(ilev) = SQRT(SUM((flux_up2(:,ilev)-flux_up1(:,ilev))**2)/Ncol)
    flux_up_bias(ilev) = SUM(flux_up2(:,ilev)-flux_up1(:,ilev))/Ncol
    flux_net_rms(ilev) = SQRT(SUM((flux_net2(:,ilev)-flux_net1(:,ilev))**2)/Ncol)
    flux_net_bias(ilev) = SUM(flux_net2(:,ilev)-flux_net1(:,ilev))/Ncol
  ENDDO

  OPEN (UNIT=1, FILE=profoutfile, STATUS='UNKNOWN')
  WRITE (1,'(A,I2,A)') '! Broadband flux comparison over ',Ncol,' atmospheres:'
  WRITE (1,'(A,A)') '!  Input 1: ', TRIM(fluxfile1)
  WRITE (1,'(A,A)') '!  Input 2: ', TRIM(fluxfile2)
  WRITE (1,'(A)') '!             RMS differences             Biases (2 - 1)'
  WRITE (1,'(A)') '! Pres      Fdn     Fup     Fnet      Fdn      Fup      Fnet'
  DO ilev = 1, Nlev
    WRITE (1,'(1X,F7.2,1X,3(1X,F7.4),1X,3(1X,F8.4))') preslevavg(ilev), &
         flux_dn_rms(ilev), flux_up_rms(ilev), flux_net_rms(ilev), &
         flux_dn_bias(ilev), flux_up_bias(ilev), flux_net_bias(ilev)
  ENDDO
  CLOSE (1)

  IF (compband == 0) THEN
    ALLOCATE (band_flux_dn_rms(Nband), band_flux_dn_bias(Nband))
    ALLOCATE (band_flux_up_rms(Nband), band_flux_up_bias(Nband))
    ALLOCATE (band_flux_net_rms(Nband), band_flux_net_bias(Nband))

     ! Do the band flux comparisons
    DO ibnd = 1, Nband
      band_flux_dn_rms(ibnd) = SQRT(SUM((band_flux_dn2(ibnd,:,:)-band_flux_dn1(ibnd,:,:))**2)/(Ncol*Nlev))
      band_flux_dn_bias(ibnd) = SUM(band_flux_dn2(ibnd,:,:)-band_flux_dn1(ibnd,:,:))/(Ncol*Nlev)
      band_flux_up_rms(ibnd) = SQRT(SUM((band_flux_up2(ibnd,:,:)-band_flux_up1(ibnd,:,:))**2)/(Ncol*Nlev))
      band_flux_up_bias(ibnd) = SUM(band_flux_up2(ibnd,:,:)-band_flux_up1(ibnd,:,:))/(Ncol*Nlev)
      band_flux_net_rms(ibnd) = SQRT(SUM((band_flux_net2(ibnd,:,:)-band_flux_net1(ibnd,:,:))**2)/(Ncol*Nlev))
      band_flux_net_bias(ibnd) = SUM(band_flux_net2(ibnd,:,:)-band_flux_net1(ibnd,:,:))/(Ncol*Nlev)
    ENDDO

    OPEN (UNIT=2, FILE=bandoutfile, STATUS='UNKNOWN')
    WRITE (2,'(A,I2,A,I2,A)') '! Band flux comparison over ',Ncol,' atmospheres and ',Nlev,' levels:'
    WRITE (2,'(A,A)') '!  Input 1: ', TRIM(fluxfile1)
    WRITE (2,'(A,A)') '!  Input 2: ', TRIM(fluxfile2)
    WRITE (2,'(A)') '!                      RMS differences            Biases (2 - 1)'
    WRITE (2,'(A)') '!   Band_waveno     Fdn     Fup     Fnet      Fdn      Fup      Fnet'
    DO ibnd = 1, Nband
      WRITE (2,'(I2,2(1X,F6.0),1X,3(1X,F7.4),1X,3(1X,F8.4))') ibnd, band_waveno1(ibnd,:), &
           band_flux_dn_rms(ibnd), band_flux_up_rms(ibnd), band_flux_net_rms(ibnd), &
           band_flux_dn_bias(ibnd), band_flux_up_bias(ibnd), band_flux_net_bias(ibnd)
    ENDDO
    CLOSE (2)

  ELSE IF (compband > 0 .AND. compband <= Nband) THEN

    ALLOCATE (band_flux_dn_rms(Nlev), band_flux_dn_bias(Nlev))
    ALLOCATE (band_flux_up_rms(Nlev), band_flux_up_bias(Nlev))
    ALLOCATE (band_flux_net_rms(Nlev), band_flux_net_bias(Nlev))

     ! Do the band flux comparisons for one band
    ibnd = compband
    DO ilev = 1, Nlev
      band_flux_dn_rms(ilev) = SQRT(SUM((band_flux_dn2(ibnd,:,ilev)-band_flux_dn1(ibnd,:,ilev))**2)/Ncol)
      band_flux_dn_bias(ilev) = SUM(band_flux_dn2(ibnd,:,ilev)-band_flux_dn1(ibnd,:,ilev))/Ncol
      band_flux_up_rms(ilev) = SQRT(SUM((band_flux_up2(ibnd,:,ilev)-band_flux_up1(ibnd,:,ilev))**2)/Ncol)
      band_flux_up_bias(ilev) = SUM(band_flux_up2(ibnd,:,ilev)-band_flux_up1(ibnd,:,ilev))/Ncol
      band_flux_net_rms(ilev) = SQRT(SUM((band_flux_net2(ibnd,:,ilev)-band_flux_net1(ibnd,:,ilev))**2)/Ncol)
      band_flux_net_bias(ilev) = SUM(band_flux_net2(ibnd,:,ilev)-band_flux_net1(ibnd,:,ilev))/Ncol
    ENDDO

    OPEN (UNIT=2, FILE=bandoutfile, STATUS='UNKNOWN')
    WRITE (2,'(A,I2,A,I2,A)') '! Band flux comparison over ',Ncol,' atmospheres for band ', compband,':'
    WRITE (2,'(A,A)') '!  Input 1: ', TRIM(fluxfile1)
    WRITE (2,'(A,A)') '!  Input 2: ', TRIM(fluxfile2)
    WRITE (2,'(A)') '!                RMS differences                 Biases (2 - 1)'
    WRITE (2,'(A)') '! Pres       Fdn       Fup       Fnet        Fdn        Fup        Fnet'
    DO ilev = 1, Nlev
      WRITE (2,'(1X,F7.2,1X,3(1X,F9.6),1X,3(1X,F10.6))') preslevavg(ilev), &
           band_flux_dn_rms(ilev),  band_flux_up_rms(ilev),  band_flux_net_rms(ilev), &
           band_flux_dn_bias(ilev), band_flux_up_bias(ilev), band_flux_net_bias(ilev)
    ENDDO
    CLOSE (2)
  ENDIF

  DEALLOCATE (band_flux_dn_rms, band_flux_dn_bias)
  DEALLOCATE (band_flux_up_rms, band_flux_up_bias)
  DEALLOCATE (band_flux_net_rms, band_flux_net_bias)
  DEALLOCATE (flux_dn_rms, flux_up_rms, flux_net_rms)
  DEALLOCATE (flux_dn_bias, flux_up_bias, flux_net_bias)
  DEALLOCATE (preslevavg, preslev1, preslev2)
  DEALLOCATE (band_flux_dn2, band_flux_up2, band_flux_net2)
  DEALLOCATE (band_flux_dn1, band_flux_up1, band_flux_net1)
  DEALLOCATE (flux_dn2, flux_up2, flux_net2)
  DEALLOCATE (flux_dn1, flux_up1, flux_net1)
END PROGRAM compare_netcdf_fluxes

MODULE MIE_GAMMA_WATER

CONTAINS

SUBROUTINE LIQCLOUD_MIE_CALCS (WAVELEN1, WAVELEN2, AVGFLAG, DELTAWAVE, &
                               ALPHA, NREFF, REFF, MAXRADIUS, MAXLEG, &
                               EXTINCT, SSALB, NLEG, LEGCOEF)
! Uses Mie theory to calculate scattering properties of specified 
! effective radii for gamma size distributions of spherical water droplets.
! The program provides the index of refraction depending on wavelength.
! The scattering properties may be averaged over the desired spectral 
! range with Planck function weighting or use a weighted central wavelength
! with an averaged index of refraction.  The output volume extinction
! coefficients are for a liquid water content of 1 g/m^3.  The phase 
! functions in the output are represented with Legendre series coefficients.
!
!   Input parameters:
! WAVELEN1, WAVELEN2  Wavelength range (micron) WAVELEN1<=WAVELEN2
! AVGFLAG             Average Mie properties over wavelength or use Planck
!                       weighted refractive index at center wavelength (A or C)
! DELTAWAVE           Wavelength spacing for averaging (micron)
! ALPHA               Gamma size distribution shape parameter
! NREFF               Number of input effective radii
! REFF(:)             List of effective radii (micron)
! MAXRADIUS           Maxium particle radius in size distributions (micron)
! MAXLEG              Maximum order of the Legendre phase function series
!   Output parameters:
! EXTINCT             Extinction coefficient (km^-1) for each effective radius
! SSALB               Single scattering albedo for each r_eff
! NLEG                Number of Legendre coefficients for each r_eff
! LEGCOEF             Legendre coefficients for each r_eff; 
!                       LEGCOEF(0,:)=1.0.  LEGCOEF(1,:)=3*asymmetry_param

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NREFF, MAXLEG
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE
  REAL,    INTENT(IN) :: REFF(NREFF), ALPHA, MAXRADIUS
  CHARACTER(LEN=1), INTENT(IN) :: AVGFLAG
  INTEGER, INTENT(OUT) :: NLEG(NREFF)
  REAL,    INTENT(OUT) :: EXTINCT(NREFF), SSALB(NREFF), LEGCOEF(0:MAXLEG,NREFF)
   ! 
  CHARACTER(LEN=1) :: DISTFLAG
  INTEGER :: NSIZE, MAXLEG2, I, J, L, NL
  REAL    :: PI, WAVELENCEN, XMAX, SCATTER, PARDENS
  COMPLEX :: RINDEX
  INTEGER, ALLOCATABLE :: NLEG1(:)
  REAL, ALLOCATABLE :: RADII(:), ND(:)
  REAL, ALLOCATABLE :: QEXT(:), QSCA(:)
  REAL, ALLOCATABLE :: EXTINCT1(:), SCATTER1(:), LEGEN1(:,:)


   ! Particle size distribution type: G = Gamma, L = Lognormal
  DISTFLAG = 'G'
  PARDENS = 1.0  ! for liquid water
  
   ! Calculate the maximum size parameter and the max number of Legendre terms
  CALL GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
  PI = ACOS(-1.0)
!  IF (AVGFLAG == 'A') THEN
!    XMAX = 2*PI*MAXRADIUS/WAVELEN1
!  ELSE
!    XMAX = 2*PI*MAXRADIUS/WAVELENCEN
!  ENDIF  
!   MAXLEG2 = NINT(2*(XMAX + 4.0*XMAX**0.3334 + 2))
!  IF (MAXLEG < MAXLEG2) THEN
!    PRINT *, 'Input MAXLEG is smaller than maximum possible needed:',MAXLEG,MAXLEG2
!    STOP
!  ENDIF

   ! Get the average index of refraction for water or ice
  CALL GET_REFRACT_INDEX ('W', WAVELEN1, WAVELEN2, RINDEX)

   ! Figure the number of radii there will be
  CALL GET_NSIZE (MINVAL(REFF(:)), MAXRADIUS, WAVELENCEN, NSIZE)

   ! Allocate all the arrays here
  ALLOCATE (RADII(NSIZE), ND(NSIZE), NLEG1(NSIZE))
  ALLOCATE (EXTINCT1(NSIZE), SCATTER1(NSIZE), LEGEN1(0:MAXLEG,NSIZE))

   ! Make up the discrete particle radii to use
  CALL GET_SIZES (MINVAL(REFF(:)), MAXRADIUS, WAVELENCEN, NSIZE, RADII)

   ! Do the Mie computations for each radius, which may involve several
   !   Mie calculation over the wavelength integration
  CALL COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, 'W', &
                              WAVELENCEN, RINDEX, NSIZE, RADII, MAXLEG, &
                              EXTINCT1, SCATTER1, NLEG1, LEGEN1)

  ! Loop over the number of output effective radii
  DO I = 1, NREFF
     ! Calculate the discrete size number concentrations (ND), which vary
     !   according to a truncated gamma or lognormal distribution,
     !   that gives the desired effective radius (REFF) and LWC (1 g/m^3).
    CALL MAKE_SIZE_DIST (DISTFLAG, PARDENS, NSIZE, RADII, REFF(I), ALPHA, ND)

     ! Sum the scattering properties over the discrete size distribution
    EXTINCT(I) = 0.0
    SCATTER = 0.0
    LEGCOEF(:,I) = 0.0
    NL = 1
    DO J = 1, NSIZE
      EXTINCT(I) = EXTINCT(I) + ND(J)*EXTINCT1(J)
      SCATTER = SCATTER + ND(J)*SCATTER1(J)
      NL = MAX(NL,NLEG1(J))
      LEGCOEF(0:NL,I) = LEGCOEF(0:NL,I) + ND(J)*LEGEN1(0:NL,J)
    ENDDO
    DO L = 0, NL
      LEGCOEF(L,I) = LEGCOEF(L,I)/SCATTER
      IF (LEGCOEF(L,I) .GT. 0.5E-5) NLEG(I) = L
    ENDDO
    IF (ABS(LEGCOEF(0,I)-1.0) > 0.0001) THEN
      PRINT *,'Phase function not normalized for Reff=',REFF(I),LEGCOEF(0,I)
      STOP
    ENDIF 
    IF (EXTINCT(I) > 0.0) THEN
      SSALB(I) = SCATTER/EXTINCT(I)
    ENDIF
    EXTINCT(I) = 0.001*EXTINCT(I)

  ENDDO  ! end of effective radius loop

  DEALLOCATE (ND, EXTINCT1, SCATTER1, NLEG1, LEGEN1)
END SUBROUTINE LIQCLOUD_MIE_CALCS 



SUBROUTINE GET_NSIZE (SRETAB, MAXRADIUS, WAVELEN, NSIZE)
 ! Calculates the number of radii for which the Mie computation will be run.
 ! The formula and spacing in size parameter can be changed to trade
 ! off size distribution integration accuracy vs. computer time.
  IMPLICIT NONE
  REAL,    INTENT(IN)  :: SRETAB, MAXRADIUS, WAVELEN
  INTEGER, INTENT(OUT) :: NSIZE
  REAL    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD

  TWOPI = 2.0*ACOS(-1.0)
  RADMIN = 0.02*SRETAB
  RAD = RADMIN
  NSIZE = 1
  DO WHILE (RAD < MAXRADIUS)
    X = TWOPI*RAD/WAVELEN
    DELX = MAX(0.01,0.03*X**0.5)    ! coarser spacing at large size parameters
!    DELX = 0.1                     ! One alternative method
    DELRAD = MIN(MAXRADIUS/100,DELX*WAVELEN/TWOPI)
    RAD = RAD + DELRAD
    NSIZE = NSIZE + 1
  ENDDO
END SUBROUTINE GET_NSIZE


SUBROUTINE GET_SIZES (SRETAB, MAXRADIUS, WAVELEN, NSIZE, RADII)
 ! Calculates the radii for which the Mie computation will be run and
 ! from which all the size distributions will be computed.
 ! The formula and spacing in size parameter can be changed to trade
 ! off size distribution integration accuracy vs. computer time.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE
  REAL,    INTENT(IN) :: SRETAB, MAXRADIUS, WAVELEN
  REAL,    INTENT(OUT) :: RADII(NSIZE)
  INTEGER :: N
  REAL    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD

  TWOPI = 2.0*ACOS(-1.0)
  RAD = 0.02*SRETAB
  RADII(1) = RAD
  DO N = 2, NSIZE
    X = TWOPI*RAD/WAVELEN
    DELX = MAX(0.01,0.03*X**0.5)    ! coarser spacing at large size parameters
!    DELX = 0.1                     ! One alternative method
    DELRAD = MIN(MAXRADIUS/100,DELX*WAVELEN/TWOPI)
    RAD = RAD + DELRAD
    RADII(N) = RAD
  ENDDO
END SUBROUTINE GET_SIZES



SUBROUTINE GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
!  Returns the Planck weighted center wavelength averaged over the 
! wavelength interval (WAVELEN1 < WAVELEN2 [microns]).  A solar
! blackbody temperature of 5800 K is used for the Planck weighting
! if the average wavelength is less than 3 microns, no Planck weighting
! is done for an average wavelength between 3 and 5 microns, and a 
! blackbody temperature of 270 K is done for an average wavelength 
! greater than 5 microns.
  IMPLICIT NONE
  REAL, INTENT(IN)  :: WAVELEN1, WAVELEN2
  REAL, INTENT(OUT) :: WAVELENCEN
  REAL :: WAVECEN, DELWAVE, WAVE, BBTEMP, PLANCK, SUMP, SUMW

  IF (WAVELEN1 == WAVELEN2) THEN
    WAVELENCEN = WAVELEN1  
  ELSE
    WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
    IF (WAVECEN < 3.0) THEN
      BBTEMP = 5800.0
    ELSE IF (WAVECEN > 5.0) THEN
      BBTEMP = 270.0
    ELSE
      BBTEMP = -1.0
      PLANCK = 1.0
    ENDIF 
    DELWAVE = MIN(WAVECEN/100.,0.1*ABS(WAVELEN2-WAVELEN1))
    SUMP = 0.0
    SUMW = 0.0
    WAVE = WAVELEN1   
    DO WHILE (WAVE <= WAVELEN2)
      IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
      SUMP = SUMP + PLANCK
      SUMW = SUMW + PLANCK*WAVE
      WAVE = WAVE + DELWAVE
    ENDDO
    WAVELENCEN = 0.001*NINT(1000*SUMW/SUMP)
  ENDIF
END SUBROUTINE GET_CENTER_WAVELEN




SUBROUTINE GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2, RINDEX)
 ! Returns the index of refraction for water or ice averaged over
 ! the wavelength interval (WAVELEN1 < WAVELEN2 [microns]).   The
 ! averaging is done at 0.05 micron intervals and is weighted by
 ! a Planck function at 5800 K temperature for central wavelengths
 ! less than 3 microns, a flat function between 3 and 5 microns, and
 ! 270 K Planck function beyond 5 microns.  
 ! The index of refraction is at -30 C for ice and +10 C for water
 ! (the temperature dependence is important in the microwave).
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN) :: PARTYPE
  REAL, INTENT(IN) :: WAVELEN1, WAVELEN2
  COMPLEX, INTENT(OUT) :: RINDEX
  REAL :: WAVECEN, WAVECUT, DELWAVE, WAVE, BBTEMP, PLANCK
  REAL :: MRE, MIM, SUMP, SUMMR, SUMMI, A

  WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
  IF (WAVECEN < 3.0) THEN
    BBTEMP = 5800.0
  ELSE IF (WAVECEN > 5.0) THEN
    BBTEMP = 270.0
  ELSE
    BBTEMP = -1.0
    PLANCK = 1.0
  ENDIF 
  DELWAVE = MIN(WAVECEN/100.,0.1*ABS(WAVELEN2-WAVELEN1))
  DELWAVE = MAX(DELWAVE,WAVECEN*1.0E-5)
  SUMP = 0.0
  SUMMR = 0.0
  SUMMI = 0.0
  WAVE = WAVELEN1
  DO WHILE (WAVE <= WAVELEN2)
    IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
    SUMP = SUMP + PLANCK
    IF (PARTYPE == 'I') THEN
!      CALL REFICE (0, WAVE, 243.0, MRE, MIM, A, A)
      stop 'Ice index of refraction not available'
    ELSE
      CALL REFWAT (0, WAVE, 283.0, MRE, MIM, A, A)
    ENDIF
    SUMMR = SUMMR + PLANCK*MRE
    SUMMI = SUMMI + PLANCK*MIM
    WAVE = WAVE + DELWAVE
  ENDDO
  MRE = SUMMR/SUMP
  MIM = SUMMI/SUMP
  RINDEX = CMPLX(MRE,-MIM)
END SUBROUTINE GET_REFRACT_INDEX


 

SUBROUTINE COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, &
                                  PARTYPE, WAVELENCEN, RINDEX, NSIZE, RADII, &
                                  MAXLEG, EXTINCT1, SCATTER1, NLEG1, LEGEN1)
 ! Does a Mie computation for each particle radius in RADII and returns the
 ! optical properties in arrays EXTINCT1, SCATTER1, NLEG1, and LEGEN1.
 ! For AVGFLAG='C' the computation is done at a single wavelength (WAVELENCEN),
 ! using the input index of refraction (RINDEX).  For AVGFLAG='A' an
 ! integration of the Mie properties over wavelength is performed for
 ! each radius.  For each wavelength, with spacing DELTAWAVE, the water 
 ! or ice (depending on PARTYPE) index of refraction is obtained and
 ! used in the Mie computation for that wavelength, and the Mie optical
 ! properties are averaged with Planck function weighting (blackbody
 ! temperature depends on wavelength).  The Legendre coefficients are
 ! returned with the product of the phase function times the scattering
 ! coefficient.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE, MAXLEG
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE, WAVELENCEN
  REAL,    INTENT(IN) :: RADII(NSIZE)
  COMPLEX, INTENT(IN) :: RINDEX
  CHARACTER(LEN=1), INTENT(IN) :: AVGFLAG, PARTYPE
  INTEGER, INTENT(OUT) :: NLEG1(NSIZE)
  REAL,    INTENT(OUT) :: EXTINCT1(NSIZE), SCATTER1(NSIZE)
  REAL,    INTENT(OUT) :: LEGEN1(0:MAXLEG,NSIZE)
  INTEGER :: I, NL
  REAL    :: WAVECEN, WAVE, BBTEMP, PLANCK, SUMP, A
  REAL    :: MRE, MIM, EXT, SCAT, LEG(0:MAXLEG)
  COMPLEX :: REFIND

  IF (AVGFLAG == 'C') THEN
     ! For using one central wavelength: just call Mie routine for each radius
    DO I = 1, NSIZE
      CALL MIE_ONE (WAVELENCEN, RINDEX, RADII(I), MAXLEG, &
                    EXTINCT1(I), SCATTER1(I), NLEG1(I), LEGEN1(:,I) )
    ENDDO

  ELSE
     ! For averaging over wavelength range:
    WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
    IF (WAVECEN < 3.0) THEN
      BBTEMP = 5800.0
    ELSE IF (WAVECEN > 5.0) THEN
      BBTEMP = 270.0
    ELSE
      BBTEMP = -1.0
      PLANCK = 1.0
    ENDIF 
    EXTINCT1(:) = 0.0
    SCATTER1(:) = 0.0
    NLEG1(:) = 1
    LEGEN1(:,:) = 0.0
    SUMP = 0.0
    WAVE = WAVELEN1
    DO WHILE (WAVE <= WAVELEN2)   ! Loop over the wavelengths
      IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
      SUMP = SUMP + PLANCK
      IF (PARTYPE == 'I') THEN   ! Get the index of refraction of water or ice
!        CALL REFICE (0, WAVE, 243.0, MRE, MIM, A, A)
        stop 'Ice index of refraction not available'
      ELSE
        CALL REFWAT (0, WAVE, 283.0, MRE, MIM, A, A)
      ENDIF
      REFIND = CMPLX(MRE,-MIM)
      DO I = 1, NSIZE
        CALL MIE_ONE (WAVE, REFIND, RADII(I), MAXLEG, EXT, SCAT, NL, LEG)
        EXTINCT1(I) = EXTINCT1(I) + PLANCK*EXT
        SCATTER1(I) = SCATTER1(I) + PLANCK*SCAT
        NLEG1(I) = MAX(NLEG1(I),NL)
        LEGEN1(0:NL,I) = LEGEN1(0:NL,I) + PLANCK*LEG(0:NL)
      ENDDO
      WAVE = WAVE + DELTAWAVE
    ENDDO
    EXTINCT1(:) = EXTINCT1(:)/SUMP
    SCATTER1(:) = SCATTER1(:)/SUMP
    LEGEN1(:,:) = LEGEN1(:,:)/SUMP
  ENDIF
END SUBROUTINE COMPUTE_MIE_ALL_SIZES





SUBROUTINE MAKE_SIZE_DIST (DISTFLAG, PARDENS, NSIZE, RADII, REFF, ALPHA, ND)
 ! Calculates the number concentrations (ND in cm^-3) for the NSIZE
 ! discrete particle radii (micron) of a gamma or lognormal size distribution
 ! with an effective radius of REFF (micron), gamma shape parameter or 
 ! lognormal standard deviation of ALPHA, and mass content of 1 g/m^3.
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NSIZE
  REAL,    INTENT(IN)  :: RADII(NSIZE), REFF, ALPHA, PARDENS
  REAL,    INTENT(OUT) :: ND(NSIZE)
  CHARACTER(LEN=1), INTENT(IN) :: DISTFLAG
  REAL, PARAMETER :: TOL=0.001  ! fractional tolerance in achieving Reff
  INTEGER :: I
  REAL    :: TRUERE, F, REHI, RELO, REMID

   ! See if the true effective radius is already close enough
  CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, REFF, NSIZE,RADII, ND, TRUERE)
  IF (ABS(TRUERE-REFF) < TOL*REFF) RETURN
  F = REFF/TRUERE

  IF (TRUERE < REFF) THEN
    ! Find Reff that gives true Reff above desired value
    RELO = REFF
    REHI = F*REFF
    I = 0
    TRUERE = REFF/F
    DO WHILE (TRUERE <= REFF .AND. I < 8)
      REHI = F*REHI
      I = I + 1
      CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, REHI, NSIZE,RADII, ND, TRUERE)
    ENDDO
    IF (TRUERE <= REFF) THEN
      PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,TRUERE
      STOP
    ENDIF
  ELSE
    ! Find Reff that gives true Reff below desired value
    REHI = REFF
    RELO = F*REFF
    I = 0
    TRUERE = REFF/F
    DO WHILE (TRUERE >= REFF .AND. I < 8)
      RELO = F*RELO
      I = I + 1
      CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, RELO, NSIZE,RADII, ND, TRUERE)
    ENDDO
    IF (TRUERE >= REFF) THEN
      PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,TRUERE
      STOP
    ENDIF
  ENDIF
  ! Do bisection to get correct effective radius
  DO WHILE (ABS(TRUERE-REFF) > TOL*REFF)
    REMID = 0.5*(RELO+REHI)
    CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, REMID, NSIZE,RADII, ND, TRUERE)
    IF (TRUERE < REFF) THEN
      RELO = REMID
    ELSE
      REHI = REMID
    ENDIF
  ENDDO  
END SUBROUTINE MAKE_SIZE_DIST



SUBROUTINE DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, RE, NSIZE, RADII, &
                         ND, TRUERE)
 ! For the input effective radius (RE) [um], returns the number concentrations 
 ! ND [cm^-3] and the calculated effective radius TRUERE [um] for a 
 ! gamma or lognormal size distribution with mass content of 1 g/m^3.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE
  REAL,    INTENT(IN) :: PARDENS, ALPHA, RE, RADII(NSIZE)
  CHARACTER(LEN=1), INTENT(IN) :: DISTFLAG
  REAL,    INTENT(OUT) :: ND(NSIZE), TRUERE
  INTEGER :: J
  REAL    :: DENS, PI, A, B, LWC, R, DELR, SUM2, SUM3

  PI = ACOS(-1.0)
  IF (DISTFLAG == 'G') THEN
    B = (ALPHA+3)/RE
    A = 1.E6/( (4*PI/3.)*PARDENS *B**(-ALPHA-4) *EXP(GAMMLN(ALPHA+4.)) )
  ELSE
    B = RE*EXP(-2.5*ALPHA**2)
    A = 1.E6/( (4*PI/3.)*PARDENS *SQRT(2*PI)*ALPHA * B**3 *EXP(4.5*ALPHA**2) )
  ENDIF
  LWC = 0.0
  SUM2 = 0.0
  SUM3 = 0.0
  DO J = 1, NSIZE
    R = RADII(J)
    DELR = SQRT(RADII(J)*RADII(MIN(NSIZE,J+1))) &
         - SQRT(RADII(J)*RADII(MAX(1,J-1)))
    IF (DISTFLAG == 'G') THEN
      ND(J) = A* R**ALPHA *EXP(-B*R) *DELR
    ELSE
      ND(J) = (A/R)*EXP(-0.5*(LOG(R/B))**2/ALPHA**2) *DELR
    ENDIF
    LWC = LWC + 1.0E-6*PARDENS*ND(J)*(4*PI/3)*R**3
    SUM2 = SUM2 + ND(J)*R**2
    SUM3 = SUM3 + ND(J)*R**3
  ENDDO
  ND(:) = (1.0/LWC)*ND(:)
  TRUERE = SUM3/SUM2
END SUBROUTINE DO_SIZE_DIST



SUBROUTINE MIE_ONE (WAVELENGTH, MINDEX, RADIUS, MAXLEG, &
                    EXTINCTION, SCATTER, NLEG, LEGEN)
 ! Computes the Mie scattering properties for a single homogeneous
 ! sphere of radius RADIUS.  The phase function times the scattering
 ! coefficient is returned as Legendre series coefficients.
  IMPLICIT NONE
  REAL,    INTENT(IN) :: WAVELENGTH, RADIUS
  COMPLEX, INTENT(IN) :: MINDEX
  INTEGER, INTENT(IN) :: MAXLEG
  REAL,    INTENT(OUT) :: EXTINCTION, SCATTER, LEGEN(0:)
  INTEGER, INTENT(OUT) :: NLEG
  INTEGER :: NMIE, NQUAD
  INTEGER :: I, L
  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  REAL(DP) :: X, PI
  REAL(DP) :: QEXT, QSCAT, GEOMAREA
  REAL(DP), ALLOCATABLE :: MU(:), WTS(:)
  REAL(DP), ALLOCATABLE :: COEF1(:)
  REAL(DP) :: P1, PL, PL1, PL2
  COMPLEX(DP) :: MSPHERE
  COMPLEX(DP), ALLOCATABLE :: A(:), B(:)
      
  PI = ACOS(-1.0D0)      
  X = 2.0D0*PI*RADIUS/WAVELENGTH
  GEOMAREA = PI*RADIUS**2
  MSPHERE = MINDEX

   ! Calculate the number of terms in the Mie series (Nstop)
  NMIE = X + 4.0*X**0.3334 + 2
  ALLOCATE (A(NMIE), B(NMIE))

   ! Compute the An's and Bn's and the cross sections
  CALL MIECALC (NMIE, X, MSPHERE, A, B)
  CALL MIECROSS (NMIE, X, A, B, QEXT, QSCAT)
  EXTINCTION = GEOMAREA*QEXT
  SCATTER = GEOMAREA*QSCAT

   !  Get the Gauss-Legendre quadrature abscissas and weights
  NLEG = MIN(MAXLEG,2*NMIE)
  NQUAD  = (NLEG + 2*NMIE + 2)/2
  ALLOCATE (MU(NQUAD), WTS(NQUAD), COEF1(0:NLEG))
  CALL GAUSQUAD (NQUAD, MU, WTS)

   ! Calculate the phase function for quadrature angles and then
   ! integrate over angle to get the Legendre coefficients. 
  COEF1(:) = 0.0
  DO I = 1, NQUAD
    CALL MIEANGLE (NMIE, A, B, MU(I), P1)
    PL1 = 1.0
    PL = 1.0
    DO L = 0, NLEG
      IF (L > 0)  PL = (2*L-1)*MU(I)*PL1/L - (L-1)*PL2/L
      COEF1(L) = COEF1(L) + PL*P1*WTS(I)
      PL2 = PL1
      PL1 = PL
    ENDDO
  ENDDO
  DO L = 0, NLEG
    LEGEN(L) = (2*L+1)/2.0 *(WAVELENGTH**2/PI) *COEF1(L)
  ENDDO
       
END SUBROUTINE MIE_ONE



SUBROUTINE MIECALC (NTERMS, X, MN, A, B)
 ! MIECALC calculates the complex Mie coefficients An and Bn
 ! given the dimensionless size parameter X and the complex
 ! index of refraction (Mre,Mim).  The number of terms calculated
 ! is input in NTERMS.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NTERMS
  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  REAL(DP), INTENT(IN) :: X
  COMPLEX(DP), INTENT(IN) :: MN
  COMPLEX(DP), INTENT(OUT) :: A(NTERMS), B(NTERMS)
  INTEGER :: NSTOP, N, NN
  REAL(DP) :: PSIN, PSIM, CHIN, CHIM, TMP
  REAL(DP) :: DCOS, DSIN
  COMPLEX(DP) :: M, Y, D(NTERMS+15), XIN, XIM, CTMP
 
   ! Generate the Dn's by down recurrence  D = d(log(PSI(y)))/dy
  M = CONJG(MN)
  Y = M*X
  NN = NTERMS + 15
  D(NN) = CMPLX (0.0D0, 0.0D0)
  DO N = NN, 2, -1
      D(N-1) = N/Y - 1.0/ (D(N) + N/Y)
  ENDDO
 
   ! Generate the PSIn's and XIn'S by upward recurrence
   ! and calculate the An's and Bn's from them.
   !  (PSIN = PSI(n), PSIM = PSI(n-1), same for CHI)
  PSIM = DCOS(X)
  PSIN = DSIN(X)
  CHIM = -DSIN(X)
  CHIN = DCOS(X)
  DO N = 1, NTERMS
    TMP = PSIN
    PSIN = (2*N-1)/X *PSIN - PSIM
    PSIM = TMP
    TMP = CHIN
    CHIN = (2*N-1)/X *CHIN - CHIM
    CHIM = TMP
    XIN = CMPLX (PSIN, -CHIN)
    XIM = CMPLX (PSIM, -CHIM)
    CTMP = D(N)/M + N/X
    A(N) = (CTMP*PSIN - PSIM) / (CTMP*XIN - XIM)
    CTMP = M*D(N) + N/X
    B(N) = (CTMP*PSIN - PSIM) / (CTMP*XIN - XIM)
  ENDDO
END SUBROUTINE MIECALC


SUBROUTINE MIECROSS (NTERMS, X, A, B, QEXT, QSCAT)
 ! MIECROSS calculates the extinction and scattering efficiencies 
 ! given the Mie coefficients An and Bn and the size parameter X.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NTERMS
  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  REAL(DP), INTENT(IN) :: X
  COMPLEX(DP), INTENT(IN) :: A(NTERMS), B(NTERMS)
  REAL(DP), INTENT(OUT) :: QEXT, QSCAT
  INTEGER :: N
  REAL(DP) :: SUM1, SUM2, DREAL
 
  SUM1 = 0.0D0
  SUM2 = 0.0D0
  DO N = 1, NTERMS
    SUM1 = SUM1 + (2*N+1)*( REAL(A(N)) + REAL(B(N)) )
    SUM2 = SUM2 + (2*N+1)*( REAL(A(N)*CONJG(A(N))) &
                          + REAL(B(N)*CONJG(B(N))) )
  ENDDO
  QEXT = 2.0D0/X**2 * SUM1
  QSCAT = 2.0D0/X**2 * SUM2
ENDSUBROUTINE MIECROSS


SUBROUTINE MIEANGLE (NTERMS, A, B, MU, P1)
 ! MIEANGLE calculates the intensity scattering matrix element
 ! (P1) for a particular value of MU (cos(theta)) from the
 ! Mie coefficients An's and Bn's.  The matrix element is for the
 ! Stokes intensity vector (I) and is calculated from the
 ! complex scattering amplitudes S1 and S2.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NTERMS
  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  REAL(DP), INTENT(IN) :: MU
  COMPLEX(DP), INTENT(IN) :: A(NTERMS), B(NTERMS)
  REAL(DP), INTENT(OUT) :: P1
  INTEGER :: N
  REAL(DP) :: TMP, PIN, PIM, TAUN, C
  COMPLEX(DP) :: S1, S2
 
  S1 = CMPLX(0.0D0,0.0D0)
  S2 = CMPLX(0.0D0,0.0D0)
    ! Sum up the series using the An's and Bn's
  PIN = 1.0
  PIM = 0.0
  DO N = 1, NTERMS
    TAUN = N*MU*PIN - (N+1)*PIM
     ! Calculate the scattering functions at +mu and -mu
     ! using the PIn's and the TAUn's.
    C = (2*N+1.0D0) / (N*(N+1.0D0))
    S1 = S1 + C*( A(N)*PIN + B(N)*TAUN)
    S2 = S2 + C*( B(N)*PIN + A(N)*TAUN)
      ! Calculate the angular function PIn by up recurrence
    TMP = PIN
    PIN = ( (2*N+1)*MU*PIN - (N+1)*PIM ) / N
    PIM = TMP
  ENDDO
   ! Calculate the first Stokes parameter scattering matrix element
  P1 = 0.5*( ABS(S2)**2 + ABS(S1)**2 )
END SUBROUTINE MIEANGLE

 

SUBROUTINE GAUSQUAD (N, XA, WT)
 ! Generates the abscissas (X) and weights (W) for an N point
 ! Gauss-Legendre quadrature.  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  REAL(DP), INTENT(OUT) :: XA(N), WT(N)
  INTEGER :: K, I, J, L
  REAL(DP) :: X, XP, PL, PL1, PL2, DPL, TINY=3.0D-13

  K = (N+1)/2
  DO J = 1, K
    X = COS(3.141592654*(J-.25)/(N+.5))
    XP = X + 1.0
    I = 0
    DO WHILE (ABS(X-XP) > TINY .AND. I < 10)
      PL1 = 1
      PL = X
      DO L = 2, N
        PL2 = PL1
        PL1 = PL
        PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
      ENDDO
      DPL = N*(X*PL-PL1)/(X*X-1)
      XP = X
      X = XP - PL/DPL
      I = I + 1
    ENDDO
    XA(J)     = -X
    XA(N-J+1) = X
    WT(J  )   = 2.0D0/((1.0D0-X*X)*DPL*DPL)
    WT(N-J+1) = WT(J)
  ENDDO
END SUBROUTINE GAUSQUAD


REAL FUNCTION GAMMLN(XIN)
  ! Returns the natural log of the Gamma function of XIN.
  IMPLICIT NONE
  REAL, INTENT(IN) :: XIN
  INTEGER :: J
  INTEGER, PARAMETER :: DP=KIND(1.0D0)
  REAL(DP) :: COF(6) = (/ 76.18009173D0,-86.50532033D0,24.01409822D0, &
                                -1.231739516D0,.120858003D-2,-.536382D-5 /)
  REAL(DP) :: STP=2.50662827465D0, X, TMP, SER

  X=XIN-1.0D0
  TMP=X+5.5D0
  TMP=(X+0.5D0)*LOG(TMP)-TMP
  SER=1.0D0
  DO J = 1, 6
    X=X+1.0D0
    SER=SER+COF(J)/X
  ENDDO
  GAMMLN=TMP+LOG(STP*SER)
END FUNCTION GAMMLN



SUBROUTINE REFWAT(IUNIT,XLAM,T,RN,CN,ABSIND,ABSCOF)
  !  DEFINES WAVELENGTH DEPENDENT COMPLEX INDEX OF REFRACTION FOR WATER
  !  ALLOWABLE WAVELENGTH RANGE EXTENDS FROM .2 MICRONS TO 10 CM
  !  TEMPERATURE DEPENDENCE ONLY CONSIDERED BEYOND 0.1 CM
  !  ERIC A. SMITH
  !  DEPT OF ATMOSPHERIC SCIENCE
  !  COLORADO STATE UNIVERSITY
  !  FORT COLLINS,CO  80523
  !  TEL   303-491-8533
  !
  !  REFERENCES

  !  0.2 UM - 0.69 UM

  !  HALE,G., AND M. QUERRY,1972.
  !  OPTICAL CONSTANTS OF WATER IN THE 200 NM TO 200 UM WAVELENGTH REGI
  !  APPLIED OPTICS,12,3,555-563.

  !  0.69 UM - 2.0 UM

  !  PALMER,K.F., AND D. WILLIAMS,1974.
  !  OPTICAL PROPERTIES OF WATER IN THE NEAR INFRARED.
  !  JOURNAL OF THE OPTICAL SOCIETY OF AMERICA,64,8,1107-1110.

  !  2.0 UM - 1000.0 UM

  !  DOWNING,H.D., AND D. WILLIAMS,1975.
  !  OPTICAL CONSTANTS OF WATER IN THE INFRARED.
  !  JOURNAL OF GEOPHYSICAL REVIEW,80,12,1656-1661.

  !  1.0 MM - 10.0 CM

  !  RAY,P.S.,1972.
  !  BROADBAND COMPLEX REFRACTIVE INDICES OF ICE AND WATER.
  !  APPLIED OPTICS,11,8,1836-1844.

  !  INPUT PARAMETERS

  !  IUNIT = 0 FOR WAVELENGTH SPECIFIED IN MICRONS
  !        = 1 FOR WAVELENGTH SPECIFIED IN MILLIMETERS
  !        = 2 FOR WAVELENGTH SPECIFIED IN CENTIMETERS
  !        = 3 FOR WAVELENGTH SPECIFIED IN INVERSE CENTIMETERS ( WAVE N
  !  XLAM = WAVELENGTH ( MICRONS OR MM OR CM OR CM**-1 )
  !  T = TEMPERATURE ( DEGREES KELVIN )

  !  OUTPUT PARAMETERS

  !  RN = REAL PORTION ( SCATTERING )
  !  CN = COMPLEX PORTION ( ABSORPTION )
  !  ABSIND = ABSORPTIVE INDEX ( CN/RN )
  !  ABSCOF = ABORPTION COEFFICIENT ( 4*PI*CN/XLAM )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: IUNIT
  REAL,    INTENT(IN) :: XLAM, T
  REAL,    INTENT(OUT) :: RN,CN,ABSIND,ABSCOF

  INTEGER :: NUMWAT, I, I1, I2
  REAL    :: WL, WLCEN, BET, DEL, GAM, ABSFUNC
  REAL    :: WLMIN, WLMAX, CUTWAT, PI
  REAL    :: FAC, TC, T1, T2, XL, SIGMA, ALPHA
  REAL    :: ES, E00, XLAMS, TERM, SINT, COST, XLRAT, POWTRM, DENOM, ER, EI
  COMPLEX E,M
  REAL :: WLTABW(518),RNTABW(518),CNTABW(518)

  DATA NUMWAT/518/
  DATA WLMIN,WLMAX/0.2,100000.0/
  DATA CUTWAT/1000.0/
  DATA WLTABW(1:66) / &
         .20000,    .22500,    .25000,    .27500,    .30000,    .32500, &
         .35001,    .37500,    .40000,    .42501,    .45000,    .47499, &
         .50000,    .52499,    .54999,    .57501,    .59999,    .62500, &
         .64998,    .67499,    .68966,    .70175,    .71429,    .72464, &
         .73529,    .74627,    .75188,    .75758,    .76923,    .78125, &
         .79365,    .80645,    .81301,    .81967,    .83333,    .84746, &
         .86207,    .87719,    .89286,    .90909,    .92593,    .93458, &
         .94340,    .95238,    .96154,    .97276,    .98039,    .99010, &
        1.00000,   1.01010,   1.02041,   1.03093,   1.04167,   1.05263, &
        1.06952,   1.08696,   1.09890,   1.11111,   1.12360,   1.13636, &
        1.14943,   1.16279,   1.17647,   1.19048,   1.20482,   1.21951/
      DATA WLTABW(67:132) / &
        1.23457,   1.25000,   1.26582,   1.28205,   1.29870,   1.31579, &
        1.33333,   1.35135,   1.36986,   1.38889,   1.40845,   1.42857, &
        1.44300,   1.47059,   1.49254,   1.51515,   1.53846,   1.56250, &
        1.58730,   1.61290,   1.63934,   1.66667,   1.69492,   1.72414, &
        1.75439,   1.78571,   1.80180,   1.81818,   1.85185,   1.88679, &
        1.92678,   1.96078,   2.00000,   2.02020,   2.04082,   2.06186, &
        2.08333,   2.10526,   2.12766,   2.15054,   2.17391,   2.19780, &
        2.22222,   2.24719,   2.27273,   2.29885,   2.32558,   2.35294, &
        2.38095,   2.40964,   2.43902,   2.46914,   2.50000,   2.50627, &
        2.51256,   2.51889,   2.52525,   2.53165,   2.53807,   2.54453, &
        2.55102,   2.55754,   2.56410,   2.57069,   2.57732,   2.58398/
      DATA WLTABW(133:198) / &
        2.59067,   2.59740,   2.60417,   2.61097,   2.61780,   2.62467, &
        2.63158,   2.63852,   2.64550,   2.65252,   2.65957,   2.66667, &
        2.67380,   2.68097,   2.68817,   2.69542,   2.70270,   2.71003, &
        2.71739,   2.72480,   2.73224,   2.73973,   2.74725,   2.75482, &
        2.76243,   2.77008,   2.77778,   2.78552,   2.79330,   2.80112, &
        2.80899,   2.81690,   2.82486,   2.83286,   2.84091,   2.84900, &
        2.85714,   2.86533,   2.87356,   2.88184,   2.89017,   2.89855, &
        2.90698,   2.91545,   2.92398,   2.93255,   2.94118,   2.94985, &
        2.95858,   2.96736,   2.97619,   2.98507,   2.99401,   3.00300, &
        3.01205,   3.02115,   3.03030,   3.03951,   3.04878,   3.05810, &
        3.06748,   3.07692,   3.08642,   3.09598,   3.10559,   3.11526/
      DATA WLTABW(199:264) / &
        3.12500,   3.13480,   3.14465,   3.15457,   3.16456,   3.17460, &
        3.18471,   3.19489,   3.20513,   3.21543,   3.22581,   3.23625, &
        3.24675,   3.25733,   3.26797,   3.27869,   3.28947,   3.30033, &
        3.31126,   3.32226,   3.33333,   3.34448,   3.35570,   3.36700, &
        3.37838,   3.38983,   3.40136,   3.41297,   3.42466,   3.43643, &
        3.44828,   3.46021,   3.47222,   3.48432,   3.49650,   3.50877, &
        3.52113,   3.53357,   3.54610,   3.55872,   3.57143,   3.58423, &
        3.59712,   3.61011,   3.62319,   3.63636,   3.64964,   3.66300, &
        3.67647,   3.69004,   3.70370,   3.71747,   3.73134,   3.74532, &
        3.75940,   3.77358,   3.78788,   3.80228,   3.81679,   3.83142, &
        3.84615,   3.86100,   3.87597,   3.89105,   3.90625,   3.92157/
      DATA WLTABW(265:330) / &
        3.93701,   3.95257,   3.96825,   3.98406,   4.00000,   4.01606, &
        4.03226,   4.04858,   4.06504,   4.08163,   4.09836,   4.11523, &
        4.13223,   4.14938,   4.16667,   4.18410,   4.20168,   4.21941, &
        4.23729,   4.25532,   4.27350,   4.29185,   4.31034,   4.32900, &
        4.34783,   4.36681,   4.38596,   4.40529,   4.42478,   4.44444, &
        4.46429,   4.48430,   4.50450,   4.52489,   4.54545,   4.56621, &
        4.58716,   4.60829,   4.62963,   4.65116,   4.67290,   4.69484, &
        4.71698,   4.73934,   4.76190,   4.78469,   4.80769,   4.83092, &
        4.85437,   4.87805,   4.90196,   4.92611,   4.95050,   4.97512, &
        5.00000,   5.02513,   5.05051,   5.07614,   5.10204,   5.12821, &
        5.15464,   5.18135,   5.20833,   5.23560,   5.26316,   5.29101/ 
      DATA WLTABW(331:396) / &
        5.31915,   5.34759,   5.37634,   5.40541,   5.43478,   5.46448, &
        5.49451,   5.52486,   5.55556,   5.58659,   5.61798,   5.64972, &
        5.68182,   5.71429,   5.74713,   5.78035,   5.81395,   5.84795, &
        5.88235,   5.91716,   5.95238,   5.98802,   6.02410,   6.06061, &
        6.09756,   6.13497,   6.17284,   6.21118,   6.25000,   6.28931, &
        6.32911,   6.36943,   6.41026,   6.45161,   6.49351,   6.53595, &
        6.57895,   6.62252,   6.66667,   6.71141,   6.75676,   6.80272, &
        6.84932,   6.89655,   6.94444,   6.99301,   7.04225,   7.09220, &
        7.14286,   7.19424,   7.24638,   7.29927,   7.35294,   7.40741, &
        7.46269,   7.51880,   7.57576,   7.63359,   7.69231,   7.75194, &
        7.81250,   7.87402,   7.93651,   8.00000,   8.06452,   8.13008/
      DATA WLTABW(397:462) / &
        8.19672,   8.26446,   8.33333,   8.40336,   8.47458,   8.54701, &
        8.62069,   8.69565,   8.77193,   8.84956,   8.92857,   9.00901, &
        9.09091,   9.17431,   9.25926,   9.34579,   9.43396,   9.52381, &
        9.61538,   9.70874,   9.80392,   9.90099,  10.00000,  10.10101, &
       10.20408,  10.30928,  10.41667,  10.52632,  10.63830,  10.75269, &
       10.86957,  10.98901,  11.11111,  11.23596,  11.36364,  11.49425, &
       11.62791,  11.76471,  11.90476,  12.04819,  12.19512,  12.34568, &
       12.50000,  12.65823,  12.82051,  12.98701,  13.15789,  13.33333, &
       13.51351,  13.69863,  13.88889,  14.08451,  14.28571,  14.49275, &
       14.70588,  14.92537,  15.15152,  15.38462,  15.62500,  15.87302, &
       16.12903,  16.39344,  16.66667,  16.94915,  17.24138,  17.54386/ 
      DATA WLTABW(463:518) / &
       17.85714,  18.18182,  18.51852,  18.86792,  19.23077,  19.60784, &
       20.00000,  20.40816,  20.83333,  21.27660,  21.73913,  22.22222, &
       22.72727,  23.25581,  23.80952,  24.39024,  25.00000,  25.64103, &
       26.31579,  27.02703,  27.77778,  28.57143,  29.41176,  30.30303, &
       31.25000,  32.25806,  33.33333,  34.48276,  35.71429,  37.03704, &
       38.46154,  40.00000,  41.66667,  43.47826,  45.45455,  47.61905, &
       50.00000,  52.63158,  55.55556,  58.82353,  62.50000,  66.66667, &
       71.42857,  76.92308,  83.33333,  90.90909, 100.00000, 111.11111, &
      125.00000, 142.85714, 166.66667, 200.00000, 250.00000, 333.33333, &
      500.00000,1000.00000/ 
      DATA RNTABW(1:66) / &
     1.396,1.373,1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337, &
     1.336,1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.332,1.332, &
     1.332,1.332,1.332,1.332,1.332,1.332,1.331,1.331,1.331,1.331,1.331, &
     1.330,1.330,1.330,1.330,1.330,1.329,1.329,1.329,1.329,1.329,1.328, &
     1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328, &
     1.327,1.327,1.327,1.327,1.326,1.326,1.326,1.326,1.325,1.325,1.325/
      DATA RNTABW(67:132) / &
     1.325,1.325,1.324,1.324,1.324,1.324,1.323,1.323,1.323,1.322,1.322, &
     1.321,1.321,1.321,1.320,1.320,1.319,1.319,1.318,1.318,1.317,1.316, &
     1.315,1.314,1.314,1.313,1.312,1.312,1.311,1.310,1.309,1.307,1.306, &
     1.301,1.301,1.300,1.298,1.298,1.296,1.295,1.294,1.293,1.291,1.289, &
     1.287,1.285,1.282,1.280,1.277,1.274,1.270,1.265,1.261,1.260,1.259, &
     1.257,1.256,1.255,1.254,1.252,1.250,1.249,1.247,1.246,1.243,1.241/ 
      DATA RNTABW(133:198) / &
     1.240,1.238,1.235,1.232,1.230,1.227,1.224,1.221,1.218,1.214,1.210, &
     1.205,1.200,1.195,1.191,1.185,1.179,1.172,1.166,1.157,1.149,1.144, &
     1.139,1.138,1.138,1.139,1.141,1.144,1.149,1.154,1.158,1.161,1.165, &
     1.171,1.177,1.183,1.191,1.199,1.212,1.220,1.233,1.246,1.258,1.271, &
     1.282,1.293,1.305,1.317,1.329,1.342,1.353,1.364,1.376,1.386,1.398, &
     1.407,1.417,1.426,1.434,1.442,1.450,1.457,1.465,1.471,1.476,1.480/ 
      DATA RNTABW(199:264) / &
     1.483,1.486,1.487,1.487,1.487,1.486,1.485,1.482,1.479,1.477,1.474, &
     1.472,1.467,1.464,1.461,1.457,1.454,1.451,1.448,1.444,1.441,1.437, &
     1.434,1.431,1.427,1.425,1.421,1.418,1.415,1.413,1.410,1.407,1.405, &
     1.403,1.400,1.398,1.396,1.394,1.392,1.390,1.388,1.387,1.385,1.383, &
     1.382,1.379,1.378,1.377,1.375,1.374,1.372,1.371,1.370,1.369,1.367, &
     1.366,1.365,1.363,1.361,1.361,1.360,1.358,1.358,1.357,1.355,1.354/ 
      DATA RNTABW(265:330) / &
     1.353,1.352,1.351,1.350,1.349,1.348,1.348,1.347,1.346,1.345,1.344, &
     1.344,1.343,1.342,1.341,1.340,1.340,1.338,1.337,1.337,1.335,1.334, &
     1.334,1.333,1.332,1.332,1.331,1.330,1.330,1.330,1.329,1.329,1.329, &
     1.328,1.328,1.327,1.327,1.327,1.327,1.327,1.326,1.326,1.326,1.325, &
     1.325,1.325,1.325,1.325,1.325,1.324,1.324,1.323,1.322,1.322,1.321, &
     1.320,1.319,1.318,1.318,1.317,1.316,1.314,1.313,1.311,1.310,1.308/ 
      DATA RNTABW(331:396) / &
     1.306,1.304,1.302,1.299,1.297,1.294,1.291,1.288,1.285,1.282,1.278, &
     1.275,1.271,1.267,1.262,1.256,1.251,1.247,1.242,1.241,1.241,1.247, &
     1.265,1.289,1.311,1.332,1.349,1.354,1.356,1.354,1.350,1.345,1.341, &
     1.337,1.333,1.330,1.326,1.324,1.322,1.320,1.319,1.318,1.317,1.316, &
     1.315,1.314,1.313,1.311,1.310,1.309,1.308,1.307,1.306,1.305,1.303, &
     1.302,1.301,1.300,1.298,1.296,1.295,1.294,1.293,1.291,1.288,1.286/ 
      DATA RNTABW(397:462) / &
     1.285,1.283,1.281,1.279,1.276,1.274,1.271,1.269,1.267,1.264,1.261, &
     1.259,1.256,1.253,1.249,1.246,1.242,1.238,1.234,1.230,1.224,1.220, &
     1.214,1.208,1.202,1.194,1.189,1.181,1.174,1.168,1.162,1.156,1.149, &
     1.143,1.139,1.135,1.132,1.132,1.131,1.132,1.130,1.130,1.134,1.138, &
     1.142,1.157,1.171,1.182,1.189,1.201,1.213,1.223,1.236,1.249,1.264, &
     1.277,1.289,1.303,1.313,1.324,1.335,1.348,1.361,1.372,1.385,1.396/ 
      DATA RNTABW(463:518) / &
     1.407,1.419,1.431,1.441,1.451,1.462,1.470,1.480,1.488,1.496,1.504, &
     1.510,1.515,1.521,1.527,1.532,1.537,1.541,1.545,1.549,1.552,1.552, &
     1.552,1.550,1.546,1.543,1.541,1.539,1.537,1.534,1.532,1.529,1.525, &
     1.528,1.542,1.567,1.600,1.640,1.689,1.746,1.801,1.848,1.890,1.929, &
     1.960,1.982,1.997,2.000,2.010,2.020,2.040,2.070,2.110,2.150,2.225, &
     2.481/ 
      DATA CNTABW(1:66) / &
     1.1000E-07,4.9000E-08,3.4000E-08,2.4000E-08,1.6000E-08,1.1000E-08, &
     6.5000E-09,3.5000E-09,1.9000E-09,1.3000E-09,1.0000E-09,9.4000E-10, &
     1.0000E-09,1.3000E-09,2.0000E-09,3.6000E-09,1.1000E-08,1.4000E-08, &
     1.6000E-08,2.2000E-08,2.7000E-08,3.8000E-08,5.6000E-08,7.7300E-08, &
     1.3900E-07,1.6300E-07,1.6800E-07,1.6400E-07,1.5400E-07,1.4300E-07, &
     1.3300E-07,1.2500E-07,1.2400E-07,1.3000E-07,2.0400E-07,2.6100E-07, &
     2.9400E-07,3.5300E-07,4.3300E-07,5.4300E-07,8.7700E-07,1.1800E-06, &
     1.6100E-06,2.4400E-06,3.6000E-06,3.9800E-06,3.9200E-06,3.7000E-06, &
     3.3100E-06,2.8200E-06,2.3100E-06,1.9000E-06,1.5700E-06,1.3700E-06, &
     1.2600E-06,1.4400E-06,1.6800E-06,2.0500E-06,2.8900E-06,4.9600E-06, &
     8.8700E-06,1.0900E-05,1.1500E-05,1.1800E-05,1.2000E-05,1.1800E-05/ 
      DATA CNTABW(67:132) / &
     1.1500E-05,1.1000E-05,1.0800E-05,1.1500E-05,1.3800E-05,1.7500E-05, &
     2.3900E-05,4.1600E-05,5.9400E-05,1.0100E-04,2.4100E-04,3.5200E-04, &
     3.6400E-04,3.3400E-04,2.5800E-04,1.8800E-04,1.4800E-04,1.2000E-04, &
     1.0200E-04,8.7300E-05,7.9200E-05,7.4900E-05,7.6200E-05,8.5500E-05, &
     1.0600E-04,1.3000E-04,1.3600E-04,1.3700E-04,1.5900E-04,8.6300E-04, &
     1.9000E-03,1.7000E-03,1.1000E-03,9.0000E-04,7.3100E-04,6.1700E-04, &
     5.1400E-04,4.5200E-04,4.0000E-04,3.5900E-04,3.4100E-04,3.3800E-04, &
     3.4500E-04,3.7600E-04,4.1600E-04,4.6500E-04,5.4200E-04,6.5200E-04, &
     7.9200E-04,9.6800E-04,1.2300E-03,1.5600E-03,1.9000E-03,1.9500E-03, &
     2.0000E-03,2.0500E-03,2.0700E-03,2.1000E-03,2.1200E-03,2.1500E-03, &
     2.1900E-03,2.2400E-03,2.2700E-03,2.3100E-03,2.3400E-03,2.3900E-03/ 
      DATA CNTABW(133:198) / &
     2.4300E-03,2.4800E-03,2.5700E-03,2.7000E-03,2.9800E-03,3.3000E-03, &
     4.0200E-03,4.3700E-03,4.8200E-03,5.3600E-03,6.2700E-03,7.3200E-03, &
     8.5500E-03,1.0500E-02,1.2700E-02,1.4500E-02,1.6400E-02,1.8600E-02, &
     2.0500E-02,2.8200E-02,3.8000E-02,4.6200E-02,5.4800E-02,6.4900E-02, &
     7.4400E-02,8.3600E-02,9.2700E-02,1.0200E-01,1.1200E-01,1.2100E-01, &
     1.3100E-01,1.4200E-01,1.5400E-01,1.6700E-01,1.8000E-01,1.9400E-01, &
     2.0600E-01,2.1800E-01,2.2900E-01,2.3900E-01,2.4900E-01,2.5800E-01, &
     2.6500E-01,2.7100E-01,2.7600E-01,2.8000E-01,2.8100E-01,2.8200E-01, &
     2.8200E-01,2.7900E-01,2.7600E-01,2.7200E-01,2.6700E-01,2.6200E-01, &
     2.5500E-01,2.5000E-01,2.4300E-01,2.3600E-01,2.2800E-01,2.2000E-01, &
     2.1200E-01,2.0400E-01,1.9500E-01,1.8300E-01,1.7300E-01,1.6300E-01/ 
      DATA CNTABW(199:264) / &
     1.5300E-01,1.4400E-01,1.3400E-01,1.2500E-01,1.1700E-01,1.1000E-01, &
     9.9400E-02,9.2000E-02,8.5500E-02,7.8500E-02,7.1600E-02,6.5300E-02, &
     6.0000E-02,5.5000E-02,5.0400E-02,4.6200E-02,4.2200E-02,3.8500E-02, &
     3.4800E-02,3.1500E-02,2.9700E-02,2.7900E-02,2.6200E-02,2.5000E-02, &
     2.2900E-02,2.1000E-02,1.9300E-02,1.7700E-02,1.6300E-02,1.5100E-02, &
     1.3800E-02,1.2800E-02,1.1800E-02,1.1000E-02,1.0100E-02,9.4100E-03, &
     8.6600E-03,8.0700E-03,7.3700E-03,6.8300E-03,6.2500E-03,5.7900E-03, &
     5.3800E-03,5.0600E-03,4.7300E-03,4.4900E-03,4.2400E-03,4.0500E-03, &
     3.8900E-03,3.7600E-03,3.6300E-03,3.5500E-03,3.4700E-03,3.4000E-03, &
     3.3500E-03,3.3600E-03,3.3500E-03,3.3900E-03,3.4000E-03,3.4800E-03, &
     3.5200E-03,3.6300E-03,3.7000E-03,3.7800E-03,3.8900E-03,3.9900E-03/ 
      DATA CNTABW(265:330) / &
     4.1000E-03,4.2200E-03,4.3300E-03,4.5000E-03,4.6500E-03,4.7900E-03, &
     4.9400E-03,5.1200E-03,5.3100E-03,5.4900E-03,5.6800E-03,5.8600E-03, &
     6.0800E-03,6.3100E-03,6.5300E-03,6.7300E-03,6.9600E-03,7.2200E-03, &
     7.4900E-03,7.7900E-03,8.0600E-03,8.3300E-03,8.6400E-03,8.9600E-03, &
     9.2700E-03,9.6600E-03,1.0000E-02,1.0400E-02,1.0800E-02,1.1200E-02, &
     1.1700E-02,1.2200E-02,1.2600E-02,1.3100E-02,1.3600E-02,1.4000E-02, &
     1.4500E-02,1.4900E-02,1.5200E-02,1.5400E-02,1.5600E-02,1.5700E-02, &
     1.5700E-02,1.5700E-02,1.5500E-02,1.5300E-02,1.5100E-02,1.4800E-02, &
     1.4600E-02,1.4300E-02,1.4000E-02,1.3700E-02,1.3300E-02,1.2900E-02, &
     1.2600E-02,1.2200E-02,1.1800E-02,1.1500E-02,1.1000E-02,1.0800E-02, &
     1.0500E-02,1.0300E-02,1.0100E-02,1.0000E-02,9.9300E-03,9.9000E-03/ 
      DATA CNTABW(331:396) / &
     9.9500E-03,1.0000E-02,1.0200E-02,1.0400E-02,1.0700E-02,1.1000E-02, &
     1.1500E-02,1.2000E-02,1.2800E-02,1.3800E-02,1.5000E-02,1.6600E-02, &
     1.8500E-02,2.0500E-02,2.4200E-02,2.9300E-02,3.3200E-02,4.2900E-02, &
     5.4400E-02,6.8800E-02,8.4000E-02,1.0210E-01,1.1700E-01,1.3000E-01, &
     1.3200E-01,1.2400E-01,1.0600E-01,8.8000E-02,7.4000E-02,6.1800E-02, &
     5.3500E-02,4.8400E-02,4.4700E-02,4.2000E-02,3.9800E-02,3.8300E-02, &
     3.7300E-02,3.7000E-02,3.6600E-02,3.6300E-02,3.6000E-02,3.5700E-02, &
     3.5500E-02,3.5200E-02,3.5000E-02,3.4700E-02,3.4600E-02,3.4300E-02, &
     3.4200E-02,3.4200E-02,3.4200E-02,3.4300E-02,3.4200E-02,3.4200E-02, &
     3.4200E-02,3.4200E-02,3.4200E-02,3.4400E-02,3.4500E-02,3.4600E-02, &
     3.4900E-02,3.5100E-02,3.5100E-02,3.5100E-02,3.5200E-02,3.5600E-02/ 
      DATA CNTABW(397:462) / &
     3.5900E-02,3.6100E-02,3.6200E-02,3.6600E-02,3.7000E-02,3.7400E-02, &
     3.7800E-02,3.8300E-02,3.8700E-02,3.9200E-02,3.9800E-02,4.0500E-02, &
     4.1100E-02,4.1700E-02,4.2400E-02,4.3400E-02,4.4300E-02,4.5300E-02, &
     4.6700E-02,4.8100E-02,4.9700E-02,5.1500E-02,5.3400E-02,5.5700E-02, &
     5.8900E-02,6.2200E-02,6.6100E-02,7.0700E-02,7.6400E-02,8.2800E-02, &
     8.9800E-02,9.7300E-02,1.0700E-01,1.1800E-01,1.3000E-01,1.4400E-01, &
     1.5900E-01,1.7600E-01,1.9200E-01,2.0800E-01,2.2600E-01,2.4300E-01, &
     2.6000E-01,2.7700E-01,2.9200E-01,3.0500E-01,3.1700E-01,3.2800E-01, &
     3.3800E-01,3.4700E-01,3.5600E-01,3.6500E-01,3.7300E-01,3.7900E-01, &
     3.8600E-01,3.9200E-01,3.9700E-01,4.0300E-01,4.0800E-01,4.1200E-01, &
     4.1700E-01,4.2000E-01,4.2300E-01,4.2500E-01,4.2700E-01,4.2800E-01/ 
      DATA CNTABW(463:518) / &
     4.2700E-01,4.2700E-01,4.2600E-01,4.2500E-01,4.2300E-01,4.2100E-01, &
     4.1800E-01,4.1500E-01,4.1100E-01,4.0800E-01,4.0400E-01,4.0100E-01, &
     3.9700E-01,3.9400E-01,3.9000E-01,3.8600E-01,3.8200E-01,3.7700E-01, &
     3.7200E-01,3.6800E-01,3.6300E-01,3.5900E-01,3.5600E-01,3.5200E-01, &
     3.5300E-01,3.5700E-01,3.6100E-01,3.6800E-01,3.7500E-01,3.8500E-01, &
     3.9800E-01,4.1400E-01,4.3600E-01,4.6900E-01,5.0500E-01,5.3900E-01, &
     5.7100E-01,5.9700E-01,6.1800E-01,6.2900E-01,6.2200E-01,6.0800E-01, &
     5.9300E-01,5.7700E-01,5.5700E-01,5.3200E-01,5.0700E-01,4.8700E-01, &
     4.6600E-01,4.5000E-01,4.4400E-01,4.3800E-01,4.6000E-01,5.2700E-01, &
     7.1800E-01,8.4657E-01/ 
      DATA PI/3.14159265/

  !  ZERO PARAMETERS
  RN=0.0
  CN=0.0
  ABSIND=0.0
  ABSCOF=0.0

  !  CONVERT WAVELENGTH TO MICRONS
  WL=XLAM
  IF(IUNIT.EQ.1)WL=1000*WL
  IF(IUNIT.EQ.2)WL=10000*WL
  IF(IUNIT.EQ.3)WL=10000*(1.0/WL)
  IF(WL.LT.WLMIN.OR.WL.GT.WLMAX)RETURN

  !  REGION FROM 0.2 MICRON TO 1000.0 MICRON  -  TABLE LOOKUP
      IF (WL.GT.CUTWAT) GO TO 3
      DO 1 I=2,NUMWAT
      IF(WL.GT.WLTABW(I))GO TO 1
      I1=I-1
      I2=I
      GO TO 2
 1    CONTINUE
      I1=NUMWAT-1
      I2=NUMWAT
 2    FAC=(WL-WLTABW(I1))/(WLTABW(I2)-WLTABW(I1))
      RN=RNTABW(I1)+FAC*(RNTABW(I2)-RNTABW(I1))
      CN=CNTABW(I1)+FAC*(CNTABW(I2)-CNTABW(I1))
      GO TO 5

  !  REGION FROM 0.1 CM TO 10 CM
  !  EXTENSION OF DEBYE THEOREY BASED ON THE WORK OF
  !     COLE,K.S.,AND R.H.COLE,1941.JOUR.CHEM.PHYS.,9,P 341.

  !  DEFINE TEMPERATURE TERMS AND WAVELENGTH IN CM
 3    TC=T-273.15
      T1=TC+273.0
      T2=TC-25.0
      XL=WL/10000.0

  !  DEFINE FREQUENCY INDEPENDENT CONDUCTIVITY(SIGMA) AND
  !  SPREAD PARAMETER(ALPHA)
  !  IN CLASSICAL DEBYE THEOREY THESE TERMS ARE ZERO

  !  SIGMA GIVEN BY SAXTON,J.A.,1949.WIRELESS ENGINEER,26,P 288.
  !  ALPHA GIVEN BY RAY ( EQUATION 7B )
      SIGMA=12.5664E8
      ALPHA=-16.8129/T1+0.0609265

  !  DEFINE STATIC DIELECTRIC CONSTANT(ES) - RAY EQN 4
  !         HIGH FREQUENCY DIELECTRIC CONSTANT(E00) - RAY EQN 7A
  !         RELAXTION WAVELENGTH IN CM(XLAMS) - RAY EQN 7C

  !  TEMPERATURE DEPENDENCE OF ES GIVEN BY
  !     WYMAN,J.,AND E.N.INGALLS,1938.JOUR.AM.CHEM.SOC.,60,P 1182.
      ES=78.54*(1.0-4.579E-3*T2+1.19E-5*T2*T2-2.8E-8*T2*T2*T2)
      E00=5.27137+0.0216474*TC-0.00131198*TC*TC
      XLAMS=0.00033836*EXP(2513.98/T1)

  !  CALCULATE EXPRESSIONS USED FOR DIELECTRIC CONSTANT
      TERM=PI*ALPHA/2
      SINT=SIN(TERM)
      COST=COS(TERM)
      XLRAT=XLAMS/XL
      POWTRM=XLRAT**(1-ALPHA)
      DENOM=1.0+2*POWTRM*SINT+XLRAT**(2.0*(1-ALPHA))

  !  CALCULATION OF DIELECTRIC CONSTANT
  !  REAL PART - RAY EQN 5
      ER=E00+(ES-E00)*(1.0+POWTRM*SINT)/DENOM
  !  IMAGINARY PART OR LOSS TERM - RAY EQN 6
      EI=(SIGMA*XL/18.8496E10)+(ES-E00)*POWTRM*COST/DENOM
  !  COMPLEX PERMITTIVITY
      E=CMPLX(ER,-EI)

  !  COMPLEX INDEX OF REFRACTION - RAY EQN 1
      M=CSQRT(E)
      RN=REAL(M)
      CN=-AIMAG(M)

  !  CORRECTION TO IMAGINARY INDEX TO ACCOUNT FOR THE
  !  REMAINING ABSORPTION BANDS - RAY EQN 8(TABLE 2)
      IF (WL.GT.3000.0) GO TO 5
      CN = CN + 0.39*EXP(-ABS(ALOG10(WL/17.0)/0.45)**1.3) &
              + 0.41*EXP(-ABS(ALOG10(WL/62.0)/0.35)**1.7) &
              + 0.25*EXP(-ABS(ALOG10(WL/300.)/0.47)**3.0)

  !  ABSORPTIVE QUANTITIES
 5    ABSIND=CN/RN
      ABSCOF=4.0*PI*CN/WL
END SUBROUTINE REFWAT


END MODULE MIE_GAMMA_WATER


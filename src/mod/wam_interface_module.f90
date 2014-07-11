MODULE WAM_INTERFACE_MODULE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!          THIS MODULE COLLECTS PROCEDURES, WHICH ARE USED IN THE WAM MODEL    !
!          TO COMPUTE PARAMETERS FROM SPECTRA AND TO INTERPOLATE SPECTRA.      !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE, ONLY:  &
&         AKI                     !! WAVE NUMBER FROM FREQUENY AND DEPTH

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FRE_DIR_MODULE, ONLY: ML, KL, FR, CO, TH, DELTH, SINTH, COSTH,         &
&                             DF, DF_FR, DF_FR2,                               &
&                             DFIM, DFIMOFR, DFIM_FR, DFIM_FR2,                &
&                             MO_TAIL, MM1_TAIL, MP1_TAIL, MP2_TAIL,           &
&                             TFAK, TFAC_ST

USE WAM_GENERAL_MODULE, ONLY: G, PI, ZPI, DKMAX, ZPISQRT,                      &
&                             BF2MAX, BF2MIN, C4MAX, C4MIN

IMPLICIT NONE

PRIVATE

REAL      :: EMIN = 1.0E-12    !! REPLACES THE INTRINSIC TINY


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     D.  GENERIC INTERFACES (THIS MODULE CONTAINS THE PROCEDURES).            !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE CHARNOCK_PAR                     !! CHARNOCK PARAMETER.
MODULE PROCEDURE CHARNOCK_PAR
END INTERFACE
PUBLIC CHARNOCK_PAR

INTERFACE FEMEAN                           !! MEAN FREQUENCY.
   MODULE PROCEDURE FEMEAN
END INTERFACE
PUBLIC FEMEAN

INTERFACE INTSPEC                          !! INTERPOLATES SPECTRA.
   MODULE PROCEDURE INTSPEC
END INTERFACE
PUBLIC INTSPEC

INTERFACE MEAN_DIRECTION                   !! MEAN DIRECTION AND SPREAD
   MODULE PROCEDURE MEAN_DIRECTION_1       !! SCALAR VERSION
   MODULE PROCEDURE MEAN_DIRECTION_B       !! VECTOR VERSION
END INTERFACE
PUBLIC MEAN_DIRECTION

INTERFACE MEANSQS                          !! COMPUTATES MEAN SQUARE SLOPE.
   MODULE PROCEDURE MEANSQS
END INTERFACE
PUBLIC MEANSQS

INTERFACE PEAK_PERIOD                      !! COMPUTATES PEAK PERIOD.
   MODULE PROCEDURE PEAK_PERIOD_1          !! SCALAR VERSION
   MODULE PROCEDURE PEAK_PERIOD_B          !! VECTOR VERSION
END INTERFACE
PUBLIC PEAK_PERIOD

INTERFACE ROTSPEC                          !! ROTATE A SPECTRUM.
   MODULE PROCEDURE ROTSPEC
END INTERFACE
PUBLIC ROTSPEC

INTERFACE STRSPEC                          !! STRETCH A SPECTRUM.
   MODULE PROCEDURE STRSPEC
END INTERFACE
PUBLIC STRSPEC

INTERFACE STOKES_DRIFT                     !! COMPUTATES STOKES DRIFT.
   MODULE PROCEDURE STOKES_DRIFT
END INTERFACE
PUBLIC STOKES_DRIFT

INTERFACE TM1_TM2_PERIODS                  !! COMPUTATES TM1 AND/OR TM2 PERIODS.
   MODULE  PROCEDURE TM1_TM2_PERIODS_1     !! SCALAR VERSION
   MODULE  PROCEDURE TM1_TM2_PERIODS_B     !! VECTOR VERSION
END INTERFACE
PUBLIC TM1_TM2_PERIODS

INTERFACE TOTAL_ENERGY                     !! COMPUTES TOTAL ENERGY.
   MODULE  PROCEDURE TOTAL_ENERGY_1        !! SCALAR VERSION
   MODULE  PROCEDURE TOTAL_ENERGY_B        !! VECTOR VERSION
END INTERFACE
PUBLIC TOTAL_ENERGY

INTERFACE COS2_SPR                         !! COSINE SQUARE SPREAD.
   MODULE  PROCEDURE COS2_SPR_1            !! SCALAR VERSION
   MODULE  PROCEDURE COS2_SPR_B            !! VECTOR VERSION
END INTERFACE
PUBLIC COS2_SPR

INTERFACE WM1_WM2_WAVENUMBER               !! WM1 AND/OR WM2 WAVENUMBER.
   MODULE  PROCEDURE WM1_WM2_WAVENUMBER_1  !! SCALAR VERSION
   MODULE  PROCEDURE WM1_WM2_WAVENUMBER_B  !! VECTOR VERSION
END INTERFACE
PUBLIC WM1_WM2_WAVENUMBER

INTERFACE KURTOSIS                         !! KURTOSIS, BENJAMIN-FEIR INDEX
   MODULE  PROCEDURE KURTOSIS              !! GODA'S PEAKEDNESS PARAMETER AND
END INTERFACE                              !! NORMALIZED MAXIMUM WAVE HEIGHT.
PUBLIC KURTOSIS

INTERFACE PEAK_FREQ                        !! PARAMETERS AT PEAK OF SPECTRUM.
   MODULE  PROCEDURE PEAK_FREQ
END INTERFACE
PUBLIC PEAK_FREQ

INTERFACE H_MAX                            !! EXPECTED MAXIMUM WAVE HEIGHT.
   MODULE  PROCEDURE H_MAX
END INTERFACE
PUBLIC H_MAX

INTERFACE TRANSF2                          !! NARROW BAND LIMIT BENJAMIN-FEIR
   MODULE  PROCEDURE TRANSF2               !! INDEX FOR THE FINITE DEPTH CASE.
END INTERFACE
PUBLIC TRANSF2

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CHARNOCK_PAR (USTAR, Z0, BETA)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   CHARNOCK_PAR - DETERMINES THE CHARNOCK PARAMETER.                          !
!                                                                              !
!     P.JANSSEN      KNMI/ECMWF  JANUARY 1992                                  !
!     J.BIDLOT       ECMWF       FEBRUARY 1996  MESSAGE PASSING                !
!     J.BIDLOT       ECMWF       AUGUST 2008  REMOVE MESSAGE PASSING           !
!                                                                              !
!     PURPOSE.                                                                 !
!                                                                              !
!       COMPUTES THE CHARNOCK PARAMETER FOR ATMOSPHERIC MODEL.                 !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)            :: USTAR(:)      !! FRICTION VELOCITY IN M/S.
REAL,    INTENT(IN)            :: Z0(:)         !! ROUGHNESS LENGTH IN M.
REAL,    INTENT(OUT)           :: BETA(:)       !! CHARNOCK PARAMETER.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: EPSUS = 1.0E-6

! ---------------------------------------------------------------------------- !

BETA(:) = G*Z0(:)/MAX(USTAR(:)**2,EPSUS)

END SUBROUTINE CHARNOCK_PAR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FEMEAN (F, EMEAN, FM, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   FEMEAN - COMPUTATION OF MEAN FREQUENCY.                                    !
!                                                                              !
!     S.D. HASSELMANN                                                          !
!     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER                              !
!     H. GUNTHER     GKSS         DECEMBER 2001    FT90                        !
!                                                                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)            :: F(:,:,:)      !! BLOCK OF SPECTRA.
REAL,    INTENT(IN)            :: EMEAN(:)      !! TOTAL ENERGY.
REAL,    INTENT(OUT)           :: FM   (:)      !! MEAN FREQUENCY.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK (:,:,:)  !! INTERATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL     :: TEMP2(SIZE(F,1),SIZE(F,3))
INTEGER :: IJ, M

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP2 = SUM(F,DIM=2, MASK=MASK)
ELSE
   TEMP2 = SUM(F,DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !

DO IJ = 1,SIZE(F,1)
   FM(IJ) = MM1_TAIL*TEMP2(IJ,ML)  !! TAIL ENERGY
END DO

DO M = 1,ML
   DO IJ = 1,SIZE(F,1)
      FM(IJ) = FM(IJ) + TEMP2(IJ,M)*DFIMOFR(M)
   END DO
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. NORMALIZE.                                                            !
!        ----------                                                            !

DO IJ = 1, SIZE(F,1)
   FM(IJ) = EMEAN(IJ)/MAX(FM(IJ),EMIN)
END DO

END SUBROUTINE FEMEAN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE INTSPEC (DEL12, DEL1L,  F1, FMEAN1, EMEAN1, THETM1,                 &
&                   F2, FMEAN2, EMEAN2, THETM2, FL, FMEAN, EMEAN, THETM )

! ---------------------------------------------------------------------------- !
!                                                                              !
!   INTSPEC  -  INTERPOLATION OF SPECTRA.                                      !
!                                                                              !
!     SUSANNE HASSELMANN  MPI        JUNE 1990.                                !
!     H. GUNTHER          GKSS/ECMWF JAN. 1991   MODIFIED FOR CYCLE_4          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       INTERPOLATION OF SPECTRA.                                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       ROTATE SPECTRA ACCORDING TO MEAN OF MEAN ANGLES, TRANSFORM             !
!       FREQUENCIES ACCORDING TO MEAN OF MEAN FREQUENCIES ,ADJUST ENERGY       !
!       ACCORDCING TO MEAN OF TOTAL ENERGY AND INTERPOLATE RESULTING           !
!       SPECTRA.                                                               !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       K.HASSELMANN, 1990,                                                    !
!          INTERPOLATION OF WAVE SPECTRA. WAM NOTE 6/6/90.                     !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)  :: DEL12         !! DISTANCE SPECTRUM 2 - SPECTRUM 1.
REAL,    INTENT(IN)  :: DEL1L         !! DISTANCE SPECTRUM L - SPECTRUM 1.
REAL,    INTENT(IN)  :: F1(:,:)       !! SPECTRUM 1.
REAL,    INTENT(IN)  :: FMEAN1        !! MEAN FREQUENCY OF F1.
REAL,    INTENT(IN)  :: EMEAN1        !! MEAN ENERGY OF F1.
REAL,    INTENT(IN)  :: THETM1        !! MEAN DIRECTION OF F1.
REAL,    INTENT(IN)  :: F2(:,:)       !! SPECTRUM 2.
REAL,    INTENT(IN)  :: FMEAN2        !! MEAN FREQUENCY OF F2.
REAL,    INTENT(IN)  :: EMEAN2        !! MEAN ENERGY OF F2.
REAL,    INTENT(IN)  :: THETM2        !! MEAN DIRECTION OF F2.
REAL,    INTENT(OUT) :: FL(: ,:)      !! INTEPOLATED SPECTRUM.
REAL,    INTENT(OUT) :: FMEAN         !! INTEPOLATED MEAN FREQUENCY.
REAL,    INTENT(OUT) :: EMEAN         !! INTEPOLATED MEAN ENERGY.
REAL,    INTENT(OUT) :: THETM         !! INTEPOLATED MEAN DIRECTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL :: GW1, GW2
REAL :: F_L(SIZE(F1,1),SIZE(F1,2)), F_R(SIZE(F1,1),SIZE(F1,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTERPOLATION WEIGHTS.                                                !
!        ----------------------                                                !

GW2 = DEL1L/DEL12
GW1 = 1. - GW2

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LEFT WEIGHT OR ENERGY OF LEFT SPECTRUM IS ZERO.                       !
!        -----------------------------------------------                       !

IF (ABS(GW1).LT.EPSILON(1.) .OR. EMEAN1.LT.EPSILON(1.)) THEN
   FL = GW2*F2
   EMEAN = GW2*EMEAN2
   FMEAN = GW2*FMEAN2
   THETM = GW2*THETM2
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. RIGHT WEIGHT OR ENERGY OF RIGHT SPECTRUM IS ZERO.                     !
!        -------------------------------------------------                     !

IF (ABS(GW2).LT.EPSILON(1.) .OR. EMEAN2.LT.EPSILON(1.)) THEN
   FL = GW1*F1
   EMEAN = GW1*EMEAN1
   FMEAN = GW1*FMEAN1
   THETM = GW1*THETM1
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ENERGY AND WEIGHTS OF BOTH SPECTRA ARE GT ZERO.                       !
!        -----------------------------------------------                       !

!     3.1 INTERPOLATE MEAN VALUES.                                             !

EMEAN = GW1*EMEAN1+GW2*EMEAN2
FMEAN = GW1*FMEAN1+GW2*FMEAN2
THETM = ATAN2 (GW1*SIN(THETM1)+GW2*SIN(THETM2),GW1*COS(THETM1)+GW2*COS(THETM2))

!     3.2 ADJUST LEFT SPECTRUM TO MEAN VALUES.                                 !

CALL ROTSPEC (F1, FL, THETM-THETM1)        !! ROTATE.
CALL STRSPEC (FL, F_L, FMEAN1/FMEAN)       !! STRETCH.
GW1 = GW1*EMEAN/EMEAN1                     !! ADJUST ENERGY.

!    3.3 ADJUST RIGHT SPECTRUM TO MEAN VALUES.                                 !

CALL ROTSPEC (F2, FL, THETM-THETM2)        !! ROTATE.
CALL STRSPEC (FL, F_R, FMEAN2/FMEAN)       !! STRETCH.
GW2 = GW2*EMEAN/EMEAN2                     !! ADJUST ENERGY.

!      3.4 LINEAR INTERPOLATION TO NEW SPECTRA.                                !

FL = GW1*F_L + GW2*F_R

END SUBROUTINE INTSPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MEAN_DIRECTION_B (F3, THQ, SPREAD, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!  MEAN_DIRECTION_B - COMPUTATION OF MEAN WAVE DIRECTION FOR BLOCK.            !
!                                                                              !
!     S.D. HASSELMANN                                                          !
!     OPTIMIZED BY L. ZAMBRESKY                                                !
!     MODIFIED FOR K-MODEL BY C.SCHNEGGENBURGER                                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE MEAN WAVE DIRECTION FROM ENERGY DENSITY AT EACH GRID POINT. !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGRATION OF SPECTRUM TIMES SIN AND COS OVER DIRECTION.              !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLE                                                       !

REAL,    INTENT(IN)            :: F3(:,:,:)   !! BLOCK OF DENSITY SPECTRA.
REAL,    INTENT(OUT), OPTIONAL :: THQ(:)      !! MEAN DIRECTION [RAD].
REAL,    INTENT(OUT), OPTIONAL :: SPREAD(:)   !! MEAN SPREAD [RAD].
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLE                                                           !

INTEGER         :: K, M, IJ
REAL            :: SI(1:SIZE(F3,1)), CI(1:SIZE(F3,1)), TEMP(1:SIZE(F3,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALISE SIN AND COS ARRAYS.                                        !
!        ------------------------------                                        !

SI = 0.
CI = 0.

IF (PRESENT(SPREAD)) SPREAD = 0.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   DO K = 1,KL
      TEMP = 0.
      DO M = 1,ML
         WHERE (MASK(:,K,M)) TEMP(:) = TEMP(:) + F3(:,K,M)*DFIM(M)
      END DO
      DO IJ = 1,SIZE(F3,1)
         SI(IJ) = SI(IJ) + SINTH(K)*TEMP(IJ)
         CI(IJ) = CI(IJ) + COSTH(K)*TEMP(IJ)
      END DO
      IF (PRESENT(SPREAD)) SPREAD(:) = SPREAD(:) + TEMP(:)
   END DO
ELSE
   DO K = 1,KL
      DO IJ = 1,SIZE(F3,1)
         TEMP(IJ) = F3(IJ,K,1)*DFIM(1)
      END DO
      DO M = 2,ML
         DO IJ = 1,SIZE(F3,1)
            TEMP(IJ) = TEMP(IJ) + F3(IJ,K,M)*DFIM(M)
         END DO
      END DO

      DO IJ = 1,SIZE(F3,1)
         SI(IJ) = SI(IJ) + SINTH(K)*TEMP(IJ)
         CI(IJ) = CI(IJ) + COSTH(K)*TEMP(IJ)
      END DO
      IF (PRESENT(SPREAD)) SPREAD(:) = SPREAD(:) + TEMP(:)
   END DO
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTE MEAN DIRECTION.                                               !
!        -----------------------                                               !

IF (PRESENT(THQ)) THEN
   WHERE (CI.EQ.0.) CI = 0.1E-30
   THQ = ATAN2(SI,CI)
   WHERE (THQ.LT.0.) THQ = THQ + ZPI
   WHERE (THQ.GT.ZPI-0.001) THQ = 0.
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTE MEAN SPREAD.                                                  !
!        --------------------                                                  !

IF (PRESENT(SPREAD)) THEN
   WHERE (ABS(CI) .LT. 0.1E-15) CI = SIGN(0.1E-15,CI)
   WHERE (ABS(SI) .LT. 0.1E-15) SI = SIGN(0.1E-15,SI)
   WHERE (ABS(SPREAD) .LT. 0.1E-15) SPREAD = SIGN(0.1E-15,SPREAD)
   SPREAD = 2.*(1.-SQRT(SI**2 + CI**2)/SPREAD)
   WHERE (SPREAD.LE.0)
      SPREAD = TINY(1.)
   ELSEWHERE
      SPREAD = SQRT(SPREAD)
   END WHERE
END IF

END SUBROUTINE MEAN_DIRECTION_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MEAN_DIRECTION_1 (F3, THQ, SPREAD)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MEAN_DIRECTION_1 - COMPUTATION OF MEAN WAVE DIRECTION ONE SPECTRUM.        !
!                                                                              !
!     S.D. HASSELMANN                                                          !
!     OPTIMIZED BY L. ZAMBRESKY                                                !
!     MODIFIED FOR K-MODEL BY C.SCHNEGGENBURGER                                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE MEAN WAVE DIRECTION FROM ONE SPECTRUM.                      !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGRATION OF SPECTRUM TIMES SIN AND COS OVER DIRECTION.              !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLE                                                       !

REAL, INTENT(IN)            :: F3(:,:)  !! DENSITY SPECTRUM.
REAL, INTENT(OUT), OPTIONAL :: THQ      !! MEAN DIRECTION [RAD].
REAL, INTENT(OUT), OPTIONAL :: SPREAD   !! MEAN SPREAD [RAD].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLE                                                           !

INTEGER  :: K, M
REAL     :: SI, CI
REAL     :: TEMP(1:SIZE(F3,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.                            !
!        ------------------------------------------                            !

DO K = 1,KL
   TEMP(K) = F3(K,1)*DFIM(1)
END DO

DO M = 2,ML
   DO K = 1,KL
      TEMP(K) = TEMP(K) + F3(K,M)*DFIM(M)
   END DO
END DO

SI = TEMP(1)*SINTH(1)
CI = TEMP(1)*COSTH(1)
DO K = 2,KL
   SI = SI + TEMP(K)*SINTH(K)
   CI = CI + TEMP(K)*COSTH(K)
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN DIRECTION.                                               !
!        -----------------------                                               !

IF (PRESENT(THQ)) THEN
   IF (CI.EQ.0.) CI = 0.1E-30
   THQ = ATAN2(SI,CI)
   IF (THQ.LT.0.) THQ = THQ + ZPI
   IF (THQ.GT.ZPI-0.001) THQ = 0.
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTE MEAN SPREAD.                                                  !
!        --------------------                                                  !

IF (PRESENT(SPREAD)) THEN
   SPREAD = SUM(TEMP)
   IF (ABS(CI) .LT. 0.1E-15) CI = SIGN(0.1E-15,CI)
   IF (ABS(SI) .LT. 0.1E-15) SI = SIGN(0.1E-15,SI)
   IF (ABS(SPREAD) .LT. 0.1E-15) SPREAD = SIGN(0.1E-15,SPREAD)
   SPREAD = 2.*(1.-SQRT(SI**2 + CI**2)/SPREAD)
   IF (SPREAD.LE.0) THEN
      SPREAD = TINY(1.)
   ELSE
      SPREAD = SQRT(SPREAD)
   END IF
END IF

END SUBROUTINE MEAN_DIRECTION_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MEANSQS (F, USTAR, FM, SM, INDEP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MEANSQS - COMPUTATION OF MEAN SQUARE SLOPE.
!                                                                              !
!     P.A.E.M. JANSSEN
!     J. BIDLOT  ECMWF  FEBRUARY 1996  MESSAGE PASSING
!     H. GUNTHER   HZG  JUNE 2012                                              !
!                                                                              !
!    PURPOSE.
!     --------
!                                                                              !
!       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT.
!                                                                              !
!     METHOD.
!     -------
!                                                                              !
!       NONE.
!                                                                              !
!     EXTERNALS.
!     ----------
!                                                                              !
!       NONE.
!                                                                              !
!     REFERENCE.
!     ----------
!                                                                              !
!       NONE.
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !

REAL,    INTENT(IN)            :: F(:,:,:)    !! BLOCK OF DENSITY SPECTRA.
REAL,    INTENT(IN)            :: USTAR(:)    !! FRICTION VELOCITIES [M/S].
REAL,    INTENT(IN)            :: FM(:)       !! MEAN WAVE FREQUENCY.
REAL,    INTENT(OUT)           :: SM(:)       !! MEAN SQUARE SLOPE.
INTEGER, INTENT(IN), OPTIONAL  :: INDEP(:)    !! DEPTH TABLE INDEX.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !

REAL, PARAMETER :: XLAMBDAC = 0.0628
REAL, PARAMETER :: SURFT    = 0.000075
REAL, PARAMETER :: EPSMIN = 0.1E-32

INTEGER :: M
REAL :: FS, XKC, FC, CONST1, CONST2

REAL, DIMENSION(SIZE(F,3))           ::  FD
REAL, DIMENSION(SIZE(F,1),SIZE(F,3)) ::  TEMP   !! FREQUENCY SPECTRA


! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

TEMP(:,:) = SUM(F(:,:,:), DIM=2)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALISE MEAN SLOPE ARRAY WITH TAIL FACTOR.
!        ---------------------------------------------                         !

FS     = FR(ML)
XKC    = ZPI/XLAMBDAC
FC     = SQRT(G*XKC+SURFT*XKC**3)/ZPI
CONST1 = LOG(FC/FS)
CONST2 = CONST1*ZPI**4*FS**5/G**2*DELTH

SM(:)    = MAX(CONST2*TEMP(:,ML), EPSMIN)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

IF (PRESENT(INDEP)) THEN     !! SHALLOW WATER INTEGRATION.

   DO M = 1,ML
      SM(:) = SM(:) + DFIM(M)*TFAK(INDEP(:),M)**2*TEMP(:,M)
   END DO

ELSE                         !! DEEP WATER INTEGRATION.

   FD(1:ML) = DFIM(1:ML)*(ZPI*FR(1:ML))**4/G**2
   DO M=1,ML
      SM(:) = SM(:) + FD(M)*TEMP(:,M)
   END DO

END IF

END SUBROUTINE MEANSQS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PEAK_PERIOD_B (F, PEAKP, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    PEAK_PERIOD_B - COMPUTATES PEAK PERIOD (VECTOR VERSION).                  !
!                                                                              !
!     H. GUNTHER      ECMWF            DECEMBER 1989                           !
!     (CODE REMOVED FROM SUB. FEMEAN)                                          !
!     H. GUNTHER      GKSS            FEBRUARY 2002  CHANGED TO PERIOD.        !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE PEAK PERIOD AT EACH GRID POINT.                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE FREQUENCY INDEX OF THE 1-D SPECTRA ARE COMPUTED AND                !
!       CONVERTED TO PERIODS.                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)            :: F(:,:,:)    !! BLOCK OF SPECTRA.
REAL,    INTENT(OUT)           :: PEAKP(:)    !! PEAK PERIODS.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: IJ
INTEGER  :: IPEAK(SIZE(F,1))
REAL     :: EED1D(SIZE(F,1),SIZE(F,3))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE 1-D SPECTRUM (WITHOUT DELTA THETA).                           !
!        -------------------------------------------                           !

IF (PRESENT(MASK)) THEN
   EED1D = SUM(F, DIM=2, MASK=MASK)
ELSE
   EED1D = SUM(F, DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. DEFINE PEAK INDEX.                                                    !
!        ------------------                                                    !

DO IJ = 1,SIZE(F,1)
   IPEAK(IJ:IJ) = MAXLOC(EED1D(IJ,:))
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CALCULATE PEAK PERIOD FROM PEAK INDEX.                                !
!        --------------------------------------                                !

PEAKP = 1./FR(IPEAK)
WHERE (IPEAK.EQ.1) PEAKP = 1.

END SUBROUTINE PEAK_PERIOD_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PEAK_PERIOD_1 (F, PEAKP, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    PEAK_PERIOD_1 - COMPUTATES PEAK PERIOD (SCALAR VERSION).                  !
!                                                                              !
!     H. GUNTHER      ECMWF            DECEMBER 1989                           !
!     (CODE REMOVED FROM SUB. FEMEAN)                                          !
!     H. GUNTHER      GKSS            FEBRUARY 2002  CHANGED TO PERIOD.        !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE PEAK PERIOD.                                                   !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE FREQUENCY INDEX OF THE 1-D SPECTRUM IS COMPUTED AND                !
!       CONVERTED TO PERIOD.                                                   !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)            :: F(:,:)     !! SPECTRUM.
REAL,    INTENT(OUT)           :: PEAKP      !! PEAK PERIOD.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:)  !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: IPEAK(1:1)
REAL     :: EED1D(SIZE(F,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE 1-D SPECTRUM (WITHOUT DELTA THETA).                           !
!        -------------------------------------------                           !

IF (PRESENT(MASK)) THEN
   EED1D = SUM(F, DIM=1, MASK=MASK)
ELSE
   EED1D = SUM(F, DIM=1)
END IF


! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. DEFINE PEAK INDEX.                                                    !
!        ------------------                                                    !

IPEAK(1:1) = MAXLOC(EED1D)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CALCULATE PEAK PERIOD FROM PEAK INDEX.                                !
!        --------------------------------------                                !

PEAKP = 1./FR(IPEAK(1))
IF (IPEAK(1).EQ.1) PEAKP = 1.

END SUBROUTINE PEAK_PERIOD_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE ROTSPEC (F_IN, F_OUT, RTHET)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   ROTSPEC - ROUTINE TO ROTATE THE SPECTRUM.                                  !
!                                                                              !
!     EVA BAUER      MPI  HAMBURG    MAY 1990.                                 !
!     H. GUNTHER          GKSS/ECMWF JAN. 1991   MODIFIED FOR CYCLE_4          !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO ROTATE THE SPECTRUM.                                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLE.                                                      !
!     -------------------                                                      !

REAL,    INTENT(IN)  :: F_IN(:,:)      !! SPECTRUM TO BE ROTATED.
REAL,    INTENT(OUT) :: F_OUT(:,:)     !! ROTATED SPECTRUM.
REAL,    INTENT(IN) ::  RTHET          !! TURNING ANGLE [RAD], CLOCKWISE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER         :: INC
REAL            :: ADIF, BDIF

! ---------------------------------------------------------------------------- !

ADIF = RTHET * REAL(SIZE(F_IN,1)) / ZPI
INC = -FLOOR(ADIF)
ADIF = ADIF + REAL(INC)
BDIF = 1. - ADIF

F_OUT = BDIF * CSHIFT(F_IN, INC) + ADIF * CSHIFT(F_IN, INC-1)

END SUBROUTINE ROTSPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE STRSPEC (F_IN, F_OUT, GAMMA)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   STRSPEC - ROUTINE TO STRETCH A SPECTRUM.                                   !
!                                                                              !
!      EVA BAUER      MPI  HAMBURG    MAY 1990.                                !
!      H. GUNTHER     GKSS/ECMWF      JAN 1991  MODIFIED FOR CYCLE_4.          !
!      H. GUNTHER     GKSS            JAN 2002  FT90.                          !
!                                               ERROR FOR SHIFT TO HIGER       !
!                                               FREQUENCIES CORRECTED.         !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)  :: F_IN (:,:)    !! INPUT SPECTRUM.
REAL,    INTENT(OUT) :: F_OUT(:,:)    !! OUTPUT SPECTRUM.
REAL,    INTENT(IN)  :: GAMMA         !! STRETCHING PARAMETER.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER         :: ML, INC
REAL            :: ADIF, BDIF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALIZATION.                                                       !
!        ---------------                                                       !

IF (GAMMA.EQ.1.0) THEN
   F_OUT = F_IN
   RETURN
END IF

F_OUT = 0.0
ML = SIZE(F_IN,2)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. DETERMINE ACROSS HOW MANY FREQUENCY BINS THE  STRETCHING IS ACTING    !
!        AND THE INTERPOLATION WEIGHTS.                                        !
!        --------------------------------------------------------------------  !

INC = FLOOR(LOG10(GAMMA)/LOG10(CO))
IF (ABS(INC).GE.ML-1) RETURN      !! ENERGY IS SHIFTED OUT OF FREQUENCY RANGE

ADIF = (CO -GAMMA*CO**(-INC))/(CO-1.)
BDIF = 1. - ADIF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. STRECH SPECTRUM.                                                      !
!        ----------------                                                      !

IF (INC.GE.0) THEN

!     3.1 SHIFT TO LOWER FREQUENCIES.                                          !

   F_OUT(:,1:ML-INC-1) = ADIF*F_IN(:,1+INC:ML-1) + BDIF*F_IN(:,2+INC:ML)
   F_OUT(:,ML-INC)     = ADIF*F_IN(:,ML)
ELSE

!      3.2 SHIFT TO HIGHER FREQUENCIES.                                        !

   F_OUT(:,1-INC:ML) = ADIF*F_IN(:,1:ML+INC) + BDIF*F_IN(:,2:ML+INC+1)
   F_OUT(:,-INC)     = BDIF*F_IN(:,1)
END IF

END SUBROUTINE STRSPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE STOKES_DRIFT (F3, UST, VST, IN)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   STOKES_DRIFT - COMPUTES STOKES DRIFT FROM SPECTRA FOR DEEP WATER.          !
!                                                                              !
!     M. REISTAD AND O SAETRA     DNMI     AUGUST 1997                         !
!     H. GUNTHER    GKSS          NOVEMBER 2010  FT90.                         !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       STOKES DRIFT FIELDS ARE CREATED FROM SPECTRUM FIELDS.                  !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!        METHOD PROPOSED BY Kern E. Kenton                                     !
!                                                                              !
!     REFERENCES.                                                              !
!      -----------                                                             !
!                                                                              !
!       JGR, Vol 74 NO 28, 1969                                                !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !
    
REAL,              INTENT(IN)  :: F3(:,:,:)  !! BLOCK OF SPECTRA.
REAL,              INTENT(OUT) :: UST(:)     !! U COMPONENTS OF STOCKES DRIFT. 
REAL,              INTENT(OUT) :: VST(:)     !! V COMPONENTS OF STOCKES DRIFT.
INTEGER, OPTIONAL, INTENT(IN)  :: IN(:)      !! SHALLOW WATER TABLE INDEX

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: FAK_DEEP = 16.*PI**3/G
REAL            :: FAK, TAILFAC
REAL            :: FAKT(1:SIZE(F3,1))
REAL            :: SI(1:SIZE(F3,1)), CI(1:SIZE(F3,1))
INTEGER         :: IJ, K, M

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIAL.                                                              !
!        --------                                                              !

UST(:) = 0.
VST(:) = 0.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.                            !
!        ------------------------------------------                            !

IF (PRESENT(IN)) THEN     !! SHALLOW WATER
   FREQ_LOOP1: DO M = 1,ML
      FAKT(:)  = TFAC_ST(IN(:),M) * DFIM(M)
      DO IJ=1,SIZE(F3,1)
         SI(IJ) = F3(IJ,1,M)*SINTH(1)
         CI(IJ) = F3(IJ,1,M)*COSTH(1)
      END DO
      DIR_LOOP1: DO K=2,KL
         DO IJ=1,SIZE(F3,1)   
            SI(IJ) = SI(IJ) + F3(IJ,K,M)*SINTH(K)
            CI(IJ) = CI(IJ) + F3(IJ,K,M)*COSTH(K)
         END DO
      END DO DIR_LOOP1

      DO IJ=1,SIZE(F3,1)   
         SI(IJ) = FAKT(IJ) * SI(IJ) 
         CI(IJ) = FAKT(IJ) * CI(IJ)
         UST(IJ) = UST(IJ) + SI(IJ)
         VST(IJ) = VST(IJ) + CI(IJ)
      END DO
   END DO FREQ_LOOP1

ELSE                        !! DEEP WATER

   FREQ_LOOP2: DO M = 1,ML
      FAK   = FAK_DEEP * FR(M)**3 * DFIM(M)
      DO IJ=1,SIZE(F3,1)   
         SI(IJ) = F3(IJ,1,M)*SINTH(1)
         CI(IJ) = F3(IJ,1,M)*COSTH(1)
      END DO
      DIR_LOOP2: DO K=2,KL
         DO IJ=1,SIZE(F3,1)   
            SI(IJ) = SI(IJ) + F3(IJ,K,M)*SINTH(K)
            CI(IJ) = CI(IJ) + F3(IJ,K,M)*COSTH(K)
         END DO
      END DO DIR_LOOP2

      DO IJ=1,SIZE(F3,1)   
         SI(IJ) = FAK * SI(IJ) 
         CI(IJ) = FAK * CI(IJ)
         UST(IJ) = UST(IJ) + SI(IJ)
         VST(IJ) = VST(IJ) + CI(IJ)
      END DO
   END DO FREQ_LOOP2
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ADD CONTRIBUTION FROM TAIL.                                           !
!        ---------------------------                                           !

TAILFAC = FR(ML)**2 / (DF(ML)*(FR(ML)+0.5*DF(ML)))

UST(:) = UST(:) + TAILFAC*SI(:)
VST(:) = VST(:) + TAILFAC*CI(:)

END SUBROUTINE STOKES_DRIFT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TM1_TM2_PERIODS_B (F, EMEAN, TM1, TM2, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TM1_TM2_PERIODS_B - COMPUTES TM1 AND/OR TM2 PERIODS (VECTOR VESION).       !
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE TM1 AND TM2 PERIODS.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN )           :: F(:,:,:)    !! BLOCK OF SPECTRA.
REAL,    INTENT(IN )           :: EMEAN(:)    !! TOTAL ENERGY [M*M].
REAL,    INTENT(OUT), OPTIONAL :: TM1(:)      !! TM1 PERIOD [S].
REAL,    INTENT(OUT), OPTIONAL :: TM2(:)      !! TM2 PERIOD [S].
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, M
REAL    :: TEMP(SIZE(F,1),SIZE(F,3))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F, DIM=2, MASK=MASK)
ELSE
   TEMP = SUM(F, DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TM1 PERIOD.                                                           !
!        -----------                                                           !

IF (PRESENT(TM1)) THEN
   TM1(:) = MP1_TAIL * TEMP(:,ML)     !! TAIL
   DO M = 1,ML
      DO IJ = 1,SIZE(TEMP,1)
         TM1(IJ) = TM1(IJ) + TEMP(IJ,M)*DFIM_FR(M)
      END DO
   END DO

   WHERE (EMEAN.GT.EMIN)                          !! NORMALIZE WITH ENERGY.
      TM1 = EMEAN/TM1
   ELSEWHERE
      TM1 = 1.
   END WHERE
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. TM2 PERIOD.                                                           !
!        -----------                                                           !

IF (PRESENT(TM2)) THEN
   TM2(:) = MP2_TAIL * TEMP(:,ML)     !! TAIL
   DO M = 1,ML
      DO IJ = 1,SIZE(TEMP,1)
         TM2(IJ) = TM2(IJ) + TEMP(IJ,M)*DFIM_FR2(M)
      END DO
   END DO

   WHERE (EMEAN.GT.EMIN)                          !! NORMALIZE WITH ENERGY.
      TM2 = SQRT(EMEAN/TM2)
   ELSEWHERE
      TM2 = 1.
   END WHERE
END IF

END SUBROUTINE TM1_TM2_PERIODS_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TM1_TM2_PERIODS_1 (F, EMEAN, TM1, TM2, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TM1_TM2_PERIODS_1 - COMPUTES TM1 AND/OR TM2 PERIODS (SCALAR VESION).       !
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE TM1 AND TM2 PERIODS.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN )           :: F(:,:)    !! BLOCK OF SPECTRA.
REAL,    INTENT(IN )           :: EMEAN     !! TOTAL ENERGY [M*M].
REAL,    INTENT(OUT), OPTIONAL :: TM1       !! TM1 PERIOD [S].
REAL,    INTENT(OUT), OPTIONAL :: TM2       !! TM2 PERIOD [S].
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: M
REAL    :: TEMP(SIZE(F,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F, DIM=1, MASK=MASK)
ELSE
   TEMP = SUM(F, DIM=1)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TM1 PERIOD.                                                           !
!        -----------                                                           !

IF (PRESENT(TM1)) THEN
   TM1 = MP1_TAIL * TEMP(ML)     !! TAIL
   DO M = 1,ML
      TM1 = TM1 + TEMP(M)*DFIM_FR(M)
   END DO

   IF (EMEAN.GT.EMIN) THEN                    !! NORMALIZE WITH TOTAL ENERGY.
      TM1 = EMEAN/TM1
   ELSE
      TM1 = 1.
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. TM2 PERIOD.                                                           !
!        -----------                                                           !

IF (PRESENT(TM2)) THEN
   TM2 = MP2_TAIL * TEMP(ML)     !! TAIL
   DO M = 1,ML
      TM2 = TM2 + TEMP(M)*DFIM_FR2(M)
   END DO

   IF (EMEAN.GT.EMIN) THEN                     !! NORMALIZE WITH TOTAL ENERGY.
      TM2 = SQRT(EMEAN/TM2)
   ELSE
      TM2 = 1.
   END IF
END IF

END SUBROUTINE TM1_TM2_PERIODS_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TOTAL_ENERGY_B (F3, EMEAN, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TOTAL_ENERGY_B - COMPUTES TOTAL ENERGY (VECTOR VERSION).                   !
!                                                                              !
!     S.D. HASSELMANN.                                                         !
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER                               !
!     H. GUENTHER     GKSS   DECEMBER 2001  FT90                               !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGRATION OVER DIRECTION AND FREQUENCY. A TAIL CORRECTION IS ADDED.  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)            :: F3(:,:,:)    !! BLOCK OF SPECTRA.
REAL,    INTENT(OUT)           :: EMEAN(:)     !! TOTAL ENERGY.
LOGICAL, INTENT(IN), OPTIONAL  :: MASK(:,:,:)  !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: IJ, M
REAL     :: TEMP(SIZE(F3,1),SIZE(F3,3))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTION (WITHOUT DELTH).                             !
!        -----------------------------------------                             !

IF (PRESENT(MASK)) THEN
   TEMP(:,:) = SUM(F3(:,:,:), DIM=2, MASK=MASK(:,:,:))
ELSE
   TEMP(:,:) = SUM(F3(:,:,:), DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !

EMEAN(:) = MO_TAIL*TEMP(:,ML)    !! TAIL ENERGY

DO M = 1,ML
   DO IJ = 1,SIZE(F3,1)
      EMEAN(IJ) = EMEAN(IJ) + TEMP(IJ,M)*DFIM(M)
   END DO
END DO

EMEAN(:) = MAX(EMEAN(:), EMIN)

END SUBROUTINE TOTAL_ENERGY_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TOTAL_ENERGY_1 (F3, EMEAN, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TOTAL_ENERGY_1 - COMPUTES TOTAL ENERGY (SCALAR VERSION).                   !
!                                                                              !
!     S.D. HASSELMANN.                                                         !
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER                               !
!     H. GUENTHER     GKSS   DECEMBER 2001  FT90                               !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGRATION OVER DIRECTION AND FREQUENCY. A TAIL CORRECTION IS ADDED.  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)           :: F3(:,:)    !! SPECTRUM.
REAL,    INTENT(OUT)          :: EMEAN      !! TOTAL ENERGY.
LOGICAL, INTENT(IN), OPTIONAL :: MASK(:,:)  !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: M
REAL    :: TEMP(SIZE(F3,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTION (WITHOUT DELTH).                             !
!        -----------------------------------------                             !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F3 ,DIM=1,MASK=MASK)
ELSE
   TEMP = SUM(F3, DIM=1)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !

EMEAN = MO_TAIL*TEMP(ML)    !! TAIL ENERGY

DO M = 1,ML
   EMEAN = EMEAN + TEMP(M)*DFIM(M)
END DO

EMEAN = MAX(EMEAN, EMIN)

END SUBROUTINE TOTAL_ENERGY_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE COS2_SPR_1 (TH, THES, ST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     COS2_SPR - ROUTINE TO COMPUTE SPREADING FACTOR (SCALAR VERSION).         !
!                                                                              !
!     SUSANNE HASSELMANN  JULY 1986.                                           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF COS**2 SPREADING FUNCTION.                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN) :: TH(:)       !! DIRECTIONS.
REAL,    INTENT(IN) :: THES        !! MEAN WAVE DIRECTION.
REAL,    INTENT(OUT) :: ST(:)      !! SPREADING FUNCTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: ZDP=2./PI

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COSINE SQUARE SPREAD.                                                         !
!     ------------------------                                                         !

ST(:) = MAX(0. ,COS(TH(:)-THES))
ST = ZDP*ST**2
WHERE (ST.LT.0.1E-08) ST = 0.

END SUBROUTINE COS2_SPR_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE COS2_SPR_B (TH, THES, ST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     COS2_SPR - ROUTINE TO COMPUTE SPREADING FACTOR (VECTOR VERSION).         !
!                                                                              !
!     SUSANNE HASSELMANN  JULY 1986.                                           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF COS**2 SPREADING FUNCTION.                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN) :: TH(:)       !! DIRECTIONS.
REAL,    INTENT(IN) :: THES(:)     !! MEAN WAVE DIRECTIONS.
REAL,    INTENT(OUT) :: ST(:,:)    !! SPREADING FUNCTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: ZDP=2./PI
INTEGER         :: K

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COSINE SQUARE SPREAD.                                                         !
!     ------------------------                                                         !

DO K = 1,SIZE(TH)
   ST(:,K) = MAX(0. ,COS(TH(K)-THES(:)))
END DO
ST = ZDP*ST**2
WHERE (ST.LT.0.1E-08) ST = 0.

END SUBROUTINE COS2_SPR_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WM1_WM2_WAVENUMBER_B (F, EMEAN, WM1, WM2, IN, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   WM1_WM2_WAVENUMBER_B - COMPUTES WM1 AND/OR WM2 WAVENUMBERS (VECTOR VESION).!
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE  WM1 AND/OR WM2 WAVENUMBERS                                    !
!          WM1 IS SQRT(1/K)*F INTGRATION                                       !
!          WM2 IS SQRT(K)*F INTGRATION                                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN )           :: F(:,:,:)    !! BLOCK OF SPECTRA.
REAL,    INTENT(IN )           :: EMEAN(:)    !! TOTAL ENERGY [M*M].
REAL,    INTENT(OUT), OPTIONAL :: WM1(:)      !! WM1 WAVENUMBER [M].
REAL,    INTENT(OUT), OPTIONAL :: WM2(:)      !! WM2 WAVENUMBER [M].
INTEGER, INTENT(IN),  OPTIONAL :: IN   (:)    !! DEPTH TABLE INDEX.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, M
REAL    :: DEL2
REAL    :: TEMP(SIZE(F,1),SIZE(F,3)), TEMP2(SIZE(F,1),SIZE(F,3))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F, DIM=2, MASK=MASK)
ELSE
   TEMP = SUM(F, DIM=2)
END IF

IF (PRESENT(IN)) THEN
   TEMP2(:,:) = SQRT(TFAK(IN(:),:))
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM1.                                         !
!        -----------------------------                                         !

IF (PRESENT(WM1)) THEN
   DEL2 = SQRT(G)/ZPI
   DO IJ = 1,SIZE(F,1)
      WM1(IJ) = MM1_TAIL*DEL2*TEMP(IJ,ML)   !! TAIL.
   END DO

   IF (PRESENT(IN)) THEN
      DO M = 1,ML
         DO IJ = 1,SIZE(F,1)
            WM1(IJ) = WM1(IJ) + TEMP(IJ,M) / TEMP2(IJ,M) * DFIM(M)
         END DO
      END DO
   ELSE
      DO M = 1,ML
         DO IJ = 1,SIZE(TEMP,1)
            WM1(IJ) = WM1(IJ) + TEMP(IJ,M)*DFIMOFR(M)*DEL2
         END DO
      END DO
   END IF

   WHERE (EMEAN.GT.EMIN)
      WM1 = (EMEAN/WM1)**2                        !! NORMALIZE WITH ENERGY.
   ELSEWHERE
      WM1 = 1.
   END WHERE
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM2.                                         !
!        -----------------------------                                         !

IF (PRESENT(WM2)) THEN
   DEL2 = ZPI/SQRT(G)
   DO IJ = 1,SIZE(TEMP,1)
      WM2(IJ) = MP1_TAIL*DEL2*TEMP(IJ,ML)   !! ADD TAIL.
   END DO

   IF (PRESENT(IN)) THEN
      DO M = 1,ML
         DO IJ = 1,SIZE(TEMP,1)
            WM2(IJ) = WM2(IJ) + TEMP(IJ,M) * TEMP2(IJ,M) * DFIM(M)
         END DO
      END DO
   ELSE
      DO M = 1,ML
         DO IJ = 1,SIZE(TEMP,1)
            WM2(IJ) = WM2(IJ) + TEMP(IJ,M)*DFIM_FR(M)*DEL2
         END DO
      END DO
   END IF

   WHERE (EMEAN.GT.EMIN)
      WM2 = (WM2/EMEAN)**2                        !! NORMALIZE WITH ENERGY.
   ELSEWHERE
      WM2 = 1.
   END WHERE
END IF

END SUBROUTINE WM1_WM2_WAVENUMBER_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WM1_WM2_WAVENUMBER_1 (F, EMEAN, WM1, WM2, IN, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   WM1_WM2_WAVENUMBER_1 - COMPUTES WM1 AND/OR WM2 WAVENUMBER (SCALAR VESION). !
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE  WM1 AND/OR WM2 WAVENUMBERS                                    !
!          WM1 IS SQRT(1/K)*F INTGRATION                                       !
!          WM2 IS SQRT(K)*F INTGRATION                                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN )           :: F(:,:)    !! SPECTRUM.
REAL,    INTENT(IN )           :: EMEAN     !! TOTAL ENERGY [M*M].
REAL,    INTENT(OUT), OPTIONAL :: WM1       !! WM1 WAVENUMBER [M].
REAL,    INTENT(OUT), OPTIONAL :: WM2       !! WM2 WAVENUMBER [M].
INTEGER, INTENT(IN),  OPTIONAL :: IN        !! DEPTH TABLE INDEX.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: M
REAL    :: DEL2
REAL    :: TEMP(SIZE(F,2)), TEMP2(SIZE(F,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F, DIM=1, MASK=MASK)
ELSE
   TEMP = SUM(F, DIM=1)
END IF

IF (PRESENT(IN)) THEN
   TEMP2(:) = SQRT(TFAK(IN,:))
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM1.                                         !
!        -----------------------------                                         !

IF (PRESENT(WM1)) THEN
   DEL2 = SQRT(G)/ZPI
   WM1 = MM1_TAIL*DEL2*TEMP(ML)   !! ADD TAIL.

   IF (PRESENT(IN)) THEN
      DO M = 1,ML
         WM1 = WM1 + TEMP(M) / TEMP2(M) * DFIM(M)
      END DO
   ELSE
      DO M = 1,ML
         WM1 = WM1 + TEMP(M)*DFIMOFR(M)*DEL2
      END DO
   END IF

   IF (EMEAN.GT.EMIN) THEN
      WM1 = (EMEAN/WM1)**2                            !! NORMALIZE WITH ENERGY.
   ELSE
      WM1 = 1.
   END IF  
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM2.                                         !
!        -----------------------------                                         !

IF (PRESENT(WM2)) THEN
   DEL2 = ZPI/SQRT(G)
   WM2 = MP1_TAIL*DEL2*TEMP(ML)   !! ADD TAIL.

   IF (PRESENT(IN)) THEN
      DO M = 1,ML
         WM2 = WM2 + TEMP(M) * TEMP2(M) * DFIM(M)
      END DO
   ELSE
      DO M = 1,ML
         WM2 = WM2 + TEMP(M)*DFIM_FR(M)*DEL2
      END DO
   END IF

   IF (EMEAN.GT.EMIN) THEN
      WM2 = (WM2/EMEAN)**2                            !! NORMALIZE WITH ENERGY.
   ELSE
      WM2 = 1.
   END IF  
END IF

END SUBROUTINE WM1_WM2_WAVENUMBER_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE KURTOSIS (F3, DEPTH, C4, BF2, QP, FP, THMAX, HMAX, TMAX)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    KURTOSIS   DETERMINES KURTOSIS, BENJAMIN-FEIR INDEX
!               GODA'S PEAKEDNESS PARAMETER AND NORMALIZED
!               MAXIMUM WAVE HEIGHT.
!                                                                              !
!     PETER JANSSEN       JULY 2007.
!                                                                              !
!     PURPOSE.
!     --------
!                                                                              !
!           DETERMINATION OF KURTOSIS, B-F INDEX , GODA QP AND
!                              HMAX AND TMAX
!                                                                              !
!     METHOD.
!     -------
!            USED NUMERICAL SIMULATIONS OF THE NONLINEAR SCHROEDINGER
!            EQUATION TO OBTAIN THE DEPENDENCE OF KURTOSIS ON THE
!            SQUARE OF THE BENJAMIN-FEIR INDEX, BF2, AND THE WIDTH OF
!            THE DIRECTIONAL DISTRIBUTION, SIG_TH. FITTING THE NUMERICAL
!            RESULTS ONE FINDS (ONORATO, MORI AND JANSSEN, 2007)
!
!              C4 = A_MORI/SIG_TH * BF2 * PI/(3*SQRT(3))
!
!            WHERE A_MORI=0.031. NOTE THAT FOR NARROW SPECTRA, I.E.
!            SIG_TH => A_MORI, THE ONE-DIMENSIONAL RESULT, PREVIOUSLY
!            IMPLEMENTED, IS OBTAINED.
!
!            IN ADDITION, IT IS WELL-KNOWN THAT WHEN A WAVE GROUP
!            APPROACHES SHALLOW WATER THE NONLINEAR FOCUSSING
!            EFFECT DISAPPEARS FOR K*D=1.363 AND WHEN K*D BECOMES
!            LESS THAN THIS CRITICAL VALUE EVEN DEFOCUSSING OCCURS.
!            AS A CONSEQUENCE, FOR K*D>1.363 C4 > 0 WHILE IN THE
!            OPPOSITE CASE C4 < 0. HERE, FOLLOWING THE WORK
!            OF JANSSEN & ONORATO (2007) WE HAVE PARAMETRIZED THIS SHALLOW
!            WATER EFFECT BY MEANS OF THE FUNCTION TRANSF2.
!                                                                              !
!     EXTERNALS.
!     ----------
!             AKI,TRANSF2,H_MAX,PEAK_FREQ
!                                                                              !
!     REFERENCES:
!     ---------
!                                                                              !
!     NONLINEAR FOUR WAVE INTERACTIONS AND FREAK WAVES
!     PETER A.E.M. JANSSEN, JPO, 863-884, APRIL 2003

!     ON KURTOSIS AND OCCURRENCE PROBABILITY OF FREAK WAVES
!     NOBUHITO MORI AND PETER A.E.M. JANSSEN, JPO, 1471-1483, JULY 2006
!                                                                              !
!     THE INTERMEDIATE WATER DEPTH LIMIT OF THE ZAKHAROV EQUATION AND
!     CONSEQUENCES FOR WAVE PREDICTION
!     PETER A.E.M. JANSSEN AND MIGUEL ONORATO, ACCEPTED FOR PUBLICATION
!     IN JPO, 2007

!     EFFECTS OF DIRECTIONALITY ON THE NONLINEAR FOCUSSING IN RANDOM SEAS.
!     MIGUEL ONORATO, NOBUHITO MORI AND PETER A.E.M. JANSSEN
!     IN PREPARATION, 2007
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !
!                                                                              !
REAL,    INTENT(IN)            :: F3(:,:,:)  !! 2-DIMENSIONAL SPECTRA.
REAL,    INTENT(IN)            :: DEPTH(:)   !! WATER DEPTH.
REAL,    INTENT(OUT)           :: C4(:)      !! KURTOSIS.
REAL,    INTENT(OUT)           :: BF2(:)     !! BENJAMIN-FEIR INDEX.
REAL,    INTENT(OUT)           :: QP(:)      !! GODA'S QUALITY FACTOR.
REAL,    INTENT(OUT)           :: FP(:)      !! PEAK FREQUENCY.
REAL,    INTENT(OUT)           :: THMAX(:)   !! PEAK DIRECTION.
REAL,    INTENT(OUT)           :: HMAX(:)    !! MAXIMUM WAVE HEIGHT.
REAL,    INTENT(OUT)           :: TMAX(:)    !! MAXIMUM WAVE PERIOD.

! REAL, EXTERNAL :: AKI, TRANSF2

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: DURATION = 10800.
REAL, PARAMETER :: ZA_MORI = 0.062
REAL, PARAMETER :: ZALPHA = 0.4

INTEGER :: IJ, M, K

INTEGER ::  NW(SIZE(F3,1))    !! NUMBER OF WAVES.

REAL :: ZSUM3, ZJ , ZEPS, ZNU, HS, SIG_OM1

REAL :: F1D(SIZE(F3,1),SIZE(F3,3))    !! FREQUENCY SPECTRA
REAL :: XMAX(SIZE(F3,1))              !! MAX FREQ SPECTRA
REAL :: SIG_OM(SIZE(F3,1))            !! RELATIVE WIDTH IN FREQUENCY

REAL :: F1A(SIZE(F3,1),SIZE(F3,2))    !! DIRECTIONAL SPECTRA
REAL :: YMAX(SIZE(F3,1))              !! MAX DIRECTIONAL SPECTRUM
REAL :: SIG_TH(SIZE(F3,1))            !! RELATIVE WIDTH IN DIRECTION


REAL, DIMENSION(SIZE(F3,1)) :: ETA2,SUM0,SUM1, SUM2, SUM4, XKP, EPS
REAL, DIMENSION(SIZE(F3,3)) :: FAC4

ZSUM3  = PI/(3.*SQRT(3.))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DETERMINE ONE-D FREQUENCY AND ANGULAR SPECTRUM.
!        -----------------------------------------------

DO M = 1, ML
   FAC4(M) = 2.*DF_FR(M)
END DO

F1D = SUM(F3, DIM=2)*DELTH  !! 1-D Frequency Spectra

DO K = 1, KL
   M = 1
   F1A(:,K) = F3(:,K,M)*DF(M)
   DO M = 2,ML
      F1A(:,K) = F1A(:,K)+F3(:,K,M)*DF(M)
   END DO
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. DETERMINE PEAK PARAMETERS FP,SIG_OM AND SIG_TH.                       !
!        -----------------------------------------------                       !

CALL PEAK_FREQ (F1D, F1A, XMAX, FP, SIG_OM, YMAX, THMAX, SIG_TH)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. DETERMINE BF2, QP AND C4.                                             !
!        -------------------------                                             !


SUM0(:) = F1D(:,1)*DF(1)
SUM1(:) = F1D(:,1)*DF_FR(1)
SUM2(:) = F1D(:,1)*DF_FR2(1)

DO M = 2, ML
   DO IJ = 1, SIZE(F3,1)
      SUM0(IJ) = SUM0(IJ) + F1D(IJ,M)*DF(M)       !! ETOT
      SUM1(IJ) = SUM1(IJ) + F1D(IJ,M)*DF_FR(M)    !! F*E   =1/TM1
      SUM2(IJ) = SUM2(IJ) + F1D(IJ,M)*DF_FR2(M)   !! F*F*E =1/TM2
   END DO
END DO

ETA2(:)  = 0.
SUM4(:)  = 0.

DO M = 1, ML
   DO IJ = 1, SIZE(F3,1)
      IF (F1D(IJ,M) .LE. ZALPHA*XMAX(IJ)) CYCLE
      ETA2(IJ) = ETA2(IJ)+F1D(IJ,M)*DF(M)
      SUM4(IJ) = SUM4(IJ)+F1D(IJ,M)**2*FAC4(M)
   END DO
END DO

DO IJ = 1, SIZE(F3,1)
   IF (ETA2(IJ).GT.0.) THEN
      QP(IJ)     = SUM4(IJ)/ETA2(IJ)**2
      SIG_OM1    = 1./(QP(IJ)*ZPISQRT)
      SIG_OM(IJ) = MIN(SIG_OM(IJ),SIG_OM1)
      XKP(IJ)    = AKI(ZPI*FP(IJ),DEPTH(IJ))
      EPS(IJ)    = XKP(IJ)*SQRT(SUM0(IJ))

      BF2(IJ) = 2.*EPS(IJ)**2/SIG_OM(IJ)**2*TRANSF2(XKP(IJ),DEPTH(IJ))
      BF2(IJ) = MAX(MIN(BF2(IJ),BF2MAX),BF2MIN)
   ELSE
      BF2(IJ)    = 0.
      EPS(IJ)    = 0.
      SIG_OM(IJ) = 0.
      QP(IJ)     = 0.
   ENDIF
END DO

DO IJ = 1, SIZE(F3,1)
   IF (ETA2(IJ).GT.0.) THEN
      SIG_TH(IJ) = MAX(SIG_TH(IJ),ZA_MORI)
      ZJ = ZA_MORI/SIG_TH(IJ)*ZSUM3
      C4(IJ) = ZJ*BF2(IJ)+8.*EPS(IJ)**2
      C4(IJ) = MAX(MIN(C4(IJ),C4MAX),C4MIN)
   ELSE
      C4(IJ)  = 0.
   ENDIF
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. DETERMINE HMAX AND TMAX.
!       (EXPECTED IN RECORD OF LENGTH GIVEN BY DURATION).
!       -------------------------------------------------

DO IJ = 1, SIZE(F3,1)
   IF (FP(IJ).EQ.0.) THEN
      NW(IJ) = NINT(0.1*DURATION)
   ELSE
      NW(IJ) = NINT(DURATION*FP(IJ))
   ENDIF
END DO

CALL H_MAX (C4, NW, HMAX)

DO IJ = 1, SIZE(F3,1)
   IF (SUM1(IJ).GT.0. .AND. HMAX(IJ).GT.0.) THEN
      ZNU = SUM0(IJ)*SUM2(IJ)/SUM1(IJ)**2 - 1.
      ZNU = SQRT(MAX(0.,ZNU))
      ZEPS = ZNU/(SQRT(2.)*HMAX(IJ))
      TMAX(IJ) = 1.+0.5*ZEPS**2+0.75*ZEPS**4
      TMAX(IJ) = TMAX(IJ)*SUM0(IJ)/SUM1(IJ)
   ELSE
      TMAX(IJ) = 0.
   ENDIF
END DO

DO IJ = 1, SIZE(F3,1)
   IF (SUM0(IJ).GT.0.) THEN
      HS = 4.*SQRT(SUM0(IJ))
      HMAX(IJ) = HMAX(IJ)*HS
   ELSE
      HMAX(IJ) = 0.
   ENDIF
END DO

END SUBROUTINE KURTOSIS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PEAK_FREQ (F1D, F1A, XMAX, FP, SIG_OM, YMAX, THMAX, SIG_TH)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     PEAK_FREQ   DETERMINES PARAMETERS AT PEAK OF SPECTRUM.                   !
!                                                                              !
!     PETER JANSSEN                                                            !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!              DETERMINATION OF PEAK PARAMETERS.                               !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!              NONE                                                            !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !
!              NONE                                                            !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !
!                                                                              !
REAL,    INTENT(IN)    :: F1D(:,:)   !! FREQUENCY SPECTRA
REAL,    INTENT(IN)    :: F1A(:,:)   !! DIRECTIONAL SPECTRA
REAL,    INTENT(OUT)   :: XMAX(:)    !! MAX FREQ SPECTRA
REAL,    INTENT(OUT)   :: FP(:)      !! PEAK FREQUENCY (HZ)
REAL,    INTENT(OUT)   :: SIG_OM(:)  !! RELATIVE WIDTH IN FREQUENCY
REAL,    INTENT(OUT)   :: YMAX(:)    !! MAX DIRECTIONAL SPECTRUM
REAL,    INTENT(OUT)   :: THMAX(:)   !! PEAK DIRECTION (RAD)
REAL,    INTENT(OUT)   :: SIG_TH(:)  !! RELATIVE WIDTH IN DIRECTION

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, K, KPLUS, KMIN

INTEGER ::  MMAX(SIZE(F1D,1))  !! FREQUENCY INDICES OF MAX FREQ SPECTRA.
INTEGER ::  KMAX(SIZE(F1A,1))  !! DIRECTION INDICES OF MAX DIRECTIONAL SPECTRA.

REAL :: A, B, C, XP1, XM1, F1P1, F1M1, SIG_TH2
REAL :: SUM1(SIZE(F1A,1))
REAL :: SUM2(SIZE(F1A,1))
REAL :: SIG_TH1(SIZE(F1A,1))

! ---------------------------------------------------------------------------- !
!
!     1. DETERMINE F_PEAK USING QUADRATIC FIT.
!        -------------------------------------

!     MAX OF 1-D SPECTRUM

MMAX(:) = MAXLOC(F1D(:,:), DIM=2)
FP(:) = FR(MMAX(:))
DO IJ = 1, SIZE(F1D,1)
   XMAX(IJ) = F1D(IJ,MMAX(IJ))
END DO

SIG_OM(:) = 1.0

DO IJ = 1, SIZE(F1D,1)
   IF (XMAX(IJ).LE.0. .OR. MMAX(IJ).GE.ML .OR. MMAX(IJ).LE.1) CYCLE
   XP1 = FR(MMAX(IJ)+1)-FP(IJ)
   XM1 = FR(MMAX(IJ)-1)-FP(IJ)
   A = XMAX(IJ)
   F1P1 = (F1D(IJ,MMAX(IJ)+1) - A)/XP1
   F1M1 = (F1D(IJ,MMAX(IJ)-1) - A)/XM1

   C = (F1M1-F1P1)/(XM1-XP1)

   IF (C.GE.0.) CYCLE
   B = (XM1*F1P1-XP1*F1M1)/(XM1-XP1)
   FP(IJ) = FP(IJ)-B/(2.*C)
   XMAX(IJ) = A-B**2/(4.*C)
   SIG_OM(IJ) = SQRT((B**2-4.*A*C)/(8.*C**2))/FP(IJ)
END DO

! ---------------------------------------------------------------------------- !

!     2. DETERMINE YMAX, THMAX AND SIG_TH USING QUADRATIC FIT.
!     --------------------------------------------------------

KMAX(:) = MAXLOC(F1A(:,:), DIM=2)
THMAX(:) = TH(KMAX(:))
DO IJ = 1, SIZE(F1A,1)
   YMAX(IJ) = F1A(IJ,KMAX(IJ))
END DO

SIG_TH1(:) = 1.0

DO IJ = 1, SIZE(F1A,1)
   IF (YMAX(IJ).LE.0.) CYCLE
   XP1 = DELTH
   XM1 = -DELTH
   KPLUS = KMAX(IJ)+1
   IF (KPLUS.GT.KL) KPLUS = KPLUS - KL
   KMIN = KMAX(IJ)-1
   IF (KMIN.LT.1) KMIN = KMIN + KL
   A  = YMAX(IJ)
   F1P1 = (F1A(IJ,KPLUS) - A)/XP1
   F1M1 = (F1A(IJ,KMIN)  - A)/XM1

   C = (F1M1-F1P1)/(XM1-XP1)

   IF (C.GE.0.) CYCLE
   B = (XM1*F1P1-XP1*F1M1)/(XM1-XP1)
   THMAX(IJ) = THMAX(IJ)-B/(2.*C)
   YMAX(IJ) = A-B**2/(4.*C)
   SIG_TH1(IJ) = SQRT((B**2-4.*A*C)/(8.*C**2))
END DO

SUM1(:) = 0.
SUM2(:) = 0.

DO K = 1, KL
   DO IJ = 1, SIZE(F1A,1)
      IF (COS(TH(K)-THMAX(IJ)).LT.0.5) CYCLE
      SUM1(IJ) = SUM1(IJ) + F1A(IJ,K)
      SUM2(IJ) = SUM2(IJ) + COS(TH(K)-THMAX(IJ))*F1A(IJ,K)
   END DO
END DO

DO IJ = 1, SIZE(F1A,1)
   IF (SUM1(IJ).GT.0.) THEN
       SIG_TH2 = SQRT(2.*(1.-SUM2(IJ)/SUM1(IJ)))
   ELSE
      SIG_TH2 = 1.
   ENDIF
   SIG_TH(IJ) = MIN(SIG_TH1(IJ),SIG_TH2)
ENDDO

END SUBROUTINE PEAK_FREQ

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE H_MAX (C4, NW, HMAX_N)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   H_MAX   DETERMINE EXPECTED MAXIMUM WAVE HEIGHT.
!                                                                              !
!     P.A.E.M. JANSSEN, JANUARY 2008
!                                                                              !
!     PURPOSE:
!     -------
!                                                                              !
!     DETERMINE EXPECTED MAXIMUM WAVE HEIGHT ACCORDING TO
!                 /
!         <H_MAX>=| dh h p(h)
!                 /
!     WHERE p(h) IS DERIVED IN MORI AND JANSSEN(2006)
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !
!                                                                              !
REAL,    INTENT(IN)   :: C4(:)     !! KURTOSIS(<ETA>^4/(3<ETA^2>^2)-1)
INTEGER, INTENT(IN)   :: NW(:)     !! NUMBER OF WAVES
REAL,    INTENT(OUT)  :: HMAX_N(:) !! MAXIMUM WAVE HEIGHT NORMALIZED WITH
!! SIGNIFICANT WAVE HEIGHT H_S.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: GAMMA_E = 0.57721566  ! EULER CONSTANT
REAL, PARAMETER :: CONST   = 0.5*(GAMMA_E**2+PI**2/6.)

INTEGER :: IJ
REAL    :: C4M, Z0, B_Z, ARG

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DETERMINE EXPECTED MAXIMUM WAVE HEIGHT.
!        ---------------------------------------

DO IJ = 1, SIZE(C4)
   C4M = MIN (C4MAX, C4(IJ))
   C4M = MAX (C4MIN, C4M)

   Z0 = 0.5*LOG(REAL(NW(IJ)))
   B_Z = 2.*Z0*(Z0-1.)
   ARG = MAX(1.+C4M*(B_Z-GAMMA_E*(1.-2.*Z0)-CONST), 0.1)
   HMAX_N(IJ) = SQRT(Z0 + 0.5*(GAMMA_E+LOG(ARG)))
END DO

END SUBROUTINE H_MAX

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

REAL FUNCTION TRANSF2 (XK, D)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    TRANSF2   DETERMINE NARROW BAND LIMIT BENJAMIN-FEIR INDEX FOR             !
!              THE FINITE DEPTH CASE.
!                                                                              !
!     PETER JANSSEN       JULY 2007.
!                                                                              !
!     PURPOSE.
!     --------
!                                                                              !
!           DETERMINATION OF THENARROW BAND LIMIT BENJAMIN-FEIR INDEX FOR      !
!           THE FINITE DEPTH CASE.                                             !
!
!           BF**2 = (2 S^2)/SIG_OM^2) . TRANSF2(XK,D)
!
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !
!                                                                              !
REAL,    INTENT(IN)    :: XK   !! WAVE NUMBER
REAL,    INTENT(IN)    :: D    !! DEPTH

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL ::  X, T_0, OM, C_0, C_S, V_G, D2OM, XNL_1, XNL_2, T_NL

! ---------------------------------------------------------------------------- !
!
!     1. DETERMINE TRANSFER FUNCTION.
!        ----------------------------

X   = MIN(XK*D, DKMAX)
X   = MAX(X, 0.5)
T_0 = TANH(X)
OM  = SQRT(G*XK*T_0)
C_0 = OM/XK
C_S = SQRT(G*D)
V_G = 0.5*C_0*(1.+2.*X/SINH(2.*X))
D2OM = V_G**2-C_S**2*(1.-T_0**2)*(1.-X*T_0)

XNL_1 = (9.*T_0**4-10.*T_0**2+9.)/(8.*T_0**3)
XNL_2 = ((2.*V_G-0.5*C_0)**2/(G*D-V_G**2)+1.)/X

T_NL = XK**3*(XNL_1-XNL_2)
TRANSF2 = (V_G/C_0)**2*G*T_NL/XK**4/D2OM

END FUNCTION TRANSF2

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_INTERFACE_MODULE

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
&                             MO_TAIL, MM1_TAIL, MP1_TAIL, MP2_TAIL

USE WAM_TABLES_MODULE,  ONLY: TFAK, TFAC_ST

USE WAM_GENERAL_MODULE, ONLY: G, PI, ZPI, DKMAX, ZPISQRT, GAMMA_E,             &
&                             BF2MAX, BF2MIN, C4MAX, C4MIN

USE WAM_OUTPUT_SET_UP_MODULE, ONLY: WMDUR, WMDX, WMDY                          !! WAM-MAX

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

PRIVATE

REAL      :: EMIN = 1.0E-12    !! REPLACES THE INTRINSIC TINY

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

INTERFACE WAMAX                            !! EXPECTED MAXIMUM WAVE HEIGHT.    !! WAM-MAX
   MODULE  PROCEDURE WAMAX                                                     !! WAM-MAX
END INTERFACE                                                                  !! WAM-MAX
PUBLIC WAMAX                                                                   !! WAM-MAX

INTERFACE TRANSF                           !! NARROW BAND LIMIT BENJAMIN-FEIR
   MODULE  PROCEDURE TRANSF                !! INDEX FOR THE FINITE DEPTH CASE.
END INTERFACE
PUBLIC TRANSF

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
!
!==>   deactivate tail correction since it's not always possible !
!
CONST1 = 0.*LOG(FC/FS)
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
!     1. COSINE SQUARE SPREAD.                                                 !
!     ------------------------                                                 !

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

REAL,    INTENT(IN)  :: TH(:)      !! DIRECTIONS.
REAL,    INTENT(IN)  :: THES(:)    !! MEAN WAVE DIRECTIONS.
REAL,    INTENT(OUT) :: ST(:,:)    !! SPREADING FUNCTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: ZDP=2./PI
INTEGER         :: K

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COSINE SQUARE SPREAD.                                                 !
!     ------------------------                                                 !

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
!    KURTOSIS   DETERMINES SKEWNESS, KURTOSIS, BENJAMIN-FEIR INDEX
!               GODA'S PEAKEDNESS PARAMETER AND NORMALIZED
!               MAXIMUM WAVE HEIGHT.
!                                                                              !
!     PETER JANSSEN       JULY 2007.
!     PETER JANSSEN       APRIL 2014. INCLUDES SKEWNESS EFFECTS, WHILE
!                         NUMBER OF WAVES IS DETERMINED BY EWING (1973)
!                                                                              !
!     PURPOSE.
!     --------
!                                                                              !
!           DETERMINATION OF SKEWNESS, KURTOSIS, B-F INDEX , GODA QP AND
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
!
!            APRIL 2014:
!            ----------
!
!            1) FORMULATE DYNAMIC KURTOSIS C4 IN TERMS OF THE DIMENSIONLESS
!               NUMBER R = 0.5*(SIG_TH/SIG_OM)**2 (SEE TECH MEMO 588).
!               RESULTS IN GOOD AGREEMENT WITH  MORI FORMULATION FOR
!                  C4 = XJ*BFI**2
!               WITH
!                  XJ = C4_CONST/SQRT(1+7.*R),
!               WHERE
!                  C4_CONST = PI/(3.*SQRT(3.))
!            2) BOUND WAVES HAVE ALSO FINITE SKEWNESS WHICH AFFECTS ENVELOPE
!               WAVE HEIGHT DISTRIBUTION (SEE MEMO ON UPDATES)
!            3) FOR MAXIMUM WAVE HEIGHT THE NUMBER NW OF INDEPENDENT EVENTS
!               IS REQUIRED. THUS FAR NUMBER OF INDEPENDENT EVENTS WAS OBTAINED
!               BY SAMPLING WITH THE PEAK FREQUENCY. THIS IS NOT CORRECT. IT
!               IS MORE APPROPRIATE TO USE AS ESTIMATE OF EVENTS THE NUMBER OF
!               WAVE GROUPS IN THE TIME SERIES. THIS IS OBTAINED USING JOHN
!               EWINGS (1973) RESULT ON THE LENGTH OF A WAVE GROUP REFERRED TO
!               SIGNIFICANT WAVE HEIGHT.
!     ON THE EXTENSION OF THE FREAK WAVE WARNING SYSTEM AND ITS VERIFICATION.
!     PETER A.E.M JANSSEN AND J.-R. BIDLOT, ECMWF TECH MEMO 588.

!     FURTHER UPDATES TO THE FREAK WAVE WARNING SYSTEM. PETER A.E.M.
!     JANSSEN AND J.-R. BIDLOT, TO BE PUBLISHED AS ECMWF TECH MEMO
!
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
! REAL,    INTENT(OUT)           :: C3(:)      !! SKEWNESS.
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

REAL :: ZKAPPA30, ZKAPPA03, ZKAPPA12, ZKAPPA21, ZKAPPA40, ZKAPPA04, ZKAPPA22, ZKAPPA4
REAL :: ZR, ZARG, ZC4_CONST
REAL :: ZK, ZF_WG

REAL :: C3(SIZE(F3,1))      !! SKEWNESS.

REAL :: F1D(SIZE(F3,1),SIZE(F3,3))    !! FREQUENCY SPECTRA
REAL :: XMAX(SIZE(F3,1))              !! MAX FREQ SPECTRA
REAL :: SIG_OM(SIZE(F3,1))            !! RELATIVE WIDTH IN FREQUENCY

REAL :: F1A(SIZE(F3,1),SIZE(F3,2))    !! DIRECTIONAL SPECTRA
REAL :: YMAX(SIZE(F3,1))              !! MAX DIRECTIONAL SPECTRUM
REAL :: SIG_TH(SIZE(F3,1))            !! RELATIVE WIDTH IN DIRECTION


REAL, DIMENSION(SIZE(F3,1)) :: ETA2,SUM0,SUM1, SUM2, SUM4, XKP, EPS
REAL, DIMENSION(SIZE(F3,3)) :: FAC4

ZSUM3  = PI/(3.*SQRT(3.))
ZC4_CONST = PI/(3.*SQRT(3.))

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

!! DO IJ = 1, SIZE(F3,1)
!!    IF (ETA2(IJ).GT.0.) THEN
!!       SIG_TH(IJ) = MAX(SIG_TH(IJ),ZA_MORI)
!!       ZJ = ZA_MORI/SIG_TH(IJ)*ZSUM3
!!       C4(IJ) = ZJ*BF2(IJ)+8.*EPS(IJ)**2
!!       C4(IJ) = MAX(MIN(C4(IJ),C4MAX),C4MIN)
!!    ELSE
!!       C4(IJ)  = 0.
!!    ENDIF
!! END DO



DO IJ = 1, SIZE(F3,1)
   IF (ETA2(IJ).GT.0.) THEN
!
!           SKEWNESS
!
      ZKAPPA30 = 3.*EPS(IJ)
      ZKAPPA03 = 0.
      ZKAPPA12 = ZKAPPA30/3.
      ZKAPPA21 = 0.
!
!           BOUND-WAVE PART OF KURTOSIS
!

      ZKAPPA40 = 18.*EPS(IJ)**2
      ZKAPPA04 = 0.
      ZKAPPA22 = ZKAPPA40/6.
      ZKAPPA4  = ZKAPPA40+ZKAPPA04+2.*ZKAPPA22
!
!           DYNAMIC PART OF KURTOSIS
!
      ZR = 0.5*(SIG_TH(IJ)/SIG_OM(IJ))**2
      ZJ = ZC4_CONST/SQRT(1+7.*ZR)

      C4(IJ) = ZJ*BF2(IJ)+ZKAPPA4/8.
      C4(IJ) = MAX(MIN(C4(IJ),C4MAX),C4MIN)

      ZARG = 5.*(ZKAPPA30**2+ZKAPPA03**2)+9*(ZKAPPA21**2+                    &
&            ZKAPPA12**2)+6.*(ZKAPPA30*ZKAPPA12+ZKAPPA03*ZKAPPA21)

      C3(IJ) = SQRT(ZARG/72.)
   ELSE
      C4(IJ) = 0.
      C3(IJ) = 0.
   ENDIF
ENDDO


! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. DETERMINE HMAX AND TMAX.
!       (EXPECTED IN RECORD OF LENGTH GIVEN BY DURATION).
!       -------------------------------------------------
!
!     WAVE GROUP FREQUENCY F_WG FOLLOWS FROM EWING (1973) AND CRAMER AND
!     LEADBETTER (1967), WHERE K = 2  CORRESPONDS TO SIGNIFICANT HEIGHT LEVEL.

ZK = 2.

DO IJ = 1, SIZE(F3,1)
   IF (ETA2(IJ).GT.0. .AND. FP(IJ).GT.0.) THEN
      ZF_WG  = ZK*SIG_OM(IJ)*FP(IJ)*SQRT(ZPI)
      NW(IJ) = NINT(DURATION*ZF_WG)
   ELSE
      NW(IJ) = 0
   ENDIF
ENDDO

CALL H_MAX (C3, C4, NW, HMAX)

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

SUBROUTINE H_MAX (C3, C4, NW, HMAX_N)

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
REAL,    INTENT(IN)   :: C3(:)     !! SKEWNESS
REAL,    INTENT(IN)   :: C4(:)     !! KURTOSIS(<ETA>^4/(3<ETA^2>^2)-1)
INTEGER, INTENT(IN)   :: NW(:)     !! NUMBER OF WAVES
REAL,    INTENT(OUT)  :: HMAX_N(:) !! MAXIMUM WAVE HEIGHT NORMALIZED WITH
                                   !! SIGNIFICANT WAVE HEIGHT H_S.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: GAMMA_E = 0.57721566  !! EULER CONSTANT
REAL, PARAMETER :: ZETA_3  = 1.20206     !! RIEMANN ZETA FUNCTION

!   G1, G2, G3 ARE FIRST, SECOND AND THIRD DERIVATIVE OF GAMMA FUNCTION AT Z=1.
REAL, PARAMETER :: G1 = -GAMMA_E
REAL, PARAMETER :: G2 = (GAMMA_E**2+PI**2/6.)
REAL, PARAMETER :: G3 = -2.*ZETA_3-GAMMA_E*PI**2/2.-GAMMA_E**3

INTEGER :: IJ
REAL    :: C3M, C4M, Z0, ALPHA, C_Z, D_Z, E_Z, BETA, ARG, Z

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DETERMINE EXPECTED MAXIMUM WAVE HEIGHT.
!        ---------------------------------------

DO IJ = 1, SIZE(C3)
   IF (NW(IJ).GT.0) THEN
      C3M = C3(IJ)
      C4M = MIN (C4MAX, C4(IJ))
      C4M = MAX (C4MIN, C4M)

      Z0 = 0.5*LOG(REAL(NW(IJ)))
      ALPHA = 2.*Z0*(Z0-1.)+(1.-2.*Z0)*G1+0.5*G2

      C_Z = Z0*(4.*Z0**2-12.*Z0+6.)-(6.*Z0**2-12.*Z0+3.)*G1
      D_Z = 3.*(Z0-1.)*G2
      E_Z = -0.5*G3
      BETA = C_Z+D_Z+E_Z

      ARG = 1.+C4M*ALPHA+C3M**2*BETA
      ARG = MAX(ARG,0.1)
      Z = Z0+0.5*(GAMMA_E+LOG(ARG))

      HMAX_N(IJ) = SQRT(Z)
   ELSE
      HMAX_N(IJ) = 0.
   ENDIF
END DO

END SUBROUTINE H_MAX

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WAMAX (F, DEPTH, THMAX, CMAX_F, HMAX_N, CMAX_ST, HMAX_ST, MASK)               !! WAM-MAX
                                                                                         !! WAM-MAX
 ! ---------------------------------------------------------------------------- !        !! WAM-MAX
 !                                                                              !        !! WAM-MAX
 !   WAMAX   DETERMINE EXPECTED MAXIMUM CREST AND CREST-TO-TROUGH WAVE HEIGHTS. !        !! WAM-MAX
 !                                                                              !        !! WAM-MAX
 !     FRANCESCO BARBARIOL  CNR-ISMAR   adapted from WW3 5.16 (09/2018)         !        !! WAM-MAX
 !     PAOLO PEZZUTTO       CNR-ISMAR   min of autocov golden search (03/2019)  !        !! WAM-MAX
 !                                                                              !        !! WAM-MAX
 !     PURPOSE:                                                                 !        !! WAM-MAX
 !     -------                                                                  !        !! WAM-MAX
 !                                                                              !        !! WAM-MAX
 !           DETERMINE EXPECTED MAXIMUM CREST AND CREST-TO-TROUGH WAVE HEIGHTS  !        !! WAM-MAX
 !           ACCORDING TO TWO DIFFERENT STATISTICAL APPROACHES:                 !        !! WAM-MAX
 !           1. TIME EXTREMES:                                                  !        !! WAM-MAX
 !              1a. CREST H., FORRISTALL (2000), NONLINEAR 2ND ORDER            !        !! WAM-MAX
 !              1b. WAVE H., NAESS (1985), LINEAR (ARBITRARY BANDWIDTH)         !        !! WAM-MAX
 !           2. SPACE-TIME EXTREMES:                                            !        !! WAM-MAX
 !              2a. CREST H., BENETAZZO ET AL. (2015), NONLINEAR 2ND ORDER      !        !! WAM-MAX
 !              2b. WAVE H., FEDELE (2012), BOCCOTTI (2000), LINEAR             !        !! WAM-MAX
 !                                                                              !        !! WAM-MAX
 ! ---------------------------------------------------------------------------- !        !! WAM-MAX
 !                                                                              !        !! WAM-MAX
 !     INTERFACE VARIABLES.                                                     !        !! WAM-MAX
 !     --------------------                                                     !        !! WAM-MAX
 !                                                                              !        !! WAM-MAX
REAL,    INTENT(IN)    :: F(:,:,:)     !! BLOCK OF SPECTRA                               !! WAM-MAX
REAL,    INTENT(IN)    :: DEPTH(:)     !! DEPTH                                          !! WAM-MAX
REAL,    INTENT(IN)    :: THMAX(:)     !! PEAK DIRECTION                                 !! WAM-MAX
REAL,    INTENT(OUT)   :: CMAX_F(:)    !! MAXIMUM CREST H.- TIME (FORRISTALL)            !! WAM-MAX
REAL,    INTENT(OUT)   :: HMAX_N(:)    !! MAXIMUM WAVE H.- TIME (NAESS)                  !! WAM-MAX
REAL,    INTENT(OUT)   :: CMAX_ST(:)   !! MAXIMUM CREST H.- SPACE-TIME (STQD)            !! WAM-MAX
REAL,    INTENT(OUT)   :: HMAX_ST(:)   !! MAXIMUM WAVE H.- SPACE-TIME (STQD)             !! WAM-MAX
LOGICAL, INTENT(IN),  OPTIONAL   :: MASK(:,:,:)    !! INTEGRATION MASK                   !! WAM-MAX
 ! ---------------------------------------------------------------------------- !        !! WAM-MAX
 !                                                                              !        !! WAM-MAX
 !     LOCAL VARIABLES.                                                         !        !! WAM-MAX
 !     ----------------                                                         !        !! WAM-MAX
                                                                                         !! WAM-MAX
! REAL, PARAMETER :: GAMMA_E = 0.57721566  !! EULER CONSTANT (moved as global)           !! WAM-MAX
                                                                                         !! WAM-MAX
INTEGER :: IJ, T, M, K, MAXIT                                                            !! WAM-MAX
REAL    :: Z0, NW, ALFA, BETA, URSN, STEEP, WNUM1, PHIST, XK, AXYT, N3, N2, N1, HS       !! WAM-MAX
REAL    :: LX(SIZE(F,1)), LY(SIZE(F,1)), AXT(SIZE(F,1)), AYT(SIZE(F,1)), AXY(SIZE(F,1))  !! WAM-MAX
REAL    :: TEMP_X(SIZE(F,1),SIZE(F,3)), TEMP_Y(SIZE(F,1),SIZE(F,3))                      !! WAM-MAX
REAL    :: TEMP_X2(SIZE(F,1),SIZE(F,3)), TEMP_Y2(SIZE(F,1),SIZE(F,3))                    !! WAM-MAX
REAL    :: TEMP_XY(SIZE(F,1),SIZE(F,3)), NI(SIZE(F,1)), MU(SIZE(F,1))                    !! WAM-MAX
REAL    :: CX(SIZE(F,1),SIZE(F,2)), CY(SIZE(F,1),SIZE(F,2))                              !! WAM-MAX
REAL    :: CCX(SIZE(F,1),SIZE(F,2),SIZE(F,3)), CCY(SIZE(F,1),SIZE(F,2),SIZE(F,3))        !! WAM-MAX
REAL    :: TEMP(SIZE(F,1),SIZE(F,3))                                                     !! WAM-MAX
REAL    :: ACF(SIZE(F,1)), TLG(SIZE(F,1))                                                !! WAM-MAX
REAL    :: T1(SIZE(F,1)), T2(SIZE(F,1)), EMEAN(SIZE(F,1))                                !! WAM-MAX
REAL    :: TLGS(4), ACFS(4)                                                              !! WAM-MAX
REAL, PARAMETER   :: GRR = (1.+SQRT(5.))/2. !! GOLDEN RATIO                              !! WAM-MAX
REAL    :: TOL                            ! GOLDEN SEARCH TOLERANCE                      !! WAM-MAX

! ---------------------------------------------------------------------------- !         !! WAM-MAX
!                                                                              !         !! WAM-MAX
!     0. COMPUTE SPECTRAL PARAMETER BY INTEGRATION OF SPECTRA.                 !         !! WAM-MAX
!        --------------------------------------------------------------------- !         !! WAM-MAX

! INITIALIZE VARIABLES FOR INTEGRATION                                                   !! WAM-MAX

ACF = 0.                                                                                 !! WAM-MAX
T1 = 0.                                                                                  !! WAM-MAX
T2 = 0.                                                                                  !! WAM-MAX
LX = 0.                                                                                  !! WAM-MAX
LY = 0.                                                                                  !! WAM-MAX
AXT = 0.                                                                                 !! WAM-MAX
AXY = 0.                                                                                 !! WAM-MAX
AYT = 0.                                                                                 !! WAM-MAX
EMEAN = 0.                                                                               !! WAM-MAX

! COMPUTE MEAN PERIODS (WITHOUT TAIL CONTRIBUTION)                                       !! WAM-MAX
! COMPUTE MEAN WAVE AND CREST LENGTH (WITHOUT TAIL CONTRIBUTION, WRT PEAK DIR.)          !! WAM-MAX
! COMPUTE IRREGULARITY PARAMETERS (WITHOUT TAIL CONTRIBUTION, WRT PEAK DIR.)             !! WAM-MAX

DO K = 1,KL                                                                              !! WAM-MAX
   DO IJ = 1,SIZE(TEMP,1)                                                                !! WAM-MAX
      CX(IJ,K) = COSTH(K)*COS(THMAX(IJ)/180.*PI)+SINTH(K)*SIN(THMAX(IJ)/180.*PI)         !! WAM-MAX
      CY(IJ,K) = SINTH(K)*COS(THMAX(IJ)/180.*PI)-COSTH(K)*SIN(THMAX(IJ)/180.*PI)         !! WAM-MAX
   END DO                                                                                !! WAM-MAX
END DO                                                                                   !! WAM-MAX
CCX(1:SIZE(F,1),1:KL,1:ML) = SPREAD(CX,3,ML)                                             !! WAM-MAX
CCY(1:SIZE(F,1),1:KL,1:ML) = SPREAD(CY,3,ML)                                             !! WAM-MAX

IF (PRESENT(MASK)) THEN                                                                  !! WAM-MAX
   TEMP = SUM(F, DIM=2, MASK=MASK)                                                       !! WAM-MAX
   TEMP_X = SUM(F*CCX, DIM=2, MASK=MASK)                                                 !! WAM-MAX
   TEMP_Y = SUM(F*CCY, DIM=2, MASK=MASK)                                                 !! WAM-MAX
   TEMP_X2 = SUM(F*CCX**2, DIM=2, MASK=MASK)                                             !! WAM-MAX
   TEMP_Y2 = SUM(F*CCY**2, DIM=2, MASK=MASK)                                             !! WAM-MAX
   TEMP_XY = SUM(F*CCX*CCY, DIM=2, MASK=MASK)                                            !! WAM-MAX
ELSE                                                                                     !! WAM-MAX
   TEMP = SUM(F, DIM=2)                                                                  !! WAM-MAX
   TEMP_X = SUM(F*CCX, DIM=2)                                                            !! WAM-MAX
   TEMP_Y = SUM(F*CCY, DIM=2)                                                            !! WAM-MAX
   TEMP_X2 = SUM(F*CCX**2, DIM=2)                                                        !! WAM-MAX
   TEMP_Y2 = SUM(F*CCY**2, DIM=2)                                                        !! WAM-MAX
   TEMP_XY = SUM(F*CCX*CCY, DIM=2)                                                       !! WAM-MAX
END IF                                                                                   !! WAM-MAX

DO M = 1,ML                                                                              !! WAM-MAX
   DO IJ = 1,SIZE(TEMP_X2,1)                                                             !! WAM-MAX
      XK = AKI(2.*PI*FR(M),DEPTH(IJ))                                                    !! WAM-MAX
      EMEAN(IJ) = EMEAN(IJ) + TEMP(IJ,M)*DFIM(M)                                         !! WAM-MAX
      T1(IJ) = T1(IJ) + TEMP(IJ,M)*DFIM_FR(M)                                            !! WAM-MAX
      T2(IJ) = T2(IJ) + TEMP(IJ,M)*DFIM_FR2(M)                                           !! WAM-MAX
      LX(IJ) = LX(IJ) + TEMP_X2(IJ,M)*XK**2*DFIM(M)                                      !! WAM-MAX
      LY(IJ) = LY(IJ) + TEMP_Y2(IJ,M)*XK**2*DFIM(M)                                      !! WAM-MAX
      AXY(IJ) = AXY(IJ) + TEMP_XY(IJ,M)*XK**2*DFIM(M)                                    !! WAM-MAX
      AXT(IJ) = AXT(IJ) + TEMP_X(IJ,M)*XK*(2.*PI*FR(M))*DFIM(M)                          !! WAM-MAX
      AYT(IJ) = AYT(IJ) + TEMP_Y(IJ,M)*XK*(2.*PI*FR(M))*DFIM(M)                          !! WAM-MAX
   END DO                                                                                !! WAM-MAX
END DO                                                                                   !! WAM-MAX

WHERE (EMEAN.GT.EMIN)                                                                    !! WAM-MAX
   AXY = AXY/SQRT(LX*LY)                                                                 !! WAM-MAX
   AXT = AXT/(2.*PI*SQRT(LX*T2))                                                         !! WAM-MAX
   AYT = AYT/(2.*PI*SQRT(LY*T2))                                                         !! WAM-MAX
   LX = 2.*PI*SQRT(EMEAN/LX)                                                             !! WAM-MAX
   LY = 2.*PI*SQRT(EMEAN/LY)                                                             !! WAM-MAX
   NI = SQRT(EMEAN*T2/T1**2 - 1)                                                         !! WAM-MAX
   MU = (2.*PI*T1)**2*EMEAN**(-3./2.)*(1.-NI+NI**2)/G                                    !! WAM-MAX
   T1 = EMEAN/T1                                                                         !! WAM-MAX
   T2 = SQRT(EMEAN/T2)                                                                   !! WAM-MAX
ELSEWHERE                                                                                !! WAM-MAX
   AXY = 0.                                                                              !! WAM-MAX
   AXT = 0.                                                                              !! WAM-MAX
   AYT = 0.                                                                              !! WAM-MAX
   LX = 1.                                                                               !! WAM-MAX
   LY = 1.                                                                               !! WAM-MAX
   NI = 0.                                                                               !! WAM-MAX
   MU = 0.                                                                               !! WAM-MAX
   T1 = 1.                                                                               !! WAM-MAX
   T2 = 1.                                                                               !! WAM-MAX
END WHERE                                                                                !! WAM-MAX
!                                                                                        !! WAM-MAX
! ---------------------------------------------------------------------------- !         !! WAM-MAX
!                                                                              !         !! WAM-MAX
!     0. MIN OF AUTOCOVARIANCE FUNCTION                                        !         !! WAM-MAX
!        VIA GOLDEN RATIO SEARCH                                               !         !! WAM-MAX
!        --------------------------------------------------------------------- !         !! WAM-MAX
MAXIT = 10                                                                               !! WAM-MAX
TOL = 0.01                                                                               !! WAM-MAX
DO IJ = 1,SIZE(TEMP,1)                                                                   !! WAM-MAX
   TLGS(1) = 0.3*T2(IJ)                                                                  !! WAM-MAX
   TLGS(4) = 1.3*T2(IJ)                                                                  !! WAM-MAX
   TLGS(2) = TLGS(4) - (TLGS(4) - TLGS(1))/GRR                                           !! WAM-MAX
   TLGS(3) = TLGS(1) + (TLGS(4) - TLGS(1))/GRR                                           !! WAM-MAX
   ACFS(1) = SUM( COS( 2.*PI*FR(:)*TLGS(1) ) * TEMP(IJ,:) * DFIM(:) )                    !! WAM-MAX
   ACFS(4) = SUM( COS( 2.*PI*FR(:)*TLGS(4) ) * TEMP(IJ,:) * DFIM(:) )                    !! WAM-MAX
   ACFS(2) = SUM( COS( 2.*PI*FR(:)*TLGS(2) ) * TEMP(IJ,:) * DFIM(:) )                    !! WAM-MAX
   ACFS(3) = SUM( COS( 2.*PI*FR(:)*TLGS(3) ) * TEMP(IJ,:) * DFIM(:) )                    !! WAM-MAX
   DO T = 1,MAXIT                                                                        !! WAM-MAX
      IF (ACFS(2) .LT. ACFS(3)) THEN                                                     !! WAM-MAX
         ACF(IJ) = ACFS(2)                                                               !! WAM-MAX
         TLG = TLGS(2)                                                                   !! WAM-MAX
         TLGS(4) = TLGS(3)                                                               !! WAM-MAX
         ACFS(4) = ACFS(3)                                                               !! WAM-MAX
         TLGS(3) = TLGS(2)                                                               !! WAM-MAX
         ACFS(3) = ACFS(2)                                                               !! WAM-MAX
         TLGS(2) = TLGS(4) - (TLGS(4) - TLGS(1))/GRR                                     !! WAM-MAX
         ACFS(2) = SUM( COS( 2.*PI*FR(:)*TLGS(2) ) * TEMP(IJ,:) * DFIM(:) )              !! WAM-MAX
      ELSE                                                                               !! WAM-MAX
         ACF(IJ) = ACFS(3)                                                               !! WAM-MAX
         TLGS(1) = TLGS(2)                                                               !! WAM-MAX
         ACFS(1) = ACFS(2)                                                               !! WAM-MAX
         TLGS(2) = TLGS(3)                                                               !! WAM-MAX
         ACFS(2) = ACFS(3)                                                               !! WAM-MAX
         TLG = TLGS(3)                                                                   !! WAM-MAX
         TLGS(3) = TLGS(1) + (TLGS(4) - TLGS(1))/GRR                                     !! WAM-MAX
         ACFS(3) = SUM( COS( 2.*PI*FR(:)*TLGS(3) ) * TEMP(IJ,:) * DFIM(:) )              !! WAM-MAX
      END IF                                                                             !! WAM-MAX
      IF ( (ABS(TLGS(4)-TLGS(1))) .LT. (TOL*(ABS(TLGS(2))+ABS(TLGS(3)))) ) THEN          !! WAM-MAX
         EXIT                                                                            !! WAM-MAX
      END IF                                                                             !! WAM-MAX
   END DO                                                                                !! WAM-MAX
END DO                                                                                   !! WAM-MAX
!                                                                                        !! WAM-MAX
! ---------------------------------------------------------------------------- !         !! WAM-MAX
!                                                                              !         !! WAM-MAX
!     1. COMPUTE OUTPUT VARIABLES                                              !         !! WAM-MAX
!        --------------------------------------------------------------------- !         !! WAM-MAX
DO IJ = 1,SIZE(DEPTH)                                                                    !! WAM-MAX
! ---------------------------------------------------------------------------- !         !! WAM-MAX
!                                                                              !         !! WAM-MAX
!         1a. DETERMINE EXPECTED MAXIMUM CREST HEIGHT - TIME (FORRISTALL).     !         !! WAM-MAX
!        --------------------------------------------------------------------- !         !! WAM-MAX
   IF (T1(IJ).GT.0 .AND. T2(IJ).GT.0 .AND. DEPTH(IJ).GT.0 .AND. EMEAN(IJ).GT.0) THEN     !! WAM-MAX
      HS = 4*SQRT(EMEAN(IJ))                                                             !! WAM-MAX
      WNUM1 = AKI(2.*PI/T1(IJ),DEPTH(IJ))                                                !! WAM-MAX
      STEEP = 2.*PI*HS/G/T1(IJ)**2                                                       !! WAM-MAX
      URSN = HS/WNUM1**2/DEPTH(IJ)**3                                                    !! WAM-MAX
      ALFA = 0.3536+0.2568*STEEP+0.08*URSN                                               !! WAM-MAX
      BETA = 2.-1.7912*STEEP-0.5302*URSN+0.284*URSN**2                                   !! WAM-MAX
      NW = WMDUR/T2(IJ)                                                                  !! WAM-MAX
      Z0 = LOG(REAL(NW))                                                                 !! WAM-MAX
      CMAX_F(IJ) = ALFA*Z0**(1./BETA)*(1.+GAMMA_E/(BETA*Z0))*HS                          !! WAM-MAX
! ---------------------------------------------------------------------------- !         !! WAM-MAX
!                                                                              !         !! WAM-MAX
!         1b. DETERMINE EXPECTED MAXIMUM WAVE HEIGHT - TIME (NAESS).           !         !! WAM-MAX
!        --------------------------------------------------------------------- !         !! WAM-MAX
      PHIST = ACF(IJ)/EMEAN(IJ)                                                          !! WAM-MAX
      HMAX_N(IJ) = 0.5*SQRT(1.-PHIST)*SQRT(Z0)*(1.+0.5*GAMMA_E/Z0)*HS                    !! WAM-MAX
! ---------------------------------------------------------------------------- !         !! WAM-MAX
!                                                                              !         !! WAM-MAX
!         2a. DETERMINE EXPECTED MAXIMUM CREST HEIGHT - SPACE-TIME (STQD).     !         !! WAM-MAX
!        --------------------------------------------------------------------- !         !! WAM-MAX

      AXYT = SQRT(1+2*AXT(IJ)*AXY(IJ)*AYT(IJ)-AXT(IJ)**2-AXY(IJ)**2-AYT(IJ)**2)          !! WAM-MAX
      N3 = 2.*PI*WMDX*WMDY*WMDUR/LX(IJ)/LY(IJ)/T2(IJ)*AXYT                               !! WAM-MAX
      N2 = SQRT(2.*PI)*(WMDX*WMDUR/LX(IJ)/T2(IJ)*SQRT(1-AXT(IJ)**2) &                    !! WAM-MAX
&          +WMDX*WMDY/LX(IJ)/LY(IJ)*SQRT(1-AXY(IJ)**2) &                                 !! WAM-MAX
&          +WMDY*WMDUR/LY(IJ)/T2(IJ)*SQRT(1-AYT(IJ)**2))                                 !! WAM-MAX
      N1 = WMDX/LX(IJ) + WMDY/LY(IJ) + WMDUR/T2(IJ)                                      !! WAM-MAX

      IF (WMDX.NE.0 .AND. WMDY.NE.0 .AND. WMDUR.NE.0) THEN                               !! WAM-MAX
         Z0 = SQRT(2.*LOG(REAL(N3))+2.*LOG(2.*LOG(REAL(N3))+2.*LOG(2.*LOG(REAL(N3)))))   !! WAM-MAX
      ELSE IF (WMDUR.EQ.0) THEN                                                          !! WAM-MAX
         Z0 = SQRT(2.*LOG(REAL(N2))+LOG(2.*LOG(REAL(N2))+LOG(2.*LOG(REAL(N2)))))         !! WAM-MAX
      ELSE IF (WMDX.EQ.0 .AND. WMDY.EQ.0) THEN                                           !! WAM-MAX
         Z0 = SQRT(2.*LOG(REAL(N1)))                                                     !! WAM-MAX
      END IF                                                                             !! WAM-MAX
      CMAX_ST(IJ) = ((Z0+0.5*MU(IJ)*Z0**2)+GAMMA_E*((1.+MU(IJ)*Z0) &                     !! WAM-MAX
&                   *(Z0-(2.*N3*Z0+N2)/(N3*Z0**2+N2*Z0+N1))**(-1))) *SQRT(EMEAN(IJ))     !! WAM-MAX
! ---------------------------------------------------------------------------- !         !! WAM-MAX
!                                                                              !         !! WAM-MAX
!         2b. DETERMINE EXPECTED MAXIMUM WAVE HEIGHT - SPACE-TIME (STQD).      !         !! WAM-MAX
!        --------------------------------------------------------------------- !         !! WAM-MAX

      HMAX_ST(IJ) = (Z0+GAMMA_E*(Z0-(2.*N3*Z0+N2)/(N3*Z0**2+N2*Z0+N1))**(-1))* &         !! WAM-MAX
&                   SQRT(2.*(1.-PHIST)) *SQRT(EMEAN(IJ))                                 !! WAM-MAX
   ELSE                                                                                  !! WAM-MAX
      CMAX_F(IJ) = 0.                                                                    !! WAM-MAX
      HMAX_N(IJ) = 0.                                                                    !! WAM-MAX
      CMAX_ST(IJ) = 0.                                                                   !! WAM-MAX
      HMAX_ST(IJ) = 0.                                                                   !! WAM-MAX
   END IF                                                                                !! WAM-MAX
END DO                                                                                   !! WAM-MAX
!close (87)
                                                                                         !! WAM-MAX
END SUBROUTINE WAMAX                                                                     !! WAM-MAX

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

REAL FUNCTION TRANSF (XK, D)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    TRANSF   DETERMINE NARROW BAND LIMIT BENJAMIN-FEIR INDEX FOR              !
!              THE FINITE DEPTH CASE.
!                                                                              !
!     PETER JANSSEN       JULY 2005.
!                                                                              !
!     PURPOSE.
!     --------
!                                                                              !
!           DETERMINATION OF THE NARROW BAND LIMIT BENJAMIN-FEIR INDEX FOR     !
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

REAL ::  EPS, X, T_0, OM, C_0, V_G, DV_G, XNL_1, XNL_2, XNL

EPS = 0.0001

! ---------------------------------------------------------------------------- !
!
!     1. DETERMINE TRANSFER FUNCTION.
!        ----------------------------

IF (D.LT.999. .AND. D.GT.0.) THEN
   X   = XK*D
   IF (X .GT. DKMAX) THEN
      TRANSF = 1.
   ELSE
      T_0 = TANH(X)
      OM  = SQRT(G*XK*T_0)
      C_0 = OM/XK
      IF (X .LT. EPS) THEN
         V_G = 0.5*C_0
         V_G = C_0
      ELSE
         V_G = 0.5*C_0*(1.+2.*X/SINH(2.*X))
      ENDIF
      DV_G = (T_0-X*(1.-T_0**2))**2+4.*X**2*T_0**2*(1.-T_0**2)

      XNL_1 = (9.*T_0**4-10.*T_0**2+9.)/(8.*T_0**3)
      XNL_2 = ((2.*V_G-0.5*C_0)**2/(G*D-V_G**2)+1.)/X

      XNL = XNL_1-XNL_2
      TRANSF = XNL**2/(DV_G*T_0**8)
   ENDIF
ELSE
   TRANSF = 1.
ENDIF

END FUNCTION TRANSF

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
!           DETERMINATION OF THE NARROW BAND LIMIT BENJAMIN-FEIR INDEX FOR     !
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
X   = MAX(X, 0.05)
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

END MODULE WAM_INTERFACE_MODULE

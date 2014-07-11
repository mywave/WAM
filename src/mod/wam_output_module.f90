MODULE WAM_OUTPUT_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS: ARRAYS NECESSARY FOR GRIDDED FIELDS OF PARAMTERS.    !
!                         ALL PROCEEDURES TO COMPUTE AND WRITE OUTPUT.         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1,                  &  !! TERMINATE PROCESSING.
&       OPEN_FILE,               &  !! OPENS A FILE.
&       INCDATE,                 &  !! UPDATE DATE TIME GROUP.
&       PRINT_ARRAY,             &  !! PRINTER OUTPUT OF AN ARRAY.
&       PRINT_SPECTRUM,          &  !! PRINT A SPECTRUM.
&       difdate                     !! difference between two dates in seconds

USE WAM_INTERFACE_MODULE, ONLY:  &
&       FEMEAN,                  &  !! COMPUTATION OF MEAN FREQUENCY.
&       MEAN_DIRECTION,          &  !! COMPUTATION OF MEAN DIRECTION AND SPREAD.
&       PEAK_PERIOD,             &  !! COMPUTATION OF PEAK PERIOD.
&       TM1_TM2_PERIODS,         &  !! COMPUTATION OF TM1 AND/OR TM2 PERIOD.
&       TOTAL_ENERGY,            &  !! COMPUTATION OF TOTAL ENERGY.
&       MEANSQS,                 &  !! COMPUTATION OF TMEAN SQUARE SLOPE.
&       CHARNOCK_PAR,            &  !! COMPUTATION OF CHARNOCK_PARAMETER,
&       KURTOSIS                    !! COMPUTATION OF KURTOSIS,
                                    !! BENJAMIN-FEIR INDEX,
                                    !! GODA'S PEAKEDNESS PARAMETER,
                                    !! PEAK FREQUENCY (INTERPOLATED),
                                    !! PEAK DIRECTION (INTERPOLATED),
                                    !! NORMALIZED MAXIMUM WAVE HEIGHT,
                                    !! MAXIMUM WAVE PERIOD

USE WAM_OUTPUT_SET_UP_MODULE, ONLY:  &
&       SAVE_OUTPUT_FILES,           & !! CLOSES AND OPENS OUTPUT FILES.
&       UPDATE_OUTPUT_TIME             !! UPDATES OUTPUT TIMES.

USE WAM_ICE_MODULE,    ONLY:     & 
&       PUT_ICE                     !! PUTS ICE INDICATOR INTO DATA FIELD.

USE WAM_TOPO_MODULE,   ONLY:     & 
&       PUT_DRY                     !! PUTS DRY INDICATOR INTO DATA FIELD.

USE WAM_SOURCE_MODULE, ONLY:     & 
&       AIRSEA                      !! EVALUATE USTAR.

use wam_mpi_comp_module, only:   &
&       mpi_gather_block,        &
&       mpi_gather_grid,         &
&       mpi_gather_spp

use wam_special_module, only:    &
&      create_ready_file_name

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE,ONLY: G, ZPI, DEG, RCHAR

USE WAM_FRE_DIR_MODULE,ONLY: KL, ML, CO, FR, C, DFIM, TH, COSTH, SINTH, TFAK

USE WAM_GRID_MODULE,   ONLY: NX, NY, NLON_RG, AMOWEP, AMOSOP, AMOEAP, AMONOP,  &
&                             L_S_MASK, NSEA, ZDELLO

USE WAM_FILE_MODULE,   ONLY: IU06, ITEST, IU20, IU25, FILE20, FILE25,          &
&                            area, iu67

USE WAM_MODEL_MODULE,  ONLY: FL3, U10, UDIR, USTAR, TAUW, Z0, DEPTH, INDEP, U, V

USE WAM_TIMOPT_MODULE, ONLY: IDELPRO, CDTPRO,                                  &
&                            SHALLOW_RUN, REFRACTION_C_RUN, cdatea

USE WAM_OUTPUT_SET_UP_MODULE, ONLY:                                            &
&       CDTINTT, CDTSPT, IDELINT, IDELSPT, NOUTT, COUTT, CDT_OUT,              &
&       NOUT_P, FFLAG_P, FFLAG20, PFLAG_P, PFLAG20, CFLAG_P, CFLAG20,          &
&       TITL_P, SCAL_P,                                                        &
&       NOUT_S, FFLAG_S, FFLAG25, PFLAG_S, PFLAG25, CFLAG_S, CFLAG25,          &
&       TITL_S, NOUTP, OUTLAT, OUTLONG, NAME, IJAR,                            &
&       ready_outf, owpath

USE WAM_ICE_MODULE,    ONLY: ICE_RUN
USE WAM_TOPO_MODULE,   ONLY: N_DRY

use wam_mpi_module,    only: irank, nijs, nijl, ipfgtbl, i_out_par, i_out_spec

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
include 'mpif.h'
PRIVATE

INTEGER :: I, ishift

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. BLOCKED INTEGRATED PARAMETER TOTAL SPECTRUM.                          !
!        --------------------------------------------                          !

REAL, ALLOCATABLE, DIMENSION(:,:) :: BLOCK  !! FIRST INDEX COUNTS SEAPOINTS
                                            !! SECOND INDEX COUNTS PARAMETER.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE MODEL_OUTPUT_CONTROL        !! CONTROLS MODEL OUTPUT.
   MODULE PROCEDURE MODEL_OUTPUT_CONTROL
END INTERFACE
PUBLIC MODEL_OUTPUT_CONTROL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE COMPUTE_OUTPUT_PARAMETER    !! COMPUTES OUTPUT PARAMETER.
   MODULE PROCEDURE COMPUTE_OUTPUT_PARAMETER
END INTERFACE
PRIVATE  COMPUTE_OUTPUT_PARAMETER

INTERFACE SIGMA_TO_OMEGA              !! MAP SPECTRUM FROM SIGMA TO OMEGA SPACE.
   MODULE PROCEDURE SIGMA_TO_OMEGA
END INTERFACE
PRIVATE SIGMA_TO_OMEGA

INTERFACE SWELL_SEPARATION            !! SWELL SEPARATION AND INTEGRATED
   MODULE PROCEDURE SWELL_SEPARATION  !! PARAMETER OF SEA AND SWELL.
END INTERFACE
PRIVATE SWELL_SEPARATION

INTERFACE WRITE_MODEL_OUTPUT           !! WRITE MODEL OUTPUT.
   MODULE PROCEDURE WRITE_MODEL_OUTPUT
END INTERFACE
PRIVATE WRITE_MODEL_OUTPUT

INTERFACE WRITE_INT_PAR_OUTPUT        !! WRITE OUTPUT OF INTEGRATED PARAMETER.
   MODULE PROCEDURE WRITE_INT_PAR_OUTPUT
END INTERFACE
PRIVATE  WRITE_INT_PAR_OUTPUT

INTERFACE WRITE_SPECTRA_OUTPUT        !! WRITE OUTPUT OF SPECTRA.
   MODULE PROCEDURE WRITE_SPECTRA_OUTPUT
END INTERFACE
PRIVATE WRITE_SPECTRA_OUTPUT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MODEL_OUTPUT_CONTROL (fl3, iu20, iu25)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MODEL_OUTPUT_CONTROL - CONTROLS MODEL OUTPUT.                              !
!                                                                              !
!     H. GUNTHER         GKSS/ECMWF         JUNE 1990                          !
!                                                                              !
!    PURPOSE.                                                                  !
!     --------                                                                 !
!                                                                              !
!       CONTROL OUTPUT OF WAVE, WIND, DEPTH AND CURRENT FIELDS.                !
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
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL, INTENT(IN)    :: FL3(:,:,:)  !! BLOCK OF SPECTRA.
INTEGER, INTENT(IN) :: IU20        !! PARAMETER UNIT NUMBER.
INTEGER, INTENT(IN) :: IU25        !! SPECTRA UNIT NUMBER.

!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, ALLOCATABLE  :: FL1(:,:,:)  !! BLOCK OF SPECTRA.
REAL, ALLOCATABLE  :: FL (:,:,:)  !! SWELL SPECTRA.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. ALLOCATE PARAMETER ARRAYS.                                            !
!        --------------------------                                            !

! ARRAY FOR INTEGRATED PARAMETERS

IF (ANY(CFLAG_P(:))) THEN
   ALLOCATE(BLOCK(nijs:nijl,1:NOUT_P))
   BLOCK = -1.
END IF

!  ARRAY FOR SWELL-SEA SEPARATION

IF (ANY(CFLAG_P(17:32)).OR.ANY(CFLAG_S(2:NOUT_S)) ) THEN
   ALLOCATE (FL(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)))
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. PROCESS OUTPUT.                                                       !
!        ---------------                                                       !

IF (REFRACTION_C_RUN) THEN
   IF (.NOT.ALLOCATED(FL1)) ALLOCATE(FL1(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)))
   CALL SIGMA_TO_OMEGA (FL3, FL1, u(nijs:nijl), v(nijs:nijl),                  &
&                      indep(nijs:nijl))
   CALL COMPUTE_OUTPUT_PARAMETER (FL1, FL)
   CALL WRITE_MODEL_OUTPUT (FL1, FL, IU20, IU25)
ELSE
   CALL COMPUTE_OUTPUT_PARAMETER (FL3, FL)
   CALL WRITE_MODEL_OUTPUT (FL3, FL, IU20, IU25)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. DEALLOCATE ARRAYS.                                                    !
!        ------------------                                                    !

IF (ALLOCATED (BLOCK) ) DEALLOCATE (BLOCK)
IF (ALLOCATED (FL   ) ) DEALLOCATE (FL   )
IF (ALLOCATED (FL1  ) ) DEALLOCATE (FL1  )

END SUBROUTINE MODEL_OUTPUT_CONTROL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SIGMA_TO_OMEGA (F3, F1, u, v, indep)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SIGMA_TO_OMEGA - TRANSFORMATION OF SPECTRA FROM SIGMA TO OMEGA.            !
!                                                                              !
!     S.D.HASSELMANN      MPI            1.1.91                                !
!     H. GUNTHER          GKSS/ECMWF     1.2.01  MODIFIED FOR CYCLE_4          !
!     H. GUNTHER          HZG/ECMWF     1.11.10  NEW INTERPOLATION.            !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFORMATION OF A MOVING COORDINATE SYSTEM TO AN ABSOLUTE            !
!       COORDINATE SYSTEM.                                                     !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       ABSOLUTE FREQUENCIES (OMEGA) ARE COMPUTED FROM THE MOVING FREQUENCIES  !
!       (SIGMA) AND  THE SIGMA SPECTRUN IS TRANSFORMED TO THE OMEGA SPECTRUM   !
!       BY MULTIPLICATION WITH THE DELTA_SIGMAN/DELT_OMEGA.                    !
!       THE OMEGA SPECTRUM IS THAN INTERPOLATED TO ABSOLUTE OUTPUT FREQUENCIES !
!       BY RE-DISTRIBUTING THE ENGERY IN OVERLAPPING FREQUENCY INTERVALLS.     !
!       THE ABSOLUTE OUTPUT FREQUENCIES ARE THE MOVING INTPUT FREQUENCIES.     !
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
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: F3(:,:,:)     !! BLOCK OF SPECTRA.
REAL,    INTENT(OUT)   :: F1(:,:,:)     !! BLOCK OF SPECTRA.
REAL,    INTENT(IN)    :: U(:)          !! U COMPONENT OF CURRENT.
REAL,    INTENT(IN)    :: V(:)          !! V COMPONENT OF CURRENT.
INTEGER, INTENT(IN)    :: INDEP(:)      !! DEPTH INDEX

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: PI2G = ZPI/G

INTEGER  :: K, KN, M1, M2,  LL, IJ
REAL     :: FRR, F1L, F1R, F2L, F2R

REAL     :: OM(SIZE(F3,3))        !! ABSOLUTE FREQUENCIES.
REAL     :: F1I(1:SIZE(F3,3),1:2) !! INTERVALL BOARDER OF RELATIVE FREQUENCIES.
REAL     :: DF1I(1:SIZE(F3,3))    !! INTERVALL SIZE OF RELATIVE FREQUENCIES.
REAL     :: WAVN(SIZE(F3,3))      !! WAVE NUMBERS OF RELATIVE FREQUENCIES.
REAL     :: F2I(1:SIZE(F3,3),1:2) !! INTERVALL BOARDER OF ABSOLUTE FREQUENCIES.
REAL     :: E_IN(SIZE(F3,3))      !! SPECTRUM AT ABSOLUTE FREQUENCIES.
REAL     :: E_OUT(SIZE(F3,3))     !! SPECTRUM AT OUTPUT ABSOLUTE FREQUENCIES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. INITIAL OUTPUT ARRAY WITH ZERO.                                       !
!        -------------------------------                                       !

F1 = 0.0
F1I(:,:)  = INTERVAL_BOARDERS (FR)   !! BOARDER FREQUENCIES.
DF1I(:)   = F1I(:,2) - F1I(:,1)      !! FREQUENCY INTERVALL SIZES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. LOOP OVER GRID POINTS.                                                !
!        ----------------------                                                !

POINT: DO IJ = 1,SIZE(F3,1)

!   IF (ABS(U(IJ)).LT.0.001 .AND. ABS(V(IJ)).LT.0.001) THEN
!      F1(IJ,:,:) = F3(IJ,:,:)
!      CYCLE POINT
!   END IF

!     1.1 WAVE NUMBER FOR ALL RELATIVE FREQUENCIES.

   IF (SHALLOW_RUN) THEN
      WAVN(:) = TFAK(INDEP(IJ),:)/ZPI
   ELSE
      WAVN(:) = PI2G*FR(:)*FR(:)
   END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LOOP OVER DIRECTIONS.                                                 !
!        ---------------------                                                 !

   DIR: DO K = 1, KL

!     2.1 ABSOLUTE FREQUENCIES FOR A FIXED DIRECTION.

      OM(:) = FR(:) + WAVN(:)*(COSTH(K)*V(IJ) + SINTH(K)*U(IJ))
      LL = COUNT(OM.GT.0.)         !! INDEX OF LAST POSITVE ABSOLUTE FREQUENCY.
      F2I = INTERVAL_BOARDERS (OM) !! ABSOLUTE FREQUENCY INTERVALL BOARDES.

!     2.2 ENERGY DENSITIES AT ABSOLUTE FREQUENCIES FOR A FIXED DIRECTION.

     WHERE (F2I(1:ML,2).NE.F2I(1:ML,1)) 
         E_IN(:) = F3(IJ,K,:) * DF1I(:)/ (F2I(1:ML,2)-F2I(1:ML,1))
      ELSEWHERE
         E_IN(:) = 0.
      END WHERE
      
!     2.3 INTEPOLATE TO NEW ABSOLUTE FREQUENCIES.

      IF (LL.GT.0) THEN       !! FOR POSITIVE OMEGAS: KEEP DIRECTION.
         E_OUT(:) = 0.
         DO M1 = 1,ML         !! LOOP OVER NEW COORDINATES.
            F1L = F1I(M1,1)
            F1R = F1I(M1,2)
            DO M2 = 1,LL      !! LOOP OVER OLD COORDINATES.
               F2L = F2I(M2,1)
               F2R = F2I(M2,2)
               FRR = MAX(0.,MIN(F2R,F1R)-MAX(F2L,F1L))  !! COMMON AREA.
               E_OUT(M1) = E_OUT(M1) + E_IN(M2) * FRR
            END DO
            F1(IJ,K,M1) = F1(IJ,K,M1) + E_OUT(M1)/DF1I(M1)
         END DO
      END IF

      IF (LL.LT.ML) THEN      !! FOR NEGATIVE OMEGAS: TURN DIRECTION. 
         KN = MOD(K+KL/2-1,KL) + 1
         E_OUT(:) = 0.
         DO M1 = 1,ML         !! LOOP OVER NEW COORDINATES.
            F1L = F1I(M1,1)
            F1R = F1I(M1,2)
            DO M2 = LL+1,ML   !! LOOP OVER OLD COORDINATES.
               F2L = F2I(M2,1)
               F2R = F2I(M2,2)
               FRR = MAX(0.,MIN(F2R,F1R)-MAX(F2L,F1L))  !! COMMON AREA.
               E_OUT(M1) = E_OUT(M1) + E_IN(M2) * FRR
            END DO
            F1(IJ,KN,M1) = F1(IJ,KN,M1) + E_OUT(M1)/DF1I(M1)
         END DO
      END IF
      
   END DO DIR
END DO POINT

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

   FUNCTION INTERVAL_BOARDERS (A)  RESULT(B)

   REAL, INTENT(IN) :: A(:)   !! INTERVAL CENTER
   REAL :: B(1:SIZE(A),1:2)   !! INTERVAL BOARDERS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. INITIAL OUTPUT ARRAY WITH ZERO.                                       !
!        -------------------------------                                       !

   INTEGER :: ML
   REAL    :: X(1:SIZE(A))   !! INTERVAL BOARDERS
   ML = SIZE(A)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTERVAL BOARDERS.                                                    !
!        ------------------                                                    !

   B(1,1) = MAX(0., A(1) - 0.5*(A(2)-A(1)))
   B(2:ML,1) = 0.5*(A(1:ML-1)+A(2:ML))
   B(1:ML-1,2) = B(2:ML,1)
   B(ML,2) =  A(ML) + 0.5*(A(ML)-A(ML-1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. IF NEIGHBOUR CENTRES CHANGE SIGN, FIX ZERO AS INTERVALL BOARDER.      !
!        -----------------------------------------------------------------     !

   WHERE (A(1:ML-1)*A(2:ML).LT.0.)
      B(1:ML-1,2) = 0.
      B(2:ML,1) = 0.
   END WHERE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. POSITIVE INTERVALL BOARDERS AND RIGHT BOARDERS GREATER THAN LEFT ONES.!
!        ----------------------------------------------------------------------!

   B = ABS(B)
   WHERE (B(1:ML,1).GT.B(1:ML,2))
      X = B(1:ML,1)
      B(1:ML,1) = B(1:ML,2)
      B(1:ML,2) = X
   END WHERE

   END FUNCTION INTERVAL_BOARDERS

END SUBROUTINE SIGMA_TO_OMEGA

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE COMPUTE_OUTPUT_PARAMETER (FL3, FL)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    COMPUTE_OUTPUT_PARAMETER - CALCULATES THE MODEL OUTPUT.                   !
!                                                                              !
!     H. GUNTHER         GKSS/ECMWF         JUNE 1990                          !
!                                                                              !
!    PURPOSE.                                                                  !
!     --------                                                                 !
!                                                                              !
!       COMPUTE  REQUESTED OUTPUT PARAMETERS OF WAVES AND WINDS.               !
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
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)  :: FL3(:,:,:)  !! BLOCK OF SPECTRA.
REAL,    INTENT(OUT) :: FL (:,:,:)  !! BLOCK OF SWELL SPECTRA.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL  :: TAU(NIJS:NIJL)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE OUTPUT PARAMETERS.                                            !
!        --------------------------                                            !

BLOCK = 0.

IF (MAXVAL(USTAR(:)).EQ.0.) THEN !! EVALUATE FRICTION VELOCITY
  CALL AIRSEA  (U10, TAUW, USTAR, Z0)
END IF
IF (CFLAG_P(1)) BLOCK(NIJS:NIJL,1) = U10(:)
IF (CFLAG_P(2)) BLOCK(NIJS:NIJL,2) = DEG * UDIR(:)
IF (CFLAG_P(3)) BLOCK(NIJS:NIJL,3) = USTAR(:)

IF (CFLAG_P(4).OR.CFLAG_P(16)) THEN
   TAU(:) = USTAR(:)**2 +0.0001
   IF (CFLAG_P(4) ) BLOCK(NIJS:NIJL,4) = TAU(:)/(U10(:)**2+0.01)
   IF (CFLAG_P(16)) BLOCK(NIJS:NIJL,16) = TAUW(:)/TAU(:)
END IF

IF (CFLAG_P(5)) CALL CHARNOCK_PAR (USTAR, Z0, BLOCK(NIJS:NIJL,5))

IF (CFLAG_P(6)) THEN
   IF (SHALLOW_RUN) THEN
      BLOCK(NIJS:NIJL,6) = DEPTH(NIJS:NIJL)
   ELSE
      BLOCK(NIJS:NIJL,6) = 999.
   END IF
END IF

IF (CFLAG_P(7)) THEN
   BLOCK(NIJS:NIJL,7) = SQRT(U(NIJS:NIJL)**2+V(NIJS:NIJL)**2)
END IF
IF (CFLAG_P(8)) THEN  
   BLOCK(NIJS:NIJL,8) = ATAN2(U(NIJS:NIJL),V(NIJS:NIJL))*DEG
   WHERE (BLOCK(NIJS:NIJL,8).LT.0.) BLOCK(NIJS:NIJL,8) = BLOCK(NIJS:NIJL,8)+360.
END IF

IF (CFLAG_P(9).OR.ANY(CFLAG_P(11:13)).OR.CFLAG_P(39)) THEN
   CALL TOTAL_ENERGY (FL3, BLOCK(NIJS:NIJL,9))
   IF (CFLAG_P(11).OR.CFLAG_P(39)) THEN
      CALL FEMEAN (FL3,  BLOCK(NIJS:NIJL,9), FM=BLOCK(NIJS:NIJL,11))
      IF (CFLAG_P(40)) THEN
         IF (SHALLOW_RUN) THEN
            CALL MEANSQS (FL3, USTAR, FM=BLOCK(NIJS:NIJL,11),                 &
&                         SM=BLOCK(NIJS:NIJL,40), INDEP = INDEP(NIJS:NIJL))
         ELSE
            CALL MEANSQS (FL3, USTAR, FM=BLOCK(NIJS:NIJL,11),                 &
&                           SM=BLOCK(NIJS:NIJL,40))
         END IF
      END IF
      IF (CFLAG_P(11)) BLOCK(NIJS:NIJL,11) = 1./BLOCK(NIJS:NIJL,11)
   END IF
END IF

IF (CFLAG_P(12).AND.CFLAG_P(13)) THEN
   CALL TM1_TM2_PERIODS (FL3, BLOCK(NIJS:NIJL,9), TM1=BLOCK(NIJS:NIJL,12),   &
&                                TM2=BLOCK(NIJS:NIJL,13))
ELSE IF (CFLAG_P(12)) THEN
   CALL TM1_TM2_PERIODS (FL3, BLOCK(NIJS:NIJL,9), TM1=BLOCK(NIJS:NIJL,12))
ELSE IF (CFLAG_P(13)) THEN
   CALL TM1_TM2_PERIODS (FL3, BLOCK(NIJS:NIJL,9), TM2=BLOCK(NIJS:NIJL,13))
END IF

IF (CFLAG_P(9)) BLOCK(NIJS:NIJL,9) = 4.*SQRT(BLOCK(NIJS:NIJL,9))

IF (CFLAG_P(10)) CALL PEAK_PERIOD (FL3, BLOCK(NIJS:NIJL,10))

IF (CFLAG_P(14).AND.CFLAG_P(15)) THEN
   CALL MEAN_DIRECTION (FL3, BLOCK(NIJS:NIJL,14), BLOCK(NIJS:NIJL,15))
ELSE IF (CFLAG_P(14)) THEN
   CALL MEAN_DIRECTION (FL3, THQ=BLOCK(NIJS:NIJL,14))
ELSE IF (CFLAG_P(15)) THEN
   CALL MEAN_DIRECTION (FL3, SPREAD=BLOCK(NIJS:NIJL,15))
END IF
IF (CFLAG_P(14)) BLOCK(NIJS:NIJL,14) = BLOCK(NIJS:NIJL,14)*DEG
IF (CFLAG_P(15)) BLOCK(NIJS:NIJL,15) = BLOCK(NIJS:NIJL,15)*DEG


IF (ANY(CFLAG_P(33:39))) THEN
   CALL KURTOSIS (FL3, DEPTH(NIJS:NIJL), BLOCK(NIJS:NIJL,34),                  &
&                                        BLOCK(NIJS:NIJL,35),                  &
&                                        BLOCK(NIJS:NIJL,33),                  &
&                                        BLOCK(NIJS:NIJL,38),                  &
&                                        BLOCK(NIJS:NIJL,39),                  &
&                                        BLOCK(NIJS:NIJL,36),                  &
&                                        BLOCK(NIJS:NIJL,37))
END IF
IF (CFLAG_P(39)) BLOCK(NIJS:NIJL,39) = DEG * BLOCK(NIJS:NIJL,39)

IF (ITEST.GE.3)   WRITE(IU06,*)                                                &
& '      SUB. COMPUTE_OUTPUT_PARAMETER: INTEGRATED PARAMETERS COMPUTED'


! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. WINDSEA SWELL SEPARATION.                                             !
!        -------------------------                                             !

IF (ANY(CFLAG_P(17:32)).OR.ANY(CFLAG_S(2:NOUT_S)) ) THEN
   CALL SWELL_SEPARATION (FL3, FL)
   IF (ITEST.GE.3) THEN
      WRITE(IU06,*) '      SUB. COMPUTE_OUTPUT_PARAMETER: SWELL /SEA SEPARATION DONE'
   END IF
END IF

END SUBROUTINE COMPUTE_OUTPUT_PARAMETER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SWELL_SEPARATION (FL3, FL1)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SWELL_SEPARATION - COMPUTES THE SWELL SPECTRUM AND INTEGRATED PARAMETER    !
!                      OF SWELL AND WINDSEA.                                   !
!                                                                              !
!     P.LIONELLO     FEBRUARY 87                                               !
!                                                                              !
!     L.ZAMBRESKY    NOVEMBER 87   GKSS/ECMWF   OPTIMIZED SUB.                 !
!     H. GUENTHER    FEBRUARY 2002   GKSS       FT90 AND INTEGRATED PARAMETERS !
!     A. Behrens     December 2003 MSC/GKSS     Message passing                !
!     E. Myklebust   February 2005              MPI parallelization
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO SEPARATE THE SWELL FROM THE WIND INTERACTING SEA AND COMPUTE        !
!       INTEGRATED PARAMETER FOR BOTH SYSTEMS.                                 !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE WAVES WHICH DO NOT INTERACT WITH THE WIND ARE CONSIDERED SWELL.    !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)  :: FL3(:,:,:)        !! BLOCK OF SPECTRA.
REAL,    INTENT(OUT) :: fl1(:,:,:)        !! BLOCK OF SWELL SPECTRA.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: FRIC = 28.
LOGICAL :: SWELL(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)) !! TRUE FOR SWELL ENERGY.

INTEGER  :: K, M
REAL     :: CM
REAL     :: SIS(SIZE(FL3,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. THE SWELL DISTRIBUTION IS COMPUTED.                                   !
!        -----------------------------------                                   !

FL1 = 0.
DO M = 1,ML
   CM = FRIC/C(M)
   DO K=1,KL
      SIS(:) = 1.2*USTAR(:)*COS(TH(K)-UDIR(:))
      WHERE (cm*sis(:)<1.0) fl1(:,k,m) = fl3(:,k,m)
      SWELL(:,k,m) =  cm*sis(:)<1.0
   END DO
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTATION INTEGRATED PARAMETER FOR WINDSEA.                         !
!        ---------------------------------------------                         !

IF (CFLAG_P(17).OR.ANY(CFLAG_P(19:21))) THEN
   CALL TOTAL_ENERGY (FL3, BLOCK(NIJS:NIJL,17), MASK=.NOT.SWELL)
   IF (CFLAG_P(19)) THEN
      CALL FEMEAN (FL3, BLOCK(NIJS:NIJL,17), FM=BLOCK(NIJS:NIJL,19),           &
&                  MASK=.NOT.SWELL)
      BLOCK(NIJS:NIJL,19) = 1./BLOCK(NIJS:NIJL,19)
   END IF
   IF (CFLAG_P(20).AND.CFLAG_P(21)) THEN
      CALL TM1_TM2_PERIODS (FL3, BLOCK(NIJS:NIJL,17), TM1=BLOCK(NIJS:NIJL,20), &
&                           TM2=BLOCK(NIJS:NIJL,21), MASK=.NOT.SWELL)
   ELSE IF (CFLAG_P(20)) THEN
      CALL TM1_TM2_PERIODS (FL3, BLOCK(NIJS:NIJL,17), TM1=BLOCK(NIJS:NIJL,20), &
&                           MASK=.NOT.SWELL)
   ELSE IF (CFLAG_P(21)) THEN
      CALL TM1_TM2_PERIODS (FL3, BLOCK(NIJS:NIJL,17), TM2=BLOCK(NIJS:NIJL,21), &
&                           MASK=.NOT.SWELL)
   END IF
   IF (CFLAG_P(17)) BLOCK(NIJS:NIJL,17) = 4.*SQRT(BLOCK(NIJS:NIJL,17))
END IF

IF (CFLAG_P(18)) CALL PEAK_PERIOD (FL3,BLOCK(NIJS:NIJL,18), MASK=.NOT.SWELL)

IF (CFLAG_P(22).AND.CFLAG_P(23)) THEN
   CALL MEAN_DIRECTION (FL3, THQ=BLOCK(NIJS:NIJL,22),                          &
&                       SPREAD=BLOCK(NIJS:NIJL,23), MASK=.NOT.SWELL)
ELSE IF (CFLAG_P(22)) THEN
   CALL MEAN_DIRECTION (FL3, THQ=BLOCK(NIJS:NIJL,22), MASK=.NOT.SWELL)
ELSE IF (CFLAG_P(23)) THEN
   CALL MEAN_DIRECTION (FL3, SPREAD=BLOCK(NIJS:NIJL,23), MASK=.NOT.SWELL)
END IF
IF (CFLAG_P(22)) BLOCK(NIJS:NIJL,22) = BLOCK(NIJS:NIJL,22)*DEG
IF (CFLAG_P(23)) BLOCK(NIJS:NIJL,23) = BLOCK(NIJS:NIJL,23)*DEG

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTATION INTEGRATED PARAMETER FOR SWELL.                           !
!        -------------------------------------------                           !

IF (CFLAG_P(25).OR.ANY(CFLAG_P(27:29))) THEN
   CALL TOTAL_ENERGY (FL3, BLOCK(NIJS:NIJL,25), MASK=SWELL)
   IF (CFLAG_P(27)) THEN
      CALL FEMEAN (FL3,  BLOCK(NIJS:NIJL,25), FM=BLOCK(NIJS:NIJL,27),          &
&                  MASK=SWELL)
      BLOCK(NIJS:NIJL,27) = 1./BLOCK(NIJS:NIJL,27)
   END IF
   IF (CFLAG_P(28).AND.CFLAG_P(29)) THEN
      CALL TM1_TM2_PERIODS (FL3, BLOCK(NIJS:NIJL,25), TM1=BLOCK(NIJS:NIJL,28), &
&                           TM2=BLOCK(NIJS:NIJL,29), MASK=SWELL)
   ELSE IF (CFLAG_P(28)) THEN
      CALL TM1_TM2_PERIODS (FL3, BLOCK(NIJS:NIJL,25), TM1=BLOCK(NIJS:NIJL,28), &
&                           MASK=SWELL)
   ELSE IF (CFLAG_P(29)) THEN
      CALL TM1_TM2_PERIODS (FL3, BLOCK(NIJS:NIJL,25), TM2=BLOCK(NIJS:NIJL,29), &
&                           MASK=SWELL)
   END IF
   IF (CFLAG_P(25)) BLOCK(NIJS:NIJL,25) = 4.*SQRT( BLOCK(NIJS:NIJL,25))
END IF

IF (CFLAG_P(26)) CALL PEAK_PERIOD (FL3, BLOCK(NIJS:NIJL,26), MASK=SWELL)

IF (CFLAG_P(30).AND.CFLAG_P(31)) THEN
   CALL MEAN_DIRECTION (FL3, THQ=BLOCK(NIJS:NIJL,30),                          &
&                            SPREAD=BLOCK(NIJS:NIJL,31), MASK=SWELL)
ELSE IF (CFLAG_P(30)) THEN
   CALL MEAN_DIRECTION (FL3, THQ=BLOCK(NIJS:NIJL,30), MASK=SWELL)
ELSE IF (CFLAG_P(31)) THEN
   CALL MEAN_DIRECTION (FL3, SPREAD=BLOCK(NIJS:NIJL,31), MASK=SWELL)
END IF
IF (CFLAG_P(30)) BLOCK(NIJS:NIJL,30) = BLOCK(NIJS:NIJL,30)*DEG
IF (CFLAG_P(31)) BLOCK(NIJS:NIJL,31) = BLOCK(NIJS:NIJL,31)*DEG

END SUBROUTINE SWELL_SEPARATION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WRITE_MODEL_OUTPUT (FL3, FL, IU20, IU25)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    WRITE_MODEL_OUTPUT - WRITES THE MODEL OUTPUT.                             !
!                                                                              !
!     H. GUNTHER         GKSS/ECMWF         JUNE 1990                          !
!                                                                              !
!    PURPOSE.                                                                  !
!     --------                                                                 !
!                                                                              !
!       WRITE REQUESTED OUTPUT TO PRINTER AND FILE.                            !
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
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN) :: FL3(:,:,:)  !! BLOCK OF SPECTRA.
REAL,    INTENT(IN) :: FL (:,:,:)  !! BLOCK OF SWELL SPECTRA.
INTEGER, INTENT(IN) :: IU20        !! PARAMETER UNIT NUMBER.
INTEGER, INTENT(IN) :: IU25        !! SPECTRA UNIT NUMBER.

character (len= 12) :: filename
character (len=128) :: ready_fname
integer :: ilen

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. OUTPUT OF SPECTRA AT SELECTED GRID POINTS.                            !
!        ------------------------------------------                            !

IF (CDTSPT.EQ.CDTPRO) THEN
   CALL WRITE_SPECTRA_OUTPUT (FL3, FL, IU25)
   IF (ITEST.GE.3) WRITE(IU06,*)                                               &
&     '      SUB. MODEL_OUTPUT: SUB WRITE_SPECTRA_OUTPUT DONE'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. OUTPUT OF INTEGRATED PARAMETERS.                                      !
!        --------------------------------                                      !

IF (CDTINTT.EQ.CDTPRO) THEN
   IF (ITEST.GE.3) WRITE(IU06,*)                                               &
&    '      SUB. MODEL_OUTPUT: SUB WRITE_INT_PAR_OUTPUT START'
   CALL WRITE_INT_PAR_OUTPUT (IU20)
   IF (ITEST.GE.3) WRITE(IU06,*)                                               &
&    '      SUB. MODEL_OUTPUT: SUB WRITE_INT_PAR_OUTPUT DONE'
!
!==> write ready file for output of integrated parameters
!
   if (ready_outf) then
      call difdate (cdatea, cdtintt, ishift)
      ishift = ishift/3600
      if (ishift<0) ishift = ishift+12
      call create_ready_file_name (ishift, area//'_', filename)
      ready_fname = trim(owpath)//'/'//filename
      ilen = len_trim(ready_fname)
      write (iu06,*) ' +++ ready-file written : '//ready_fname(1:ilen)
      open (iu67,file=ready_fname(1:ilen))
      write (iu67,'(a)') 'ready'
      close (iu67)
   endif
END IF

END SUBROUTINE WRITE_MODEL_OUTPUT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WRITE_INT_PAR_OUTPUT (IU20)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      WRITE_INT_PAR_OUTPUT  - OUTPUT OF INTEGRATED FIELDS.                    !
!                                                                              !
!     P.JANSSEN      KNMI                                                      !
!                                                                              !
!     P.LIONELLO  FEB. 87                                                      !
!                 OUTPUT OF SWELL ENERGY ,MEAN SWELL DIRECTION ,               !
!                 MEAN SWELL FREQUENCY AND MEAN WIND-SEA WAVE                  !
!                 DIRECTION AT ALL ACTIVE GRID POINTS                          !
!                                                                              !
!     H.GUNTHER   GKSS/ECMWF     NOVEMBER 1989                                 !
!     H.GUNTHER   GKSS           FEBRUARY 2002                                 !
!     A.Behrens   MSC            November 2003   (Message passing)             !
!     E. Myklebust               November 2004   MPI parallelization           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       SAVE INTEGRATED PARAMETER IN GRID ARRAYS.                              !
!                                                                              !
!     externals :                                                              !
!     -----------                                                              !
!                                                                              !
!       mpi_gather_block                                                       !
!       mpi_gather_grid                                                        !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN) :: IU20         !! PARAMETER UNIT NUMBER.
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLE.                                                          !
!     ---------------                                                          !

REAL, PARAMETER       :: ZMISS     =  -999.  !! MISSING VALUE (LAND)
REAL, PARAMETER       :: ZMISS_ICE =  -997.  !! MISSING VALUE (ICE)
REAL, PARAMETER       :: ZMISS_DRY =  -998.  !! MISSING VALUE (DRY)
REAL, PARAMETER       :: ZMISS_PRI =  -995.  !! MISSING VALUE for Print

INTEGER                          :: IP, ierr
REAL,ALLOCATABLE, DIMENSION(:,:) :: GRID        !! GRIDDED PARAMETER FIELD.
REAL,ALLOCATABLE, DIMENSION(:)   :: BLOCK_TOTAL !! FULL PARAMETER FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. WRITE INTEGRATED PARAMETER TO FILE HEADER.                            !
!        ------------------------------------------                            !

IF (FFLAG20) THEN
   if (irank==i_out_par) then
      WRITE (IU20) CDTPRO, REAL(NX), REAL(NY), AMOWEP, AMOSOP, AMOEAP, AMONOP, &
&                  cdatea
      WRITE (IU20) NLON_RG, ZDELLO
      WRITE (IU20) FFLAG_P(1:NOUT_P)
   endif
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.  STORE INTEGRATED PARAMETERS OF TOTAL SEA IN GRID ARRAYS.             !
!         --------------------------------------------------------             !

IF (irank==i_out_par) THEN
   ALLOCATE (BLOCK_TOTAL(1:nsea))
   ALLOCATE (GRID(1:NX,1:NY))
END IF

DO IP = 1,NOUT_P
   IF (.NOT. FFLAG_P(IP) .AND. .NOT. PFLAG_P(IP) ) CYCLE

!     2.1 GATHER PARAMETER BLOCKS ON ONE PROCESSOR FOR EACH BLOCK.             !
!         --------------------------------------------------------             !

   if (irank==i_out_par) then
      CALL mpi_gather_block(i_out_par, BLOCK(:,IP), BLOCK_TOTAL)
   else
      CALL mpi_gather_block(i_out_par, BLOCK(:,IP))
   end if
   call mpi_barrier (mpi_comm_world, ierr)

!     2.1 INSERT ICE AND DRY POINT.                                            !
!         -------------------------                                            !

   if (irank==i_out_par) then

      IF ( ICE_RUN .OR. N_DRY.GT.0 ) THEN
         IF (IP.EQ.3 .OR.IP.EQ.4 .OR.IP.GT.8) THEN
            IF (ICE_RUN)    CALL PUT_ICE (BLOCK_TOTAL(:), ZMISS_ICE)
            IF (N_DRY.GT.0) CALL PUT_DRY (BLOCK_TOTAL(:), ZMISS_DRY)
         END IF
         IF (IP.EQ.5) THEN
            IF (ICE_RUN)    CALL PUT_ICE (BLOCK_TOTAL(:), RCHAR)
            IF (N_DRY.GT.0) CALL PUT_DRY (BLOCK_TOTAL(:), RCHAR)
END IF
     END IF

!     2.2 MAKE GRID FIELD.                                                     !
!        -----------------                                                     !

      GRID = UNPACK (BLOCK_TOTAL(:), L_S_MASK, ZMISS)

!     2.3 WRITE OUTPUT.                                                        !
!         -------------                                                        !

      IF (FFLAG_P(IP)) WRITE (IU20) GRID
      IF (PFLAG_P(IP)) THEN
         IF (IP.EQ.5) GRID = MIN(GRID, 999.)
         CALL PRINT_ARRAY (IU06, CDTPRO, TITL_P(IP), GRID,                     &
&                        AMOWEP, AMOSOP, AMOEAP, AMONOP, SCAL_P(IP),ZMISS_PRI)
      END IF
   END IF
   
END DO

IF (ALLOCATED(BLOCK_TOTAL)) DEALLOCATE(BLOCK_TOTAL)
IF (ALLOCATED(GRID)       ) DEALLOCATE(GRID)

END SUBROUTINE WRITE_INT_PAR_OUTPUT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WRITE_SPECTRA_OUTPUT (FL3, FL1, IU25)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   WRITE_SPECTRA_OUTPUT -  OUTPUT OF SPECTRA.                                 !
!                                                                              !
!      S.D.HASSELMANN.                                                         !
!      P.JANSSEN       KNMI         FEBRUARY 1986                              !
!                                                                              !
!      P. LIONELLO                  FEBRUARY 1987                              !
!                  OUTPUT OF SWELL 2-D DISTRIBUTION , SWELL WAVE HEIGHT        !
!                  MEAN SWELL FREQUENCY , MEAN SWELL DIRECTION ,WIND-SE        !
!                  WAVE HEIGTH AND WIND SEA WAVE DIRECTION.                    !
!                                                                              !
!      L. ZAMBRESKY   GKSS/ECMWF    JULY 89                                    !
!                 VARIANCE ENERGY SUMMED OVER FREQUENCY AND DIRECTION.         !
!                                                                              !
!      H. GUNTHER     ECMWF         OCTOBER 1990                               !
!                 EXTENDED HEADER AND OUTPUT OF 1D SPECTRA IN PRINT.           !
!                                                                              !
!     E. Myklebust                     November 2004 MPI parallelization       !
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
!     INTERFACE VARIABLE.                                                      !
!     -------------------                                                      !

REAL,    INTENT(IN) :: FL3(:,:,:)  !! Block OF SPECTRA.
REAL,    INTENT(IN) :: FL1(:,:,:)  !! BLOCK OF SWELL SPECTRA.
INTEGER, INTENT(IN) :: IU25        !! SPECTRA UNIT NUMBER.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

real, allocatable, dimension (:,:,:,:) :: flpts
real, allocatable, dimension (:,:)     :: Block_sp
real, allocatable, dimension (:,:)     :: SPEC

CHARACTER (LEN=40) :: TEXT

INTEGER  :: IJ, nspfld, itag, ierr
INTEGER  :: XLON, XLAT
REAL     :: XANG, XFRE, TH0

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. LOOP OVER OUTPUT POINTS.                                              !
!        ------------------------                                              !

if (noutp>=1) then
   if (any(cflag_s(2:nout_s))) then
      nspfld = 2
   else
      nspfld = 1
   endif
   itag   = 25
   if (irank==i_out_spec) then
      allocate (flpts(kl, ml, noutp, nspfld))
      allocate (block_sp(noutp,NOUT_P))
      allocate (spec(KL,ML))

      call mpi_gather_spp (i_out_spec, itag, nspfld, NOUT_P, fl3, fl1, BLOCK,  &
&                       flpts, Block_sp)
   else
      call mpi_gather_spp (i_out_spec, itag, nspfld, NOUT_P, fl3, fl1, BLOCK)
   endif
   call mpi_barrier (mpi_comm_world, ierr)
endif
        
if (irank==i_out_spec) then
   DO ij = 1,NOUTP
      XLON = OUTLONG(ij)
      XLAT = OUTLAT (ij)

!     1.1 FILE OUTPUT.                                                         !
!         ------------                                                         !

      IF (FFLAG25) THEN
         XANG = REAL(KL)
         XFRE = REAL(ML)
         TH0 = TH(1)*DEG
         WRITE(IU25) XLON, XLAT, CDTPRO, XANG, XFRE, TH0, FR(1), CO
         WRITE(IU25) FFLAG_S(1:NOUT_S)
         WRITE(IU25) Block_sp(IJ,1), Block_sp(IJ,2), Block_sp(IJ,3),           &
&                    Block_sp(IJ,6), Block_sp(IJ,7), Block_sp(IJ,8)

!     1.1.1 TOTAL SPECTRUM.                                                    !
!           ---------------                                                    !

         IF (FFLAG_S(1)) THEN
            WRITE(IU25) Block_sp(IJ,9:15)
            write (iu25) flpts(1:kl,1:ml,ij,1)
         END IF

!     1.1.2 SEA OUTPUT.                                                        !
!           -----------                                                        !

         IF (FFLAG_S(2)) THEN
            WRITE (IU25) Block_sp(IJ,17:23)
            write (iu25) flpts(1:kl,1:ml,ij,1)-flpts(1:kl,1:ml,ij,2)
         END IF

!     1.1.3 SWELL OUTPUT.                                                      !
!           -------------                                                      !

         IF (FFLAG_S(3)) THEN
            WRITE(IU25) Block_sp(IJ,25:31)
            write (iu25) flpts(1:kl,1:ml,ij,2)
         END IF
      END IF

!     1.2 PRINT OUTPUT.                                                        !
!         -------------                                                        !

!     1.1.2 PRINT SPECTRUM.                                                    !
!           ---------------                                                    !

      IF (PFLAG25) THEN
         IF (PFLAG_S(1)) THEN
            spec(1:kl,1:ml) = flpts(1:kl,1:ml,ij,1)
            TEXT(1:40) = TITL_S(1)(1:20)//NAME(ij)
            CALL PRINT_SPECTRUM (IU06, CDTPRO, XLON, XLAT, TEXT, FR, TH, SPEC, &
&                 Block_sp(IJ, 1), Block_sp(IJ, 2), Block_sp(IJ, 3),           &
&                 Block_sp(IJ, 6), Block_sp(IJ, 7), Block_sp(IJ, 8),           &
&                 Block_sp(IJ, 9), Block_sp(IJ,10), Block_sp(IJ,11),           &
&                 Block_sp(IJ,12), Block_sp(IJ,13), Block_sp(IJ,14),           &
&                 Block_sp(IJ,15))
         END IF


!     1.2.2 PRINT OUTPUT SEA SPECTRUM.                                         !
!           --------------------------                                         !

         IF (PFLAG_S(2)) THEN
            spec(1:kl,1:ml) = flpts(1:kl,1:ml,ij,1)-flpts(1:kl,1:ml,ij,2)
            TEXT(1:40) = TITL_S(2)(1:20)//NAME(ij)
            CALL PRINT_SPECTRUM (IU06, CDTPRO, XLON, XLAT, TEXT, FR, TH, SPEC, &
&                 Block_sp(IJ, 1), Block_sp(IJ, 2), Block_sp(IJ, 3),           &
&                 Block_sp(IJ, 6), Block_sp(IJ, 7), Block_sp(IJ, 8),           &
&                 Block_sp(IJ,17), Block_sp(IJ,18), Block_sp(IJ,19),           &
&                 Block_sp(IJ,20), Block_sp(IJ,21), Block_sp(IJ,22),           &
&                 Block_sp(IJ,23))
         END IF


!     1.2.4 PRINT SWELL SPECTRUM.                                              !
!           ---------------------                                              !

         IF (PFLAG_S(3)) THEN
            spec(1:kl,1:ml) = flpts(1:kl,1:ml,ij,2)
            TEXT(1:40) = TITL_S(3)(1:20)//NAME(ij)
            CALL PRINT_SPECTRUM (IU06, CDTPRO, XLON, XLAT, TEXT, FR, TH, SPEC, &
&                 Block_sp(IJ, 1), Block_sp(IJ, 2), Block_sp(IJ, 3),           &
&                 Block_sp(IJ, 6), Block_sp(IJ, 7), Block_sp(IJ, 8),           &
&                 Block_sp(IJ,25), Block_sp(IJ,26), Block_sp(IJ,27),           &
&                 Block_sp(IJ,28), Block_sp(IJ,29), Block_sp(IJ,30),           &
&                 Block_sp(IJ,31))
         END IF
      END IF
   END DO

   IF (ALLOCATED(flpts   )) DEALLOCATE (flpts)
   IF (ALLOCATED(block_sp)) DEALLOCATE(block_sp)
   IF (ALLOCATED(spec    )) DEALLOCATE(spec    )

END IF

END SUBROUTINE WRITE_SPECTRA_OUTPUT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_OUTPUT_MODULE

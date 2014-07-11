MODULE WAM_SOURCE_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL VARIABLES, CONSTANTS AND PROCEDURESN FOR THE      !
!   SOURCE FUNCTION INTEGRATION.                                               !
!                                                                              !
!    JUNE 2005:                                                                !
!       Changes introduced as documented in:                                   !
!       (Jean Bidlot, Peter Janssen and Saleh Abdalla: A revised formulation   !
!       for ocean wave dissipation in CY29R1. ECMWF MEMORANDUM RESEARCH        !
!       DEPARTMENT:April 7, 2005 File: R60.9/JB/0516                           !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_INTERFACE_MODULE, ONLY:  &
&       FEMEAN,                  &  !! COMPUTATION OF MEAN FREQUENCY.
&       TOTAL_ENERGY,            &  !! COMPUTATION OF TOTAL ENERGY.
&       TM1_TM2_PERIODS,         &  !! COMPUTATION OF MEAN WAVENUMBER.
&       WM1_WM2_WAVENUMBER          !! COMPUTATION OF MEAN WAVENUMBER.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE, ONLY: G, PI, ZPI, DEG, XEPS, XNLEV, XINVEPS, ALPHA,    &
&                             BETAMAX, XKAPPA, ZALP
USE WAM_FRE_DIR_MODULE, ONLY: KL, ML, FR, CO, TH, DELTH, COSTH, SINTH, DFIM,   &
&                             C, TFAK, FMIN
USE WAM_TIMOPT_MODULE,  ONLY: IDELT, SHALLOW_RUN, WAVE_BREAKING_RUN,           &
&                             PHILLIPS_RUN
USE WAM_FILE_MODULE,    ONLY: IU06, ITEST

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

! ---------------------------------------------------------------------------- !
!
!    1. INDICES AND WEIGHTS USED TO COMPUTE THE NONLINEAR TRANSFER RATE.
!       ----------------------------------------------------------------

INTEGER,ALLOCATABLE, DIMENSION(:)   :: IKP    !! FREQUENCY INDEX FOR STORING 
                                              !! ENERGY TRANSFER INCREMENTS
                                              !! INTO BINS, WAVE NO. 3.
INTEGER,ALLOCATABLE, DIMENSION(:)   :: IKP1   !! IKP+1. 
INTEGER,ALLOCATABLE, DIMENSION(:)   :: IKM    !! FREQUENCY INDEX FOR STORING
                                              !! ENERGY TRANSFER INCREMENTS
                                              !! INTO BINS, WAVE NO. 4. 
INTEGER,ALLOCATABLE, DIMENSION(:)   :: IKM1   !! IKM+1
INTEGER,ALLOCATABLE, DIMENSION(:,:) :: K1W    !! ANGULAR INDEX FOR STORING
                                              !! ENERGY TRANSFER INCREMENTS
                                              !! INTO BINS, WAVE NO. 3.
INTEGER,ALLOCATABLE, DIMENSION(:,:) :: K2W    !! ANGULAR INDEX FOR STORING
                                              !! ENERGY TRANSFER INCREMENTS
                                              !! INTO BINS, WAVE NO. 4.
INTEGER,ALLOCATABLE, DIMENSION(:,:) :: K11W   !! K1W(.,1)-1, K1W(.,2)+1.
INTEGER,ALLOCATABLE, DIMENSION(:,:) :: K21W   !! K2W(.,1)+1, K2W(.,2)-1.

REAL,   ALLOCATABLE, DIMENSION(:)   :: AF11   !! WEIGHTS FOR APPROXIMATION
                                              !! OF NONL. TRANSFER (ONE 
                                              !! TERM ONLY SET TO 3000). 
                                              !! MULTIPLIED BY FREQ. **11.
REAL,   ALLOCATABLE, DIMENSION(:)   :: FKLAP  !! WEIGHT IN FREQUENCY GRID  
                                              !! FOR INTERPOLATION, WAVE  
                                              !! NO. 3 ("1+LAMBDA" TERM)
REAL,   ALLOCATABLE, DIMENSION(:)   :: FKLAP1 !! 1-FKLAP.
REAL,   ALLOCATABLE, DIMENSION(:)   :: FKLAM  !! WEIGHT IN FREQUENCY GRID  
                                              !! FOR INTERPOLATION, WAVE 
                                              !! NO. 4 ("1-LAMBDA" TERM).
REAL,   ALLOCATABLE, DIMENSION(:)   :: FKLAM1 !! 

REAL    :: ACL1        !! WEIGHT IN ANGULAR GRID FOR INTERPOLATION,
                          !! WAVE NO. 3 ("1+LAMBDA" TERM).
REAL    :: ACL2        !! WEIGHT IN ANGULAR GRID FOR INTERPOLATION,
                          !! WAVE NO. 4 ("1-LAMBDA" TERM).
REAL    :: CL11        !! 1.-ACL1.
REAL    :: CL21        !! 1.-ACL2.
REAL    :: DAL1        !! 1./ACL1.
REAL    :: DAL2        !! 1./ACL2.
REAL    :: FRH(30)     !! TAIL FREQUENCY RATIO **5

! ---------------------------------------------------------------------------- !
!
!    2. TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
!       --------------------------------------------

INTEGER, PARAMETER :: ITAUMAX = 200 !! TABLE DIMENSION For WAVE STRESS.
INTEGER, PARAMETER :: JUMAX   = 200 !! TABLE DIMENSION FOR U10.
INTEGER, PARAMETER :: IUSTAR  = 100 !! TABLE DIMENSION FOR USTAR.
INTEGER, PARAMETER :: IALPHA  = 200 !! TABLE DIMENSION FOR ALPHA.
 
REAL,    PARAMETER :: UMAX = 50.   !! MAXIMUM WIND SPEED IN STRESS TABLE.
REAL,    PARAMETER :: USTARM = 5.  !! MAXIMUM FRICTION VELOCITY IN STRESS TABLE.

REAL :: TAUT(0:ITAUMAX,0:JUMAX)     !! STRESS TABLE.
REAL :: DELTAUW                     !! WAVE STRESS INCREMENT.
REAL :: DELU                        !! WIND INCREMENT.
REAL :: TAUHFT(0:IUSTAR,0:IALPHA)   !! HIGH FREQUENCY STRESS TABLE.
REAL :: DELUST                      !! USTAR INCREMENT.
REAL :: DELALP                      !! ALPHA INCREMENT.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE IMPLSCH                  !! IMPLICIT INTEGRATION OF SOURCE FUNCTIONS.
   MODULE PROCEDURE IMPLSCH
END INTERFACE
PUBLIC IMPLSCH

INTERFACE PREPARE_SOURCE           !! PRE-COMPUTES MODULE CONSTANTS.
   MODULE PROCEDURE PREPARE_SOURCE
END INTERFACE
PUBLIC PREPARE_SOURCE

INTERFACE PRINT_SOURCE_STATUS      !! PRINTS THE STATUS OF THIS MODULE.
   MODULE PROCEDURE PRINT_SOURCE_STATUS
END INTERFACE
PUBLIC PRINT_SOURCE_STATUS

INTERFACE AIRSEA              !! DETERMINE TOTAL STRESS IN SURFACE LAYER.
   MODULE PROCEDURE AIRSEA
END INTERFACE
PUBLIC AIRSEA

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE NLWEIGT             !! COMPUTATION OF INDEX ARRAYS AND WEIGHTS
   MODULE PROCEDURE NLWEIGT   !! FOR THE NONLINEAR TRANSFER RATE.
END INTERFACE
PRIVATE NLWEIGT

INTERFACE SBOTTOM             !! COMPUTATION OF BOTTOM FRICTION.
   MODULE PROCEDURE SBOTTOM
END INTERFACE
PRIVATE SBOTTOM

INTERFACE SDISSIP             !! COMPUTATION OF DISSIPATION SOURCE FUNCTION.
   MODULE PROCEDURE SDISSIP
END INTERFACE
PRIVATE SDISSIP

INTERFACE SFBRK               !! COMPUTATION OF DEPTH_INDUCED WAVE BREAKING.
   MODULE PROCEDURE SFBRK
END INTERFACE
PRIVATE SFBRK

INTERFACE CMPQB               !! COMPUTATION OF OF WAVE BREAKING SOUREC FKT.
   MODULE PROCEDURE CMPQB
END INTERFACE
PRIVATE CMPQB

INTERFACE SINPUT              !! COMPUTATION OF INPUT SOURCE FUNCTION.
   MODULE PROCEDURE SINPUT
END INTERFACE
PRIVATE SINPUT

INTERFACE SNONLIN             !! COMPUTATION OF NONLINEAR TRANSFER RATE.
   MODULE PROCEDURE SNONLIN
END INTERFACE
PRIVATE SNONLIN
 
INTERFACE SOURCE_PHILLIPS     !! COMPUTATION OF PHILLIP'S INPUT SOURCE FUNCTION.
   MODULE PROCEDURE SOURCE_PHILLIPS
END INTERFACE
PRIVATE SOURCE_PHILLIPS

INTERFACE STRESS              !! COMPUTATION OF TOTAL STRESS TABLE.
   MODULE PROCEDURE STRESS
END INTERFACE
PRIVATE STRESS

INTERFACE STRESSO             !! COMPUTATION OF WAVE STRESS.
   MODULE PROCEDURE STRESSO
END INTERFACE
PRIVATE STRESSO

INTERFACE TAUHF               !! COMPUTATION OF HIGH-FREQUENCY STRESS TABLE.
   MODULE PROCEDURE TAUHF
END INTERFACE
PRIVATE TAUHF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE IMPLSCH (FL3, U10, UDIR, TAUW, USTAR, Z0, DEPTH, INDEP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   IMPLSCH - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE FUNCTIONS.        !
!                                                                              !
!     S.D.HASSELMANN.  MPI                                                     !
!     H. GUENTHER AND L. ZAMBRESKY  OPTIMIZATION PERFORMED.                    !
!     H. GUENTHER      GKSS/ECMWF   OCTOBER 1989  NEW WIND FIELD               !
!                                                 INTERFACE AND                !
!                                                 TIME COUNTING                !
!     P.A.E.M. JANSSEN KNMI         AUGUST  1990  COUPLED MODEL                !
!     H. GUENTHER      GKSS/ECMWF   JUNE    1991  NEW SEPARATION OF            !
!                                                  DIAG- AND PROGNOSTIC        !
!                                                  PART OF SPECTRUM.           !
!     P.A.E.M. JANSSEN ECMWF        FEBRUARY 1995  ADD MINIMUM VALUE (FMIN).   !
!     H. GUENTHER      GKSS         JANUARY  2002  FT90                        !
!     ROOP LALBEHARRY  MSC/ARMN     APRIL    2003
!                      PHILLIPS SOURCE
!     ERIK MYKLEBUST                FEBRUARY 2005 OPTIMIZATION
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS             !
!       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.             !
!       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH                  !
!       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.                                 !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS FN+1=FN+DELT*(SN+SN+1)/2.,  !
!       WHERE SN IS THE TOTAL SOURCE FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF   !
!       - ONLY THE DIAGONAL TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE COMPUTED, !
!         THE NONDIAGONAL TERMS ARE NEGLIGIBLE.                                !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",                !
!       30/6/85 (UNPUBLISHED NOTE)                                             !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(INOUT) :: FL3(:,:,:)     !! FREQUENCY SPECTRUM.
REAL,    INTENT(IN)    :: U10   (:)      !! WIND SPEED [M/S].
REAL,    INTENT(IN)    :: UDIR  (:)      !! WIND DIRECTION [RAD].
REAL,    INTENT(INOUT) :: TAUW  (:)      !! WAVE STRESS IN (M/S)**2 
REAL,    INTENT(OUT)   :: USTAR (:)      !! FRICTION VELOCITY [M/S].
REAL,    INTENT(INOUT) :: Z0    (:)      !! ROUGHNESS LENGTH [M].
REAL,    INTENT(IN)    :: DEPTH (:)      !! WATER DEPTH [M].
INTEGER, INTENT(IN)    :: INDEP (:)      !! DEPTH TABLE INDEX.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !

INTEGER :: K, M, KL, ML, IJ
REAL    :: DELT, FPMH, CM

REAL    :: FL(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3))  !! DIAGONAL MATRIX OF
                                                    !! FUNCTIONAL DERIVATIVE.
REAL    :: SL(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3))  !! TOTAL SOURCE FUNCTION.

INTEGER :: MIJ(SIZE(FL3,1))  !! LAST FREQUENCY INDEX OF PROGNOSTIC PART.
REAL    :: EMEAN(SIZE(FL3,1))  !! TOTAL ENERGY
REAL    :: FMEAN(SIZE(FL3,1))  !! MEAN FREQUENCY
REAL    :: AKMEAN(SIZE(FL3,1)) !! MEAN WAVENUMBER BASED ON SQRT(1/K)-MOMENT
REAL    :: GADIAG(SIZE(FL3,1)), TEMP(SIZE(FL3,1),SIZE(FL3,3))
REAL    :: DELFL(SIZE(FL3,3))
LOGICAL :: MASK(SIZE(FL3,1),SIZE(FL3,3))
LOGICAL :: LLWS(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3))

REAL    :: F1MEAN(SIZE(FL3,1))    !! MEAN FREQUENCY BASED ON F-MOMENT
REAL    :: XKMEAN(SIZE(FL3,1))    !! MEAN WAVENUMBER BASED ON SQRT(K)-MOMENT
REAL    :: EMEANWS(SIZE(FL3,1))   !! TOTAL WINDSEA ENERGY
REAL    :: FMEANWS(SIZE(FL3,1))   !! MEAN WINDSEA FREQUENCY
    
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CALCULATE ROUGHNESS LENGTH AND FRICTION VELOCITIES.                   !
!        ---------------------------------------------------                   !

KL = SIZE(FL3,2)
ML = SIZE(FL3,3)

CALL AIRSEA (U10, TAUW, USTAR, Z0)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN PARAMETERS.                                              !
!        ------------------------                                              !

CALL TOTAL_ENERGY (FL3, EMEAN)
CALL FEMEAN (FL3, EMEAN, FMEAN)
CALL TM1_TM2_PERIODS (FL3, EMEAN, TM1=F1MEAN)
F1MEAN = 1./F1MEAN
IF (SHALLOW_RUN) THEN
   CALL WM1_WM2_WAVENUMBER (FL3, EMEAN, WM1=AKMEAN, WM2=XKMEAN, IN=INDEP)
ELSE
   CALL WM1_WM2_WAVENUMBER (FL3, EMEAN, WM1=AKMEAN, WM2=XKMEAN)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTATION OF SOURCE FUNCTIONS AND DERIVATIVES.                      !
!        ------------------------------------------------                      !

SL = 0.     !! SOURCE FUNCTION
FL = 0.     !! DERIVATIVE ARRAY

CALL SINPUT  (FL3, SL, FL, USTAR, UDIR, Z0, INDEP, LLWS)
IF (PHILLIPS_RUN) THEN
   call source_phillips (sl, ustar, udir, depth, indep)
END IF
CALL STRESSO (FL3, SL, USTAR, UDIR, Z0, TAUW)

CALL SNONLIN (FL3, SL, FL, DEPTH, AKMEAN)
CALL SDISSIP (FL3, SL, FL, EMEAN, F1MEAN, XKMEAN, INDEP)

IF (SHALLOW_RUN) THEN
   CALL SBOTTOM (FL3, SL, FL, DEPTH, INDEP)
   IF (WAVE_BREAKING_RUN) CALL SFBRK (FL3, SL, FL, EMEAN, FMEAN, DEPTH)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTATION OF NEW SPECTRA.                                           !
!        ---------------------------                                           !
!                                                                              !
!       INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE             !
!       FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.                    !

!     4.1 REEVALUATE WAVE STRESS AND MEAN WINDSEA PARAMETER.                   !
!        --------------------------------------------------                    !

CALL AIRSEA (U10, TAUW, USTAR, Z0)

DO K=1,KL
   TEMP(:,1) = USTAR*COS(TH(K)-UDIR)
   DO M=1,ML
      CM = 28./C(M)
      WHERE (CM*TEMP(:,1).GE.1.)  LLWS(:,K,M) =.TRUE.
   END DO
END DO
CALL TOTAL_ENERGY (FL3, EMEANWS, LLWS)
CALL FEMEAN (FL3, EMEANWS, FMEANWS, LLWS)

!     4.2 INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE           !
!         FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.                  !
!         ----------------------------------------------------------           !

DELT = IDELT
DELFL = 3.0E-07*G/FR**4*DELT
DO M=1,ML
   TEMP(:,2) = USTAR*DELFL(M)*MAX(FMEANWS,FMEAN)
   DO K=1,KL
      TEMP(:,1) = DELT*SL(:,K,M)/ MAX(1.,1.-DELT*FL(:,K,M))
      TEMP(:,3) = MIN(ABS(TEMP(:,1)),TEMP(:,2))
      FL3(:,K,M) = FL3(:,K,M)+SIGN(TEMP(:,3),TEMP(:,1))
   END DO
END DO
 
FL3 = MAX (FL3, FMIN)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.                 !
!        -----------------------------------------------------                 !
!                                                                              !
!    5.1 COMPUTE MEAN PARAMETERS.                                              !
!        ------------------------                                              !

CALL TOTAL_ENERGY (FL3, EMEAN)
CALL FEMEAN (FL3, EMEAN, FMEAN)
CALL TOTAL_ENERGY (FL3, EMEANWS, LLWS)
CALL FEMEAN (FL3, EMEANWS, FMEANWS, LLWS)

!     5.2 COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.         !
!         FREQUENCIES LE 2.5*MAX(FMEANWS , FMEAN).                             !
!         ------------------------------------------------------------         !

FPMH = 2.5/FR(1)
MIJ = INT(LOG10(MAX(FMEANWS,FMEAN)*FPMH)*24.1589)+1
MIJ = MIN(MIJ,ML)

!     5.3 COMPUTE TAIL ENERGY RATIOS.                                          !
!         ---------------------------                                          !

DELFL = (1./FR)**5.
GADIAG = FR(MIJ)**5.
DO M = 1,ML
   TEMP(:,M) = GADIAG*DELFL(M)
END DO

!     5.4 MERGE TAIL INTO SPECTRA.                                             !
!         ------------------------                                             !

MASK = .TRUE.
DO IJ = 1,SIZE(FL3,1)
   MASK(IJ,1:MIJ(IJ)) = .FALSE.
END DO

DO K =1 ,KL
   DO IJ=1,SIZE(FL3,1)
      GADIAG(IJ) = FL3(IJ,K,MIJ(IJ))
   END DO
   DO M = 1,ML
      WHERE (MASK(:,M)) FL3(:,K,M) = GADIAG * TEMP(:,M)
   END DO
END DO
FL3 = MAX (FL3, FMIN)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. REEVALUATE WIND INPUT SOURCE TERM AND WAVE STRESS.                    !
!        --------------------------------------------------                    !

SL = 0.     !! SOURCE FUNCTION
FL = 0.     !! DERIVATIVE ARRAY
CALL SINPUT  (FL3, SL, FL, USTAR, UDIR, Z0, INDEP, LLWS)
CALL STRESSO (FL3, SL, USTAR, UDIR, Z0, TAUW)
CALL AIRSEA  (U10, TAUW, USTAR, Z0)

END SUBROUTINE IMPLSCH

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_SOURCE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PREPARE_SOURCE - ROUTINE TO PREPARE WAM SOURCE MODULE.                     !
!                                                                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO COMPUTE ALL VARIABLES IN WAM SOURCE MODULE.                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       DONE BY CALLS TO MANY SUBS.                                            !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. WEIGHT OF NON-LINEAR INTERACTION.                                     !
!        ---------------------------------                                     !

CALL NLWEIGT
IF (ITEST.GE.2) WRITE (IU06,*) '    SUB. PREPARE_SOURCE: NLWEIGT DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. STRESS TABLES.                                                        !
!        --------------                                                        !

CALL STRESS
IF (ITEST.GE.2) WRITE (IU06,*) '    SUB. PREPARE_SOURCE: STRESS DONE'
CALL TAUHF (FR(ML))
IF (ITEST.GE.2) WRITE (IU06,*) '    SUB. PREPARE_SOURCE: TAUHF DONE'

END SUBROUTINE PREPARE_SOURCE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_SOURCE_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PRINT_SOURCE_STATUS - PRINT STATUS OF WAM_SOURCE_MODULE.                   !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       MAKE A PRINTER OUTPUT OF THE DATA STORED IN WAN_SOURCE MODULE.         !
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
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER      :: K, M, KH

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NON LINEAR INTERACTION PARAMETERS.                                    !
!        ----------------------------------                                    !

IF (.NOT. ALLOCATED(IKP)) RETURN

WRITE(IU06,'(/,'' -------------------------------------------------'')')
WRITE(IU06,*)'        NON LINEAR INTERACTION PARAMETERS' 
WRITE(IU06,'(  '' -------------------------------------------------'')')
WRITE(IU06,'(/,''  FREQUENCY ARRAYS'')')
WRITE(IU06,'(''     ACL1       ACL2       CL11       CL21   '',                &
&            ''    DAL1       DAL2'')')
WRITE(IU06,'(1X,6F11.8)') ACL1, ACL2, CL11, CL21, DAL1, DAL2
WRITE(IU06,*) ' '
WRITE(IU06,'(''  M   IKP IKP1  IKM IKM1   FKLAP       FKLAP1 '',               &
&            ''   FKLAM       FKLAM1     AF11'')')

DO M = 1,SIZE(IKP)
   WRITE(IU06,'(1X,I2,4I5,4F11.8,E11.3)') M, IKP(M), IKP1(M), IKM(M), IKM1(M), &
&                           FKLAP(M), FKLAP1(M), FKLAM(M), FKLAM1(M), AF11(M)
END DO

WRITE(IU06,'(/,''  ANGULAR ARRAYS'')')
WRITE(IU06,'(''   |--------KH = 1----------||--------KH = 2----------|'')')
WRITE(IU06,'(''  K   K1W   K2W  K11W  K21W   K1W   K2W  K11W  K21W'')')
DO  K= 1,SIZE(K1W,1)
   WRITE(IU06,'(1X,I2,8I6)') K,(K1W(K,KH), K2W(K,KH), K11W(K,KH),              &
&                            K21W(K,KH),KH=1,2)
END DO

WRITE(IU06,'(/,''  TAIL ARRAY FRH'')')
WRITE(IU06,'(1X,8F10.7)') FRH(1:30)

! ---------------------------------------------------------------------------- !
!
!    2. TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
!       --------------------------------------------

WRITE(IU06,'(/,'' -------------------------------------------------'')')
WRITE(IU06,*)'       TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS' 
WRITE(IU06,'(  '' -------------------------------------------------'')')
WRITE(IU06,*) ' '
WRITE(IU06,*) ' WAVE STRESS TABLE LENGTH .: ', ITAUMAX
WRITE(IU06,*) ' WAVE STRESS INCREMENT ....: ', DELTAUW 
WRITE(IU06,*) ' WIND TABLE LENGTH ........: ', JUMAX
WRITE(IU06,*) ' WIND INCREMENT ...........: ', DELU 
WRITE(IU06,*) ' USTAR TABLE LENGTH .......: ', IUSTAR
WRITE(IU06,*) ' USTAR INCREMENT ..........: ', DELUST 
WRITE(IU06,*) ' ALPHA TABLE LENGTH .......: ', IALPHA
WRITE(IU06,*) ' ALPHA INCREMENT ..........: ', DELALP

END SUBROUTINE PRINT_SOURCE_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE AIRSEA (U10, TAUW, US, Z0)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   AIRSEA - DETERMINE TOTAL STRESS IN SURFACE LAYER.                          !
!                                                                              !
!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE TOTAL STRESS.                                                  !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       USE TABLE TAUT(TAUW,U) AND LINEAR INTERPOLATION.                       !
!                                                                              !
!     REFERENCE.                                                               !
!     ---------                                                                !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: U10(:)   !! WINDSPEEDS U10.
REAL,    INTENT(IN)    :: TAUW(:)  !! WAVE STRESSES.
REAL,    INTENT(OUT)   :: US(:)    !! SURFACE STRESSES.
REAL,    INTENT(INOUT) :: Z0(:)    !! ROUGHNESS LENGTH.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: IJ, I, J
REAL     :: XI, DELI1, DELI2, XJ, DELJ1, DELJ2

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DETERMINE TOTAL STRESS FROM TABLE.                                    !
!        ----------------------------------                                    !

DO IJ = 1,SIZE(U10)
   XI      = SQRT(TAUW(IJ))/DELTAUW
   I       = MIN ( ITAUMAX-1, INT(XI) )
   DELI1   = MIN(1.,XI - REAL(I))
   DELI2   = 1. - DELI1
   XJ      = U10(IJ)/DELU
   J       = MIN ( JUMAX-1, INT(XJ) )
   DELJ1   = MIN(1.,XJ - REAL(J))
   DELJ2   = 1. - DELJ1
   US(IJ)  = (TAUT(I,J  )*DELI2 + TAUT(I+1,J  )*DELI1)*DELJ2                   &
&          + (TAUT(I,J+1)*DELI2 + TAUT(I+1,J+1)*DELI1)*DELJ1
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. DETERMINE ROUGHNESS LENGTH.                                           !
!        ---------------------------                                           !

where (us==0.)
   z0 = 100.
elsewhere
   z0 = MIN(U10/US,100.0)
end where
Z0  = XNLEV*EXP(-XKAPPA*Z0)

END SUBROUTINE AIRSEA

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE NLWEIGT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   NLWEIGT - COMPUTATION OF INDEX ARRAYS AND WEIGHTS FOR THE COMPUTATION OF   !
!             THE NONLINEAR TRANSFER RATE.                                     !
!                                                                              !
!     SUSANNE HASSELMANN JUNE 86.                                              !
!                                                                              !
!     H. GUNTHER   ECMWF/GKSS  DECEMBER 90 - CYCLE_4 MODIFICATIONS.            !
!                                            4 FREQUENCIES ADDED.              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF PARAMETERS USED IN DISCRETE INTERACTION                 !
!       PARAMETERIZATION OF NONLINEAR TRANSFER.                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     INTERNALS.                                                               !
!     ----------                                                               !
!                                                                              !
!       JAFU      - FUNCTION FOR COMPUTATION OF ANGULAR INDICES OF K(F,THET).  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       S. HASSELMANN AND K. HASSELMANN, JPO, 1985 B.                          !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

!     *PARAMETER*  FOR DISCRETE APPROXIMATION OF NONLINEAR TRANSFER.           !

REAL, PARAMETER :: ALAMD   = 0.25     !! LAMBDA
REAL, PARAMETER :: CON     = 3000.    !! WEIGHT FOR DISCRETE APPROXIMATION OF
                                      !! NONLINEAR TRANSFER
REAL, PARAMETER :: DELPHI1 = -11.48   !!
REAL, PARAMETER :: DELPHI2 = 33.56    !!

INTEGER :: KLP1, IC, KH, KLH, K, KS, ICL1, ICL2, ISG, K1, K11, K2, K21
INTEGER :: M, IKN, I, IE

REAL    :: DELTHA, CL1, CL2, AL11, AL12, CH, CL1H, CL2H
REAL    :: F1P1, FRG, FLP, FLM, FKP, FKM
REAL    :: FRLON(2*ML+2)
INTEGER :: JA1(KL,2), JA2(KL,2)
INTEGER :: JAFU
EXTERNAL   JAFU

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. ALLOCATE ARRAYS.                                                      !
!        ----------------                                                      !

IF (.NOT.ALLOCATED (IKP ))  ALLOCATE (IKP (ML+4))
IF (.NOT.ALLOCATED (IKP1))  ALLOCATE (IKP1(ML+4))
IF (.NOT.ALLOCATED (IKM ))  ALLOCATE (IKM (ML+4))
IF (.NOT.ALLOCATED (IKM1))  ALLOCATE (IKM1(ML+4))
IF (.NOT.ALLOCATED (K1W ))  ALLOCATE (K1W (KL,2))
IF (.NOT.ALLOCATED (K2W ))  ALLOCATE (K2W (KL,2))
IF (.NOT.ALLOCATED (K11W))  ALLOCATE (K11W(KL,2))
IF (.NOT.ALLOCATED (K21W))  ALLOCATE (K21W(KL,2))
IF (.NOT.ALLOCATED (AF11))  ALLOCATE (AF11(ML+4))
IF (.NOT.ALLOCATED (FKLAP ))  ALLOCATE (FKLAP(ML+4))
IF (.NOT.ALLOCATED (FKLAP1))  ALLOCATE (FKLAP1(ML+4))
IF (.NOT.ALLOCATED (FKLAM ))  ALLOCATE (FKLAM(ML+4))
IF (.NOT.ALLOCATED (FKLAM1))  ALLOCATE (FKLAM1(ML+4))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTATION FOR ANGULAR GRID.                                         !
!        -----------------------------                                         !

DELTHA = DELTH*DEG
CL1 = DELPHI1/DELTHA
CL2 = DELPHI2/DELTHA

!     1.1 COMPUTATION OF INDICES OF ANGULAR CELL.                              !
!         ---------------------------------------                              !

KLP1 = KL+1
IC = 1

DO KH = 1,2
   KLH = KL
   IF (KH.EQ.2) KLH=KLP1
   DO K = 1,KLH
      KS = K
      IF (KH.GT.1) KS=KLP1-K+1
      IF (KS.GT.KL) CYCLE
      CH = IC*CL1
      JA1(KS,KH) = JAFU(CH,K,KL)
      CH = IC*CL2
      JA2(KS,KH) = JAFU(CH,K,KL)
   END DO
   IC = -1
END DO

!     1.2 COMPUTATION OF ANGULAR WEIGHTS.                                      !
!         -------------------------------                                      !

ICL1 = CL1
CL1  = CL1-ICL1
ICL2 = CL2
CL2  = CL2-ICL2
ACL1 = ABS(CL1)
ACL2 = ABS(CL2)
CL11 = 1.-ACL1
CL21 = 1.-ACL2
AL11 = (1.+ALAMD)**4
AL12 = (1.-ALAMD)**4
DAL1 = 1./AL11
DAL2 = 1./AL12

!     1.3 COMPUTATION OF ANGULAR INDICES.                                      !
!         -------------------------------                                      !

ISG = 1
DO KH = 1,2
   CL1H = ISG*CL1
   CL2H = ISG*CL2
   DO K = 1,KL
      KS = K
      IF (KH.EQ.2) KS = KL-K+2
      IF(K.EQ.1) KS = 1
      K1 = JA1(K,KH)
      K1W(KS,KH) = K1
      IF (CL1H.LT.0.) THEN
         K11 = K1-1
         IF (K11.LT.1) K11 = KL
      ELSE
         K11 = K1+1
         IF (K11.GT.KL) K11 = 1
      END IF
      K11W(KS,KH) = K11
      K2 = JA2(K,KH)
      K2W(KS,KH) = K2
      IF (CL2H.LT.0) THEN
         K21 = K2-1
         IF(K21.LT.1) K21 = KL
      ELSE
         K21 = K2+1
         IF (K21.GT.KL) K21 = 1
      END IF
      K21W(KS,KH) = K21
   END DO
   ISG = -1
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTATION FOR FREQUENCY GRID.                                       !
!        -------------------------------                                       !

FRLON(1:ML) = FR(1:ML)

DO M = ML+1,2*ML+2
   FRLON(M) = CO*FRLON(M-1)
END DO

F1P1 = ALOG10(CO)
DO M = 1,ML+4
   FRG = FRLON(M)
   AF11(M) = CON * FRG**11
   FLP = FRG*(1.+ALAMD)
   FLM = FRG*(1.-ALAMD)
   IKN = INT(ALOG10(1.+ALAMD)/F1P1+.000001)
   IKN = M+IKN
   IKP(M) = IKN
   FKP = FRLON(IKP(M))
   IKP1(M) = IKP(M)+1
   FKLAP(M) = (FLP-FKP)/(FRLON(IKP1(M))-FKP)
   FKLAP1(M) = 1.-FKLAP(M)
   IF (FRLON(1).GE.FLM) THEN
      IKM(M) = 1
      IKM1(M) = 1
      FKLAM(M) = 0.
      FKLAM1(M) = 0.
   ELSE
      IKN = INT(ALOG10(1.-ALAMD)/F1P1+.0000001)
      IKN = M+IKN-1
      IF (IKN.LT.1) IKN = 1
      IKM(M) = IKN
      FKM = FRLON(IKM(M))
      IKM1(M) = IKM(M)+1
      FKLAM(M) = (FLM-FKM)/(FRLON(IKM1(M))-FKM)
      FKLAM1(M) = 1.-FKLAM(M)
   END IF
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTE TAIL FREQUENCY RATIOS.                                        !
!        ------------------------------                                        !

IE = MIN(30,ML+3)
DO I = 1,IE
   M = ML+I-1
   FRH(I) = (FRLON(ML)/FRLON(M))**5
END DO

! ---------------------------------------------------------------------------- !

RETURN

END SUBROUTINE NLWEIGT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SBOTTOM (F, SL, FL, DEPTH, INDEP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SBOTTOM - COMPUTATION OF BOTTOM FRICTION.                                  !
!                                                                              !
!     G.J.KOMEN AND Q.D.GAO                                                    !
!     OPTIMIZED BY L.F. ZAMBRESKY                                              !
!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90                            !
!     E. MYKLEBUST        FEBRUARY 2005       OPTIMIZATION                     !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF BOTTOM FRICTION DISSIPATION                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       SEE REFERENCES.                                                        !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)                 !
!       BOUWS AND KOMEN, JPO 13(1983)1653-1658                                 !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL, INTENT(IN)     :: F (:, :, :)  !! SPECTRUM.
REAL, INTENT(INOUT)  :: SL(:, :, :)  !! TOTAL SOURCE FUNCTION ARRAY
REAL, INTENT(INOUT)  :: FL(:, :, :)  !! DIAGONAL MATRIX OF FUNC. DERIVATIVE.
REAL, INTENT(IN)     :: DEPTH(:)     !! WATER DEPTH
INTEGER, INTENT(IN)  :: INDEP(:)     !! DEPTH TABLE INDEX.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER  :: CONST = -2.0*0.038/G

INTEGER :: M, K
REAL    :: WAV(SIZE(F,1))
REAL    :: SBO(SIZE(F,1))

! ---------------------------------------------------------------------------- !

FRE: DO M = 1,SIZE(F,3)
    WAV = TFAK(INDEP,M)
    SBO = MIN (2.* DEPTH*WAV ,50.)
    SBO = CONST*WAV/SINH(SBO)

    DIR: DO K = 1,SIZE(F,2)
       SL(:,K,M) = SL(:,K,M) + SBO*F(:,K,M)
       FL(:,K,M) = FL(:,K,M) + SBO
   END DO DIR
END DO FRE

END SUBROUTINE SBOTTOM

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SDISSIP (F, SL, FL, EMEAN, FMEAN, AKMEAN, INDEP)

! ---------------------------------------------------------------------------- !
!                                                                              !
! *** *SDISSIP* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.                  !
!                                                                              !
!     S.D.HASSELMANN.                                                          !
!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN                        !
!     OPTIMIZATION : L. ZAMBRESKY                                              !
!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90                            !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO          !
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE          !
!       OF DISSIPATION SOURCE FUNCTION.                                        !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       SEE REFERENCES.                                                        !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE             !
!          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL, INTENT(IN)     :: F (:, :, :)  !! SPECTRUM.
REAL, INTENT(INOUT)  :: SL(:, :, :)  !! TOTAL SOURCE FUNCTION ARRAY
REAL, INTENT(INOUT)  :: FL(:, :, :)  !! DIAGONAL MATRIX OF FUNCTIONAL
REAL, INTENT(IN)     :: EMEAN (:)    !! TOTAL ENERGY
REAL, INTENT(IN)     :: FMEAN (:)    !! MEAN FREQUENCY BASED ON 1. MOMENT
REAL, INTENT(IN)     :: AKMEAN(:)    !! MEAN WAVE NUMBER BASED ON SQRT(K) MOMENT
INTEGER, INTENT(IN)  :: INDEP (:)    !! DEPTH TABLE INDEX

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: CDIS = 2.1          !! DISSIPATION CONSTANT
REAL, PARAMETER :: CONSD = -CDIS*ZPI**9/G**4
REAL, PARAMETER :: CONSS = -CDIS*ZPI
REAL, PARAMETER :: DELTA = 0.6         !! WEIGHT LINEAR, QUADRATIC PART.

INTEGER :: K, M
REAL   :: TEMP1(SIZE(F,1)), SDS(SIZE(F,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE        !
!        FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.                          !
!        --------------------------------------------------------------        !

IF (SHALLOW_RUN) THEN

   SDS = CONSS*FMEAN*EMEAN**2*AKMEAN**4
   FRES: DO M = 1,SIZE(F,3)
      TEMP1 = TFAK(INDEP,M)/AKMEAN
      TEMP1 = SDS * ((1.-DELTA)*TEMP1 +  DELTA*TEMP1**2)
      DIRS: DO K = 1,SIZE(F,2)
         SL(:,K,M) = SL(:,K,M) + TEMP1*F(:,K,M)
         FL(:,K,M) = FL(:,K,M) + TEMP1
      END DO DIRS
   END DO FRES

ELSE

   SDS = CONSD*EMEAN**2*FMEAN**9
   FRED: DO M = 1,SIZE(F,3)
      TEMP1 = (FR(M)/FMEAN)**2
      TEMP1 = SDS * ((1.-DELTA)*TEMP1 + DELTA*TEMP1**2)
      DIRD: DO K = 1,SIZE(F,2)
         SL(:,K,M) = SL(:,K,M) + TEMP1*F(:,K,M)
         FL(:,K,M) = FL(:,K,M) + TEMP1
      END DO DIRD
   END DO FRED

END IF

END SUBROUTINE SDISSIP

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SFBRK (F, SL, FL, EMEAN, FMEAN, DEPTH)

! ---------------------------------------------------------------------------- !
!
!     WEIMIN LUO, POL, MAY 1996, COMPUTATION OF WAVE BREAKING
!
!     PURPOSE
!     -------
!
!     COMPUTE DISSIPATION DUE TO DEPTH-INDUCED WAVE BREAKING
!
!     METHOD.
!     -------
!
!       SEE REFERENCES.
!
!     REFERENCES.
!     -----------
!
!     BATTJES & JANSSEN (1978)
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL, INTENT(IN)     :: F (:, :, :)  !! SPECTRUM.
REAL, INTENT(INOUT)  :: SL(:, :, :)  !! TOTAL SOURCE FUNCTION ARRAY
REAL, INTENT(INOUT)  :: FL(:, :, :)  !! DIAGONAL MATRIX OF FUNC. DERIVATIVE.
REAL, INTENT(IN)     :: EMEAN (:)    !! TOTAL ENERGY
REAL, INTENT(IN)     :: FMEAN(:)     !! MEAN FREQUENCY
REAL, INTENT(IN)     :: DEPTH(:)     !! WATER DEPTH

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: GAMD  = 0.8
REAL, PARAMETER :: ALPHA = 1.0

INTEGER :: M, K
REAL :: HMAX(SIZE(F,1)), HRMS(SIZE(F,1)), QB(SIZE(F,1)), BB(SIZE(F,1))
REAL :: SBR(SIZE(F,1)),DSBR(SIZE(F,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!   1. compute total dissipation rate according to Battjes-Janssen
!      -----------------------------------------------------------

HMAX = GAMD*DEPTH              !! compute Hmax
HRMS = SQRT(8.*EMEAN)          !! compute Hrms

CALL CMPQB(HRMS, HMAX, QB)

QB = MIN(1.,QB)
SBR = -ALPHA*2.*FMEAN

BB = (HRMS/HMAX)**2
WHERE (BB.LE.1.) SBR = SBR*QB/BB

WHERE (BB .LT. 1. .AND. ABS(BB - QB) .GT. 0.)
   DSBR = SBR * (1. - QB) / (BB - QB)
ELSEWHERE
   DSBR = 0.
ENDWHERE

DO M = 1,SIZE(F,3)
  DO K = 1,SIZE(F,2)
      SL(:,K,M) = SL(:,K,M)+SBR*F(:,K,M)
      FL(:,K,M) = FL(:,K,M)+DSBR
   END DO
END DO

END SUBROUTINE SFBRK

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CMPQB (HRMS, HMAX, QB)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   CMPQB - 
!
!     ADDED BY WEIMIN LUO, POL, MAY 1996
!     BASED ON THE CODE OF G. Ph. van Vledder, Delft Hydraulics
!
!  1. Purpose
!
!     Compute fraction of breaking waves for use in
!     SFBRK wave breaking dissipation function
!
!  2. Method
!
!     Newton-Raphson implementation of
!
!     1 - QB
!     ------ = - (HRMS/HMAX)^2
!     ln(QB)
!
!  3. Parameter list
!
!  4. Subroutines used
!
!  5. Error messages
!
!  6. Remarks
!
!     If HRMS > HMAX --> QB = 1
!
!  7. Structure
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

IMPLICIT NONE

REAL    , INTENT(IN)  :: HRMS(:)   !! Root mean square wave height
REAL    , INTENT(IN)  :: HMAX(:)   !! Maximum wave height
REAL    , INTENT(OUT) :: QB(:)     !! Fraction of breaking waves

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL    :: B(SIZE(HRMS)), Q0(SIZE(HRMS))
INTEGER :: I

! ---------------------------------------------------------------------------- !
!                                                                              !

QB = 1.
I = COUNT((HRMS.LT.HMAX))
IF (I.EQ.0) RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !

B = HRMS/HMAX
Q0 = 0.
WHERE (B.GE.0.5) Q0 = (2.*B-1.)**2

! ---------------------------------------------------------------------------- !
!                                                                              !

WHERE (HRMS.LT.HMAX)
   Q0 = Q0 - B**2*(Q0 - exp((Q0-1.)/B**2))/(B**2-exp((Q0-1.)/B**2))
END WHERE

WHERE (HRMS.LT.HMAX) QB = Q0

END SUBROUTINE CMPQB

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SINPUT (F, SL, FL, USTAR, UDIR, Z0, INDEP, LLWS)

! ---------------------------------------------------------------------------- !
!                                                                              !
! *** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.                         !
!                                                                              !
!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990                             !
!                                                                              !
!     OPTIMIZED BY : H. GUENTHER                                               !
!                                                                              !
!     MODIFIED BY : J-R BIDLOT NOVEMBER 1995                                   !
!                   J-R BIDLOT FEBRUARY 1996-97                                !
!                   J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT          !
!                                              BY ONLY USING NEW               !
!                                              STRESS AND ROUGHNESS.           !
!                                                                              !
!     PURPOSE.                                                                 !
!     ---------                                                                !
!                                                                              !
!       COMPUTE INPUT SOURCE FUNCTION AND THE FUNCTIONAL DERIVATIVE OF INPUT   !
!       SOURCE FUNCTION.                                                       !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       SEE REFERENCE.                                                         !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       P. JANSSEN, J.P.O., 1989.                                              !
!       P. JANSSEN, J.P.O., 1991.                                              !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: F (:, :, :)    !! SPECTRUM.
REAL,    INTENT(INOUT) :: SL(:, :, :)    !! TOTAL SOURCE FUNCTION ARRAY
REAL,    INTENT(INOUT) :: FL(:, :, :)    !! DIAGONAL MATRIX OF FUNCTIONAL
                                         !! DERIVATIVE
REAL,    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY.
REAL,    INTENT(IN)    :: UDIR (:)       !! WIND DIRECTION.
REAL,    INTENT(IN)    :: Z0   (:)       !! ROUGHNESS LENGTH.
INTEGER, INTENT(IN)    :: INDEP(:)       !! DEPTH TABLE INDEX.
LOGICAL, INTENT(OUT)   :: LLWS(:, :, :)  !! TRUE WHERE SINPUT IS POSITIVE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: K, M,IJ
REAL  :: X, ZLOG, ZLOG2X, CMT

REAL  :: UCN(SIZE(F,1)), ZCN(SIZE(F,1)), UFAC(SIZE(F,1)), CM(SIZE(F,1))
REAL  :: TEMP(SIZE(F,1),SIZE(F,2))
REAL  :: FAC(SIZE(F,3)), CONST(SIZE(F,3))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. PRECALCULATED ANGULAR DEPENDENCE.                                     !
!        ---------------------------------                                     !

DO K = 1,SIZE(F,2)
   TEMP(:,K) = COS(TH(K)-UDIR)
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. PRECALCULATED FREQUENCY DEPENDENCE.                                   !
!        -----------------------------------                                   !

FAC = ZPI*FR
CONST = FAC*XEPS*BETAMAX/XKAPPA**2
LLWS = .FALSE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. LOOP OVER FREQUENCIES.                                                !
!        ----------------------                                                !

FRE: DO M = 1,SIZE(F,3)

   IF (SHALLOW_RUN) THEN
      CM = TFAK(INDEP,M)/FAC(M)   !! INVERSE OF PHASE VELOCITY
      UCN = USTAR*CM + ZALP
      ZCN = LOG(G*Z0*CM**2)
   ELSE
      CMT = FAC(M)/G              !! INVERSE OF PHASE VELOCITY
      UCN = USTAR*CMT + ZALP
      ZCN = LOG(G*Z0*CMT**2)
   END IF

!     3.1 LOOP OVER DIRECTIONS.                                                !
!         ---------------------                                                !

   DIR: DO K = 1,SIZE(F,2)
      UFAC = 0.
      DO IJ = 1,SIZE(F,1)
         IF (TEMP(IJ,K).LE.0.01) CYCLE
         X    = TEMP(IJ,K)*UCN(IJ)
         ZLOG = ZCN(IJ) + XKAPPA/X
         IF (ZLOG.GE.0.) CYCLE
         ZLOG2X = ZLOG*ZLOG*X
         UFAC(IJ) = CONST(M)*EXP(ZLOG)*ZLOG2X*ZLOG2X
         LLWS(IJ,K,M) = .TRUE.
      END DO

!     3.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.                     !
!         ------------------------------------------------                     !

      SL(:,K,M) = UFAC*F(:,K,M)
      FL(:,K,M) = UFAC
   END DO DIR
END DO FRE

END SUBROUTINE SINPUT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SNONLIN (F, SL, FL, DEPTH, AKMEAN)

! ---------------------------------------------------------------------------- !
!                                                                              !
! *** *SNONLIN* - COMPUTATION OF NONLINEAR TRANSFER RATE AND ITS               !
! ***             FUNCTIONAL DERIVATIVE (DIAGONAL TERMS ONLY) AND              !
! ***             ADDITION TO CORRESPONDING NET EXPRESSIONS.                   !
!                                                                              !
!     S.D. HASSELMANN.  MPI                                                    !
!                                                                              !
!     G. KOMEN, P. JANSSEN   KNMI             MODIFIED TO SHALLOW WATER        !
!     H. GUENTHER, L. ZAMBRESKY               OPTIMIZED                        !
!     H. GUENTHER       GKSS/ECMWF  JUNE 1991 INTERACTIONS BETWEEN DIAG-       !
!                                             AND PROGNOSTIC PART.             !
!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90                            !
!     E. MYKLEBUST        FEBRUARY 2005       OPTIMIZATION                     !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       SEE ABOVE.                                                             !
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

REAL, INTENT(IN)    :: F (:,:,:)     !! SPECTRA.
REAL, INTENT(INOUT) :: SL(:,:,:)     !! TOTAL SOURCE FUNCTION ARRAY.
REAL, INTENT(INOUT) :: FL(:,:,:)     !! DIAGONAL MATRIX OF
REAL, INTENT(IN)    :: DEPTH (:)     !! WATER DEPTH.
REAL, INTENT(IN)    :: AKMEAN (:)    !! MEAN WAVE NUMBER.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: MC, MP, MP1, MM, MM1, IC, IP, IP1, IM, IM1, KH, K, K1, K2, K11, K21
INTEGER :: IJ, ML

REAL    :: FFACP, FFACP1, FFACM1, FTAIL, FKLAMP, FKLAMP1, GW1, GW2, GW3, GW4
REAL    :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAPA2, FKLAPB2, FKLAP12, FKLAP22
REAL    :: FKLAMM, FKLAMM1, GW5, GW6, GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2
REAL    :: FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
REAL    :: SAP, SAM, FIJ, FAD1, FAD2, FCEN

REAL    :: FTEMP(SIZE(F,1)), AD
REAL    :: DELAD, DELAP, DELAM, ENH(SIZE(F,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. SHALLOW WATER INITIALISATION.                                         !
!        -----------------------------                                         !

ML = SIZE(F,3)

IF (SHALLOW_RUN) THEN
   DO IJ = 1,SIZE(F,1)
      ENH(IJ) = MAX(0.75*DEPTH(IJ)*AKMEAN(IJ) , .5)
      ENH(IJ) = 1. + (5.5/ENH(IJ)) * (1.-.833*ENH(IJ)) * EXP(-1.25*ENH(IJ))
   END DO
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. FREQUENCY LOOP.                                                       !
!        ---------------                                                       !

FRE: DO MC = 1,ML+4
        MP  = IKP (MC)
        MP1 = IKP1(MC)
        MM  = IKM (MC)
        MM1 = IKM1(MC)
        FFACP  = 1.
        FFACP1 = 1.
        FFACM1 = 1.
        FTAIL  = 1.
        IC  = MC
        IP  = MP
        IP1 = MP1
        IM  = MM
        IM1 = MM1
        IF (IP1.GT.ML) THEN
          FFACP1 = FRH(IP1-ML+1)
          IP1 = ML
          IF (IP .GT.ML) THEN
            FFACP  = FRH(IP -ML+1)
            IP  = ML
            IF (IC .GT.ML) THEN
              FTAIL  = FRH(IC -ML+1)
              IC  = ML
              IF (IM1.GT.ML) THEN
                FFACM1 = FRH(IM1-ML+1)
                IM1 = ML
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        FKLAMP  = FKLAP(MC)
        FKLAMP1 = FKLAP1(MC)
        GW2 = FKLAMP1*FFACP*DAL1
        GW1 = GW2*CL11
        GW2 = GW2*ACL1
        GW4 = FKLAMP*FFACP1*DAL1
        GW3 = GW4*CL11
        GW4 = GW4*ACL1
        FKLAMPA = FKLAMP*CL11
        FKLAMPB = FKLAMP*ACL1
        FKLAMP2 = FKLAMP1*ACL1
        FKLAMP1 = FKLAMP1*CL11
        FKLAPA2 = FKLAMPA**2
        FKLAPB2 = FKLAMPB**2
        FKLAP12 = FKLAMP1**2
        FKLAP22 = FKLAMP2**2
        FKLAMM  = FKLAM(MC)
        FKLAMM1 = FKLAM1(MC)
        GW6 = FKLAMM1*DAL2
        GW5 = GW6*CL21
        GW6 = GW6*ACL2
        GW8 = FKLAMM*FFACM1*DAL2
        GW7 = GW8*CL21
        GW8 = GW8*ACL2
        FKLAMMA = FKLAMM*CL21
        FKLAMMB = FKLAMM*ACL2
        FKLAMM2 = FKLAMM1*ACL2
        FKLAMM1 = FKLAMM1*CL21
        FKLAMA2 = FKLAMMA**2
        FKLAMB2 = FKLAMMB**2
        FKLAM12 = FKLAMM1**2
        FKLAM22 = FKLAMM2**2

   IF (SHALLOW_RUN) THEN
      DO IJ = 1,SIZE(F,1)
         FTEMP(IJ) = AF11(MC)*ENH(IJ)
      END DO
   ELSE
      DO IJ = 1,SIZE(F,1)
         FTEMP(IJ) = AF11(MC) 
      END DO
   END IF
   
   IF (MC.GT.4) THEN! UNTIL 7.
      IF (MM1.LE.ML) THEN! UNTIL 6.
         IF (MC .LE.ML) THEN! UNTIL 5.
            IF (MP .LE.ML) THEN! UNTIL 4.
               IF (MP1.LE.ML) THEN! UNTIL 3.
               !     2.1 LOOP FOR ANLULAR SYMMETRY.                            !
                  MIR2: DO KH = 1,2
               !     2.1.1   ANGULAR LOOP.                                     !
                     DIR2: DO K = 1,SIZE(F,2)
                        K1  = K1W (K,KH)
                        K2  = K2W (K,KH)
                        K11 = K11W(K,KH)
                        K21 = K21W(K,KH)
               !     2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND     !
               !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.         !
               !             ----------------------------------------------    !

                        DO IJ = 1,SIZE(F,1)
                           SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )  &
&                              + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
                           SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )  &
&                              + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
                           FIJ = F(IJ,K  ,IC )*FTAIL
                           FAD1 = FIJ*(SAP+SAM)
                           FAD2 = FAD1-2.*SAP*SAM
                           FAD1 = FAD1+FAD2
                           FCEN = FTEMP(IJ)*FIJ
                           AD = FAD2*FCEN
                           DELAD = FAD1*FTEMP(IJ)
                           DELAP = (FIJ-2.*SAM)*DAL1*FCEN
                           DELAM = (FIJ-2.*SAP)*DAL2*FCEN
                           SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD*FKLAMM1
                           SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD*FKLAMM2
                           FL(IJ,K2 ,MM ) = FL(IJ,K2 ,MM ) + DELAM*FKLAM12
                           FL(IJ,K21,MM ) = FL(IJ,K21,MM ) + DELAM*FKLAM22

                           SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD*FKLAMMA
                           SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD*FKLAMMB
                           FL(IJ,K2 ,MM1) = FL(IJ,K2 ,MM1) + DELAM*FKLAMA2
                           FL(IJ,K21,MM1) = FL(IJ,K21,MM1) + DELAM*FKLAMB2

                           SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.*AD
                           FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD

                           SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD*FKLAMP1
                           SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD*FKLAMP2
                           FL(IJ,K1 ,MP ) = FL(IJ,K1 ,MP ) + DELAP*FKLAP12
                           FL(IJ,K11,MP ) = FL(IJ,K11,MP ) + DELAP*FKLAP22

                           SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1) + AD*FKLAMPA
                           SL(IJ,K11,MP1) = SL(IJ,K11,MP1) + AD*FKLAMPB
                           FL(IJ,K1 ,MP1) = FL(IJ,K1 ,MP1) + DELAP*FKLAPA2
                           FL(IJ,K11,MP1) = FL(IJ,K11,MP1) + DELAP*FKLAPB2
                        END DO
                     END DO DIR2!!  BRANCH BACK TO 2.1.1 FOR NEXT DIRECTION.
                  END DO MIR2   !!  BRANCH BACK TO 2.1 FOR MIRROR INTERACTIONS.
               ELSE!IF (MP1.LE.ML) THEN
                 !     3.1 LOOP FOR ANLULAR SYMMETRY.                          !
                  MIR3: DO KH = 1,2
                 !     3.1.1   ANGULAR LOOP.                                   !
                     DIR3: DO K = 1,SIZE(F,2)
                        K1  = K1W (K,KH)
                        K2  = K2W (K,KH)
                        K11 = K11W(K,KH)
                        K21 = K21W(K,KH)
                 !     3.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND   !
                 !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.       !
                 !             ----------------------------------------------  !

                        DO IJ = 1,SIZE(F,1)
                           SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )  &
&                              + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
                           SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )  &
&                              + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
                           FIJ = F(IJ,K  ,IC )*FTAIL
                           FAD1 = FIJ*(SAP+SAM)
                           FAD2 = FAD1-2.*SAP*SAM
                           FAD1 = FAD1+FAD2
                           FCEN = FTEMP(IJ)*FIJ
                           AD = FAD2*FCEN
                           DELAD = FAD1*FTEMP(IJ)
                           DELAP = (FIJ-2.*SAM)*DAL1*FCEN
                           DELAM = (FIJ-2.*SAP)*DAL2*FCEN
                           SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD*FKLAMM1
                           SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD*FKLAMM2
                           FL(IJ,K2 ,MM ) = FL(IJ,K2 ,MM ) + DELAM*FKLAM12
                           FL(IJ,K21,MM ) = FL(IJ,K21,MM ) + DELAM*FKLAM22

                           SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD*FKLAMMA
                           SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD*FKLAMMB
                           FL(IJ,K2 ,MM1) = FL(IJ,K2 ,MM1) + DELAM*FKLAMA2
                           FL(IJ,K21,MM1) = FL(IJ,K21,MM1) + DELAM*FKLAMB2

                           SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.*AD
                           FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD

                           SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD*FKLAMP1
                           SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD*FKLAMP2
                           FL(IJ,K1 ,MP ) = FL(IJ,K1 ,MP ) + DELAP*FKLAP12
                           FL(IJ,K11,MP ) = FL(IJ,K11,MP ) + DELAP*FKLAP22
                        END DO
                     END DO DIR3!  BRANCH BACK TO 3.1.1 FOR NEXT DIRECTION.
                  END DO MIR3   !  BRANCH BACK TO 3.1 FOR MIRROR INTERACTIONS.
               END IF
            ELSE!IF (MP .LE.ML) THEN
          !     4.1 LOOP FOR ANLULAR SYMMETRY.                                 !
               MIR4: DO KH = 1,2
          !     4.1.1   ANGULAR LOOP.                                          !
                  DIR4: DO K = 1,SIZE(F,2)
                     K1  = K1W (K,KH)
                     K2  = K2W (K,KH)
                     K11 = K11W(K,KH)
                     K21 = K21W(K,KH)
          !     4.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND          !
          !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.              !
          !             ----------------------------------------------         !

                     DO IJ = 1,SIZE(F,1)
                        SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )     &
&                           + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
                        SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )     &
&                           + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
                        FIJ = F(IJ,K  ,IC )*FTAIL
                        FAD1 = FIJ*(SAP+SAM)
                        FAD2 = FAD1-2.*SAP*SAM
                        FAD1 = FAD1+FAD2
                        FCEN = FTEMP(IJ)*FIJ
                        AD = FAD2*FCEN
                        DELAD = FAD1*FTEMP(IJ)
                        DELAP = (FIJ-2.*SAM)*DAL1*FCEN
                        DELAM = (FIJ-2.*SAP)*DAL2*FCEN
                        SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD*FKLAMM1
                        SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD*FKLAMM2
                        FL(IJ,K2 ,MM ) = FL(IJ,K2 ,MM ) + DELAM*FKLAM12
                        FL(IJ,K21,MM ) = FL(IJ,K21,MM ) + DELAM*FKLAM22

                        SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD*FKLAMMA
                        SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD*FKLAMMB
                        FL(IJ,K2 ,MM1) = FL(IJ,K2 ,MM1) + DELAM*FKLAMA2
                        FL(IJ,K21,MM1) = FL(IJ,K21,MM1) + DELAM*FKLAMB2

                        SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.*AD
                        FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD
                     END DO
                  END DO DIR4  !! BRANCH BACK TO 4.1.1 FOR NEXT DIRECTION.
               END DO MIR4     !! BRANCH BACK TO 4.1 FOR MIRROR INTERACTIONS.
            END IF
         ELSE!IF (MC .LE.ML) THEN
           !     5.1 LOOP FOR ANLULAR SYMMETRY.                                !
            MIR5: DO KH = 1,2
           !     5.1.1   ANGULAR LOOP.                                         !
               DIR5: DO K = 1,SIZE(F,2)
                  K1  = K1W (K,KH)
                  K2  = K2W (K,KH)
                  K11 = K11W(K,KH)
                  K21 = K21W(K,KH)
           !     5.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND         !
           !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.             !
           !             ----------------------------------------------        !

                  DO IJ = 1,SIZE(F,1)
                     SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )        &
&                        + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
                     SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )        &
&                        + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
                     FIJ = F(IJ,K  ,IC )*FTAIL
                     FAD1 = FIJ*(SAP+SAM)
                     FAD2 = FAD1-2.*SAP*SAM
                     FAD1 = FAD1+FAD2
                     FCEN = FTEMP(IJ)*FIJ
                     AD = FAD2*FCEN
                     DELAD = FAD1*FTEMP(IJ)
                     DELAP = (FIJ-2.*SAM)*DAL1*FCEN
                     DELAM = (FIJ-2.*SAP)*DAL2*FCEN
                     SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD*FKLAMM1
                     SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD*FKLAMM2
                     FL(IJ,K2 ,MM ) = FL(IJ,K2 ,MM ) + DELAM*FKLAM12
                     FL(IJ,K21,MM ) = FL(IJ,K21,MM ) + DELAM*FKLAM22

                     SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD*FKLAMMA
                     SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD*FKLAMMB
                     FL(IJ,K2 ,MM1) = FL(IJ,K2 ,MM1) + DELAM*FKLAMA2
                     FL(IJ,K21,MM1) = FL(IJ,K21,MM1) + DELAM*FKLAMB2
                  END DO
               END DO DIR5    !! BRANCH BACK TO 5.1.1 FOR NEXT DIRECTION.
            END DO MIR5       !! BRANCH BACK TO 5.1 FOR MIRROR INTERACTIONS.
         END IF
      ELSE!IF (MM1.LE.ML) THEN
         !     6.1 LOOP FOR ANLULAR SYMMETRY.                                  !
         MIR6: DO KH = 1,2
            !     6.1.1   ANGULAR LOOP.                                        !
            DIR6: DO K = 1,SIZE(F,2)
               K1  = K1W (K,KH)
               K2  = K2W (K,KH)
               K11 = K11W(K,KH)
               K21 = K21W(K,KH)
               !     6.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND     !
               !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.         !
               !             ----------------------------------------------    !

              DO IJ = 1,SIZE(F,1)
                  SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )           &
&                     + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
                  SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )           &
&                     + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
                  FIJ = F(IJ,K  ,IC )*FTAIL
                  FAD1 = FIJ*(SAP+SAM)
                  FAD2 = FAD1-2.*SAP*SAM
                  FAD1 = FAD1+FAD2
                  FCEN = FTEMP(IJ)*FIJ
                  AD = FAD2*FCEN
                  DELAD = FAD1*FTEMP(IJ)
                  DELAP = (FIJ-2.*SAM)*DAL1*FCEN
                  DELAM = (FIJ-2.*SAP)*DAL2*FCEN
                  SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD*FKLAMM1
                  SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD*FKLAMM2
                  FL(IJ,K2 ,MM ) = FL(IJ,K2 ,MM ) + DELAM*FKLAM12
                  FL(IJ,K21,MM ) = FL(IJ,K21,MM ) + DELAM*FKLAM22
               END DO
            END DO DIR6       !! BRANCH BACK TO 6.1.1 FOR NEXT DIRECTION.
         END DO MIR6          !! BRANCH BACK TO 6.1 FOR MIRROR INTERACTIONS.
      END IF
   ELSE!IF (MC.GT.4) THEN
      !     7.1 LOOP FOR ANLULAR SYMMETRY.                                     !
      MIR7: DO KH = 1,2
         !     7.1.1   ANGULAR LOOP.                                           !
         DIR7: DO K = 1,SIZE(F,2)
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)
            !     6.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND     !
            !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.         !
            !             ----------------------------------------------    !

            DO IJ = 1,SIZE(F,1)
               SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )              &
&                  + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
               FIJ = F(IJ,K,IC)
               FAD2 = FIJ*SAP
               FAD1 = 2.*FAD2
               FCEN = FTEMP(IJ)*FIJ
               AD = FAD2*FCEN
               DELAD = FAD1*FTEMP(IJ)
               DELAP = FIJ*DAL1*FCEN
               SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.*AD
               SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD*FKLAMP1
               SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD*FKLAMP2
               SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1) + AD*FKLAMPA
               SL(IJ,K11,MP1) = SL(IJ,K11,MP1) + AD*FKLAMPB
               FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD
               FL(IJ,K1 ,MP ) = FL(IJ,K1 ,MP ) + DELAP*FKLAP12
               FL(IJ,K11,MP ) = FL(IJ,K11,MP ) + DELAP*FKLAP22
               FL(IJ,K1 ,MP1) = FL(IJ,K1 ,MP1) + DELAP*FKLAPA2
             FL(IJ,K11,MP1) = FL(IJ,K11,MP1) + DELAP*FKLAPB2
            END DO
         END DO DIR7         !! BRANCH BACK TO 7.1.1 FOR NEXT DIRECTION.
      END DO MIR7            !! BRANCH BACK TO 7.1 FOR MIRROR INTERACTIONS.
   END IF

END DO FRE                  !! BRANCH BACK TO 2. FOR NEXT FREQUENCY.
END SUBROUTINE SNONLIN
 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SOURCE_PHILLIPS (SL, USTAR, UDIR, DEPTH, INDEP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SOURCE_PHILLIPS - COMPUTATION OF PHILLIPS INPUT.                           !
!                                                                              !
!     ROOP LALBEHARRY          ARMN/MSC        NOVEMBER 2003
!                                                                              !
!*    PURPOSE.                                                                 !
!     ---------                                                                !
!                                                                              !
!       COMPUTE PHILLIPS INPUT SOURCE FUNCTION FACTOR.                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       SEE REFERENCE.                                                         !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       O. PHILLIPS 1966, CAVALERI AND RIZZOLI 1981, H. TOLMAN 1992            !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!*    INTERFACE VARIABLE

REAL,   INTENT(INOUT)        :: SL(:,:,:)  !! TOTAL SOURCE FUNCTION.
REAL,   INTENT(IN)           :: USTAR(:)   !! FRICTION VELOCITY
REAL,   INTENT(IN)           :: UDIR(:)    !! WIND DIRECTION
REAL,   INTENT(IN)           :: DEPTH(:)   !! DEPTH
INTEGER, INTENT(IN)          :: INDEP(:)   !! DEPTH TABLE INDEX.

! ---------------------------------------------------------------------------- !
!
!*    LOCAL VARIABLE

INTEGER :: K, M, IJ
 
REAL    :: CONST1
REAL    :: TEMP(1:SIZE(SL,1)), FPM(1:SIZE(SL,1)), KD
REAL    :: TPHOLD(1:SIZE(SL,1),1:SIZE(SL,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    0. CONSTANTS.                                                            !
!        ----------                                                            !

CONST1 = 28/(0.13*G)
DO IJ = 1,SIZE(SL,1)
   FPM(IJ) = MAX(0.0000001, CONST1*USTAR(IJ))
END DO

CONST1 = (80.*16.0*XEPS**2)/(0.5*3.0*G**2)

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    1. PRECALCULATED ANGULAR DEPENDENCE.                                     !
!        ---------------------------------                                     !

DO K = 1,SIZE(SL,2)
   DO IJ = 1,SIZE(SL,1)
      TPHOLD(IJ,K) = MAX(0.,COS(TH(K)-UDIR(IJ)))
      TPHOLD(IJ,K) = (USTAR(IJ)*TPHOLD(IJ,K))**4
   END DO
END DO

IF (SHALLOW_RUN) THEN

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    2. SHALLOW WATER                                                         !
!        -------------                                                         !

   DO M = 1,SIZE(SL,3)
      DO IJ = 1,SIZE(SL,1)
        KD = MIN(TFAK(INDEP(IJ),M)*DEPTH(IJ),40.)
        TEMP(IJ) = EXP(-(FR(M)*FPM(IJ))**(-4))
        TEMP(IJ)=CONST1*TEMP(IJ)/(1.0+2.0*KD/SINH(2.*KD))
      END DO

      DO K = 1,SIZE(SL,2)
         DO IJ = 1,SIZE(SL,1)
            SL(IJ,K,M) = SL(IJ,K,M) + TEMP(IJ) * TPHOLD(IJ,K)
         END DO
      END DO
   END DO

ELSE

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    3. DEEP WATER (INDEP NOT DEFINED FOR DEEP WATER)                         !
!        ----------                                                            !

   DO M = 1,SIZE(SL,3)
      DO IJ = 1,SIZE(SL,1)
         TEMP(IJ) = CONST1*EXP(-(FR(M)*FPM(IJ))**(-4))
      END DO

      DO K=1,SIZE(SL,2)
         DO IJ = 1,SIZE(SL,1)
            SL(IJ,K,M) = SL(IJ,K,M) + TEMP(IJ) * TPHOLD(IJ,K)
         END DO
      END DO
   END DO

END IF

END SUBROUTINE SOURCE_PHILLIPS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE STRESS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   STRESS - COMPUTATION OF TOTAL STRESS.                                      !
!                                                                              !
!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990                             !
!     J. BIDLOT           ECMWF     SEPTEMBER 1996 : REMOVE Z0 DUE TO          !
!                                   VISCOSITY AND ADD ZERO STRESS FOR          !
!                                   ZERO WIND.                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     ---------                                                                !
!                                                                              !
!       TO GENERATE STRESS TABLE TAU(TAUW,U10).                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       A STEADY STATE WIND PROFILE IS ASSUMED.                                !
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH                  !
!                                                                              !
!                  Z1=Z0/SQRT(1-TAUW/TAU)                                      !
!                                                                              !
!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-                  !
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.                            !
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.           !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL,    PARAMETER :: XM    = 0.50    !! POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
INTEGER, PARAMETER :: NITER = 10      !! NUMBER OF ITERATIONS TO OBTAIN
                                      !! TOTAL STRESS
REAL,    PARAMETER :: EPS1  = 0.00001 !! SMALL NUMBER TO MAKE SURE THAT A
                                      !! SOLUTION IS OBTAINED IN ITERATION
                                      !! WITH TAU>TAUW.

INTEGER :: I, J, ITER
REAL    :: TAUWMAX, ZTAUW, UTOP, CDRAG, WCD, USTOLD, TAUOLD, X, UST, Z0
REAL    :: F, DELF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALISE CONSTANTS.                                                 !
!        ---------------------                                                 !

TAUWMAX = SQRT(USTARM)
DELU    = UMAX/REAL(JUMAX)
DELTAUW = TAUWMAX/REAL(ITAUMAX)
CDRAG   = 0.0012875
WCD     = SQRT(CDRAG)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2 DETERMINE STRESS.                                                      !
!       -----------------                                                      !

DO I = 0,ITAUMAX
   ZTAUW   = (REAL(I)*DELTAUW)**2
   DO J = 0,JUMAX
      UTOP    = REAL(J)*DELU
      USTOLD  = UTOP*WCD
      TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)

      DO ITER = 1,NITER
         X    = ZTAUW/TAUOLD
         UST  = SQRT(TAUOLD)
         Z0   = ALPHA*TAUOLD/(G)/(1.-X)**XM
         F    = UST-XKAPPA*UTOP/(LOG(XNLEV/Z0))
         DELF = 1.-XKAPPA*UTOP/(LOG(XNLEV/Z0))**2*2./UST*(1.-(XM+1)*X)/(1.-X)
         UST  = UST-F/DELF
         TAUOLD =  MAX(UST**2., ZTAUW+EPS1)
      END DO
      TAUT(I,J)  = SQRT(TAUOLD)
   END DO
END DO  !! END LOOP OVER INDICES OF TAU-TABLE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. FORCE ZERO WIND TO HAVE ZERO STRESS                                   !

TAUT(0:ITAUMAX,0)=0.0

END SUBROUTINE STRESS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE STRESSO (F, SL, USTAR, UDIR, Z0, TAUW)

! ---------------------------------------------------------------------------- !
!                                                                              !
! *** *STRESSO* - COMPUTATION OF WAVE STRESS.                                  !
!                                                                              !
!     H. GUNTHER      GKSS/ECMWF  NOVEMBER  1989 CODE MOVED FROM SINPUT.       !
!     P.A.E.M. JANSSEN      KNMI  AUGUST    1990                               !
!     J. BIDLOT             ECMWF FEBRUARY  1996-97                            !
!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90                            !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY AND DIRECTIONS. !
!       BECAUSE ARRAY *SL* IS USED, ONLY THE INPUT SOURCE HAS TO BE STORED IN  !
!       *SL* (CALL FIRST SINPUT, THEN STRESSO, AND THEN THE REST OF THE SOURCE !
!       FUNCTIONS)                                                             !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       R SNYDER ET AL,1981.                                                   !
!       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.                  !
!       P. JANSSEN, JPO, 1985                                                  !
!                                                                              !

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: F(:,:,:)       !! SPECTRUM.
REAL,    INTENT(IN)    :: SL(:,:,:)      !! INPUT SOURCE FUNCTION.
REAL,    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY.
REAL,    INTENT(IN)    :: UDIR(:)        !! WIND DIRECTION.
REAL,    INTENT(IN)    :: Z0(:)          !! ROUGHNESS LENGTH.
REAL,    INTENT(OUT)   :: TAUW(:)        !! WAVE STRESS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: I, J, K, M, KL, ML, IJ

REAL    :: CONST0, UST
REAL    :: XI, DELI1, DELI2, XJ, DELJ1, DELJ2

REAL    :: CONSTF(SIZE(F,3))
REAL    :: TEMP(SIZE(F,1)), XSTRESS(SIZE(F,1)), YSTRESS(SIZE(F,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. PRECOMPUTE FREQUENCY SCALING.                                         !
!        -----------------------------                                         !

KL = SIZE(F,2)
ML = SIZE(F,3)

CONST0  = DELTH*(ZPI)**4*FR(ML)**5/G**2
CONSTF =ZPI*XINVEPS*FR*DFIM

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE WAVE STRESS OF BLOCK.                                         !
!        -----------------------------                                         !
!                                                                              !
!     2.1 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.       !
!         --------------------------------------------------------------       !

XSTRESS = 0.
YSTRESS = 0.
DO M = 1,ML
   DO K = 1,KL
   DO IJ = 1,SIZE(F,1)
      XSTRESS(IJ) = XSTRESS(IJ) + SL(IJ,K,M)*CONSTF(M)*SINTH(K)
      YSTRESS(IJ) = YSTRESS(IJ) + SL(IJ,K,M)*CONSTF(M)*COSTH(K)
   END DO
   END DO
END DO

!     2.2 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.                     !
!     ----------------------------------------------------                     !

TEMP = 0.
DO K = 1,KL
   DO IJ = 1,SIZE(F,1)
      TEMP(IJ)  = TEMP(IJ) + F(IJ,K,ML)*MAX(COS(TH(K)-UDIR(IJ)),0.)**3
   END DO
END DO

DO IJ = 1,SIZE(F,1)
   UST   = MAX(USTAR(IJ),0.000001)
   XI    = MIN(REAL(IUSTAR), UST/DELUST)
   I     = MAX (0, MIN (IUSTAR-1, INT(XI)))
   DELI1 = MIN (1. ,XI-REAL(I))
   DELI2   = 1. - DELI1

   XJ    = MIN(REAL(IALPHA), (G*Z0(IJ)/UST**2-ALPHA)/DELALP)
   J     = MAX (0, MIN (IALPHA-1, INT(XJ)))
   DELJ1 = MIN (1. , XJ-REAL(J))
   DELJ2   = 1. - DELJ1

   TAUW(IJ) = (TAUHFT(I,J  )*DELI2 + TAUHFT(I+1,J  )*DELI1)*DELJ2              &
&           + (TAUHFT(I,J+1)*DELI2 + TAUHFT(I+1,J+1)*DELI1)*DELJ1
   TAUW(IJ) = CONST0*TEMP(IJ)*USTAR(IJ)**2*TAUW(IJ)
   XSTRESS(IJ) = XSTRESS(IJ) + TAUW(IJ)*SIN(UDIR(IJ))
   YSTRESS(IJ) = YSTRESS(IJ) + TAUW(IJ)*COS(UDIR(IJ))
   TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)
END DO


END SUBROUTINE STRESSO

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TAUHF (FRMAX)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TAUHF - COMPUTATION OF HIGH-FREQUENCY STRESS.                              !
!                                                                              !
!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     ---------                                                                !
!                                                                              !
!       COMPUTE HIGH-FREQUENCY WAVE STRESS                                     !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.                             !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL, INTENT(IN)   :: FRMAX               !! LAST MODEL FREQUENCY FR(ML).

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER, PARAMETER :: JTOT = 250

INTEGER :: L, K, J
REAL    :: ALPHAM, CONST1, OMEGAC, X0, UST, Z0,                        &
&          OMEGACC, YC, DELY, Y, OMEGA, CM, ZX, ZARG, ZMU, ZLOG, ZBETA

REAL :: W(1:JTOT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. PRELIMINARY CALCULATIONS.                                             !
!        -------------------------                                             !

ALPHAM = 20.*ALPHA
DELUST = USTARM/REAL(IUSTAR)
DELALP = ALPHAM/REAL(IALPHA)

CONST1 = BETAMAX/XKAPPA**2
OMEGAC = ZPI*FRMAX

TAUHFT(0:IUSTAR,0:IALPHA) = 0.

W = 1.
W(1) = 0.5
W(JTOT) = 0.5

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.                      !
!        ------------------------------------------------                      !

X0 = 0.05
DO L = 0,IALPHA
   DO K = 0,IUSTAR
      UST      = MAX(REAL(K)*DELUST,0.000001)
      Z0       = UST**2*(ALPHA+REAL(L)*DELALP)/G
      OMEGACC  = MAX(OMEGAC,X0*G/UST)
      YC       = OMEGACC*SQRT(Z0/G)
      DELY     = MAX((1.-YC)/REAL(JTOT),0.)
      DO J = 1,JTOT
         Y        = YC+REAL(J-1)*DELY
         OMEGA    = Y*SQRT(G/Z0)
         CM       = G/OMEGA
         ZX       = UST/CM +ZALP
         ZARG     = MIN(XKAPPA/ZX,20.)
         ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)
         ZLOG         = MIN(LOG(ZMU),0.)
         ZBETA        = CONST1*ZMU*ZLOG**4
         TAUHFT(K,L)  = TAUHFT(K,L)+W(J)*ZBETA/Y*DELY
      END DO
   END DO
END DO

END SUBROUTINE TAUHF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_SOURCE_MODULE

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
&       WM1_WM2_WAVENUMBER,      &  !! COMPUTATION OF MEAN WAVENUMBER.
&       TRANSF                      !! NARROW BAND LIMIT BENJAMIN-FEIR INDEX
                                    !! FOR THE FINITE DEPTH.
USE WAM_GENERAL_MODULE, ONLY:    &
&       AKI,                     &  !! WAVE NUMBER FROM FREQUENY AND DEPTH
&       ABORT1                      !! TERMINATE PROCESSING.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE, ONLY: G, PI, ZPI, DEG, ROAIR, RNUAIR, RNUAIRM, ROWATER,&
&                             RAD, XEPS, XNLEV, XINVEPS, XKAPPA,               &
&                             ALPHA,     BETAMAX,    ZALP,    TAUWSHELTER,     &
&                             SDSBR, ISDSDTH , ISB, IPSAT, SSDSC2, SSDSC4,     &
&                             SSDSC6,  MICHE, SSDSC3, SSDSBRF1, BRKPBCOEF,     &
&                             SSDSC5

USE WAM_FRE_DIR_MODULE, ONLY: KL, ML, FR, CO, TH, DELTH, COSTH, SINTH, DFIM,   &
&                             C, FMIN, FR5, FRM5, RHOWG_DFIM, INV_LOG_CO
USE WAM_TIMOPT_MODULE,  ONLY: IDELT, SHALLOW_RUN, IPHYS, WAVE_BREAKING_RUN,    &
&                             PHILLIPS_RUN, ISNONLIN, LCFLX
USE WAM_FILE_MODULE,    ONLY: IU06, ITEST
USE WAM_TABLES_MODULE,  ONLY: FLMINFR, NDEPTH, TFAK, TCGOND, T_TAIL, EPS1,     &
&                             JUMAX, DELU  

USE WAM_FLUX_MODULE,    ONLY: PHIOC, PHIAW, TAUOC_X, TAUOC_Y,                  &
&                             PHIBOT, TAUBOT_X, TAUBOT_Y
use wam_mpi_module,     only: NIJS, NIJL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

REAL, PARAMETER :: GAMD  = 0.8  !! Parameter of depth limited wave height

! ---------------------------------------------------------------------------- !
!
!    1. INDICES AND WEIGHTS USED TO COMPUTE THE NONLINEAR TRANSFER RATE.
!       ----------------------------------------------------------------

INTEGER, PARAMETER :: NINL = 5  !! SIZE OF INLCOEF
INTEGER, PARAMETER :: NRNL = 25 !! SIZE OF RNLCOEF

INTEGER :: KFRH     !! SIZE OF FRH
INTEGER :: MFRSTLW  !! INDEX OF FIRST EXTRA LOW FREQUENCY FOR SNL
INTEGER :: MLSTHG   !! INDEX OF LAST EXTRA HIGH FREQUENCY FOR SNL

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
INTEGER,ALLOCATABLE, DIMENSION(:,:) :: INLCOEF  !! STORES ALL FREQUENCY DEPENDENT
                                                !! INDICES FOUND IN SNL

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

REAL, ALLOCATABLE, DIMENSION(:)     :: FRH    !! TAIL FREQUENCY RATION **5
REAL, ALLOCATABLE, DIMENSION(:)     :: FTRF   !! FRONT TAIL REDUCTION FACTOR 
                                              !! USED TO A SPECTRAL TAIL IN FRONT
                                              !! OF THE FIRST DISCRETISED FREQUENCY
REAL, ALLOCATABLE, DIMENSION(:,:)   :: RNLCOEF !! STORES ALL FREQUENCY DEPENDENT
                                               !! COEFFICIENT FOUND IN SNL

REAL    :: ACL1        !! WEIGHT IN ANGULAR GRID FOR INTERPOLATION,
                          !! WAVE NO. 3 ("1+LAMBDA" TERM).
REAL    :: ACL2        !! WEIGHT IN ANGULAR GRID FOR INTERPOLATION,
                          !! WAVE NO. 4 ("1-LAMBDA" TERM).
REAL    :: CL11        !! 1.-ACL1.
REAL    :: CL21        !! 1.-ACL2.
REAL    :: DAL1        !! 1./ACL1.
REAL    :: DAL2        !! 1./ACL2.

! ---------------------------------------------------------------------------- !
!                                                                              !
!    2. NONLINEAR TRANSFER FUNCTION COEFFICIENTS FOR SHALLOW WATER.            !
!       -----------------------------------------------------------            !

REAL,   ALLOCATABLE, DIMENSION(:,:) :: ENH

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. INTEGRATION WEIGHT FOR HIGH FREQ STRESS AND ENERGY (SUB.TAU_PHI_HF).   !
!       --------------------------------------------------------------         !

INTEGER, PARAMETER :: JTOT_TAUHF=19  !! DIMENSION OF WTAUHF. IT MUST BE ODD !!!
REAL               :: WTAUHF(JTOT_TAUHF) !! INTEGRATION WEIGHT FOR TAU_PHI_HF
REAL               :: X0TAUHF        !! LOWEST LIMIT FOR INTEGRATION IN TAU_PHI_HF: X0 *(G/USTAR)

! ---------------------------------------------------------------------------- !
!
!    4. TABLE FOR SINPUT_ARD AND SDISSIP_ARD.
!       -------------------------------------

! FOR SINPUT_ARD
INTEGER, PARAMETER :: IAB=200       !! DIMENSION OF SWELLFT.
REAL               :: SWELLFT(IAB)  !! FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS.

! FOR SDISSIP_ARD
INTEGER :: NSDSNTH  !! NUMBER OF DIRECTIONS TO COMPUTE THE SPECTRAL SATURATION
INTEGER :: NDIKCUMUL !! INTEGER DIFFERENCE IN FREQUENCY BANDS.
INTEGER, ALLOCATABLE :: INDICESSAT(:,:)
REAL, ALLOCATABLE :: SATWEIGHTS(:,:)
REAL, ALLOCATABLE :: CUMULW(:,:,:,:)

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

INTERFACE MAKE_SHALLOW_SNL    !! COMPUTES THE NONLINEAR TRANSFER FUNCTION ENH.
   MODULE PROCEDURE MAKE_SHALLOW_SNL
END INTERFACE
PUBLIC MAKE_SHALLOW_SNL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE NLWEIGT             !! COMPUTATION OF INDEX ARRAYS AND WEIGHTS
   MODULE PROCEDURE NLWEIGT   !! FOR THE NONLINEAR TRANSFER RATE.
END INTERFACE
PRIVATE NLWEIGT

INTERFACE INIT_SNONLIN        !! COMPUTATION OF INDEX ARRAYS AND WEIGHTS
   MODULE PROCEDURE INIT_SNONLIN   !! FOR THE NONLINEAR TRANSFER RATE.
END INTERFACE
PRIVATE INIT_SNONLIN

INTERFACE SBOTTOM             !! COMPUTATION OF BOTTOM FRICTION.
   MODULE PROCEDURE SBOTTOM
END INTERFACE
PRIVATE SBOTTOM

INTERFACE SDISSIP             !! COMPUTATION OF DISSIPATION SOURCE FUNCTION.
   MODULE PROCEDURE SDISSIP
END INTERFACE
PRIVATE SDISSIP

INTERFACE SDISSIP_ARD             !! COMPUTATION OF DISSIPATION SOURCE FUNCTION.
   MODULE PROCEDURE SDISSIP_ARD
END INTERFACE
PRIVATE SDISSIP_ARD

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

INTERFACE SINPUT_ARD              !! COMPUTATION OF INPUT SOURCE FUNCTION (ECMWF (CY46R1).
   MODULE PROCEDURE SINPUT_ARD
END INTERFACE
PRIVATE SINPUT_ARD

INTERFACE WSIGSTAR                 !! COMPUTES STANDARD DEVIATION OF USTAR.
   MODULE PROCEDURE  WSIGSTAR
END INTERFACE
PRIVATE WSIGSTAR 

INTERFACE SNONLIN             !! COMPUTATION OF NONLINEAR TRANSFER RATE.
   MODULE PROCEDURE SNONLIN
END INTERFACE
PRIVATE SNONLIN
 
INTERFACE SOURCE_PHILLIPS     !! COMPUTATION OF PHILLIP'S INPUT SOURCE FUNCTION.
   MODULE PROCEDURE SOURCE_PHILLIPS
END INTERFACE
PRIVATE SOURCE_PHILLIPS

INTERFACE STRESSO             !! COMPUTATION OF WAVE STRESS.
   MODULE PROCEDURE STRESSO
END INTERFACE
PRIVATE STRESSO

INTERFACE FRCUTINDEX         !!  Frequency cut index
   MODULE PROCEDURE FRCUTINDEX
END INTERFACE
PRIVATE FRCUTINDEX

INTERFACE IMPHFTAIL         !!  Diagnostic high frequency tail
   MODULE PROCEDURE IMPHFTAIL
END INTERFACE
PRIVATE IMPHFTAIL

INTERFACE SDEPTHLIM         !!  LIMITS THE SPECTRAL VARIANCE
   MODULE PROCEDURE SDEPTHLIM
END INTERFACE
PRIVATE SDEPTHLIM

INTERFACE TABU_SWELLFT        !! FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS
   MODULE PROCEDURE TABU_SWELLFT
END INTERFACE
PRIVATE TABU_SWELLFT

INTERFACE INIT_SDISSP_ARD     !! INITIALISATION FOR SDISS_ARD
   MODULE PROCEDURE INIT_SDISSP_ARD 
END INTERFACE
PRIVATE INIT_SDISSP_ARD

INTERFACE TAU_PHI_HF !! COMPUTATION OF HIGH-FREQUENCY STRESS AND
                     !! HIGH-FREQUENCY ENERGY FLUX.
   MODULE PROCEDURE TAU_PHI_HF
END INTERFACE
PRIVATE TAU_PHI_HF

INTERFACE INIT_X0TAUHF        !! INITIALISATION FOR TAU_PHI_HF
   MODULE PROCEDURE INIT_X0TAUHF
END INTERFACE
PRIVATE INIT_X0TAUHF

INTERFACE WNFLUXES            !! DETERMINE SURFACE WAVE FLUXES.
   MODULE PROCEDURE WNFLUXES
END INTERFACE
PRIVATE WNFLUXES

INTERFACE KERKEI              !! Computes zeroth order Kelvin function Ker and Kei
   MODULE PROCEDURE KERKEI
END INTERFACE
PRIVATE KERKEI

INTERFACE KZEONE              !! MODIFIED BESSEL FUNCTIONS OF THE SECOND KIND,K0 AND K1
   MODULE PROCEDURE KZEONE
END INTERFACE
PRIVATE KZEONE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE IMPLSCH (FL3, U10, UDIR, TAUW, USTAR, Z0, ROAIRN, WSTAR, DEPTH, INDEP)

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
!     ROOP LALBEHARRY  MSC/ARMN     APRIL    2003 PHILLIPS SOURCE              !
!     ERIK MYKLEBUST                FEBRUARY 2005 OPTIMIZATION                 !
!     H. GUENTHER      GKSS         JANUARY  2014 Wave breaking                !
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
REAL,    INTENT(IN)    :: ROAIRN(:)      !! SURFACE AIR DENSITY [kg/m**3].
REAL,    INTENT(IN)    :: WSTAR (:)      !! CONVECTIVE VELOCITY SCALE [m/s]
REAL,    INTENT(INOUT) :: TAUW  (:)      !! WAVE STRESS IN (M/S)**2 
REAL,    INTENT(OUT)   :: USTAR (:)      !! FRICTION VELOCITY [M/S].
REAL,    INTENT(INOUT) :: Z0    (:)      !! ROUGHNESS LENGTH [M].
REAL,    INTENT(IN)    :: DEPTH (:)      !! WATER DEPTH [M].
INTEGER, INTENT(IN)    :: INDEP (:)      !! DEPTH TABLE INDEX.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !

INTEGER :: K, M, KL, ML
REAL    :: DELT

REAL    :: FL(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3))  !! DIAGONAL MATRIX OF
                                                    !! FUNCTIONAL DERIVATIVE.
REAL    :: SL(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3))  !! TOTAL SOURCE FUNCTION.
REAL    :: SL_BOT(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)) !! BOTTOM SOURCE FUNCTION.

INTEGER :: MIJ(SIZE(FL3,1))    !! LAST FREQUENCY INDEX OF PROGNOSTIC PART.
INTEGER :: JU(SIZE(FL3,1))     !! U10 TABLE INDEX.
REAL    :: EMEAN(SIZE(FL3,1))  !! TOTAL ENERGY
REAL    :: FMEAN(SIZE(FL3,1))  !! MEAN FREQUENCY
REAL    :: AKMEAN(SIZE(FL3,1)) !! MEAN WAVENUMBER BASED ON SQRT(1/K)-MOMENT
REAL    :: TEMP(SIZE(FL3,1),SIZE(FL3,3))
REAL    :: DELFL(SIZE(FL3,3))
LOGICAL :: LLWS(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3))

REAL    :: F1MEAN(SIZE(FL3,1))    !! MEAN FREQUENCY BASED ON F-MOMENT
REAL    :: XKMEAN(SIZE(FL3,1))    !! MEAN WAVENUMBER BASED ON SQRT(K)-MOMENT
REAL    :: EMEANWS(SIZE(FL3,1))   !! TOTAL WINDSEA ENERGY
REAL    :: FMEANWS(SIZE(FL3,1))   !! MEAN WINDSEA FREQUENCY
REAL    :: SPRD(SIZE(FL3,1),SIZE(FL3,2)) !! SPREAD FUNCTIONS.

REAL    :: SPOS(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)) !! POSITIVE SINPUT
REAL    :: SMIN(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)) !! NEGATIVE SINPUT
REAL    :: SSOURCE(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)) !! SOURCE TERMS CONTRIBUTING
                                                        !! SURFACE WAVE FLUXES. 

! MODULATION OF SOURCE TERM BY IMPLICIT FACTOR IN THE CALCULATION
! OF THE SURFACE WAVE FLUXES To THE OCEANS.
LOGICAL, PARAMETER :: LLIMPFLX=.FALSE.
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CALCULATE ROUGHNESS LENGTH AND FRICTION VELOCITIES.                   !
!        ---------------------------------------------------                   !

KL = SIZE(FL3,2)
ML = SIZE(FL3,3)

DO  K = 1,KL
   SPRD(:,K) = MAX(0.,COS(TH(K)-UDIR(:)))**2   !! cosine spreading.
END DO

JU(:) = MIN(JUMAX, MAX(NINT(U10(:)/DELU),1))   !! u10 tabel index.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN PARAMETERS.                                              !
!        ------------------------                                              !

!     2.1 Total wave energy.                                                   !
!         ------------------                                                   !

CALL TOTAL_ENERGY (FL3, EMEAN)

!     2.2 Reduce wave energy if larger than depth limited wave height.         !
!         ------------------------------------------------------------         !

IF (WAVE_BREAKING_RUN) THEN
  CALL SDEPTHLIM(DEPTH, EMEAN, FL3) 
ENDIF

!     2.3 Mean frequencies and wave numbers.                                   !
!         ---------------------------------                                    !

CALL FEMEAN (FL3, EMEAN, FMEAN)
CALL TM1_TM2_PERIODS (FL3, EMEAN, TM1=F1MEAN)
F1MEAN(:) = 1./F1MEAN(:)
IF (SHALLOW_RUN) THEN
   CALL WM1_WM2_WAVENUMBER (FL3, EMEAN, WM1=AKMEAN, WM2=XKMEAN, IN=INDEP)
ELSE
   CALL WM1_WM2_WAVENUMBER (FL3, EMEAN, WM1=AKMEAN, WM2=XKMEAN)
END IF

! -------------------------------------------------------------------------!
!                                                                          !
!     3. COMPUTATION OF SOURCE FUNCTIONS AND DERIVATIVES.                  !
!        ------------------------------------------------                  !
CALL AIRSEA (U10, TAUW, USTAR, Z0)

IF (IPHYS .EQ. 1 ) THEN
   CALL SINPUT_ARD (FL3, SL, SPOS, FL, USTAR, UDIR, Z0, ROAIRN, WSTAR,     &
&                   INDEP, LLWS)
ELSE
   CALL SINPUT     (FL3, SL, SPOS, FL, USTAR, UDIR, Z0, ROAIRN, WSTAR,     &
&                   INDEP, LLWS)
ENDIF

CALL TOTAL_ENERGY (FL3, EMEANWS, LLWS)
CALL FEMEAN (FL3, EMEANWS, FMEANWS, LLWS)
CALL FRCUTINDEX (FMEAN, FMEANWS, USTAR, MIJ)

CALL STRESSO (FL3, SPOS, USTAR, UDIR, Z0, MIJ, TAUW, PHIAW, INDEP)

! re-evalute the input
CALL AIRSEA (U10, TAUW, USTAR, Z0)
CALL IMPHFTAIL (MIJ, INDEP, FL3)

IF (IPHYS .EQ. 1 ) THEN
   CALL SINPUT_ARD (FL3, SL, SPOS, FL, USTAR, UDIR, Z0, ROAIRN, WSTAR,     &
&                   INDEP, LLWS)
   IF (LCFLX) SMIN(:,:,:) = SL(:,:,:) - SPOS(:,:,:)
ELSE
   CALL SINPUT     (FL3, SL, SPOS, FL, USTAR, UDIR, Z0, ROAIRN, WSTAR,     &
&                   INDEP, LLWS)
   IF (LCFLX) SMIN(:,:,:) = SL(:,:,:) - SPOS(:,:,:)
ENDIF

CALL STRESSO (FL3, SPOS, USTAR, UDIR, Z0, MIJ, TAUW, PHIAW, INDEP)


IF (IPHYS .EQ. 1 ) THEN
  CALL SDISSIP_ARD (FL3, SL, FL, USTAR, UDIR, ROAIRN, INDEP)
ELSE
  CALL SDISSIP     (FL3, SL, FL, EMEAN, F1MEAN, XKMEAN, INDEP)
ENDIF

CALL SNONLIN (FL3, SL, FL, DEPTH, AKMEAN)

IF (SHALLOW_RUN) THEN
   SL_BOT = 0.
   CALL SBOTTOM (FL3, SL_BOT, FL, DEPTH, INDEP)
END IF

! COMPUTE WAVE SURFACE FLUXES
IF (LCFLX) THEN
!!!!!!  SL must only contain contributions contributed to surface fluxes into the oceans
   IF (LLIMPFLX) THEN     !!   MODULATE SL BY IMPLICIT FACTOR
      DELT = IDELT
      DO M=1,ML
         DO K=1,KL
            TEMP(:,1) = 1.0/ MAX(1.,1.-DELT*FL(:,K,M))
            SSOURCE(:,K,M) = (SL(:,K,M)-SMIN(:,K,M))*TEMP(:,1)
         END DO
      END DO
      CALL WNFLUXES (SSOURCE, SL_BOT, USTAR, UDIR, PHIAW, DEPTH, INDEP, MIJ)

   ELSE       !!  NO MODULATION OF SL BY IMPLICIT FACTOR

      SSOURCE(:,:,:) = SL(:,:,:)-SMIN(:,:,:)
      CALL WNFLUXES (SSOURCE, SL_BOT, USTAR, UDIR, PHIAW, DEPTH, INDEP, MIJ)
   ENDIF
ENDIF

IF (SHALLOW_RUN) SL(:,:,:) = SL(:,:,:) + SL_BOT(:,:,:)
IF (WAVE_BREAKING_RUN) THEN
   CALL SFBRK (FL3, SL, FL, EMEAN, FMEAN, DEPTH)
END IF

IF (PHILLIPS_RUN) THEN
   call source_phillips (sl, ustar, udir, depth, indep)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTATION OF NEW SPECTRA.                                           !
!        ---------------------------                                           !
!                                                                              !
!       INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE             !
!       FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.                    !

!     4.2 INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE           !
!         FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.                  !
!         ----------------------------------------------------------           !

DELT = IDELT
DELFL = 5.0E-07*G/FR**4*DELT
DO M=1,ML
   TEMP(:,2) = USTAR*DELFL(M)*MAX(FMEANWS,FMEAN)
   DO K=1,KL
      TEMP(:,1) = DELT*SL(:,K,M)/ MAX(1.,1.-DELT*FL(:,K,M))
      TEMP(:,3) = MIN(ABS(TEMP(:,1)),TEMP(:,2))
      FL3(:,K,M) = FL3(:,K,M)+SIGN(TEMP(:,3),TEMP(:,1))
      FL3(:,K,M) = MAX(FL3(:,K,M),FLMINFR(JU(:),M)*SPRD(:,K))
   END DO
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.                 !
!        -----------------------------------------------------                 !

!    5.1 COMPUTE MEAN PARAMETERS.                                              !
!        ------------------------                                              !

CALL TOTAL_ENERGY (FL3, EMEAN)
CALL FEMEAN (FL3, EMEAN, FMEAN)
CALL TOTAL_ENERGY (FL3, EMEANWS, LLWS)
CALL FEMEAN (FL3, EMEANWS, FMEANWS, LLWS)

!     5.2 COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.         !
!         ------------------------------------------------------------         !

CALL FRCUTINDEX (FMEAN, FMEANWS, USTAR, MIJ)


!     5.3 COMPUTE TAIL ENERGY RATIOS AND MERGE TAIL INTO SPECTRA.              !
!         -------------------------------------------------------              !

CALL IMPHFTAIL (MIJ, INDEP, FL3)

END SUBROUTINE IMPLSCH

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_SOURCE (DEPTH)

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
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: DEPTH (:)      !! WATER DEPTH [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. WEIGHT OF NON-LINEAR INTERACTION.                                     !
!        ---------------------------------                                     !

CALL NLWEIGT
IF (ITEST.GE.2) WRITE (IU06,*) '    SUB. PREPARE_SOURCE: NLWEIGT DONE'

CALL INIT_SNONLIN
IF (ITEST.GE.2) WRITE (IU06,*) '    SUB. PREPARE_SOURCE: INIT_SNONLIN'

IF (SHALLOW_RUN) THEN
   CALL MAKE_SHALLOW_SNL (DEPTH)
   IF (ITEST.GE.2) THEN
      WRITE (IU06,*) '    SUB. PREPARE_SOURCE: MAKE_SHALLOW_SNL DONE '
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TABLES AND PRE_COMPUTED CONSTANT FOR WAM PHYSICS.                     !
!        -------------------------------------------------                     !

IF (IPHYS.EQ.1) THEN
   CALL TABU_SWELLFT
   IF (ITEST.GE.2) WRITE (IU06,*) '    SUB. PREPARE_SOURCE: TABU_SWELLFT DONE'

   CALL INIT_SDISSP_ARD
   IF (ITEST.GE.2) WRITE (IU06,*) '    SUB. PREPARE_SOURCE: INIT_SDISSP_ARD DONE'
ENDIF

CALL INIT_X0TAUHF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ALLOCATE ARRAYS FOR FLUXES.                                           !
!        ---------------------------                                           !

IF (ALLOCATED(PHIOC)) DEALLOCATE(PHIOC)
ALLOCATE (PHIOC(1:NIJL-NIJS+1))
PHIOC(:) = 0.0
IF (ALLOCATED(PHIAW)) DEALLOCATE(PHIAW)
ALLOCATE (PHIAW(1:NIJL-NIJS+1))
PHIAW(:) = 0.0
IF (ALLOCATED(TAUOC_X)) DEALLOCATE(TAUOC_X)
ALLOCATE (TAUOC_X(1:NIJL-NIJS+1))
TAUOC_X(:) = 0.0
IF (ALLOCATED(TAUOC_Y)) DEALLOCATE(TAUOC_Y)
ALLOCATE (TAUOC_Y(1:NIJL-NIJS+1))
TAUOC_Y(:) = 0.0
IF (ALLOCATED(PHIBOT)) DEALLOCATE(PHIBOT)
ALLOCATE (PHIBOT(1:NIJL-NIJS+1))
PHIBOT(:) = 0.0
IF (ALLOCATED(TAUBOT_X)) DEALLOCATE(TAUBOT_X)
ALLOCATE (TAUBOT_X(1:NIJL-NIJS+1))
TAUBOT_X(:) = 0.0
IF (ALLOCATED(TAUBOT_Y)) DEALLOCATE(TAUBOT_Y)
ALLOCATE (TAUBOT_Y(1:NIJL-NIJS+1))
TAUBOT_Y(:) = 0.0

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

DO M = MFRSTLW,MLSTHG
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
WRITE(IU06,'(1X,8F10.7)') FRH(1:KFRH)

END SUBROUTINE PRINT_SOURCE_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE AIRSEA (UTOP, TAUW, USTAR, Z0)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   AIRSEA - - COMPUTATION OF TOTAL STRESS AND ROUGHNESS LENGTH SCALE.         !
!                                                                              !
!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE TOTAL STRESS.                                                  !
!                                                                              !
!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)+RNUAIRM/USTAR

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: UTOP(:)       !! WIND SPEED AT REFERENCE LEVEL XNLEV.
REAL,    INTENT(IN)    :: TAUW(:)       !! WAVE STRESS.
REAL,    INTENT(OUT)   :: USTAR(:)      !! FRICTION VELOCITY.
REAL,    INTENT(OUT)   :: Z0(:)         !! ROUGHNESS LENGTH.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------   

INTEGER, PARAMETER :: NITER=15

REAL, PARAMETER :: TWOXMP1=3.0
REAL, PARAMETER :: EPSUS = 1.0E-6

!     *ACD*       COEFFICIENTS FOR SIMPLE CD(U10) RELATION
!     *BCD*       CD = ACD + BCD*U10
REAL, PARAMETER :: ACD=8.0E-4
REAL, PARAMETER :: BCD=8.0E-5

INTEGER :: IJ, ITER

! REAL :: ALPHA
REAL :: XLOGXL, ALPHAOG, XKUTOP, XOLOGZ0
REAL :: USTOLD, TAUOLD, TAUNEW, X, F, DELF
REAL :: USTM1, Z0TOT, Z0CH, Z0VIS, ZZ

! ----------------------------------------------------------------------

XLOGXL=LOG(XNLEV)
ALPHAOG=ALPHA/G

DO IJ=1,SIZE(USTAR)
   XKUTOP = XKAPPA*UTOP(IJ)
   USTOLD = UTOP(IJ)*SQRT(ACD+BCD*UTOP(IJ))
   TAUOLD = MAX(USTOLD**2,TAUW(IJ)+EPS1)
   USTAR(IJ) = SQRT(TAUOLD)
   USTM1 = 1.0/MAX(USTAR(IJ),EPSUS)

   DO ITER=1,NITER
      X = TAUW(IJ)/TAUOLD
!      Z0CH = ALPHAOG*TAUOLD/SQRT(1.0-X)
      Z0CH = ALPHAOG*TAUOLD/SQRT(MAX(1.0-X,EPS1))
      Z0VIS = RNUAIRM*USTM1
      Z0TOT = Z0CH+Z0VIS

      XOLOGZ0= 1.0/(XLOGXL-LOG(Z0TOT))
      F = USTAR(IJ)-XKUTOP*XOLOGZ0
      ZZ = USTM1*(Z0CH*(2.0-TWOXMP1*X)/(1.0-X)-Z0VIS)/Z0TOT
      DELF= 1.0-XKUTOP*XOLOGZ0**2*ZZ

      USTAR(IJ) = USTAR(IJ)-F/DELF
      TAUNEW = MAX(USTAR(IJ)**2,TAUW(IJ)+EPS1)
      USTAR(IJ) = SQRT(TAUNEW)
      IF (TAUNEW.EQ.TAUOLD) EXIT
      USTM1 = 1.0/MAX(USTAR(IJ),EPSUS)
      TAUOLD = TAUNEW
   ENDDO

   Z0(IJ)=Z0CH

ENDDO

END SUBROUTINE AIRSEA

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MAKE_SHALLOW_SNL (DEPTH)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      MAKE_SHALLOW_SNL - COMPUTE THE NONLINEAR TRANSFER FUNCTION COEFFICIENTS !
!                         FOR SHALLOW WATER.                                   !
!                                                                              !
!      P. JANSSEN     ECMWF  JUNE 2005                                         !
!      H. GUNTHER     HZG    JANUARY 2015  CYCLE_4.5.4                         !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!          NONE                                                                !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!     -----------                                                              !

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,       ONLY: IU06, ITEST

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: DEPTH (:)      !! WATER DEPTH [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, M
REAL    :: D, OM, XK

REAL, PARAMETER  :: ENH_MAX=10. !! MAXIMUM COEFFICIENT.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.  THE FIRST CALL TO WAVEMDL PERFORMS INITIALIZATION.                   !
!         --------------------------------------------------                   !

IF (.NOT.ALLOCATED(ENH)) ALLOCATE(ENH(SIZE(DEPTH),ML+4))

IF (.NOT.ALLOCATED(ENH)) ALLOCATE(ENH(SIZE(DEPTH),MLSTHG))

IF (ISNONLIN.NE.0) THEN
   DO M = 1,ML
      DO IJ = 1, SIZE(DEPTH)
         D = DEPTH(IJ)
         OM = ZPI*FR(M)
         XK = AKI(OM,D)
         ENH(IJ,M) = MIN(ENH_MAX,TRANSF(XK,D))
      END DO
    END DO
!       NOTE THAT FR IS NOT DEFINED FOR M>ML.
    DO M = ML+1, MLSTHG
       DO IJ = 1, SIZE(DEPTH)
          D = DEPTH(IJ)
          OM = ZPI*FR(ML)*CO**(M-ML)
          XK = AKI(OM,D)
          ENH(IJ,M) = MIN(ENH_MAX,TRANSF(XK,D))
       END DO
    END DO
ENDIF

END SUBROUTINE MAKE_SHALLOW_SNL

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
!REAL, PARAMETER :: DELPHI1 = -11.48   !!
!REAL, PARAMETER :: DELPHI2 = 33.56    !!

INTEGER :: KLP1, IC, KH, KLH, K, KS, ISG, K1, K11, K2, K21
INTEGER :: M, IKN, I, ISP, ISM
INTEGER, ALLOCATABLE :: JA1(:,:)
INTEGER, ALLOCATABLE :: JA2(:,:)

REAL    :: DELPHI1
REAL    :: DELPHI2
REAL    :: DELTHA, CL1, CL2, AL11, AL12, CH, CL1H, CL2H
REAL    :: F1P1, FRG, FLP, FLM, FKP, FKM, XF, COSTH3, COSTH4
REAL,    ALLOCATABLE :: FRLON(:)

INTEGER, EXTERNAL :: JAFU

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. ALLOCATE ARRAYS.                                                      !
!        ----------------                                                      !

F1P1 = LOG10(CO)
ISP = INT(LOG10(1.+ALAMD)/F1P1+.000001)
ISM = FLOOR(LOG10(1.-ALAMD)/F1P1+.0000001)

MFRSTLW = 1+ISM
MLSTHG = ML-ISM

KFRH=-ISM+ISP+2

ALLOCATE(JA1(KL,2))
ALLOCATE(JA2(KL,2))
ALLOCATE(FRLON(MFRSTLW:ML+KFRH))

IF (.NOT.ALLOCATED (IKP ))  ALLOCATE (IKP (MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (IKP1))  ALLOCATE (IKP1(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (IKM ))  ALLOCATE (IKM (MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (IKM1))  ALLOCATE (IKM1(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (K1W ))  ALLOCATE (K1W (KL,2))
IF (.NOT.ALLOCATED (K2W ))  ALLOCATE (K2W (KL,2))
IF (.NOT.ALLOCATED (K11W))  ALLOCATE (K11W(KL,2))
IF (.NOT.ALLOCATED (K21W))  ALLOCATE (K21W(KL,2))
IF (.NOT.ALLOCATED (AF11))  ALLOCATE (AF11(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (FKLAP ))  ALLOCATE (FKLAP(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (FKLAP1))  ALLOCATE (FKLAP1(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (FKLAM ))  ALLOCATE (FKLAM(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (FKLAM1))  ALLOCATE (FKLAM1(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (FRH))  ALLOCATE(FRH(KFRH))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTATION FOR ANGULAR GRID.                                         !
!        -----------------------------                                         !

XF      = ((1.+ALAMD)/(1.-ALAMD))**4
COSTH3  = (1.+2.*ALAMD+2.*ALAMD**3)/(1.+ALAMD)**2
DELPHI1 = -180./PI*ACOS(COSTH3)
COSTH4  = SQRT(1.-XF+XF*COSTH3**2)
DELPHI2 = 180./PI*ACOS(COSTH4)

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

CL1  = CL1-INT(CL1)
CL2  = CL2-INT(CL2)
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

DO M=0,MFRSTLW,-1
   FRLON(M)=FRLON(M+1)/CO
ENDDO
DO M=ML+1,ML+KFRH
   FRLON(M) = CO*FRLON(M-1)
ENDDO

DO M = MFRSTLW,MLSTHG
   FRG = FRLON(M)
   AF11(M) = CON * FRG**11
   FLP = FRG*(1.+ALAMD)
   FLM = FRG*(1.-ALAMD)
   IKN = M+ISP
   IKP(M) = IKN
   FKP = FRLON(IKP(M))
   IKP1(M) = IKP(M)+1
   FKLAP(M) = (FLP-FKP)/(FRLON(IKP1(M))-FKP)
   FKLAP1(M) = 1.-FKLAP(M)
   IF (FRLON(MFRSTLW).GE.FLM) THEN
      IKM(M) = 1
      IKM1(M) = 1
      FKLAM(M) = 0.
      FKLAM1(M) = 0.
   ELSE
      IKN = M+ISM
      IKM(M) = IKN
      FKM = FRLON(IKM(M))
      IKM1(M) = IKM(M)+1
      FKLAM(M) = (FLM-FKM)/(FRLON(IKM1(M))-FKM)
      FKLAM1(M) = 1.-FKLAM(M)
      IF (IKN.LT.MFRSTLW) THEN
         IKM(M) = 1
         FKLAM1(M) = 0.
      ENDIF
   END IF
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTE TAIL FREQUENCY RATIOS.                                        !
!        ------------------------------                                        !

DO I=1,KFRH
   M = ML+I-1
   FRH(I) = (FRLON(ML)/FRLON(M))**5
END DO

END SUBROUTINE NLWEIGT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE INIT_SNONLIN

! ---------------------------------------------------------------------------- !
!                                                                              !
!    INIT_SNONLIN - INITIALISE ALL FREQUENCY DEPENDENT ARRAYS USED BY SNONLIN  !
!                                                                              !
!     J. BIDLOT   ECMWF  MAY 2012                                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       USED TO BE IN SNONLIN BUT NOW IT IS ONLY COMPUTED ONCE.                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
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

INTEGER :: ICOUNT, IRCOUNT
INTEGER :: MC, MP, MP1, MM, MM1, IC, IP, IP1, IM , IM1, ITEMP

REAL :: ALPH, FRR
REAL :: FFACP, FFACP1, FFACM, FFACM1, FTAIL, FKLAMP, FKLAMP1
REAL :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAPA2, FKLAPB2
REAL :: FKLAP12, FKLAP22, FKLAMM, FKLAMM1, FKLAMMA, FKLAMMB
REAL :: FKLAMM2, FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
REAL :: GW1, GW2, GW3, GW4, GW5, GW6, GW7, GW8

!     INLINE FUNCTION (PIERSON-MOSKOWITZ SMOOTH CUT-OFF)
!     X == FR(1)/FREQUENCY
REAL :: EPMMA, X
EPMMA(X) = EXP(-MIN(1.25*X**4,50.))*(X**5)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. FRONT SPECTRAL TAIL REDUCTION COEFFICIENTS

IF(.NOT.ALLOCATED(FTRF)) ALLOCATE(FTRF(MFRSTLW:1))
ALPH = 1./EPMMA(1.)
FRR = 1.
DO MC=1,MFRSTLW,-1
   FTRF(MC)=ALPH*EPMMA(FRR)
   FRR=FRR*CO
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. WORK ARRAYS STORING THE DIFFERENT INDICES AND COEFFICIENTS

IF(.NOT.ALLOCATED(INLCOEF)) ALLOCATE(INLCOEF(NINL,1:MLSTHG))
IF(.NOT.ALLOCATED(RNLCOEF)) ALLOCATE(RNLCOEF(NRNL,1:MLSTHG))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. FREQUENCY LOOP.
!        ---------------

DO MC=1,MLSTHG
   MP  = IKP (MC)
   MP1 = IKP1(MC)
   MM  = IKM (MC)
   MM1 = IKM1(MC)
   FFACP  = 1.
   FFACP1 = 1.
   FFACM  = 1.
   FFACM1 = 1.
   FTAIL  = 1.
   IC  = MC
   IP  = MP
   IP1 = MP1
   IM  = MM
   IM1 = MM1
!       LOW FREQUENCY FRONT TAIL
   IF (IM.LT.1) THEN
      FFACM = FTRF(IM)
      IM = 1
      IF (IM1.LT.1) THEN
         FFACM1 = FTRF(IM1)
         IM1 = 1
      ENDIF
   ENDIF
!       HIGH FREQUENCY TAIL
   IF (IP1.GT.ML) THEN
! Quick fix from Deborah
      ITEMP=IP1-ML+1
      IF(ITEMP .GT. SIZE(FRH))THEN
         ITEMP=SIZE(FRH)
      ENDIF
!         FFACP1 = FRH(IP1-ML+1)
      FFACP1 = FRH(ITEMP)

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

   ICOUNT=1
   INLCOEF(ICOUNT,MC) = IC
   ICOUNT=ICOUNT+1
   INLCOEF(ICOUNT,MC) = IP
   ICOUNT=ICOUNT+1
   INLCOEF(ICOUNT,MC) = IP1
   ICOUNT=ICOUNT+1
   INLCOEF(ICOUNT,MC) = IM
   ICOUNT=ICOUNT+1
   INLCOEF(ICOUNT,MC) = IM1

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
   IRCOUNT=1
   RNLCOEF(IRCOUNT,MC) = FTAIL
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW1
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW3
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW4
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMPA
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMPB
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMP2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMP1
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAPA2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAPB2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAP12
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAP22

   FKLAMM  = FKLAM(MC)
   FKLAMM1 = FKLAM1(MC)
   GW6 = FKLAMM1*FFACM*DAL2
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
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW5
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW6
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW7
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW8
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMMA
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMMB
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMM2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMM1
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMA2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMB2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAM12
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAM22

ENDDO

IF(ICOUNT.NE.NINL) THEN
   WRITE(IU06,*) '*************************************'
   WRITE(IU06,*) 'ERROR IN INISNONLIN : ICOUNT NE NINL'
   WRITE(IU06,*) 'ICOUNT= ',ICOUNT
   WRITE(IU06,*) 'NINL= ',NINL
   WRITE(IU06,*) '*************************************'
   CALL ABORT1
ENDIF
IF(IRCOUNT.NE.NRNL) THEN
   WRITE(IU06,*) '*************************************'
   WRITE(IU06,*) 'ERROR IN INISNONLIN : IRCOUNT NE NRNL'
   WRITE(IU06,*) 'IRCOUNT= ',IRCOUNT
   WRITE(IU06,*) 'NRNL= ',NRNL
   WRITE(IU06,*) '*************************************'
   CALL ABORT1
ENDIF

END SUBROUTINE INIT_SNONLIN

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
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL             !
!     H. GUENTHER GKSS  FEBRUARY 2002       FT 90                              !
!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON AKMEAN          !
!                                       AND FMEAN.                             !
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

REAL, INTENT(IN)     :: F (:, :, :) !! SPECTRUM.
REAL, INTENT(INOUT)  :: SL(:, :, :) !! TOTAL SOURCE FUNCTION ARRAY
REAL, INTENT(INOUT)  :: FL(:, :, :) !! DIAGONAL MATRIX OF FUNCTIONAL
REAL, INTENT(IN)     :: EMEAN (:)   !! TOTAL ENERGY
REAL, INTENT(IN)     :: FMEAN (:)   !! MEAN FREQUENCY BASED ON 1. MOMENT
REAL, INTENT(IN)     :: AKMEAN(:)   !! MEAN WAVE NUMBER BASED ON SQRT(K) MOMENT
INTEGER, INTENT(IN)  :: INDEP (:)   !! DEPTH TABLE INDEX

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: CDIS = 1.33
REAL, PARAMETER :: DELTA = 0.5
! REAL, PARAMETER :: CDIS = 2.1          !! DISSIPATION CONSTANT
REAL, PARAMETER :: CONSD = -CDIS*ZPI**9/G**4
REAL, PARAMETER :: CONSS = -CDIS*ZPI
! REAL, PARAMETER :: DELTA = 0.6         !! WEIGHT LINEAR, QUADRATIC PART.


INTEGER :: K, M,IJ
REAL ::  FAC, SDISS

REAL, DIMENSION(SIZE(F,1)) :: TEMP1
REAL, DIMENSION(SIZE(F,1)) :: SDS
REAL, DIMENSION(SIZE(F,1)) :: CM

! ---------------------------------------------------------------------------- !
!                                                                              !
!    1.  ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE        !
!        FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.                          !
!        --------------------------------------------------------------        !

IF (SHALLOW_RUN) THEN

   SDS(:) = CONSS*FMEAN(:)*EMEAN(:)**2*AKMEAN(:)**4
   FRES: DO M = 1,SIZE(F,3)
      FAC = ZPI*FR(M)
      TEMP1(:) = TFAK(INDEP(:),M)/AKMEAN(:)
      TEMP1(:) = SDS(:) * ((1.-DELTA)*TEMP1(:) +  DELTA*TEMP1(:)**2)
      CM(:) = TFAK(INDEP(:),M)/FAC

      DIRS: DO K = 1,SIZE(F,2)
         DO IJ = 1, SIZE(F,1)
            SDISS = TEMP1(IJ)*F(IJ,K,M)
            SL(IJ,K,M) = SL(IJ,K,M) + SDISS
            FL(IJ,K,M) = FL(IJ,K,M) + TEMP1(IJ)
         END DO
      END DO DIRS
   END DO FRES

ELSE

   SDS(:) = CONSD*EMEAN(:)**2*FMEAN(:)**9
   FRED: DO M = 1,SIZE(F,3)
      TEMP1(:) = (FR(M)/FMEAN(:))**2
      TEMP1(:) = SDS(:) * ((1.-DELTA)*TEMP1(:) + DELTA*TEMP1(:)**2)
      CM(1)  = ZPI*FR(M)/G

      DIRD: DO K = 1,SIZE(F,2)
         DO IJ = 1, SIZE(F,1)
            SDISS = TEMP1(IJ)*F(IJ,K,M)
            SL(IJ,K,M) = SL(IJ,K,M) + SDISS
            FL(IJ,K,M) = FL(IJ,K,M) + TEMP1(IJ)
         END DO
      END DO DIRD
   END DO FRED

END IF

END SUBROUTINE SDISSIP

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SDISSIP_ARD (F, SL, FL, USTAR, UDIR, ROAIRN, INDEP)

! ---------------------------------------------------------------------------- !
!**** *SDISSIP_ARD* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: F (:, :, :)    !! SPECTRUM.
REAL,    INTENT(OUT)   :: SL(:, :, :)    !! TOTAL SOURCE FUNCTION ARRAY
REAL,    INTENT(OUT)   :: FL(:, :, :)    !! DIAGONAL MATRIX OF FUNCTIONAL
                                         !! DERIVATIVE
REAL,    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY.
REAL,    INTENT(IN)    :: UDIR (:)       !! WIND DIRECTION.
REAL,    INTENT(IN)    :: ROAIRN(:)      !! AIR DENSITY
INTEGER, INTENT(IN)    :: INDEP(:)       !! DEPTH TABLE INDEX.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------

INTEGER :: IJ, K, M, M2, K2, KK, KLD
INTEGER, DIMENSION(SIZE(F,2)) :: KKD

REAL :: XK(SIZE(F,1),SIZE(F,3))
REAL :: TPIINV, TPIINVH, TMP01, TMP03
REAL :: EPSR
REAL :: ROG
REAL :: SSDSC6M1
REAL :: SIG(SIZE(F,3))
REAL :: SSDSC2_SIG(SIZE(F,3))
REAL :: FACTURB(SIZE(F,1))
REAL :: FACSAT(SIZE(F,1),SIZE(F,3))
REAL :: FACWTRB(SIZE(F,1),SIZE(F,3))
REAL :: TEMP1(SIZE(F,1),SIZE(F,3))
REAL :: BTH0(SIZE(F,1),SIZE(F,3))  !saturation spectrum 
REAL :: BTH(SIZE(F,1),SIZE(F,2),SIZE(F,3))  !saturation spectrum 

REAL, DIMENSION(SIZE(F,1),SIZE(F,2),SIZE(F,3)) :: D

!!! the following 3 arrays are only used when  LLSSDSC3
!!! is not used, there should be a way to save the memory
REAL, DIMENSION(SIZE(F,1),SIZE(F,2),SIZE(F,3)) :: SCUMUL 
REAL, DIMENSION(SIZE(F,1),SIZE(F,2),SIZE(F,3)) :: RENEWALFREQ
REAL, DIMENSION(SIZE(F,1),0:SIZE(F,2)/2,SIZE(F,3)) :: WCUMUL

LOGICAL :: LLSSDSC3,  LLSSDSC5

! ----------------------------------------------------------------------

! INITIALISATION

EPSR=SQRT(SDSBR)

TPIINV = 1.0/ZPI
TPIINVH= 0.5*TPIINV

KLD=SIZE(F,2)/2

ROG = ROWATER*G

LLSSDSC3=(SSDSC3.NE.0.0)

TMP03 = 1.0/(SDSBR*MICHE)

LLSSDSC5=(SSDSC5.NE.0.0)

SSDSC6M1=1.-SSDSC6

DO M = 1,SIZE(F,3)
  SIG(M) = ZPI*FR(M)
  SSDSC2_SIG(M)=SSDSC2*SIG(M)
END DO

IF (SHALLOW_RUN) THEN
  DO M = 1,SIZE(F,3)
    DO IJ = 1,SIZE(F,1)
       XK(IJ,M) = TFAK(INDEP(IJ),M)
       FACSAT(IJ,M) = XK(IJ,M)**3*TPIINV*TCGOND(INDEP(IJ),M)
    ENDDO
  ENDDO
ELSE
  DO M= 1,SIZE(F,3)
    DO IJ = 1,SIZE(F,1)
      XK(IJ,M) = (SIG(M)**2)/G
      FACSAT(IJ,M) = XK(IJ,M)**3*TPIINVH*G/SIG(M)
    ENDDO
  ENDDO
ENDIF


! COMPUTE SATURATION SPECTRUM
BTH0(:,:) = 0.0
BTH(:,:,:)=0.0

DO M = 1,SIZE(F,3)
  DO K = 1,SIZE(F,2)
    ! integrates in directional sector
    DO K2 = 1,NSDSNTH*2+1
      KK=INDICESSAT(K,K2)
      DO IJ = 1,SIZE(F,1)
         BTH(IJ,K,M) = BTH(IJ,K,M) + SATWEIGHTS(K,K2)*F(IJ,KK,M)
      ENDDO
    ENDDO
    DO IJ = 1,SIZE(F,1)
      BTH(IJ,K,M)=BTH(IJ,K,M)*FACSAT(IJ,M)
      BTH0(IJ,M)=MAX(BTH0(IJ,M),BTH(IJ,K,M))
    ENDDO
  ENDDO
ENDDO


! SATURATION TERM

DO  M = 1,SIZE(F,3)
  DO IJ = 1,SIZE(F,1)
    TEMP1(IJ,M)=SSDSC6*(MAX(0.,BTH0(IJ,M)*TMP03-SSDSC4))**IPSAT
  ENDDO
ENDDO

DO  M = 1,SIZE(F,3)
  DO K = 1,SIZE(F,2)
    DO IJ = 1,SIZE(F,1)
      D(IJ,K,M)= SSDSC2_SIG(M)*(TEMP1(IJ,M)+SSDSC6M1*(MAX(0.,BTH(IJ,K,M)*TMP03-SSDSC4))**IPSAT)
    ENDDO
  ENDDO
ENDDO


! CUMULATIVE TERM
IF (LLSSDSC3) THEN

  DO M2 = 1,SIZE(F,3)-NDIKCUMUL
    DO IJ = 1,SIZE(F,1)
      IF(BTH0(IJ,M2).GT.SDSBR) THEN
        TEMP1(IJ,M2)=1.0
      ELSE
        TEMP1(IJ,M2)=0.0
      ENDIF
    ENDDO
  ENDDO
  DO M2 = 1,SIZE(F,3)-NDIKCUMUL
    DO K2 = 1,SIZE(F,2)
      DO IJ = 1,SIZE(F,1)
        SCUMUL(IJ,K2,M2)=TEMP1(IJ,M2)*(MAX(SQRT(BTH(IJ,K2,M2))-EPSR,0.0))**2
      ENDDO
    ENDDO
  ENDDO

  DO M = NDIKCUMUL+1,SIZE(F,3)
    DO K = 1,SIZE(F,2)
      DO IJ = 1,SIZE(F,1)
        RENEWALFREQ(IJ,K,M)=0.0
      ENDDO
    ENDDO
  ENDDO


  DO M = NDIKCUMUL+1,SIZE(F,3)

    DO M2 = 1,M-NDIKCUMUL
      DO KK = 0,KLD
        DO IJ = 1,SIZE(F,1)
          WCUMUL(IJ,KK,M2)=CUMULW(INDEP(IJ),KK,M2,M)
        ENDDO
      ENDDO
    ENDDO

    DO K = 1,SIZE(F,2)
      ! Correction of saturation level for shallow-water kinematics
      ! Cumulative effect based on lambda   (breaking probability is
      ! the expected rate of sweeping by larger breaking waves)

      DO K2 = 1,SIZE(F,2)
        KKD(K2)=ABS(K2-K)
        IF(KKD(K2).GT.KLD) KKD(K2)=KKD(K2)-KLD
      ENDDO

      DO M2 = 1,M-NDIKCUMUL
        DO K2 = 1,SIZE(F,2)
          KK=KKD(K2)
          DO IJ = 1,SIZE(F,1)
          ! Integrates over frequencies M2 and directions K2 to 
          ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
          RENEWALFREQ(IJ,K,M)=RENEWALFREQ(IJ,K,M)+ WCUMUL(IJ,KK,M2)*SCUMUL(IJ,K2,M2)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO M = NDIKCUMUL+1,SIZE(F,3)
    DO K = 1,SIZE(F,2)
      DO IJ = 1,SIZE(F,1)
        D(IJ,K,M)= D(IJ,K,M) + RENEWALFREQ(IJ,K,M)
      ENDDO
    ENDDO
  ENDDO
ENDIF  ! LLSSDSC3


!     WAVE-TURBULENCE INTERACTION TERM
IF (LLSSDSC5) THEN
  TMP01 = 2.*SSDSC5/ROG
  DO IJ = 1,SIZE(F,1)
    FACTURB(IJ) = TMP01*ROAIRN(IJ)*USTAR(IJ)*USTAR(IJ)
   ENDDO
   DO M= 1, SIZE(F,3)
     DO IJ = 1,SIZE(F,1)
       FACWTRB(IJ,M) = SIG(M)*XK(IJ,M)*FACTURB(IJ)
     ENDDO
   DO K = 1,SIZE(F,2)
     DO IJ = 1,SIZE(F,1)
       D(IJ,K,M)= D(IJ,K,M)- FACWTRB(IJ,M)*COS(UDIR(IJ)-TH(K))
     ENDDO
   ENDDO
  ENDDO
ENDIF


! ADD ALL CONTRIBUTIONS TO SOURCE TERM
DO  M= 1, SIZE(F,3)
  DO K= 1, SIZE(F,2)
    DO IJ = 1,SIZE(F,1)
      SL(IJ,K,M) = SL(IJ,K,M)+D(IJ,K,M)*F(IJ,K,M)
      FL(IJ,K,M) = FL(IJ,K,M)+D(IJ,K,M)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE SDISSIP_ARD

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

REAL, PARAMETER :: ALPHA = 1.0

INTEGER :: M, K
REAL :: QB(SIZE(F,1)), BB(SIZE(F,1))
REAL :: SBR(SIZE(F,1)),DSBR(SIZE(F,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!   1. compute total dissipation rate according to Battjes-Janssen
!      -----------------------------------------------------------

BB = 8.*EMEAN/(GAMD*DEPTH)**2 !! (Hrms / Hmax)**2

CALL CMPQB (BB, QB)           !! fraction of breaking waves

QB = MIN(1.,QB)
SBR = -ALPHA*2.*FMEAN

WHERE (BB.LE.1.) SBR = SBR*QB/BB

WHERE (BB .LT. 1. .AND. ABS(BB - QB) .GT. 0.)
   DSBR = SBR * (1. - QB) / (BB - QB)
ELSEWHERE
   DSBR = 0.
ENDWHERE

DO M = 1,SIZE(F,3)
  DO K = 1,SIZE(F,2)
     SL(:,K,M) = SL(:,K,M) + SBR*F(:,K,M)
     FL(:,K,M) = FL(:,K,M) + DSBR
   END DO
END DO

END SUBROUTINE SFBRK

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CMPQB (BB, QB)

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

REAL    , INTENT(IN)  :: BB(:)     !! (RMS wave height/Maximum wave height)**2 !
REAL    , INTENT(OUT) :: QB(:)     !! Fraction of breaking waves

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL    :: Q0(SIZE(BB))

! ---------------------------------------------------------------------------- !
!                                                                              !

Q0 = 0.
WHERE (BB.GE.0.25) Q0 = (2.*SQRT(BB)-1.)**2

! ---------------------------------------------------------------------------- !
!                                                                              !

WHERE (BB.LT.1.)
   QB = Q0 - BB*(Q0 - exp((Q0-1.)/BB))/(BB-exp((Q0-1.)/BB))
ELSEWHERE
    QB = 1.
END WHERE

END SUBROUTINE CMPQB

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SINPUT (F, SL, SPOS, FL, USTAR, UDIR, Z0, ROAIRN, WSTAR, INDEP, LLWS)

! ---------------------------------------------------------------------------- !
!                                                                              !
! *** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.                         !
!                                                                              !
!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990                             !
!                                                                              !
!     OPTIMIZED BY : H. GUENTHER                                               !
!                                                                              !
!     MODIFIED BY :                                                            !
!       J-R BIDLOT NOVEMBER 1995                                               !
!       J-R BIDLOT FEBRUARY 1996-97                                            !
!       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY USING NEW    !
!                                              STRESS AND ROUGHNESS.           !
!       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR DENSITY         !
!                                 AND STABILITY-DEPENDENT WIND GUSTINESS       !
!       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE        !
!                                      RUNNING FASTER THAN THE WIND.           !
!       J-R BIDLOT JANUARY 2013: SHALLOW WATER FORMULATION.                    !
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
REAL,    INTENT(OUT)   :: SPOS(:, :, :)  !! POSITIVE SOURCE FUNCTION ARRAY
REAL,    INTENT(INOUT) :: FL(:, :, :)    !! DIAGONAL MATRIX OF FUNCTIONAL
                                         !! DERIVATIVE
REAL,    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY.
REAL,    INTENT(IN)    :: UDIR (:)       !! WIND DIRECTION.
REAL,    INTENT(IN)    :: Z0   (:)       !! ROUGHNESS LENGTH.
REAL,    INTENT(IN)    :: ROAIRN(:)      !! AIR DENSITY
REAL,    INTENT(IN)    :: WSTAR(:)       !! FREE CONVECTION VELOCITY SCLALE
INTEGER, INTENT(IN)    :: INDEP(:)       !! DEPTH TABLE INDEX.
LOGICAL, INTENT(OUT)   :: LLWS(:, :, :)  !! TRUE WHERE SINPUT IS POSITIVE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER, PARAMETER :: NSIN = 2
INTEGER :: K, M, IJ, ISIN
REAL  :: X, ZLOG, ZLOG2X, CONST3

REAL, DIMENSION(NSIN)      :: WSIN
REAL, DIMENSION(SIZE(F,3)) :: FAC, CONST
REAL, DIMENSION(SIZE(F,1)) :: CM
REAL, DIMENSION(SIZE(F,1)) :: XK
REAL, DIMENSION(SIZE(F,1)) :: SIG_N
REAL, DIMENSION(SIZE(F,1)) :: CNSN
REAL, DIMENSION(SIZE(F,1)) :: UFAC1
REAL, DIMENSION(SIZE(F,1)) :: UFAC2
REAL, DIMENSION(SIZE(F,1),NSIN) :: SIGDEV, USN, Z0N, UCN, ZCN
REAL, DIMENSION(SIZE(F,1),NSIN) :: XV1D

REAL, DIMENSION(SIZE(F,1),SIZE(F,2)) :: TEMP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. PRECALCULATED ANGULAR DEPENDENCE.                                     !
!        ---------------------------------                                     !

DO K = 1,SIZE(F,2)
   TEMP(:,K) = COS(TH(K)-UDIR(:))
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. PRECALCULATED FREQUENCY DEPENDENCE.                                   !
!        -----------------------------------                                   !

FAC(:) = ZPI*FR(:)
CONST(:) = FAC(:)*BETAMAX/XKAPPA**2/ROWATER
CONST3   = 2.*XKAPPA/BETAMAX*XKAPPA**2  ! SEE IDAMPING

LLWS(:,:,:) = .FALSE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.                         !
!        ---------------------------------------------                         !

IF (NSIN.GT.1) CALL WSIGSTAR (USTAR, Z0, WSTAR, SIG_N)

IF (NSIN.EQ.1) THEN
   WSIN(1) = 1.0
   SIGDEV(:,1) = 1.0
ELSE IF (NSIN.EQ.2) THEN
   WSIN(1) = 0.5
   WSIN(2) = 0.5
   SIGDEV(:,1) = 1.0-SIG_N(:)
   SIGDEV(:,2) = 1.0+SIG_N(:)
ELSE
   WRITE (IU06,*) '**************************************'
   WRITE (IU06,*) '*    FATAL ERROR                     *'
   WRITE (IU06,*) '*    ===========                     *'
   WRITE (IU06,*) '* IN SINPUT_JAN: NSIN > 2            *'
   WRITE (IU06,*) '* NSIN = ', NSIN
   WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
   WRITE (IU06,*) '*                                    *'
   WRITE (IU06,*) '**************************************'
   CALL ABORT1
ENDIF

IF (NSIN.EQ.1) THEN
   USN(:,1) = USTAR(:)
   Z0N(:,1) = Z0(:)
ELSE
   DO ISIN=1,NSIN
      USN(:,ISIN) = USTAR(:)*SIGDEV(:,ISIN)
      Z0N(:,ISIN) = Z0(:)
   ENDDO
ENDIF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. LOOP OVER FREQUENCIES.                                                !
!        ----------------------                                                !

FRE: DO M = 1,SIZE(F,3)

   IF (SHALLOW_RUN) THEN
      XK(:) = TFAK(INDEP(:),M)
      CM(:) = XK(:)/FAC(M)         !! INVERSE OF PHASE VELOCITY
      CNSN(:) = CONST(M)*FAC(M)**2/(G*XK(:))*ROAIRN(:)
   ELSE
      XK(:) = FAC(M)**2/G
      CM(:) = FAC(M)/G              !! INVERSE OF PHASE VELOCITY
      CNSN(:) = CONST(M)*ROAIRN(:)
   END IF

   DO ISIN = 1, NSIN
      UCN(:,ISIN) = USN(:,ISIN)*CM(:) + ZALP
      ZCN(:,ISIN) = LOG(XK(:)*Z0N(:,ISIN))
      XV1D(:,ISIN) = -1./(USN(:,ISIN)/XKAPPA*ZCN(:,ISIN)*CM(:))
   END DO

!     3.1 LOOP OVER DIRECTIONS.                                                !
!         ---------------------                                                !

   DIR: DO K = 1,SIZE(F,2)
      UFAC1(:) = 0.
      DO ISIN=1,NSIN
         DO IJ = 1,SIZE(F,1)
            IF (TEMP(IJ,K).GT.0.01) THEN
               X    = TEMP(IJ,K)*UCN(IJ,ISIN)
               ZLOG = ZCN(IJ,ISIN) + XKAPPA/X
               IF (ZLOG.LT.0.) THEN
                  ZLOG2X = ZLOG*ZLOG*X
                  UFAC1(IJ) = UFAC1(IJ) + WSIN(ISIN)*EXP(ZLOG)*ZLOG2X*ZLOG2X
                  LLWS(IJ,K,M) = .TRUE.
               END IF
            END IF
         END DO
      END DO

!     3.2 SWELL DAMPING.                                                       !
!         --------------                                                       !

      UFAC2(:) = CONST3*(TEMP(:,K)-XV1D(:,1))*UCN(:,1)**2
      IF (NSIN.EQ.2) THEN
         UFAC2(:) = WSIN(1)*UFAC2(:)+                                          &
&                   WSIN(2)*CONST3*(TEMP(:,K)-XV1D(:,2))*UCN(:,2)**2
      END IF

!     3.3 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.                     !
!         ------------------------------------------------                     !

      SPOS(:,K,M) = CNSN(:)*UFAC1(:)
      FL(:,K,M) = SPOS(:,K,M)+CNSN(:)*UFAC2(:)
      SPOS(:,K,M) = SPOS(:,K,M)*F(:,K,M)
      SL(:,K,M) = FL(:,K,M)*F(:,K,M)

   END DO DIR
END DO FRE

END SUBROUTINE SINPUT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SINPUT_ARD (F,SL,SPOS,FL,USTAR,UDIR,Z0,ROAIRN,WSTAR,INDEP,LLWS)

! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY :
!       J-R BIDLOT NOVEMBER 1995
!       J-R BIDLOT FEBRUARY 1996-97
!       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
!       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
!       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
!                                  USING NEW STRESS AND ROUGHNESS.
!       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
!                                 DENSITY AND STABILITY-DEPENDENT
!                                 WIND GUSTINESS
!       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE
!                                      RUNNING FASTER THAN THE WIND.
!       J-R BIDLOT JANUARY 2013: SHALLOW WATER FORMULATION.


!       L. AOUF    March 2011 : USE OF NEW DISSIPATION DEVELOPED BY ARDHUIN ET AL.2010

!       JEAN BIDLOT : ADAPTED TO ECWAM AND WAM.

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION.
!
!       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
!       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
!       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
!       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
!       FINDS:
!
!             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
!
!       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
!       LEVEL.

!     METHOD.                                                                  !
!     ------- 

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       WSIGSTAR.

!     REFERENCE.
!     ----------

!       P. JANSSEN, J.P.O., 1989.
!       P. JANSSEN, J.P.O., 1991
!       ARDHUIN ET AL. 2010

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: F (:, :, :)    !! SPECTRUM.
REAL,    INTENT(OUT)   :: SL(:, :, :)    !! TOTAL SOURCE FUNCTION ARRAY
REAL,    INTENT(OUT)   :: SPOS(:, :, :)  !! POSITIVE SOURCE FUNCTION ARRAY
REAL,    INTENT(OUT)   :: FL(:, :, :)    !! DIAGONAL MATRIX OF FUNCTIONAL
                                         !! DERIVATIVE
REAL,    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY.
REAL,    INTENT(IN)    :: UDIR (:)       !! WIND DIRECTION.
REAL,    INTENT(IN)    :: Z0   (:)       !! ROUGHNESS LENGTH.
REAL,    INTENT(IN)    :: ROAIRN(:)      !! AIR DENSITY
REAL,    INTENT(IN)    :: WSTAR(:)       !! FREE CONVECTION VELOCITY SCLALE 
INTEGER, INTENT(IN)    :: INDEP(:)       !! DEPTH TABLE INDEX.
LOGICAL, INTENT(OUT)   :: LLWS(:, :, :)  !! TRUE WHERE SINPUT IS POSITIVE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------   

INTEGER :: IJ,K,M,IND,IGST

INTEGER , PARAMETER :: NGST=2

!     TEST471 (modified SWELLF):
REAL, PARAMETER :: SWELLF = 0.66 ! controls the turbulent swell dissipation
REAL, PARAMETER :: SWELLF2 = -0.018
REAL, PARAMETER :: SWELLF3 = 0.022
REAL, PARAMETER :: SWELLF4 = 1.5E05
REAL, PARAMETER :: SWELLF5 = 1.2  ! controls the viscous swell dissipation
REAL, PARAMETER :: SWELLF6 = 1.0
REAL, PARAMETER :: SWELLF7 = 3.6E05
REAL, PARAMETER :: Z0RAT = 0.04
REAL, PARAMETER :: Z0TUBMAX = 0.0005

REAL, PARAMETER :: ABMIN = 0.3
REAL, PARAMETER :: ABMAX = 8.0 


REAL :: ROG
REAL :: AVG_GST, ABS_TAUWSHELTER 
REAL :: CONST1
REAL :: ZLOG, ZLOG2X
REAL :: XI,X,DELI1,DELI2
REAL :: FU,FUD,NU_AIR,SMOOTH, HFTSWELLF6,Z0TUB
REAL :: FAC_NU_AIR,FACM1_NU_AIR
REAL :: DELABM1
REAL :: TAUPX,TAUPY
REAL :: DSTAB2
REAL, DIMENSION(SIZE(F,3)) :: CONST, SIG, SIGM1, SIG2, COEF, COEF5, DFIM_SIG2
REAL, DIMENSION(SIZE(F,1)) :: CONSTF, CONST11, CONST22
REAL, DIMENSION(SIZE(F,1)) :: Z0VIS, Z0NOZ, FWW
REAL, DIMENSION(SIZE(F,1)) :: PVISC, PTURB
REAL, DIMENSION(SIZE(F,1)) :: ZCN
REAL, DIMENSION(SIZE(F,1)) :: SIG_N, UORBT, AORB, TEMP, RE, RE_C, ZORB
REAL, DIMENSION(SIZE(F,1)) :: CNSN
REAL, DIMENSION(SIZE(F,1)) :: FLP_AVG, SLP_AVG
REAL, DIMENSION(SIZE(F,1)) :: GZ0, ROGOROAIR, ROAIRN_PVISC
REAL, DIMENSION(SIZE(F,1),NGST) :: XSTRESS, YSTRESS, FLP, SLP
REAL, DIMENSION(SIZE(F,1),NGST) :: USG2, TAUX, TAUY, USTP, USDIRP, UCN
REAL, DIMENSION(SIZE(F,1),NGST) :: UCNZALPD
REAL, DIMENSION(SIZE(F,1),SIZE(F,3)) :: CM, XK, SH
REAL, DIMENSION(SIZE(F,1)) :: DSTAB1, TEMP1, TEMP2
REAL, DIMENSION(SIZE(F,1),SIZE(F,2),NGST) :: COSLP, UFAC, DSTAB

! ----------------------------------------------------------------------

      ROG = ROWATER*G
      AVG_GST = 1.0/NGST
      CONST1  = BETAMAX/XKAPPA**2 /ROWATER
      NU_AIR = RNUAIR
      FAC_NU_AIR= RNUAIRM
      FACM1_NU_AIR=4.0/NU_AIR

      FU=ABS(SWELLF3)
      FUD=SWELLF2
      DELABM1= REAL(IAB)/(ABMAX-ABMIN)

      ABS_TAUWSHELTER=ABS(TAUWSHELTER)

!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      CALL WSIGSTAR (USTAR, Z0, WSTAR, SIG_N)

! ----------------------------------------------------------------------

      DO IJ = 1,SIZE(F,1)
        UORBT(IJ) = TINY(1.0) 
        AORB(IJ) = TINY(1.0) 
      ENDDO

      DO M = 1,SIZE(F,3)
        SIG(M) = ZPI*FR(M)
        SIGM1(M) = 1.0/SIG(M)
        SIG2(M) = SIG(M)**2
        DFIM_SIG2(M)=DFIM(M)*SIG2(M)
        CONST(M)=SIG(M)*CONST1
        COEF(M) =-SWELLF*16.*SIG2(M)/(G*ROWATER)
        COEF5(M)=-SWELLF5*2.*SQRT(2.*NU_AIR*SIG(M))/ROWATER

        K=1
        DO IJ = 1,SIZE(F,1)
          TEMP(IJ) = F(IJ,K,M)
        ENDDO
        DO K = 2,SIZE(F,2)
          DO IJ = 1,SIZE(F,1)
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ = 1,SIZE(F,1)
          UORBT(IJ) = UORBT(IJ)+DFIM_SIG2(M)*TEMP(IJ)
          AORB(IJ) = AORB(IJ)+DFIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

      DO IJ = 1,SIZE(F,1)
        UORBT(IJ) = 2.0*SQRT(UORBT(IJ))  ! this is the significant orbital amplitude
        AORB(IJ)  = 2.0*SQRT(AORB(IJ))   ! this 1/2 Hs
        RE(IJ)    = FACM1_NU_AIR*UORBT(IJ)*AORB(IJ) ! this is the Reynolds number 
        Z0VIS(IJ) = FAC_NU_AIR/MAX(USTAR(IJ),0.0001)
        Z0TUB = Z0RAT*MIN(Z0TUBMAX,Z0(IJ))
        Z0NOZ(IJ) = MAX(Z0VIS(IJ),Z0TUB)
        ZORB(IJ)  = AORB(IJ)/Z0NOZ(IJ)

! conpute fww
        XI=(LOG10(MAX(ZORB(IJ),3.0))-ABMIN)*DELABM1
        IND  = MIN (IAB-1, INT(XI))
        DELI1= MIN (1.0 ,XI-FLOAT(IND))
        DELI2= 1.0 - DELI1
        FWW(IJ)= SWELLFT(IND)*DELI2+SWELLFT(IND+1)*DELI1
        TEMP2(IJ) = FWW(IJ)*UORBT(IJ)
      ENDDO

! Define the critical Reynolds number
      IF( SWELLF6 .EQ. 1.0) THEN
        DO IJ = 1,SIZE(F,1)
          RE_C(IJ) = SWELLF4
        ENDDO
      ELSE
        HFTSWELLF6=1.0-SWELLF6
        DO IJ = 1,SIZE(F,1)
          RE_C(IJ) = SWELLF4*(2.0/AORB(IJ))**HFTSWELLF6
        ENDDO
      ENDIF

! Swell damping weight between viscous and turbulent boundary layer
      IF (SWELLF7.GT.0.0) THEN
        DO IJ = 1,SIZE(F,1)
          SMOOTH=0.5*TANH((RE(IJ)-RE_C(IJ))/SWELLF7)
          PTURB(IJ)=0.5+SMOOTH
          PVISC(IJ)=0.5-SMOOTH
        ENDDO
      ELSE
        DO IJ = 1,SIZE(F,1)
          IF (RE(IJ).LE.RE_C(IJ)) THEN
            PTURB(IJ)=0.0
            PVISC(IJ)=0.5
          ELSE
            PTURB(IJ)=0.5
            PVISC(IJ)=0.0
          ENDIF
        ENDDO
      ENDIF


! Initialisation

      DO IGST=1,NGST
        DO IJ = 1,SIZE(F,1)
          XSTRESS(IJ,IGST)=0.0
          YSTRESS(IJ,IGST)=0.0
        ENDDO
      ENDDO
      DO IJ = 1,SIZE(F,1)
        USG2(IJ,1)=(USTAR(IJ)*(1.0+SIG_N(IJ)))**2
        USG2(IJ,2)=(USTAR(IJ)*(1.0-SIG_N(IJ)))**2
      ENDDO

      DO IGST=1,NGST
        DO IJ = 1,SIZE(F,1)
          TAUX(IJ,IGST)=USG2(IJ,IGST)*SIN(UDIR(IJ))
          TAUY(IJ,IGST)=USG2(IJ,IGST)*COS(UDIR(IJ))
        ENDDO
      ENDDO

      DO IJ = 1,SIZE(F,1)
        GZ0(IJ) = G*Z0(IJ)
      ENDDO
      DO IJ = 1,SIZE(F,1)
        ROGOROAIR(IJ) = ROG/MAX(ROAIRN(IJ),1.0)
      ENDDO

      DO IJ = 1,SIZE(F,1)
        ROAIRN_PVISC(IJ) = ROAIRN(IJ)*PVISC(IJ)
      ENDDO


!     INVERSE OF PHASE VELOCITIES AND WAVE NUMBER.
      IF (SHALLOW_RUN) THEN
        DO M = 1,SIZE(F,3)
          DO IJ = 1,SIZE(F,1)
            XK(IJ,M) = TFAK(INDEP(IJ),M)
          ENDDO
        ENDDO
        DO M = 1,SIZE(F,3)
          DO IJ = 1,SIZE(F,1)
            CM(IJ,M) = XK(IJ,M)*SIGM1(M)
!!            SH(IJ,M) = SIG2(M)/(G*XK(IJ,M)) ! tanh(kh)
            SH(IJ,M) = 1.0
          ENDDO
        ENDDO
      ELSE
        DO M = 1,SIZE(F,3)
          DO IJ = 1,SIZE(F,1)
            CM(IJ,M) = SIG(M)/G
            XK(IJ,M) = SIG2(M)/G
            SH(IJ,M) = 1.0
          ENDDO
        ENDDO
      ENDIF

!*    2. MAIN LOOP OVER FREQUENCIES.
!        ---------------------------

      DO M = 1,SIZE(F,3)

        DO IGST=1,NGST
          DO IJ = 1,SIZE(F,1)
            TAUPX=TAUX(IJ,IGST)-ABS_TAUWSHELTER*XSTRESS(IJ,IGST)
            TAUPY=TAUY(IJ,IGST)-ABS_TAUWSHELTER*YSTRESS(IJ,IGST)
            USTP(IJ,IGST)=(TAUPX**2+TAUPY**2)**0.25
            USDIRP(IJ,IGST)=ATAN2(TAUPX,TAUPY)
          ENDDO
        ENDDO

        DO IJ = 1,SIZE(F,1)
          CONSTF(IJ) = ROGOROAIR(IJ)*CM(IJ,M)*DFIM(M)
        ENDDO


        DO IGST=1,NGST
          DO K = 1,SIZE(F,2)
            DO IJ = 1,SIZE(F,1)
              COSLP(IJ,K,IGST) = COS(TH(K)-USDIRP(IJ,IGST))
            ENDDO
          ENDDO
        ENDDO

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IGST=1,NGST
          DO IJ = 1,SIZE(F,1)
            UCN(IJ,IGST) = USTP(IJ,IGST)*CM(IJ,M)
            UCNZALPD(IJ,IGST) = XKAPPA/(UCN(IJ,IGST) + ZALP)
          ENDDO
        ENDDO
        DO IJ = 1,SIZE(F,1)
          ZCN(IJ) = LOG(XK(IJ,M)*Z0(IJ))
          CNSN(IJ) = CONST(M)*SH(IJ,M)*ROAIRN(IJ)
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K = 1,SIZE(F,2)
          DO IJ = 1,SIZE(F,1)
            LLWS(IJ,K,M)=.FALSE.
          ENDDO
        ENDDO

        DO IGST=1,NGST
          DO K = 1,SIZE(F,2)
            DO IJ = 1,SIZE(F,1)
              IF (COSLP(IJ,K,IGST).GT.0.01) THEN
                X    = COSLP(IJ,K,IGST)*UCN(IJ,IGST)
                ZLOG = ZCN(IJ) + UCNZALPD(IJ,IGST)/COSLP(IJ,K,IGST)
                IF (ZLOG.LT.0.0) THEN
                  ZLOG2X=ZLOG*ZLOG*X
                  UFAC(IJ,K,IGST) = EXP(ZLOG)*ZLOG2X*ZLOG2X
                  LLWS(IJ,K,M)=.TRUE.
                ELSE
                  UFAC(IJ,K,IGST)=0.0
                ENDIF
              ELSE
                UFAC(IJ,K,IGST)=0.0
              ENDIF
            ENDDO
          ENDDO
        ENDDO

!       SWELL DAMPING:

        DO IJ = 1,SIZE(F,1)
          DSTAB1(IJ) = COEF5(M)*ROAIRN_PVISC(IJ)*XK(IJ,M)
          TEMP1(IJ) = COEF(M)*ROAIRN(IJ)
        ENDDO

        DO IGST=1,NGST
          DO K = 1,SIZE(F,2)
            DO IJ = 1,SIZE(F,1)
              DSTAB2 = TEMP1(IJ)*(TEMP2(IJ)+(FU+FUD*COSLP(IJ,K,IGST))*USTP(IJ,IGST))
              DSTAB(IJ,K,IGST) = DSTAB1(IJ)+PTURB(IJ)*DSTAB2
            ENDDO
          ENDDO
        ENDDO


!*    2.2 UPDATE THE SHELTERING STRESS,
!         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ---------------------------------------------------------

        DO K = 1,SIZE(F,2)
          DO IJ = 1,SIZE(F,1)
            CONST11(IJ)=CONSTF(IJ)*SINTH(K)
            CONST22(IJ)=CONSTF(IJ)*COSTH(K)
          ENDDO

          DO IGST=1,NGST
            DO IJ = 1,SIZE(F,1)
              ! SLP: only the positive contributions
              SLP(IJ,IGST) = CNSN(IJ)*UFAC(IJ,K,IGST)
              FLP(IJ,IGST) = SLP(IJ,IGST)+DSTAB(IJ,K,IGST)
            ENDDO
          ENDDO

          DO IGST=1,NGST
            DO IJ = 1,SIZE(F,1)
              SLP(IJ,IGST) = SLP(IJ,IGST)*F(IJ,K,M)
              XSTRESS(IJ,IGST)=XSTRESS(IJ,IGST)+SLP(IJ,IGST)*CONST11(IJ)
              YSTRESS(IJ,IGST)=YSTRESS(IJ,IGST)+SLP(IJ,IGST)*CONST22(IJ)
            ENDDO
          ENDDO

          IGST=1
            DO IJ = 1,SIZE(F,1)
              SLP_AVG(IJ) = SLP(IJ,IGST)
              FLP_AVG(IJ) = FLP(IJ,IGST)
            ENDDO
          DO IGST=2,NGST
            DO IJ = 1,SIZE(F,1)
              SLP_AVG(IJ) = SLP_AVG(IJ)+SLP(IJ,IGST)
              FLP_AVG(IJ) = FLP_AVG(IJ)+FLP(IJ,IGST)
            ENDDO
          ENDDO

          DO IJ = 1,SIZE(F,1)
            SPOS(IJ,K,M) = AVG_GST*SLP_AVG(IJ)
          ENDDO

          DO IJ = 1,SIZE(F,1)
            FL(IJ,K,M) = AVG_GST*FLP_AVG(IJ)
            SL(IJ,K,M) = FL(IJ,K,M)*F(IJ,K,M)
          ENDDO

        ENDDO

      ENDDO ! END LOOP OVER FREQUENCIES


END SUBROUTINE SINPUT_ARD

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
!     P. JANSSEN  ECMWF  JUNE 2005       IMPROVED SCALING IN SHALLOW WATER     !
!     J. BIDLOT   ECMWF  AUGUST 2006     KEEP THE OLD FORMULATION              !
!                                        UNDER A SWITCH (ISNONLIN = 0 for OLD  !
!                                                                 = 1 for NEW  !
!                                        BE AWARE THAT THE OLD FORMULATION     !
!                                        REQUIRES THE MEAN WAVE NUMBER AKMEAN. !
!     J. BIDLOT   ECMWF  JANUARY 2012    ADD EXTENSION TO LOW FREQUENCIES      !
!                                        OPTIMISATION FOR IBM.                 !
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

REAL, INTENT(IN)    :: F (:,:,:)    !! SPECTRA.
REAL, INTENT(INOUT) :: SL(:,:,:)    !! TOTAL SOURCE FUNCTION ARRAY.
REAL, INTENT(INOUT) :: FL(:,:,:)    !! DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
REAL, INTENT(IN)    :: DEPTH (:)    !! WATER DEPTH.
REAL, INTENT(IN)    :: AKMEAN (:)   !! MEAN WAVE NUMBER.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: MC, MP, MP1, MM, MM1, IC, IP, IP1, IM, IM1, KH, K, K1, K2, K11, K21
INTEGER :: IJ
INTEGER :: MFR1STFR, MFRLSTFR

REAL    :: FTAIL, FKLAMP, FKLAMP1, GW1, GW2, GW3, GW4
REAL    :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAPA2, FKLAPB2, FKLAP12, FKLAP22
REAL    :: FKLAMM, FKLAMM1, GW5, GW6, GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2
REAL    :: FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
REAL    :: SAP, SAM, FIJ, FAD1, FAD2, FCEN

REAL, DIMENSION(SIZE(F,1)) :: FTEMP
REAL, DIMENSION(SIZE(F,1)) :: ENHFR
REAL, DIMENSION(SIZE(F,1)) :: AD
REAL, DIMENSION(SIZE(F,1)) :: DELAD
REAL, DIMENSION(SIZE(F,1)) :: DELAP
REAL, DIMENSION(SIZE(F,1)) :: DELAM

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. SHALLOW WATER INITIALISATION.                                         !
!        -----------------------------                                         !

IF (SHALLOW_RUN) THEN
   IF (ISNONLIN.EQ.0) THEN
      DO IJ = 1,SIZE(F,1)
         ENHFR(IJ) = MAX(0.75*DEPTH(IJ)*AKMEAN(IJ) , .5)
         ENHFR(IJ) = 1. + (5.5/ENHFR(IJ)) * (1.-.833*ENHFR(IJ))               &
&                    * EXP(-1.25*ENHFR(IJ))
      END DO
      DO MC=1,MLSTHG
         ENH(:,MC) = ENHFR(:)
      END DO
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. FREQUENCY LOOP.                                                       !
!        ---------------                                                       !

MFR1STFR = -MFRSTLW+1
MFRLSTFR = ML-KFRH+MFR1STFR

FRE: DO MC = 1,MLSTHG
   MP  = IKP (MC)
   MP1 = IKP1(MC)
   MM  = IKM (MC)
   MM1 = IKM1(MC)
   IC  = INLCOEF(1,MC)
   IP  = INLCOEF(2,MC)
   IP1 = INLCOEF(3,MC)
   IM  = INLCOEF(4,MC)
   IM1 = INLCOEF(5,MC)

   FTAIL  = RNLCOEF(1,MC)

   FKLAMP  = FKLAP(MC)
   FKLAMP1 = FKLAP1(MC)
   GW1 = RNLCOEF(2,MC)
   GW2 = RNLCOEF(3,MC)
   GW3 = RNLCOEF(4,MC)
   GW4 = RNLCOEF(5,MC)
   FKLAMPA = RNLCOEF(6,MC)
   FKLAMPB = RNLCOEF(7,MC)
   FKLAMP2 = RNLCOEF(8,MC)
   FKLAMP1 = RNLCOEF(9,MC)
   FKLAPA2 = RNLCOEF(10,MC)
   FKLAPB2 = RNLCOEF(11,MC)
   FKLAP12 = RNLCOEF(12,MC)
   FKLAP22 = RNLCOEF(13,MC)

   FKLAMM  = FKLAM(MC)
   FKLAMM1 = FKLAM1(MC)
   GW5 = RNLCOEF(14,MC)
   GW6 = RNLCOEF(15,MC)
   GW7 = RNLCOEF(16,MC)
   GW8 = RNLCOEF(17,MC)
   FKLAMMA = RNLCOEF(18,MC)
   FKLAMMB = RNLCOEF(19,MC)
   FKLAMM2 = RNLCOEF(20,MC)
   FKLAMM1 = RNLCOEF(21,MC)
   FKLAMA2 = RNLCOEF(22,MC)
   FKLAMB2 = RNLCOEF(23,MC)
   FKLAM12 = RNLCOEF(24,MC)
   FKLAM22 = RNLCOEF(25,MC)

   IF (SHALLOW_RUN) THEN
      FTEMP(:) = AF11(MC)*ENH(:,MC)
   ELSE
      FTEMP(:) = AF11(MC)
   ENDIF

   IF (MC.GT.MFR1STFR .AND. MC.LT.MFRLSTFR ) THEN
!       MC is within the fully resolved spectral domain

      DO KH=1,2
         DO K=1,KL
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)

!*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
!*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!             ----------------------------------------------

            DO IJ = 1,SIZE(F,1)
               SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )                     &
&                  + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
               SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )                     &
&                  + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
!!!! not needed ftail always=1.                FIJ = F(IJ,K  ,IC )*FTAIL
               FIJ = F(IJ,K  ,IC )
               FAD1 = FIJ*(SAP+SAM)
               FAD2 = FAD1-2.*SAP*SAM
               FAD1 = FAD1+FAD2
               FCEN = FTEMP(IJ)*FIJ
               AD(IJ) = FAD2*FCEN
               DELAD(IJ) = FAD1*FTEMP(IJ)
               DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FCEN
               DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FCEN
            ENDDO

            SL(:,K  ,MC ) = SL(:,K  ,MC ) - 2.*AD(:)
            FL(:,K  ,MC ) = FL(:,K  ,MC ) - 2.*DELAD(:)
            SL(:,K2 ,MM ) = SL(:,K2 ,MM ) + AD(:)*FKLAMM1
            FL(:,K2 ,MM ) = FL(:,K2 ,MM ) + DELAM(:)*FKLAM12
            SL(:,K21,MM ) = SL(:,K21,MM ) + AD(:)*FKLAMM2
            FL(:,K21,MM ) = FL(:,K21,MM ) + DELAM(:)*FKLAM22
            SL(:,K2 ,MM1) = SL(:,K2 ,MM1) + AD(:)*FKLAMMA
            FL(:,K2 ,MM1) = FL(:,K2 ,MM1) + DELAM(:)*FKLAMA2
            SL(:,K21,MM1) = SL(:,K21,MM1) + AD(:)*FKLAMMB
            FL(:,K21,MM1) = FL(:,K21,MM1) + DELAM(:)*FKLAMB2
            SL(:,K1 ,MP ) = SL(:,K1 ,MP ) + AD(:)*FKLAMP1
            FL(:,K1 ,MP ) = FL(:,K1 ,MP ) + DELAP(:)*FKLAP12
            SL(:,K11,MP ) = SL(:,K11,MP ) + AD(:)*FKLAMP2
            FL(:,K11,MP ) = FL(:,K11,MP ) + DELAP(:)*FKLAP22
            SL(:,K1 ,MP1) = SL(:,K1 ,MP1) + AD(:)*FKLAMPA
            FL(:,K1 ,MP1) = FL(:,K1 ,MP1) + DELAP(:)*FKLAPA2
            SL(:,K11,MP1) = SL(:,K11,MP1) + AD(:)*FKLAMPB
            FL(:,K11,MP1) = FL(:,K11,MP1) + DELAP(:)*FKLAPB2
         ENDDO
      ENDDO

   ELSEIF (MC.GE.MFRLSTFR ) THEN
      DO KH=1,2
         DO K=1,KL
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)

            DO IJ = 1,SIZE(F,1)
               SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )                     &
&                  + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
               SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )                     &
&                  + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
               FIJ = F(IJ,K  ,IC )*FTAIL
               FAD1 = FIJ*(SAP+SAM)
               FAD2 = FAD1-2.*SAP*SAM
               FAD1 = FAD1+FAD2
               FCEN = FTEMP(IJ)*FIJ
               AD(IJ) = FAD2*FCEN
               DELAD(IJ) = FAD1*FTEMP(IJ)
               DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FCEN
               DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FCEN
            ENDDO

            SL(:,K2 ,MM ) = SL(:,K2 ,MM ) + AD(:)*FKLAMM1
            FL(:,K2 ,MM ) = FL(:,K2 ,MM ) + DELAM(:)*FKLAM12
            SL(:,K21,MM ) = SL(:,K21,MM ) + AD(:)*FKLAMM2
            FL(:,K21,MM ) = FL(:,K21,MM ) + DELAM(:)*FKLAM22

            IF (MM1.LE.ML) THEN
               SL(:,K2 ,MM1) = SL(:,K2 ,MM1) + AD(:)*FKLAMMA
               FL(:,K2 ,MM1) = FL(:,K2 ,MM1) + DELAM(:)*FKLAMA2
               SL(:,K21,MM1) = SL(:,K21,MM1) + AD(:)*FKLAMMB
               FL(:,K21,MM1) = FL(:,K21,MM1) + DELAM(:)*FKLAMB2

               IF (MC .LE.ML) THEN
                  SL(:,K  ,MC ) = SL(:,K  ,MC ) - 2.*AD(:)
                  FL(:,K  ,MC ) = FL(:,K  ,MC ) - 2.*DELAD(:)

                  IF (MP .LE.ML) THEN
                     SL(:,K1 ,MP ) = SL(:,K1 ,MP ) + AD(:)*FKLAMP1
                     FL(:,K1 ,MP ) = FL(:,K1 ,MP ) + DELAP(:)*FKLAP12
                     SL(:,K11,MP ) = SL(:,K11,MP ) + AD(:)*FKLAMP2
                     FL(:,K11,MP ) = FL(:,K11,MP ) + DELAP(:)*FKLAP22

                     IF (MP1.LE.ML) THEN
                        SL(:,K1 ,MP1) = SL(:,K1 ,MP1) + AD(:)*FKLAMPA
                        FL(:,K1 ,MP1) = FL(:,K1 ,MP1) + DELAP(:)*FKLAPA2
                        SL(:,K11,MP1) = SL(:,K11,MP1) + AD(:)*FKLAMPB
                        FL(:,K11,MP1) = FL(:,K11,MP1) + DELAP(:)*FKLAPB2
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

   ELSE

      DO KH=1,2
         DO K=1,KL
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)

            DO IJ = 1,SIZE(F,1)
               SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )     &
&                  + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
               SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )     &
&                  + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
               FIJ = F(IJ,K  ,IC )*FTAIL
               FAD1 = FIJ*(SAP+SAM)
               FAD2 = FAD1-2.*SAP*SAM
               FAD1 = FAD1+FAD2
               FCEN = FTEMP(IJ)*FIJ
               AD(IJ) = FAD2*FCEN
               DELAD(IJ) = FAD1*FTEMP(IJ)
               DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FCEN
               DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FCEN
            ENDDO

            IF (MM1.GE.1) THEN
               SL(:,K2 ,MM1) = SL(:,K2 ,MM1) + AD(:)*FKLAMMA
               FL(:,K2 ,MM1) = FL(:,K2 ,MM1) + DELAM(:)*FKLAMA2
               SL(:,K21,MM1) = SL(:,K21,MM1) + AD(:)*FKLAMMB
               FL(:,K21,MM1) = FL(:,K21,MM1) + DELAM(:)*FKLAMB2
            ENDIF

            SL(:,K  ,MC ) = SL(:,K  ,MC ) - 2.*AD(:)
            FL(:,K  ,MC ) = FL(:,K  ,MC ) - 2.*DELAD(:)
            SL(:,K1 ,MP ) = SL(:,K1 ,MP ) + AD(:)*FKLAMP1
            FL(:,K1 ,MP ) = FL(:,K1 ,MP ) + DELAP(:)*FKLAP12
            SL(:,K11,MP ) = SL(:,K11,MP ) + AD(:)*FKLAMP2
            FL(:,K11,MP ) = FL(:,K11,MP ) + DELAP(:)*FKLAP22
            SL(:,K1 ,MP1) = SL(:,K1 ,MP1) + AD(:)*FKLAMPA
            FL(:,K1 ,MP1) = FL(:,K1 ,MP1) + DELAP(:)*FKLAPA2
            SL(:,K11,MP1) = SL(:,K11,MP1) + AD(:)*FKLAMPB
            FL(:,K11,MP1) = FL(:,K11,MP1) + DELAP(:)*FKLAPB2
         ENDDO
      ENDDO

   ENDIF
END DO FRE                  !! BRANCH BACK FOR NEXT FREQUENCY.

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

SUBROUTINE STRESSO (F, SL, USTAR, UDIR, Z0, MIJ, TAUW, PHIAW, INDEP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    STRESSO - COMPUTATION OF WAVE STRESS.                                     !
!                                                                              !
!     H. GUNTHER      GKSS/ECMWF  NOVEMBER  1989 CODE MOVED FROM SINPUT.       !
!     P.A.E.M. JANSSEN      KNMI  AUGUST    1990                               !
!     J. BIDLOT             ECMWF FEBRUARY  1996-97                            !
!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90                            !
!     J. BIDLOT             ECMWF           2007  ADD MIJ                      !
!     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS         !
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

REAL,    INTENT(IN)    :: F(:,:,:)       !! WAVE SPECTRUM.
REAL,    INTENT(IN)    :: SL(:,:,:)      !! INPUT SOURCE FUNCTION.
REAL,    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY.
REAL,    INTENT(IN)    :: UDIR(:)        !! WIND DIRECTION.
REAL,    INTENT(IN)    :: Z0(:)          !! ROUGHNESS LENGTH.
INTEGER, INTENT(IN)    :: MIJ(:)         !! LAST FREQUENCY INDEX OF THE
                                         !! PROGNOSTIC RANGE.
REAL,    INTENT(OUT)   :: TAUW(:)        !! WAVE STRESS.
REAL,    INTENT(OUT)   :: PHIAW(:)       !! ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
                                         !! OVER THE FULL FREQUENCY RANGE.

INTEGER, INTENT(IN)    :: INDEP(:)       !! DEPTH TABLE INDEX.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: K, M, KL, ML, IJ

REAL :: GM1, COSW, CONST, SINPLUS

REAL, DIMENSION(SIZE(F,1)) :: UST2
REAL, DIMENSION(SIZE(F,1)) :: TAUHF, TEMP1, CONST1, XSTRESS, YSTRESS
REAL, DIMENSION(SIZE(F,1)) :: PHIHF, TEMP2, CONST2
REAL, DIMENSION(SIZE(F,1)) :: SUMT, SUMX, SUMY
REAL, DIMENSION(SIZE(F,1)) :: TAU1, PHI1, XLEVTAIL
REAL, DIMENSION(SIZE(F,3)) :: SIG, SIGM1
REAL, DIMENSION(SIZE(F,1),SIZE(F,3)) :: CM, RHOWGDFTH
REAL, DIMENSION(SIZE(F,1)) :: CMRHOWGDFTH

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. PRECOMPUTE FREQUENCY SCALING.                                         !
!        -----------------------------                                         !

KL = SIZE(F,2)
ML = SIZE(F,3)

XLEVTAIL(:) = 0.0
UST2(:) = ustar(:)**2

GM1 = 1.0/G
CONST = DELTH*(ZPI)**4*GM1

DO IJ = 1,SIZE(F,1)
   CONST1(IJ)  = CONST*FR5(MIJ(IJ))*GM1
   CONST2(IJ)  = ROAIR*CONST*FR5(MIJ(IJ))
ENDDO

!     INVERSE OF PHASE VELOCITIES

IF (SHALLOW_RUN) THEN
    DO M = 1,SIZE(F,3)
       SIGM1(M) = 1.0/(ZPI*FR(M))
       DO IJ = 1,SIZE(F,1)
          CM(IJ,M) = TFAK(INDEP(IJ),M)*SIGM1(M)
       ENDDO
    ENDDO
ELSE
   DO M = 1,SIZE(F,3)
      SIG(M) = ZPI*FR(M)
      DO IJ = 1,SIZE(F,1)
         CM(IJ,M) = SIG(M)*GM1
      ENDDO
    ENDDO
ENDIF

DO IJ = 1,SIZE(USTAR)
   DO M=1,MIJ(IJ)
      RHOWGDFTH(IJ,M) = RHOWG_DFIM(M)
   ENDDO
   IF(MIJ(IJ).NE.ML) RHOWGDFTH(IJ,MIJ(IJ))=0.5*RHOWGDFTH(IJ,MIJ(IJ))
   DO M=MIJ(IJ)+1,ML
      RHOWGDFTH(IJ,M) = 0.0
   ENDDO
ENDDO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE WAVE STRESS OF BLOCK.                                         !
!        -----------------------------                                         !
!                                                                              !
!     2.1 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.       !
!         --------------------------------------------------------------       !

PHIAW(:)   = 0.0
XSTRESS(:) = 0.0
YSTRESS(:) = 0.0

DO M=1,MAXVAL(MIJ(:)) !! THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
   K=1
   DO IJ = 1,SIZE(F,1)
      SINPLUS = MAX (SL(IJ,K,M),0.)
      SUMT(IJ) = SINPLUS
      SUMX(IJ) = SINPLUS*SINTH(K)
      SUMY(IJ) = SINPLUS*COSTH(K)
   ENDDO
   DO K = 2,SIZE(F,2)
      DO IJ = 1,SIZE(F,1)
         SINPLUS = MAX (SL(IJ,K,M),0.)
         SUMT(IJ) = SUMT(IJ) + SINPLUS
         SUMX(IJ) = SUMX(IJ) + SINPLUS*SINTH(K)
         SUMY(IJ) = SUMY(IJ) + SINPLUS*COSTH(K)
      ENDDO
   ENDDO
   DO IJ = 1,SIZE(F,1)
      PHIAW(IJ)   =  PHIAW(IJ) + SUMT(IJ)*RHOWGDFTH(IJ,M)
      CMRHOWGDFTH(IJ) = CM(IJ,M)*RHOWGDFTH(IJ,M)
      XSTRESS(IJ) = XSTRESS(IJ) + SUMX(IJ)*CMRHOWGDFTH(IJ)
      YSTRESS(IJ) = YSTRESS(IJ) + SUMY(IJ)*CMRHOWGDFTH(IJ)
   ENDDO
ENDDO

DO IJ = 1,SIZE(F,1)
   XSTRESS(IJ) = XSTRESS(IJ)/MAX(ROAIR,1.)
   YSTRESS(IJ) = YSTRESS(IJ)/MAX(ROAIR,1.)
ENDDO

!     2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!     ----------------------------------------------------

K=1
DO IJ = 1,SIZE(F,1)
   COSW     = MAX(COS(TH(K)-UDIR(IJ)),0.)
   TEMP2(IJ) = F(IJ,K,MIJ(IJ))*COSW**2
   TEMP1(IJ) = TEMP2(IJ)*COSW
ENDDO

DO K=2,KL
   DO IJ = 1,SIZE(F,1)
      COSW     = MAX(COS(TH(K)-UDIR(IJ)),0.)
      TEMP1(IJ) = TEMP1(IJ)+F(IJ,K,MIJ(IJ))*COSW**3
      TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,MIJ(IJ))*COSW**2
   ENDDO
ENDDO

CALL TAU_PHI_HF(MIJ, USTAR, Z0, XLEVTAIL, TAU1, PHI1)

DO IJ = 1,SIZE(F,1)
   TAUHF(IJ) = CONST1(IJ)*TEMP1(IJ)*TAU1(IJ)
   PHIHF(IJ) = CONST2(IJ)*TEMP2(IJ)*PHI1(IJ)
ENDDO

DO IJ = 1,SIZE(F,1)
   PHIAW(IJ)   = PHIAW(IJ)   + PHIHF(IJ)
   XSTRESS(IJ) = XSTRESS(IJ) + TAUHF(IJ)*SIN(UDIR(IJ))
   YSTRESS(IJ) = YSTRESS(IJ) + TAUHF(IJ)*COS(UDIR(IJ))
   TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)

   TAUW(IJ) = MIN(TAUW(IJ),UST2(IJ)-EPS1)
   TAUW(IJ) = MAX(TAUW(IJ),0.)
END DO

END SUBROUTINE STRESSO

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TABU_SWELLFT

!**** *TABU_SWELLFT* - FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS

!     FABRICE ARDHUIN  IFREMER  2013

!*    PURPOSE.
!     --------
!     TO ESTIMATE FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS

!     METHOD.
!     -------
!       TABULATION ON KELVIN FUNCTIONS.

!     EXTERNALS.
!     -----------

!     KERKEI  (zeroth order Kelvin function Ker and Kei)

! ----------------------------------------------------------------------

      INTEGER, PARAMETER :: NITER=100
      REAL,    PARAMETER :: ABMIN=0.3
      REAL,    PARAMETER :: ABMAX=8.0, KAPPA=0.40
!     VARIABLE.   TYPE.     PURPOSE.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
! ----------------------------------------------------------------------
      INTEGER :: I,ITER
      REAL :: DELAB
      REAL :: KER, KEI
      REAL :: ABR,ABRLOG,L10,FACT,FSUBW,FSUBWMEMO,DZETA0,DZETA0MEMO

! ----------------------------------------------------------------------

      DZETA0 = 0.0
!
      DELAB   = (ABMAX-ABMIN)/REAL(IAB)
      L10=ALOG(10.0)
      DO I=1,IAB
         ABRLOG=ABMIN+REAL(I)*DELAB
         ABR=EXP(ABRLOG*L10)
         FACT=1/ABR/(21.2*KAPPA)
         FSUBW=0.05
         DO ITER=1,NITER
            FSUBWMEMO=FSUBW
            DZETA0MEMO=DZETA0
            DZETA0=FACT*FSUBW**(-0.5)
            CALL KERKEI(2.0*SQRT(DZETA0),KER,KEI)
            FSUBW=0.08/(KER**2+KEI**2)
            FSUBW=0.5*(FSUBWMEMO+FSUBW)
            DZETA0=0.5*(DZETA0MEMO+DZETA0)
         ENDDO   
         SWELLFT(I)=FSUBW
      ENDDO

END SUBROUTINE TABU_SWELLFT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE INIT_SDISSP_ARD

! ----------------------------------------------------------------------

!**** *INIT_SDISSP_ARD* - INITIALISATION FOR SDISS_ARD

!     FABRICE ARDHUIN  IFREMER  2013

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *INIT_SDISSP_ARD

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1

! ----------------------------------------------------------------------


 INTEGER :: JD, K, M, I_INT, J_INT, M2, KK, KLD

 REAL :: TPIINV, TMP02
 REAL :: DELTH_TRUNC, DELTH_LOC

 REAL :: XLOGDFRTH
 REAL :: BRLAMBDA
 REAL, DIMENSION(0:KL/2) :: COSDTH
 REAL, DIMENSION(ML) :: SIG, C_, C_C, CGM1, DSIP, TRPZ_DSIP 

! ----------------------------------------------------------------------

      TPIINV = 1.0/ZPI

      KLD=KL/2

      XLOGDFRTH=LOG(CO)*DELTH

!     COMPUTE SATWEIGHTS

!     l(k,th)=1/(2*pi²)= the breaking crest density
      BRLAMBDA=BRKPBCOEF/(2.0*ZPI**2)

      TMP02 = SSDSC3*BRLAMBDA

      NSDSNTH  = MIN(NINT(ISDSDTH*RAD/(DELTH)),KLD-1)
      DELTH_TRUNC=(TH(1)+ISDSDTH*RAD)-(TH(1+NSDSNTH)-0.5*DELTH)
      DELTH_TRUNC=MAX(0.,MIN(DELTH_TRUNC,DELTH))

      IF(ALLOCATED(INDICESSAT)) DEALLOCATE(INDICESSAT)
      ALLOCATE(INDICESSAT(KL,NSDSNTH*2+1))
      IF(ALLOCATED(SATWEIGHTS)) DEALLOCATE(SATWEIGHTS)
      ALLOCATE(SATWEIGHTS(KL,NSDSNTH*2+1))

      DO K=1,KL
        DO I_INT=K-NSDSNTH, K+NSDSNTH
          J_INT=I_INT
          IF (I_INT.LT.1)  J_INT=I_INT+KL
          IF (I_INT.GT.KL) J_INT=I_INT-KL
          INDICESSAT(K,I_INT-(K-NSDSNTH)+1)=J_INT

          IF(I_INT.EQ.K-NSDSNTH .OR. I_INT.EQ.K+NSDSNTH) THEN
            DELTH_LOC=DELTH_TRUNC
          ELSE
            DELTH_LOC=DELTH
          ENDIF
          SATWEIGHTS(K,I_INT-(K-NSDSNTH)+1)=DELTH_LOC*COS(TH(K)-TH(J_INT))**ISB
        END DO
      END DO

!     COMPUTE CUMULW (only if needed)
      IF (SSDSC3.NE.0.0) THEN
        IF(ALLOCATED(CUMULW)) DEALLOCATE(CUMULW)
        ALLOCATE(CUMULW(NDEPTH,0:KLD,ML,ML))

!       NDIKCUMUL is the  integer difference in frequency bands
!       between the "large breakers" and short "wiped-out waves"
!!! wrong !!???        NDIKCUMUL = NINT(SSDSBRF1/(CO-1.))
        NDIKCUMUL = NINT(-LOG(SSDSBRF1)/LOG(CO))

        DO KK=0,KLD
          COSDTH(KK)=COS(KK*DELTH)
        ENDDO

        DO M=1,ML
          SIG(M) = ZPI*FR(M)
        ENDDO

        DO JD=1,NDEPTH

          IF (SHALLOW_RUN) THEN
            DO M=1,ML
              C_(M)=SIG(M)/TFAK(JD,M)
              CGM1(M)=1.0/TCGOND(JD,M)
            ENDDO
          ELSE
            DO M=1,ML
              C_(M)=G/SIG(M)  ! Valid in deep water only !!!!!!!!!!!!
              CGM1(M)=2.0/C_(M) ! deep water !
            ENDDO
          ENDIF

          DO M=1,ML
            C_C(M)=C_(M)*C_(M)
            DSIP(M)=TMP02*SIG(M)*XLOGDFRTH*CGM1(M) !  coef*dtheta*dk = coef*dtheta*dsigma/cg
          ENDDO

          DO M=NDIKCUMUL+1,ML

            IF(M-NDIKCUMUL.GE.3) THEN
              TRPZ_DSIP(1)=0.5*DSIP(1)
              DO M2=2,M-NDIKCUMUL-1
                TRPZ_DSIP(M2)=DSIP(M2)
              ENDDO
              TRPZ_DSIP(M-NDIKCUMUL)=0.5*DSIP(M-NDIKCUMUL)
            ELSE
              DO M2=1,M-NDIKCUMUL
                TRPZ_DSIP(M2)=DSIP(M2)
              ENDDO
            ENDIF

            DO M2=1,M-NDIKCUMUL
              DO KK=0,KLD
                CUMULW(JD,KK,M2,M)=SQRT(ABS(C_C(M)+C_C(M2)-2.0*C_(M)*C_(M2)*COSDTH(KK)))*TRPZ_DSIP(M2)
              ENDDO 
            ENDDO
          ENDDO

        ENDDO ! JD

      ENDIF

 END SUBROUTINE INIT_SDISSP_ARD

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FRCUTINDEX (FM, FMWS, USTAR, MIJ)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   FRCUTINDEX - RETURNS THE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!
!     METHOD.
!     -------
!                                                                              !
!     COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.

!!! be aware that if this is NOT used, for iphys=1, the cumulative dissipation has to be
!!! re-activated !!!
!                                                                              !
!                                                                              !
!     EXTERNALS.
!     ---------
!                                                                              !
!     REFERENCE.
!     ----------
!                                                                              !
! ----------------------------------------------------------------------

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.
!     --------------------
!                                                                              !
REAL,    INTENT(IN)    :: FM(:)          !! MEAN FREQUENCY
REAL,    INTENT(IN)    :: FMWS(:)        !! MEAN FREQUENCY OF WINDSEA
REAL,    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY
INTEGER, INTENT(OUT)   :: MIJ(:)         !! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !
!                                                                              !
INTEGER :: IJ

REAL, PARAMETER :: EPSUS = 1.0E-6
REAL, PARAMETER :: FRIC = 28.0
REAL, PARAMETER :: TAILFACTOR = 2.5
REAL, PARAMETER :: TAILFACTOR_PM = 3.0

REAL :: FPMH, FPPM, FM2, FPM, FPM4

! -----------------------------------------------------------------------------!
!                                                                              !
!*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
!*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
!*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
!     ------------------------------------------------------------

FPMH = TAILFACTOR/FR(1)
FPPM = TAILFACTOR_PM*G/(FRIC*ZPI*FR(1))

DO IJ = 1,SIZE(USTAR)
    FM2 = MAX(FMWS(IJ),FM(IJ))*FPMH
    FPM = FPPM/MAX(USTAR(IJ),EPSUS)
    FPM4 = MAX(FM2,FPM)
    MIJ(IJ) = NINT(LOG10(FPM4)*INV_LOG_CO)+1
    MIJ(IJ) = MIN(MAX(1,MIJ(IJ)),ML)
ENDDO

END SUBROUTINE FRCUTINDEX

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE IMPHFTAIL (MIJ, INDEP, FL3) 

! ----------------------------------------------------------------------

!**** *IMPHFTAIL* - IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM

!*    PURPOSE.
!     --------

!     IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM ABOVE FREQUENCY INDEX MIJ

!     METHOD.
!     -------

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)    :: MIJ(:)         !! LAST FREQUENCY INDEX OF THE
                                         !! PROGNOSTIC RANGE.
INTEGER, INTENT(IN)    :: INDEP(:)       !! DEPTH TABLE INDEX.
REAL,    INTENT(INOUT) :: FL3 (:, :, :)  !! SPECTRUM.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------   

INTEGER :: NFRE
INTEGER :: IJ, K, M

REAL :: AKM1, TFAC
REAL, DIMENSION(SIZE(FL3,1),SIZE(FL3,3)) :: TEMP2

! ----------------------------------------------------------------------

!*    DIAGNOSTIC TAIL.
!     ----------------

NFRE = SIZE(FL3,3)
IF (SHALLOW_RUN) THEN
  DO IJ = 1,SIZE(FL3,1)
    DO M=MIJ(IJ),NFRE
      AKM1 = 1./TFAK(INDEP(IJ),M)
      TEMP2(IJ,M) = AKM1**3/TCGOND(INDEP(IJ),M)
    ENDDO
  ENDDO
ELSE
  DO IJ = 1,SIZE(FL3,1)
    DO M=MIJ(IJ),NFRE
      TEMP2(IJ,M) = FRM5(M)
    ENDDO
  ENDDO
ENDIF

DO IJ = 1,SIZE(FL3,1)
  DO M=MIJ(IJ)+1,NFRE
    TEMP2(IJ,M) = TEMP2(IJ,M)/TEMP2(IJ,MIJ(IJ))
  ENDDO
ENDDO

!*    MERGE TAIL INTO SPECTRA.
!     ------------------------

DO K=1,SIZE(FL3,2)
  DO IJ = 1,SIZE(FL3,1)
    TFAC = FL3(IJ,K,MIJ(IJ))
    DO M=MIJ(IJ)+1,NFRE
      FL3(IJ,K,M) = TEMP2(IJ,M)*TFAC
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE IMPHFTAIL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SDEPTHLIM(DEPTH, EMEAN, FL3)

! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------
!     LIMITS THE SPECTRAL VARIANCE SUCH THAT THE TOTAL VARIANCE
!     DOES NOT EXCEED THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: DEPTH (:)      !! WATER DEPTH [M].
REAL,    INTENT(INOUT) :: EMEAN (:)      !! SPECTRAL VARIANCE 
REAL,    INTENT(INOUT) :: FL3(:,:,:)     !! SPECTRUM.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !

INTEGER :: K, M
REAL, DIMENSION(SIZE(FL3,1)) :: EM

! ----------------------------------------------------------------------

EM(:)=MIN((0.25*GAMD*DEPTH(:))**2/EMEAN(:),1.0)
DO M = 1,SIZE(FL3,3)
  DO K = 1,SIZE(FL3,2)
     FL3(:,K,M) = FL3(:,K,M)*EM(:)
  ENDDO
ENDDO
EMEAN(:) = EM(:)*EMEAN(:)

END SUBROUTINE SDEPTHLIM

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE INIT_X0TAUHF

! ----------------------------------------------------------------------

!**** *INIT_X0TAUHF* -

!*    PURPOSE.
!     ---------

!     INITIALISATION FOR TAU_PHI_HF


!**   INTERFACE.
!     ----------

!       *CALL* *INIT_X0TAUHF

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------   

INTEGER :: J

REAL :: CONST1, X0, FF, F, DF

! ---------------------------------------------------------------------------- !


!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

!     find lowest limit for integration X0 *(G/USTAR)
!     ALPHA*X0**2*EXP(XKAPPA/(X0+ZALP))=1
      X0=0.005
      DO J=1,30
        FF=EXP(XKAPPA/(X0+ZALP))
        F=ALPHA*X0**2*FF-1.0
        IF (F.EQ.0.0) EXIT
        DF=ALPHA*FF*(2.0*X0-XKAPPA*(X0/(X0+ZALP))**2)
        X0=X0-F/DF
      ENDDO
      X0TAUHF=X0

      CONST1 = (BETAMAX/XKAPPA**2)/3.0

      ! Simpson Integration weights (JTOT_TAUHF must be odd) !
      WTAUHF(1)=CONST1
      DO J=2,JTOT_TAUHF-1,2
        WTAUHF(J)=4.0*CONST1
        WTAUHF(J+1)=2.0*CONST1
      ENDDO
      WTAUHF(JTOT_TAUHF)=CONST1

END SUBROUTINE INIT_X0TAUHF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TAU_PHI_HF(MIJ, USTAR, Z0, XLEVTAIL, TAUHF, PHIHF)

! ----------------------------------------------------------------------

!**** *TAU_PHI_HF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
!                                   HIGH-FREQUENCY ENERGY FLUX.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
!     JEAN BIDLOT  ECMWF  JANUARY 2017

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS AND ENERGY FLUX
!       FOR BOTH ECMWF PHYSICS AND METREO FRANCE PHYSICS.


!     METHOD.
!     -------

!       IT NEEDS A CALL TO INIT_X0TAUHF TO INITIALISE
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ---------------------------------------------------------------------------- !

!
!     INTERFACE VARIABLES.
!
!     --------------------
!
INTEGER, INTENT(IN)    :: MIJ(:)         !! LAST FREQUENCY INDEX OF THE
                                         !! PROGNOSTIC RANGE.
REAL,    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY
REAL,    INTENT(IN)    :: Z0(:)          !! ROUGHNESS LENGTH.
REAL,    INTENT(IN)    :: XLEVTAIL(:)    !! TAIL LEVEL
REAL,    INTENT(OUT)   :: TAUHF(:)       !! HIGH-FREQUENCY STRESS 
REAL,    INTENT(OUT)   :: PHIHF(:)       !! HIGH-FREQUENCY ENERGY FLUX INTO OCEAN 

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: J, IJ

REAL, PARAMETER :: ZSUP = 0.0  !  LOG(1.)
! REAL :: BETAMAX, ZALP, ALPHA
REAL :: OMEGA, OMEGAC, OMEGACC
REAL :: X0G, UST, UST0, TAUW, TAUW0
REAL :: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
REAL :: DELZ, ZINF
REAL :: FNC2, SQRTZ0OG, SQRTGZ0, GM1, GZ0, XLOGGZ0

! ----------------------------------------------------------------------

      GM1= 1.0/G

!     See INIT_X0TAUHF
      X0G=X0TAUHF*G

!*    COMPUTE THE INTEGRALS
!     ---------------------

      DO IJ=1,SIZE(USTAR)
        OMEGAC    = ZPI*FR(MIJ(IJ))
        UST0      = USTAR(IJ)
        TAUW0     = UST0**2
        GZ0       = G*Z0(IJ)
        XLOGGZ0   = LOG(GZ0)
        OMEGACC   = MAX(OMEGAC,X0G/UST0)

        SQRTZ0OG  = SQRT(Z0(IJ)*GM1)
        SQRTGZ0   = 1.0/SQRTZ0OG
        YC        = OMEGACC*SQRTZ0OG
        ZINF      = LOG(YC)
        DELZ      = MAX((ZSUP-ZINF)/REAL(JTOT_TAUHF-1),0.0)

        TAUHF(IJ)= 0.0
        PHIHF(IJ)= 0.0

        TAUW     = TAUW0
        UST      = UST0
        ! Intergrals are integrated following a change of variable : Z=LOG(Y)
        DO J=1,JTOT_TAUHF
          Y         = EXP(ZINF+REAL(J-1)*DELZ)
          OMEGA     = Y*SQRTGZ0
          CM1       = OMEGA*GM1
          ZX        = UST*CM1 +ZALP
          ZARG      = XKAPPA/ZX
          ZLOG      = XLOGGZ0+2.0*LOG(CM1)+ZARG 
          ZLOG      = MIN(ZLOG,0.0)
          ZBETA     = EXP(ZLOG)*ZLOG**4
          FNC2      = ZBETA*TAUW*WTAUHF(J)*DELZ
          TAUW      = MAX(TAUW-XLEVTAIL(IJ)*FNC2,0.0)
          UST       = SQRT(TAUW)
          TAUHF(IJ) = TAUHF(IJ) + FNC2 
          PHIHF(IJ) = PHIHF(IJ) + FNC2/Y
        ENDDO
        TAUHF(IJ) = TAUHF(IJ)
        PHIHF(IJ) = SQRTZ0OG*PHIHF(IJ)

      ENDDO

END SUBROUTINE TAU_PHI_HF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WNFLUXES (SL, SL_BOT, USTAR, UDIR, PHIAW, DEPTH, INDEP, MIJ)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)    :: SL (:, :, :)    !!  CONTRIBUTION OF ALL SOURCE
                                          !!  TERMS ACTING ON SURFACE
                                          !!  MOMENTUM AND ENERGY FLUXES.
REAL,    INTENT(IN)    :: SL_BOT (:, :, :)  !!  CONTRIBUTION OF ALL SOURCE
                                            !!  TERMS ACTING ON BOTTOM
                                            !!  MOMENTUM AND ENERGY FLUXES.
REAL,    INTENT(IN)    :: USTAR(:)        !! FRICTION VELOCITY.
REAL,    INTENT(IN)    :: UDIR (:)        !! WIND DIRECTION.
REAL,    INTENT(IN)    :: PHIAW(:)        !! ENERGY FLUX FROM WIND INTO WAVES
                                          !! OVER THE FULL FREQUENCY RANGE.
REAL,    INTENT(IN)    :: DEPTH(:)        !! ENERGY FLUX FROM WIND INTO WAVES
INTEGER, INTENT(IN)    :: INDEP(:)        !! DEPTH TABLE INDEX.
INTEGER, INTENT(IN)    :: MIJ(:)          !! LAST FREQUENCY INDEX OF THE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------

INTEGER :: IJ, K, M

REAL :: ROG, TAU, SCDFM
REAL, DIMENSION(SIZE(SL,1)) :: CM

REAL, DIMENSION(SIZE(SL,1),SIZE(SL,3)) :: CONSTFM

! ---------------------------------------------------------------------------- !
!                                                                              !
!    7.2 IF THE CALCULATION FOR THE FLUXES ARE PERFORMED, INITIAL VARIABLES.   !
!        ------------------------------------------------------------------    !

   ROG = ROWATER*G
   PHIOC(:) = PHIAW(:)
   DO IJ = 1,SIZE(SL,1)
      TAU         = ROAIR*USTAR(IJ)**2
      TAUOC_X(IJ) = TAU*SIN(UDIR(IJ))
      TAUOC_Y(IJ) = TAU*COS(UDIR(IJ))
   END DO
   PHIBOT(:) = 0.   !! BOTTOM ENERGY FLUX TO OCEAN.
   TAUBOT_X(:) = 0. !! BOTTOM MOMENTUM FLUX INTO OCEAN.
   TAUBOT_Y(:) = 0. !! BOTTOM MOMENTUM FLUX INTO OCEAN.
   CONSTFM(:,:) = 0.
   SCDFM = 0.5*DELTH*(1.-1./CO)
   DO IJ=1,SIZE(SL,1)
      DO M=1,MIJ(IJ)-1            !! CONSTFM is only defined up to M=MIJ(IJ)
         CONSTFM(IJ,M) = ROG*DFIM(M)
      ENDDO
      CONSTFM(IJ,MIJ(IJ)) = ROG*SCDFM*FR(MIJ(IJ))
   ENDDO

   IF (SHALLOW_RUN) THEN
      DO IJ = 1, SIZE(SL,1)
         FRES: DO M = 1,MIJ(IJ)
            CM(:) = TFAK(INDEP(:),M)/(ZPI*FR(M))
            DIRS: DO K = 1,KL
               PHIOC(IJ) = PHIOC(IJ)-SL(IJ,K,M)*CONSTFM(IJ,M)
               TAUOC_X(IJ) = TAUOC_X(IJ)-CM(IJ)*SL(IJ,K,M)*CONSTFM(IJ,M)*SINTH(K)
               TAUOC_Y(IJ) = TAUOC_Y(IJ)-CM(IJ)*SL(IJ,K,M)*CONSTFM(IJ,M)*COSTH(K)
               PHIBOT(IJ) = PHIBOT(IJ)-SL_BOT(IJ,K,M)*CONSTFM(IJ,M)
               TAUBOT_X(IJ) = TAUBOT_X(IJ)-CM(IJ)*SL_BOT(IJ,K,M)*CONSTFM(IJ,M)*SINTH(K)
               TAUBOT_Y(IJ) = TAUBOT_Y(IJ)-CM(IJ)*SL_BOT(IJ,K,M)*CONSTFM(IJ,M)*COSTH(K)
            END DO DIRS
         END DO FRES
      END DO
   ELSE
      DO IJ = 1, SIZE(SL,1)
         FRED: DO M = 1,MIJ(IJ)
            CM(1)  = ZPI*FR(M)/G
            DIRD: DO K = 1,KL
               PHIOC(IJ) = PHIOC(IJ)-SL(IJ,K,M)*CONSTFM(IJ,M)
               TAUOC_X(IJ) = TAUOC_X(IJ)-CM(1)*SL(IJ,K,M)*CONSTFM(IJ,M)*SINTH(K)
               TAUOC_Y(IJ) = TAUOC_Y(IJ)-CM(1)*SL(IJ,K,M)*CONSTFM(IJ,M)*COSTH(K)
            END DO DIRD
         END DO FRED
      END DO
   ENDIF

END SUBROUTINE WNFLUXES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE KERKEI (X, KER, KEI)

!**********************************************************************
! Computes the values of the zeroth order Kelvin function Ker and Kei
! These functions are used to determine the friction factor fw as a
! function of the bottom roughness length assuming a linear profile
! of eddy viscosity (See Grant and Madsen, 1979)
!**********************************************************************

      REAL :: ZR, ZI, CYR, CYI, CYR1, CYI1

      REAL :: X, KER, KEI

      ZR = X*0.50*SQRT(2.0)
      ZI = ZR
      CALL KZEONE (ZR, ZI, CYR, CYI, CYR1, CYI1)
      KER = CYR/EXP(ZR)
      KEI = CYI/EXP(ZR)
END SUBROUTINE KERKEI

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE KZEONE(X, Y, RE0, IM0, RE1, IM1)
!  June 1999 adaptation to CRESTb, all tests on range of (x,y) have been
!  bypassed, we implicitly expect X to be positive or |x,y| non zero
! 
! THE VARIABLES X AND Y ARE THE REAL AND IMAGINARY PARTS OF
! THE ARGUMENT OF THE FIRST TWO MODIFIED BESSEL FUNCTIONS
! OF THE SECOND KIND,K0 AND K1.  RE0,IM0,RE1 AND IM1 GIVE
! THE REAL AND IMAGINARY PARTS OF EXP(X)*K0 AND EXP(X)*K1,
! RESPECTIVELY.  ALTHOUGH THE REAL NOTATION USED IN THIS
! SUBROUTINE MAY SEEM INELEGANT WHEN COMPARED WITH THE
! COMPLEX NOTATION THAT FORTRAN ALLOWS, THIS VERSION RUNS
! ABOUT 30 PERCENT FASTER THAN ONE WRITTEN USING COMPLEX
! VARIABLES.
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IMPLICIT NONE
      REAL :: X, Y, X2, Y2, RE0, IM0, RE1, IM1,               &
     &        R1, R2, T1, T2, P1, P2, RTERM, ITERM, L

      REAL, PARAMETER, DIMENSION(8) :: EXSQ =                 &
     &   (/ 0.5641003087264E0,  0.4120286874989E0,            &
     &      0.1584889157959E0,  0.3078003387255E-1,           &
     &      0.2778068842913E-2, 0.1000044412325E-3,           &
     &      0.1059115547711E-5, 0.1522475804254E-8 /)

      REAL, PARAMETER, DIMENSION(8) :: TSQ =                  &
     &   (/ 0.0E0,              3.19303633920635E-1,          &
     &      1.29075862295915E0, 2.95837445869665E0,           &
     &      5.40903159724444E0, 8.80407957805676E0,           &
     &      1.34685357432515E1, 2.02499163658709E1 /)

      INTEGER :: N, M, K, LL
! THE ARRAYS TSQ AND EXSQ CONTAIN THE SQUARE OF THE
! ABSCISSAS AND THE WEIGHT FACTORS USED IN THE GAUSS-
! HERMITE QUADRATURE.
      R2 = X*X + Y*Y
      IF (R2.GE.1.96E2) GO TO 50
      IF (R2.GE.1.849E1) GO TO 30
! THIS SECTION CALCULATES THE FUNCTIONS USING THE SERIES
! EXPANSIONS
      X2 = X/2.0E0
      Y2 = Y/2.0E0
      P1 = X2*X2
      P2 = Y2*Y2
      T1 = -(LOG(P1+P2)/2.0E0+0.5772156649015329E0)
! THE CONSTANT IN THE PRECEDING STATEMENT IS EULER*S
! CONSTANT
      T2 = -ATAN2(Y,X)
      X2 = P1 - P2
      Y2 = X*Y2
      RTERM = 1.0E0
      ITERM = 0.0E0
      RE0 = T1
      IM0 = T2
      T1 = T1 + 0.5E0
      RE1 = T1
      IM1 = T2
      P2 = SQRT(R2)
      L = 2.106E0*P2 + 4.4E0
      IF (P2.LT.8.0E-1) L = 2.129E0*P2 + 4.0E0
      LL=NINT(L)
      DO N=1,LL
        P1 = N
        P2 = N*N
        R1 = RTERM
        RTERM = (R1*X2-ITERM*Y2)/P2
        ITERM = (R1*Y2+ITERM*X2)/P2
        T1 = T1 + 0.5E0/P1
        RE0 = RE0 + T1*RTERM - T2*ITERM
        IM0 = IM0 + T1*ITERM + T2*RTERM
        P1 = P1 + 1.0E0
        T1 = T1 + 0.5E0/P1
        RE1 = RE1 + (T1*RTERM-T2*ITERM)/P1
        IM1 = IM1 + (T1*ITERM+T2*RTERM)/P1
      END DO
      R1 = X/R2 - 0.5E0*(X*RE1-Y*IM1)
      R2 = -Y/R2 - 0.5E0*(X*IM1+Y*RE1)
      P1 = EXP(X)
      RE0 = P1*RE0
      IM0 = P1*IM0
      RE1 = P1*R1
      IM1 = P1*R2
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE INTEGRAL
! REPRESENTATION, EQN 3, EVALUATED WITH 15 POINT GAUSS-
! HERMITE QUADRATURE
   30 X2 = 2.0E0*X
      Y2 = 2.0E0*Y
      R1 = Y2*Y2
      P1 = SQRT(X2*X2+R1)
      P2 = SQRT(P1+X2)
      T1 = EXSQ(1)/(2.0E0*P1)
      RE0 = T1*P2
      IM0 = T1/P2
      RE1 = 0.0E0
      IM1 = 0.0E0
      DO N=2,8
        T2 = X2 + TSQ(N)
        P1 = SQRT(T2*T2+R1)
        P2 = SQRT(P1+T2)
        T1 = EXSQ(N)/P1
        RE0 = RE0 + T1*P2
        IM0 = IM0 + T1/P2
        T1 = EXSQ(N)*TSQ(N)
        RE1 = RE1 + T1*P2
        IM1 = IM1 + T1/P2
      END DO
      T2 = -Y2*IM0
      RE1 = RE1/R2
      R2 = Y2*IM1/R2
      RTERM = 1.41421356237309E0*COS(Y)
      ITERM = -1.41421356237309E0*SIN(Y)
! THE CONSTANT IN THE PREVIOUS STATEMENTS IS,OF COURSE,
! SQRT(2.0).
      IM0 = RE0*ITERM + T2*RTERM
      RE0 = RE0*RTERM - T2*ITERM
      T1 = RE1*RTERM - R2*ITERM
      T2 = RE1*ITERM + R2*RTERM
      RE1 = T1*X + T2*Y
      IM1 = -T1*Y + T2*X
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE
! ASYMPTOTIC EXPANSIONS
   50 RTERM = 1.0E0
      ITERM = 0.0E0
      RE0 = 1.0E0
      IM0 = 0.0E0
      RE1 = 1.0E0
      IM1 = 0.0E0
      P1 = 8.0E0*R2
      P2 = SQRT(R2)
      L = 3.91E0+8.12E1/P2
      LL=NINT(L)
      R1 = 1.0E0
      R2 = 1.0E0
      M = -8
      K = 3
      DO N=1,LL
        M = M + 8
        K = K - M
        R1 = FLOAT(K-4)*R1
        R2 = FLOAT(K)*R2
        T1 = FLOAT(N)*P1
        T2 = RTERM
        RTERM = (T2*X+ITERM*Y)/T1
        ITERM = (-T2*Y+ITERM*X)/T1
        RE0 = RE0 + R1*RTERM
        IM0 = IM0 + R1*ITERM
        RE1 = RE1 + R2*RTERM
        IM1 = IM1 + R2*ITERM
      END DO
      T1 = SQRT(P2+X)
      T2 = -Y/T1
      P1 = 8.86226925452758E-1/P2
! THIS CONSTANT IS SQRT(PI)/2.0, WITH PI=3.14159...
      RTERM = P1*COS(Y)
      ITERM = -P1*SIN(Y)
      R1 = RE0*RTERM - IM0*ITERM
      R2 = RE0*ITERM + IM0*RTERM
      RE0 = T1*R1 - T2*R2
      IM0 = T1*R2 + T2*R1
      R1 = RE1*RTERM - IM1*ITERM
      R2 = RE1*ITERM + IM1*RTERM
      RE1 = T1*R1 - T2*R2
      IM1 = T1*R2 + T2*R1

END SUBROUTINE KZEONE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WSIGSTAR (USTAR, Z0, WSTAR, SIG_N)

! ----------------------------------------------------------------------

!**** *WSIGSTAR* - COMPUTATION OF THE RELATIVE STANDARD DEVIATION OF USTAR.

!*    PURPOSE.
!     ---------

!     COMPUTES THE STANDARD DEVIATION OF USTAR DUE TO SMALL SCALE GUSTINESS
!     RELATIVE TO USTAR

!     METHOD.
!     -------

!     USE PANOFSKY (1991) TO EXPRESS THE STANDARD DEVIATION OF U10 IN TERMS
!     USTAR AND (Zi/L) THE MONIN-OBOKHOV LENGTH (Zi THE INVERSION HEIGHT).
!     (but with the background gustiness set to 0.)
!     and USTAR=SQRT(Cd)*U10 to DERIVE THE STANDARD DEVIATION OF USTAR.
!     WITH CD=A+B*U10 (see below).

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     REFERENCE.
!     ----------

!     SEE SECTION 3.2.1 OF THE WAM DOCUMENTATION.
!     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)

! ---------------------------------------------------------------------------- !

!
!     INTERFACE VARIABLES.
!
!     --------------------
!

REAL,    INTENT(IN)    :: USTAR(:) !! FRICTION VELOCITY
REAL,    INTENT(IN)    :: Z0(:)    !! ROUGHNESS LENGTH.
REAL,    INTENT(IN)    :: WSTAR(:) !! FREE CONVECTION VELOCITY SCALE.
REAL,    INTENT(OUT)   :: SIG_N(:) !! RELATIVE STANDARD DEVIATION OF USTAR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: WSPMIN = 1.0 ! MINIMUM WIND SPEED
REAL, PARAMETER :: BG_GUST = 0.0 ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
REAL, PARAMETER :: ONETHIRD = 1.0/3.0
REAL, PARAMETER :: LOG10 = LOG(10.0)
REAL, PARAMETER :: C1 = 1.03E-3
REAL, PARAMETER :: C2 = 0.04E-3
REAL, PARAMETER :: P1 = 1.48
REAL, PARAMETER :: P2 = -0.21

INTEGER :: IJ

REAL :: U10, C_D, DC_DDU, SIG_CONV
REAL :: XKAPPAD
REAL :: U10M1, C2U10P1, U10P2

! ---------------------------------------------------------------------------- !

      XKAPPAD=1.0/XKAPPA

      DO IJ=1,SIZE(USTAR)
!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
        U10 = USTAR(IJ)*XKAPPAD*(LOG10-LOG(Z0(IJ)))
        U10 = MAX(U10,WSPMIN)
        U10M1=1.0/U10
        C2U10P1=C2*U10**P1
        U10P2=U10**P2
        C_D = (C1 + C2U10P1)*U10P2
        DC_DDU = (P2*C1+(P1+P2)*C2U10P1)*U10P2*U10M1
        SIG_CONV = 1.0 + 0.5*U10/C_D*DC_DDU
        SIG_N(IJ) = MIN(0.1, SIG_CONV * U10M1*(BG_GUST*USTAR(IJ)**3 + &
     &                  0.5*XKAPPA*WSTAR(IJ)**3)**ONETHIRD )
      ENDDO

END SUBROUTINE WSIGSTAR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !


END MODULE WAM_SOURCE_MODULE

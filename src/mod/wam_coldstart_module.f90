MODULE WAM_COLDSTART_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL DATA WHICH ARE USED FOR A COLD START.             !
!   ALL PROCEDURES ARE INCLUDED TO COMPUTE THE WAM MODEL COLDSTART FILEDS.     !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE, ONLY: &
&       ABORT1                     !! ABORT PROCESSING.
USE WAM_WIND_MODULE,    ONLY: &
&       WAM_WIND                   !! READ AND PREPARE WINDS.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE, ONLY: G, PI, ZPI, DEG, RAD

USE WAM_FRE_DIR_MODULE, ONLY: FR, TH, ML, KL

USE WAM_GRID_MODULE,    ONLY: NSEA, DELPHI

USE WAM_MODEL_MODULE,   ONLY: FL3, U10, UDIR

USE WAM_TIMOPT_MODULE,  ONLY: CDATEA, CDATEE, CDTPRO, CDTSOU,                 &
&                             CDA, CDTA, CDCA, COLDSTART

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST

USE WAM_WIND_MODULE,    ONLY: IDELWO, IDELWI

use wam_mpi_module,     only: nijs, nijl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.0 START OPTION.                                                        !
!         -------------                                                        !

INTEGER :: IOPTI = -1   !! =  0 WIND INDEPENDENT INITIAL VALUES.                              
                        !! =  1 WIND DEPENDENT INITIAL VALUES AND                             
                        !!      ENERGY EQUAL ZERO IF WINDSPEED IS ZERO                        
                        !! =  2 WIND DEPENDENT INITIAL VALUES AND                             
                        !!      ENERGY COMPUTED FROM GIVEN PARAMETERS IF                      
                        !!      WINDSPEED IS ZERO. 
! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.0 JONSWAP PARAMETERS AND FETCH.                                        !
!         -----------------------------                                        !

REAL, ALLOCATABLE, DIMENSION(:) :: FP    !! PEAK FREQUENCY IN A BLOCK [HZ].
REAL, ALLOCATABLE, DIMENSION(:) :: ALPHJ !! ALPHA PARAMETER IN A BLOCK.
REAL, ALLOCATABLE, DIMENSION(:) :: THES  !! MEAN DIRECTION IN A BLOCK [RAD].

REAL :: FM = 0.      !! PEAK FREQUENCY [HZ].
REAL :: ALFA = 0.    !! ALPHA PARAMETER.
REAL :: GAMMA = 0.   !! OVERSHOOT FACTOR.
REAL :: SA = 0.      !! LEFT PEAK WIDTH.
REAL :: SB = 0.      !! RIGHT PEAK WIDTH.
REAL :: THETAQ = 0.  !! MEAN DIRECTION [RAD].

REAL :: FETCH = 0.   !! FETCH [M].

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE PREPARE_COLDSTART           !! PREPARES COLDSTART FIELDS FOR WAMODEL.
   MODULE PROCEDURE PREPARE_COLDSTART
END INTERFACE
PUBLIC PREPARE_COLDSTART

INTERFACE PRINT_COLDSTART_STATUS      !! PRINTS COLDSTART MODULE.
   MODULE PROCEDURE PRINT_COLDSTART_STATUS
END INTERFACE
PUBLIC PRINT_COLDSTART_STATUS

INTERFACE SET_C_START_PAR             !! TRANSFER INITIAL PARAMETER TO MODULE.
   MODULE PROCEDURE SET_C_START_PAR
END INTERFACE
PUBLIC SET_C_START_PAR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE FETCH_LAW              !! COMPUTE JONSWAP PARAMETERS FROM FETCH LAW.
   MODULE PROCEDURE FETCH_LAW
END INTERFACE
PRIVATE FETCH_LAW

INTERFACE SPECTRA                !! COMPUTATION OF 2-D SPECTRA.
   MODULE PROCEDURE SPECTRA
END INTERFACE
PRIVATE SPECTRA

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_COLDSTART

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PREPARE_COLDSTART - PREPARES COLDSTART FIELDS FOR WAMODEL.                 !
!                                                                              !
!      H. GUNTHER    ECMWF     MAY 1990                                        !
!      H. GUNTHER    ECMWF     DECEMBER 90  MODIFIED FOR CYCLE_4.              !
!      H. GUNTHER    GKSS      DECEMBER 2000  FT90.                            !
!      A. Behrens    MSC/GKSS  December 2003  Message passing                  !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO GENERATE WAMODEL START FIELDS.                                      !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!    REFERENCE.                                                                !
!    ----------                                                                !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: ij

real, allocatable, dimension (:) :: rudir
real, allocatable, dimension (:) :: ru10

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DEFINE MODEL DATES.                                                   !
!        -------------------                                                   !

CDTPRO  = CDATEA  !! PROPAGATION
CDTSOU  = CDATEA  !! SOURCE FUNCTION
CDTA    = ' '     !! DEPTH
CDCA    = ' '     !! CURRENTS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. PREPARE FIRST WINDFIELD.                                              !
!        ------------------------                                              !

allocate (ru10(nsea), rudir(nsea))
CDA = CDATEA
CALL WAM_WIND (ru10, rudir, CDA)

u10(nijs:nijl)  = ru10(nijs:nijl)
udir(nijs:nijl) = rudir(nijs:nijl)

deallocate (ru10, rudir)

IF (ITEST.GE.3) THEN
      WRITE (IU06,*) '      SUB. PREPARE_COLDSTART: WAM_WIND DONE'
      WRITE (IU06,*) '        FIRST WIND FIELD PROCESSED'
      WRITE (IU06,*) '        DATE IS ................ CDA = ', CDA
      WRITE (IU06,*) '        FIELD IS SAVED IN WAM_MODEL_MODULE'
      WRITE (IU06,*) '  '
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. INITIAL PARAMETERS DUE TO OPTION.                                     !
!        ---------------------------------                                     !

IF (ALLOCATED(FP   )) DEALLOCATE(FP)
IF (ALLOCATED(ALPHJ)) DEALLOCATE(ALPHJ)
IF (ALLOCATED(THES )) DEALLOCATE(THES)
   
ALLOCATE (FP(nijs:nijl), ALPHJ(nijs:nijl), THES(nijs:nijl))

IF (IOPTI.EQ.0) THEN
   FP = FM
   ALPHJ = ALFA
   THES = THETAQ
ELSE IF (IOPTI.EQ.1) THEN
   FP = FM
   ALPHJ = 0.
   do ij=nijs,nijl
      THES(ij)  = UDIR(ij)
   enddo
ELSE IF (IOPTI.EQ.2) THEN
   FP = FM
   ALPHJ = ALFA
   do ij= nijs,nijl
      if (u10(ij)>0.1E-08) then
         thes(ij) = udir(ij)
      else
         thes(ij) = thetaq
      endif
   enddo
END IF

IF (IOPTI.NE.0) THEN
   IF (FETCH.LT.0.1E-5) FETCH = 0.5*DELPHI
   CALL FETCH_LAW (FETCH, FM, U10)
   IF (ITEST.GE.3) THEN
      WRITE (IU06,*) '      SUB. PREPARE_COLDSTART: FETCH_LAW DONE'
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTE INITIAL SPECTRA FROM PARAMETERS.                              !
!        ----------------------------------------                              !

CALL SPECTRA (FL3)
IF (ITEST.GE.3) WRITE (IU06,*) '      SUB. PREPARE_COLDSTART: SPECTRA DONE'

DEALLOCATE (FP, ALPHJ, THES)

END SUBROUTINE  PREPARE_COLDSTART

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_COLDSTART_STATUS

WRITE(IU06,*) '  '
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '              MODEL COLDSTART STATUS:'
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '  '

IF (.NOT.COLDSTART) THEN
   WRITE (IU06,*) ' THIS IS A RESTART RUN:'
   WRITE (IU06,*) ' THEREFORE THE DATA STORED IN THIS MODULE ARE NOT USED.'   
   RETURN
END IF
WRITE(IU06,*) '  '
IF (IOPTI.EQ.0) THEN
   WRITE (IU06,*) ' COLD START VALUES ARE COMPUTED FROM JONSWAP PARAMETERS.'
   WRITE (IU06,*) '    PHILLIPS  PARAMETER ALPHA = ', ALFA
   WRITE (IU06,*) '    PEAK FREQUENCY         FM = ', FM, 'HZ'
   WRITE (IU06,*) '    OVERSHOOT FACTOR    GAMMA = ', GAMMA
   WRITE (IU06,*) '    LEFT PEAK WIDTH        SA = ', SA
   WRITE (IU06,*) '    RIGHT PEAK WIDTH       SB = ', SB
   WRITE (IU06,*) '    MEAN DIRECTION     THETAQ = ', THETAQ*DEG, 'DEG'
ELSE IF (IOPTI.EQ.1) THEN
   WRITE (IU06,*) ' COLD START VALUES ARE COMPUTED FROM FETCH LAW.'
   WRITE (IU06,*) '    FETCH IS            FETCH = ', FETCH, 'M'
   WRITE (IU06,*) '    MAXIMUM PEAK FREQUENCY FM = ', FM, 'HZ'
   WRITE (IU06,*) '    WAVE ENERGY IS ZERO IN CALM WIND AREAS.'
ELSE IF (IOPTI.EQ.2) THEN
   WRITE (IU06,*) ' COLD START VALUES ARE COMPUTED FROM FETCH LAW.'
   WRITE (IU06,*) '    FETCH IS            FETCH = ', FETCH, 'M'
   WRITE (IU06,*) '    MAXIMUM PEAK FREQUENCY FM = ', FM, 'HZ'
   WRITE (IU06,*) '    JONSWAP PARAMETERS USED IN CALM WIND AREAS:'
   WRITE (IU06,*) '    PHILLIPS  PARAMETER ALPHA = ', ALFA
   WRITE (IU06,*) '    PEAK FREQUENCY         FM = ', FM, 'HZ'
   WRITE (IU06,*) '    OVERSHOOT FACTOR    GAMMA = ', GAMMA
   WRITE (IU06,*) '    LEFT PEAK WIDTH        SA = ', SA
   WRITE (IU06,*) '    RIGHT PEAK WIDTH       SB = ', SB
   WRITE (IU06,*) '    MEAN DIRECTION     THETAQ = ', THETAQ*DEG, 'DEG'
ELSE
   WRITE (IU06,*) ' COLD START MODULE IS NOT INITIALIZED.'
END IF
WRITE(IU06,*) '  '

END SUBROUTINE PRINT_COLDSTART_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_C_START_PAR (OPTION, ALPHA, FPEAK, GA, SIGA, SIGB, DIR, FET)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)           :: OPTION  !! MODEL START OPTION.
REAL,    INTENT(IN)           :: ALPHA   !! ALPHA PARAMETER.
REAL,    INTENT(IN)           :: FPEAK   !! PEAK FREQUENCY [HZ].
REAL,    INTENT(IN)           :: GA      !! OVERSHOOT FACTOR.
REAL,    INTENT(IN)           :: SIGA    !! LEFT PEAK WIDTH.
REAL,    INTENT(IN)           :: SIGB    !! RIGTH PEAK WIDTH.
REAL,    INTENT(IN)           :: DIR     !! MEAN WAVE DIRECTION [DEG].
REAL,    INTENT(IN), OPTIONAL :: FET     !! FETCH LENGTH [M].

! ---------------------------------------------------------------------------- !

IOPTI  = OPTION

IF (IOPTI.LT.0 .OR. IOPTI.GT.2) THEN
   WRITE(IU06,*) ' ************************************************'
   WRITE(IU06,*) ' *                                              *'
   WRITE(IU06,*) ' *     FATAL ERROR IN SUB. SET_START_OPTION     *'
   WRITE(IU06,*) ' *     ====================================     *'
   WRITE(IU06,*) ' * OPTION FOR INITIAL VALUES IS OUT OF RANGE.   *'
   WRITE(IU06,*) ' * OPTION IS       IOPTI = ', IOPTI
   WRITE(IU06,*) ' * POSSIBLE COLD START OPTIONS ARE: 0, 1, OR 2  *'
   WRITE(IU06,*) ' *                                              *'
   WRITE(IU06,*) ' *      PROGRAM ABORTS  PROGRAM ABORTS          *'
   WRITE(IU06,*) ' *                                              *'
   WRITE(IU06,*) ' ************************************************'
   CALL ABORT1
END IF

ALFA    = ALPHA
FM      = FPEAK
GAMMA   = GA
SA      = SIGA
SB      = SIGB
THETAQ  = DIR*RAD
IF (PRESENT(FET)) THEN
   FETCH = FET
ELSE
   FETCH = 0.
END IF

END SUBROUTINE SET_C_START_PAR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FETCH_LAW (FETCH, FPMAX, U10)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   FETCH_LAW - COMPUTE JONSWAP PARAMETERS FROM FETCH LAW.                     !
!                                                                              !
!     S. HASSELMANN  - JULY 1990                                               !
!     H. GUNTHER     - DECEMBER 1990   MODIFIED FOR CYCLE_4.                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTES FOR EACH GRID POINT THE PEAK FREQUENCY FROM A FETCH LAW       !
!       AND THE JONSWAP ALPHA FROM THE ALPHA NY RELATION.                      !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       FP = A * (G*FETCH/U_10**2)**D    A = 2.84                              !
!       FP = MAX [FP, 0.13]              D = -3./10.                           !
!       FP = MIN [FP, FRMAX*U_10/G]                                            !
!       ALPHJ = B * FP**2/3              B = 0.033                             !
!       ALPHJ = MAX [ALPHJ, 0.0081]                                            !
!       FP = G/U_10*FP                                                         !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       K.HASSELMAN,D.B.ROOS,P.MUELLER AND W.SWELL                             !
!          A PARAMETRIC WAVE PREDICTION MODEL                                  !
!          JOURNAL OF PHSICAL OCEANOGRAPHY, VOL. 6, NO. 2, MARCH 1976.         !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN) :: FETCH          !! FETCH TO BE USED (METRES).
REAL,    INTENT(IN) :: FPMAX          !! MAXIMUM PEAK FREQUENCY (HERTZ).
REAL,    INTENT(IN) :: U10(:)         !! MODULUS OF WIND VELOCITY [M/S].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: A = 2.84,  D = -(3./10.) !! PEAKFREQUENCY FETCH LAW CONSTANTS
REAL, PARAMETER :: B = 0.033, E = 2./3.     !! ALPHA-PEAKFREQUENCY LAW CONSTANTS

REAL :: UG(SIZE(U10))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE VALUES FROM FETCH LAWS.                                       !
!        -------------------------------                                       !

WHERE (U10 .GT. 0.1E-08)
   UG = G/U10
   FP = MAX(0.13, A*((G*FETCH)/(U10**2))**D)
   FP = MIN(FP, FPMAX/UG)
   ALPHJ = MAX(0.0081, B * FP**E)
   FP = FP*UG
END WHERE

END SUBROUTINE FETCH_LAW

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SPECTRA (FL3)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SPECTRA  - COMPUTATION OF 2-D SPECTRA FOR ONE BLOCK.                       !
!                                                                              !
!     S. HASSELMANN  - JULY 1990                                               !
!     H. GUNTHER     - DECEMBER 1990   MODIFIED FOR CYCLE_4.                   !
!     H. GUNTHER     - DECEMBER 2001   TMA SCALING.                            !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       INITIALISATION OF A BLOCK BY 2-D SPECTRA.                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       1-D JONSWAP SPECTRA AND COSINE**2 SPREADING FUNCTIONS ARE              !
!       COMPUTED FROM GIVEN WINDS AND PARAMETERS AT EACH GRID POINT.           !
!       THE 1-D SPECTRA ARE SPREAD OVER THE DIRECTIONS BY MULTIPLICATION       !
!       WITH THE SPREADING FUNCTION.                                           !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       K.HASSELMAN,D.B.ROSS,P.MUELLER AND W.SWELL                             !
!          A PARAMETRIC WAVE PREDICTION MODEL                                  !
!          JOURNAL OF PHYSICAL OCEANOGRAPHY, VOL. 6, NO. 2, MARCH 1976         !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL, INTENT(OUT) :: FL3(:,:,:)      !! block of 2-d spectr.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: ZDP=2./PI
REAL, PARAMETER :: PIRHF = G**2/(ZPI**4)  !! PM SPECTRUM FACTOR

INTEGER :: M, K
REAL    :: FRH, SIGMA(SIZE(FL3,1)), ARG(SIZE(FL3,1))
REAL    :: ST(SIZE(FL3,1),SIZE(FL3,2)),                          &
&          ET(SIZE(FL3,1),SIZE(FL3,3))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE JONSWAP SPECTRUM.                                             !
!        -------------------------                                             !

ET = 0.
DO M = 1,ML
   FRH = FR(M)
   ARG = 1.25*(FP/FRH)**4
   WHERE (ARG.LT.99.) ET(:,M) = ALPHJ*PIRHF/FRH**5*EXP(-ARG)
   SIGMA = SA
   WHERE (FRH.GT.FP) SIGMA = SB
   ARG = .5*((FRH-FP) / (SIGMA*FP))**2
   
   WHERE (ARG.LT.99.) ET(:,M) = ET(:,M)*exp(log(GAMMA)*EXP(-ARG))
!AB   WHERE (ARG.LT.99.) ET(:,M) = ET(:,M)*GAMMA**EXP(-ARG)
    
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTATION OF SPREADING FUNCTION.                                    !
!        ----------------------------------                                    !

DO K = 1,KL
   ST(:,K) = MAX(0. ,COS(TH(K)-THES))
END DO
ST = ZDP*ST**2
WHERE (ST.LT.0.1E-08) ST = 0.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTATION OF 2-D SPECTRUM.                                          !
!        ----------------------------                                          !

DO M = 1,ML
   DO K = 1,KL
      FL3(:,K,M) = ET(:,M) * ST(:,K)
   END DO
END DO

END SUBROUTINE SPECTRA

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_COLDSTART_MODULE

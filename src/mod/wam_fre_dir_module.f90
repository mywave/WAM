MODULE WAM_FRE_DIR_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS THE SET-UP OF THE FREQUENCY-DIRECTION GRID AND        !
!   THE PRECOMPUTED FREQUENY DIRECTION DEPENDENT CONSTANTS OF THE WAM MODEL.   !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1                      !! TERMINATES PROCESSING.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST
USE WAM_GENERAL_MODULE, ONLY: G, PI, ZPI, DEG, R

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!    1. FREQUENCY AND DIRECTION GRID.                                          !

INTEGER            :: ML = -1        !! NUMBER OF FREQUENCIES
INTEGER            :: KL = -1        !! NUMBER OF DIRECTIONS
REAL               :: DELTH = -1     !! ANGULAR INCREMENT OF SPECTRUM [RAD].
REAL               :: DELTR = -1     !! DELTH TIMES RADIUS OF EARTH [RAD*M].
REAL,    PARAMETER :: CO = 1.1       !! FREQUENCY RATIO.
REAL,    SAVE      :: INV_LOG_CO     !! 1./LOG10(FREQUENCY RATIO).

INTEGER, PARAMETER :: EX_TAIL = -5   !! TAIL FREQUENCY EXPONENT.
REAL,    SAVE      :: MO_TAIL        !! MOMENT  O TAIL FACTOR.
REAL,    SAVE      :: MM1_TAIL       !! MOMENT -1 TAIL FACTOR.
REAL,    SAVE      :: MP1_TAIL       !! MOMENT +1 TAIL FACTOR.
REAL,    SAVE      :: MP2_TAIL       !! MOMENT +2 TAIL FACTOR.

REAL,    SAVE      :: COEF4 = 3.0E-07 !! COEFFICIENT USED TO COMPUTE THE 
                                      !! SPECTRAL LIMITER IN IMPLSCH.
REAL,    SAVE      :: FMIN            !! MINIMUM ENERGY DENSITY = HS = 7CM

REAL,ALLOCATABLE :: FR(:)        !! FREQUENCIES [HZ].
REAL,ALLOCATABLE :: DF(:)        !! FREQUENCY INTERVAL [HZ].
REAL,ALLOCATABLE :: DF_FR(:)     !! DF*FR [HZ**2].
REAL,ALLOCATABLE :: DF_FR2(:)    !! DF*FR*FR [HZ**3].
REAL,ALLOCATABLE :: DFIM_FR(:)   !! DFIM*FR [HZ*HZ*RAD].
REAL,ALLOCATABLE :: DFIM_FR2(:)  !! DFIM*FR**2  [HZ**3*RAD].
REAL,ALLOCATABLE :: FR5(:)       !! FR(M)**5
REAL,ALLOCATABLE :: FRM5(:)      !! 1./FR(M))**5
REAL,ALLOCATABLE :: DFIM(:)      !! FREQUENCY INTERVAL*DIRECTION INTER.
REAL,ALLOCATABLE :: DFIMOFR(:)   !! DFIM/FR
REAL,ALLOCATABLE :: DFFR(:)      !! DFIM*FR
REAL,ALLOCATABLE :: DFFR2(:)     !! DFIM*FR**2 

REAL,ALLOCATABLE, DIMENSION(:) :: GOM   !! DEEP WATER GROUP VELOCITIES [M/S].
REAL,ALLOCATABLE, DIMENSION(:) :: C     !! DEEP WATER PHASE VELOCITIES [M/S].
REAL,ALLOCATABLE, DIMENSION(:) :: TH    !! DIRECTIONS IN RADIANS.
REAL,ALLOCATABLE, DIMENSION(:) :: COSTH !! COS OF DIRECTION.
REAL,ALLOCATABLE, DIMENSION(:) :: SINTH !! SIN OF DIRECTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!    2. SHALLOW WATER TABLES.                                                  !

INTEGER, PARAMETER :: NDEPTH = 69   !! LENGTH OF SHALLOW WATER TABLES.
REAL,    PARAMETER :: DEPTHA = 1.0  !! MINIMUM DEPTH FOR TABLES [M].
REAL,    PARAMETER :: DEPTHD = 1.1  !! DEPTH RATIO.

REAL,   ALLOCATABLE, DIMENSION(:,:) :: TCGOND  !! SHALLOW WATER GROUP VELOCITY.
REAL,   ALLOCATABLE, DIMENSION(:,:) :: TFAK    !! SHALLOW WATER WAVE NUMBER.
REAL,   ALLOCATABLE, DIMENSION(:,:) :: TSIHKD  !! TABLE FOR OMEGA/SINH(2KD).
REAL,   ALLOCATABLE, DIMENSION(:,:) :: TFAC_ST !! TABLE FOR 2*G*K**2/
!                                                           (OMEGA*TANH(2KD)).

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE SET_FRE_DIR                !! TRANSFERS FREQUENY AND DIRECTION
   MODULE PROCEDURE SET_FRE_DIR      !! GRID DEFINITIONS TO MODULE
END INTERFACE
PUBLIC SET_FRE_DIR

INTERFACE PRINT_FRE_DIR_STATUS       !! PRINTS WAM_FRE_DIR_MODULE DATA.
   MODULE PROCEDURE PRINT_FRE_DIR_STATUS
END INTERFACE
PUBLIC PRINT_FRE_DIR_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE MAKE_FRE_DIR               !! COMPUTE FREQUENCY-DIRECTION ARRAYS.
   MODULE PROCEDURE MAKE_FRE_DIR
END INTERFACE
PRIVATE MAKE_FRE_DIR

INTERFACE MAKE_SHALLOW_TABLES        !! COMPUTE TABLES FOR SHALLOW WATER.
   MODULE PROCEDURE MAKE_SHALLOW_TABLES
END INTERFACE
PRIVATE MAKE_SHALLOW_TABLES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_FRE_DIR (N_DIR, N_FRE, FR1)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_FRE_DIR - ROUTINE TO PREPARE WAM_FRE_DIR MODULE.                       !
!                                                                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO COMPUTE ALL VARAIABLES IN WAM_FRE_DIR MODULE.                       !
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

INTEGER , INTENT(IN) :: N_DIR   !! NUMBER OF DIRECTIONS.
INTEGER , INTENT(IN) :: N_FRE   !! NUMBER OF FREQUENCIES.
REAL    , INTENT(IN) :: FR1     !! FIRST FREQUENCY [HZ]

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLEAR FRE_DIR GRID.                                                   !
!        -------------------                                                   !

KL = -1
ML = -1

IF (ALLOCATED (FR)      ) DEALLOCATE (FR)
IF (ALLOCATED (GOM)     ) DEALLOCATE (GOM)
IF (ALLOCATED (C)       ) DEALLOCATE (C)
IF (ALLOCATED (TH)      ) DEALLOCATE (TH)
IF (ALLOCATED (COSTH)   ) DEALLOCATE (COSTH)
IF (ALLOCATED (SINTH)   ) DEALLOCATE (SINTH)
IF (ALLOCATED (DF)      ) DEALLOCATE (DF)
IF (ALLOCATED (DF_FR)   ) DEALLOCATE (DF_FR)
IF (ALLOCATED (DF_FR2)  ) DEALLOCATE (DF_FR2)

IF (ALLOCATED (DFIM)    ) DEALLOCATE (DFIM)
IF (ALLOCATED (DFIMOFR) ) DEALLOCATE (DFIMOFR)
IF (ALLOCATED (DFIM_FR) ) DEALLOCATE (DFIM_FR)
IF (ALLOCATED (DFIM_FR2)) DEALLOCATE (DFIM_FR2)

IF (ALLOCATED (FR5)     ) DEALLOCATE (FR5)
IF (ALLOCATED (FRM5)    ) DEALLOCATE (FRM5)

IF (ALLOCATED (TCGOND)  ) DEALLOCATE (TCGOND)
IF (ALLOCATED (TFAK)    ) DEALLOCATE (TFAK)
IF (ALLOCATED (TSIHKD)  ) DEALLOCATE (TSIHKD)
IF (ALLOCATED (TFAC_ST) ) DEALLOCATE (TFAC_ST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK INPUT DATA.                                                     !
!        -----------------                                                     !

IF (N_DIR.LE.0 .OR. N_FRE.LE.0 .OR. FR1.LE.0.) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *          FATAL ERROR IN SUB. SET_FRE_DIR               *'
   WRITE (IU06,*) ' *          ===============================               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' * FREQUENCY DIRECTION GRID IS NOT DEFINED.               *'
   WRITE (IU06,*) ' * NUMBER OF DIRECTIONS IS N_DIR = ', N_DIR
   WRITE (IU06,*) ' * NUMBER OF FREUENCIES IS N_FRE = ', N_FRE
   WRITE (IU06,*) ' * FIRST FREQUENCY [HZ]: FR1                              *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COPY FRE_DIR GRID DEFINITIONS.                                        !
!        ------------------------------                                        !

ML = N_FRE
KL = N_DIR
ALLOCATE (FR(ML))
FR(1) = FR1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. FREQUENCY-DIRECTION GRID.                                             !
!        -------------------------                                             !

CALL MAKE_FRE_DIR
IF (ITEST.GT.1) WRITE (IU06,*) ' SUB. MAKE_FRE_DIR DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. SHALLOW WATER TABLES.                                                 !
!        ---------------------                                                 !

CALL MAKE_SHALLOW_TABLES
IF (ITEST.GT.1) WRITE (IU06,*) ' SUB. MAKE_SHALLOW_TABLES DONE'

END SUBROUTINE SET_FRE_DIR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_FRE_DIR_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PRINT_FRE_DIRT_STATUS - PRINT STATUS OF WAM_FRE_DIR_MODULE.                !
!                                                                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000                           !
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

INTEGER, PARAMETER :: NAN  = 10 !! STEPS FOR SHALLOW WATER TABLE PRINT.
INTEGER :: K, I

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. FREQUENCY AND DIRECTION GRID.                                         !
!        -----------------------------                                         !

WRITE (IU06,'(/,'' ----------------------------------------'')')
WRITE (IU06,'(  ''  FREQUENCY AND DIRECTION MODULE STATUS'')')
WRITE (IU06,'(  '' ----------------------------------------'')')
WRITE (IU06,*)' '
IF (KL.GT.0 .AND. ML.GT.0) THEN
   WRITE (IU06,'('' NUMBER OF DIRECTIONS  IS  KL = '',I3)') KL
   WRITE (IU06,'('' NUMBER OF FREQUENCIES IS  ML = '',I3)') ML
   WRITE (IU06,'('' FIRST FREQUENCY IS     FR(1) = '',F10.8,'' [HZ]'')') FR(1)
ELSE
   WRITE (IU06,*) '  INPUT DATA TO MODULE ARE NOT DEFINED'
END IF

IF (ALLOCATED(DFIM)) THEN
   WRITE (IU06,'('' MODEL FREQUENCIES IN HERTZ:'')')
   WRITE (IU06,'(1X,13F10.5)') FR(1:ML)
   WRITE (IU06,'(/,'' MODEL DIRECTIONS IN DEGREE (CLOCKWISE FROM NORTH):'')')
   WRITE (IU06,'(1X,13F10.5)') TH(1:KL)*DEG
   IF (ITEST.GT.0) THEN
      WRITE (IU06,'(/,'' MODEL FREQUENCY INTERVALLS TIMES DIRECTION'',         &
&              '' INTERVALL IN HERTZ*RADIENS'')')
      WRITE (IU06,'(1X,13F10.5)') DFIM(1:ML)
      WRITE (IU06,'(/,'' MODEL DEEP WATER GROUPVELOCITY IN M/S:'')')
      WRITE (IU06,'(1X,13F10.5)') GOM(1:ML)
      WRITE (IU06,'(/,'' MODEL DEEP WATER PHASEVELOCITY IN M/S:'')')
      WRITE (IU06,'(1X,13F10.5)') C(1:ML)
   END IF
ELSE
   WRITE (IU06,*) '  MODULE DATA ARE NOT PREPARED'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. SHALLOW WATER TABLES.                                                 !
!        ---------------------                                                 !

WRITE (IU06,'(/,'' ----------------------------------------'')')
WRITE (IU06,'(  ''            SHALLOW WATER TABLES'')')
WRITE (IU06,'(  '' ----------------------------------------'',/)')
K = MAX(NDEPTH/NAN,1)
WRITE (IU06,'(''  LOGARITHMIC DEPTH FROM: DEPTHA = '',F5.1,'' TO DEPTHE  = '', &
&             F5.1, ''IN STEPS OF DEPTHD = '',F5.1)')                          &
&             DEPTHA, DEPTHA*DEPTHD**(NDEPTH-1), DEPTHD
IF (ALLOCATED(TSIHKD)) THEN
   IF (ITEST.GT.0) THEN
      WRITE (IU06,'(''  PRINTED IN STEPS OF '',I3,'' ENTRIES'',/)') K
      DO I = 1,NDEPTH,K
         WRITE (IU06,'('' DEPTH = '',F7.1,'' METRES '')') DEPTHA*DEPTHD**(I-1)
         WRITE (IU06,'('' GROUP VELOCITY IN METRES/SECOND'')')
         WRITE (IU06,'(1X,13F10.5)') TCGOND(I,1:ML)
         WRITE (IU06,'('' WAVE NUMBER IN 1./METRES'')')
         WRITE (IU06,'(1X,13F10.5)') TFAK(I,1:ML)
         WRITE (IU06,'('' OMEGA/SINH(2KD) IN 1./SECOND'')')
         WRITE (IU06,'(1X,13F10.5)') TSIHKD(I,1:ML)
         WRITE (IU06,'('' 2*G*K**2/(OMEGA*TANH(2KD)) IN 1./(METRE*SECOND)'')')
         WRITE (IU06,'(1X,13F10.5)') TFAC_ST(I,1:ML)
      END DO
   END IF
ELSE
   WRITE (IU06,*) '  MODULE DATA ARE NOT PREPARED'
END IF

END SUBROUTINE PRINT_FRE_DIR_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVAT MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MAKE_FRE_DIR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: M, K
REAL    :: CO1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE FREQUENCIES.                                                  !
!        --------------------                                                  !

INV_LOG_CO = 1./LOG10(CO)

DO M = 2,ML
   FR(M) = CO*FR(M-1)
END DO

IF (.NOT.ALLOCATED (FR5 )) ALLOCATE (FR5 (1:ML)) !! FR**5.
IF (.NOT.ALLOCATED (FRM5)) ALLOCATE (FRM5(1:ML)) !! 1. / FR**5.

FR5(:) = FR(:)**5
FRM5(:) = 1/FR5(:)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTE DEEP WATER GROUP VELOCITIES.                                  !
!        ------------------------------------                                  !

IF (.NOT.ALLOCATED (GOM))  ALLOCATE (GOM(1:ML))
GOM = G/(2.*ZPI*FR)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTE PHASE VELOCITY IN DEEP WATER.                                 !
!         -------------------------------------                                !

IF (.NOT.ALLOCATED (C))  ALLOCATE (C(1:ML))
C = G/(ZPI*FR)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. COMPUTATION OF DIRECTION BANDWIDTH.                                   !
!        -----------------------------------                                   !

DELTH = ZPI/REAL(KL)
DELTR = DELTH*R

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. COMPUTATION OF DIRECTIONS.                                            !
!        --------------------------                                            !

IF (.NOT.ALLOCATED (TH))  ALLOCATE (TH(1:KL))
DO K = 1,KL
   TH(K) = REAL(K-1)*DELTH + 0.5*DELTH
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. COMPUTATION SIN AND COS OF DIRECTIONS.                                !
!        --------------------------------------                                !

IF (.NOT.ALLOCATED (COSTH))  ALLOCATE (COSTH(1:KL))
IF (.NOT.ALLOCATED (SINTH))  ALLOCATE (SINTH(1:KL))
COSTH = COS(TH)
SINTH = SIN(TH)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     8. COMPUTATION FREQUENCY DIRECTION INTERVALS AND AREAS.                  !
!        ----------------------------------------------------                  !

IF (.NOT.ALLOCATED (DF)   )  ALLOCATE (DF(1:ML))     !! FREQUENCY INTERVAL.
IF (.NOT.ALLOCATED (DF_FR))  ALLOCATE (DF_FR(1:ML))  !! DF*FR.
IF (.NOT.ALLOCATED (DF_FR2)) ALLOCATE (DF_FR2(1:ML)) !! DF*FR*FR.

CO1 = 0.5*(CO-1.)
DF(1) = CO1*FR(1)
DF(2:ML-1) = CO1 * (FR(2:ML-1)+FR(1:ML-2))
DF(ML) = CO1*FR(ML-1)

DF_FR(:)  = DF(:)*FR(:)
DF_FR2(:) = DF_FR(:)*FR(:)

IF (.NOT.ALLOCATED (DFIM)    )  ALLOCATE (DFIM(1:ML))     !! DF*DELTH.
IF (.NOT.ALLOCATED (DFIMOFR) )  ALLOCATE (DFIMOFR(1:ML))  !! DFIM/FR.
IF (.NOT.ALLOCATED (DFIM_FR) )  ALLOCATE (DFIM_FR(1:ML))  !! DFIM*FR.
IF (.NOT.ALLOCATED (DFIM_FR2))  ALLOCATE (DFIM_FR2(1:ML)) !! DFIM*FR**2.

DFIM(:)     = DF(:)*DELTH                !! MO  INTEGRATION WEIGHTS.
DFIMOFR(:)  = DFIM(:)/FR(:)              !! M-1 INTEGRATION WEIGHTS.
DFIM_FR(:)  = DF_FR(:)*DELTH             !! M+1 INTEGRATION WEIGHTS.
DFIM_FR2(:) = DF_FR2(:)*DELTH            !! M+2 INTEGRATION WEIGHTS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     9. COMPUTATION AIL FACTOR.                                               !
!        ------------------------                                              !

MO_TAIL  = -DELTH/ REAL(EX_TAIL+1)*FR(ML)      !! MO  TAIL FACTOR.
MM1_TAIL = -DELTH/ REAL(EX_TAIL)               !! M-1 TAIL FACTOR.
MP1_TAIL = -DELTH/ REAL(EX_TAIL+2)*FR(ML)**2   !! M+1 TAIL FACTOR.
MP2_TAIL = -DELTH/ REAL(EX_TAIL+3)*FR(ML)**3   !! M+2 TAIL FACTOR.

! ---------------------------------------------------------------------------- !
!                                                                              !
!    10. MINIMUM ENEGRY DENSITY (HS = 7CM).                                    !
!        ----------------------------------                                    !

FMIN = 0.07**2 /(16.*(FR(ML)-FR(1))*ZPI)

END SUBROUTINE MAKE_FRE_DIR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MAKE_SHALLOW_TABLES

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MAKE_SHALLOW_TABLES - ROUTINE TO COMPUTE TABLES USED FOR SHALLOW WATER.    !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE TABLES USED FOR SHALLOW WATER.                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!      TABLES FOR GROUP VELOCITY, WAVE NUMBER AND OMEGA/SINH(2KD) ARE COMPUTED !
!      AT ALL FREQUENCIES AND FOR A DEPTH TABLE OF LENGTH NDEPTH, STARTING AT  !
!      DEPTHA METERS AND INCREMENTED BY DEPTHD METRES.                         !
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

INTEGER        :: M, JD
REAL           :: GH, OM, AD, AK, AKD

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. ALLOCATE TABLE ARRAYS.                                                !
!        ----------------------                                                !

IF (.NOT.ALLOCATED (TFAK))    ALLOCATE (TFAK  (NDEPTH,ML))
IF (.NOT.ALLOCATED (TCGOND))  ALLOCATE (TCGOND(NDEPTH,ML))
IF (.NOT.ALLOCATED (TSIHKD))  ALLOCATE (TSIHKD(NDEPTH,ML))
IF (.NOT.ALLOCATED (TFAC_ST)) ALLOCATE (TFAC_ST(NDEPTH,ML))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. GROUP VELOCITY AND WAVE NUMBER.                                       !
!        -------------------------------                                       !

GH = G/(4.*PI)
DO M = 1,ML                             !! LOOP OVER FREQUENCIES.
   OM=ZPI*FR(M)
   DO JD = 1,NDEPTH                     !! LOOP OVER DEPTH.
      AD = DEPTHA*DEPTHD**(JD-1)
      AK = AKI(OM,AD)
      TFAK(JD,M) = AK
      AKD = AK*AD
      IF (AKD.LE.10.0) THEN
         TCGOND(JD,M) = 0.5*SQRT(G*TANH(AKD)/AK) * (1.0+2.0*AKD/SINH(2.*AKD))
         TSIHKD(JD,M) = OM/SINH(2.*AKD)
         TFAC_ST(JD,M) = 2.*G*AK**2/(OM*TANH(2.*AKD))
      ELSE
         TCGOND(JD,M) = GH/FR(M)
         TSIHKD(JD,M) = 0.
         TFAC_ST(JD,M) = 2./G*OM**3
      END IF
   END DO
END DO

! ---------------------------------------------------------------------------- !

RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTERNAL FUNCTION.                                                    !
!        ------------------                                                    !

CONTAINS

REAL FUNCTION AKI (OM, BETA)

   ! ------------------------------------------------------------------------- !
   !                                                                           !
   !   AKI - FUNCTION TO COMPUTE WAVE NUMBER.                                  !
   !                                                                           !
   !     G. KOMEN, P. JANSSEN   KNMI        01/06/1986                         !
   !                                                                           !
   !     PURPOSE.                                                              !
   !     -------                                                               !
   !                                                                           !
   !       WAVE NUMBER AS FUNCTION OF CIRCULAR FREQUENCY AND WATER DEPTH.      !
   !                                                                           !
   !     METHOD.                                                               !
   !     -------                                                               !
   !                                                                           !
   !       NEWTONS METHOD TO SOLVE THE DISPERSION RELATION IN SHALLOW WATER.   !
   !                                                                           !
   !     REFERENCE.                                                            !
   !     ----------                                                            !
   !                                                                           !
   !       NONE.                                                               !
   !                                                                           !
   ! ------------------------------------------------------------------------- !

   REAL, INTENT(IN) :: OM    !! CIRCULAR FREQUENCY.
   REAL, INTENT(IN) :: BETA  !! WATER DEPTH.

   ! ------------------------------------------------------------------------- !
   !                                                                           !
   !     LOCAL VARIABLES.                                                      !
   !     ----------------                                                      !

   REAL, PARAMETER :: EBS = 0.0001  !! RELATIVE ERROR LIMIT OF NEWTON'S METHOD.

   REAL :: AKP, BO, TH, STH

   ! ------------------------------------------------------------------------- !
   !                                                                           !
   !     1. START WITH MAXIMUM FROM DEEP AND EXTREM SHALLOW WATER WAVE NUMBER. !
   !        ------------------------------------------------------------------ !

   AKI   = MAX ( OM**2/(4.*G), OM/(2.*SQRT(G*BETA)) )

   ! ------------------------------------------------------------------------- !
   !                                                                           !
   !     2. ITERATION LOOP.                                                    !
   !        ---------------                                                    !

   AKP = 10000.
   DO WHILE (ABS(AKP-AKI).GT.EBS*AKI)
      BO = BETA*AKI
      IF (BO.GT.40.) THEN
         AKI = OM**2/G
         EXIT
      ELSE
         AKP = AKI
         TH = G*AKI*TANH(BO)
         STH = SQRT(TH)
         AKI = AKI+(OM-STH)*STH*2./(TH/AKI+G*BO/COSH(BO)**2)
      END IF
   END DO

   END FUNCTION AKI

END SUBROUTINE MAKE_SHALLOW_TABLES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_FRE_DIR_MODULE

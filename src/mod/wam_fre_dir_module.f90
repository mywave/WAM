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
USE WAM_GENERAL_MODULE, ONLY: G, PI, ZPI, DEG, R, ROWATER

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
REAL,ALLOCATABLE :: FRM5(:)      !! 1./FR(M)**5
REAL,ALLOCATABLE :: DFIM(:)      !! FREQUENCY INTERVAL*DIRECTION INTER.
REAL,ALLOCATABLE :: DFIMOFR(:)   !! DFIM/FR
REAL,ALLOCATABLE :: DFFR(:)      !! DFIM*FR
REAL,ALLOCATABLE :: DFFR2(:)     !! DFIM*FR**2 
REAL,ALLOCATABLE :: RHOWG_DFIM(:)!! ROWATER*G*DELTH*LOG(CO)

REAL,ALLOCATABLE, DIMENSION(:) :: GOM   !! DEEP WATER GROUP VELOCITIES [M/S].
REAL,ALLOCATABLE, DIMENSION(:) :: C     !! DEEP WATER PHASE VELOCITIES [M/S].
REAL,ALLOCATABLE, DIMENSION(:) :: TH    !! DIRECTIONS IN RADIANS.
REAL,ALLOCATABLE, DIMENSION(:) :: COSTH !! COS OF DIRECTION.
REAL,ALLOCATABLE, DIMENSION(:) :: SINTH !! SIN OF DIRECTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INDEX OF NEIGHTBOUR TERMS IN PROPAGATION.                             !
!        -----------------------------------------                             !

INTEGER, PUBLIC, ALLOCATABLE :: MPM(:,:)    !! INDEX FOR DIRECTION TERMS.
INTEGER, PUBLIC, ALLOCATABLE :: KPM(:,:)    !! INDEX FOR FREQUENCY TERMS.
INTEGER, PUBLIC, ALLOCATABLE :: JXO(:,:)    !! INDEX FOR EAST-WEST TERMS.
INTEGER, PUBLIC, ALLOCATABLE :: JYO(:,:)    !! INDEX FOR NORTH-SOUTH TERMS.

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

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_FRE_DIR (N_DIR, N_FRE, FR1, IREF)

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
REAL    , INTENT(IN) :: FR1     !! REFERENCE FREQUENCY [HZ].
INTEGER , INTENT(IN) :: IREF    !! FREQUENCY BIN NUMBER OF REFERENCE FREQU.

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
IF (ALLOCATED (RHOWG_DFIM)) DEALLOCATE (RHOWG_DFIM)

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
FR(1) = CO**(-IREF+1)*FR1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. FREQUENCY-DIRECTION GRID.                                             !
!        -------------------------                                             !

CALL MAKE_FRE_DIR
IF (ITEST.GT.1) WRITE (IU06,*) ' SUB. MAKE_FRE_DIR DONE'

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
IF (.NOT.ALLOCATED (RHOWG_DFIM)) ALLOCATE (RHOWG_DFIM(1:ML))

DFIM(:)     = DF(:)*DELTH                !! MO  INTEGRATION WEIGHTS.
DFIMOFR(:)  = DFIM(:)/FR(:)              !! M-1 INTEGRATION WEIGHTS.
DFIM_FR(:)  = DF_FR(:)*DELTH             !! M+1 INTEGRATION WEIGHTS.
DFIM_FR2(:) = DF_FR2(:)*DELTH            !! M+2 INTEGRATION WEIGHTS.
RHOWG_DFIM(:)=ROWATER*G*DELTH*LOG(CO)*FR(:) !! MOMENTUM AND ENERGY FLUX WEIGHTS.
RHOWG_DFIM(1)=0.5*RHOWG_DFIM(1)             !! TRAPEZOIDAL INTEGRATION
RHOWG_DFIM(ML)=0.5*RHOWG_DFIM(ML)           !! WITH CHANGE OF VARIABLE x=LOG(f)
                                            !! HENCE df = f dx,
                                            !! dx=LOG(FR(n+1))-LOG(FR(n)) = LOG(CO)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     9. COMPUTATION TAIL FACTOR.                                              !
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

! ---------------------------------------------------------------------------- !
!                                                                              !
!     11. INDEX OF NEIGHTBOUR TERMS.                                           !
!        ---------------------------                                           !

IF (.NOT. ALLOCATED(MPM)) ALLOCATE(MPM(ML,-1:1))  !! FREQUENCY NEIGHTBOURS.
DO M=1,ML
   MPM(M,-1)= MAX(1,M-1)
   MPM(M,0) = M
   MPM(M,1) = MIN(ML,M+1)
ENDDO

IF(.NOT. ALLOCATED(KPM)) ALLOCATE(KPM(KL,-1:1))  !! DIRECTION NEIGHTBOURS.
IF(.NOT. ALLOCATED(JXO)) ALLOCATE(JXO(KL,2))
IF(.NOT. ALLOCATED(JYO)) ALLOCATE(JYO(KL,2))

DO K=1,KL

   KPM(K,-1) = K-1
   IF (KPM(K,-1).LT.1) KPM(K,-1) = KL
   KPM(K,0) = K
   KPM(K,1) = K+1
   IF (KPM(K,1).GT.KL) KPM(K,1) = 1

   IF (COSTH(K).GE.0.) THEN
      JYO(K,1)=1
      JYO(K,2)=2
   ELSE
      JYO(K,1)=2
      JYO(K,2)=1
   ENDIF

   IF (SINTH(K).GE.0.) THEN
      JXO(K,1)=1
      JXO(K,2)=2
   ELSE
      JXO(K,1)=2
      JXO(K,2)=1
   ENDIF
ENDDO

END SUBROUTINE MAKE_FRE_DIR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_FRE_DIR_MODULE

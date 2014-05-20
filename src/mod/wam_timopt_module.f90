MODULE WAM_TIMOPT_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL MODEL TIMES, TIMESTEPS AND OPTIONS                !
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

USE WAM_FILE_MODULE,     ONLY: IU06

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. TIME STEPS.                                                           !
!        -----------                                                           !

INTEGER  :: IDELPRO     = -1 !! TIMESTEP WAM PROPAGATION IN SECONDS.
INTEGER  :: IDELT       = -1 !! TIMESTEP SOURCE FUNCTION IN SECONDS.
INTEGER  :: IDEL_WAM    = -1 !! TIMESTEP WAMODEL CALLS IN SECONDS.
integer  :: ifcst       =  0 !! ready file time in days and hours

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. OPTIONS.                                                              !
!        ---------                                                             !

LOGICAL :: COLDSTART         !! TRUE IF COLDSTART
LOGICAL :: SPHERICAL_RUN     !! TRUE: SPHERICAL  PROPAGATION,  
                             !! FALSE: CARTESIAN PROPAGATION.
LOGICAL :: SHALLOW_RUN       !! TRUE:  SHALLOW WATER MODEL, 
                             !! FALSE:  DEEP WATER MODEL. 
LOGICAL :: REFRACTION_D_RUN  !! TRUE: DEPTH REFRACTION ON.                                            
LOGICAL :: REFRACTION_C_RUN  !! TRUE: CURRENT REFRACTION ON.                                            
LOGICAL :: WAVE_BREAKING_RUN !! TRUE: WAVE BREAKING ON.                                          
LOGICAL :: TOPO_RUN          !! TRUE: INSTATIONARY WATER DEPTH.
LOGICAL :: CURRENT_RUN       !! TRUE: INSTATIONARY CURRENTS.
LOGICAL :: PHILLIPS_RUN      !! TRUE: PHILLIPS SOURCE TERM 0N.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. TIME STATUS OF INTEGRATION.                                           !
!        ---------------------------                                           !

CHARACTER (LEN=14) :: CDATEA      =' ' !! START DATE OF RUN  (VVYYMMDDHHMMSS).
CHARACTER (LEN=14) :: CDATEE      =' ' !! END DATE OF RUN (VVYYMMDDHHMMSS).
CHARACTER (LEN=14) :: CDTPRO      =' ' !! END DATE OF PROPAGATION.
CHARACTER (LEN=14) :: CDTSOU      =' ' !! END DATE OF SOURCE INTEGRATION.
CHARACTER (LEN=14) :: CDA         =' ' !! DATE OF WINDFIELD USED IN WAM.
CHARACTER (LEN=14) :: CDATEWO     =' ' !! DATE OF NEXT WIND FIELD TO BE USED.
CHARACTER (LEN=14) :: CDTA        =' ' !! DATE OF DEPTH FIELD USED IN WAM.
CHARACTER (LEN=14) :: CD_TOPO_NEW =' ' !! DATE OF NEXT DEPTH FIELD TO BE USED.
CHARACTER (LEN=14) :: CDCA        =' ' !! DATE OF CURRENT FIELD USED IN WAM.
CHARACTER (LEN=14) :: CD_CURR_NEW =' ' !! DATE OF NEXT CURRENT FIELD TO BE USED.
character (len=14) :: cdtstop     =' ' !! full 2-d-spectral output up to cdtstop

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE SET_MODEL_OPTION             !! SETS MODEL OPTIONS.
   MODULE PROCEDURE SET_MODEL_OPTION
END INTERFACE
PUBLIC SET_MODEL_OPTION

INTERFACE SET_INTEGRATION_TIMESTEPS    !! SETS MODEL TIMESTEPS.
   MODULE PROCEDURE SET_INTEGRATION_TIMESTEPS
END INTERFACE
PUBLIC SET_INTEGRATION_TIMESTEPS

INTERFACE SET_INTEGRATION_PERIOD      !! SETS MODEL TIMESTEPS.
   MODULE PROCEDURE SET_INTEGRATION_PERIOD
END INTERFACE
PUBLIC SET_INTEGRATION_PERIOD

INTERFACE PRINT_TIMOPT_STATUS          !! PRINTS MODULE STATUS.
   MODULE PROCEDURE PRINT_TIMOPT_STATUS
END INTERFACE
PUBLIC PRINT_TIMOPT_STATUS

INTERFACE SET_START_OPTION             !! SETS START OPTION.
   MODULE PROCEDURE SET_START_OPTION
END INTERFACE
PUBLIC SET_START_OPTION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_MODEL_OPTION (SPHERICAL, SHALLOW, REFRACTION_D, REFRACTION_C,   &
&                            WAVE_BREAKING, PHILLIPS)

LOGICAL, OPTIONAL, INTENT(IN) :: SPHERICAL      !! PROPAGATION FLAG
                                                !! TRUE:  SPHERICAL COORDINATES
                                                !! FALSE: CARTESIAN COORDINATES.
LOGICAL, OPTIONAL, INTENT(IN) :: SHALLOW        !! SHALLOW WATER MODEL FLAG
                                                !! TRUE:  SHALLOW WATER MODEL.
                                                !! FALSE: DEEP WATER MODEL
LOGICAL, OPTIONAL, INTENT(IN) :: REFRACTION_D   !! DEPTH REFRACTION OPTION.
                                                !! FALSE: NO DEPTH REFRACTION.
                                                !! TRUE:  DEPTH REFRACTION.
LOGICAL, OPTIONAL, INTENT(IN) :: REFRACTION_C   !! CURRENT REFRACTION OPTION.
                                                !! TRUE:  CURRENT REFRACTION.
                                                !! FALSE: NO CURRENT REFRACTION.
LOGICAL, OPTIONAL, INTENT(IN) :: WAVE_BREAKING  !! WAVE BREAKING OPTION.
                                                !! FALSE:  NO BREAKING.
                                                !! TRUE:   BREAKING ACTIVE.
LOGICAL, OPTIONAL, INTENT(IN) :: PHILLIPS       !! PHILLIPS SOURCE.
                                                !! FALSE:  OFF.
                                                !! TRUE:   ON.
                                    
SPHERICAL_RUN     = .TRUE.
SHALLOW_RUN       = .TRUE.
REFRACTION_D_RUN  = .FALSE.
REFRACTION_C_RUN  = .FALSE.
WAVE_BREAKING_RUN = .FALSE.
PHILLIPS_RUN      = .FALSE.

IF (PRESENT(SPHERICAL    )) SPHERICAL_RUN     = SPHERICAL
IF (PRESENT(SHALLOW      )) SHALLOW_RUN       = SHALLOW
IF (PRESENT(REFRACTION_D )) REFRACTION_D_RUN  = REFRACTION_D
IF (PRESENT(REFRACTION_C )) REFRACTION_C_RUN  = REFRACTION_C
IF (PRESENT(WAVE_BREAKING)) WAVE_BREAKING_RUN = WAVE_BREAKING
IF (PRESENT(PHILLIPS     )) PHILLIPS_RUN = PHILLIPS

IF (SHALLOW_RUN .AND.REFRACTION_C_RUN .AND. .NOT.REFRACTION_D_RUN) THEN
   WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' +        WARNING ERROR SUB.SET_MODEL_OPTION.        +'
   WRITE (IU06,*) ' +        ===================================        +'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' + A SHALLOW WATER MODEL RUN WITH CURRENT REFRACTION +'
   WRITE (IU06,*) ' + AND WITHOUT DEPTH REFRACTION IS REQUESTED.        +'
   WRITE (IU06,*) ' + THIS IS NOT IMPLEMENTED.                          +'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' +               MODEL CONTINUES                     +'
   WRITE (IU06,*) ' +       WITH DEPTH AND CURRENT REFRACTION           +'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++++++'
   REFRACTION_D_RUN = .TRUE.
END IF

IF (.NOT.SHALLOW_RUN .AND. REFRACTION_D_RUN) THEN
   WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' +        WARNING ERROR SUB.SET_MODEL_OPTION.        +'
   WRITE (IU06,*) ' +        ===================================        +'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' + A DEEP WATER MODEL RUN WITH DEPTH REFRACTION      +'
   WRITE (IU06,*) ' + IS REQUESTED.                                     +'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' +               MODEL CONTINUES                     +'
   WRITE (IU06,*) ' +           WITHOUT DEPTH REFRACTION                +'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++++++'
   REFRACTION_D_RUN = .FALSE.
END IF
IF (.NOT.SHALLOW_RUN .AND. WAVE_BREAKING_RUN) THEN
   WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' +        WARNING ERROR SUB.SET_MODEL_OPTION.        +'
   WRITE (IU06,*) ' +        ===================================        +'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' + A DEEP WATER MODEL RUN WITH WAVE BREAKING         +'
   WRITE (IU06,*) ' + IS REQUESTED.                                     +'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' +               MODEL CONTINUES                     +'
   WRITE (IU06,*) ' +             WITHOUT WAVE BREAKING                 +'
   WRITE (IU06,*) ' +                                                   +'
   WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WAVE_BREAKING_RUN = .FALSE.
END IF

END SUBROUTINE SET_MODEL_OPTION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_INTEGRATION_TIMESTEPS (P, S)

INTEGER, INTENT(IN)   :: P  !! PROPAGATION TIMSTEP [S].
INTEGER, INTENT(IN)   :: S  !! SOURCE FUNCTION TIMESTEP [S].

IDELPRO = P
IDELT   = S
IF (S.LE.0) IDELT = IDELPRO

END SUBROUTINE SET_INTEGRATION_TIMESTEPS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_INTEGRATION_PERIOD (B, E)

CHARACTER (LEN=14), INTENT(IN)   :: B   !! START DATE OF MODEL RUN.
CHARACTER (LEN=14), INTENT(IN)   :: E   !! END DATE OF MODEL RUN.

CDATEA = B
CDATEE = E

END SUBROUTINE SET_INTEGRATION_PERIOD

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_TIMOPT_STATUS

WRITE(IU06,*) '  '
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '              MODEL TIMOPT STATUS:'
WRITE(IU06,*) ' ------------------------------------------------- '

WRITE(IU06,*) '  '
WRITE(IU06,*) ' START DATE (YYYYMMDDHHMMSS) IS.....: ', CDATEA
WRITE(IU06,*) ' END   DATE (YYYYMMDDHHMMSS) IS.....: ', CDATEE
WRITE(IU06,*) ' DATE OF PROPAGATION INTEGRATION IS.: ', CDTPRO
WRITE(IU06,*) ' DATE OF SOURCE INTEGRATION IS......: ', CDTSOU
WRITE(IU06,*) ' DATE OF WINDFIELD IS...............: ', CDA
WRITE(IU06,*) ' DATE FOR NEXT WIND FIELD IS........: ', CDATEWO
WRITE(IU06,*) ' DATE OF DEPTH FIELD IS.............: ', CDTA
WRITE(IU06,*) ' DATE FOR NEXT DEPTH FIELD IS.......: ', CD_TOPO_NEW
WRITE(IU06,*) ' DATE OF CURRENT FIELD IS...........: ', CDCA
WRITE(IU06,*) ' DATE FOR NEXT CURRENT FIELD IS.....: ', CD_CURR_NEW

WRITE(IU06,*) '  '
WRITE(IU06,*) ' MODEL TIME STEPS:'
WRITE(IU06,*) ' SOURCE TERM INTEGRATION TIME STEP..: ',IDELT,' SECS'
WRITE(IU06,*) ' PROPAGATION TIME STEP .............: ',IDELPRO,' SECS'
WRITE(IU06,*) ' WAMODEL TIME STEP .................: ',IDEL_WAM,' SECS'
WRITE(IU06,*) '  '
WRITE(IU06,*) ' MODEL OPTIONS:'
IF (SPHERICAL_RUN) THEN
   WRITE(IU06,*) ' PROPAGATION GRID SPHERICAL LAT/LON COORDINATES'
ELSE
   WRITE(IU06,*) ' PROPAGATION GRID CARTESIAN COORDINATES'
END IF
IF (SHALLOW_RUN) THEN
   WRITE(IU06,*) ' THIS IS A SHALLOW WATER RUN '
ELSE
   WRITE(IU06,*) ' THIS IS A DEEP WATER RUN '
END IF
IF (REFRACTION_D_RUN) THEN
   WRITE(IU06,*) ' MODEL RUNS WITH DEPTH REFRACTION'
ELSE
   WRITE(IU06,*) ' MODEL RUNS WITHOUT DEPTH REFRACTION'
END IF
IF (REFRACTION_C_RUN) THEN
   WRITE(IU06,*) ' MODEL RUNS WITH CURRENT REFRACTION'
ELSE
   WRITE(IU06,*) ' MODEL RUNS WITHOUT CURRENT REFRACTION'
END IF

IF (WAVE_BREAKING_RUN) THEN
   WRITE(IU06,*) ' MODEL RUNS WITH  WAVE BREAKING'
ELSE 
   WRITE(IU06,*) ' MODEL RUNS WITHOUT WAVE BREAKING'
END IF

IF (PHILLIPS_RUN) THEN
   WRITE(IU06,*) ' MODEL RUNS WITH PHILLIPS SOURCE TERM'
ELSE 
   WRITE(IU06,*) ' MODEL RUNS WITHOUT PHILLIPS SOURCE TERM'
END IF

IF (COLDSTART) THEN
   WRITE (IU06,*) ' INITIAL VALUES ARE PROCESSED DUE TO OPTION (COLD START).'
ELSE 
   WRITE (IU06,*) ' INITIAL VALUES FROM A PREVIOUS MODEL RUN (HOT START).'
END IF
WRITE(IU06,*) '  '

END SUBROUTINE PRINT_TIMOPT_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_START_OPTION (OPTION)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

LOGICAL, INTENT(IN)  :: OPTION   !! MODEL START OPTION.

! ---------------------------------------------------------------------------- !

COLDSTART = OPTION

END SUBROUTINE SET_START_OPTION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_TIMOPT_MODULE

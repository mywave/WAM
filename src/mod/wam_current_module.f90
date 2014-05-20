MODULE WAM_CURRENT_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE STORES THE CURRENT INPUT GRID SPECFICATIONS AND                !
!   THE CURRENT FIELDS, WHICH ARE PASSED TO THE WAM-MODELL.                    !
!   IT CONTAINS ALL PROCEDURES NESSECARY FOR THE PROCESSING.                   !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE          !! COORDINATE PROCEDURES

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1,                  & !! TERMINATES PROCESSING.
&       INCDATE                    !! INCREMENTS DATE TIME GROUP.

USE WAM_GRID_MODULE, ONLY:       &
&       INTERPOLATION_TO_GRID,   & !! INTERPOLATE TO WAM POINTS.
&       EQUAL_TO_M_GRID            !! COMPARES CURRENT GRID TO MODEL GRID.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST
USE WAM_GENERAL_MODULE, ONLY: RAD
USE WAM_GRID_MODULE,    ONLY: NSEA, L_S_MASK
USE WAM_MODEL_MODULE,   ONLY: U, V
USE WAM_TIMOPT_MODULE,  ONLY: CURRENT_RUN, CDCA, CDATEE, IDEL_WAM

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INPUT CURRENT GRID SPECFICATIONS, DATE AND CURRENT FIELDS.            !
!        ----------------------------------------------------------            !

INTEGER :: NX_IN = -1    !! NUMBER OF LONGITUDES.
INTEGER :: NY_IN = -1    !! NUMBER OF LATITUDES.
LOGICAL :: PER = .FALSE. !! .TRUE. IF GRID IS PERIODICAL.
INTEGER :: CODE_IN = 0   !! CURRENT CODE: 
                         !! 1 = SPEED AND DIRECTION, 
                         !! OTHERWISE: COMPONENTS
INTEGER :: DX_IN         !! STEPSIZE BETWEEN LONGITUDES [M_SEC].
INTEGER :: DY_IN         !! STEPSIZE BETWEEN LATITUDES [M_SEC].
INTEGER :: SOUTH_IN      !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER :: NORTH_IN      !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER :: WEST_IN       !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER :: EAST_IN       !! EAST LONGITUDE OF GRID [M_SEC].
LOGICAL :: EQUAL_GRID =.FALSE. !! .TRUE. IF CURRENT GRID IS EQUAL TO MODEL GRID.

CHARACTER (LEN=14) :: CD_READ = ' ' !! DATE OF LAST DATA READ FROM INPUT.
REAL, ALLOCATABLE, DIMENSION(:,:) :: U_IN  !! U-COMP. OF CURRENTS [M/S].
REAL, ALLOCATABLE, DIMENSION(:,:) :: V_IN  !! V-COMP. OF CURRENTS [M/S].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CURRENT TIMESTEPS.                                                    !
!        ------------------                                                    !

INTEGER, PUBLIC   :: IDELCI = -1  !! INPUT CURRENT TIMESTEP IN SECONDS.
INTEGER, PUBLIC   :: IDELCO = -1  !! OUTPUT CURRENT TIMESTEP IN SECONDS
                                  !! EQUAL TO INPUT TIMESTEP INTO WAMODEL.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CURRENT FIELDS PREPARED FOR WAM-MODEL.                                !
!        (FIRST INDEX IS POINTS, SECOND IS TIME)                               !
!        ----------------------------------------------                        !

INTEGER  :: M_STORE = 0    !! NUMBER OF CURRENT FIELDS STORED.

REAL,               ALLOCATABLE, DIMENSION(:,:) :: U_STORE   !! U-COMP.
REAL,               ALLOCATABLE, DIMENSION(:,:) :: V_STORE   !! V-COMP.
CHARACTER (LEN=14), ALLOCATABLE, DIMENSION(:)   :: CD_STORE  !! DATE/TIME.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE GET_CURRENT                   !! GETS CURRENTS FORM THIS MODULE.
   MODULE PROCEDURE GET_CURRENT
END INTERFACE
PUBLIC GET_CURRENT

INTERFACE PREPARE_CURRENT               !! PREPARES CURRENT DATA FOR WAM MODEL.
   MODULE PROCEDURE PREPARE_CURRENT
END INTERFACE
PUBLIC PREPARE_CURRENT

INTERFACE PRINT_CURRENT_STATUS          !! PRINTS CURRENT STATUS.
   MODULE PROCEDURE PRINT_CURRENT_STATUS
END INTERFACE
PUBLIC PRINT_CURRENT_STATUS

INTERFACE SET_CURRENT_FIELD              !! SETS A CURRENT FIELD.
   MODULE PROCEDURE SET_CURRENT_FIELD
END INTERFACE
PUBLIC SET_CURRENT_FIELD

INTERFACE SET_CURRENT_HEADER             !! SETS CURRENT HEADER.
   MODULE PROCEDURE SET_CURRENT_HEADER_C !! CHARACTER VERSION
   MODULE PROCEDURE SET_CURRENT_HEADER_D !! DEGREE VERSION
   MODULE PROCEDURE SET_CURRENT_HEADER_M !! M_SEC VERSION
END INTERFACE
PUBLIC SET_CURRENT_HEADER

INTERFACE SET_CURRENT_TIMESTEPS          !! SETS CURRENT TIMESTEPS.
   MODULE PROCEDURE SET_CURRENT_TIMESTEPS
END INTERFACE
PUBLIC SET_CURRENT_TIMESTEPS

INTERFACE WAM_CURRENT     !! READS AND TRANSFORMS INPUT CURRENTS TO WAM POINTS.
   MODULE PROCEDURE WAM_CURRENT
END INTERFACE
PUBLIC WAM_CURRENT

INTERFACE
   SUBROUTINE READ_CURRENT_INPUT       !! READS A CURRENT FIELD
   END SUBROUTINE READ_CURRENT_INPUT
END INTERFACE 

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE NOTIM           !! STEERING SUB IF TIME INTERPOLATION IS NOT WANTED.
   MODULE PROCEDURE NOTIM
END INTERFACE
PRIVATE NOTIM

INTERFACE TIMIN           !! STEERING SUB IF TIME INTERPOLATION IS WANTED.
   MODULE PROCEDURE TIMIN
END INTERFACE
PRIVATE TIMIN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE GET_CURRENT (CD_NEW)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   GET_CURRENT - TRANSFERS NEW CURRENT DATA TO WAM MODEL.                     !
!                                                                              !
!     H. GUENTHER  GKSS    DECEMBER 2009                                       !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER CURRENTS FROM WAM_CURRENT_MODULE TO THE WAM_MODEL_MODULE.     !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       COPY NEW CURRENTS.                                                     !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=14),INTENT(INOUT) :: CD_NEW !! DATE OF A NEW CURRENT FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER   :: IT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NEW CURRENTS FROM TRANSFER ARRAYS.                                    !
!        ----------------------------------                                    !

CDCA  = ' '
DO IT = 1, M_STORE
   IF ( CD_STORE(IT).GT.CD_NEW) THEN
      U =  U_STORE(:,IT)
      V =  V_STORE(:,IT)
      CDCA = CD_STORE(IT)
      EXIT
   END IF
END DO

IF (CDCA .EQ. ' ') THEN
   WRITE(IU06,*) ' *******************************************'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *   FATAL ERROR IN SUB. GET_CURRENT       *'
   WRITE(IU06,*) ' *   ===============================       *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' * CURRENT FIELD WAS NOT FOUND             *'
   WRITE(IU06,*) ' *  CHANGE DATE IS ............. CD_NEW:', CD_NEW
   WRITE(IU06,*) ' *  LAST DATE STORED IN MODULE IS .... :', CD_STORE(M_STORE)
   WRITE(IU06,*) ' *  NO. OF DATES IN MODULE IS M_STORE = ', M_STORE
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *******************************************'
   CALL ABORT1
END IF

CALL INCDATE(CD_NEW,IDELCO)

IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '  '
   WRITE (IU06,*) '   SUB. GET_CURRENT: NEW CURRENT FIELDS PASSED TO WAM MODEL'
   WRITE (IU06,*) '     DATE OF CURRENT FIELD IS ......... CDCA = ', CDCA
   WRITE (IU06,*) '     DATE FOR NEXT CURRENT CHANGE IS  CD_NEW = ', CD_NEW
   WRITE (IU06,*) '  '
END IF

END SUBROUTINE GET_CURRENT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_CURRENT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PREPARE_CURRENT - PREPARES CURRENT DATA FOR WAVE MODEL.                    !
!                                                                              !
!     H. GUNTHER      GKSS  DECEMBER 2009                                      !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       EVALUATE CURRENT AT WAVE MODEL GRID POINTS.                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INPUT CURRENT FIELDSCAN BE COMPONENTS OR MAGNITUDE AND DIRECTION.      !
!       THE INPUT FIELDS HAVE TO BE ON A LAT /LONG GRID.                       !
!       SEE SUB READ_CURRENT FOR FORMATS AND HEADER INFORMATION, WHICH HAVE    !
!       TO BE GIVEN TO THE PROGRAM.                                            !
!                                                                              !
!       A DOUBLE LINEAR INTERPOLATION IN SPACE IS PERFORMED ONTO THE MODEL     !
!       GRID.                                                                  !
!       IF THE CURRENT OUTPUT TIMSTEP IS LESS THAN THE INPUT TIMESTEP          !
!       A LINEAR INTERPOLATION IN TIME IS PERFORMED.                           !
!       ALL CURRENT FIELDS ARE STORED IN WAM_CURRENT_MODULE.                   !
!                                                                              !
!     REFERENCE.                                                               !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

CHARACTER (LEN=14)            :: CD_START, CD_END 

! ---------------------------------------------------------------------------- !
!                                                                              !
!                                                                              !
!     1. BEGIN AND END DATES OF CURRENT FIELDS TO BE PROCESSED.                !
!        ------------------------------------------------------                !

CD_START = CDCA
CALL INCDATE (CD_START,IDELCO)
CD_END = CDCA
CALL INCDATE (CD_END, IDEL_WAM)
IF (CD_END.GE.CDATEE) CD_END = CDATEE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. NUMBER OF CURRENT FIELDS TO BE GENERATED.                             !
!        -----------------------------------------                             !

 M_STORE = IDEL_WAM/IDELCO

IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '  '
   WRITE (IU06,*) '   SUB. PREPARE_CURRENT: CURRENT REQUEST'
   WRITE (IU06,*) '     START OF PERIOD IS ............ CD_START = ', CD_START
   WRITE (IU06,*) '     END   OF PERIOD IS .............. CD_END = ', CD_END
   WRITE (IU06,*) '     CURRENT INPUT TIME STEP ......... IDELCI = ', IDELCI
   WRITE (IU06,*) '     CURRENT OUTPUT TIME STEP ........ IDELCO = ', IDELCO
   WRITE (IU06,*) '     NUMBER FIELDS TO BE GENERATED IS M_STORE = ', M_STORE
   WRITE (IU06,*) '     FIELDS ARE SAVED IN WAM_CURRENT_MODULE'
END IF

IF ( ALLOCATED(CD_STORE)) then
   if (M_STORE.NE.SIZE(CD_STORE)) DEALLOCATE (CD_STORE, U_STORE,  V_STORE)
endif
IF (.NOT. ALLOCATED(U_STORE) ) ALLOCATE (U_STORE(1:NSEA, 1:M_STORE))
IF (.NOT. ALLOCATED(V_STORE) ) ALLOCATE (V_STORE(1:NSEA, 1:M_STORE))
IF (.NOT. ALLOCATED(CD_STORE)) ALLOCATE (CD_STORE(1:M_STORE))
CD_STORE  = ' '

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. PROCESS CURRENT FIELDS.                                               !
!        -----------------------                                               !

IF (IDELCO.GE.IDELCI) THEN

   IF (ITEST.GE.2) WRITE (IU06,*) '     NO TIME INTERPOLATION'
   CALL NOTIM (CD_START, CD_END)

ELSE

   IF (ITEST.GE.2) WRITE (IU06,*) '     TIME INTERPOLATION'
   CALL TIMIN (CD_START, CD_END)

END IF

END SUBROUTINE PREPARE_CURRENT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_CURRENT_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PRINT_CURRENT_STATUS - PRINTS STATUS OF THE CURRENT MODULE.                !
!                                                                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2009                           !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       MAKE A PRINTER OUTPUT OF THE STATUS OF THE CURRENT MODULE.             !
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

character (len=len_coor) :: formtext
    
WRITE (IU06,*) '  '
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '              CURRENT MODULE STATUS:'
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '  '

IF (.NOT. CURRENT_RUN) THEN
   WRITE (IU06,*) ' CURRENTS ARE STATIONARY. '
   RETURN
END IF

WRITE (IU06,*) ' CURRENT INPUT TIME STEP ....................: ', IDELCI,' SECS'
WRITE (IU06,*) ' CURRENT OUTPUT TIME STEP....................: ', IDELCO,' SECS'
IF (IDELCO.GT.0 .AND. IDELCI.GT.0) THEN
   IF (IDELCO.GE.IDELCI) THEN
      WRITE(IU06,*) ' CURRENT FIELDS ARE NOT INTERPOLATED IN TIME'
   ELSE
      WRITE(IU06,*) ' CURRENT FIELDS ARE INTERPOLATED IN TIME'
   END IF
END IF
WRITE (IU06,*) '  '
WRITE (IU06,*) ' DATE OF LAST FIELD GIVEN TO MODULE..........: ', CD_READ
WRITE (IU06,*) ' DATE OF LAST FIELD GIVEN TO WAVE MODEL......: ', CDCA
WRITE (IU06,*) '  '

IF (NX_IN.GT.0 .AND. NY_IN.GT.0) THEN
   WRITE (IU06,'('' CURRENT INPUT GRID SPECIFICATION ARE: '')')
   WRITE (IU06,*) ' NUMBER OF LONGITUDE ..........: ', NX_IN
   formtext = write_coor_text (west_in)
   WRITE (IU06,*) ' WESTERN MOST LONGITUDE........: ', formtext
   formtext = write_coor_text (east_in)
   WRITE (IU06,*) ' EASTERN MOST LONGITUDE........: ', formtext
   formtext = write_coor_text (dx_in)
   WRITE (IU06,*) ' LONGITUDE INCREMENT ..........: ', formtext
   WRITE (IU06,*) ' NUMBER OF LATITUDE ...........: ', NY_IN
   formtext = write_coor_text (south_in)
   WRITE (IU06,*) ' SOUTHERN MOST LATITUDE .......: ', formtext
   formtext = write_coor_text (north_in)
   WRITE (IU06,*) ' NORTHERN MOST LATITUDE .......: ', formtext
   formtext = write_coor_text (dy_in)
   WRITE (IU06,*) ' LATITUDE INCREMENT ...........: ', formtext
   IF (PER) THEN
      WRITE (IU06,*) 'THE GRID IS EAST-WEST PERIODIC'
   ELSE
      WRITE (IU06,*) 'THE GRID IS NOT EAST-WEST PERIODIC'
   END IF
   IF (EQUAL_GRID) THEN
      WRITE (IU06,*) 'GRID IS IDENTICAL TO MODEL GRID: NO SPACE INTERPOLATION'
   ELSE
      WRITE (IU06,*) 'GRID IS NOT IDENTICAL TO MODEL GRID: SPACE INTERPOLATION'
   END IF
   IF (CODE_IN.EQ.1) THEN
      WRITE (IU06,*) ' CURRENT INPUTS ARE SPEED AND DIRECTION'
   ELSE
      WRITE (IU06,*) ' CURRENT INPUTS ARE COMPONENTS'
   END IF
   
   WRITE (IU06,*) '  '

ELSE
   WRITE (IU06,*) ' CURRENT INPUT GRID IS NOT DEFINED'
END IF
WRITE (IU06,*) '  '

IF ( M_STORE.GT.0 .AND. ANY(CD_STORE.NE. ' ')) THEN
   WRITE (IU06,*) ' NUMBER OF CURRENT FIELDS STORED IN MODULE ..: ', M_STORE
   WRITE (IU06,*) ' DATES OF CURRENT FIELDS ARE:'
   WRITE (IU06,'(5(3X,A14,2X))')  CD_STORE
ELSE
   WRITE (IU06,*) ' CURRENT FIELDS ARE NOT STORED IN MODULE'
END IF
WRITE (IU06,*) '  '

END SUBROUTINE PRINT_CURRENT_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_CURRENT_FIELD (CDT, U_MAP, V_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_CURRENT_FIELD - SET INPUT CURRENT FIELD IN WAM_CURRENT MODULE.         !
!                                                                              !
!     H. GUENTHER  GKSS  DECEMBER 2009                                         !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER CURRENTS TO WAM_CURRENT_MODULE.                               !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY NEW CURRENTS. IF CURRENTS ARE GIVEN AS      !
!       MAGNITUDE AND DIRECTION, THEY ARE CHANGED TO VECTOR COMPONENTS.        !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=14),INTENT(IN) :: CDT        !! DATE/TIME OF CURRENT FIELD.
REAL,              INTENT(IN) :: U_MAP(:,:) !! U-COMP. OR SPEED [M/S].
REAL,              INTENT(IN) :: V_MAP(:,:) !! V-COMP. [M/S] OR DIRECTION [DEG]

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COPY DATE.                                                            !
!        ----------                                                            !

CD_READ = CDT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK CONSISTENCY.                                                    !
!        ------------------                                                    !

IF (NX_IN.NE.SIZE(U_MAP,1) .OR. NY_IN.NE.SIZE(U_MAP,2) .OR.                    &
&   NX_IN.NE.SIZE(V_MAP,1) .OR. NY_IN.NE.SIZE(V_MAP,2) ) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *         FATAL  ERROR IN SUB. SET_CURRENT_FIELD         *'
   WRITE (IU06,*) ' *         ======================================         *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' * GRID SPECIFICATIONS ARE NOT CONSISTENT                 *'
   WRITE (IU06,*) ' * OR CURRENT HEADER IS NOT DEFINED.                      *'
   WRITE (IU06,*) ' * (SET_CURRENT_HEADER HAS TO BE PROCESSED BEFORE)        *'
   WRITE (IU06,*) ' * GRID SIZES AS DEFINED IN MODULE ARE:                   *'
   WRITE (IU06,*) ' * NO. OF LONGITUDES    NX_IN = ', NX_IN
   WRITE (IU06,*) ' * NO. OF LATITUDE      NY_IN = ', NY_IN
   WRITE (IU06,*) ' * DIMENSIONS OF U CURRENT ARRAY ARE :                    *'
   WRITE (IU06,*) ' * NO. OF LONGITUDES 1. DIMENSION = ', SIZE(U_MAP,1)
   WRITE (IU06,*) ' * NO. OF LATITUDE   2. DIMENSION = ', SIZE(U_MAP,2)
   WRITE (IU06,*) ' * DIMENSIONS OF V CURRENT ARRAY ARE :                    *'
   WRITE (IU06,*) ' * NO. OF LONGITUDES 1. DIMENSION = ', SIZE(V_MAP,1)
   WRITE (IU06,*) ' * NO. OF LATITUDE   2. DIMENSION = ', SIZE(V_MAP,2)
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COPY CURRENT FIELDS.                                                  !
!        --------------------                                                  !

IF (ALLOCATED(U_IN)) DEALLOCATE(U_IN)
IF (ALLOCATED(V_IN)) DEALLOCATE(V_IN)
ALLOCATE (U_IN(NX_IN,NY_IN))
ALLOCATE (V_IN(NX_IN,NY_IN))

IF (CODE_IN.EQ.1) THEN
   U_IN = U_MAP*SIN(V_MAP*RAD)
   V_IN = U_MAP*COS(V_MAP*RAD)
ELSE
   U_IN = U_MAP
   V_IN = V_MAP
END IF

END SUBROUTINE SET_CURRENT_FIELD

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_CURRENT_HEADER_C (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,       &
&                                N_LON, N_LAT, CODE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_CURRENT_HEADER - SET INPUT GRID FOR CURRENT FIELDS.                    !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2011                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS TO WAM_CURRENT_MODULE.             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY DEFINITIONS.                                !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=LEN_COOR), INTENT(IN)           :: WEST  !! WEST LONGITUDE.
CHARACTER (LEN=LEN_COOR), INTENT(IN)           :: SOUTH !! SOUTH LATITUDE.
CHARACTER (LEN=LEN_COOR), INTENT(IN), OPTIONAL :: EAST  !! EAST LONGITUDE.
CHARACTER (LEN=LEN_COOR), INTENT(IN), OPTIONAL :: NORTH !! NORTH LATITUDE.
CHARACTER (LEN=LEN_COOR), INTENT(IN), OPTIONAL :: D_LON !! LONGITUDE INCREMENT.
CHARACTER (LEN=LEN_COOR), INTENT(IN), OPTIONAL :: D_LAT !! LATITUDE INCREMENT.
INTEGER,                  INTENT(IN), OPTIONAL :: N_LON !! NUMBER OF LONGITUDES.
INTEGER,                  INTENT(IN), OPTIONAL :: N_LAT !! NUMBER OF LATITUDES.
INTEGER, INTENT(IN), OPTIONAL :: CODE    !! CURRENT CODE:
                                         !! 1 = SPEED AND DIRECTION.
                                         !! OTHERWISE: COMPONENTS (DEFAULT).

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: WEST_CO    !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER :: SOUTH_CO   !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER :: EAST_CO    !! EAST LONGITUDE OF GRID [M_SEC].
INTEGER :: NORTH_CO   !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER :: D_LON_CO   !! LONGITUDE INCREMENT [M_SEC].
INTEGER :: D_LAT_CO   !! LATITUDE INCREMENT [M_SEC].
INTEGER :: N_LON_CO   !! NUMBER OF LONGITUDES.
INTEGER :: N_LAT_CO   !! NUMBER OF LATITUDES.

! ---------------------------------------------------------------------------- !
! 
!    1. RE-FORMAT INPUT PARAMETERS.
!       --------------------------

WEST_CO  =  READ_COOR_TEXT(WEST)
SOUTH_CO =  READ_COOR_TEXT(SOUTH)
IF (PRESENT(EAST)) THEN
   EAST_CO =  READ_COOR_TEXT(EAST)
ELSE
   EAST_CO = COOR_UNDEF
END IF
IF (PRESENT(NORTH)) THEN
   NORTH_CO =  READ_COOR_TEXT(NORTH)
ELSE
   NORTH_CO = COOR_UNDEF
END IF
IF (PRESENT(D_LON)) THEN
   D_LON_CO =  READ_COOR_TEXT(D_LON)
ELSE
   D_LON_CO = COOR_UNDEF
END IF
IF (PRESENT(D_LAT)) THEN
   D_LAT_CO =  READ_COOR_TEXT(D_LAT)
ELSE
   D_LAT_CO = COOR_UNDEF
END IF

IF (PRESENT(N_LON)) THEN
   N_LON_CO =  N_LON
ELSE
   N_LON_CO = -1
END IF

IF (PRESENT(N_LAT)) THEN
   N_LAT_CO =  N_LAT
ELSE
   N_LAT_CO = -1
END IF

! ---------------------------------------------------------------------------- !
! 
!    2. TRANSFER INPUT PARAMETERS.
!       --------------------------

IF (PRESENT(CODE)) THEN
   CALL SET_CURRENT_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,      &
&                             NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,  &
&                             N_LON=N_LON_CO, N_LAT=N_LAT_CO, CODE=CODE)
ELSE
   CALL SET_CURRENT_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,      &
&                             NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,  &
&                             N_LON=N_LON_CO, N_LAT=N_LAT_CO)
END IF

END SUBROUTINE SET_CURRENT_HEADER_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_CURRENT_HEADER_D (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,       &
&                                N_LON, N_LAT, CODE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_CURRENT_HEADER - SET INPUT GRID FOR CURRENT FIELDS.                    !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2011                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS TO WAM_CURRENT_MODULE.             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY DEFINITIONS.                                !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL (KIND=KIND_D), INTENT(IN)           :: WEST   !! WEST LONG. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN)           :: SOUTH  !! SOUTH LAT. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: EAST   !! EAST LONG. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: NORTH  !! NORTH LAT. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: D_LON  !! LONG. INC. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: D_LAT  !! LAT.  INC. OF GRID [DEG].
INTEGER,            INTENT(IN), OPTIONAL :: N_LON  !! NUMBER OF LONG. IN GRID.
INTEGER,            INTENT(IN), OPTIONAL :: N_LAT  !! NUMBER OF LAT. IN GRID.
INTEGER,            INTENT(IN), OPTIONAL :: CODE   !! CURRENT CODE:
                                                   !! 1 = SPEED AND DIRECTION.
                                                   !! OTHERWISE: COMPONENTS (DEFAULT).

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: WEST_CO    !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER :: SOUTH_CO   !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER :: EAST_CO    !! EAST LONGITUDE OF GRID [M_SEC].
INTEGER :: NORTH_CO   !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER :: D_LON_CO   !! LONGITUDE INCREMENT [M_SEC].
INTEGER :: D_LAT_CO   !! LATITUDE INCREMENT [M_SEC].
INTEGER :: N_LON_CO   !! NUMBER OF LONGITUDES.
INTEGER :: N_LAT_CO   !! NUMBER OF LATITUDES.

! ---------------------------------------------------------------------------- !
! 
!    1. RE-FORMAT INPUT PARAMETERS.
!       --------------------------

WEST_CO  =  DEG_TO_M_SEC(WEST)
SOUTH_CO =  DEG_TO_M_SEC(SOUTH)
IF (PRESENT(EAST)) THEN
   EAST_CO =  DEG_TO_M_SEC(EAST)
ELSE
   EAST_CO = COOR_UNDEF
END IF
IF (PRESENT(NORTH)) THEN
   NORTH_CO =  DEG_TO_M_SEC(NORTH)
ELSE
   NORTH_CO = COOR_UNDEF
END IF
IF (PRESENT(D_LON)) THEN
   D_LON_CO =  DEG_TO_M_SEC(D_LON)
ELSE
   D_LON_CO = COOR_UNDEF
END IF
IF (PRESENT(D_LAT)) THEN
   D_LAT_CO =  DEG_TO_M_SEC(D_LAT)
ELSE
   D_LAT_CO = COOR_UNDEF
END IF

IF (PRESENT(N_LON)) THEN
   N_LON_CO =  N_LON
ELSE
   N_LON_CO = -1
END IF

IF (PRESENT(N_LAT)) THEN
   N_LAT_CO =  N_LAT
ELSE
   N_LAT_CO = -1
END IF

! ---------------------------------------------------------------------------- !
! 
!    2. TRANSFER INPUT PARAMETERS.
!       --------------------------

IF (PRESENT(CODE)) THEN
   CALL SET_CURRENT_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,      &
&                             NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,  &
&                             N_LON=N_LON_CO, N_LAT=N_LAT_CO, CODE=CODE)
ELSE
   CALL SET_CURRENT_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,      &
&                             NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,  &
&                             N_LON=N_LON_CO, N_LAT=N_LAT_CO)
END IF

END SUBROUTINE SET_CURRENT_HEADER_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_CURRENT_HEADER_M (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,       &
&                                N_LON, N_LAT, CODE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_CURRENT_HEADER - SET INPUT GRID FOR CURRENT FIELDS.                    !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS TO WAM_CURRENT_MODULE.             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY DEFINITIONS.                                !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,  INTENT(IN)           :: WEST    !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER,  INTENT(IN)           :: SOUTH   !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER,  INTENT(IN), OPTIONAL :: EAST    !! EAST LONGITUDE OF GRID [M_SEC].
INTEGER,  INTENT(IN), OPTIONAL :: NORTH   !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER,  INTENT(IN), OPTIONAL :: D_LON   !! LONGITUDE INCREMENT [M_SEC].
INTEGER,  INTENT(IN), OPTIONAL :: D_LAT   !! LATITUDE INCREMENT [M_SEC].
INTEGER,  INTENT(IN), OPTIONAL :: N_LON   !! NUMBER OF LONGITUDES.
INTEGER,  INTENT(IN), OPTIONAL :: N_LAT   !! NUMBER OF LATITUDES.
INTEGER,  INTENT(IN), OPTIONAL :: CODE    !! CURRENT CODE:
                                          !! 1 = SPEED AND DIRECTION.
                                          !! OTHERWISE: COMPONENTS (DEFAULT).
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

LOGICAL  :: ERROR                         !! ERROR FLAG

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLEAR GRID DEFINITIONS.                                               !
!        -----------------------                                               !

NX_IN  =-1      !! NUMBER OF COLUMNES IN CURRENT INPUT GRID.
NY_IN  =-1      !! NUMBER OF ROWS     IN CURRENT INPUT GRID.
PER   =.FALSE.  !! .TRUE. IF PERIODIC GRID.
CODE_IN = 0     !! CURRENT CODE: 1 = SPEED AND DIRECTION, OTHERWISE: COMPONENTS
DX_IN    = COOR_UNDEF   !! STEPSIZE BETWEEN LONGITUDES.
DY_IN    = COOR_UNDEF   !! STEPSIZE BETWEEN LATITUDES.
SOUTH_IN = COOR_UNDEF   !! MOST SOUTHERN LATITUDE.
NORTH_IN = COOR_UNDEF   !! MOST NORTHERN LATITUDE.
WEST_IN  = COOR_UNDEF   !! LEFT MOST LONGITUDE.
EAST_IN  = COOR_UNDEF   !! RIGHT MOST LONGITUDE.
EQUAL_GRID = .FALSE.    !! INPUT GRID IS NOT EQUAL TO MODEL GRID.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COPY GRID DEFINITIONS AND CONVERT TO M_SEC.                           !
!        ---------------------- --------------------                           !

WEST_IN  = WEST
SOUTH_IN = SOUTH
IF (PRESENT(EAST )) EAST_IN  = EAST
IF (PRESENT(NORTH)) NORTH_IN = NORTH
IF (PRESENT(D_LON)) DX_IN    = D_LON
IF (PRESENT(D_LAT)) DY_IN    = D_LAT
IF (PRESENT(N_LON)) NX_IN    = N_LON
IF (PRESENT(N_LAT)) NY_IN    = N_LAT
IF (PRESENT(CODE))  CODE_IN  = CODE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CHECK GRID.                                                           !
!        -----------                                                           !

CALL CHECK_GRID_DEFINITION (WEST_IN, SOUTH_IN, EAST_IN, NORTH_IN,              &
&                           DX_IN, DY_IN, NX_IN, NY_IN, ERROR)

IF (ERROR) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *        FATAL ERROR IN SUB. SET_CURRENT_HEADER          *'
   WRITE (IU06,*) ' *        ======================================          *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *       A CURRENT GRID COULD NOT BE DEFINED.             *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF


PER = PERIODIC (WEST_IN, EAST_IN, DX_IN, NX_IN)   !! PERIODIC?
EQUAL_GRID = EQUAL_TO_M_GRID (NX_IN, NY_IN, DX_IN, DY_IN,                      &
&                             WEST_IN, SOUTH_IN, EAST_IN, NORTH_IN)

END SUBROUTINE SET_CURRENT_HEADER_M

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_CURRENT_TIMESTEPS (IN, OUT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_CURRENT_TIMESTEPS - SET CURRENT TIMESTEPS IN WAM_CURRENT MODULE.       !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE CURRENT TIMESTEPS TO WAM_CURRENT_MODULE.                  !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY DEFINITIONS.                                !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)           :: IN    !! CURRENT INPUT TIME STEP.
INTEGER, INTENT(IN), OPTIONAL :: OUT   !! CURRENT OUTPUT TIME STEP.

! ---------------------------------------------------------------------------- !

IDELCI = MAX (IN,0)
IF (PRESENT(OUT)) THEN
   IDELCO = OUT
ELSE
   IDELCO = IN
END IF

IF (IDELCI.LE.0) THEN    !! STATIONARY CURRENT.
   IDELCO = 0
   RETURN
END IF

IF (IDELCO.LE.0) THEN    !! CURRENS ARE NOT INTERPOLATED.
   IDELCO = IDELCI
END IF

!         CHECK CURRENT INPUT AND CURRENT OUTPUT TIMESTEP.                     !
!         ------------------------------------------------                     !

IF (IDELCI.LT.IDELCO) THEN
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+       WARNING ERROR IN SUB. SET_CURRENT_TIMESTEPS       +'
   WRITE(IU06,*) '+       ===========================================       +'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+ CURRENT INPUT TIME STEP IS LESS THAN CURRENT OUTPUT STEP+'
   WRITE(IU06,*) '+ CURRENT INPUT TIMESTEP    IN = ', IDELCI, ' SECONDS'
   WRITE(IU06,*) '+ CURRENT OUTPUT TIMESTEP  OUT = ', IDELCO, ' SECONDS'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+ CURRENT INPUT CHANGED TO CURRENT OUTPUT TIME STEP       +'
   WRITE(IU06,*) '+ THE MODEL MAY IGNORE CURRENT FIELDS                     +'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF
IF (MOD(IDELCI,IDELCO).NE.0) THEN
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+       WARNING ERROR IN SUB. SET_CURRENT_TIMESTEPS       +'
   WRITE(IU06,*) '+       ===========================================       +'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+ CURRENT INPUT AND CURRENT OUTPUT TIME STEPS DO NOT HAVE +'
   WRITE(IU06,*) '+ INTEGER RATIO.                                          +'
   WRITE(IU06,*) '+ CURRENT INPUT TIMESTEP    IN = ', IDELCI, ' SECONDS'
   WRITE(IU06,*) '+ CURRENT OUTPUT TIMESTEP  OUT = ', IDELCO, ' SECONDS'
   WRITE(IU06,*) '+ OUTPUT TIMESTEP IS CHANGED TO NEAREST MULTIPLE OF THE   +'
   WRITE(IU06,*) '+ INPUT TIME STEP.                                        +'
   IDELCO = IDELCI/MAX(NINT(REAL(IDELCI)/REAL(IDELCO)),1)
   WRITE(IU06,*) '+ NEW OUTPUT TIMESTEP   IS: ',IDELCO, ' SECONDS'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

END SUBROUTINE SET_CURRENT_TIMESTEPS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WAM_CURRENT (US, DS, CD_START)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   WAM_CURRENT - ROUTINE TO READ AND PROCESS ONE CURRENT FIELD.               !
!                                                                              !
!     H. GUNTHER      GKSS    DECEMBER 2009                                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       READ, INTERPOLATE AND CONVERT INPUT CURRENTS TO WAM CURRENTS.          !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       READ A CURRENT FIELD FROM THE CURRENT FILE (SEARCH FOR IT) AND         !
!       INTERPOLATED IT TO THE WAVE MODEL SEA POINTS.                          !
!       THE INPUT GRID MUST BE ON A LATITUDE/LONGITUDE GRID EITHER PERIODIC    !
!       OR NON PERIODIC.                                                       !
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

REAL,          INTENT(OUT)    :: US(:)    !! OUTPUT CURRENT (U-COMPONENT).
REAL,          INTENT(OUT)    :: DS(:)    !! OUTPUT CURRENT (V-COMPONENT).
CHARACTER (LEN=14),INTENT(IN) :: CD_START !! DATE OF FIELD TO BE LOOKED FOR.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. READ CURRENT DATA AND CHECK DATE.                                        !
!        ------------------------------                                        !

DO
   CALL READ_CURRENT_INPUT

   IF (CD_READ.EQ.CD_START) EXIT
   IF (CD_READ.GT.CD_START) THEN
         WRITE (IU06,*) ' ********************************************'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' *        FATAL ERROR SUB. WAM_CURRENT      *'
         WRITE (IU06,*) ' *        ============================      *'
         WRITE (IU06,*) ' * CURRENT DATE READ IS LATER THAN EXPECTED *'
         WRITE (IU06,*) ' * DATE READ IS       CD_READ = ', CD_READ
         WRITE (IU06,*) ' * DATE EXPECTED IS  CD_START = ', CD_START
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' *     PROGRAM ABORTS  PROGRAM ABORTS       *'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' ********************************************'
         CALL ABORT1
   END IF
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTERPOLATE CURRENT FIELD TO GRID.                                    !
!        ----------------------------------                                    !


US = 0.  !! INITIALISE CURRENT ARRAY WITH ZERO.
DS = 0.
IF (EQUAL_GRID) THEN
   US = PACK(U_IN, L_S_MASK)
   DS = PACK(V_IN, L_S_MASK)
ELSE
   CALL INTERPOLATION_TO_GRID (IU06, PER, DX_IN, DY_IN, WEST_IN, SOUTH_IN,     &
&                              U_IN, US, V_IN, DS)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. TEST OUTPUT.                                                          !
!        ------------                                                          !

IF (ITEST.GE.3) THEN
   IJ = MIN (10,SIZE(US))
   WRITE (IU06,*) ' '
   WRITE (IU06,*) '      SUB. WAM_CURRENT: CURRENT FIELDS CONVERTED TO GRID'
   WRITE (IU06,*) ' '
   WRITE (IU06,*) ' US(1:10) = ', US(1:IJ)
   WRITE (IU06,*) ' DS(1:10) = ', DS(1:IJ)
END IF

END SUBROUTINE WAM_CURRENT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE NOTIM (CD_START, CD_END)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   NOTIM - STEERING SUB IF TIME INTERPOLATION IS NOT WANTED.                  !
!                                                                              !
!     H. GUNTHER    ECMWF    MAY 1990     MODIFIED FOR SUB VERSION.            !
!     H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.                !
!     H. GUNTHER      GKSS    SEPTEMBER 2000   FT90.                           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       NOTIM: NO TIME INTERPOLATION: PROCESS CURRENT FIELDS.                  !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NO TIME INTERPOLATION:                                                 !
!       CURRENT FIELDS ARE PROCESSED EVERY IDELCI SECONDS (U,V),               !
!       THE CURRENTS ARE INTERPOLATED IN SPACE ONLY.                           !
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

CHARACTER (LEN=14), INTENT(IN)  :: CD_START  !! DATE OF FIRST CURRENT FIELD.
CHARACTER (LEN=14), INTENT(IN)  :: CD_END    !! DATE OF LAST CURRENT FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER            :: MP
CHARACTER (LEN=14) :: CDTWIH
REAL               :: US(1:NSEA)  !! OUTPUT CURRENT FIELD ARRAY (U-COMPONENT).
REAL               :: DS(1:NSEA)  !! OUTPUT CURRENT FIELD ARRAY (V-COMPONENT).

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. LOOP OVER OUTPUT CURRENT TIMES.                                       !
!        -------------------------------                                       !

CDTWIH = CD_START
MP = 0
DO WHILE (CDTWIH.LE.CD_END)

!     1.1 READ ONE CURRENT FIELD AND TRANSFORM TO GRID.                        !
!         ---------------------------------------------                        !

   CALL WAM_CURRENT (US, DS, CDTWIH)
   MP = MP + 1

!     1.2 SAVE IN MODULE WAM_CCURRENT.                                         !
!         ----------------------------                                         !

    U_STORE(:,MP)  = US
    V_STORE(:,MP)  = DS
    CD_STORE(MP)   = CDTWIH

   IF (ITEST.GE.3) THEN
      WRITE(IU06,*) '       SUB. NOTIM: NEW CURRENT FILES AT CDTWIH = ', CDTWIH
   END IF

!     1.3 UPDATE CURRENT FIELD REQUEST TIME.                                   !
!         ----------------------------------                                   !

   CALL INCDATE (CDTWIH,IDELCO)
END DO

END SUBROUTINE NOTIM

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TIMIN (CD_START, CD_END)

! ---------------------------------------------------------------------------- !
!                                                                              !
!                                                                              !
!   TIMIN - STEERING MODULE IF TIME INTERPOLATION IS WANTED.                   !
!                                                                              !
!     H. GUNTHER    ECMWF    MAY 1990         MODIFIED FOR SUB VERSION.        !
!     H. GUNTHER    ECMWF    DECEMBER 90      MODIFIED FOR CYCLE_4.            !
!     H. GUNTHER    GKSS     SEPTEMBER 2000   FT90.                            !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TIME INTERPOLATION: PROCESS CURRENT FIELDS.                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CURRENT FIELDS ARE READ IN EVERY IDELCI SECONDS,                       !
!       BI-LINEAR INTERPOLATED IN SPACE, AND TRANSFORMED TO WAM GRID.          !
!       COMONENTS ARE INTERPOLATED LINEARLY IN TIME.                           !
!       CURRENT FIELDS AND SAVED IN WAM_CURRENT_MODULE.                        !
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

CHARACTER (LEN=14), INTENT(IN)  :: CD_START  !! DATE OF FIRST CURRENT FIELD.
CHARACTER (LEN=14), INTENT(IN)  :: CD_END    !! DATE OF LAST CURRENT FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER            :: MP, NTS, N
REAL               :: DEL
CHARACTER (LEN=14) :: CDT1, CDT2, CDTH

REAL               :: US1(1:NSEA)  !! OUTPUT CURRENT FIELD ARRAY (U-COMPONENT).
REAL               :: DS1(1:NSEA)  !! OUTPUT CURRENT FIELD ARRAY (V-COMPONENT).
REAL               :: US2(1:NSEA)  !! OUTPUT CURRENT FIELD ARRAY (U-COMPONENT).
REAL               :: DS2(1:NSEA)  !! OUTPUT CURRENT FIELD ARRAY (V-COMPONENT).

! ---------------------------------------------------------------------------- !
!                                                                              !
!    1. INITIALIZE TIME COUNTER AND FIRST CURRENT FIELD USED FOR INTERPOLATION.!
!       -----------------------------------------------------------------------!

CDT1 = CDCA
US1 = U
DS1 = V
CDTH = CDT1
CALL INCDATE (CDTH,IDELCO)
IF (CD_START.NE.CDTH) THEN
   WRITE(IU06,*) ' *******************************************'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *        FATAL ERROR IN --TIMIN--         *'
   WRITE(IU06,*) ' *        =========================        *'
   WRITE(IU06,*) ' * DATES DO NOT MATCH.                     *'
   WRITE(IU06,*) ' * START DATE FOR CURRENT IS CD_START = ', CD_START
   WRITE(IU06,*) ' * LAST DATE SAVED IN MODULE IS   CDT1 = ', CDT1
   WRITE(IU06,*) ' * PROCESSING WILL BE ABORTED              *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *******************************************'
   CALL ABORT1
ENDIF

CDT2 = CDT1
CALL INCDATE(CDT2,IDELCI)
NTS = IDELCI/IDELCO
DEL = REAL(IDELCO)/REAL(IDELCI)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LOOP OVER INPUT CURRENTFIELDS.                                        !
!        ------------------------------                                        !

MP  = 0
DO

!     2.1 READ ONE CURRENT FIELD AND TRANSFORM TO GRID.                        !
!         ---------------------------------------------                        !

   CDT2 = CDT1
   CALL INCDATE(CDT2,IDELCI)
   CALL WAM_CURRENT (US2, DS2, CDT2)

!     2.2 INTERPOLATE IN TIME AND SAVE CURRENT FIELDS.                         !
!         --------------------------------------------                         !

   CDTH = CDT1
   DO N = 1,NTS
      MP = MP + 1
      CALL INCDATE(CDTH,IDELCO)
       CD_STORE(MP) = CDTH
       U_STORE(:,MP) = US1 + REAL(N)*DEL*(US2-US1)
       V_STORE(:,MP) = DS1 + REAL(N)*DEL*(DS2-DS1)
   END DO

   IF (ITEST.GE.3) THEN
      WRITE(IU06,*) '       SUB. TIMIN: CURRENT FIELDS FOR CD_STORE = ',       &
&                                           CD_STORE(MP-NTS+1:MP)
   END IF

!     2.3 UPDATE CURRENT FIELD REQUEST TIME AND READ NEXT IF REQUESTED.        !
!         -------------------------------------------------------------        !

   US1 = US2
   DS1 = DS2
   CDT1 = CDT2
   CALL INCDATE (CDTH,IDELCI)
   IF (CDTH.GT.CD_END) EXIT
END DO

END SUBROUTINE TIMIN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_CURRENT_MODULE

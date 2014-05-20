SUBROUTINE READ_TOPO_INPUT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READ_TOPO_INPUT - ROUTINE TO READ TOPO FIELDS.                             !
!                                                                              !
!     HEINZ GUNTHER    GKSS    DECEMBER 2009                                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO READ A TOPO FIELD AND TRANSFER IT TO THE WAM_TOPO_MODULE.           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        FORMATTED READ FROM UNIT IU08, FILE08.                                !
!                                                                              !
!       FILE08, WHICH IS DEFINED IN THE USER INPUT, IS ASSIGNED TO IU08.       !
!       THE FILE MUST STORE:                                                   !
!           1. RECORD: THE TOPO DATA HEADER.                                   !
!           FOLLOWING RECORDS: THE TOPO DATA MATRIX.                           !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_TOPO_HEADER TO THE WAM_TOPO_MODULE:                           !
!          INTEGER :: N_LON      !! NUMBER OF LONGITUDES IN GRID.              !
!          INTEGER :: N_LAT      !! NUMBER OF LATITUDES IN GRID.               !
!          REAL*8  :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].          !
!          REAL*8  :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].         !
!          REAL*8  :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: NORTH      !! NORTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: WEST       !! WEST LONGITUDE OF GRID [DEG].              !
!          REAL*8  :: EAST       !! EAST LONGITUDE OF GRID [DEG].              !
!          INTEGER :: CODE       !! TOPO CODE: 1 = TOTAL WATER DEPTH           !
!                                   OTHERWISE ELEVATION OVER NN.               !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_TOPO_FIELD TO THE WAM_TOPO_MODULE:                            !
!          CHARACTER (LEN=14) :: CDTTOR     !! DATE/TIME OF TOPO FIELD.        !
!          REAL               :: D_MAP(:,:) !! TOPO MAP [M].                   !
!                                                                              !
!       THE TOPO DATA MUST BE ON A REGULAR LATITUDE-LONGITUDE GRID ARRANGED    !
!       FROM  WEST TO EAST AND FROM SOUTH TO NORTH, WHICH IS:                  !
!       IN THE ARRAY "D_MAP(I,K)" THE CORNER POINTS ARE:                       !
!                 (    1,    1 ) <==> SOUTH WEST                               !
!                 (N_LON,    1 ) <==> SOUTH EAST                               !
!                 (    1, N_LAT) <==> NORTH WEST                               !
!                 (N_LON, N_LAT) <==> NORTH EAST                               !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1                     !! TERMINATES PROCESSING.

USE WAM_TOPO_MODULE,       ONLY: &
&       SET_TOPO_HEADER,         & !! SETS TOPO HEADER
&       SET_TOPO_FIELD             !! SETS TOPO FIELD

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE, ONLY: IU06, ITEST, IU08, FILE08

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER, PARAMETER :: KIND_D = 8

INTEGER               :: CODE       !! TOPO CODE: 1 = TOTAL WATER DEPTH,
                                    !!      OTHERWISE ELEVATION OVER NN.
INTEGER, SAVE         :: N_LON      !! NUMBER OF LONGITUDES IN GRID.
INTEGER, SAVE         :: N_LAT      !! NUMBER OF LATITUDES IN GRID.
REAL (KIND=KIND_D)    :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: NORTH      !! NORTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: WEST       !! WEST LONGITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: EAST       !! EAST LONGITUDE OF GRID [DEG].
REAL,    ALLOCATABLE  :: D_MAP(:,:) !! TOPOGRAPHY DATA [M].
CHARACTER (LEN=14)    :: CDTTOR     !! DATE/TIME OF TOPO FIELD

LOGICAL, SAVE  :: FRSTIME = .TRUE.
INTEGER        :: LEN, IOS, J


! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. FOR FIRST CALL: OPEN FILE AND READ HEADER.                            !
!        ------------------------------------------                            !

IF (FRSTIME) THEN
   LEN = LEN_TRIM(FILE08)
   OPEN (UNIT=IU08,FILE=FILE08(1:LEN),FORM='FORMATTED',STATUS='OLD',IOSTAT=IOS)
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_TOPO_INPUT        *'
      WRITE (IU06,*) ' *       ===================================        *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * TOPO INPUT FILE COULD NOT BE OPENED              *'
      WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
      WRITE (IU06,*) ' *    FILE NAME IS  FILE08 = ', TRIM(FILE08)
      WRITE (IU06,*) ' *    UNIT IS         IU08 = ', IU08
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   END IF

   READ (IU08, '(6F10.5,3I3)', IOSTAT=IOS) SOUTH, NORTH, WEST, EAST,          &
&                                          D_LON, D_LAT, N_LAT, N_LON, CODE
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_TOPO_INPUT        *'
      WRITE (IU06,*) ' *       ===================================        *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * READ ERROR ON TOPO FILE.                         *'
      WRITE (IU06,*) ' * FILE HEADER EXPECTED                             *'
      WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
      WRITE (IU06,*) ' *    FILE NAME IS  FILE08 = ', TRIM(FILE08)
      WRITE (IU06,*) ' *    UNIT IS         IU08 = ', IU08
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   END IF

   CALL SET_TOPO_HEADER (WEST=WEST,   SOUTH=SOUTH,    &
&                        EAST=EAST,   NORTH=NORTH,    &
&                        D_LON=D_LON, D_LAT=D_LAT,    &
&                        N_LON=N_LON, N_LAT=N_LAT,    &
&                        CODE=CODE)
   FRSTIME = .FALSE.
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. ALLOCATE TOPO INPUT ARRAYS.                                           !
!        ---------------------------                                           !

IF (.NOT.ALLOCATED(D_MAP) ) ALLOCATE(D_MAP(N_LON,N_LAT))

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. READ TOPO FIELD.                                                       !
!       -------------------                                                    !

READ(IU08,'(A14)',IOSTAT=IOS) CDTTOR
IF (IOS.NE.0) THEN
   WRITE (IU06,*) ' ****************************************************'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_TOPO_INPUT        *'
   WRITE (IU06,*) ' *       ===================================        *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' * READ ERROR ON TOPO FILE.                         *'
   WRITE (IU06,*) ' * DATE/TIME OF TOPOFIELD EXPECTED                  *'
   WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
   WRITE (IU06,*) ' *    FILE NAME IS  FILE08 = ', TRIM(FILE08)
   WRITE (IU06,*) ' *    UNIT IS         IU08 = ', IU08
   WRITE (IU06,*) ' *    RECORD IS          J = ', J
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' ****************************************************'
   CALL ABORT1
END IF

DO J = 1,N_LAT
   READ(IU08,'(8F10.3)',IOSTAT=IOS) D_MAP(1:N_LON,J)
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_TOPO_INPUT        *'
      WRITE (IU06,*) ' *       ===================================        *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * READ ERROR ON TOPO FILE.                         *'
      WRITE (IU06,*) ' * U - COMPONENTS EXPECTED                          *'
      WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
      WRITE (IU06,*) ' *    FILE NAME IS  FILE08 = ', TRIM(FILE08)
      WRITE (IU06,*) ' *    UNIT IS         IU08 = ', IU08
      WRITE (IU06,*) ' *    RECORD IS          J = ', J
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   END IF
END DO

CALL SET_TOPO_FIELD (CDTTOR, D_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    4. WRITE TEST OUTPUT AND DEALLOCATE ARRAYS.                               !
!       ----------------------------------------                               !

IF (ITEST.GT.1) THEN
   WRITE(IU06,*) ' '
   WRITE(IU06,*) '     SUB. READ_TOPO_INPUT -  TOPO FIELD FOR CDTTOR = ',CDTTOR
   WRITE(IU06,'(5X,20F6.0)') D_MAP(1:MIN(20,N_LON),1:MIN(5,N_LAT))
END IF

DEALLOCATE(D_MAP)

END SUBROUTINE READ_TOPO_INPUT

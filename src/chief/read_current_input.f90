SUBROUTINE READ_CURRENT_INPUT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READ_CURRENT_INPUT - READ CURRENT FIELDS.                                  !
!                                                                              !
!     HEINZ GUNTHER    GKSS    DECEMBER 2009                                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO READ A CURRENT FIELD AND TRANSFER IT TO THE WAM_CURRENT_MODULE.     !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       FORMATTED READ FROM UNIT IU09, FILE09.                                 !
!                                                                              !
!       FILE09, WHICH IS DEFINED IN THE USER INPUT, IS ASSIGNED TO IU09.       !
!       THE FILE MUST STORE:                                                   !
!           1. RECORD: THE CURRENT DATA HEADER.                                !
!           FOLLOWING RECORDS: THE CURRENT DATA MATRIX.                        !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_CURRENT_HEADER TO THE WAM_CURRENT_MODULE:                     !
!          INTEGER :: N_LON      !! NUMBER OF LONGITUDES IN GRID.              !
!          INTEGER :: N_LAT      !! NUMBER OF LATITUDES IN GRID.               !
!          REAL*8  :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].          !
!          REAL*8  :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].         !
!          REAL*8  :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: NORTH      !! NORTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: WEST       !! WEST LONGITUDE OF GRID [DEG].              !
!          REAL*8  :: EAST       !! EAST LONGITUDE OF GRID [DEG].              !
!          INTEGER :: CODE       !! CURRENT CODE:                              !
!                                !! 1 = SPEED (U-MAP) AND DIRECTION (V_MAP),   !
!                                !! OTHERWISE: COMPONENTS                      !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_CURRENT_FIELD TO THE WAM_CURRENT_MODULE:                      !
!          CHARACTER (LEN=14) :: CDTCUR     !! DATE/TIME OF CURRENT FIELD.     !
!          REAL               :: U_MAP(:,:) !! U COMPONENT OF CURRENT MAP [M/S]!
!          REAL               :: V_MAP(:,:) !! V COMPONENT OF CURRENT MAP [M/S]!
!                                                                              !
!       THE CURRENTS MUST BE ON A REGULAR LATITUDE-LONGITUDE GRID ARRANGED     !
!       FROM  WEST TO EAST AND FROM SOUTH TO NORTH, WHICH IS:                  !
!       IN THE ARRAYS "U_MAP(I,K)" AND  "V_MAP(I,K)" THE CORNER POINTS ARE:    !
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
&       ABORT1                       !! TERMINATES PROCESSING.

USE WAM_CURRENT_MODULE,   ONLY:  &
&       SET_CURRENT_HEADER,         & !! SETS CURRENT HEADER
&       SET_CURRENT_FIELD             !! SETS CURRENT FIELD

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY:  IU06, ITEST, IU09, FILE09

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER, PARAMETER :: KIND_D = 8

INTEGER               :: CODE       !! CURRENT CODE:
                                    !! 1 = SPEED (U-MAP) AND DIRECTION (V_MAP),
                                    !! OTHERWISE: COMPONENTS
INTEGER, SAVE         :: N_LON      !! NUMBER OF LONGITUDES IN GRID.
INTEGER, SAVE         :: N_LAT      !! NUMBER OF LATITUDES IN GRID.
REAL (KIND=KIND_D)    :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: NORTH      !! NORTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: WEST       !! WEST LONGITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: EAST       !! EAST LONGITUDE OF GRID [DEG].
REAL,    ALLOCATABLE  :: U_MAP(:,:) !! 1. COMPONENT OF CURRENT MAP [M/S].
REAL,    ALLOCATABLE  :: V_MAP(:,:) !! 2. COMPONENT OF CURRENT MAP [M/S].
CHARACTER (LEN=14)    :: CDTCUR     !! DATE/TIME OF CURRENT FIELD

LOGICAL, SAVE  :: FRSTIME = .TRUE.
INTEGER        :: LEN, IOS, J

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. FOR FIRST CALL: OPEN FILE AND READ HEADER.                            !
!        ------------------------------------------                            !

IF (FRSTIME) THEN
   LEN = LEN_TRIM(FILE09)
   OPEN (UNIT=IU09,FILE=FILE09(1:LEN),FORM='FORMATTED',STATUS='OLD',IOSTAT=IOS)
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READ_CURRENT_INPUT      *'
      WRITE (IU06,*) ' *      ======================================      *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * CURRENT INPUT FILE COULD NOT BE OPENED           *'
      WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
      WRITE (IU06,*) ' *    FILE NAME IS  FILE09 = ', TRIM(FILE09)
      WRITE (IU06,*) ' *    UNIT IS         IU09 = ', IU09
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   END IF

   READ (IU09, '(6F10.5,3I3)', IOSTAT=IOS) SOUTH, NORTH, WEST, EAST,          &
&                                          D_LON, D_LAT, N_LAT, N_LON, CODE
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READ_CURRENT_INPUT      *'
      WRITE (IU06,*) ' *      ======================================      *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * READ ERROR ON CURRENT FILE.                      *'
      WRITE (IU06,*) ' * FILE HEADER EXPECTED                             *'
      WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
      WRITE (IU06,*) ' *    FILE NAME IS  FILE09 = ', TRIM(FILE09)
      WRITE (IU06,*) ' *    UNIT IS         IU09 = ', IU09
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   END IF

   CALL SET_CURRENT_HEADER (WEST=WEST,   SOUTH=SOUTH,    &
&                           EAST=EAST,   NORTH=NORTH,    &
&                           D_LON=D_LON, D_LAT=D_LAT,    &
&                           N_LON=N_LON, N_LAT=N_LAT,    &
&                           CODE=CODE)
   FRSTIME = .FALSE.
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. ALLOCATE CURRENT INPUT ARRAYS.                                        !
!        ------------------------------                                        !

IF (.NOT.ALLOCATED(U_MAP)) ALLOCATE(U_MAP(N_LON,N_LAT))
IF (.NOT.ALLOCATED(V_MAP)) ALLOCATE(V_MAP(N_LON,N_LAT))

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. READ CURRENT FIELD.                                                    !
!       -------------------                                                    !

READ(IU09,'(A14)',IOSTAT=IOS) CDTCUR
IF (IOS.NE.0) THEN
   WRITE (IU06,*) ' ****************************************************'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READ_CURRENT_INPUT      *'
   WRITE (IU06,*) ' *      ======================================      *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' * READ ERROR ON CURRENT FILE.                      *'
   WRITE (IU06,*) ' * DATE/TIME OF WINDFIELD EXPECTED                  *'
   WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
   WRITE (IU06,*) ' *    FILE NAME IS  FILE09 = ', TRIM(FILE09)
   WRITE (IU06,*) ' *    UNIT IS         IU09 = ', IU09
   WRITE (IU06,*) ' *    RECORD IS          J = ', J
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' ****************************************************'
   CALL ABORT1
END IF

DO J = 1,N_LAT
   READ(IU09,'(8E10.3)',IOSTAT=IOS) U_MAP(1:N_LON,J)
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READ_CURRENT_INPUT      *'
      WRITE (IU06,*) ' *      ======================================      *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * READ ERROR ON CURRENT FILE.                      *'
      WRITE (IU06,*) ' * U - COMPONENTS EXPECTED                          *'
      WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
      WRITE (IU06,*) ' *    FILE NAME IS  FILE09 = ', TRIM(FILE09)
      WRITE (IU06,*) ' *    UNIT IS         IU09 = ', IU09
      WRITE (IU06,*) ' *    RECORD IS          J = ', J
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   END IF
END DO
DO J = 1,N_LAT
   READ(IU09,'(8E10.3)',IOSTAT=IOS) V_MAP(1:N_LON,J)
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READ_CURRENT_INPUT      *'
      WRITE (IU06,*) ' *      ======================================      *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * READ ERROR ON CURRENT FILE.                      *'
      WRITE (IU06,*) ' * V - COMPONENTS EXPECTED                          *'
      WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
      WRITE (IU06,*) ' *    FILE NAME IS  FILE09 = ', TRIM(FILE09)
      WRITE (IU06,*) ' *    UNIT IS         IU09 = ', IU09
      WRITE (IU06,*) ' *    RECORD IS          J = ', J
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   END IF
END DO

CALL SET_CURRENT_FIELD (CDTCUR, U_MAP, V_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    4. WRITE TEST OUTPUT AND DEALLOCATE ARRAYS.                               !
!       ----------------------------------------                               !

IF (ITEST.GT.1) THEN
   WRITE(IU06,*) ' '
   WRITE(IU06,*) '     SUB. READ_CURRENT_INPUT: FIELD FOR CDTCUR = ',CDTCUR
   WRITE(IU06,'(5X,20F6.2)') U_MAP(1:MIN(20,N_LON),1:MIN(5,N_LAT))
   WRITE(IU06,'(5X,20F6.2)') V_MAP(1:MIN(20,N_LON),1:MIN(5,N_LAT))
END IF

DEALLOCATE (U_MAP, V_MAP)

END SUBROUTINE READ_CURRENT_INPUT

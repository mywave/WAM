SUBROUTINE READ_ICE_INPUT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READ_ICE_INPUT - READ AN ICE MAP.                                          !
!                                                                              !
!     HEINZ GUNTHER    GKSS    JANUARY 1995                                    !
!     ERIK MYKLEBUST           NOVEMBER 2004                                   !
!                                                                              !
!     PURPOSE                                                                  !
!     -------                                                                  !
!                                                                              !
!       TO READ AN ICE MAP.                                                    !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        FORMATTED READ FROM UNIT IU03, FILE03.                                !
!                                                                              !
!       FILE03, WHICH IS DEFINED IN THE USER INPUT, IS ASSIGNED TO IU03.       !
!       THE FILE MUST STORE:                                                   !
!           1. RECORD: THE ICE DATA HEADER.                                    !
!           2. RECORD: THE ICE DATA TIME.                                      !
!           FOLLOWING RECORDS: THE ICE DATA MATRIX.                            !
!           (RECORD 2 AND FOLLOWING RECORDS ARE REPEATED FOR NEXT ICE FIELD.)  !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_ICE_HEADER TO THE WAM_WIND_MODULE:                            !
!          INTEGER :: N_LON      !! NUMBER OF LONGITUDES IN GRID.              !
!          INTEGER :: N_LAT      !! NUMBER OF LATITUDES IN GRID.               !
!          REAL*8  :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].          !
!          REAL*8  :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].         !
!          REAL*8  :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: NORTH      !! NORTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: WEST       !! WEST LONGITUDE OF GRID [DEG].              !
!          REAL*8  :: EAST       !! EAST LONGITUDE OF GRID [DEG].              !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_ICE TO THE WAM_ICE_MODULE:                                    !
!          CHARACTER (LEN=14) :: CDATE     !! DATE/TIME OFICE FIELD.           !
!          INTEGER            :: ICE_GRID(:,:) !! ICE MAP.                     !
!                                                                              !
!       THE ICE MAP MUST BE ON A REGULAR LATITUDE-LONGITUDE GRID ARRANGED      !
!       FROM  WEST TO EAST AND FROM SOUTH TO NORTH, WHICH IS:                  !
!                       (1,  1 ) <==> SOUTH WEST                               !
!                       (NX, 1 ) <==> SOUTH EAST                               !
!                       (1,  NY) <==> NORTH WEST                               !
!                       (NX, NY) <==> NORTH EAST                               !
!                                                                              !
!        ICE POINTS MUST BE EQUAL TO 1 IN ICE_GRID.                            !
!        THE SAME SUB. SET_ICE CAN BE USED, IF ICE_GRID IS DEFINED AS          !
!        REAL OR LOGICAL. IN THIS CASE ICE POINTS ARE WHERE                    !
!           NINT(ICE_GRID) = 1  OR ICE_GRID =.TRUE. , RESPECTIVELY.            !
!                                                                              !
!     REFERENCES                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1                     !! TERMINATES PROCESSING.

USE WAM_ICE_MODULE,       ONLY:  &
&       SET_ICE,                 & !! ICE INPUT INTO MODULE.
&       SET_ICE_HEADER             !! ICE INPUT HEADER INTO MODULE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY: IU06, ITEST, IU03, FILE03

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER, PARAMETER :: KIND_D = 8

INTEGER, SAVE         :: NX_ICE        !! ICE MAP DIMENSIONS
INTEGER, SAVE         :: NY_ICE        !! ICE MAP DIMENSIONS
REAL (KIND=KIND_D)    :: D_LAT         !! LATITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: D_LON         !! LONGITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: SOUTH         !! SOUTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: NORTH         !! NORTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: WEST          !! WEST LONGITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: EAST          !! EAST LONGITUDE OF GRID [DEG].
INTEGER, ALLOCATABLE  :: ICE_GRID(:,:) !! ICE MAP 
CHARACTER (LEN=14)    :: CDATE         !! ICE DATE

LOGICAL, SAVE         :: FRSTIME = .TRUE.

INTEGER               :: K, IOS
LOGICAL, SAVE         :: FORMATTED

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. OPEN ICE DATA FILE.                                                   !
!        -------------------                                                   !

IF (FRSTIME) THEN
   FORMATTED = .TRUE.
   IOS = 0
   OPEN (UNIT=IU03, FILE=TRIM(FILE03), FORM='FORMATTED', STATUS='OLD',         &
&        IOSTAT=IOS)

   READ (IU03,*,IOSTAT=IOS) D_LAT, D_LON, SOUTH, NORTH, WEST, EAST           
!   READ (IU03,'(6F10.5)',IOSTAT=IOS) D_LAT, D_LON, SOUTH, NORTH, WEST, EAST

   IF (IOS.NE.0) THEN       !! NOT ASCII FILE. TRY UNFORMATTED
      WRITE(IU06,*) 'FORMATTED READING OF FILE ', TRIM(FILE03), ' FAILED'
      WRITE(IU06,*) 'TRYING UNFORMATTED READING'
      CLOSE(UNIT=IU03, STATUS='KEEP')
      OPEN (UNIT=IU03, FILE=TRIM(FILE03), FORM='UNFORMATTED', STATUS='OLD',    &
&           IOSTAT=IOS)
      IF (IOS.NE.0) THEN
         WRITE(IU06,*) ' *******************************************'
         WRITE(IU06,*) ' *                                         *'
         WRITE(IU06,*) ' *   FATAL ERROR IN SUB. READ_ICE_INPUT    *'
         WRITE(IU06,*) ' *   ==================================    *'
         WRITE(IU06,*) ' *                                         *'
         WRITE(IU06,*) ' * OPEN ERROR IN ICE FILE                  *'
         WRITE(IU06,*) ' * ERROR CODE IS IOSTAT = ', IOS
         WRITE(IU06,*) ' * UNIT IS         IU03 = ', IU03
         WRITE(IU06,*) ' * UNIT IS       FILE03 = ', TRIM(FILE03)
         WRITE(IU06,*) ' *                                         *'
         WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
         WRITE(IU06,*) ' *                                         *'
         WRITE(IU06,*) ' *******************************************'
         CALL ABORT1
      END IF
      FORMATTED = .FALSE.
      READ (IU03,IOSTAT=IOS) D_LAT, D_LON, SOUTH, NORTH, WEST, EAST
      IF (IOS.NE.0) THEN
         WRITE(IU06,*) ' *******************************************'
         WRITE(IU06,*) ' *                                         *'
         WRITE(IU06,*) ' *   FATAL ERROR IN SUB. READ_ICE_INPUT    *'
         WRITE(IU06,*) ' *   ==================================    *'
         WRITE(IU06,*) ' *                                         *'
         WRITE(IU06,*) ' * READ ERROR IN ICE GRID HEADER RECORD    *'
         WRITE(IU06,*) ' * ERROR CODE IS IOSTAT = ', IOS
         WRITE(IU06,*) ' * UNIT IS         IU03 = ', IU03
         WRITE(IU06,*) ' * UNIT IS       FILE03 = ', TRIM(FILE03)
         WRITE(IU06,*) ' *                                         *'
         WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
         WRITE(IU06,*) ' *                                         *'
         WRITE(IU06,*) ' *******************************************'
         CALL ABORT1
      END IF
   END IF
   CALL SET_ICE_HEADER (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT)

   FRSTIME = .FALSE.
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. READ ICE DATA HEADER.                                                 !
!        ---------------------                                                 !

IF (FORMATTED) THEN
   READ (UNIT=IU03, FMT='(A14,2I10)', IOSTAT=IOS) CDATE, NX_ICE, NY_ICE
ELSE
   READ (UNIT=IU03, IOSTAT=IOS) CDATE, NX_ICE, NY_ICE
END IF
IF (IOS.NE.0) THEN
   WRITE(IU06,*) ' *******************************************'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *   FATAL ERROR IN SUB. READ_ICE_INPUT    *'
   WRITE(IU06,*) ' *   ==================================    *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' * READ ERROR IN ICE GRID HEADER RECORD    *'
   WRITE(IU06,*) ' * ERROR CODE IS IOSTAT = ', IOS
   WRITE(IU06,*) ' * UNIT IS         IU03 = ', IU03
   WRITE(IU06,*) ' * UNIT IS       FILE03 = ', TRIM(FILE03)
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *******************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. READ ICE DATA.                                                        !
!        --------------                                                        !

IF (.NOT.ALLOCATED(ICE_GRID)) ALLOCATE (ICE_GRID(1:NX_ICE,1:NY_ICE))

DO K=1,NY_ICE
   IF (FORMATTED) THEN
      READ (UNIT=IU03, FMT='(80I1)', IOSTAT=IOS) ICE_GRID(1:NX_ICE,K)
   ELSE
      READ (UNIT=IU03, IOSTAT=IOS) ICE_GRID(1:NX_ICE,K)
   END IF
   IF (IOS.NE.0) THEN
      WRITE(IU06,*) ' *******************************************'
      WRITE(IU06,*) ' *                                         *'
      WRITE(IU06,*) ' *   FATAL ERROR IN SUB. READ_ICE_INPUT    *'
      WRITE(IU06,*) ' *   ==================================    *'
      WRITE(IU06,*) ' *                                         *'
      WRITE(IU06,*) ' * READ ERROR ON ICE GRID                  *'
      WRITE(IU06,*) ' * ERROR CODE IS IOSTAT = ', IOS
      WRITE(IU06,*) ' * GRID LINE  IS      K = ', K
      WRITE(IU06,*) ' * UNIT IS         IU03 = ', IU03
      WRITE(IU06,*) ' * UNIT IS       FILE03 = ', TRIM(FILE03)
      WRITE(IU06,*) ' *                                         *'
      WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
      WRITE(IU06,*) ' *                                         *'
      WRITE(IU06,*) ' *******************************************'
      CALL ABORT1
   END IF
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. TRANSFER ICE DATA TO MODULE.                                          !
!        ----------------------------                                          !

CALL SET_ICE (CDATE, ICE_GRID)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      5. WRITE TEST OUTPUT AND DEALLOCATE ARRAY.                              !
!         ---------------------------------------                              !

IF (ITEST.GE.1) THEN
   WRITE(IU06,*) ' '
   WRITE(IU06,*) '     SUB. READ_ICE_INPUT: FIELD FOR CDATE = ',CDATE
END IF   

IF (ALLOCATED(ICE_GRID)) DEALLOCATE (ICE_GRID)

END SUBROUTINE READ_ICE_INPUT

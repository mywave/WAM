SUBROUTINE READ_TOPOGRAPHY

! ---------------------------------------------------------------------------- !
!                                                                              !
!    READ_TOPOGRAPHY - READ TOPOGRAPHY.                                        !
!                                                                              !
!     H. GUNTHER      GKSS  JANUARY 2002                                       !
!     A. Behrens      HZG   January 2014  Topography real                      !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO READ A TOPOGRAPHY FILE AND TRANSFER THE DATA TO THE PREPROC MODULE. !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       FILE08, WHICH IS DEFINED IN THE USER INPUT, IS ASSIGNED TO IU08.       !
!       THE FILE MUST STORE:                                                   !
!           1. RECORD: THE TOPOGRAPHIC DATA HEADER.                            !
!           FOLLOWING RECORDS: THE DEPTH DATA MATRIX.                          !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY SUB. SET_TOPOGRAPHY !
!       TO THE WAM_GRID_MODULE:                                                !
!         INTEGER    :: N_LON      !! NUMBER OF LONGITUDES.                    !
!         INTEGER    :: N_LAT      !! NUMBER OF LATITUDES.                     !
!         REAL*8     :: D_LAT      !! LATITUDE INCREMENT.                      !
!         REAL*8     :: D_LON      !! LONGITUDE INCREMENT.                     !
!         REAL*8     :: SOUTH      !! SOUTH LATITUDE.                          !
!         REAL*8     :: NORTH      !! NORTH LATITUDE.                          !
!         REAL*8     :: WEST       !! WEST LONGITUDE.                          !
!         REAL*8     :: EAST       !! EAST LONGITUDE.                          !
!         REAL       :: D_MAP(:,:) !! WATER DEPTH [M].                         !
!                                                                              !
!      ALL INCREMENTS, LATITUDES AND LONGITUDES MUST BE REAL*8 IN DEGREES      !
!      OR  INTEGER IN M_SEC, OR A CHARCTER STRING                              !
!                                                                              !
!      THE TOPOGRAPHY MUST BE ON A REGULAR LATITUDE-LONGITUDE GRID ARRANGED    !
!      FROM  WEST TO EAST AND FROM SOUTH TO NORTH, WHICH IS                    !
!      THE DEPTH ARRAY "D_MAP(I,K)" MUST BE ORDERED AS                         !
!                 (    1,    1 ) <==> SOUTH WEST                               !
!                 (N_LON,    1 ) <==> SOUTH EAST                               !
!                 (    1, N_LAT) <==> NORTH WEST                               !
!                 (N_LON, N_LAT) <==> NORTH EAST                               !
!       POSITIVE VALUES ARE SEA DEPTHS AND NEGATIVE VALUES ARE LAND.           !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !

USE WAM_COORDINATE_MODULE           !! COORDINATE PROCEDURES

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1                      !! TERMINATES PROCESSING.

USE PREPROC_MODULE,       ONLY:  &
&       SET_TOPOGRAPHY              !! TRANSFERS DEPTH DATA TO MODULE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,   ONLY: IU06, IU08, FILE08

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER               :: N_LON      !! NUMBER OF LONGITUDES IN GRID.
INTEGER               :: N_LAT      !! NUMBER OF LATITUDES IN GRID.
REAL (KIND=KIND_D)    :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: NORTH      !! NORTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: WEST       !! WEST LONGITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: EAST       !! EAST LONGITUDE OF GRID [DEG].
REAL, ALLOCATABLE     :: D_MAP(:,:) !! WATER DEPTH [M].

INTEGER                  :: I, K, L, IMAX, IOS, IA, IE
CHARACTER*120            :: LINE
CHARACTER*15             :: IFORM
CHARACTER*1, ALLOCATABLE :: AX(:)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. OPEN TOPOGRAPHY DATA FILE.                                            !
!        --------------------------                                            !

L = LEN_TRIM(FILE08)
OPEN (UNIT=IU08, FILE=FILE08(1:L), FORM='FORMATTED', STATUS='OLD', IOSTAT=IOS)
IF (IOS.NE.0) THEN
   WRITE (IU06,*) ' *****************************************************'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' *     FATAL  ERROR IN SUB. READ_TOPOGRAPHY          *'
   WRITE (IU06,*) ' *     =====================================         *'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' * ERROR WHEN OPENING TOPOGRAPHY_FILE                *'
   WRITE (IU06,*) ' * FILE NAME IS            FILE08 = ', FILE08(1:L)
   WRITE (IU06,*) ' * ASSIGNED TO UNIT          IU08 = ', IU08
   WRITE (IU06,*) ' * ERROR CODE             IOSTAT  = ', IOS
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' *      PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' *****************************************************'
   CALL ABORT1
ELSE

   WRITE(IU06,*) ' SUB. READ TOPOGRAPHY: FILE CONNECTED TO UNIT =', IU08,      &
&                ' FILE NAME IS: ',FILE08(1:L)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. READ HEADER.                                                          !
!        ------------                                                          !

READ (IU08,*,IOSTAT=IOS) D_LAT, D_LON, SOUTH, NORTH, WEST, EAST
IF (IOS.NE.0) THEN
   WRITE (IU06,*) ' *****************************************************'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' *     FATAL  ERROR IN SUB. READ_TOPOGRAPHY          *'
   WRITE (IU06,*) ' *     ====================================          *'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' * READ ERROR ON TOPOGRAPHY_FILE                     *'
   WRITE (IU06,*) ' * FILE NAME IS            FILE08 = ', FILE08(1:L)
   WRITE (IU06,*) ' * ASSIGNED TO UNIT          IU08 = ', IU08
   WRITE (IU06,*) ' * PROGRAM TRIES TO READ HEADER INFORMATION          *'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' *      PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' *****************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. PROCESS HEADER.                                                       !
!        ---------------                                                       !

CALL ADJUST (WEST, EAST)                      !!  WEST < EAST
N_LON = NINT((EAST-WEST)/D_LON+1.0)           !! NO. OF LONGITUDES IN FILE
N_LAT = NINT((NORTH-SOUTH)/D_LAT+1.0)         !! NO. OF LATITUDES IN FILE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. READ THE TOPOGRAPHY.                                                  !
!        --------------------                                                  !

ALLOCATE (D_MAP(N_LON,N_LAT), AX(N_LON))        !! INPUT ARRAYS.
IMAX = (N_LON+11)/12        !! NUMBER OF RECORDS OF INPUT GRID PER LATITUDE.

READ (IU08,'(A)',IOSTAT=IOS) LINE
BACKSPACE(IU08)
if(INDEX(LINE,'.').gt.0) THEN
   IFORM='(12(f8.2,A1))'    !! read topography in m and cm
ELSE
   IFORM='(12(f5.0,A1))'    !! read topography in m only
END IF

DO K = 1,N_LAT
   DO I = 1,IMAX
      IA = 12*(I-1)+1
      IE = MIN(12*I,N_LON)
      READ (IU08,trim(IFORM)) (D_MAP(L,K),AX(L),L=IA,IE)
      IF (IOS.NE.0) THEN
         WRITE (IU06,*) ' *****************************************************'
         WRITE (IU06,*) ' *                                                   *'
         WRITE (IU06,*) ' *     FATAL  ERROR IN SUB. READ_TOPOGRAPHY          *'
         WRITE (IU06,*) ' *     ====================================          *'
         WRITE (IU06,*) ' *                                                   *'
         WRITE (IU06,*) ' * READ ERROR ON TOPOGRAPHY_FILE                     *'
         WRITE (IU06,*) ' * FILE NAME IS            FILE08 = ', FILE08
         WRITE (IU06,*) ' * ASSIGNED TO UNIT          IU08 = ', IU08
         WRITE (IU06,*) ' * LATITUDE  NO. IS             K = ', K
         WRITE (IU06,*) ' * LONGITUDE SECTION IS         I = ', I
         WRITE (IU06,*) ' * LONGITUDE NO. IS FROM IA = ', IA,' TO IE = ', IE
         WRITE (IU06,*) ' *                                                   *'
         WRITE (IU06,*) ' *      PROGRAM ABORTS  PROGRAM ABORTS               *'
         WRITE (IU06,*) ' *                                                   *'
         WRITE (IU06,*) ' *****************************************************'
         CALL ABORT1
      END IF
   END DO

   WHERE (AX(1:N_LON).EQ.'E'.AND.D_MAP(1:N_LON,K).LT.0)                        &
&                                      D_MAP(1:N_LON,K) = -D_MAP(1:N_LON,K)
   WHERE (AX(1:N_LON).EQ.'D'.AND.D_MAP(1:N_LON,K).GT.0)                        &
&                                      D_MAP(1:N_LON,K) = -D_MAP(1:N_LON,K)
   WHERE (AX(1:N_LON).EQ.'E'.AND.D_MAP(1:N_LON,K).EQ.0) D_MAP(1:N_LON,K) = 1
   WHERE (AX(1:N_LON).EQ.'D'.AND.D_MAP(1:N_LON,K).EQ.0) D_MAP(1:N_LON,K) = -1

END DO

D_MAP = -D_MAP   !! SEA DEPTH MUST BE POSITIVE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. TRANSFER TO MODEULE.                                                  !
!        --------------------                                                  !

CALL SET_TOPOGRAPHY (N_LON, N_LAT, D_LON, D_LAT, SOUTH, NORTH, WEST, EAST,     &
&                    REAL(D_MAP))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. DEALLOCATE ARRAYS AND CLOSE FILE.                                     !
!        ---------------------------------                                     !

DEALLOCATE (D_MAP, AX)
CLOSE (UNIT=IU08)

END SUBROUTINE READ_TOPOGRAPHY

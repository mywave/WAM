SUBROUTINE READ_TOPOGRAPHY

! ---------------------------------------------------------------------------- !
!                                                                              !
!    READ_TOPOGRAPHY - READ TOPOGRAPHY.   Formatted BSH Topography             !
!                                                                              !
!     H. GUNTHER      GKSS  JANUARY 2002                                       !
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
real, ALLOCATABLE  :: D_MAP(:,:)    !! WATER DEPTH [M].

INTEGER                  :: I, K, L, IOS

Real :: xlat,xlon, Tief
Integer ::iin, kin, NUM

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
READ (IU08,*,IOSTAT=IOS)
READ (IU08,*,IOSTAT=IOS)
READ (IU08,*,IOSTAT=IOS)
READ (IU08,*,IOSTAT=IOS)
READ (IU08,*,IOSTAT=IOS)

D_LAT = 0.008333
D_LON = 0.013889

NORTH = 56.445835
WEST = 6.173611
N_LON = 630           !! NO. OF LONGITUDES IN FILE
N_LAT = 387           !! NO. OF LATITUDES IN FILE
EAST  =  WEST+Real(N_LON-1)*D_LON
SOUTH = NORTH-Real(N_LAT-1)*D_LAT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. PROCESS HEADER.                                                       !
!        ---------------                                                       !

CALL ADJUST (WEST, EAST)                      !!  WEST < EAST

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. READ THE TOPOGRAPHY.                                                  !
!        --------------------                                                  !

ALLOCATE (D_MAP(N_LON,N_LAT))        !! INPUT ARRAYS.
D_MAP = -999.           !! NUMBER OF RECORDS OF INPUT GRID PER LATITUDE.

IOS = 0
DO While (Ios.eq.0)
      READ (IU08,*,IOSTAT=IOS) xlat, xlon, Tief, Kin, iin, NUM
      i = iin
      k = N_LAT+1-kin
      D_MAP(i,k) = tief
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. TRANSFER TO MODEULE.                                                  !
!        --------------------                                                  !

CALL SET_TOPOGRAPHY (N_LON, N_LAT, D_LON, D_LAT, SOUTH, NORTH, WEST, EAST,     &
&                    D_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. DEALLOCATE ARRAYS AND CLOSE FILE.                                     !
!        ---------------------------------                                     !

DEALLOCATE (D_MAP)
CLOSE (UNIT=IU08)

END SUBROUTINE READ_TOPOGRAPHY

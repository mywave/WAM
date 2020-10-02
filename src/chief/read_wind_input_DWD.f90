SUBROUTINE READ_WIND_INPUT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READ_WIND_INPUT - ROUTINE TO READ WINDFIELDS.                              !
!                                                                              !
!     HEINZ GUNTHER    GKSS    JANUARY 2001                                    !
!     ERIK MYKLEBUST           NOVEMBER 2004                                   !
!     Arno Behrens     GKSS    February 2007  (DWD wind fields)                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO READ A WIND FIELD AND TRANSFER IT TO THE WAM_WIND_MODULE.           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        FORMATTED READ FROM UNIT IU01, FILE01.                                !
!                                                                              !
!       FILE01, WHICH IS DEFINED IN THE USER INPUT, IS ASSIGNED TO IU01.       !
!       THE FILE MUST STORE:                                                   !
!           1. RECORD: THE WIND DATA HEADER.                                   !
!           2. RECORD: THE WIND DATA TIME.                                     !
!           FOLLOWING RECORDS: THE WIND DATA MATRIX.                           !
!           (RECORD 2 AND FOLLOWING RECORDS ARE REPEATED FOR NEXT WIND FIELD.) !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_WIND_HEADER TO THE WAM_WIND_MODULE:                           !
!          INTEGER :: N_LON      !! NUMBER OF LONGITUDES IN GRID.              !
!          INTEGER :: N_LAT      !! NUMBER OF LATITUDES IN GRID.               !
!          REAL*8  :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].          !
!          REAL*8  :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].         !
!          REAL*8  :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: NORTH      !! NORTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: WEST       !! WEST LONGITUDE OF GRID [DEG].              !
!          REAL*8  :: EAST       !! EAST LONGITUDE OF GRID [DEG].
!          INTEGER :: ICODE      !! WIND CODE: 1= USTAR; 2= USTRESS; 3= U10    !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_WIND_FIELD TO THE WAM_WIND_MODULE:                            !
!          CHARACTER (LEN=14) :: CDTWIR     !! DATE/TIME OF WIND FIELD.        !
!          REAL               :: U_MAP(:,:) !! U COMPONENT OF WIND MAP [M/S].  !
!          REAL               :: V_MAP(:,:) !! V COMPONENT OF WIND MAP [M/S].  !
!                                                                              !
!       THE WINDS MUST BE ON A REGULAR LATITUDE-LONGITUDE GRID ARRANGED        !
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

USE WAM_GENERAL_MODULE,    ONLY: &
&       ABORT1,                  & !! TERMINATES PROCESSING.
&       incdate,                 & !! calculate new date
&       difdate                    !! difference between two dates in seconds

USE WAM_WIND_MODULE,       ONLY: &
&       SET_WIND_HEADER,         & !! SETS WIND HEADER 
&       SET_WIND_FIELD,          & !! SETS WIND FIELD 
&       PRINT_WIND_STATUS          !! PRINTS WIND MODULE STATUS
 
use wam_special_module,    only: &
&       chready                    !! wait for wind/ice files

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST, IU01, FILE01
use wam_timopt_module,  only: ifcst, cdatea, cda
use wam_wind_module,    only: idelwi
use wam_special_module, only: readyf

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER, PARAMETER    :: KIND_D = 8
INTEGER, SAVE         :: ICODE = 3  !! WIND CODE: 1= USTAR; 2= USTRESS; 3= U10
INTEGER, SAVE         :: N_LON      !! NUMBER OF LONGITUDES IN GRID.
INTEGER, SAVE         :: N_LAT      !! NUMBER OF LATITUDES IN GRID.
INTEGER, SAVE         :: ix         !! String position of filename
REAL (KIND=KIND_D)    :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: NORTH      !! NORTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: WEST       !! WEST LONGITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: EAST       !! EAST LONGITUDE OF GRID [DEG].
REAL,    ALLOCATABLE  :: U_MAP(:,:) !! 1. COMPONENT OF WIND MAP [M/S].
REAL,    ALLOCATABLE  :: V_MAP(:,:) !! 2. COMPONENT OF WIND MAP [M/S].
CHARACTER (LEN=14)    :: CDTWIR     !! DATE/TIME OF WIND FIELD

LOGICAL, SAVE  :: FRSTIME = .TRUE.
INTEGER        :: LEN, IOS, J

LOGICAL, SAVE  :: FORMATTED = .TRUE.
character (len=14), save :: chelp

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. FOR FIRST CALL: OPEN FILE AND READ HEADER.                            !
!        ------------------------------------------                            !

if (frstime) then
   chelp = '00000000000000'
   ix=index(file01,'/',.TRUE.)          !! find the last occurrence of ./.
   chelp(1:10)=file01(ix+2:ix+11)       !! extract the filename
else
   call incdate (chelp,idelwi)
   if (readyf) then
      call difdate (cdatea, cda, ifcst)
      ifcst = (ifcst+idelwi)/3600
      if (ifcst<0) ifcst = ifcst+12
      call chready (ifcst)              !! check ready file
   endif
   file01(ix+2:ix+11)=chelp(1:10)       !! rebuild the filename
endif
LEN = LEN_TRIM(FILE01)
OPEN (UNIT=IU01,FILE=FILE01(1:LEN),FORM='FORMATTED',STATUS='OLD',IOSTAT=IOS)
READ (IU01,*, IOSTAT=IOS) SOUTH,NORTH,WEST,EAST,D_LON,D_LAT,N_LON,N_LAT,ICODE
!READ (IU01, '(6F10.5,3I4)', IOSTAT=IOS) SOUTH, NORTH, WEST, EAST,             &
!&                                       D_LON, D_LAT, N_LON, N_LAT, ICODE
IF (IOS.NE.0) THEN
   WRITE (IU06,*) 'FORMATTED READING FROM FILE ', FILE01(1:LEN), ' FAILED'
   WRITE (IU06,*) ' ****************************************************'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_WIND_INPUT        *'
   WRITE (IU06,*) ' *       ===================================        *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' * WIND INPUT FILE COULD NOT BE OPENED              *'
   WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
   WRITE (IU06,*) ' *    FILE NAME IS  FILE01 = ', FILE01(1:LEN)
   WRITE (IU06,*) ' *    UNIT IS         IU01 = ', IU01
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' ****************************************************'
   CALL ABORT1
END IF
   
if (frstime) then
   CALL SET_WIND_HEADER (WEST=WEST,   SOUTH=SOUTH,    &
&                        EAST=EAST,   NORTH=NORTH,    &
&                        D_LON=D_LON, D_LAT=D_LAT,    &
&                        N_LON=N_LON, N_LAT=N_LAT,    &
&                        CODE=ICODE)
   IF (ITEST.GT.0) CALL PRINT_WIND_STATUS
   frstime = .false.
endif

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. ALLOCATE WIND INPUT ARRAYS.                                           !
!        ---------------------------                                           !

IF (.NOT.ALLOCATED(U_MAP) ) ALLOCATE(U_MAP(N_LON,N_LAT))
IF (.NOT.ALLOCATED(V_MAP) ) ALLOCATE(V_MAP(N_LON,N_LAT))

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. READ WIND FIELD.                                                       !
!       -------------------                                                    !

CALL READ_WIND_FORMATTED
CALL SET_WIND_FIELD (CDTWIR, U_MAP, V_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    4. WRITE TEST OUTPUT AND DEALLOCATE ARRAYS.                               !
!       ----------------------------------------                               !

IF (ITEST.GT.1) THEN
   WRITE(IU06,*) ' READ_WIND_INPUT -  WIND FIELD FOR THE CDTWIR = ', CDTWIR
   WRITE(IU06,'(1X,24F6.1)') U_MAP(1:MIN(24,N_LON),1:MIN(5,N_LAT))
   WRITE(IU06,*) ' '
   WRITE(IU06,'(1X,24F6.1)') V_MAP(1:MIN(24,N_LON),1:MIN(5,N_LAT))
END IF

DEALLOCATE(U_MAP)
DEALLOCATE(V_MAP)

close (iu01)
RETURN

CONTAINS

  SUBROUTINE READ_WIND_FORMATTED

!==> DWD - wind fields

    IMPLICIT NONE
    cdtwir = '00000000000000'
    j = 1
    read (iu01,'(3a2,a4)',iostat=ios) cdtwir(9:10), cdtwir(7:8),            &
&                                     cdtwir(5:6), cdtwir(1:4)
    IF (IOS.NE.0) THEN
       WRITE (IU06,*) ' ****************************************************'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_WIND_INPUT        *'
       WRITE (IU06,*) ' *       ===================================        *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' * READ ERROR ON WIND FILE.                         *'
       WRITE (IU06,*) ' * DATE/TIME OF WINDFIELD EXPECTED                  *'
       WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
       WRITE (IU06,*) ' *    FILE NAME IS  FILE01 = ', FILE01(1:LEN)
       WRITE (IU06,*) ' *    UNIT IS         IU01 = ', IU01
       WRITE (IU06,*) ' *    RECORD IS          J = ', J
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' ****************************************************'
       CALL ABORT1
    END IF
    READ(IU01,'(13F6.1)',IOSTAT=IOS) U_MAP
    IF (IOS.NE.0) THEN
       WRITE (IU06,*) ' ****************************************************'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_WIND_INPUT        *'
       WRITE (IU06,*) ' *       ===================================        *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' * READ ERROR ON WIND FILE.                         *'
       WRITE (IU06,*) ' * U - COMPONENTS EXPECTED                          *'
       WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
       WRITE (IU06,*) ' *    FILE NAME IS  FILE01 = ', FILE01(1:LEN)
       WRITE (IU06,*) ' *    UNIT IS         IU01 = ', IU01
       WRITE (IU06,*) ' *    RECORD IS          J = ', J
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' ****************************************************'
       CALL ABORT1
    END IF
    read (iu01,'(1x)')
    READ(IU01,'(13F6.1)',IOSTAT=IOS) V_MAP
    IF (IOS.NE.0) THEN
       WRITE (IU06,*) ' ****************************************************'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_WIND_INPUT        *'
       WRITE (IU06,*) ' *       ===================================        *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' * READ ERROR ON WIND FILE.                         *'
       WRITE (IU06,*) ' * V - COMPONENTS EXPECTED                          *'
       WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
       WRITE (IU06,*) ' *    FILE NAME IS  FILE01 = ', FILE01(1:LEN)
       WRITE (IU06,*) ' *    UNIT IS         IU01 = ', IU01
       WRITE (IU06,*) ' *    RECORD IS          J = ', J
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' ****************************************************'
       CALL ABORT1
    END IF
    
  END SUBROUTINE READ_WIND_FORMATTED

END SUBROUTINE READ_WIND_INPUT

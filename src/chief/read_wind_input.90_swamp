SUBROUTINE READ_WIND_INPUT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READ_WIND_INPUT - ROUTINE TO READ WINDFIELDS.                              !
!                                                                              !
!     HEINZ GUNTHER    GKSS    JANUARY 2001                                    !
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
!          REAL*8  :: EAST       !! EAST LONGITUDE OF GRID [DEG].              !
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

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1                     !! TERMINATES PROCESSING.

USE WAM_WIND_MODULE,       ONLY: &
&       SET_WIND_HEADER,         & !! SETS WIND HEADER 
&       SET_WIND_FIELD,          & !! SETS WIND FIELD 
&       PRINT_WIND_STATUS          !! PRINTS WIND MODULE STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE, ONLY: IU06, ITEST, IU01, FILE01

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !


INTEGER, PARAMETER :: KIND_D = 8

INTEGER, SAVE         :: ICODE = 3  !! WIND CODE: 1= USTAR; 2= USTRESS; 3= U10
INTEGER, SAVE         :: N_LON      !! NUMBER OF LONGITUDES IN GRID.
INTEGER, SAVE         :: N_LAT      !! NUMBER OF LATITUDES IN GRID.
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
INTEGER        :: I
INTEGER        :: IDATEC, ITIMEC, IYYYY, IMM, IDD, IHH, IMI, ISS, IPS, IDAT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. FOR FIRST CALL: OPEN FILE AND READ HEADER.                            !
!        ------------------------------------------                            !

IF (FRSTIME) THEN
   LEN = LEN_TRIM(FILE01)
   OPEN (UNIT=IU01,FILE=FILE01(1:LEN),FORM='FORMATTED',STATUS='OLD',IOSTAT=IOS)
   READ (IU01, '(6F10.5,3I3)', IOSTAT=IOS) SOUTH, NORTH, WEST, EAST,         &
&                                         D_LON, D_LAT, N_LAT, N_LON, ICODE
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) 'FORMATTED READING FROM FILE ', FILE01(1:LEN), ' FAILED'
      WRITE (IU06,*) 'TRYING UNFORMATTED READING'
      CLOSE(UNIT=IU01, STATUS='KEEP')
      OPEN (UNIT=IU01,FILE=FILE01(1:LEN),FORM='UNFORMATTED',STATUS='OLD',IOSTAT=IOS)
      IF (IOS.NE.0) THEN
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
      FORMATTED = .FALSE.

      READ (IU01,IOSTAT=IOS) SOUTH, NORTH, WEST, EAST, D_LON, D_LAT, N_LAT,    &
           &                         N_LON
      IF (IOS.NE.0) THEN
         WRITE (IU06,*) ' ****************************************************'
         WRITE (IU06,*) ' *                                                  *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_WIND_INPUT        *'
         WRITE (IU06,*) ' *       ===================================        *'
         WRITE (IU06,*) ' *                                                  *'
         WRITE (IU06,*) ' * READ ERROR ON WIND FILE.                         *'
         WRITE (IU06,*) ' * FILE HEADER EXPECTED                             *'
         WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
         WRITE (IU06,*) ' *    FILE NAME IS  FILE01 = ', FILE01(1:LEN)
         WRITE (IU06,*) ' *    UNIT IS         IU01 = ', IU01
         WRITE (IU06,*) ' *                                                  *'
         WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
         WRITE (IU06,*) ' *                                                  *'
         WRITE (IU06,*) ' ****************************************************'
         CALL ABORT1
      END IF
   END IF
   
   CALL SET_WIND_HEADER (WEST=WEST,   SOUTH=SOUTH,    &
&                        EAST=EAST,   NORTH=NORTH,    &
&                        D_LON=D_LON, D_LAT=D_LAT,    &
&                        N_LON=N_LON, N_LAT=N_LAT,    &
&                        CODE=ICODE)
   IF (ITEST.GT.0) CALL PRINT_WIND_STATUS
   FRSTIME = .FALSE.
END IF

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

IF (FORMATTED) THEN
   CALL READ_WIND_FORMATTED
ELSE
   CALL READ_WIND_UNFORMATTED
END IF
CALL SET_WIND_FIELD (CDTWIR, U_MAP, V_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    4. WRITE TEST OUTPUT AND DEALLOCATE ARRAYS.                               !
!       ----------------------------------------                               !

IF (ITEST.GT.1) THEN
   WRITE(IU06,*) ' READ_WIND_INPUT -  WIND FIELD FOR THE CDTWIR = ', CDTWIR
   WRITE(IU06,'(1X,24F5.2)') U_MAP(1:MIN(24,N_LON),1:MIN(5,N_LAT))
   WRITE(IU06,*) ' '
   WRITE(IU06,'(1X,24F5.2)') V_MAP(1:MIN(24,N_LON),1:MIN(5,N_LAT))
END IF

DEALLOCATE(U_MAP)
DEALLOCATE(V_MAP)

RETURN

CONTAINS

  SUBROUTINE READ_WIND_FORMATTED

    IMPLICIT NONE
    READ(IU01,'(A14)',IOSTAT=IOS) CDTWIR
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
    
    DO J = 1,N_LAT
       READ(IU01,'(8F10.3)',IOSTAT=IOS) U_MAP(1:N_LON,J)
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
    END DO
    DO J = 1,N_LAT
       READ(IU01,'(8F10.3)',IOSTAT=IOS) V_MAP(1:N_LON,J)
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
    END DO
    
  END SUBROUTINE READ_WIND_FORMATTED


  SUBROUTINE READ_WIND_UNFORMATTED

    IMPLICIT NONE
    READ (IU01, IOSTAT=IOS) IDATEC,ITIMEC
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
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' ****************************************************'
       CALL ABORT1
    END IF
    !
200 FORMAT(I4,5I2.2)
    IPS  = MOD(ITIMEC,100)
    IDAT = ITIMEC/100
    ISS  = MOD(IDAT,100)
    IDAT = IDAT/100
    IMI  = MOD(IDAT,100)
    IHH  = IDAT/100
    IDD  = MOD(IDATEC,100)
    IDAT = IDATEC/100
    IMM  = MOD(IDAT,100)
    IYYYY  = IDAT/100
    WRITE(CDTWIR,200) IYYYY,IMM,IDD,IHH,IMI,ISS
    WRITE(IU06,*)'WIND INPUT DATE = ',CDTWIR
    !
    
    READ(IU01,IOSTAT=IOS) ((U_MAP(I,J),I=1,N_LON),J=1,N_LAT)
    IF (IOS.NE.0) THEN
       WRITE (IU06,*) ' ****************************************************'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_WIND_INPUT        *'
       WRITE (IU06,*) ' *       ===================================        *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' * READ ERROR ON WIND FILE.                         *'
       WRITE (IU06,*) ' * U - COMPONENTS EXPECTED                          *'
       WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT      = ', IOS
       WRITE (IU06,*) ' *    FILE NAME IS  FILE01      = ', FILE01(1:LEN)
       WRITE (IU06,*) ' *    UNIT IS         IU01      = ', IU01
       WRITE (IU06,*) ' *    WIND INPUT DATE IS CDTWIR = ', CDTWIR
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' ****************************************************'
       CALL ABORT1
    END IF
    READ(IU01,IOSTAT=IOS) ((V_MAP(I,J),I=1,N_LON),J=1,N_LAT)
    IF (IOS.NE.0) THEN
       WRITE (IU06,*) ' ****************************************************'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READ_WIND_INPUT        *'
       WRITE (IU06,*) ' *       ===================================        *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' * READ ERROR ON WIND FILE.                         *'
       WRITE (IU06,*) ' * V - COMPONENTS EXPECTED                          *'
       WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT      = ', IOS
       WRITE (IU06,*) ' *    FILE NAME IS  FILE01      = ', FILE01(1:LEN)
       WRITE (IU06,*) ' *    UNIT IS         IU01      = ', IU01
       WRITE (IU06,*) ' *    WIND INPUT DATE IS CDTWIR = ', CDTWIR
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
       WRITE (IU06,*) ' *                                                  *'
       WRITE (IU06,*) ' ****************************************************'
       CALL ABORT1
    END IF
    
  END SUBROUTINE READ_WIND_UNFORMATTED

END SUBROUTINE READ_WIND_INPUT

MODULE WAM_PRINT_USER_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL DEFAULT WAM MODEL SETTINGS, WHICH CAN BE          !
!   CONTROLLED OR WHICH MUST BE DEFINED BY THE USER IN DIFFERNT VERSIONS OF    !
!   USER INPUT.                                                                !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE          !! COORDINATE PROCEDURES

USE WAM_PRINT_MODULE,         ONLY: &
&       SET_OUTPUT_PERIOD,          & !! SET OUTPUT PERIOD. 
&       SET_OUTPUT_TIMES,           & !! SET OUTPUT TIMES.
&       SET_OUTPUT_SITES,           & !! SET OUTPUT SITES FOR SPECTRA.
&       SET_TEST_OPTION,            & !! SETS TEST OPTION.
&       SET_INPUT_FILE,             & !! INTEGRATED DATA FILE (UNFORM. OUTPUT).
&       SET_PARAMETER_OUTPUT_FLAGS, & 
&       SET_SPECTRA_OUTPUT_FLAGS,   & 
&       SET_INTERPOLATION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,         ONLY: IU05, IU06
   
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

! ---------------------------------------------------------------------------- !

CHARACTER (LEN=14) :: START_DATE        !! START DATE OF RUN.
CHARACTER (LEN=14) :: END_DATE          !! END DATE OF RUN.
INTEGER            :: OUTPUT_TIMESTEP
CHARACTER (LEN=1)  :: OUTPUT_TIMESTEP_UNIT

! ---------------------------------------------------------------------------- !

INTEGER, PARAMETER :: MOUTT = 20
CHARACTER (LEN=14), DIMENSION(MOUTT) :: COUTT  !! SPECIFIED OUTPUT TIMES.

! ---------------------------------------------------------------------------- !

INTEGER :: ITEST               !! TEST OUTPUT UP TO LEVEL.

! ---------------------------------------------------------------------------- !

INTEGER            :: INPUT_FILE_UNIT
CHARACTER (LEN=80) :: INPUT_FILE_NAME
CHARACTER (LEN=14) :: INPUT_FILE_DATE
INTEGER            :: INPUT_FILE_TIMESTEP
CHARACTER (LEN=1)  :: INPUT_FILE_TIMESTEP_UNIT

! ---------------------------------------------------------------------------- !

INTEGER, PARAMETER         :: NOUT_P = 70
LOGICAL, DIMENSION(NOUT_P) :: CFLAG_P         !! PARAMETER OUTPUT FLAG.

! ---------------------------------------------------------------------------- !

INTEGER, PARAMETER         :: NOUT_S = 4
LOGICAL, DIMENSION(NOUT_S) :: CFLAG_S         !! SPECTRA OUTPUT FLAG. 

! ---------------------------------------------------------------------------- !

INTEGER, PARAMETER                        :: MOUTP   = 200
CHARACTER(LEN=LEN_COOR), DIMENSION(MOUTP) :: OUTLAT 
CHARACTER(LEN=LEN_COOR), DIMENSION(MOUTP) :: OUTLONG   
CHARACTER (LEN=20),      DIMENSION(MOUTP) :: NAME

! ---------------------------------------------------------------------------- !

LOGICAL :: REGULAR    !! INTERPOLATION OF FIELDS TO REGULAR GRID

! ---------------------------------------------------------------------------- !

NAMELIST /PRINT_NAMELIST/                                                      &
&       START_DATE,                 END_DATE,                                  &
&       ITEST,                                                                 &
&       OUTPUT_TIMESTEP,  OUTPUT_TIMESTEP_UNIT,                                &
&       INPUT_FILE_NAME,  INPUT_FILE_DATE, INPUT_FILE_UNIT,                    &
&       INPUT_FILE_TIMESTEP,  INPUT_FILE_TIMESTEP_UNIT,                        &
&       COUTT,                                                                 &
&       CFLAG_P, CFLAG_S,                                                      &
&       OUTLAT,   OUTLONG,   NAME,                                             &
&       REGULAR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE CLEAR_PRINT_USER_MODULE         !! SETS ALL DATA TO DEFAULT VALUES.
   MODULE PROCEDURE CLEAR_PRINT_USER_MODULE
END INTERFACE
PUBLIC CLEAR_PRINT_USER_MODULE

INTERFACE PRINT_PRINT_NAMELIST           !! PRINTS PRINT NAMELIST.
   MODULE PROCEDURE PRINT_PRINT_NAMELIST
END INTERFACE
PUBLIC PRINT_PRINT_NAMELIST

INTERFACE READ_PRINT_NAMELIST            !! READS PRINT NAMELIST.
   MODULE PROCEDURE READ_PRINT_NAMELIST
END INTERFACE
PUBLIC READ_PRINT_NAMELIST

INTERFACE SET_PRINT_USER_PARAMETER       !! TRANSFERS USER PARAMETERS TO MODULE.
   MODULE PROCEDURE SET_PRINT_USER_PARAMETER
END INTERFACE
PUBLIC SET_PRINT_USER_PARAMETER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E. PRIVATE  INTERFACES.                                                  !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTERFACE CHANGE_TO_SECONDS            !! CHANGES HOURS OR MINUTES TO SECONDS.
   MODULE PROCEDURE CHANGE_TO_SECONDS
END INTERFACE
PRIVATE CHANGE_TO_SECONDS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CLEAR_PRINT_USER_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   CLEAR_PRINT_USER_MODULE - SETS ALL MODULE DATA TO DEFAULT VALUES.          !
!                                                                              !
!       H. GUNTHER   HZG     DECEMBER 2010                                     !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       CLEARS THE PRINT_USER_MODULE.                                          !
!                                                                              !
!     METHOD.                                                                  !
!     --------                                                                 !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. SET ALL DATA TO DEFAULT.                                              !
!        ------------------------                                              !

START_DATE = ' '  !! START DATE OF RUN.
END_DATE   = ' '  !! END   DATE OF RUN.

! ---------------------------------------------------------------------------- !

ITEST                = 0        !! TEST OUTPUT UP TO LEVEL.

! ---------------------------------------------------------------------------- !

OUTPUT_TIMESTEP          = 1
OUTPUT_TIMESTEP_UNIT     = 'H'
INPUT_FILE_NAME          = ' '
INPUT_FILE_DATE          = ' '
INPUT_FILE_TIMESTEP      = 24            
INPUT_FILE_TIMESTEP_UNIT = 'H'            
INPUT_FILE_UNIT          = 20


COUTT = ' '           !! SPECIFIED OUTPUT TIMES.

CFLAG_P     = .TRUE.  !! PARAMETER FILE OUTPUT FLAG.

CFLAG_S     = .TRUE.  !! SPECTRA FILE OUTPUT FLAG.

OUTLAT      = ' '       !! LATITUDES OF OUTPUT SITES.
OUTLONG     = ' '       !! LONGITUDES OF OUTPUT SITES.
NAME        = ' '        !! OUTPUT SITES NAMES.

REGULAR     = .TRUE.  !! INTERPOLATION TO REGULAR GRID

END SUBROUTINE CLEAR_PRINT_USER_MODULE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE READ_PRINT_NAMELIST (FILE, IOS)

! ---------------------------------------------------------------------------- !
!                                                                              !
!       H. GUNTHER   HZG    DECEMBER 2010                                      !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTEGER, INTENT(IN)  :: FILE   !! > 0 READ FROM FILE ASSIGNED TO IU05
                               !! ELSE READ FROM STANDARD INPUT
INTEGER, INTENT(OUT) :: IOS    !! = 0 SUCCESSFULLY READ
                               !! ELSE READ ERRROR
			       
! ---------------------------------------------------------------------------- !

IOS = 0

IF (FILE .GT. 0) THEN
   READ (UNIT=IU05, NML=PRINT_NAMELIST, IOSTAT=IOS)
ELSE
   READ (*, NML=PRINT_NAMELIST,IOSTAT=IOS)
END IF

IF (IOS.EQ.0) THEN
   WRITE (IU06,*) '     SUB. READ_PRINT_NAMELIST SUCCESSFULLY COMPLETED. ' 
END IF

END SUBROUTINE READ_PRINT_NAMELIST

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_PRINT_NAMELIST

WRITE (IU06,*) '  '
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '              PRINT_USER_MODULE STATUS:'
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '  '
WRITE (IU06,NML=PRINT_NAMELIST)
write (iu06,*) '  '

END SUBROUTINE PRINT_PRINT_NAMELIST

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_PRINT_USER_PARAMETER

! ---------------------------------------------------------------------------- !

if (itest>5) then
   CALL PRINT_PRINT_NAMELIST
endif

! ---------------------------------------------------------------------------- !

CALL CHANGE_TO_SECONDS (OUTPUT_TIMESTEP, OUTPUT_TIMESTEP_UNIT)
CALL SET_OUTPUT_PERIOD (B=START_DATE, E=END_DATE, STEP=OUTPUT_TIMESTEP)
CALL SET_OUTPUT_TIMES (TIME=COUTT)

CALL CHANGE_TO_SECONDS (INPUT_FILE_TIMESTEP,INPUT_FILE_TIMESTEP_UNIT)
CALL SET_INPUT_FILE (NAME=INPUT_FILE_NAME, DATE=INPUT_FILE_DATE,               &
&                    STEP=INPUT_FILE_TIMESTEP, UNIT=INPUT_FILE_UNIT) 

CALL SET_TEST_OPTION   (TEST=ITEST)

CALL SET_PARAMETER_OUTPUT_FLAGS (PF=CFLAG_P)
CALL SET_SPECTRA_OUTPUT_FLAGS   (PF=CFLAG_S)

CALL SET_OUTPUT_SITES (LONG=OUTLONG, LAT=OUTLAT, NA=NAME)

CALL SET_INTERPOLATION (REG=REGULAR)

END SUBROUTINE SET_PRINT_USER_PARAMETER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CHANGE_TO_SECONDS (TIME, UNIT) 

INTEGER,           INTENT(INOUT) :: TIME
CHARACTER (LEN=1), INTENT(INOUT) :: UNIT

IF (UNIT.EQ.'M' .OR. UNIT.EQ.'m') TIME = TIME*60.
IF (UNIT.EQ.'H' .OR. UNIT.EQ.'h') TIME = TIME*3600.
UNIT = 'S'

END SUBROUTINE CHANGE_TO_SECONDS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_PRINT_USER_MODULE

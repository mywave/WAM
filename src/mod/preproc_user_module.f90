MODULE PREPROC_USER_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL DEFAULT PREPROC SETTINGS, WHICH CAN BE            !
!   CONTROLLED OR WHICH MUST BE DEFINED BY THE USER IN DIFFERNT VERSIONS OF    !
!   USER INPUT.                                                                !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE          !! COORDINATE PROCEDURES

USE WAM_FRE_DIR_MODULE,  ONLY:  &
&       SET_FRE_DIR                !! DEFINE FREQUENCY DIRECTION GRID.

USE WAM_NEST_MODULE, ONLY:      &
&       SET_BOUNDARY_OPTION,    &  !! DEFINES BOUNDARY OPTION
&       SET_NEST                   !! DEFINES A NEST.

USE PREPROC_MODULE,      ONLY:  &
&       SET_GRID_DEF,           &  !! DEFINE GRID MODEL GRID AREA.
&       SET_GRID_CORRECTIONS,   &  !! TRANSFER DEPTH CORRECTIONS.
&       SET_HEADER                 !! DEFINE MODEL HEADER.

USE WAM_FILE_MODULE,     ONLY:  &
&       SET_TEST_OPTION,        &  !! DEFINE TEST OPTION.
&       SET_TOPO_FILE,          &  !! DEFINE DEPTH DATA FILE NAME.
&       SET_PREPROC_FILE,       &  !! DEFINE PREPROC OUTPUT FILE NAME.
&       SET_C_PREPROC_FILE         !! DEFINE COARSE GRID PREPROC FILE NAME.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,    ONLY: IU05, FILE05, IU06

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
! PRIVATE

! ---------------------------------------------------------------------------- !

CHARACTER (LEN=80) :: HEADER            !! HEADER OF MODEL RUN.
INTEGER            :: ITEST             !! TEST OUTPUT OPTION.

! ---------------------------------------------------------------------------- !
! 
!    2. FREQUENCY-DIRECTION GRID.
!       -------------------------

INTEGER            :: KL                !! NUMBER OF DIRECTIONS.
INTEGER            :: ML                !! NUMBER OF FREQUENCIES.
REAL               :: FR1               !! FIRST FREQUENCY [HZ].

! ---------------------------------------------------------------------------- !
! 
!    3. BASIC MODEL GRID.
!       -----------------

LOGICAL            :: REDUCED_GRID      !! REDUCED GRID OPTION.
INTEGER            :: NX                !! NUMBER OF LONGITUDES.
INTEGER            :: NY                !! NUMBER OF LATITUDES.
CHARACTER(LEN=LEN_COOR) :: XDELLA       !! LONGITUDE INCREMENT.
CHARACTER(LEN=LEN_COOR) :: XDELLO       !! LATITUDE  INCREMENT.
CHARACTER(LEN=LEN_COOR) :: AMOSOP       !! SOUTH LATITUDE.
CHARACTER(LEN=LEN_COOR) :: AMONOP       !! NORTH LATITUDE.
CHARACTER(LEN=LEN_COOR) :: AMOWEP       !! WEST LONGITUDE.
CHARACTER(LEN=LEN_COOR) :: AMOEAP       !! EAST LONGITUDE.
REAL               :: LAND              !! DEPTH >= LAND ARE SEAPOINTS [M].

! ---------------------------------------------------------------------------- !
! 
!    3. CORRECTION AREAS FOR MODEL GRID.
!       --------------------------------

INTEGER, PARAMETER :: NOUT  = 80        !! MAX. NO. OF DEPTH CORRECTION AREAS.
CHARACTER(LEN=LEN_COOR) :: XOUTS(NOUT)       !! S - LATITUDE OF AREA.
CHARACTER(LEN=LEN_COOR) :: XOUTN(NOUT)       !! N - LATITUDE OF AREA.
CHARACTER(LEN=LEN_COOR) :: XOUTW(NOUT)       !! W - LONGITUIDE OF AREA.
CHARACTER(LEN=LEN_COOR) :: XOUTE(NOUT)       !! E - LONGITUIDE OF AREA.
REAL                    :: XOUTD(NOUT)       !! DEPTH OF AREA [M].

! ---------------------------------------------------------------------------- !
! 
!    4. NEST DEFINITIONS.
!       -----------------

INTEGER, PARAMETER :: N_NEST = 20        !! MAX. NO. OF NESTS.
CHARACTER(LEN=LEN_COOR) :: AMOSOC(N_NEST)     !! SOUTH LATITUDE OF NEST.
CHARACTER(LEN=LEN_COOR) :: AMONOC(N_NEST)     !! NORTH LATITUDE OF NEST.
CHARACTER(LEN=LEN_COOR) :: AMOWEC(N_NEST)     !! WEST LONGITUDE OF NEST.
CHARACTER(LEN=LEN_COOR) :: AMOEAC(N_NEST)     !! EAST LONGITUDE OF NEST.
CHARACTER*20       :: NEST_NAME(N_NEST)  !! NAME OF NEST.
integer            :: nestcode(n_nest)   !! boundary values in ascii or binary
INTEGER            :: PREPROC_C_INPUT_FILE_UNIT
CHARACTER (LEN=80) :: PREPROC_C_INPUT_FILE_NAME 

! ---------------------------------------------------------------------------- !
! 
!    5. TOPOGRAPHY FILE.
!       ----------------

INTEGER            :: TOPO_INPUT_FILE_UNIT
CHARACTER (LEN=80) :: TOPO_INPUT_FILE_NAME

! ---------------------------------------------------------------------------- !
! 
!    6. PREPROC OUTPUT FILE.
!       --------------------

INTEGER            :: PREPROC_OUTPUT_FILE_UNIT
CHARACTER (LEN=80) :: PREPROC_OUTPUT_FILE_NAME

! ---------------------------------------------------------------------------- !
! 
!    7. PREPROC NAMELIST.
!       -----------------

NAMELIST /PREPROC_NAMELIST/                                                    &
&       HEADER,   ITEST,                                                       &      
&       KL,       ML,      FR1,                                                &
&       REDUCED_GRID,                                                          &
&       NX,       NY,      XDELLA,  XDELLO,                                    &
&       AMOSOP,   AMONOP,  AMOWEP,  AMOEAP,                                    &
&       LAND,                                                                  &
&       XOUTS,    XOUTN,   XOUTW,   XOUTE,   XOUTD,                            &
&       AMOSOC,   AMONOC,  AMOWEC,  AMOEAC,  NEST_NAME, nestcode,              &
&       PREPROC_C_INPUT_FILE_UNIT,  PREPROC_C_INPUT_FILE_NAME,                 & 
&       TOPO_INPUT_FILE_UNIT,       TOPO_INPUT_FILE_NAME,                      &
&       PREPROC_OUTPUT_FILE_UNIT,   PREPROC_OUTPUT_FILE_NAME 

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE CLEAR_PREPROC_USER_MODULE        !! SETS ALL DATA TO DEFAULT VALUES.
   MODULE PROCEDURE CLEAR_PREPROC_USER_MODULE
END INTERFACE
PUBLIC CLEAR_PREPROC_USER_MODULE

INTERFACE PRINT_PREPROC_NAMELIST           !! PRINTS PREPROC NAMELIST.
   MODULE PROCEDURE PRINT_PREPROC_NAMELIST
END INTERFACE
PUBLIC PRINT_PREPROC_NAMELIST

INTERFACE READ_PREPROC_NAMELIST            !! READS PREPROC NAMELIST.
   MODULE PROCEDURE READ_PREPROC_NAMELIST
END INTERFACE
PUBLIC READ_PREPROC_NAMELIST

INTERFACE SET_PREPROC_USER_PARAMETER       !! TRANSFERS ALL DATA TO MODULES.
   MODULE PROCEDURE SET_PREPROC_USER_PARAMETER
END INTERFACE
PUBLIC SET_PREPROC_USER_PARAMETER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CLEAR_PREPROC_USER_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   CLEAR_PREPROC_USER_MODULE - SETS ALL MODULE DATA TO DEFAULT VALUES.        !
!                                                                              !
!       H. GUNTHER   GKSS    DECEMBER 2009                                     !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       CLEARS THE PREPROC_USER_MODULE.                                        !
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

HEADER = ' '         !! HEADER OF MODEL RUN.
ITEST  = 0           !! TEST OUTPUT OPTION.

! ---------------------------------------------------------------------------- !

KL     = 24          !! NUMBER OF DIRECTIONS.
ML     = 25          !! NUMBER OF FREQUENCIES.
FR1    = .04177248   !! FIRST FREQUENCY [HZ].

! ---------------------------------------------------------------------------- !

REDUCED_GRID = .FALSE.   !! REDUCED GRID OPTION.
NX     = -1              !! NUMBER OF LONGITUDES.
NY     = -1              !! NUMBER OF LATITUDES.
XDELLA = ' '             !! LONGITUDE INCREMENT.
XDELLO = ' '             !! LATITUDE  INCREMENT.
AMOSOP = ' '             !! SOUTH LATITUDE.
AMONOP = ' '             !! NORTH LATITUDE.
AMOWEP = ' '             !! WEST LONGITUDE.
AMOEAP = ' '             !! EAST LONGITUDE.
LAND    = 0.             !! DEPTH >= LAND ARE SEAPOINTS [M].

! ---------------------------------------------------------------------------- !

XOUTS = ' '        !! S - LATITUDE OF AREA [DEG].
XOUTN = ' '        !! N - LATITUDE OF AREA [DEG].
XOUTW = ' '        !! W - LONGITUIDE OF AREA [DEG].
XOUTE = ' '        !! E - LONGITUIDE OF AREA [DEG].
XOUTD = -999.      !! DEPTH IN AREA [M].

! ---------------------------------------------------------------------------- !

AMOSOC = ' '       !! SOUTH LATITUDE OF NESTS [DEG].
AMONOC = ' '       !! NORTH LATITUDE OF NESTS [DEG].
AMOWEC = ' '       !! WEST LONGITUDE OF NESTS [DEG].
AMOEAC = ' '       !! EAST LONGITUDE OF NESTS [DEG].
NEST_NAME = ' '    !! NAME OF NESTS.
nestcode  = 0      !! boundary values in ascii or binary
PREPROC_C_INPUT_FILE_UNIT = 10
PREPROC_C_INPUT_FILE_NAME = ' ' 

! ---------------------------------------------------------------------------- !

TOPO_INPUT_FILE_UNIT = 1
TOPO_INPUT_FILE_NAME = ' '

! ---------------------------------------------------------------------------- !

PREPROC_OUTPUT_FILE_UNIT = 7
PREPROC_OUTPUT_FILE_NAME = 'Grid_Info'

END SUBROUTINE CLEAR_PREPROC_USER_MODULE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE READ_PREPROC_NAMELIST (FILE, IOS)


! ---------------------------------------------------------------------------- !
!                                                                              !
!    Arno Behrens    MSC/GKSS    January 2004                                  !
!                    module for namelist management, namelist read shifted     !
!                    to chief for message passing purposes,                    !
!                    namelist replaces the common user input file              !
!                                                                              !
!    Erik Myklebust             November 2004                                  !
!                    reading of namelist can now optionally be from file       !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTEGER, INTENT(IN)  :: FILE   !! > 0 READ FROM FILE ASSIGNED TO IU05
                               !! ELSE READ FROM STANDARD INPUT
INTEGER, INTENT(OUT) :: IOS    !! = 0 SUCCESSFULLY READ
                               !! ELSE READ ERRROR
			       
! ---------------------------------------------------------------------------- !

IOS = 0

IF (FILE .GT. 0) THEN
   READ (UNIT=IU05, NML=PREPROC_NAMELIST, IOSTAT=IOS)
ELSE
   READ (*, NML=PREPROC_NAMELIST,IOSTAT=IOS)
END IF

IF (IOS.EQ.0) THEN
   WRITE (IU06,*) '     SUB. READ_PREPROC_NAMELIST SUCCESSFULLY COMPLETED. ' 
END IF

END SUBROUTINE READ_PREPROC_NAMELIST

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_PREPROC_NAMELIST

WRITE (IU06,*) '  '
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '             PREPROC_USER_MODULE STATUS:'
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '  '
WRITE (IU06,NML=PREPROC_NAMELIST)

END SUBROUTINE PRINT_PREPROC_NAMELIST

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_PREPROC_USER_PARAMETER

! ---------------------------------------------------------------------------- !

if (itest>5) then
   CALL PRINT_PREPROC_NAMELIST
endif

! ---------------------------------------------------------------------------- !

CALL SET_HEADER (HEADER)
CALL SET_TEST_OPTION (TEST=ITEST)

CALL SET_FRE_DIR (N_DIR=KL, N_FRE=ML, FR1=FR1)
  
CALL SET_GRID_DEF (N_LON=NX, N_LAT=NY, D_LON=XDELLO, D_LAT=XDELLA,             &
&                 SOUTH=AMOSOP, NORTH=AMONOP, WEST=AMOWEP, EAST=AMOEAP,        &
&                  LAND=LAND, R_GRID=REDUCED_GRID)

CALL SET_GRID_CORRECTIONS (SOUTH=XOUTS, NORTH=XOUTN, WEST=XOUTW, EAST=XOUTE,   &
&                          D_COR=XOUTD)

CALL SET_NEST (SOUTH=AMOSOC, NORTH=AMONOC, WEST=AMOWEC, EAST=AMOEAC,           &
&              NAME=NEST_NAME, ncode=nestcode)

CALL SET_C_PREPROC_FILE (NAME=PREPROC_C_INPUT_FILE_NAME,                       &
&                        UNIT=PREPROC_C_INPUT_FILE_UNIT) 

CALL SET_BOUNDARY_OPTION (C=COUNT(AMOSOC.NE.' ').GT.0,                         &
&                         F=PREPROC_C_INPUT_FILE_NAME .NE. ' ')

! ---------------------------------------------------------------------------- !

CALL SET_TOPO_FILE (NAME=TOPO_INPUT_FILE_NAME,                                 &
&                   UNIT=TOPO_INPUT_FILE_UNIT) 

! ---------------------------------------------------------------------------- !

CALL SET_PREPROC_FILE (NAME=PREPROC_OUTPUT_FILE_NAME,                          &
&                      UNIT=PREPROC_OUTPUT_FILE_UNIT) 

END SUBROUTINE SET_PREPROC_USER_PARAMETER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE PREPROC_USER_MODULE

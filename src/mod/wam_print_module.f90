MODULE WAM_PRINT_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS: OUTPUT TIMES, FLAGS, AND ARRAYS NECESSARY FOR GRIDDED!
!                         FIELDS OF PARAMTERS.                                 !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE          !! COORDINATE PROCEDURES

USE WAM_GENERAL_MODULE, ONLY:  &   !! TERMINATES PROCESSING. 
&        ABORT1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,         ONLY: IU06, ITEST
   
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.  NEXT OUTPUT TIMES AND TIME INCREMENTS.                               !
!        ---------------------------------------                               !

CHARACTER (LEN=14) :: CDATEA    !! START DATE OF PRINT OUPUT  (YYMMDDHHMM).
CHARACTER (LEN=14) :: CDATEE    !! END DATE OF PRINT OUPUT (YYMMDDHHMM).
INTEGER            :: IDELDO    !! PRINT OUTPUT TIMESTEP IN SECONDS.
CHARACTER (LEN=14) :: CDTINTT   !! DATE OF DATA IN MODULE (YYMMDDHHMM).


! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. FIXED OUTPUT TIMES.                                                   !
!        -------------------                                                   !

INTEGER            :: NOUTT = 0
CHARACTER (LEN=14), DIMENSION(:), ALLOCATABLE :: COUTT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. OUTPUT POINTS.                                                        !
!        --------------                                                        !

INTEGER                            :: NOUTP = 0 !! NUMBER OF OUTPUT POINTS.
INTEGER, DIMENSION(:), ALLOCATABLE :: OUTLAT    !! LATITUDE OF POINTS [M_SEC].
INTEGER, DIMENSION(:), ALLOCATABLE :: OUTLONG   !! LONGITUDE OF POINTS [M_SEC].
CHARACTER (LEN=20), DIMENSION(:), ALLOCATABLE :: NAME(:)  !! OUTPUT SITE NAMES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. PARAMETER FLAGS, TITLE AND SCALING FACTORS.                           !
!        -------------------------------------------                           !

INTEGER, PARAMETER         :: NOUT_P = 40
LOGICAL, DIMENSION(NOUT_P) :: CFLAG_P         !! FILE PARAMETER OUTPUT FLAG.
LOGICAL, DIMENSION(NOUT_P) :: PFLAG_P         !! FLAG OF DATA IN MODULE.

CHARACTER(LEN=60), DIMENSION(NOUT_P) :: TITL_P = (/                  &
& ' WIND SPEED U10 ( METRES/SECOND )                           ',    &   !!  1
& ' WIND DIRECTION ( DEGREE FROM NORTH TO )                    ',    &   !!  2
& ' FRICTION VELOCITY ( METRES/SECOND )                        ',    &   !!  3
& ' DRAG COEFFICIENT ( PROMILLE )                              ',    &   !!  4
& ' CHARNOCK PARAMETER                                         ',    &   !!  5
& ' WATER DEPTH (METRES) (DEEPER THAN 999M ARE PRINTED AS 999) ',    &   !!  6
& ' CURRENT SPEED ( METRES/SECOND )                            ',    &   !!  7
& ' CURRENT DIRECTION ( DEGREE FROM NORTH TO )                 ',    &   !!  8
& ' SIGNIFICANT WAVE HEIGHT ( METRES )                         ',    &   !!  9
& ' WAVE PEAK PERIOD ( SECONDS )                               ',    &   !! 10
& ' WAVE MEAN PERIOD (SECONDS )                                ',    &   !! 11
& ' WAVE TM1 PERIOD ( SECONDS )                                ',    &   !! 12
& ' WAVE TM2 PERIOD ( SECONDS )                                ',    &   !! 13
& ' WAVE DIRECTION ( DEGREE FROM NORTH TO )                    ',    &   !! 14
& ' DIRECTIONAL SPREAD ( DEGREES )                             ',    &   !! 15
& ' NORMALISED WAVE STRESS ( % )                               ',    &   !! 16
& ' SEA SIGNIFICANT WAVE HEIGHT ( METRES )                     ',    &   !! 17
& ' SEA PEAK PERIOD ( SECONDS )                                ',    &   !! 18
& ' SEA MEAN PERIOD ( SECONDS )                                ',    &   !! 19
& ' SEA TM1 PERIOD ( SECONDS )                                 ',    &   !! 20
& ' SEA TM2 PERIOD (  SECONDS )                                ',    &   !! 21
& ' SEA DIRECTION ( DEGREE FROM NORTH TO )                     ',    &   !! 22
& ' SEA DIRECTIONAL SPREAD ( DEGREES )                         ',    &   !! 23
& ' DUMMY                                                      ',    &   !! 24
& ' SWELL SIGNIFICANT WAVE HEIGHT ( METRES )                   ',    &   !! 25
& ' SWELL PEAK PERIOD ( SECONDS )                              ',    &   !! 26
& ' SWELL MEAN PERIOD ( SECONDS )                              ',    &   !! 27
& ' SWELL TM1 PERIOD ( SECONDS )                               ',    &   !! 28
& ' SWELL TM2 PERIOD ( SECONDS )                               ',    &   !! 29
& ' SWELL DIRECTION ( DEGREE FROM NORTH TO )                   ',    &   !! 30
& ' SWELL DIRECTIONAL SPREAD ( DEGREES )                       ',    &   !! 31
& ' DUMMY                                                      ',    &   !! 32
& ' GODA PEAKEDNESS PARAMETER                                  ',    &   !! 33
& ' KURTOSIS                                                   ',    &   !! 34
& ' BENJAMIN-FEIR INDEX                                        ',    &   !! 35
& ' NORMALIZED MAXIMUM WAVE HEIGHT                             ',    &   !! 36
& ' MAXIMUM WAVE PERIOD ( SECONDS )                            ',    &   !! 37
& ' PEAK FREQUENCY (INTERPOLATED) ( HZ )                       ',    &   !! 38
& ' PEAK DIRECTION ( DEGREE FROM NORTH TO )                    ',    &   !! 39
& ' MEAN SQUARE SLOPE                                          '/)       !! 40

REAL, PARAMETER, DIMENSION(NOUT_P) :: SCAL_P = (/                              &
&                     10.            ,    &   !!  1
&                      1.            ,    &   !!  2
&                    100.            ,    &   !!  3
&                  10000.            ,    &   !!  4
&                  10000.            ,    &   !!  5
&                      1.            ,    &   !!  6
&                    100.            ,    &   !!  7
&                      1.            ,    &   !!  8
&                     10.            ,    &   !!  9
&                     10.            ,    &   !! 10
&                     10.            ,    &   !! 11
&                     10.            ,    &   !! 12
&                     10.            ,    &   !! 13
&                      1.            ,    &   !! 14
&                      1.            ,    &   !! 15
&                    100.            ,    &   !! 16
&                     10.            ,    &   !! 17
&                     10.            ,    &   !! 18
&                     10.            ,    &   !! 19
&                     10.            ,    &   !! 20
&                     10.            ,    &   !! 21
&                      1.            ,    &   !! 22
&                      1.            ,    &   !! 23
&                      1.            ,    &   !! 24
&                     10.            ,    &   !! 25
&                     10.            ,    &   !! 26
&                     10.            ,    &   !! 27
&                     10.            ,    &   !! 28
&                     10.            ,    &   !! 29
&                      1.            ,    &   !! 30
&                      1.            ,    &   !! 31
&                      1.            ,    &   !! 32
&                     10.            ,    &   !! 33
&                    100.            ,    &   !! 34
&                     10.            ,    &   !! 35
&                     10.            ,    &   !! 36
&                     10.            ,    &   !! 37
&                   1000.            ,    &   !! 38
&                      1.            ,    &   !! 39
&                   1000.            /)       !! 40

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. SPECTRA FLAGS AND TITLE.                                              !
!        ------------------------                                              !

INTEGER, PARAMETER         :: NOUT_S = 4
LOGICAL, DIMENSION(NOUT_S) :: CFLAG_S         !! SPECTRA OUTPUT FLAG.
LOGICAL, DIMENSION(NOUT_S) :: PFLAG_S         !! FLAG OF DATA IN MODULE.

CHARACTER(LEN=60), DIMENSION(NOUT_S) :: TITL_S = (/                  &
& ' SPECTRUM                                                   ',    &   !!  1
& ' SEA SPECTRUM                                               ',    &   !!  2
& ' SWELL SPECTRUM                                             ',    &   !!  3
& ' DUMMY                                                      '/)       !!  4

! ---------------------------------------------------------------------------- !
!                                                                              !
!     8. RADIATION FLAGS, TITLE AND SCALING FACTORS.                           !
!        -------------------------------------------                           !

INTEGER, PARAMETER         :: NOUT_R = 8
LOGICAL, DIMENSION(NOUT_R) :: CFLAG_R         !! SPECTRA OUTPUT FLAG.
LOGICAL, DIMENSION(NOUT_R) :: PFLAG_R         !! FLAG OF DATA IN MODULE.

CHARACTER(LEN=60), DIMENSION(NOUT_R) :: TITL_R = (/                  &
& ' RADIATION STRESS TENSOR SXX ( KG/S/S )                     ',    &   !!  1
& ' RADIATION STRESS TENSOR SYY ( KG/S/S )                     ',    &   !!  2
& ' RADIATION STRESS TENSOR SXY ( KG/S/S )                     ',    &   !!  3
& ' DUMMY                                                      ',    &   !!  4
& ' X-COMP. WAVE FORCE PER SURFACE UNIT ( N/M/M )              ',    &   !!  5
& ' Y-COMP. WAVE FORCE PER SURFACE UNIT ( N/M/M )              ',    &   !!  6
& ' X-COMP. STOKES DRIFT ( M/S )                               ',    &   !!  7
& ' Y-COMP. STOKES DRIFT ( M/S )                               '/)       !!  8

REAL, PARAMETER, DIMENSION(NOUT_R) :: SCAL_R = (/                    &
&                      0.            ,    &   !!  1
&                      0.            ,    &   !!  2
&                      0.            ,    &   !!  3
&                      1.            ,    &   !!  4
&                      0.            ,    &   !!  5
&                      0.            ,    &   !!  6
&                      0.            ,    &   !!  7
&                      0.            /)       !!  8

! ---------------------------------------------------------------------------- !
!                                                                              !
!     9. INTERPOLATION OPTION.                                                 !
!        ---------------------                                                 !

LOGICAL  :: REGULAR  !! IF TRUE, REDUCED GRIDS ARE INTERPOLATED TO REGULAR ONES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!    10. FILE NAME INFORMATION.                                                !
!        ----------------------                                                !

INTEGER            :: IU01 = 1          !! INPUT UNIT FOR INTEGRATED PARAMETER.
CHARACTER (LEN=80) :: FILE01 = 'MAP'    !! FILE IDENTIFIER.
CHARACTER (LEN=14) :: CDTFILE           !! DATE OF FILE.
INTEGER            :: IDFILE            !! FILE DATE INCREMENT IN SECONDS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!   11. GENERAL GRID INFORMATION.                                             !
!        -------------------------                                             !

INTEGER              :: NX         !! NUMBER OF LONGITUDES IN GRID.
INTEGER              :: NY         !! NUMBER OF LATITUDES  IN GRID.
INTEGER, ALLOCATABLE :: NLON_RG(:) !! GRID INCREMENT FOR LONGITUDES [M_SEC].
LOGICAL              :: PER        !! = .TRUE. IF GRID IS PERIODIC.
INTEGER              :: AMOWEP     !! MOST WESTERN LONGITUDE IN GRID [M_SEC].
INTEGER              :: AMOSOP     !! MOST SOUTHERN LATITUDE IN GRID [M_SEC].
INTEGER              :: AMOEAP     !! MOST EASTERN LONGITUDE IN GRID [M_SEC].
INTEGER              :: AMONOP     !! MOST NORTHERN LATITUDE IN GRID [M_SEC].
INTEGER              :: XDELLA     !! GRID INCREMENT FOR LATITUDE [[M_SEC].
INTEGER              :: XDELLO     !! GRID INCREMENT FOR LONGITUDE AT EQUATOR [M_SEC].
INTEGER, ALLOCATABLE :: ZDELLO(:)  !! GRID INCREMENT FOR LONGITUDES [M_SEC].

! ---------------------------------------------------------------------------- !
!                                                                              !
!    12. GRIDDED INTEGRATED PARAMETER.                                         !
!        -----------------------------                                         !

REAL, DIMENSION(:,:,:), ALLOCATABLE :: GRID !! GRIDDED MODEL OUTPUT PARAMETERS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!    13. GENERAL SPECTRA INFORMATION.                                          !
!        ----------------------------                                          !

INTEGER                          :: KL        !! NUMBER OF DIRECTIONS.
INTEGER                          :: ML        !! NUMBER OF FREQUENCIES.
REAL                             :: CO        !! LOGARTHMIC FREQUENCY INCREMENT.
REAL, DIMENSION(:), ALLOCATABLE  :: FR        !! FREQUENCIES [HZ].
REAL, DIMENSION(:), ALLOCATABLE  :: THETA     !! DIRECTIONS [DEG].
INTEGER                          :: SPEC_LAT  !! LATITUDE OF SPECTRUM [M_SEC].
INTEGER                          :: SPEC_LON  !! LONGITUDE ODF SPECTRUM [M_SEC].
CHARACTER (LEN=14)               :: SPEC_DATE !! DATE OF SPECTRUM.

REAL, DIMENSION(:,:), ALLOCATABLE :: SPEC  !! SPECTRUM [M*M/HZ].
REAL               :: U10         !! WIND SPEED U10 [M/S].
REAL               :: UDIR        !! WIND DIRECTION [DEG].
REAL               :: US          !! FRICTION VELOCITY (M/S).
REAL               :: DEPTH       !! WATER DEPTH [M].
REAL               :: CSPEED      !! CURRENT SPEED [M/S].
REAL               :: CDIR        !! CURRENT DIRECTION [M/S].
REAL               :: HS          !! SIG. WAVE HEIGHT [M].
REAL               :: PPER        !! PEAK PERIOD [S].
REAL               :: MPER        !! MEAN PERIOD [S].
REAL               :: TM1         !! TM1 PERIOD [S].
REAL               :: TM2         !! TM2 PERIOD [S].
REAL               :: MDIR        !! MEAN DIRECTION [DEG].
REAL               :: SPRE        !! MEAN SPREAD [DEG].
REAL               :: TAUW        !! NORMALISED WAVE STRESS [%].

REAL, DIMENSION(:,:), ALLOCATABLE :: SPEC_SEA  !! SEA SPECTRUM [M*M/HZ].
REAL               :: HS_SEA      !! SEA  SIG. WAVE HEIGHT [M].
REAL               :: PPER_SEA    !! SEA PEAK PERIOD [S].
REAL               :: MPER_SEA    !! SEA MEAN PERIOD [S].
REAL               :: TM1_SEA     !! SEA TM1 PERIOD [S].
REAL               :: TM2_SEA     !! SEA TM2 PERIOD [S].
REAL               :: MDIR_SEA    !! SEA MEAN DIRECTION [DEG].
REAL               :: SPRE_SEA    !! SEA MEAN SPREAD [DEG].

REAL, DIMENSION(:,:), ALLOCATABLE :: SPEC_SWELL  !! SWELL SPECTRUM [M*M/HZ].
REAL               :: HS_SWELL    !! SWELL SIG. WAVE HEIGHT [M].
REAL               :: PPER_SWELL  !! SWELL PEAK PERIOD [S].
REAL               :: MPER_SWELL  !! SWELL MEAN PERIOD [S].
REAL               :: TM1_SWELL   !! SWELL TM1 PERIOD [S].
REAL               :: TM2_SWELL   !! SWELL TM2 PERIOD [S].
REAL               :: MDIR_SWELL  !! SWELL MEAN DIRECTION [DEG].
REAL               :: SPRE_SWELL  !! SWELL MEAN SPREAD [DEG].

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE SET_TEST_OPTION
   MODULE PROCEDURE SET_TEST_OPTION
END INTERFACE
PUBLIC SET_TEST_OPTION

INTERFACE SET_INPUT_FILE
   MODULE PROCEDURE SET_INPUT_FILE
END INTERFACE
PUBLIC SET_INPUT_FILE

INTERFACE SET_INTERPOLATION
   MODULE PROCEDURE SET_INTERPOLATION
END INTERFACE
PUBLIC SET_INTERPOLATION

INTERFACE SET_OUTPUT_PERIOD
   MODULE PROCEDURE SET_OUTPUT_PERIOD
END INTERFACE
PUBLIC SET_OUTPUT_PERIOD

INTERFACE SET_OUTPUT_SITES
   MODULE PROCEDURE SET_OUTPUT_SITES_C    !! CHARACTER VERSION
   MODULE PROCEDURE SET_OUTPUT_SITES_D    !! DEGREE VERSION
   MODULE PROCEDURE SET_OUTPUT_SITES_M    !! M_SEC VERSION
END INTERFACE
PUBLIC SET_OUTPUT_SITES

INTERFACE SET_OUTPUT_TIMES
   MODULE PROCEDURE SET_OUTPUT_TIMES
END INTERFACE
PUBLIC SET_OUTPUT_TIMES

INTERFACE SET_PARAMETER_OUTPUT_FLAGS
   MODULE PROCEDURE SET_PARAMETER_OUTPUT_FLAGS
END INTERFACE
PUBLIC SET_PARAMETER_OUTPUT_FLAGS

INTERFACE SET_RADIATION_OUTPUT_FLAGS
   MODULE PROCEDURE SET_RADIATION_OUTPUT_FLAGS
END INTERFACE
PUBLIC SET_RADIATION_OUTPUT_FLAGS

INTERFACE SET_SPECTRA_OUTPUT_FLAGS
   MODULE PROCEDURE SET_SPECTRA_OUTPUT_FLAGS
END INTERFACE
PUBLIC SET_SPECTRA_OUTPUT_FLAGS

INTERFACE PRINT_GRID_USER
   MODULE PROCEDURE PRINT_GRID_USER
END INTERFACE
PUBLIC PRINT_GRID_USER

INTERFACE PRINT_TIME_USER
   MODULE PROCEDURE PRINT_TIME_USER
END INTERFACE
PUBLIC PRINT_TIME_USER

INTERFACE PRINT_RADIATION_USER
   MODULE PROCEDURE PRINT_RADIATION_USER
END INTERFACE
PUBLIC PRINT_RADIATION_USER

INTERFACE PRINT_SPECTRA_USER
   MODULE PROCEDURE PRINT_SPECTRA_USER
END INTERFACE
PUBLIC PRINT_SPECTRA_USER

INTERFACE PREPARE_EXTRACTION
   MODULE PROCEDURE PREPARE_EXTRACTION
END INTERFACE
PUBLIC PREPARE_EXTRACTION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_TEST_OPTION (TEST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,   INTENT(IN) :: TEST      !! TEST LEVEL

! ---------------------------------------------------------------------------- !

ITEST = MAX(TEST,0)

END SUBROUTINE SET_TEST_OPTION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_INPUT_FILE (NAME, DATE, STEP, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=*),  INTENT(IN)    :: NAME  !! FILE NAME
CHARACTER (LEN=14), INTENT(IN)    :: DATE  !! FILE DATE
INTEGER,            INTENT(IN)    :: STEP  !! FILE TIMESTEP (SEC)
INTEGER,            INTENT(IN)    :: UNIT  !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

INTEGER  :: LENT

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE01 = NAME(1:LENT)
CDTFILE = DATE
IDFILE = STEP
IU01 = UNIT

END SUBROUTINE SET_INPUT_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_INTERPOLATION (REG)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

LOGICAL,  INTENT(IN)    :: REG  !! FILE NAME

! ---------------------------------------------------------------------------- !

REGULAR = REG

END SUBROUTINE SET_INTERPOLATION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_OUTPUT_PERIOD (B, E, STEP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=14), INTENT(IN)   :: B     !! START DATE OF OUTPUT.
CHARACTER (LEN=14), INTENT(IN)   :: E     !! END DATE OF OUTPUT.
INTEGER,            INTENT(IN)   :: STEP  !! OUTPUT TIMESTEP (SEC)

! ---------------------------------------------------------------------------- !

CDATEA = B
CDATEE = E
IDELDO = STEP

IF (CDATEE.LT.CDATEA) THEN
   WRITE(IU06,*) '*********************************************'
   WRITE(IU06,*) '*                                           *'
   WRITE(IU06,*) '*   FATAL ERROR IN SUB. SET_OUTPUT_PERIOD   *'
   WRITE(IU06,*) '*   =====================================   *'
   WRITE(IU06,*) '*                                           *'
   WRITE(IU06,*) '* END DATE IS BEFORE START DATE             *'
   WRITE(IU06,*) '* START DATE = ', CDATEA
   WRITE(IU06,*) '* END  DATE  = ', CDATEE
   WRITE(IU06,*) '*                                           *'
   WRITE(IU06,*) '* CORRECT USER INPUT                        *'
   WRITE(IU06,*) '*                                           *'
   WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.         *'
   WRITE(IU06,*) '*********************************************'
   CALL ABORT1
END IF

END SUBROUTINE SET_OUTPUT_PERIOD

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_OUTPUT_TIMES (TIME)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER*14, INTENT(IN) :: TIME(:) !! OUTPUT TIMES.

! ---------------------------------------------------------------------------- !

NOUTT = COUNT(TIME.NE.' ')
IF (NOUTT.GT.0) THEN
   ALLOCATE (COUTT(NOUTT))
   COUTT = TIME(1:NOUTT)
END IF

END SUBROUTINE SET_OUTPUT_TIMES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_OUTPUT_SITES_C (LONG, LAT, NA)

CHARACTER (LEN=LEN_COOR), INTENT(IN) :: LONG(:) !! LONGITUDES.
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: LAT(:)  !! LATITUDES.
CHARACTER (LEN=*), INTENT(IN) :: NA(:)          !! SITE NAME.

! ---------------------------------------------------------------------------- !
! 
!    1. RE-FORMAT INPUT PARAMETERS AND CALL SET_OUTPUT_SITES_M.
!       -------------------------------------------------------

CALL SET_OUTPUT_SITES_M ( READ_COOR_TEXT(LONG),  READ_COOR_TEXT(LAT), NA)

END SUBROUTINE SET_OUTPUT_SITES_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_OUTPUT_SITES_D (LONG, LAT, NA)

REAL (KIND=KIND_D), INTENT(IN) :: LONG(:) !! LONGITUDES.
REAL (KIND=KIND_D), INTENT(IN) :: LAT(:)  !! LATITUDES.
CHARACTER (LEN=*),  INTENT(IN) :: NA(:)   !! SITE NAME.

! ---------------------------------------------------------------------------- !
! 
!    1. RE-FORMAT INPUT PARAMETERS AND CALL SET_OUTPUT_SITES_M.
!       -------------------------------------------------------

CALL SET_OUTPUT_SITES_M (DEG_TO_M_SEC(LONG), DEG_TO_M_SEC(LAT), NA)

END SUBROUTINE SET_OUTPUT_SITES_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_OUTPUT_SITES_M (LONG, LAT, NA)

INTEGER,           INTENT(IN) :: LONG(:) !! LONGITUDES.
INTEGER,           INTENT(IN) :: LAT(:)  !! LATITUDES.
CHARACTER (LEN=*), INTENT(IN) :: NA(:)   !! SITE NAME.

NOUTP = COUNT (LONG.NE.COOR_UNDEF .AND. LAT.NE.COOR_UNDEF)
IF (NOUTP.GT.0) THEN
   ALLOCATE (OUTLONG(NOUTP))
   ALLOCATE (OUTLAT (NOUTP))
   ALLOCATE (NAME   (NOUTP))
   OUTLONG(1:NOUTP) = LONG(1:NOUTP)
   OUTLAT (1:NOUTP) = LAT (1:NOUTP)
   NAME   (1:NOUTP) = NA  (1:NOUTP)
END IF

END SUBROUTINE SET_OUTPUT_SITES_M

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_PARAMETER_OUTPUT_FLAGS (PF)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

LOGICAL, INTENT(IN) :: PF(:)   !! PRINTER FLAGS.

! ---------------------------------------------------------------------------- !

IF (SIZE(PF).NE.NOUT_P) THEN
   WRITE(IU06,*) '*  PROGRAM NEEDS ', NOUT_P,' FLAGS FOR PARAMETER OUTPUT *'
   WRITE(IU06,*) '*  NUMBER OF PRINTER FLAGS IS : ', SIZE(PF)
   CALL ABORT1
END IF

CFLAG_P = PF

CFLAG_P(24) = .FALSE.
CFLAG_P(32) = .FALSE.

IF (.NOT.ANY(CFLAG_P(:))) THEN
   WRITE(IU06,*) '*****************************************************'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*   FATAL ERROR IN SUB. SET_PARAMETER_OUTPUT_FLAGS  *'
   WRITE(IU06,*) '*   ==============================================  *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*  ALL FLAGS FOR PARAMETER OUTPUT ARE .FALSE.       *'
   WRITE(IU06,*) '*  CORRECT USER INPUT                               *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*          PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*****************************************************'
   CALL ABORT1
END IF

END SUBROUTINE SET_PARAMETER_OUTPUT_FLAGS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_RADIATION_OUTPUT_FLAGS (PF)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

LOGICAL, INTENT(IN) :: PF(:)   !! PRINTER FLAGS.

! ---------------------------------------------------------------------------- !

IF (SIZE(PF).NE.NOUT_R) THEN
   WRITE(IU06,*) '*  PROGRAM NEEDS ', NOUT_R,' FLAGS FOR RADIATION OUTPUT *'
   WRITE(IU06,*) '*  NUMBER OF PRINTER FLAGS IS : ', SIZE(PF)
   CALL ABORT1
END IF

CFLAG_R = PF

CFLAG_R(4)  = .FALSE.    !! CORRECT DUMMY PARAMETER.

IF (.NOT.ANY(CFLAG_R(:))) THEN
   WRITE(IU06,*) '*****************************************************'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*   FATAL ERROR IN SUB. SET_RADIATION_OUTPUT_FLAGS  *'
   WRITE(IU06,*) '*   ==============================================  *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*  ALL FLAGS FOR RADIATION OUTPUT ARE .FALSE.       *'
   WRITE(IU06,*) '*  CORRECT USER INPUT                               *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*          PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*****************************************************'
   CALL ABORT1
END IF

END SUBROUTINE SET_RADIATION_OUTPUT_FLAGS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_SPECTRA_OUTPUT_FLAGS (PF)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

LOGICAL, INTENT(IN) :: PF(:)   !! PRINTER FLAGS.

! ---------------------------------------------------------------------------- !

IF (SIZE(PF).NE.NOUT_S) THEN
   WRITE(IU06,*) '*  PROGRAM NEEDS ', NOUT_S,' FLAGS FOR SPECTRA OUTPUT *'
   WRITE(IU06,*) '*  NUMBER OF PRINTER FLAGS IS : ', SIZE(PF)
   CALL ABORT1
END IF

CFLAG_S = PF

CFLAG_S(4) = .FALSE.    !! CORRECT DUMMY PARAMETER.

IF (.NOT.ANY(CFLAG_S(:))) THEN
   WRITE(IU06,*) '*****************************************************'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*   FATAL ERROR IN SUB. SET_SPECTRA_OUTPUT_FLAGS    *'
   WRITE(IU06,*) '*   ============================================    *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*  ALL FLAGS FOR SPECTRA OUTPUT ARE .FALSE.         *'
   WRITE(IU06,*) '*  CORRECT USER INPUT                               *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*          PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*****************************************************'
   CALL ABORT1
END IF


END SUBROUTINE SET_SPECTRA_OUTPUT_FLAGS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_GRID_USER

INTEGER  :: I

WRITE(IU06,'(''1'')')
WRITE(IU06,*) ' USER INPUT PROG. PRINT_GRID_FILE:'
WRITE(IU06,*) '  '
IF (NOUTT.EQ.0) THEN
   WRITE(IU06,*) ' START  DATE (FORMAT:YYYYMMDDHHMMSS) : ',CDATEA,             &
&                ' END DATE :', CDATEE
   WRITE(IU06,*) '  '
   WRITE(IU06,*) ' OUTPUT EVERY ',IDELDO, ' SECONDS'
ELSE
   WRITE(IU06,*) ' GRIDS ARE PRINTED AT:'
   DO I = 1, NOUTT
      WRITE(IU06,'(5(1X,A14),/)') COUTT(I)
   END DO
END IF
WRITE(IU06,*) '  '
WRITE(IU06,*) ' INPUT FILE HANDLING:'
WRITE(IU06,*) ' FILE ID IS ..................... ', FILE01
WRITE(IU06,*) ' THE FIRST FILE DATE IS ......... ', CDTFILE
IF (IDFILE.GT.0) THEN
   WRITE(IU06,*) ' A NEW FILE WILL BE FETCHED EVERY ', IDFILE, ' SECONDS'
ELSE
   WRITE(IU06,*) ' A NEW FILE WILL NOT BE FETCHED'
END IF
WRITE(IU06,*) '  '
WRITE(IU06,*) ' LIST OF OUTPUTS TO BE GENERATED:'
WRITE(IU06,*) '  '
DO I=1,NOUT_P
   IF (CFLAG_P(I))  WRITE(IU06,'(1X,A50)') TITL_P(I)
END DO
WRITE(IU06,*) '  '
IF (REGULAR) THEN
   WRITE(IU06,*) ' REDUCED GRIDS ARE INTERPOLATED TO REGULAR GRIDS'
ELSE
   WRITE(IU06,*) ' REDUCED GRIDS ARE NOT INTERPOLATED TO REGULAR GRIDS'
END IF
WRITE(IU06,*) '  '

END SUBROUTINE PRINT_GRID_USER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_RADIATION_USER

INTEGER  :: I

WRITE(IU06,'(''1'')')
WRITE(IU06,*) ' USER INPUT PROG. PRINT_RADIATION_FILE:'
WRITE(IU06,*) '  '
IF (NOUTT.EQ.0) THEN
   WRITE(IU06,*) ' START  DATE (FORMAT:YYYYMMDDHHMMSS) : ',CDATEA,             &
&                ' END DATE :', CDATEE
   WRITE(IU06,*) '  '
   WRITE(IU06,*) ' OUTPUT EVERY ',IDELDO, ' SECONDS'
ELSE
   WRITE(IU06,*) ' GRIDS ARE PRINTED AT:'
   DO I = 1, NOUTT
      WRITE(IU06,'(5(1X,A14),/)') COUTT(I)
   END DO
END IF
WRITE(IU06,*) '  '
WRITE(IU06,*) ' INPUT FILE HANDLING:'
WRITE(IU06,*) ' FILE ID IS ..................... ', FILE01
WRITE(IU06,*) ' THE FIRST FILE DATE IS ......... ', CDTFILE
IF (IDFILE.GT.0) THEN
   WRITE(IU06,*) ' A NEW FILE WILL BE FETCHED EVERY ', IDFILE, ' SECONDS'
ELSE
   WRITE(IU06,*) ' A NEW FILE WILL NOT BE FETCHED'
END IF
WRITE(IU06,*) '  '
WRITE(IU06,*) ' LIST OF OUTPUTS TO BE GENERATED:'
WRITE(IU06,*) '  '
DO I=1,NOUT_R
   IF (CFLAG_R(I))  WRITE(IU06,'(1X,A50)') TITL_R(I)
END DO
WRITE(IU06,*) '  '
IF (REGULAR) THEN
   WRITE(IU06,*) ' REDUCED GRIDS ARE INTERPOLATED TO REGULAR GRIDS'
ELSE
   WRITE(IU06,*) ' REDUCED GRIDS ARE NOT INTERPOLATED TO REGULAR GRIDS'
END IF
WRITE(IU06,*) '  '

END SUBROUTINE PRINT_RADIATION_USER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_TIME_USER

INTEGER  :: I
character (len=len_coor) :: ftext1, ftext2

WRITE(IU06,'(''1'')')
WRITE(IU06,*) ' USER INPUT PROG. PRINT_TIME:'
WRITE(IU06,*) '  '
WRITE(IU06,*) ' START  DATE (FORMAT:YYYYMMDDHHMMSS) : ',CDATEA,                &
&             ' END DATE :', CDATEE
WRITE(IU06,*) '  '
WRITE(IU06,*) ' OUTPUT EVERY ',IDELDO, ' SECONDS'
WRITE(IU06,*) '  '
WRITE(IU06,*) ' INPUT FILE HANDLING:'
WRITE(IU06,*) ' FILE ID IS ..................... ', FILE01
WRITE(IU06,*) ' THE FIRST FILE DATE IS ......... ', CDTFILE
IF (IDFILE.GT.0) THEN
   WRITE(IU06,*) ' A NEW FILE WILL BE FETCHED EVERY ', IDFILE, ' SECONDS'
ELSE
   WRITE(IU06,*) ' A NEW FILE WILL NOT BE FETCHED'
END IF
WRITE(IU06,*) '  '
WRITE(IU06,*) ' LIST OF OUTPUT SITES TO BE PROCESSED:'
WRITE(IU06,*) '  '
WRITE(IU06,*) ' TOTAL NUMBER OF SITES IS.........', NOUTP
WRITE(IU06,*) '  '
WRITE(IU06,'('' |   LONGITUDE   |    LATITUDE   |       SITE NAME      |'')')
WRITE(IU06,'('' |---------------|---------------|----------------------|'')') 
DO I = 1,NOUTP
   ftext1 = write_coor_text (outlong(i))
   ftext2 = write_coor_text (outlat(i))
   WRITE(IU06,'('' | '',A,'' | '',A,'' | '',A20,'' | '')')                     &
&     ftext1, ftext2, name(i)
END DO
WRITE(IU06,*) '  '

END SUBROUTINE PRINT_TIME_USER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_SPECTRA_USER

INTEGER  :: I
character (len=len_coor) :: ftext1, ftext2

WRITE(IU06,'(''1'')')
WRITE(IU06,*) ' USER INPUT PROG. PRINT_SPECTRA:'
WRITE(IU06,*) '  '
IF (NOUTT.EQ.0) THEN
   WRITE(IU06,*) ' START  DATE (FORMAT:YYYYMMDDHHMMSS) : ',CDATEA,             &
&                ' END DATE :', CDATEE
   WRITE(IU06,*) '  '
   WRITE(IU06,*) ' OUTPUT EVERY ',IDELDO, ' SECONDS'
ELSE
   WRITE(IU06,*) ' SPECTRA ARE PRINTED AT:'
   DO I = 1, NOUTT
      WRITE(IU06,'(5(1X,A14),/)') COUTT(I)
   END DO
END IF
WRITE(IU06,*) '  '
WRITE(IU06,*) ' INPUT FILE HANDLING:'
WRITE(IU06,*) ' FILE ID IS ..................... ', FILE01
WRITE(IU06,*) ' THE FIRST FILE DATE IS ......... ', CDTFILE
IF (IDFILE.GT.0) THEN
   WRITE(IU06,*) ' A NEW FILE WILL BE FETCHED EVERY ', IDFILE, ' SECONDS'
ELSE
   WRITE(IU06,*) ' A NEW FILE WILL NOT BE FETCHED'
END IF
WRITE(IU06,*) '  '
WRITE(IU06,*) ' LIST OF OUTPUT SITES TO BE PROCESSED:'
WRITE(IU06,*) '  '
WRITE(IU06,*) ' TOTAL NUMBER OF SITES IS.........', NOUTP
WRITE(IU06,*) '  '
WRITE(IU06,'('' |   LONGITUDE   |    LATITUDE   |       SITE NAME      |'')')
WRITE(IU06,'('' |---------------|---------------|----------------------|'')') 
DO I = 1,NOUTP
   ftext1 = write_coor_text (outlong(i))
   ftext2 = write_coor_text (outlat(i))
   WRITE(IU06,'('' | '',A,'' | '',A,'' | '',A20,'' | '')')                     &
&     ftext1, ftext2, name(i)
END DO

WRITE(IU06,*) '  '
WRITE(IU06,*) ' LIST OF OUTPUTS TO BE GENERATED:'
WRITE(IU06,*) '  '
DO I=1,NOUT_S
   IF (CFLAG_S(I))  WRITE(IU06,'(1X,A50)') TITL_S(I)
END DO
WRITE(IU06,*) '  '

END SUBROUTINE PRINT_SPECTRA_USER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_EXTRACTION (LONGITUDE, LATITUDE, I, K)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    PREPARE_EXTRACTION -  COMPUTE GRID INICES FROM LATITUES AND LANGITUDES    !
!                                                                              !
!      H. GUNTHER       GKSS      DECEMBER 1998                                !
!      H. GUNTHER       HZG       DECEMBER 2010      RE-ORGANISED              !
!                                                                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE GRID INICES FROM LATITUES AND LANGITUDES.                   !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NEAREST NEIGHBOUR. INPUT LATITUDES AND LONGITUDES ARE CHANGED TO       !
!       TO THE COORDINATES OF THE NEAREST GRID POINT.                          !
!       IF A POINT IS NOT IN THE GRID THE FIRST GRID INDEX IS -1.              !
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

INTEGER, INTENT(INOUT) :: LONGITUDE(:)   !! LONGITUTES (DEGREE)
INTEGER, INTENT(INOUT) :: LATITUDE(:)    !! LATITUDES (DEGREE)
INTEGER, INTENT(OUT)   :: I(:)           !! FIRST GRID INDICES.
INTEGER, INTENT(OUT)   :: K(:)           !! SECOND GRID INDICES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CHECK GRID.                                                           !
!        -----------                                                           !

IF (NX.GT.1) THEN
   XDELLO = (AMOEAP-AMOWEP)/(NX-1)
ELSE
   XDELLO = 1.
END IF
IF (NY.GT.1) THEN
   XDELLA = (AMONOP-AMOSOP)/(NY-1)
ELSE
   XDELLA = 1.
END IF

PER = PERIODIC(AMOWEP, AMOEAP, XDELLO, NX) !! IS GRID EAST-WEST PERIODIC?

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE INDICES.                                                      !
!        ----------------                                                      !

K(:) = NINT(REAL(LATITUDE(:)-AMOSOP)/REAL(XDELLA)+1.)
WHERE (K(:).GE.1.AND.K(:).LE.NY)
   I(:) = NINT(REAL(MOD(LONGITUDE(:)-AMOWEP+2*M_S_PER,M_S_PER))                &
&            / REAL(ZDELLO(K(:)))+1.)
ELSEWHERE
   I(:) = -1
END WHERE
   
WHERE (I(:).EQ.NLON_RG(K(:))+1 .AND. PER) I(:) = 1
WHERE (I(:).EQ.0               .AND. PER) I(:) = NLON_RG(K(:))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. MARK GRIDPOINTS OUTSIDE OF GRID BY -1.                                !
!        OR RE-COMPUTE LAT. AND LONG. OF OUTPUT GRID POINTS.                   !
!        ---------------------------------------------------                   !

WHERE (I(:).LT.1.OR.I(:).GT.NLON_RG(K(:)).OR.K(:).LT.1.OR.K(:).GT.NY) 
   I(:) = -1
ELSEWHERE
   LONGITUDE(:) = AMOWEP+(I(:)-1)*ZDELLO(K(:))
   LATITUDE(:)  = AMOSOP+(K(:)-1)*XDELLA
END WHERE

END SUBROUTINE PREPARE_EXTRACTION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_PRINT_MODULE

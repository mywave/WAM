MODULE WAM_OUTPUT_SET_UP_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS: OUTPUT TIMES, FLAGS.                                 !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE          !! COORDINATE TYPE AND PROCEDURES

USE WAM_GENERAL_MODULE, ONLY:  &
&       ABORT1,                &  !! TERMINATE PROCESSING.
&       DIFDATE,               &  !! COMPUTES TIME DIFFERENCE.
&       OPEN_FILE,             &  !! OPENS A FILE.
&       INCDATE                   !! UPDATE DATE TIME GROUP.

USE WAM_GRID_MODULE,    ONLY:  & 
&       FIND_SEA_POINT            !! FIND SEA POINT NUMBER.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST, IU20, FILE20, IU25, FILE25

USE WAM_TIMOPT_MODULE,  ONLY: CDATEA, CDATEE, IDELPRO, CDTPRO,                 &
&                             SHALLOW_RUN, REFRACTION_C_RUN, COLDSTART

use wam_special_module, only: ispec2d, ispecode

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

INTEGER :: I, lent

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NEXT OUTPUT TIMES AND TIME INCREMENTS.                                !
!        --------------------------------------                                !

CHARACTER (LEN=14) :: CDTINTT     !! NEXT DATE TO WRITE INTEG. PARA. (TOTAL).
CHARACTER (LEN=14) :: CDTSPT      !! NEXT DATE TO WRITE SPECTRA (TOTAL).
INTEGER            :: IDELINT = 0 !! INTEG. PARAMETER OUTPUT TIMESTEP IN SECONDS.
INTEGER            :: IDELSPT = 0 !! SPECTRA OUTPUT TIMESTEP IN SECONDS.


! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. OUTPUT FILE TIMESTEP AND DATE.                                        !
!        ------------------------------                                        !

INTEGER            :: IDEL_OUT = -1    !! TIMESTEP TO SAVE OUTPUT FILES.
CHARACTER (LEN=14) :: CDT_OUT  = ' '   !! NEXT DATE TO SAVE OUTPUT FILES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. OUTPUT TIMES AS DEFINED IN INPUT FILE.                                !
!        --------------------------------------                                !

INTEGER                         :: NOUTT = 0 !! NUMBER OF OUTPUT TIMES.
CHARACTER( LEN=14), ALLOCATABLE :: COUTT(:)  !! OUTPUT TIMES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. OUTPUT FLAGS FOR PARAMETER AND SPECTRA.                               !
!        ---------------------------------------                               !

INTEGER, PARAMETER :: NOUT_P = 40        !! NUMBER OF INTEGRATED PARAMETER.

LOGICAL, DIMENSION(NOUT_P) :: FFLAG_P    !! FILE OUTPUT FLAG.
LOGICAL, DIMENSION(NOUT_P) :: PFLAG_P    !! PRINTER OUTPUT FLAG.
LOGICAL, DIMENSION(NOUT_P) :: CFLAG_P    !! COMPUTATION FLAG.

LOGICAL :: FFLAG20 = .FALSE. !! .TRUE. IF FIELDS ARE WRITTEN TO FILE20.
LOGICAL :: PFLAG20 = .FALSE. !! .TRUE. IF FIELDS ARE PRINTED.
LOGICAL :: CFLAG20 = .FALSE. !! .TRUE. IF ANY COMPUTATION OF FIELDS.

INTEGER, PARAMETER :: NOUT_S = 4         !! NUMBER OF SPECTRA TYPES.
integer, parameter :: npout  = 36
   
LOGICAL, DIMENSION(NOUT_S) :: FFLAG_S    !! FILE OUTPUT FLAG. 
LOGICAL, DIMENSION(NOUT_S) :: PFLAG_S    !! PRINTER OUTPUT FLAG.
LOGICAL, DIMENSION(NOUT_S) :: CFLAG_S    !! COMPUTATION FLAG.

LOGICAL :: FFLAG25 = .FALSE. !! .TRUE. IF ANY SPECTRA ARE WRITTEN TO FILE25.
LOGICAL :: PFLAG25 = .FALSE. !! .TRUE. IF ANY SPECTRA ARE PRINTED.
LOGICAL :: CFLAG25 = .FALSE. !! .TRUE. IF ANY COMPUTATION OF SPECTRA.
logical :: ready_outf

character (len=128) :: owpath

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. TITLE FOR OUTPUT PARAMETER AND SPECTRA.                               !
!        ---------------------------------------                               !

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

CHARACTER(LEN=60), DIMENSION(NOUT_S) :: TITL_S = (/                  &
& ' SPECTRUM                                                   ',    &   !!  1
& ' SEA SPECTRUM                                               ',    &   !!  2
& ' SWELL SPECTRUM                                             ',    &   !!  3
& ' DUMMY                                                      '/)       !!  4

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. SCALING FACTORS FOR OUTPUT PARAMETER.                                 !
!        -------------------------------------                                 !

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
!      7. OUTPUT POINTS FOR SPECTRA.                                           !
!         --------------------------                                           !

INTEGER                         :: NOUTP = 0  !! NUMBER OF OUTPUT POINTS.
INTEGER, ALLOCATABLE            :: OUTLAT (:) !! LATITUDE OF POINTS [M_SEC].
INTEGER, ALLOCATABLE            :: OUTLONG(:) !! LONGITUDE OF POINTS [M_SEC].
CHARACTER (LEN=20), ALLOCATABLE :: NAME(:)    !! NAMES OF OUTPUT SITES.

INTEGER, ALLOCATABLE            :: IJAR(:)    !! GRIDPOINT NUMBER OF OUTPUT POINT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE SET_OUTPUT_FILE_STEP
   MODULE PROCEDURE SET_OUTPUT_FILE_STEP
END INTERFACE
PUBLIC SET_OUTPUT_FILE_STEP

INTERFACE SET_PARAMETER_OUTPUT_FLAGS    !! SETS FLAGS FOR PARAMETER OUTPUT.
   MODULE PROCEDURE SET_PARAMETER_OUTPUT_FLAGS
END INTERFACE
PUBLIC SET_PARAMETER_OUTPUT_FLAGS

INTERFACE SET_SPECTRA_OUTPUT_FLAGS      !! SETS FLAGS FOR SPECTRA OUTPUT.
   MODULE PROCEDURE SET_SPECTRA_OUTPUT_FLAGS
END INTERFACE
PUBLIC SET_SPECTRA_OUTPUT_FLAGS

INTERFACE SET_OUTPUT_SITES              !! SETS SITES FOR OUTPUT.
   MODULE PROCEDURE SET_OUTPUT_SITES_C  !! CHARACTER VERSION
   MODULE PROCEDURE SET_OUTPUT_SITES_D  !! DEGREE VERSION
   MODULE PROCEDURE SET_OUTPUT_SITES_M  !! M_SEC VERSION
END INTERFACE
PUBLIC SET_OUTPUT_SITES

INTERFACE SET_OUTPUT_TIMES              !! SETS TIMES FOR OUTPUT.
   MODULE PROCEDURE SET_OUTPUT_TIMES_F  !! FIXED OUTPUT TIMES.
   MODULE PROCEDURE SET_OUTPUT_TIMES_I  !! OUTPUT INCREMENTS.
END INTERFACE
PUBLIC SET_OUTPUT_TIMES

interface set_spectral_code             !! sets code for spectral output
  module procedure set_spectral_code    !! and hours od 2d spectral output
end interface
public set_spectral_code

interface set_ready_outfile_flag        !! sets ready outfile flag
  module procedure set_ready_outfile_flag
end interface
public set_ready_outfile_flag

interface set_ready_outfile_directory   !! sets ready outfile directory
  module procedure set_ready_outfile_directory 
end interface
public set_ready_outfile_directory

INTERFACE PREPARE_OUTPUT                !! PREPARES OUTPUT.
   MODULE PROCEDURE PREPARE_OUTPUT
END INTERFACE
PUBLIC PREPARE_OUTPUT

INTERFACE PRINT_OUTPUT_STATUS           !! PRINT OUTPUT SETTING.
   MODULE PROCEDURE PRINT_OUTPUT_STATUS
END INTERFACE
PUBLIC PRINT_OUTPUT_STATUS

INTERFACE UPDATE_OUTPUT_TIME            !! UPDATES OUTPUT TIMES.
   MODULE PROCEDURE UPDATE_OUTPUT_TIME
END INTERFACE
PUBLIC UPDATE_OUTPUT_TIME

INTERFACE SAVE_OUTPUT_FILES             !! SAVES AND OPENS OUTPUT FILES.
   MODULE PROCEDURE SAVE_OUTPUT_FILES 
END INTERFACE
PUBLIC SAVE_OUTPUT_FILES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE MAKE_OUTPUT_SITES           !! INDEX OF OUTPUT POINTS.
   MODULE PROCEDURE MAKE_OUTPUT_SITES
END INTERFACE
PRIVATE MAKE_OUTPUT_SITES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_OUTPUT_FILE_STEP (STEP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN) :: STEP    !! TIME INCREMENT TO SAVE FILES.

! ---------------------------------------------------------------------------- !

IF (STEP.GT.0) THEN
   IDEL_OUT = STEP
ELSE
   IDEL_OUT = 0
END IF

END SUBROUTINE SET_OUTPUT_FILE_STEP

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine set_spectral_code (ihour, icode)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     interface variables                                                      !
!     -------------------                                                      !

integer, intent(in) :: ihour, icode

! ---------------------------------------------------------------------------- !

if (icode>0) then
   ispecode = icode
else
   ispecode = 0
endif
if (ihour>0) then
   ispec2d = ihour
else
   ispec2d = 0
endif

end subroutine set_spectral_code

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine set_ready_outfile_flag (op)

logical, intent(in) :: op       !! ready output file flag

ready_outf = op

end subroutine set_ready_outfile_flag

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine set_ready_outfile_directory (oname)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

character (len=*), intent(in) :: oname     !! full path of ready file directory

! ---------------------------------------------------------------------------- !

lent = len_trim (oname)
if (lent>=0) owpath = oname(1:lent)

end subroutine set_ready_outfile_directory
   
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_PARAMETER_OUTPUT_FLAGS (PF, FF)

LOGICAL, INTENT(IN) :: PF(:)   !! PRINTER FLAGS.
LOGICAL, INTENT(IN) :: FF(:)   !! FILE FLAGS.

IF (SIZE(PF).NE.NOUT_P .OR. SIZE(FF).NE.NOUT_P) THEN
   WRITE(IU06,*) '*  PROGRAM NEEDS ', NOUT_P,' FLAGS FOR PARAMETER OUTPUT *'
   WRITE(IU06,*) '*  NUMBER OF PRINTER FLAGS IS : ', SIZE(PF)
   WRITE(IU06,*) '*  NUMBER OF FILE    FLAGS IS : ', SIZE(FF)
   CALL ABORT1
END IF

PFLAG_P = PF
FFLAG_P = FF

FFLAG_P(24) = .FALSE.    !! CORRECT DUMMY PARAMETER.
FFLAG_P(32) = .FALSE.

PFLAG_P(24) = .FALSE.
PFLAG_P(32) = .FALSE.

CFLAG_P = FFLAG_P.OR.PFLAG_P
FFLAG20 = ANY(FFLAG_P(:))
PFLAG20 = ANY(PFLAG_P(:))
CFLAG20 = FFLAG20.OR.PFLAG20

END SUBROUTINE SET_PARAMETER_OUTPUT_FLAGS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_SPECTRA_OUTPUT_FLAGS (PF, FF)

LOGICAL, INTENT(IN) :: PF(:)   !! PRINTER FLAGS.
LOGICAL, INTENT(IN) :: FF(:)   !! FILE FLAGS.

IF (SIZE(PF).NE.NOUT_S .OR. SIZE(FF).NE.NOUT_S) THEN
   WRITE(IU06,*) '*  PROGRAM NEEDS ', NOUT_S,' FLAGS FOR SPECTRA OUTPUT *'
   WRITE(IU06,*) '*  NUMBER OF PRINTER FLAGS IS : ', SIZE(PF)
   WRITE(IU06,*) '*  NUMBER OF FILE    FLAGS IS : ', SIZE(FF)
   CALL ABORT1
END IF

PFLAG_S = PF
FFLAG_S = FF

PFLAG_S(4) = .FALSE.    !! CORRECT DUMMY PARAMETER.
FFLAG_S(4) = .FALSE.    !! CORRECT DUMMY PARAMETER.

CFLAG_S = FFLAG_S.OR.PFLAG_S

FFLAG25 = ANY(FFLAG_S(:))
PFLAG25 = ANY(PFLAG_S(:))
CFLAG25 = FFLAG25.OR.PFLAG25

END SUBROUTINE SET_SPECTRA_OUTPUT_FLAGS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_OUTPUT_SITES_C (LONG, LAT, NA)

CHARACTER (LEN=LEN_COOR), INTENT(IN) :: LONG(:) !! LONGITUDES.
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: LAT(:)  !! LATITUDES.
CHARACTER (LEN=*), INTENT(IN) :: NA(:)   !! SITE NAME.

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

SUBROUTINE SET_OUTPUT_TIMES_F (TIME)

CHARACTER*14, INTENT(IN) :: TIME(:) !! OUTPUT TIMES.

NOUTT = COUNT(TIME.NE.' ')
IF (NOUTT.GT.0) THEN
   ALLOCATE (COUTT(NOUTT))
   COUTT = TIME(1:NOUTT)
END IF

END SUBROUTINE SET_OUTPUT_TIMES_F

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_OUTPUT_TIMES_I (INT, SPE)

INTEGER,      INTENT(IN) :: INT     !! OUTPUT TIME INCREMENT FOR PARAMETER.
INTEGER,      INTENT(IN) :: SPE     !! OUTPUT TIME INCREMENT FOR SPECTRA.

IDELINT = INT
IDELSPT = SPE

END SUBROUTINE SET_OUTPUT_TIMES_I

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_OUTPUT

! ---------------------------------------------------------------------------- !

INTEGER :: J, NW, ISHIFT, IW(1:NOUTT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CHECK CONSISTENCY FOR PARAMETER OUTPUT AND MODEL OPTIONS.             !
!        ---------------------------------------------------------             !

IF (.NOT.SHALLOW_RUN) THEN
   FFLAG_P(6) = .FALSE.
   PFLAG_P(6) = .FALSE.
END IF
IF (.NOT.REFRACTION_C_RUN) THEN
   FFLAG_P(7:8) = .FALSE.
   PFLAG_P(7:8) = .FALSE.
END IF

CFLAG_P(:) = FFLAG_P(:).OR.PFLAG_P(:)
FFLAG20 = ANY(FFLAG_P(:))
PFLAG20 = ANY(PFLAG_P(:))
CFLAG20 = FFLAG20.OR.PFLAG20

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK CONSISTENCY FOR SPECTRA OUTPUT.                                 !
!        -------------------------------------                                 !

IF (CFLAG25 .AND. NOUTP.GT.0) CALL MAKE_OUTPUT_SITES

IF (NOUTP.LE.0 .AND. CFLAG25)  THEN
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) ' +                                             +'
   WRITE(IU06,*) ' +    WARNING ERROR IN SUB. PREPARE_OUTPUT     +'
   WRITE(IU06,*) ' +    ====================================     +'
   WRITE(IU06,*) ' +                                             +'
   WRITE(IU06,*) ' + NUMBER OF OUTPUT SITES FOR SPECTRA          +'
   WRITE(IU06,*) ' + FOUND IN THE GRID IS ZERO.                  +'
   WRITE(IU06,*) ' + SPECTRA OUTPUT FLAGS CHANGED TO NO OUTPUT.  +'
   WRITE(IU06,*) ' +                                             +'
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++'
   NOUTP = 0
   FFLAG_S(:) = .FALSE.
   PFLAG_S(:) = .FALSE.
   CFLAG_S(:) = FFLAG_S(:).OR.PFLAG_S(:)
   FFLAG25 = ANY(FFLAG_S(:))
   PFLAG25 = ANY(PFLAG_S(:))
   CFLAG25 = FFLAG25.OR.PFLAG25
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CHECK CONSISTENCY OUTPUT TIMES AND MODEL TIMESTEPS.                   !
!        ---------------------------------------------------                   !

IF (NOUTT.GT.0) THEN
   IW = 0
   DO J = 1,NOUTT
      IF (COUTT(J).LT.CDATEA .OR. COUTT(J).GT.CDATEE) THEN
         IW(J) = J
         CYCLE
      END IF
      CALL DIFDATE (CDATEA, COUTT(J), ISHIFT)
      IF (ISHIFT.LT.0 .OR. MOD(ISHIFT,IDELPRO).NE.0) IW(J) = J
   END DO
   NW = COUNT(IW.NE.0)
   IF (NW.GT.0) THEN
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +     WARNING ERROR IN SUB. PREPARE_OUTPUT      +'
      WRITE(IU06,*) ' +     ====================================      +'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + THE FOLLOWING OUTPUT DATES ARE NOT AT THE END +'
      WRITE(IU06,*) ' + OF A PROPAGATION TIMESTEP  OR                 +'
      WRITE(IU06,*) ' + ARE NOT IN THE MODEL INTEGRATION PERIOD.      +'
      WRITE(IU06,*) ' + THE FOLLOWING DATES WILL BE IGNORED:          +'
      WRITE(IU06,*) ' +                                               +'
      DO J = 1,NOUTT
         IF (IW(J).EQ.J) WRITE(IU06,'(5X,A14)') COUTT(J)
      END DO
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'

      IF (NW.LT.NOUTT) THEN
         COUTT = PACK(COUTT,IW.EQ.0)
         NOUTT = NOUTT-NW
      ELSE IF (NW.GE.NOUTT) THEN
         WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
         WRITE(IU06,*) ' +                                               +'
         WRITE(IU06,*) ' +     WARNING ERROR IN SUB. PREPARE_OUTPUT      +'
         WRITE(IU06,*) ' +     ====================================      +'
         WRITE(IU06,*) ' +                                               +'
         WRITE(IU06,*) ' + ALL OUTPUT DATES ARE NOT AT THE END OF A      +'
         WRITE(IU06,*) ' + PROPAGATION TIMESTEP  OR                      +'
         WRITE(IU06,*) ' + ARE NOT IN THE MODEL INTEGRATION PERIOD.      +'
         WRITE(IU06,*) ' + PROGRAM WILL TRY WITH OUTPUT INCREMENTS.      +'
         WRITE(IU06,*) ' +                                               +'
         WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
         NOUTT = 0
      END IF
   END IF
END IF

IF (NOUTT.LE.0) THEN
   IF (CFLAG20 .AND. IDELINT.LT.IDELPRO) THEN
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +     WARNING ERROR IN SUB. PREPARE_OUTPUT      +'
      WRITE(IU06,*) ' +     ====================================      +'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT TIMESTEP FOR INTEGRATED PARAMETER IS   +'
      WRITE(IU06,*) ' + LESS THAN:                                   +'
      WRITE(IU06,*) ' + PROPAGATION     TIMESTEP IS: ',IDELPRO, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT          TIMESTEP IS: ',IDELINT, ' SECONDS'
      IDELINT = IDELPRO
      WRITE(IU06,*) ' + OUTPUT TIMESTEP  CHANGED TO: ',IDELINT, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
   ELSE IF (CFLAG20 .AND. MOD(IDELINT,IDELPRO).NE.0) THEN
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +     WARNING ERROR IN SUB. PREPARE_OUTPUT      +'
      WRITE(IU06,*) ' +     ====================================      +'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT TIMESTEP FOR INTEGRATED PARAMETER IS   +'
      WRITE(IU06,*) ' + NOT AN INTEGER MULTIPLE OF:                   +'
      WRITE(IU06,*) ' + PROPAGATION     TIMESTEP IS: ',IDELPRO, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT          TIMESTEP IS: ',IDELINT, ' SECONDS'
      IDELINT = MAX(NINT(REAL(IDELINT)/REAL(IDELPRO)),1) * IDELPRO
      WRITE(IU06,*) ' + OUTPUT TIMESTEP  CHANGED TO: ',IDELINT, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
   ELSE IF (.NOT.CFLAG20) THEN
      IDELINT = 0
   END IF

   IF (CFLAG25 .AND. IDELSPT.LT.IDELPRO) THEN
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +     WARNING ERROR IN SUB. PREPARE_OUTPUT      +'
      WRITE(IU06,*) ' +     ====================================      +'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT TIMESTEP FOR SPECTRA IS LESS THAN:     +'
      WRITE(IU06,*) ' + PROPAGATION     TIMESTEP IS: ',IDELPRO, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT          TIMESTEP IS: ',IDELSPT, ' SECONDS'
      IDELSPT = IDELPRO
      WRITE(IU06,*) ' + OUTPUT TIMESTEP  CHANGED TO: ',IDELSPT, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
   ELSE IF (CFLAG25 .AND. MOD(IDELSPT,IDELPRO).NE.0) THEN
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +     WARNING ERROR IN SUB. PREPARE_OUTPUT      +'
      WRITE(IU06,*) ' +     ====================================      +'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT TIMESTEP FOR SPECTRA IS NOT AN         +'
      WRITE(IU06,*) ' + INTEGER MULTIPLE OF:                          +'
      WRITE(IU06,*) ' + PROPAGATION     TIMESTEP IS: ',IDELPRO, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT          TIMESTEP IS: ',IDELINT, ' SECONDS'
      IDELSPT = MAX(NINT(REAL(IDELSPT)/REAL(IDELPRO)),1)  * IDELPRO
      WRITE(IU06,*) ' + OUTPUT TIMESTEP  CHANGED TO: ',IDELSPT, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
   ELSE IF (.NOT.CFLAG25) THEN
      IDELSPT = 0
   END IF
END IF

! IF SPECTRA OUTPUT COMPUTE SOME INTEGRATED PARAMETERS.

CFLAG_P( 1: 3) = CFLAG_P( 1: 3) .OR. CFLAG25     !! WIND
CFLAG_P(    6) = CFLAG_P(    6) .OR. CFLAG25     !! DEPTH
CFLAG_P( 7: 8) = CFLAG_P( 7: 8) .OR. CFLAG25     !! CURRENT
CFLAG_P( 9:15) = CFLAG_P( 9:15) .OR. CFLAG_S(1)  !! PARATETER OF TOTAL SPECTRUM
CFLAG_P(17:23) = CFLAG_P(17:23) .OR. CFLAG_S(2)  !! PARATETER OF SEA   SPECTRUM
CFLAG_P(25:31) = CFLAG_P(25:31) .OR. CFLAG_S(3)  !! PARATETER OF SWELL SPECTRUM

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. INITIALIZE OUTPUT TIME VARIABLES.                                     !
!        ---------------------------------                                     !

CDTINTT = ' '
CDTSPT  = ' '
CDT_OUT = ' ' 

IF (.NOT.CFLAG20 .AND. .NOT.CFLAG25) THEN
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) ' +                                               +'
   WRITE(IU06,*) ' +     WARNING ERROR IN SUB. PREPARE_OUTPUT      +'
   WRITE(IU06,*) ' +     ====================================      +'
   WRITE(IU06,*) ' +                                               +'
   WRITE(IU06,*) ' +    THIS IS A MODEL RUN WITHOUT OUTPUT OF      +'
   WRITE(IU06,*) ' +         INTEGRATED WAVE PARAMETERS            +'
   WRITE(IU06,*) ' +          AND WITHOUT OUTPUT OF                +'
   WRITE(IU06,*) ' +               WAVE SPECTRA                    +'
   WRITE(IU06,*) ' +                                               +'
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
   RETURN
END IF
   
IF (NOUTT.LE.0) THEN 
   IF (CFLAG20) THEN
      CDTINTT = CDATEA
      IF (COLDSTART) CALL INCDATE(CDTINTT ,-IDELINT)
   END IF
   IF (CFLAG25) THEN
      CDTSPT = CDATEA
      IF (COLDSTART) CALL INCDATE(CDTSPT ,-IDELSPT)
   END IF
ELSE
   IF (CFLAG20) THEN
      IF (.NOT.COLDSTART) CDTINTT = CDATEA
   END IF
   IF (CFLAG25) THEN
      IF (.NOT.COLDSTART) CDTSPT = CDATEA
   END IF
END IF
CALL UPDATE_OUTPUT_TIME

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5.  OPEN FIRST OUTPUT FILES.                                             !
!        -------------------------                                             !

IF (NOUTT.GT.0 .AND. IDEL_OUT.GT.0) THEN
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) ' +                                               +'
   WRITE(IU06,*) ' +     WARNING ERROR IN SUB. PREPARE_OUTPUT      +'
   WRITE(IU06,*) ' +     ====================================      +'
   WRITE(IU06,*) ' +                                               +'
   WRITE(IU06,*) ' + SPECIAL OUTPUT DATES ARE DEFINED AND          +'
   WRITE(IU06,*) ' + OUTPUT FILES SAVE TIMESTEP IS POSITIVE.       +'
   WRITE(IU06,*) ' + OUTPUT FILES SAVE    TIMESTEP IS: ',IDEL_OUT, ' SECONDS'
   IDEL_OUT = 0
   WRITE(IU06,*) ' + OUTPUT FILES SAVE CHANGED TO    : ',IDEL_OUT, ' SECONDS'
   WRITE(IU06,*) ' +  (THE FILE DATE IS THE END OF RUN DATE)       +'
   WRITE(IU06,*) ' +                                               +'
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

IF (IDEL_OUT.GT.0) THEN
   IF (IDEL_OUT.LT.MAX(IDELINT,IDELSPT)) THEN
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +     WARNING ERROR IN SUB. PREPARE_OUTPUT      +'
      WRITE(IU06,*) ' +     ====================================      +'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT FILES SAVE TIMESTEP IS LESS THAN:      +'
      WRITE(IU06,*) ' + INTEGRATED PARAMETER TIMESTEP IS: ',IDELINT, ' SECONDS'
      WRITE(IU06,*) ' + SPECTRA              TIMESTEP IS: ',IDELSPT, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT FILES SAVE    TIMESTEP IS: ',IDEL_OUT, ' SECONDS'
      IDEL_OUT = MAX(IDELINT,IDELSPT)
      WRITE(IU06,*) ' + OUTPUT FILES SAVE CHANGED TO    : ',IDEL_OUT, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
   ELSE IF (MOD(IDEL_OUT,MAX(IDELINT,IDELSPT)).NE.0) THEN
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +     WARNING ERROR IN SUB. PREPARE_OUTPUT      +'
      WRITE(IU06,*) ' +     ====================================      +'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT FILES SAVE TIMESTEP IS NOT AN INTEGER  +'
      WRITE(IU06,*) ' + MULTIPLE OF THE MAXIMUM OF:                   +'
      WRITE(IU06,*) ' + INTEGRATED PARAMETER TIMESTEP IS: ',IDELINT, ' SECONDS'
      WRITE(IU06,*) ' + SPECTRA              TIMESTEP IS: ',IDELSPT, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' + OUTPUT FILES SAVE    TIMESTEP IS: ',IDEL_OUT, ' SECONDS'
      IDEL_OUT = MAX(NINT(REAL(IDEL_OUT)/REAL(MAX(IDELINT,IDELSPT))),1)        &
&               * MAX(IDELINT,IDELSPT)
      WRITE(IU06,*) ' + OUTPUT FILES SAVE CHANGED TO    : ',IDEL_OUT, ' SECONDS'
      WRITE(IU06,*) ' +                                               +'
      WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++'
   END IF
   CDT_OUT = CDTPRO
   IF (COLDSTART) CALL INCDATE(CDT_OUT ,-IDEL_OUT)
ELSE
   CDT_OUT = CDATEE
END IF

CALL SAVE_OUTPUT_FILES (IU20, FILE20, IU25, FILE25)
CALL INCDATE(CDT_OUT ,IDEL_OUT)

END SUBROUTINE PREPARE_OUTPUT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_OUTPUT_STATUS

INTEGER :: I
character (len=len_coor) :: ftext1, ftext2

WRITE(IU06,*) '  '
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '              MODEL OUTPUT SELECTION:'
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '  '

IF (NOUTT.GT.0) THEN
   WRITE(IU06,*) ' NUMBER OF OUTPUT TIMES IS NOUTT = ', NOUTT
   WRITE(IU06,*) ' OUTPUT WILL BE PROCESSED AT:'
   WRITE(IU06,'(60(2X,A14))') (COUTT(I),I=1,NOUTT)
   WRITE(IU06,*) '  '
END IF

IF (CFLAG20) THEN
   I = LEN_TRIM(FILE20)
   WRITE(IU06,*) ' TO PRINTER AND/OR FILE: ', FILE20(1:I),'YYYYMMDDHHMMSS'
   IF (NOUTT.EQ.0) THEN
      WRITE(IU06,*) ' INTEGRATED PARAMETER OUTPUT EVERY : ', IDELINT, ' SECONDS'
      WRITE(IU06,*) ' NEXT OUTPUT DATE IS...............: ', CDTINTT
   END IF
   IF (IDEL_OUT.GT.0) THEN
      WRITE(IU06,*) ' OUTPUT FILES ARE SAVED EVERY......:', IDEL_OUT, ' SECONDS'
      WRITE(IU06,*) ' NEXT DATE TO SAVE FILES IS........: ', CDT_OUT
   ELSE
      WRITE(IU06,*) ' ALL OUTPUT TIMES ARE SAVED THE SAME FILE'
      WRITE(IU06,*) ' THE FILE DATE IS THE END OF RUN DATE: ', CDT_OUT
   END IF
   WRITE(IU06,*) '                              F = FALSE   T = TRUE '
   WRITE(IU06,*) '                                                           ',&
&                                                        'PRINTER     UNIT'
   DO I=1,NOUT_P
      WRITE(IU06,*) TITL_P(I),'....', PFLAG_P( I),'......', FFLAG_P( I)
   END DO
ELSE
   WRITE(IU06,*) '  OUTPUT OF INTEGRATED PARAMETERS IS NOT REQUESTED '
END IF
WRITE(IU06,*) '  '

IF  (CFLAG25) THEN
   I = LEN_TRIM(FILE25)
   WRITE(IU06,*) ' TO PRINTER AND/OR FILE : ', FILE25(1:I),'YYYYMMDDHHMMSS'
   IF (NOUTT.EQ.0) THEN
      WRITE(IU06,*) ' SPECTRA OUTPUT EVERY .............: ', IDELSPT, ' SECONDS'
      WRITE(IU06,*) ' NEXT OUTPUT DATE IS...............: ', CDTSPT
   END IF
   IF (IDEL_OUT.GT.0) THEN
      WRITE(IU06,*) ' OUTPUT FILES ARE SAVED EVERY......:', IDEL_OUT, ' SECONDS'
      WRITE(IU06,*) ' NEXT DATE TO SAVE FILES IS........: ', CDT_OUT
   ELSE
      WRITE(IU06,*) ' ALL OUTPUT TIMES ARE SAVED THE SAME FILE'
      WRITE(IU06,*) ' THE FILE DATE IS THE END OF RUN DATE: ', CDT_OUT
   END IF
   WRITE(IU06,*) '                                                           ',&
&                                                        'PRINTER     UNIT'
   DO I=1,NOUT_S
      WRITE(IU06,*) TITL_S(I),'....', PFLAG_S( I),'......', FFLAG_S( I)
   END DO

   WRITE(IU06,*) '  '
   WRITE(IU06,*) ' OUTPUT SITES FOR SPECTRA:'
   WRITE(IU06,*) ' TOTAL NUMBER OF SITES IS..........:', NOUTP
   WRITE(IU06,*) '  '
   WRITE(IU06,'('' |   LONGITUDE   |    LATITUDE   |       SITE NAME      |'', &
&                  ''    POINT   |'')')
   WRITE(IU06,'('' |---------------|---------------|----------------------|'', &
&                  ''------------|'')')
   DO I = 1,NOUTP
      ftext1 = write_coor_text (outlong(i))
      ftext2 = write_coor_text (outlat(i))
      WRITE(IU06,'('' | '',A,'' | '',A,'' | '',A20,'' | '',I10,'' |'')')       &
&        ftext1, ftext2, name(i), ijar(i)
   END DO
ELSE
   WRITE(IU06,*) '  OUTPUT OF SPECTRA IS NOT REQUESTED '
   WRITE(IU06,*) '  OUTPUT SITES OR PARAMETER WERE NOT SPECIFIED'
   NOUTP = 0
END IF
WRITE(IU06,*) '  '

END SUBROUTINE PRINT_OUTPUT_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MAKE_OUTPUT_SITES

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MAKE_OUTPUT_SITES - ROUTINE TO COMPUTE OUTPUT INDICES.                     !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO COMPUTES THE INDICES OF SPECTRA OUTPUT POINTS.                      !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE LATITUDE AND LOGITUDE ARE CONVERTED TO INDICES.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: IO, NG
character (len=len_coor) :: ftext1, ftext2

LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NO OUTPUT POINTS SPECIFIED.                                           !
!        ---------------------------                                           !

IF (NOUTP.EQ.0) RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. SEARCH BLOCK NUMBER AND SEA POINT NUMBER.                             !
!        -----------------------------------------                             !

IF (.NOT.ALLOCATED(IJAR)) ALLOCATE (IJAR(1:NOUTP))

CALL FIND_SEA_POINT (NOUTP, OUTLAT, OUTLONG, IJAR)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. REMOVE OUTPUT POINTS WHICH ARE NOT IN GRID.                           !
!        -------------------------------------------                           !

IF (.NOT.ALLOCATED(MASK)) ALLOCATE(MASK(1:NOUTP))

MASK(1:NOUTP) = (IJAR(1:NOUTP).GT.0)
NG = COUNT(MASK(1:NOUTP))

IF (NG.LT.NOUTP) THEN
   WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                             +'
   WRITE (IU06,*) ' +  WARNING ERROR FROM SUB. MAKE_OUTPUT_SITES  +'
   WRITE (IU06,*) ' +  =========================================  +'
   WRITE (IU06,*) ' +                                             +'
   WRITE (IU06,*) ' + A SEAPOINT WAS NOT FOUND FOR ', NOUTP-NG
   WRITE (IU06,*) ' + OUTPUT SITE REQUESTS.                       +'
   WRITE (IU06,*) ' + THE FOLLOWING SITES ARE IGNORED.            +'
   DO IO = 1,NOUTP
      ftext1 = write_coor_text (outlong(io))
      ftext2 = write_coor_text (outlat(io))
      IF (IJAR(IO).LE.0) WRITE(IU06,'(4X,I5,2(2X,A),2X,A20)')                  &
&                        IO, ftext1, ftext2, name(io)
   END DO
   WRITE (IU06,*) ' + NUMBER OF OUTPUT POINTS IS NGOUT = ', NG
   WRITE (IU06,*) ' +                                             +'
   WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++'

   IJAR(1:NG) = PACK (IJAR(1:NOUTP), MASK(1:NOUTP))
   NAME(1:NG) = PACK (NAME(1:NOUTP), MASK(1:NOUTP))
   OUTLONG(1:NG) = PACK (OUTLONG(1:NOUTP), MASK(1:NOUTP))
   OUTLAT(1:NG)  = PACK (OUTLAT(1:NOUTP), MASK(1:NOUTP))
   NOUTP = NG
END IF

DEALLOCATE (MASK)

END SUBROUTINE MAKE_OUTPUT_SITES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE UPDATE_OUTPUT_TIME

INTEGER :: J, J1

IF (CFLAG20 .OR. CFLAG25) THEN
   IF (NOUTT.GT.0) THEN
      J1 = 0
      DO J = 1,NOUTT
         IF (CFLAG20.AND.CDTINTT.GE.COUTT(J)) CYCLE
         IF (CFLAG25.AND.CDTSPT.GE.COUTT(J)) CYCLE
         IF (J1.EQ.0) THEN
            J1 = J
         ELSE
            IF (COUTT(J).LT.COUTT(J1)) J1 = J
         END IF
      END DO
      IF (J1.NE.0) THEN
         IF (CFLAG20) CDTINTT = COUTT(J1)
         IF (CFLAG25) CDTSPT  = COUTT(J1)
      ELSE
         CDTINTT = ' '
         CDTSPT  = ' ' 
      END IF
   ELSE
      IF (CFLAG20 .AND. CDTPRO.GE.CDTINTT) CALL INCDATE (CDTINTT,IDELINT)
      IF (CFLAG25 .AND. CDTPRO.GE.CDTSPT) CALL INCDATE (CDTSPT,IDELSPT)
   END IF
END IF

END SUBROUTINE UPDATE_OUTPUT_TIME

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SAVE_OUTPUT_FILES (IU_PA, FILE_PA, IU_SP, FILE_SP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SAVE_OUTPUT_FILES - GETS/SAVES FILES FROM/TO MASS STORAGE.                 !
!                                                                              !
!     H. GUNTHER    GKSS/ECMWF    OCTOBER 1989                                 !
!     P. JANSSEN    KNMI          OCTOBER 1990   YMP-MODIFICATION              !
!     H. GUNTHER    GKSS/ECMWF    OCTOBER 1990   NEW FILE NAMES.               !
!     H. GUNTHER    GKSS          NOVEMBER 1999  NEW DATES AND FT90.           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       GETS OR SAVES FILES FROM / TO MASS STORAGE.                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!        OLD FILES ARE CLOSED AND                                              !
!        THE NEW FILES ARE OPENED WITH SUB. OPEN_FILE.                         !
!                                                                              !
!                                                                              !
!     REFERENCES.                                                              !
!      -----------                                                             !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLE.                                                      !

INTEGER,           INTENT(IN) :: IU_PA    !! PARAMETER FILE UNIT
CHARACTER (LEN=*), INTENT(IN) :: FILE_PA  !! PARAMETER FILE NAME
INTEGER,           INTENT(IN) :: IU_SP    !! SPECTRA FILE UNIT
CHARACTER (LEN=*), INTENT(IN) :: FILE_SP  !! SPECTRA FILE NAME
   
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLE                                                           !

INTEGER :: IFAIL                   !! OPEN ERROR CODE
character (len=14) :: cdt_new

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLOSE OLD FILES.                                                      !
!     -------------------                                                      !

IF (FFLAG20) CLOSE (UNIT=IU_PA, STATUS ="KEEP")     !! INTEGRATED PARAMETER FILE.
IF (FFLAG25) CLOSE (UNIT=IU_SP, STATUS ="KEEP")     !! SPECTRA FILE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. OPEN A NEW FILE IF MODEL DATE IS BEFORE END DATE.                     !
!        -------------------------------------------------                     !

IF (CDTPRO.LT.CDATEE) THEN
   CDT_NEW = CDT_OUT
   CALL INCDATE(CDT_NEW, IDEL_OUT)
   IF (FFLAG20) THEN                              !! INTEGRATED PARAMETER FILE.
      CALL OPEN_FILE (IU06, IU_PA, FILE_PA, CDT_NEW, 'UNKNOWN', IFAIL)
      IF (IFAIL.NE.0) CALL ABORT1
   END IF

   IF (FFLAG25) THEN                              !! SPECTRA FILE.
      CALL OPEN_FILE (IU06, IU_SP, FILE_SP, CDT_NEW, 'UNKNOWN', IFAIL)
      IF (IFAIL.NE.0) CALL ABORT1
   END IF
END IF

END SUBROUTINE SAVE_OUTPUT_FILES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_OUTPUT_SET_UP_MODULE

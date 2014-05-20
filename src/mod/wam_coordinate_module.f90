MODULE WAM_COORDINATE_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE PROVIDES THE ROUTINES AND OPERATORS NECESSARY FOR SHERICAL     !
!   COORDINATES REPRESENTATION IN THE WAM MODEL.                               !
!                                                                              !
!   MODEL SECONDS (M_SEC) ARE DEFINED AS                                       !
!           NEAREST INTEGER OF (SECONDS * SCALING FACTOR).                     !
!   THE SCALING FACTOR IS FRAC_SEC = 100                                       !
!                                                                              !
!   ROUTINES ARE PROVIDED TO                                                   !
!           CORRECT WEST EAST BOARDERS,                                        !
!           CHECK WEST EAST PERIODIC,                                          !
!           READ A COORDINATE FROM TEXT STRING AND CONVERT TO M_SEC,           !   
!           CONVERT DEGREES TO M_SEC,                                          !   
!           CONVERT M_SEC TO DEGREES,                                          !   
!           WRITE M_SEC TO A TEXT STRING.                                      !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,    ONLY: IU06

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B.  MODULE DATA.                                                         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

IMPLICIT NONE
PRIVATE

INTEGER, PARAMETER :: KIND_D = 8

INTEGER, PARAMETER :: COOR_UNDEF = 1          !! UNDEFINED MODEL COORDINATE

INTEGER, PARAMETER :: FRAC_SEC = 100          !! M_SEC = SECONDS * FRAC_SEC.
INTEGER, PARAMETER :: M_MINUTE = 60*FRAC_SEC  !! M_SEC = MINUTES * M_MINUTE.
INTEGER, PARAMETER :: M_DEGREE = 60*M_MINUTE  !! M_SEC = DEGREES * M_DEGREE.
INTEGER, PARAMETER :: M_S_PER  = 360*M_DEGREE  !! M_SEC PER CYCLE OF LONGITUDES.

REAL (KIND=KIND_D), PARAMETER :: FRAC_SEC_R = 100.            !! M_SEC = SECONDS * FRAC_SEC.
REAL (KIND=KIND_D), PARAMETER :: M_MINUTE_R = 60.*FRAC_SEC_R  !! M_SEC = MINUTES * M_MINUTE.
REAL (KIND=KIND_D), PARAMETER :: M_DEGREE_R = 60.*M_MINUTE_R  !! M_SEC = DEGREES * M_DEGREE.
REAL (KIND=KIND_D), PARAMETER :: M_S_PER_R = 360.*M_DEGREE_R  !! M_SEC PER CYCLE OF LONGITUDES.

INTEGER, PARAMETER :: LEN_COOR = 13   !! LENGTH OF COORDINATE TEXT STRING

PUBLIC :: COOR_UNDEF, FRAC_SEC, M_MINUTE, M_DEGREE, M_S_PER, LEN_COOR
PUBLIC :: FRAC_SEC_R, M_MINUTE_R, M_DEGREE_R, M_S_PER_R
PUBLIC :: KIND_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTERFACE ADJUST                      !! CORRECT WEST EAST BOARDERS.
   MODULE PROCEDURE ADJUST_D_S        !! DEGREE SCALAR VERSION.
   MODULE PROCEDURE ADJUST_S          !! M_SEC  SCALAR VERSION.
   MODULE PROCEDURE ADJUST_D_A        !! DEGREE ARRAY VERSION. 
   MODULE PROCEDURE ADJUST_A          !! M_SEC  ARRAY VERSION. 
END INTERFACE
PUBLIC ADJUST

INTERFACE CHECK_GRID_DEFINITION
   MODULE PROCEDURE CHECK_GRID_DEFINITION
END INTERFACE
PUBLIC CHECK_GRID_DEFINITION

INTERFACE DEG_TO_M_SEC                !! CONVERTS DEGREES TO M_SEC
   MODULE PROCEDURE DEG_TO_M_SEC_S    !! SCALAR VERSION. 
   MODULE PROCEDURE DEG_TO_M_SEC_A    !! ARRAY VERSION.
END INTERFACE
PUBLIC DEG_TO_M_SEC 

INTERFACE M_SEC_TO_DEG                !! CONVERTS M_SEC TO DEGREES
   MODULE PROCEDURE M_SEC_TO_DEG_S    !! SCALAR VERSION. 
   MODULE PROCEDURE M_SEC_TO_DEG_A    !! ARRAY VERSION.
END INTERFACE
PUBLIC M_SEC_TO_DEG

INTERFACE PERIODIC                    !! WEST EAST PERIODIC.
   MODULE PROCEDURE PERIODIC
END INTERFACE
PUBLIC PERIODIC

INTERFACE READ_COOR_TEXT              !! M_SEC FROM A TEXT STRING.
   MODULE PROCEDURE READ_COOR_TEXT_S  !! SCALAR VERSION. 
   MODULE PROCEDURE READ_COOR_TEXT_A  !! ARRAY VERSION.
END INTERFACE
PUBLIC READ_COOR_TEXT

INTERFACE WRITE_COOR_TEXT             !! WRITE M_SEC TO A TEXT STRING.
   MODULE PROCEDURE WRITE_COOR_TEXT   
END INTERFACE
PUBLIC WRITE_COOR_TEXT

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE ADJUST_D_S (WEST, EAST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   ADJUST_S - ROUTINE TO CORRECT BORDERS OF INTERVALS (SCALAR VERSION).       !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       ADJUSTS INTERVAL BORDERS GIVEN IN DEGREE.                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE INTERVAL BORDERS ARE CHANGED TO FULLFILL:                          !
!         0. .LE. EAST  .AND. EAST .LT. 360. .AND. WEST .LE. EAST              !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

REAL (KIND=KIND_D), INTENT(INOUT) :: WEST  !! LEFT INTERVAL BORDER  [DEG].
REAL (KIND=KIND_D), INTENT(INOUT) :: EAST  !! RIGHT INTERVAL BORDER [DEG].

! ---------------------------------------------------------------------------- !
!                                                                              !
!  1. CORRECT BORDERS.                                                         !
!     ----------------                                                         !

IF (WEST.LT.0.) THEN
   WEST = WEST+360.
   IF (WEST.LT.0.) WEST = WEST+360.
END IF
IF (WEST.GE.360.) THEN
   WEST = WEST-360.
   IF (WEST.GE.360.) WEST = WEST-360.
END IF
IF (EAST.LT.0.) THEN
   EAST = EAST+360.
   IF (EAST.LT.0.) EAST = EAST+360.
END IF
IF (EAST.GE.360.) THEN
   EAST = EAST-360.
   IF (EAST.GE.360.) EAST = EAST-360.
END IF
IF (WEST.GT.EAST) WEST = WEST-360.

END SUBROUTINE ADJUST_D_S

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE ADJUST_S (WEST, EAST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   ADJUST_S - ROUTINE TO CORRECT BORDERS OF INTERVALS (SCALAR VERSION).       !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       ADJUSTS INTERVAL BORDERS GIVEN IN M_SEC.                               !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE INTERVAL BORDERS ARE CHANGED TO FULLFILL:                          !
!         0. .LE. EAST  .AND. EAST .LT. 360. .AND. WEST .LE. EAST              !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

INTEGER, INTENT(INOUT) :: WEST  !! LEFT INTERVAL BORDER  [M_SEC].
INTEGER, INTENT(INOUT) :: EAST  !! RIGHT INTERVAL BORDER [M_SEC].

! ---------------------------------------------------------------------------- !
!                                                                              !
!  1. CORRECT BORDERS.                                                         !
!     ----------------                                                         !

IF (WEST.LT.0) THEN
   WEST = WEST+M_S_PER
   IF (WEST.LT.0) WEST = WEST+M_S_PER
END IF
IF (WEST.GE.M_S_PER) THEN
   WEST = WEST-M_S_PER
   IF (WEST.GE.M_S_PER) WEST = WEST-M_S_PER
END IF
IF (EAST.LT.0) THEN
   EAST = EAST+M_S_PER
   IF (EAST.LT.0) EAST = EAST+M_S_PER
END IF
IF (EAST.GE.M_S_PER) THEN
   EAST = EAST-M_S_PER
   IF (EAST.GE.M_S_PER) EAST = EAST-M_S_PER
END IF
IF (WEST.GT.EAST) WEST = WEST-M_S_PER

END SUBROUTINE ADJUST_S

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE ADJUST_D_A (WEST, EAST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   ADJUST_V - ROUTINE TO CORRECT BORDERS OF INTERVALS (ARRAY VERSION).        !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       ADJUSTS INTERVAL BORDERS GIVEN IN DEGREE.                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE INTERVAL BORDERS ARE CHANGED TO FULLFILL:                          !
!         0. .LE. EAST  .AND. EAST .LT. 360. .AND. WEST .LE. EAST              !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

REAL (KIND=KIND_D), INTENT(INOUT) :: WEST(:) !! LEFT INTERVAL BORDER  [DEG].
REAL (KIND=KIND_D), INTENT(INOUT) :: EAST(:) !! RIGHT INTERVAL BORDER [DEG].

! ---------------------------------------------------------------------------- !
!                                                                              !
!  1. CORRECT BORDERS.                                                         !
!     ----------------                                                         !

WHERE (WEST.LT.0.) WEST = WEST+360.
WHERE (WEST.LT.0.) WEST = WEST+360.
WHERE (WEST.GE.360.) WEST = WEST-360.
WHERE (WEST.GE.360.) WEST = WEST-360.

WHERE (EAST.LT.0.) EAST = EAST+360.
WHERE (EAST.LT.0.) EAST = EAST+360.
WHERE (EAST.GE.360.) EAST = EAST-360.
WHERE (EAST.GE.360.) EAST = EAST-360.

WHERE (WEST.GT.EAST) WEST = WEST-360.

END SUBROUTINE ADJUST_D_A

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE ADJUST_A (WEST, EAST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   ADJUST_A - ROUTINE TO CORRECT BORDERS OF INTERVALS (ARRAY VERSION).        !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!     H.GUNTHER            HZG        JANUARY 2011     M_SEC                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       ADJUSTS INTERVAL BORDERS GIVEN IN DEGREE.                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE INTERVAL BORDERS ARE CHANGED TO FULLFILL:                          !
!         0. .LE. EAST  .AND. EAST .LT. 360. .AND. WEST .LE. EAST              !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

INTEGER, INTENT(INOUT) :: WEST(:)  !! LEFT INTERVAL BORDER  [M_SEC].
INTEGER, INTENT(INOUT) :: EAST(:)  !! RIGHT INTERVAL BORDER [M_SEC].

! ---------------------------------------------------------------------------- !
!                                                                              !
!  1. CORRECT BORDERS.                                                         !
!     ----------------                                                         !

WHERE (WEST.LT.0) WEST = WEST+M_S_PER
WHERE (WEST.LT.0) WEST = WEST+M_S_PER
WHERE (WEST.GE.M_S_PER) WEST = WEST-M_S_PER
WHERE (WEST.GE.M_S_PER) WEST = WEST-M_S_PER

WHERE (EAST.LT.0) EAST = EAST+M_S_PER
WHERE (EAST.LT.0) EAST = EAST+M_S_PER
WHERE (EAST.GE.M_S_PER) EAST = EAST-M_S_PER
WHERE (EAST.GE.M_S_PER) EAST = EAST-M_S_PER

WHERE (WEST.GT.EAST) WEST = WEST-M_S_PER

END SUBROUTINE ADJUST_A

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CHECK_GRID_DEFINITION (WEST, SOUTH, EAST, NORTH,                    &
&                                 DX, DY, NX, NY, ERROR)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   CHECK_GRID_DEFINITION - CHECKS THE GRID COORDINATES.                       !
!                                                                              !
!     H. GUENTHER  HZG  MAY 2002                                               !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       GRID DEFINITIONS ARE CHECKED FOR CONSISTENCY AND COMPLETNESS           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       COMPUTE MISSING VALUES IF POSSIBLE AND CHECK CONSISTENCY.              !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,  INTENT(INOUT) :: WEST    !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER,  INTENT(INOUT) :: SOUTH   !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER,  INTENT(INOUT) :: EAST    !! EAST LONGITUDE OF GRID [M_SEC].
INTEGER,  INTENT(INOUT) :: NORTH   !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER,  INTENT(INOUT) :: DX      !! LONGITUDE INCREMENT [M_SEC].
INTEGER,  INTENT(INOUT) :: DY      !! LATITUDE INCREMENT [M_SEC].
INTEGER,  INTENT(INOUT) :: NX      !! NUMBER OF LONGITUDES.
INTEGER,  INTENT(INOUT) :: NY      !! NUMBER OF LATITUDES.
LOGICAL,  INTENT(OUT)   :: ERROR   !! .TRUE. IF ERROR DETECTED.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     local VARIABLES.                                                         !
!     ----------------                                                         !

character (len=len_coor) :: formtext

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE MISSING VALUES.                                               !
!        -----------------------                                               !

ERROR = .FALSE.

IF (WEST.EQ.COOR_UNDEF) ERROR = .TRUE.
IF (WEST.NE.COOR_UNDEF .AND. EAST.NE.COOR_UNDEF) CALL ADJUST (WEST, EAST)

IF (EAST.EQ.COOR_UNDEF) THEN
   IF (DX.GT.COOR_UNDEF .AND. NX.GT.1) THEN
      EAST = WEST + (NX-1)*DX !! EAST LONGITUDE OF GRID [M_SEC].
      CALL ADJUST (WEST, EAST)
   ELSE
      ERROR = .TRUE.
   END IF
ELSE IF (DX.EQ.COOR_UNDEF) THEN
   IF (NX.GT.1) THEN
      DX = (EAST-WEST+(NX-1)/2)/(NX-1) !! LONGITUDE INCREMENT [M_SEC].
   ELSE
      ERROR = .TRUE.
   END IF
ELSE IF (NX.EQ.-1) THEN
   IF (DX.GT.COOR_UNDEF) THEN
      NX = (EAST-WEST+DX/2)/DX+1   !! NO. OF LONGITUDES.
   ELSE
      ERROR = .TRUE.
   END IF
END IF

IF (SOUTH.EQ.COOR_UNDEF) ERROR = .TRUE.

IF (NORTH.EQ.COOR_UNDEF) THEN
   IF (DY.GT.COOR_UNDEF .AND. NY.GT.1) THEN
      NORTH = SOUTH + (NY-1)*DY !! EAST LONGITUDE OF GRID [M_SEC].
   ELSE
      ERROR = .TRUE.
   END IF
ELSE IF (DY.EQ.COOR_UNDEF) THEN
   IF (NY.GT.1) THEN
      DY = (NORTH-SOUTH+(NY-1)/2)/(NY-1) !! LONGITUDE INCREMENT [M_SEC].
   ELSE
      ERROR = .TRUE.
   END IF
ELSE IF (NY.EQ.-1) THEN
   IF (DY.GT.COOR_UNDEF) THEN
      NY = (NORTH-SOUTH+DY/2)/DY+1   !! NO. OF LONGITUDES.
   ELSE
      ERROR = .TRUE.
   END IF
END IF

IF (ERROR) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *        FATAL ERROR IN SUB. CHECK_GRID_DEFINITION       *'
   WRITE (IU06,*) ' *        =========================================       *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' * GRID SPECIFICATINOS ARE INSUFFICIENT TO SET-UP A GRID. *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' * USER PROVIDED GRID DEFINITIONS ARE:                    *'
   formtext = write_coor_text (coor_undef)
   WRITE (IU06,*) ' *  (-1 AND ', formtext,                                    &
&                                          ' INDICATED UNDEFINED VALUES)     *'
   WRITE (IU06,*) ' *                                                        *'
   formtext = write_coor_text(west)
   WRITE (IU06,*) ' * LONGITUDE             WEST = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (WEST), 'deg'
   formtext = write_coor_text(east)
   WRITE (IU06,*) ' * LONGITUDE             EAST = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (EAST), 'deg'
   formtext = write_coor_text(dx)
   WRITE (IU06,*) ' * LONGITUDE INCREMENT     DX = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (DX), 'deg'
   WRITE (IU06,*) ' * NO. OF LONGITUDES       NX = ', NX
   formtext = write_coor_text(south)
   WRITE (IU06,*) ' * LATITUDE             SOUTH = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (SOUTH), 'deg'
   formtext = write_coor_text(north)
   WRITE (IU06,*) ' * LATITUDE             NORTH = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (NORTH), 'deg'
   formtext = write_coor_text(dy)
   WRITE (IU06,*) ' * LATITUDE  INCREMENT     DY = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (DY), 'deg'
   WRITE (IU06,*) ' * NO. OF LATITUDE         NY = ', NY
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK CONSISTENCY.                                                    !
!        ------------------                                                    !

IF (ABS(SOUTH - NORTH + (NY-1)*DY).GT. NY/2) THEN
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' +      WARNING ERROR IN SUB. CHECK_GRID_DEFINITION       +'
   WRITE (IU06,*) ' +      ===========================================       +'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' + LATITUDE SPECIFICATINOS ARE NOT CONSISTENT WITHIN AN   +'
   WRITE (IU06,*) ' + ACCURACY OF 1/', FRAC_SEC, ' SECONDS'
   formtext = write_coor_text(south)
   WRITE (IU06,*) ' * LATITUDE             SOUTH = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (SOUTH), 'deg'
   formtext = write_coor_text(north)
   WRITE (IU06,*) ' * LATITUDE             NORTH = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (NORTH), 'deg'
   formtext = write_coor_text(dy)
   WRITE (IU06,*) ' * LATITUDE  INCREMENT     DY = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (DY), 'deg'
   WRITE (IU06,*) ' * NO. OF LATITUDE         NY = ', NY
   WRITE (IU06,*) ' + LATITUDE INCREMENT CHANGED:                            +'
    DY = (NORTH-SOUTH+(NY-1)/2)/(NY-1) 
   WRITE (IU06,*) ' +                                                        +'
   formtext = write_coor_text(dy)
   WRITE (IU06,*) ' + NEW LATITUDE INCREMENT  DY = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (DY), 'deg'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' +            PROGRAM WILL CONTINUE                       +'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

IF (ABS(WEST- EAST + (NX-1)*DX) .GT. NX/2 ) THEN
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' +      WARNING ERROR IN SUB. CHECK_GRID_DEFINITION       +'
   WRITE (IU06,*) ' +      ===========================================       +'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' + LONGITUDE SPECIFICATINOS ARE NOT CONSISTENT WITHIN AN  +'
   WRITE (IU06,*) ' + ACCURACY OF 1/', FRAC_SEC, ' SECONDS'
   formtext = write_coor_text(west)
   WRITE (IU06,*) ' * LONGITUDE             WEST = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (WEST), 'deg'
   formtext = write_coor_text(east)
   WRITE (IU06,*) ' * LONGITUDE             EAST = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (EAST), 'deg'
   formtext = write_coor_text(dx)
   WRITE (IU06,*) ' * LONGITUDE INCREMENT     DX = ',  formtext,               &
&                                      ' = ',M_SEC_TO_DEG (DX), 'deg'
   WRITE (IU06,*) ' * NO. OF LONGITUDES       NX = ', NX
   WRITE (IU06,*) ' + LONGITUDE INCREMENT CHANGED:                           +'
   DX = (EAST-WEST+(NX-1)/2)/(NX-1)
   WRITE (IU06,*) ' +                                                        +'
   formtext = write_coor_text(dx)
   WRITE (IU06,*) ' + NEW LONGITUDE INCREMENT DX = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (DX), 'deg'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' +            PROGRAM WILL CONTINUE                       +'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

ERROR = (ABS(SOUTH - NORTH + (NY-1)*DY) .GT. NY/2) .OR.                             &
&       (ABS(WEST- EAST + (NX-1)*DX) .GT. NX/2 ) 

END SUBROUTINE CHECK_GRID_DEFINITION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTEGER FUNCTION DEG_TO_M_SEC_S (DEGREE) RESULT(M_SEC)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   DEG_TO_M_SEC - CONVERTS DEGREES INTO M_SEC.   (SCALAR VERSION)             !
!                                                                              !
!     H.GUNTHER            HZG        JANUARY 2011                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!      TO CONVERT REAL DEGREES INTO INTEGER M_SEC.                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

REAL (KIND=KIND_D), INTENT(IN) :: DEGREE   !! INPUT COORDINATE IN [DEG].

! ---------------------------------------------------------------------------- !

M_SEC = NINT(DEGREE*M_DEGREE_R)

END FUNCTION DEG_TO_M_SEC_S

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

FUNCTION DEG_TO_M_SEC_A (DEGREE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   DEG_TO_M_SEC - CONVERTS DEGREES INTO M_SEC.  (ARRAY VERSION)               !
!                                                                              !
!     H.GUNTHER            HZG        JANUARY 2011                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!      TO CONVERT REAL DEGREES INTO INTEGER M_SEC.                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

REAL (KIND=KIND_D), INTENT(IN) :: DEGREE(:)  !! INPUT COORDINATES IN [DEG].

INTEGER ::  DEG_TO_M_SEC_A(SIZE(DEGREE))     !! OUTPUT COORDINATES IN [M_SEC].

! ---------------------------------------------------------------------------- !

DEG_TO_M_SEC_A = NINT(DEGREE*M_DEGREE_R)

END FUNCTION DEG_TO_M_SEC_A

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

REAL (KIND=KIND_D) FUNCTION M_SEC_TO_DEG_S (M_SEC) RESULT(DEGREE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   M_SEC_TO_DEG - CONVERTS M_SEC INTO DEGREES.   (SCALAR VERSION)             !
!                                                                              !
!     H.GUNTHER            HZG        JANUARY 2011                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!      TO CONVERT INTEGER M_SEC INTO REAL DEGREES.                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

INTEGER, INTENT(IN) :: M_SEC   !! INPUT COORDINATE IN [M_SEC].

! ---------------------------------------------------------------------------- !

DEGREE = DBLE(M_SEC)/M_DEGREE_R

END FUNCTION M_SEC_TO_DEG_S

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

FUNCTION M_SEC_TO_DEG_A (M_SEC)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   M_SEC_TO_DEG - CONVERTS M_SEC INTO DEGREES.  (ARRAY VERSION)               !
!                                                                              !
!     H.GUNTHER            HZG        JANUARY 2011                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!      TO CONVERT INTEGER M_SEC INTO REAL*8 DEGREES.                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

INTEGER, INTENT(IN) :: M_SEC(:)                !! INPUT COORDINATES IN [M_SEC].

REAL (KIND=KIND_D) ::  M_SEC_TO_DEG_A(SIZE(M_SEC)) !! OUTPUT RESULT IN [DEG]

! ---------------------------------------------------------------------------- !

M_SEC_TO_DEG_A = DBLE(M_SEC)/M_DEGREE_R

END FUNCTION M_SEC_TO_DEG_A

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

LOGICAL FUNCTION PERIODIC (WEST, EAST, DEL, N) RESULT(L)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PERIODIC - CHECKS IF GRID IS WEST-EAST PERIODIC.                           !
!                                                                              !
!     H.GUNTHER            HZG        JANUARY 2011                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!      TO CHECK IF GRID IS WEST-EAST PERIODIC.                                 !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE RESULT IS TRUE IF  WEST + M_S_PER = EAST + DEL                     !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

INTEGER, INTENT(IN) :: WEST  !! WESTERN COORDINATE IN [M_SEC].
INTEGER, INTENT(IN) :: EAST  !! EASTERN COORDINATE IN [M_SEC].
INTEGER, INTENT(IN) :: DEL   !! LONGITUDE STEP IN [M_SEC].
INTEGER, INTENT(IN) :: N     !! NUMBER OF LONGITUDES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!  1. CORRECT BORDERS.                                                         !
!     ----------------                                                         !

L = ABS(WEST - EAST - DEL + M_S_PER) .LT. N/2


END FUNCTION PERIODIC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTEGER FUNCTION READ_COOR_TEXT_S (TEXT) RESULT (M_SEC)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     READ_COOR_TEXT_S - READ A COORDINATE FROM A CHARACTER STRING.            !
!                       (SCALAR VERSION)                                       !
!                                                                              !
!     HEINZ GUNTHER    HZG    JANUARY 2011                                     !
!                                                                              !
!     PURPOSE                                                                  !
!     -------                                                                  !
!                                                                              !
!     TO CONVERT A COORDINATE GIVEN AS CHARACTER STRING INTO MODEL             !
!     REPRESENTATION.                                                          !
!                                                                              !
!     METHOD                                                                   !
!     ------                                                                   !
!                                                                              !
!     READ A COORDINATE FROM CHARACTER VARIABLE FROMATTED AS                   !
!         +DDD:MM:SS.SS OR  -DDD:MM:SS.SS                                      !
!     OR                                                                       !
!         '(F13.8)'  REAL*8 DEGREES.                                             !
!                                                                              !
!     THE RESULT IS THE COORDINATE IN INTEGER [M_SEC].                         !        
!                                                                              !
!                                                                              !
!                                                                              !
!     REFERENCES                                                               !
!     ----------                                                               !
!                                                                              !
!        NONE.                                                                 !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=LEN_COOR), INTENT(IN)  :: TEXT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IOS 
INTEGER :: DEGREE, MINUTE
REAL (KIND=KIND_D) :: SECOND

! ---------------------------------------------------------------------------- !

IF (TEXT.EQ.' ') THEN
   M_SEC = COOR_UNDEF
   RETURN
END IF

READ (TEXT,'(F13.8)',IOSTAT=IOS) SECOND
IF (IOS.EQ.0) THEN
   M_SEC = NINT(SECOND*M_DEGREE_R)
ELSE
   READ (TEXT,'(1X,I3,1X,I2,1X,F5.2)',IOSTAT=IOS) DEGREE, MINUTE, SECOND
   IF (IOS.EQ.0) THEN
      M_SEC = DEGREE*M_DEGREE + MINUTE*M_MINUTE + NINT(SECOND*FRAC_SEC_R)
      IF (TEXT(1:1).EQ.'-') M_SEC= -M_SEC
   ELSE 
      M_SEC = COOR_UNDEF   
   END IF
END IF

END FUNCTION READ_COOR_TEXT_S

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

FUNCTION READ_COOR_TEXT_A (TEXT) 

! ---------------------------------------------------------------------------- !
!                                                                              !
!     READ_COOR_TEXT_A - READ COORDINATES FROM A CHARACTER STRINGS.            !
!                       (ARRAY VERSION)                                        !
!                                                                              !
!     HEINZ GUNTHER    HZG    JANUARY 2011                                     !
!                                                                              !
!     PURPOSE                                                                  !
!     -------                                                                  !
!                                                                              !
!     TO CONVERT AN ARRAY OF COORDINATES GIVEN AS CHARACTER STRINGS INTO MODEL !
!     REPRESENTATIONS (M_SEC).                                                 !
!                                                                              !
!     METHOD                                                                   !
!     ------                                                                   !
!                                                                              !
!     READS COORDINATES FROM CHARACTER VARIABLE FROMATTED AS                   !
!         +DDD:MM:SS.SS OR  -DDD:MM:SS.SS                                      !
!     OR                                                                       !
!         '(F13.8)'  REAL*8 DEGREES.                                           !
!                                                                              !
!     THE RESULT ARE THE COORDINATES IN INTEGER [M_SEC].                         !        
!                                                                              !
!                                                                              !
!                                                                              !
!     REFERENCES                                                               !
!     ----------                                                               !
!                                                                              !
!        NONE.                                                                 !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=LEN_COOR),           INTENT(IN)  :: TEXT(:)

INTEGER :: READ_COOR_TEXT_A(SIZE(TEXT))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: I

! ---------------------------------------------------------------------------- !

DO I=1,SIZE(TEXT)
   READ_COOR_TEXT_A(I) =  READ_COOR_TEXT_S (TEXT(I))
END DO

END FUNCTION READ_COOR_TEXT_A

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CHARACTER (LEN=LEN_COOR) FUNCTION WRITE_COOR_TEXT (M_SEC) RESULT (TEXT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     TEXT_WRITE_COOR - WRITES A MODEL COORDINATE TO A CHARACTER STRING.       !
!                                                                              !
!     HEINZ GUNTHER    HZG    JANUARY 2011                                     !
!                                                                              !
!     PURPOSE                                                                  !
!     -------                                                                  !
!                                                                              !
!     TO CONVERT A MODEL COORDINATE GIVEN AS M_SEC INTO A TEXT STRING.         !
!                                                                              !
!     METHOD                                                                   !
!     ------                                                                   !
!                                                                              !
!     WRITES A COORDINATE TO A CHARACTER VARIABLE FROMATTED AS                 !
!         +DDD:MM:SS.SS OR  -DDD:MM:SS.SS                                      !
!                                                                              !
!     REFERENCES                                                               !
!     ----------                                                               !
!                                                                              !
!        NONE.                                                                 !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)  :: M_SEC

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER           :: DEGREE, MINUTE, SEC_H
REAL              :: SECOND
CHARACTER (LEN=1) :: SIG

! ---------------------------------------------------------------------------- !

IF (M_SEC.LT.0) THEN
   SIG = '-'
   SEC_H = -M_SEC
ELSE
   SIG = '+'
   SEC_H = M_SEC
END IF
DEGREE = SEC_H/M_DEGREE
SEC_H = SEC_H-DEGREE*M_DEGREE
MINUTE = SEC_H/M_MINUTE
SEC_H = SEC_H-MINUTE*M_MINUTE
SECOND = REAL(SEC_H)/REAL(FRAC_SEC)

WRITE (TEXT,'(A1,I3.3,A1,I2.2,A1,F5.2)') SIG, DEGREE, ':', MINUTE, ':', SECOND
IF (TEXT(9:9) .EQ. ' ') TEXT(9:9)= '0'

END FUNCTION WRITE_COOR_TEXT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_COORDINATE_MODULE

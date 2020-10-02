MODULE WAM_FILE_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL STANDARD UNITS, FILE NAMES, AND FILE IDENTIFIER   !
!   USED BY PREPROC, PRESET AND CHIEF.                                         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

INTEGER, PRIVATE :: LENT

! ---------------------------------------------------------------------------- !
!                                                                              !
!    1.0 PRINTER OUTPUT UNIT, FILE NAME AND TEST FLAGS.                        !
!        ----------------------------------------------                        !

INTEGER            :: IU06   = 66            !! UNIT FOR PRINTER OUTPUT.
CHARACTER (LEN=80) :: FILE06 = 'WAM_Prot'

INTEGER            :: ITEST = 0  !! TEST OUTPUT LEVEL:
                                 !!   .LE. 0  NO OUTPUT
                                 !!   .GE. I  OUTPUT UP TO LEVEL I

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.0 USER INPUT FILE UNIT AND NAME.                                       !
!         ------------------------------                                       !

INTEGER            :: IU05   = 55           !! INPUT OF USER DATA.
CHARACTER (LEN=80) :: FILE05 = 'WAM_User'   !!   (SEE SUB READ_WAM_USER).

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3.0 INPUT AND OUTPUT UNIT AND FILE NAMES.                                 !
!        -------------------------------------                                 !
!                                                                              !
!         FILE UNITS AND NAMES CAN BE REDEFINED BY SET_XXX_FILE SUBROUTINES.   !

INTEGER            :: IU01   =  1           !! INPUT OF WIND FIELDS.
                                            !! (SEE SUB READ_WIND_INPUT).
CHARACTER (LEN=80) :: FILE01 = 'WIND_INPUT' !! WIND INPUT FILE NAME

INTEGER            :: IU02   =  2           !! INPUT OF BOUDARY VALUES
                                            !! FROM A PREVIOUS COARSE GRID WAM
CHARACTER (LEN=80) :: FILE02 = ' '          !! BOUNDARY INPUT FILE IDENTIFIER

INTEGER            :: IU03   =  3           !! INPUT OF AN ICE FIELD
                                            !! (SEE SUB READ_ICE_INPUT).
CHARACTER (LEN=80) :: FILE03 = ' '          !! ICE INPUT FILE NAME

INTEGER            :: IU07   =  7           !! PREPROC OUTPUT FILE.
CHARACTER (LEN=80) :: FILE07 = 'Grid_Info'  !!  (INPUT TO CHIEF).

INTEGER            :: IU08   = 8            !! INPUT OF TOPOGRAPHIC DATA.
CHARACTER (LEN=80) :: FILE08 = ' '          !!   (SEE SUB READ_TOPOGRAPHY).

INTEGER            :: IU09   = 9            !! INPUT OF CURRENT DATA.
CHARACTER (LEN=80) :: FILE09 = ' '          !!   (SEE SUB READ_CURRENT).

INTEGER            :: IU10   = 10           !! OUTPUT FILE FROM THE COARSE.
                                            !! GRID PREPROC.
CHARACTER (LEN=80) :: FILE10 = ' '          !!   (SEE SUB MAKE_FINE_BOUNDARY).

! ---------------------------------------------------------------------------- !
!                                                                              !
!    4.0 STANDARD OUTPUT UNITS AND FILE IDENTIFIER.                            !
!        ------------------------------------------                            !
!                                                                              !
!         FILE UNITS AND NAMES CAN BE REDEFINED BY SET_XXX_FILE SUBROUTINES.   !
!         THE FULL FILE NAME IS THE FILE IDENTIFIER EXTENDED BY THE FILE DATE. !
!         SEE OPEN_FILE IN WAM_GENERAL_MODULE FOR DETAILS.                     !

INTEGER            :: IU17   = 17       !! RESTART FILE.
CHARACTER (LEN=80) :: FILE17 = 'BLS'    !! RESTART FILE IDENTIFIER.

INTEGER            :: IU19 = 70         !! OUTPUT UNIT FOR BOUNDARY VALUES IF
                                        !! THIS IS A FINE GRID RUN.
CHARACTER (LEN=80) :: FILE19 = 'CBO'    !! FILE IDENTIFIER

INTEGER            :: IU20 = 20         !! OUTPUT UNIT FOR INTEGRATED PARAMETER
CHARACTER (LEN=80) :: FILE20 = 'MAP'    !! FILE IDENTIFIER

INTEGER            :: IU25 = 25         !! OUTPUT UNIT FOR SPECTRA AT CERTAIN
                                        !! GRID POINTS.
CHARACTER (LEN=80) :: FILE25 = 'OUT'    !! FILE IDENTIFIER

INTEGER            :: IU27 = 27         !! OUTPUT UNIT FOR RADIATION STRESS.
CHARACTER*80       :: FILE27 = 'RAD'    !! FILE IDENTIFIER

integer            :: iu67 = 67         !! ready file for output of integrated
                                        !! parameters
character (len=128) :: wpath
character (len=  3) :: area
				       
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE SET_USER_FILE
   MODULE PROCEDURE SET_USER_FILE
END INTERFACE
PUBLIC SET_USER_FILE

INTERFACE SET_TEST_OPTION
   MODULE PROCEDURE SET_TEST_OPTION
END INTERFACE
PUBLIC SET_TEST_OPTION

INTERFACE SET_PROTOCOL_FILE
   MODULE PROCEDURE SET_PROTOCOL_FILE
END INTERFACE
PUBLIC SET_PROTOCOL_FILE

INTERFACE SET_SPECTRA_FILE
   MODULE PROCEDURE SET_SPECTRA_FILE
END INTERFACE
PUBLIC SET_SPECTRA_FILE

INTERFACE SET_MAP_FILE
   MODULE PROCEDURE SET_MAP_FILE
END INTERFACE
PUBLIC SET_MAP_FILE

INTERFACE SET_RESTART_FILE
   MODULE PROCEDURE SET_RESTART_FILE
END INTERFACE
PUBLIC SET_RESTART_FILE

INTERFACE SET_ICE_FILE
   MODULE PROCEDURE SET_ICE_FILE
END INTERFACE
PUBLIC SET_ICE_FILE

INTERFACE SET_PREPROC_FILE
   MODULE PROCEDURE SET_PREPROC_FILE
END INTERFACE
PUBLIC SET_PREPROC_FILE

INTERFACE SET_B_INPUT_FILE
   MODULE PROCEDURE SET_B_INPUT_FILE
END INTERFACE
PUBLIC SET_B_INPUT_FILE

INTERFACE SET_B_OUTPUT_FILE
   MODULE PROCEDURE SET_B_OUTPUT_FILE
END INTERFACE
PUBLIC SET_B_OUTPUT_FILE

INTERFACE SET_WIND_FILE
   MODULE PROCEDURE SET_WIND_FILE
END INTERFACE
PUBLIC SET_WIND_FILE

INTERFACE SET_TOPO_FILE
   MODULE PROCEDURE SET_TOPO_FILE
END INTERFACE
PUBLIC SET_TOPO_FILE

INTERFACE SET_CURRENT_FILE
   MODULE PROCEDURE SET_CURRENT_FILE
END INTERFACE
PUBLIC SET_CURRENT_FILE

INTERFACE SET_C_PREPROC_FILE
   MODULE PROCEDURE SET_C_PREPROC_FILE
END INTERFACE
PUBLIC SET_C_PREPROC_FILE

INTERFACE PRINT_FILE_STATUS
   MODULE PROCEDURE PRINT_FILE_STATUS
END INTERFACE
PUBLIC PRINT_FILE_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

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

SUBROUTINE SET_USER_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE05 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU05 = UNIT

END SUBROUTINE SET_USER_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_PROTOCOL_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE06 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU06 = UNIT

END SUBROUTINE SET_PROTOCOL_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_SPECTRA_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE25 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU25 = UNIT

END SUBROUTINE SET_SPECTRA_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_MAP_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE20 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU20 = UNIT

END SUBROUTINE SET_MAP_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_RESTART_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE17 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU17 = UNIT

END SUBROUTINE SET_RESTART_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_ICE_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE03 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU03 = UNIT

END SUBROUTINE SET_ICE_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_PREPROC_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE07 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU07 = UNIT

END SUBROUTINE SET_PREPROC_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_B_INPUT_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE02 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU02 = UNIT

END SUBROUTINE SET_B_INPUT_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_B_OUTPUT_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE19 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU19 = UNIT

END SUBROUTINE SET_B_OUTPUT_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_WIND_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE01 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU01 = UNIT

END SUBROUTINE SET_WIND_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_TOPO_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE08 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU08 = UNIT

END SUBROUTINE SET_TOPO_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_CURRENT_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE09 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU09 = UNIT

END SUBROUTINE SET_CURRENT_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_C_PREPROC_FILE (NAME, UNIT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER, INTENT(IN)           :: NAME*(*)  !! FILE NAME
INTEGER,   INTENT(IN), OPTIONAL :: UNIT      !! UNIT NUMBER

! ---------------------------------------------------------------------------- !

LENT = LEN_TRIM(NAME)
IF (LENT.GT.0) FILE10 = NAME(1:LENT)
IF (PRESENT(UNIT)) IU10 = UNIT

END SUBROUTINE SET_C_PREPROC_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_FILE_STATUS

LOGICAL       :: OPND

WRITE(IU06,*) '  '
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '               FILE MODULE STATUS:'
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '  '
WRITE(IU06,*) ' STANDARD FILES: '
WRITE(IU06,*) '  '

LENT = LEN_TRIM(FILE05)
WRITE (IU06,'(''  STANDARD USER INPUT FILE ... UNIT:'',I3,'', NAME: '',A)')    &
&                                                         IU05, FILE05(1:LENT)

INQUIRE (UNIT=IU05, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED'
END IF

LENT = LEN_TRIM(FILE06)
WRITE (IU06,'(''  PROTOCOLL FILE ............. UNIT:'',I3,'', NAME: '',A)')    &
&                                                         IU06, FILE06(1:LENT)
INQUIRE (UNIT=IU06, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED'
END IF

WRITE(IU06,*) '  '
WRITE(IU06,*) ' INPUT FILES: '
WRITE(IU06,*) '  '

LENT = LEN_TRIM(FILE01)
WRITE (IU06,'(''  WIND DATA FILE ............. UNIT:'',I3,'', NAME: '',A)')    &
&                                                         IU01, FILE01(1:LENT)
INQUIRE (UNIT=IU01, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED'
END IF

LENT = LEN_TRIM(FILE02)
WRITE (IU06,'(''  BOUNDARY INPUT SPECTRA FILE  UNIT:'',I3,'', ID .: '',A)')    &
&                                                         IU02, FILE02(1:LENT)
INQUIRE (UNIT=IU02, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED'
END IF

LENT = LEN_TRIM(FILE03)
WRITE (IU06,'(''  ICE DATA FILE .............. UNIT:'',I3,'', NAME: '',A)')    &
&                                                         IU03, FILE03(1:LENT)
INQUIRE (UNIT=IU03, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED'
END IF

LENT = LEN_TRIM(FILE07)
WRITE (IU06,'(''  PREPROC FILE ............... UNIT:'',I3,'', NAME: '',A)')    &
&                                                         IU07, FILE07(1:LENT)
INQUIRE (UNIT=IU07, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED'
END IF

LENT = LEN_TRIM(FILE08)
WRITE (IU06,'(''  TOPO DATA FILE ............. UNIT:'',I3,'', NAME: '',A)')    &
&                                                         IU08, FILE08(1:LENT)
INQUIRE (UNIT=IU08, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED'
END IF

LENT = LEN_TRIM(FILE09)
WRITE (IU06,'(''  CURRENT DATA FILE .......... UNIT:'',I3,'', NAME: '',A)')    &
&                                                         IU09, FILE09(1:LENT)
INQUIRE (UNIT=IU09, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED'
END IF

WRITE(IU06,*) '  '
WRITE(IU06,*) ' RESTART FILE: '
WRITE(IU06,*) '  '

LENT = LEN_TRIM(FILE17)
WRITE (IU06,'(''  RESTART FILE ............... UNIT:'',I3,'', ID .: '',A)')    &
&                                                         IU17, FILE17(1:LENT)
INQUIRE (UNIT=IU17, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) ' A FILE IS PRESENTLY ASSIGNED.'
ELSE
   WRITE(IU06,*) '  A FILE IS PRESENTLY NOT ASSIGNED.'
END IF

WRITE(IU06,*) '  '
WRITE(IU06,*) ' BOUNDARY OUTPUT SPECTRA FILE: '
WRITE(IU06,*) '  '

LENT = LEN_TRIM(FILE19)
WRITE (IU06,'(''  BOUNDARY OUTPUT SPECTRA FILE UNIT:'',I3,'', ID .: '',A)')    &
&                                                         IU19, FILE19(1:LENT)
INQUIRE (UNIT=IU19, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED.'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED.'
END IF

WRITE(IU06,*) '  '
WRITE(IU06,*) ' OUTPUT FILES: '
WRITE(IU06,*) '  '

LENT = LEN_TRIM(FILE20)
WRITE (IU06,'(''  INTEGRATED PARAMETER FILE .. UNIT:'',I3,'', ID .: '',A)')    &
&                                                         IU20, FILE20(1:LENT)
INQUIRE (UNIT=IU20, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED.'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED.'
END IF

LENT = LEN_TRIM(FILE25)
WRITE (IU06,'(''  SPECTRA FILE ............... UNIT:'',I3,'', ID .: '',A)')    &
&                                                         IU25, FILE25(1:LENT)
INQUIRE (UNIT=IU25, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED.'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED.'
END IF

WRITE(IU06,*) '  '
WRITE(IU06,*) ' RADIATION OUTPUT FILE: '
WRITE(IU06,*) '  '

LENT = LEN_TRIM(FILE27)
WRITE (IU06,'(''  RADIATION OUTPUT FILE ...... UNIT:'',I3,'', ID .: '',A)')    &
&                                                         IU27, FILE27(1:LENT)
INQUIRE (UNIT=IU27, OPENED=OPND)
IF (OPND) THEN
   WRITE(IU06,*) '  FILE IS ASSIGNED.'
ELSE
   WRITE(IU06,*) '  FILE IS NOT ASSIGNED.'
END IF

WRITE(IU06,*) '  '
WRITE(IU06,*) ' TEST OUTPUT LEVEL IS ..............: ', ITEST
WRITE(IU06,*) '  '

END SUBROUTINE PRINT_FILE_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_FILE_MODULE

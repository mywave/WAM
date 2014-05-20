SUBROUTINE READ_PREPROC_USER

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READ_PREPROC_USER - ROUTINE TO READ USER INPUT FOR PREPROC.                !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO READ USER INPUT OF PROGRAM PREPROC AND TRANSFER INFORMATION         !
!       TO MODULES.                                                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       FILE05 IS ASSIGNED TO IU05, DATA ARE READ AND TRANSFERED TO MODULE     !
!       BY CALLS OF SET_XXX SUBROUTINES.                                       !
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

USE WAM_GENERAL_MODULE,  ONLY:  &
&       ABORT1                    !! TERMINATES PROCESSING.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,   ONLY: IU06, IU05, FILE05

USE PREPROC_USER_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

IMPLICIT NONE

CHARACTER (LEN=80) :: LINE
INTEGER            :: L, IOS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. OPEN FILE AND READ HEADER.                                            !
!        --------------------------                                            !

CALL CLEAR_PREPROC_USER_MODULE

L = LEN_TRIM(FILE05)
IOS = 0
OPEN (UNIT=IU05, FILE=FILE05(1:L), FORM='FORMATTED', STATUS='OLD', IOSTAT=IOS)

IF (IOS.EQ.0) THEN
   CALL READ_PREPROC_NAMELIST (1, IOS)
   IF (IOS.EQ.0) THEN
      CALL SET_PREPROC_USER_PARAMETER
      CLOSE (UNIT=IU05, STATUS="KEEP")
      RETURN
   ELSE
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +   INFORMATION FROM SUB. READ_PREPROC_USER        +'
      WRITE (IU06,*) ' +   =======================================        +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' + THE PREPROC USER FILE COULD BE OPENED.           +'
      WRITE (IU06,*) ' + BUT DOES NOT CONTAIN THE PREPROC NAMELIST.       +'
      WRITE (IU06,*) ' + PREPROC USER FILE NAME IS FILE05 = ', TRIM(FILE05)  
      WRITE (IU06,*) ' +         UNIT IS         IU05 = ', IU05
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' *   PROGRAM TRIES FOR FIXED FORMATTED TEXT FILE    +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      REWIND (UNIT=IU05)
   END IF
ELSE
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +     INFORMATION FROM SUB. READ_PREPROC_USER      +'
   WRITE (IU06,*) ' +     =======================================      +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + THE WAM USER FILE COULD NOT BE OPENED.           +'
   WRITE (IU06,*) ' + PREPROC USER FILE NAME IS FILE05 = ', TRIM(FILE05)  
   WRITE (IU06,*) ' +         UNIT IS         IU05 = ', IU05
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' *  PROGRAM TRIES FOR NAMELIST FROM STANDARD INPUT  +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'

   CALL READ_PREPROC_NAMELIST (0, IOS)
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *    FATAL ERROR IN SUB. READ_PREPROC_NAMELIST     *'
      WRITE (IU06,*) ' *    =========================================     *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * READ ERROR ON NAMLIST FROM STANDARD INPUT        *'  
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   END IF
   CALL SET_PREPROC_USER_PARAMETER
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. GRID NAME.                                                            !
!        -----------                                                           !

CALL F_NEW_DATA
HEADER = LINE(2:80)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. FREQUENCY AND DIRECTION GRID DEFINITIONS.                             !
!        -----------------------------------------                             !

CALL F_NEW_DATA
IF (LINE( 2: 6).NE.' ') READ (LINE( 2: 6),'(I5)',IOSTAT=IOS) ML
IF (IOS.NE.0) CALL ERROR_MESSAGE ('ML')
IF (LINE( 8:12).NE.' ') READ (LINE( 8:12),'(I5)',IOSTAT=IOS) KL
IF (IOS.NE.0) CALL ERROR_MESSAGE ('KL')
IF (LINE(14:23).NE.' ') READ (LINE(14:23),'(F10.8)',IOSTAT=IOS) FR1
IF (IOS.NE.0) CALL ERROR_MESSAGE('FR1')

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. OUTPUT GRID DEFINITIONS.                                              !
!        ------------------------                                              !

CALL F_NEW_DATA
REDUCED_GRID = (SCAN(LINE( 2:11),'T').GT.0 .OR. SCAN(LINE( 2:11),'t').GT.0) 

CALL F_NEW_DATA
IF (LINE( 2:14).NE.' ') AMOSOP = LINE( 2:14)
IF (LINE(16:28).NE.' ') AMONOP = LINE(16:28)
IF (LINE(30:42).NE.' ') AMOWEP = LINE(30:42)
IF (LINE(44:56).NE.' ') AMOEAP = LINE(44:56)

CALL F_NEW_DATA

IF (LINE(2:14).NE.' ') THEN
   READ (LINE(2:14),'(I13)',IOSTAT=IOS)  NX
   IF (IOS.NE.0) CALL ERROR_MESSAGE ('OUTPUT GRID DEFINITIONS: NX')
END IF
IF (LINE(16:28).NE.' ') THEN
   READ (LINE(16:28),'(I13)',IOSTAT=IOS)  NY
   IF (IOS.NE.0) CALL ERROR_MESSAGE ('OUTPUT GRID DEFINITIONS: NY')
END IF
IF (LINE(30:42).NE.' ') XDELLA = LINE(30:42)
IF (LINE(44:56).NE.' ') XDELLO = LINE(44:56)

CALL F_NEW_DATA
READ (LINE,'(1X,F20.3)',IOSTAT=IOS) LAND
IF (IOS.NE.0) CALL ERROR_MESSAGE ('LAND')

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. OUTPUT GRID CORRECTIONS.                                              !
!        ------------------------                                              !

L = 0
DO
   CALL F_NEW_DATA
   IF (LINE(2:4).EQ.'END') EXIT
   L = L + 1
   IF (L.LE.NOUT) THEN
      IF (LINE( 2:14).NE.' ') XOUTS(L) = LINE( 2:14)
      IF (LINE(16:28).NE.' ') XOUTN(L) = LINE(16:28)
      IF (LINE(30:42).NE.' ') XOUTW(L) = LINE(30:42)
      IF (LINE(44:56).NE.' ') XOUTE(L) = LINE(44:56)
      READ (LINE(58:67),'(F10.2)',IOSTAT=IOS) XOUTD(L)
      IF (IOS.NE.0) CALL ERROR_MESSAGE ('OUTPUT GRID CORRECTIONS')
   END IF
END DO

IF (L.GT.NOUT) THEN
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +     WARNING ERROR IN SUB. READ_PREPROC_USER      +'
   WRITE (IU06,*) ' +     =======================================      +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + NUMBER OF AREAS TO BE CORRECTED EXCEEDS          +'
   WRITE (IU06,*) ' + DIMENSION NOUT = ', NOUT
   WRITE (IU06,*) ' + NOUT IS DEFINED IN PREPROC_USER_MODULE           +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + THE FIRST NOUT AREAS ARE ONLY USED.              +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. TEST OPTION.                                                          !
!        ------------                                                          !

CALL F_NEW_DATA
READ (LINE,'(1X,I8)',IOSTAT=IOS) ITEST
IF (IOS.NE.0) CALL ERROR_MESSAGE ('ITEST')

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. NESTED GRID INFORMATION.                                              !
!        ------------------------                                              !

!     6.1 COARSE GRID OPTION.                                                  !

L = 0
DO
   CALL F_NEW_DATA
   IF (LINE(2:4).EQ.'END') EXIT
   L = L+1
   IF (L.LE.N_NEST) THEN
      IF (LINE( 2:14).NE.' ') AMOSOC(L) = LINE( 2:14)
      IF (LINE(16:28).NE.' ') AMONOC(L) = LINE(16:28)
      IF (LINE(30:42).NE.' ') AMOWEC(L) = LINE(30:42)
      IF (LINE(44:56).NE.' ') AMOEAC(L) = LINE(44:56)
      IF (LINE(58:77).NE.' ') NEST_NAME(L) = LINE(58:77)

      READ (LINE(78:80),'(i3)', IOSTAT=IOS) nestcode(l)
      IF (IOS.NE.0) CALL ERROR_MESSAGE ('NEST DEFINITIONS')
   END IF
END DO
IF (L.GT.N_NEST) THEN
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +     WARNING ERROR IN SUB. READ_PREPROC_USER      +'
   WRITE (IU06,*) ' +     =======================================      +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + NUMBER OF NESTS EXCEEDS DIMENSION N_NEST = ', N_NEST
   WRITE (IU06,*) ' + N_NEST IS DEFINED IN PREPROC_USER_MODULE.        +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + THE FIRST N_NEST AREAS ARE ONLY USED.            +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

CALL F_NEW_DATA

!     6.2 FINE GRID OPTION.                                                  !

PREPROC_C_INPUT_FILE_NAME  = TRIM(LINE(2:80)) !! COARSE PREPROC FILE NAME

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. INPUT AND OUTPUT FILE NAMES.                                          !
!        ----------------------------                                          !

CALL F_NEW_DATA
TOPO_INPUT_FILE_NAME = TRIM(LINE(2:80))     !! DEPTH DATA FILE NAME

CALL F_NEW_DATA
PREPROC_OUTPUT_FILE_NAME = TRIM(LINE(2:80))  !! PREPROC OUTPUT FILE NAME

! ---------------------------------------------------------------------------- !
!                                                                              !
!     8. RETURN.                                                               !
!        -------                                                               !

WRITE (IU06,*) '     SUB. READ_PREPROC_USER SUCCESSFULLY COMPLETED. ' 

CALL SET_PREPROC_USER_PARAMETER

CLOSE (UNIT=IU05)

RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     9. INCLUDED FUNCTIONS.                                                   !
!        -------------------                                                   !

CONTAINS

   SUBROUTINE F_NEW_DATA      !! FIND A NEW RECORD STARTING WITHOUT 'C'

   IOS = 0
   LINE(1:1) = 'C'
   DO WHILE (LINE(1:1).EQ.'C' .OR. LINE(1:1).EQ.'c')
      READ (IU05, '(A)',IOSTAT=IOS) LINE
      IF (IOS.NE.0) THEN
         WRITE (IU06,*) ' ****************************************************'
         WRITE (IU06,*) ' *                                                  *'
         WRITE (IU06,*) ' *     FATAL ERROR IN SUB. READ_PREPROC_USER        *'
         WRITE (IU06,*) ' *     =====================================        *'
         WRITE (IU06,*) ' *                                                  *'
         WRITE (IU06,*) ' * READ ERROR ON INPUT FILE:                        *'
         WRITE (IU06,*) ' * LAST LINE READ IS     LINE = ', LINE
         WRITE (IU06,*) ' * ERROR NO. IS        IOSTAT = ', IOS
         WRITE (IU06,*) ' *                                                  *'
         WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
         WRITE (IU06,*) ' *                                                  *'
         WRITE (IU06,*) ' ****************************************************'
         CALL ABORT1
      END IF
   END DO

   END SUBROUTINE F_NEW_DATA

   SUBROUTINE ERROR_MESSAGE (MESSAGE)
   CHARACTER (LEN=*) :: MESSAGE
      WRITE(IU06,*) ' *********************************************'
      WRITE(IU06,*) ' *                                           *'
      WRITE(IU06,*) ' *   FATAL ERROR IN SUB. READ_PREPROC_USER   *'
      WRITE(IU06,*) ' *   =====================================   *'
      WRITE(IU06,*) ' *  READ ERROR ON CHARACTER STRING           *'
      WRITE(IU06,*) ' *                     IOSTAT = ', IOS
      WRITE(IU06,*) ' *  CHARACTER STRING IS  LINE = ', LINE
      WRITE(IU06,*) ' *  PROGRAM WANTS TO READ: ', MESSAGE
      WRITE(IU06,*) ' *                                           *'
      WRITE(IU06,*) ' *    PROGRAM ABORTS  PROGRAM ABORTS         *'
      WRITE(IU06,*) ' *                                           *'
      WRITE(IU06,*) '*********************************************'
      CALL ABORT1
   END SUBROUTINE ERROR_MESSAGE

END SUBROUTINE READ_PREPROC_USER

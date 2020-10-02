SUBROUTINE READ_SPECTRA_USER

! ---------------------------------------------------------------------------- !
!                                                                              !
!    READ_SPECTRA_USER - ROUTINE TO READ USER INPUT OF PROG PRINT_SPECTRA.     !
!                                                                              !
!     H. GUNTHER     GKSS/ECMWF  NOVEMBER 1989                                 !
!                    HZG         DECEMBER 2010      RE-ORGANISED               !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        READ USER INPUT CONCERNING PERIOD OF INTEREST, TIMESTEPS AND OPTIONS. !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        ALL USER PARAMETERS, WHICH CAN OR MUST BE DEFINED BY THE USER, ARE    !
!        COMBINED IN THE PRINT_NAMELIST DEFINED IN  WAM_PRINT_USER_MODULE.     !
!        DEFAULT VALUES ARE DEFINED BY SUB. CLEAR_PRINT_USER_MODULE, WHICH IS  !
!        CONTAINED IN WAM_PRINT_USER_MODULE.                                   !
!        THE SUB. TRIES TO OPEN THE STANDARD USER INPUT FILE.                  !
!        IF THE FILES EXISTS IT FIRST TRIES TO READ THE NAMELIST, AND          !
!        IF THE READING FAILS A FAORMATTED READ IS USED.                       !
!        IF THE FILES DOES NOT EXIST, THE NAMELIST IS READ FROM STANDARD INPUT.!
!                                                                              !
!        USER INFORMATION IN THE FORMATTED FILE IS READ WITH THE               !
!        PRESUMPTIONS THAT:                                                    !
!         1. EVERY LINE STARTING WITH 'C' OR 'c' IS A COMMENT LINE             !
!         2. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'              !
!         3. IF VALUES ARE NOT SPECIFIED DEFAULT VALUES WILL BE USED.          !
!         4. DEFAULT VALUES ARE DEFINED BY SUB. CLEAR_PRINT_USER_MODULE        !
!            CONTAINED IN WAM_PRINT_USER_MODULE.                               !
!                                                                              !
!        FINALLY ALL USER PARAMETERS ARE TRANSFERED TO THE WAM_PINT_MODULE BY  !
!        SUB: SET_PRINT_USER_PARAMETER                                         !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!     -----------                                                              !

USE WAM_GENERAL_MODULE,     ONLY:  &
&       ABORT1                          !! TERMINATES PROCESSING.

USE WAM_PRINT_USER_MODULE,  ONLY:  &
&       CLEAR_PRINT_USER_MODULE,   &    !! DEFINES DEFAULT USER PARAMETERS.
&       READ_PRINT_NAMELIST,       &    !! READS THE NAMELIST
&       SET_PRINT_USER_PARAMETER        !! TRANSFERS USER PARAMETERS TO MODULE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

USE WAM_FILE_MODULE, ONLY: IU05, FILE05, IU06

USE WAM_PRINT_USER_MODULE, ONLY: START_DATE, END_DATE,                         &
&                                OUTPUT_TIMESTEP, OUTPUT_TIMESTEP_UNIT,        &
&                                MOUTT, COUTT, NOUT_S, CFLAG_S,                &
&                                INPUT_FILE_NAME, INPUT_FILE_DATE,             &
&                                INPUT_FILE_TIMESTEP, INPUT_FILE_TIMESTEP_UNIT,&
&                                MOUTP, OUTLAT, OUTLONG, NAME

! ---------------------------------------------------------------------------- !

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER                         :: I, IOS, LEN, NOUTP, NOUTT
CHARACTER (LEN=80)              :: LINE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. OPEN USER INPUT FILE.                                                 !
!        ---------------------                                                 !

CALL CLEAR_PRINT_USER_MODULE

LEN = LEN_TRIM(FILE05)
IOS = 0
OPEN (UNIT=IU05, FILE=FILE05(1:LEN), FORM='FORMATTED', STATUS='OLD', IOSTAT=IOS)

IF (IOS.EQ.0) THEN
   CALL READ_PRINT_NAMELIST (1, IOS)
   IF (IOS.EQ.0) THEN
      CALL SET_PRINT_USER_PARAMETER
      CLOSE (UNIT=IU05, STATUS="KEEP")
      RETURN
   ELSE
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +     INFORMATION FROM SUB. READ_SPECTRA_USER      +'
      WRITE (IU06,*) ' +     =======================================      +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' + THE SPECTRA USER FILE COULD BE OPENED.           +'
      WRITE (IU06,*) ' + BUT DOES NOT CONTAIN THE PRINT NAMELIST.         +'
      WRITE (IU06,*) ' + GRID USER FILE NAME IS FILE05 = ', TRIM(FILE05)
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
   WRITE (IU06,*) ' +     INFORMATION FROM SUB. READ_SPECTRA_USER      +'
   WRITE (IU06,*) ' +     =======================================      +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + THE PRINT USER FILE COULD NOT BE OPENED.         +'
   WRITE (IU06,*) ' + SPECTRA USER FILE NAME IS FILE05 = ', TRIM(FILE05)
   WRITE (IU06,*) ' +         UNIT IS         IU05 = ', IU05
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' *  PROGRAM TRIES FOR NAMELIST IN STANDARD INPUT    +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'

   CALL READ_PRINT_NAMELIST (0, IOS)
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *     FATAL ERROR IN SUB. READ_SPECTRA_NAMELIST    *'
      WRITE (IU06,*) ' *     =========================================    *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * READ ERROR ON NAMLIST FROM STANDARD INPUT        *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   END IF
   CALL SET_PRINT_USER_PARAMETER
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. PRINT PERIOD AND TIMESTEP.                                            !
!        --------------------------                                            !

CALL F_NEW_DATA
START_DATE = LINE(2:15)
END_DATE   = LINE(18:31)

IF (LINE(32:43).NE.' ') THEN
   READ(LINE(34:40),'(I7 )',IOSTAT=IOS) OUTPUT_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('OUTPUT_TIMESTEP')
   IF (LINE(43:43) .EQ.'S'.OR.LINE(43:43) .EQ.'s') OUTPUT_TIMESTEP_UNIT = 'S'
   IF (LINE(43:43) .EQ.'M'.OR.LINE(43:43) .EQ.'m') OUTPUT_TIMESTEP_UNIT = 'M'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. MODEL OUTPUT AT FIXED TIMES.                                          !
!        ----------------------------                                          !

I = 0
DO
   CALL F_NEW_DATA
   IF (LINE(2:4).EQ.'END') EXIT
   I = I+1
   IF (I.LE.MOUTT) COUTT(I) = LINE( 2:15)
   IF (LINE(18:20).EQ.' ') CYCLE
   I = I+1
   IF (I.LE.MOUTT) COUTT(I) = LINE(18:31)
   IF (LINE(34:36).EQ.' ') CYCLE
   I = I+1
   IF (I.LE.MOUTT) COUTT(I) = LINE(34:47)
   IF (LINE(50:52).EQ.' ') CYCLE
   I = I+1
   IF (I.LE.MOUTT) COUTT(I) = LINE(50:63)
END DO

IF (I.GT.MOUTT) THEN
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                               +'
   WRITE(IU06,*) '+    WARNING ERROR IN SUB. READ_SPECTRA_USER    +'
   WRITE(IU06,*) '+    =======================================    +'
   WRITE(IU06,*) '+                                               +'
   WRITE(IU06,*) '+ NUMBER OF OUTPUT TIMES IN INPUT EXCEEDS       +'
   WRITE(IU06,*) '+ DIMENSION MOUTT                = ', MOUTT
   WRITE(IU06,*) '+ NUMBER OF TIMES INPUT IS     I = ', I
   WRITE(IU06,*) '+ CHANGE PARAMETER MOUTT IN PRINT_USER_MODULE   +'
   WRITE(IU06,*) '+                                               +'
   WRITE(IU06,*) '+   PROGRAM WILL IGNORE THE LAST OUTPUT TIMES   +'
   WRITE(IU06,*) '+                                               +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++'
   I = MOUTT
END IF
NOUTT = I

! ---------------------------------------------------------------------------- !
!                                                                              !
!    4. MODEL OUTPUT SITES.                                                    !
!       -------------------                                                    !

I = 0
DO
   CALL F_NEW_DATA                           !! SITES FOR SPECTRA.
   IF (LINE(2:4).EQ.'END') EXIT
   I = I + 1
   IF (I.LE.MOUTP) THEN
      OUTLONG(I) = LINE( 2:14)
      OUTLAT(I)  = LINE(16:28)
      NAME(I) = LINE(30:49)
   END IF
END DO

IF (I.GT.MOUTP) THEN
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                             +'
   WRITE(IU06,*) '+    WARNING ERROR IN SUB. READ_SPECTRA_USER  +'
   WRITE(IU06,*) '+    =======================================  +'
   WRITE(IU06,*) '+                                             +'
   WRITE(IU06,*) '+ NUMBER OF OUTPUT SITES IN INPUT EXCEEDS     +'
   WRITE(IU06,*) '+ DIMENSION MOUTP                = ', MOUTP
   WRITE(IU06,*) '+ NUMBER OF SITES INPUT IS I     = ', I
   WRITE(IU06,*) '+ CHANGE PARAMETER MOUTP IN PRINT_USER_MODULE +'
   WRITE(IU06,*) '+                                             +'
   WRITE(IU06,*) '+ PROGRAM WILL IGNORE THE LAST OUTPUT SITES   +'
   WRITE(IU06,*) '+                                             +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
   I = MOUTP
END IF
NOUTP = I

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. MODEL OUTPUT SELECTION.                                               !
!        -----------------------                                               !

DO I= 1,NOUT_S,2
   CALL F_NEW_DATA
   CFLAG_S(  I) = .NOT. (LINE( 2: 2).EQ.'N' .OR. LINE( 2: 2).EQ.'n')
   CFLAG_S(I+1) = .NOT. (LINE(38:38).EQ.'N' .OR. LINE(38:38).EQ.'n')
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. FIRST INPUT DATE AND INCREMENT.                                       !
!        -------------------------------                                       !

CALL F_NEW_DATA
INPUT_FILE_DATE = LINE( 2:15)

IF (LINE(16:27).NE.' ') THEN
   READ(LINE(18:24),'(I7 )',IOSTAT=IOS) INPUT_FILE_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('INPUT_FILE_TIMESTEP')
   IF (LINE(27:27).EQ.'S'.OR.LINE(27:27).EQ.'s') INPUT_FILE_TIMESTEP_UNIT = 'S'
   IF (LINE(27:27).EQ.'M'.OR.LINE(27:27).EQ.'m') INPUT_FILE_TIMESTEP_UNIT = 'M'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. INPUT FILENAME.                                                       !
!        ---------------                                                       !

CALL F_NEW_DATA
LEN = LEN_TRIM(LINE)
IF (LEN.GT.1) INPUT_FILE_NAME = LINE(2:LEN)  !! INPUT FILE (UNFORM. INPUT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     8. TRANSFER TO WAM_PRINT_MODULE.                                         !
!        -----------------------------                                         !

CALL SET_PRINT_USER_PARAMETER

RETURN

! ---------------------------------------------------------------------------- !

CONTAINS

! ---------------------------------------------------------------------------- !
!                                                                              !
!    18. INCLUDED FUNCTIONS.                                                   !
!        -------------------                                                   !
!                                                                              !
   SUBROUTINE F_NEW_DATA        !! FIND A NEW RECORD STARTING WITHOUT 'C'

   LINE(1:1) = 'C'
   DO WHILE (LINE(1:1).EQ.'C' .OR. LINE(1:1).EQ.'c')
      READ (IU05, '(A)',IOSTAT=IOS) LINE

      IF (IOS.EQ.0) CYCLE
      WRITE(IU06,*) ' ***********************************************'
      WRITE(IU06,*) ' *                                             *'
      WRITE(IU06,*) ' *   FATAL ERROR IN SUB. READ_SPECTRA_USER     *'
      WRITE(IU06,*) ' *   =====================================     *'
      WRITE(IU06,*) ' * READ ERROR ON INPUT FILE:                   *'
      WRITE(IU06,*) ' * LAST LINE READ IS     LINE = ', LINE
      WRITE(IU06,*) ' * ERROR NO. IS        IOSTAT = ', IOS
      WRITE(IU06,*) ' *                                             *'
      WRITE(IU06,*) ' *      PROGRAM ABORTS  PROGRAM ABORTS         *'
      WRITE(IU06,*) ' *                                             *'
      WRITE(IU06,*) ' ***********************************************'
      CALL ABORT1
   END DO

   END SUBROUTINE F_NEW_DATA

   SUBROUTINE ERROR_MESSAGE (MESSAGE)  !! READ ERROR MESSAGE
   CHARACTER (LEN=*) :: MESSAGE
      WRITE(IU06,*) ' *********************************************'
      WRITE(IU06,*) ' *                                           *'
      WRITE(IU06,*) ' *   FATAL ERROR IN SUB. READ_SPECTRA_USER   *'
      WRITE(IU06,*) ' *   =====================================   *'
      WRITE(IU06,*) ' *  READ ERROR ON CHARACTER STRING           *'
      WRITE(IU06,*) ' *                     IOSTAT = ', IOS
      WRITE(IU06,*) ' *  CHARACTER STRING IS  LINE = ', LINE
      WRITE(IU06,*) ' *  PROGRAM WANTS TO READ: ', MESSAGE
      WRITE(IU06,*) ' *                                           *'
      WRITE(IU06,*) ' *    PROGRAM ABORTS  PROGRAM ABORTS         *'
      WRITE(IU06,*) ' *                                           *'
      WRITE(IU06,*) '**********************************************'
      CALL ABORT1
   END SUBROUTINE ERROR_MESSAGE

END SUBROUTINE READ_SPECTRA_USER

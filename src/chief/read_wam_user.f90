SUBROUTINE READ_WAM_USER

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READ_WAM_USER - ROUTINE TO READ THE WAM USER INPUT.                        !
!                                                                              !
!     H. GUNTHER        GKSS/ECMWF     NOVEMBER 1989                           !
!     H. GUNTHER        GKSS           NOVEMBER 2001                           !
!     ROOP LALBEHARRY   MSC/ARMN       DECEMBER 2003                           !
!                       NAMELIST user input                                    !
!     ERIK MYKLEBUST                   NOVEMBER 2004                           !
!     H. GUNTHER        GKSS           JANUARY 2010   CYCLE 4.5.3              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       READS USER INPUT CONCERNING PERIOD OF INTEREST,TIMESTEPS,              !
!       MODEL OPTIONS, AND FILE NAMES ETC                                      !
!       IS DONE TOO.                                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        USER INFORMATION IS BEING READ WITH THE PRESUMPTIONS THAT:            !
!         1. EVERY LINE STARTING WITH 'C' IS A COMMENT LINE                    !
!         2. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'              !
!            (RIGHT-JUSTIFIED)                                                 !
!         THE USER INPUT SPECIFICATIONS ARE TRANSFERED TO DIFFERENT MODULES    !
!         BY CALLS OF SET_XXXX SUBROUTINES.                                    !
!                                                                              !
!         READING FROM NAMELIST IS TRIED, IF FAILED READ FROM Wam_User         !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!     -----------                                                              !

USE WAM_GENERAL_MODULE,   ONLY:    &
&       ABORT1                       !! TERMINATE PROCESSING.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY: IU05, FILE05, IU06

USE WAM_USER_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

IMPLICIT NONE

character (len=128) :: line
INTEGER             :: IOS, I, L

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. OPEN USER INPUT FILE.                                                 !
!        ---------------------                                                 !

CALL CLEAR_WAM_USER_MODULE

L = LEN_TRIM(FILE05)
IOS = 0
OPEN (UNIT=IU05, FILE=FILE05(1:L), FORM='FORMATTED', STATUS='OLD', IOSTAT=IOS)

IF (IOS.EQ.0) THEN
   CALL READ_WAM_NAMELIST (1, IOS)
   IF (IOS.EQ.0) THEN
      CALL SET_WAM_USER_PARAMETER
      CLOSE (UNIT=IU05, STATUS="KEEP")
      RETURN
   ELSE
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +     INFORMATION FROM SUB. READ_WAM_USER          +'
      WRITE (IU06,*) ' +     ===================================          +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' + THE WAM USER FILE COULD BE OPENED.               +'
      WRITE (IU06,*) ' + BUT DOES NOT CONTAIN THE WAM NAMELIST.           +'
      WRITE (IU06,*) ' + WAM USER FILE NAME IS FILE05 = ', TRIM(FILE05)  
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
   WRITE (IU06,*) ' +     INFORMATION FROM SUB. READ_WAM_USER          +'
   WRITE (IU06,*) ' +     ===================================          +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + THE WAM USER FILE COULD NOT BE OPENED.           +'
   WRITE (IU06,*) ' + WAM USER FILE NAME IS FILE05 = ', TRIM(FILE05)  
   WRITE (IU06,*) ' +         UNIT IS         IU05 = ', IU05
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' *  PROGRAM TRIES FOR NAMELIST IN STANDARD INPUT    +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'

   CALL READ_WAM_NAMELIST (1, IOS)
   IF (IOS.NE.0) THEN
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *     FATAL ERROR IN SUB. READ_WAM_NAMELIST        *'
      WRITE (IU06,*) ' *     =====================================        *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * READ ERROR ON NAMLIST FROM STANDARD INPUT        *'  
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   END IF
   CALL SET_WAM_USER_PARAMETER
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. MODEL PERIOD.                                                         !
!        -------------                                                         !

CALL F_NEW_DATA
START_DATE = LINE(2:15)
END_DATE   = LINE(18:31)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. MODEL START OPTIONS.                                                  !
!        --------------------                                                  !

CALL F_NEW_DATA
COLDSTART = .NOT.(SCAN(LINE(2:11),'F').GT.0 .OR. SCAN(LINE(2:11),'f').GT.0 )

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COLDSTART PARAMETERS.                                                 !
!        ---------------------                                                 !

CALL F_NEW_DATA
IF (LINE(2:9).NE.' ') THEN
   READ (LINE(2:9),'(I8)',IOSTAT=IOS) IOPTI
   IF (IOS.NE.0) CALL ERROR_MESSAGE('IOPTI')
END IF

CALL F_NEW_DATA
IF (LINE( 2:11).NE.' ') THEN
   READ (LINE( 2:11),'(F10.5)',IOSTAT=IOS) ALPHA
   IF (IOS.NE.0) CALL ERROR_MESSAGE('ALPHA')
END IF
IF (LINE(14:23).NE.' ') THEN
   READ (LINE(14:23),'(F10.5)',IOSTAT=IOS) FM
   IF (IOS.NE.0) CALL ERROR_MESSAGE('FM')
END IF
IF (LINE(26:35).NE.' ') THEN
   READ (LINE(26:35),'(F10.5)',IOSTAT=IOS) GAMMA
   IF (IOS.NE.0) CALL ERROR_MESSAGE('GAMMA')
END IF
IF (LINE(38:47).NE.' ') THEN
   READ (LINE(38:47),'(F10.5)',IOSTAT=IOS) SIGMA_A
   IF (IOS.NE.0) CALL ERROR_MESSAGE('SIGMA_A')
END IF
IF (LINE(50:59).NE.' ') THEN
   READ (LINE(50:59),'(F10.5)',IOSTAT=IOS) SIGMA_B
   IF (IOS.NE.0) CALL ERROR_MESSAGE('SIGMA_B')
END IF
IF (LINE(62:71).NE.' ') THEN
   READ (LINE(62:71),'(F10.5)',IOSTAT=IOS) THETAQ
   IF (IOS.NE.0) CALL ERROR_MESSAGE('THETAQ')
END IF

CALL F_NEW_DATA
IF (LINE(2:11).NE.' ') THEN
   READ (LINE(2:11),'(F10.1)',IOSTAT=IOS) FETCH
   IF (IOS.NE.0) CALL ERROR_MESSAGE('FETCH')
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. MODEL OPTIONS.                                                        !
!        --------------                                                        !

CALL F_NEW_DATA
SPHERICAL_RUN     = .NOT.(SCAN(LINE( 2: 8),'F').GT.0 .OR. SCAN(LINE( 2: 8),'f').GT.0) 
SHALLOW_RUN       = .NOT.(SCAN(LINE(11:17),'F').GT.0 .OR. SCAN(LINE(11:17),'f').GT.0)
REFRACTION_D_RUN  = (SCAN(LINE(20:26),'T').GT.0 .OR. SCAN(LINE(20:26),'t').GT.0)
REFRACTION_C_RUN  = (SCAN(LINE(29:35),'T').GT.0 .OR. SCAN(LINE(29:35),'t').GT.0)
WAVE_BREAKING_RUN = (SCAN(LINE(38:44),'T').GT.0 .OR. SCAN(LINE(38:44),'t').GT.0)
PHILLIPS_RUN      = (SCAN(LINE(47:53),'T').GT.0 .OR. SCAN(LINE(47:53),'t').GT.0)
IF (LINE(56:62).NE.' ') THEN
   READ(LINE(56:62),'(I7)',IOSTAT=IOS) ITEST
   IF (IOS.NE.0) CALL ERROR_MESSAGE('ITEST')
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. PROPAGATION AND SOURCE TIMESTEP.                                      !
!        --------------------------------                                      !

CALL F_NEW_DATA
READ(LINE( 2: 8),'(I7)',IOSTAT=IOS) PROPAGATION_TIMESTEP
IF (IOS.NE.0) CALL ERROR_MESSAGE('PROPAGATION_TIMESTEP')
IF (LINE(10:10) .EQ.'H'.OR.LINE(10:10) .EQ.'h') PROPAGATION_TIMESTEP_UNIT = 'H'
IF (LINE(10:10) .EQ.'M'.OR.LINE(10:10) .EQ.'m') PROPAGATION_TIMESTEP_UNIT = 'M'

IF (LINE(13:19).NE.' ') THEN
   READ(LINE(13:19),'(I7)',IOSTAT=IOS) SOURCE_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('SOURCE_TIMESTEP')
   IF (LINE(21:21) .EQ.'H'.OR.LINE(21:21) .EQ.'h') SOURCE_TIMESTEP_UNIT = 'H'
   IF (LINE(21:21) .EQ.'M'.OR.LINE(21:21) .EQ.'m') SOURCE_TIMESTEP_UNIT = 'M'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. RESTART.                                                              !
!        --------                                                              !

CALL F_NEW_DATA
IF (LINE( 2: 8).NE.' ') THEN
   READ(LINE( 2: 8),'(I7)',IOSTAT=IOS) RESTART_SAVE_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('RESTART_SAVE_TIMESTEP')
   IF (LINE(10:10) .EQ.'H'.OR.LINE(10:10) .EQ.'h') RESTART_SAVE_TIMESTEP_UNIT = 'H'
   IF (LINE(10:10) .EQ.'M'.OR.LINE(10:10) .EQ.'m') RESTART_SAVE_TIMESTEP_UNIT = 'M'
END IF
! ---------------------------------------------------------------------------- !
!                                                                              !
!     8. COARSE OR FINE GRID.                                                  !
!        --------------------                                                  !

CALL F_NEW_DATA
COARSE_GRID_RUN = (SCAN(LINE( 2: 8),'T').GT.0 .OR. SCAN(LINE( 2: 8),'t').GT.0)
FINE_GRID_RUN   = (SCAN(LINE(11:17),'T').GT.0 .OR. SCAN(LINE(11:17),'t').GT.0)

CALL F_NEW_DATA
IF (COARSE_GRID_RUN) THEN
   IF (LINE( 2: 8).NE.' ') THEN
      READ(LINE( 2: 8),'(I7)',IOSTAT=IOS) COARSE_OUTPUT_TIMESTEP
      IF (IOS.NE.0) CALL ERROR_MESSAGE('COARSE_OUTPUT_TIMESTEP')
      IF (LINE(10:10) .EQ.'H'.OR.LINE(10:10) .EQ.'h') COARSE_OUTPUT_TIMESTEP_UNIT = 'H'
      IF (LINE(10:10) .EQ.'M'.OR.LINE(10:10) .EQ.'m') COARSE_OUTPUT_TIMESTEP_UNIT = 'M'
   END IF

   IF (LINE(13:19).NE.' ') THEN
      READ(LINE(13:19),'(I7)',IOSTAT=IOS) COARSE_FILE_SAVE_TIMESTEP
      IF (IOS.NE.0) CALL ERROR_MESSAGE('COARSE_FILE_SAVE_TIMESTEP')
      IF (LINE(21:21) .EQ.'H'.OR.LINE(21:21) .EQ.'h') COARSE_FILE_SAVE_TIMESTEP_UNIT = 'H'
      IF (LINE(21:21) .EQ.'M'.OR.LINE(21:21) .EQ.'m') COARSE_FILE_SAVE_TIMESTEP_UNIT = 'M'
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     9. MODEL INPUT TIMESTEPS.                                                !
!        ----------------------                                                !

CALL F_NEW_DATA
IF (LINE( 2: 8).NE.' ') THEN
   READ(LINE( 2: 8),'(I7)',IOSTAT=IOS) WIND_OUTPUT_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('WIND_OUTPUT_TIMESTEP')
   IF (LINE(10:10) .EQ.'H'.OR.LINE(10:10) .EQ.'h') WIND_OUTPUT_TIMESTEP_UNIT = 'H'
   IF (LINE(10:10) .EQ.'M'.OR.LINE(10:10) .EQ.'m') WIND_OUTPUT_TIMESTEP_UNIT = 'M'
END IF

READ(LINE(13:19),'(I7)',IOSTAT=IOS) WIND_INPUT_TIMESTEP
IF (IOS.NE.0) CALL ERROR_MESSAGE('WIND_INPUT_TIMESTEP')
IF (LINE(21:21) .EQ.'H'.OR.LINE(21:21) .EQ.'h') WIND_INPUT_TIMESTEP_UNIT = 'H'
IF (LINE(21:21) .EQ.'M'.OR.LINE(21:21) .EQ.'m') WIND_INPUT_TIMESTEP_UNIT = 'M'

IF (LINE(24:30).NE.' ') THEN
   READ(LINE(24:30),'(I7)',IOSTAT=IOS) TOPO_OUTPUT_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('TOPO_OUTPUT_TIMESTEP')
   IF (LINE(32:32) .EQ.'H'.OR.LINE(32:32) .EQ.'h') TOPO_OUTPUT_TIMESTEP_UNIT = 'H'
   IF (LINE(32:32) .EQ.'M'.OR.LINE(32:32) .EQ.'m') TOPO_OUTPUT_TIMESTEP_UNIT = 'M'
END IF

IF (LINE(35:41).NE.' ') THEN
   READ(LINE(35:41),'(I7)',IOSTAT=IOS) TOPO_INPUT_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('TOPO_INPUT_TIMESTEPP')
   IF (LINE(43:43) .EQ.'H'.OR.LINE(43:43) .EQ.'h') TOPO_INPUT_TIMESTEP_UNIT = 'H'
   IF (LINE(43:43) .EQ.'M'.OR.LINE(43:43) .EQ.'m') TOPO_INPUT_TIMESTEP_UNIT = 'M'
END IF
IF (LINE(46:52).NE.' ') THEN
   READ(LINE(46:52),'(I7)',IOSTAT=IOS) CURRENT_OUTPUT_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('CURRENT_OUTPUT_TIMESTEP')
   IF (LINE(54:54) .EQ.'H'.OR.LINE(54:54) .EQ.'h') CURRENT_OUTPUT_TIMESTEP_UNIT = 'H'
   IF (LINE(54:54) .EQ.'M'.OR.LINE(54:54) .EQ.'m') CURRENT_OUTPUT_TIMESTEP_UNIT = 'M'
END IF

IF (LINE(57:63).NE.' ') THEN
   READ(LINE(57:63),'(I7)',IOSTAT=IOS) CURRENT_INPUT_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('CURRENT_INPUT_TIMESTEP')
   IF (LINE(65:65) .EQ.'H'.OR.LINE(65:65) .EQ.'h') CURRENT_INPUT_TIMESTEP_UNIT = 'H'
   IF (LINE(65:65) .EQ.'M'.OR.LINE(65:65) .EQ.'m') CURRENT_INPUT_TIMESTEP_UNIT = 'M'
END IF

IF (LINE(68:74).NE.' ') THEN
   READ(LINE(68:74),'(I7)',IOSTAT=IOS) ICE_INPUT_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('ICE_INPUT_TIMESTEP')
   IF (LINE(76:76) .EQ.'H'.OR.LINE(76:76) .EQ.'h') ICE_INPUT_TIMESTEP_UNIT = 'H'
   IF (LINE(76:76) .EQ.'M'.OR.LINE(76:76) .EQ.'m') ICE_INPUT_TIMESTEP_UNIT = 'M'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!    10. MODEL OUTPUT TIMESTEPS.                                               !
!        -----------------------                                               !

CALL F_NEW_DATA
IF (LINE( 2: 8).NE.' ') THEN
   READ(LINE( 2: 8),'(I7)',IOSTAT=IOS) PARAMETER_OUTPUT_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('PARAMETER_OUTPUT_TIMESTEP')
   IF (LINE(10:10) .EQ.'S'.OR.LINE(10:10) .EQ.'s') PARAMETER_OUTPUT_TIMESTEP_UNIT = 'S'
   IF (LINE(10:10) .EQ.'M'.OR.LINE(10:10) .EQ.'m') PARAMETER_OUTPUT_TIMESTEP_UNIT = 'M'
END IF

IF (LINE(13:19).NE.' ') THEN
   READ(LINE(13:19),'(I7)',IOSTAT=IOS) SPECTRA_OUTPUT_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('SPECTRA_OUTPUT_TIMESTEP')
   IF (LINE(21:21) .EQ.'S'.OR.LINE(21:21) .EQ.'s') SPECTRA_OUTPUT_TIMESTEP_UNIT = 'S'
   IF (LINE(21:21) .EQ.'M'.OR.LINE(21:21) .EQ.'m') SPECTRA_OUTPUT_TIMESTEP_UNIT = 'M'
END IF

IF (LINE(24:30).NE.' ') THEN
   READ(LINE(24:30),'(I7)',IOSTAT=IOS) OUTPUT_FILE_SAVE_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('OUTPUT_FILE_SAVE_TIMESTEP')
   IF (LINE(32:32) .EQ.'S'.OR.LINE(32:32) .EQ.'s') OUTPUT_FILE_SAVE_TIMESTEP_UNIT = 'S'
   IF (LINE(32:32) .EQ.'M'.OR.LINE(32:32) .EQ.'m') OUTPUT_FILE_SAVE_TIMESTEP_UNIT = 'M'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!    11. MODEL OUTPUT AT FIXED TIMES.                                          !
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
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                           +'
   WRITE(IU06,*) '+    WARNING ERROR IN SUB. READ_WAM_USER    +'
   WRITE(IU06,*) '+    ===================================    +'
   WRITE(IU06,*) '+                                           +'
   WRITE(IU06,*) '+ NUMBER OF OUTPUT TIMES IN INPUT EXCEEDS   +'
   WRITE(IU06,*) '+ DIMENSION MOUTT                = ', MOUTT
   WRITE(IU06,*) '+ NUMBER OF TIMES INPUT IS     I = ', I
   WRITE(IU06,*) '+ CHANGE PARAMETER MOUTT IN WAM_USER_MODULE +'
   WRITE(IU06,*) '+                                           +'
   WRITE(IU06,*) '+ PROGRAM WILL IGNORE THE LAST OUTPUT TIMES +'
   WRITE(IU06,*) '+                                           +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
   I = MOUTT
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!    12. MODEL OUTPUT SELECTION.                                               !
!        -----------------------                                               !

DO I=1,NOUT_P,2                            !! INTEGRATED PARAMETERS.
   CALL F_NEW_DATA
   PFLAG_P(  I) = LINE( 2: 2).EQ.'T' .OR. LINE( 2: 2).EQ.'t'
   FFLAG_P(  I) = .NOT. (LINE( 4: 4).EQ.'F' .OR. LINE( 4: 4).EQ.'f')
   PFLAG_P(I+1) = LINE(40:40).EQ.'T' .OR. LINE(40:40).EQ.'t'
   FFLAG_P(I+1) = .NOT. (LINE(42:42).EQ.'F' .OR. LINE(42:42).EQ.'f')
END DO

DO I=1,NOUT_S,2                              !! SPECTRA.
   CALL F_NEW_DATA
   PFLAG_S(  I) = LINE( 2: 2).EQ.'T' .OR. LINE( 2: 2).EQ.'t'
   FFLAG_S(  I) = .NOT. (LINE( 4: 4).EQ.'F' .OR. LINE( 4: 4).EQ.'f')
   PFLAG_S(I+1) = LINE(40:40).EQ.'T' .OR. LINE(40:40).EQ.'t'
   FFLAG_S(I+1) = .NOT. (LINE(42:42).EQ.'F' .OR. LINE(42:42).EQ.'f')
END DO
   
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
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                           +'
   WRITE(IU06,*) '+    WARNING ERROR IN SUB. READ_WAM_USER    +'
   WRITE(IU06,*) '+    ===================================    +'
   WRITE(IU06,*) '+                                           +'
   WRITE(IU06,*) '+ NUMBER OF OUTPUT SITES IN INPUT EXCEEDS   +'
   WRITE(IU06,*) '+ DIMENSION MOUTP                = ', MOUTP
   WRITE(IU06,*) '+ NUMBER OF SITES INPUT IS I     = ', I
   WRITE(IU06,*) '+ CHANGE PARAMETER MOUTP IN WAM_USER_MODULE +'
   WRITE(IU06,*) '+                                           +'
   WRITE(IU06,*) '+ PROGRAM WILL IGNORE THE LAST OUTPUT SITES +'
   WRITE(IU06,*) '+                                           +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
   I = MOUTP
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!    13. RADIATION STRESS.                                                     !
!        -----------------                                                     !

CALL F_NEW_DATA

IF (LINE( 2: 8).NE.' ') THEN
   READ(LINE( 2: 8),'(I7)',IOSTAT=IOS) RADIATION_OUTPUT_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('RADIATION_OUTPUT_TIMESTEP')
   IF (LINE(10:10) .EQ.'H'.OR.LINE(10:10) .EQ.'h') RADIATION_OUTPUT_TIMESTEP_UNIT = 'H'
   IF (LINE(10:10) .EQ.'M'.OR.LINE(10:10) .EQ.'m') RADIATION_OUTPUT_TIMESTEP_UNIT = 'M'
END IF

IF (LINE(13:19).NE.' ') THEN
   READ(LINE(13:19),'(I7)',IOSTAT=IOS) RADIATION_FILE_TIMESTEP
   IF (IOS.NE.0) CALL ERROR_MESSAGE('RADIATION_FILE_TIMESTEP')
   IF (LINE(21:21) .EQ.'H'.OR.LINE(21:21) .EQ.'h') RADIATION_FILE_TIMESTEP_UNIT = 'H'
   IF (LINE(21:21) .EQ.'M'.OR.LINE(21:21) .EQ.'m') RADIATION_FILE_TIMESTEP_UNIT = 'M'
END IF

DO I = 1,NOUT_R,2
   CALL F_NEW_DATA
   PFLAG_R(  I) = LINE( 2: 2).EQ.'T' .OR. LINE( 2: 2).EQ.'t'
   FFLAG_R(  I) = .NOT. (LINE( 4: 4).EQ.'F' .OR. LINE( 4: 4).EQ.'f')
   PFLAG_R(I+1) = LINE(40:40).EQ.'T' .OR. LINE(40:40).EQ.'t'
   FFLAG_R(I+1) = .NOT. (LINE(42:42).EQ.'F' .OR. LINE(42:42).EQ.'f')
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!    14. MODEL FILES.                                                          !
!        ------------                                                          !

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') WIND_INPUT_FILE_NAME = LINE(2:80)

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') FINE_INPUT_FILE_NAME = LINE(2:80)

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') ICE_INPUT_FILE_NAME = LINE(2:80)

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') PREPROC_OUTPUT_FILE_NAME = LINE(2:80)

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') TOPO_INPUT_FILE_NAME = LINE(2:80)

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') CURRENT_INPUT_FILE_NAME = LINE(2:80)

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') RESTART_FILE_NAME = LINE(2:80)

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') COARSE_OUTPUT_FILE_NAME = LINE(2:80)

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') PARAMETER_OUTPUT_FILE_NAME = LINE(2:80)

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') SPECTRA_OUTPUT_FILE_NAME = LINE(2:80)

CALL F_NEW_DATA
IF ( LINE(2:80).NE.' ') RADIATION_OUTPUT_FILE_NAME = LINE(2:80)
 
call f_new_data
if ( line(2: 4)/=' ') model_area = line(2:4)
 
call f_new_data
ready_file_flag = line(2:2)=='T'.or.line(2:2)=='t'

call f_new_data
if ( line(2:128)/=' ') ready_file_directory = line(2:128)

call f_new_data
ready_outfile_flag = line(2:2)=='T'.or.line(2:2)=='t'

call f_new_data
if ( line(2:128)/=' ') ready_outfile_directory = line(2:128)

call f_new_data
read (line( 2: 8),'(i7)',IOSTAT=IOS) hours_2d_spectra
IF (IOS.NE.0) CALL ERROR_MESSAGE('hours_2d_spectra')
read (line(11:17),'(i7)',IOSTAT=IOS) spectral_code
IF (IOS.NE.0) CALL ERROR_MESSAGE('spectral_code')

! ---------------------------------------------------------------------------- !
!
!    15. ASSIMILATION OPTIONS.
!
     
CALL F_NEW_DATA
IF (SCAN(LINE(2:8),'1').GT.0) THEN
   assimilation_flag = 1
   READ(LINE(11:17),'(F7.3)',IOSTAT=IOS) influence_radius
   IF (IOS.NE.0) CALL ERROR_MESSAGE('influence_radius')
   READ(LINE(20:26),'(F7.3)',IOSTAT=IOS) observation_scatter
   IF (IOS.NE.0) CALL ERROR_MESSAGE('observation_scatter')
   READ(LINE(29:35),'(F7.3)',IOSTAT=IOS) model_scatter
   IF (IOS.NE.0) CALL ERROR_MESSAGE('model_scatter')

   CALL F_NEW_DATA
   READ(LINE(2:15),'(A14)',IOSTAT=IOS) assimilation_start_date
   IF (IOS.NE.0) CALL ERROR_MESSAGE('assimilation_start_date')
   READ(LINE(18:31),'(A14)',IOSTAT=IOS) assimilation_end_date
   IF (IOS.NE.0) CALL ERROR_MESSAGE('assimilation_end_date')
   READ(LINE(34:40),'(I7)',IOSTAT=IOS) assimilation_time_step
   IF (IOS.NE.0) CALL ERROR_MESSAGE('assimilation_time_step')
   assimilation_time_step_unit = LINE(42:42)

   CALL F_NEW_DATA
   first_guess_output_flag = (SCAN(LINE(2:10),'T').GT.0                         &
&                             .OR. SCAN(LINE(2:10),'t').GT.0 )

   CALL F_NEW_DATA
   observation_filename    = line(2:80)    !! OBSERVATION INPUT FILE IDENTIFIER

   CALL F_NEW_DATA
   first_guess_ip_filename = line(2:80)    !! INTEGRATED DATA FILE (UNFORM. OUTPUT)

   CALL F_NEW_DATA
   first_guess_sp_filename = line(2:80)    !! SPECTRA DATA FILE (UNFORM. OUTPUT)
ELSE
   assimilation_flag = 0
END IF
    
! ---------------------------------------------------------------------------- !
!                                                                              !
!    16. TRANSFER USER PARAMETER INTO MODULES.                                 !
!        -------------------------------------                                 !

CALL SET_WAM_USER_PARAMETER

! ---------------------------------------------------------------------------- !
!
!    17. CLOSE INPUT FILE AND RETURN.                                          !
!        ----------------------------                                          !

CLOSE (UNIT=IU05, STATUS="KEEP")
write (iu06,*) ' '
WRITE (IU06,*) '   SUB. READ_WAM_USER SUCCESSFULLY COMPLETED. ' 

RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !
!    18. INCLUDED FUNCTIONS.                                                   !
!        -------------------                                                   !

CONTAINS

   SUBROUTINE F_NEW_DATA        !! FIND A NEW RECORD STARTING WITHOUT 'C'

   LINE(1:1) = 'C'
   DO WHILE (LINE(1:1).EQ.'C' .OR. LINE(1:1).EQ.'c')
      READ (IU05, '(A)',IOSTAT=IOS) LINE

      IF (IOS.EQ.0) CYCLE
      WRITE(IU06,*) ' ***********************************************'
      WRITE(IU06,*) ' *                                             *'
      WRITE(IU06,*) ' *     FATAL ERROR IN SUB. READ_WAM_USER       *'
      WRITE(IU06,*) ' *     ====================================    *'
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

   SUBROUTINE ERROR_MESSAGE (MESSAGE)
   CHARACTER (LEN=*) :: MESSAGE
      WRITE(IU06,*) ' *********************************************'
      WRITE(IU06,*) ' *                                           *'
      WRITE(IU06,*) ' *     FATAL ERROR IN SUB. READ_WAM_USER     *'
      WRITE(IU06,*) ' *     =================================     *'
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

END SUBROUTINE READ_WAM_USER

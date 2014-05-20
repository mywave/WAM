SUBROUTINE INITMDL

! ---------------------------------------------------------------------------- !
!                                                                              !
!    INITMDL - INITIALIZES THE WAM MODEL.                                      !
!                                                                              !
!     L. ZAMBRESKY   GKSS/ECMWF    JULY 1988                                   !
!                                                                              !
!     MODIFIED BY:   H. GUNTHER    NOVEMBER 1989                               !
!     H. GUNTHER      GKSS         OCTOBER 2000  FT90                          !
!     A. Behrens      MSC/ARMN     October 2003  MPI parallelization
!     E. Myklebust                 November 2004 MPI parallelization
!     H. GUNTHER      GKSS         JANUARY 2010  CYCLE_4.5.3                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       INITIALIZE THE WAM MODEL.                                              !
!                                                                              !
!     INTERFACE.                                                               !
!     ----------                                                               !
!                                                                              !
!          ---- INPUT/OUTPUT UNITS ---                                         !
!                                                                              !
!      THE PROGRAM OPENS AUTOMATICALLY THE FOLLOWING FILES, WHICH ARE          !
!      DEFINED IN "WAM_FILE_MODULE.f90":                                       !
!                                                                              !
!       UNIT = IU05 = 55  FILE05 = 'WAM_User' TO READ USER INPUT FILE.         !
!       UNIT = IU06 = 66  FILE06 = 'WAM_Prot' TO WRITE A PROTOCOL.             !
!                                                                              !
!      ALL OTHER FILE NAMES AND THE UNITS ARE PRE-DEFINED IN WAM_FILE_MODULE   !
!      TOO.                                                                    !
!      THE NAMES OF INPUT AND OUTPUT FILES CAN BE CHANGED IN THE WAM USER      !
!      INPUT FILE.                                                             !
!                                                                              !
!      THE PROGRAM USES OPEN TO ASSIGN FILES.                                  !
!      MODEL OUTPUT FILES ARE EXTENTED BY A DATE/TIME.                         !
!      FOR DETAILS OF THE FILE NAME CONVENTION OF THESE FILES SEE              !
!      IN MODULE WAM_GENERAL_MODULE SUB. OPEN_FILE.                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!          THIS ROUTINE INITIALISES THE WAVEMODEL:                             !
!            -  READS THE USER INPUT FILE,                                     !
!            -  READS THE DATA PRECOMPUTED BY PROG. PREPROC,                   !
!            -  READS THE RECOVERY FILES FOR A HOT START OR                    !
!            -  GENERATES AN INITIAL FIELD FOR A COLD START,                   !
!            -  DOES SOME GENERAL BOOKEEPING REGARDING                         !
!               DATES, INTEGRATION TIME STEPS AND OUTPUT TIME STEPS.           !
!            -  PREPARES THE FIRST ICE DATA.                                   !
!            -  PREPARES PROPAGATION.                                          !
!            -  PERFORMS A CFL CHECK.                                          !
!            -  PREPARES SOURCE FUNCTIONS.                                     !            
!            -  PREPARES OUTPUT.                                               !
!            -  PREPARES BOUNDARY PROCESSING.                                  !
!            -  OPENS THE FIRST RESULT FILES.                                  !
!                                                                              !
!     REFERENCE                                                                !
!     ---------                                                                !
!                                                                              !
!          A MORE DETAILED DISCUSSION MAY BE FOUND IN SUB WAMODEL.             !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!     -----------                                                              !

USE WAM_BOUNDARY_MODULE,     ONLY: &
&       PREPARE_BOUNDARY             !! PREPARES BOUNDARY PROCESSING.

USE WAM_GENERAL_MODULE,      ONLY: &
&       ABORT1,                    & !! TERMINATES PROCESSING.
&       INCDATE,                   & !! UPDATE DATE TIME GROUP.
&       READ_WAM_USER                !! READS USER INPUT.

USE WAM_INITIAL_MODULE,      ONLY: &
&       PREPARE_START,             & !! PREPARES START FIELDS.
&       READ_PREPROC_FILE            !! READS PREPROC OUTPUT FILE.

USE WAM_OUTPUT_SET_UP_MODULE,ONLY: &
&       PREPARE_OUTPUT,            & !! PREPARES OUTPUT.
&       SAVE_OUTPUT_FILES,         & !! CLOSES AND OPENS OUTPUT FILES.
&       UPDATE_OUTPUT_TIME           !! UPDATES OUTPUT TIMES.

USE WAM_OUTPUT_MODULE,       ONLY: &
&       MODEL_OUTPUT_CONTROL         !! CONTROLS MODEL OUTPUT.

USE WAM_PROPAGATION_MODULE,  ONLY: &
&       PREPARE_PROPAGATION          !! PREPARES PROPAGATION, DOES CFL CHECK.

USE WAM_RADIATION_MODULE,    ONLY: &
&       PREPARE_RADIATION            !! PREPARES RADIATION MODULE.

USE WAM_RESTART_MODULE,      ONLY: & 
&       PREPARE_RESTART_FILE         !! PREPARES RESTART_FILE.

USE WAM_SOURCE_MODULE,       ONLY: & 
&       PREPARE_SOURCE               !! PREPARES SOURCE FUNCTIONS.

use wam_mpi_comp_module,     only: &
&       mpi_decomp

use wam_assi_set_up_module,  only: &
&       prepare_assimilation         !! prepares the data assimilation

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE, ONLY: IU06, FILE06, ITEST, IU20, IU25, FILE20, FILE25
USE WAM_NEST_MODULE, ONLY: COARSE, FINE

USE WAM_OUTPUT_SET_UP_MODULE, ONLY: CDTINTT, CDTSPT, CDT_OUT, IDEL_OUT
USE WAM_TIMOPT_MODULE,        ONLY: CDTPRO

use wam_grid_module,        only: nsea, klat, klon, nx
use wam_mpi_module,         only: ninf
use wam_model_module,       only: fl3
use wam_assi_set_up_module, only: iassi

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INPUT OF USER PARAMETER.                                              !
!        ------------------------                                              !

CALL READ_WAM_USER
IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. INITMDL: READ_WAM_USER DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. READ PREPROC OUTPUT.                                                  !
!        ---------------------                                                 !

CALL READ_PREPROC_FILE
IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. INITMDL: READ_PREPROC_FILE DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.5 Decomposition of grid domain among processes.
!         ---------------------------------------------

call mpi_decomp 
if (itest>=2)  write (iu06,*) '   mpi_decomp: mpi_decomp'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. PREPARE START SPECTRA WINDS, TOPO AND CURRENTS.                       !
!        -----------------------------------------------                       !

WHERE (KLAT .EQ. 0) KLAT = NINF-1  !! UPDATE KLAT AND KLON FOR DECOMPOSED GRID.
WHERE (KLON .EQ. 0) KLON = NINF-1

CALL PREPARE_START
IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. INITMDL: PREPARE_START DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. PREPARE PROPAGATION AND PERFORM CFL CHECK.                            !
!         ------------------------------------------                           !

CALL PREPARE_PROPAGATION
IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. INITMDL: PREPARE_PROPAGATION DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. PREPARE SOURCE FUNCTIONS.                                             !
!        -------------------------                                             !

CALL PREPARE_SOURCE
IF (ITEST.GE.2) WRITE (IU06,*) '   SUB. INITMDL: PREPARE_SOURCE DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. PREPARE OUTPUT FOR INTEGRATED PARAMETER AND/OR SPECTRA.               !
!        -------------------------------------------------------               !

CALL PREPARE_OUTPUT   
IF (ITEST.GE.2)  WRITE(IU06,*) '   SUB. INITMDL: PREPARE_OUTPUT DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. PREPARE BOUNDARY VALUE HANDLING.                                      !
!        --------------------------------                                      !

IF (COARSE .OR. FINE) THEN
   CALL PREPARE_BOUNDARY
   IF (ITEST.GE.2)  WRITE(IU06,*) '   SUB. INITMDL: PREPARE_BOUNDARY DONE'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     8. INITIALIZE RESTART TIME.                                              !
!        ------------------------                                              !

CALL PREPARE_RESTART_FILE
IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '    SUB. INITMDL: PREPARE_RESTART_FILE DONE '
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     9. INITIALIZE RADIATION COMPUTATIONS AND OUTPUT.                         !
!        ---------------------------------------------                         !

CALL PREPARE_RADIATION (FL3)
IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '    SUB. INITMDL: PREPARE_RADIATION DONE '
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!    10. PREPARE ASSIMILATION.                                                 !
!

if (iassi==1) then
   call prepare_assimilation
   if (itest>=2) write (iu06,*) '    sub. initmdl: prepare_assimilation done'
endif

! ---------------------------------------------------------------------------- !
!                                                                              !
!    11. DO OUTPUT OF INITIAL FIELD.                                           !
!        ---------------------------                                           !

IF (CDTINTT.EQ.CDTPRO .OR. CDTSPT.EQ.CDTPRO) THEN
   CALL MODEL_OUTPUT_CONTROL (FL3, iu20, iu25)
   IF (CDT_OUT.EQ.CDTPRO) THEN
      CALL SAVE_OUTPUT_FILES (IU20, FILE20, IU25, FILE25)
      CALL INCDATE(CDT_OUT, IDEL_OUT)
   END IF
   CALL UPDATE_OUTPUT_TIME                          !! UPDATE OUTPUT TIMES.
   IF (ITEST.GE.2) WRITE(IU06,*) '    SUB. INITMDL: MODEL_OUTPUT_CONTROL DONE'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!    12. PRINT MODULE STATUS.                                                  !
!        --------------------                                                  !

WRITE (IU06,*) ' '
WRITE (IU06,*) '        MODULE STATUS AFTER INITIALISATION'
WRITE (IU06,*) '        ----------------------------------'
WRITE (IU06,*) ' '

CALL PRINT_WAM_STATUS

WRITE(IU06,*) '  '
WRITE(IU06,*) '         END OF WAM MODEL INITIALISATION'
WRITE(IU06,*) '         -------------------------------'
WRITE(IU06,*) '  '

END SUBROUTINE INITMDL

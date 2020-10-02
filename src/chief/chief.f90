PROGRAM CHIEF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     CHIEF - SUPERVISES WAVE MODEL EXECUTION.                                 !
!                                                                              !
!     LIANA ZAMBRESKY      GKSS/ECMWF  JUNE 1989                               !
!     H. GUNTHER           ECMWF       JUNE 1990  MODIFIED FOR CYCLE_4.        !
!     H. GUNTHER           GKSS        JANUARY 2010  CYCLE_4.5.3               !
!     A. Behrens           MSC/GKSS    January 2004  MPI parallelization       !
!     E. Myklebust                     November 2004 MPI parallelization       !
!     H. GUNTHER           GKSS        JANUARY 2010  CYCLE_4.5.3               !
!                                                                              !
!    PURPOSE.                                                                  !
!    --------                                                                  !
!                                                                              !
!       THIS PROGRAM SUPERVISES THE EXECUTION OF THE WAM MODEL.                !
!                                                                              !
!     INTERFACE.                                                               !
!     ----------                                                               !
!                                                                              !
!       IN ORDER FOR THE WAM MODEL TO EXECUTE,                                 !
!       IT NEEDS THE FOLLOWING INPUT FILES:                                    !
!                                                                              !
!       1. AN UNFORMATTED FILE CREATED BY PREPROC,                             !
!                                                                              !
!       2. A USER INPUT FILE,                                                  !
!                                                                              !
!       3. A WIND INPUT FILE,                                                  !
!                                                                              !
!       AND OPTIONAL:                                                          !
!                                                                              !
!       4. A RESTART FILE CREATED BY BY A PREVIOUS WAM,                        !
!                                                                              !
!       5. A BOUNDARY VALUE INPUT FILES CREATED BY A COARSE GRID WAM,          !
!                                                                              !
!       6. AN ICE INPUT FILE,                                                  !
!                                                                              !
!       7. A TOPOGRAPHY FILE WITH TIME DEPENDENT WATER DEPTH FIELDS,           !
!                                                                              !
!       8. A CURRENT FILE WITH TIME DEPENDENT CURRENT FIELDS.                  !
!                                                                              !
!     LIBRARIES.                                                               !
!     ----------                                                               !
!                                                                              !
!         NONE.                                                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       STARTING WITH THE INITIAL SEA STATE FIELD, THE WAM MODEL INTEGRATES    !
!       FORWARD IN TIME, DRIVEN BY WINDS.                                      !
!       SEA STATE AND RESULT FILES ARE SAVED IN REGULAR INTERVALLS.            !
!       THE SEA STATE FILE SERVE AS THE INITIAL CONDITION FOR A RESTART.       !
!                                                                              !
!       EACH CALL OF THE SUB WAVEMDL INTEGRATES FORWARD IN TIME BY             !
!       ONE WIND INPUT OR ONE PROPAGATION OR ONE SOURCE FUNCTION TIMESTEP,     !
!       WHAT EVER IS LONGER.                                                   !
!       IN THE FIRST CALL TO WAVEMDL AN INITIALIZATION IS  DONE IN ADDITION.   !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       EACH MODULE IS OF ITSELF THOROUGHLY DOCUMENTED.                        !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!     -----------                                                              !

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1,                  &  !! TERMINATES PROCESSING.
&       WAVEMDL                     !! SUPERVISES THE OVERALL FLOW THROUGH
                                    !! THE MAIN MODULES: INITMDL, PREWIND
                                    !! AND WAMODEL.

use wam_mpi_comp_module,     only: &
&       expand_string                !! prepares output of individual pu.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_TIMOPT_MODULE, ONLY: CDATEE, CDTPRO
USE WAM_FILE_MODULE,   ONLY: IU06, FILE06
use wam_mpi_module,    only: pelocal, petotal, nprevious, nnext,               &
&                            irank, extime, comtime

IMPLICIT NONE 
INCLUDE 'mpif.h'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

real (kind=kind(1d0)) :: time0, time
integer  :: ierr
character (len=80), dimension (1) :: logfilename

! ---------------------------------------------------------------------------- !
!
call MPI_INIT (ierr)
TIME0 = MPI_WTIME()

CDTPRO = ' '
CDATEE = '99999999999900'
    
! ---------------------------------------------------------------------------- !
!
!*    1. Initialize MPI
!        --------------

CALL MPI_COMM_RANK (MPI_COMM_WORLD, pelocal, ierr)
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, petotal, ierr)

irank = pelocal+1
nprevious = irank-1
if (irank==petotal) then
   nnext = 0
else
   nnext = irank+1
endif

iu06 = 66
IF (petotal.GT.1) THEN
   logfilename(1) ='logfile.%p'
   call expand_string (pelocal,petotal,0,0,logfilename,1)
   open (iu06, file=logfilename(1),status='unknown')
ELSE
   OPEN (UNIT=IU06, FILE=FILE06, FORM="FORMATTED", STATUS="UNKNOWN")
END IF     

write (iu06,*) ' +++ This is a parallel run using MPI '
write (iu06,*) ' +++ rank of local process      : ', pelocal
write (iu06,*) ' +++ total number of processors : ', petotal
write (iu06,*) 

! ---------------------------------------------------------------------------- !
!
!*    2. CALLS TO WAVEMDL UNTIL MODEL DATE REACHES END DATE. EACH CALL 
!*       INTEGRATES ONE WIND INPUT TIMESTEP, OR ONE PROPAGATION TIMESTEP, 
!*       OR ONE SOURCE FUNCTION TIMESTEP WHAT EVER IS LONGER.
!        ----------------------------------------------------------------

DO WHILE (CDTPRO<CDATEE)
   CALL WAVEMDL
END DO
     
TIME = MPI_WTIME()-TIME0
call MPI_finalize (ierr)
if (ierr==0) then
   write (iu06,*)
   write (iu06,*) ' +++ MPI successfully finalized ! '
else
   write (iu06,*) ' +++ error : finalize MPI ! '
endif
 
! ---------------------------------------------------------------------------- !
!
!*    3.  TERMINATE PROTOCOL.
!        --------------------

WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++'
WRITE (IU06,*) ' + TOTAL USER TIME IN SECONDS    +'
WRITE (IU06,*) ' + ', TIME
WRITE (IU06,*) ' +                               +'
WRITE (IU06,*) ' + COMMUNICATION TIME IN SECONDS +'
WRITE (IU06,*) ' + ', NINT(extime+comtime)
WRITE (IU06,*) ' +                               +'
WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++'

WRITE (*,*) ' Chief all done '

STOP
END PROGRAM CHIEF

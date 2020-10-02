MODULE WAM_RESTART_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL DATA WHICH ARE USED BY THE PRESET PROGRAM.        !
!   ALL PROCEDURES ARE INCLUDED TO COMPUTE THE WAM MODEL COLDSTART FILE,       !
!   TO COMPUTE SAVE AND CONNECT RESTART FILES TO THE WAM MODEL.                !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1,                  & !! TERMINATES PROCESSING.
&       INCDATE,                 & !! UPDATE DATE/TIME GROUP.
&       OPEN_FILE                  !! OPENS A FILE.

use wam_mpi_comp_module, only:   &
&       mpi_gather_fl,           &
&       mpi_gather_block,        &
&       mpi_exchng

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_MODEL_MODULE,   ONLY: FL3, U10, UDIR, TAUW, DEPTH, U, V

USE WAM_TIMOPT_MODULE,  ONLY: CDATEA, CDATEE, CDTPRO, CDTSOU, CDA, CDTA, CDCA, & 
&                             IDEL_WAM, cdtstop, L_DECOMP

USE WAM_FILE_MODULE,    ONLY: IU06, IU17, FILE17

use wam_fre_dir_module, only: ml, kl
use wam_grid_module,    only: nsea
use wam_mpi_module,     only: irank, nijs, nijl, ninf, nsup, i_out_restart,    &
&                             NGBTOPE, NTOPEMAX, NTOPELST, NTOPE, IJTOPE,      &
&                             NGBFROMPE, NFROMPEMAX, NFROMPELST, NFROMPE,      &
&                             NIJSTART,IJ2NEWIJ
use wam_special_module, only: ispec2d, ispecode
  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
include 'mpif.h'
PRIVATE

! ---------------------------------------------------------------------------- !
!                                                                              !
!    1.0 RESTART OPTION, TIMESTEP AND DATE.                                    !
!        ----------------------------------                                    !

INTEGER            :: IDEL_RES = -1    !! > O : TIMESTEP TO SAVE RESTART FILE.
                                       !! = 0 : FILE IS SAVED AT END OF RUN,
                                       !! < 0 : RESTART FILE IS NOT SAVED.
CHARACTER (LEN=14) :: CDT_RES  = ' '   !! NEXT DATE TO SAVE RESTART FILES.

PUBLIC CDT_RES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE CONNECT_RESTART           !! CONNECT RESTART FILE TO THE WAM MODEL.
   MODULE PROCEDURE CONNECT_RESTART
END INTERFACE
PUBLIC CONNECT_RESTART

INTERFACE SAVE_RESTART_FILE          !! SAVE RESTART FILE FOR THE WAM MODEL.
   MODULE PROCEDURE SAVE_RESTART_FILE
END INTERFACE
PUBLIC SAVE_RESTART_FILE

INTERFACE SET_RESTART_FILE_STEP
   MODULE PROCEDURE SET_RESTART_FILE_STEP
END INTERFACE
PUBLIC SET_RESTART_FILE_STEP

INTERFACE PREPARE_RESTART_FILE
   MODULE PROCEDURE PREPARE_RESTART_FILE
END INTERFACE
PUBLIC PREPARE_RESTART_FILE

INTERFACE PRINT_RESTART_STATUS
   MODULE PROCEDURE PRINT_RESTART_STATUS
END INTERFACE
PUBLIC PRINT_RESTART_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CONNECT_RESTART

! ---------------------------------------------------------------------------- !
!                                                                              !
!   CONNECT_RESTART - CONNECT RESTART FILE TO THE WAM MODEL.                   !
!                                                                              !
!       H. GUNTHER   GKSS      JULY 2001      F90                              !
!       A. Behrens   MSC/GKSS  December 2003  Message passing                  !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO CONNECT THE RESTART FILES TO THE WAM MODEL.                         !
!                                                                              !
!     METHOD.                                                                  !
!     --------                                                                 !
!                                                                              !
!       THE RESTART FILES FOR THE WAM MODEL ARE OPEND AND READ.                !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!                                                                              !
!     define local allocatable arrays
!     -------------------------------

real, allocatable, dimension (:,:,:) :: rfl

!     local variables
!     ----------------                                                         !

integer            :: ifail, NSEA_R, ios
character (len=14) :: zero = ' '
logical            :: unformatted = .true.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. OPEN RESTART FILE.                                                    !
!        ------------------                                                    !

if (unformatted) then
   CALL OPEN_FILE (IU06, IU17, FILE17, CDATEA, 'OLD', IFAIL)
   IF (IFAIL.NE.0) CALL ABORT1
endif

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. READ RESTART (binary or ascii code)                                   !
!        -----------------------------------                                   !

allocate (rfl(1:nsea,1:kl,1:ml))

ios = 0
if (unformatted) then
   READ (IU17,iostat=ios)  NSEA_R, CDTPRO, CDTSOU, CDA, CDTA, CDCA
   if (ios/=0) then
      close (iu17, status='keep')
      unformatted = .false.
      write (iu06,*) ' UNFORMATTED READING FROM FILE ', trim(FILE17),         &
&                    ' FAILED'
      write (iu06,*) ' TRYING FORMATTED READING'
   endif
endif
if (.not.unformatted) then
   call open_file (iu06, iu17, file17, cdatea, 'old', ifail, 'formatted')
   if (ifail/=0) call abort1
   read (iu17,*) nsea_r, cdtpro, cdtsou, cda, cdta, cdca
endif

IF (NSEA_R.NE.NSEA) THEN
   WRITE(IU06,*) ' *****************************************************'
   WRITE(IU06,*) ' *                                                   *'
   WRITE(IU06,*) ' *        FATAL ERROR IN SUB. CONNECT_RESTART        *'
   WRITE(IU06,*) ' *        ===================================        *'
   WRITE(IU06,*) ' *                                                   *'
   WRITE(IU06,*) ' * RESTART FIELDS ARE INCONSISTENT WITH MODEL GRID.  *'
   WRITE(IU06,*) ' * NO. OF SEA POINTS IN MODEL GRID IS       NSEA = ', NSEA
   WRITE(IU06,*) ' * NO. OF SEA POINTS IN RESTART FIELDS IS NSEA_R = ', NSEA_R
   WRITE(IU06,*) ' *                                                   *'
   WRITE(IU06,*) ' *           PROGRAM ABORTS     PROGRAM ABORTS       *'
   WRITE(IU06,*) ' *                                                   *'
   WRITE(IU06,*) ' *****************************************************'
   CALL ABORT1
END IF

if (unformatted) then                       !! binary code
   IF (L_DECOMP) THEN
     read (iu17) rfl(1:nsea,1,1)
     u10(:)     = rfl(nijs:nijl,1,1)
     read (iu17) rfl(1:nsea,1,1)
     udir(:)    = rfl(nijs:nijl,1,1)
     read (iu17) rfl(1:nsea,1,1)
     tauw(:)    = rfl(nijs:nijl,1,1)
     read (iu17) rfl(1:nsea,:,:)
     fl3(:,:,:) = rfl(nijs:nijl,:,:)
   ELSE
     read (iu17) rfl(IJ2NEWIJ(1:nsea),1,1)
     u10(:)     = rfl(nijs:nijl,1,1)
     read (iu17) rfl(IJ2NEWIJ(1:nsea),1,1)
     udir(:)    = rfl(nijs:nijl,1,1)
     read (iu17) rfl(IJ2NEWIJ(1:nsea),1,1)
     tauw(:)    = rfl(nijs:nijl,1,1)
     read (iu17) rfl(IJ2NEWIJ(1:nsea),:,:)
     fl3(:,:,:) = rfl(nijs:nijl,:,:)
   END IF
   IF (CDTA.NE.' ') then
      IF (L_DECOMP) THEN
         READ (IU17) rfl(1:nsea,1,1)
         depth(ninf:nsup) = rfl(ninf:nsup,1,1)
      ELSE
        READ (IU17) rfl(IJ2NEWIJ(1:nsea),1,1)
	    depth(nijs:nijl)=rfl(nijs:nijl,1,1)
	    call mpi_exchng(depth)
      END IF
   end if
   IF (CDCA.NE.' ') then
      IF (L_DECOMP) THEN
         READ (IU17) rfl(1:nsea,1,1), rfl(1:nsea,2,1)
         u(ninf:nsup) = rfl(ninf:nsup,1,1)
         v(ninf:nsup) = rfl(ninf:nsup,2,1)
      ELSE
         READ (IU17) rfl(IJ2NEWIJ(1:nsea),1,1), rfl(IJ2NEWIJ(1:nsea),2,1)
	     u(nijs:nijl)=rfl(nijs:nijl,1,1)
	     call mpi_exchng(u)
	     v(nijs:nijl)=rfl(nijs:nijl,2,1)
	     call mpi_exchng(v)
      END IF
   end if
else                           !! ascii code

   IF (L_DECOMP) THEN
      read (iu17,*) rfl(1:nsea,1,1)
      u10(:)     = rfl(nijs:nijl,1,1)
      read (iu17,*) rfl(1:nsea,1,1)
      udir(:)    = rfl(nijs:nijl,1,1)
      read (iu17,*) rfl(1:nsea,1,1)
      tauw(:)    = rfl(nijs:nijl,1,1)
      read (iu17,*) rfl(1:nsea,:,:)
      fl3(:,:,:) = rfl(nijs:nijl,:,:)
   ELSE
      read (iu17,*) rfl(IJ2NEWIJ(1:nsea),1,1)
      u10(:)     = rfl(nijs:nijl,1,1)
      read (iu17,*) rfl(IJ2NEWIJ(1:nsea),1,1)
      udir(:)    = rfl(nijs:nijl,1,1)
      read (iu17,*) rfl(IJ2NEWIJ(1:nsea),1,1)
      tauw(:)    = rfl(nijs:nijl,1,1)
      read (iu17,*) rfl(IJ2NEWIJ(1:nsea),:,:)
      fl3(:,:,:) = rfl(nijs:nijl,:,:)
   END IF
   if (cdta/='xxxxxxxxxxxxxx') then
      IF (L_DECOMP) THEN
         READ (IU17,*) rfl(1:nsea,1,1)
         depth(ninf:nsup) = rfl(ninf:nsup,1,1)
      else
         READ (IU17,*) rfl(IJ2NEWIJ(1:nsea),1,1)
	     depth(nijs:nijl)=rfl(nijs:nijl,1,1)
	     call mpi_exchng(depth)
      end if
   else
      cdta = ' '
   end if
   if (cdca/='xxxxxxxxxxxxxx') then
      IF (L_DECOMP) THEN
         READ (IU17,*) rfl(1:nsea,1,1), rfl(1:nsea,2,1)
         u(ninf:nsup) = rfl(ninf:nsup,1,1)
         v(ninf:nsup) = rfl(ninf:nsup,2,1)
     else
        READ (IU17,*) rfl(IJ2NEWIJ(1:nsea),1,1), rfl(IJ2NEWIJ(1:nsea),2,1)
	    u(nijs:nijl)=rfl(nijs:nijl,1,1)
	    call mpi_exchng(u)
	    v(nijs:nijl)=rfl(nijs:nijl,2,1)
	    call mpi_exchng(v)
     end if
   else
     cdca = ' '
   end if
endif
deallocate(rfl)

CLOSE (UNIT=IU17, STATUS="KEEP")

WRITE(IU06,*) ' '
WRITE(IU06,*) ' SUB. CONNECT_RESTART: RESTART FILE READ'
WRITE(IU06,*) ' '
WRITE(IU06,*) ' PROPAGATION DATE IS .............. CDTPRO  = ', CDTPRO
WRITE(IU06,*) ' SOURCE FUNCTION DATE IS .......... CDTSOU  = ', CDTSOU
WRITE(IU06,*) ' WIND FIELD DATE IS .................. CDA  = ', CDA
IF (CDTA.NE.' ') THEN
   WRITE(IU06,*) ' DEPTH FIELD DATE IS ................ CDTA  = ', CDTA
ELSE
   WRITE(IU06,*) ' DEPTH FIELD WAS NOT STORED IN RESTART FILE '
END IF
IF (CDCA.NE.' ') THEN
   WRITE(IU06,*) ' CURRENT FIELD DATE IS .............. CDCA  = ', CDCA
ELSE
   WRITE(IU06,*) ' CURRENT FIELD WAS NOT STORED IN RESTART FILE '
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. CHECK RESTART TIME.                                                   !
!        -------------------                                                   !

IF (CDTPRO.NE.CDATEA) THEN
   WRITE(IU06,*) ' ************************************************'
   WRITE(IU06,*) ' *                                              *'
   WRITE(IU06,*) ' *      FATAL ERROR IN SUB. CONNECT_RESTART     *'
   WRITE(IU06,*) ' *      ===================================     *'
   WRITE(IU06,*) ' *                                              *'
   WRITE(IU06,*) ' * START DATE AND RESTART FIELD DO NOT MATCH    *'
   WRITE(IU06,*) ' * START DATE OF RUN       IS CDATEA = ', CDATEA
   WRITE(IU06,*) ' * START DATE FROM RESTART IS CDTPRO = ', CDTPRO
   WRITE(IU06,*) ' *                                              *'
   WRITE(IU06,*) ' *      PROGRAM ABORTS     PROGRAM ABORTS       *'
   WRITE(IU06,*) ' *                                              *'
   WRITE(IU06,*) ' ************************************************'
   CALL ABORT1
END IF

END SUBROUTINE CONNECT_RESTART

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_RESTART_FILE

IF (IDEL_RES.GT.0) THEN

   IF (MOD(IDEL_RES,IDEL_WAM).NE.0) THEN
      WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(IU06,*) '+                                                       +'
      WRITE(IU06,*) '+      WARNING ERROR IN SUB. PREPARE_RESTART_FILE       +'
      WRITE(IU06,*) '+      ==========================================       +'
      WRITE(IU06,*) '+                                                       +'
      WRITE(IU06,*) '+ A RESTART FILE IS REQUESTED EVERY                     +'
      WRITE(IU06,*) '+                  IDEL_RES  :',IDEL_RES,' SECONDS'
      WRITE(IU06,*) '+                                                       +'
      WRITE(IU06,*) '+     IDEL_RES MUST BE AN INTEGER MULTIPLE OF           +'
      WRITE(IU06,*) '+             THE MAXIMUM OF THE                        +'
      WRITE(IU06,*) '+   WIND INPUT,  DEPTH INPUT,  CURRENT INPUT,           +'
      WRITE(IU06,*) '+   SOURCE FUNCTION, AND PROPAGATION TIMESTEP:          +'
      WRITE(IU06,*) '+                  IDEL_WAM  :',IDEL_WAM,' SECONDS'
      WRITE(IU06,*) '+                                                       +'
      WRITE(IU06,*) '+ RESTART TIMESTEP IS CHANGED TO NEAREST MULTIPLE       +'
      IDEL_RES = MAX(NINT(REAL(IDEL_RES)/REAL(IDEL_WAM)),1) * IDEL_WAM
      WRITE(IU06,*) '+   NEW RESTART TIMESTEP     : ',IDEL_RES, ' SECONDS'
      WRITE(IU06,*) '+                                                       +'
      WRITE(IU06,*) '+                  PROGRAM CONTINUES                    +'
      WRITE(IU06,*) '+                                                       +'
      WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   END IF

   CDT_RES = CDTPRO
   CALL INCDATE(CDT_RES, IDEL_RES)

ELSE IF (IDEL_RES.EQ.0) THEN
   CDT_RES = CDATEE
ELSE
   CDT_RES = ' '
END IF

if (ispec2d>0) then
   cdtstop = cdatea
   call incdate (cdtstop,ispec2d*3600)  !! last date for full spectral output
else
   cdtstop = cdatee
endif

END SUBROUTINE PREPARE_RESTART_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_RESTART_STATUS

WRITE(IU06,*) '  '
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '               RESTART MODULE STATUS:'
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '  '

IF (IDEL_RES.GT.0) THEN
   WRITE (IU06,*) ' RESTART FILE IS SAVED EVERY......: ', IDEL_RES,' SECONDS'
   WRITE (IU06,*) ' NEXT DATE TO SAVE FILE IS........: ', CDT_RES
   if (ispec2d/=0) then
      write (iu06,*) ' FULL SPECTRAL OUTPUT FOR.........  ', ispec2d, ' HOURS'
      write (iu06,*) ' LAST DATE FOR RESTART FILE.......: ', cdtstop
   endif
   if (ispecode==1) then
      write (iu06,*) ' SPECTRAL DATA IS WRITTEN IN ASCII CODE'
   else
      write (iu06,*) ' SPECTRAL DATA IS WRITTEN IN BINARY CODE'
   endif
ELSE IF (IDEL_RES.EQ.0) THEN
   WRITE(IU06,*) ' RESTART FILE IS SAVED AT END OF RUN.'
   if (ispecode==1) then
      write (iu06,*) ' SPECTRAL DATA IS WRITTEN IN ASCII CODE'
   else
      write (iu06,*) ' SPECTRAL DATA IS WRITTEN IN BINARY CODE'
   endif
ELSE
   WRITE(IU06,*) ' RESTART FILE IS NOT SAVED.'
END IF

END SUBROUTINE PRINT_RESTART_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SAVE_RESTART_FILE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SAVE_RESTART_FILE  - TO SAVE RESTART FILE FOR THE WAM MODEL.               !
!                                                                              !
!                                                                              !
!     H. GUENTHER   GKSS      FEBRUARY 2002   FT 90                            !
!     A. Behrens    MSC/GKSS  December 2003   Message passing                  !
!     E. Myklebust            November 2004   MPI parallelization              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        COPY ALL DATA FIELDS AND VARIABLES NESSECCARY FOR A RESTART OF THE    !
!        WAM MODEL TO FILE.                                                    !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!        A FILE IS OPENED AND THE FOLLOWING ACTUAL INFORMATION                 !
!        IS COPIED OUT OF THE MODEL RUN AT THE TIME OF THE RESTART:            !
!          - THE WIND INFORMATION TO GATHER WITH THE TIME COUNTERS,            !
!          - THE SPECTRA AT ALL GRID POINTS,                                   !
!        THE FILE NAME CONVENTION IS DEFINED IN SUB. GFILE.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!        NONE                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     define local allocatable arrays
!     -------------------------------

real, allocatable, dimension (:,:,:) :: rfl
real, allocatable, dimension (:)     :: ru10, rudir, rtauw, ru, rv, rdepth

!     local variables
!     ---------------

integer :: ifail = 0
integer :: ierr
     
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. OPEN FILE AND WRITE OUT.                                              !
!        ------------------------                                              !

if (irank==i_out_restart) then
   allocate (rfl(1:nsea,1:kl,1:ml))
   allocate (ru10(1:nsea), rudir(1:nsea), rtauw(1:nsea))

   call mpi_gather_fl (i_out_restart,121,fl3,rfl)
   CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   CALL mpi_gather_block(i_out_restart, u10, ru10)
   CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   CALL mpi_gather_block(i_out_restart, udir, rudir)
   CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   CALL mpi_gather_block(i_out_restart, tauw, rtauw)
   CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   IF (CDTA.NE.' ') then
      allocate (rdepth(1:nsea))
      CALL mpi_gather_block(i_out_restart, depth(nijs:nijl), rdepth)
      CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   end if
   IF (CDCA.NE.' ') then
      allocate (ru(1:nsea), rv(1:nsea))
      CALL mpi_gather_block(i_out_restart, u(nijs:nijl), ru)
      CALL mpi_barrier(MPI_COMM_WORLD,ierr)
      CALL mpi_gather_block(i_out_restart, v(nijs:nijl), rv)
      CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   end if
   REWIND IU17
   if (ispecode==1) then                                   !! asci code
      call open_file (iu06, iu17, file17, cdtpro, 'unknown', ifail, 'FORMATTED')
      if (ifail/=0) call abort1
      if (cdta==' ') cdta = 'xxxxxxxxxxxxxx'
      if (cdca==' ') cdca = 'xxxxxxxxxxxxxx'
      write (iu17,'(i8,5(2x,a14))') nsea, cdtpro, cdtsou, cda, cdta, cdca

      IF (L_DECOMP) THEN
         write (iu17,*) ru10
         write (iu17,*) rudir
         write (iu17,*) rtauw
         write (iu17,*) rfl
         if (cdta/='xxxxxxxxxxxxxx') write (iu17,*) rdepth
         if (cdca/='xxxxxxxxxxxxxx') write (iu17,*) ru, rv
      ELSE
         write (iu17,*) ru10(IJ2NEWIJ(1:NSEA))
         write (iu17,*) rudir(IJ2NEWIJ(1:NSEA))
         write (iu17,*) rtauw(IJ2NEWIJ(1:NSEA))
         write (iu17,*) rfl(IJ2NEWIJ(1:NSEA),:,:)
         if (cdta/='xxxxxxxxxxxxxx') write (iu17,*) rdepth(IJ2NEWIJ(1:NSEA))
         if (cdca/='xxxxxxxxxxxxxx') write (iu17,*) ru(IJ2NEWIJ(1:NSEA)), rv(IJ2NEWIJ(1:NSEA))
      ENDIF
      if (cdta=='xxxxxxxxxxxxxx') cdta = ' '
      if (cdca=='xxxxxxxxxxxxxx') cdca = ' '

   else                                                    !! binary code
      call open_file (iu06, iu17, file17, cdtpro, 'unknown', ifail)
      if (ifail/=0) call abort1

      IF (L_DECOMP) THEN
         WRITE (IU17) NSEA, CDTPRO, CDTSOU, CDA, CDTA, CDCA
         WRITE (IU17) ru10
         WRITE (IU17) rudir
         WRITE (IU17) rtauw
         WRITE (IU17) rfl
         IF (CDTA.NE.' ') WRITE (IU17) rdepth
         IF (CDCA.NE.' ') WRITE (IU17) rU, rV
      ELSE
         WRITE (IU17) NSEA, CDTPRO, CDTSOU, CDA, CDTA, CDCA
         WRITE (IU17) ru10(IJ2NEWIJ(1:NSEA))
         WRITE (IU17) rudir(IJ2NEWIJ(1:NSEA))
         WRITE (IU17) rtauw(IJ2NEWIJ(1:NSEA))
         WRITE (IU17) rfl(IJ2NEWIJ(1:NSEA),:,:)
         IF (CDTA.NE.' ') WRITE (IU17) rdepth(IJ2NEWIJ(1:NSEA))
         IF (CDCA.NE.' ') WRITE (IU17) rU(IJ2NEWIJ(1:NSEA)), rV(IJ2NEWIJ(1:NSEA))
      ENDIF
   endif
   CLOSE (UNIT=IU17, STATUS="KEEP")
   if (allocated(rdepth)) deallocate(rdepth)
   if (allocated(rfl)) deallocate(rfl)
   if (allocated(ru10)) deallocate(ru10)
   if (allocated(rudir)) deallocate(rudir)
   if (allocated(rtauw)) deallocate(rtauw)
   if (allocated(ru)) deallocate(ru)
   if (allocated(rv)) deallocate(rv)
else
   CALL mpi_gather_fl(i_out_restart,121,fl3)
   CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   CALL mpi_gather_block(i_out_restart, u10)
   CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   CALL mpi_gather_block(i_out_restart, udir)
   CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   CALL mpi_gather_block(i_out_restart, tauw)
   CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   IF (CDTA.NE.' ') then
      CALL mpi_gather_block(i_out_restart, depth(nijs:nijl))
      CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   end if
   IF (CDCA.NE.' ') then
      CALL mpi_gather_block(i_out_restart, u(nijs:nijl))
      CALL mpi_barrier(MPI_COMM_WORLD,ierr)
      CALL mpi_gather_block(i_out_restart, v(nijs:nijl))
      CALL mpi_barrier(MPI_COMM_WORLD,ierr)
   end if
endif

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. PRINT A MESSAGE.                                                      !
!        ----------------                                                      !

WRITE(IU06,*) ' '
WRITE(IU06,*) ' SUB. SAVE_RESTART: RESTART FILE SAVED'
WRITE(IU06,*) ' '
WRITE(IU06,*) ' PROPAGATION DATE IS .............. CDTPRO  = ', CDTPRO
WRITE(IU06,*) ' SOURCE FUNCTION DATE IS .......... CDTSOU  = ', CDTSOU
WRITE(IU06,*) ' WIND FIELD DATE IS .................. CDA  = ', CDA
IF (CDTA.NE.' ') THEN
   WRITE(IU06,*) ' DEPTH FIELD DATE IS ................ CDTA  = ', CDTA
END IF
IF (CDCA.NE.' ') THEN
   WRITE(IU06,*) ' CURRENT FIELD DATE IS .............. CDCA  = ', CDCA
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. NEXT DATE TO WRITE A RESTART FILE.                                    !
!        ----------------------------------                                    !

IF (CDT_RES.LT.CDATEE) THEN
   CALL INCDATE(CDT_RES,IDEL_RES)
   WRITE(IU06,*) ' NEXT RESTART FILE IS SAVED AT ... CDT_RES  = ', CDT_RES
END IF

END SUBROUTINE SAVE_RESTART_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_RESTART_FILE_STEP (STEP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, OPTIONAL, INTENT(IN) :: STEP  !! > O : TIMESTEP TO SAVE RESTART FILE.
                                       !! = 0 : FILE IS SAVED AT END OF RUN,
                                       !! < 0 : RESTART FILE IS NOT SAVED.

! ---------------------------------------------------------------------------- !

IF (PRESENT(STEP)) THEN
   IDEL_RES = STEP
ELSE
   IDEL_RES = -1
END IF

END SUBROUTINE SET_RESTART_FILE_STEP

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_RESTART_MODULE

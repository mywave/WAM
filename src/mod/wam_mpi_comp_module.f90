MODULE WAM_MPI_COMP_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   This module contains all required MPI subroutines                          !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1,                  &   !! TERMINATES PROCESSING.
&       PRINT_ARRAY,             &    !! PRINTS AN ARRAY
&       SORTI,                   &
&       SORTINI

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

use WAM_COORDINATE_MODULE

use wam_file_module,          only: iu06, itest
use wam_fre_dir_module,       only: ml, kl
use wam_grid_module,          only: nx, ny, nsea, NLON_RG, IXLG, KXLT,         &
&                                   AMOWEP, AMOEAP, AMOSOP, AMONOP, IPER,      &
&                                   XDELLO, ZDELLO, KLON, KLAT, DEPTH_B,       &
&                                   dellam, COSPH, WLAT, IFROMIJ, KFROMIJ,     &
&                                   OBSLAT, OBSLON
use wam_output_set_up_module, only: nout_P, noutp, ijar, npout
use wam_nest_module,          only: n_nest, max_nest, nbounc, ijarc
use wam_mpi_module,           only: petotal, irank, nstart, nend, nlen,        &
&                                   klentop, klenbot, mpmaxlength,             &
&                                   nnext, nprevious, ninf, nsup, nijs, nijl,  &
&                                   npoi,                                      &
&                                   ipfgtbl, i_out_par, i_out_spec,            &
&                                   i_out_b_spec, i_out_restart,               &
&                                   noutp_ga, ijar_ga, ngou_ga,                &
&                                   nbounc_ga, ijarc_ga, ngouc_ga,             &
&                                   extime, comtime,                           &
&                                   NGBTOPE, NTOPEMAX, NTOPELST, NTOPE, IJTOPE,&
&                                   NGBFROMPE, NFROMPEMAX, NFROMPELST, NFROMPE,&
&                                   NIJSTART, IJ2NEWIJ
USE WAM_TIMOPT_MODULE,        ONLY: L_DECOMP, L_OBSTRUCTION

implicit none
include 'mpif.h'

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

PRIVATE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE expand_string              !! EXPAND STRING BY ADDING THE PROCESS NUMBER,
   MODULE PROCEDURE expand_string    !! AND/OR THE TOTAL NUMBER OF PROCESSES,    
END INTERFACE                        !! AND/OR THE TIME STEP NUMBER,
PUBLIC expand_string                 !! AND/OR THE INDEX OF THE STRING

INTERFACE mpi_crtbl                  !! creates a table associating integrated
   MODULE PROCEDURE mpi_crtbl        !! output parameter with a specific process
END INTERFACE
PUBLIC mpi_crtbl

INTERFACE mpi_decomp                 !! developes an even decomposition of the
   MODULE PROCEDURE mpi_decomp       !! grid domain among the available processes
END INTERFACE
PUBLIC mpi_decomp

INTERFACE mpi_exchng                 !! exchanges messages between the process
   MODULE PROCEDURE mpi_exchng       !! pelocal and the previous and next one
   MODULE PROCEDURE mpi_exchng_V     !! pelocal and the previous and next one
END INTERFACE
PUBLIC mpi_exchng

INTERFACE mpi_gather_block           !! gather scalar block data field
   MODULE PROCEDURE mpi_gather_block !! onto a single process
END INTERFACE
PUBLIC mpi_gather_block

INTERFACE mpi_gather_bound           !! gathers spectrum onto a single process
   MODULE PROCEDURE mpi_gather_bound !! for ouput of boundary values
END INTERFACE
PUBLIC mpi_gather_bound

INTERFACE mpi_gather_fl              !! gather spectral field fl onto a 
   MODULE PROCEDURE mpi_gather_fl    !! single process
END INTERFACE
PUBLIC mpi_gather_fl

INTERFACE mpi_gather_grid            !! gather grid data field from the 
   MODULE PROCEDURE mpi_gather_grid  !! process isend onto the process irecv
END INTERFACE
PUBLIC mpi_gather_grid

INTERFACE mpi_gather_spp             !! gathers spectrum and scalar fields onto 
   MODULE PROCEDURE mpi_gather_spp   !! a single process for output 
END INTERFACE                        !! at selected grid points
PUBLIC mpi_gather_spp

interface mpi_gather_oifl
   module procedure mpi_gather_oifl_real     !! broadcast of scalar real or
   module procedure mpi_gather_oifl_logical  !! logical grid data contained
end interface                        !! in array field stored on process isend
public mpi_gather_oifl               !! onto all processes 

interface mpi_gather_cfl
   module procedure mpi_gather_cfl  !! broadcast of number contained
end interface                       !! in array field stored on process irank
public mpi_gather_cfl               !! onto all processes 


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

subroutine expand_string (                                             &
&          myproc,                                                     &    ! %p
&          nproc,                                                      &    ! %n
&          timestep,                                                   &    ! %t
&          max_timestep,                                               &
&          s,                                                          &    ! %s
&          n)
!
! ---------------------------------------------------------------------------- !
!
!     S. SAARINEM   ECMWF     MAY 1996      MESSAGE PASSING
!     A. Behrens    MSC/GKSS  October 2003  f90
!
!*    PURPOSE.
!     --------
!
!     EXPAND STRING BY ADDING AN EXTENSION WHICH IS FUNCTION OF
!     THE PROCESS NUMBER, AND/OR THE TOTAL NUMBER OF PROCESSES,
!     AND/OR THE TIME STEP NUMBER, AND/OR THE INDEX OF THE STRING
!
!**   INTERFACE.
!     ----------
!
!     CALL EXPAND_STRING (MYPROC,NPROC,TIMESTEP,MAX_TIMESTEP,S,N)
!
!     *MYPROC*       INTEGER  PROCESS NUMBER (PE)
!     *NPROC*        INTEGER  TOTAL NUMBER OF PE'S
!     *TIMESTEP*     INTEGER  TIME STEP NUMBER (OR ANY OTHER INDEX)
!     *MAX_TIMESTEP* INTEGER  MAXIMUM TIME STEP NUMBER (OR ANY OTHER INDEX)
!     *S*            CHARACTER OR CHARACTER ARRAY THAT WILL BE EXPANDED
!     *N*            INTEGER  DIMENSION OF S
!
!
!     METHOD.
!     -------
!     THE STRING S MUST CONTAIN THE FOLLOWING SUB-STRING THAT WILL THAN BE
!     SUBSTITUTED BY THE CORRESPONDING EXPANDED SUB STRING, THAT STARTS WITH
!     % AND p, n, t, or s
!
!     %p   -->  myproc in character mode, left justified
!     %n   -->  nproc  in character mode, left justified
!     %t   -->  timestep in character mode, left justified
!     %s   -->  the array index of S, left justified
!
!     e.g.   s='name_%n.%p'
!            call(myproc,nproc,0,0,s,1)
!
!     will produce, if nproc=11
!
!     on PE1,  s = name_11.1
!     on PE2,  s = name_11.2
!     on PE3,  s = name_11.3
!
!     ...
!
!     on PE10,  s = name_11.10
!     on PE11,  s = name_11.11
!
!
!     EXTERNALS.
!     ---------
!
!      none
!
!     REFERENCE.
!     ----------
!
!      none
!
! ---------------------------------------------------------------------------- !

integer, intent(in)  :: myproc, nproc
integer, intent(in)  :: timestep, max_timestep
integer, intent(in)  :: n
      
character (len=*), intent(inout) :: s(:)
    
!
! === END OF INTERFACE BLOCK ===

character (len=len(s))   :: t
character (len=2*len(s)) :: tt
character (len=6), dimension (4) :: fmt

integer :: i, j, jj, loc_p, len_t
integer, dimension (4) :: ndigs, num
    
if (n < 1) return              !! Setup output formats
num(1) = myproc
num(2) = max(nproc,myproc)
num(3) = n
num(4) = max(max_timestep,timestep)

do j=1,4                       !! Count number of digits in each integer
   ndigs(j) = 1
   if (num(j) /= 0) then
      ndigs(j) = 1 + log10(dble(abs(num(j))))
      if (num(j) < 0) ndigs(j) = ndigs(j) + 1 ! Room for minus sign
   endif
   ndigs(j) = min(ndigs(j),9)   ! Max 9 digits supported; i.e. '999 999 99
   write(fmt(j),'("(i",i1,")")') ndigs(j)
enddo
!
!*    Expand fields '%s', '%p', '%n' and '%t' with their values
!
!*    A special treatment with the sequence numbering

if (n>1) then
   loc_p = index(s(1),'%s')
   if (loc_p > 0) then
      s(2:) = s(1)
   endif
endif
do i=1,n
   t = adjustl(s(i))//' '
   loc_p = index(t,'%')
   if (loc_p > 0) then
      len_t = len_trim(t)
      j = loc_p
      tt(:j-1) = t(:j-1)
      tt(j:) = ' '
      jj = j-1
      do while (j <= len_t)
         if (t(j:j) == '%') then
            j = j + 1
            if (j <= len_t) then
               select case ( t(j:j) )
                 case ( 'p' )   ! myproc
                    write(tt(jj+1:jj+ndigs(1)),fmt(1)) myproc
                    jj = jj + ndigs(1)
                 case ( 'n' )   ! nproc
                    write(tt(jj+1:jj+ndigs(2)),fmt(2)) nproc
                    jj = jj + ndigs(2)
                 case ( 's' )   ! sequence number i=[1..n]
                    write(tt(jj+1:jj+ndigs(3)),fmt(3)) i
                    jj = jj + ndigs(3)
                 case ( 't' )   ! timestep
                    write(tt(jj+1:jj+ndigs(4)),fmt(4)) timestep
                    jj = jj + ndigs(4)
                 case default
                    tt(jj+1:jj+2) = '%'//t(j:j)
                    jj = jj + 2
               end select
            else
               tt(jj+1:jj+1) = '%'
               jj = jj + 1
            endif
         else
            tt(jj+1:jj+1) = t(j:j)
            jj = jj + 1
         endif
         j = j + 1
      enddo
      t = adjustl(tt)
!
!*   Get also rid of any blanks in the middle of the string
!
      len_t = len_trim(t)
      j = 1
      do while (j < len_t)
         if (t(j:j) == ' ') then
            t(j:) = t(j+1:)
            len_t = len_trim(t)
         else
            j = j + 1
         endif
      enddo
   endif
   s(i) = t
enddo
     
end subroutine expand_string

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_crtbl
     
! ---------------------------------------------------------------------------- !
!                                                                              !
!     mpi_crtbl - creates a table for input / output processors                !
!                                                                              !
!     A. Behrens   MSC/GKSS  November 2003   MPI parallelization (RPN_COMM)    !
!     E. Myklebust           February 2005   MPI parallelization               !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     purpose :
!     ---------
!
!     The subroutine creates a table that associates each integrated output
!     parameter with a specific process. This table is used then to spread
!     the retrieval of those parameter fields in grid format on as many
!     processing units as possible
!
!     method :
!     --------
!
!     assign a process to a parameter and start again with the first
!     one when run out of processing units.
!
!     externals :
!     -----------
!     
!       none
!
!     reference :
!     -----------
!
!       none
!
! ---------------------------------------------------------------------------- !
!
!     local variables :
!     -----------------
      
integer :: iflag

! ---------------------------------------------------------------------------- !
! 
!     Create table mapping integrated output parameter with a process rank.
!     ---------------------------------------------------------------------

comtime = MPI_WTIME()-comtime

if (allocated(ipfgtbl)) deallocate (ipfgtbl)
allocate (ipfgtbl(npout+1))

do iflag=1,npout
   ipfgtbl(iflag) = petotal             !! output integrated parameters
enddo

ipfgtbl(npout+1) = MAX(petotal-1,1)

i_out_par = petotal              !! process number for parameter output
i_out_spec = 1                   !! process number for spectra output
i_out_b_spec = MAX(petotal-1,1)  !! process number for boundary spectra output
i_out_restart = MAX(petotal-1,1) !! process number for restart field output

write (iu06,*)
write (iu06,*) '    sub. mpi_crtbl: Input / output processors'
write (iu06,*) ' process number for :'
write (iu06,*) ' parameter output           i_out_par     : ', i_out_par
write (iu06,*) ' spectra output             i_out_spec    : ', i_out_spec
write (iu06,*) ' boundary spectra output    i_out_b_spec  : ', i_out_b_spec
write (iu06,*) ' restart field input/output i_out_restart : ', i_out_restart

write (iu06,*) '    sub. mpi_crtbl: processors for output parameter computations'
write (iu06,*) ' ipfgtbl : ', ipfgtbl
write (iu06,*)

comtime = MPI_WTIME()-comtime
     
end subroutine mpi_crtbl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_decomp
    
! ---------------------------------------------------------------------------- !
!                                                                              !
!   mpi_decomp - developes an even decomposition of the grid domain among      !
!                the available processes                                       !
!                                                                              !
!    A. Behrens   MSC/GKSS      October 2003     MPI parallelization           !
!    E. Myklebust               November 2004    MPI parallelization           !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!    purpose :
!    ---------
!      
!    Determination of an even decomposition of the grid arrays for application
!    on a computer with distributed memory using the MPI message passing 
!    protocol for the exchange of information across the different processes.
!    Furthermore the length of the messages transferred between different
!    processes is computed.
!    It is assumed that each grid sub domain will send or receive messages
!    from the previous and the next process only.
!
!    method :
!    --------
!
!    Since the grid domain is mapped onto a one dimensional array of following
!    increasing latitude lines, sub grid domains are segments of that long
!    array. The length of those segments is determined by distributing the
!    total number of active sea points as even as possible.
!    The advection scheme uses neighboring grid points in the 2-d grid only
!    (see figure). Therefore for a given segment (process) information from 
!    previous grid points are required from a few grid points only starting
!    from the corresponding grid point in the 2-d grid just below (on the 
!    previous latitude line) of the starting point. This previous grid point
!    is a valid sea point - if not, then find the next one which is a valid
!    sea point but does not belong to the segment in question.
!    The same procedure is valid for the following grid points. The first
!    segment does not require any information from a previous segment.
!
!    figure : 2-d grid
!    -----------------
!
!       +    +    +    +    +    +    +    +
!       +    +    +    +    +    +    +    +
!       +    +    +    I    I    I    I    I
!       I    I    I    *    *    *    *    *
!       *    *    *    *    *    *    *    *
!       *    *    *    *    *    *    *    *
!       *    *    *    *    *    *    *    *
!       *    *    *    *    *    *    *    *
!       *    *    *    *    *    *    I    I
!       I    I    I    I    I    I    +    +
!       +    +    +    +    +    +    +    +
!       +    +    +    +    +    +    +    +
!
!    figure legend :
!    ---------------
!
!    *  grid points of the selected segment (process)
!    I  grid points of the previous or following segment for which 
!       information is necessary if the advection scheme is computed
!       at the * points only
!    +  other grid points that do not affect the * points
!
!
!    references :
!    ------------
!
!    none
!
! ---------------------------------------------------------------------------- !
!
!    externals :
!    -----------
!     
!    none
!
! ---------------------------------------------------------------------------- !
!
!    local variables :
!    -----------------
  
integer :: nmean, nrest, npts, ip, INBNGH
integer :: KMNOP, KMSOP, ij, NXDECOMP, NYDECOMP, NYCUT, IPROC, icount,NIJ, NAREA,IPR,IAR,NTOT
integer :: ISTAGGER, IIPER, KLATBOT,KLATTOP, KXLAT, NLONGMAX, IXLONMIN, IXLONMAX, IX, JSN
integer :: JC, KMIN, JCS, JCM, IIL, IC, ICL
integer :: MAXPERMLEN, MXPRLEN, IH, NH, JH
integer :: NLENHALO_MAX !! MAXIMUM NUMBER OF SEA POINTS CONTRIBUTING TO PROCESSES.

REAL :: XDELLOINV, STAGGER, XLON
INTEGER, ALLOCATABLE, DIMENSION(:) :: NSTART1D,NEND1D
INTEGER, ALLOCATABLE, DIMENSION(:) :: NEWIJ2IJ
INTEGER, ALLOCATABLE, DIMENSION(:) :: INDEX, IJNDEX ,NTOTSUB
INTEGER, ALLOCATABLE, DIMENSION(:) :: KSTART1, KEND1, NLON, ILON, KDUM
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IXLON
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IJFROMPE, IPROCFROM
INTEGER, ALLOCATABLE, DIMENSION(:) :: ITEMP,IJHALO
INTEGER, ALLOCATABLE, DIMENSION(:) :: NLENHALO   !! NUMBER OF SEA POINTS CONTRIBUTING TO PROCESSES.
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: KTEMP


INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: KDUM2
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KDUM3
REAL,    ALLOCATABLE, DIMENSION(:)     :: RDUM
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: RDUM2
REAL,    ALLOCATABLE, DIMENSION(:,:,:) :: R_HELP

CHARACTER (LEN=1)   :: CHAR
CHARACTER (LEN=100) :: TITL
CHARACTER (LEN=14)  :: ZERO = ' '

CHARACTER (LEN=1), DIMENSION(NX,NY) :: CGRID

comtime = MPI_WTIME()-comtime

! ---------------------------------------------------------------------------- !
!                                                                              !
!   1. Allocate arrary.                                                        !
!      ----------------                                                        !

allocate (nstart(petotal), nend(petotal), klenbot(petotal), klentop(petotal))
allocate (nlen(petotal))

IF (ALLOCATED(NFROMPE)) DEALLOCATE(NFROMPE)
ALLOCATE (NFROMPE(petotal))
IF (ALLOCATED(NTOPE)) DEALLOCATE(NTOPE)
ALLOCATE (NTOPE(petotal))
IF (ALLOCATED(NIJSTART)) DEALLOCATE(NIJSTART)
ALLOCATE (NIJSTART(petotal))

! ---------------------------------------------------------------------------- !
!                                                                              !
!   2. find the index of first (south) and last (north) latitude               !
!      containing sea points.                                                  !
!      ---------------------------------------------------------               !

KMNOP = MAXVAL(KXLT(1:NSEA))
KMSOP = MINVAL(KXLT(1:NSEA))

! ---------------------------------------------------------------------------- !
!                                                                              !
!   3. find the number of points per process, the start and the end index      !
!      ------------------------------------------------------------------      !

IF (L_DECOMP) THEN              !! 1D DECOMPOSITION ONLY (old)
   NXDECOMP = 1
   NYDECOMP = petotal
   NYCUT = petotal
ELSE                        !! 2D DECOMPOSITION (new)
   IF(petotal.EQ.1) THEN
      NXDECOMP=1
      NYDECOMP=1
      NYCUT=1
      L_DECOMP = .true.
   ELSEIF(petotal.EQ.2) THEN
      NXDECOMP=2
      NYDECOMP=1
      NYCUT=1
   ELSE

!   3.1 find whether petotal can be expressed as 2*i**2 i=1,2,3,...
!       because in that case
!         petotal = NXDECOMP*NYCUT + (NYDECOMP-NYCUT)*(NXDECOMP-1)  (1)
!       is satisfied with NXDECOMP=2*NYDECOMP and NYCUT=NYDECOMP
!       which is a perfect subdivision into identical squares

      IPROC=0
      ICOUNT=0
      DO WHILE (IPROC.LT.petotal)
         ICOUNT = ICOUNT+1
         IPROC = 2*ICOUNT**2
      ENDDO

      IF(IPROC.EQ.petotal) THEN
         NYDECOMP = INT(SQRT(FLOAT(petotal)/2))
         NXDECOMP = 2*NYDECOMP
         NYCUT = NYDECOMP
      ELSE

!   3.2 if the even decomposition into squares is not possible
!       start with the following approximation for NYDECOMP
!       found by setting NXDECOMP=2*NYDECOMP into (1) and
!       play around NXDECOMP=2*NYDECOMP,NYDECOMP,-1 and
!       NYCUT=NYDECOMP,1,-1 until a solution to (1) is reached.

         IPROC=0
         NYDECOMP = INT(SQRT(FLOAT(petotal)/2))+1
         DO NXDECOMP = 2*NYDECOMP,NYDECOMP,-1
            DO NYCUT = NYDECOMP,1,-1
               IPROC = NYDECOMP*(NXDECOMP-1)+ NYCUT
               IF (IPROC.EQ.petotal) EXIT
            ENDDO
            IF (IPROC.EQ.petotal) EXIT
         ENDDO
         IF (IPROC.NE.petotal) THEN
            WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!!'
            CALL ABORT1
         ENDIF
      ENDIF
   ENDIF
ENDIF

! ---------------------------------------------------------------------------- !
!                                                                              !
!   4. FIRST 1-D DECOMPOSITION IN LATITUDINAL BANDS

IF (ALLOCATED(NSTART1D)) DEALLOCATE(NSTART1D)
ALLOCATE(NSTART1D(NYDECOMP))
IF (ALLOCATED(NEND1D)) DEALLOCATE(NEND1D)
ALLOCATE(NEND1D(NYDECOMP))

IF(NYCUT.EQ.NYDECOMP) THEN

!   4.1 if the number of subareas per latitunal bands is the same
!       in all bands then the number of sea points in each band will
!       be determined to be as even as possible

   NMEAN=NSEA/NYDECOMP
   NREST=NSEA-NMEAN*NYDECOMP
     
   NSTART1D(1)=1
   IF (NREST.GT.0) THEN
      NPTS = NMEAN+1
      NREST = NREST-1
   ELSE
      NPTS=NMEAN
   ENDIF
   NEND1D(1)=NSTART1D(1)+NPTS-1

   DO IP = 2,NYDECOMP
      NSTART1D(IP) = NSTART1D(IP-1)+NPTS
      IF (NREST.GT.0) THEN
         NPTS = NMEAN+1
         NREST = NREST-1
      ELSE
         NPTS = NMEAN
      ENDIF
      NEND1D(IP) = NSTART1D(IP)+NPTS-1
   ENDDO

ELSE

!   4.2 if the number of subareas per latitunal bands is not the same
!       then the number of sea points per latitudinal bands will be
!       determined in such a way that number of points for the bands with
!       the least subareas (the top (nydecomp-nycut) bands) will
!       roughly scale like (nxdecomp-1)/nxdecomp the number of points
!       in the remaining bottom nycut bands

   NMEAN = NXDECOMP*NSEA/((NXDECOMP-1)*NYDECOMP+NYCUT)
      
   NSTART1D(1) = 1
   NPTS = NMEAN
   NEND1D(1) = NSTART1D(1)+NPTS-1

   DO IP = 2,NYCUT
      NSTART1D(IP) = NSTART1D(IP-1)+NPTS
      NPTS = NMEAN
      NEND1D(IP) = NSTART1D(IP)+NPTS-1
   ENDDO

   NMEAN =(NSEA-NEND1D(NYCUT))/(NYDECOMP-NYCUT)
   NREST =(NSEA-NEND1D(NYCUT))-NMEAN*(NYDECOMP-NYCUT)

   DO IP = NYCUT+1,NYDECOMP
      NSTART1D(IP)=NSTART1D(IP-1)+NPTS
      IF (NREST.GT.0) THEN
         NPTS=NMEAN+1
         NREST=NREST-1
      ELSE
         NPTS=NMEAN
      ENDIF
      NEND1D(IP)=NSTART1D(IP)+NPTS-1
   ENDDO
ENDIF

! ---------------------------------------------------------------------------- !
!                                                                              !
!   5. SECOND 1-D DECOMPOSITION IN EACH LATITUDINAL BAND

IF (L_DECOMP .OR. petotal.EQ.1) THEN
   DO IP=1,NYDECOMP
      NSTART(IP) = NSTART1D(IP)
      NEND(IP) = NEND1D(IP)
   ENDDO

ELSE

   ALLOCATE(NEWIJ2IJ(0:NSEA))
   IF(ALLOCATED(IJ2NEWIJ)) DEALLOCATE(IJ2NEWIJ)
   ALLOCATE(IJ2NEWIJ(0:NSEA))
   NEWIJ2IJ(0)=0
   IJ2NEWIJ(0)=0

   XDELLOINV=1.0/M_SEC_TO_DEG(XDELLO)
   IF (IPER) THEN
      IIPER = 1
   ELSE
      IIPER = 0
   END IF
   STAGGER=0.5*M_SEC_TO_DEG(AMOEAP-AMOWEP+IIPER*XDELLO)/NXDECOMP
   STAGGER=FLOAT(NINT(100*STAGGER))/100.
   ISTAGGER=NINT(STAGGER*XDELLOINV)

   IPROC = 0
   NIJ = 0
   DO IPR = 1,NYDECOMP                    !! loop over latitude bands

!   5.1 find no. sea points in all subareas in lat. band

      IPROC = IPROC+1
      NSTART(IPROC) = NIJ+1
      NTOT = NEND1D(IPR)-NSTART1D(IPR)+1  !! no. of points in latitude band

      IF (IPR.LE.NYCUT) THEN              !! no. of subarea in latitude band
         NAREA = NXDECOMP
      ELSE
         NAREA = NXDECOMP-1
      ENDIF

      IF (ALLOCATED(NTOTSUB)) DEALLOCATE(NTOTSUB)
      ALLOCATE(NTOTSUB(NAREA))          !! no. points in subareas in lat. band

      NMEAN = NTOT/NAREA
      NREST = NTOT-NMEAN*NAREA
      DO IAR=1,NAREA
         IF (NREST.GT.0) THEN
            NTOTSUB(IAR) = NMEAN+1
            NREST = NREST-1
         ELSE
            NTOTSUB(IAR)= NMEAN
         ENDIF
      ENDDO

!   5.2 find start and end sea point no. of each latitude

      KLATBOT = KXLT(NSTART1D(IPR))      !! bottom index lat. in band
      KLATTOP = KXLT(NEND1D(IPR))        !! top index at. in band
      ALLOCATE(KSTART1(KLATBOT:KLATTOP)) !! start seapoint no. of latitude
      KSTART1 = 0
      ALLOCATE(KEND1(KLATBOT:KLATTOP))   !! end sea point no. of latitude
      KEND1 = 0
      KXLAT = KLATBOT
      KSTART1(KXLAT) = NSTART1D(IPR)
      DO IJ = NSTART1D(IPR)+1,NEND1D(IPR)
         IF (KXLAT.LT.KXLT(IJ)) THEN
            KXLAT = KXLT(IJ)
            KSTART1(KXLAT) = IJ
            KEND1(KXLAT-1) = IJ-1
         ENDIF
      ENDDO
      KEND1(KLATTOP) = NEND1D(IPR)

!   5.3 find start and end sea point no. of each latitude


      ALLOCATE(NLON(KLATBOT:KLATTOP))   !! no. of longitudes on each latitude
      NLON(KLATBOT:KLATTOP)=0

      NLONGMAX=0
      DO KXLAT=KLATBOT,KLATTOP
         NLONGMAX = MAX(KEND1(KXLAT)-KSTART1(KXLAT)+1,NLONGMAX)
      ENDDO

      ALLOCATE(IXLON(NLONGMAX,KLATBOT:KLATTOP))  !! index of longitudes on each
                                                 !! latitude in a regular grid
      IXLONMAX=INT(M_SEC_TO_DEG(AMOWEP)*XDELLOINV)-1

      KXLAT=KLATBOT
      DO IJ=NSTART1D(IPR),NEND1D(IPR)
         IF(KXLAT.LT.KXLT(IJ)) THEN
            KXLAT = KXLT(IJ)
         ENDIF
         NLON(KXLAT) = NLON(KXLAT)+1
         IX = IXLG(IJ)
         JSN= KXLT(IJ)
         XLON = M_SEC_TO_DEG(AMOWEP+(IX-1)*ZDELLO(JSN))
         XLON = FLOAT(NINT(100*XLON))/100.
         IXLON(NLON(KXLAT),KXLAT) = NINT(XLON*XDELLOINV)
         IXLONMAX = MAX(IXLONMAX,IXLON(NLON(KXLAT),KXLAT))
      ENDDO

!   5.4 find sea point numbers along longitudes in a latitude band.

      ALLOCATE (IJNDEX(NTOT))   !! sea point numbers along longitudes

      ALLOCATE(ILON(KLATBOT:KLATTOP))
      ILON(KLATBOT:KLATTOP) = 1

      JC=0
      KMIN = KLATBOT
      DO WHILE(KMIN.GT.0)
         IXLONMIN = IXLONMAX+1
         KMIN=0
         DO KXLAT = KLATBOT,KLATTOP
            IF (ILON(KXLAT).LE.NLON(KXLAT)) THEN
               IF (IXLON(ILON(KXLAT),KXLAT).LT.IXLONMIN) THEN
                  KMIN=KXLAT
                  IXLONMIN=IXLON(ILON(KXLAT),KXLAT)
               ENDIF
            ENDIF
         ENDDO
         IF (KMIN.GT.0) THEN
            IJ = KSTART1(KMIN)+ILON(KMIN)-1
            JC = JC+1
            IJNDEX(JC) = IJ
            ILON(KMIN) = ILON(KMIN)+1
         ENDIF
      ENDDO

!   5.5 find which points belong to a subarea

      JCS=1
      IF(MOD(IPR,2).EQ.0) THEN
!         staggering
            JCM=1
            DO KXLAT=KLATBOT,KLATTOP
              IIL=1
              DO WHILE (IXLON(MIN(IIL,NLON(KXLAT)),KXLAT).LT.ISTAGGER  &
&                  .AND.   IIL.LE.NLON(KXLAT) .AND.   NLON(KXLAT).GT.0 )
                IIL=IIL+1
                JCM=JCM+1 
              ENDDO
            ENDDO
      ELSE
            JCM=1
      ENDIF

      IAR=1
      IC=0
      DO JC = JCM,NTOT
         NIJ = NIJ+1
         IC = IC+1
         IF (IC.EQ.NTOTSUB(IAR)) THEN
            NEND(IPROC) = NIJ
         ELSEIF (IC.GT.NTOTSUB(IAR)) THEN
              IC = 1
              IAR = IAR+1
              IPROC = IPROC+1
              NSTART(IPROC) = NIJ
         ENDIF
         IJ = IJNDEX(JC)
         NEWIJ2IJ(NIJ) = IJ
         IJ2NEWIJ(IJ) = NIJ
      ENDDO

      DO JC=JCS,JCM-1
            NIJ=NIJ+1
            IC=IC+1 
            IF(IC.EQ.NTOTSUB(IAR)) THEN
              NEND(IPROC)=NIJ
            ELSEIF(IC.GT.NTOTSUB(IAR)) THEN
              IC=1
              IAR=IAR+1
              IPROC=IPROC+1
              NSTART(IPROC)=NIJ
            ENDIF
            IJ=IJNDEX(JC)
            NEWIJ2IJ(NIJ)=IJ
            IJ2NEWIJ(IJ)=NIJ
      ENDDO

      DEALLOCATE(KSTART1)
      DEALLOCATE(KEND1)
      DEALLOCATE(NLON)
      DEALLOCATE(ILON)
      DEALLOCATE(IXLON)
      DEALLOCATE(IJNDEX)
      DEALLOCATE(NTOTSUB)

   ENDDO !! Loop over latidude bands

! ---------------------------------------------------------------------------- !
!
!   5.6 RELABEL KLAT, KLON, DEPTH, IXLG, KXLT, WLAT, OBSLAT, AND OBSLON

   DO ICL = 1,2
      DO IC = 1,2
         DO IJ = 1,NSEA
            IF(KLAT(IJ,IC,ICL).GT.0) KLAT(IJ,IC,ICL) = IJ2NEWIJ(KLAT(IJ,IC,ICL))
         ENDDO
      ENDDO
   ENDDO
   DO IC=1,2
      DO IJ = 1,NSEA
         IF(KLON(IJ,IC).GT.0) KLON(IJ,IC) = IJ2NEWIJ(KLON(IJ,IC))
      ENDDO
   ENDDO

   ALLOCATE(KDUM(NSEA))
   DO ICL = 1,2
      DO IC = 1,2
         KDUM(1:NSEA) = KLAT(NEWIJ2IJ(1:NSEA),IC,ICL)
         KLAT(1:NSEA,IC,ICL)=KDUM(1:NSEA)
      ENDDO
   ENDDO
   DO IC = 1,2
      KDUM(1:NSEA)=KLON(NEWIJ2IJ(1:NSEA),IC)
      KLON(1:NSEA,IC)=KDUM(1:NSEA)
   ENDDO

   KDUM(1:NSEA)=IXLG(NEWIJ2IJ(1:NSEA))
   IXLG(1:NSEA)=KDUM(1:NSEA)

   KDUM(1:NSEA)=KXLT(NEWIJ2IJ(1:NSEA))
   KXLT(1:NSEA)=KDUM(1:NSEA)

   DEALLOCATE(KDUM)

   ALLOCATE(RDUM(NSEA))
   RDUM(1:NSEA) = DEPTH_B(NEWIJ2IJ(1:NSEA))
   DEPTH_B(1:NSEA) = RDUM(1:NSEA)

   DO IC = 1,2
      RDUM(1:NSEA) = WLAT(NEWIJ2IJ(1:NSEA),IC)
      WLAT(1:NSEA,IC) = RDUM(1:NSEA)
   ENDDO

   IF (L_OBSTRUCTION) THEN
      DO ICL = 1,ML
         DO IC = 1,2
            RDUM(1:NSEA) = OBSLAT(NEWIJ2IJ(1:NSEA),IC,ICL)
            OBSLAT(1:NSEA,IC,ICL) = RDUM(1:NSEA)
         ENDDO
      ENDDO

      DO ICL = 1,ML
         DO IC = 1,2
            RDUM(1:NSEA) = OBSLON(NEWIJ2IJ(1:NSEA),IC,ICL)
            OBSLON(1:NSEA,IC,ICL) = RDUM(1:NSEA)
         ENDDO
      ENDDO
   ENDIF
   DEALLOCATE(RDUM)

END IF  !!   END IF L_DECOMP

DEALLOCATE(NSTART1D)
DEALLOCATE(NEND1D)

! ---------------------------------------------------------------------------- !
!
!     6. DETERMINE THE LENGTH OF THE MESSAGE THAT WILL BE EXCHANGED
!        BETWEEN NEIGHBORING SUB GRID DOMAINS
!        -----------------------------------------------------------

nlen(:) = nend(:)-nstart(:)+1   !! number of points on a process.
mpmaxlength = maxval(nlen(:))   !! maximum number of points on a process.


MAXPERMLEN = 2*MAX(mpmaxlength,NX)+12
MXPRLEN = 6*MAXPERMLEN

!     6.1 FIND THE NUMBER OF SEA POINT AND THE SEA POINT NUMBERS 
!         CONTRIBUTING TO PROCESSES

ALLOCATE(IJFROMPE(MAXPERMLEN,petotal))  !! SEA POINT NUMBERS CONTRIBUTING TO PROCESSES
ALLOCATE(NLENHALO(petotal))             !! NUMBER OF SEA POINT CONTRIBUTING TO PROCESSES

IJFROMPE(:,:) = 0
NLENHALO_MAX = 0

ALLOCATE(ITEMP(MXPRLEN))  !! SEA POINT NUMBERS CONTRIBUTING TO ONE PROCESS
ALLOCATE(INDEX(MXPRLEN))

DO IP = 1,petotal

!       CONTRIBUTION TO PROCESS IP FROM OTHER PROCESSORS

   NH = 0          !! NUMBER OF POINTS CONTRIBUTING TO PROCESS IP
   ITEMP(:) = 0    !! SEA POINT NUMBERS CONTRIBUTING TO PROCESS IP

   DO IC = 1,2
      DO IJ = NSTART(IP),NEND(IP)

         IF ( KLON(IJ,IC).GT.0 .AND. KLON(IJ,IC).LE.NSEA .AND.     &
&            (KLON(IJ,IC).LT.NSTART(IP) .OR. KLON(IJ,IC).GT.NEND(IP) ) ) THEN
            NH = NH + 1
            IF (NH.GT.MXPRLEN) THEN
               WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!'
               WRITE(IU06,*) 'MXPRLEN TOO SMALL !!!'
               CALL ABORT1
            ENDIF
            ITEMP(NH) = KLON(IJ,IC)
         ENDIF
      ENDDO ! END DO ON IJ
   ENDDO ! END DO ON IC

   DO ICL = 1,2
      DO IC = 1,2
         DO IJ = NSTART(IP),NEND(IP)

            IF ( KLAT(IJ,IC,ICL).GT.0 .AND. KLAT(IJ,IC,ICL).LE.NSEA .AND.     &
&               (KLAT(IJ,IC,ICL).LT.NSTART(IP).OR. KLAT(IJ,IC,ICL).GT.NEND(IP) ) ) THEN
                NH=NH+1
                IF(NH.GT.MXPRLEN) THEN
                  WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!'
                  WRITE(IU06,*) 'MXPRLEN TOO SMALL !!!'
                  CALL ABORT1
                ENDIF
                ITEMP(NH) = KLAT(IJ,IC,ICL)
             ENDIF

         ENDDO ! END DO ON IJ
      ENDDO ! END DO ON IC
   ENDDO ! END DO ON ICL

   IF(NH.GT.1) THEN                        !! sort points if more than one.
      CALL SORTINI(ITEMP(1:NH),INDEX(1:NH))
      CALL SORTI  (ITEMP(1:NH),INDEX(1:NH))
   ENDIF

   JH = 0
   IF (NH.GT.0) THEN
      JH = 1
      IJFROMPE(JH,IP) = ITEMP(1)
      DO IH = 2,NH
         IF (ITEMP(IH).GT.ITEMP(IH-1)) THEN      !! remove douple points
            JH = JH+1
            IF (JH.GT.MAXPERMLEN) THEN
               WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!'
               WRITE(IU06,*) 'MAXPERMLEN TOO SMALL !!!'
               WRITE(IU06,*) 'IH = ',IH,' JH = ',JH
               CALL ABORT1
            ENDIF
            IJFROMPE(JH,IP) = ITEMP(IH)
         ENDIF
      ENDDO
   ENDIF
   NLENHALO(IP) = JH
END DO        !! loop over processors

DEALLOCATE(ITEMP)
DEALLOCATE(INDEX)

! ---------------------------------------------------------------------------- !
!
!     7. RESIZE ARRAYS BASED ON MAXIMUM SIZE OF HALO.
!        --------------------------------------------

NLENHALO_MAX = MAXVAL(NLENHALO(:))
NLENHALO_MAX = MAX(1,NLENHALO_MAX)

ALLOCATE(KTEMP(NLENHALO_MAX,petotal))

DO IP=1,petotal
   DO IH=1,NLENHALO_MAX
      KTEMP(IH,IP )= IJFROMPE(IH,IP)
   ENDDO
ENDDO

DEALLOCATE(IJFROMPE)
ALLOCATE(IJFROMPE(NLENHALO_MAX,petotal))
DO IP = 1,petotal
   DO IH = 1,NLENHALO_MAX
      IJFROMPE(IH,IP) = KTEMP(IH,IP)
  ENDDO
ENDDO

DEALLOCATE(KTEMP)

! ---------------------------------------------------------------------------- !
!
!     8. DETERMINE IPROC FROM, KLENBOT and KLENTOP.
!        ------------------------------------------

ALLOCATE(IPROCFROM(NLENHALO_MAX,petotal))

IPROCFROM(:,:) = petotal+1   !! process numbers sending data

DO IP=1,petotal
   KLENBOT(IP)=0    !! number of points sent by the previous sub grid domains
   KLENTOP(IP)=0    !! number of points sent by the next sub grid domains
   DO IH = 1,NLENHALO(IP)
      DO IPROC=1,petotal
         IF ( IJFROMPE(IH,IP).GE.NSTART(IPROC) .AND.   &
&             IJFROMPE(IH,IP).LE. NEND(IPROC) ) THEN
            IPROCFROM(IH,IP) = IPROC
            IF (IPROC.LT.IP) THEN
               KLENBOT(IP) = KLENBOT(IP)+1
            ELSEIF(IPROC.GT.IP) THEN
               KLENTOP(IP) = KLENTOP(IP)+1
            ENDIF
            EXIT
         ENDIF
      ENDDO
   ENDDO
ENDDO

! ---------------------------------------------------------------------------- !
!
!     9. DETERMINE INFROMATION FOR THE LOCAL PROCESS.
!        --------------------------------------------

nijs = nstart(irank)   !! start index  of the sub grid domain of local process
nijl = nend(irank)     !! end index of the sub grid domain of local process
ninf = nstart(irank)-klenbot(irank)  !! first index used for propagation.
nsup = nend(irank)+klentop(irank)    !! last  index used for propagation.

NTOPE(1:petotal) = 0  !! number of points send from local pe to other pe's.

DO IP = 1,petotal
   DO IH = 1,NLENHALO(IP)
      IF (IPROCFROM(IH,IP).EQ.IRANK) THEN
         NTOPE(IP) = NTOPE(IP)+1
      ENDIF
   ENDDO
ENDDO
NTOPEMAX = MAXVAL(NTOPE) !! total number of points send by local pe.

NFROMPE(1:petotal)=0     !! number of points received by the local pe 
                         !! from other pe's.

DO IH=1,NLENHALO(IRANK)
   NFROMPE(IPROCFROM(IH,IRANK)) = NFROMPE(IPROCFROM(IH,IRANK))+1
ENDDO
NFROMPEMAX = MAXVAL(NFROMPE)   !! total number of points received
                               !! by the local pe.

NGBTOPE = 0    !! total number of pe's to which information will be
               !! sent by the local pe.
DO IP = 1,petotal
   IF (NTOPE(IP).GT.0) NGBTOPE = NGBTOPE+1
ENDDO

IF (ALLOCATED(NTOPELST)) DEALLOCATE(NTOPELST)
ALLOCATE(NTOPELST(NGBTOPE))    !! list of pe's to which information will
                               !! be sent by the local pe.
INBNGH = 0
DO IP = 1,petotal
   IF (NTOPE(IP).GT.0) THEN
      INBNGH = INBNGH+1
      NTOPELST(INBNGH) = IP
   ENDIF
ENDDO

NGBFROMPE = 0 !! total number of pe's from which information will be
              !! received by the local pe.
DO IP = 1,petotal
   IF (NFROMPE(IP).GT.0) NGBFROMPE = NGBFROMPE+1
ENDDO

IF (ALLOCATED(NFROMPELST)) DEALLOCATE(NFROMPELST)
ALLOCATE(NFROMPELST(MAX(1,NGBFROMPE)))  !! list of pe's from which
                                        !! information will be received
INBNGH=0
DO IP = 1,petotal
   IF(NFROMPE(IP).GT.0) THEN
     INBNGH = INBNGH+1
     NFROMPELST(INBNGH) = IP
   ENDIF
ENDDO

!     DETERMINE WHICH SEA POINTS NEED TO BE SEND TO THE OTHER PE'S

IF (ALLOCATED(IJTOPE)) DEALLOCATE(IJTOPE)
ALLOCATE (IJTOPE(NTOPEMAX,petotal))  !! sea point indices to be send
                                     !! by the local pe to other pe's.

IJTOPE(:,:) = NINF-1

DO IP=1,petotal
   JH = 0
   DO IH = 1,NLENHALO(IP)
      IF (IPROCFROM(IH,IP).EQ.IRANK) THEN
         JH = JH+1
         IJTOPE(JH,IP) = IJFROMPE(IH,IP)
      ENDIF
   ENDDO
ENDDO


ALLOCATE (IJHALO(MAX(1,NLENHALO(IRANK))))

DO IH=1,NLENHALO(IRANK)
   IF (IPROCFROM(IH,IRANK).LT.IRANK) THEN
      IJHALO(IH) = NINF+IH-1
   ELSEIF(IPROCFROM(IH,IRANK).GT.IRANK) THEN
      IJHALO(IH) = NEND(IRANK)+IH-KLENBOT(IRANK)
   ENDIF
ENDDO

!     CHANGE THE LOCAL ADDRESSING OF KLAT, KLON, AND DEPTH
!     FOR POINTS IN THE HALO.
!     NOTE THAT THIS IMPLIES THAT THESE ARRAYS ARE LOCAL
!     BECAUSE THEY ARE DIFFERENT IN THE HALO REGIONS

DO IC=1,2
   DO IJ = NSTART(IRANK),NEND(IRANK)
      DO IH=1,NLENHALO(IRANK)
         IF(KLON(IJ,IC) .EQ. IJFROMPE(IH,IRANK)) THEN
            KLON(IJ,IC) = IJHALO(IH)
            EXIT
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO ICL=1,2
   DO IC=1,2
      DO IJ=NSTART(IRANK),NEND(IRANK)
         DO IH=1,NLENHALO(IRANK)
            IF(KLAT(IJ,IC,ICL).EQ.IJFROMPE(IH,IRANK)) THEN
               KLAT(IJ,IC,ICL) = IJHALO(IH)
               EXIT
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

ALLOCATE(RDUM(MAX(1,NLENHALO(IRANK))))
DO IH = 1,NLENHALO(IRANK)
   IJ = IJFROMPE(IH,IRANK)
   RDUM(IH) = DEPTH_B(IJ)
ENDDO
DO IH = 1,NLENHALO(IRANK)
   IJ = IJHALO(IH)
   DEPTH_B(IJ) = RDUM(IH)
ENDDO
DEALLOCATE(RDUM)

!     CREATE IFROMIJ, JFROMIJ
!     !!!! IT IS ONLY DEFINED FOR GRID POINTS ON A GIVEN PE AND THEIR HALO !!!!

IF (ALLOCATED(IFROMIJ)) DEALLOCATE(IFROMIJ)
ALLOCATE(IFROMIJ(NINF-1:NSUP))
IF(ALLOCATED(KFROMIJ)) DEALLOCATE(KFROMIJ)
ALLOCATE(KFROMIJ(NINF-1:NSUP))

IFROMIJ(NINF-1)=0
KFROMIJ(NINF-1)=0

IF (petotal.GT.1) THEN
   DO IJ = NSTART(IRANK),NEND(IRANK)
      IFROMIJ(IJ) = IXLG(IJ)
      KFROMIJ(IJ) = KXLT(IJ)
   ENDDO
!       LOCAL HALO
   DO IH = 1,NLENHALO(IRANK)
      IJ = IJFROMPE(IH,IRANK)
      IFROMIJ(IJHALO(IH)) = IXLG(IJ)
      KFROMIJ(IJHALO(IH)) = KXLT(IJ)
   ENDDO
ELSE
   DO IJ = NINF, NSUP
      IFROMIJ(IJ) = IXLG(IJ)
      KFROMIJ(IJ) = KXLT(IJ)
   ENDDO
ENDIF


! ---------------------------------------------------------------------------- !
!
!     FIND INDEX OF THE FIRST HALO POINT FROM A GIVEN PE IN THE LOCAL
!      HALO BUFFERS WHICH PADS BOTH ENDS OF THE LOCAL GRID POINT 1-D ARRAY

DO IP = 1,petotal
   NIJSTART(IP) = NINF-1
ENDDO
IF (petotal.GT.1) THEN
   IF (IPROCFROM(1,IRANK).LT.IRANK) THEN
      NIJSTART(IPROCFROM(1,IRANK)) = NINF
   ELSEIF (IPROCFROM(1,IRANK).GT.IRANK) THEN
      NIJSTART(IPROCFROM(1,IRANK)) = NEND(IRANK)+1
   ENDIF
   DO IH = 2,NLENHALO(IRANK)
      IF (IPROCFROM(IH,IRANK).NE.IPROCFROM(IH-1,IRANK))THEN
         IF (IPROCFROM(IH,IRANK).LT.IRANK) THEN
            NIJSTART(IPROCFROM(IH,IRANK)) = NINF+IH-1
         ELSEIF (IPROCFROM(IH,IRANK).GT.IRANK) THEN
            NIJSTART(IPROCFROM(IH,IRANK)) = NEND(IRANK)+IH-KLENBOT(IRANK)
         ENDIF
      ENDIF
   ENDDO
ENDIF

DEALLOCATE(IJFROMPE)
DEALLOCATE(IPROCFROM)
DEALLOCATE(NLENHALO)
DEALLOCATE(IJHALO)

! ---------------------------------------------------------------------------- !
!
!     4. KEEP THE PART OF KLAT,KLON, DEPTH ,WLAT, WHICH IS NECESSARY
!        ----------------------------------------------------------


IF (petotal.GT.1) THEN

   ALLOCATE(KDUM3(NSTART(IRANK):NEND(IRANK),2,2))

   DO ICL=1,2
      DO IC=1,2
         DO IJ=NSTART(IRANK),NEND(IRANK)
            KDUM3(IJ,IC,ICL) = KLAT(IJ,IC,ICL)
         ENDDO
      ENDDO
   ENDDO

   DEALLOCATE(KLAT)
   ALLOCATE(KLAT(NSTART(IRANK):NEND(IRANK),2,2))

   DO ICL=1,2
      DO IC=1,2
         DO IJ=NSTART(IRANK),NEND(IRANK)
            KLAT(IJ,IC,ICL) = KDUM3(IJ,IC,ICL)
         ENDDO
      ENDDO
   ENDDO
   DEALLOCATE(KDUM3)

   ALLOCATE(KDUM2(NSTART(IRANK):NEND(IRANK),2))
   DO IC=1,2
      DO IJ=NSTART(IRANK),NEND(IRANK)
         KDUM2(IJ,IC) = KLON(IJ,IC)
      ENDDO
   ENDDO

   DEALLOCATE(KLON)
   ALLOCATE(KLON(NSTART(IRANK):NEND(IRANK),2))

   DO IC=1,2
      DO IJ=NSTART(IRANK),NEND(IRANK)
         KLON(IJ,IC) = KDUM2(IJ,IC)
      ENDDO
   ENDDO
   DEALLOCATE(KDUM2)

   ALLOCATE(RDUM(NINF:NSUP))
   DO IJ=NINF,NSUP
      RDUM(IJ)=DEPTH_B(IJ)
   ENDDO
   DEALLOCATE(DEPTH_B)
   ALLOCATE(DEPTH_B(NINF:NSUP))
   DO IJ = NINF,NSUP
      DEPTH_B(IJ)=RDUM(IJ)
   ENDDO
   DEALLOCATE(RDUM)

   ALLOCATE(RDUM2(NINF:NSUP,2))
   DO IC = 1,2
      DO IJ = NINF,NSUP
         RDUM2(IJ,IC) = WLAT(IJ,IC)
      ENDDO
   ENDDO
   DEALLOCATE(WLAT)
   ALLOCATE(WLAT(NINF:NSUP,2))
   DO IC = 1,2
      DO IJ = NINF,NSUP
         WLAT(IJ,IC) = RDUM2(IJ,IC)
      ENDDO
   ENDDO
   DEALLOCATE(RDUM2)

   IF (L_OBSTRUCTION) THEN
      ALLOCATE (R_HELP (ninf:nsup,1:2,1:ML))
      R_HELP (ninf:nsup,1:2,1:ML) = OBSLAT(ninf:nsup,1:2,1:ML)
      DEALLOCATE (OBSLAT)
      ALLOCATE (OBSLAT(ninf:nsup,1:2,1:ML))
      OBSLAT(ninf:nsup,1:2,1:ML) = R_HELP (ninf:nsup,1:2,1:ML)

      R_HELP (ninf:nsup,1:2,1:ML) = OBSLON(ninf:nsup,1:2,1:ML)
      DEALLOCATE (OBSLON)
      ALLOCATE (OBSLON(ninf:nsup,1:2,1:ML))
      OBSLON(ninf:nsup,1:2,1:ML) = R_HELP (ninf:nsup,1:2,1:ML)

      IF (ALLOCATED(R_HELP)) DEALLOCATE (R_HELP)
   END IF


ENDIF

!     5. MODIFY KLAT AND KLON SUCH THAT POINT INDICES FOR LAND IS
!        NINF-1.
!        ---------------------------------------------------------

IF (petotal.GT.1) THEN
   DO ICL=1,2
      DO IC=1,2
         DO IJ = NSTART(IRANK),NEND(IRANK)
            IF(KLAT(IJ,IC,ICL).EQ.0) KLAT(IJ,IC,ICL) = NINF-1
         ENDDO
      ENDDO
   ENDDO
   DO IC=1,2
      DO IJ = NSTART(IRANK),NEND(IRANK)
         IF(KLON(IJ,IC).EQ.0) KLON(IJ,IC) = NINF-1
      ENDDO
   ENDDO
ENDIF

npoi = nlen(irank)                   !! number of points on local process.

! ---------------------------------------------------------------------------- !
!  
!    4. Test print out 
!       --------------
     
if (itest>=0) then
   write (iu06,*) '   mpi_decomp: Sea point distribution to processors'
   do ip=1,petotal
      write (iu06,*)
      write (iu06,*) ' process number                          : ', ip
      write (iu06,*) ' first sea point number           nstart : ', nstart(ip)
      write (iu06,*) ' last sea point number              nend : ', nend(ip)
      write (iu06,*) ' number of sea points               nlen : ', nlen (ip)
      write (iu06,*) ' points from previous process,   klenbot : ', klenbot(ip)
      write (iu06,*) ' points from next     process,   klentop : ', klentop(ip)
      write (iu06,*) ' -------------------------------------------------------'
   enddo
   write (iu06,*)
   write (iu06,*)
   write (iu06,*) ' local process number                irank : ', irank
   write (iu06,*) ' first sea point number               nijs : ', nijs
   write (iu06,*) ' last sea point number                nijl : ', nijl
   write (iu06,*) ' number of sea points                 npoi : ', npoi
   write (iu06,*) ' first index used for propagation     ninf : ', ninf
   write (iu06,*) ' last  index used for propagation     nsup : ', nsup
   write (iu06,*) ' number of processes recieving data from local pe ngbtope: ', NGBTOPE
   write (iu06,*) ' process numbers recieving data    ntopelst: ', NTOPELST

   write (iu06,*) ' number of points send to other processors    ntope : ', NTOPE
   DO ij=1,NGBTOPE
     ip = NTOPELST(ij)
    write (iu06,*) ' sea point indices to be send by the local pe  IP IJTOPE : ', IP, IJTOPE(1:NTOPE(IP),ip)
   end do

   write (iu06,*) ' number of processes sending data ngbfrompe: ', NGBFROMPE
   write (iu06,*) ' process numbers sending data    nfrompelst: ', NFROMPELST
   write (iu06,*) ' number of points received to other processors   nfrompe : ', NFROMPE

DO ij=1,NGBFROMPE
   ip = NFROMPELST(ij)
   write (iu06,*) ' INDEX OF THE FIRST HALO POINT  from pe IP  NIJSTART : ',  ip, NIJSTART(ip)
end do
   write (iu06,*) ' ----------------------------------------------------------'
endif

if (itest>=5) then
   write (iu06,*)
   cgrid = 'L'
   do ip = 1,petotal
      write (char,'(i1)')  mod(ip,10)
      do ij = nstart(ip),nend(ip)
         cgrid(ixlg(ij),kxlt(ij)) = char
      enddo
   enddo
   TITL = '   mpi_decomp: Sea point distribution to processors'
   CALL PRINT_ARRAY (IU06, ZERO, TITL, CGRID,                                  &
&                       AMOWEP, AMOSOP, AMOEAP, AMONOP, NG_R=NLON_RG)

endif

comtime = MPI_WTIME()-comtime

! ---------------------------------------------------------------------------- !
!  
!     5. Create table for input / output processors.
!        ------------------------------------------- 
    
call mpi_crtbl
if (itest>=2)  write (iu06,*) '    sub. impi_decomp: mpi_crtbl done'

end subroutine mpi_decomp

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_exchng (fl)

! ---------------------------------------------------------------------------- !
!                                                                              !
!**  *mpi_exchng* - exchanges messages between the process pelocal and the     !
!**                 previous and next one.                                     !
!                                                                              !
!     A. Behrens   MSC/GKSS  November 2003   MPI parallelization (RPN_COMM)    !
!     E. Myklebust           November 2004   MPI parallelization               !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!*    purpose :
!     ---------
!    
!     This subroutine exchanges messages between one processes and all its 
!     neighours using the MPI message passing protocol.
!
!*    method :
!     --------
!
!     first send all messages then collect
!     them accordingly with a blocking receive.
!
!     references :
!     ------------
!
!     chapter 4 of :
!     Using MPI, Portable parallel programming with the message passing
!     interface. W.Cropp, E.Lusk, A.Skjellum, MIT Press 1995.
!
!     externals :
!     -----------
!
!       MPI_send
!       MPI_recv
!       MPI_wait
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

real, dimension (ninf:nsup,kl,ml), intent(inout) :: fl !! frequency spectra

! ---------------------------------------------------------------------------- !
!
!*    local variables :
!     -----------------

integer :: ir, iproc, ingb, ij, ih, kcount, m, k, ktag, ierr(1)
integer :: istatus(MPI_STATUS_SIZE,ngbtope+ngbfrompe)
integer :: nbufmax
integer, dimension(ngbtope+ngbfrompe) :: ireq

real, allocatable :: zcombufs(:,:)
real, allocatable :: zcombufr(:,:)

! ---------------------------------------------------------------------------- !

IF (petotal.LE.1) RETURN

extime = MPI_WTIME()-extime
ktag = 1

nbufmax = MAX(ntopemax,nfrompemax)*kl*ml
allocate (zcombufs(nbufmax,ngbtope))
allocate (zcombufr(nbufmax,ngbfrompe))

! ---------------------------------------------------------------------------- !
!
!      1. Pack send buffers for ngbtope neighbourin pe's
!        -----------------------------------------------

do ingb = 1, ngbtope
   iproc = ntopelst(ingb)
   kcount = 0
   do m = 1, ml
      do k = 1, kl
         do ih = 1, ntope(iproc)
            ij = ijtope(ih,iproc)
            kcount =  kcount +1
            zcombufs(kcount,ingb) = fl(ij,k,m)
         enddo
      enddo
   enddo
enddo

! ---------------------------------------------------------------------------- !
!
!      2. send and receive all messages.
!         -----------------------------

ir = 0

DO ingb = 1, ngbfrompe
   ir = ir + 1
   iproc = nfrompelst(ingb)
   kcount = ml*kl*nfrompe(iproc)
   call MPI_irecv(zcombufr(1,ingb), kcount, MPI_REAL, iproc-1, KTAG,            &
&                MPI_COMM_WORLD, ireq(ir), ierr)
   if (ierr(1)/=0) then
      write (iu06,*) ' +++ error: Sub. mpi_MPEXCHNG'
      write (iu06,*) ' +++ error: MPI_recv, ierr = ', ierr
      write (iu06,*) ' receiving from iproc = ', iproc
      call abort1
   endif
enddo

DO ingb = 1, ngbtope
   ir = ir + 1
   iproc = ntopelst(ingb)
   kcount = ml*kl*ntope(iproc)
   call MPI_iSend(zcombufs(1,ingb), kcount, MPI_REAL, iproc-1, KTAG,             &
&                MPI_COMM_WORLD, ireq(ir), ierr)
   if (ierr(1)/=0) then
      write (iu06,*) ' +++ error: Sub. mpi_MPEXCHNG'
      write (iu06,*) ' +++ error: MPI_send, ierr = ', ierr
      write (iu06,*) ' send to iproc = ', iproc
      call abort1
   endif
enddo

! ---------------------------------------------------------------------------- !
!
!      3. wait untill all messages are completed.
!         --------------------------------------

call MPI_Waitall(ir, ireq(1:ir),istatus, ierr)

if (ierr(1)/=0) then
   write (iu06,*) ' +++ error: Sub. mpi_exchng'
   write (iu06,*) ' +++ error: MPI_Wait, ierr = ', ierr
   call abort1
endif

! ---------------------------------------------------------------------------- !
!
!      4. Decode received buffers.
!        -------------------------

do ingb = 1, ngbfrompe
   iproc = nfrompelst(ingb)
   kcount = 0
   do m = 1, ml
      do k = 1,kl
         do ih = 1, nfrompe(iproc)
            ij = nijstart(iproc)+ih-1
            kcount =  kcount +1
            fl(ij,k,m) = zcombufr(kcount,ingb)
         enddo
      enddo
   enddo
enddo

deallocate (zcombufs)
deallocate (zcombufr)

extime = MPI_WTIME()-extime

end subroutine mpi_exchng

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_exchng_v (fl)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     mpi_exchng_v - exchanges messages between the process pelocal and the    !
!                    previous and next one.                                    !
!                                                                              !
!     A. Behrens   MSC/GKSS  November 2003   MPI parallelization (RPN_COMM)    !
!     E. Myklebust           November 2004   MPI parallelization               !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!*    purpose :
!     ---------
!
!     This subroutine exchanges messages between a process and its previous
!     and its next process using the MPI message passing protocol. For each
!     process 2 "send" and 2 "receive" have to be performed, except for the
!     first and the last process which require one of each only.
!
!*    method :
!     --------
!
!     using non-blocking send - first send all messages then collect
!     them accordingly with a blocking receive.
!
!     references :
!     ------------
!
!     chapter 4 of :
!     Using MPI, Portable parallel programming with the message passing
!     interface. W.Cropp, E.Lusk, A.Skjellum, MIT Press 1995.
!
!     externals :
!     -----------
!
!       MPI_send
!       MPI_recv
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

real, dimension (ninf:nsup), intent(inout) :: fl !! PARAMETER FIELD

! ---------------------------------------------------------------------------- !
!
!*    local variables :
!     -----------------

integer :: ir, iproc, ingb, ij, ih, kcount, ktag, ierr(1)
integer :: istatus(MPI_STATUS_SIZE,ngbtope+ngbfrompe)
integer :: nbufmax
integer, dimension(ngbtope+ngbfrompe) :: ireq

real, allocatable :: zcombufs(:,:)
real, allocatable :: zcombufr(:,:)

! ---------------------------------------------------------------------------- !

IF (petotal.LE.1) RETURN

extime = MPI_WTIME()-extime

nbufmax = MAX(ntopemax,nfrompemax)
allocate (zcombufs(nbufmax,ngbtope))
allocate (zcombufr(nbufmax,ngbfrompe))

! ---------------------------------------------------------------------------- !
!
!      1. Pack send buffers for ngbtope neighbourin pe's
!        -----------------------------------------------

do ingb = 1, ngbtope
   iproc = ntopelst(ingb)
   kcount = 0
   do ih = 1, ntope(iproc)
      ij = ijtope(ih,iproc)
      kcount =  kcount +1
      zcombufs(kcount,ingb) = fl(ij)
   enddo
enddo

! ---------------------------------------------------------------------------- !
!
!      2. send and receive all messages.
!         -----------------------------

ir = 0

ktag=0	!MPI_ANY_TAG
DO ingb = 1, ngbfrompe
   ir = ir + 1
   iproc = nfrompelst(ingb)
   kcount = nfrompe(iproc)
   call MPI_irecv(zcombufr(1,ingb), kcount, MPI_REAL, iproc-1, 0,             &
&                MPI_COMM_WORLD,ireq(ir), ierr)
   if (ierr(1)/=0) then
      write (iu06,*) ' +++ error: Sub. mpi_MPEXCHNG_V'
      write (iu06,*) ' +++ error: MPI_recv, ierr = ', ierr
      write (iu06,*) ' receiving from iproc = ', iproc
      call abort1
   endif
enddo

DO ingb = 1, ngbtope
   ir = ir + 1
   iproc = ntopelst(ingb)
   kcount = ntope(iproc)
   call MPI_iSend(zcombufs(1,ingb), kcount, MPI_REAL, iproc-1, 0,             &
&                MPI_COMM_WORLD,ireq(ir), ierr)
   if (ierr(1)/=0) then
      write (iu06,*) ' +++ error: Sub. mpi_MPEXCHNG_V'
      write (iu06,*) ' +++ error: MPI_send, ierr = ', ierr
      write (iu06,*) ' send to iproc = ', iproc
      call abort1
   endif
enddo
    
! ---------------------------------------------------------------------------- !
!
!      3. wait untill all messages are completed.
!         --------------------------------------

call MPI_Waitall(ir, ireq(1:ir),istatus, ierr)

if (ierr(1)/=0) then
   write (iu06,*) ' +++ error: Sub. mpi_exchng_V'
   write (iu06,*) ' +++ error: MPI_Wait, ierr = ', ierr
   call abort1
endif

! ---------------------------------------------------------------------------- !
!
!      4. Decode received buffers.
!        -------------------------

do ingb = 1, ngbfrompe
   iproc = nfrompelst(ingb)
   kcount = 0
   do ih = 1, nfrompe(iproc)
      ij = nijstart(iproc)+ih-1
      kcount =  kcount +1
      fl(ij) = zcombufr(kcount,ingb)
   enddo
enddo

deallocate (zcombufs)
deallocate (zcombufr)

extime = MPI_WTIME()-extime

end subroutine mpi_exchng_v

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_gather_block (irecv, field, rfield)
      
! ---------------------------------------------------------------------------- !
!                                                                              !
!     mpi_gather_block* - gather scalar block data field onto a single process !
!                                                                              !
!     A. Behrens   MSC/GKSS  November 2003   MPI parallelization (RPN_COMM)    !
!     E. Myklebust           November 2004   MPI parallelization               !
!     H. Bockelmann DKRZ     December 2011   use MPI_Gatherv instead of send/recv
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     purpose :
!     ---------
!
!     gather scalar block data contained in the array field, distributed
!     across the different processes onto the single process irecv
!
!     method :
!     --------
!
!     MPI_Gatherv with displacement by nstart
!
!     externals :
!     -----------
!    
!      MPI_Gatherv
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

integer,                    intent(in) :: irecv !! process rank receiving the grid field
real, dimension(nijs:nijl), intent(in) :: field !! containing the part be gathered
real, optional, dimension(1:nsea), intent(inout) :: rfield !! the gathered field
    
! ---------------------------------------------------------------------------- !
!
!*    local variables :
!     -----------------

integer :: ierr,i
integer, dimension(petotal) :: displ     ! HB TODO: dies koennte im mpi_module liegen

! ---------------------------------------------------------------------------- !

comtime = MPI_WTIME()-comtime

! ---------------------------------------------------------------------------- !
!
!    1.0 default action if field gathering is not required.
!         -------------------------------------------------

if (irecv==0.or.petotal==1) then
  comtime = MPI_WTIME()-comtime
  IF (.NOT. PRESENT(rfield)) THEN
    WRITE(iu06,*) ' +++ error: Sub. mpi_gather_block'
    WRITE(iu06,*) ' +++ error: gathering process i = ', irank,' does not posess rfield'
    CALL abort1
  END IF
  rfield = field
  return

else

  do i = 1,petotal
    displ(i) = nstart(i)-1
  end do

  CALL MPI_Gatherv(field, nlen(irank), MPI_REAL, rfield, nlen, displ,          &
                   MPI_REAL, irecv-1, MPI_COMM_WORLD, ierr)
  IF (ierr/=0) THEN
    write (iu06,*) ' +++ error: SUB. mpi_gather_block'
    write (iu06,*) ' +++ error: MPI-task ',irank,' -> MPI_Gatherv ierr = ',ierr
    call abort1
  END IF

  comtime = MPI_WTIME()-comtime
end if

end subroutine mpi_gather_block

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_gather_bound (irecv, itag, fl3, depth, fbc, depthbc)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      mpi_gather_bound - gathers spectrum onto a single process for output    !
!                         of boundary values                                   !
!                                                                              !
!     A. Behrens   MSC/GKSS  December 2003   MPI parallelization (RPN_COMM)    !
!     E. Myklebust           November 2004   MPI parallelization               !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     purpose :
!     ---------
!
!     The subroutine gathers the spectrum distributed across the different
!     processes onto the single process irecv for selected grid points.
!
!     method :
!     --------
!
!     MPI send of a message buffer containing the 2d-specctrum at all
!     selected grid points available at a local process to the process
!     corresponding to irecv where the message buffers are received and
!     reordered.
!
!     externals :
!     -----------
!
!      MPI_send
!      MPI_recv
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

integer, intent(in) :: irecv !! process rank receiving the grid field
integer, intent(in) :: itag  !! tag to differentiate calls to this subroutine

real, dimension (nijs:nijl,kl,ml),       intent(in)  :: fl3 !! block of spectra
real, dimension (ninf:nsup),             intent(in)  :: depth !! block of depth
real, dimension (max_nest,kl,ml,n_nest), intent(out) :: fbc !! spectra at
                                                    !! the selected grid points   
real, dimension (max_nest,n_nest), intent(out) :: depthbc !! depth at
                                                    !! the selected grid points

! ---------------------------------------------------------------------------- !
!
!*    local variables :
!     -----------------

integer :: maxlength, ngou, ij, k, m, ip, i, ij1
integer :: kcount, ierr, istatus(MPI_STATUS_SIZE),ir,isc
integer,allocatable,dimension(:)  :: ireq
real, allocatable, dimension (:)  :: zcombufs
real, allocatable, dimension (:,:):: zcombufr

! ---------------------------------------------------------------------------- !

comtime = MPI_WTIME()-comtime
maxlength = (1+kl*ml)*max_nest*n_nest

! ---------------------------------------------------------------------------- !
!
!*    1.0 default action if field gathering is not required
!         -------------------------------------------------

if (irecv==0.or.petotal==1) then
   do i = 1,n_nest
      do ngou=1,nbounc(i)
         ij = ijarc(ngou,i)
         depthbc(ngou,i) = depth(ij)
         do m=1,ml
            do k=1,kl
               fbc(ngou,k,m,i) = fl3(ij,k,m)
            enddo
         enddo
      enddo
   enddo

else if (irank/=irecv) then
   allocate(zcombufs(maxlength))

!*    1.1 send to the process that gathers the whole field
!         ------------------------------------------------
!
!*    load communication buffer
!
   kcount = 0
   do i = 1,n_nest
      do ngou=1,nbounc_ga(irank,i)
         ij = ijarc_ga(ngou,irank,i)
         kcount = kcount+1
         zcombufs(kcount) = depth(ij)
            do m=1,ml
               do k=1,kl
                  kcount = kcount+1
                  zcombufs(kcount) = fl3(ij,k,m)
               enddo
            enddo
      enddo   
   enddo   
!
!*    send buffer
!
   if(kcount>0)then
      call MPI_Send(zcombufs, kcount, MPI_REAL, irecv-1, itag,                 &
&                MPI_COMM_WORLD, ierr)

      if (ierr/=0) then
	 write (iu06,*) ' +++ error: Sub. mpi_gather_bound'
	 write (iu06,*) ' +++ error: MPI_send, ierr = ', ierr
	 call abort1
      endif
   endif

   deallocate (zcombufs)
else

!*    1.2.1  receive contribution to the field from the other processes
!            ----------------------------------------------------------
!
   allocate(zcombufr(maxlength,petotal))
   allocate(ireq(petotal))
   isc=0
   ireq=MPI_REQUEST_NULL
   do ip=1,petotal 
      if (ip==irecv) cycle
      if (sum(nbounc_ga(ip,:))==0) cycle
      isc=isc+1

      call MPI_iRecv(zcombufr(1,ip), maxlength, MPI_REAL, ip-1, itag, &
&                   MPI_COMM_WORLD, ireq(ip), ierr)
      
      if (ierr/=0) then
         write (iu06,*) ' +++ error: Sub. mpi_gather_bound'
         write (iu06,*) ' +++ error: MPI_irecv, ierr = ', ierr
         call abort1
      endif

   enddo

!     1.2.2  contribution from receiving process
!            -----------------------------------

   do i = 1,n_nest
      do ngou=1,nbounc_ga(irank,i)
         ij  = ngouc_ga(ngou,irank,i)
         ij1 = ijarc_ga(ngou,irank,i)
         depthbc(ij,i) = depth(ij1)
         do m=1,ml
            do k=1,kl
               fbc(ij,k,m,i) = fl3(ij1,k,m)
            enddo
         enddo
      enddo
   enddo

!*    1.2.1  receive contribution to the field from the other processes
!            ----------------------------------------------------------
!
   do ip=1,isc
      call MPI_Waitany(petotal, ireq, ir, istatus, ierr)
      if (ierr/=0) then
         write (iu06,*) ' +++ error: Sub. mpi_gather_bound'
         write (iu06,*) ' +++ error: MPI_waitany, ierr = '
         call abort1
      endif
!
!*    decode buffer
!
      kcount = 0
      do i = 1,n_nest
         do ngou=1,nbounc_ga(ir,I)
            ij = ngouc_ga(ngou,ir,i)
            kcount = kcount+1
            depthbc(ngou,i) = zcombufr(kcount,ir)
            do m=1,ml
               do k=1,kl
                  kcount = kcount+1
                  fbc(ij,k,m,i) = zcombufr(kcount,ir)
               enddo
            enddo
         enddo
      enddo
   enddo

   deallocate(ireq)
   deallocate (zcombufr)
endif
    
comtime = MPI_WTIME()-comtime

end subroutine mpi_gather_bound

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_gather_fl (irecv, itag, fl, rfl)
     
! ---------------------------------------------------------------------------- !
!                                                                              !
!     mpi_gather_fl - gather spectral field fl onto a single process           !
!                                                                              !
!     A. Behrens   MSC/GKSS  December 2003  MPI parallelization (RPN_COMM)     !
!     E. Myklebust           November 2004  MPI parallelization                !
!     H. Bockelmann DKRZ     December 2011  modified arguments to avoid unnec. !
!                                           tmp-arrays on sender site          !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     purpose :
!     ---------
!
!     gather array fl which is distributed across the different
!     processes onto the single process irecv
!
!     method :
!     --------
!     
!     MPI send of array fl to the process corresponding to irecv
!     for all processes except for the process corresponding to irecv
!     where it is received.
!
!     externals :
!     -----------
!  
!      MPI_send
!      MPI_recv
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

integer, intent(in) :: irecv !! process rank receiving the grid field
integer, intent(in) :: itag  !! tag to differentiate calls to this subroutine
real, dimension(nijs:nijl,kl,ml), intent(in) :: fl !! containing the part 
                                                   !! of the spectrum
real, optional, dimension(1:nsea,kl,ml), intent(inout) :: rfl !! the gathered 
                                                              !! spectrum
! ---------------------------------------------------------------------------- !
!
!     local variables :
!     -----------------
   
real, allocatable, dimension (:) :: zcombuf
integer :: kcount, mplength, len, ij, m, k, ip
integer :: ierr, istatus(MPI_STATUS_SIZE)
    
! ---------------------------------------------------------------------------- !

comtime = MPI_WTIME()-comtime
mplength = mpmaxlength*kl*ml

! ---------------------------------------------------------------------------- !
!
!     1.0 default action if gathering is not required
!         -------------------------------------------

if (irecv==0.or.petotal==1) then
   comtime = MPI_WTIME()-comtime
  IF (.NOT. PRESENT(rfl)) THEN
    WRITE(iu06,*) ' +++ error: Sub. mpi_gather_fl'
    WRITE(iu06,*) ' +++ error: gathering process i = ', irank,' does not posess rfl'
    CALL abort1
  END IF
  rfl = fl
  return

else if (irank/=irecv) then
      
!     1.1 send to the process that gathers the whole field
!         ------------------------------------------------
       
   mplength = mpmaxlength*kl*ml
   allocate (zcombuf(mplength))
   len = nlen(irank)*kl*ml
   kcount = 0
   do m=1,ml
      do k=1,kl
         do ij=nstart(irank),nend(irank)
            kcount = kcount+1
            zcombuf(kcount) = fl(ij,k,m)
         enddo
      enddo
   enddo
!
!*    send contribution to receiving pe
!
   call MPI_send(zcombuf,len, MPI_REAL, irecv-1, itag, MPI_COMM_WORLD, ierr)

   if (ierr<0) then
      write (iu06,*) ' +++ error: Sub. mpi_gather_fl'
      write (iu06,*) ' +++ error: MPI_send, ierr = ', ierr
      call abort1
   endif
    
else

!     1.2 receive contribution to the field from the other process
!         --------------------------------------------------------

!     check that collecting array is present

  IF (.NOT. PRESENT(rfl)) THEN
    WRITE(iu06,*) ' +++ error: Sub. mpi_gather_fl'
    WRITE(iu06,*) ' +++ error: gathering process i = ', irank,' does not posess rfl'
    CALL abort1
  END IF
!
!     get contribution from other pe's
!
   mplength = mpmaxlength*kl*ml
   allocate (zcombuf(mplength))

   do ip=1,petotal

      IF (ip==irecv) THEN
     
         rfl(nijs:nijl,:,:) = fl(:,:,:)  !! decode own part
      ELSE
   
         len = nlen(ip)*kl*ml
   
         CALL MPI_recv(zcombuf, len, MPI_REAL, ip-1, itag, MPI_COMM_WORLD,  &
&                      istatus, ierr)
         IF (ierr/=0) then
            write (iu06,*) ' +++ error: Sub. mpi_gather_fl'
            write (iu06,*) ' +++ error: MPI_recv, ierr = ', ierr
            call abort1
         END IF

!    decode buffer

         kcount = 0
         do m=1,ml
            do k=1,kl
               do ij=nstart(ip),nend(ip)
                  kcount = kcount+1
                  rfl(ij,k,m) = zcombuf(kcount)
               enddo
            enddo
         enddo
      endif
   enddo

endif

deallocate (zcombuf)
comtime = MPI_WTIME()-comtime

end subroutine mpi_gather_fl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_gather_grid (isend, irecv, itag, field)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    mpi_gather_grid - gather grid data field from the process isend           !
!                      onto the process irecv                                  !
!                                                                              !
!     A. Behrens   MSC/GKSS  November 2003   MPI parallelization (RPN_COMM)    !
!     E. Myklebust           November 2004   MPI parallelization               !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     purpose :
!     ---------
!
!     transmit grid (2-d) data contained in the array field that might be
!     on a different processes to the receiving process
!
!     method :
!     --------
!
!     MPI send to process irecv from process isend
!
!     externals :
!     -----------
!        
!       MPI_send
!       MPI_recv
!   
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

integer, intent(in) :: isend !! process rank sending the grid field
integer, intent(in) :: irecv !! process rank receiving the grid field
integer, intent(in) :: itag  !! tag to differentiate calls to this subroutine

real, dimension (nx,ny), intent(inout) :: field !! array containing the field

! ---------------------------------------------------------------------------- !
!
!     local variables :
!     -----------------

integer :: len, ierr, ij, j, i
integer :: istatus(MPI_STATUS_SIZE)

real, dimension (nx*ny) :: zcombuf

! ---------------------------------------------------------------------------- !

comtime = MPI_WTIME()-comtime
len = nx*ny

! ---------------------------------------------------------------------------- !
!
!     1.0 default action if field gathering is not required
!         -------------------------------------------------
     
if (isend==irecv.or.petotal==1) then
   comtime = MPI_WTIME()-comtime
   return

else if (irank==isend) then

!     1.1 send field from process isend to process irecv
!         -----------------------------------------------
!
!*    load communication buffer
!
   ij = 0
   do j=1,ny
      do i=1,nx
         ij = ij+1
         zcombuf(ij) = field(i,j)
      enddo
   enddo
      
   call MPI_Send(zcombuf, len, MPI_REAL, irecv-1, itag, MPI_COMM_WORLD, ierr)

   if (ierr/=0) then
      write (iu06,*) ' +++ error: Sub. mpi_gather_grid'
      write (iu06,*) ' +++ error: MPI_send, ierr = ', ierr
      call abort1
   endif
    
else if (irank==irecv) then

!     1.2 receive field from process isend on process irecv
!         -------------------------------------------------

   call MPI_recv(zcombuf, len, MPI_REAL, isend-1, itag, MPI_COMM_WORLD,        &
&                istatus, ierr)
    
   if (ierr/=0) then
      write (iu06,*) ' +++ error: Sub. mpi_gather_grid'
      write (iu06,*) ' +++ error: MPI_recv, ierr = ', ierr
      call abort1
   endif
!
!*    decode communication buffer
!
   ij = 0
   do j=1,ny
      do i=1,nx
         ij = ij+1
         field(i,j) = zcombuf(ij)
      enddo
   enddo
endif
    
comtime = MPI_WTIME()-comtime

end subroutine mpi_gather_grid

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_gather_spp (irecv, itag, nspfld, nscfld, fl3, fl1, block,       &
&                          flpts, BLOCK_sp)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    mpi_gather_spp - gathers spectrum and scalar fields onto a single         !
!                     process for output at selected grid points               !
!                                                                              !
!     A. Behrens   MSC/GKSS  December 2003   MPI parallelization (RPN_COMM)    !
!     E. Myklebust           November 2004   MPI parallelization               !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     purpose :
!     ---------
!
!     gather nspfld spectra and nscfld scalar parameter distributed across the
!     different processes onto the single process irecv for selected grid 
!     points.
!
!     method :
!     --------
!
!     MPI send of a message buffer containing all the fields at all
!     selected grid points held by that process to the process corresponding
!     to irecv for all processes,  except for the process corresponding to
!     irecv where the message buffers will be received and reordered.
!
!     externals :
!     -----------
!     
!      MPI_send
!      MPI_recv
!      MPI_WTIME
!      MPI_Get_Count
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!     -----------                                                              !
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN) :: irecv  !! rank of the process onto which field 
                              !! is collected
INTEGER, INTENT(IN) :: itag   !! tag associated with a particular call to the 
                              !! subroutine.T his is necessary to differentiate
                              !! the different calls. The use of a separate
                              !! tag for each specific call should ensure that
                              !! no conflict arises between different messages.
INTEGER, INTENT(IN) :: nspfld !! number of spectra
INTEGER, INTENT(IN) :: nscfld !! number of scalar fields

real, dimension (nijs:nijl,kl,ml),    INTENT(IN)  :: fl3   !! spectra
real, dimension (nijs:nijl,kl,ml),    INTENT(IN)  :: fl1   !! swell spectra
real, dimension (nijs:nijl,nscfld),   INTENT(IN)  :: block !! INTEGRATED PARAMETER
real, optional, dimension (kl,ml,noutp,nspfld), INTENT(OUT) :: flpts 
                              !! spectra and parameter at the selected 
                              !! grid points
real, optional, dimension (noutp,nscfld), INTENT(OUT)  :: block_sp !! INTEGRATED PARAMETER

! ---------------------------------------------------------------------------- !
!
!     local variables.
!     ----------------

integer :: istatus(MPI_STATUS_SIZE)
integer :: ij, ij1, kcount, ierr, ip, rcount, klml
integer :: ngou, maxlength
real, allocatable, dimension (:) :: zcombuf

! ---------------------------------------------------------------------------- !

comtime = MPI_WTIME()-comtime

maxlength = (nscfld+nspfld*kl*ml)*noutp
allocate (zcombuf(maxlength))

! ---------------------------------------------------------------------------- !
!
!     1.0 default action if no field gathering
!         ------------------------------------

if (irecv==0.or.petotal==1) then
   do ngou=1,noutp
      ij = ijar(ngou)
      flpts(:,:,ngou,1) = fl3(ij,:,:)
      block_sp(ngou,:) = block(ij,:)
   enddo
   If (nspfld==2) then
      do ngou=1,noutp
         ij = ijar(ngou)
         flpts(:,:,ngou,2) = fl1(ij,:,:)
      enddo
   end if

else if (irank/=irecv .and. noutp_ga(irank).gt.0) then

!     1.1 send to the process that gathers the whole field
!         ------------------------------------------------
!
!     load communication buffer

   klml= kl*ml
   kcount = 0
   do ngou=1,noutp_ga(irank)
      ij = ijar_ga(ngou,irank)

      zcombuf(kcount+1:kcount+nscfld) = BLOCK(ij,1:nscfld)
      kcount = kcount + nscfld

      zcombuf(kcount+1:kcount+klml) = RESHAPE(fl3(ij,:,:), (/klml/))
      kcount = kcount+ml*kl
 
      If (nspfld==2) then
         zcombuf(kcount+1:kcount+klml) = RESHAPE(fl1(ij,:,:), (/klml/))
         kcount = kcount+klml
      endif
   enddo
!
!     send buffer
!
   call MPI_send(zcombuf, kcount, MPI_REAL, irecv-1, itag,            &
&                MPI_COMM_WORLD, ierr)

   if (ierr/=0) then
      write (iu06,*) ' +++ error : MPI_send, ierr = ', ierr
      call abort1
   endif

else if (irank==irecv) then

!     1.2.1  receive contribution to the field from the other processes
!            ----------------------------------------------------------

   klml= kl*ml
   Process: do ip = 1,petotal
      if (ip==irecv .or. noutp_ga(ip).eq.0) cycle Process

      call MPI_recv(zcombuf, maxlength, MPI_REAL, ip-1, itag,           &
&                   MPI_COMM_WORLD, istatus, ierr)
      if (ierr/=0) then
         write (iu06,*) ' +++ error: Sub. mpi_gather_spp'
         write (iu06,*) ' +++ error: MPI_recv, ierr = ', ierr
        call abort1
      endif

!    decode buffer

      kcount = 0
      do ngou=1,noutp_ga(ip)
         ij = ngou_ga(ngou,ip)

         block_sp(ij,1:nscfld) = zcombuf(kcount+1:kcount+nscfld)
         kcount = kcount+nscfld
 
         flpts(:,:,ij,1) = RESHAPE(zcombuf(kcount+1:kcount+klml), (/kl,ml/))
         kcount = kcount+klml

         If (nspfld==2) then
           flpts(:,:,ij,2) = RESHAPE(zcombuf(kcount+1:kcount+klml), (/kl,ml/))
           kcount = kcount+klml
         endif
      enddo
      call MPI_Get_Count(istatus, MPI_REAL, rcount, ierr)
      if (kcount/=rcount) then
         write (iu06,*) ' +++ receive-error im mpi_gather_spp : '
         write (iu06,*) ' +++ wrong message length'
         call abort1
      endif     
   enddo Process

!     1.2.2  contribution from receiving process
!            -----------------------------------
!
   do ngou=1,noutp_ga(irank)
     ij  = ngou_ga(ngou,irank)
     ij1 = ijar_ga(ngou,irank)
     block_sp(ij,:) = block(ij1,:)
     flpts(:,:,ij,1) = fl3(ij1,:,:)
     If (nspfld==2)  flpts(:,:,ij,2) = fl1(ij1,:,:)
   enddo
endif

deallocate (zcombuf)
comtime = MPI_WTIME()-comtime

end subroutine mpi_gather_spp
    
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_gather_oifl_real (isend, field)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     mpi_gather_oifl - send scalar grid data field onto all processes         !
!                                                                              !
!     A. Behrens   GKSS     September 2010                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     purpose :
!     ---------
!
!     broadcast of scalar grid data contained in array field stored            !
!     on process isend onto all processes                                      !
!
!     method :
!     --------
!
!     mpi broadcast of array field to all processes
!
!     externals :
!     -----------
!
!       mpi_broadcast
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     interface variables
!     -------------------

integer, intent(in) :: isend           !! process rank sending the grid field

real, dimension (nx,ny), intent(inout) :: field !! array containing the field

! ---------------------------------------------------------------------------- !
!
!     local variables.
!     ----------------

integer :: len, ierr, ij, j, i

real, dimension (nx*ny) :: zcombuf

! ---------------------------------------------------------------------------- !
!
!     1.0 broadcast field from process isend.
!         -----------------------------------

comtime = MPI_WTIME()-comtime

if (petotal/=1) then

   len = nx*ny
!
!        load communication buffer

   ij = 0
   do j=1,ny
      do i=1,nx
         ij = ij+1
         zcombuf(ij) = field(i,j)
      enddo
   enddo

   call MPI_Bcast (zcombuf, len, MPI_real, isend-1, MPI_COMM_WORLD, ierr)

   if (ierr/=0)  then
      write (iu06,*) ' +++ error: Sub. mpi_gather_oifl'
      write (iu06,*) ' +++ error: MPI_Bcast, ierr = ', ierr
      call abort1
   endif
!
!*    1.2 decode communication buffer on all processes
!         --------------------------------------------
!
   ij = 0
   do j=1,ny
      do i=1,nx
         ij = ij+1
         field(i,j) = zcombuf(ij)
      enddo
   enddo
endif

comtime = MPI_WTIME()-comtime
return
end subroutine mpi_gather_oifl_real

    
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_gather_oifl_logical (isend, field)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     mpi_gather_oifl - send logical grid data field onto all processes         !
!                                                                              !
!     A. Behrens   GKSS     September 2010                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     purpose :
!     ---------
!
!     broadcast of scalar grid data contained in array field stored            !
!     on process isend onto all processes                                      !
!
!     method :
!     --------
!
!     mpi broadcast of array field to all processes
!
!     externals :
!     -----------
!
!       mpi_broadcast
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     interface variables
!     -------------------

integer, intent(in) :: isend           !! process rank sending the grid field

logical, dimension (nx,ny), intent(inout) :: field !! array containing the field

! ---------------------------------------------------------------------------- !
!
!     local variables.
!     ----------------

integer :: len, ierr, ij, j, i

logical, dimension (nx*ny) :: zcombuf

! ---------------------------------------------------------------------------- !
!
!     1.0 broadcast field from process isend.
!         -----------------------------------

comtime = MPI_WTIME()-comtime

if (petotal/=1) then

   len = nx*ny
!
!        load communication buffer

   ij = 0
   do j=1,ny
      do i=1,nx
         ij = ij+1
         zcombuf(ij) = field(i,j)
      enddo
   enddo

   call MPI_Bcast (zcombuf, len, MPI_logical, isend-1, MPI_COMM_WORLD, ierr)

   if (ierr/=0)  then
      write (iu06,*) ' +++ error: Sub. mpi_gather_oifl'
      write (iu06,*) ' +++ error: MPI_Bcast, ierr = ', ierr
      call abort1
   endif
!
!*    1.2 decode communication buffer on all processes
!         --------------------------------------------
!
   ij = 0
   do j=1,ny
      do i=1,nx
         ij = ij+1
         field(i,j) = zcombuf(ij)
      enddo
   enddo
endif

comtime = MPI_WTIME()-comtime
return
end subroutine mpi_gather_oifl_logical

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine mpi_gather_cfl (field)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     mpi_gather_cfl - send a number onto all processes                        !
!                                                                              !
!     A. Behrens   GKSS     September 2010                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     purpose :
!     ---------
!
!     broadcast of a number contained in array field(irank) stored             !
!     on process irank onto all processes                                      !
!
!     method :
!     --------
!
!     mpi broadcast of array field to all processes
!
!     externals :
!     -----------
!
!       mpi_broadcast
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     interface variables
!     -------------------

real,    intent(inout) :: field(petotal) !! array containing the field

! ---------------------------------------------------------------------------- !
!
!     local variables.
!     ----------------

integer :: len, ierr, isend
real, dimension (1:1) :: zcombuf

! ---------------------------------------------------------------------------- !
!
!     1.0 broadcast field from process isend.
!         -----------------------------------

comtime = MPI_WTIME()-comtime

if (petotal/=1) then

   len = 1
   do isend=1,petotal
      zcombuf(1) = field (irank)
      call MPI_Bcast (zcombuf, len, MPI_real, isend-1,MPI_COMM_WORLD, ierr)

      if (ierr/=0)  then
         write (iu06,*) ' +++ error: Sub. mpi_gather_cfl'
         write (iu06,*) ' +++ error: MPI_Bcast, ierr = ', ierr
         call abort1
      endif
      field(isend) = zcombuf(1)
   enddo
endif

comtime = MPI_WTIME()-comtime

end subroutine mpi_gather_cfl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_MPI_COMP_MODULE

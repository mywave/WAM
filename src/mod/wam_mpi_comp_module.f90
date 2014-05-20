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
&       ABORT1                        !! TERMINATES PROCESSING.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

use wam_file_module,          only: iu06, itest
use wam_fre_dir_module,       only: ml, kl
use wam_grid_module,          only: nx, ny, nsea, klat
use wam_output_set_up_module, only: nout_P, cflag_P, noutp, ijar, npout
use wam_nest_module,          only: n_nest, max_nest, nbounc, ijarc
use wam_mpi_module,           only: petotal, irank, nstart, nend, nlen,        &
&                                   klentop, klenbot, mpmaxlength,             &
&                                   nnext, nprevious, ninf, nsup, nijs, nijl,  &
&                                   ipfgtbl, i_out_par, i_out_spec, i_out_rad, &
&                                   i_out_b_spec, i_out_restart,               &
&                                   noutp_ga, ijar_ga, ngou_ga,                &
&                                   extime, comtime

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
END INTERFACE
PUBLIC mpi_exchng

INTERFACE mpi_exchng_V               !! exchanges messages between the process
   MODULE PROCEDURE mpi_exchng_V     !! pelocal and the previous and next one
END INTERFACE                        !! PARAMETER
PUBLIC mpi_exchng_V

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
   module procedure mpi_gather_oifl  !! broadcast of scalar grid data contained
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
i_out_rad = petotal              !! process number for radiation stress output
i_out_b_spec = MAX(petotal-1,1)  !! process number for boundary spectra output
i_out_restart = MAX(petotal-1,1) !! process number for restart field output

write (iu06,*)
write (iu06,*) '    sub. mpi_crtbl: Input / output processors'
write (iu06,*) ' process number for :'
write (iu06,*) ' parameter output           i_out_par     : ', i_out_par
write (iu06,*) ' spectra output             i_out_spec    : ', i_out_spec
write (iu06,*) ' radiation stress output    i_out_rad     : ', i_out_spec 
write (iu06,*) ' boundary spectra output    i_out_b_spec  : ', i_out_rad  
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
  
integer :: nmean, nrest, npts, ip, i, isb, iet, n
    
! ---------------------------------------------------------------------------- !
!  
!   1. Allocate arrary 
!      ------------------------------------------------------------------

allocate (nstart(petotal), nend(petotal), klenbot(petotal), klentop(petotal))
allocate (nlen(petotal))

! ---------------------------------------------------------------------------- !
!  
!   2. find the number of points per process, the start and the end index
!      ------------------------------------------------------------------

comtime = MPI_WTIME()-comtime
nmean = nsea/petotal
nrest = nsea-nmean*petotal

nstart(1) = 1               !! index of the first point of the sub grid domain
if (nrest>0) then
   npts = nmean+1
   nrest = nrest-1
else
   npts = nmean
endif
nend(1) = nstart(1)+npts-1  !! index of the last point of the sub grid domain
 
do ip=2,petotal
   nstart(ip) = nstart(ip-1)+npts
   if (nrest>0) then
      npts = nmean+1
      nrest = nrest-1
   else
      npts = nmean
   endif
   nend(ip) = nstart(ip)+npts-1
enddo

nijs = nstart(irank)   !! start index  of the sub grid domain of local process  
nijl = nend(irank)     !! end index of the sub grid domain of local process   

! ---------------------------------------------------------------------------- !
!  
!    3. determine the length of the message that will be 
!       exchanged between neighbouring sub grid domains
!       ------------------------------------------------
   
klenbot(1) = 0     !! length of the message sent by the previous sub grid domain
do ip=2,petotal
   i = -1
20 i = i+1
   isb = klat(nstart(ip)+i,1)
   if (isb<=0) goto 20
   klenbot(ip) = max (0,nstart(ip)-isb)
enddo
     
klentop(petotal) = 0  !! length of the message sent by the next sub grid domain
do ip=1,petotal-1
   i = -1
40 i = i+1
   iet = klat(nend(ip)-i,2)
   if (iet<=0) goto 40
   klentop(ip) = max (0,iet-nend(ip))
enddo

nlen(:) = nend(:)-nstart(:)+1   !! number of points on a process.  
mpmaxlength = maxval(nlen(:))   !! maximum number of points on a process.
n = minval(nlen(:))             !! minimum number of points on a process.

ninf = nstart(irank)-klenbot(irank)  !! first index used for propagation
nsup = nend(irank)+klentop(irank)    !! last  index used for propagation

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
   write (iu06,*) ' local process number               irank : ', irank
   write (iu06,*) ' first sea point number              nijs : ', nijs
   write (iu06,*) ' last sea point number               nijl : ', nijl
   write (iu06,*) ' first index used for propagation    ninf : ', ninf
   write (iu06,*) ' last  index used for propagation    nsup : ', nsup
   write (iu06,*) ' -------------------------------------------------------'
endif

if (n<nx) then
   write (iu06,*)
   write (iu06,*) ' *******************************************************'
   write (iu06,*) ' ***  Error : not enough active grid points for one PE'
   write (iu06,*) ' ***  Grid points per PE : ', n, '   NX = ', nx
   write (iu06,*) ' ***  Number of grid points must be greater than NX !'
   write (iu06,*) ' ***  Please reduce number of processing units '
   write (iu06,*) ' *******************************************************'
   write (iu06,*)
   call abort1
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

real, dimension (ninf:nsup,kl,ml), intent(inout) :: fl !! frequency spectra

! ---------------------------------------------------------------------------- !
!
!*    local variables :
!     -----------------

logical, save :: frstime=.true.
integer, save :: mpitype_sn, mpitype_sp,  mpitype_rn, mpitype_rp
integer, save :: requests(4), nrequests
integer       :: istatus(MPI_STATUS_SIZE,4), ierr

! ---------------------------------------------------------------------------- !

extime = MPI_WTIME()-extime

if (frstime) then
   call exchng_init()
   frstime = .false.
end if
       
call MPI_StartAll(nrequests,requests,ierr)
if (ierr/=0) then
   write (iu06,*) ' +++ error: Sub. mpi_exchng'
   write (iu06,*) ' +++ error: MPI_Start_All, ierr = ', ierr
   call abort1
endif

!
!  It might be efficient to put some work here that does not 
!  depend on the exchanged data.
!

call MPI_WaitAll(nrequests,requests,istatus,ierr)
if (ierr/=0) then
   write (iu06,*) ' +++ error: Sub. mpi_exchng'
   write (iu06,*) ' +++ error: MPI_Wait_All, ierr = ', ierr
   call abort1
endif

extime = MPI_WTIME()-extime
return

! ---------------------------------------------------------------------------- !

contains

! ---------------------------------------------------------------------------- !

subroutine exchng_init

integer :: ist, isb, ierrs(12), ktag

ktag = 1

!*    0.3 find index of the start or end point of the different messages
!     ------------------------------------------------------------------

if (nnext/=0) then
   ist = nend(irank)-klenbot(irank+1)+1
endif
isb = nstart(irank)-klenbot(irank)

!*    1.0 Initialize communication
!

nrequests = 0
ierrs=0

if (nnext-1>0.and.nnext-1<petotal) then
!
!*     Send to next process
!
   call MPI_Type_Vector(ml*kl, klenbot(nnext), nsup-ninf+1, MPI_REAL, &
        & mpitype_sn, ierrs(1))
   call MPI_Type_Commit(mpitype_sn,ierrs(3))
   call MPI_Send_Init(fl(ist,1,1),1,mpitype_sn,nnext-1,ktag, &
        &  MPI_COMM_WORLD,requests(nrequests+1),ierrs(5))
!
!*     Receive from next process
!
   call MPI_Type_Vector(ml*kl, klentop(irank), nsup-ninf+1, &
        & MPI_REAL, mpitype_rn, ierrs(2))
   call MPI_Type_Commit(mpitype_rn,ierrs(4))
   call MPI_Recv_Init(fl(nend(irank)+1,1,1),1,mpitype_rn,nnext-1,  &
&                  ktag+1,MPI_COMM_WORLD,requests(nrequests+2),ierrs(6))
   nrequests = nrequests+2
end if

if (nprevious-1>=0) then
!
!*     Send to previous process
!
   call MPI_Type_Vector(ml*kl, klentop(nprevious), nsup-ninf+1, &
        & MPI_REAL, mpitype_sp, ierrs(7))
   call MPI_Type_Commit(mpitype_sp,ierrs(9))
   call MPI_SEND_Init(fl(nstart(irank),1,1),1,mpitype_sp,nprevious-1, &
&                  ktag+1,MPI_COMM_WORLD,requests(nrequests+1),ierrs(11))
!
!*     Receive from previous process
!
   call MPI_Type_Vector(ml*kl, klenbot(irank), nsup-ninf+1, &
        &MPI_REAL, mpitype_rp, ierrs(8))
   call MPI_Type_Commit(mpitype_rp,ierrs(10))
   call MPI_Recv_Init(fl(isb,1,1),1,mpitype_rp,nprevious-1,ktag, &
        & MPI_COMM_WORLD,requests(nrequests+2),ierrs(12))
   nrequests = nrequests+2
end if

if (any(ierrs/=0)) then
   write (iu06,*) ' +++ error: Sub. exchng_init'
   write (iu06,*) ' +++ error: initializing communications, ierrs = ', ierrs
   call abort1
endif

return
end subroutine exchng_init

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
!     local variables :
!     -----------------

integer, save :: requests(4), nrequests
integer       :: istatus(MPI_STATUS_SIZE,4), ierr

! ---------------------------------------------------------------------------- !

extime = MPI_WTIME()-extime

call exchng_init_v()

call MPI_StartAll(nrequests,requests,ierr)
if (ierr/=0) then
   write (iu06,*) ' +++ error: Sub. mpi_exchng_v'
   write (iu06,*) ' +++ error: MPI_Start_All, ierr = ', ierr
   call abort1
endif

!  It might be efficient to put some work here that does not
!  depend on the exchanged data.

call MPI_WaitAll(nrequests,requests,istatus,ierr)
if (ierr/=0) then
   write (iu06,*) ' +++ error: Sub. mpi_exchng_v'
   write (iu06,*) ' +++ error: MPI_Wait_All, ierr = ', ierr
   call abort1
endif

extime = MPI_WTIME()-extime
return

! ---------------------------------------------------------------------------- !

contains

! ---------------------------------------------------------------------------- !

subroutine exchng_init_v

integer :: ist, isb, ierrs(12), ktag

ktag = 400

!     0.3 find index of the start or end point of the different messages
!     ------------------------------------------------------------------

if (nnext/=0) then
   ist = nend(irank)-klenbot(irank+1)+1
endif
isb = nstart(irank)-klenbot(irank)

!     1.0 Initialize communication
!

nrequests = 0
ierrs=0

if (nnext-1>0.and.nnext-1<petotal) then

!      Send to next process

   call MPI_Send_Init (fl(ist), klenbot(nnext), MPI_REAL, nnext-1, ktag,       &
&       MPI_COMM_WORLD, requests(nrequests+1), ierrs(5))

!      Receive from next process

   call MPI_Recv_Init (fl(nend(irank)+1),klentop(irank),MPI_REAL,nnext-1,      &
&       ktag+1, MPI_COMM_WORLD, requests(nrequests+2), ierrs(6))
   nrequests = nrequests+2
end if

if (nprevious-1>=0) then

!     Send to previous process

   call MPI_SEND_Init (fl(nstart(irank)), klentop(nprevious), MPI_REAL,        &
&       nprevious-1, ktag+1, MPI_COMM_WORLD, requests(nrequests+1), ierrs(11))

!     Receive from previous process

   call MPI_Recv_Init (fl(isb), klenbot(irank), MPI_REAL, nprevious-1, ktag,   &
&       MPI_COMM_WORLD, requests(nrequests+2), ierrs(12))
   nrequests = nrequests+2
end if

if (any(ierrs/=0)) then
   write (iu06,*) ' +++ error: Sub. exchng_init_v'
   write (iu06,*) ' +++ error: initializing communications, ierrs = ', ierrs
   call abort1
endif

end subroutine exchng_init_v

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

subroutine mpi_gather_bound (irecv, itag, fl3, fbc)

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
real, dimension (max_nest,kl,ml,n_nest), intent(out) :: fbc !! spectra at 
                                                    !! the selected grid points   

! ---------------------------------------------------------------------------- !
!
!*    local variables :
!     -----------------

integer :: maxlength, ngou, ij, k, m, ip, i
integer :: kcount, ierr, istatus(MPI_STATUS_SIZE)
real, allocatable, dimension (:)  :: zcombuf

! ---------------------------------------------------------------------------- !

comtime = MPI_WTIME()-comtime
maxlength = kl*ml*max_nest*n_nest
allocate (zcombuf(maxlength))

! ---------------------------------------------------------------------------- !
!
!*    1.0 default action if field gathering is not required
!         -------------------------------------------------

if (irecv==0.or.petotal==1) then
   do i = 1,n_nest
      do ngou=1,nbounc(i)
         ij = ijarc(ngou,i)
         do m=1,ml
            do k=1,kl
               fbc(ngou,k,m,i) = fl3(ij,k,m)
            enddo
         enddo
      enddo
   enddo

else if (irank/=irecv) then

!*    1.1 send to the process that gathers the whole field
!         ------------------------------------------------
!
!*    load communication buffer
!
   kcount = 0
   do i = 1,n_nest
      do ngou=1,nbounc(i)
         ij = ijarc(ngou,i)
         if(ij>=nijs.and.ij<=nijl) then
            do m=1,ml
               do k=1,kl
                  kcount = kcount+1
                  zcombuf(kcount) = fl3(ij,k,m)
               enddo
            enddo
         endif
      enddo   
   enddo   
!
!*    send buffer
!
  call MPI_Send(zcombuf, kcount, MPI_REAL, irecv-1, itag,                 &
&                MPI_COMM_WORLD, ierr)

   if (ierr/=0) then
      write (iu06,*) ' +++ error: Sub. mpi_gather_bound'
      write (iu06,*) ' +++ error: MPI_send, ierr = ', ierr
      call abort1
   endif

else

!*    1.2.1  receive contribution to the field from the other processes
!            ----------------------------------------------------------
!
   do ip=1,petotal 
      if (ip==irecv) cycle
   
      call MPI_Recv(zcombuf, maxlength, MPI_REAL, ip-1, itag, &
&                   MPI_COMM_WORLD, istatus, ierr)
      
      if (ierr/=0) then
         write (iu06,*) ' +++ error: Sub. mpi_gather_bound'
         write (iu06,*) ' +++ error: MPI_recv, ierr = ', ierr
         call abort1
      endif
!
!*    decode buffer
!
      kcount = 0
      do i = 1,n_nest
         do ngou=1,nbounc(I)
            ij = ijarc(ngou,i)
            if (ij>=nstart(ip).and.ij<=nend(ip)) then
               do m=1,ml
                  do k=1,kl
                     kcount = kcount+1
                     fbc(ngou,k,m,i) = zcombuf(kcount)
                  enddo
               enddo
            endif
         enddo
      enddo
   enddo

!     1.2.2  contribution from receiving process
!            -----------------------------------

   do i = 1,n_nest
      do ngou=1,nbounc(i)
         ij = ijarc(ngou,i)
         if (ij>=nijs.and.ij<=nijl) then
            do m=1,ml
               do k=1,kl
                  fbc(ngou,k,m,i) = fl3(ij,k,m)
               enddo
            enddo
         endif
      enddo
   enddo
endif
    
deallocate (zcombuf)
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

Logical, save :: frstime = .true.

! ---------------------------------------------------------------------------- !

comtime = MPI_WTIME()-comtime

maxlength = (nscfld+nspfld*kl*ml)*noutp
allocate (zcombuf(maxlength))

! ---------------------------------------------------------------------------- !
!
!     1.0 default action if no field gathering
!         ------------------------------------

If (frstime) then
   if (irecv/=0.and.petotal/=1) then
      If (allocated (noutp_ga)) deallocate (noutp_ga)
      If (allocated (ijar_ga)) deallocate (ijar_ga)
      If (allocated (ngou_ga )) deallocate (ngou_ga)
      allocate (noutp_ga(petotal))
      allocate (ijar_ga (noutp,petotal))
      allocate (ngou_ga (noutp,petotal))
      noutp_ga = 0
      ijar_ga = 0
      ngou_ga = 0
      do ngou = 1,noutp
         ij = ijar(ngou)
         do ip = 1, petotal
            if (ij>=nstart(ip).and.ij<=nend(ip)) then
               noutp_ga(ip) = noutp_ga(ip) + 1
               ijar_ga(noutp_ga(ip),ip) = ij
               ngou_ga(noutp_ga(ip),ip) = ngou
               exit
           end if
         end do
      end do
      IF (sum(noutp_ga(:)).NE.noutp) THEN
         write (iu06,*) ' +++ error: Sub. mpi_gather_spp'
         write (iu06,*) ' +++ decomposion error for output spectra'
         write (iu06,*) ' +++ ', noutp_ga
         write (iu06,*) ' +++ ', sum(noutp_ga(:)), noutp
         call abort1
      end if
   end if   
   frstime = .false.
end if
   
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

subroutine mpi_gather_oifl (isend, field)

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
end subroutine mpi_gather_oifl

    
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

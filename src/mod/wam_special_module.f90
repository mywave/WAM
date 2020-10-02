module wam_special_module

! ---------------------------------------------------------------------------- !
!                                                                              !
!   This module contains a couple of special options and extensions            !
!                                                                              !
!   Arno Behrens      GKSS     August 2010                                     !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!      a.  externals                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

use wam_general_module,   only:  &
&       abort1                      !! terminates processing
    
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!      b.  variables from other modules                                        !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
use wam_file_module,  only: iu06, wpath, area
use wam_mpi_module,   only: irank
 
implicit none
include 'mpif.h'

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!      c.  module variables                                                    !
!                                                                              !
! ---------------------------------------------------------------------------- !

integer :: ispec2d      !! store 2-d spectra up to n hours, 0: standard 
integer :: ispecode     !! 2-d spectra in ascii (1) or in binary code
logical :: readyf       !! wait for ready files ?

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!      d.  public interfaces                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
interface chready                             !! wait for wind input
   module procedure chready
end interface
public chready
  
interface create_ready_file_name              !! build filename
   module procedure create_ready_file_name
end interface
public create_ready_file_name

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

contains

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!      f.  public module procedures                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine chready (ifcst)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   chready - waiting for ready files from atmospherical model to              !
!             start wind input                                                 !
!                                                                              !
!     U. Schaettler      DWD        September 1999                             !
!     A. Behrens         GKSS       August    2010                             !
!                                                                              !
! ---------------------------------------------------------------------------- !
!

character (len=128) :: ytrans_in, yready
character (len= 60) :: yerrmsg
character (len= 12) :: filename
   
logical :: lzexist
integer :: ilen, my_cart_id, ifcst, ierr, nlen
integer :: nmaxwait, nincwait, nwait

! ---------------------------------------------------------------------------- !
!**   If required, check whether ready-files are available                     !
! ---------------------------------------------------------------------------- !

ilen = len_trim (wpath)
ytrans_in = wpath(1:ilen)
my_cart_id = irank-1

nmaxwait = 2400
nincwait = 30

if (ytrans_in /= '   ') then
   if (my_cart_id == 0) then
      nwait   = 0
      lzexist = .false.
      if (area=='GSM'.or.area=='GWM') then
         call create_ready_file_name (ifcst, 'GME_', filename)
      else if (area=='ESM'.or.area=='EWM') then
         call create_ready_file_name (ifcst, 'LMA_', filename)
      else if (area=='KSM'.or.area=='KWM') then
         call create_ready_file_name (ifcst, 'LMA_', filename)
      else
         call create_ready_file_name (ifcst, area//'_', filename)
      endif
      write (iu06,*) ' +++ ready-file ',filename,' has been found !'
      ilen = len_trim (ytrans_in)
      nlen = len_trim (filename)
      yready = ytrans_in(1:ilen)//'/ready/'//filename(1:nlen)

      do while ( (.not. lzexist) .and. (nwait < nmaxwait) )
         inquire (file=yready, exist=lzexist)
         if (.not. lzexist) then

!! ==> wait nincwait seconds and try again  !!
     
            write (iu06,*) ' +++ ready-file not available: ',                  &
&                            yready(1:len_trim(yready))
            write (iu06,*) ' +++ sleep ', nincwait,' seconds'
            call sleep (nincwait)
            nwait = nwait + nincwait
         endif
      enddo
   endif

!! ==> Synchronize the processes again      !!
    
   call mpi_barrier (mpi_comm_world, ierr)

!!  ==> Distribute lzexist to all nodes     !!
 
   call mpi_bcast (lzexist,1,mpi_logical,0,mpi_comm_world, ierr)

   if (.not. lzexist) then
      yerrmsg  = ' *** ERROR:  ready-file not available'
      write (iu06,*) yerrmsg
      call abort1
   endif

endif
end subroutine chready
    
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
 
subroutine create_ready_file_name (ifcst, stbeg, xfile)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   create_ready_file_name                                                     !
!                                                                              !
!   Arno Behrens   GKSS    June 1999 (August 2010 : f90)                       !
!                                                                              !
! ---------------------------------------------------------------------------- !

character (len=12) :: xfile
character (len= 4) :: stbeg
integer :: ifcst, ifcd, ifch
 
xfile(1:4) = stbeg
if (ifcst<24) then
   ifcd = 0
   ifch = ifcst
else
   ifcd = ifcst/24
   ifch = ifcst - (ifcd*24)
endif
write (xfile(5:6),'(i2.2)') ifcd
write (xfile(7:8),'(i2.2)') ifch
xfile(9:12) = '0000'

return
end subroutine create_ready_file_name

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

end module wam_special_module

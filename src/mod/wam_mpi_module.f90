module wam_mpi_module

! ---------------------------------------------------------------------------- !
!                                                                              !
!   this module contains required parameters for the MPI parallelization       !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!   Arno Behrens   MSC/ARMN  October 2003    MPI-parallelization               !
!   Erik Myklebust          November 2004    MPI parallelization               !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!   A. module variables.                                                       !
!                                                                              !
! ---------------------------------------------------------------------------- !

implicit none

! ---------------------------------------------------------------------------- !
!                                                                              !

integer :: pelocal       !!  local process number
integer :: petotal       !!  total number of processors            
integer :: irank         !!  rank of local process (pelocal+1)
integer :: mpmaxlength   !!  maximum length of scalar block message
integer :: nprevious     !!  rank of previous process
integer :: nnext         !!  rank of next process
integer :: ninf          !!  smallest grid point index used by the local process
integer :: nsup          !!  largest grid point index used by the local process
integer :: nijs          !!  index of the first point of the sub grid domain
                         !!  of local process  = nstart(irank)
integer :: nijl          !!  index of the last point of the sub grid domain
                         !!  of local process  = nend(irank)
 
real (kind=kind(1d0)) :: extime      !!  time spent in mpi_exchange
real (kind=kind(1d0)) :: comtime     !!  time spent for other communications

! ---------------------------------------------------------------------------- !
!                                                                              !
!      model decomposition details 
!      ---------------------------

integer, allocatable, dimension (:) :: nstart     !! index of the first point
                                                  !! of the sub grid domain
integer, allocatable, dimension (:) :: nend       !! index of the last point
                                                  !! of the sub grid domain
integer, allocatable, dimension (:) :: nlen       !! number of points in 
                                                  !! the sub grid domain
integer, allocatable, dimension (:) :: klentop    !! length of the message
                                                  !! sent by the next sub
                                                  !! grid domain
integer, allocatable, dimension (:) :: klenbot    !! length of the message
                                                  !! sent by the previous
                                                  !! sub grid domain

! ---------------------------------------------------------------------------- !
!                                                                              !
!         MPI PROCESS NUMBER FOR OUTPUTS.                                      !
!         -------------------------------                                      !
   
integer :: i_out_par      !! process numbers for parameter output
integer :: i_out_spec     !! process numbers for spectra output
integer :: i_out_rad      !! process numbers for radiation stress output
integer :: i_out_b_spec   !! process numbers forboundary spectra output
integer :: i_out_restart  !! process numbers for restart field output

integer, allocatable, dimension (:) :: ipfgtbl !! process numbers for each parameter

! ---------------------------------------------------------------------------- !
!                                                                              !
!         spectra output points on processes.                                  !
!         -----------------------------------                                  !
   
Integer, allocatable, dimension (:)   :: noutp_ga !! no. of points on each pe.
Integer, allocatable, dimension (:,:) :: ijar_ga  !! sea point numbers on each pe.
Integer, allocatable, dimension (:,:) :: ngou_ga  !! output numbers on each pe.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

end module wam_mpi_module

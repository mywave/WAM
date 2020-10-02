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

integer :: petotal       !!  total number of processors
integer :: pelocal       !!  local process number
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
integer :: npoi          !!  number of points of the sub grid domain
                         !!  of local process  = nlen(irank)

real (kind=kind(1d0)) :: extime      !!  time spent in mpi_exchange
real (kind=kind(1d0)) :: comtime=0.  !!  time spent for other communications

! ---------------------------------------------------------------------------- !
!                                                                              !
!      model decomposition details for all processes.                          !
!      ----------------------------------------------                          !

integer, allocatable, dimension (:) :: nstart     !! index of the first point
                                                  !! of the sub grid domain
integer, allocatable, dimension (:) :: nend       !! index of the last point
                                                  !! of the sub grid domain
integer, allocatable, dimension (:) :: nlen       !! number of points in 
                                                  !! the sub grid domain
integer, allocatable, dimension (:) :: klentop    !! number of points
                                                  !! sent by the next sub
                                                  !! grid domains
integer, allocatable, dimension (:) :: klenbot    !! number of points
                                                  !! sent by the previous
                                                  !! sub grid domains

integer, allocatable, dimension (:) :: IJ2NEWIJ   !! TRANSFORMS THE OLD INDEX IJ
                                                  !! INTO THE NEW IJ AS
                                                  !! DETERMINED IN THE CASE OF
                                                  !! 2-D DECOMPOSITION

! ---------------------------------------------------------------------------- !
!                                                                              !
!      model decomposition details for local process.                          !
!      ----------------------------------------------                          !

integer :: NGBTOPE=0    !! total number of pe's to which information will be
                        !! sent from the local pe.
integer :: NTOPEMAX=0   !! total number of points send by local pe.
integer, allocatable, dimension (:) :: NTOPELST   !! list of pe's to which
                                                  !! information will be sent
integer, allocatable, dimension (:) :: NTOPE      !! number of points send
                                                  !! to other processors
integer, allocatable, dimension (:,:) :: IJTOPE   !! sea point indices to be
                                                  !! send by the local pe 
                                                  !! to other pe's.


integer :: NGBFROMPE=0  !! total number of pe's from which information will be
                        !! received by the local pe.
integer :: NFROMPEMAX=0 !! total number  of points received by the local pe.
integer, allocatable, dimension (:) :: NFROMPELST !! list of pe's from which
                                                  !! information will be received
integer, allocatable, dimension (:) :: NFROMPE    !! number of points recieved
                                                  !! from other processors

integer, allocatable, dimension (:) :: NIJSTART   !! INDEX OF THE FIRST HALO
                                                  !! POINT FROM A GIVEN PE
                                                  !! IN THE LOCAL HALO BUFFERS
                                                  !! WHICH PADS BOTH ENDS OF THE
                                                  !! LOCAL GRID POINT 1-D ARRAY

! ---------------------------------------------------------------------------- !
!                                                                              !
!         MPI PROCESS NUMBER FOR OUTPUTS.                                      !
!         -------------------------------                                      !
   
integer :: i_out_par      !! process numbers for parameter output
integer :: i_out_spec     !! process numbers for spectra output
integer :: i_out_b_spec   !! process numbers for boundary spectra output
integer :: i_out_restart  !! process numbers for restart field output

integer, allocatable, dimension (:) :: ipfgtbl !! process numbers for each parameter

! ---------------------------------------------------------------------------- !
!                                                                              !
!         spectra output points on processes.                                  !
!         -----------------------------------                                  !
   
integer, allocatable, dimension (:)   :: noutp_ga !! no. of points on each pe.
integer, allocatable, dimension (:,:) :: ijar_ga  !! sea point numbers on each pe.
integer, allocatable, dimension (:,:) :: ngou_ga  !! output numbers on each pe.

! ---------------------------------------------------------------------------- !
!                                                                              !
!         Boundary output points on processes.                                 !
!         ------------------------------------                                 !

integer, allocatable, dimension (:,:)   :: nbounc_ga !! no. of points on each pe.
integer, allocatable, dimension (:,:,:) :: ijarc_ga  !! sea point numbers on each pe.
integer, allocatable, dimension (:,:,:) :: ngouc_ga  !! output numbers on each pe.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

end module wam_mpi_module

program make_netcdf

! ---------------------------------------------------------------------------- !
!                                                                              !    
!     *make_netcdf* - change output from binary format to NetCDF    ---------- !
!                                                                              !
!     A. Behrens    Helmholtz-Zentrum Geesthacht       March 2011   ---------- !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     PURPOSE.
!     --------
!
!       POSTPROCESSING OF WAM MODEL INTEGRATED DATA.
!
!     INTERFACE.
!     ----------
!
!       *PROGRAM* *PRINT_GRID_FILE*
!          IU01     INPUT UNIT WAVE AND WIND FIELDS (WAMODEL IU20).
!          IU05     USER INPUT.
!          IU06     PRINTER OUTPUT.
!
!     EXTERNALS
!     ---------
!
!       *ABORT1*         - TERMINATES PROCESSING.
!       *INCDATE*        - INCREMENTS DATE-TIME-GROUP
!       *READ_GRID_USER* - READS IN USER INPUT.
!
!     METHOD.
!     -------
!
!       THIS PROGRAM TAKES THE  WAM MODEL OUTPUTS AS INPUT AND
!       PRINTS FIELDS OF INTEGRATED DATA.
!
!     REFERENCE.
!     ----------
!
!        NONE.
!
! ---------------------------------------------------------------------------- !
!
!      EXTERNALS.
!     -----------

use wam_coordinate_module, only : m_sec_to_deg

USE WAM_GENERAL_MODULE, ONLY:  &
&       INCDATE,               &  !! UPDATES A DATE/TIME GROUP.
&       OPEN_FILE,             &  !! OPEN A FILE.
&       reduced_to_regular,    &  !! transformation from a reduced to a regular grid
&       abort1

use wam_netcdf_module

! ---------------------------------------------------------------------------- !
!
!     INTERFACE VARIABLES

USE WAM_PRINT_MODULE, ONLY: ITEST,                                             &
&                           CDATEA, CDATEE, IDELDO,                            &
&                           IU01, FILE01, CDTFILE, IDFILE,                     &
&                           NX, NY, GRID

USE WAM_OUTPUT_SET_UP_MODULE,ONLY: CDTINTT, NOUTT, COUTT,                      &
&                                  NOUT_P, PFLAG_P, CFLAG_P,                   &
&                                  TITL_P, SCAL_P, idelint

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!
!*    LOCAL VARIABLE

LOGICAL            :: IEOF=.false.
character (len=14) :: ihh
character (len=17) :: xfile
integer            :: i, ifail, tstep, ios
LOGICAL,SAVE       :: FRSTIME = .TRUE.
integer            :: iu05 = 55
integer            :: iu06 = 6
character (len=80) :: file05 = 'NETCDF_User'
character (len=80) :: file06 = 'pnetcdf_prot'

integer, allocatable, dimension (:)  :: nlon_rg
integer, allocatable, dimension (:)  :: zdello
real, allocatable, dimension (:)     :: rzdello
real, allocatable, dimension (:,:,:) :: reg_grid
    
character (len=14) :: cstart
integer :: iwest, ieast, isouth, inorth, iredu
real    :: dnx, dny, amowep, amosop, amoeap, amonop, xdella, xdello

! ---------------------------------------------------------------------------- !
!
!*    1. INITALISATION.
!        --------------

namelist /nlnetcdf/ cdatea, cdatee, ideldo, cflag_p, cdtfile, idfile,          &
&                   file01, xdello, xdella, iredu

OPEN (UNIT=IU05, FILE=FILE05, FORM="FORMATTED", STATUS="OLD")
OPEN (UNIT=IU06, FILE=FILE06, FORM="FORMATTED", STATUS="UNKNOWN")
     
!*    1.1 read namelist user input
!         ------------------------

read (iu05,nlnetcdf,iostat=ios)
if (ios==0) then
   write (iu06,*) ' +++ read namelist successfully !'
else
   write (iu06,*) ' +++ read namelist error !'
   call abort1
endif
ideldo  = ideldo*3600
idfile  = idfile*3600
idelint = ideldo
if (iredu==1) then
   write (iu06,*) ' +++ grid is a reduced one !'
else
   write (iu06,*) ' +++ grid is a regular one !'
endif

!*    1.2 FIRST WAVEMODEL OUTPUT FILE DATE AND PRINT DATE.
!         ------------------------------------------------

IF (NOUTT.GT.0) THEN
   CDATEE = '  '                                                             
   CDATEA = COUTT(1)                                                      
   DO I = 1,NOUTT                                                      
      IF (COUTT(I).LT.CDATEA) CDATEA = COUTT(I)                                       
      IF (COUTT(I).GT.CDATEE) CDATEE = COUTT(I)                                       
   END DO                                                              
END IF
CDTINTT = '  '

! ---------------------------------------------------------------------------- !
!
!     2. LOOP OVER OUTPUT FILES.
!        -----------------------

tstep = 0
FILES: DO

!     2.1 FETCH FILE.
!         -----------

   CALL OPEN_FILE (IU06, IU01, FILE01, CDTFILE, 'OLD', IFAIL)
   IF (IFAIL.NE.0) STOP

!     2.2  LOOP OVER OUTPUT TIMES.
!          -----------------------
    
    TIMES: DO   
    
!     2.2.1 READ IN WIND AND WAVE FIELDS.
!           -----------------------------

      read (iu01,err=999,end=999,iostat=ios) cdtintt, dnx, dny,              &
&                                  iwest, isouth, ieast, inorth, cstart    
      amowep = m_sec_to_deg (iwest)
      amoeap = m_sec_to_deg (ieast)
      amosop = m_sec_to_deg (isouth)
      amonop = m_sec_to_deg (inorth)
      nx = nint (dnx)
      ny = nint (dny)
      if (.not.allocated(nlon_rg))  allocate(nlon_rg(1:ny))
      if (.not.allocated(zdello ))  allocate(zdello(1:ny))
      if (.not.allocated(rzdello))  allocate(rzdello(1:ny))

      read (iu01,err=999,end=999) nlon_rg, zdello
      if (iredu==1) then
         rzdello = m_sec_to_deg (zdello)
      endif
      read (iu01,err=999,end=999) pflag_p(1:nout_p)
      if (.not.allocated(grid)) allocate(grid(nx,ny,nout_p))
      if (.not.allocated(reg_grid)) allocate(reg_grid(nx,ny,nout_p))

      do i=1,nout_p
         if (pflag_p(i)) then
            read (iu01,err=999,end=999) grid(:,:,i)
            if (iredu==1) then
               call reduced_to_regular (grid(:,:,i), reg_grid(:,:,i),        &
&                                       nlon_rg, rzdello, xdello)
            endif
         endif
      enddo

!     2.2.2 OUTPUT TIME FOUND?
!           ------------------

      IF (CDTINTT.LT.CDATEA) CYCLE TIMES
      DO WHILE (CDTINTT.GT.CDATEA)
         CALL NEXT_OUTPUT_TIME
         IF (CDATEA.GT.CDATEE) EXIT FILES
         IF (CDTINTT.LT.CDATEA) CYCLE TIMES
      END DO

!     2.2.3 DO OUTPUT OF REQUESTED FIELDS.
!           ------------------------------

      WRITE (IU06,*) ' '
      IF (FRSTIME) THEN
 
!*       open NetCDF file
!        ----------------

         xfile( 1: 4) = 'WAVE'
         xfile( 5:14) = cdatea(1:10)
         xfile(15:17) = '.nc'
         write (iu06,*) ' +++'
         write (iu06,*) ' +++ NetCDF file has been opened - name is: ', trim(xfile)
         write (iu06,*) ' +++'
         call wknco (xfile, cdatea, cflag_p, amowep, amosop, xdella, xdello)
         write (iu06,*) ' +++ idelint, nx, ny, amosop, amowep, xdella, xdello : ',  &
&                             idelint, nx, ny, amosop, amowep, xdella, xdello
         DO I = 1,NOUT_P
           IF (.NOT.PFLAG_P(I) .AND. CFLAG_P(I)) THEN
              WRITE(IU06,*) TITL_P(I), 'IS NOT STORED IN FILE'
            END IF
         END DO     
         FRSTIME = .FALSE.
         WRITE (IU06,*) ' ' 
      END IF
 
      do i = 1,nout_p
         if (pflag_p(i).and.cflag_p(i)) then
            if (iredu==1) then
               call wkncw (i,reg_grid(:,:,i),cdatea,tstep)
            else
               call wkncw (i,grid(:,:,i),cdatea,tstep)
            endif
         endif
      enddo

!*    2.2.7 NEXT OUTPUT TIME.
!           -----------------

      CALL NEXT_OUTPUT_TIME
      tstep = tstep + 1
      IF (CDATEA.GT.CDATEE) EXIT FILES
      if (cdatea>cdtfile) exit times
   END DO TIMES
   
   CALL INCDATE (CDTFILE, IDFILE)       !! INCREMENT DATE FOR THE NEXT FILE.
   if (cdtfile>cdatee) exit files
   CLOSE (UNIT=IU01, STATUS='KEEP')     !! CLOSE OLD FILE
END DO FILES
 
call wkncc                              !! close NetCDF file
999 stop

! ---------------------------------------------------------------------------- !

CONTAINS 

SUBROUTINE NEXT_OUTPUT_TIME
IF (NOUTT.EQ.0) THEN
   CALL INCDATE (CDATEA,IDELDO)
ELSE
   IHH = '99999999999999'
   DO I=1,NOUTT
      IF (COUTT(I).GT.CDATEA .AND. COUTT(I).LT.IHH) IHH = COUTT(I)
   END DO
   CDATEA = IHH
ENDIF
RETURN
END SUBROUTINE NEXT_OUTPUT_TIME

end program make_netcdf

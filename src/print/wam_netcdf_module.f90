!
!########################################################
!#                                                      #
!# NetCDF Outputmodule                                  #
!#                                                      #
!#      Wolfgang Koch  GKSS    July 2009                #
!#      Arno Behrens   GKSS    October 2010             #
!#      Arno Behrens   HZG     June 2014                #
!#                             additional parameters    #
!#                                                      #
!########################################################
!
module wam_netcdf_module

use wam_file_module,          only: iu06
use wam_print_module,         only: nx, ny, amosop, amowep, xdello, xdella
use wam_output_set_up_module, only: idelint

implicit none
private
integer, parameter :: nf=70             !! maximum number of fields !! WAM-MAX
integer, save :: nc = -1
logical, save :: fc = .true.            !! erster Aufruf
integer, save :: no = -1                !! Offset beim schreiben
integer, save :: ncid(0:nf+3) = -1
character (len=14), save :: fsd

real*8, allocatable, dimension (:) :: lat
real*8, allocatable, dimension (:) :: lon

public wknco, wkncw, wkncc, pf, nx, ny, idelint

contains

!##################################################
!# open NetCDF file                               #
!#                                                #
!#  name   filename                               #
!#  n      number of longitudes                   #
!#  m      number of latitudes                    #
!#  l      maximum number of depth layers         #
!#  lon    longitudes                             #
!#  lat    latitudes                              #
!#  ncid   NetCDF file ID-number                  #
!#  sd     date [yyyymmddhh]                      #
!#  dt     time step [s]                          #
!#  flg    parameter switch                       #
!##################################################

subroutine wknco (name, sd, cflag_p, amowep, amosop, xdella, xdello)
use netcdf

implicit none
character (len=17), intent(in) :: name
character (len=14), intent(in) :: sd
logical, intent(in) :: cflag_p(:)
real :: amowep, amosop, xdella, xdello

character (len=33) :: tua
character (len=21) :: tda
character (len=100), dimension (4,nf) :: vl

integer, dimension (3)   :: diid
integer, dimension (0:2) :: did
integer, dimension (nf)  :: flg, ty
integer :: n, m, l, i, loop

real*8 :: dtor, dt
real*4, dimension (nf) :: fw, vmin, vmax

allocate (lat(ny), lon(nx))
dt = real(idelint)
lat(1) = amosop
do loop=2,ny
   lat(loop) = lat(loop-1)+xdella    !! latitudes
enddo
lon(1) = amowep
do loop=2,nx
   lon(loop) = lon(loop-1)+xdello    !! longitudes
enddo

IF (ncid(0)<0) THEN
   fsd = sd
   vl  = ''                          !! integrated parameters
   vl(1, 1) = 'ff'
   vl(2, 1) = 'wind_speed'
   vl(3, 1) = 'Wind speed'
   vl(4, 1) = 'm s-1'
   vmin(1)  = 0.
   vmax(1)  = 40.

   vl(1, 2) = 'dd'
   vl(2, 2) = 'wind_to_direction'
   vl(3, 2) = 'Wind direction'
   vl(4, 2) = 'degree'
   vmin(2)  = 0.
   vmax(2)  = 360.

   vl(1, 3) = 'FV'
   vl(2, 3) = 'FV'
   vl(3, 3) = 'friction velocity'
   vl(4, 3) = 'm s-1'
   vmin(3)  = 0.
   vmax(3)  = 3.

   vl(1, 4) = 'DC'
   vl(2, 4) = 'DC'
   vl(3, 4) = 'drag coefficient'
   vl(4, 4) = ''
   vmin(4)  = 0.
   vmax(4)  = 0.01

   vl(1, 5) = 'CP'
   vl(2, 5) = 'Charnock_parameter'
   vl(3, 5) = 'Charnock parameter'
   vl(4, 5) = ''
   vmin(5)  = 0.
   vmax(5)  = 0.1

   vl(1, 9) = 'hs'
   vl(2, 9) = 'sea_surface_wave_significant_height'
   vl(3, 9) = 'Total significant wave height'
   vl(4, 9) = 'm'
   vmin(9)  = 0.
   vmax(9)  = 20.

   vl(1,10) = 'tp'
   vl(2,10) = 'sea_surface_wave_peak_period_from_variance_spectral_density'
   vl(3,10) = 'Total peak period'
   vl(4,10) = 's'
   vmin(10) = 1.
   vmax(10) = 30.

   vl(1,11) = 'tmp'
   vl(2,11) = 'sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment'
   vl(3,11) = 'Total mean period'
   vl(4,11) = 's'
   vmin(11) = 1.
   vmax(11) = 20.

   vl(1,12) = 'tm1'
   vl(2,12) = 'sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment'
   vl(3,12) = 'Total m1-period'
   vl(4,12) = 's'
   vmin(12) = 1.
   vmax(12) = 20.

   vl(1,13) = 'tm2'
   vl(2,13) = 'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment'
   vl(3,13) = 'Total m2-period'
   vl(4,13) = 's'
   vmin(13) = 1.
   vmax(13) = 20.

   vl(1,14) = 'thq'
   vl(2,14) = 'sea_surface_wave_to_direction'
   vl(3,14) = 'Total mean wave direction'
   vl(4,14) = 'degree'
   vmin(14) = 0.
   vmax(14) = 360.

   vl(1,15) = 'ds'
   vl(2,15) = 'sea_surface_wave_directional_spread'
   vl(3,15) = 'Total directional spreed'
   vl(4,15) = 'degree'
   vmin(15) = 0.
   vmax(15) = 120.

   vl(1,16) = 'NWS'
   vl(2,16) = 'normalised_wave_stress'
   vl(3,16) = 'normalised_wave_stress'
   vl(4,16) = ''
   vmin(16) = 0.
   vmax(16) = 2.

   vl(1,17) = 'hs_sea'
   vl(2,17) = 'sea_surface_wind_wave_significant_height'
   vl(3,17) = 'Sea significant wave height'
   vl(4,17) = 'm'
   vmin(17) = 0.
   vmax(17) = 20.

   vl(1,18) = 'tp_sea'
   vl(2,18) = 'sea_surface_wind_wave_peak_period_from_variance_spectral_density'
   vl(3,18) = 'Sea peak period'
   vl(4,18) = 's'
   vmin(18) = 1.
   vmax(18) = 20.

   vl(1,19) = 'tmp_sea'
   vl(2,19) = 'sea_surface_wind_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment'
   vl(3,19) = 'Sea mean period'
   vl(4,19) = 's'
   vmin(19) = 1.
   vmax(19) = 20.

   vl(1,20) = 'tm1_sea'
   vl(2,20) = 'sea_surface_wind_wave_mean_period_from_variance_spectral_density_first_frequency_moment'
   vl(3,20) = 'Sea m1-period'
   vl(4,20) = 's'
   vmin(20) = 1.
   vmax(20) = 20.

   vl(1,21) = 'tm2_sea'
   vl(2,21) = 'sea_surface_wind_wave_mean_period_from_variance_spectral_density_second_frequency_moment'
   vl(3,21) = 'Sea m2-period'
   vl(4,21) = 's'
   vmin(21) = 1.
   vmax(21) = 20.

   vl(1,22) = 'thq_sea'
   vl(2,22) = 'sea_surface_wind_wave_to_direction'
   vl(3,22) = 'Sea mean wave direction'
   vl(4,22) = 'degree'
   vmin(22) = 0.
   vmax(22) = 360.

   vl(1,23) = 'ds_sea'
   vl(2,23) = 'sea_surface_wind_wave_directional_spread'
   vl(3,23) = 'Sea directional spreed'
   vl(4,23) = 'degree'
   vmin(23) = 0.
   vmax(23) = 120.

   vl(1,25) = 'hs_swell'
   vl(2,25) = 'sea_surface_swell_wave_significant_height'
   vl(3,25) = 'Swell significant wave height'
   vl(4,25) = 'm'
   vmin(25) = 0.
   vmax(25) = 20.

   vl(1,26) = 'tp_swell'
   vl(2,26) = 'sea_surface_swell_wave_peak_period_from_variance_spectral_density'
   vl(3,26) = 'Swell peak period'
   vl(4,26) = 's'
   vmin(26) = 1.
   vmax(26) = 30.

   vl(1,27) = 'tmp_swell'
   vl(2,27) = 'sea_surface_swell_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment'
   vl(3,27) = 'Swell mean period'
   vl(4,27) = 's'
   vmin(27) = 1.
   vmax(27) = 25.

   vl(1,28) = 'tm1_swell'
   vl(2,28) = 'sea_surface_swell_wave_mean_period_from_variance_spectral_density_first_frequency_moment'
   vl(3,28) = 'Swell m1-period'
   vl(4,28) = 's'
   vmin(28) = 1.
   vmax(28) = 25.

   vl(1,29) = 'tm2_swell'
   vl(2,29) = 'sea_surface_swell_wave_mean_period_from_variance_spectral_density_second_frequency_moment'
   vl(3,29) = 'Swell tm2-period'
   vl(4,29) = 's'
   vmin(29) = 1.
   vmax(29) = 25.

   vl(1,30) = 'thw_swell'
   vl(2,30) = 'sea_surface_swell_wave_to_direction'
   vl(3,30) = 'Swell mean wave direction'
   vl(4,30) = 'degree'
   vmin(30) = 0.
   vmax(30) = 360.

   vl(1,31) = 'ds_swell'
   vl(2,31) = 'sea_surface_swell_wave_directional_spread'
   vl(3,31) = 'Swell directional spreed'
   vl(4,31) = 'degree'
   vmin(31) = 0.
   vmax(31) = 120.

   vl(1,32) = 'Z0'
   vl(2,32) = 'roughness length'
   vl(3,32) = 'roughness length'
   vl(4,32) = 'm'
   vmin(32) = 0.
   vmax(32) = 100.

   vl(1,33) = 'goda'
   vl(2,33) = 'goda_peakness_parameter'
   vl(3,33) = 'goda peakness parameter'
   vl(4,33) = ''
   vmin(33) = 0.
   vmax(33) = 15.

   vl(1,34) = 'kurtosis'
   vl(2,34) = 'kurtosis'
   vl(3,34) = 'kurtosis'
   vl(4,34) = ''
   vmin(34) = 0.
   vmax(34) = 1.

   vl(1,35) = 'BFI'
   vl(2,35) = 'Benjamin_Feir_index'
   vl(3,35) = 'Benjamin Feir index'
   vl(4,35) = ''
   vmin(35) = -10.
   vmax(35) = 10.

   vl(1,36) = 'mHs'
   vl(2,36) = 'normalized_maximum_wave_height'
   vl(3,36) = 'normalized maximum wave height'
   vl(4,36) = 'm'
   vmin(36) = 0.
   vmax(36) = 40.

   vl(1,37) = 'mwp'
   vl(2,37) = 'maximum_wave_period'
   vl(3,37) = 'maximum wave period'
   vl(4,37) = 's'
   vmin(37) = 1.
   vmax(37) = 30.

   vl(1,38) = 'TpI'
   vl(2,38) = 'interpolated_peak_frequency'
   vl(3,38) = 'interpolated peak frequency'
   vl(4,38) = 's'
   vmin(38) = 0.
   vmax(38) = 0.8

   vl(1,39) = 'Pdir'
   vl(2,39) = 'peak_direction'
   vl(3,39) = 'peak direction'
   vl(4,39) = 'degree'
   vmin(39) = 0.
   vmax(39) = 360.

   vl(1,40) = 'msqs'
   vl(2,40) = 'mean_square_slope'
   vl(3,40) = 'mean square slope'
   vl(4,40) = ''
   vmin(40) = 0.
   vmax(40) = 1.

   vl(1,41) = 'fshs'
   vl(2,41) = 'first swell significant wave height'
   vl(3,41) = 'first swell significant wave height'
   vl(4,41) = 'm'
   vmin(41) = 0.
   vmax(41) = 20.

   vl(1,42) = 'fstm1'
   vl(2,42) = 'first swell mean period'
   vl(3,42) = 'first swell mean period'
   vl(4,42) = 'm/s'
   vmin(42) = 1.
   vmax(42) = 25.

   vl(1,43) = 'fsdir'
   vl(2,43) = 'first swell direction'
   vl(3,43) = 'first swell direction'
   vl(4,43) = 'degree'
   vmin(43) = 0.
   vmax(43) = 360.

   vl(1,44) = 'sshs'
   vl(2,44) = 'second swell significant wave height'
   vl(3,44) = 'second swell significant wave height'
   vl(4,44) = 'm'
   vmin(44) = 0.
   vmax(44) = 20.

   vl(1,45) = 'sstm1'
   vl(2,45) = 'second swell mean period'
   vl(3,45) = 'second swell mean period'
   vl(4,45) = 'm/s'
   vmin(45) = 1.
   vmax(45) = 25.

   vl(1,46) = 'ssdir'
   vl(2,46) = 'second swell direction'
   vl(3,46) = 'second swell direction'
   vl(4,46) = 'degree'
   vmin(46) = 0.
   vmax(46) = 360.

   vl(1,47) = 'tshs'
   vl(2,47) = 'third swell significant wave height'
   vl(3,47) = 'third swell significant wave height'
   vl(4,47) = 'm'
   vmin(47) = 0.
   vmax(47) = 20.

   vl(1,48) = 'tstm1'
   vl(2,48) = 'third swell mean period'
   vl(3,48) = 'third swell mean period'
   vl(4,48) = 'm/s'
   vmin(48) = 1.
   vmax(48) = 25.

   vl(1,49) = 'tsdir'
   vl(2,49) = 'third swell direction'
   vl(3,49) = 'third swell direction'
   vl(4,49) = 'degree'
   vmin(49) = 0.
   vmax(49) = 360.

   vl(1,51) = 'radsxx'
   vl(2,51) = 'radiation stress tensor sxx'
   vl(3,51) = 'radiation stress tensor sxx'
   vl(4,51) = 'kg/s/s'
   vmin(51) = -1.e07
   vmax(51) = 1.e07

   vl(1,52) = 'radsyy'
   vl(2,52) = 'radiation stress tensor syy'
   vl(3,52) = 'radiation stress tensor syy'
   vl(4,52) = 'kg/s/s'
   vmin(52) = -1.e07
   vmax(52) = 1.e07

   vl(1,53) = 'radsxy'
   vl(2,53) = 'radiation stress tensor sxy'
   vl(3,53) = 'radiation stress tensor sxy'
   vl(4,53) = 'kg/s/s'
   vmin(53) = -1.e07
   vmax(53) = 1.e07

   vl(1,55) = 'wforx'
   vl(2,55) = 'x-comp. wave force per surface unit'
   vl(3,55) = 'x-comp. wave force per surface unit'
   vl(4,55) = 'n/m/m'
   vmin(55) = -1.
   vmax(55) = 1.

   vl(1,56) = 'wfory'
   vl(2,56) = 'y-comp. wave force per surface unit'
   vl(3,56) = 'y-comp. wave force per surface unit'
   vl(4,56) = 'n/m/m'
   vmin(56) = -1.
   vmax(56) = 1.

   vl(1,57) = 'sdx'
   vl(2,57) = 'x-comp. Stokes drift'
   vl(3,57) = 'x-comp. Stokes drift'
   vl(4,57) = 'm/s'
   vmin(57) = -1.
   vmax(57) = 1.

   vl(1,58) = 'sdy'
   vl(2,58) = 'y-comp. Stokes drift'
   vl(3,58) = 'y-comp. Stokes drift'
   vl(4,58) = 'm/s'
   vmin(58) = -1.
   vmax(58) = 1.

   vl(1,59) = 'phioc'
   vl(2,59) = 'normalized energy flux to ocean'
   vl(3,59) = 'normalized energy flux to ocean'
   vl(4,59) = 'kg/s/s'
   vmin(59) = -10.
   vmax(59) = 30.

   vl(1,60) = 'phiaw'
   vl(2,60) = 'normalized energy flux from wind to waves'
   vl(3,60) = 'normalized energy flux from wind to waves'
   vl(4,60) = 'kg/s/s'
   vmin(60) = 0.
   vmax(60) = 30.

   vl(1,61) = 'tauocx'
   vl(2,61) = 'x-comp. normalized momentum flux into ocean'
   vl(3,61) = 'x-comp. normalized momentum flux into ocean'
   vl(4,61) = 'kg/s/s'
   vmin(61) = -5.
   vmax(61) = 5.

   vl(1,62) = 'tauocy'
   vl(2,62) = 'y-comp. normalized momentum flux into ocean'
   vl(3,62) = 'y-comp. normalized momentum flux into ocean'
   vl(4,62) = 'kg/s/s'
   vmin(62) = -5.
   vmax(62) = 5.

   vl(1,63) = 'phibot'
   vl(2,63) = 'energy flux from waves to bottom'
   vl(3,63) = 'energy flux from waves to bottom'
   vl(4,63) = 'kg/s/s'
   vmin(63) = -1.
   vmax(63) = 10.

   vl(1,64) = 'taubot_x'
   vl(2,64) = 'x-comp. momentum flux from waves into bottom'
   vl(3,64) = 'x-comp. momentum flux from waves into bottom'
   vl(4,64) = 'kg/s/s'
   vmin(64) = -0.1
   vmax(64) = 0.1

   vl(1,65) = 'taubot_y'
   vl(2,65) = 'y-comp. momentum flux from waves into bottom'
   vl(3,65) = 'y-comp. momentum flux from waves into bottom'
   vl(4,65) = 'kg/s/s'
   vmin(65) = -0.1
   vmax(65) = 0.1
      
   vl(1,67) = 'cmax_F'                                      !! WAM-MAX
   vl(2,67) = 'exp max crest (Forristall)'                  !! WAM-MAX
   vl(3,67) = 'exp maximum crest height as per Forristall ()'    !! WAM-MAX
   vl(4,67) = 'm'                                           !! WAM-MAX
   vmin(67) = 0.                                            !! WAM-MAX
   vmax(67) = 30.                                           !! WAM-MAX
                                                            !! WAM-MAX
   vl(1,68) = 'Hmax_N'                                      !! WAM-MAX
   vl(2,68) = 'exp max height (Naess)'                      !! WAM-MAX
   vl(3,68) = 'exp maximum wave height as per Naess ()'     !! WAM-MAX
   vl(4,68) = 'm'                                           !! WAM-MAX
   vmin(68) = 0.                                            !! WAM-MAX
   vmax(68) = 50.                                           !! WAM-MAX
                                                            !! WAM-MAX
   vl(1,69) = 'cmax_st'                                     !! WAM-MAX
   vl(2,69) = 'CRMAX (SPACE-TIME (STQD))'                   !! WAM-MAX
   vl(3,69) = 'maximum crest height - space-time (stqd)'    !! WAM-MAX
   vl(4,69) = 'm'                                           !! WAM-MAX
   vmin(69) = 0.                                            !! WAM-MAX
   vmax(69) = 30.                                           !! WAM-MAX
                                                            !! WAM-MAX
   vl(1,70) = 'hmax_st'                                     !! WAM-MAX
   vl(2,70) = 'HMAX (SPACE-TIME (STQD))'                    !! WAM-MAX
   vl(3,70) = 'maximum wave height - space-time (stqd)'     !! WAM-MAX
   vl(4,70) = 'm'                                           !! WAM-MAX
   vmin(70) = 0.                                            !! WAM-MAX
   vmax(70) = 50.                                           !! WAM-MAX
    
   flg = 1                   !! prepare table for required parameters only
   do loop=1,nf
      if (.not.cflag_p(loop)) then
         flg([loop]) = 0
      endif
   enddo

   ty = NF90_FLOAT           !! type of parameter
   fw = -9999999.            !! dummy (zmiss)
   did = [3,2,1]
   i = dt
   WRITE(tda,'("0000-00-00 (",2(i2.2,":"),i2.2,")")')i/3600,MOD(i/60,60),MOD(i,60)
   tua="seconds since "//sd(1:4)//"-"//sd(5:6)//"-"//sd(7:8)//" "//sd(9:10)//":"//sd(11:12)//":"//sd(13:14)

   CALL Pf(NF90_CREATE(name,ior(NF90_CLOBBER,NF90_SHARE),ncid(0)))
   CALL Pf(NF90_DEF_DIM(ncid(0),'time',NF90_UNLIMITED,diid(1)))
   CALL Pf(NF90_DEF_DIM(ncid(0),'lat',ny,diid(2)))
   CALL Pf(NF90_DEF_DIM(ncid(0),'lon',nx,diid(3)))
   CALL Pf(NF90_DEF_VAR(ncid(0),'lat',NF90_FLOAT,[diid(2)],ncid(nf+1)))
   CALL Pf(NF90_DEF_VAR(ncid(0),'lon',NF90_FLOAT,[diid(3)],ncid(nf+2)))
   CALL Pf(NF90_DEF_VAR(ncid(0),'time',NF90_INT,[diid(1)],ncid(nf+3)))
   call pf(nf90_put_att(ncid(0),ncid(nf+1),'units',"degrees_north"))
   call pf(nf90_put_att(ncid(0),ncid(nf+2),'units',"degrees_east"))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'delta_t',tda))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'units',tua))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'dt',i))
   call pf(nf90_put_att(ncid(0),ncid(nf+3),'standard_name',"time"))
   call pf(nf90_put_att(ncid(0),ncid(nf+3),'calendar',"standard"))

   DO i=1,nf
      IF (flg(i)>0) THEN
         CALL Pf(NF90_DEF_VAR(ncid(0),vl(1,i),ty(i),diid(did),ncid(i)))
         CALL Pf(NF90_PUT_ATT(ncid(0),ncid(i),'_FillValue',fw(i)))
         call pf(nf90_put_att(ncid(0),ncid(i),'long_name',vl(3,i)))
         call pf(nf90_put_att(ncid(0),ncid(i),'standard_name',vl(2,i)))
         call pf(nf90_put_att(ncid(0),ncid(i),'units',vl(4,i)))
         call pf(nf90_put_att(ncid(0),ncid(i),'valid_min',vmin(i)))
         call pf(nf90_put_att(ncid(0),ncid(i),'valid_max',vmax(i)))
      ENDIF
   ENDDO
   CALL Pf(NF90_ENDDEF(ncid(0)))
   CALL Pf(NF90_PUT_VAR(ncid(0),ncid(nf+1),lat))
   CALL Pf(NF90_PUT_VAR(ncid(0),ncid(nf+2),lon))
ELSE
   write (iu06,*) ' +++ NetCDF-file is open already'
   STOP           ' +++ Problem occurs during NetCDF-Output handling '
ENDIF
END subroutine wknco

!########################################################
!#      write a field to NetCDF file                    #
!#                                                      #
!#      nid     pointer to NetCDF file ID-number        #
!#      vid     pointer to NetCDF parameter ID-number   #
!########################################################

subroutine wkncw (ip,grid,ad,t)
use netcdf

implicit none
real,    intent(in)    :: grid(:,:)   !! grid of parameter
integer, intent(in)    :: ip          !! parameter number
integer, dimension (3) :: sta

character (len=14), intent(in) :: ad
integer :: vid, j, t, tda

if (fc) then
   no = 1-t                           !! t+no=1 : first output
   fc = .false.
endif
tda = idelint
vid = ip
if (ncid(vid)>0) then
   sta = [1,1,t+no]
   CALL Pf(NF90_PUT_VAR(ncid(0),ncid(vid),grid,sta))
   CALL Pf(NF90_PUT_VAR(ncid(0),ncid(nf+3),t*tda,[t+no]))
else
   write (iu06,*) ' +++ parameter ',vid,' is not available in NetCDF-file'
endif
end subroutine wkncw

!##############################
!#      close NetCDF-file     #
!##############################

subroutine wkncc
use netcdf
if (ncid(0)>=0) CALL Pf (NF90_CLOSE(ncid(0)))
ncid = -1
end subroutine wkncc

!################################################
!#   pf      write a NetCDF-error message       #
!#   en      NetCDF error number                #
!################################################

subroutine pf(en)
use netcdf
integer :: en
if (en/=0) write (iu06,*) ' +++ Error : ', NF90_STRERROR(en)
end subroutine pf
end module wam_netcdf_module

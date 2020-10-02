MODULE WAM_USER_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL DEFAULT WAM MODEL SETTINGS, WHICH CAN BE          !
!   CONTROLLED OR WHICH MUST BE DEFINED BY THE USER IN DIFFERNT VERSIONS OF    !
!   USER INPUT.                                                                !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE             !! COORDINATE PROCEDURES

USE WAM_OUTPUT_SET_UP_MODULE, ONLY: &
&       SET_OUTPUT_FILE_STEP,       & !! SETS OUTPUT FILES TIMESTEP.
&       SET_PARAMETER_OUTPUT_FLAGS, & !! SET PARAMETER OUTPUT FLAGS.
&       SET_SPECTRA_OUTPUT_FLAGS,   & !! SET SPECTRA OUTPUT FLAGS.
&       SET_OUTPUT_SITES,           & !! SET OUTPUT SITES FOR SPECTRA.
&       SET_OUTPUT_TIMES,           & !! SET OUTPUT TIMES.
&       set_spectral_code,          & !! set spectral code and output time period
&       set_wammax_options,         & !! set wam-max time/space-time options       !! WAM-MAX
&       set_ready_outfile_flag,     & !! ready files for integrated parameters ?
&       set_ready_outfile_directory   !! set full path name of output ready file
                                      !! diectory (integrated parameters)

USE WAM_GENERAL_MODULE, ONLY:       & !! Sets Betamax
&       SET_GENERAL_MODULE

USE WAM_WIND_MODULE,      ONLY:     &
&       SET_WIND_TIMESTEPS,         & !! SET WIND TIMESTEPS.
&       set_ready_file_flag,        & !! wait for ready files ?
&       set_ready_file_directory      !! set full pathname of ready file directory

USE WAM_TOPO_MODULE,      ONLY:     &
&       SET_TOPO_TIMESTEPS            !! SET TOPO TIMESTEPS.

USE WAM_CURRENT_MODULE,   ONLY:     &
&       SET_CURRENT_TIMESTEPS         !! SET CURRENT TIMESTEPS.

USE WAM_ICE_MODULE,       ONLY:     &
&       SET_ICE_TIMESTEP              !! SET ICE TIMESTEPS.

USE WAM_NEST_MODULE,  ONLY:         &
&       SET_BOUNDARY_OPTION           !! SET THE BOUNDARY OPTIONS.

USE WAM_BOUNDARY_MODULE,  ONLY:     &
&       SET_BOUNDARY_OUTPUT_TIMESTEPS

USE WAM_TIMOPT_MODULE,    ONLY:     &
&       SET_INTEGRATION_PERIOD,     & !! SET INTEGRATION PERIOD. 
&       SET_INTEGRATION_TIMESTEPS,  & !! SET INTEGRATION TIMESTEPS.
&       SET_MODEL_OPTION,           & !! SET MODEL OPTIONS.
&       SET_START_OPTION              !! SET START OPTION.

USE WAM_COLDSTART_MODULE,    ONLY:  &
&       SET_C_START_PAR               !! SETS PARAMETER FOR INITIAL SPECTRA.

USE WAM_FILE_MODULE,      ONLY:     &
&       SET_TEST_OPTION,            & !! SETS TEST OPTION.
&       SET_WIND_FILE,              & !! WIND DATA FILE (FORM. INPUT).
&       SET_B_INPUT_FILE,           & !! BOUNDARY VALUE INPUT FILE IDENTIFIER.
&       SET_ICE_FILE,               & !! ICE DATA FILE (FORM. INPUT).
&       SET_TOPO_FILE,              & !! DEPTH DATA FILE (FORM. INPUT).
&       SET_CURRENT_FILE,           & !! CURRENT DATA FILE (FORM. INPUT).
&       SET_PREPROC_FILE,           & !! GRID DATA FILE (UNFORM. INPUT).
&       SET_RESTART_FILE,           & !! RESTART FILE (UNFORM. INPUT/OUTPUT).
&       SET_B_OUTPUT_FILE,          & !! BOUNDARY DATA FILE (UNFORM. OUTPUT).
&       SET_MAP_FILE,               & !! INTEGRATED DATA FILE (UNFORM. OUTPUT).
&       SET_SPECTRA_FILE              !! SPECTRA DATA FILE (UNFORM. OUTPUT).

USE WAM_RESTART_MODULE,   ONLY:     &
&       SET_RESTART_FILE_STEP         !! SETS RESTART FILE TIMESTEP.

USE WAM_ASSI_SET_UP_MODULE, ONLY:   &
&       SET_ASSIMILATION_OPTION,    & !! SETS ASSIMILATION OPTIONS.
&       SET_ASSIMILATION_OUTPUT,    & !! SETS ASSIMILATION OUTPUT.
&       SET_ASSIMILATION_PERIOD,    & !! SETS ASSIMILATION PERIOD AND TIMESTEP.
&       SET_OBSERVATION_FILE,       & !! SETS OBSERVATION FILE IDENTIFIER.
&       SET_ASSI_MAP_FILE,          & !! INTEGRATED DATA FILE (UNFORM. OUTPUT)
&       SET_ASSI_SPECTRA_FILE         !! SPECTRA DATA FILE (UNFORM. OUTPUT)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,    ONLY: IU05, FILE05, IU06
   
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

! ---------------------------------------------------------------------------- !

CHARACTER (LEN=14) :: START_DATE = ' '  !! START DATE OF RUN.
CHARACTER (LEN=14) :: END_DATE   = ' '  !! END DATE OF RUN.

! ---------------------------------------------------------------------------- !

LOGICAL :: COLDSTART     !! TRUE = COLDSTART, FALSE = HOT START.

INTEGER :: IOPTI         !! COLD START OPTION:
                         !! =  0 WIND INDEPENDENT INITIAL VALUES.
                         !! =  1 WIND DEPENDENT INITIAL VALUES AND
                         !!      ENERGY EQUAL ZERO IF WINDSPEED IS ZERO
                         !! =  2 WIND DEPENDENT INITIAL VALUES AND
                         !!      ENERGY COMPUTED FROM GIVEN PARAMETERS IF
                         !!      WINDSPEED IS ZERO.

REAL :: ALPHA     !! PHILLIPS' PARAMETER  (NOT USED IF IOPTI = 1)
REAL :: FM        !! PEAK FREQUENCY (HZ) AND/OR MAXIMUM FREQUENCY
REAL :: GAMMA     !! OVERSHOOT FACTOR
REAL :: SIGMA_A   !! LEFT PEAK WIDTH
REAL :: SIGMA_B   !! RIGHT PEAK WIDTH
REAL :: THETAQ    !! WAVE DIRECTION (DEG) (NOT USED IF IOPTI = 1)
REAL :: FETCH     !! FETCH IN METRES (IF ZERO THEN 0.5 OF THE
                  ! LATITUDE INCREMENT IS USED.)

! ---------------------------------------------------------------------------- !

LOGICAL :: SPHERICAL_RUN       !! TRUE: SPHERICAL  PROPAGATION,
                               !! FALSE: CARTESIAN PROPAGATION.
LOGICAL :: SHALLOW_RUN         !! TRUE:  SHALLOW WATER MODEL,
                               !! FALSE:  DEEP WATER MODEL.
LOGICAL :: REFRACTION_D_RUN    !! TRUE: DEPTH REFRACTION ON.
LOGICAL :: REFRACTION_C_RUN    !! TRUE: CURRENT REFRACTION ON.
LOGICAL :: L_OBSTRUCTION       !! FALSE: NO REDUCTION DUE TO SUB-GRID FEATURES.
LOGICAL :: L_DECOMP            !! TRUE: IF MPI RUN, THEN 1D DECOMPOSITION OF GRID.
INTEGER :: IPHYS               !! PHYSICS PARAMETERISATION FOR INPUT AND
                               !! OPEN OCEAN DISSIPATION
                               !!  0 : ECMWF CY45R1
                               !!  1 : ECMWF CY46R1, based on Ardhuin et al. 2010
LOGICAL :: WAVE_BREAKING_RUN   !! TRUE: WAVE BREAKING ON.
LOGICAL :: PHILLIPS_RUN        !! TRUE:  PHILLIPS SOURCE ON.
INTEGER :: ISNONLIN            !! = 0 OLD DEPTH SCALING
                               !! ≠ 0 NEW DEPTH SCALING
INTEGER :: ITEST               !! TEST OUTPUT UP TO LEVEL.

REAL    :: BETAMAX             !! PARAMETER FOR WIND INPUT

! ---------------------------------------------------------------------------- !

INTEGER            :: PROPAGATION_TIMESTEP
CHARACTER (LEN=1)  :: PROPAGATION_TIMESTEP_UNIT
INTEGER            :: SOURCE_TIMESTEP
CHARACTER (LEN=1)  :: SOURCE_TIMESTEP_UNIT

! ---------------------------------------------------------------------------- !

INTEGER            :: RESTART_SAVE_TIMESTEP
CHARACTER (LEN=1)  :: RESTART_SAVE_TIMESTEP_UNIT
INTEGER            :: RESTART_FILE_UNIT
CHARACTER (LEN=80) :: RESTART_FILE_NAME

! ---------------------------------------------------------------------------- !

LOGICAL            :: COARSE_GRID_RUN             !! TRUE: COARSE GRID MODEL.
INTEGER            :: COARSE_OUTPUT_TIMESTEP
CHARACTER (LEN=1)  :: COARSE_OUTPUT_TIMESTEP_UNIT
INTEGER            :: COARSE_FILE_SAVE_TIMESTEP
CHARACTER (LEN=1)  :: COARSE_FILE_SAVE_TIMESTEP_UNIT
INTEGER            :: COARSE_OUTPUT_FILE_UNIT
CHARACTER (LEN=80) :: COARSE_OUTPUT_FILE_NAME

LOGICAL            :: FINE_GRID_RUN               !! TRUE: FINE GRID MODEL.
INTEGER            :: FINE_INPUT_FILE_UNIT
CHARACTER (LEN=80) :: FINE_INPUT_FILE_NAME

! ---------------------------------------------------------------------------- !

INTEGER            :: WIND_INPUT_TIMESTEP
CHARACTER (LEN=1)  :: WIND_INPUT_TIMESTEP_UNIT
INTEGER            :: WIND_OUTPUT_TIMESTEP
CHARACTER (LEN=1)  :: WIND_OUTPUT_TIMESTEP_UNIT
INTEGER            :: WIND_INPUT_FILE_UNIT
CHARACTER (LEN=80) :: WIND_INPUT_FILE_NAME

! ---------------------------------------------------------------------------- !

INTEGER            :: TOPO_INPUT_TIMESTEP
CHARACTER (LEN=1)  :: TOPO_INPUT_TIMESTEP_UNIT
INTEGER            :: TOPO_OUTPUT_TIMESTEP
CHARACTER (LEN=1)  :: TOPO_OUTPUT_TIMESTEP_UNIT
INTEGER            :: TOPO_INPUT_FILE_UNIT
CHARACTER (LEN=80) :: TOPO_INPUT_FILE_NAME

! ---------------------------------------------------------------------------- !

INTEGER            :: CURRENT_INPUT_TIMESTEP
CHARACTER (LEN=1)  :: CURRENT_INPUT_TIMESTEP_UNIT
INTEGER            :: CURRENT_OUTPUT_TIMESTEP
CHARACTER (LEN=1)  :: CURRENT_OUTPUT_TIMESTEP_UNIT
INTEGER            :: CURRENT_INPUT_FILE_UNIT
CHARACTER (LEN=80) :: CURRENT_INPUT_FILE_NAME

! ---------------------------------------------------------------------------- !

INTEGER            :: ICE_INPUT_TIMESTEP
CHARACTER (LEN=1)  :: ICE_INPUT_TIMESTEP_UNIT
INTEGER            :: ICE_INPUT_FILE_UNIT
CHARACTER (LEN=80) :: ICE_INPUT_FILE_NAME

! ---------------------------------------------------------------------------- !

INTEGER            :: PARAMETER_OUTPUT_TIMESTEP
CHARACTER (LEN=1)  :: PARAMETER_OUTPUT_TIMESTEP_UNIT
INTEGER            :: PARAMETER_OUTPUT_FILE_UNIT
CHARACTER (LEN=80) :: PARAMETER_OUTPUT_FILE_NAME

INTEGER            :: SPECTRA_OUTPUT_TIMESTEP
CHARACTER (LEN=1)  :: SPECTRA_OUTPUT_TIMESTEP_UNIT
INTEGER            :: SPECTRA_OUTPUT_FILE_UNIT
CHARACTER (LEN=80) :: SPECTRA_OUTPUT_FILE_NAME

INTEGER            :: OUTPUT_FILE_SAVE_TIMESTEP
CHARACTER (LEN=1)  :: OUTPUT_FILE_SAVE_TIMESTEP_UNIT

! ---------------------------------------------------------------------------- !

INTEGER, PARAMETER :: MOUTT = 100
CHARACTER (LEN=14), DIMENSION(MOUTT) :: COUTT  !! SPECIFIED OUTPUT TIMES.

! ---------------------------------------------------------------------------- !

INTEGER, PARAMETER         :: NOUT_P = 70
LOGICAL, DIMENSION(NOUT_P) :: FFLAG_P         !! FILE PARAMETER OUTPUT FLAG.
LOGICAL, DIMENSION(NOUT_P) :: PFLAG_P         !! PRINTER PARAMETER OUTPUT FLAG.
logical :: orientation_of_directions          !! coming from or going to ?

! ---------------------------------------------------------------------------- !

INTEGER, PARAMETER         :: NOUT_S = 4
LOGICAL, DIMENSION(NOUT_S) :: FFLAG_S         !! FILE SPECTRA OUTPUT FLAG. 
LOGICAL, DIMENSION(NOUT_S) :: PFLAG_S         !! PRINTER SPECTRA  OUTPUT FLAG.

INTEGER, PARAMETER         :: MOUTP   = 200
CHARACTER(LEN=LEN_COOR), DIMENSION(MOUTP) :: OUTLAT
CHARACTER(LEN=LEN_COOR), DIMENSION(MOUTP) :: OUTLONG
    
CHARACTER (LEN=20), DIMENSION(MOUTP) :: NAME

INTEGER, PARAMETER         :: NOUT_R = 12
LOGICAL, DIMENSION(NOUT_R) :: FFLAG_R         !! FILE OUTPUT FLAG.
LOGICAL, DIMENSION(NOUT_R) :: PFLAG_R         !! PRINTER OUTPUT FLAG.

! ---------------------------------------------------------------------------- !

INTEGER            :: PREPROC_OUTPUT_FILE_UNIT
CHARACTER (LEN=80) :: PREPROC_OUTPUT_FILE_NAME

! ---------------------------------------------------------------------------- !

logical             :: ready_file_flag, ready_outfile_flag
character (len=128) :: ready_file_directory, ready_outfile_directory
character (len=  3) :: model_area
integer             :: spectral_code
integer             :: hours_2d_spectra

! ---------------------------------------------------------------------------- !

integer             :: assimilation_flag
character (len= 14) :: assimilation_start_date
character (len= 14) :: assimilation_end_date
integer             :: assimilation_time_step
CHARACTER (LEN=1)   :: assimilation_time_step_UNIT



INTEGER             :: observation_file_unit
character (len= 80) :: observation_filename
logical             :: first_guess_output_flag
INTEGER             :: first_guess_ip_file_unit
character (len= 80) :: first_guess_ip_filename
INTEGER             :: first_guess_sp_file_unit
character (len= 80) :: first_guess_sp_filename

real                :: influence_radius
real                :: observation_scatter
real                :: model_scatter

real                :: wammax_dur                                          !! WAM-MAX
character (len=  5) :: wammax_dur_unit                                     !! WAM-MAX
real                :: wammax_dx                                           !! WAM-MAX
character (len=  5) :: wammax_dx_unit                                      !! WAM-MAX
real                :: wammax_dy                                           !! WAM-MAX
character (len=  5) :: wammax_dy_unit                                      !! WAM-MAX

! ---------------------------------------------------------------------------- !

NAMELIST /WAM_NAMELIST/                                                        &
&       START_DATE,                 END_DATE,                                  &
&       COLDSTART,                                                             &
&       IOPTI, ALPHA, FM, GAMMA, SIGMA_A, SIGMA_B, THETAQ, FETCH,              &
&       SPHERICAL_RUN,              SHALLOW_RUN,                               &
&       REFRACTION_D_RUN,           REFRACTION_C_RUN,                          &
&       L_OBSTRUCTION,              L_DECOMP,                                  &
&       BETAMAX,                    IPHYS,                                     &
&       WAVE_BREAKING_RUN,          PHILLIPS_RUN,                              &
&       ISNONLIN,                   ITEST,                                     &
&       PROPAGATION_TIMESTEP,       PROPAGATION_TIMESTEP_UNIT,                 &
&       SOURCE_TIMESTEP,            SOURCE_TIMESTEP_UNIT,                      &
&       RESTART_SAVE_TIMESTEP,      RESTART_SAVE_TIMESTEP_UNIT,                &
&       RESTART_FILE_UNIT,          RESTART_FILE_NAME,                         &
&       COARSE_GRID_RUN,                                                       &
&       COARSE_OUTPUT_TIMESTEP,     COARSE_OUTPUT_TIMESTEP_UNIT,               &
&       COARSE_FILE_SAVE_TIMESTEP,  COARSE_FILE_SAVE_TIMESTEP_UNIT,            &
&       COARSE_OUTPUT_FILE_UNIT,    COARSE_OUTPUT_FILE_NAME,                   &
&       FINE_GRID_RUN,                                                         &
&       FINE_INPUT_FILE_UNIT,       FINE_INPUT_FILE_NAME,                      &
&       WIND_INPUT_TIMESTEP,        WIND_INPUT_TIMESTEP_UNIT,                  &
&       WIND_OUTPUT_TIMESTEP,       WIND_OUTPUT_TIMESTEP_UNIT,                 &
&       WIND_INPUT_FILE_UNIT,       WIND_INPUT_FILE_NAME,                      &
&       TOPO_INPUT_TIMESTEP,        TOPO_INPUT_TIMESTEP_UNIT,                  &
&       TOPO_OUTPUT_TIMESTEP,       TOPO_OUTPUT_TIMESTEP_UNIT,                 &
&       TOPO_INPUT_FILE_UNIT,       TOPO_INPUT_FILE_NAME,                      &
&       CURRENT_INPUT_TIMESTEP,     CURRENT_INPUT_TIMESTEP_UNIT,               &
&       CURRENT_OUTPUT_TIMESTEP,    CURRENT_OUTPUT_TIMESTEP_UNIT,              &
&       CURRENT_INPUT_FILE_UNIT,    CURRENT_INPUT_FILE_NAME,                   &
&       ICE_INPUT_TIMESTEP,         ICE_INPUT_TIMESTEP_UNIT,                   &
&       ICE_INPUT_FILE_UNIT,        ICE_INPUT_FILE_NAME,                       &
&       PARAMETER_OUTPUT_TIMESTEP,  PARAMETER_OUTPUT_TIMESTEP_UNIT,            &
&       PARAMETER_OUTPUT_FILE_UNIT, PARAMETER_OUTPUT_FILE_NAME,                &
&       SPECTRA_OUTPUT_TIMESTEP,    SPECTRA_OUTPUT_TIMESTEP_UNIT,              &
&       SPECTRA_OUTPUT_FILE_UNIT,   SPECTRA_OUTPUT_FILE_NAME,                  &
&       OUTPUT_FILE_SAVE_TIMESTEP,  OUTPUT_FILE_SAVE_TIMESTEP_UNIT,            &
&       COUTT,                                                                 &
&       FFLAG_P,  PFLAG_P,  FFLAG_S,  PFLAG_S, orientation_of_directions,      &
&       OUTLAT,   OUTLONG,   NAME,                                             &
&       spectral_code,              hours_2d_spectra,                          &
&       ready_file_flag,            ready_file_directory,   model_area,        &
&       ready_outfile_flag,         ready_outfile_directory,                   &
&       PREPROC_OUTPUT_FILE_UNIT,   PREPROC_OUTPUT_FILE_NAME,                  &
&       wammax_dur,                 wammax_dx,                                 &  !! WAM-MAX
&       wammax_dy,                                                             &  !! WAM-MAX
&       assimilation_flag,          influence_radius,                          &
&       observation_scatter,        model_scatter,                             &
&       assimilation_start_date,    assimilation_end_date,                     &
&       assimilation_time_step,     first_guess_output_flag,                   &
&       observation_file_unit,      observation_filename,                      &
&       first_guess_ip_file_unit,   first_guess_ip_filename,                   &
&       first_guess_sp_file_unit,   first_guess_sp_filename

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE CLEAR_WAM_USER_MODULE         !! SETS ALL DATA TO DEFAULT VALUES.
   MODULE PROCEDURE CLEAR_WAM_USER_MODULE
END INTERFACE
PUBLIC CLEAR_WAM_USER_MODULE

INTERFACE PRINT_WAM_NAMELIST           !! PRINTS WAM NAMELIST.
   MODULE PROCEDURE PRINT_WAM_NAMELIST
END INTERFACE
PUBLIC PRINT_WAM_NAMELIST

INTERFACE READ_WAM_NAMELIST            !! READS WAM NAMELIST.
   MODULE PROCEDURE READ_WAM_NAMELIST
END INTERFACE
PUBLIC READ_WAM_NAMELIST

INTERFACE SET_WAM_USER_PARAMETER       !! TRANSFERS ALL DATA TO MODULES.
   MODULE PROCEDURE SET_WAM_USER_PARAMETER
END INTERFACE
PUBLIC SET_WAM_USER_PARAMETER

INTERFACE CHANGE_TO_SECONDS            !! CHANGES HOURS OR MINUTES TO SECONDS.
   MODULE PROCEDURE CHANGE_TO_SECONDS
END INTERFACE
PRIVATE CHANGE_TO_SECONDS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CLEAR_WAM_USER_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   CLEAR_WAM_USER_MODULE - SETS ALL MODULE DATA TO DEFAULT VALUES.            !
!                                                                              !
!       H. GUNTHER   GKSS    DECEMBER 2009                                     !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       CLEARS THE WAM_USER_MODULE.                                            !
!                                                                              !
!     METHOD.                                                                  !
!     --------                                                                 !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. SET ALL DATA TO DEFAULT.                                              !
!        ------------------------                                              !

START_DATE = ' '  !! START DATE OF RUN.
END_DATE   = ' '  !! END   DATE OF RUN.

! ---------------------------------------------------------------------------- !

COLDSTART = .TRUE.  !! TRUE: COLDSTART, FALSE: HOT START.

IOPTI     = 1       !! COLD START OPTION:
                    !! =  0 WIND INDEPENDENT INITIAL VALUES.
                    !! =  1 WIND DEPENDENT INITIAL VALUES AND
                    !!      ENERGY EQUAL ZERO IF WINDSPEED IS ZERO
                    !! =  2 WIND DEPENDENT INITIAL VALUES AND
                    !!      ENERGY COMPUTED FROM GIVEN PARAMETERS IF
                    !!      WINDSPEED IS ZERO.
ALPHA     = 0.018   !! PHILLIPS' PARAMETER  (NOT USED IF IOPTI = 1).
FM        = 0.2     !! PEAK FREQUENCY (HZ) AND/OR MAXIMUM FREQUENCY.
GAMMA     = 3.0     !! OVERSHOOT FACTOR.
SIGMA_A   = 0.07    !! LEFT PEAK WIDTH.
SIGMA_B   = 0.09    !! RIGHT PEAK WIDTH.
THETAQ    = 0.0     !! WAVE DIRECTION (DEG) (NOT USED IF IOPTI = 1).
FETCH     = 30000.  !! FETCH IN METRES (IF ZERO THEN 0.5 OF THE
                    !! LATITUDE INCREMENT IS USED.).

! ---------------------------------------------------------------------------- !

SPHERICAL_RUN        = .TRUE.   !! TRUE:  SPHERICAL  PROPAGATION,
                                !! FALSE: CARTESIAN PROPAGATION.
SHALLOW_RUN          = .TRUE.   !! TRUE:  SHALLOW WATER MODEL,
                                !! FALSE: DEEP WATER MODEL.
REFRACTION_D_RUN     = .FALSE.  !! TRUE:  DEPTH REFRACTION ON.
REFRACTION_C_RUN     = .FALSE.  !! TRUE:  CURRENT REFRACTION ON.
L_OBSTRUCTION        = .FALSE.  !! FALSE: NO REDUCTION DUE TO SUB-GRID FEATURES
L_DECOMP             = .TRUE.   !! TRUE: 1D GRID DECOMPOSITION.
IPHYS                = 0        !! Janssen physics
WAVE_BREAKING_RUN    = .FALSE.  !! TRUE:  WAVE BREAKING ON.
PHILLIPS_RUN         = .FALSE.  !! TRUE:  PHILLIPS SOURCE ON.
ISNONLIN             = 0        !! OLD DEPTH SCALING FOR SNL.
ITEST                = 0        !! TEST OUTPUT UP TO LEVEL.
BETAMAX              = 1.20     !! PARAMETER FOR WIND INPUT.

! ---------------------------------------------------------------------------- !

PROPAGATION_TIMESTEP       = 0    !! PROPAGATION TIMESTEP.
PROPAGATION_TIMESTEP_UNIT  = 'S'  !! PROPAGATION TIMESTEP UNIT.
SOURCE_TIMESTEP            = 0    !! SOURCE TIME STEP.
SOURCE_TIMESTEP_UNIT       = 'S'  !! SOURCE TIME STEP UNIT.

! ---------------------------------------------------------------------------- !

RESTART_SAVE_TIMESTEP      =  0    !! RESTART FILE SAVE TIMESTEP.
RESTART_SAVE_TIMESTEP_UNIT = 'H'   !! RESTART FILE SAVE TIMESTEP UNIT.
RESTART_FILE_UNIT          = 17    !! RESTART FILE UNIT.
RESTART_FILE_NAME          = 'BLS' !! RESTART FILE IDENTIFIER.

! ---------------------------------------------------------------------------- !

COARSE_GRID_RUN                = .FALSE. !! TRUE: COARSE GRID MODEL.
COARSE_OUTPUT_TIMESTEP         = 0       !! OUTPUT TIMESTEP.
COARSE_OUTPUT_TIMESTEP_UNIT    = 'S'
COARSE_FILE_SAVE_TIMESTEP      = 0       !! TIMESTEP TO SAVE BOUNDARY FILE.
COARSE_FILE_SAVE_TIMESTEP_UNIT = 'S'
COARSE_OUTPUT_FILE_UNIT        = 70      !! COARSE SPECTRA OUTPUT FILE UNIT.
COARSE_OUTPUT_FILE_NAME        = 'CBO'   !! COARSE SPECTRA OUTPUT FILE IDENTIFIER.

FINE_GRID_RUN                  = .FALSE. !! TRUE: FINE GRID MODEL.
FINE_INPUT_FILE_UNIT           = 2       !! FINE SPECTRA INPUT FILE UNIT.
FINE_INPUT_FILE_NAME           = 'CBO'   !! FINE SPECTRA INPUT FILE IDENTIFIER.

! ---------------------------------------------------------------------------- !

WIND_INPUT_TIMESTEP         = 0
WIND_INPUT_TIMESTEP_UNIT    = 'S'
WIND_OUTPUT_TIMESTEP        = 0
WIND_OUTPUT_TIMESTEP_UNIT   = 'S'
WIND_INPUT_FILE_UNIT        = 1
WIND_INPUT_FILE_NAME        = ' '

! ---------------------------------------------------------------------------- !

TOPO_INPUT_TIMESTEP          = 0
TOPO_INPUT_TIMESTEP_UNIT     = 'S'
TOPO_OUTPUT_TIMESTEP         = 0
TOPO_OUTPUT_TIMESTEP_UNIT    = 'S'
TOPO_INPUT_FILE_UNIT         = 8
TOPO_INPUT_FILE_NAME         = ' '

! ---------------------------------------------------------------------------- !

CURRENT_INPUT_TIMESTEP       = 0
CURRENT_INPUT_TIMESTEP_UNIT  = 'S'
CURRENT_OUTPUT_TIMESTEP      = 0
CURRENT_OUTPUT_TIMESTEP_UNIT = 'S'
CURRENT_INPUT_FILE_UNIT      = 9
CURRENT_INPUT_FILE_NAME      = ' '

! ---------------------------------------------------------------------------- !

ICE_INPUT_TIMESTEP           = 0
ICE_INPUT_TIMESTEP_UNIT      = 'S'
ICE_INPUT_FILE_UNIT          = 3
ICE_INPUT_FILE_NAME          = ' '

! ---------------------------------------------------------------------------- !

PARAMETER_OUTPUT_TIMESTEP      = 1
PARAMETER_OUTPUT_TIMESTEP_UNIT = 'H'
PARAMETER_OUTPUT_FILE_UNIT     = 20
PARAMETER_OUTPUT_FILE_NAME     = 'MAP'
SPECTRA_OUTPUT_TIMESTEP        = 1
SPECTRA_OUTPUT_TIMESTEP_UNIT   = 'H'
SPECTRA_OUTPUT_FILE_UNIT       = 25
SPECTRA_OUTPUT_FILE_NAME       = 'OUT'
OUTPUT_FILE_SAVE_TIMESTEP      = 24
OUTPUT_FILE_SAVE_TIMESTEP_UNIT = 'H'


COUTT = ' '           !! SPECIFIED OUTPUT TIMES.

FFLAG_P     = .TRUE.  !! PARAMETER FILE OUTPUT FLAG.
PFLAG_P     = .FALSE. !! PARAMETER PRINTER OUTPUT FLAG.
orientation_of_directions = .true.   !! coming from or going to (default)

FFLAG_S     = .TRUE.  !! SPECTRA FILE OUTPUT FLAG.
PFLAG_S     = .FALSE. !! SPECTRA PRINTER OUTPUT FLAG.

OUTLAT      = ' '     !! LATITUDES OF OUTPUT SITES.
OUTLONG     = ' '     !! LONGITUDES OF OUTPUT SITES.
NAME        = ' '     !! OUTPUT SITES NAMES.

! ---------------------------------------------------------------------------- !

PREPROC_OUTPUT_FILE_UNIT = 7
PREPROC_OUTPUT_FILE_NAME = 'Grid_info'

! ---------------------------------------------------------------------------- !

ready_file_flag         = .false.
ready_file_directory    = 'wind'
model_area              = 'GSM'
ready_outfile_flag      = .false.
ready_outfile_directory = 'ready'
spectral_code           = -1
hours_2d_spectra        = -1

! ---------------------------------------------------------------------------- !
   
assimilation_flag           = 0      !! assimilation flag (1: assimilation)
assimilation_start_date     = ' '
assimilation_end_date       = ' '
assimilation_time_step      = 3
assimilation_time_step_unit = 'H'
influence_radius            = 3.0
observation_scatter         = 0.5
model_scatter               = 0.5
first_guess_output_flag     = .false.
observation_file_unit       = 80
observation_filename        = 'OBS'
first_guess_ip_file_unit    = 30
first_guess_ip_filename     = 'MAPFG'
first_guess_sp_file_unit    = 35
first_guess_sp_filename     = 'OUTFG'

! ---------------------------------------------------------------------------- ! !! WAM-MAX

wammax_dur       = 1200.0                                                        !! WAM-MAX
wammax_dx        = 1000.0                                                        !! WAM-MAX
wammax_dy        = 1000.0                                                        !! WAM-MAX
wammax_dur_unit  = 's'                                                           !! WAM-MAX
wammax_dx_unit   = 'm'                                                           !! WAM-MAX
wammax_dy_unit   = 'm'                                                           !! WAM-MAX

END SUBROUTINE CLEAR_WAM_USER_MODULE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE READ_WAM_NAMELIST (FILE, IOS)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    Arno Behrens    MSC/GKSS   January 2004                                   !
!                    module for namelist management, namelist read shifted     !
!                    to chief for message passing purposes,                    !
!                    namelist replaces the common user input file              !
!                                                                              !
!    Erik Myklebust             November 2004                                  !
!                    reading of namelist can now optionally be from file       !
!                                                                              !
!    Arno Behrens    GKSS       August 2010                                    !
!                    adding data assimilation parameters to namelist           !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTEGER, INTENT(IN)  :: FILE   !! > 0 READ FROM FILE ASSIGNED TO IU05
                               !! ELSE READ FROM STANDARD INPUT
INTEGER, INTENT(OUT) :: IOS    !! = 0 SUCCESSFULLY READ
                               !! ELSE READ ERRROR
			       
! ---------------------------------------------------------------------------- !

IOS = 0

IF (FILE .GT. 0) THEN
   READ (UNIT=IU05, NML=WAM_NAMELIST, IOSTAT=IOS)
ELSE
   READ (*, NML=WAM_NAMELIST,IOSTAT=IOS)
END IF

IF (IOS.EQ.0) THEN
   WRITE (IU06,*) '     SUB. READ_WAM_NAMELIST SUCCESSFULLY COMPLETED. '
END IF

END SUBROUTINE READ_WAM_NAMELIST

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_WAM_NAMELIST

WRITE (IU06,*) '  '
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '              WAM_USER_MODULE STATUS:'
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '  '
WRITE (IU06,NML=WAM_NAMELIST)
write (iu06,*) '  '

END SUBROUTINE PRINT_WAM_NAMELIST

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_WAM_USER_PARAMETER

! ---------------------------------------------------------------------------- !

if (itest>5) then
   CALL PRINT_WAM_NAMELIST
endif

! ---------------------------------------------------------------------------- !

CALL SET_INTEGRATION_PERIOD (B=START_DATE, E=END_DATE)
CALL SET_START_OPTION (COLDSTART)
CALL SET_C_START_PAR (IOPTI, ALPHA, FM, GAMMA, SIGMA_A, SIGMA_B, THETAQ, FETCH)
CALL SET_MODEL_OPTION (SPHERICAL     = SPHERICAL_RUN,                          &
&                      SHALLOW       = SHALLOW_RUN,                            &
&                      REFRACTION_D  = REFRACTION_D_RUN,                       &
&                      REFRACTION_C  = REFRACTION_C_RUN,                       &
&                      R_FACTOR      = L_OBSTRUCTION,                          &
&                      L_DEC         = L_DECOMP,                               &
&                      IWAMPHYS      = IPHYS,                                  &
&                      WAVE_BREAKING = WAVE_BREAKING_RUN,                      &
&                      PHILLIPS      = PHILLIPS_RUN,                           &
&                      INONLIN       = ISNONLIN)
CALL SET_TEST_OPTION (TEST=ITEST)

CALL SET_GENERAL_MODULE (B_MAX = BETAMAX)

! ---------------------------------------------------------------------------- !

CALL CHANGE_TO_SECONDS (PROPAGATION_TIMESTEP, PROPAGATION_TIMESTEP_UNIT)
CALL CHANGE_TO_SECONDS (SOURCE_TIMESTEP, SOURCE_TIMESTEP_UNIT)
CALL SET_INTEGRATION_TIMESTEPS (P=PROPAGATION_TIMESTEP,                        &
&                               S=SOURCE_TIMESTEP)

! ---------------------------------------------------------------------------- !

CALL CHANGE_TO_SECONDS (RESTART_SAVE_TIMESTEP, RESTART_SAVE_TIMESTEP_UNIT)
CALL SET_RESTART_FILE_STEP (STEP=RESTART_SAVE_TIMESTEP )
CALL SET_RESTART_FILE (NAME=RESTART_FILE_NAME,                                 &
&                      UNIT=RESTART_FILE_UNIT)

! ---------------------------------------------------------------------------- !

CALL SET_BOUNDARY_OPTION (C=COARSE_GRID_RUN,                                   &
&                         F=FINE_GRID_RUN)

CALL CHANGE_TO_SECONDS (COARSE_OUTPUT_TIMESTEP, COARSE_OUTPUT_TIMESTEP_UNIT)
CALL CHANGE_TO_SECONDS (COARSE_FILE_SAVE_TIMESTEP,                             &
&                       COARSE_FILE_SAVE_TIMESTEP_UNIT)
CALL SET_BOUNDARY_OUTPUT_TIMESTEPS (STEP_OUTPUT=COARSE_OUTPUT_TIMESTEP,        &
&                                   STEP_FILE=COARSE_FILE_SAVE_TIMESTEP)
CALL SET_B_OUTPUT_FILE (NAME=COARSE_OUTPUT_FILE_NAME,                          &
&                       UNIT=COARSE_OUTPUT_FILE_UNIT)
CALL SET_B_INPUT_FILE  (NAME=FINE_INPUT_FILE_NAME,                             &
&                       UNIT=FINE_INPUT_FILE_UNIT)

! ---------------------------------------------------------------------------- !

CALL CHANGE_TO_SECONDS (WIND_INPUT_TIMESTEP, WIND_INPUT_TIMESTEP_UNIT)
CALL CHANGE_TO_SECONDS (WIND_OUTPUT_TIMESTEP, WIND_OUTPUT_TIMESTEP_UNIT)
CALL SET_WIND_TIMESTEPS (IN=WIND_INPUT_TIMESTEP,                               &
&                        OUT=WIND_OUTPUT_TIMESTEP)
CALL SET_WIND_FILE (NAME=WIND_INPUT_FILE_NAME,                                 &
&                    UNIT=WIND_INPUT_FILE_UNIT)

! ---------------------------------------------------------------------------- !

CALL CHANGE_TO_SECONDS (TOPO_INPUT_TIMESTEP, TOPO_INPUT_TIMESTEP_UNIT)
CALL CHANGE_TO_SECONDS (TOPO_OUTPUT_TIMESTEP, TOPO_OUTPUT_TIMESTEP_UNIT)
CALL SET_TOPO_TIMESTEPS (IN=TOPO_INPUT_TIMESTEP,                               &
&                        OUT=TOPO_OUTPUT_TIMESTEP)
CALL SET_TOPO_FILE (NAME=TOPO_INPUT_FILE_NAME,                                 &
&                   UNIT=TOPO_INPUT_FILE_UNIT)

! ---------------------------------------------------------------------------- !

CALL CHANGE_TO_SECONDS (CURRENT_INPUT_TIMESTEP, CURRENT_INPUT_TIMESTEP_UNIT)
CALL CHANGE_TO_SECONDS (CURRENT_OUTPUT_TIMESTEP, CURRENT_OUTPUT_TIMESTEP_UNIT)
CALL SET_CURRENT_TIMESTEPS (IN=CURRENT_INPUT_TIMESTEP,                         &
&                           OUT=CURRENT_OUTPUT_TIMESTEP)
CALL SET_CURRENT_FILE (NAME=CURRENT_INPUT_FILE_NAME,                           &
&                      UNIT=CURRENT_INPUT_FILE_UNIT)

! ---------------------------------------------------------------------------- !

CALL CHANGE_TO_SECONDS (ICE_INPUT_TIMESTEP, ICE_INPUT_TIMESTEP_UNIT)
CALL SET_ICE_TIMESTEP (IN=ICE_INPUT_TIMESTEP)
CALL SET_ICE_FILE (NAME=ICE_INPUT_FILE_NAME,                                   &
&                  UNIT=ICE_INPUT_FILE_UNIT)

! ---------------------------------------------------------------------------- !

CALL CHANGE_TO_SECONDS (PARAMETER_OUTPUT_TIMESTEP,                             &
&                       PARAMETER_OUTPUT_TIMESTEP_UNIT)
CALL CHANGE_TO_SECONDS (SPECTRA_OUTPUT_TIMESTEP,SPECTRA_OUTPUT_TIMESTEP_UNIT)
CALL SET_OUTPUT_TIMES (INT=PARAMETER_OUTPUT_TIMESTEP,                          &
&                      SPE=SPECTRA_OUTPUT_TIMESTEP)

CALL CHANGE_TO_SECONDS (OUTPUT_FILE_SAVE_TIMESTEP,                             &
&                       OUTPUT_FILE_SAVE_TIMESTEP_UNIT)
CALL SET_OUTPUT_FILE_STEP (STEP=OUTPUT_FILE_SAVE_TIMESTEP)

CALL SET_MAP_FILE (NAME=PARAMETER_OUTPUT_FILE_NAME,                            &
&                  UNIT=PARAMETER_OUTPUT_FILE_UNIT) 
CALL SET_SPECTRA_FILE (NAME=SPECTRA_OUTPUT_FILE_NAME,                          &
&                      UNIT=SPECTRA_OUTPUT_FILE_UNIT)

CALL SET_OUTPUT_TIMES (TIME=COUTT)

CALL SET_PARAMETER_OUTPUT_FLAGS (PF=PFLAG_P, FF=FFLAG_P, od=                   &
&                                orientation_of_directions)
CALL SET_SPECTRA_OUTPUT_FLAGS (PF=PFLAG_S, FF=FFLAG_S)

CALL SET_OUTPUT_SITES (LONG=OUTLONG, LAT=OUTLAT, NA=NAME)

! ---------------------------------------------------------------------------- !

CALL SET_PREPROC_FILE (NAME=PREPROC_OUTPUT_FILE_NAME,                          &
&                      UNIT=PREPROC_OUTPUT_FILE_UNIT)

! ---------------------------------------------------------------------------- !

call set_ready_file_flag (p=ready_file_flag)
call set_ready_file_directory (name=ready_file_directory, ar=model_area)
call set_ready_outfile_flag (op=ready_outfile_flag)
call set_ready_outfile_directory (oname=ready_outfile_directory)
call set_spectral_code (ihour=hours_2d_spectra, icode=spectral_code)

! ---------------------------------------------------------------------------- ! !! WAM-MAX
!                                                                                !! WAM-MAX
call set_wammax_options ( dt=wammax_dur, dx=wammax_dx, dy=wammax_dy)             !! WAM-MAX
!                                                                                !! WAM-MAX
! ---------------------------------------------------------------------------- !

call set_assimilation_option (as=assimilation_flag, ra=influence_radius,       &
&                             so=observation_scatter, sm=model_scatter)
CALL CHANGE_TO_SECONDS (assimilation_time_step,                                &
&                       assimilation_time_step_UNIT)
call set_assimilation_period (b=assimilation_start_date,                       &
&                             e=assimilation_end_date,                         &
&                             step=assimilation_time_step)
call set_assimilation_output (p=first_guess_output_flag)
call set_observation_file  (name=observation_filename,                         &
&                           unit=observation_file_unit)
call set_assi_map_file     (name=first_guess_ip_filename,                      &
&                           unit=first_guess_ip_file_unit)
call set_assi_spectra_file (name=first_guess_sp_filename,                      &
&                           unit=first_guess_sp_file_unit)

END SUBROUTINE SET_WAM_USER_PARAMETER


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CHANGE_TO_SECONDS (TIME, UNIT)

INTEGER,           INTENT(INOUT) :: TIME
CHARACTER (LEN=1), INTENT(INOUT) :: UNIT

IF (UNIT.EQ.'M' .OR. UNIT.EQ.'m') TIME = TIME*60.
IF (UNIT.EQ.'H' .OR. UNIT.EQ.'h') TIME = TIME*3600.
UNIT = 'S'

END SUBROUTINE CHANGE_TO_SECONDS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_USER_MODULE

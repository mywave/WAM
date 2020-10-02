MODULE PREPROC_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL VARIABLES AND CONSTANTS NECESSARY FOR THE         !
!   PREPROC PROGRAM. ALL PROCEDURES ARE INCLUDED TO COMPUTE THE INFORMATION    !
!   STORED IN WAM_CONST_MODULE, WAM_TABLE_MODULE WAM_NEST_MODULE AND           !
!   WAM_GRID_MODULE.                                                           !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE          !! COORDINATE TYPE AND PROCEDURES

USE WAM_GENERAL_MODULE,   ONLY:  &
&       AKI,                     &
&       ABORT1,                  & !! TERMINATES PROCESSING.
&       PRINT_ARRAY,             & !! PRINTS AN ARRAY.
&       PRINT_GENERAL_MODULE

USE WAM_GRID_MODULE,      ONLY:  &
&       EQUAL_TO_M_GRID,         & !! COMPARES TWO GRIDS.
&       PRINT_GRID_STATUS          !! PRINTS THE MODULE INFORMATION.

USE WAM_FRE_DIR_MODULE,   ONLY:  &
&       PRINT_FRE_DIR_STATUS       !! PRINTS THE MODULE INFORMATION.

USE WAM_NEST_MODULE,  ONLY:      &
&       PREPARE_BOUNDARY_nest,   & !! MAKES NEST BOUNDARIES.
&       PRINT_NEST_STATUS          !! PRINTS THE MODULE INFORMATION.

USE WAM_TABLES_MODULE,    ONLY:  &
&       PRINT_TABLES_STATUS        !! PRINTS THE MODULE INFORMATION.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST, IU07, FILE07

USE WAM_GENERAL_MODULE, ONLY: ZPI, DEG, RAD, CIRC

USE WAM_FRE_DIR_MODULE, ONLY: KL, ML, FR, CO, TH, DELTH, DELTR, COSTH, SINTH,  &
&                             GOM, C, INV_LOG_CO,                              &
&                             DF, DF_FR, DF_FR2,                               &
&                             DFIM, DFIMOFR, DFIM_FR, DFIM_FR2, FR5, FRM5,     &
&                             RHOWG_DFIM,                                      &
&                             FMIN, MO_TAIL, MM1_TAIL, MP1_TAIL, MP2_TAIL,     &
&                             MPM, KPM, JXO, JYO

USE WAM_GRID_MODULE,    ONLY: HEADER, NX, NY, NSEA, NLON_RG, IPER,             &
&                             AMOWEP, AMOSOP, AMOEAP, AMONOP,                  &
&                             XDELLA, XDELLO, DELLAM, ZDELLO, DELPHI,          &
&                             SINPH, COSPH, DEPTH_B, KLAT, KLON, WLAT,         &
&                             IXLG, KXLT, L_S_MASK, ONE_POINT, REDUCED_GRID,   &
&                             OBSLAT, OBSLON

USE WAM_NEST_MODULE,    ONLY: N_NEST, MAX_NEST, N_NAME, n_code,                &
&                             NBOUNC, IJARC, BLATC, BLNGC, DLAMAC, DPHIAC,     &
&                             N_SOUTH, N_NORTH, N_EAST, N_WEST, N_ZDEL,        &
&                             NBINP, NBOUNF, C_NAME, BLNGF, BLATF,             &
&                             IJARF, IBFL, IBFR, BFW

USE WAM_TABLES_MODULE,  ONLY: NDEPTH, DEPTHA, DEPTHD, DEPTHE,                  &
&                             FLMINFR, TCGOND, TFAK, TSIHKD, TFAC_ST, T_TAIL,  &
&                             DELU

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

LOGICAL :: SET_STATUS    = .FALSE. !! .TRUE. IF USER INPUT IS DEFINED.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TOPOGRAPHY INPUT DATA AS PROVIDED BY USER.                            !
!        ------------------------------------------                            !

INTEGER              :: NX_T = -1     !! NUMBER OF LONGITUDES.
INTEGER              :: NY_T = -1     !! NUMBER OF LATITUDES.
LOGICAL              :: PER_T         !! .TRUE. IF GRID IS PERIODICAL.
INTEGER              :: DY_T          !! LATITUDE INCREMENT [M_SEC].
INTEGER              :: DX_T          !! LONGITUDE INCREMENT [M_SEC].
INTEGER              :: SOUTH_T       !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER              :: NORTH_T       !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER              :: WEST_T        !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER              :: EAST_T        !! EAST LONGITUDE OF GRID [M_SEC].
LOGICAL :: EQUAL_GRID =.FALSE.        !! .TRUE. IF GRID IS EQUAL TO MODEL GRID.
LOGICAL              :: L_INTERPOL_T  !! INTERPOLATION OPTION FOR MODEL DEPTH.
LOGICAL, PUBLIC      :: L_OBSTRUCTION_T   !! OBSTRUCTION FACTORS DUE TO SUB-GRID FEATURES.
REAL,    ALLOCATABLE :: GRID_IN(:,:)  !! WATER DEPTH [M].
REAL                 :: LAND_LIMIT    !! DEPTH <= LAND_LIMIT ARE NOT ACTIVE

INTEGER, ALLOCATABLE, DIMENSION(:) :: ALON        !! longitudes in INPUT GRID
INTEGER, ALLOCATABLE, DIMENSION(:) :: ALAT        !! latitudes  in INPUT GRID

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. GRID CORRECTION AREAS.                                                !
!        ----------------------                                                !

INTEGER              :: N_CA = 0      !! NO. OF CORRECTION AREAS
INTEGER, ALLOCATABLE :: SOUTH_CA(:)   !! S - LATITUDE OF AREA
INTEGER, ALLOCATABLE :: NORTH_CA(:)   !! N - LATITUDE OF AREA
INTEGER, ALLOCATABLE :: WEST_CA (:)   !! W - LONGITUIDE OF AREA
INTEGER, ALLOCATABLE :: EAST_CA (:)   !! E - LONGITUIDE OF AREA
REAL,    ALLOCATABLE :: DEPTH_CA(:)   !! DEPTH OF AREA

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. GRIDDED MODEL DEPTH DATA.                                             !
!        -------------------------                                             !

REAL,    ALLOCATABLE, DIMENSION(:,:) :: GRD             !! GRIDDED TOPOGRAPHY
REAL,    ALLOCATABLE, DIMENSION(:,:) :: PERCENTLAND     !! percent of land points in WAM grid cell
REAL,    ALLOCATABLE, DIMENSION(:,:) :: PERCENTSHALLOW  !! percent of shallow points in WAM grid cell
INTEGER, ALLOCATABLE, DIMENSION(:)   :: XLAT            !! LATITUDES OF GRID.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !


INTERFACE SET_GRID_CORRECTIONS       !! TRANSFERS GRID CORRECTIONS TO MODULE.
   MODULE PROCEDURE SET_GRID_CORRECTIONS_C  !! TEXT COORDINATES
   MODULE PROCEDURE SET_GRID_CORRECTIONS_D  !! DEGREE COORDINATES
   MODULE PROCEDURE SET_GRID_CORRECTIONS_M  !! M_SEC
END INTERFACE
PUBLIC SET_GRID_CORRECTIONS

INTERFACE SET_GRID_DEF               !! TRANSFERS MODEL GRID DEFINITIONS 
   MODULE PROCEDURE SET_GRID_DEF_C     !! TEXT COORDINATES
   MODULE PROCEDURE SET_GRID_DEF_D     !! DEGREE COORDINATES
   MODULE PROCEDURE SET_GRID_DEF_M     !! M_SEC
END INTERFACE
PUBLIC SET_GRID_DEF

INTERFACE SET_HEADER                 !! TRANSFERS MODEL HEADER 
   MODULE PROCEDURE SET_HEADER
END INTERFACE
PUBLIC SET_HEADER

INTERFACE SET_TOPOGRAPHY             !! TRANSFERS DEPTH DATA TO MODULE.
   MODULE PROCEDURE SET_TOPOGRAPHY_C   !! TEXT COORDINATES
   MODULE PROCEDURE SET_TOPOGRAPHY_D   !! DEGREE COORDINATES
   MODULE PROCEDURE SET_TOPOGRAPHY_M   !! M_SEC
END INTERFACE
PUBLIC SET_TOPOGRAPHY

INTERFACE PREPARE_CONST              !! COMPUTES WAM_CONST_MODULE DATA.
   MODULE PROCEDURE PREPARE_CONST
END INTERFACE
PUBLIC PREPARE_CONST

INTERFACE PRINT_PREPROC_STATUS        !! PRINTS PREPROC STATUS.
   MODULE PROCEDURE PRINT_PREPROC_STATUS
END INTERFACE
PUBLIC PRINT_PREPROC_STATUS

INTERFACE WRITE_PREPROC_FILE          !! WRITES PREPROC OUTPUT FILE.
   MODULE PROCEDURE WRITE_PREPROC_FILE
END INTERFACE
PUBLIC WRITE_PREPROC_FILE

INTERFACE
   SUBROUTINE READ_PREPROC_USER        !! READ USER INPUT FOR PREPROC.
   END SUBROUTINE READ_PREPROC_USER

   SUBROUTINE READ_TOPOGRAPHY          !! READ TOPOGRAPHY INPUT FILE
   END SUBROUTINE READ_TOPOGRAPHY
END INTERFACE
PUBLIC READ_PREPROC_USER, READ_TOPOGRAPHY

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE ADJUST_DEPTH                  !! ADJUST DEPTH.
   MODULE PROCEDURE ADJUST_DEPTH
END INTERFACE
PRIVATE ADJUST_DEPTH

INTERFACE CLOSEST_GP_LAT                !! FINDS CLOSEST AND SECOND CLOSEST GRID
   MODULE PROCEDURE CLOSEST_GP_LAT      !! POINT AT N OR S LATITUDE NEIGHTBOUR.
END INTERFACE
PRIVATE CLOSEST_GP_LAT

INTERFACE COUNT_POINTS_E_W             !! ROUTINE TO COUNT POINTS.
   MODULE PROCEDURE COUNT_POINTS_E_W
END INTERFACE
PRIVATE COUNT_POINTS_E_W

INTERFACE COUNT_POINTS_N_S             !! ROUTINE  TO COUNT POINTS..
   MODULE PROCEDURE COUNT_POINTS_N_S
END INTERFACE
PRIVATE COUNT_POINTS_N_S

INTERFACE CREATE_OBSTRUCTIONS          !! CREATES OBSTRUCTION DUE TO SUB-GRID FEATRURES.
MODULE PROCEDURE CREATE_OBSTRUCTIONS
END INTERFACE
PRIVATE CREATE_OBSTRUCTIONS

INTERFACE MEAN_DEPTH                   !! ROUTINE TO ARRANGE WAMODEL GRID.
MODULE PROCEDURE MEAN_DEPTH
END INTERFACE
PRIVATE MEAN_DEPTH

INTERFACE PREPARE_GRID                 !! ROUTINE TO ARRANGE WAMODEL GRID.
   MODULE PROCEDURE PREPARE_GRID
END INTERFACE
PRIVATE PREPARE_GRID

INTERFACE SUBGRID_TOPOGRAPHY           !! ARRANGE SUBGRID TOPOGRAPHY.
   MODULE PROCEDURE SUBGRID_TOPOGRAPHY
END INTERFACE
PRIVATE SUBGRID_TOPOGRAPHY

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_GRID_CORRECTIONS_C (SOUTH, NORTH, WEST, EAST, D_COR)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=LEN_COOR), INTENT(IN)  :: SOUTH(:)  !! SOUTH LATITUDES.
CHARACTER (LEN=LEN_COOR), INTENT(IN)  :: NORTH(:)  !! NORTH LATITUDES.
CHARACTER (LEN=LEN_COOR), INTENT(IN)  :: WEST(:)   !! WEST LONGITUDES.
CHARACTER (LEN=LEN_COOR), INTENT(IN)  :: EAST(:)   !! EAST LONGITUDES.
REAL,                     INTENT(IN)  :: D_COR(:)  !! DEPTH IN CORR. AREAS [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CONVERT TO M_SEC AND CALL SET_GRID_CORRECTIONS_M.                     !
!        ------------------------------------------------                      !

CALL SET_GRID_CORRECTIONS_M (SOUTH=READ_COOR_TEXT(SOUTH),                      &
&                            NORTH=READ_COOR_TEXT(NORTH),                      &
&                            WEST=READ_COOR_TEXT(WEST),                        &
&                            EAST=READ_COOR_TEXT(EAST),                        &
&                            D_COR=D_COR)

END SUBROUTINE SET_GRID_CORRECTIONS_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_GRID_CORRECTIONS_D (SOUTH, NORTH, WEST, EAST, D_COR)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL (KIND=KIND_D), INTENT(IN)  :: SOUTH(:) !! SOUTH LAT. OF CORR. AREAS [DEG].
REAL (KIND=KIND_D), INTENT(IN)  :: NORTH(:) !! NORTH LAT. OF CORR. AREAS [DEG].
REAL (KIND=KIND_D), INTENT(IN)  :: WEST(:)  !! WEST LONG. OF CORR. AREAS [DEG].
REAL (KIND=KIND_D), INTENT(IN)  :: EAST(:)  !! EAST LONG. OF CORR. AREAS [DEG].
REAL,               INTENT(IN)  :: D_COR(:) !! DEPTH IN CORR. AREAS [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CONVERT TO M_SEC AND CALL SET_GRID_CORRECTIONS_M.                     !
!        ------------------------------------------------                      !

CALL SET_GRID_CORRECTIONS_M (SOUTH=DEG_TO_M_SEC(SOUTH),                        &
&                            NORTH=DEG_TO_M_SEC(NORTH),                        &
&                            WEST=DEG_TO_M_SEC(WEST),                          &
&                            EAST=DEG_TO_M_SEC(EAST),                          &
&                            D_COR=D_COR)

END SUBROUTINE SET_GRID_CORRECTIONS_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_GRID_CORRECTIONS_M (SOUTH, NORTH, WEST, EAST, D_COR)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)  :: SOUTH(:)   !! SOUTH LATITUDES OF CORR. AREAS [M_SEC].
INTEGER, INTENT(IN)  :: NORTH(:)   !! NORTH LATITUDES OF CORR. AREAS [M_SEC].
INTEGER, INTENT(IN)  :: WEST(:)    !! WEST LONGITUDES OF CORR. AREAS [M_SEC].
INTEGER, INTENT(IN)  :: EAST(:)    !! EAST LONGITUDES OF CORR. AREAS [M_SEC].
REAL,    INTENT(IN)  :: D_COR(:)   !! DEPTH IN CORR. AREAS [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLEAR CORRECTION AREAS.                                               !
!        -----------------------                                               !

N_CA = 0
IF (ALLOCATED(SOUTH_CA)) DEALLOCATE(SOUTH_CA)
IF (ALLOCATED(NORTH_CA)) DEALLOCATE(NORTH_CA)
IF (ALLOCATED(WEST_CA )) DEALLOCATE(WEST_CA)
IF (ALLOCATED(EAST_CA )) DEALLOCATE(EAST_CA)
IF (ALLOCATED(DEPTH_CA)) DEALLOCATE(DEPTH_CA)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COPY CORRECTION AREAS DEFINITIONS.                                    !
!        ----------------------------------                                    !

N_CA = COUNT (SOUTH.NE.COOR_UNDEF)
IF (N_CA.GT.0) THEN
   ALLOCATE(SOUTH_CA(1:N_CA))
   ALLOCATE(NORTH_CA(1:N_CA))
   ALLOCATE(WEST_CA(1:N_CA))
   ALLOCATE(EAST_CA(1:N_CA))
   ALLOCATE(DEPTH_CA(1:N_CA))

   SOUTH_CA(1:N_CA) = SOUTH(1:N_CA)
   NORTH_CA(1:N_CA) = NORTH(1:N_CA)
   WEST_CA(1:N_CA)  = WEST (1:N_CA)
   EAST_CA(1:N_CA)  = EAST (1:N_CA)
   DEPTH_CA(1:N_CA) = D_COR(1:N_CA)

   CALL ADJUST (WEST_CA, EAST_CA)
END IF

END SUBROUTINE SET_GRID_CORRECTIONS_M

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_GRID_DEF_C (N_LON, N_LAT, D_LON, D_LAT,                         &
&                          SOUTH, NORTH, WEST, EAST,                           &
&                          LAND, R_GRID, L_INTERPOL, L_OBSTRUCTION)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,                  INTENT(IN) :: N_LON   !! NUMBER OF LONGITUDES IN GRID.
INTEGER,                  INTENT(IN) :: N_LAT   !! NUMBER OF LATITUDES IN GRID.
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: D_LON   !! LONGITUDE INCREMENT OF GRID.
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: D_LAT   !! LATITUDE INCREMENT OF GRID.
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: SOUTH   !! SOUTH LATITUDE OF GRID.
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: NORTH   !! NORTH LATITUDE OF GRID.
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: WEST    !! WEST LONGITUDE OF GRID.
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: EAST    !! EAST LONGITUDE OF GRID.
REAL,      INTENT(IN) :: LAND     !! DEPTH >= LAND ARE SEAPOINTS [M].
LOGICAL,   INTENT(IN) :: R_GRID   !! REDUCED GRID OPTION.
LOGICAL,   INTENT(IN) :: L_INTERPOL     !! INTERPOLATION OPTION FOR MODEL DEPTH.
LOGICAL,   INTENT(IN) :: L_OBSTRUCTION !! REDUCTION FACTORS DUE TO SUB-GRID FEATURES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CONVERT TO M_SEC AND CALL SET_GRID_CORRECTIONS_M.                     !
!        -------------------------------------------------                     !

CALL SET_GRID_DEF_M (N_LON, N_LAT,                                             &
&                    READ_COOR_TEXT(D_LON), READ_COOR_TEXT(D_LAT),             &
&                    READ_COOR_TEXT(SOUTH), READ_COOR_TEXT(NORTH),             &
&                    READ_COOR_TEXT(WEST),  READ_COOR_TEXT(EAST),              &
&                    LAND, R_GRID, L_INTERPOL, L_OBSTRUCTION)

END SUBROUTINE SET_GRID_DEF_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_GRID_DEF_D (N_LON, N_LAT, D_LON, D_LAT,                         &
&                          SOUTH, NORTH, WEST, EAST,                           &
&                          LAND, R_GRID, L_INTERPOL, L_OBSTRUCTION)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,            INTENT(IN) :: N_LON   !! NUMBER OF LONGITUDES IN GRID.
INTEGER,            INTENT(IN) :: N_LAT   !! NUMBER OF LATITUDES IN GRID.
REAL (KIND=KIND_D), INTENT(IN) :: D_LON   !! LONGITUDE INCREMENT OF GRID.
REAL (KIND=KIND_D), INTENT(IN) :: D_LAT   !! LATITUDE INCREMENT OF GRID.
REAL (KIND=KIND_D), INTENT(IN) :: SOUTH   !! SOUTH LATITUDE OF GRID.
REAL (KIND=KIND_D), INTENT(IN) :: NORTH   !! NORTH LATITUDE OF GRID.
REAL (KIND=KIND_D), INTENT(IN) :: WEST    !! WEST LONGITUDE OF GRID.
REAL (KIND=KIND_D), INTENT(IN) :: EAST    !! EAST LONGITUDE OF GRID.
REAL,      INTENT(IN) :: LAND     !! DEPTH >= LAND ARE SEAPOINTS [M].
LOGICAL,   INTENT(IN) :: R_GRID   !! REDUCED GRID OPTION.
LOGICAL,   INTENT(IN) :: L_INTERPOL     !! INTERPOLATION OPTION FOR MODEL DEPTH.
LOGICAL,   INTENT(IN) :: L_OBSTRUCTION !! REDUCTION FACTORS DUE TO SUB-GRID FEATURES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CONVERT TO M_SEC AND CALL SET_GRID_CORRECTIONS_M.                     !
!        -------------------------------------------------                     !

CALL SET_GRID_DEF_M (N_LON, N_LAT,                                             &
&                    DEG_TO_M_SEC(D_LON), DEG_TO_M_SEC(D_LAT),                 &
&                    DEG_TO_M_SEC(SOUTH), DEG_TO_M_SEC(NORTH),                 &
&                    DEG_TO_M_SEC(WEST),  DEG_TO_M_SEC(EAST),                  &
&                    LAND, R_GRID, L_INTERPOL, L_OBSTRUCTION)

END SUBROUTINE SET_GRID_DEF_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_GRID_DEF_M (N_LON, N_LAT, D_LON, D_LAT,                         &
&                          SOUTH, NORTH, WEST, EAST,                           &
&                          LAND, R_GRID, L_INTERPOL, L_OBSTRUCTION)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,           INTENT(IN) :: N_LON   !! NUMBER OF LONGITUDES IN GRID.
INTEGER,           INTENT(IN) :: N_LAT   !! NUMBER OF LATITUDES IN GRID.
INTEGER,           INTENT(IN) :: D_LON   !! LONGITUDE INCREMENT OF GRID.
INTEGER,           INTENT(IN) :: D_LAT   !! LATITUDE INCREMENT OF GRID.
INTEGER,           INTENT(IN) :: SOUTH   !! SOUTH LATITUDE OF GRID.
INTEGER,           INTENT(IN) :: NORTH   !! NORTH LATITUDE OF GRID.
INTEGER,           INTENT(IN) :: WEST    !! WEST LONGITUDE OF GRID.
INTEGER,           INTENT(IN) :: EAST    !! EAST LONGITUDE OF GRID.
REAL,      INTENT(IN) :: LAND          !! DEPTH >= LAND ARE SEAPOINTS [M].
LOGICAL,   INTENT(IN) :: R_GRID        !! REDUCED GRID OPTION.
LOGICAL,   INTENT(IN) :: L_INTERPOL     !! INTERPOLATION OPTION FOR MODEL DEPTH.
LOGICAL,   INTENT(IN) :: L_OBSTRUCTION !! REDUCTION FACTORS DUE TO SUB-GRID FEATURES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

LOGICAL  :: ERROR = .FALSE.              !! ERROR FLAG
character (len=len_coor) :: formtext

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLEAR GRID DEFINITIONS.                                               !
!        -----------------------                                               !

IPER         = .FALSE. !! PERIODIC GRID

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COPY GRID DEFINITIONS AND CONVERT TO M_SEC.                           !
!        -------------------------------------------                           !

L_INTERPOL_T = L_INTERPOL
L_OBSTRUCTION_T =  L_OBSTRUCTION
AMOWEP = WEST
AMOSOP = SOUTH
AMOEAP = EAST
AMONOP = NORTH
XDELLO = D_LON
XDELLA = D_LAT
NX = N_LON
NY = N_LAT

REDUCED_GRID = R_GRID
LAND_LIMIT = LAND

IF (NX .EQ.-1 .AND. NY .EQ.-1 .AND.                                      &
&   XDELLO.EQ.COOR_UNDEF .AND. XDELLA.EQ.COOR_UNDEF .AND.                      &
&   AMOSOP.EQ.COOR_UNDEF .AND. AMONOP.EQ.COOR_UNDEF .AND.                      &
&   AMOWEP.EQ.COOR_UNDEF .AND. AMOEAP.EQ.COOR_UNDEF ) THEN
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' +         INFORMATION FROM SUB. SET_GRID_DEF             +'
   WRITE (IU06,*) ' +         ==================================             +'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' +     ALL MODEL GRID SPECIFICATIONS ARE UNDEFINED.       +'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' +            PROGRAM WILL CONTINUE                       +'
   WRITE (IU06,*) ' +  USING THE DEFINITIONS FROM THE TOPOGRAPHY INPUT FILE  +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   SET_STATUS = (NX.GT.0 .AND. NY.GT.0 .AND. NX_T.GT.0 .AND. NY_T.GT.0)
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MISSING VALUES.                                               !
!        ----------------------                                                !

ERROR = .FALSE.
IF (AMOWEP.EQ.COOR_UNDEF) ERROR = .TRUE.
IF (AMOSOP.EQ.COOR_UNDEF) ERROR = .TRUE.

IF (.NOT.ERROR) THEN
   IF (NX.EQ.1 .OR. XDELLO.EQ.0 .OR. AMOEAP.EQ.AMOWEP) THEN
      NX = 1
      XDELLO = M_S_PER
      AMOEAP = AMOWEP
    END IF
    IF (NY.EQ.1 .OR. XDELLA.EQ.0 .OR. AMOSOP.EQ.AMONOP) THEN
      NY = 1
      XDELLA = 0
      AMONOP = AMOSOP
    END IF
END IF

CALL CHECK_GRID_DEFINITION (AMOWEP, AMOSOP, AMOEAP, AMONOP,                   &
&                           XDELLO, XDELLA, NX, NY, ERROR)

IF (ERROR) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           FATAL ERROR IN SUB. SET_GRID_DEF             *'
   WRITE (IU06,*) ' *           ================================             *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *         A MODEL GRID COULD NOT BE DEFINED.             *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF

IPER = PERIODIC (AMOWEP, AMOEAP, XDELLO, NX)   !! PERIODIC ?
ONE_POINT = NX.EQ.1 .AND. NY.EQ.1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CHECK CONSISTENCY.                                                    !
!        ------------------                                                    !

IF ((NX.LE.1 .OR. NY.LE.1) .AND. (NX.NE.1 .OR. NY.LT.1)) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           FATAL ERROR IN SUB. SET_GRID_DEF             *'
   WRITE (IU06,*) ' *           ================================             *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *  A MODEL GRID MUST HAVE                                *'
   WRITE (IU06,*) ' *  MORE THAN ONE POINT ON BOTH GRID AXES                 *'
   WRITE (IU06,*) ' *                OR                                      *'
   WRITE (IU06,*) ' *  ONE POINT ON BOTH GRID AXES (PROPAGATION IS NOT DONE) *'
   WRITE (IU06,*) ' *                OR                                      *'
   WRITE (IU06,*) ' *  ONE POINT ON LONGITUDE AXIS AND MORE THAN ONE POINT   *'
   WRITE (IU06,*) ' *  ON LATITUDE AXIS AND THE LONGITUDE INCREMENT MUST     *'
   WRITE (IU06,*) ' *  BE 360 DEG.  (PERIODIC GRID IN LONGITUDE)             *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' * USER PROVIDED GRID DEFINITIONS ARE:                    *'
   formtext = write_coor_text (coor_undef)
   WRITE (IU06,*) ' *  (-1 AND ', formtext,                                    &
&                                          ' INDICATED UNDEFINED VALUES)     *'
   WRITE (IU06,*) ' *                                                        *'
   formtext = write_coor_text (amowep)
   WRITE (IU06,*) ' * LONGITUDE             WEST = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (amowep), 'deg'
   formtext = write_coor_text (amoeap)
   WRITE (IU06,*) ' * LONGITUDE             EAST = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (amoeap), 'deg'
   formtext = write_coor_text (xdello)
   WRITE (IU06,*) ' * LONGITUDE INCREMENT  D_LON = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (xdello), 'deg'
   WRITE (IU06,*) ' * NO. OF LONGITUDES    N_LON = ', NX
   formtext = write_coor_text (amosop)
   WRITE (IU06,*) ' * LATITUDE             SOUTH = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (amosop), 'deg'
   formtext = write_coor_text (amonop)
   WRITE (IU06,*) ' * LATITUDE             NORTH = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (amonop), 'deg'
   formtext = write_coor_text (xdella)
   WRITE (IU06,*) ' * LATITUDE  INCREMENT  D_LAT = ', formtext,                &
&                                      ' = ',M_SEC_TO_DEG (xdella), 'deg'
   WRITE (IU06,*) ' * NO. OF LATITUDE      N_LAT = ', NY
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. INPUT MODULE DATA ARE DEFINED.                                        !
!        ------------------------------                                        !

IF (L_INTERPOL_T .AND. L_OBSTRUCTION_T) THEN

   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' +         WARNING FROM SUB. SET_GRID_DEF                 +'
   WRITE (IU06,*) ' +         ==============================                 +'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' + MODEL DEPTH IS NEAREST NEIGHBOUR DEPTH OF INPUT',      &
&                     ' TOPOGRAPHY:    L_INTERPOL = ', L_INTERPOL_T
   WRITE (IU06,*) ' + OBSTRUCTION FACTORS DUE TO SUB-GRID FEATURES ARE',      &
&                     ' REQUESTED: L_OBSTRUCTION = ', L_OBSTRUCTION
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' +                                                        +'
   WRITE (IU06,*) ' +            PROGRAM WILL CONTINUE                       +'
   WRITE (IU06,*) ' +         WITHOUT OBSTRUCTION FACTORS                    +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   L_OBSTRUCTION_T = .FALSE.
END IF
! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. INPUT MODULE DATA ARE DEFINED.                                        !
!        ------------------------------                                        !

SET_STATUS = (NX.GT.0 .AND. NY.GT.0 .AND. NX_T.GT.0 .AND. NY_T.GT.0)

END SUBROUTINE SET_GRID_DEF_M

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_HEADER (TEXT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=*), INTENT(IN)  :: TEXT      !! MODEL HEADER.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COPY HEADER DEFINITION.                                               !
!        -----------------------                                               !

HEADER = TEXT(1:MIN(LEN_TRIM(TEXT), LEN(HEADER)))

END SUBROUTINE SET_HEADER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_TOPOGRAPHY_D (N_LON, N_LAT, D_LON, D_LAT,                       &
&                            SOUTH, NORTH, WEST, EAST, D_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,            INTENT(IN) :: N_LON    !! NUMBER OF LONGITUDES IN GRID.
INTEGER,            INTENT(IN) :: N_LAT    !! NUMBER OF LATITUDES IN GRID.
REAL (KIND=KIND_D), INTENT(IN) :: D_LON    !! LONGITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN) :: D_LAT    !! LATITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN) :: SOUTH    !! SOUTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN) :: NORTH    !! NORTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN) :: WEST     !! WEST LONGITUDE OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN) :: EAST     !! EAST LONGITUDE OF GRID [DEG].
REAL,               INTENT(IN) :: D_MAP(1:N_LON,1:N_LAT) !! WATER DEPTH [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CONVERT TO M_SEC AND CALL SET_TOPOGRAPHY_M.                           !
!        -------------------------------------------                           !

CALL SET_TOPOGRAPHY_M (N_LON, N_LAT,                                           & 
&                      DEG_TO_M_SEC(D_LON),  DEG_TO_M_SEC(D_LAT),              &
&                      DEG_TO_M_SEC(SOUTH),  DEG_TO_M_SEC(NORTH),              &
&                      DEG_TO_M_SEC(WEST),   DEG_TO_M_SEC(EAST),               &
&                      D_MAP)

END SUBROUTINE SET_TOPOGRAPHY_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_TOPOGRAPHY_M (N_LON, N_LAT, D_LON, D_LAT,                       &
&                            SOUTH, NORTH, WEST, EAST, D_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN) :: N_LON      !! NUMBER OF LONGITUDES IN GRID.
INTEGER, INTENT(IN) :: N_LAT      !! NUMBER OF LATITUDES IN GRID.
INTEGER, INTENT(IN) :: D_LON      !! LONGITUDE INCREMENT OF GRID [M_SEC].
INTEGER, INTENT(IN) :: D_LAT      !! LATITUDE INCREMENT OF GRID [M_SEC].
INTEGER, INTENT(IN) :: SOUTH      !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER, INTENT(IN) :: NORTH      !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER, INTENT(IN) :: WEST       !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER, INTENT(IN) :: EAST       !! EAST LONGITUDE OF GRID [M_SEC].
REAL,    INTENT(IN)  :: D_MAP(1:N_LON,1:N_LAT) !! WATER DEPTH [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL    :: LAND
LOGICAL :: R_GRID
LOGICAL :: L_OBSTRUCTION
LOGICAL :: L_INTERPOL
LOGICAL :: ERROR
character (len=len_coor) :: formtext

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLEAR INPUT TOPOGRAPHY.                                               !
!        -----------------------                                               !

NX_T = -1
NY_T = -1
IF (ALLOCATED(GRID_IN)) DEALLOCATE(GRID_IN)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COPY GRID DEFINITIONS.                                                !
!        ----------------------                                                !

DY_T    = D_LAT                     !! LATITUDE INCREMENT.
DX_T    = D_LON                     !! LONGITUDE INCREMENT.
SOUTH_T = SOUTH                     !! SOUTH LATITUDE OF GRID.
NORTH_T = NORTH                     !! NORTH LATITUDE OF GRID.
WEST_T  = WEST                      !! WEST LONGITUDE OF GRID.
EAST_T  = EAST                      !! EAST LONGITUDE OF GRID.
CALL ADJUST (WEST_T, EAST_T)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CHECK CONSISTENCY.                                                    !
!        ------------------                                                    !

CALL CHECK_GRID_DEFINITION (WEST_T, SOUTH_T, EAST_T, NORTH_T,                  &
&                           DX_T, DY_T, NX_T, NY_T, ERROR)

IF (ERROR) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *          FATAL ERROR IN SUB. SET_TOPOGRAPHY            *'
   WRITE (IU06,*) ' *          =================================             *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *        A TOPOGRAPY INPUT GRID COULD NOT BE DEFINED.    *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF

IF (N_LON.NE.NX_T .OR. N_LAT.NE.NY_T .OR.                                      &
&   N_LON.NE.SIZE(D_MAP,1) .OR. N_LAT.NE.SIZE(D_MAP,2)) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *         FATAL  ERROR IN SUB. SET_TOPOGRAPHY            *'
   WRITE (IU06,*) ' *         ===================================            *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' * GRID SPECIFICATIONS ARE NOT CONSISTENT                 *'
   WRITE (IU06,*) ' * USER PROVIDED GRID DEFINITIONS ARE:                    *'
   formtext = write_coor_text(WEST)
   WRITE (IU06,*) ' * LONGITUDE            WEST = ', formtext,                 &
&                                      ' = ',M_SEC_TO_DEG (WEST), 'deg'
   formtext = write_coor_text(EAST)
   WRITE (IU06,*) ' * LONGITUDE            EAST = ', formtext,                 &
&                                      ' = ',M_SEC_TO_DEG (EAST), 'deg'
   formtext = write_coor_text(D_LON)
   WRITE (IU06,*) ' * LONGITUDE INCREMENT D_LON = ', formtext,                 &
&                                      ' = ',M_SEC_TO_DEG (D_LON), 'deg'
   WRITE (IU06,*) ' * NO. OF LONGITUDES   N_LON = ', N_LON
   formtext = write_coor_text(NORTH)
   WRITE (IU06,*) ' * LATITUDE            NORTH = ', formtext,                 &
&                                      ' = ',M_SEC_TO_DEG (NORTH), 'deg'
   formtext = write_coor_text(SOUTH)
   WRITE (IU06,*) ' * LATITUDE            SOUTH = ', formtext,                 &
&                                      ' = ',M_SEC_TO_DEG (SOUTH), 'deg'
   formtext = write_coor_text(D_LAT)
   WRITE (IU06,*) ' * LATITUDE  INCREMENT D_LAT = ', formtext,                 &
&                                      ' = ',M_SEC_TO_DEG (D_LAT), 'deg'
   WRITE (IU06,*) ' * NO. OF LATITUDE     N_LAT = ', N_LAT
   WRITE (IU06,*) ' * COMPUTED GRID SIZES ARE:                               *'
   WRITE (IU06,*) ' * NO. OF LONGITUDES    NX_T = ', NX_T
   WRITE (IU06,*) ' * NO. OF LATITUDE      NY_T = ', NY_T
   WRITE (IU06,*) ' * DIMENSIONS OF DEPTH MAP ARRAY ARE :                    *'
   WRITE (IU06,*) ' * NO. OF LONGITUDES 1. DIMENSION = ', SIZE(D_MAP,1)
   WRITE (IU06,*) ' * NO. OF LATITUDE   2. DIMENSION = ', SIZE(D_MAP,2)
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   NX_T = -1
   NY_T = -1
   CALL ABORT1
END IF

PER_T = PERIODIC (WEST_T, EAST_T, DX_T, NX_T)   !! PERIODIC ?

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COPY DEPTH FIELD.                                                     !
!        -----------------                                                     !

ALLOCATE (GRID_IN(NX_T,NY_T))
GRID_IN = D_MAP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. INPUT MODULE DATA ARE DEFINED.                                        !
!        ------------------------------                                        !

SET_STATUS = (NX.GT.0 .AND. NY.GT.0 .AND. NX_T.GT.0 .AND. NY_T.GT.0)
IF (.NOT.SET_STATUS) THEN
   LAND   = LAND_LIMIT
   R_GRID = REDUCED_GRID
   L_INTERPOL = L_INTERPOL_T
   L_OBSTRUCTION = .FALSE.
   CALL SET_GRID_DEF (N_LON, N_LAT, D_LON, D_LAT, SOUTH, NORTH, WEST, EAST,    &
&                     LAND=LAND, R_GRID=R_GRID, L_INTERPOL=L_INTERPOL,         &
&                     L_OBSTRUCTION=L_OBSTRUCTION)
END IF

EQUAL_GRID = EQUAL_TO_M_GRID (NX_T, NY_T, DX_T, DY_T,                          &
&                             WEST_T, SOUTH_T, EAST_T, NORTH_T) 

END SUBROUTINE SET_TOPOGRAPHY_M

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_TOPOGRAPHY_C (N_LON, N_LAT, D_LON, D_LAT,                       &
&                            SOUTH, NORTH, WEST, EAST, D_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,                  INTENT(IN) :: N_LON !! NUMBER OF LONGITUDES IN GRID.
INTEGER,                  INTENT(IN) :: N_LAT !! NUMBER OF LATITUDES IN GRID.
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: D_LON !! LONGITUDE INCREMENT OF GRID [DMS].
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: D_LAT !! LATITUDE INCREMENT OF GRID [DMS].
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: SOUTH !! SOUTH LATITUDE OF GRID [DMS].
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: NORTH !! NORTH LATITUDE OF GRID [DMS].
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: WEST  !! WEST LONGITUDE OF GRID [DMS].
CHARACTER (LEN=LEN_COOR), INTENT(IN) :: EAST  !! EAST LONGITUDE OF GRID [DMS].
REAL,             INTENT(IN) :: D_MAP(1:N_LON,1:N_LAT) !! WATER DEPTH [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CONVERT TO M_SEC AND CALL SET_TOPOGRAPHY_M.                           !
!        -------------------------------------------                           !

CALL SET_TOPOGRAPHY_M (N_LON, N_LAT,                                           & 
&                      READ_COOR_TEXT(D_LON), READ_COOR_TEXT(D_LAT),           &
&                      READ_COOR_TEXT(SOUTH), READ_COOR_TEXT(NORTH),           &
&                      READ_COOR_TEXT(WEST),  READ_COOR_TEXT(EAST),            &
&                      D_MAP)

END SUBROUTINE SET_TOPOGRAPHY_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_CONST

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PREPARE_CONST - ROUTINE TO PREPARE WAM CONST MODULE.                       !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO COMPUTE ALL VARAIABLES IN WAM CONST MODULE FROM THE USER INPUT.     !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       DONE BY CALLS TO MANY SUBS.                                            !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. GENERATE MODEL GRID.                                                  !
!        --------------------                                                  !

CALL PREPARE_GRID
IF (ITEST.GT.1) WRITE (IU06,*) ' SUB. PREPARE_GRID DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE NEST INFORMATION.                                             !
!        -------------------------                                             !

CALL PREPARE_BOUNDARY_nest
IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. PREPARE_BOUNDARY DONE'

END SUBROUTINE PREPARE_CONST

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_PREPROC_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PRINT_PREPROC_STATUS - PRINT STATUS OF PREPROC.                            !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       MAKE A PRINTER OUTPUT OF THE PREPROC RESULTS, WHICH ARE SAVED IN       !
!       WAM_CONST_MODULE.                                                      !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTEGER             :: I
INTEGER             :: IGRID(NX,NY)  !! ARRAY FOR GRIDDED PRINT OUTPUT.
CHARACTER (LEN=1)   :: C_GRID(NX,NY) !! ARRAY FOR GRIDDED PRINT OUTPUT.
CHARACTER (LEN=100) :: TITL
CHARACTER (LEN=14)  :: ZERO = ' '
character (len=len_coor) :: formtext, ftext1, ftext2, ftext3, ftext4

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. FREQUENCY DIRECTION.                                                  !
!        --------------------                                                  !

CALL PRINT_FRE_DIR_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INPUT GRID DEFINITIONS.                                               !
!        -----------------------                                               !

WRITE (IU06,'(/,'' ----------------------------------------'')')
WRITE (IU06,'(  ''             INPUT TOPOGRAPHY '')')
WRITE (IU06,'(  '' ----------------------------------------'')')
WRITE(IU06,*) ' '

IF (NX_T.GT.0 .AND. NY_T.GT.0) THEN
   WRITE (IU06,'('' TOPOGRAPHY INPUT GRID: '')')
   WRITE (IU06,'('' NUMBER OF LONGITUDES IS      NX_T = '',I5)') NX_T
   formtext = write_coor_text(west_t)
   WRITE (IU06,'('' MOST WESTERN LONGITUDE IS  WEST_T = '',A)') formtext
   formtext = write_coor_text(east_t)
   WRITE (IU06,'('' MOST EASTERN LONGITUDE IS  EAST_T = '',A)') formtext
   formtext = write_coor_text(dx_t)
   WRITE (IU06,'('' LONGITUDE INCREMENT IS       DX_T = '',A)') formtext
   WRITE (IU06,'('' NUMBER OF LATITUDES IS       NY_T = '',I5)') NY_T
   formtext = write_coor_text(south_t)
   WRITE (IU06,'('' MOST SOUTHERN LATITUDE IS SOUTH_T = '',A)') formtext
   formtext = write_coor_text(north_t)
   WRITE (IU06,'('' MOST NORTHERN LATITUDE IS NORTH_T = '',A)') formtext
   formtext = write_coor_text(dy_t)
   WRITE (IU06,'('' LATITUDE INCREMENT IS        DY_T = '',A)') formtext
   IF (PER_T) THEN
      WRITE (IU06,*) 'THE GRID IS EAST-WEST PERIODIC'
   ELSE
      WRITE (IU06,*) 'THE GRID IS NOT EAST-WEST PERIODIC'
   END IF
   IF (EQUAL_GRID) THEN
      WRITE (IU06,*) 'GRID IS IDENTICAL TO MODEL GRID: NO SPACE INTERPOLATION'
   ELSE
      WRITE (IU06,*) 'GRID IS NOT IDENTICAL TO MODEL GRID: SPACE INTERPOLATION'
   END IF
ELSE
   WRITE (IU06,*) ' TOPOGRAPHY INPUT GRID IS NOT DEFINED'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. GRID CORRECTIONS.                                                     !
!        -----------------                                                     !

IF (N_CA.GT.0) THEN
   WRITE (IU06,'(/4X,'' AREAS TO BE CORRECTED IN OUTPUT GRID'',                &
&               /,4X,''  NO.   SOUTHERN LAT   NORTHERN LAT   WESTERN LONG '',  &
&                    ''  EASTERN LONG     DEPTH'')')
   DO I = 1,N_CA
     ftext1 = write_coor_text (south_ca(i))
     ftext2 = write_coor_text (north_ca(i))
     ftext3 = write_coor_text (west_ca(i))
     ftext4 = write_coor_text (east_ca(i))
     WRITE (IU06,'(4X,I5,4(2X,A),F10.2 )') I, ftext1, ftext2, ftext3,          &
&           ftext4, DEPTH_CA(I)
   END DO
ELSE
   WRITE (IU06,*) ' CORRECTION AREAS ARE NOT DEFINED'
END IF

WRITE (IU06,'('' DRY LAND POINTS UP TO  LAND_LIMIT = '',F10.2)') LAND_LIMIT
IF (L_INTERPOL_T) THEN
   WRITE (IU06,*) ' MODEL DEPTH IS NEAREST NEIGHBOUR DEPTH OF INPUT TOPOGRAPHY'
ELSE
   WRITE (IU06,*) ' MODEL DEPTH IS MEAN OF INPUT TOPOGRAPHY DEPTH IN GRID CELL'
END IF
IF (L_OBSTRUCTION_T) THEN
   WRITE (IU06,*) ' OBSTRUCTION FACTORS DUE TO SUBGRID FEATURES ARE COMPUTED'
ELSE
   WRITE (IU06,*) ' OBSTRUCTION FACTORS DUE TO SUBGRID FEATURES ARE NOT COMPUTED'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. BASIC MODEL GRID.                                                     !
!        -----------------                                                     !

CALL PRINT_GRID_STATUS

IF (NSEA.GT.0) THEN
   IF (ITEST.GT.3) THEN
      WRITE (IU06,*) ' '
      TITL = 'LAND - SEA MASK'
      C_GRID = 'L'
      WHERE (L_S_MASK) C_GRID = 'S'
      CALL PRINT_ARRAY (IU06, ZERO, TITL, C_GRID,                             &
&                       AMOWEP, AMOSOP, AMOEAP, AMONOP,NG_R=NLON_RG)
   END IF
   IF (ITEST.GT.9) THEN
      WRITE (IU06,*) ' '
      IGRID=UNPACK((/(I,I=1,NSEA)/), L_S_MASK,-999)
      IGRID=MOD(IGRID,1000)
      TITL = 'SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,       &
&                       AMONOP, NG_R=NLON_RG)
      IGRID=UNPACK(KLON(:,1), L_S_MASK,-999)
      IGRID=MOD(IGRID,1000)
      TITL = 'WEST NEIGHBOUR SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,       &
&                       AMONOP, NG_R=NLON_RG)
      WRITE (IU06,*) ' '
      IGRID=UNPACK(KLON(:,2), L_S_MASK,-999)
      IGRID=MOD(IGRID,1000)
      TITL = 'EAST NEIGHBOUR SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,       &
&                       AMONOP, NG_R=NLON_RG)
      WRITE (IU06,*) ' '
      IGRID=UNPACK(KLAT(:,1,1), L_S_MASK,-999)
      IGRID=MOD(IGRID,1000)
      TITL = 'NEAREST SOUTH NEIGHBOUR SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,       &
&                       AMONOP, NG_R=NLON_RG)
      WRITE (IU06,*) ' '
      IGRID=UNPACK(KLAT(:,1,2), L_S_MASK,-999)
      IGRID=MOD(IGRID,1000)
      TITL = 'SECOND NEAREST SOUTH NEIGHBOUR SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,       &
&                       AMONOP, NG_R=NLON_RG)
      WRITE (IU06,*) ' '
      IGRID=UNPACK(KLAT(:,2,1), L_S_MASK,-999)
      IGRID=MOD(IGRID,1000)
      TITL = 'NEAREST NORTH NEIGHBOUR SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,       &
&                       AMONOP, NG_R=NLON_RG)
      WRITE (IU06,*) ' '
      IGRID=UNPACK(KLAT(:,2,2), L_S_MASK,-999)
      IGRID=MOD(IGRID,1000)
      TITL = 'SECOND NEAREST NORTH NEIGHBOUR SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,       &
&                       AMONOP, NG_R=NLON_RG)
      WRITE (IU06,*) ' '
      IGRID=NINT(UNPACK(WLAT(:,1), L_S_MASK,-999.)*100.)
      TITL = 'WEIGHT OF CLOSEST SOUTH NEIGHBOUR POINT (*100)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,       &
&                       AMONOP, NG_R=NLON_RG)
      WRITE (IU06,*) ' '
      IGRID=NINT(UNPACK(WLAT(:,2), L_S_MASK,-999.)*100.)
      TITL = 'WEIGHT OF CLOSEST NORTH NEIGHBOUR POINT (*100)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,       &
&                       AMONOP, NG_R=NLON_RG)
      WRITE (IU06,*) ' '

      IF (L_OBSTRUCTION_T) THEN
         IGRID=NINT(UNPACK(OBSLAT (1:NSEA,1,1), L_S_MASK,-999.)*100.)
         TITL = 'OBSTRUCTION (*100) FOR S-N ADVECTION FOR FREQUENCY 1'
         CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,    &
&                          AMONOP, NG_R=NLON_RG)
         WRITE (IU06,*) ' '
         IGRID=NINT(UNPACK(OBSLAT (1:NSEA,2,1), L_S_MASK,-999.)*100.)
         TITL = 'OBSTRUCTION (*100) FOR N-S ADVECTION FOR FREQUENCY 1'
         CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,    &
&                          AMONOP, NG_R=NLON_RG)
         WRITE (IU06,*) ' '
         IGRID=NINT(UNPACK(OBSLON (1:NSEA,1,1), L_S_MASK,-999.)*100.)
         TITL = 'OBSTRUCTION (*100) FOR W-E ADVECTION FOR FREQUENCY 1'
         CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,    &
&                          AMONOP, NG_R=NLON_RG)
         WRITE (IU06,*) ' '
         IGRID=NINT(UNPACK(OBSLON (1:NSEA,2,1), L_S_MASK,-999.)*100.)
         TITL = 'OBSTRUCTION (*100) FOR E-W ADVECTION FOR FREQUENCY 1'
         CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP,    &
&                          AMONOP, NG_R=NLON_RG)
         WRITE (IU06,*) ' '
      END IF
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. NEST INFORMATION.                                                     !
!        -----------------                                                     !

CALL PRINT_NEST_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. GENERAL INFORMATION.                                                  !
!        --------------------                                                  !

CALL PRINT_GENERAL_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. TABLES INFORMATION.                                                   !
!        -------------------                                                   !

CALL PRINT_TABLES_STATUS

END SUBROUTINE PRINT_PREPROC_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WRITE_PREPROC_FILE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   WRITE_PREPROC_FILE - ROUTINE TO WRITE PREPROC OUTPUT TO FILE               !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO WRITE OUT THE COMPUTED CONSTANTS WHICH ARE STORED IN MODULE         !
!       WAM_CONST_MODULE.                                                      !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       UNFORMATTED WRITE AS SPECIFIED TO UNIT = IU07.                         !
!       FILENAME IS 'FILE07' AS DEFINED IN THE USER INPUT                      !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER      :: LEN, I

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. OPEN FILES.                                                           !
!        -----------                                                           !

LEN = LEN_TRIM(FILE07)
OPEN (UNIT=IU07, FILE=FILE07(1:LEN), FORM='UNFORMATTED', STATUS='UNKNOWN')

WRITE(IU07) HEADER

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. WRITE COARSE GRID BOUNDARY OUTPUT INFORMATION.                        !
!        ----------------------------------------------                        !

WRITE (IU07) N_NEST, MAX_NEST
DO I=1,N_NEST
   WRITE (IU07) NBOUNC(I), N_NAME(I), n_code(i)
   IF (NBOUNC(I).GT.0) THEN
      WRITE(IU07) IJARC(1:NBOUNC(I),I)
      WRITE(IU07) XDELLO, XDELLA, N_SOUTH(I), N_NORTH(I), N_EAST(I), N_WEST(I),&
&              BLNGC(1:NBOUNC(I),I), BLATC(1:NBOUNC(I),I), N_ZDEL(1:NBOUNC(I),I)
   END IF
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. WRITE FINE GRID BOUNDARY INPUT INFORMATION.                           !
!        -------------------------------------------                           !

WRITE(IU07) NBOUNF, NBINP, C_NAME
IF (NBOUNF.GT.0) THEN
   WRITE(IU07) BLNGF(1:NBOUNF), BLATF(1:NBOUNF), IJARF(1:NBOUNF),              &
&              IBFL(1:NBOUNF), IBFR(1:NBOUNF), BFW(1:NBOUNF)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. WRITE FREQUENCY DIRECTION GRID.                                       !
!        ------------------------------                                        !

WRITE (IU07) ML, KL
WRITE (IU07) FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH, INV_LOG_CO,     &
&            DF, DF_FR, DF_FR2, DFIM, DFIMOFR, DFIM_FR, DFIM_FR2, FR5, FRM5,   &
&            RHOWG_DFIM,                                                       &
&            FMIN, MO_TAIL, MM1_TAIL, MP1_TAIL, MP2_TAIL
WRITE (IU07) MPM, KPM, JXO, JYO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. WRITE GRID INFORMATION.                                               !
!        -----------------------                                               !

WRITE (IU07) NX, NY, NSEA, IPER, ONE_POINT, REDUCED_GRID, L_OBSTRUCTION_T
WRITE (IU07) NLON_RG
WRITE (IU07) DELPHI, DELLAM, SINPH, COSPH, AMOWEP, AMOSOP, AMOEAP, AMONOP,     &
&            XDELLA, XDELLO, ZDELLO
WRITE (IU07) IXLG, KXLT, L_S_MASK
WRITE (IU07) KLAT, KLON, WLAT, DEPTH_B
IF (L_OBSTRUCTION_T) THEN
   WRITE (IU07) OBSLAT, OBSLON
END IF
! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. WRITE TABLES.                                                         !
!        -------------                                                         !

WRITE (IU07) NDEPTH, DEPTHA, DEPTHD, DEPTHE
WRITE (IU07) FLMINFR, TCGOND, TFAK, TSIHKD, TFAC_ST, T_TAIL
WRITE (IU07) DELU

! ---------------------------------------------------------------------------- !
!                                                                              !
!    10. CLOSE FILES.                                                          !
!        ------------                                                          !

CLOSE (UNIT=IU07, STATUS="KEEP")

END SUBROUTINE WRITE_PREPROC_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVAT MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE ADJUST_DEPTH

! ---------------------------------------------------------------------------- !
!                                                                              !
!   ADJUST_DEPTH - ADJUST DEPTH IN AREAS.                                      !
!                                                                              !
!     S. HASSELMANN     MPIFM           1/6/86.                                !
!                                                                              !
!     MODIFIED BY       H. GUNTHER      1/4/90  -  REARANGEMENT OF CODE.       !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO CHANGE THE MODEL DEPTH IN USER SPECIFIED AREAS.                     !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: I, K, L
INTEGER :: XLAT, XLON

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. MANUAL ADJUSTMENT OF TOPOGRAPHY.                                      !
!        --------------------------------                                      !

IF (N_CA.GT.0) THEN
   DO K = 1,NY
      XLAT = AMOSOP+(K-1)*XDELLA
      DO I = 1,NLON_RG(K)
         XLON = AMOWEP+(I-1)*ZDELLO(K)
         IF (XLON.GE.M_S_PER) XLON = XLON-M_S_PER
         DO L = 1,N_CA
            IF (XLON.LT.WEST_CA(L)) XLON = XLON+M_S_PER
            IF (XLON.GT.EAST_CA(L)) XLON = XLON-M_S_PER
            IF (XLON.GE.WEST_CA(L) .AND. XLAT.GE.SOUTH_CA(L) .AND.             &
&               XLON.LE.EAST_CA(L) .AND. XLAT.LE.NORTH_CA(L))                  &
&                       GRD(I,K) = DEPTH_CA(L)
         END DO
      END DO
   END DO
END IF

END SUBROUTINE ADJUST_DEPTH

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CLOSEST_GP_LAT (I, K, IP, XH, KH, IP1, IP2, WH)

! ---------------------------------------------------------------------------- !
!                                                                              !
!  CLOSEST_GP_LAT - ROUTINE TO FIND THE CLOSEST AND SECOND CLOSEST GRID POINT  !
!                   ON LATITUDE NORTH OR SOUTH AND THE INTEPOLATION WEIGHT.    !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     J. BIDLOT            ECMWF       APRIL 2000: add second closest          !
!                                                  grid points.                !
!     J. BIDLOT            ECMWF       CLOSEST AND SECOND CLOSEST GRID         !
!                                      POINT FOR THE ROTATED CELL.             !
!     H.GUNTHER            HZG         CODE TAKEN FROM SUB. UBUF               !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO ARRANGE NEIGHBOUR GRID POINT INDICES FOR A GIVEN SEA POINT ON       !
!       LATITUDE NORTH OR SOUTH AND THE INTEPOLATION WEIGHT.                   !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE INDICES OF THE NEXT POINTS ON LAT. AND LONG. ARE                   !
!       COMPUTED. ZERO INDICATES A LAND POINT IS NEIGHBOUR.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)  :: I      !! LONGITUDE INDEX OF GRID POINT.
INTEGER, INTENT(IN)  :: K      !! LATITUTE INDEX OF GRID POINT.
INTEGER, INTENT(IN)  :: IP     !! SEA POINT NUMBER OF GRID POINT.
REAL,    INTENT(IN)  :: XH     !! LONGITUDE OF GRID POINT.
INTEGER, INTENT(IN)  :: KH     !! LATITUDE INDEX NORTH (K+1) OR SOUTH (K-1).
INTEGER, INTENT(OUT) :: IP1    !! NUMBER OF CLOSEST SEA POINT TO IP
!! ON LATITUDE KH.
INTEGER, INTENT(OUT) :: IP2    !! NUMBER OF SECOND CLOSEST SEA POINT TO IP
!! ON LATITUDE KH.
REAL,    INTENT(OUT) :: WH     !! INTERPOLATION WEIGHT FOR CLOSEST POINT.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: IC, ICA, ICE, ICS, IPH, IPH1, IPH2
REAL     :: D3, D5, XP, D4, D6, Z

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIAL.                                                              !
!        --------                                                              !

WH  = 1.               !! INTERPOLATION WEIGHT
IP1 = 0                !! NUMBER OF CLOSEST SEA POINT
IP2 = 0                !! NUMBER OF SECOND CLOSEST SEA POINT
IF (KH.LT.K) THEN      !! SOUTH LATITUDE
   IF (KH.LT.1) RETURN
   ICA = IP            !! SEA POINT COUNTER START
   ICE = 1             !! SEA POINT COUNTER END
   ICS = -1            !! SEA POINT COUNTER STEP
ELSE                   !! NORTH LATITUDE
   IF (KH.GT.NY) RETURN
   ICA = IP            !! SEA POINT COUNTER START
   ICE = NSEA          !! SEA POINT COUNTER END
   ICS = +1            !! SEA POINT COUNTER STEP
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CLOSEST AND SECOND CLOSEST GRID POINT.                                !
!        --------------------------------------                                !

IF (NLON_RG(KH).EQ.NLON_RG(K)) THEN
   IPH1 = I
   IPH2 = I
ELSE
   IPH1 = FLOOR(XH)
   IPH2 = CEILING(XH)
   IF (IPH1.EQ.IPH2) THEN
      IPH1 = IPH1+1
      IPH2 = IPH1
   ELSE
      IF (IPH2.EQ.NINT(XH)) THEN
         IPH = IPH1
         IPH1 = IPH2+1
         IPH2 = IPH+1
      ELSE
         IPH1 = IPH1+1
         IPH2 = IPH2+1
      END IF
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. INTERPOLATION WEIGHT FOR CLOSEST POINT.                               !
!        ---------------------------------------                               !

IF (IPH1.NE.IPH2) THEN
   Z = REAL(ZDELLO(KH))/REAL(ZDELLO(K))
   D3 = XH*Z - 0.5
   D5 = XH*Z + 0.5
   XP = REAL(IPH1-1)*Z
   D4 = XP-0.5*Z
   D6 = XP+0.5*Z
   IF (D4.LT.D3 .OR. D6.GT.D5) THEN
      WH = MIN(1.,(MIN(D5,D6)-MAX(D3, D4)))
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. SEAPOINT NUMMER OF FOR CLOSEST POINT.                                 !
!        -------------------------------------                                 !

IF (IPER) THEN
   IF (IPH1.LT.1) THEN
      IPH1 = IPH1 + NLON_RG(KH)
   END IF
   IF (IPH1.GT.NLON_RG(KH)) THEN
      IPH1 = IPH1 - NLON_RG(KH)
   END IF
END IF
IF (IPH1.GE.1.AND.IPH1.LE.NLON_RG(KH)) THEN
   IF (L_S_MASK(IPH1,KH)) THEN
      DO IC = ICA,ICE,ICS
         IF (IXLG(IC).EQ.IPH1 .AND.KXLT(IC).EQ.KH) EXIT
      END DO
      IP1 = IC
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. SEAPOINT NUMMER OF FOR SECOND CLOSEST POINT.                          !
!        --------------------------------------------                          !

IF (IPH1.NE.IPH2) THEN
   IF (IPER) THEN
      IF (IPH2.LT.1) THEN
         IPH2 = IPH2 + NLON_RG(KH)
      END IF
      IF (IPH2.GT.NLON_RG(KH)) THEN
         IPH2 = IPH2 - NLON_RG(KH)
      END IF
   END IF
   IF (IPH2.GE.1.AND.IPH2.LE.NLON_RG(KH)) THEN
      IF (L_S_MASK(IPH2,KH)) THEN
         DO IC = ICA,ICE,ICS
            IF (IXLG(IC).EQ.IPH2 .AND.KXLT(IC).EQ.KH) EXIT
         END DO
         IP2 = IC
      END IF
   END IF
ELSE
   IP2 = IP1
END IF

END SUBROUTINE CLOSEST_GP_LAT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !


SUBROUTINE COUNT_POINTS_E_W (NX_TL, NX_TR, NY_TB, NY_TT,                       &
&                        IBLOCKDPT, ITHRSHOLD, IEXCLTHRSHOLD, GRD,             &
&                        PERCENTSHALLOW, PSHALLOWTRHS, PERCENTLAND, PLANDTRHS, &
&                        NOBSTRCT, NTOTPTS)
INTEGER, INTENT(IN) :: NX_TL, NX_TR, NY_TB, NY_TT
INTEGER, INTENT(IN) :: IBLOCKDPT, ITHRSHOLD, IEXCLTHRSHOLD
REAL,    INTENT(IN) ::  GRD
REAL,    INTENT(IN) ::  PERCENTSHALLOW, PSHALLOWTRHS, PERCENTLAND, PLANDTRHS

INTEGER, INTENT(OUT) :: NOBSTRCT, NTOTPTS

INTEGER :: II, I, J, NIOBSLON, NX_TH, IREINF, NBLOCKLAND

LOGICAL :: LLAND, LREALLAND, LNSW, L1ST

IREINF= 1
NBLOCKLAND = 0
NOBSTRCT = 0
NX_TH = NX_TR
IF (NX_TR.LT.NX_TL) NX_TH = NX_TR + NX_T-1
NTOTPTS = (NY_TT-NY_TB+1)*(NX_TH-NX_TL+1)

DO J = NY_TB, NY_TT
   NIOBSLON=0
   LLAND=.FALSE.
   LREALLAND=.FALSE.

   DO II = NX_TL, NX_TH
      I = II
      IF (II.GT.NX_T) I = II-NX_T+1
      IF (GRID_IN(I,J).LE.IBLOCKDPT ) THEN
         IF (GRID_IN(I,J).LE.0 ) LREALLAND = .TRUE.
         LLAND =.TRUE.
         NIOBSLON = NIOBSLON+1

!         IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!         IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA (SEE BELOW)
!         LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!         -----------------------------------------------------------

      ELSEIF (GRID_IN(I,J) .LE. ITHRSHOLD .AND.  &
&                GRD.GT.IEXCLTHRSHOLD) THEN

!         IF SEA ABOVE THE THRESHOLD THEN ONLY THAT GRID POINTS BLOCKS
!         ------------------------------------------------------------
         NIOBSLON = NIOBSLON + 1
      ENDIF
   ENDDO

   IF (LLAND) THEN
      IF (LREALLAND) THEN
         LNSW=.TRUE.
         IF (GRID_IN(NX_TL,J).LE.IBLOCKDPT) THEN
            L1ST=.TRUE.
         ELSE
            L1ST=.FALSE.
         ENDIF
         DO II = NX_TL+1,NX_TH
            I = II
            IF (II.GT.NX_T) I = II-NX_T+1
            IF ( ((GRID_IN(I,J).LE.IBLOCKDPT) .NEQV.L1ST)  &
&                                      .AND. LNSW ) THEN
               LNSW=.FALSE.
            ENDIF
            IF (((GRID_IN(I,J).LE.IBLOCKDPT).EQV.L1ST)  &
&                                    .AND. .NOT. LNSW ) THEN

!                           LAND IS BLOCKING
               NIOBSLON=IREINF*(NX_TH-NX_TL+1)
               NBLOCKLAND=NBLOCKLAND+1
               EXIT
            ENDIF
         ENDDO
         IF (LNSW) NIOBSLON = NX_TH-NX_TL+1
      ELSE
         IF (PERCENTSHALLOW.GT.PSHALLOWTRHS) THEN
!                         mostly shallow, do not enhance obstruction
            NIOBSLON = NX_TH-NX_TL+1
         ELSEIF (PERCENTLAND.LT.PLANDTRHS) THEN
!                         does not contain too much land
            NIOBSLON = IREINF*(NX_TH-NX_TL+1)
            NBLOCKLAND = NBLOCKLAND+1
         ELSE
            NIOBSLON = 0
         ENDIF
      ENDIF     !! LREALLAND
   ENDIF     !! LLAND

   NOBSTRCT = NOBSTRCT + NIOBSLON
ENDDO

END SUBROUTINE COUNT_POINTS_E_W

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE COUNT_POINTS_N_S (NX_TL, NX_TR, NY_TB, NY_TT,                       &
&                        IBLOCKDPT, ITHRSHOLD, IEXCLTHRSHOLD, GRD,        &
&                        PERCENTSHALLOW, PSHALLOWTRHS, PERCENTLAND, PLANDTRHS, &
&                        NOBSTRCT, NTOTPTS)
INTEGER, INTENT(IN)  :: NX_TL, NX_TR, NY_TB, NY_TT
INTEGER, INTENT(IN)  :: IBLOCKDPT, ITHRSHOLD, IEXCLTHRSHOLD
REAL,    INTENT(IN) ::  GRD
REAL,    INTENT(IN) ::  PERCENTSHALLOW, PSHALLOWTRHS, PERCENTLAND, PLANDTRHS

INTEGER, INTENT(OUT) :: NOBSTRCT, NTOTPTS

INTEGER :: II, I, J, NIOBSLON, NX_TH, IREINF, NBLOCKLAND

LOGICAL :: LLAND, LREALLAND, LNSW, L1ST

IREINF= 1
NBLOCKLAND = 0
NOBSTRCT = 0
NX_TH = NX_TR
IF (NX_TR.LT.NX_TL) NX_TH = NX_TR + NX_T-1
NTOTPTS = (NY_TT-NY_TB+1)*(NX_TH-NX_TL+1)

DO II = NX_TL, NX_TH
   I = II
   IF (II.GT.NX_T) I = II-NX_T+1
   NIOBSLON=0
   LLAND=.FALSE.
   LREALLAND=.FALSE.

   DO J = NY_TB, NY_TT
      IF (GRID_IN(I,J).LE.IBLOCKDPT ) THEN
         IF (GRID_IN(I,J).LE.0 ) LREALLAND = .TRUE.
         LLAND =.TRUE.
         NIOBSLON = NIOBSLON+1

!         IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!         IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA (SEE BELOW)
!         LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!         -----------------------------------------------------------

      ELSEIF (GRID_IN(I,J) .LE. ITHRSHOLD .AND.  &
&                GRD.GT.IEXCLTHRSHOLD) THEN

!         IF SEA ABOVE THE THRESHOLD THEN ONLY THAT GRID POINTS BLOCKS
!         ------------------------------------------------------------

         NIOBSLON = NIOBSLON + 1
      ENDIF
   ENDDO

   IF (LLAND) THEN
      IF (LREALLAND) THEN
         LNSW=.TRUE.
         IF (GRID_IN(I,NY_TT).LE.IBLOCKDPT) THEN
            L1ST=.TRUE.
         ELSE
            L1ST=.FALSE.
         ENDIF
         DO J = NY_TT-1,NY_TB,-1
            IF ( ((GRID_IN(I,J).LE.IBLOCKDPT) .NEQV.L1ST)  &
&                                      .AND. LNSW ) THEN
               LNSW=.FALSE.
            ENDIF
            IF (((GRID_IN(I,J).LE.IBLOCKDPT).EQV.L1ST)  &
&                                    .AND. .NOT. LNSW ) THEN

!                           LAND IS BLOCKING
               NIOBSLON=IREINF*(NY_TT-NY_TB+1)
               NBLOCKLAND=NBLOCKLAND+1
               EXIT
            ENDIF
         ENDDO
         IF (LNSW) NIOBSLON = NY_TT-NY_TB+1
      ELSE
         IF (PERCENTSHALLOW.GT.PSHALLOWTRHS) THEN
!                         mostly shallow, do not enhance obstruction
            NIOBSLON=NY_TT-NY_TB+1
         ELSEIF (PERCENTLAND.LT.PLANDTRHS) THEN
!                         does not contain too much land
            NIOBSLON = IREINF*(NY_TT-NY_TB+1)
            NBLOCKLAND = NBLOCKLAND+1
         ELSE
            NIOBSLON = 0
         ENDIF
      ENDIF     !! LREALLAND
   ENDIF     !! LLAND

   NOBSTRCT = NOBSTRCT + NIOBSLON
ENDDO

END SUBROUTINE COUNT_POINTS_N_S

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CREATE_OBSTRUCTIONS

! ---------------------------------------------------------------------------- !
!

INTEGER :: M, IS, K, IX
REAL    :: OMEGA, XKDMAX, DEPTH, XX
INTEGER :: KT, KB, NY_TB, NY_TT, NX_TL, NX_TR, IREINF
INTEGER :: XLATT, XLATB, XLONL, XLONR
INTEGER :: STEPT, STEPB, STEP_LON, STEP_LAT
INTEGER :: NOBSTRCT, NTOTPTS

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ITHRSHOLD     !! Threshold depth for blocking
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IEXCLTHRSHOLD !! exceeds threshold (no blocking)
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IBLOCKDPT     !! Blocking depth
REAL, ALLOCATABLE, DIMENSION(:,:) :: HOBSLAT
REAL, ALLOCATABLE, DIMENSION(:,:) :: HOBSLON

REAL :: PLANDTRHS
REAL :: PSHALLOWTRHS
REAL (KIND=KIND_D) :: X_HELP

! ---------------------------------------------------------------------------- !
!

X_HELP = 2*DY_T
IF (XDELLA .LE. X_HELP ) THEN !! COMPUTE OBSTRUCTIONS ONLY WHEN IT IS MEANINGFUL
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +     WARNING ERROR SUB. CREATE_OBSTRUCTIONS.      +'
   WRITE (IU06,*) ' +     =======================================      +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + THE REQUESTED RESOLUTION IS SMALL ENOUGH WITH    +'
   WRITE (IU06,*) ' + RESPECT TO THE INPUT BATHYMETRY DATA.            +'
   WRITE (IU06,*) ' + NO OBSTRUCTIONS WILL BE COMPUTED.                +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +               MODEL CONTINUES                    +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
   L_OBSTRUCTION_T = .FALSE.
   RETURN
END IF

PLANDTRHS = 0.3
PSHALLOWTRHS = 0.8

!     IREINF IS USED TO REINFORCE LAND OBSTRUCTIONS FOR SMALL GRID SPACING.

X_HELP = 0.5
IF (XDELLA.LE.DEG_TO_M_SEC(X_HELP)) THEN
   IREINF = 2
ELSE
   IREINF = 1
ENDIF

WRITE (*,*) ' IREINF', IREINF

ALLOCATE(HOBSLAT(NX,NY))
ALLOCATE(HOBSLON(NX,NY))

HOBSLAT(:,:) = 1.
HOBSLON(:,:) = 1.

ALLOCATE(ITHRSHOLD(NX,NY))
ALLOCATE(IBLOCKDPT(NX,NY))
ALLOCATE(IEXCLTHRSHOLD(NX,NY))

FREQ: DO M = 1,ML

! ---------------------------------------------------------------------------- !
!
!       COMPUTE THE THRESHOLD AT WHICH THE WAVES ARE OBSTRUCTED BY THE BOTTOM,
!       EXCEPT IF WAM DEPTH OF THE SAME ORDER OF MAGITUDE.
!       ALSO COMPUTE THE DEPTH THAT IS CONSIDERED TO BE FULLY BLOCKING
!       AS IF IT WAS LAND.

   OMEGA  = ZPI*FR(M)
   XKDMAX = 2.
   DO K=1,NY
      DO IX=1,NLON_RG(K)
         IF (GRD(IX,K).GT.0.) THEN
            DEPTH = GRD(IX,K)
            XX = XKDMAX/AKI(OMEGA,DEPTH)
            ITHRSHOLD(IX,K) = NINT(XX)
            IEXCLTHRSHOLD(IX,K) = MIN(10*ITHRSHOLD(IX,K),998)
            IBLOCKDPT(IX,K) = INT(0.05*XX)
         END IF
     END DO
   END DO

! ---------------------------------------------------------------------------- !
!
!       NORTH-SOUTH OBSTRUCTIONS
!       -----------------------
!       IS=1 is for the south-north advection
!       IS=2 is for the north-south advection

   DO IS = 1,2
      DO K = 1,NY
         IF (IS.EQ.1) THEN
               KT = K
               KB = K-1
               STEPT = -DY_T
               STEPB = 0
               STEP_LON = ZDELLO(K)/2
         ELSE
               KT = K+1
               KB = K
               STEPT = 0
               STEPB = DY_T
               STEP_LON = ZDELLO(K)/2
         END IF
         XLATT = XLAT(KT) + STEPT
         XLATB = XLAT(KB) + STEPB
         NY_TT = NINT(REAL(XLATT-SOUTH_T)/REAL(DY_T)) + 1
         NY_TT = MAX(1,MIN(NY_TT,NY_T))
         NY_TB = NINT(REAL(XLATB-SOUTH_T)/REAL(DY_T)) + 1
         NY_TB = MAX(1,MIN(NY_TB,NY_T))

         DO IX = 1,NLON_RG(K)
            IF (GRD(IX,K).LE.0) CYCLE
            XLONL = AMOWEP + (IX-1)*ZDELLO(K) - STEP_LON
            IF (XLONL.GT.EAST_T) XLONL = XLONL - M_S_PER
            XLONR = AMOWEP + (IX-1)*ZDELLO(K) + STEP_LON
            IF (XLONR.GT.EAST_T) XLONR = XLONR - M_S_PER
            NX_TL = NINT(REAL(XLONL - WEST_T)/REAL(DX_T)) + 1
            NX_TR = NINT(REAL(XLONR - WEST_T)/REAL(DX_T)) + 1
            CALL COUNT_POINTS_N_S (NX_TL, NX_TR, NY_TB, NY_TT,                 &
&                                  IBLOCKDPT(IX,K), ITHRSHOLD(IX,K),           &
&                                  IEXCLTHRSHOLD(IX,K), GRD(IX,K),             &
&                                  PERCENTSHALLOW(IX,K), PSHALLOWTRHS,         &
&                                  PERCENTLAND(IX,K), PLANDTRHS,               &
&                                  NOBSTRCT, NTOTPTS)

            HOBSLAT(IX,K) = (1.-FLOAT(NOBSTRCT)/NTOTPTS)

         ENDDO           !! LONGITUDES
      ENDDO              !! LATTITUDES
      OBSLAT (1:NSEA,IS,M) = PACK (HOBSLAT, L_S_MASK)
   ENDDO                 !! NORTH-SOUTH OBSTRUCTIONS


! ---------------------------------------------------------------------------- !
!
!       EAST-WEST OBSTRUCTIONS
!       -----------------------
!       IS=1 is for the west-east advection
!       IS=2 is for the east-west advection

   DO IS=1,2
      DO K=1,NY
         STEP_LAT = NINT(REAL(XDELLA)/2.)
         XLATT = XLAT(K) + STEP_LAT
         XLATB = XLAT(K) - STEP_LAT
         NY_TT = NINT(REAL(XLATT-SOUTH_T)/REAL(DY_T)) + 1
         NY_TT = MAX(1,MIN(NY_TT,NY_T))
         NY_TB = NINT(REAL(XLATB-SOUTH_T)/REAL(DY_T)) + 1
         NY_TB = MAX(1,MIN(NY_TB,NY_T))

         DO IX=1,NLON_RG(K)
            IF (GRD(IX,K).LE.0) CYCLE
            IF(IS.EQ.1) THEN
                  XLONL = AMOWEP + (IX-2)*ZDELLO(K)
                  XLONR = AMOWEP + (IX-1)*ZDELLO(K) - DX_T
            ELSE
                  XLONL = AMOWEP + (IX-1)*ZDELLO(K) + DX_T
                  XLONR = AMOWEP + IX*ZDELLO(K)
            ENDIF
            IF (XLONL.GT.EAST_T) XLONL = XLONL - M_S_PER
            IF (XLONR.GT.EAST_T) XLONR = XLONR - M_S_PER

            NX_TL = NINT(REAL(XLONL - WEST_T)/REAL(DX_T)) + 1
            NX_TR = NINT(REAL(XLONR - WEST_T)/REAL(DX_T)) + 1

            CALL COUNT_POINTS_E_W (NX_TL, NX_TR, NY_TB, NY_TT,                 &
&                                  IBLOCKDPT(IX,K), ITHRSHOLD(IX,K),           &
&                                  IEXCLTHRSHOLD(IX,K), GRD(IX,K),        &
&                                  PERCENTSHALLOW(IX,K), PSHALLOWTRHS,         &
&                                  PERCENTLAND(IX,K), PLANDTRHS, &
&                                  NOBSTRCT, NTOTPTS)

            HOBSLON(IX,K)= (1.-FLOAT(NOBSTRCT)/NTOTPTS)

         ENDDO           !! LONGITUDES
      ENDDO              !! LATTITUDES
      OBSLON(1:NSEA,IS,M) = PACK (HOBSLON, L_S_MASK)
   ENDDO                 !! EAST-WEST OBSTRUCTIONS

END DO FREQ

END SUBROUTINE CREATE_OBSTRUCTIONS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MEAN_DEPTH

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MEAN_DEPTH - Average depth.                                                !
!                                                                              !
!     JEAN BIDLOT JUNE 2007                                                    !
!                                                                              !
!     MODIFIED BY       H. GUNTHER      1/4/2015  -  REARANGEMENT OF CODE.     !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO CHANGE THE MODEL DEPTH IN USER SPECIFIED AREAS.                     !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: K, J, JJ, IX, I, II, IK

INTEGER :: NLANDCENTREPM      !! HALF CENTER SIZE OF A GRID BOX
INTEGER :: NLANDCENTREMAX     !! Number of points in HALF CENTER SIZE OF A GRID BOX
INTEGER :: NLANDCENTRE        !! Number of Land points in the box centre.
INTEGER :: NIM, NIP           !! LONGITUDE LEFT AND RIGHT AVARAGE LIMITS
INTEGER :: NJM, NJP           !! LATITUDE BOTTOM AND TOP AVARAGE LIMITS
INTEGER :: NSEA               !! counts sea points
REAL    :: SEA                !! cumulates depth
INTEGER :: NLAND              !! counts land points
REAL    :: XLAND              !! cumulates land depth
INTEGER :: NSEASH             !! counts shallow points (depth less than 500m)
REAL    :: SEASH              !! cumulates depth at shallow points
REAL    :: XLON

! ---------------------------------------------------------------------------- !

ALLOCATE (PERCENTLAND(NX,NY))
ALLOCATE (PERCENTSHALLOW(NX,NY))

NLANDCENTREPM  = (NINT(0.2*REAL(XDELLA)/REAL(DY_T))-1)/2 !! HALF CENTER SIZE OF A GRID BOX
NLANDCENTREPM  = MAX(NLANDCENTREPM,1)
NLANDCENTREMAX = (2*NLANDCENTREPM+1)**2

NJM = INT(0.5*REAL(XDELLA)/REAL(DY_T))            !! LATITUDE BOTTOM AND TOP AVARAGE LIMITS
NJP = NJM
ALLOCATE (ALON(NX_T))
ALLOCATE (ALAT(NY_T))
ALLOCATE (XLAT(0:NY+1))

DO I = 1, NX_T
   ALON(I) = WEST_T +(I-1)*DX_T
END DO
DO J = 1, NY_T
   ALAT(J)  = SOUTH_T-(J-1)*DY_T
END DO
DO K = 0,NY+1
   XLAT(K) = AMOSOP + (K-1)*XDELLA
END DO

!        WE ASSUME THAT WAMGRID IS ALWAYS WITHIN ETOPO2

DO K=1,NY

!        DETERMINE CLOSEST ETO2 J INDEX TO WAM POINT

   DO J = NY_T-1,1,-1
      IF (ALAT(J+1).LT.XLAT(K) .AND. XLAT(K).LE.ALAT(J) ) EXIT
   ENDDO

   J = NINT(REAL(XLAT(K)-SOUTH_T)/REAL(DY_T)) + 1

   IF (J.GT.NY_T .OR. J.LT.1) THEN
      WRITE(*,*) 'PROBLEM WITH J !!!'
      CALL ABORT1
   ENDIF

   DO IX = 1,NLON_RG(K)
      XLON = AMOWEP + (IX-1)*ZDELLO(K)
      IF (XLON.GE.EAST_T) XLON = XLON - M_S_PER     ! ETOPO2 STARTS AT -180

!          DETERMINE CLOSEST ETOPO2 I INDEX TO WAM POINT

      I = NINT(REAL(XLON-WEST_T)/REAL(DX_T)) + 1

!          AVERAGE OVER LAND AND SEA SEPARATELY AROUND POINT I,J

      NIM = INT(0.5*REAL(ZDELLO(K))/REAL(DX_T))
      NIP = NIM

      NSEA   = 0     ! counts sea points
      SEA    = 0.    ! cumulates depth
      NLAND  = 0     ! counts land points
      XLAND  = 0.    ! cumulates land depth
      NSEASH = 0     ! counts shallow points (depth less than 500m)
      SEASH  = 0.    ! cumulates depth at shallow points

      DO JJ = J-NJM, J+NJP
         IF (JJ.GE.1 .AND. JJ.LE.NY_T) THEN
            DO II = I-NIM,I+NIP
               IK = II
               IF (II.LT.1)    IK = NX_T-1 + II
               IF (II.GT.NX_T) IK = II - NX_T+1
               IF (GRID_IN(IK,JJ).GE.LAND_LIMIT) THEN
                  NSEA = NSEA+1
                  SEA = SEA + MIN(999.,GRID_IN(IK,JJ))  ! IN WAM 999M IS THE MAXIMUM DEPTH

                  IF (GRID_IN(IK,JJ).LT.500.) THEN      ! FIND SHALLOWER AREAS
                    NSEASH = NSEASH+1
                    SEASH = SEASH + GRID_IN(IK,JJ)
                  ENDIF

               ELSE
                  NLAND = NLAND+1
                  XLAND = XLAND+GRID_IN(IK,JJ)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

!          SEARCH FOR LAND AT THE CENTER OF THE GRID BOX

      NLANDCENTRE = 0
      DO JJ = J-NLANDCENTREPM, J+NLANDCENTREPM
         IF (JJ.GE.1 .AND. JJ.LE.NY_T) THEN
            DO II = I-NLANDCENTREPM, I+NLANDCENTREPM
               IK = II
               IF (II.LT.1)    IK = NX_T-1+II
               IF (II.GT.NX_T) IK = II-NX_T+1
               IF (GRID_IN(IK,JJ).LT.0.) NLANDCENTRE = NLANDCENTRE+1
            ENDDO
         ENDIF
      ENDDO

!   IF 40% OR MORE LAND, THEN AVERAGE OVER LAND POINTS OR THE CENTER OF THE GRID BOX IS LAND.
!          ELSE AVERAGE OVER SEA POINTS

      PERCENTLAND(IX,K) = FLOAT(NLAND)/FLOAT(NLAND+NSEA)
      IF (PERCENTLAND(IX,K).GT.0.60 .OR. NLANDCENTRE.GE.NLANDCENTREMAX ) THEN
         GRD(IX,K) = XLAND/NLAND
      ELSE

!   IF THERE IS A PERCENTAGE OF SHALLOWER POINTS THEN THE AVERAGE IS TAKEN OVER THOSE POINTS ALONE.

         PERCENTSHALLOW(IX,K) = FLOAT(NSEASH)/FLOAT(NSEA)
         IF (PERCENTSHALLOW(IX,K).GE.0.3) THEN
            GRD(IX,K) = SEASH/NSEASH
         ELSE 
            GRD(IX,K) = SEA/NSEA
         ENDIF
      ENDIF
   ENDDO
ENDDO

GRD = MAX(GRD,-999.)
GRD = MIN(GRD, 999.)

END SUBROUTINE MEAN_DEPTH

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_GRID

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PREPARE_GRID - ROUTINE TO ARRANGE WAMODEL GRID.                            !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO ARRANGE WAMODEL GRID FOR A GIVEN AREA AND COMPUTE VARIOUS           !
!       MODEL CONSTANTS.                                                       !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE MODEL GRID AREA IS EXTRACTED FROM THE INPUT TOPOGRAGHY.            !
!       THE NEAREST GRID POINT IS TAKEN FOR DEPTH INTERPOLATION.               !
!       LAND POINTS ARE REMOVED FROM THE GRID.                                 !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: I, K, KH, IP
REAL    :: XLAT, XH

character (len=len_coor) :: formtext

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLEAR MODEL GRID AND CHECK MODULE STATUS.                             !
!        -----------------------------------------                             !

NSEA = -1
IF (ALLOCATED (SINPH))  DEALLOCATE (SINPH)
IF (ALLOCATED (COSPH))  DEALLOCATE (COSPH)
IF (ALLOCATED (DEPTH_B))  DEALLOCATE (DEPTH_B)
IF (ALLOCATED (L_S_MASK))  DEALLOCATE (L_S_MASK)

IF (ALLOCATED (IXLG))    DEALLOCATE (IXLG)
IF (ALLOCATED (KXLT))    DEALLOCATE (KXLT)
IF (ALLOCATED (KLAT))    DEALLOCATE (KLAT)
IF (ALLOCATED (KLON))    DEALLOCATE (KLON)
IF (ALLOCATED (WLAT))    DEALLOCATE (WLAT)
IF (ALLOCATED (GRD))     DEALLOCATE (GRD)
IF (ALLOCATED (NLON_RG)) DEALLOCATE (NLON_RG)
IF (ALLOCATED (ZDELLO))  DEALLOCATE (ZDELLO)
IF (ALLOCATED (DELLAM))  DEALLOCATE (DELLAM)

IF (.NOT.SET_STATUS) THEN
   WRITE (IU06,*) ' ********************************************************'
   WRITE (IU06,*) ' *                                                      *'
   WRITE (IU06,*) ' *          FATAL  ERROR IN SUB. PREPARE_GRID           *'
   WRITE (IU06,*) ' *          =================================           *'
   IF (NX_T .EQ. -1 .OR. NY_T .EQ. -1) THEN
      WRITE (IU06,*) ' *                                                      *'
      WRITE (IU06,*) ' * TOPOGRAPHY INPUT DATA ARE NOT DEFINED IN GRID MODULE.*'
      WRITE (IU06,*) ' * DATA MUST BE DEFINED BY SUB. SET_TOPOGRAPHY.         *'
   END IF
   IF (NX .EQ. -1 .OR. NY .EQ. -1) THEN
      WRITE (IU06,*) ' *                                                      *'
      WRITE (IU06,*) ' * MODEL GRID DIMENISIONS ARE NOT DEFINED GRID MODULE.  *'
      WRITE (IU06,*) ' * DATA MUST BE DEFINED BY SUB. SET_GRID_DEF.           *'
   END IF
   WRITE (IU06,*) ' *                                                      *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS             *'
   WRITE (IU06,*) ' *                                                      *'
   WRITE (IU06,*) ' ********************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!    2. INITIALISE SIN AND COS OF LATITUDES.                                   !
!       ------------------------------------                                   !

ALLOCATE (SINPH(1:NY))
ALLOCATE (COSPH(1:NY))

DO K = 1,NY
   XLAT = REAL(AMOSOP + (K-1)*XDELLA)/REAL(M_DEGREE)*RAD
   SINPH(K)   = SIN(XLAT)
   COSPH(K)   = COS(XLAT)
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. INITIALISE GRID INCREMENTS.                                            !
!       ---------------------------                                            !

DELPHI = REAL(XDELLA)/REAL(M_DEGREE)*CIRC/360.  !! LATITUDE INCREMENT IN  METRES.

ALLOCATE (NLON_RG(NY))

IF (REDUCED_GRID) THEN
   NLON_RG(:) = NINT(NX*COSPH(:)/MAXVAL(COSPH(:)))
   WHERE (MOD(NLON_RG(:),2).EQ.1) NLON_RG(:) = NLON_RG(:)+1
   NLON_RG(:) = MIN(NLON_RG(:),NX)
ELSE
   NLON_RG(:) = NX
END IF

IF (REDUCED_GRID .AND. MINVAL(NLON_RG(:)).LE.0) THEN
   WRITE (IU06,*) '********************************************************'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*          FATAL ERROR IN SUB. PREPARE_GRID            *'
   WRITE (IU06,*) '*          ================================            *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*  A REDUCED GRID CANNOT BE SET UP.                    *'
   WRITE (IU06,*) '*  THE INCREASED LONGITUDE INCREMENT IS GREATER THAN   *'
   WRITE (IU06,*) '*  THE TOTAL LATITUDE LENGTH IN THE GRID AREA.         *'
   WRITE (IU06,*) '*  THE LATITUDE HAS NO GRID POINT (NLON < 2).          *'
   WRITE (IU06,*) '*                                                      *'

   WRITE (IU06,*) '  NO. |    LATITUDE   |   NLON |'
   WRITE (IU06,*) '------|---------------|--------|'
   DO I = NY, 1, -1
      formtext = write_coor_text (AMOSOP + (I-1)*XDELLA)
      WRITE (IU06,'(I6,'' | '',A,'' | '',I6,'' | '')')                         &
&             I, formtext, NLON_RG(I)  
   END DO

   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*   REDUCE THE GRID SOUTH-NORTH EXTENT OR              *'
   WRITE (IU06,*) '*   DECREASE THE BASIC LONGITUDE INCREMENT.            *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*              THE PROGRAM ABORTS                      *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '********************************************************'
   CALL ABORT1
ENDIF

ALLOCATE (ZDELLO(NY))
ALLOCATE (DELLAM(NY))

IF (IPER) THEN             !! LONG. INCREMENTS IN (M_SEC)
   ZDELLO(:) = NINT(REAL(AMOEAP-AMOWEP+XDELLO)/REAL(NLON_RG(:)))
ELSE
   ZDELLO(:) = NINT(REAL(AMOEAP-AMOWEP)/REAL(NLON_RG(:)-1))
END IF

DELLAM(:) = REAL(ZDELLO(:))/REAL(M_DEGREE)*CIRC/360.  !! LONG. INCREMENTS IN (M)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    4. ARRANGE THE TOPOGRAPHY ON REQUESTED MODEL AREA AND RESOLUTION.         !
!       -------------------------------------------------------------          !

IF (.NOT.ALLOCATED(GRD)) ALLOCATE(GRD(NX,NY))

GRD(1:NX,1:NY) = -999.

IF (EQUAL_GRID) THEN
   GRD(1:NX,1:NY) = GRID_IN(1:NX_T,1:NY_T)
   IF (L_OBSTRUCTION_T) THEN
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +       WARNING ERROR SUB. PREPARE_GRID.           +'
      WRITE (IU06,*) ' +       ================================           +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' + THE REQUESTED MODEL GRID AND                     +'
      WRITE (IU06,*) ' + THE INPUT BATHYMETRY DATA GRID ARE IDENTICAL.    +'
      WRITE (IU06,*) ' + NO OBSTRUCTIONS CAN BE COMPUTED.                 +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +      MODEL CONTINUES WITHOUT OBSTRUCTIONS        +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      L_OBSTRUCTION_T = .FALSE.
   END IF

   CALL ADJUST_DEPTH
   IF (ITEST.GT.1) WRITE (IU06,*) '   SUB. PREPARE_GRID: ADJUST_DEPTH DONE'

ELSE

   IF (L_INTERPOL_T) THEN
      CALL SUBGRID_TOPOGRAPHY
      IF (ITEST.GT.1) WRITE (IU06,*) '   SUB. PREPARE_GRID: SUBGRID_TOPOGRAPHY DONE '
   ELSE
      CALL MEAN_DEPTH
      IF (ITEST.GT.1) WRITE (IU06,*) '   SUB. PREPARE_GRID: MEAN_DEPTH DONE '
   END IF

   CALL ADJUST_DEPTH
   IF (ITEST.GT.1) WRITE (IU06,*) '   SUB. PREPARE_GRID: ADJUST_DEPTH DONE'
END IF

WHERE (GRD.LT.LAND_LIMIT) GRD = -999.  !! MARK LAND BY -999

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. GENERATE LAND-SEA MASK AND COUNT NUMBER OF SEA POINTS.                !
!        ------------------------------------------------------                !

ALLOCATE (L_S_MASK(1:NX,1:NY))
L_S_MASK = GRD.GE.LAND_LIMIT
NSEA = COUNT (L_S_MASK)

IF (NSEA.LE.0) THEN
   WRITE (IU06,*) '********************************************************'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*          FATAL ERROR IN SUB. PREPARE_GRID            *'
   WRITE (IU06,*) '*          ================================            *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*  THE MODEL GRID DOES NOT CONTAIN SEA POINTS.         *'
   WRITE (IU06,*) '*  NUMBER OF SEA POINTS IS NSEA = ', NSEA
   formtext = write_coor_text (amowep)
   WRITE (IU06,*) ' + LONGITUDE               WEST = ', formtext
   formtext = write_coor_text (amoeap)
   WRITE (IU06,*) ' + LONGITUDE               EAST = ', formtext
   formtext = write_coor_text (xdello)
   WRITE (IU06,*) ' + LONGITUDE INCREMENT    D_LON = ', formtext
   WRITE (IU06,*) ' + NO. OF LONGITUDES      N_LON = ', NX
   formtext = write_coor_text (amosop)
   WRITE (IU06,*) ' + LATITUDE               SOUTH = ', formtext
   formtext = write_coor_text (amonop)
   WRITE (IU06,*) ' + LATITUDE               NORTH = ', formtext
   formtext = write_coor_text (xdella)
   WRITE (IU06,*) ' + LATITUDE  INCREMENT    D_LAT = ', formtext
   WRITE (IU06,*) ' + NO. OF LATITUDE        N_LAT = ', NY
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*  DOES INPUT DEPTH DATA GRID INCLUDE MODEL GRID?      *'
   WRITE (IU06,*) '*              THE PROGRAM ABORTS                      *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '********************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. COMPUTE ARRAYS TO MAP POINT INDEX TO GRID POINT INDEX.                !
!        ------------------------------------------------------                !

ALLOCATE (IXLG (1:NSEA))
ALLOCATE (KXLT (1:NSEA))

IP = 0
DO K = 1,NY
   DO I = 1,NX
      IF (L_S_MASK(I,K)) THEN
         IP = IP+1
         IXLG(IP) = I
         KXLT(IP) = K
      END IF
   END DO
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. REMOVE LAND POINT FROM DEPTH ARRAY.                                   !
!        -----------------------------------                                   !

ALLOCATE (DEPTH_B(1:NSEA))
DEPTH_B = PACK (GRD, L_S_MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     8. COMPUTE INDICES OF NEIGHBOUR SEA POINTS AND WEIGHT IN ADVECTION SCHEME!
!        FOR CLOSEST GRIDPOINT IN NORTH-SOUTH DIRECTION.                       !
!        ----------------------------------------------------------------------!

ALLOCATE (KLAT(1:NSEA,1:2,1:2))       !! CLOSEST AND SECOND CLOSEST POINTS 
                                      !! ON SOUTH AND NORTH LATITUDE.
ALLOCATE (KLON(1:NSEA,1:2))           !! NEXT WEST AND EAST POINTS ON LATITUDE.
KLON = 0
ALLOCATE (WLAT(1:NSEA,1:2))           !! WEIGHT IN ADVECTION FOR N-S DIRECTION.

IP = 0
DO IP = 1, NSEA
   I = IXLG(IP)
   K = KXLT(IP)

!     8.1 WEST LONGITUDE NEIGHBOURS.                                           !

   IF (I.GT.1) THEN
      IF (L_S_MASK(I-1,K)) KLON(IP,1) = IP-1
   ELSE
      IF (IPER.AND.L_S_MASK(NLON_RG(K),K)) KLON(IP,1) = IP+COUNT(L_S_MASK(2:NX,K))
   END IF

!     8.2 EAST LONGITUDE NEIGHBOURS.                                           !

   IF (I.LT.NLON_RG(K)) THEN
      IF (L_S_MASK(I+1,K)) KLON(IP,2) = IP+1
   ELSE
      IF (IPER.AND.L_S_MASK(1,K)) KLON(IP,2) = IP+1-COUNT(L_S_MASK(1:NX,K))
   END IF

!     8.3 CLOSEST AND SECOND CLOSEST SOUTH LATITUDE NEIGHBOURS.                !

   IF (K.GT.1) THEN
      KH = K-1
      XH = REAL(I-1)*ZDELLO(K)/ZDELLO(KH)     ! LONGITUDE OF GRID POINT
      CALL CLOSEST_GP_LAT (I, K, IP, XH, KH, KLAT(IP,1,1), KLAT(IP,1,2), WLAT(IP,1))
   ELSE
      KLAT(IP,1,1) = 0
      KLAT(IP,1,2) = 0
      WLAT(IP,1)   = 1
   ENDIF

!     8.4 CLOSEST AND SECOND CLOSEST NORTH LATITUDE NEIGHBOURS.                !

   IF (K.LT.NY) THEN
      KH = K+1
      XH = REAL(I-1)*ZDELLO(K)/ZDELLO(KH)     ! LONGITUDE OF GRID POINT
      CALL CLOSEST_GP_LAT (I, K, IP, XH, KH, KLAT(IP,2,1), KLAT(IP,2,2), WLAT(IP,2))
   ELSE
      KLAT(IP,2,1) = 0
      KLAT(IP,2,2) = 0
      WLAT(IP,2)   = 1.
   ENDIF
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     11. CREATE OBSTRUCTIONS.                                                 !
!        ---------------------                                                 !

IF (L_OBSTRUCTION_T) THEN
   ALLOCATE (OBSLAT (NSEA,2,ML))
   ALLOCATE (OBSLON (NSEA,2,ML))
   CALL CREATE_OBSTRUCTIONS
   IF (ITEST.GT.1) WRITE(IU06,*) '   SUB. OBSTRUCTION DONE'
ENDIF
IF (.NOT.L_OBSTRUCTION_T) THEN
   IF (ALLOCATED(OBSLAT)) DEALLOCATE(OBSLAT)
   IF (ALLOCATED(OBSLON)) DEALLOCATE(OBSLON)
ENDIF

END SUBROUTINE PREPARE_GRID

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SUBGRID_TOPOGRAPHY

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SUBGRID_TOPOGRAPHY - ARRANGE SUBGRID TOPOGRAPHY.                           !
!                                                                              !
!     S. HASSELMANN     MPIFM           1/6/86.                                !
!                                                                              !
!     MODIFIED BY       H. GUNTHER      1/4/90  -  REARANGEMENT OF CODE.       !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO CONVERT THE INPUT GRID TO THE MODEL GRID.                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE TOPOGRAPHIC DATA ARE PUT ON THE REQUESTED SUBGRID LAT-LONG         !
!       RESOLUTION, ALWAYS USING THE NEAREST POINT.                            !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: I, K
INTEGER :: IH(NX), KH(NY)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. STORING TOPOGRAPHIC DATA AT LONGITUDES AND LATITUDES OF GRID AREA.    !
!        ------------------------------------------------------------------    !

DO K = 1,NY
   KH(K) = NINT( REAL(AMOSOP + (K-1)*XDELLA - SOUTH_T)/REAL(DY_T) ) + 1
END DO
IF (MINVAL(KH).LT.1 .OR. MAXVAL(KH).GT.NY_T) THEN
   WRITE (IU06,*) ' *****************************************************'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' *     FATAL  ERROR IN SUB. SUBGRID_TOPOGRAPHY       *'
   WRITE (IU06,*) ' *     =======================================       *'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' * INPUT TOPOGRAPHY DOES NOT FIT TO REQUESTED GRID.  *'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' *      PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                   *'
   WRITE (IU06,*) ' *****************************************************'
   CALL PRINT_PREPROC_STATUS
   CALL ABORT1
END IF

DO K=1,NY
   DO I=1,NLON_RG(K)
      IF (PER_T) THEN
         IH(I) = NINT( REAL(AMOWEP + (I-1)*ZDELLO(K)+M_S_PER - WEST_T)/REAL(DX_T) )
         IH(I) = MOD(IH(I)+NX_T,NX_T) + 1
      ELSE
         IH(I) = NINT( REAL(AMOWEP + (I-1)*ZDELLO(K) - WEST_T)/REAL(DX_T) ) + 1
      END IF
   END DO
   IF (MINVAL(IH(1:NLON_RG(K))).LT.1 .OR. MAXVAL(IH(1:NLON_RG(K))).GT.NX_T) THEN
      WRITE (IU06,*) ' *****************************************************'
      WRITE (IU06,*) ' *                                                   *'
      WRITE (IU06,*) ' *     FATAL  ERROR IN SUB. SUBGRID_TOPOGRAPHY       *'
      WRITE (IU06,*) ' *     =======================================       *'
      WRITE (IU06,*) ' *                                                   *'
      WRITE (IU06,*) ' * INPUT TOPOGRAPHY DOES NOT FIT TO REQUESTED GRID.  *'
      WRITE (IU06,*) ' *                                                   *'
      WRITE (IU06,*) ' *      PROGRAM ABORTS  PROGRAM ABORTS               *'
      WRITE (IU06,*) ' *                                                   *'
      WRITE (IU06,*) ' *****************************************************'
      CALL PRINT_PREPROC_STATUS
      CALL ABORT1
   END IF
   
   GRD(1:NLON_RG(K),K) = GRID_IN(IH(1:NLON_RG(K)),KH(K))
END DO

END SUBROUTINE SUBGRID_TOPOGRAPHY

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE PREPROC_MODULE

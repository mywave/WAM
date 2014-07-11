MODULE PREPROC_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL VARIABLES AND CONSTANTS NECESSARY FOR THE         !
!   PREPROC PROGRAM. ALL PROCEDURES ARE INCLUDED TO COMPUTE THE INFORMATION    !
!   STORED IN WAM_CONNST_MODULE, WAM_NEST_MODULE AND WAM_GRID_MODULE.          !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE          !! COORDINATE TYPE AND PROCEDURES

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1,                  & !! TERMINATES PROCESSING.
&       PRINT_ARRAY                !! PRINTS AN ARRAY.

USE WAM_GRID_MODULE,      ONLY:  &
&       EQUAL_TO_M_GRID,         & !! COMPARES TWO GRIDS.
&       PRINT_GRID_STATUS          !! PRINTS THE MODULE INFORMATION.

USE WAM_FRE_DIR_MODULE,   ONLY:  &
&       PRINT_FRE_DIR_STATUS       !! PRINTS THE MODULE INFORMATION.

USE WAM_NEST_MODULE,  ONLY:      &
&       PREPARE_BOUNDARY_nest,   & !! MAKES NEST BOUNDARIES.
&       PRINT_NEST_STATUS          !! PRINTS THE MODULE INFORMATION.

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
&                             FMIN, MO_TAIL, MM1_TAIL, MP1_TAIL, MP2_TAIL,     &
&                             NDEPTH, TCGOND, TSIHKD, TFAK, TFAC_ST

USE WAM_GRID_MODULE,    ONLY: HEADER, NX, NY, NSEA, NLON_RG, IPER,             &
&                             AMOWEP, AMOSOP, AMOEAP, AMONOP,                  &
&                             XDELLA, XDELLO, DELLAM, ZDELLO, DELPHI,          &
&                             SINPH, COSPH, DEPTH_B, KLAT, KLON, IXLG, KXLT,   &
&                             L_S_MASK, ONE_POINT, REDUCED_GRID

USE WAM_NEST_MODULE,    ONLY: N_NEST, MAX_NEST, N_NAME, n_code,                &
&                             NBOUNC, IJARC, BLATC, BLNGC, DLAMAC, DPHIAC,     &
&                             N_SOUTH, N_NORTH, N_EAST, N_WEST, N_ZDEL,        &
&                             NBINP, NBOUNF, C_NAME, BLNGF, BLATF,             &
&                             IJARF, IBFL, IBFR, BFW

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
REAL,    ALLOCATABLE :: GRID_IN(:,:)  !! WATER DEPTH [M].
REAL                 :: LAND_LIMIT    !! DEPTH <= LAND_LIMIT ARE NOT ACTIVE

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

REAL, ALLOCATABLE, DIMENSION(:,:) :: GRD   !! GRIDDED TOPOGRAPHY

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

INTERFACE PREPARE_GRID               !! ROUTINE TO ARRANGE WAMODEL GRID.
   MODULE PROCEDURE PREPARE_GRID
END INTERFACE
PRIVATE PREPARE_GRID

INTERFACE SUBGRID_TOPOGRAPHY         !! ARRANGE SUBGRID TOPOGRAPHY.
   MODULE PROCEDURE SUBGRID_TOPOGRAPHY
END INTERFACE
PRIVATE SUBGRID_TOPOGRAPHY

INTERFACE ADJUST_DEPTH               !! ADJUST DEPTH.
   MODULE PROCEDURE ADJUST_DEPTH
END INTERFACE
PRIVATE ADJUST_DEPTH

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
&                          LAND, R_GRID)

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
REAL,    OPTIONAL,        INTENT(IN) :: LAND    !! DEPTH >= LAND ARE SEAPOINTS [M].
LOGICAL, OPTIONAL,        INTENT(IN) :: R_GRID  !! REDUCED GRID OPTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CONVERT TO M_SEC AND CALL SET_GRID_CORRECTIONS_M.                     !
!        -------------------------------------------------                     !

IF (PRESENT(LAND) .AND. PRESENT(R_GRID)) THEN
   CALL SET_GRID_DEF_M (N_LON, N_LAT,                                          & 
&                       READ_COOR_TEXT(D_LON), READ_COOR_TEXT(D_LAT),          &
&                       READ_COOR_TEXT(SOUTH), READ_COOR_TEXT(NORTH),          &
&                       READ_COOR_TEXT(WEST),  READ_COOR_TEXT(EAST),           &
&                       LAND, R_GRID)
ELSE IF (PRESENT(LAND)) THEN
   CALL SET_GRID_DEF_M (N_LON, N_LAT,                                          & 
&                       READ_COOR_TEXT(D_LON), READ_COOR_TEXT(D_LAT),          &
&                       READ_COOR_TEXT(SOUTH), READ_COOR_TEXT(NORTH),          &
&                       READ_COOR_TEXT(WEST),  READ_COOR_TEXT(EAST),           &
&                       LAND=LAND)
ELSE IF (PRESENT(R_GRID)) THEN
   CALL SET_GRID_DEF_M (N_LON, N_LAT,                                          & 
&                       READ_COOR_TEXT(D_LON), READ_COOR_TEXT(D_LAT),          &
&                       READ_COOR_TEXT(SOUTH), READ_COOR_TEXT(NORTH),          &
&                       READ_COOR_TEXT(WEST),  READ_COOR_TEXT(EAST),           &
&                       R_GRID=R_GRID)
ELSE
   CALL SET_GRID_DEF_M (N_LON, N_LAT,                                          & 
&                       READ_COOR_TEXT(D_LON), READ_COOR_TEXT(D_LAT),          &
&                       READ_COOR_TEXT(SOUTH), READ_COOR_TEXT(NORTH),          &
&                       READ_COOR_TEXT(WEST),  READ_COOR_TEXT(EAST))
END IF

END SUBROUTINE SET_GRID_DEF_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_GRID_DEF_D (N_LON, N_LAT, D_LON, D_LAT,                         &
&                          SOUTH, NORTH, WEST, EAST,                           &
&                          LAND, R_GRID)

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
REAL,     OPTIONAL, INTENT(IN) :: LAND    !! DEPTH >= LAND ARE SEAPOINTS [M].
LOGICAL,  OPTIONAL, INTENT(IN) :: R_GRID  !! REDUCED GRID OPTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CONVERT TO M_SEC AND CALL SET_GRID_CORRECTIONS_M.                     !
!        -------------------------------------------------                     !

IF (PRESENT(LAND) .AND. PRESENT(R_GRID)) THEN
   CALL SET_GRID_DEF_M (N_LON, N_LAT,                                          & 
&                       DEG_TO_M_SEC(D_LON), DEG_TO_M_SEC(D_LAT),              &
&                       DEG_TO_M_SEC(SOUTH), DEG_TO_M_SEC(NORTH),              &
&                       DEG_TO_M_SEC(WEST),  DEG_TO_M_SEC(EAST),               &
&                       LAND, R_GRID)
ELSE IF (PRESENT(LAND)) THEN
   CALL SET_GRID_DEF_M (N_LON, N_LAT,                                          & 
&                       DEG_TO_M_SEC(D_LON), DEG_TO_M_SEC(D_LAT),              &
&                       DEG_TO_M_SEC(SOUTH), DEG_TO_M_SEC(NORTH),              &
&                       DEG_TO_M_SEC(WEST),  DEG_TO_M_SEC(EAST),               &
&                       LAND=LAND)
ELSE IF (PRESENT(R_GRID)) THEN
   CALL SET_GRID_DEF_M (N_LON, N_LAT,                                          & 
&                       DEG_TO_M_SEC(D_LON), DEG_TO_M_SEC(D_LAT),              &
&                       DEG_TO_M_SEC(SOUTH), DEG_TO_M_SEC(NORTH),              &
&                       DEG_TO_M_SEC(WEST),  DEG_TO_M_SEC(EAST),               &
&                       R_GRID=R_GRID)
ELSE
   CALL SET_GRID_DEF_M (N_LON, N_LAT,                                          & 
&                       DEG_TO_M_SEC(D_LON), DEG_TO_M_SEC(D_LAT),              &
&                       DEG_TO_M_SEC(SOUTH), DEG_TO_M_SEC(NORTH),              &
&                       DEG_TO_M_SEC(WEST),  DEG_TO_M_SEC(EAST))
END IF

END SUBROUTINE SET_GRID_DEF_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_GRID_DEF_M (N_LON, N_LAT, D_LON, D_LAT,                         &
&                          SOUTH, NORTH, WEST, EAST,                           &
&                          LAND, R_GRID)

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
REAL,    OPTIONAL, INTENT(IN) :: LAND    !! DEPTH >= LAND ARE SEAPOINTS [M].
LOGICAL, OPTIONAL, INTENT(IN) :: R_GRID  !! REDUCED GRID OPTION.


! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

LOGICAL  :: ERROR = .FALSE.              !! ERROR FLAG
character (len=len_coor) :: formtext

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLEAR GRID DEFINITIONS.                                               !
!        ----------------------_                                               !

NX     = -1     !! NUMBER OF COLUMNES IN TOPO INPUT GRID.
NY     = -1     !! NUMBER OF ROWS      IN TOPO INPUT GRID.
XDELLO = COOR_UNDEF     !! STEPSIZE BETWEEN LONGITUDES.
XDELLA = COOR_UNDEF     !! STEPSIZE BETWEEN LATITUDES.
AMOSOP = COOR_UNDEF     !! MOST SOUTHERN LATITUDE.
AMONOP = COOR_UNDEF     !! MOST NORTHERN LATITUDE.
AMOWEP = COOR_UNDEF     !! LEFT MOST LONGITUDE.
AMOEAP = COOR_UNDEF     !! RIGHT MOST LONGITUDE.
LAND_LIMIT   = 0.
REDUCED_GRID = .FALSE.
ONE_POINT    = .FALSE.
IPER         = .FALSE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COPY GRID DEFINITIONS AND CONVERT TO M_SEC.                           !
!        -------------------------------------------                           !

AMOWEP = WEST
AMOSOP = SOUTH
AMOEAP = EAST
AMONOP = NORTH
XDELLO = D_LON
XDELLA = D_LAT

IF (PRESENT(R_GRID)) REDUCED_GRID = R_GRID

IF (PRESENT(LAND)) LAND_LIMIT = LAND

NX = N_LON
NY = N_LAT

IF (N_LON .EQ.-1 .AND. N_LAT .EQ.-1 .AND.                                      &
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
      XDELLA = M_S_PER
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
!        -------------------------------------------                            !

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
   CALL SET_GRID_DEF (N_LON, N_LAT, D_LON, D_LAT, SOUTH, NORTH, WEST, EAST,    &
&                     LAND=LAND, R_GRID=R_GRID)
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
      DO I=1,NY
          C_GRID(NLON_RG(I)+1:NX,I) = ' '
      END DO
      
      CALL PRINT_ARRAY (IU06, ZERO, TITL, C_GRID,                             &
&                       AMOWEP, AMOSOP, AMOEAP, AMONOP)
   END IF
   IF (ITEST.GT.9) THEN
      WRITE (IU06,*) ' '
      IGRID=UNPACK((/(I,I=1,NSEA)/), L_S_MASK,-99)
      IGRID=MOD(IGRID,1000)
      TITL = 'SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP, AMONOP)
      IGRID=UNPACK(KLON(:,1), L_S_MASK,-99)
      IGRID=MOD(IGRID,1000)
      TITL = 'WEST NEIGHBOUR SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP, AMONOP)
      WRITE (IU06,*) ' '
      IGRID=UNPACK(KLON(:,2), L_S_MASK,-99)
      IGRID=MOD(IGRID,1000)
      TITL = 'EAST NEIGHBOUR SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP, AMONOP)
      WRITE (IU06,*) ' '
      IGRID=UNPACK(KLAT(:,1), L_S_MASK,-99)
      IGRID=MOD(IGRID,1000)
      TITL = 'SOUTH NEIGHBOUR SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP, AMONOP)
      WRITE (IU06,*) ' '
      IGRID=UNPACK(KLAT(:,2), L_S_MASK,-99)
      IGRID=MOD(IGRID,1000)
      TITL = 'NORTH NEIGHBOUR SEA POINT NUMBER (MODULO 1000)'
      CALL PRINT_ARRAY (IU06, ZERO, TITL,IGRID, AMOWEP, AMOSOP, AMOEAP, AMONOP)
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. NEST INFORMATION.                                                     !
!        -----------------                                                     !

CALL PRINT_NEST_STATUS

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
&            FMIN, MO_TAIL, MM1_TAIL, MP1_TAIL, MP2_TAIL

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. WRITE GRID INFORMATION.                                               !
!        -----------------------                                               !

WRITE (IU07) NX, NY, NSEA, IPER, ONE_POINT, REDUCED_GRID
WRITE (IU07) NLON_RG
WRITE (IU07) DELPHI, DELLAM, SINPH, COSPH, AMOWEP, AMOSOP, AMOEAP, AMONOP,     &
&            XDELLA, XDELLO, ZDELLO
WRITE (IU07) IXLG, KXLT, L_S_MASK
WRITE (IU07) KLAT, KLON, DEPTH_B

! ---------------------------------------------------------------------------- !
!                                                                              !
!     8. WRITE SHALLOW WATER TABLES.                                           !
!        ---------------------------                                           !

WRITE (IU07) TCGOND, TFAK, TSIHKD, TFAC_ST

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

INTEGER :: I, K, IP, IH
REAL    :: XLAT

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
ALLOCATE (ZDELLO(NY))
ALLOCATE (DELLAM(NY))

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

IF (IPER) THEN             !! LONG. INCREMENTS IN (M_SEC)
   ZDELLO(:) = NINT(REAL(AMOEAP-AMOWEP+XDELLO)/REAL(NLON_RG(:)))
ELSE
   ZDELLO(:) = NINT(REAL(AMOEAP-AMOWEP)/REAL(NLON_RG(:)-1))
END IF

DELLAM(:) = REAL(ZDELLO(:))/REAL(M_DEGREE)*CIRC/360.  !! LONG. INCREMENTS IN (M)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    4. ARRANGE THE TOPOGRAPHYON REQUESTED MODEL AREA AND RESOLUTION.          !
!       -------------------------------------------------------------          !

WHERE (GRID_IN.LT.LAND_LIMIT) GRID_IN = -999.  !! MARK LAND BY -999
CALL SUBGRID_TOPOGRAPHY
CALL ADJUST_DEPTH

IF (ITEST.GT.1) WRITE (IU06,*) ' SUB SUBGRID_TOPOGRAPHY DONE'

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

IXLG = 0
KXLT = 0

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
!     8. COMPUTE INDICES OF NEIGHBOUR SEA POINTS.                              !
!        ----------------------------------------                              !

ALLOCATE (KLAT(1:NSEA,1:2))   !! NEXT POINTS ON LATITUDE.
ALLOCATE (KLON(1:NSEA,1:2))   !! NEXT POINTS ON LONGITUDE.
KLAT = 0
KLON = 0

IP = 0
DO K = 1,NY
   DO I = 1,NX
      IF ( .NOT.L_S_MASK(I,K) ) CYCLE
      IP = IP+1

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

!     8.3 SOUTH LATITUDE NEIGHBOURS.                                           !

      IF (K.GT.1) THEN
         IH = NINT (REAL(I-1)*ZDELLO(K)/ZDELLO(K-1)) + 1
         IF (L_S_MASK(IH,K-1)) THEN
            KLAT(IP,1) = IP+1-COUNT(L_S_MASK(1:I,K))-COUNT(L_S_MASK(IH:NX,K-1))
         END IF
      END IF

!     8.4 NORTH LATITUDE NEIGHBOURS.                                           !

      IF (K.LT.NY) THEN
         IH = NINT(REAL(I-1)*ZDELLO(K)/ZDELLO(K+1)) + 1
         IF (L_S_MASK(IH,K+1)) THEN
            KLAT(IP,2) = IP-1+COUNT(L_S_MASK(I:NX,K))+COUNT(L_S_MASK(1:IH,K+1))
         END IF
      END IF
   END DO
END DO

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

IF (.NOT.ALLOCATED(GRD)) ALLOCATE(GRD(NX,NY))

GRD(1:NX,1:NY) = -999.

IF (EQUAL_GRID) THEN
   GRD(1:NX,1:NY) = GRID_IN(1:NX_T,1:NY_T)
   RETURN
END IF


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

END MODULE PREPROC_MODULE

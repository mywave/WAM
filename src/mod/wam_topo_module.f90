MODULE WAM_TOPO_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE STORES THE DEPTH INPUT GRID SPECFICATIONS AND THE MODEL DEPTH  !
!   FIELDS WHICH ARE PASSED WAM-MODELL.                                        !
!   IT CONTAINS ALL PROCEDURES NESSECARY FOR THE PROCESSING.                   !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE          !! COORDINATE PROCEDURES

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1,                  & !! TERMINATES PROCESSING.
&       INCDATE                    !! INCREMENTS DATE TIME GROUP.

USE WAM_GRID_MODULE,      ONLY:  &
&       INTERPOLATION_TO_GRID,   & !! INTERPOLATE TO WAM POINTS.
&       EQUAL_TO_M_GRID            !! COMPARES WIND GRID TO MODEL GRID.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST
USE WAM_FRE_DIR_MODULE, ONLY: KL, ML, NDEPTH, DEPTHA, DEPTHD
USE WAM_GRID_MODULE,    ONLY: NSEA, DEPTH_B, L_S_MASK
USE WAM_MODEL_MODULE,   ONLY: DEPTH, INDEP
USE WAM_TIMOPT_MODULE,  ONLY: TOPO_RUN, CDTA, CDATEE, IDEL_WAM

use wam_mpi_module,     only: nijs, nijl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INPUT TOPO GRID SPECFICATIONS, DATE AND TOPO FIELD.                   !
!        ---------------------------------------------------                   !

INTEGER :: NX_IN = -1    !! NUMBER OF LONGITUDES.
INTEGER :: NY_IN = -1    !! NUMBER OF LATITUDES.
LOGICAL :: PER = .FALSE. !! .TRUE. IF GRID IS PERIODICAL.
INTEGER :: CODE_IN = 1   !! TOPO CODE: 1 = TOTAL WATER DEPTH
                         !! OTHERWISE  ELEVATION OVER NN
INTEGER :: DX_IN         !! STEPSIZE BETWEEN LONGITUDES [M_SEC].
INTEGER :: DY_IN         !! STEPSIZE BETWEEN LATITUDES [M_SEC].
INTEGER :: SOUTH_IN      !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER :: NORTH_IN      !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER :: WEST_IN       !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER :: EAST_IN       !! EAST LONGITUDE OF GRID [M_SEC].
LOGICAL :: EQUAL_GRID =.FALSE. !! .TRUE. IF TOPO GRID IS EQUAL TO MODEL GRID.

CHARACTER (LEN=14) :: CD_READ = ' '  !! DATE OF LAST DATA READ FROM INPUT.
REAL, ALLOCATABLE, DIMENSION(:,:) :: TOPO_IN  !! INPUT TOPO DATA [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TOPO TIMESTEPS AND DRY POINTS FLAG.                                   !
!        -----------------------------------                                    !

INTEGER, PUBLIC      :: N_DRY   = 0       !! NUMBER OF DRY POINTS.
INTEGER, ALLOCATABLE :: IJ_DRY(:)         !! INDEX OF DRY POINTS.

INTEGER, PUBLIC   :: IDELTI = -1  !! INPUT TOPO TIMESTEP IN SECONDS.
INTEGER, PUBLIC   :: IDELTO = -1  !! OUTPUT TOPO TIMESTEP IN SECONDS
                                  !! EQUAL TO INPUT TIMESTEP INTO WAMODEL.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. TOPO FIELDS PREPARED BY PREPARE_DEPTH FOR WAM-MODEL.                  !
!        (FIRST INDEX IS POINTS, SECOND IS TIME)                               !
!        ----------------------------------------------                        !

INTEGER  :: M_STORE = 0    !! NUMBER OF DEPTH FIELDS STORED IN TRANSFER ARRAYS.

REAL,               ALLOCATABLE, DIMENSION(:,:) :: D_STORE   !! DEPTH DATA [M].
CHARACTER (LEN=14), ALLOCATABLE, DIMENSION(:)   :: CD_STORE  !! DATE/TIME.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE GET_TOPO                       !! GETS DEPTH FORM THIS MODULE.
   MODULE PROCEDURE GET_TOPO
END INTERFACE
PUBLIC  :: GET_TOPO

INTERFACE PREPARE_TOPO                   !! PREPARES DEPTH DATA FOR WAVE MODEL.
   MODULE PROCEDURE PREPARE_TOPO
END INTERFACE
PUBLIC :: PREPARE_TOPO

INTERFACE PRINT_TOPO_STATUS              !! PRINTS MODULE STATUS.
   MODULE PROCEDURE PRINT_TOPO_STATUS
END INTERFACE
PUBLIC PRINT_TOPO_STATUS

INTERFACE PUT_DRY                        !! MARKS DRY POINTS
   MODULE PROCEDURE PUT_DRY_PAR          !! IN PARAMETER ARRAY.
   MODULE PROCEDURE PUT_DRY_PAR_B        !! IN PARAMETER BLOCK.
   MODULE PROCEDURE PUT_DRY_SPEC         !! IN A SPECTRA BLOCK.
END INTERFACE
PUBLIC PUT_DRY

INTERFACE SET_TOPO_FIELD                 !! SETS TOPO FIELD.
   MODULE PROCEDURE SET_TOPO_FIELD
END INTERFACE
PUBLIC SET_TOPO_FIELD

INTERFACE SET_TOPO_HEADER                !! SETS TOPO HEADER.
   MODULE PROCEDURE SET_TOPO_HEADER_C    !! CHARACTER VERSION
   MODULE PROCEDURE SET_TOPO_HEADER_D    !! DEGREE VERSION
   MODULE PROCEDURE SET_TOPO_HEADER_M    !! M_SEC VERSION
END INTERFACE
PUBLIC SET_TOPO_HEADER

INTERFACE SET_TOPO_TIMESTEPS             !! SETS TOPO TIMESTEPS.
   MODULE PROCEDURE SET_TOPO_TIMESTEPS
END INTERFACE
PUBLIC SET_TOPO_TIMESTEPS

INTERFACE FIND_DRY_POINTS                !! FINDS DRY POINTS.
   MODULE PROCEDURE FIND_DRY_POINTS
END INTERFACE
PUBLIC FIND_DRY_POINTS

INTERFACE WAM_TOPO        !! READS AND TRANSFORMS INPUT TOPOS TO WAM POINTS.
   MODULE PROCEDURE WAM_TOPO
END INTERFACE
PUBLIC WAM_TOPO

INTERFACE
   SUBROUTINE READ_TOPO_INPUT            !! READ ONE TOPO FIELD.
   END SUBROUTINE READ_TOPO_INPUT         
END INTERFACE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE NOTIM           !! STEERING SUB IF TIME INTERPOLATION IS NOT WANTED.
   MODULE PROCEDURE NOTIM
END INTERFACE
PRIVATE NOTIM

INTERFACE TIMIN           !! STEERING SUB IF TIME INTERPOLATION IS WANTED.
   MODULE PROCEDURE TIMIN
END INTERFACE
PRIVATE TIMIN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE GET_TOPO (CD_NEW)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   GET_TOPO - TRANSFERS NEW DEPTH DATA TO WAM MODEL.                          !
!                                                                              !
!    H. GUENTHER  GKSS  DECEMBER 2009                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER DEPTH FROM WAM_TOPO_MODULE TO THE WAM_MODEL_MODULE.           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       COPY NEW DEPTH.                                                        !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=14),INTENT(INOUT) :: CD_NEW !! DATE TO USE A NEW TOPO FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER   :: IT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NEW TOPO FROM TRANSFER ARRAYS.                                        !
!        ------------------------------                                        !

CDTA = ' '
DO IT = 1, M_STORE
   IF (CD_STORE(IT).GT.CD_NEW) THEN
      DEPTH = D_STORE(:,IT)
      CDTA  = CD_STORE(IT)
      EXIT
   END IF
END DO
IF (CDTA .EQ. ' ') THEN
   WRITE(IU06,*) ' *******************************************'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *   FATAL ERROR IN SUB. GET_TOPO          *'
   WRITE(IU06,*) ' *   ============================          *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' * DEPTH FIELD WAS NOT FOUND               *'
   WRITE(IU06,*) ' *  CHANGE DATE IS ............ CD_NEW: ', CD_NEW
   WRITE(IU06,*) ' *  LAST DATE IN MODULE IS ...........: ', CD_STORE(M_STORE)
   WRITE(IU06,*) ' *  NO. OF DATES IN MODULE IS M_STORE = ', M_STORE
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *******************************************'
   CALL ABORT1
END IF

CALL FIND_DRY_POINTS

CALL INCDATE(CD_NEW,IDELTO)

IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '  '
   WRITE (IU06,*) '   SUB. GET_TOPO: NEW DEPTH FIELD PASSED TO WAM MODEL'
   WRITE (IU06,*) '     DATE OF DEPTH FIELD IS ......... CDTA = ', CDTA
   WRITE (IU06,*) '     NUMBER OF DRY POINTS IS ....... N_DRY = ', N_DRY
   WRITE (IU06,*) '     DATE FOR NEXT DEPTH CHANGE IS  CD_NEW = ', CD_NEW
   WRITE (IU06,*) '  '
END IF

END SUBROUTINE GET_TOPO

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_TOPO

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PREPARE_TOPO - PREPARES DEPTH DATA FOR WAVE MODEL.                         !
!                                                                              !
!    H. GUENTHER  GKSS  DECEMBER 2009                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       EVALUATE DEPTH DATA AT WAVE MODEL GRID POINTS.                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INPUT DEPTH FIELDS, WHICH CAN BE TOTAL WATER DEPTH OR SURFACE          !
!       ELEVATION OVER NN ARE TRANSFORMED INTO WATER DEPTH FOR THE WAM MODEL.  !
!       THE INPUT FIELDS HAVE TO BE ON A LAT /LONG GRID.                       !
!       SEE SUB READ_TOPO_INPUT FOR FORMATS AND HEADER INFORMATION,            !
!       WHICH HAVE TO BE GIVEN TO THE PROGRAM.                                 !
!                                                                              !
!       A DOUBLE LINEAR INTERPOLATION IN SPACE IS PERFORMED ONTO THE MODEL     !
!       GRID.                                                                  !
!       IF THE TOPO OUTPUT TIMSTEP IS LESS THAN THE INPUT TIMESTEP             !
!       A LINEAR INTERPOLATION IN TIME IS PERFORMED.                           !
!       ALL DEPTH FIELDS ARE STORED IN WAM_TOPO_MODULE.                        !
!                                                                              !
!     REFERENCE.                                                               !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

CHARACTER (LEN=14)            :: CD_START, CD_END 

! ---------------------------------------------------------------------------- !
!                                                                              !
!                                                                              !
!     1. BEGIN AND END DATES OF DEPTH FIELDS TO BE PROCESSED.                  !
!        ----------------------------------------------------                  !

CD_START = CDTA
CALL INCDATE (CD_START,IDELTO)
CD_END = CDTA
CALL INCDATE (CD_END,IDEL_WAM)
IF (CD_END.GE.CDATEE) CD_END = CDATEE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. NUMBER OF DEPTH FIELDS TO BE GENERATED.                               !
!        ----------------------------------------                              !

M_STORE = IDEL_WAM/IDELTO

IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '  '
   WRITE (IU06,*) '   SUB. PREPARE_TOPO: DEPTH REQUEST'
   WRITE (IU06,*) '     START OF PERIOD IS ............ CD_START = ', CD_START
   WRITE (IU06,*) '     END   OF PERIOD IS .............. CD_END = ', CD_END
   WRITE (IU06,*) '     TOPO INPUT TIME STEP ............ IDELTI = ', IDELTI
   WRITE (IU06,*) '     TOPO OUTPUT TIME STEP ........... IDELTO = ', IDELTO
   WRITE (IU06,*) '     NUMBER FIELDS TO BE GENERATED IS M_STORE = ', M_STORE
   WRITE (IU06,*) '     FIELDS ARE SAVED IN WAM_TOPO_MODULE'
END IF

IF ( ALLOCATED(CD_STORE))  then
   if (M_STORE.NE.SIZE(CD_STORE)) DEALLOCATE (CD_STORE, D_STORE)
endif
IF (.NOT. ALLOCATED(D_STORE) ) ALLOCATE (D_STORE(1:NSEA,M_STORE))
IF (.NOT. ALLOCATED(CD_STORE)) ALLOCATE (CD_STORE(M_STORE))
CD_STORE  = ' '

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. PROCESS TOPO FIELDS.                                                  !
!        --------------------                                                  !

IF (IDELTO.GE.IDELTI) THEN

   IF (ITEST.GE.2) WRITE (IU06,*) '     NO TIME INTERPOLATION'
   CALL NOTIM (CD_START, CD_END)

ELSE

   IF (ITEST.GE.2) WRITE (IU06,*) '     TIME INTERPOLATION'
   CALL TIMIN (CD_START, CD_END)

END IF

END SUBROUTINE PREPARE_TOPO

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FIND_DRY_POINTS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   FIND_DRY_POINTS - PREPARES TABLES FOR DRY POINTS.                          !
!                                                                              !
!    H. GUENTHER  GKSS  DECEMBER 2009                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       EVALUATE DEPTH DATA AT WAVE MODEL GRID POINTS.                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       GRID POINTS WITH WATER DEPTH LESS THAN MIN_DEPTH ARE DEFIEND AS DRY.   !
!       SEA POINT NUMBERS OF DRY POINTS ARE STORED IN AN INDEX ARRAY.          !
!       THE DEPTH IS CHANGED TO THE MIN_DEPTH.                                 !
!                                                                              !
!     REFERENCE.                                                               !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, N
REAL :: MIN_DEPTH

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DRY SEA POINTS.                                                       !
!        ---------------                                                       !

MIN_DEPTH = DEPTHA

N_DRY = COUNT(DEPTH.LT.MIN_DEPTH)
IF (ALLOCATED(IJ_DRY)) DEALLOCATE(IJ_DRY)
IF (N_DRY.GT.0) THEN
   ALLOCATE (IJ_DRY(1:N_DRY))
   N = 0
   DO IJ = 1,NSEA
      IF (DEPTH(IJ).LT.MIN_DEPTH) THEN
         N = N+1
         IJ_DRY(N) = IJ
         DEPTH(IJ) = 0.9*MIN_DEPTH
      END IF
   END DO
END IF
IF (ITEST.GE.3) THEN
   WRITE(IU06,*) '         SUB. FIND_DRY_POINTS: NUMBER OF DRY POINTS IS',    &
&                                                       ' N_DRY = ', N_DRY
END IF

END SUBROUTINE FIND_DRY_POINTS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_TOPO_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PRINT_TOPO_STATUS - PRINTS STATUS OF THE TOPO MODULE.                      !
!                                                                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2009                           !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       MAKE A PRINTER OUTPUT OF THE STATUS OF THE TOPO MODULE.                !
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

character (len=len_coor) :: formtext
   
WRITE (IU06,*) '  '
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '              TOPO MODULE STATUS:'
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '  '

IF (.NOT. TOPO_RUN) THEN
   WRITE (IU06,*) ' WATER DEPTH IS STATIONARY. '
   WRITE (IU06,*) ' NO. OF DRY POINTS IS N_DRY = ', N_DRY
   RETURN
END IF

WRITE (IU06,*) ' MODEL TOPO INPUT TIME STEP .................: ', IDELTI,' SECS'
WRITE (IU06,*) ' MODEL TOPO OUTPUT TIME STEP.................: ', IDELTO,' SECS'
IF (IDELTO.GT.0 .AND. IDELTI.GT.0) THEN
   IF (IDELTO.GE.IDELTI) THEN
      WRITE(IU06,*) ' TOPO FIELDS ARE NOT INTERPOLATED IN TIME'
   ELSE
      WRITE(IU06,*) ' TOPO FIELDS ARE INTERPOLATED IN TIME'
   END IF
END IF
WRITE (IU06,*) '  '
WRITE (IU06,*) ' DATE OF LAST TOPO FIELD GIVEN TO MODULE....: ', CD_READ
WRITE (IU06,*) ' DATE OF LAST TOPO FIELD GIVEN TO WAM MODEL.: ', CDTA
WRITE (IU06,*) '  '

IF (NX_IN.GT.0 .AND. NY_IN.GT.0) THEN
   WRITE (IU06,*) ' TOPO INPUT GRID SPECIFICATION ARE:'
   WRITE (IU06,*) ' NUMBER OF LONGITUDE ..........: ', NX_IN
   formtext = write_coor_text (west_in)
   WRITE (IU06,*) ' WESTERN MOST LONGITUDE........: ', formtext
   formtext = write_coor_text (east_in)
   WRITE (IU06,*) ' EASTERN MOST LONGITUDE........: ', formtext
   formtext = write_coor_text (dx_in)
   WRITE (IU06,*) ' LONGITUDE INCREMENT ..........: ', formtext
   WRITE (IU06,*) ' NUMBER OF LATITUDE ...........: ', NY_IN
   formtext = write_coor_text (south_in)
   WRITE (IU06,*) ' SOUTHERN MOST LATITUDE .......: ', formtext
   formtext = write_coor_text (north_in)
   WRITE (IU06,*) ' NORTHERN MOST LATITUDE .......: ', formtext
   formtext = write_coor_text (dy_in)
   WRITE (IU06,*) ' LATITUDE INCREMENT ...........: ', formtext
   IF (PER) THEN
      WRITE (IU06,*) ' THE GRID IS PERIODIC IN EAST-WEST DIRECTION'
   ELSE
      WRITE (IU06,*) ' THE GRID IS NOT PERIODIC IN EAST-WEST DIRECTION'
   END IF
   IF (CODE_IN.EQ.1) THEN
      WRITE (IU06,*) ' TOPO INPUTS ARE STOTAL WATER DEPTH'
   ELSE
      WRITE (IU06,*) ' TOPO INPUTS ARE ELEVATION OVER BASIC DEPTH'
   END IF
   WRITE (IU06,*) '  '
ELSE
   WRITE (IU06,*) ' TOPO INPUT GRID IS NOT DEFINED'
END IF

WRITE (IU06,*) ' NO. OF DRY POINTS IS........................: ', N_DRY

IF ( M_STORE.GT.0 .AND. ANY(CD_STORE.NE. ' ')) THEN
   WRITE (IU06,*) ' NUMBER OF DEPTH FIELDS STORED IN MODULE ....: ', M_STORE
   WRITE (IU06,*) ' DATES OF DEPTH FIELDS ARE:'
   WRITE (IU06,'(5(3X,A14,2X))')  CD_STORE
ELSE
   WRITE (IU06,*) ' DEPTH FIELDS ARE NOT STORED IN MODULE'
END IF
WRITE (IU06,*) '  '

END SUBROUTINE PRINT_TOPO_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PUT_DRY_PAR (PAR, INDICATOR)

REAL, INTENT(INOUT) :: PAR(:)      !! BLOCK OF PARAMETER.
REAL, INTENT(IN)    :: INDICATOR   !! VALUE TO BE INSERTED AT DRY POINTS.

IF (N_DRY.GT.0) PAR(IJ_DRY) = INDICATOR

END SUBROUTINE PUT_DRY_PAR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PUT_DRY_PAR_B (PAR, na, nb, INDICATOR)

integer, intent(in) :: na          !! first index
integer, intent(in) :: nb          !! last index
REAL, INTENT(INOUT) :: PAR(na:nb)  !! BLOCK OF PARAMETER.
REAL, INTENT(IN)    :: INDICATOR   !! VALUE TO BE INSERTED AT DRY POINTS.

INTEGER :: IJ

DO IJ = 1,N_DRY
   IF (IJ_DRY(IJ).GE.na .AND. IJ_DRY(IJ).LE.nb) PAR(IJ_DRY(IJ)) = INDICATOR
END DO

END SUBROUTINE PUT_DRY_PAR_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PUT_DRY_SPEC (FL3, INDICATOR)

REAL, INTENT(INOUT) :: FL3(nijs:,:,:)   !! BLOCK OF SPECTRA.
REAL, INTENT(IN)    :: INDICATOR        !! VALUE TO BE INSERTED AT DRY POINTS.

INTEGER :: IJ

DO IJ = 1, N_DRY
   IF (IJ_DRY(IJ).GE.nijs .AND. IJ_DRY(IJ).LE.nijl)                            &
&        FL3(IJ_DRY(IJ),:,:) = INDICATOR
END DO

END SUBROUTINE PUT_DRY_SPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_TOPO_FIELD (CDT, D_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_TOPO_FIELD - SET INPUT TOPO FIELD IN WAM_TOPO MODULE.                  !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER DEPTH TO WAM_TOPO_MODULE.                                     !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY NEW DEPTH. IF DEPTH IS GIVEN AS SURFACE     !
!       ELEVATION OVER NN THE BASIC DEPTH IS USED TO CALCULATE THE ACTUALL     !
!       WATER DEPTH TO BE USED BY THE MODEL.                                   !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=14),INTENT(IN) :: CDT        !! DATE/TIME OF TOPO FIELD.
REAL,              INTENT(IN) :: D_MAP(:,:) !! TOPO FIELD  [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COPY DATE.                                                            !
!        ----------                                                            !

CD_READ = CDT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK CONSISTENCY.                                                    !
!        ------------------                                                    !

IF (NX_IN.NE.SIZE(D_MAP,1) .OR. NY_IN.NE.SIZE(D_MAP,2) ) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *          FATAL  ERROR IN SUB. SET_TOPO_FIELD           *'
   WRITE (IU06,*) ' *          ===================================           *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' * GRID SPECIFICATIONS ARE NOT CONSISTENT                 *'
   WRITE (IU06,*) ' * OR TOPO HEADER IS NOT DEFINED:                         *'
   WRITE (IU06,*) ' * GRID SIZES AS DEFINED IN MODULE ARE:                   *'
   WRITE (IU06,*) ' * NO. OF LONGITUDES    NX_IN = ', NX_IN
   WRITE (IU06,*) ' * NO. OF LATITUDE      NY_IN = ', NY_IN
   WRITE (IU06,*) ' * DIMENSIONS OF TOPO ARRAY ARE :                         *'
   WRITE (IU06,*) ' * NO. OF LONGITUDES 1. DIMENSION = ', SIZE(D_MAP,1)
   WRITE (IU06,*) ' * NO. OF LATITUDE   2. DIMENSION = ', SIZE(D_MAP,2)
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COPY TOPO FIELDS.                                                     !
!        -----------------                                                     !

IF (ALLOCATED(TOPO_IN)) DEALLOCATE(TOPO_IN)
ALLOCATE (TOPO_IN(NX_IN,NY_IN))

TOPO_IN = D_MAP

END SUBROUTINE SET_TOPO_FIELD

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_TOPO_HEADER_C (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,          &
&                             N_LON, N_LAT, CODE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_TOPO_HEADER - SET INPUT GRID FOR TOPO FIELDS IN WAM_TOPO MODULE.       !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 20011                                         !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS FOR TOPO FIELDS TO WAM_TOPO_MODULE.!
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY DEFINITIONS.                                !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=LEN_COOR), INTENT(IN)           :: WEST  !! WEST LONGITUDE.
CHARACTER (LEN=LEN_COOR), INTENT(IN)           :: SOUTH !! SOUTH LATITUDE.
CHARACTER (LEN=LEN_COOR), INTENT(IN), OPTIONAL :: EAST  !! EAST LONGITUDE.
CHARACTER (LEN=LEN_COOR), INTENT(IN), OPTIONAL :: NORTH !! NORTH LATITUDE.
CHARACTER (LEN=LEN_COOR), INTENT(IN), OPTIONAL :: D_LON !! LONGITUDE INCREMENT.
CHARACTER (LEN=LEN_COOR), INTENT(IN), OPTIONAL :: D_LAT !! LATITUDE INCREMENT.
INTEGER,                  INTENT(IN), OPTIONAL :: N_LON !! NUMBER OF LONGITUDES.
INTEGER,                  INTENT(IN), OPTIONAL :: N_LAT !! NUMBER OF LATITUDES.
INTEGER, INTENT(IN), OPTIONAL :: CODE  !! 1 = TOTAL WATER DEPTH (DEFAULT)
                                       !! OTHERWISE ELEVATION OVER NN.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: WEST_CO    !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER :: SOUTH_CO   !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER :: EAST_CO    !! EAST LONGITUDE OF GRID [M_SEC].
INTEGER :: NORTH_CO   !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER :: D_LON_CO   !! LONGITUDE INCREMENT [M_SEC].
INTEGER :: D_LAT_CO   !! LATITUDE INCREMENT [M_SEC].
INTEGER :: N_LON_CO   !! NUMBER OF LONGITUDES.
INTEGER :: N_LAT_CO   !! NUMBER OF LATITUDES.

! ---------------------------------------------------------------------------- !
! 
!    1. RE-FORMAT INPUT PARAMETERS.
!       --------------------------

WEST_CO  =  READ_COOR_TEXT(WEST)
SOUTH_CO =  READ_COOR_TEXT(SOUTH)
IF (PRESENT(EAST)) THEN
   EAST_CO =  READ_COOR_TEXT(EAST)
ELSE
   EAST_CO = COOR_UNDEF
END IF
IF (PRESENT(NORTH)) THEN
   NORTH_CO =  READ_COOR_TEXT(NORTH)
ELSE
   NORTH_CO = COOR_UNDEF
END IF
IF (PRESENT(D_LON)) THEN
   D_LON_CO =  READ_COOR_TEXT(D_LON)
ELSE
   D_LON_CO = COOR_UNDEF
END IF
IF (PRESENT(D_LAT)) THEN
   D_LAT_CO =  READ_COOR_TEXT(D_LAT)
ELSE
   D_LAT_CO = COOR_UNDEF
END IF

IF (PRESENT(N_LON)) THEN
   N_LON_CO =  N_LON
ELSE
   N_LON_CO = -1
END IF

IF (PRESENT(N_LAT)) THEN
   N_LAT_CO =  N_LAT
ELSE
   N_LAT_CO = -1
END IF

! ---------------------------------------------------------------------------- !
! 
!    2. TRANSFER INPUT PARAMETERS.
!       --------------------------

IF (PRESENT(CODE)) THEN
   CALL SET_TOPO_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,         &
&                          NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,     &
&                          N_LON=N_LON_CO, N_LAT=N_LAT_CO, CODE=CODE)
ELSE
   CALL SET_TOPO_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,         &
&                          NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,     &
&                          N_LON=N_LON_CO, N_LAT=N_LAT_CO)
END IF

END SUBROUTINE SET_TOPO_HEADER_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_TOPO_HEADER_D (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,          &
&                             N_LON, N_LAT, CODE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_TOPO_HEADER - SET INPUT GRID FOR TOPO FIELDS IN WAM_TOPO MODULE.       !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2011                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS FOR TOPO FIELDS TO WAM_TOPO_MODULE.!
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY DEFINITIONS.                                !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL (KIND=KIND_D), INTENT(IN)           :: WEST   !! WEST LONG. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN)           :: SOUTH  !! SOUTH LAT. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: EAST   !! EAST LONG. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: NORTH  !! NORTH LAT.F GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: D_LON  !! LONG. INC. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: D_LAT  !! LAT.  INC. OF GRID [DEG].
INTEGER,            INTENT(IN), OPTIONAL :: N_LON  !! NUMBER OF LONGITUDES.
INTEGER,            INTENT(IN), OPTIONAL :: N_LAT  !! NUMBER OF LATITUDES.
INTEGER,            INTENT(IN), OPTIONAL :: CODE   !! 1 = TOTAL WATER DEPTH (DEFAULT)
                                                   !! OTHERWISE ELEVATION OVER NN.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: WEST_CO    !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER :: SOUTH_CO   !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER :: EAST_CO    !! EAST LONGITUDE OF GRID [M_SEC].
INTEGER :: NORTH_CO   !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER :: D_LON_CO   !! LONGITUDE INCREMENT [M_SEC].
INTEGER :: D_LAT_CO   !! LATITUDE INCREMENT [M_SEC].
INTEGER :: N_LON_CO   !! NUMBER OF LONGITUDES.
INTEGER :: N_LAT_CO   !! NUMBER OF LATITUDES.

! ---------------------------------------------------------------------------- !
! 
!    1. RE-FORMAT INPUT PARAMETERS.
!       --------------------------

WEST_CO  =  DEG_TO_M_SEC(WEST)
SOUTH_CO =  DEG_TO_M_SEC(SOUTH)
IF (PRESENT(EAST)) THEN
   EAST_CO =  DEG_TO_M_SEC(EAST)
ELSE
   EAST_CO = COOR_UNDEF
END IF
IF (PRESENT(NORTH)) THEN
   NORTH_CO =  DEG_TO_M_SEC(NORTH)
ELSE
   NORTH_CO = COOR_UNDEF
END IF
IF (PRESENT(D_LON)) THEN
   D_LON_CO =  DEG_TO_M_SEC(D_LON)
ELSE
   D_LON_CO = COOR_UNDEF
END IF
IF (PRESENT(D_LAT)) THEN
   D_LAT_CO =  DEG_TO_M_SEC(D_LAT)
ELSE
   D_LAT_CO = COOR_UNDEF
END IF

IF (PRESENT(N_LON)) THEN
   N_LON_CO =  N_LON
ELSE
   N_LON_CO = -1
END IF

IF (PRESENT(N_LAT)) THEN
   N_LAT_CO =  N_LAT
ELSE
   N_LAT_CO = -1
END IF

! ---------------------------------------------------------------------------- !
! 
!    2. TRANSFER INPUT PARAMETERS.
!       --------------------------

IF (PRESENT(CODE)) THEN
   CALL SET_TOPO_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,         &
&                          NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,     &
&                          N_LON=N_LON_CO, N_LAT=N_LAT_CO, CODE=CODE)
ELSE
   CALL SET_TOPO_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,         &
&                          NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,     &
&                          N_LON=N_LON_CO, N_LAT=N_LAT_CO)
END IF

END SUBROUTINE SET_TOPO_HEADER_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_TOPO_HEADER_M (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,          &
&                             N_LON, N_LAT, CODE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_TOPO_HEADER - SET INPUT GRID FOR TOPO FIELDS IN WAM_TOPO MODULE.       !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS FOR TOPO FIELDS TO WAM_TOPO_MODULE.!
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY DEFINITIONS.                                !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,  INTENT(IN)           :: WEST    !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER,  INTENT(IN)           :: SOUTH   !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER,  INTENT(IN), OPTIONAL :: EAST    !! EAST LONGITUDE OF GRID [M_SEC].
INTEGER,  INTENT(IN), OPTIONAL :: NORTH   !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER,  INTENT(IN), OPTIONAL :: D_LON   !! LONGITUDE INCREMENT [M_SEC].
INTEGER,  INTENT(IN), OPTIONAL :: D_LAT   !! LATITUDE INCREMENT [M_SEC].
INTEGER,  INTENT(IN), OPTIONAL :: N_LON   !! NUMBER OF LONGITUDES.
INTEGER,  INTENT(IN), OPTIONAL :: N_LAT   !! NUMBER OF LATITUDES.
INTEGER,  INTENT(IN), OPTIONAL :: CODE    !! 1 = TOTAL WATER DEPTH (DEFAULT)
                                          !! OTHERWISE ELEVATION OVER NN.
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

LOGICAL  :: ERROR                         !! ERROR FLAG

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLEAR GRID DEFINITIONS.                                               !
!        -----------------------                                               !

NX_IN    = -1           !! NUMBER OF COLUMNES IN TOPO INPUT GRID.
NY_IN    = -1           !! NUMBER OF ROWS      IN TOPO INPUT GRID.
PER      =.FALSE.       !! .TRUE. IF PERIODIC GRID.
CODE_IN  = 1            !! TOPO CODE IS TOTAL WATER DEPTH
DX_IN    = COOR_UNDEF   !! STEPSIZE BETWEEN LONGITUDES.
DY_IN    = COOR_UNDEF   !! STEPSIZE BETWEEN LATITUDES.
SOUTH_IN = COOR_UNDEF   !! MOST SOUTHERN LATITUDE.
NORTH_IN = COOR_UNDEF   !! MOST NORTHERN LATITUDE.
WEST_IN  = COOR_UNDEF   !! LEFT MOST LONGITUDE.
EAST_IN  = COOR_UNDEF   !! RIGHT MOST LONGITUDE.
EQUAL_GRID = .FALSE.    !! INPUT GRID IS NOT EQUAL TO MODEL GRID.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COPY GRID DEFINITION.                                                 !
!        ---------------------                                                 !

WEST_IN  = WEST
SOUTH_IN = SOUTH
IF (PRESENT(EAST)) EAST_IN  = EAST
IF (PRESENT(NORTH)) NORTH_IN = NORTH
IF (PRESENT(D_LON)) DX_IN    = D_LON
IF (PRESENT(D_LAT)) DY_IN    = D_LAT
IF (PRESENT(N_LON)) NX_IN    = N_LON
IF (PRESENT(N_LAT)) NY_IN    = N_LAT
IF (PRESENT(CODE))  CODE_IN = CODE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CHECK GRID.                                                           !
!        -----------                                                           !

CALL CHECK_GRID_DEFINITION (WEST_IN, SOUTH_IN, EAST_IN, NORTH_IN,              &
&                           DX_IN, DY_IN, NX_IN, NY_IN, ERROR)

IF (ERROR) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *        FATAL ERROR IN SUB. SET_TOPO_HEADER             *'
   WRITE (IU06,*) ' *        ===================================             *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *       A TOPO GRID COULD NOT BE DEFINED.                *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF


PER = PERIODIC (WEST_IN, EAST_IN, DX_IN, NX_IN)   !! PERIODIC?
EQUAL_GRID = EQUAL_TO_M_GRID (NX_IN, NY_IN, DX_IN, DY_IN,                      &
&                             WEST_IN, SOUTH_IN, EAST_IN, NORTH_IN)

END SUBROUTINE SET_TOPO_HEADER_M

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_TOPO_TIMESTEPS (IN, OUT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_TOPO_TIMESTEPS - SET TOPO TIMESTEPS IN WAM_TOPO_MODULE.                !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE TOPO TIMESTEPS TO WAM_TOPO_MODULE.                        !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY DEFINITIONS.                                !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)           :: IN    !! TOPO INPUT TIME STEP.
INTEGER, INTENT(IN), OPTIONAL :: OUT   !! TOPO OUTPUT TIME STEP.

! ---------------------------------------------------------------------------- !

IDELTI = MAX (IN,0)
IF (PRESENT(OUT)) THEN
   IDELTO = OUT
ELSE
   IDELTO = IN
END IF

IF (IDELTI.LE.0) THEN    !! STATIONARY DEPTH.
   IDELTO = 0
   RETURN
END IF

IF (IDELTO.LE.0) THEN    !! DEPTH IS NOT INTERPOLATED.
   IDELTO = IDELTI
END IF

!         CHECK TOPO INPUT AND TOPO OUTPUT TIMESTEP.                           !
!         ------------------------------------------                           !

IF (IDELTI.LT.IDELTO) THEN
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+       WARNING ERROR IN SUB. SET_TOPO_TIMESTEPS          +'
   WRITE(IU06,*) '+       ========================================          +'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+ TOPO INPUT TIME STEP IS LESS THAN TOPO OUTPUT STEP      +'
   WRITE(IU06,*) '+ TOPO INPUT TIMESTEP    IN = ', IDELTI, ' SECONDS'
   WRITE(IU06,*) '+ TOPO OUTPUT TIMESTEP  OUT = ', IDELTO, ' SECONDS'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+ TOPO INPUT CHANGED TO TOPO OUTPUT TIMESTEP              +'
   WRITE(IU06,*) '+ THE MODEL MAY IGNORE DEPTH FIELDS                       +'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

IF (MOD(IDELTI,IDELTO).NE.0) THEN
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+       WARNING ERROR IN SUB. SET_TOPO_TIMESTEPS          +'
   WRITE(IU06,*) '+       ========================================          +'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+ TOPO INPUT AND TOPO OUTPUT TIME STEPS DO NOT HAVE       +'
   WRITE(IU06,*) '+ INTEGER RATIO.                                          +'
   WRITE(IU06,*) '+ TOPO INPUT TIMESTEP    IN = ', IDELTI, ' SECONDS'
   WRITE(IU06,*) '+ TOPO OUTPUT TIMESTEP  OUT = ', IDELTO, ' SECONDS'
   WRITE(IU06,*) '+ OUTPUT TIMESTEP IS CHANGED TO NEAREST MULTIPLE OF THE   +'
   WRITE(IU06,*) '+ INPUT TIME STEP.                                        +'
   IDELTO = IDELTI/MAX(NINT(REAL(IDELTI)/REAL(IDELTO)),1)
   WRITE(IU06,*) '+ NEW OUTPUT TIMESTEP   IS: ',IDELTO, ' SECONDS'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

END SUBROUTINE SET_TOPO_TIMESTEPS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WAM_TOPO (DT, CD_START)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   WAM_TOPO - ROUTINE TO READ AND PROCESS ONE TOPO FIELD.                     !
!                                                                              !
!     H. GUNTHER      GKSS    DECEMBER 2009                                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       READ, INTERPOLATE AND CONVERT INPUT TOPOGRAPHY FIELDS TO WAM DEPTHS.   !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        READ A TOPO FIELD FROM THE TOPO FILE (SEARCH FOR IT),                 !
!       INTERPOLATED IT TO THE WAVE MODEL SEA POINTS AND                       !
!        CALCULATES THE TOTAL WATER DEPTH IF SURFACE ELEVATION OVER NN.        !
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

REAL,          INTENT(OUT)    :: DT(:)    !! DEPTH (M).
CHARACTER (LEN=14),INTENT(IN) :: CD_START !! DATE OF FIELD TO BE LOOKED FOR.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. READ TOPO DATA AND CHECK DATE.                                        !
!        ------------------------------                                        !

DO
   CALL READ_TOPO_INPUT

   IF (CD_READ.EQ.CD_START) EXIT
   IF (CD_READ.GT.CD_START) THEN
         WRITE (IU06,*) ' ******************************************'
         WRITE (IU06,*) ' *                                        *'
         WRITE (IU06,*) ' *      FATAL ERROR SUB. WAM_TOPO         *'
         WRITE (IU06,*) ' *      =========================         *'
         WRITE (IU06,*) ' * TOPO DATE READ IS LATER THAN EXPECTED  *'
         WRITE (IU06,*) ' * DATE READ IS      CD_READ = ', CD_READ
         WRITE (IU06,*) ' * DATE EXPECTED IS CD_START = ', CD_START
         WRITE (IU06,*) ' *                                        *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
         WRITE (IU06,*) ' *                                        *'
         WRITE (IU06,*) ' ******************************************'
         CALL ABORT1
   END IF
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTERPOLATE TOPO FIELD TO GRID.                                       !
!        -------------------------------                                       !

DT = 0.  !! INITIALISE TOPO ARRAY WITH ZERO.
IF (EQUAL_GRID) THEN
   DT = PACK(TOPO_IN, L_S_MASK)
ELSE
   CALL INTERPOLATION_TO_GRID (IU06, PER, DX_IN, DY_IN, WEST_IN, SOUTH_IN,     &
&                              TOPO_IN, DT)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. PROCESS TOPO ACCORDING TO TYPE                                        !
!        NOTHING TO DO FOR TOTAL WATER DEPTH (CODE_IN = 1).                    !
!        ------------------------------------------------                      !

IF (CODE_IN.NE.1)  DT = DT + DEPTH_B  !! INPUT IS SURFACE ELEVATION OVER NN.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. TEST OUTPUT OF WAVE MODEL BLOCKS                                      !
!        ---------------------------------                                     !

IF (ITEST.GE.3) THEN
   IJ = MIN (10,SIZE(DT))
   WRITE (IU06,*) ' '
   WRITE (IU06,*) '      SUB. WAM_TOPO: TOPO FIELDS CONVERTED TO MODEL GRID'
   WRITE (IU06,*) ' '
   WRITE (IU06,*) ' DT(1:10) = ', DT(1:IJ)
END IF

END SUBROUTINE WAM_TOPO

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVAT MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE NOTIM (CD_START, CD_END)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   NOTIM - STEERING SUB IF TIME INTERPOLATION IS NOT WANTED.                  !
!                                                                              !
!     H. GUNTHER    ECMWF    MAY 1990     MODIFIED FOR SUB VERSION.            !
!     H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.                !
!     H. GUNTHER      GKSS    SEPTEMBER 2000   FT90.                           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       NOTIM NO TIME INTERPOLATION: PROCESSES TOPO FIELDS.                    !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NO TIME INTERPOLATION:                                                 !
!       TOPO FIELDS ARE PROCESSED EVERY IDELTI SECONDS,                        !
!       THE TOPO IS INTERPOLATED IN SPACE ONLY.                                !
!       DEPTH FIELDS AND SAVED IN WAM_TOPO_MODULE.                             !
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

CHARACTER (LEN=14), INTENT(IN)  :: CD_START  !! DATE OF FIRST TOPO FIELD.
CHARACTER (LEN=14), INTENT(IN)  :: CD_END    !! DATE OF LAST TOPO FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER            :: MP
CHARACTER (LEN=14) :: CDTTOH
REAL               :: DT(1:NSEA)  !! OUTPUT TOPO FIELD ARRAY.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. LOOP OVER OUTPUT DEPTH TIMES.                                         !
!        -----------------------------                                         !

CDTTOH = CD_START
MP = 0
DO WHILE (CDTTOH.LE.CD_END)

!     1.1 READ ONE TOPO FIELD AND TRANSFORM TO GRID.                           !
!         --------------------------------------------                         !

   CALL WAM_TOPO (DT, CDTTOH)
   MP = MP + 1

!     1.2 SAVE IN MODULE WAM_TOPO.                                             !
!         ------------------------                                             !

   D_STORE(:,MP)  = DT
   CD_STORE(MP) = CDTTOH

   IF (ITEST.GE.3) THEN
      WRITE(IU06,*) '       SUB. NOTIM: NEW TOPO FILES AT CDTTOH = ', CDTTOH
   END IF

!     1.3 UPDATE TOPO FIELD REQUEST TIME.                                      !
!         -------------------------------                                      !

   CALL INCDATE (CDTTOH,IDELTO)
END DO

END SUBROUTINE NOTIM

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TIMIN (CD_START, CD_END)

! ---------------------------------------------------------------------------- !
!                                                                              !
!                                                                              !
!   TIMIN - STEERING MODULE IF TIME INTERPOLATION IS WANTED.                   !
!                                                                              !
!     H. GUNTHER    ECMWF    MAY 1990         MODIFIED FOR SUB VERSION.        !
!     H. GUNTHER    ECMWF    DECEMBER 90      MODIFIED FOR CYCLE_4.            !
!     H. GUNTHER    GKSS     SEPTEMBER 2000   FT90.                            !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TIME INTERPOLATION: PROCESSES TOPO FIELDS.                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       TOPO FIELDS ARE READ IN EVERY IDELTI SECONDS,                          !
!       BI-LINEAR INTERPOLATED IN SPACE, AND TRANSFORMED TO WAM GRID.          !
!       TOTAL DEPTH IS INTERPOLATED LINEARLY IN TIME.                          !
!       DEPTH FIELDS AND SAVED IN WAM_TOPO_MODULE.                             !
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

CHARACTER (LEN=14), INTENT(IN)  :: CD_START  !! DATE OF FIRST TOPO FIELD.
CHARACTER (LEN=14), INTENT(IN)  :: CD_END    !! DATE OF LAST TOPO FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER            :: MP, NTS, N
REAL               :: DEL
CHARACTER (LEN=14) :: CDT1, CDT2, CDTH

REAL               :: DT1(1:NSEA)  !! OUTPUT TOPO FIELD ARRAY.
REAL               :: DT2(1:NSEA)  !! OUTPUT TOPO FIELD ARRAY.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALIZE TIME COUNTER AND FIRST DEPTH FIELD USED FOR INTERPOLATION. !
!        --------------------------------------------------------------------- !

CDT1 = CDTA
DT1  = DEPTH
CDTH = CDT1
CALL INCDATE (CDTH,IDELTO)
IF (CD_START.NE.CDTH) THEN
   WRITE(IU06,*) ' *******************************************'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *        FATAL ERROR IN --TIMIN--         *'
   WRITE(IU06,*) ' *        =========================        *'
   WRITE(IU06,*) ' * DATES DO NOT MATCH.                     *'
   WRITE(IU06,*) ' * START DATE FOR TOPO IS     CD_START = ', CD_START
   WRITE(IU06,*) ' * LAST DATE SAVED IN MODULE IS   CDT1 = ', CDT1
   WRITE(IU06,*) ' * PROCESSING WILL BE ABORTED              *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *******************************************'
   CALL ABORT1
END IF

CDT2 = CDT1
CALL INCDATE(CDT2,IDELTI)
NTS = IDELTI/IDELTO
DEL = REAL(IDELTO)/REAL(IDELTI)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LOOP OVER INPUT TOPO FIELDS.                                          !
!        ----------------------------                                          !

MP = 0
DO

!     2.1 READ ONE TOPO FIELD AND TRANSFORM TO GRID.                           !
!         ------------------------------------------                           !

   CDT2 = CDT1
   CALL INCDATE(CDT2,IDELTI)
   CALL WAM_TOPO (DT2, CDT2)

!     2.2 INTERPOLATE AND SAVE DEPTH FIELDS.                                   !
!         ----------------------------------                                   !

   CDTH = CDT1
   DO N = 1,NTS
      MP = MP + 1
      CALL INCDATE(CDTH,IDELTO)
      CD_STORE(MP) = CDTH
      D_STORE(:,MP) = DT1 + REAL(N)*DEL*(DT2-DT1)
   END DO
   IF (ITEST.GE.3) THEN
      WRITE(IU06,*) '       SUB. TIMIN: TOPO FIELDS FOR CD_STORE = ',          &
&                                           CD_STORE(MP-NTS+1:MP)
   END IF

!     2.3 UPDATE TOPO FIELD REQUEST TIME AND READ NEXT IF REQUESTED.           !
!         ----------------------------------------------------------           !

   DT1  = DT2
   CDT1 = CDT2
   CALL INCDATE (CDTH,IDELTI)
   IF (CDTH.GT.CD_END) EXIT
END DO

END SUBROUTINE TIMIN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_TOPO_MODULE

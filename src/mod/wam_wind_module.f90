MODULE WAM_WIND_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE STORES THE WIND INPUT GRID SPECFICATIONS, THE WIND FIELDS      !
!   WHICH ARE PASSED FROM PREWIND TO WAM-MODELL.                               !
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

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST, wpath, area
USE WAM_GENERAL_MODULE, ONLY: PI, ZPI, G, RAD, ROAIR, XKAPPA
USE WAM_GRID_MODULE,    ONLY: NSEA, L_S_MASK
USE WAM_MODEL_MODULE,   ONLY: U10, UDIR
USE WAM_TIMOPT_MODULE,  ONLY: CDA, CDATEE, IDEL_WAM
use wam_mpi_module,     only: nijs, nijl
use wam_special_module, only: readyf

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INPUT WIND GRID SPECFICATIONS DATE AND WIND FIELDS.                   !
!        ---------------------------------------------------                   !

INTEGER :: NX_IN  =-1      !! NUMBER OF LONGITUDES.
INTEGER :: NY_IN  =-1      !! NUMBER OF LATITUDES.
LOGICAL :: PER =.FALSE.    !! .TRUE. IF PERIODIC GRID.
INTEGER :: CODE_IN = 3     !! WIND CODE:
                           !! 1 = USTAR;  2 = USTRESS; 3 = U10
INTEGER :: DX_IN  =-1.     !! STEPSIZE BETWEEN LONGITUDES [M_SEC].
INTEGER :: DY_IN  =-1.     !! STEPSIZE BETWEEN LATITUDES [M_SEC].
INTEGER :: SOUTH_IN =-1.   !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER :: NORTH_IN =-1.   !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER :: WEST_IN =-1.    !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER :: EAST_IN =-1.    !! EAST LONGITUDE OF GRID [M_SEC].
LOGICAL :: EQUAL_GRID =.FALSE. !! .TRUE. IF WIND GRID IS EQUAL TO MODEL GRID.
    
CHARACTER (LEN= 14) :: CD_READ =' '!! DATE OF LAST DATA READ FROM INPUT.
   
REAL, ALLOCATABLE, DIMENSION(:,:)  :: U_IN  !! W-E WIND COMPONENT.
REAL, ALLOCATABLE, DIMENSION(:,:)  :: V_IN  !! S-N WIND COMPONENT.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. WIND TIMESTEPS.                                                       !
!        ---------------                                                       !

INTEGER, PUBLIC   :: IDELWI = -1  !! INPUT WIND TIMESTEP INTO IN SECONDS.
INTEGER, PUBLIC   :: IDELWO = -1  !! OUTPUT WIND TIMESTEP IN SECONDS
                                  !! EQUAL TO INPUT TIMESTEP INTO WAMODEL.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. WIND FIELDS PREPARED BY PREWIND FOR WAM-MODEL.                        !
!        (FIRST INDEX IS POINTS, SECOND IS TIME)                               !
!        ----------------------------------------------                        !

INTEGER  :: M_STORE = 0     !! NUMBER OF WINDFIELDS.

REAL,               ALLOCATABLE, DIMENSION(:,:) :: U_STORE  !! WIND SPEEDS.
REAL,               ALLOCATABLE, DIMENSION(:,:) :: D_STORE  !! WIND DIRECTIONS.
CHARACTER (LEN=14), ALLOCATABLE, DIMENSION(:)   :: CD_STORE !! DATE/TIME.

integer :: lent

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE GET_WIND                       !! GETS WINDS FORM THIS MODULE.
   MODULE PROCEDURE GET_WIND
END INTERFACE
PUBLIC GET_WIND

INTERFACE PREPARE_WIND                   !! PREPARES WIND DATA FOR WAVE MODEL.
   MODULE PROCEDURE PREPARE_WIND
END INTERFACE
PUBLIC PREPARE_WIND

INTERFACE PRINT_WIND_STATUS              !! PRINTS WIND STATUS.
   MODULE PROCEDURE PRINT_WIND_STATUS
END INTERFACE
PUBLIC PRINT_WIND_STATUS

INTERFACE SET_WIND_FIELD                 !! SETS WIND FIELD.
   MODULE PROCEDURE SET_WIND_FIELD
END INTERFACE
PUBLIC SET_WIND_FIELD

INTERFACE SET_WIND_HEADER                !! SETS WIND HEADER.
   MODULE PROCEDURE SET_WIND_HEADER_C    !! CHARACTER VERSION
   MODULE PROCEDURE SET_WIND_HEADER_D    !! DEGREE VERSION
   MODULE PROCEDURE SET_WIND_HEADER_M    !! M_SEC VERSION
END INTERFACE
PUBLIC SET_WIND_HEADER

INTERFACE SET_WIND_TIMESTEPS             !! SETS WIND TIMESTEPS.
   MODULE PROCEDURE SET_WIND_TIMESTEPS
END INTERFACE
PUBLIC SET_WIND_TIMESTEPS

interface set_ready_file_flag            !! wait for ready files ?
   module procedure set_ready_file_flag
end interface
public set_ready_file_flag
    
interface set_ready_file_directory       !! sets full path of ready file directory
   module procedure set_ready_file_directory
end interface
public set_ready_file_directory

INTERFACE WAM_WIND       !! READS AND TRANSFORMS INPUT WINDS TO WAM POINTS.
   MODULE PROCEDURE WAM_WIND
END INTERFACE
PUBLIC WAM_WIND

INTERFACE
   SUBROUTINE READ_WIND_INPUT            !! READS A WIND FIELD
   END SUBROUTINE READ_WIND_INPUT
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

SUBROUTINE GET_WIND (CD_NEW)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   GET_WIND - TRANSFERS NEW WIND DATA TO WAM MODEL.                           !
!                                                                              !
!     P.A.E.M. JANSSEN  KNMI/ECMWF  SEPTEMBER 1994                             !
!     A. Behrens   MSC/GKSS  December 2003  Message passing                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER WINDS FROM WAM_WIND_MODULE TO THE WAM_MODEL_MODULE.           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       COPY NEW WINDS.                                                        !
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

CHARACTER (LEN=14),INTENT(INOUT) :: CD_NEW !! DATE TO USE OF NEW WIND FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER   :: IT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NEW WINDS FROM TRANSFER ARRAYS.                                       !
!        -------------------------------                                       !

CDA  = ' '
DO IT = 1, M_STORE
   IF (CD_STORE(IT).GT.CD_NEW) THEN
      U10(:)  = U_STORE(:,IT)
      UDIR(:) = D_STORE(:,IT)
      CDA  = CD_STORE(IT)
      EXIT
   END IF
END DO

IF (CDA .EQ. ' ') THEN
   WRITE(IU06,*) ' *******************************************'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *   FATAL ERROR IN SUB. GET_WIND          *'
   WRITE(IU06,*) ' *   ============================          *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' * WIND FIELD WAS NOT FOUND                *'
   WRITE(IU06,*) ' *  CHANGE DATE IS ............ CD_NEW: ', CD_NEW
   WRITE(IU06,*) ' *  LAST DATE IN MODULE IS ...........: ', CD_STORE(M_STORE)
   WRITE(IU06,*) ' *  NO. OF DATES IN MODULE IS M_STORE = ', M_STORE
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *******************************************'
   CALL ABORT1
END IF

CALL INCDATE(CD_NEW,IDELWO)

IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '  '
   WRITE (IU06,*) '   SUB. GET_WIND: NEW WIND FIELDS PASSED TO WAM MODEL'
   WRITE (IU06,*) '     DATE OF WIND FIELD IS .......... CDA = ', CDA
   WRITE (IU06,*) '     DATE FOR NEXT WIND CHANGE IS  CD_NEW = ', CD_NEW
   WRITE (IU06,*) '  '
END IF

END SUBROUTINE GET_WIND

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_WIND

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PREPARE_WIND - PREPARES WIND DATA FOR WAVE MODEL.                          !
!                                                                              !
!     P.GROENWOUD     DELFT HYDRAULICS LABORATORY  OKTOBER 1986                !
!                                                                              !
!     E. BAUER        MPI         FEB 1987   VERSION FOR CDC 205 HAMBURG.      !
!                                                                              !
!     S. HASSELMANN   MPI         MAY 1987   COMBINED CDC 205 AND CRAY         !
!                                            VERSIONS.                         !
!     W. BRUEGGEMANN  MPI      AUGUST 1988   SIMPLIFIED PROGRAM.               !
!                                                                              !
!     L. ZAMBRESKY    ECMWF      JUNE 1988   MODIFIED EXTENSIVELY FOR          !
!                                            COUPLING TO SPECTRAL MODEL.       !
!                                                                              !
!     H. GUNTHER      ECMWF      JUNE 1990   MODIFIED FOR CYCLE_4.             !
!     H. GUNTHER      GKSS  SEPTEMBER 2000   FT90.                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       EVALUATE WIND DATA AT WAVE MODEL GRID POINTS.                          !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INPUT WIND FIELDS WHICH CAN BE COMPONENTS OR MAGNITUDE AND DIRECTION   !
!       OF USTAR, U10, USTRESS ARE TRANSFORMED TO FRICTION VELOCITIES.         !
!       THE INPUT FIELDS HAVE TO BE ON A LAT /LONG GRID.                       !
!       SEE SUB READ_WIND FOR FORMATS AND HEADER INFORMATION, WHICH HAVE TO BE !
!       GIVEN TO THE PROGRAM.                                                  !
!                                                                              !
!       A DOUBLE LINEAR INTERPOLATION IN SPACE IS PERFORMED ONTO THE MODEL     !
!       BLOCKS.                                                                !
!       IF THE WIND OUTPUT TIMSTEP IS LESS THAN THE INPUT TIMESTEP             !
!       A LINEAR INTERPOLATION IN TIME IS PERFORMED.                           !
!       ALL WIND FIELDS ARE STORED IN WAM_WIND_MODULE.                         !
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

CHARACTER (LEN=14) :: CD_START, CD_END 

! ---------------------------------------------------------------------------- !
!                                                                              !
!                                                                              !
!     1. BEGIN AND END DATES OF WIND FIELDS TO BE PROCESSED.                   !
!        ---------------------------------------------------                   !

CD_START = CDA
CALL INCDATE (CD_START,IDELWO)
CD_END = CDA
CALL INCDATE (CD_END, IDEL_WAM)
IF (CD_END.GE.CDATEE) CD_END = CDATEE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. NUMBER OF WIND FIELDS TO BE GENERATED.                                !
!        --------------------------------------                                !

M_STORE = IDEL_WAM/IDELWO

IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '  '
   WRITE (IU06,*) '   SUB. PREPARE_WIND: WIND REQUEST'
   WRITE (IU06,*) '     START OF PERIOD IS ............ CD_START = ', CD_START
   WRITE (IU06,*) '     END   OF PERIOD IS .............. CD_END = ', CD_END
   WRITE (IU06,*) '     WIND INPUT TIME STEP ............ IDELWI = ', IDELWI
   WRITE (IU06,*) '     WIND OUTPUT TIME STEP ........... IDELWO = ', IDELWO
   WRITE (IU06,*) '     NUMBER FIELDS TO BE GENERATED IS M_STORE = ', M_STORE
   WRITE (IU06,*) '     FIELDS ARE SAVED IN WAM_WIND_MODULE'
END IF

IF (ALLOCATED(CD_STORE))  then
   if (M_STORE.NE.SIZE(CD_STORE)) DEALLOCATE (CD_STORE, U_STORE, D_STORE)
endif
IF (.NOT. ALLOCATED(U_STORE) ) ALLOCATE (U_STORE(nijs:nijl,1:M_STORE))
IF (.NOT. ALLOCATED(D_STORE) ) ALLOCATE (D_STORE(nijs:nijl,1:M_STORE))
IF (.NOT. ALLOCATED(CD_STORE)) ALLOCATE (CD_STORE(1:M_STORE))
CD_STORE = ' '

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. PROCESS WIND FIELDS.                                                  !
!        --------------------                                                  !

IF (IDELWO.GE.IDELWI) THEN

   IF (ITEST.GE.2) WRITE (IU06,*) '     NO TIME INTERPOLATION'
   CALL NOTIM (CD_START, CD_END)

ELSE

   IF (ITEST.GE.2) WRITE (IU06,*) '     TIME INTERPOLATION'
   CALL TIMIN (CD_START, CD_END)

END IF

END SUBROUTINE PREPARE_WIND

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_WIND_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PRINT_WIND_STATUS - PRINTS STATUS OF THE WIND MODULE.                      !
!                                                                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2009                           !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       MAKE A PRINTER OUTPUT OF THE STATUS OF THE WIND MODULE.                !
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
WRITE (IU06,*) '              WIND MODULE STATUS:'
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '  '
WRITE (IU06,*) ' MODEL WIND INPUT TIME STEP .................: ', IDELWI,' SECS'
WRITE (IU06,*) ' MODEL WIND OUTPUT TIME STEP.................: ', IDELWO,' SECS'
IF (IDELWO.GT.0 .AND. IDELWI.GT.0) THEN
   IF (IDELWO.GE.IDELWI) THEN
      WRITE(IU06,*) ' WIND FIELDS ARE NOT INTERPOLATED IN TIME'
   ELSE
      WRITE(IU06,*) ' WIND FIELDS ARE INTERPOLATED IN TIME'
   END IF
END IF
WRITE (IU06,*) '  '
WRITE (IU06,*) ' DATE OF LAST WIND FIELD GIVEN TO MODULE.....: ', CD_READ
WRITE (IU06,*) ' DATE OF LAST WIND FIELD GIVEN TO WAVE MODEL.: ', CDA
WRITE (IU06,*) '  '

IF (NX_IN.GT.0 .AND. NY_IN.GT.0) THEN
   WRITE (IU06,'('' WIND INPUT GRID SPECIFICATION ARE: '')')
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
   WRITE (IU06,*) ' WIND INPUT GRID IS NOT DEFINED'
END IF
WRITE (IU06,*) ' WIND CODE: 1 USTAR;  2 USTRESS; 3 U10.......: ', CODE_IN
WRITE (IU06,*) '  '

IF ( M_STORE.GT.0 .AND. ANY(CD_STORE.NE.' ')) THEN
   WRITE (IU06,*) ' NUMBER OF WIND FIELDS STORED IN MODULE .....: ', M_STORE
   WRITE (IU06,*) ' DATES OF WIND FIELDS ARE:'
   WRITE (IU06,'(5(3X,A14,2X))')  CD_STORE
ELSE
   WRITE (IU06,*) ' WIND FIELDS ARE NOT STORED IN MODULE'
END IF
WRITE (IU06,*) '  '

END SUBROUTINE PRINT_WIND_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_WIND_FIELD (CDT, U_MAP, V_MAP, CODE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_WIND_FIELD - SET INPUT WIND FIELD IN WAM_WIND MODULE.                  !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER WINDS TO WAM_WIND_MODULE.                                     !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       CHECK CONSISTENCY AND COPY NEW WINDS. IF WINDS ARE GIVEN AS MAGNITUDE  !
!       AND DIRECTION, THEY ARE CHANGED TO VECTOR COMPONENTS.                  !
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

CHARACTER (LEN=14),INTENT(IN) :: CDT        !! DATE/TIME OF WIND FIELD.
REAL,              INTENT(IN) :: U_MAP(:,:) !! U-COMP. OR SPEED OF WINDS [M/S].
REAL,              INTENT(IN) :: V_MAP(:,:) !! V-COMP. [M/S] OR DIRECTION [DEG]
INTEGER, INTENT(IN), OPTIONAL :: CODE       !! 1 = SPEED AND DIRECTION
                                            !! OTHERWISE: COMPONENTS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COPY DATE.                                                            !
!        ----------                                                            !

CD_READ = CDT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK CONSISTENCY.                                                    !
!        ------------------                                                    !

IF (NX_IN.NE.SIZE(U_MAP,1) .OR. NY_IN.NE.SIZE(U_MAP,2) .OR.                    &
&   NX_IN.NE.SIZE(V_MAP,1) .OR. NY_IN.NE.SIZE(V_MAP,2) ) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *          FATAL  ERROR IN SUB. SET_WIND_FIELD           *'
   WRITE (IU06,*) ' *          ===================================           *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' * GRID SPECIFICATIONS ARE NOT CONSISTENT                 *'
   WRITE (IU06,*) ' * OR WIND HEADER IS NOT DEFINED:                         *'
   WRITE (IU06,*) ' * GRID SIZES AS DEFINED IN MODULE ARE:                   *'
   WRITE (IU06,*) ' * NO. OF LONGITUDES    NX_IN = ', NX_IN
   WRITE (IU06,*) ' * NO. OF LATITUDE      NY_IN = ', NY_IN
   WRITE (IU06,*) ' * DIMENSIONS OF U WIND ARRAY ARE :                       *'
   WRITE (IU06,*) ' * NO. OF LONGITUDES 1. DIMENSION = ', SIZE(U_MAP,1)
   WRITE (IU06,*) ' * NO. OF LATITUDE   2. DIMENSION = ', SIZE(U_MAP,2)
   WRITE (IU06,*) ' * DIMENSIONS OF V WIND ARRAY ARE :                       *'
   WRITE (IU06,*) ' * NO. OF LONGITUDES 1. DIMENSION = ', SIZE(V_MAP,1)
   WRITE (IU06,*) ' * NO. OF LATITUDE   2. DIMENSION = ', SIZE(V_MAP,2)
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COPY WIND FIELDS.                                                     !
!        -----------------                                                     !

IF (ALLOCATED(U_IN)) DEALLOCATE(U_IN)
IF (ALLOCATED(V_IN)) DEALLOCATE(V_IN)
ALLOCATE (U_IN(NX_IN,NY_IN))
ALLOCATE (V_IN(NX_IN,NY_IN))

IF (PRESENT(CODE)) THEN
   IF (CODE.EQ.1) THEN
      U_IN = U_MAP*SIN(V_MAP*RAD)
      V_IN = U_MAP*COS(V_MAP*RAD)
   ELSE
      U_IN = U_MAP
      V_IN = V_MAP
   END IF
ELSE
   U_IN = U_MAP
   V_IN = V_MAP
END IF

END SUBROUTINE SET_WIND_FIELD

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine set_ready_file_flag (p)

logical, intent(in) :: p       !! ready file flag

readyf = p

end subroutine set_ready_file_flag

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
 
subroutine set_ready_file_directory (name, ar)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

character (len=*), intent(in) :: name     !! full path of ready file directory
character (len=*), intent(in) :: ar       !! name of model_area
   
! ---------------------------------------------------------------------------- !

lent = len_trim (name)
if (lent>=0) wpath = name(1:lent)
lent = len_trim (ar)
if (lent>=0) area = ar(1:lent)

end subroutine set_ready_file_directory

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_WIND_HEADER_C (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,          &
&                             N_LON, N_LAT, CODE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_WIND_HEADER - SET INPUT GRID FOR WIND FIELDS IN WAM_WIND MODULE.       !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS FOR WIND FIELDS TO WAM_WIND_MODULE.!
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
INTEGER, INTENT(IN), OPTIONAL :: CODE    !! WIND CODE:
                                         !!   1 = USTAR,
                                         !!   2 = USTRESS,
                                         !!   3 = U10 (DEFAULT).

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
   CALL SET_WIND_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,         &
&                          NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,     &
&                          N_LON=N_LON_CO, N_LAT=N_LAT_CO, CODE=CODE)
ELSE
   CALL SET_WIND_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,         &
&                          NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,     &
&                          N_LON=N_LON_CO, N_LAT=N_LAT_CO)
END IF

END SUBROUTINE SET_WIND_HEADER_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_WIND_HEADER_D (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,          &
&                             N_LON, N_LAT, CODE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_WIND_HEADER - SET INPUT GRID FOR WIND FIELDS IN WAM_WIND MODULE.       !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS FOR WIND FIELDS TO WAM_WIND_MODULE.!
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
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: NORTH  !! NORTH LAT. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: D_LON  !! LONG. INC. OF GRID [DEG].
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: D_LAT  !! LATI. INC. OF GRID [DEG].
INTEGER,            INTENT(IN), OPTIONAL :: N_LON  !! NUMBER OF LONG. IN GRID.
INTEGER,            INTENT(IN), OPTIONAL :: N_LAT  !! NUMBER OF LAT. IN GRID.
INTEGER,            INTENT(IN), OPTIONAL :: CODE   !! WIND CODE:
                                                   !!   1 = USTAR,
                                                   !!   2 = USTRESS,
                                                   !!   3 = U10 (DEFAULT).

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
   CALL SET_WIND_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,         &
&                          NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,     &
&                          N_LON=N_LON_CO, N_LAT=N_LAT_CO, CODE=CODE)
ELSE
   CALL SET_WIND_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,         &
&                          NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,     &
&                          N_LON=N_LON_CO, N_LAT=N_LAT_CO)
END IF

END SUBROUTINE SET_WIND_HEADER_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_WIND_HEADER_M (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,          &
&                            N_LON, N_LAT, CODE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_WIND_HEADER - SET INPUT GRID FOR WIND FIELDS IN WAM_WIND MODULE.       !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS FOR WIND FIELDS TO WAM_WIND_MODULE.!
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
INTEGER,  INTENT(IN), OPTIONAL :: CODE    !! WIND CODE:
                                          !!   1 = USTAR,
                                          !!   2 = USTRESS,
                                          !!   3 = U10 (DEFAULT).

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

LOGICAL  :: ERROR                         !! ERROR FLAG

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLEAR GRID DEFINITIONS.                                               !
!        -----------------------                                               !

NX_IN  =-1              !! NUMBER OF COLUMNES IN WIND INPUT GRID.
NY_IN  =-1              !! NUMBER OF ROWS     IN WIND INPUT GRID.
PER =.FALSE.            !! .TRUE. IF PERIODIC GRID.
CODE_IN = 3             !! WIND CODE 1 = USTAR;  2 = USTRESS; 3 = U10
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
IF (PRESENT(EAST )) EAST_IN  = EAST
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
   WRITE (IU06,*) ' *         FATAL ERROR IN SUB. SET_WIND_HEADER            *'
   WRITE (IU06,*) ' *         ===================================            *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *       A WIND GRID COULD NOT BE DEFINED.                *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF


PER = PERIODIC (WEST_IN, EAST_IN, DX_IN, NX_IN)   !! PERIODIC?
EQUAL_GRID = EQUAL_TO_M_GRID (NX_IN, NY_IN, DX_IN, DY_IN,                      &
&                             WEST_IN, SOUTH_IN, EAST_IN, NORTH_IN)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CHECK WIND CODE.                                                      !
!        ----------------                                                      !

IF (CODE_IN.LT.1 .OR. CODE_IN.GT.3) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *        FATAL ERROR IN SUB. SET_WIND_HEADER             *'
   WRITE (IU06,*) ' *        ===================================             *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' * WIND CODE OUT OF RANGE. CODE MUST BE:                  *'
   WRITE (IU06,*) ' * 1 = USTAR;  2 = USTRESS; OR 3 = U10 (DEFAULT).         *'
   WRITE (IU06,*) ' * BUT CODE IS: ', CODE_IN
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF

END SUBROUTINE SET_WIND_HEADER_M

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_WIND_TIMESTEPS (IN, OUT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_WIND_TIMESTEPS - SET WIND TIMESTEPS IN WAM_WIND MODULE.                !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE WIND TIMESTEPS TO WAM_WIND_MODULE.                        !
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

INTEGER, INTENT(IN)           :: IN    !! WIND INPUT TIME STEP.
INTEGER, INTENT(IN), OPTIONAL :: OUT   !! WIND OUTPUT TIME STEP.

! ---------------------------------------------------------------------------- !

IDELWI = IN
IF (PRESENT(OUT)) THEN
   IDELWO = OUT
   IF (OUT.LE.0) IDELWO = IN
ELSE
   IDELWO = IN
END IF

!         CHECK WIND INPUT AND WIND OUTPUT TIMESTEP.                           !
!         ------------------------------------------                           !

IF (IN.LE.0 .OR. IDELWO.LE.0) THEN
   WRITE(IU06,*) '***********************************************************'
   WRITE(IU06,*) '*                                                         *'
   WRITE(IU06,*) '*       FATAL ERROR IN SUB. SET_WIND_TIMESTEPS            *'
   WRITE(IU06,*) '*       ========================================          *'
   WRITE(IU06,*) '* WIND INPUT OR OUTPUT TIME STEP IS NOT POSITIVE.         *'
   WRITE(IU06,*) '* WIND INPUT TIMESTEP    IN = ', IDELWI, ' SECONDS'
   WRITE(IU06,*) '* WIND OUTPUT TIMESTEP  OUT = ', IDELWO, ' SECONDS'
   WRITE(IU06,*) '*                                                         *'
   WRITE(IU06,*) '*        PROGRAM ABORTS         PROGRAM ABORTS            *'
   WRITE(IU06,*) '*                                                         *'
   WRITE(IU06,*) '***********************************************************'
   CALL ABORT1
END IF
IF (IDELWI.LT.IDELWO) THEN
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+       WARNING ERROR IN SUB. SET_WIND_TIMESTEPS          +'
   WRITE(IU06,*) '+       ========================================          +'
   WRITE(IU06,*) '+ WIND INPUT TIME STEP IS LESS THAN WIND OUTPUT STEP      +'
   WRITE(IU06,*) '+ WIND INPUT TIMESTEP    IN = ', IDELWI, ' SECONDS'
   WRITE(IU06,*) '+ WIND OUTPUT TIMESTEP  OUT = ', IDELWO, ' SECONDS'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+ WIND INPUT AND OUTPUT CHANGED TO THE MULIPLE OF THE     +'
   WRITE(IU06,*) '+ INPUT STEP NEAREST TO WIND OLD OUTPUT TIME STEP         +'
   WRITE(IU06,*) '+ THE MODEL MAY IGNORE WIND FIELDS                        +'
   IDELWI = IDELWI*MAX(NINT(REAL(IDELWO)/REAL(IDELWI)),1)
   IDELWO = IDELWI
   WRITE(IU06,*) '+ NEW INPUT  TIMESTEP   IS: ',IDELWI, ' SECONDS'
   WRITE(IU06,*) '+ NEW OUTPUT TIMESTEP   IS: ',IDELWO, ' SECONDS'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF
IF (MOD(IDELWI,IDELWO).NE.0) THEN
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+       WARNING ERROR IN SUB. SET_WIND_TIMESTEPS          +'
   WRITE(IU06,*) '+       ========================================          +'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+ WIND INPUT AND WIND OUTPUT TIME STEPS DO NOT HAVE       +'
   WRITE(IU06,*) '+ INTEGER RATIO.                                          +'
   WRITE(IU06,*) '+ WIND INPUT TIMESTEP    IN = ', IDELWI, ' SECONDS'
   WRITE(IU06,*) '+ WIND OUTPUT TIMESTEP  OUT = ', IDELWO, ' SECONDS'
   WRITE(IU06,*) '+                                                         +'
   IDELWO = IDELWI/MAX(NINT(REAL(IDELWI)/REAL(IDELWO)),1)
   WRITE(IU06,*) '+ OUTPUT TIMESTEP CHANGED TO: ',IDELWO, ' SECONDS'
   WRITE(IU06,*) '+                                                         +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

END SUBROUTINE SET_WIND_TIMESTEPS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WAM_WIND (US, DS, CD_START)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   WAM_WIND - ROUTINE TO READ AND PROCESS ONE WINDFIELD.                      !
!                                                                              !
!     H. GUNTHER      ECMWF   MAY 1990     MODIFIED FOR SUB VERSION.           !
!     H. GUNTHER      ECMWF   DECEMBER 90  MODIFIED FOR CYCLE_4.               !
!     H. GUNTHER      GKSS    SEPTEMBER 2000   FT90.                           !
!     H. GUNTHER      GKSS    DECEMBER 2009   RE-ORGANIZED.                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       READ, INTERPOLATE AND CONVERT INPUT WINDS TO WAM WINDS.                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       READ A WINDFIELD FROM THE WINDFILE (SEARCH FOR IT) AND                 !
!       INTERPOLATED IT TO THE WAVE MODEL SEA POINTS.                          !
!       THE INTERPOLATED VALUES ARE TRANSFORMED TO MAGNITUDE AND DIRECTION.    !
!       INPUT MAY BE WIND IN 10M HEIGHT, SURFACE WINDS OR FRICTION VELOCETIES. !
!       THE INPUT GRID MUST BE ON A LATITUDE/LONGITUDE GRID EITHER PERIODIC    !
!       OR NON PERIODIC.                                                       !
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

REAL,          INTENT(OUT)    :: US(:)    !! WIND SPEED (U10).
REAL,          INTENT(OUT)    :: DS(:)    !! DIRECTION.
CHARACTER (LEN=14),INTENT(IN) :: CD_START !! DATE OF FIELD TO BE LOOKED FOR.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: ALPHACH = 0.0185

INTEGER :: IJ
REAL    :: UU, VV, USTAR, Z0, CD

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. READ WIND DATA AND CHECK DATE.                                        !
!        ------------------------------                                        !

DO
   CALL READ_WIND_INPUT

   IF (CD_READ.EQ.CD_START) EXIT
   IF (CD_READ.GT.CD_START) THEN
         WRITE (IU06,*) ' ******************************************'
         WRITE (IU06,*) ' *                                        *'
         WRITE (IU06,*) ' *      FATAL ERROR SUB. WAM_WIND         *'
         WRITE (IU06,*) ' *      =========================         *'
         WRITE (IU06,*) ' * WIND DATE READ IS LATER THAN EXPECTED  *'
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
!     2. INTERPOLATE AND BLOCK WINDFIELD                                       !
!        -------------------------------                                       !


US = 0.  !! INITIALISE WIND ARRAY WITH ZERO.
DS = 0.
IF (EQUAL_GRID) THEN
   US = PACK(U_IN, L_S_MASK)
   DS = PACK(V_IN, L_S_MASK)
ELSE
   CALL INTERPOLATION_TO_GRID (IU06, PER, DX_IN, DY_IN, WEST_IN, SOUTH_IN,     &
&                              U_IN, US, V_IN, DS)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. TRANSFORM TO MAGNITUDE AND DIRECTION.                                 !
!         -------------------------------------                                !

DO IJ = 1,SIZE(US)
   UU = US(IJ)
   VV = DS(IJ)
   US(IJ) = SQRT(UU**2 + VV**2)
   IF (US(IJ).NE.0.) DS(IJ) = ATAN2(UU,VV)
   IF (DS(IJ).LT.0.) DS(IJ) = DS(IJ) + ZPI
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. PROCESS WINDS ACCORDING TO TYPE                                       !
!        NOTHING TO DO FOR WIND SPEED U10 (CODE_IN = 3).                       !
!        ---------------------------------------------                         !

IF (CODE_IN.EQ.1) THEN

!     3.2  INPUT IS FRICTION VELOCITY.                                         !
!          ---------------------------                                         !

   DO IJ = 1,SIZE(US)
         USTAR = MAX(0.01,US(IJ))
         Z0  = ALPHACH/G*USTAR**2
         CD  = XKAPPA/ALOG(10./Z0)
         US(IJ) = USTAR/CD
   END DO

ELSE IF (CODE_IN.EQ.2) THEN

!     3.3 INPUT WINDS ARE SURFACE STRESSES.                                    !
!         ---------------------------------                                    !
!                                                                              !
   DO IJ = 1,SIZE(US)
         USTAR = MAX (0.01, SQRT(US(IJ)/ROAIR))
         Z0  = ALPHACH/G*USTAR**2
         CD  = XKAPPA/ALOG(10./Z0)
         US(IJ) = USTAR/CD
   END DO
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. TEST OUTPUT OF WAVE MODEL BLOCKS                                      !
!        ---------------------------------                                     !

IF (ITEST.GE.3) THEN
   IJ = MIN (10,SIZE(US))
   WRITE (IU06,*) ' '
   WRITE (IU06,*) '      SUB. WAM_WIND: WINDFIELDS CONVERTED TO MODEL GRID'
   WRITE (IU06,*) ' '
   WRITE (IU06,*) ' US(1:10) = ', US(1:IJ)
   WRITE (IU06,*) ' DS(1:10) = ', DS(1:IJ)
END IF

END SUBROUTINE WAM_WIND

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
!       NOTIM: NO TIME INTERPOLATION: PROCESS WINDFIELDS.                      !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NO TIME INTERPOLATION:                                                 !
!       WINDFIELDS ARE PROCESSED EVERY IDELWI SECONDS (U,V),                   !
!       THE WINDS INTERPOLATED IN SPACE ONLY.                                  !
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

CHARACTER (LEN=14), INTENT(IN)  :: CD_START  !! DATE OF FIRST WIND FIELD.
CHARACTER (LEN=14), INTENT(IN)  :: CD_END    !! DATE OF LAST WIND FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER            :: MP
CHARACTER (LEN=14) :: CDTWIH
REAL               :: US(1:NSEA)  !! OUTPUT WIND FIELD ARRAY (U10).
REAL               :: DS(1:NSEA)  !! OUTPUT WIND FIELD ARRAY (DIRECTION).

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. LOOP OVER OUTPUT WIND TIMES.                                          !
!        ----------------------------                                          !

CDTWIH = CD_START
MP = 0
DO WHILE (CDTWIH.LE.CD_END)

!     1.1 READ ONE WIND FIELD AND TRANSFORM TO GRID.                           !
!         ------------------------------------------                           !

   CALL WAM_WIND (US, DS, CDTWIH)
   MP = MP + 1

!     1.2 SAVE IN MODULE WAM_WIND.                                             !
!         ------------------------                                             !

   U_STORE(:,MP)  = US(nijs:nijl)
   D_STORE(:,MP)  = DS(nijs:nijl)
   CD_STORE(MP) = CDTWIH

   IF (ITEST.GE.3) THEN
      WRITE(IU06,*) '       SUB. NOTIM: NEW WIND FILES AT CDTWIH = ', CDTWIH
   END IF

!     1.3 UPDATE WIND FIELD REQUEST TIME.                                      !
!         -------------------------------                                      !

   CALL INCDATE (CDTWIH,IDELWO)
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
!       TIME INTERPOLATION: PROCESS WINDFIELDS.                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       WINDFIELDS ARE READ IN EVERY IDELWI SECONDS,                           !
!       BI-LINEAR INTERPOLATED IN SPACE, AND TRANSFORMED TO WAM GRID.          !
!       MAGNITUDE AND DIRECTION ARE INTERPOLATED LINEARLY IN TIME.             !
!       WIND FIELDS AND SAVED IN WAM_WIND_MODULE.                              !
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

CHARACTER (LEN=14), INTENT(IN)  :: CD_START  !! DATE OF FIRST WIND FIELD.
CHARACTER (LEN=14), INTENT(IN)  :: CD_END    !! DATE OF LAST WIND FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER            :: MP, NTS, N
REAL               :: DEL
CHARACTER (LEN=14) :: CDT1, CDT2, CDTH

REAL               :: US1(nijs:nijl)  !! OUTPUT WIND FIELD ARRAY (U10).
REAL               :: DS1(nijs:nijl)  !! OUTPUT WIND FIELD ARRAY (DIRECTION).
REAL               :: US2(1:NSEA)     !! OUTPUT WIND FIELD ARRAY (U10).
REAL               :: DS2(1:NSEA)     !! OUTPUT WIND FIELD ARRAY (DIRECTION).

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALIZE TIME COUNTER AND FIRST WIND FIELD USED FOR INTERPOLATION.  !
!        --------------------------------------------------------------------  !

CDT1 = CDA
US1(:) = U10(:)
DS1(:) = UDIR(:)
CDTH = CDT1
CALL INCDATE (CDTH,IDELWO)
IF (CD_START.NE.CDTH) THEN
   WRITE(IU06,*) ' *******************************************'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *        FATAL ERROR IN --TIMIN--         *'
   WRITE(IU06,*) ' *        =========================        *'
   WRITE(IU06,*) ' * DATES DO NOT MATCH.                     *'
   WRITE(IU06,*) ' * START DATE FOR WIND IS     CD_START = ', CD_START
   WRITE(IU06,*) ' * LAST DATE SAVED IN MODULE IS   CDT1 = ', CDT1
   WRITE(IU06,*) ' * PROCESSING WILL BE ABORTED              *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *******************************************'
   CALL ABORT1
ENDIF

CDT2 = CDT1
CALL INCDATE(CDT2,IDELWI)
NTS = IDELWI/IDELWO
DEL = REAL(IDELWO)/REAL(IDELWI)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LOOP OVER INPUT WINDFIELDS.                                           !
!        ---------------------------                                           !

MP = 0
DO

!     2.1 READ ONE WINDFIELD AND TRANSFORM TO BLOCKS.                          !
!         -------------------------------------------                          !

   CDT2 = CDT1
   CALL INCDATE(CDT2,IDELWI)
   CALL WAM_WIND (US2, DS2, CDT2)

!     2.2 INTERPOLATE AND SAVE BLOCKED WIND FIELDS.                            !
!         -----------------------------------------                            !

   CDTH = CDT1
   DO N = 1,NTS
      MP = MP + 1
      CALL INCDATE(CDTH,IDELWO)
      CD_STORE(MP) = CDTH
      U_STORE(:,MP) = US1(:) + REAL(N)*DEL*(US2(nijs:nijl)-US1(:))
      D_STORE(:,MP) = DS2(nijs:nijl) - DS1(:)

      WHERE (ABS(D_STORE(:,MP)).GT.PI)                                         &
&           D_STORE(:,MP) = D_STORE(:,MP)-ZPI*SIGN(1.,D_STORE(:,MP))
      D_STORE(:,MP) = DS1(:) + REAL(N)*DEL*D_STORE(:,MP)
      D_STORE(:,MP) = MOD(D_STORE(:,MP)+ZPI,ZPI)
   END DO

   IF (ITEST.GE.3) THEN
      WRITE(IU06,*) '       SUB. TIMIN: WIND FIELDS FOR CD_STORE = ',          &
&                                           CD_STORE(MP-NTS+1:MP)
   END IF

!     2.3 UPDATE WIND FIELD REQUEST TIME AND READ NEXT IF REQUESTED.           !
!         ----------------------------------------------------------           !

   US1(:) = US2(nijs:nijl)
   DS1(:) = DS2(nijs:nijl)
   CDT1 = CDT2
   CALL INCDATE (CDTH,IDELWI)
   IF (CDTH.GT.CD_END) EXIT
END DO

END SUBROUTINE TIMIN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_WIND_MODULE

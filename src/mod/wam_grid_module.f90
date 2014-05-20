MODULE WAM_GRID_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL VARIABLES AND CONSTANTS TO DEFINE THE MODEL GRID. !
!   ALL PROCEDURES ARE INCLUDED TO COMPUTE THE INFORMATION.                    !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE          !! COORDINATE TYPE AND PROCEDURES

USE WAM_GENERAL_MODULE,  ONLY:  &
&       ABORT1,                 &  !! TERMINATES PROCESSING.
&       PRINT_ARRAY                !! PRINTS AN ARRAY.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

! ---------------------------------------------------------------------------- !
!
!     5. GENERAL MODEL GRID INFORMATION.
!        -------------------------------                                       !

CHARACTER (LEN=80)   :: HEADER=' ' !! HEADER OF MODEL RUN.
INTEGER              :: NX = -1    !! NUMBER OF LONGITUDES IN GRID.
INTEGER              :: NY = -1    !! NUMBER OF LATITUDES  IN GRID.
INTEGER              :: NSEA = -1  !! NUMBER OF SEA POINTS. 
LOGICAL              :: IPER       !! .TRUE. IF GRID IS PERIODIC.
LOGICAL              :: ONE_POINT  !! .TRUE. IF GRID HAS ONE POINT ONLY.
LOGICAL              :: REDUCED_GRID  !! .TRUE. IF GRID IS REDUCED.
INTEGER              :: AMOWEP     !! MOST WESTERN LONGITUDE IN GRID [M_SEC].
INTEGER              :: AMOSOP     !! MOST SOUTHERN LATITUDE IN GRID [M_SEC].
INTEGER              :: AMOEAP     !! MOST EASTERN LONGITUDE IN GRID [M_SEC].
INTEGER              :: AMONOP     !! MOST NORTHERN LATITUDE IN GRID [M_SEC].
INTEGER              :: XDELLA     !! GRID INCREMENT FOR LATITUDE [M_SEC].
INTEGER              :: XDELLO     !! GRID INCREMENT FOR LONG. AT EQUATOR [M_SEC].
REAL                 :: DELPHI     !! GRID INCREMENT FOR LATITUDE [M].
INTEGER, ALLOCATABLE :: NLON_RG(:) !! NUMBER OF GRID POINT ON LATITUDES.
INTEGER, ALLOCATABLE :: ZDELLO(:)  !! GRID INCREMENTS FOR LONG. [M_SEC].
REAL,    ALLOCATABLE :: DELLAM(:)  !! GRID INCREMENTS FOR LONG. [M].
REAL,    ALLOCATABLE :: SINPH(:)   !! SIN OF LATITUDE.
REAL,    ALLOCATABLE :: COSPH(:)   !! COS OF LATITUDE. 

LOGICAL, ALLOCATABLE :: L_S_MASK(:,:) !! .TRUE. AT SEA POINTS.

REAL,    ALLOCATABLE :: DEPTH_B(:) !! WATER DEPTH [M].
INTEGER, ALLOCATABLE :: IXLG(:)    !! LONGITUDE GRID INDEX.
INTEGER, ALLOCATABLE :: KXLT(:)    !! LATITUDE GRID INDEX.

INTEGER, ALLOCATABLE :: KLAT(:,:)  !! INDEX OF GRIDPOINT SOUTH AND NORTH
                                   !! LANDPOINTS ARE MARKED BY ZERO.
INTEGER, ALLOCATABLE :: KLON(:,:)  !! INDEX OF GRIDPOINT WEST AND EAST
                                   !! LANDPOINTS ARE MARKED BY ZERO.

PUBLIC :: HEADER, NX, NY, NSEA, NLON_RG, IPER, AMOWEP, AMOSOP, AMOEAP, AMONOP,&
&         XDELLA, XDELLO, DELPHI, ZDELLO, DELLAM, SINPH, COSPH, DEPTH_B,      &
&         IXLG, KXLT, KLAT, KLON, L_S_MASK, REDUCED_GRID, ONE_POINT 

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE PRINT_GRID_STATUS          !! PRINTS GRID MODULE STATUS.
   MODULE PROCEDURE PRINT_GRID_STATUS
END INTERFACE
PUBLIC PRINT_GRID_STATUS

INTERFACE FIND_SEA_POINT              !! FIND SEA POINT NUMBER. 
   MODULE  PROCEDURE FIND_SEA_POINT
END INTERFACE
PUBLIC FIND_SEA_POINT

INTERFACE INTERPOLATION_TO_GRID       !! INTERPOLATES TO MODEL GRID POINTS. 
   MODULE  PROCEDURE INTERPOLATION_TO_GRID
END INTERFACE
PUBLIC INTERPOLATION_TO_GRID

INTERFACE EQUAL_TO_M_GRID             !! COMPARES A GRID TO MODEL GRID. 
   MODULE  PROCEDURE EQUAL_TO_M_GRID
END INTERFACE
PUBLIC EQUAL_TO_M_GRID

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_GRID_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PRINT_GRID_STATUS - PRINT STATUS OF GRID MODULE.                           !
!                                                                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       MAKE A PRINTER OUTPUT OF THE DATA, WHICH ARE SAVED IN WAM_GRID_MODULE. !
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

INTEGER             :: I
INTEGER             :: ILAT
REAL                :: GRID(NX,NY)   !! ARRAY FOR GRIDDED PRINT OUTPUT.
CHARACTER (LEN=100) :: TITL
CHARACTER (LEN=14)  :: ZERO = ' '
character (len=len_coor) :: formtext, ftext1, ftext2

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. PRINT STATUS.                                                         !
!        -------------                                                         !

WRITE (IU06,'(/,'' ----------------------------------------'')')
WRITE (IU06,'(  ''            GRID MODULE STATUS'')')
WRITE (IU06,'(  '' ----------------------------------------'')')
WRITE (IU06,*) ' '
WRITE (IU06,*) ' GRID HEADER: ', HEADER
WRITE (IU06,*) ' '

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. MODEL GRID DEFINITIONS.                                               !
!        -----------------------                                               !

IF (NX.GT.0 .AND. NY.GT.0) THEN
   WRITE (IU06,'(/,'' BASIC MODEL GRID: '')')
   WRITE (IU06,'('' NUMBER OF LONGITUDES IS       NX = '',I5)') NX
   formtext = write_coor_text (amowep)
   WRITE (IU06,'('' MOST WESTERN LONGITUDE IS  AMOWEP = '',A)') formtext
   formtext = write_coor_text (amoeap)
   WRITE (IU06,'('' MOST EASTERN LONGITUDE IS  AMOEAP = '',A)') formtext
   formtext = write_coor_text (xdello)
   WRITE (IU06,'('' LONGITUDE INCREMENT IS     XDELLO = '',A)') formtext
   WRITE (IU06,'('' NUMBER OF LATITUDES IS        NY = '',I5)') NY
   formtext = write_coor_text (amosop)
   WRITE (IU06,'('' MOST SOUTHERN LATITUDE IS  AMOSOP = '',A)') formtext
   formtext = write_coor_text (amonop)
   WRITE (IU06,'('' MOST NORTHERN LATITUDE IS  AMONOP = '',A)') formtext
   formtext = write_coor_text (xdella)
   WRITE (IU06,'('' LATITUDE INCREMENT IS      XDELLA = '',A)') formtext
    
   IF (ONE_POINT) THEN
      WRITE (IU06,*) 'THIS A ONE POINT GRID: PROPAGATION IS NOT DONE'
   ELSE
      IF (IPER) THEN
         WRITE (IU06,*) 'THE GRID IS EAST-WEST PERIODIC'
      ELSE
         WRITE (IU06,*) 'THE GRID IS NOT EAST-WEST PERIODIC'
      END IF
   END IF
   WRITE(IU06,*) ' '
   IF (REDUCED_GRID) THEN
      WRITE (IU06,*) 'A REDUCED GRID IS SET-UP'
      WRITE(IU06,*) ' '
      
      WRITE (IU06,*) '  NO. |    LATITUDE   |   NLON |   DELTA_LON   |'
      WRITE (IU06,*) '------|---------------|--------|---------------|'
      DO I = NY, 1, -1
        ILAT = AMOSOP +REAL(I-1)*XDELLA
        ftext1 = write_coor_text (ilat)
        ftext2 = write_coor_text (zdello(i))
        WRITE (IU06,'(I6,'' | '',A,'' | '',I6,'' | '',A,'' | '')')             &
&             I, ftext1, NLON_RG(I), ftext2
      END DO
   ELSE
      WRITE (IU06,*) 'A REGULAR GRID IS SET-UP'
   END IF
ELSE
   WRITE (IU06,*) ' MODEL GRID DEFINITIONS ARE NOT DEFINED'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. OUTPUT OF DEPTH FIELD.                                                !
!        ----------------------                                                !
  
IF (NSEA.GT.0) THEN
   WRITE(IU06,*) ' '
   WRITE (IU06,'('' NUMBER OF SEAPOINTS IS       NSEA = '',I7)') NSEA
   GRID = UNPACK(MIN(DEPTH_B,999.),L_S_MASK,99999.)
   DO I=1,NY
       GRID(NLON_RG(I)+1:NX,I) = -99
   END DO
   WRITE (IU06,*) ' '
   TITL = 'BASIC WATER DEPTH (M). (DEPTH DEEPER THAN 999M ARE PRINTED AS 999)'
   CALL PRINT_ARRAY (IU06, ZERO, TITL, GRID, AMOWEP, AMOSOP, AMOEAP, AMONOP)
ELSE
   WRITE (IU06,*) '  MODULE DATA ARE NOT PREPARED'
END IF

END SUBROUTINE PRINT_GRID_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FIND_SEA_POINT (NUMBER, LATITUDE, LONGITUDE, POINT_NO)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   FIND_SEA_POINT - FIND SEA POINT NUMBER.                                    !
!                                                                              !
!     R. PORTZ     MPI         15/01/1991                                      !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       FIND SEA POINT NUMBERS FOR A GIVEN ARRAY OF LONGITUDES AND LATITUDES.  !
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
! 
!     INTERFACE VARIABLES.
!     --------------------

INTEGER, INTENT(IN)    :: NUMBER       !! NUMBER OF POINTS IN ARRAYS.
INTEGER, INTENT(IN)    :: LATITUDE(*)  !! INPUT LATITUDES.
INTEGER, INTENT(IN)    :: LONGITUDE(*) !! INPUT LONGITUDES.
INTEGER, INTENT(OUT)   :: POINT_NO(*)  !! OUTPUT SEA POINT NUMBERS.

! ---------------------------------------------------------------------------- !
!
!     LOCAL VARIABLES.
!     ----------------

INTEGER  :: IO, IOLT, IOLG

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. LOOP OVER INPUT LATITUDES, LONGITUDES.                                !
!        --------------------------------------                                !

POINT_NO(1:NUMBER) = 0
POINT: DO IO = 1,NUMBER

!     1.1 COMPUTE GRID MATRIX INDICES.                                         !
!         ----------------------------                                         !

   IOLT = NINT(REAL(LATITUDE(IO)-AMOSOP)/REAL(XDELLA)+1.0)
   IF (IOLT.LT.1.OR.IOLT.GT.NY) CYCLE POINT
   IOLG = NINT(REAL(MOD(LONGITUDE(IO)-AMOWEP+2*M_S_PER,M_S_PER))               &
&             / REAL(ZDELLO(IOLT))+1.0)
   IF (IOLG.EQ.NLON_RG(IOLT)+1 .AND. IPER) IOLG = 1

!     1.2 IF SEA POINT FIND SEA POINT NUMBER.                                  !
!         -----------------------------------                                  !

   IF (IOLG.LT.1.OR.IOLG.GT.NLON_RG(IOLT)) CYCLE POINT
   IF (.NOT. L_S_MASK(IOLG,IOLT)) CYCLE POINT
   POINT_NO(IO) = COUNT(L_S_MASK(1:IOLG,IOLT))
   IF (IOLT.GT.1) POINT_NO(IO) = POINT_NO(IO) + COUNT(L_S_MASK(1:NX,1:IOLT-1))

END DO POINT

END SUBROUTINE FIND_SEA_POINT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE INTERPOLATION_TO_GRID (IU06, IGPER, DLAM, DPHI, RLONL, RLATS,       &
&                                 U, US, V, VS)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   INTERPOLATION_TO_GRID - INTERPOLATES TO MODEL GRID POINTS.                 !
!                                                                              !
!     H. GUNTHER    GKSS  DECEMBER 2001.                                       !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        LOCATE AND INTERPOLATE IN INPUT GRID.                                 !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       DOUBLE LINEAR INTERPOLATION IN INPUT GRID. OPTIONAL A SECOND INPUT     !
!       CAN BE INTERPOLATED AT THE SAME CALL.                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.
!     --------------------

INTEGER, INTENT(IN)           :: IU06   !! PRINTER OUTPUT UNIT.
LOGICAL, INTENT(IN)           :: IGPER  !! .TRUE. = PERIODIC INPUT GRID 
                                        !! OTHERWISE NON-PERIODICAL
INTEGER, INTENT(IN)           :: DLAM   !! INPUT LONGITUDE STEP    (M_SEC).
INTEGER, INTENT(IN)           :: DPHI   !! INPUT LATITUDE  STEP    (M_SEC).
INTEGER, INTENT(IN)           :: RLATS  !! INPUT SOUTHERN LATITUDE (M_SEC).
INTEGER, INTENT(IN)           :: RLONL  !! INPUT WESTERN LONGITUDE (M_SEC).

REAL,    INTENT(IN)           :: U(:,:) !! INPUT FIELD.
REAL,    INTENT(OUT)          :: US(:)  !! SPACE INTERPOLATED OUTPUT FIELD.

REAL,    INTENT(IN), OPTIONAL :: V(:,:) !! OPTIONAL SECOND INPUT FIELD.
REAL,    INTENT(OUT),OPTIONAL :: VS(:)  !! OPTIONAL SECOND OUTPUT FIELD.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.
!     ----------------

INTEGER :: NC, NR, NSEA
INTEGER :: IJ, I1(1:SIZE(US)), I2(1:SIZE(US)), K1(1:SIZE(US)), K2(1:SIZE(US))
REAL    :: DI(1:SIZE(US)), DK(1:SIZE(US))

NC = SIZE(U,1)
NR = SIZE(U,2)
NSEA = SIZE(US)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. TRANSFORM MODEL COORDINATES TO INPUT GRID.                            !
!        ------------------------------------------                            !

I1(:) = AMOWEP + (IXLG(:)-1)*ZDELLO(KXLT(:)) - RLONL
I1(:) = MOD(I1(:)+2*M_S_PER,M_S_PER)
DI(:) = REAL(I1(:))/REAL(DLAM)+1.00001

K1(:) = AMOSOP + (KXLT(:)-1)*XDELLA - RLATS
DK(:) = REAL(K1(:))/REAL(DPHI)+1.00001

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE CORNER POINT INDICES IN INPUT GRID.                           !
!        -------------------------------------------                           !

I1  = INT(DI)
K1  = INT(DK)
K2  = MIN(NR,K1+1)
I2  = I1+1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. DISTANCES OF INTERPOLATION POINT FROM LOW LEFT CORNER POINT.          !
!        ------------------------------------------------------------          !

DI = DI-REAL(I1)
DK = DK-REAL(K1)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. CORRECTIONOF FIRST AND LAST GRID LINES (PERIODIC OR UNPERIODIC GRID). !
!        --------------------------------------------------------------------- !

IF (IGPER) THEN
   WHERE (I1.EQ.NC) I2 = 1
   WHERE (I1.EQ.0 ) I1 = NC
ELSE
   WHERE (I1.EQ.NC) I2 = NC
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. CHECK WHETHER POINTS ARE IN GRID.                                     !
!        ---------------------------------                                     !

IF (MINVAL(I1).LT.1 .OR. MAXVAL(I1).GT.NC .OR.                                 &
&   MINVAL(K1).LT.1 .OR. MAXVAL(K1).GT.NR) THEN
    WRITE(IU06,*) ' *******************************************'
    WRITE(IU06,*) ' *                                         *'
    WRITE(IU06,*) ' *  FATAL ERROR IN INTERPOLATION_TO_GRID   *'
    WRITE(IU06,*) ' *  ====================================   *'
    WRITE(IU06,*) ' * POINT IS OUTSIDE OF INPUT GRID          *'
    WRITE(IU06,*) ' * DIMENSION OF INPUT GRID IS   NC = ', NC
    WRITE(IU06,*) ' *                              NR = ', NR
    WRITE(IU06,*) ' * MIN AND MAX OF INDEX ARE                *'
    WRITE(IU06,*) ' * I1:  MIN, MAX = ', MINVAL(I1), MAXVAL(I1)
    WRITE(IU06,*) ' * I2:  MIN, MAX = ', MINVAL(I2), MAXVAL(I2)
    WRITE(IU06,*) ' * K1:  MIN, MAX = ', MINVAL(K1), MAXVAL(K1)
    WRITE(IU06,*) ' * K2:  MIN, MAX = ', MINVAL(K2), MAXVAL(K2)
    WRITE(IU06,*) ' *                                         *'
    WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS      *'
    WRITE(IU06,*) ' *                                         *'
    WRITE(IU06,*) ' *******************************************'
    CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. LINEAR INTERPOLATION.                                                 !
!        ----------------------                                                !

!     6.1 FIRST FIELD.

US = 0.
DO IJ = 1,NSEA
   US(IJ) = (U(I1(IJ),K1(IJ))*(1.-DI(IJ))+U(I2(IJ),K1(IJ))*DI(IJ))*(1.-DK(IJ)) &
&         + (U(I1(IJ),K2(IJ))*(1.-DI(IJ))+U(I2(IJ),K2(IJ))*DI(IJ))*DK(IJ)
END DO

!     6.2 OPTIONAL SECOND FIELD.

IF (PRESENT(V) .AND. PRESENT(VS)) THEN
   VS = 0.
   DO IJ = 1,NSEA
      VS(IJ) = (V(I1(IJ),K1(IJ))*(1.-DI(IJ))+V(I2(IJ),K1(IJ))*DI(IJ))         &
&                                                                *(1.-DK(IJ)) &
&            + (V(I1(IJ),K2(IJ))*(1.-DI(IJ))+V(I2(IJ),K2(IJ))*DI(IJ))*DK(IJ)
   END DO
END IF      

END SUBROUTINE INTERPOLATION_TO_GRID

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

LOGICAL FUNCTION EQUAL_TO_M_GRID (NX_IN, NY_IN, D_LON, D_LAT,                  &
&                                 WEST, SOUTH, EAST, NORTH)     RESULT(L)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   EQUAL_TO_M_GRID - COMPARES A GRID TO THE MODEL GRID.                       !
!                                                                              !
!     H. GUNTHER    GKSS  DECEMBER 2001.                                       !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        TO CHECK IF A GRID IS EQUAL TO THE MODEL GRID.                        !
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
!     INTERFACE VARIABLES.
!     --------------------

INTEGER, INTENT(IN)           :: NX_IN   !! INPUT NO. OF LONGITUDES.
INTEGER, INTENT(IN)           :: NY_IN   !! INPUT NO. OF LATITUDE.
INTEGER, INTENT(IN)           :: D_LON   !! INPUT LONGITUDE STEP    (M_SEC).
INTEGER, INTENT(IN)           :: D_LAT   !! INPUT LATITUDE  STEP    (M_SEC).
INTEGER, INTENT(IN)           :: WEST    !! INPUT WESTERN LONGITUDE (M_SEC).
INTEGER, INTENT(IN)           :: SOUTH   !! INPUT SOUTHERN LATITUDE (M_SEC).
INTEGER, INTENT(IN)           :: EAST    !! INPUT EASTERN LONGITUDE (M_SEC).
INTEGER, INTENT(IN)           :: NORTH   !! INPUT NORTHERN LATITUDE (M_SEC).

! ---------------------------------------------------------------------------- !

IF (REDUCED_GRID) THEN
   L = .FALSE.
   RETURN
END IF

L = NX_IN.EQ.NX     .AND. NY_IN.EQ.NY     .AND.                                 &
&   D_LON.EQ.XDELLO .AND. D_LAT.EQ.XDELLA .AND.                                 &
&   SOUTH.EQ.AMOSOP .AND. WEST.EQ.AMOWEP  .AND.                                 &
&   EAST.EQ.AMOEAP  .AND. NORTH.EQ.AMONOP

END FUNCTION EQUAL_TO_M_GRID

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_GRID_MODULE

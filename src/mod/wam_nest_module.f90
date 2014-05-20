MODULE WAM_NEST_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE STORES THE BOUNDARY INPUT VALUES FOR A FINE GRID RUN.          !
!   THE VALUES WERE PRODUCED BY A PREVIOUS COARSE GRID RUN.                    !
!                                                                              !
! ---------------------------------------------------------------------------- !

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  METHODS FROM BASIC MODULES.                                          !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

USE WAM_COORDINATE_MODULE           !! COORDINATE PROCEDURES

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1,                  &  !! TERMINATES PROCESSING.
&       PRINT_ARRAY                 !! PRINTS AN ARRAY.

USE WAM_GRID_MODULE,      ONLY:  &
&       FIND_SEA_POINT              !! FIND SEA POINT NUMBERS.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B.  DATA FROM BASIC MODULES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

USE WAM_FILE_MODULE,      ONLY: IU06, ITEST, IU10, FILE10

USE WAM_GRID_MODULE,      ONLY: NX, NY, NSEA, NLON_RG, XDELLA, XDELLO, ZDELLO, &
&                               AMOWEP, AMOSOP, AMOEAP, AMONOP, IPER, IXLG,    &
&                               KXLT, L_S_MASK

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C.  MODULE DATA.                                                         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

IMPLICIT NONE

CHARACTER (LEN=11) , PARAMETER :: MODULE_NAME = 'NEST_MODULE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. BOUNDARY OPTIONS.                                                     !
!        -----------------                                                     !

LOGICAL      :: COARSE = .FALSE.  !! = .TRUE. IF COARSE GRID RUN.
LOGICAL      :: FINE   = .FALSE.  !! = .TRUE. IF FINE GRID RUN.

! ---------------------------------------------------------------------------- !
!
!    2. COURSE GRID NEST ORGANIZATION.

INTEGER, SAVE                   :: N_NEST = 0 !! NUMBER OF NESTS
INTEGER, SAVE                   :: MAX_NEST=0 !! MAXIMUM NUMBER OF NEST POINTS.
INTEGER,            ALLOCATABLE :: NBOUNC(:)  !! NUMBERS OF NEST POINTS.
INTEGER,            ALLOCATABLE :: N_SOUTH(:) !! MOST SOUTHERN LAT. OF NESTS.
INTEGER,            ALLOCATABLE :: N_NORTH(:) !! MOST NORTHERN LAT. OF NESTS.
INTEGER,            ALLOCATABLE :: N_WEST(:)  !! MOST WESTERN LONG. OF NESTS.
INTEGER,            ALLOCATABLE :: N_EAST(:)  !! MOST EASTERN LONG. OF NESTS.
CHARACTER (LEN=20), ALLOCATABLE :: N_NAME(:)  !! NAMES OF NESTS.
integer,            allocatable :: n_code(:)  !! ascii or binary code


INTEGER,      ALLOCATABLE :: IJARC(:,:)  !! SEA POINT INDEX OF NEST POINTS.
INTEGER,      ALLOCATABLE :: BLATC(:,:)  !! LATITUDES  OF NEST POINTS.
INTEGER,      ALLOCATABLE :: BLNGC(:,:)  !! LONGITUDES  OF NEST POINTS.
INTEGER,      ALLOCATABLE :: N_ZDEL(:,:) !! LONGITUDE INCREMENTS. 
                                         !! FIRST INDEX IS POINT NUMBER 
                                         !! SECOND INDEX IS NEST NUMBER 
! ---------------------------------------------------------------------------- !
!
!    3. FINE GRID ORGANIZATION THE BOUNDARY POINTS.

INTEGER            :: NBOUNF = 0   !! NUMBER OF FINE GRID BOUNDARY POINTS.
INTEGER            :: NBINP  = 0   !! NUMBER OF INPUT BOUNDARY POINTS.
CHARACTER (LEN=20) :: C_NAME       !! NAME OF NEST GIVEN IN COARSE GRID.
INTEGER            :: DLAMAC = -1. !! LONGITUDE INCREMENT OF COARSE GRID (M_SEC).
INTEGER            :: DPHIAC = -1. !! LATITUDE INCREMENT OF COARSE GRID (M_SEC).

INTEGER, ALLOCATABLE :: IJ_C(:)
INTEGER, ALLOCATABLE :: LAT_C(:)
INTEGER, ALLOCATABLE :: LON_C(:)
INTEGER, ALLOCATABLE :: BLATF(:)
INTEGER, ALLOCATABLE :: BLNGF(:)
INTEGER, ALLOCATABLE :: ZDEL_C(:)

INTEGER,ALLOCATABLE, DIMENSION(:) :: IJARF !! POINT INDEX OF BOUNDARY POINT. 
INTEGER,ALLOCATABLE, DIMENSION(:) :: IBFL  !! INDEX OF LEFT 
                                           !! COARSE GRID OUTPUT POINT.
INTEGER,ALLOCATABLE, DIMENSION(:) :: IBFR  !! INDEX OF RIGHT
                                           !! COARSE GRID OUTPUT POINT.
REAL,   ALLOCATABLE, DIMENSION(:) :: BFW   !! SPACE INTERPOLATION WEIGHT.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTERFACE SET_NEST                    !! TRANSFERS NEST (COARSE) DEFINITION 
   MODULE PROCEDURE SET_NEST_C
   MODULE PROCEDURE SET_NEST_D
   MODULE PROCEDURE SET_NEST_M
END INTERFACE
PUBLIC SET_NEST

INTERFACE PREPARE_BOUNDARY_nest       !! CHECKS THE BOUNDARY OPTION.
   MODULE PROCEDURE PREPARE_BOUNDARY_nest 
END INTERFACE
PUBLIC PREPARE_BOUNDARY_nest

INTERFACE SET_BOUNDARY_OPTION         !! SETS THE BOUNDARY OPTION.
   MODULE PROCEDURE SET_BOUNDARY_OPTION
END INTERFACE
PUBLIC SET_BOUNDARY_OPTION

INTERFACE PRINT_NEST_STATUS           !! PRINTS THIS MODULE STATUS.
   MODULE PROCEDURE PRINT_NEST_STATUS 
END INTERFACE
PUBLIC PRINT_NEST_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTERFACE MAKE_NEST                  !! MAKES INFORMATION FOR ONE NEST.
   MODULE PROCEDURE MAKE_NEST
END INTERFACE
PUBLIC MAKE_NEST

INTERFACE MAKE_FINE_BOUNDARY         !! MAKE FINE GRID BOUNDARY.
   MODULE PROCEDURE MAKE_FINE_BOUNDARY
END INTERFACE
PUBLIC MAKE_FINE_BOUNDARY

INTERFACE MAKE_BOX                   !! MAKES A BOX IN THE GRID
   MODULE  PROCEDURE MAKE_BOX
END INTERFACE
PRIVATE :: MAKE_BOX

INTERFACE MINTF                      !! INTERPOLATION TABLES FOR BOUNDARY INPUT.
   MODULE PROCEDURE MINTF
END INTERFACE
PRIVATE MINTF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_NEST_C (SOUTH, NORTH, WEST, EAST, NAME, ncode)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=LEN_COOR), INTENT(IN)  :: SOUTH(:)   !! SOUTH LATITUDE OF NEST.
CHARACTER (LEN=LEN_COOR), INTENT(IN)  :: NORTH(:)   !! NORTH LATITUDE OF NEST.
CHARACTER (LEN=LEN_COOR), INTENT(IN)  :: WEST(:)    !! WEST LONGITUDE OF NEST.
CHARACTER (LEN=LEN_COOR), INTENT(IN)  :: EAST(:)    !! EAST LONGITUDE OF NEST.
CHARACTER (LEN=*),        INTENT(IN)  :: NAME(:)    !! NAME OF NEST
integer          ,        intent(in)  :: ncode(:)   !! ascii or binary code


! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CONVERT TO M_SEC AND CALL SET_NEST_M.                                 !
!        -------------------------------------                                 !

CALL SET_NEST_M (SOUTH=READ_COOR_TEXT(SOUTH), NORTH=READ_COOR_TEXT(NORTH),     &
&                WEST=READ_COOR_TEXT(WEST),  EAST=READ_COOR_TEXT(EAST),        &
&                NAME=NAME, ncode=ncode)

END SUBROUTINE SET_NEST_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_NEST_D (SOUTH, NORTH, WEST, EAST, NAME, ncode)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL (KIND=KIND_D), INTENT(IN)  :: SOUTH(:)   !! SOUTH LATITUDE OF NEST [DEG].
REAL (KIND=KIND_D), INTENT(IN)  :: NORTH(:)   !! NORTH LATITUDE OF NEST [DEG].
REAL (KIND=KIND_D), INTENT(IN)  :: WEST(:)    !! WEST LONGITUDE OF NEST [DEG].
REAL (KIND=KIND_D), INTENT(IN)  :: EAST(:)    !! EAST LONGITUDE OF NEST [DEG].
CHARACTER (LEN=*),  INTENT(IN)  :: NAME(:)    !! NAME OF NEST
integer          ,  intent(in)  :: ncode(:)   !! ascii or binary code


! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CONVERT TO M_SEC AND CALL SET_NEST_M.                                 !
!        -------------------------------------                                 !

CALL SET_NEST_M (SOUTH=DEG_TO_M_SEC(SOUTH), NORTH=DEG_TO_M_SEC(NORTH),         &
&                WEST=DEG_TO_M_SEC(WEST),   EAST=DEG_TO_M_SEC(EAST),           &
&                NAME=NAME, ncode=ncode)

END SUBROUTINE SET_NEST_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_NEST_M (SOUTH, NORTH, WEST, EAST, NAME, ncode)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,           INTENT(IN)  :: SOUTH(:)   !! SOUTH LATITUDE OF NEST [M_SEC].
INTEGER,           INTENT(IN)  :: NORTH(:)   !! NORTH LATITUDE OF NEST [M_SEC].
INTEGER,           INTENT(IN)  :: WEST(:)    !! WEST LONGITUDE OF NEST [M_SEC].
INTEGER,           INTENT(IN)  :: EAST(:)    !! EAST LONGITUDE OF NEST [M_SEC].
CHARACTER (LEN=*), INTENT(IN)  :: NAME(:)    !! NAME OF NEST
integer          , intent(in)  :: ncode(:)   !! ascii or binary code

INTEGER :: I, LEN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLEAR NEST INFORMATION IN MODULE.                                     !
!        ---------------------------------                                     !

N_NEST=0
IF (ALLOCATED(N_SOUTH)) DEALLOCATE(N_SOUTH)
IF (ALLOCATED(N_NORTH)) DEALLOCATE(N_NORTH)
IF (ALLOCATED(N_WEST)) DEALLOCATE(N_WEST)
IF (ALLOCATED(N_EAST)) DEALLOCATE(N_EAST)
IF (ALLOCATED(N_NAME)) DEALLOCATE(N_NAME)
if (allocated(n_code)) deallocate(n_code)

IF (ALLOCATED(NBOUNC)) DEALLOCATE(NBOUNC)
IF (ALLOCATED(IJARC)) DEALLOCATE(IJARC)
IF (ALLOCATED(BLATC)) DEALLOCATE(BLATC)
IF (ALLOCATED(BLNGC)) DEALLOCATE(BLNGC)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COPY BOUNDARY DEFINITIONS.                                            !
!        --------------------------                                            !

N_NEST = COUNT (SOUTH(:).NE.COOR_UNDEF)

IF (N_NEST.EQ.0) RETURN

ALLOCATE (N_SOUTH(1:N_NEST))
ALLOCATE (N_NORTH(1:N_NEST))
ALLOCATE (N_WEST(1:N_NEST))
ALLOCATE (N_EAST(1:N_NEST))
ALLOCATE (N_NAME(1:N_NEST))
allocate (n_code(1:n_nest))

N_SOUTH(1:N_NEST) = SOUTH(1:N_NEST)
N_NORTH(1:N_NEST) = NORTH(1:N_NEST)
N_WEST(1:N_NEST)  = WEST(1:N_NEST)
N_EAST(1:N_NEST)  = EAST(1:N_NEST)
DO I = 1,N_NEST
   LEN = MIN(LEN_TRIM(NAME(I)),20)
   N_NAME(I)  = NAME(I)(1:LEN)
END DO
n_code(1:n_nest) = ncode(1:n_nest)
CALL ADJUST (N_WEST, N_EAST)

END SUBROUTINE SET_NEST_M

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_BOUNDARY_OPTION (C , F)
LOGICAL, INTENT(IN), OPTIONAL :: C   !! COARSE GRID OPTION
LOGICAL, INTENT(IN), OPTIONAL :: F   !! FINE GRID OPTION

IF (PRESENT(C)) COARSE = C
IF (PRESENT(F)) FINE = F
 
END SUBROUTINE SET_BOUNDARY_OPTION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_NEST_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER             :: I, L
CHARACTER (LEN=1)   :: CGRID(NX,NY)  !! ARRAY FOR GRIDDED PRINT OUTPUT.
CHARACTER (LEN=100) :: TITL
character (len=len_coor) :: formtext1, formtext2

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. PRINT STATUS.                                                         !
!        -------------                                                         !

WRITE (IU06,*)' '
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '              NEST MODULE STATUS:'
WRITE(IU06,*) ' ------------------------------------------------- '

WRITE(IU06,*) '  '
WRITE(IU06,*) ' COARSE GRID (BOUNDARY VALUE OUTPUT): '
WRITE(IU06,*) ' ------------------------------------ '
WRITE(IU06,*) '  '

IF (COARSE) THEN
   WRITE(IU06,*) ' COARSE GRID RUN (OUTPUT OF BOUNDARY VALUES)'
   WRITE(IU06,*) ' NO. OF NESTS DEFINED IN GRID IS N_NEST = ', N_NEST
   WRITE(IU06,*) ' THE MAXIMUM NUMBER OF POINTS IN ALL NESTS IS: ',MAX_NEST

   DO L=1,N_NEST
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' INFORMATION OF NEST ',L,' NEST NAME IS : ',N_NAME(L)
      formtext1 = write_coor_text (n_north(l))
      WRITE(IU06,*) ' NORTH LATIDUTE = ', formtext1
      formtext1 = write_coor_text (n_south(l))
      WRITE(IU06,*) ' SOUTH LATIDUTE = ', formtext1
      formtext1 = write_coor_text (n_west(l))
      WRITE(IU06,*) ' WEST LONGITUDE = ', formtext1
      formtext1 = write_coor_text (n_east(l))
      WRITE(IU06,*) ' EAST LONGIDUTE = ', formtext1
      IF (ALLOCATED(NBOUNC)) THEN
         WRITE (IU06,*) ' NUMBER OF POINTS IS: ',NBOUNC(L)
         if (n_code(l)==1) then
            write (iu06,*) ' BOUNDARY VALUES WRITTEN IN +++ ASCII +++ CODE'
         else
            write (iu06,*) ' BOUNDARY VALUES WRITTEN IN +++ BINARY +++ CODE'
         endif
         IF (NBOUNC(L) .NE. 0 .AND. ITEST.GE.2) THEN
            WRITE(IU06,*) '  '
            WRITE(IU06,*) '      NO. |   LONGITUDE   |    LATITUDE   ',        &
&                         '|    POINT   |'
            WRITE(IU06,*) ' ---------|---------------|---------------',        &
&                         '|------------|'
            DO I = 1,NBOUNC(L)
               formtext1 = write_coor_text (blngc(i,l))
               formtext2 = write_coor_text (blatc(i,l))
               WRITE(IU06,'(I10,'' | '',A,'' | '',A,'' | '',I10,'' |'')')      &
&                             I, formtext1, formtext2, ijarc(i,l)
            END DO
         END IF
      ELSE
         WRITE(IU06,*) ' BOUNDARY POINTS ARE NOT INITIALIZED '
      END IF
   END DO
ELSE
   WRITE(IU06,*) ' MODEL RUNS WITHOUT BOUNDARY OUTPUT'
END IF

WRITE(IU06,*) '  '
WRITE(IU06,*) ' FINE GRID (BOUNDARY VALUE INPUT): '
WRITE(IU06,*) ' --------------------------------- '
WRITE(IU06,*) '  '
IF (FINE) THEN
   WRITE(IU06,*) ' FINE GRID RUN (INPUT OF BOUNDARY VALUES)'
   IF (NBOUNF .NE. 0) THEN
      WRITE (IU06,*) ' NUMBER OF SPECTRA FROM PREVIOUS COARSE GRID IS: ', NBINP
      WRITE (IU06,*) ' NUMBER OF BOUNDARY POINTS IS: ',NBOUNF
      WRITE (IU06,*) ' THE COARSE GRID NEST NAME IS: ',C_NAME
      WRITE (IU06,*) '  '
      IF (ITEST.GE.2) THEN
         WRITE (IU06,*) '         |------------FINE GRID INPUT POINT-',        &
&                       '---------|----RELATED COURSE GRID INDICES-----|'
         WRITE (IU06,*) '     NO. |   LONGITUDE   |    LATITUDE   |',          &
&                          ' POINT NO. |    LEFT   |    RIGHT  |    WEIGHT  |'
         WRITE (IU06,*) '---------|---------------|---------------|',          &
&                       '-----------|-----------|-----------|------------|'
         DO I = 1, NBOUNF
            formtext1 = write_coor_text (blngf(i))
            formtext2 = write_coor_text (blatf(i))
            WRITE (IU06,'(I9,'' | '',A,'' | '',A,'' | '',3(I9,'' | ''),        &
&                         F10.4,'' |'')')                                      &
&              I, formtext1, formtext2, ijarf(i), ibfl(i), ibfr(i), bfw(i)
         END DO
      END IF
   ELSE
      WRITE(IU06,*) ' BOUNDARY POINTS ARE NOT INITIALIZED '
   END IF
ELSE
   WRITE(IU06,*) ' MODEL RUNS WITHOUT BOUNDARY POINTS (FINE GRID)'
END IF
WRITE(IU06,*) '  '

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. OUTPUT OF LAND-SEA MASK.                                              !
!        ------------------------                                              !
  
IF (NSEA.GT.0 .AND. (COARSE .OR. FINE)) THEN
   WRITE (IU06,'('' NUMBER OF SEAPOINTS IS       NSEA = '',I7)') NSEA
   WRITE(IU06,*) ' '
   TITL = 'LAND-SEA MASK  (S = SEA, L = LAND)'
   CGRID ='L'
   WHERE (L_S_MASK) CGRID ='S'
   DO I=1,NY
       CGRID(NLON_RG(I)+1:NX,I) = ' '
   END DO
   IF (COARSE) THEN
      DO L=1,N_NEST
         DO I=1,NBOUNC(L)
           CGRID(IXLG(IJARC(I,L)), KXLT(IJARC(I,L))) = CHAR(64+L)
         END DO
      END DO
      L = LEN_TRIM(TITL)
      TITL = TITL(1:L-1)//', '//CHAR(65)//'-'// CHAR(64+N_NEST)//' = NESTS)'
   END IF
   IF (FINE) THEN
      DO I=1,NBOUNF
           CGRID(IXLG(IJARF(I)), KXLT(IJARF(I))) = '/'
      END DO
      L = LEN_TRIM(TITL)
      TITL = TITL(1:L-1)//', / = BOUNDARY INPUT)'
   END IF
   CALL PRINT_ARRAY (IU06, '              ', TITL, CGRID,                      &
&                    AMOWEP, AMOSOP, AMOEAP, AMONOP)
   WRITE(IU06,*) ' '
END IF

END SUBROUTINE PRINT_NEST_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_BOUNDARY_nest

INTEGER  :: NUMBER, I, N_NEST_NEW, I_WEG

IF (.NOT.COARSE .AND. .NOT.FINE) RETURN

IF (NSEA.LE.0) THEN
   WRITE (IU06,*) '********************************************************'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*        FATAL ERROR IN SUB. PREPARE_BOUNDARY          *'
   WRITE (IU06,*) '*        ====================================          *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*  THE MODEL GRID IS NOT DEFINED.                      *'
   WRITE (IU06,*) '*  NUMBER OF SEA POINTS IS NSEA = ', NSEA
   WRITE (IU06,*) '*  SUB. PREPARE_GRID HAS TO BE EXECUTED BEFORE         *'
   WRITE (IU06,*) '*  SUB. PREPARE_BOUNDARY                               *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*              THE PROGRAM ABORTS                      *'
   WRITE (IU06,*) '********************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COARSE GRIDS.                                                         !
!        -------------                                                         !

IF (COARSE) THEN
   DPHIAC = XDELLA
   DLAMAC = XDELLO

   NUMBER = MAXVAL(   (NINT(REAL(N_EAST - N_WEST)/REAL(XDELLO))                          &
&                    + NINT(REAL(N_NORTH - N_SOUTH)/REAL(XDELLA))) * 2)

   IF (ALLOCATED(NBOUNC)) DEALLOCATE(NBOUNC)
   IF (ALLOCATED(IJARC)) DEALLOCATE(IJARC)
   IF (ALLOCATED(BLATC)) DEALLOCATE(BLATC)
   IF (ALLOCATED(BLNGC)) DEALLOCATE(BLNGC)
   IF (ALLOCATED(N_ZDEL)) DEALLOCATE(N_ZDEL)

   ALLOCATE (NBOUNC(1:N_NEST))
   ALLOCATE (IJARC(1:NUMBER,1:N_NEST))
   ALLOCATE (BLATC(1:NUMBER,1:N_NEST))
   ALLOCATE (BLNGC(1:NUMBER,1:N_NEST))
   ALLOCATE (N_ZDEL(1:NUMBER,1:N_NEST))

   DO I = 1, N_NEST
      CALL MAKE_NEST (N_SOUTH(I), N_NORTH(I), N_WEST(I), N_EAST(I), N_NAME(I), &
&                    NBOUNC(I), IJARC(:,I), BLATC(:,I), BLNGC(:,I), N_ZDEL(:,I))
   END DO

   N_NEST_NEW = COUNT(NBOUNC.GT.0)
   MAX_NEST = MAXVAL(NBOUNC)
   IF (N_NEST_NEW .NE. N_NEST) THEN
      N_SOUTH(1: N_NEST_NEW) = PACK(N_SOUTH,NBOUNC.GT.0)
      N_NORTH(1: N_NEST_NEW) = PACK(N_NORTH,NBOUNC.GT.0)
      N_WEST (1: N_NEST_NEW) = PACK(N_WEST, NBOUNC.GT.0)
      N_EAST (1: N_NEST_NEW) = PACK(N_EAST, NBOUNC.GT.0)
      N_NAME (1: N_NEST_NEW) = PACK(N_NAME, NBOUNC.GT.0)
      I_WEG=0
      DO I = 1, N_NEST
         IF (NBOUNC(I).LE.0) THEN
            I_WEG = I_WEG + 1
            CYCLE
         END IF
         IF (I_WEG.EQ.0) CYCLE
         IF (I-I_WEG.LE.0) CYCLE
         IJARC(:,I-I_WEG) = IJARC(:,I)
         BLATC(:,I-I_WEG) = BLATC(:,I)
         BLNGC(:,I-I_WEG) = BLNGC(:,I)
         N_ZDEL(:,I-I_WEG) = N_ZDEL(:,I)
      END DO
      NBOUNC (1: N_NEST_NEW) = PACK(NBOUNC, NBOUNC.GT.0)
      N_NEST = N_NEST_NEW  
   END IF
   IF (N_NEST.EQ.0) COARSE =.FALSE.
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. FINE GRID.                                                            !
!        ----------                                                            !

IF (FINE) CALL MAKE_FINE_BOUNDARY

END SUBROUTINE PREPARE_BOUNDARY_nest

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MAKE_NEST (SOUTH, NORTH, WEST, EAST, NAME,                          &
&                     NUMBER, IJ_N, LAT_N, LON_N, ZDEL_N)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MAKE_NEST - MAKES COARSE GRID BOUNDARY.                                    !
!                                                                              !
!     R. PORTZ     MPI         15/01/1991                                      !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       COMPUTE ALL INFORMATION FOR COARSE GRID BOUNDARY VALUE OUTPUT.         !
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
!      INTERFACE VARIABLES.                                                    !

INTEGER           , INTENT(IN)    :: SOUTH    !! SOUTH LATITUDE OF NEST [M_SEC].
INTEGER           , INTENT(IN)    :: NORTH    !! NORTH LATITUDE OF NEST [M_SEC].
INTEGER           , INTENT(IN)    :: WEST     !! WEST LONGITUDE OF NEST [M_SEC].
INTEGER           , INTENT(IN)    :: EAST     !! EAST LONGITUDE OF NEST [M_SEC].
CHARACTER (LEN=20), INTENT(IN)    :: NAME     !! NAME OF NEST

INTEGER           , INTENT(INOUT) :: NUMBER   !! NUMBER OF NESTS POINTS.
INTEGER           , INTENT(INOUT) :: IJ_N(*)  !! SEA POINT NUMBERS.
INTEGER           , INTENT(INOUT) :: LAT_N(*) !! LATITUDES OF NEST [M_SEC].
INTEGER           , INTENT(INOUT) :: LON_N(*) !! LONGITUDES OF NEST [M_SEC].
INTEGER           , INTENT(INOUT) :: ZDEL_N(*) ! LONGITUDE STEP OF NEST [M_SEC].
 
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK

integer :: k, i, nbounew
character (len=len_coor) :: formtext

INTEGER, DIMENSION(4) :: C_NEST_I, C_NEST_K

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CHECK: IS FINE GRID IN COARSE GRID?                                   !
!        -----------------------------------                                   !

NUMBER = 1
IF (AMOSOP .GT. SOUTH .OR. NORTH .GT. AMONOP .OR. SOUTH .GE. NORTH) NUMBER = 0

IF ((.NOT.IPER) .AND.                                                          &
&   (AMOWEP.GT.WEST      .OR. WEST.GT.AMOEAP     )  .AND.                      &
&   (AMOWEP.GT.WEST+360. .OR. EAST+360..GT.AMOEAP)  .AND.                      &
&   (AMOWEP.GT.WEST-360. .OR. EAST-360..GT.AMOEAP) ) NUMBER = 0

IF (NUMBER .EQ. 0) THEN
   WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '+        WARNING ERROR IN SUB. MAKE__NEST            +'
   WRITE (IU06,*) '+        ================================            +'
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '+ ERROR IN NEST SPECIFICATIONS.                      +'
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '+ NEST NAME  IS: ', NAME
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '+ WEST, EAST, NORTH OR SOUTH BOUNDARY IS NOT         +'
   WRITE (IU06,*) '+ IN GRID AREA  OR                                   +'
   WRITE (IU06,*) '+ SOUTH IS GREATER OR EQUAL NORTH.                   +'
   WRITE (IU06,*) '+                                                    +'
   formtext = write_coor_text (amowep)
   WRITE (IU06,*) '+ GRID WEST  IS AMOWEP = ', formtext
   formtext = write_coor_text (amoeap)
   WRITE (IU06,*) '+ GRID EAST  IS AMOEAP = ', formtext
   formtext = write_coor_text (amonop)
   WRITE (IU06,*) '+ GRID NORTH IS AMONOP = ', formtext
   formtext = write_coor_text (amosop)
   WRITE (IU06,*) '+ GRID SOUTH IS AMOSOP = ', formtext
   formtext = write_coor_text (west)
   WRITE (IU06,*) '+ NEST WEST  IS   WEST = ', formtext
   formtext = write_coor_text (east)
   WRITE (IU06,*) '+ NEST EAST  IS   EAST = ', formtext
   formtext = write_coor_text (north)
   WRITE (IU06,*) '+ NEST NORTH IS  NORTH = ', formtext
   formtext = write_coor_text (south)
   WRITE (IU06,*) '+ NEST SOUTH IS  SOUTH = ', formtext
   WRITE (IU06,*) '+ NEST NAME  IS: ', NAME
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '+  NEST INFORMATION WILL NOT BE GENERATED            +'
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INITIAL.                                                              !
!        --------                                                              !

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTED THE SQUARE BOX.                                              !
!        ------------------------                                              !

NUMBER = (NINT(REAL(EAST -WEST)/REAL(XDELLO)+2)                                  &
&       + NINT(REAL(NORTH - SOUTH)/REAL(XDELLA))+2) * 2

CALL MAKE_BOX (WEST, SOUTH, EAST, NORTH, NUMBER, C_NEST_I, C_NEST_K,           &
&              LAT_N, LON_N)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. SEARCH SEA POINT NUMBER.                                              !
!        ------------------------                                              !

IJ_N(1:NUMBER)  = 0
CALL FIND_SEA_POINT (NUMBER, LAT_N, LON_N, IJ_N)

K = 0
DO I = C_NEST_I(1), C_NEST_I(2)
   K = K+1
   LAT_N(K) = SOUTH
END DO
LON_N(1) = WEST
LON_N(K) = EAST

DO I = C_NEST_K(1)+1, C_NEST_K(3)-1
   K = K+1
   LON_N(K) = WEST

   K = K+1
   LON_N(K) = EAST
END DO

LON_N(K+1) = WEST
DO I = C_NEST_I(3), C_NEST_I(4)
   K = K+1
   LAT_N(K) = NORTH
END DO
LON_N(K) = EAST

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. PACKED ALL ARRAYS.                                                    !
!        -------------------                                                   !

IF (.NOT.ALLOCATED(MASK)) ALLOCATE(MASK(1:NUMBER))

MASK = (IJ_N(1:NUMBER).GT.0)
NBOUNEW = COUNT(MASK)

LAT_N (1:NBOUNEW) = PACK (LAT_N(1:NUMBER), MASK(1:NUMBER))
LON_N (1:NBOUNEW) = PACK (LON_N(1:NUMBER), MASK(1:NUMBER))
IJ_N  (1:NBOUNEW) = PACK (IJ_N(1:NUMBER), MASK(1:NUMBER))
ZDEL_N(1:NBOUNEW) = ZDELLO(NINT(REAL(LAT_N(1:NBOUNEW) - AMOSOP)/REAL(XDELLA))+1)


NUMBER = NBOUNEW
DEALLOCATE (MASK)

END SUBROUTINE MAKE_NEST

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MAKE_FINE_BOUNDARY

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MAKE_FINE_BOUNDARY - MAKE FINE GRID BOUNDARY.                              !
!                                                                              !
!     R. PORTZ     MPI         15/01/1991                                      !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       COMPUTE ALL INFORMATION FOR FINE GRID BOUNDARY VALUE INPUT.            !
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

INTEGER              :: I, IO, NBOUNEW, LEN, N_NEST_C, MAX_NEST_C 
INTEGER              :: SOUTH_C, NORTH_C, EAST_C, WEST_C

LOGICAL, ALLOCATABLE :: MASK(:)
CHARACTER (LEN=80)   :: HEADER_COARSE

INTEGER :: F_NEST_I(4), F_NEST_K(4)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. READ INFO ABOUT COARSE GRID.                                          !
!     -------------------------------                                          !

LEN = LEN_TRIM(FILE10)
IO = 0
OPEN (UNIT=IU10, FILE=FILE10(1:LEN), FORM='UNFORMATTED', STATUS='OLD',IOSTAT=IO)
IF (IO.NE.0) THEN
   WRITE (IU06,*) '********************************************************'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*       FATAL ERROR IN SUB. MAKE_FINE_BOUNDARY         *'
   WRITE (IU06,*) '*       ======================================         *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*  A FINE GRID RUN IS REQUESTED BUT COARSE GRID        *'
   WRITE (IU06,*) '*  PREPROC OUTPUT COULD NOT BE OPENED.                 *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*  UNIT IS         IU10 = ', IU10
   WRITE (IU06,*) '*  FILE NAME IS  FILE10 = ', FILE10(1:LEN)
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*              THE PROGRAM ABORTS                      *'
   WRITE (IU06,*) '********************************************************'
   CALL ABORT1
END IF

READ (IU10) HEADER_COARSE
READ (IU10) N_NEST_C, MAX_NEST_C 
IF (N_NEST_C.EQ.0) THEN
   WRITE (IU06,*) '********************************************************'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*       FATAL ERROR IN SUB. MAKE_FINE_BOUNDARY         *'
   WRITE (IU06,*) '*       ======================================         *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*  A FINE GRID RUN IS REQUESTED BUT COARSE GRID OUTPUT *'
   WRITE (IU06,*) '*  INFORMATION IS NOT IN THE COARSE PREPROC OUTPUT     *'
   WRITE (IU06,*) '*  COARSE GRID HEADER IS: ', HEADER_COARSE
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*  UNIT IS         IU10 = ', IU10
   WRITE (IU06,*) '*  FILE NAME IS  FILE10 = ', FILE10(1:LEN)
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*              THE PROGRAM ABORTS                      *'
   WRITE (IU06,*) '********************************************************'
   CALL ABORT1
END IF

DO I=1,N_NEST_C
   READ (IU10) NBINP, C_NAME

   IF (ALLOCATED(IJ_C)) DEALLOCATE(IJ_C)
   ALLOCATE(IJ_C(1:NBINP))
   IF (ALLOCATED(LAT_C))  DEALLOCATE(LAT_C)
   ALLOCATE(LAT_C(1:NBINP))
   IF (ALLOCATED(LON_C))  DEALLOCATE(LON_C)
   ALLOCATE(LON_C(1:NBINP))
   IF (ALLOCATED(ZDEL_C))  DEALLOCATE(ZDEL_C)
   ALLOCATE(ZDEL_C(1:NBINP))

   READ (IU10) IJ_C
   READ (IU10) DLAMAC, DPHIAC, SOUTH_C, NORTH_C, EAST_C, WEST_C, LON_C, LAT_C, ZDEL_C
   IF (AMOWEP.EQ.WEST_C .AND. AMOEAP.EQ.EAST_C .AND.                                &
&      AMONOP.EQ.NORTH_C .AND. AMOSOP.EQ.SOUTH_C) EXIT
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK THE INPUT.                                                      !
!     -------------------                                                      !

!     IS THE FINE GRID THE SAME AS IN THE COURSE GRID ?                        !

IF (AMOWEP.NE.WEST_C .OR. AMOEAP.NE.EAST_C .OR.                                &
&    AMONOP.NE.NORTH_C .OR. AMOSOP.NE.SOUTH_C) THEN
   WRITE (IU06,*) '********************************************************'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*       FATAL ERROR IN SUB. MAKE_FINE_BOUNDARY         *'
   WRITE (IU06,*) '*       ======================================         *'
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*  THE GRID IS NOT IN THE COURSE GRID SET-UP           *'
   WRITE (IU06,*) '*  COARSE GRID HEADER IS: ', HEADER_COARSE
   WRITE (IU06,*) '*                                                      *'
   WRITE (IU06,*) '*              THE PROGRAM ABORTS                      *'
   WRITE (IU06,*) '********************************************************'
   WRITE (IU06,*) 'AMOWEP.NE.WEST_C :',AMOWEP,WEST_C 
   WRITE (IU06,*) 'AMOEAP.NE.EAST_C :',AMOEAP,EAST_C 
   WRITE (IU06,*) 'AMONOP.NE.NORTH_C:',AMONOP,NORTH_C
   WRITE (IU06,*) 'AMOSOP.NE.SOUTH_C:',AMOSOP,SOUTH_C
   WRITE (IU06,*) '********************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ALLOCATE ARRAYS AND INITIAL.                                          !
!        ----------------------------                                          !

NBOUNF = (NX+NY)*2-4
IF (.NOT.ALLOCATED(IJARF)) ALLOCATE(IJARF(1:NBOUNF))
IF (.NOT.ALLOCATED(IBFL))  ALLOCATE(IBFL(1:NBOUNF))
IF (.NOT.ALLOCATED(IBFR))  ALLOCATE(IBFR(1:NBOUNF))
IF (.NOT.ALLOCATED(BFW))   ALLOCATE(BFW(1:NBOUNF))
IF (.NOT.ALLOCATED(BLATF)) ALLOCATE(BLATF(1:NBOUNF))
IF (.NOT.ALLOCATED(BLNGF)) ALLOCATE(BLNGF(1:NBOUNF))

IJARF = 0
IBFL  = 0
IBFR  = 0
BFW   = 0.
BLATF = 0.
BLNGF = 0.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTED THE SQUARE BOX.                                              !
!        ------------------------                                              !

CALL MAKE_BOX (AMOWEP, AMOSOP, AMOEAP, AMONOP, NBOUNF, F_NEST_I, F_NEST_K,     &
&              BLATF, BLNGF)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. SEARCH BLOCK NUMBER AND SEA POINT NUMBER.                             !
!        -----------------------------------------                             !

CALL FIND_SEA_POINT (NBOUNF, BLATF, BLNGF, IJARF)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. MAKE INTERPOLATED ARRAYS.                                             !
!        -------------------------                                             !

CALL MINTF (NBOUNF)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. PACK ALL ARRAYS.                                                      !
!        ----------------                                                      !

IF (.NOT.ALLOCATED(MASK)) ALLOCATE(MASK(1:NBOUNF))

MASK = IJARF.GT.0
NBOUNEW = COUNT(MASK)

BLATF(1:NBOUNEW) = PACK (BLATF(1:NBOUNF), MASK(1:NBOUNF))
BLNGF(1:NBOUNEW) = PACK (BLNGF(1:NBOUNF), MASK(1:NBOUNF))
BFW  (1:NBOUNEW) = PACK (BFW  (1:NBOUNF), MASK(1:NBOUNF))
IJARF(1:NBOUNEW) = PACK (IJARF(1:NBOUNF), MASK(1:NBOUNF))
IBFR (1:NBOUNEW) = PACK (IBFR (1:NBOUNF), MASK(1:NBOUNF))
IBFL (1:NBOUNEW) = PACK (IBFL (1:NBOUNF), MASK(1:NBOUNF))

NBOUNF = NBOUNEW
DEALLOCATE (MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. DEALLOCATE ARRAYS.                                                    !
!        ------------------                                                    !

DEALLOCATE (IJ_C, LAT_C, LON_C)

END SUBROUTINE MAKE_FINE_BOUNDARY

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MAKE_BOX (WEST, SOUTH, EAST, NORTH, N_POINT, CORNER_I, CORNER_K,    &
&                    LATITUDE, LONGITUDE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MAKE_BOX - MAKE BOX IN A GRID.                                             !
!                                                                              !
!     R. PORTZ     MPI         15/01/1991                                      !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       COMPUTE LATITUDES AND LONGITUDES OF ALL GRID POINTS AT A BOX BOUNDARY. !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE BOX BOUNDARY IS SCANNED ALONG LATITUDES FROM WEST TO EAST          !
!       MOVING FROM SOUTH TO NORTH.                                            !
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

IMPLICIT NONE

INTEGER, INTENT(IN)    :: WEST                !! WESTERN  LONGITUDE OF BOX.
INTEGER, INTENT(IN)    :: SOUTH               !! SOUTHERN LATITUDE  OF BOX.
INTEGER, INTENT(IN)    :: EAST                !! EASTERN  LONGITUDE OF BOX.
INTEGER, INTENT(IN)    :: NORTH               !! NORTHERN LATITUDE  OF BOX.
INTEGER, INTENT(INOUT) :: N_POINT             !! NUMBER OF BOUNDARY POINTS.
INTEGER, INTENT(OUT)   :: CORNER_I (4)        !! I-INDICES OF CORNER POINTS.
INTEGER, INTENT(OUT)   :: CORNER_K (4)        !! K-INDICES OF CORNER POINTS.
INTEGER, INTENT(OUT)   :: LATITUDE (N_POINT)  !! LATITUDES  OF BOUNDARY POINTS.
INTEGER, INTENT(OUT)   :: LONGITUDE(N_POINT)  !! LONGITUDES OF BOUNDARY POINTS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

integer :: i, k, ls, ln, number

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. FIND NEAREST COARSE GRID POINTS OT CORNER POINTS OF FINE GRID.        !
!        --------------------------------------------------------------        !


CORNER_K(1) = NINT(REAL(SOUTH - AMOSOP)/REAL(XDELLA))+1
CORNER_K(2) = NINT(REAL(SOUTH - AMOSOP)/REAL(XDELLA))+1
CORNER_K(3) = NINT(REAL(NORTH - AMOSOP)/REAL(XDELLA))+1
CORNER_K(4) = NINT(REAL(NORTH - AMOSOP)/REAL(XDELLA))+1

CORNER_I(1) = NINT(REAL(WEST  - AMOWEP)/REAL(ZDELLO(CORNER_K(1))))+1
CORNER_I(2) = NINT(REAL(EAST  - AMOWEP)/REAL(ZDELLO(CORNER_K(2))))+1
CORNER_I(3) = NINT(REAL(WEST  - AMOWEP)/REAL(ZDELLO(CORNER_K(3))))+1
CORNER_I(4) = NINT(REAL(EAST  - AMOWEP)/REAL(ZDELLO(CORNER_K(4))))+1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTED THE SQUARE BOX.                                              !
!        ------------------------                                              !

NUMBER = (CORNER_I(2)-CORNER_I(1)) + (CORNER_I(4)-CORNER_I(3))                 &
&      + (CORNER_K(3)-CORNER_K(1)) + (CORNER_K(4)-CORNER_K(2))


IF (NUMBER.GT.N_POINT) THEN
   WRITE(IU06,*) ' *************************************************'
   WRITE(IU06,*) ' *                                               *'
   WRITE(IU06,*) ' *         FATAL ERROR IN SUB. MAKE_BOX          *'
   WRITE(IU06,*) ' *         ============================          *'
   WRITE(IU06,*) ' *  NUMBER OF BOUNDARY POINTS EXCEEDS DIMENSION. *'
   WRITE(IU06,*) ' *  DIMENSION IS      N_POINT = ', N_POINT
   WRITE(IU06,*) ' *  NUMBER OF POINTS IS         ', NUMBER
   WRITE(IU06,*) ' *                                               *'
   WRITE(IU06,*) ' *************************************************'
   CALL ABORT1
END IF
N_POINT = NUMBER

LATITUDE(1:N_POINT) = 0.
LATITUDE(1:N_POINT) = 0.

LS = CORNER_K(1) 
K = 1
LONGITUDE(1) = WEST
LATITUDE (1) = SOUTH

DO I = CORNER_I(1)+1, CORNER_I(2)-1
   K = K+1
   LATITUDE(K) = SOUTH
   LONGITUDE(K) = AMOWEP+(I-1)*ZDELLO(LS)
END DO
K = K+1
LONGITUDE(K) = EAST
LATITUDE (K) = SOUTH

DO I = CORNER_K(1)+1, CORNER_K(3)-1
   K = K+1
   LATITUDE(K) = AMOSOP+REAL(I-1)*XDELLA
   LONGITUDE(K) = WEST
   K = K+1
   LATITUDE(K) = AMOSOP+REAL(I-1)*XDELLA
   LONGITUDE(K) = EAST
END DO

LN = CORNER_K(3) 
K = K+1
LONGITUDE(K) = WEST
LATITUDE(K) = NORTH

DO I = CORNER_I(3)+1, CORNER_I(4)-1
   K = K+1
   LATITUDE(K) = NORTH
   LONGITUDE(K) = AMOWEP+(I-1)*ZDELLO(LN)
END DO
K = K+1
LONGITUDE(K) = EAST
LATITUDE (K) = NORTH


IF (K.NE.N_POINT) THEN
   WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '+        INFORMATION FROM SUB. MAKE_BOX              +'
   WRITE (IU06,*) '+        ==============================              +'
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '+ NUMBER OF BOUNDARY POINTS DO NOT MATCH             +'
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '+    NUMBER OF POINTS EXPECTED N_POINT = ', N_POINT
   WRITE (IU06,*) '+    NUMBER OF POINTS FOUND          K = ', K
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '+  NEST INFORMATION WILL NOT BE GENERATED            +'
   WRITE (IU06,*) '+                                                    +'
   WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

END SUBROUTINE MAKE_BOX

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MINTF (NBOUNF)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MINTF - MAKE INTERPOLATION TABLES FOR BOUNDARY INPUT.                      !
!                                                                              !
!     R. PORTZ     MPI         15/01/1991                                      !
!     H. GUNTHER   GKSS/ECMWF  15/01/1991                                      !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       GENERATE SPACE INTERPOLATION TABLES USED FOR BOUNDARY VALUE INPUT      !
!       INTO A FINE GRID MODEL.                                                !
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
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

IMPLICIT NONE

INTEGER, INTENT(IN) :: NBOUNF   !! NUMBER OF FINE GRID POINTS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: I, M, M1, M2, ID, ID1, ID2
REAL     :: D

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. LOOP OVER FINE GRID POINTS.                                           !
!         --------------------------                                           !

FINE: DO I = 1, NBOUNF
   M1 = 0
   M2 = 0
   ID1 = 99999999

!     1.1 FINE NEAREST COARSE GRID POINT.                                      !

   COARSE1: DO M = 1,NBINP
      IF (LAT_C(M) .NE. BLATF(I) .AND.  LON_C(M) .NE. BLNGF(I) ) THEN
          CYCLE COARSE1
      ELSE IF ((LAT_C(M) .EQ. BLATF(I) .AND.  LON_C(M) .EQ. BLNGF(I) )) THEN
         ID1 = 0
         M1 = M
         IBFL(I) = M
         IBFR(I) = M   
         BFW(I) = 1.
         CYCLE FINE        !! TAKE NEXT FINE GP
      ELSE IF (LAT_C(M) .EQ. BLATF(I)) THEN
          ID = ABS(LON_C(M) - BLNGF(I))
          IF (ID.LT.ID1) THEN   
             ID1 = ID
             M1 = M
          END IF
      ELSE IF (LON_C(M) .EQ. BLNGF(I)) THEN
          ID = ABS(LAT_C(M) - BLATF(I))
          IF (ID.LT.ID1) THEN   
             ID1 = ID
             M1 = M
          END IF
      END IF
   END DO COARSE1

!     1.2 FINE GP AND NEAREST COARSE GP HAVE DIFFERENT LONGITUDES AND LATITUDE.                                           !
     
   IF (M1.EQ.0 ) CYCLE FINE 

!     1.3 FINE GP AND NEAREST COARSE GP HAVE THE SAME LATITUDE.                                           !

   M2 = 0
   IF (LAT_C(M1) .EQ. BLATF(I)) THEN

!     1.3.1 NEAREST COARSE GP IS TO FAR AWAY.                                           !

      IF (ID1 .GT. ZDEL_C(M1)) CYCLE FINE

!     1.3.2 FINED NEAREST COARSE GP OPPOSITE TO FINE GP ON THE SAME LATITUDE.

      ID2 = NINT(1.5 * REAL(ZDEL_C(M1))) - ID1
      ID1 = LON_C(M1)-BLNGF(I)
      DO M = 1, NBINP
         IF (M.EQ.M1) CYCLE
         IF (LAT_C(M) .NE. BLATF(I)) CYCLE
         ID = LON_C(M)-BLNGF(I)
         IF (SIGN(1,ID1).EQ.SIGN(1,ID)) CYCLE
         IF (ABS(ID).LT.ID2) THEN
            ID2 = ABS(ID)
            M2 = M
         END IF
      END DO

!     1.3.3 CHECK SECOND POINT
       
      IF (M2.EQ.0) THEN
         M2 = 0
         D = REAL(ABS(ID1))/REAL(ZDEL_C(M1))
       ELSE
         ID1 = ABS(ID1)
         D = REAL(ID1)/REAL(ID1+ID2)
      END IF
       
   END IF

!     1.4 FINE GP AND NEAREST COARSE GP HAVE THE SAME LONGITUDE.                                           !

   IF (LON_C(M1) .EQ. BLNGF(I)) THEN

!     1.4.1 NEAREST COARSE GP IS TO FAR AWAY.                                           !

      IF (ID1 .GT. DPHIAC) CYCLE FINE

!     1.4.2 FINED NEAREST COARSE GP OPPOSITE TO FINE GP ON THE SAME LONGITUDE.


      ID2 = NINT(1.5 * REAL(DPHIAC)) - ID1
      ID1 = LAT_C(M1)-BLATF(I)
      DO M = 1, NBINP
         IF (M.EQ.M1) CYCLE
         IF (LON_C(M) .NE. BLNGF(I)) CYCLE
         ID = LAT_C(M)-BLATF(I)
         IF (SIGN(1,ID1).EQ.SIGN(1,ID)) CYCLE
         IF (ABS(ID).LT.ID2) THEN
             ID2 = ABS(ID)
             M2 = M
         END IF
      END DO

!     1.4.3 CHECK SECOND POINT

      IF (M2.EQ.0) THEN
         M2 = 0
         D = REAL(ABS(ID1))/REAL(DPHIAC)
      ELSE
         ID1 = ABS(ID1)
         D = REAL(ID1)/REAL(ID1+ID2)
      END IF
   END IF

!     1.5 DEFINE INTERPOLATION INDEICES AND WEIGHT.

   IF (M1.GT.M2) THEN
      IBFL(I) = M2   
      IBFR(I) = M1   
      BFW(I)  = 1. - D
   ELSE
      IBFL(I) = M1   
      IBFR(I) = M2   
      BFW(I)  = D
   END IF   

END DO FINE

END SUBROUTINE MINTF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_NEST_MODULE

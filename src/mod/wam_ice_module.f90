MODULE WAM_ICE_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE STORES THE ICE INPUT GRID SPECFICATIONS, AND THE MODEL ICE     !
!   ICE INFORMATION AND THE ICE PROCESSING PROGRAMMS.                          !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS FROM OTHER MODULES.                                        !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

USE WAM_COORDINATE_MODULE          !! COORDINATE PROCEDURES

USE WAM_GENERAL_MODULE, ONLY:  &
&       ABORT1,                &   !! TERMINATES PROCESSING.
&       INCDATE,               &   !! INCREMENTS DATE TIME GROUP.
&       PRINT_ARRAY                !! PRINTS AN ARRAY.

USE WAM_GRID_MODULE,      ONLY:  &
&       EQUAL_TO_M_GRID            !! COMPARES WIND GRID TO MODEL GRID.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST
USE WAM_GRID_MODULE,    ONLY: NX, NY, NSEA, L_S_MASK, NLON_RG, ZDELLO,         &
&                             AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLO, XDELLA
USE WAM_TIMOPT_MODULE,  ONLY: CDTPRO
USE WAM_FRE_DIR_MODULE, ONLY: ML, KL

use wam_mpi_module,     only: nijs, nijl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C.  MODULE DATA.                                                         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

IMPLICIT NONE
PRIVATE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INPUT ICE GRID SPECFICATIONS AND ICE FIELD.                           !
!        -------------------------------------------                           !

INTEGER   :: NX_IN  =-1      !! NUMBER OF LONGITUDES.
INTEGER   :: NY_IN  =-1      !! NUMBER OF LATITUDES.
LOGICAL   :: PER =.FALSE.    !! .TRUE. IF PERIODIC GRID.
INTEGER   :: DX_IN  =-1.     !! STEPSIZE BETWEEN LONGITUDES [M_SEC].
INTEGER   :: DY_IN  =-1.     !! STEPSIZE BETWEEN LATITUDES  [M_SEC].
INTEGER   :: SOUTH_IN =-1.   !! SOUTH LATITUDE OF GRID [M_SEC].
INTEGER   :: NORTH_IN =-1.   !! NORTH LATITUDE OF GRID [M_SEC].
INTEGER   :: WEST_IN =-1.    !! WEST LONGITUDE OF GRID [M_SEC].
INTEGER   :: EAST_IN =-1.    !! EAST LONGITUDE OF GRID [M_SEC].
LOGICAL   :: EQUAL_GRID =.FALSE. !! .TRUE. IF WIND GRID IS EQUAL TO MODEL GRID.

LOGICAL, ALLOCATABLE :: ICE_MAP_IN(:,:) !! ICE MAP.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. ICE POINTS.                                                           !
!        -----------                                                           !

LOGICAL              :: ICE_RUN = .FALSE. !! TRUE IF ICE IS TAKEN INTO ACCOUNT
CHARACTER (LEN=14)   :: CD_ICE      = ' ' !! DATE OF ICE FIELD USED IN WAM.
CHARACTER (LEN=14)   :: CD_ICE_NEW  = ' ' !! PROPAGATION DATE FOR NEXT ICE FIELD
INTEGER              :: N_ICE   = 0       !! NUMBER OF ICE POINTS.
INTEGER, ALLOCATABLE :: IJ_ICE(:)         !! INDEX OF ICE POINTS.

PUBLIC ICE_RUN, CD_ICE_NEW

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ICE INPUT TIMESTEP.                                                   !
!        -------------------                                                   !

INTEGER   :: IDEL_ICE_I = 0  !! INPUT ICE TIMESTEP IN SECONDS.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTERFACE GET_ICE                       !! GET ICE FIELDS FROM FILE.
   MODULE PROCEDURE GET_ICE 
END INTERFACE
PUBLIC GET_ICE

INTERFACE PRINT_ICE_STATUS              !! PRINTS STATUS OF THIS MODULE.
   MODULE PROCEDURE PRINT_ICE_STATUS
END INTERFACE
PUBLIC PRINT_ICE_STATUS

INTERFACE PUT_ICE                       !! MARKS ICE POINTS
   MODULE PROCEDURE PUT_ICE_PAR         !! IN AN INTEGRATED PARAMETER ARRAY.
   MODULE PROCEDURE PUT_ICE_PAR_B       !! IN AN INTEGRATED PARAMETER BLOCK.
   MODULE PROCEDURE PUT_ICE_SPEC        !! IN SPECTRA ARRAY.
END INTERFACE
PUBLIC PUT_ICE

INTERFACE SET_ICE                       !! SETS A NEW ICE MAP
   MODULE PROCEDURE SET_ICE_I           !! INTEGER VERSION
   MODULE PROCEDURE SET_ICE_L           !! LOGICAL VERSION
   MODULE PROCEDURE SET_ICE_R           !! REAL VERSION
END INTERFACE
PUBLIC SET_ICE 

INTERFACE SET_ICE_HEADER                !! SETS ICE HEADER.
   MODULE PROCEDURE SET_ICE_HEADER_C    !! CHARACTER VERSION
   MODULE PROCEDURE SET_ICE_HEADER_D    !! DEGREE VERSION
   MODULE PROCEDURE SET_ICE_HEADER_M    !! M_SEC VERSION
END INTERFACE
PUBLIC SET_ICE_HEADER

INTERFACE SET_ICE_TIMESTEP              !! SETS ICE TIMESTEP.
   MODULE PROCEDURE SET_ICE_TIMESTEP 
END INTERFACE
PUBLIC SET_ICE_TIMESTEP

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  EXPLICIT INTERFACES.                                                 !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTERFACE
   SUBROUTINE READ_ICE_INPUT            !! READ AN ICE MAP.
   END SUBROUTINE READ_ICE_INPUT
END INTERFACE 

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE GET_ICE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   GET_ICE - ROUTINE TO READ AND ONE ICE FIELD.                               !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        READ AN ICE FIELD FROM THE ICE FILE (SEARCH FOR IT).                  !
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
!     1. READ ICE DATA AND CHECK DATE.                                         !
!        ------------------------------                                        !

DO
   CALL READ_ICE_INPUT
   IF (IDEL_ICE_I.LE.0) THEN
       CD_ICE_NEW = '999912315959'
       EXIT
   ELSE
      IF (CD_ICE.GE.CDTPRO) THEN
         CD_ICE_NEW = CD_ICE
         CALL INCDATE(CD_ICE_NEW,IDEL_ICE_I/2)
         EXIT
      END IF
   END IF
END DO

IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '  '
   WRITE (IU06,*) '   SUB. GET_ICE: NEW ICE FIELDS PASSED TO WAM MODEL'
   WRITE (IU06,*) '     DATE OF ICE FIELD IS ... CD_ICE = ', CD_ICE
   WRITE (IU06,*) '     NEXT ICE FIELD AT .. CD_ICE_NEW = ', CD_ICE_NEW
   WRITE (IU06,*) '     NO. OF ICEPOINTS IS ..... N_ICE = ', N_ICE
   WRITE (IU06,*) '  '
END IF

END SUBROUTINE GET_ICE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_ICE_STATUS

INTEGER :: I
CHARACTER (LEN=100) :: TITL
CHARACTER (LEN=1)   :: LST(NX,NY) !! LAND SEA TABLE  L = LAND, S = SEA, I = ICE.
CHARACTER (LEN=14)  :: ZERO = ' '
CHARACTER (LEN=1)   :: PMSK(1:NSEA)
character (len=len_coor) :: formtext

WRITE (IU06,*) '  '
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '              ICE MODULE STATUS:'
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '  '
IF (ICE_RUN) THEN
   IF (NX_IN.GT.0 .AND. NY_IN.GT.0) THEN
      WRITE (IU06,'('' ICE INPUT GRID SPECIFICATION ARE: '')')
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
      WRITE (IU06,*) '  '
   ELSE
      WRITE (IU06,*) ' ICE INPUT GRID IS NOT DEFINED'
      WRITE (IU06,*) '  '
   END IF

   WRITE (IU06,*) ' DATE OF ICE FIELD IS............... CD_ICE = ', CD_ICE
   WRITE (IU06,*) ' ICE INITIALISED, NO. OF ICEPOINTS IS N_ICE = ', N_ICE
   IF (IDEL_ICE_I.LE.0) THEN
      WRITE (IU06,*) ' ICE IS KEPT CONSTANT FOR THE FULL RUN '
   ELSE
      WRITE (IU06,*) ' ICE IS UP-DATED EVERY...........IDEL_ICE_I = ',         &
&                      IDEL_ICE_I,' SECONDS'
      WRITE (IU06,*) ' NEXT ICE FIELD AT...............CD_ICE_NEW = ',         &
&                      CD_ICE_NEW
   END IF
   IF (ITEST.GT.0 .AND. N_ICE.GT.0) THEN
      PMSK =  'S'                               !! SEA POINTS
      PMSK(IJ_ICE) = 'I'                        !! ICE POINTS
      LST = UNPACK(PMSK, L_S_MASK, 'L')         !! LAND POINTS
      DO I=1,NY
          LST(NLON_RG(I)+1:NX,I) = ' '
      END DO
      WRITE (IU06,*) ' '
      TITL = ' LAND SEA MAP:  L = LAND, S = SEA,  I = ICE '
      CALL PRINT_ARRAY (IU06, ZERO, TITL, LST, AMOWEP, AMOSOP, AMOEAP, AMONOP)
   END IF
ELSE
   WRITE (IU06,*) ' ICE IS NOT INITIALIZED '
END IF
WRITE (IU06,*) '  '

END SUBROUTINE PRINT_ICE_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PUT_ICE_PAR (PAR, INDICATOR)

REAL, INTENT(INOUT) :: PAR(:)      !! BLOCK OF PARAMETER.
REAL, INTENT(IN)    :: INDICATOR   !! VALUE TO BE INSERTED AT ICE POINTS.

IF (N_ICE.GT.0) PAR(IJ_ICE) = INDICATOR

END SUBROUTINE PUT_ICE_PAR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PUT_ICE_PAR_B (PAR, NA, NE, INDICATOR)

INTEGER, INTENT(IN) :: NA
INTEGER, INTENT(IN) :: NE
REAL, INTENT(INOUT) :: PAR(NA:NE)  !! BLOCK OF PARAMETER.
REAL, INTENT(IN)    :: INDICATOR   !! VALUE TO BE INSERTED AT ICE POINTS.

INTEGER IJ

DO IJ = 1, N_ICE
   IF (IJ_ICE(IJ).GE.NA .AND. IJ_ICE(IJ).LE.NE)                            &
&        PAR(IJ_ICE(IJ)) = INDICATOR
END DO

END SUBROUTINE PUT_ICE_PAR_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PUT_ICE_SPEC (FL3, INDICATOR)

REAL, INTENT(INOUT) :: FL3(nijs:,:,:)   !! BLOCK OF SPECTRA.
REAL, INTENT(IN)    :: INDICATOR        !! VALUE TO BE INSERTED AT ICE POINTS.

INTEGER IJ

DO IJ = 1, N_ICE
   IF (IJ_ICE(IJ).GE.nijs .AND. IJ_ICE(IJ).LE.nijl)                            &
&        FL3(IJ_ICE(IJ),:,:) = INDICATOR
END DO

END SUBROUTINE PUT_ICE_SPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_ICE_I (CDT, I_GRID)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_ICE_I - TRANSFERS ICE DATA TO ICE MODULE.  (INTEGER VERSION)           !
!                                                                              !
!     HEINZ GUNTHER    GKSS    JANUARY 2002                                    !
!                                                                              !
!     PURPOSE                                                                  !
!     -------                                                                  !
!                                                                              !
!       READ AN ICE MAP AND BLOCK THE INFORMATION.                             !
!                                                                              !
!     METHOD                                                                   !
!     ------                                                                   !
!                                                                              !
!        NONE.                                                                 !
!                                                                              !
!     REFERENCES                                                               !
!     ----------                                                               !
!                                                                              !
!          NONE                                                                !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

CHARACTER (LEN=14),INTENT(IN) :: CDT          !! DATE/TIME OF ICE FIELD.
INTEGER,           INTENT(IN) :: I_GRID(:,:)  !! ICE MAP

INTEGER :: I, K, I_ICE, K_ICE, IJ, N
LOGICAL :: BLOCK(1:NSEA)
LOGICAL :: ICE_MAP(1:NX,1:NY) 

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CHECK DIMENSIONS OF ICE INFORMATION.                                  !
!        ------------------------------------                                  !

IF (SIZE(I_GRID,1).NE.NX_IN .OR. SIZE(I_GRID,2).NE.NY_IN) THEN
   WRITE(IU06,*) ' *******************************************'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *   FATAL ERROR IN SUB. SET_ICE_INPUT     *'
   WRITE(IU06,*) ' *   =================================     *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' * ICE MAP SIZE IS NOT EQUAL TO GRID       *'
   WRITE(IU06,*) ' * DEFINE WITH SUB. SET_ICE_HEADER.        *'
   WRITE(IU06,*) ' * DINENSIONS DEFINED BS HEADER IS:        *'
   WRITE(IU06,*) ' *           NX_IN  = ', NX_IN
   WRITE(IU06,*) ' *           NY_IN  = ', NY_IN
   WRITE(IU06,*) ' * ICE MAP DIMENSIONS ARE       *'
   WRITE(IU06,*) ' *   SIZE(I_GRID,1) = ', SIZE(I_GRID,1)
   WRITE(IU06,*) ' *   SIZE(I_GRID,2) = ', SIZE(I_GRID,2)
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
   WRITE(IU06,*) ' *                                         *'
   WRITE(IU06,*) ' *******************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. STORE ICE MAP AND DATE.                                               !
!        -----------------------                                               !

CD_ICE = CDT

IF (ALLOCATED(ICE_MAP_IN)) DEALLOCATE(ICE_MAP_IN)
ALLOCATE (ICE_MAP_IN(1:NX_IN,1:NY_IN))
ICE_MAP_IN = I_GRID.EQ.1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. INTERPOLATE INPUT TO MODEL GRID (NEAREST GRIDPOINT).                  !
!        ----------------------------------------------------                  !

ICE_MAP = .FALSE.

K_LOOP: DO K=1,NY
   K_ICE = NINT(REAL(AMOSOP+(K-1)*XDELLA-SOUTH_IN)/REAL(DY_IN)+1.0)
   IF (K_ICE.LT.1 .OR. K_ICE.GT. NY_IN) CYCLE K_LOOP
   I_LOOP: DO I=1,NLON_RG(K)
      IF (.NOT. L_S_MASK(I,K)) CYCLE I_LOOP
      IF (PER) THEN
         I_ICE = NINT( REAL(AMOWEP + (I-1)*ZDELLO(K)+M_S_PER - WEST_IN)/REAL(DX_IN))
         I_ICE = MOD(I_ICE+NX_IN,NX_IN) + 1
      ELSE
         I_ICE = NINT( REAL(AMOWEP + (I-1)*ZDELLO(K) - WEST_IN)/REAL(DX_IN) ) + 1
      END IF
      IF (I_ICE.LT.1 .OR. I_ICE.GT. NX_IN) CYCLE I_LOOP
      ICE_MAP(I,K) = ICE_MAP_IN(I_ICE,K_ICE)   
    END DO I_LOOP
END DO K_LOOP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. SMOOTH ICE GRID AT COAST LINES.                                       !
!        -------------------------------                                       !

K_LOOP1: DO K=2,NY-1
   I_LOOP1: DO I=2,NX-1
      IF (.NOT. L_S_MASK(I,K)) CYCLE I_LOOP1
      IF (ICE_MAP(I,K)) CYCLE I_LOOP1
      ICE_MAP(I,K) = (ICE_MAP(I+1,K).AND..NOT.L_S_MASK(I-1,K)) .OR.            &
&                    (ICE_MAP(I-1,K).AND..NOT.L_S_MASK(I+1,K)) .OR.            & 
&                    (ICE_MAP(I,K+1).AND..NOT.L_S_MASK(I,K-1)) .OR.            & 
&                    (ICE_MAP(I,K-1).AND..NOT.L_S_MASK(I,K+1)) 
    END DO I_LOOP1
END DO K_LOOP1

K = 1
I_LOOP2: DO I=2,NX-1
   IF (.NOT. L_S_MASK(I,K)) CYCLE I_LOOP2
   IF (ICE_MAP(I,K)) CYCLE I_LOOP2
   ICE_MAP(I,K) = (ICE_MAP(I+1,K).AND..NOT.L_S_MASK(I-1,K)) .OR.               &
&                 (ICE_MAP(I-1,K).AND..NOT.L_S_MASK(I+1,K)) .OR.               & 
&                  ICE_MAP(I,K+1)
END DO I_LOOP2

K = NY
I_LOOP3: DO I=2,NX-1
   IF (.NOT. L_S_MASK(I,K)) CYCLE I_LOOP3
   IF (ICE_MAP(I,K)) CYCLE I_LOOP3
   ICE_MAP(I,K) = (ICE_MAP(I+1,K).AND..NOT.L_S_MASK(I-1,K)) .OR.               &
&                 (ICE_MAP(I-1,K).AND..NOT.L_S_MASK(I+1,K)) .OR.               & 
&                  ICE_MAP(I,K-1) 
END DO I_LOOP3

IF (.NOT.PER) THEN 
   I = 1
   K_LOOP4: DO K=2,NY-1
      IF (.NOT. L_S_MASK(I,K)) CYCLE K_LOOP4
      IF (ICE_MAP(I,K)) CYCLE K_LOOP4
      ICE_MAP(I,K) =  ICE_MAP(I+1,K) .OR.                                      &
&                    (ICE_MAP(I,K+1).AND..NOT.L_S_MASK(I,K-1)) .OR.            & 
&                    (ICE_MAP(I,K-1).AND..NOT.L_S_MASK(I,K+1)) 
   END DO K_LOOP4

   I = NX
   K_LOOP5: DO K=2,NY-1
      IF (.NOT. L_S_MASK(I,K)) CYCLE K_LOOP5
      IF (ICE_MAP(I,K)) CYCLE K_LOOP5
      ICE_MAP(I,K) =  ICE_MAP(I-1,K) .OR.                                      &
&                    (ICE_MAP(I,K+1).AND..NOT.L_S_MASK(I,K-1)) .OR.            & 
&                    (ICE_MAP(I,K-1).AND..NOT.L_S_MASK(I,K+1)) 
   END DO K_LOOP5

   I=1
   K=1
   IF (L_S_MASK(I,K) .AND. .NOT.ICE_MAP(I,K)) THEN
      ICE_MAP(I,K) = ICE_MAP(I+1,K) .AND. ICE_MAP(I,K+1)  
   END IF   
   I=1
   K=NY
   IF (L_S_MASK(I,K) .AND. .NOT.ICE_MAP(I,K)) THEN
      ICE_MAP(I,K) = ICE_MAP(I+1,K) .AND. ICE_MAP(I,K-1)  
   END IF   
   I=NX
   K=1
   IF (L_S_MASK(I,K) .AND. .NOT.ICE_MAP(I,K)) THEN
      ICE_MAP(I,K) = ICE_MAP(I-1,K) .AND. ICE_MAP(I,K+1)  
   END IF   
   I=NX
   K=NY
   IF (L_S_MASK(I,K) .AND. .NOT.ICE_MAP(I,K)) THEN
      ICE_MAP(I,K) = ICE_MAP(I-1,K) .AND. ICE_MAP(I,K-1)  
   END IF   
ELSE
   I = 1
   K_LOOP6: DO K=2,NY-1
      IF (.NOT. L_S_MASK(I,K)) CYCLE K_LOOP6
      IF (ICE_MAP(I,K)) CYCLE K_LOOP6
      ICE_MAP(I,K) = (ICE_MAP(I+1,K).AND..NOT.L_S_MASK( NX,K)) .OR.            &
&                    (ICE_MAP( NX,K).AND..NOT.L_S_MASK(I+1,K)) .OR.            & 
&                    (ICE_MAP(I,K+1).AND..NOT.L_S_MASK(I,K-1)) .OR.            & 
&                    (ICE_MAP(I,K-1).AND..NOT.L_S_MASK(I,K+1)) 
   END DO K_LOOP6

   I = NX
   K_LOOP7: DO K=2,NY-1
      IF (.NOT. L_S_MASK(I,K)) CYCLE K_LOOP7
      IF (ICE_MAP(I,K)) CYCLE K_LOOP7
      ICE_MAP(I,K) = (ICE_MAP(  1,K).AND..NOT.L_S_MASK(I-1,K)) .OR.            &
&                    (ICE_MAP(I-1,K).AND..NOT.L_S_MASK(  1,K)) .OR.            & 
&                    (ICE_MAP(I,K+1).AND..NOT.L_S_MASK(I,K-1)) .OR.            & 
&                    (ICE_MAP(I,K-1).AND..NOT.L_S_MASK(I,K+1)) 
   END DO K_LOOP7
   I=1
   K=1
   IF (L_S_MASK(I,K) .AND. .NOT.ICE_MAP(I,K)) THEN
      ICE_MAP(I,K) = (ICE_MAP(I+1,K).AND..NOT.L_S_MASK( NX,K)) .OR.            &
&                    (ICE_MAP( NX,K).AND..NOT.L_S_MASK(I+1,K)) .OR.            & 
&                     ICE_MAP(I,K+1) 
   END IF   
   I=1
   K=NY
   IF (L_S_MASK(I,K) .AND. .NOT.ICE_MAP(I,K)) THEN
      ICE_MAP(I,K) = (ICE_MAP(I+1,K).AND..NOT.L_S_MASK( NX,K)) .OR.            &
&                    (ICE_MAP( NX,K).AND..NOT.L_S_MASK(I+1,K)) .OR.            & 
&                     ICE_MAP(I,K-1) 
      ICE_MAP(I,K) = ICE_MAP(I+1,K) .AND. ICE_MAP(I,K-1)  
   END IF   
   I=NX
   K=1
   IF (L_S_MASK(I,K) .AND. .NOT.ICE_MAP(I,K)) THEN
      ICE_MAP(I,K) = (ICE_MAP(  1,K).AND..NOT.L_S_MASK(I-1,K)) .OR.            &
&                    (ICE_MAP(I-1,K).AND..NOT.L_S_MASK(  1,K)) .OR.            & 
&                     ICE_MAP(I,K+1) 
   END IF   
   I=NX
   K=NY
   IF (L_S_MASK(I,K) .AND. .NOT.ICE_MAP(I,K)) THEN
      ICE_MAP(I,K) = (ICE_MAP(  1,K).AND..NOT.L_S_MASK(I-1,K)) .OR.            &
&                    (ICE_MAP(I-1,K).AND..NOT.L_S_MASK(  1,K)) .OR.            & 
&                     ICE_MAP(I,K-1) 
   END IF   
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. FIND POINT NUMBERS OF ICE POINTS.                                     !
!        ---------------------------------                                     !

BLOCK = PACK(ICE_MAP, L_S_MASK)
N_ICE = COUNT(BLOCK)
IF (ALLOCATED(IJ_ICE)) DEALLOCATE(IJ_ICE)

IF (N_ICE.GT.0) THEN
   ALLOCATE (IJ_ICE(1:N_ICE))
   N = 0
   DO IJ = 1,NSEA
      IF (BLOCK(IJ)) THEN
         N = N+1
         IJ_ICE(N) = IJ
      END IF
   END DO
END IF

END SUBROUTINE SET_ICE_I

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_ICE_L (CDT, L_GRID)  !! (LOGICAL VERSION OF SUB. SET_ICE)

CHARACTER (LEN=14),INTENT(IN) :: CDT          !! DATE/TIME OF ICE FIELD.
LOGICAL,           INTENT(IN) :: L_GRID(:,:)  !! ICE MAP

INTEGER :: I_GRID(SIZE(L_GRID,1),SIZE(L_GRID,2))

WHERE (L_GRID)
   I_GRID = 1
ELSEWHERE
   I_GRID = 0
END WHERE
CALL SET_ICE_I (CDT, I_GRID)

END SUBROUTINE SET_ICE_L

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_ICE_R (CDT, R_GRID)  !! (REAL VERSION OF SUB. SET_ICE)

CHARACTER (LEN=14),INTENT(IN) :: CDT          !! DATE/TIME OF ICE FIELD.
REAL,              INTENT(IN) :: R_GRID(:,:)  !! ICE MAP

CALL SET_ICE_I (CDT, NINT(R_GRID))

END SUBROUTINE SET_ICE_R

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_ICE_HEADER_C (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,           &
&                            N_LON, N_LAT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_ICE_HEADER - SET INPUT GRID FOR ICE FIELDS IN WAM_ICE MODULE.          !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS FOR ICE FIELDS TO WAM_ICE_MODULE.  !
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
!                                                                              !
!     0. CLEAR GRID DEFINITIONS.                                               !
!        ----------------------_                                               !

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

CALL SET_ICE_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,             &
&                      NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,         &
&                      N_LON=N_LON_CO, N_LAT=N_LAT_CO)

END SUBROUTINE SET_ICE_HEADER_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_ICE_HEADER_D (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,           &
&                            N_LON, N_LAT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_ICE_HEADER - SET INPUT GRID FOR ICE FIELDS IN WAM_ICE MODULE.          !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS FOR ICE FIELDS TO WAM_ICE_MODULE.  !
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
REAL (KIND=KIND_D), INTENT(IN), OPTIONAL :: D_LAT  !! LAT.  INC. OF GRID [DEG].
INTEGER,            INTENT(IN), OPTIONAL :: N_LON  !! NUMBER OF LONG. IN GRID.
INTEGER,            INTENT(IN), OPTIONAL :: N_LAT  !! NUMBER OF LAT. IN GRID.

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
!                                                                              !
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

CALL SET_ICE_HEADER_M (WEST=WEST_CO, SOUTH=SOUTH_CO, EAST=EAST_CO,             &
&                      NORTH=NORTH_CO, D_LON=D_LON_CO, D_LAT=D_LAT_CO,         &
&                      N_LON=N_LON_CO, N_LAT=N_LAT_CO)

END SUBROUTINE SET_ICE_HEADER_D

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_ICE_HEADER_M (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT,           &
&                            N_LON, N_LAT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_ICE_HEADER - SET INPUT GRID FOR ICE FIELDS IN WAM_ICE MODULE.          !
!                                                                              !
!     H. GUENTHER  GKSS  JANUARY 2002                                          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE INPUT GRID DEFINITIONS FOR ICE FIELDS TO WAM_ICE_MODULE.  !
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

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

LOGICAL  :: ERROR                         !! ERROR FLAG

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. CLEAR GRID DEFINITIONS.                                               !
!        ----------------------_                                               !

NX_IN    = -1           !! NUMBER OF COLUMNES IN TOPO INPUT GRID.
NY_IN    = -1           !! NUMBER OF ROWS      IN TOPO INPUT GRID.
PER      =.FALSE.       !! .TRUE. IF PERIODIC GRID.
DX_IN    = COOR_UNDEF   !! STEPSIZE BETWEEN LONGITUDES.
DY_IN    = COOR_UNDEF   !! STEPSIZE BETWEEN LATITUDES.
SOUTH_IN = COOR_UNDEF   !! MOST SOUTHERN LATITUDE.
NORTH_IN = COOR_UNDEF   !! MOST NORTHERN LATITUDE.
WEST_IN  = COOR_UNDEF   !! LEFT MOST LONGITUDE.
EAST_IN  = COOR_UNDEF   !! RIGHT MOST LONGITUDE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COPY GRID DEFINITION.                                                 !
!        ---------------------                                                 !

WEST_IN  = WEST
SOUTH_IN = SOUTH
IF (PRESENT(EAST))  EAST_IN  = EAST
IF (PRESENT(NORTH)) NORTH_IN = NORTH
IF (PRESENT(D_LON)) DX_IN    = D_LON
IF (PRESENT(D_LAT)) DY_IN    = D_LAT
IF (PRESENT(N_LON)) NX_IN   =  N_LON
IF (PRESENT(N_LAT)) NY_IN   =  N_LAT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CHECK GRID.                                                           !
!        -----------                                                           !

CALL CHECK_GRID_DEFINITION (WEST_IN, SOUTH_IN, EAST_IN, NORTH_IN,              &
&                           DX_IN, DY_IN, NX_IN, NY_IN, ERROR)

IF (ERROR) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *          FATAL ERROR IN SUB. SET_ICE_HEADER            *'
   WRITE (IU06,*) ' *          ==================================            *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *          A ICE GRID COULD NOT BE DEFINED.              *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
END IF


PER = PERIODIC (WEST_IN, EAST_IN, DX_IN, NX_IN)   !! PERIODIC?
EQUAL_GRID = EQUAL_TO_M_GRID (NX_IN, NY_IN, DX_IN, DY_IN,                      &
&                             WEST_IN, SOUTH_IN, EAST_IN, NORTH_IN)

END SUBROUTINE SET_ICE_HEADER_M

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_ICE_TIMESTEP (IN)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_ICE_TIMESTEP - SETS ICE TIMESTEP IN WAM_ICE_MODULE.                    !
!                                                                              !
!     H. GUENTHER  GKSS  DECEMBER 2009                                         !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TRANSFER THE ICE INPUT TIMESTEP TO WAM_ICE_MODULE.                     !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       COPY DEFINITION.                                                       !
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

INTEGER, INTENT(IN)           :: IN    !! ICE INPUT TIME STEP.

IDEL_ICE_I = IN

END SUBROUTINE SET_ICE_TIMESTEP

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_ICE_MODULE

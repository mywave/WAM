MODULE WAM_TABLES_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL PRE-COMPUTED TABLES USED IN THE WAM MODEL AND     !
!   ALL VARIABLES, CONSTANTS AND PROCEDURES TO COMPUTE TABLES.                 !
!   THESE TABLES ARE FOR:                                                      !
!     MINIMUM ENERGY IN SPECTRAL BINS  (U10, FREQUENCY)                        !
!     SHALLOW WATER GROUP VELOCITY     (DEPTH, FREQUENCY)                      !
!     SHALLOW WATER WAVE NUMBER        (DEPTH, FREQUENCY)                      !
!     OMEGA/SINH(2KD)                  (DEPTH, FREQUENCY)                      !
!     WAVE NUMBER**(-3)/GROUP VELOCITY (DEPTH, FREQUENCY)                      !
!                                                                              !
!    JUNE 2005:                                                                !
!       Changes introduced as documented in:                                   !
!       (Jean Bidlot, Peter Janssen and Saleh Abdalla: A revised formulation   !
!       for ocean wave dissipation in CY29R1. ECMWF MEMORANDUM RESEARCH        !
!       DEPARTMENT:April 7, 2005 File: R60.9/JB/0516                           !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE, ONLY:   &
&      AKI,                     & !! COMPUTE WAVENUMBER (DEPTH, FREQUENCY)
&      ABORT1                     !! TERMINATES PROCESSING.
USE WAM_JONSWAP_MODULE, ONLY:   &
&      FETCH_LAW,               & !! COMPUTE JONSWAP PARAMETERS FROM FETCH LAW.
&      JONSWAP                    !! COMPUTE JONSWAP SPECTRA.
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST
USE WAM_FRE_DIR_MODULE, ONLY: ML, FR, FRM5
USE WAM_GENERAL_MODULE, ONLY: G, PI, ZPI
USE WAM_GRID_MODULE,    ONLY: DELPHI

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

! ---------------------------------------------------------------------------- !
!
!    1. TABLE FOR MINIMUM ENERGY IN SPECTRAL BINS.
!       --------------------------------------------

REAL, PARAMETER :: FETCH_MAX = 25000.!! MAXIMUM FETCH USED FOR MIN ENEGRY.
REAL, PARAMETER :: FLMIN = 0.000001  !! ABSOLUTE MINIMUM ENERGY IN SPECTRAL BINS
REAL, ALLOCATABLE, DIMENSION(:,:) ::  FLMINFR !! THE MINIMUM VALUE IN SPECTRAL
                                              !! BINS FOR A GIVEN FREQUENCY.

INTEGER, PARAMETER :: JUMAX = 300   !! TABLE DIMENSION FOR U10.
REAL,    PARAMETER :: UMAX  = 75.   !! MAXIMUM WIND SPEED IN TABLE.
REAL               :: DELU  = -1.   !! WIND INCREMENT.

REAL,    PARAMETER :: EPS1 = 0.00001 !! SMALL NUMBER TO MAKE SURE THAT A
                                     !! SOLUTION IS OBTAINED IN ITERATION
                                     !! WITH TAU>TAUW.
PUBLIC :: JUMAX, DELU, UMAX, EPS1
PUBLIC :: FLMIN, FLMINFR

! ---------------------------------------------------------------------------- !
!                                                                              !
!    2. SHALLOW WATER TABLES.                                                  !

INTEGER :: NDEPTH    !! LENGTH OF SHALLOW WATER TABLES.
REAL    :: DEPTHA    !! MINIMUM DEPTH FOR TABLES [M].
REAL    :: DEPTHD    !! DEPTH RATIO.
REAL    :: DEPTHE    !! MAXIMUM DEPTH FOR TABLES [M].

REAL,   ALLOCATABLE, DIMENSION(:,:) :: TCGOND  !! SHALLOW WATER GROUP VELOCITY.
REAL,   ALLOCATABLE, DIMENSION(:,:) :: TFAK    !! SHALLOW WATER WAVE NUMBER.
REAL,   ALLOCATABLE, DIMENSION(:,:) :: TSIHKD  !! TABLE FOR OMEGA/SINH(2KD).
REAL,   ALLOCATABLE, DIMENSION(:,:) :: TFAC_ST !! TABLE FOR 2*G*K**2/
!                                                           (OMEGA*TANH(2KD)).
REAL,   ALLOCATABLE, DIMENSION(:,:) :: T_TAIL  !! TABLE FOR K**(-3)/
!                                                           GROUP VELOCITY.


PUBLIC :: NDEPTH, DEPTHA, DEPTHD, DEPTHE
PUBLIC :: TCGOND, TFAK, TSIHKD, TFAC_ST, T_TAIL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE PREPARE_TABLES                !! PRE-COMPUTES MODULE CONSTANTS.
   MODULE PROCEDURE PREPARE_TABLES
END INTERFACE
PUBLIC PREPARE_TABLES

INTERFACE PRINT_TABLES_STATUS            !! PRINTS THE STATUS OF THIS MODULE.
   MODULE PROCEDURE PRINT_TABLES_STATUS
END INTERFACE
PUBLIC PRINT_TABLES_STATUS

INTERFACE SET_DEPTH_TABLE_PARAMETER      !! DEFINE DEPTH TABLES.
   MODULE PROCEDURE SET_DEPTH_TABLE_PARAMETER
END INTERFACE
PUBLIC SET_DEPTH_TABLE_PARAMETER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE MIN_ENERGY              !! COMPUTATION OF MINMIMUM ENERGY TABLE.
   MODULE PROCEDURE MIN_ENERGY
END INTERFACE

INTERFACE MAKE_SHALLOW_TABLES     !! COMPUTE TABLES FOR SHALLOW WATER.
   MODULE PROCEDURE MAKE_SHALLOW_TABLES
END INTERFACE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_TABLES

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PREPARE_TABLES - ROUTINE TO PREPARE WAM TABLES MODULE.                     !
!                                                                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO COMPUTE ALL VARIABLES IN WAM_TABLES_MODULE.                         !
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
!     1. MINIMUM ENERGY TABLE.                                                 !
!        ---------------------                                                 !

CALL MIN_ENERGY
IF (ITEST.GE.2) WRITE (IU06,*) '   SUB.PREPARE_TABLES: MIN_ENERGY DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. SHALLOW WATER TABLES.                                                 !
!        ---------------------                                                 !

CALL MAKE_SHALLOW_TABLES
IF (ITEST.GT.1) WRITE (IU06,*) '   SUB.PREPARE_TABLES: MAKE_SHALLOW_TABLES DONE'

END SUBROUTINE PREPARE_TABLES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_TABLES_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PRINT_TABLES_STATUS - PRINT STATUS OF WAM_TABLES_MODULE.                 !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       MAKE A PRINTER OUTPUT OF THE DATA STORED IN WAM_TABLES_MODULE.         !
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
!     1. TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
!        --------------------------------------------

WRITE(IU06,'(/,'' -------------------------------------------------'')')
WRITE(IU06,*)'       WAM TABLES MODULE DATA' 
WRITE(IU06,'(  '' -------------------------------------------------'')')
WRITE(IU06,*) ' '
WRITE(IU06,*) ' DEPTH TABLES ARE PREPARED'
WRITE(IU06,*) ' DEPTH TABLE LENGTH ................: ', NDEPTH
WRITE(IU06,*) ' LOGARITHMIC DEPTH TABLE INCREMENT .: ', DEPTHD
WRITE(IU06,*) ' MINIMUM DEPTH VALUE ...............: ', DEPTHA
WRITE(IU06,*) ' MAXIMUM DEPTH VALUE ...............: ', DEPTHE

! ---------------------------------------------------------------------------- !
!
!    2. MINIMUM ENERGY TABLE.
!       ----------------------

WRITE(IU06,*) ' '
WRITE(IU06,*) ' MINIMUM ENERGY TABLE IS PREPARED'
WRITE(IU06,*) '  '
WRITE(IU06,*) ' WINDSPEED TABLE LENGTH ............: ', JUMAX
WRITE(IU06,*) ' WINDSPEED INCREMENT ...............: ', DELU
WRITE(IU06,*) ' MAXIMUM WIND SPEED ................: ', UMAX
WRITE(IU06,*) ' MAXIMUM FETCH FOR MINIMUM ENERGY ..: ', FETCH_MAX
WRITE(IU06,*) ' ABSOLUTE MINIMUM ENERGY ...........: ', FLMIN
WRITE(IU06,*) ' '

END SUBROUTINE PRINT_TABLES_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_DEPTH_TABLE_PARAMETER (N_DEPTH, DEPTH_S, DEPTH_I)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_DEPTH_TABLE_PARAMETER - DEFINES THE DEPTH TABLES.                      !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!     H.GUNTHER            GKSS       SEPTEMBER 2000   FT90                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO TRANSFER THE LENGTH, START DEPTH, AND DEPTH INCREMENT               !
!       OF THE SHALLOW WATER TABLES INTO THE WAM_TABLES_MODULE.                !
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

INTEGER,  INTENT(IN) :: N_DEPTH   !! LENGTH OF SHALLOW WATER TABLES.
REAL,     INTENT(IN) :: DEPTH_S   !! MINIMUM DEPTH FOR TABLES [M].
REAL,     INTENT(IN) :: DEPTH_I   !! DEPTH INCREMENT [%].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DEFINE VALUES FROM INPUT.                                             !
!        -------------------------                                             !

NDEPTH = N_DEPTH
DEPTHA = DEPTH_S
DEPTHD = DEPTH_I

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. DEFINE MISSING VALUES AND CHECK CONSISTENCY.                          !
!        --------------------------------------------                          !

IF (NDEPTH.LT.1 .OR. DEPTHA.LE.0. .OR. DEPTHD.LE.0.) THEN
   WRITE (IU06,*) ' **********************************************************'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *     FATAL ERROR IN SUB. SET_DEPTH_TABLE_PARAMETER      *'
   WRITE (IU06,*) ' *     =============================================      *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *     DEPTH TABLE COULD NOT BE DEFINED.                  *'
   WRITE (IU06,*) ' *     PRESENT DEFINITION VALUES ARE:                     *'
   WRITE (IU06,*) ' *       TABLE LENGTH    = ', NDEPTH
   WRITE (IU06,*) ' *       MINIMUM DEPTH   = ', DEPTHA
   WRITE (IU06,*) ' *       DEPTH INCREMENT = ', DEPTHD
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE (IU06,*) ' *                                                        *'
   WRITE (IU06,*) ' **********************************************************'
   CALL ABORT1
ELSE
   DEPTHE = DEPTHA*DEPTHD**(NDEPTH-1)
END IF

END SUBROUTINE SET_DEPTH_TABLE_PARAMETER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MIN_ENERGY

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MIN_ENERGY - COMPUTATION OF MINIMUM ENERGY IN SPECTRAL BINS.               !
!                                                                              !
!     J. BIDLOT    ECMWF                                                       !
!                                                                              !
!     PURPOSE.                                                                 !
!     ---------                                                                !
!                                                                              !
!       COMPUTATE OF A TABLE FOR MINIMUM ENERGY IN SPECTRAL BINS.              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       FOR EACH WINDSPEED U10 AS LISTED IN THE STRESS TABLE THE JONSWAP       !
!       SPECTRUM IS COMPUTED FROM FETCH LAWS AND 1% OF ITS VALUE IS STORED FOR !
!       EACH FREQUENY BIN. TO SMALL VALUES ARE REP†LAED BY A GIVEN MINIMUM.    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: GAMMA = 3.000000
REAL, PARAMETER :: SA=7.000000E-02
REAL, PARAMETER :: SB=9.000000E-02
INTEGER :: J
REAL :: FETCH
REAL, ALLOCATABLE :: U10(:), FPK(:), ALPHJO(:)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE WINDSPEEDS FROM TABLE PARAMETERS.                             !
!        -----------------------------------------                             !

ALLOCATE (U10(1:JUMAX))
DELU    = UMAX/REAL(JUMAX)
DO J = 1,JUMAX
   U10(J) = REAL(J)*DELU
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. PEAK FREQUENCIES AND ALPHA PARAMETER FROM FETCH LAW.                  !
!        ----------------------------------------------------                  !

ALLOCATE (FPK(1:JUMAX), ALPHJO(1:JUMAX))
FETCH = MIN(0.5*DELPHI,FETCH_MAX)
CALL FETCH_LAW (FETCH, FR(ML), U10, ALPHJO, FPK)
ALPHJO(:) = 0.01*ALPHJO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. JONSWAP SPECTRA.                                                      !
!        -----------------                                                     !

IF (ALLOCATED(FLMINFR))  DEALLOCATE(FLMINFR)
ALLOCATE(FLMINFR(1:JUMAX,1:ML))

CALL JONSWAP (FR, ALPHJO, GAMMA, SA, SB, FPK, FLMINFR(:,:))

FLMINFR(:,:) = MAX(FLMINFR(:,:),FLMIN)  !! AVOID TOO SMALL NUMBERS
DEALLOCATE (U10, FPK, ALPHJO)

END SUBROUTINE MIN_ENERGY

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MAKE_SHALLOW_TABLES

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MAKE_SHALLOW_TABLES - ROUTINE TO COMPUTE TABLES USED FOR SHALLOW WATER.    !
!                                                                              !
!     H.GUNTHER            ECMWF       04/04/1990                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE TABLES USED FOR SHALLOW WATER.                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!      TABLES FOR GROUP VELOCITY, WAVE NUMBER AND OMEGA/SINH(2KD) ARE COMPUTED !
!      AT ALL FREQUENCIES AND FOR A DEPTH TABLE OF LENGTH NDEPTH, STARTING AT  !
!      DEPTHA METERS AND INCREMENTED BY DEPTHD METRES.                         !
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

INTEGER        :: M, JD
REAL           :: GH, OM, AD, AK, AKD

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. ALLOCATE TABLE ARRAYS.                                                !
!        ----------------------                                                !

IF (ALLOCATED (TFAK))     DEALLOCATE (TFAK)
IF (ALLOCATED (TCGOND))   DEALLOCATE (TCGOND)
IF (ALLOCATED (TSIHKD))   DEALLOCATE (TSIHKD)
IF (ALLOCATED (TFAC_ST))  DEALLOCATE (TFAC_ST)
IF (ALLOCATED (T_TAIL))   DEALLOCATE (T_TAIL)
ALLOCATE (TFAK   (NDEPTH,ML))
ALLOCATE (TCGOND (NDEPTH,ML))
ALLOCATE (TSIHKD (NDEPTH,ML))
ALLOCATE (TFAC_ST(NDEPTH,ML))
ALLOCATE (T_TAIL (NDEPTH,ML))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. GROUP VELOCITY AND WAVE NUMBER.                                       !
!        -------------------------------                                       !

GH = G/(4.*PI)
DO M = 1,ML                             !! LOOP OVER FREQUENCIES.
   OM=ZPI*FR(M)
   DO JD = 1,NDEPTH                     !! LOOP OVER DEPTH.
      AD = DEPTHA*DEPTHD**(JD-1)
      AK = AKI(OM,AD)
      TFAK(JD,M) = AK
      AKD = AK*AD
      IF (AKD.LE.10.0) THEN
         TCGOND(JD,M) = 0.5*SQRT(G*TANH(AKD)/AK) * (1.0+2.0*AKD/SINH(2.*AKD))
         TSIHKD(JD,M) = OM/SINH(2.*AKD)
         TFAC_ST(JD,M) = 2.*G*AK**2/(OM*TANH(2.*AKD))
         T_TAIL(JD,M) = 1./(TFAK(JD,M)**3 * TCGOND(JD,M))
      ELSE
         TCGOND(JD,M) = GH/FR(M)
         TSIHKD(JD,M) = 0.
         TFAC_ST(JD,M) = 2./G*OM**3
         T_TAIL(JD,M) = 1./(TFAK(JD,M)**3 * TCGOND(JD,M))
      END IF
   END DO
END DO

END SUBROUTINE MAKE_SHALLOW_TABLES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_TABLES_MODULE

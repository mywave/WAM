MODULE WAM_SWELL_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS: ARRAYS NECESSARY FOR GRIDDED FIELDS OF PARAMTERS.    !
!                         ALL PROCEEDURES TO COMPUTE AND WRITE OUTPUT.         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_INTERFACE_MODULE, ONLY:  &
&       FEMEAN,                  &  !! COMPUTATION OF MEAN FREQUENCY.
&       MEAN_DIRECTION,          &  !! COMPUTATION OF MEAN DIRECTION AND SPREAD.
&       PEAK_PERIOD,             &  !! COMPUTATION OF PEAK PERIOD.
&       TM1_TM2_PERIODS,         &  !! COMPUTATION OF TM1 AND/OR TM2 PERIOD.
&       TOTAL_ENERGY                !! COMPUTATION OF TOTAL ENERGY.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE,ONLY: G, DEG, ZPI

USE WAM_FRE_DIR_MODULE,ONLY: KL, ML, FR, C, TH, DFIM, DFIMOFR, DFIM_FR,         &
&                            COSTH, SINTH
USE WAM_FILE_MODULE,   ONLY: IU06, ITEST

USE WAM_MODEL_MODULE,  ONLY: U10, UDIR, USTAR

USE WAM_OUTPUT_SET_UP_MODULE, ONLY: NOUT_P, CFLAG_P, orientation_of_directions

USE WAM_TABLES_MODULE,  ONLY: JUMAX, DELU, FLMINFR, FLMIN

use wam_mpi_module,    only: nijs, nijl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

include 'mpif.h'

PRIVATE
INTEGER, PARAMETER :: NTRAIN = 3         !! MAXIMUM NUMBER OF SWELL TRAINS.
REAL,    PARAMETER :: EPSMIN = 0.1E-32   !! SMALL NUMBER

REAL, ALLOCATABLE :: EMTRAIN(:,:)      !! TOTAL ENERGY OF SWELL TRAINS.
REAL, ALLOCATABLE :: THTRAIN(:,:)      !! MEAN DIRECTION OF SWELL TRAINS.
REAL, ALLOCATABLE :: PMTRAIN(:,:)      !! MEAN PERIOD (-1) OF SWELL TRAINS.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE SWELL_SEPARATION            !! SWELL SEPARATION AND INTEGRATED
   MODULE PROCEDURE SWELL_SEPARATION  !! PARAMETER OF SEA AND SWELL.
END INTERFACE
PUBLIC SWELL_SEPARATION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE PARTITIONING_SWELL            !! SWELL PARTITIONING
   MODULE PROCEDURE PARTITIONING_SWELL
END INTERFACE
PRIVATE PARTITIONING_SWELL

INTERFACE SMOOTHSARSPEC
   MODULE PROCEDURE SMOOTHSARSPEC
END INTERFACE
PRIVATE SMOOTHSARSPEC

INTERFACE FNDPRT
   MODULE PROCEDURE FNDPRT
END INTERFACE
PRIVATE FNDPRT

INTERFACE PARMEAN
   MODULE PROCEDURE PARMEAN
END INTERFACE
PRIVATE PARMEAN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SWELL_SEPARATION (FL3, FL1, BLOCK)

! ----------------------------------------------------------------------

!**** *SEPWISW* - COMPUTES THE SWELL ENERGY, THE MEAN SWELL DIRECTION,
!****             THE MEAN SWELL FREQUENCY AND THE MEAN WINDSEA DIR.

!     P.LIONELLO     FEBRUARY 87

!     L.ZAMBRESKY    NOVEMBER 87   GKSS/ECMWF   OPTIMIZED SUB.
!     J. BIDLOT      FEBRARY 1996       ECMWF   MESSAGE PASSING
!     J. BIDLOT      MAY     1999       ECMWF   ADD HIGH FREQUENCY TAIL
!     J. BIDLOT      APRIL   2000       ECMWF   ADD  EXTRA PARAMETERS
!     J. BIDLOT      DECEMBER2003       ECMWF   MOVE ALL ALLOCATION TO
!                                               *OUTBS*
!                                               TO MAKE IT CALLABLE IN AN
!                                               OPENMP LOOP.
!     L. AOUF        2013               MF      IMPLEMENTATION OF PARTITIONNING

!*    PURPOSE.
!     --------

!       TO SEPARATE THE SWELL FROM THE WIND INTERACTING SEA

!     METHOD.
!     -------

!       THE WAVES WHICH DO NOT INTERACT WITH THE WIND ARE
!       CONSIDERED SWELL.

!     EXTERNALS.
!     ----------

!       NONE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)  :: FL3(:,:,:)         !! BLOCK OF SPECTRA.
REAL,    INTENT(OUT) :: FL1(:,:,:)         !! BLOCK OF SWELL SPECTRA.
REAL,    INTENT(INOUT) :: BLOCK(NIJS:NIJL,1:NOUT_P) !! BLOCKED INTEGRATED PARAMETER
                                           !! FIRST INDEX COUNTS SEAPOINTS
                                           !! SECOND INDEX COUNTS PARAMETER.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: FRIC = 28.
INTEGER  :: K, M, IJ
REAL     :: CM
REAL,    DIMENSION(SIZE(FL3,1),KL,ML)  :: FL14
REAL,    DIMENSION(SIZE(FL3,1))    :: FLPP
REAL,    DIMENSION(SIZE(FL3,1))    :: EMWS
REAL,    DIMENSION(SIZE(FL3,1))    :: FMWS
REAL,    DIMENSION(SIZE(FL3,1))    :: FPWS
REAL,    DIMENSION(SIZE(FL3,1))    :: EMAXWS
REAL,    DIMENSION(SIZE(FL3,1))    :: DUM
REAL,    DIMENSION(SIZE(FL3,1))    :: FPMH
REAL,    DIMENSION(SIZE(FL3,1))    :: ETNP
INTEGER ,DIMENSION(SIZE(FL3,1))    :: MIJ

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. THE SWELL DISTRIBUTION IS COMPUTED.                                   !
!        -----------------------------------                                   !

EMWS(:) = EPSMIN
FMWS(:) = EPSMIN
FPWS(:) = FR(ML)
EMAXWS(:) = EPSMIN
DUM(:) = 0.

DO K=1,KL
   FLPP(:) = 1.2*USTAR(:)*COS(TH(K)-UDIR(:))
   DO M = 1,ML
      CM = FRIC/C(M)
      DO IJ=1,SIZE(FL3,1)
         IF (CM*FLPP(IJ).LT.1.0) then
            FL1(IJ,K,M) = MAX(FL3(IJ,K,M),EPSMIN)
         ELSE
            FL1(IJ,K,M) = 0.
            EMWS(IJ) = EMWS(IJ)+FL3(IJ,K,M)*DFIM(M)
            FMWS(IJ) = FMWS(IJ)+FL3(IJ,K,M)*DFIMOFR(M)
            IF (FL3(IJ,K,M).GT.EMAXWS(IJ)) THEN
               EMAXWS(IJ) = FL3(IJ,K,M)
               FPWS(IJ) = FR(M)
            END IF
         END IF
      END DO
   END DO
END DO

FMWS(:) = EMWS(:)/FMWS(:)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTATION INTEGRATED PARAMETER FOR TOTAL SWELL.                     !
!        -------------------------------------------------                     !

IF (CFLAG_P(25).OR.ANY(CFLAG_P(27:29))) THEN
   CALL TOTAL_ENERGY (FL1, BLOCK(NIJS:NIJL,25))
   IF (CFLAG_P(27)) THEN
      CALL FEMEAN (FL1,  BLOCK(NIJS:NIJL,25), FM=BLOCK(NIJS:NIJL,27))
      BLOCK(NIJS:NIJL,27) = 1./BLOCK(NIJS:NIJL,27)
   END IF
   IF (CFLAG_P(28).AND.CFLAG_P(29)) THEN
      CALL TM1_TM2_PERIODS (FL1, BLOCK(NIJS:NIJL,25), TM1=BLOCK(NIJS:NIJL,28), &
&                           TM2=BLOCK(NIJS:NIJL,29))
   ELSE IF (CFLAG_P(28)) THEN
      CALL TM1_TM2_PERIODS (FL1, BLOCK(NIJS:NIJL,25), TM1=BLOCK(NIJS:NIJL,28))
   ELSE IF (CFLAG_P(29)) THEN
      CALL TM1_TM2_PERIODS (FL1, BLOCK(NIJS:NIJL,25), TM2=BLOCK(NIJS:NIJL,29))
   END IF
END IF

IF (CFLAG_P(26)) CALL PEAK_PERIOD (FL1, BLOCK(NIJS:NIJL,26))

IF (CFLAG_P(30).AND.CFLAG_P(31)) THEN
   CALL MEAN_DIRECTION (FL1, THQ=BLOCK(NIJS:NIJL,30),                          &
&                            SPREAD=BLOCK(NIJS:NIJL,31))
ELSE IF (CFLAG_P(30)) THEN
   CALL MEAN_DIRECTION (FL1, THQ=BLOCK(NIJS:NIJL,30))
ELSE IF (CFLAG_P(31)) THEN
   CALL MEAN_DIRECTION (FL1, SPREAD=BLOCK(NIJS:NIJL,31))
END IF
if (orientation_of_directions) then
   IF (CFLAG_P(30)) BLOCK(NIJS:NIJL,30) = BLOCK(NIJS:NIJL,30)*DEG
else
   IF (CFLAG_P(30)) BLOCK(NIJS:NIJL,30) = mod(BLOCK(NIJS:NIJL,30)*DEG+180.,360.)
endif
IF (CFLAG_P(31)) BLOCK(NIJS:NIJL,31) = BLOCK(NIJS:NIJL,31)*DEG

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    4. PARTITION OF SWELL.                                                   !
!        -------------------

IF (ANY(CFLAG_P(41:49))) THEN
   IF(.NOT.ALLOCATED(EMTRAIN)) ALLOCATE(EMTRAIN(SIZE(FL3,1),NTRAIN))
   IF(.NOT.ALLOCATED(THTRAIN)) ALLOCATE(THTRAIN(SIZE(FL3,1),NTRAIN))
   IF(.NOT.ALLOCATED(PMTRAIN)) ALLOCATE(PMTRAIN(SIZE(FL3,1),NTRAIN))

   FPMH = 2.5/FR(1)
   MIJ(:) = INT(LOG10(MAX(FPWS(:),FMWS(:))*FPMH)*24.1589)+1
   MIJ(:) = MIN(MIJ(:),ML)
   CALL PARTITIONING_SWELL (FL1, MIJ, ETNP)


   IF (CFLAG_P(41)) BLOCK(NIJS:NIJL,41) = 4.*SQRT(EMTRAIN(:,1))
   IF (CFLAG_P(42)) BLOCK(NIJS:NIJL,42) = PMTRAIN(:,1)
   if (orientation_of_directions) then
      IF (CFLAG_P(43)) BLOCK(NIJS:NIJL,43) = THTRAIN(:,1)*DEG
   else
      IF (CFLAG_P(43)) BLOCK(NIJS:NIJL,43) = mod(THTRAIN(:,1)*DEG+180.,360.)
   endif
   IF (NTRAIN.GE.2) THEN
      IF (CFLAG_P(44)) BLOCK(NIJS:NIJL,44) = 4.*SQRT(EMTRAIN(:,2))
      IF (CFLAG_P(45)) BLOCK(NIJS:NIJL,45) = PMTRAIN(:,2)
      if (orientation_of_directions) then
         IF (CFLAG_P(46)) BLOCK(NIJS:NIJL,46) = THTRAIN(:,2)*DEG
      else
         IF (CFLAG_P(46)) BLOCK(NIJS:NIJL,46) = mod(THTRAIN(:,2)*DEG+180.,360.)
      endif
   END IF

   IF (NTRAIN.EQ.3) THEN
      IF (CFLAG_P(47)) BLOCK(NIJS:NIJL,47) = 4.*SQRT(EMTRAIN(:,3))
      IF (CFLAG_P(48)) BLOCK(NIJS:NIJL,48) = PMTRAIN(:,3)
      if (orientation_of_directions) then
         IF (CFLAG_P(49)) BLOCK(NIJS:NIJL,49) = THTRAIN(:,3)*DEG
      else
         IF (CFLAG_P(49)) BLOCK(NIJS:NIJL,49) = mod(THTRAIN(:,3)*DEG+180.,360.)
      endif
   END IF
   DEALLOCATE(EMTRAIN)
   DEALLOCATE(THTRAIN)
   DEALLOCATE(PMTRAIN)

ENDIF

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    6. COMPUTATION OF WIND SEA OUTPUT PARAMETERS
!        -----------------------------------------

FL14(:,:,:) = FL3(:,:,:) - FL1(:,:,:)    !! WINDSEA SPECTRA

IF (CFLAG_P(17).OR.ANY(CFLAG_P(19:21))) THEN
   CALL TOTAL_ENERGY (FL14, BLOCK(NIJS:NIJL,17))

   IF (CFLAG_P(19)) THEN
      CALL FEMEAN (FL14, BLOCK(NIJS:NIJL,17), FM=BLOCK(NIJS:NIJL,19))
      BLOCK(NIJS:NIJL,19) = 1./BLOCK(NIJS:NIJL,19)
   END IF
   IF (CFLAG_P(20).AND.CFLAG_P(21)) THEN
      CALL TM1_TM2_PERIODS (FL14, BLOCK(NIJS:NIJL,17), TM1=BLOCK(NIJS:NIJL,20),&
&                           TM2=BLOCK(NIJS:NIJL,21))
   ELSE IF (CFLAG_P(20)) THEN
      CALL TM1_TM2_PERIODS (FL14, BLOCK(NIJS:NIJL,17), TM1=BLOCK(NIJS:NIJL,20))
   ELSE IF (CFLAG_P(21)) THEN
      CALL TM1_TM2_PERIODS (FL14, BLOCK(NIJS:NIJL,17), TM2=BLOCK(NIJS:NIJL,21))
   END IF
END IF

IF (CFLAG_P(18)) CALL PEAK_PERIOD (FL14,BLOCK(NIJS:NIJL,18))

IF (CFLAG_P(22).AND.CFLAG_P(23)) THEN
   CALL MEAN_DIRECTION (FL14, THQ=BLOCK(NIJS:NIJL,22),                         &
&                       SPREAD=BLOCK(NIJS:NIJL,23))
ELSE IF (CFLAG_P(22)) THEN
   CALL MEAN_DIRECTION (FL14, THQ=BLOCK(NIJS:NIJL,22))
ELSE IF (CFLAG_P(23)) THEN
   CALL MEAN_DIRECTION (FL14, SPREAD=BLOCK(NIJS:NIJL,23))
END IF
if (orientation_of_directions) then
   IF (CFLAG_P(22)) BLOCK(NIJS:NIJL,22) = BLOCK(NIJS:NIJL,22)*DEG
else
   IF (CFLAG_P(22)) BLOCK(NIJS:NIJL,22) = mod(BLOCK(NIJS:NIJL,22)*DEG+180.,360.)
endif
IF (CFLAG_P(23)) BLOCK(NIJS:NIJL,23) = BLOCK(NIJS:NIJL,23)*DEG

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    7. RE-ASSIGN THE NON PARTITONED ENERGY TO WINDSEA                        !
!        AND REMOVE IT FROM THE SWELL.                                         !
!        -----------------------------------------------                       !

IF (ANY(CFLAG_P(41:49))) THEN
   IF (CFLAG_P(17)) BLOCK(NIJS:NIJL,17) = BLOCK(NIJS:NIJL,17)+ETNP(:)
   IF (CFLAG_P(25)) BLOCK(NIJS:NIJL,25) = MAX(BLOCK(NIJS:NIJL,25)-ETNP(:),0.0)
ENDIF

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    7. COMPUTE SIG. WAVEHEIGHT.                                              !
!        ------------------------                                              !

IF (CFLAG_P(17)) BLOCK(NIJS:NIJL,17) = 4.*SQRT(BLOCK(NIJS:NIJL,17))
IF (CFLAG_P(25)) BLOCK(NIJS:NIJL,25) = 4.*SQRT(BLOCK(NIJS:NIJL,25))

END SUBROUTINE SWELL_SEPARATION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PARTITIONING_SWELL (FL1, MIJ, ETNP)

! ----------------------------------------------------------------------

!**** *SEP3TR* -  COMPUTE ENERGY, DIRECTION AND PERIOD FOR THE 3 WAVE
!****             TRAINS: SWELL1, SWELL2, SWELL3

!     D.PETTENUZZO     MAY 2012
!     JEAN BIDLOT      JUNE 2013  SIMPLIFIED BY REMOVING THE WIND SEA PART
!                                 AND ONLY PARTITIONED THE SWELL SPECTRUM
!                                 AND ONLY LOOKING IN THE PROGNOSTIC RANGE


!*    PURPOSE.
!     --------

!     CREATED TO TEST WAVE SWELL SPECTRA PARTITIONING INTO 3 WAVE TRAINS 

!     METHOD.
!     -------

!       HANSON AND PHILLIPS 2001

!     EXTERNALS.
!     ----------

!       *FNDPRT*  - COMPUTE PARTITION MASKS
!       *PARMEAN* - COMPUTE MEAN PARAMETERS.

! ---------------------------------------------------------------------------- !
!
!     INTERFACE:
!     ----------

REAL,INTENT(IN),  DIMENSION(:,:,:) :: FL1      !! 2D SWELL SPECTRA.
INTEGER,INTENT(IN),  DIMENSION(:) :: MIJ  !! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
REAL,INTENT(OUT), DIMENSION(:) :: ETNP  !!   TOTAL ENERGY NOT HELD BY ALL PARTITIONS.
!                     THIS IS POSSIBLE BECAUSE WE SUPPLY SWELL SPECTRA
!                     IN WHICH THE WIND SEA HAS BEEN REMOVED. SOME OF
!                     THE ENERGY LEAKING OUT FROM THE WINDSEA SPECTRUM
!                     MIGHT NOT BE PARTITIONED.


! ---------------------------------------------------------------------------- !
!
!     LOCAL VARIABLES.
!     ----------------
!

! ----------------------------------------------------------------------

INTEGER, PARAMETER :: NPMAX=20
INTEGER :: IJ, M, K, IP
INTEGER :: ISORT
INTEGER :: IFL, IFH, ITHL, ITHH
INTEGER, DIMENSION(SIZE(FL1,1)) :: NPEAK
INTEGER, DIMENSION(SIZE(FL1,1)) :: MMIN, MMAX
INTEGER, DIMENSION(SIZE(FL1,1)) :: IPNOW
INTEGER, DIMENSION(SIZE(FL1,1)) :: JU
INTEGER, DIMENSION(SIZE(FL1,1),NTRAIN) :: IENERGY
INTEGER, DIMENSION(SIZE(FL1,1),NPMAX) :: NFRP, NTHP

REAL :: COSDIFF
REAL :: XJ, DELUINV
REAL :: FLLOWEST

REAL, DIMENSION(SIZE(FL1,1)) :: ENEX, SUMETRAIN
REAL, DIMENSION(SIZE(FL1,1)) :: ETOT
REAL, DIMENSION(SIZE(FL1,1)) :: ETT
REAL, DIMENSION(SIZE(FL1,1)) :: ENMAX
REAL, DIMENSION(KL,ML) :: W2
REAL, DIMENSION(SIZE(FL1,1),KL) :: SPRD
REAL, DIMENSION(SIZE(FL1,1),0:NPMAX) :: DIR, PER, ENE
REAL, DIMENSION(KL,ML,NPMAX) :: SPEC
REAL, DIMENSION(SIZE(FL1,1),KL,ML) :: FL3
REAL, DIMENSION(SIZE(FL1,1),KL,ML) :: W1
REAL, DIMENSION(SIZE(FL1,1),KL,ML) :: FLLOW

LOGICAL :: LLADDPART
LOGICAL, DIMENSION(SIZE(FL1,1),KL) :: LLCOSDIFF
LOGICAL, DIMENSION(KL,ML) :: LLW3
LOGICAL :: First = .true.
! ----------------------------------------------------------------------

!*    LOOP THROUGH THE GRID POINTS
!     ----------------------------

CALL TOTAL_ENERGY (FL1, ETOT)   !! TOTAL ENERGY IN THE SWELL SPECTRA

FL3(:,:,:) = FL1(:,:,:)
CALL SMOOTHSARSPEC (FL3)  !! SMOOTH INPUT SPECTRA
WHERE (FL1(:,:,:).LE.0.) FL3(:,:,:) = 0.   !! RE-IMPOSE THE WINDSEA MASK


EMTRAIN(:,:) = 0.
DO ISORT=1,NTRAIN
   THTRAIN(:,ISORT) = UDIR(:)
END DO
PMTRAIN(:,:) = 0.
NPEAK(:) = 0

!*    1. DETERMINATES MAXIMA POSITIONS AND NUMBER
!        ----------------------------------------

DO K = 1,KL
   DO IJ = 1,SIZE(FL1,1)
      COSDIFF=COS(TH(K)-UDIR(NIJS-1+IJ))
      LLCOSDIFF(IJ,K)=(COSDIFF.LT.-0.4)
      SPRD(IJ,K)=MAX(0.,COSDIFF)**2
   END DO
END DO

DELUINV=1./DELU
DO IJ = 1,SIZE(FL1,1)
   XJ = U10(NIJS-1+IJ)*DELUINV
   JU(IJ) =MIN(JUMAX, MAX(NINT(XJ),1))
END DO

DO M = 1,ML
   DO K = 1,KL
      DO IJ = 1,SIZE(FL1,1)
         FLLOW(IJ,K,M) = MAX(FLMINFR(JU(IJ),M)*SPRD(IJ,K),FLMIN)
      END DO
   END DO
END DO


DO IJ = 1,SIZE(FL1,1)
   OUT: DO M= 1,MIJ(IJ)-1
      IFL    = MAX ( 1 , M-1 )
      IFH    = MIN ( ML , M+1 )
      DO K=1,KL

         FLLOWEST = FLLOW(IJ,K,M)

         IF(FL3(IJ,K,M) .GT. FLLOWEST) THEN
            ITHL   = 1 + MOD(KL+K-2,KL)
            ITHH   = 1 + MOD(K,KL)
            IF ( FL3(IJ,ITHL,M   ) .GT. 0.0 .AND.                              &
&                FL3(IJ,ITHH,M   ) .GT. 0.0 .AND.                              &
&                FL3(IJ,K   ,IFL ) .GT. 0.0 .AND.                              &
&                FL3(IJ,K   ,IFH ) .GT. 0.0 .AND.                              &
&                FL3(IJ,ITHL,IFL ) .GT. 0.0 .AND.                              &
&                FL3(IJ,ITHL,IFH ) .GT. 0.0 .AND.                              &
&                FL3(IJ,ITHH,IFL ) .GT. 0.0 .AND.                              &
&                FL3(IJ,ITHH,IFH ) .GT. 0.0 ) THEN


               IF ( FL3(IJ,K,M) .GE. FL3(IJ,K   ,IFL ) .AND.                   &
&                   FL3(IJ,K,M) .GE. FL3(IJ,K   ,IFH ) .AND.                   &
&                   FL3(IJ,K,M) .GE. FL3(IJ,ITHL,IFL ) .AND.                   &
&                   FL3(IJ,K,M) .GE. FL3(IJ,ITHL,M   ) .AND.                   &
&                   FL3(IJ,K,M) .GE. FL3(IJ,ITHL,IFH ) .AND.                   &
&                   FL3(IJ,K,M) .GE. FL3(IJ,ITHH,IFL ) .AND.                   &
&                   FL3(IJ,K,M) .GE. FL3(IJ,ITHH,M   ) .AND.                   &
&                   FL3(IJ,K,M) .GE. FL3(IJ,ITHH,IFH ) ) THEN
                    NPEAK(IJ) = NPEAK(IJ) + 1
                    IF (NPEAK(IJ).GT.NPMAX) EXIT OUT
                    NFRP(IJ,NPEAK(IJ)) = M
                    NTHP(IJ,NPEAK(IJ)) = K
               END IF
            END IF
         END IF
      END DO
   END DO OUT
END DO

NPEAK(:) = MIN(NPEAK(:),NPMAX)

!*    2. GENERATE MASK FOR EACH PARTITION AND COMPUTE STATISTICS
!        -------------------------------------------------------

!     ENE must be initialised and index 0 must stay equal to 0
ENE(:,:)=0.
DIR(:,:)=0.
PER(:,:)=0.

!     POINTS BELOW THE MINIMUM LEVEL BELONG TO THE WINDSEA PART
!     AND CAN BE EXCLUDED FROM THE PARTITIONS OF THE SWELL SPECTRUM

MMIN(:) = ML
MMAX(:) = 0


DO M=1,ML
   DO K=1,KL
      DO IJ = 1,SIZE(FL1,1)
         IF(FL3(IJ,K,M).LE.FLLOW(IJ,K,M)) THEN
            W1(IJ,K,M) = 1.
         ELSE
            W1(IJ,K,M) = 0.
            MMIN(IJ) = MIN(M,MMIN(IJ))
            MMAX(IJ) = M
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO IJ = 1,SIZE(FL1,1)

   DO M=1,ML
      DO K=1,KL
         LLW3(K,M) = (W1(IJ,K,M).LT.1.0)
      ENDDO
   ENDDO

   DO IP=1,NPEAK(IJ)
      CALL FNDPRT(NTHP(IJ,IP),NFRP(IJ,IP),                   &
&                 MIJ(IJ),MMIN(IJ),MMAX(IJ),                 &
&                 FL3(IJ,1:KL,1:ML),LLW3,                    &
&                 W1(IJ,1:KL,1:ML),W2)
  
      DO M=1,ML
         DO K=1,KL
            SPEC(K,M,IP) = FL3(IJ,K,M)*W2(K,M)
         ENDDO
      ENDDO

   ENDDO

!       CHECK THAT UNASSIGNED BINS ARE IN THE WIND SECTOR, OTHERWISE THEY SHOULD BE
!       ADDED AS A PARTITION

   LLADDPART = .FALSE.
   IF (NPEAK(IJ).LT.NPMAX) THEN
      DO K = 1,KL
         DO M = 1,ML
            IF (LLCOSDIFF(IJ,K) .AND. W1(IJ,K,M).LE.0.0) THEN
                W2(K,M) = 1.
                LLADDPART=.TRUE.
            ELSE
                W2(K,M) = 0.
            ENDIF
         ENDDO
      ENDDO
   ENDIF

   IF (LLADDPART) THEN
      NPEAK(IJ)=NPEAK(IJ)+1
      IP=NPEAK(IJ)
      DO M=1,ML
         DO K=1,KL
            SPEC(K,M,IP) = FL3(IJ,K,M)*W2(K,M)
         ENDDO
      ENDDO
   ENDIF

   IF (NPEAK(IJ).GT.0) THEN
      CALL PARMEAN(SPEC(1:KL,1:ML,1:NPEAK(IJ)),                         &
&                  NPEAK(IJ), ENE(IJ,1:NPEAK(IJ)),                      &
&                  DIR(IJ,1:NPEAK(IJ)),PER(IJ,1:NPEAK(IJ)))
   ENDIF
ENDDO


!*    5. SORT PARTITIONS ACCORDING TO ENERGY AND ASSIGN THE FIRST
!        NTRAIN TRAINS WE ONLY NEED THE FIRST NTRAIN TO BE SORTED.
!        ---------------------------------------------------------

ETT(:)=ENE(:,1)
DO IP = 2,NPMAX
   ETT(:)=ETT(:)+ENE(:,IP)
ENDDO

DO ISORT = 1,NTRAIN
   IPNOW(:)=0
   ENMAX(:)=0.

   DO IP=1,NPMAX
      DO IJ = 1,SIZE(FL1,1)
         IF (ENE(IJ,IP).GT.ENMAX(IJ)) THEN
            IPNOW(IJ) = IP
            ENMAX(IJ) = ENE(IJ,IP)
         ENDIF
      ENDDO
   ENDDO

   DO IJ = 1,SIZE(FL1,1)
      EMTRAIN(IJ,ISORT)=ENE(IJ,IPNOW(IJ))
      THTRAIN(IJ,ISORT)=DIR(IJ,IPNOW(IJ))
      PMTRAIN(IJ,ISORT)=PER(IJ,IPNOW(IJ))
      ENE(IJ,IPNOW(IJ))=0.
      IENERGY(IJ,ISORT)=MIN(IPNOW(IJ),1)
   ENDDO
ENDDO


!*    6. PRESERVE TOTAL ENERGY
!        ---------------------

!     6.1 Distribute extra energy proportionally to swell trains

SUMETRAIN(:) = MAX(EMTRAIN(:,1),EPSMIN)
DO ISORT=2,NTRAIN
   SUMETRAIN(:) = SUMETRAIN(:)+EMTRAIN(:,ISORT)
ENDDO

ENEX(:) = MAX((ETT(:)-SUMETRAIN(:)),0.)/SUMETRAIN(:)

DO ISORT = 1,NTRAIN
    EMTRAIN(:,ISORT) = EMTRAIN(:,ISORT)+ENEX(:)*EMTRAIN(:,ISORT)
ENDDO

!     PREPARE OUTPUT

DO ISORT=1,NTRAIN
   DO IJ = 1,SIZE(FL1,1)
      IF(IENERGY(IJ,ISORT).EQ.0) THEN
         EMTRAIN(IJ,ISORT) = 0.0
         THTRAIN(IJ,ISORT) = UDIR(NIJS-1+IJ)
         PMTRAIN(IJ,ISORT) = 0.0
      ENDIF
   ENDDO
ENDDO

!     NON PARTITIONED ENERGY

ETNP(:) = MAX(ETOT(:)-ETT(:),0.)

END SUBROUTINE PARTITIONING_SWELL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SMOOTHSARSPEC (SPEC)
 
! ---------------------------------------------------------------------------- !

!     PURPOSE
!     -------
!     TO SMOOTH SPECTRA THAT ARE USED IN THE PARTITIONING.

!     AUTHOR
!     ------
!     J. BIDLOT  ECMWF APRIL 2001.

!     EXTERNALS
!     ---------

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE:
!     ----------

REAL, INTENT(INOUT) :: SPEC (:,:,:)   !! BLOCK OF SPECTRA

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.
!     ----------------

INTEGER :: K, M, KM , KP

REAL :: C16
REAL :: WORK(SIZE(SPEC,1),KL,ML)
REAL, DIMENSION(SIZE(SPEC,1)) ::  EM, EW, TEMP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE THE ENERGY OF THE ORIGINAL SPECTRA.                           !
!        -------------------------------------------                           !

EM(:) = EPSMIN

DO M=1,ML
   TEMP(:) = SPEC(:,1,M)
   DO K = 2,KL
      TEMP(:) = TEMP(:)+SPEC(:,K,M)
   ENDDO
   EM(:) = EM(:)+DFIM(M)*TEMP(:)
ENDDO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. SMOOTH SPECTRA.                                                       !
!        ----------------                                                      !

DO K=1,KL
   KM  = K - 1
   IF (KM .LT. 1) KM = KL
   KP = K + 1
   IF (KP.GT. KL) KP = 1

   DO M=2,ML-1
      WORK(:,K,M) = SPEC(:,KM,M-1) + 2.*SPEC(:,KM,M) +    SPEC(:,KM,M+1)      &
&              + 2.*SPEC(:,K ,M-1) + 4.*SPEC(:,K ,M) + 2.*SPEC(:,K ,M+1)      &
&              +    SPEC(:,KP,M-1) + 2.*SPEC(:,KP,M) +    SPEC(:,KP,M+1)
   ENDDO
ENDDO

C16 = 1./16

DO M = 2,ML-1
   SPEC(:,:,M) = C16 * WORK(:,:,M)
ENDDO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTE THE ENERGY OF THE SMOOTHED SPECTRA.                           !
!        -------------------------------------------                           !

EW(:) = EPSMIN

DO M = 1,ML
   TEMP(:) = SPEC(:,1,M)
   DO K = 2,KL
      TEMP(:) = TEMP(:)+SPEC(:,K,M)
   ENDDO
   EW(:) = EW(:)+DFIM(M)*TEMP(:)
ENDDO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. ENSURE THE ENERGY IS CONSERVED BETWEEN ORIGINAL AND SMOOTHED SPECTRA. !
!        --------------------------------------------------------------------- !

DO M = 1,ML
   DO K = 1,KL
      SPEC(:,K,M) = (EM(:)/EW(:))*SPEC(:,K,M)
   ENDDO
ENDDO

END SUBROUTINE SMOOTHSARSPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FNDPRT ( ITHC, IFRC, MIJ, MMIN, MMAX, SPEC, LLW3, W1, W2 )

! ---------------------------------------------------------------------------- !

!**** *FNDPRT* -  CALCULATE PARTITION MASKS

!     D.PETTENUZZO     MAY 2011
!     MODIFIED MAY 2012 TO SUIT ECMWF CODES


!*    PURPOSE.
!     --------

!     FIND ALL THE POINTS IN THE FREQUENCY DIRECTION DOMAIN
!     WHICH BELONG TO A GIVEN PARTITION IDENTIFIED BY ITS
!     PEAK

!     METHOD.
!     -------

!       STEEPEST ASCENT

!     EXTERNALS.
!     ----------

!       NONE.

! ---------------------------------------------------------------------------- !
!
!     INTERFACE:
!     ----------

INTEGER,INTENT(IN) :: ITHC    !! PEAK DISCRETE DIRECTION
INTEGER,INTENT(IN) :: IFRC    !! PEAK DISCRETE FREQUENCY
INTEGER,INTENT(IN) :: MIJ     !! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
INTEGER,INTENT(IN) :: MMIN    !! FREQUENCY INDEX BELOW WHICH SPECTRAL BINS ARE
                              !! ALREADY EXCLUDED
INTEGER,INTENT(INOUT) :: MMAX !! LAST FREQUENCY INDEX OF WHERE SPECTRAL ARE
                              !! ALREADY EXCLUDED
REAL,INTENT(IN),       DIMENSION(KL,ML) :: SPEC  !! SPECTRUM
LOGICAL,INTENT(INOUT), DIMENSION(KL,ML) :: LLW3  !! MAP OF BINS THAT ARE 
                                                 !! ALWAYS EXCLUDED (WINDSEA)
REAL,INTENT(INOUT),  DIMENSION(KL,ML) :: W1   !! OVERALL MAP OF BINS
REAL,INTENT(OUT),  DIMENSION(KL,ML) :: W2   !! MAP OF BINS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLE                                                           !

INTEGER :: M, K, NITT

INTEGER :: KLH, KK, KKMIN, KKMAX
INTEGER :: IFRL, ITHL, ITHR, IFRH
INTEGER, DIMENSION(1-KL:2*KL) :: KLOC
REAL :: HALF_SECTOR
LOGICAL :: LLCHANGE, LLADD


! ----------------------------------------------------------------------

!     ASSUME THAT ONE HAS ONLY SEARCH A SECTOR IN DIRECTION (not the full 360 degrees)
!     AROUND THE PEAK
!!!   THIS IS A STRONG ASSUMPTION ON THE DIRECTIONALITY OF SWELL !!!

HALF_SECTOR = 70.
KLH = NINT((HALF_SECTOR/360.)*KL)+1
KKMIN = ITHC-KLH
KKMAX = ITHC+KLH
DO KK = KKMIN,KKMAX
   KLOC(KK)=1+MOD(KL+KK-1,KL)
ENDDO

!*    1.  SET UP THE W2 MAP
!         -----------------

W2(:,:) = 0.

IFRL   = MAX ( 1 , IFRC-1 )
IFRH   = MIN ( ML , IFRC+1 )

DO K=ITHC-1, ITHC+1
   ITHL = 1 + MOD(KL+K-1,KL)
   DO M=IFRC-1, IFRC+1
      IFRL = MAX(1,MIN(ML,M))
      IF ( W1(ITHL,IFRL) .LE. 0.5 ) W2(ITHL,IFRL) = 0.5
   ENDDO
ENDDO

IF ( W1(ITHC,IFRC) .LT. 0.25 ) W2(ITHC,IFRC) = 1.0

!     FIND IF MORE HIGH FREQUENCY BINS HAVE BECOME EXCLUDED

OUT0: DO M = MMAX,MMIN,-1
   DO KK = KKMIN,KKMAX
      K=KLOC(KK)
      IF ( W1(K,M).LT.1.0 ) THEN
         MMAX=M
         EXIT OUT0
      ENDIF
   ENDDO
ENDDO OUT0

!*    2.  ITTERATE SEARCH
!         ---------------

NITT = 0

!     2.a Branch point

  200 CONTINUE

   NITT=NITT+1
   LLCHANGE=.FALSE.

!     2.b Determine central points

   DO M=MMIN,MIN(MIJ,MMAX)
!       by definition bins beyond M=MIJ are never extremas
!       and bins above MMAX are excluded.
      DO KK=KKMIN,KKMAX
         K=KLOC(KK)

         IF ( LLW3(K,M) ) THEN
            IF ( W2(K,M).EQ.0.5 .AND. W1(K,M).LT.0.5 ) THEN
               LLADD    = .TRUE.
               OUT1: DO ITHR=K-1, K+1
                  ITHL = 1 + MOD(KL+ITHR-1,KL)
                  DO IFRL=MAX(1,M-1), MIN(ML,M+1)
                     IF ( W2 (ITHL,IFRL).EQ.0. .AND.      &
&                       SPEC(ITHL,IFRL).GT.SPEC(K,M) ) THEN
                        LLADD = .FALSE.
                        EXIT OUT1
                     ENDIF
                  ENDDO
               ENDDO OUT1
               IF ( LLADD ) THEN
                  W2(K,M) = 1.
                  LLCHANGE = .TRUE.
               ENDIF
            ENDIF
         ENDIF

      ENDDO
   ENDDO

!     2.c Determine peripherical points

   DO M=MMIN,MMAX
      DO KK=KKMIN,KKMAX
         K=KLOC(KK)

         IF ( LLW3(K,M) .AND. W1(K,M).LT.1.0 ) THEN
            IF ( W2(K,M).EQ.0. ) THEN
               OUT2: DO ITHR=K-1, K+1
                  ITHL = 1 + MOD(KL+ITHR-1,KL)
                  DO IFRL=MAX(1,M-1), MIN(ML,M+1)
                     IF ( W2(ITHL,IFRL).EQ.1. ) THEN
                        W2(K,M) = 0.5
                        LLCHANGE = .TRUE.
                        EXIT OUT2
                     ENDIF
                  ENDDO
               ENDDO OUT2
            ENDIF
         ENDIF

      ENDDO
   ENDDO

!*    2.d Branch back ?

IF ( LLCHANGE .AND. NITT.LT.25 ) GOTO 200


!*    3   UPDATE THE OVERALL MAP

DO M = 1, ML
   DO K = 1, KL
      W1(K,M) = W1(K,M) + W2(K,M)
   ENDDO
ENDDO

END SUBROUTINE FNDPRT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PARMEAN (F, NPK, EM, THQ, MEANWP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!**** *PARMEAN* - COMPUTATION OF TOTAL ENERGY, MEAN DIRECTION,
!                  AND MEAN PERIOD BASED ON THE FIRST MOMENT
!
!      J-R BIDLOT    ECMWF     MARCH 2000
!      D PETTENUZZO  MAY 2012 MERGED 3 ROUTINES FOR PAR COMPUTATION
!                                                                              !
!*    PURPOSE.
!     --------
!
!       COMPUTE TOTAL ENERGY, MEAN DIRECTION AND MEAN PERIOD.
!                                                                              !
!     METHOD.
!     -------
!                                                                              !
!       NONE.
!                                                                              !
!     EXTERNALS.
!     ----------
!                                                                              !
!       NONE.
!                                                                              !
!     REFERENCE.
!     ----------
!                                                                              !
!       NONE.
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLE                                                       !

INTEGER, INTENT(IN)  :: NPK      !! NUMBER OF PEAKS
REAL,INTENT(IN),   DIMENSION(KL,ML,NPK) :: F  !! SPECTRA.
REAL,INTENT(OUT),  DIMENSION(NPK) :: EM       !! MEAN WAVE ENERGY
REAL,INTENT(OUT),  DIMENSION(NPK) :: THQ      !! MEAN WAVE DIRECTION
REAL,INTENT(OUT),  DIMENSION(NPK) :: MEANWP   !! MEAN PERIOD BASED ON 1 MOMENT.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLE                                                           !

INTEGER :: IPK, K, M
REAL    :: F1D(ML)
REAL    :: TEMP, SI, CI

! ---------------------------------------------------------------------------- !

DO IPK = 1,NPK

   EM(IPK) = EPSMIN
   MEANWP(IPK) = EPSMIN
   DO M=1,ML
      F1D(M) = F(1,M,IPK)
      DO K=2,KL
         F1D(M) = F1D(M)+F(K,M,IPK)
      ENDDO
      EM(IPK) = EM(IPK)+F1D(M)*DFIM(M)
      MEANWP(IPK) = MEANWP(IPK)+F1D(M)*DFIM_FR(M)
   ENDDO

   IF(EM(IPK).GT.EPSMIN) THEN
      MEANWP(IPK) = EM(IPK)/MEANWP(IPK)
   ELSE
      MEANWP(IPK) = 0.
   ENDIF

ENDDO

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    MEAN WAVE DIRECTION
!     -------------------

DO IPK=1,NPK
   SI = 0.
   CI = 0.

   DO K=1,KL
      TEMP = F(K,1,IPK)*DFIM(1)
      DO M=2,ML
         TEMP = TEMP + F(K,M,IPK)*DFIM(M)
      ENDDO
      SI = SI+SINTH(K)*TEMP
      CI = CI+COSTH(K)*TEMP
   ENDDO

   IF (CI.EQ.0.) CI = EPSMIN
   THQ(IPK) = ATAN2(SI,CI)
   IF (THQ(IPK).LT.0.) THQ(IPK) = THQ(IPK) + ZPI
ENDDO

END SUBROUTINE PARMEAN

! ---------------------------------------------------------------------------- !

END MODULE WAM_SWELL_MODULE

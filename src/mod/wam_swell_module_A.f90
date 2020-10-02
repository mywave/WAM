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

USE WAM_GENERAL_MODULE,ONLY: G, DEG

USE WAM_FRE_DIR_MODULE,ONLY: KL, ML, FR, C, TH

USE WAM_FILE_MODULE,   ONLY: IU06, ITEST

USE WAM_MODEL_MODULE,  ONLY: UDIR, USTAR

USE WAM_OUTPUT_SET_UP_MODULE, ONLY: NOUT_P, CFLAG_P, orientation_of_directions

use wam_mpi_module,    only: nijs, nijl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
include 'mpif.h'
PRIVATE

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

!**   INTERFACE.
!     ----------

!       *CALL* *SEPWISW (FL3, FL1, IJS, IJL, THWOLD,USOLD,USTAR)*
!          *FL3* - BLOCK OF SPECTRA
!          *FL1* - SWELL FLAG ARRAY
!                  = 1 IF COMPONENT IS SWELL
!                  = 0 IF COMPONENT IS SEA
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                      WIND VELOCITY.
!          *USTAR*     NEW FRICTION VELOCITY IN M/S.
!          *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
!                      VELOCITY.

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

REAL,    INTENT(IN)  :: FL3(:,:,:)        !! BLOCK OF SPECTRA.
REAL,    INTENT(OUT) :: FL1(:,:,:)        !! BLOCK OF SWELL SPECTRA.
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
REAL,DIMENSION(SIZE(FL3,1),KL,ML)  :: FL11, FL13, FL14
REAL,DIMENSION(SIZE(FL3,1))    :: FLPP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. THE SWELL DISTRIBUTION IS COMPUTED.                                   !
!        -----------------------------------                                   !

FL1 = 0.
DO K=1,KL
   FLPP(:) = 1.2*USTAR(:)*COS(TH(K)-UDIR(:))
   DO M = 1,ML
      CM = FRIC/C(M)
      WHERE (CM*FLPP(:)<1.0) FL1(:,k,m) = FL3(:,k,m)
   END DO
END DO

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
   IF (CFLAG_P(25)) BLOCK(NIJS:NIJL,25) = 4.*SQRT( BLOCK(NIJS:NIJL,25))
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

FL11(:,:,:) = 0.
FL13(:,:,:) = 0.


CALL PARTITIONING_SWELL (FL1, FL11)

FL14(:,:,:) = FL1(:,:,:) - FL11(:,:,:)

CALL PARTITIONING_SWELL (FL14, FL13)

WHERE (FL13(:,:,:).LT.0.) FL13(:,:,:) = 0.

! classement des houles primaire et secondaire

CALL TOTAL_ENERGY (FL11, BLOCK(NIJS:NIJL,41))
CALL TOTAL_ENERGY (FL13, BLOCK(NIJS:NIJL,44))

DO M=1,ML
   DO K=1,KL
      WHERE (BLOCK(NIJS:NIJL,41).LT.BLOCK(NIJS:NIJL,44))
         FLPP(:) = FL11(:,K,M)
         FL11(:,K,M) = FL13(:,K,M)
         FL13(:,K,M) = FLPP(:)
      ENDWHERE
   ENDDO
ENDDO
DO IJ = NIJS,NIJL
   IF (BLOCK(IJ,41).LT.BLOCK(IJ,44)) THEN
      CM = BLOCK(IJ,41)
      BLOCK(IJ,41) = BLOCK(IJ,44)
      BLOCK(IJ,44) = CM
   ENDIF
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. COMPUTATION INTEGRATED PARAMETER FOR FIRST SWELL.                     !
!        -------------------------------------------------                     !

IF (CFLAG_P(41).OR.CFLAG_P(42)) THEN
   IF (CFLAG_P(42)) THEN
      CALL TM1_TM2_PERIODS (FL11, BLOCK(NIJS:NIJL,41), TM1=BLOCK(NIJS:NIJL,42))
   END IF
   IF (CFLAG_P(41)) BLOCK(NIJS:NIJL,41) = 4.*SQRT( BLOCK(NIJS:NIJL,41))
END IF

IF (CFLAG_P(43)) THEN
   CALL MEAN_DIRECTION (FL11, THQ=BLOCK(NIJS:NIJL,43))
   if (orientation_of_directions) then
      BLOCK(NIJS:NIJL,43) = BLOCK(NIJS:NIJL,43)*DEG
   else
      BLOCK(NIJS:NIJL,43) = mod(BLOCK(NIJS:NIJL,43)*DEG+180.,360.)
   endif
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. COMPUTATION INTEGRATED PARAMETER FOR SECOND SWELL.                    !
!        --------------------------------------------------                    !

IF (CFLAG_P(44).OR.CFLAG_P(45)) THEN
   IF (CFLAG_P(45)) THEN
      CALL TM1_TM2_PERIODS (FL13, BLOCK(NIJS:NIJL,44), TM1=BLOCK(NIJS:NIJL,45))
   END IF
   IF (CFLAG_P(44)) BLOCK(NIJS:NIJL,44) = 4.*SQRT( BLOCK(NIJS:NIJL,44))
END IF

IF (CFLAG_P(46)) THEN
   CALL MEAN_DIRECTION (FL13, THQ=BLOCK(NIJS:NIJL,46))
   if (orientation_of_directions) then
      BLOCK(NIJS:NIJL,46) = BLOCK(NIJS:NIJL,46)*DEG
   else
      BLOCK(NIJS:NIJL,46) = mod(BLOCK(NIJS:NIJL,46)*DEG+180.,360.)
   endif
END IF

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
      CALL TM1_TM2_PERIODS (FL14, BLOCK(NIJS:NIJL,17), TM1=BLOCK(NIJS:NIJL,20), &
&                           TM2=BLOCK(NIJS:NIJL,21))
   ELSE IF (CFLAG_P(20)) THEN
      CALL TM1_TM2_PERIODS (FL14, BLOCK(NIJS:NIJL,17), TM1=BLOCK(NIJS:NIJL,20))
   ELSE IF (CFLAG_P(21)) THEN
      CALL TM1_TM2_PERIODS (FL14, BLOCK(NIJS:NIJL,17), TM2=BLOCK(NIJS:NIJL,21))
   END IF
   IF (CFLAG_P(17)) BLOCK(NIJS:NIJL,17) = 4.*SQRT(BLOCK(NIJS:NIJL,17))
END IF

IF (CFLAG_P(18)) CALL PEAK_PERIOD (FL14,BLOCK(NIJS:NIJL,18))

IF (CFLAG_P(22).AND.CFLAG_P(23)) THEN
   CALL MEAN_DIRECTION (FL14, THQ=BLOCK(NIJS:NIJL,22),                          &
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

END SUBROUTINE SWELL_SEPARATION

! ---------------------------------------------------------------------------- !

SUBROUTINE PARTITIONING_SWELL(FLH, FLH1)

! ---------------------------------------------------------------------------- !

!     PURPOSE:
!     --------


!     PARTITIONS SPECTRA.

!     SUSANNE HASSELMANN, MPI HAMBURG, 1992.
!     JUERGEN WASZKEWITZ,MPI HAMBURG, 1993 (VECTORIZATION)
!     L. AOUF 2010 ( IMPLEMENTED FOR SWELL PARTITIONING :
!                    DETERMINING THE 1ST AND 2ND SWELL)

!     METHOD:
!     -------

!     THE SPECTRAL PARTITIONING IS ACHIEVED BY FIRST RUNNING THROUGH
!     THE ARRAY OF ALL SPECTRAL GRID POINTS, DETERMINING FOR EACH
!     GRID POINT ITS HIGHEST NEIGHBOUR. EACH POINT IS ASSIGNED A PART
!     NUMBER. THE FOLLOWING CASES OCCUR:
!     - THE GRID POINT AND ITS HIGHEST NEIGHBOUR HAS NO PART NUMBER:
!       THEN BOTH ARE ASSIGNED THE SAME NEW PART NUMBER.
!     - ONE OF THEM HAS ALREADY A PART NUMBER, THE OTHER DOES NOT.
!       THEN THE OTHER IS ASSIGNED THE SAME PART NUMBER.
!     - BOTH, THE GRID POINT AND THE HIGHEST NEIGHBOUR, HAVE PART
!       NUMBERS WHICH ARE DIFFERENT, THEN PART NUMBERS ARE SAVED
!       IN A SPECIAL LIST, SUCH THAT THE PART NUMBER OF THE GRID POINT
!       IS EQUIVALENT TO THE PART NUMBER OF THE HIGHEST NEIGHBOUR.
!     AFTER ASSIGNING THE PART NUMBERS TO EACH GRID POINT, ALL
!     EQUIVALENT PART NUMBERS ARE CHANGED TO ONE UNIQUE PART NUMBER.
!
! ---------------------------------------------------------------------------- !
!
!     INTERFACE:
!     ----------

REAL,INTENT(INOUT), DIMENSION(:,:,:) :: FLH      !! 2D SPECTRA.
REAL,INTENT(OUT), DIMENSION(:,:,:) :: FLH1



! ---------------------------------------------------------------------------- !
!
!     LOCAL VARIABLES.
!     ----------------
!
INTEGER, PARAMETER :: MPART = 30

INTEGER :: NPART(SIZE(FLH,1))      !! ARRAY DIMENSIONS.
INTEGER :: PART(SIZE(FLH,1),SIZE(FLH,2),SIZE(FLH,3)) !! NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
INTEGER, DIMENSION(SIZE(FLH,1),KL*ML/2) :: PEAKANG  !! INDEX OF PEAK DIRECTION OF PARTITIONINGS.
INTEGER, DIMENSION(SIZE(FLH,1),KL*ML/2) :: PEAKFRE  !! INDEX OF PEAK FREQUENCY OF PARTITIONINGS.
INTEGER :: MAXMPART    !! MAXIMUM NUMBER OF PARTITIONINGS USED

INTEGER :: IJ, IPART, JPART , IANG, IFRE, MAXJPART   !! VLOOP INDEXES.
INTEGER :: IANGP1, IANGM1              !! "IANG PLUS 1" , "IANG MINUS 1"
INTEGER :: HANG, HFRE         !! ANGLE AND FREQUENCY INDEX OF HIGHEST NEIGHBOUR
INTEGER, DIMENSION(SIZE(FLH,1)) :: IFREINF, IFRESUP !! FREQUENCY INDEX FOR
                         !! THE LOW/HIGH FREQUENCY SPECTRAL COMPONENT <= FLMIN
INTEGER :: NPARTOVER     !! NUMBER OF TIME NPART EXCEEDS MPART
INTEGER :: EQUIPART(SIZE(FLH,1), KL*ML/2 ) !! CONTAINS EQUIVALENT PART NUMBERS
REAL :: HSPEC            !! VALUE OF HIGHEST NEIGHBOUR
REAL :: HSPEC0           !! FOR TEMPORARY USE
REAL :: FLMIN            !! MINIMUM VALUES FOR SPECTRAL COMPONENTS
REAL :: FSMALL, ZM1      !! A SMALL NUMBER

REAL  :: ZMAX
INTEGER, DIMENSION(SIZE(FLH,1))       :: IPARL

! ---------------------------------------------------------------------------- !

!     1.INTITIALIZE.
!     --------------

!      FSMALL=0.1E-13
FSMALL=0.
FLMIN =0.

FLH1(:,:,:) = FSMALL

NPART(:) = 0
PART(:,:,:) = 0

DO IJ = 1,SIZE(FLH,1)
   DO IPART = 1,KL*ML/2
      EQUIPART(IJ,IPART) = 0
      PEAKANG(IJ,IPART) = 0
      PEAKFRE(IJ,IPART) = 0
   END DO
END DO

! ---------------------------------------------------------------------------- !

!     2.  get rid of spectral values <= FLMIN for low and high frequencies

DO IANG = 1,KL
   IFREINF = 1
   DO IJ = 1,SIZE(FLH,1)
      DO WHILE (FLH(IJ,IANG,IFREINF(IJ)).LE.FLMIN .AND.                      &
&             IFREINF(IJ).LT.ML)
         IFREINF(IJ)=IFREINF(IJ)+1
      ENDDO
   ENDDO
   IFRESUP = ML
   DO IJ = 1,SIZE(FLH,1)
      IF(IFREINF(IJ).LT.ML) THEN
         DO WHILE (FLH(IJ,IANG,IFRESUP(IJ)).LE.FLMIN .AND.                   &
&                                               IFRESUP(IJ).GE.1)
            IFRESUP(IJ) = IFRESUP(IJ)-1
         ENDDO
      ENDIF
   ENDDO
   DO IJ = 1,SIZE(FLH,1)
      IF(IFREINF(IJ).EQ.ML) THEN
         DO IFRE = 1,ML
            FLH(IJ,IANG,IFRE)=FSMALL
         ENDDO
      ELSE
         DO IFRE = IFREINF(IJ)-1,1,-1
            FLH(IJ,IANG,IFRE)=FLMIN/(1+IFREINF(IJ)-IFRE)
         ENDDO
         DO IFRE = IFRESUP(IJ)+1,ML
            FLH(IJ,IANG,IFRE)=FLMIN/(IFRE-IFREINF(IJ)+1)
         ENDDO
      ENDIF
   ENDDO
ENDDO

! ---------------------------------------------------------------------------- !

!     3. PARTITION SPECTRA.
!     ---------------------

DO IANG = 1,KL
   DO IFRE = 1,ML
      DO IJ = 1,SIZE(FLH,1)

         IF( FLH(IJ,IANG,IFRE).LE.FSMALL )THEN
            PART(IJ,IANG,IFRE) = 1
         ELSE
!             FIND HIGHEST NEIGHBOUR

            HANG = IANG
            HFRE = IFRE
            HSPEC = FLH(IJ,HANG,HFRE)
            IF( IANG.NE.KL )THEN
               IANGP1 = IANG + 1
            ELSE
               IANGP1 = 1
            ENDIF
            HSPEC0 = FLH(IJ,IANGP1,IFRE)
            IF( HSPEC0.GE.HSPEC )THEN
               HANG = IANGP1
               HFRE = IFRE
               HSPEC = HSPEC0
            END IF
            IF( IANG.NE.1 )THEN
               IANGM1 = IANG - 1
            ELSE
               IANGM1 = KL
            ENDIF
            HSPEC0 = FLH(IJ,IANGM1,IFRE)
            IF( HSPEC0.GT.HSPEC )THEN
               HANG = IANGM1
               HFRE = IFRE
               HSPEC = HSPEC0
            END IF
            IF( IFRE.NE.ML )THEN
               HSPEC0 = FLH(IJ,IANG,IFRE+1)
               IF( HSPEC0.GE.HSPEC )THEN
                  HANG = IANG
                  HFRE = IFRE+1
                  HSPEC = HSPEC0
               END IF
            END IF
            IF( IFRE.NE.1 )THEN
               HSPEC0 = FLH(IJ,IANG,IFRE-1)
               IF( HSPEC0.GT.HSPEC )THEN
                  HANG = IANG
                  HFRE = IFRE-1
               END IF
            END IF

!             SET PART NUMBER TO POINT OR NEIGHBPOUR

            IF( IANG.EQ.HANG .AND. IFRE.EQ.HFRE )THEN
!               POINT IS PEAK
!               the peak should not belong to the noise level
               IF( FLH(IJ,IANG,IFRE).LE.FLMIN )THEN
                  PART(IJ,IANG,IFRE) = 1
               ELSE
                  IF( PART(IJ,IANG,IFRE).EQ.0 )THEN
                     NPART(IJ) = NPART(IJ) + 1
                     PART(IJ,IANG,IFRE) = NPART(IJ)
                  END IF
                  PEAKANG(IJ,PART(IJ,IANG,IFRE)) = IANG
                  PEAKFRE(IJ,PART(IJ,IANG,IFRE)) = IFRE
               END IF
            ELSE
!               POINT IS NO PEAK
               IF( PART(IJ,IANG,IFRE).NE.0 )THEN
                  IF( PART(IJ,HANG,HFRE).NE.0 ) THEN
                     IF (PART(IJ,IANG,IFRE) .NE. PART(IJ,HANG,HFRE)) THEN
                        EQUIPART(IJ,PART(IJ,IANG,IFRE)) = PART(IJ,HANG,HFRE)
                     ELSE
                        PEAKANG(IJ,PART(IJ,IANG,IFRE)) = HANG
                        PEAKFRE(IJ,PART(IJ,IANG,IFRE)) = HFRE
                     ENDIF
                  ELSE
                     PART(IJ,HANG,HFRE) = PART(IJ,IANG,IFRE)
                  END IF
               ELSE
                  IF( PART(IJ,HANG,HFRE).NE.0 )THEN
                     PART(IJ,IANG,IFRE) = PART(IJ,HANG,HFRE)
                  ELSE
!                   a new part should not be assigned to noise level data
                     IF( FLH(IJ,IANG,IFRE).LE.FLMIN )THEN
                        PART(IJ,IANG,IFRE) = 1
                     ELSE
                        NPART(IJ) = NPART(IJ) + 1
                        PART(IJ,IANG,IFRE) = NPART(IJ)
                        PART(IJ,HANG,HFRE) = NPART(IJ)
                     ENDIF
                  ENDIF
               END IF
            END IF
         END IF
      END DO
   END DO
END DO


!    3. ASSIGN EQUIVALENT PART NUMBERS TO ONE UNIQUE NUMBER.
!     ------------------------------------------------------

!     FIRST: ASSIGN TO EVERY PART NUMBER THIS UNIQUE NUMBER

NPARTOVER=0

DO IJ = 1,SIZE(FLH,1)
   JPART = 0
   DO IPART = 1,NPART(IJ)
      IF( EQUIPART(IJ,IPART).EQ.0 )THEN
         JPART = JPART + 1
         EQUIPART(IJ,IPART) = -JPART
         PEAKANG(IJ,JPART) = PEAKANG(IJ,IPART)
         PEAKFRE(IJ,JPART) = PEAKFRE(IJ,IPART)
      END IF
   END DO
   NPART(IJ) = JPART
   IF(NPART(IJ).GT.MPART) THEN
      NPARTOVER=NPARTOVER+1
      NPART(IJ)=MPART
   ENDIF
   IF(NPART(IJ).EQ.0) THEN
      NPART(IJ)=1
      PEAKANG(IJ,NPART(IJ))=1
      PEAKFRE(IJ,NPART(IJ))=1
   ENDIF

!       NOW NPART(IJ) CONTAINS THE NUMBER OF PARTITIONS
!       AND THE LIST EQUIPART( , ) ENDS WITH THE FIRST 0.

   IPART = 1
   DO WHILE( EQUIPART(IJ,IPART).NE.0.AND.IPART.LE.KL*ML/2 )
      JPART = IPART
      DO WHILE( JPART.GE.0 )
         IF(JPART.eq.0 .or. (JPART.eq.EQUIPART(IJ,JPART)) ) then
            WRITE(IU06,*) ' '
            WRITE(IU06,*) '  SUBROUTING PARTITIONING :'
            WRITE(IU06,*) ' '
            WRITE(IU06,*) ' !!! problem !!! '
            WRITE(IU06,*) ' IJ = ', IJ
            WRITE(IU06,*) ' JPART = ', JPART
            WRITE(IU06,*) ' EQUIPART(IJ,JPART) = ', EQUIPART(IJ,JPART)
            WRITE(IU06,*) ' '
            CALL FLUSH(IU06)
            JPART=-1
         ENDIF
         JPART = EQUIPART(IJ,JPART)
      END DO
      EQUIPART(IJ,IPART) = JPART
      IPART = IPART + 1
   END DO
END DO

IF( NPARTOVER .GT. 0 )THEN
   WRITE(IU06,*) '************************************************'
   WRITE(IU06,*) 'WARNING SUBROUTINE PARTITIONING:'
   WRITE(IU06,*) 'DIMENSION MPART MIGHT BE TOO SMALL'
   WRITE(IU06,*) 'FOR THE PARTITIONING OF SPECTRA '
   WRITE(IU06,*) 'IT WAS FOUND TO EXCEED ITS MAXIMUN VALUE'
   WRITE(IU06,*) 'MPART =', MPART,NPARTOVER,' TIMES.'
   WRITE(IU06,*) 'IF IT HAPPENS TOO OFTEN,'
   WRITE(IU06,*) 'CHANGE PARAMETER MPARTSW IN MODULE YOWSARAS'
   WRITE(IU06,*) '************************************************'
ENDIF

MAXJPART=0
! modif
MAXMPART=0
DO IJ = 1,SIZE(FLH,1)
   MAXJPART = MAX(MAXJPART,NPART(IJ))
END DO
MAXMPART = MAX(MAXMPART,MAXJPART)

WRITE(IU06,*) ' '
WRITE(IU06,*) '  SUBROUTING PARTITIONING :'
WRITE(IU06,*) '  THE MAXIMUM NUMBER OF PARTITIONING IS: ', MAXJPART
WRITE(IU06,*) '  FOR A TOTAL NUMBER OF PARTITIONING OF: ', MAXMPART
WRITE(IU06,*) ' '
CALL FLUSH(IU06)


IF( MAXJPART.GT.MPART )THEN
   WRITE(IU06,*) '************************************************'
   WRITE(IU06,*) 'ERROR! SUBROUTINE PARTITIONING:'
   WRITE(IU06,*) 'DIMENSION MPART IS TOO SMALL'
   WRITE(IU06,*) 'FOR THE PARTITIONING OF SPECTRA '
   WRITE(IU06,*) 'IT MUST BE MPART >= ',MAXJPART
   WRITE(IU06,*) 'CHANGE PARAMETER MPARTSW'
   WRITE(IU06,*) 'IN MODULE YOWSARAS'
   WRITE(IU06,*) '************************************************'
   CALL ABORT
ENDIF


!     THEN: ASSIGN TO EVERY POINT THIS UNIQUE PART NUMBER

DO IANG = 1,KL
   DO IFRE = 1,ML
      DO IJ = 1,SIZE(FLH,1)
         PART(IJ,IANG,IFRE) = MAX( 1 , -EQUIPART(IJ,PART(IJ,IANG,IFRE)))
         PART(IJ,IANG,IFRE) = MIN( MPART , PART(IJ,IANG,IFRE))
      END DO
   END DO
END DO

!     restore noise level
DO IANG = 1,KL
   DO IFRE = 1,ML
      DO IJ = 1,SIZE(FLH,1)
         FLH(IJ,IANG,IFRE)=MAX(FLH(IJ,IANG,IFRE),FLMIN)
      END DO
   END DO
END DO

DO IJ = 1,SIZE(FLH,1)
   ZM1 = FLH(IJ,PEAKANG(IJ,1),PEAKFRE(IJ,1))
   IPARL(IJ)=1
   DO IPART=1,NPART(IJ)
      ZMAX = FLH(IJ,PEAKANG(IJ,IPART),PEAKFRE(IJ,IPART))
      IF (ZMAX.GT.ZM1) THEN
         ZM1 = ZMAX
         IPARL(IJ)=IPART
      ENDIF
   ENDDO
ENDDO


DO IFRE = 1,ML
   DO IANG = 1,KL
      DO IJ = 1,SIZE(FLH,1)
         IF (PART(IJ,IANG,IFRE).EQ.IPARL(IJ)) THEN
            FLH1(IJ,IANG,IFRE)=FLH(IJ,IANG,IFRE)
         ENDIF
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE PARTITIONING_SWELL

! ---------------------------------------------------------------------------- !

END MODULE WAM_SWELL_MODULE

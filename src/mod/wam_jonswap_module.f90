MODULE WAM_JONSWAP_MODULE

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

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

USE WAM_GENERAL_MODULE, ONLY:  G, ZPI

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C.  MODULE DATA.                                                         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

IMPLICIT NONE
PRIVATE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTERFACE FETCH_LAW                  !! JONSWAP ALPHA and PEAK FREQUENCY FROM
   MODULE PROCEDURE FETCH_LAW        !! FETCH AND WINDSPEED.
END INTERFACE
PUBLIC FETCH_LAW

INTERFACE JONSWAP                    !! JONSWAP SPECTRA FROM PARAMTERS.
   MODULE PROCEDURE JONSWAP
END INTERFACE
PUBLIC JONSWAP

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FETCH_LAW (FETCH, FPMAX, U10, ALPHJ, FP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   FETCH_LAW - COMPUTE JONSWAP PARAMETERS FROM FETCH LAW.                     !
!                                                                              !
!     S. HASSELMANN  - JULY 1990                                               !
!     H. GUNTHER     - DECEMBER 1990   MODIFIED FOR CYCLE_4.                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTES FOR EACH GRID POINT THE PEAK FREQUENCY FROM A FETCH LAW       !
!       AND THE JONSWAP ALPHA FROM THE ALPHA NY RELATION.                      !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       FP = A * (G*FETCH/U_10**2)**D    A = 2.84                              !
!       FP = MAX [FP, 0.13]              D = -3./10.                           !
!       FP = MIN [FP, FRMAX*U_10/G]                                            !
!       ALPHJ = B * FP**2/3              B = 0.033                             !
!       ALPHJ = MAX [ALPHJ, 0.0081]                                            !
!       FP = G/U_10*FP                                                         !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       K.HASSELMAN,D.B.ROOS,P.MUELLER AND W.SWELL                             !
!          A PARAMETRIC WAVE PREDICTION MODEL                                  !
!          JOURNAL OF PHSICAL OCEANOGRAPHY, VOL. 6, NO. 2, MARCH 1976.         !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)  :: FETCH         !! FETCH TO BE USED (METRES).
REAL,    INTENT(IN)  :: FPMAX         !! MAXIMUM PEAK FREQUENCY (HERTZ).
REAL,    INTENT(IN)  :: U10(:)        !! MODULUS OF WIND VELOCITY [M/S].
REAL,    INTENT(OUT) :: ALPHJ(:)      !! JONSWAP ALPHA.
REAL,    INTENT(OUT) :: FP(:)         !! JONSWAP PEAK FREQUENCY [HERTZ].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: A = 2.84,  D = -(3./10.) !! PEAKFREQUENCY FETCH LAW CONSTANTS
REAL, PARAMETER :: B = 0.033, E = 2./3.     !! ALPHA-PEAKFREQUENCY LAW CONSTANTS

REAL :: UG(SIZE(U10))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE VALUES FROM FETCH LAWS.                                       !
!        -------------------------------                                       !

WHERE (U10 .GT. 0.1E-08)
   UG = G/U10
   FP = MAX(0.13, A*((G*FETCH)/(U10**2))**D)
   FP = MIN(FP, FPMAX/UG)
   ALPHJ = MAX(0.0081, B * FP**E)
   FP = FP*UG
ELSEWHERE
   ALPHJ = 0.0081
   FP = FPMAX
END WHERE

END SUBROUTINE FETCH_LAW

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE JONSWAP (FR, ALPHAJ, GAMMA, SA, SB, FP, ET)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   JONSWAP - ROUTINE TO COMPUTE THE 1-D JONSWAP SPECTRUM.                     !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!     TO COMPUTE A 1-D JONSWAP SPECTRUM.                                       !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL, INTENT(IN)  :: FR(:)       !! FREQUENCiIES.
REAL, INTENT(IN)  :: ALPHAJ(:)   !! OVERALL ENERGY LEVEL OF JONSWAP SPECTRA.
REAL, INTENT(IN)  :: GAMMA       !! OVERSHOOT FACTOR.
REAL, INTENT(IN)  :: SA          !! LEFT PEAK WIDTH.
REAL, INTENT(IN)  :: SB          !! RIGHT PEAK WIDTH.
REAL, INTENT(IN)  :: FP(:)       !! PEAK FREQUENCIES.
REAL, INTENT(OUT) :: ET(:,:)     !! JONSWAP SPECTRA.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: I, M
REAL    :: FRH, LOG_GAMMA
REAL    :: G2ZPI4FRH5M(1:SIZE(FR))
REAL    :: ARG, SAA, SBB

! ---------------------------------------------------------------------------- !

G2ZPI4FRH5M(:) = G**2/ZPI**4*FR(:)**(-5)
LOG_GAMMA = LOG(GAMMA)
SAA = 0.5/(SA*SA)
SBB = 0.5/(SB*SB)

DO M = 1,SIZE(FR)
   FRH = FR(M)
   DO I = 1,SIZE(FP)
      ARG = 1.25*(FP(I)/FRH)**4
      IF (ARG.LT.50.) THEN
         ET(I,M) = ALPHAJ(I)*G2ZPI4FRH5M(M)*EXP(-ARG)
      ELSE
         ET(I,M) = 0.
         CYCLE
      END IF
      IF (FRH.GT.FP(I)) THEN
         ARG = SBB*(FRH/FP(I)-1.)**2
      ELSE
         ARG = SAA*(FRH/FP(I)-1.)**2
      ENDIF
      IF (ARG.LT.99.) ET(I,M) = ET(I,M)*EXP(LOG_GAMMA*EXP(-ARG))
   END DO
END DO

END SUBROUTINE JONSWAP

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_JONSWAP_MODULE

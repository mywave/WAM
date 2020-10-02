MODULE WAM_GENERAL_MODULE

! ---------------------------------------------------------------------------- !
!          THIS MODULE COLLECTS ALL GLOBAL CONSTANTS AND GENERAL SUBROUTINES   !
!                       USED IN THE WAM MODEL PROGRAMS                         !
!                                                                              !
!    JUNE 2005:                                                                !
!       CHARNOCK CONSTANT CHANGED FROM ALPHA = 0.0100 TO ALPHA = 0.0095        !
!       (Jean Bidlot, Peter Janssen and Saleh Abdalla: A revised formulation   !
!       for ocean wave dissipation in CY29R1. ECMWF MEMORANDUM RESEARCH        !
!       DEPARTMENT:April 7, 2005 File: R60.9/JB/0516                           !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE, ONLY: IU06

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!
!    1. WAM MODEL GLOBAL CONSTANTS.

REAL, PARAMETER :: G = 9.806         !! ACCELERATION OF GRAVITY [M/S**2]
REAL, PARAMETER :: PI = 3.1415927    !! PI.
REAL, PARAMETER :: ZPI = 2.*PI       !! 2.* PI.
REAL, PARAMETER :: ZPISQRT = 1.7724539  !! SQRT(PI)
REAL, PARAMETER :: DEG = 180./PI     !! COVERTION FROM RADIANS TO DEGREE
REAL, PARAMETER :: RAD = PI/180.     !! COVERTION FROM DEGREE TO RADIANS 
REAL, PARAMETER :: CIRC = 40000000.  !! EARTH CIRCUMFERENCE [M].
REAL, PARAMETER :: R = CIRC/ZPI      !! EARTH RADIUS [M].
REAL, PARAMETER :: GAMMA_E = 0.57721566 !! EULER CONSTANT                      !! WAM-MAX

! ---------------------------------------------------------------------------- !
!
!    2. PARAMETERS FOR COUPLING.
!       ------------------------

REAL, PARAMETER :: ROAIR   = 1.225      !! AIR DENSITY
REAL, PARAMETER :: RNUAIR  = 1.5E-5     !! KINEMATIC AIR VISCOSITY
REAL, PARAMETER :: RNUAIRM = 0.11*RNUAIR !! KINEMATIC AIR VISCOSITY FOR MOMENTUM TRANSFER
REAL, PARAMETER :: ROWATER = 1000.      !! WATER DENSITY
REAL, PARAMETER :: XEPS = ROAIR/ROWATER
REAL, PARAMETER :: XINVEPS = 1./XEPS

REAL            :: BETAMAX = 1.20       !! PARAMETER FOR WIND INPUT (ECMWF CY45R1).
REAL            :: ZALP    = 0.0080     !! SHIFTS GROWTH CURVE (ECMWF CY45R1).
REAL            :: ALPHA   = 0.0060     !! MINIMUM CHARNOCK CONSTANT (ECMWF CY45R1).
                                        !! if LE 30 frequencies changed
                                        !! to 0.0075 in subroutine INITMDL
REAL            :: TAUWSHELTER= 0.0     !! SHELTERING COEFFICIENT  (ECMWF CY45R1)

REAL            :: BETAMAX_ARD = 1.42   !! PARAMETER FOR WIND INPUT (ECMWF CY46R1, based on Ardhuin et al. 2010).
REAL            :: ZALP_ARD= 0.0080     !! SHIFTS GROWTH CURVE (ECMWF CY46R1, based on Ardhuin et al. 2010).
REAL            :: ALPHA_ARD= 0.0062    !! MINIMUM CHARNOCK CONSTANT (ECMWF CY46R1, based on Ardhuin et al. 2010).
REAL            :: TAUWSHELTER_ARD= 0.25 !! SHELTERING COEFFICIENT (ECMWF CY46R1, based on Ardhuin et al. 2010).

REAL, PARAMETER :: XKAPPA  = 0.40       !! VON KARMAN CONSTANT.
REAL, PARAMETER :: XNLEV   = 10.0       !! WINDSPEED REF. LEVEL.
REAL, PARAMETER :: RCHAR   = 0.018      !! DEFAULT CHARNOCK VALUE FOR ICE AND
                                        !! DRY SEA.
! ---------------------------------------------------------------------------- !
!
!    3. PARAMETERS FOR BFI.
!       -------------------

REAL, PARAMETER :: DKMAX = 40.        !! MAXIMUM VALUE OF DEPTH*WAVENUMBER.
REAL, PARAMETER :: BF2MIN = -10.      !! MINIMUM VALUE ALLOWED FOR BFI SQUARED.
REAL, PARAMETER :: BF2MAX = 10.       !! MAXIMUM VALUE ALLOWED FOR BFI SQUARED.
REAL, PARAMETER :: C4MIN = -0.33      !! MINIMUM VALUE ALLOWED FOR KURTOSIS.
REAL, PARAMETER :: C4MAX = 1.         !! MAXIMUM VALUE ALLOWED FOR KURTOSIS.

! ---------------------------------------------------------------------------- !
!
!    4. PARAMETERS FOR SDISSIP_ARD ARDHUIN et al. 2010
!       ----------------------------------------------
!     Br:
REAL, PARAMETER :: SDSBR = 9.0E-4

!!     Saturation dissipation coefficient
INTEGER, PARAMETER :: ISDSDTH = 80
INTEGER, PARAMETER :: ISB=2
INTEGER, PARAMETER :: IPSAT=2
REAL, PARAMETER :: SSDSC2 = -2.2E-5
REAL, PARAMETER :: SSDSC4 = 1.0
REAL, PARAMETER :: SSDSC6 = 0.3
REAL, PARAMETER :: MICHE = 1.0

!!     Cumulative dissipation coefficient
!!! REAL, PARAMETER :: SSDSC3 = -0.40344
!!!   This is quite an expensive computation. Setting it to 0 will disable its calculation.
!!!   It was found that if the high frequency tail is prescribed (see frcutindex),
!!!   then the results are very similar.
!!!   Converserly, it will be required, in particular for the wave modified fluxes to NEMO
!!!   when the high frequency tail is not prescribed and used in all calculation
REAL, PARAMETER :: SSDSC3 = 0.0

REAL, PARAMETER :: SSDSBRF1 = 0.5
!  28.16 = 22.0 * 1.6² * 1/2 with
!  22.0 (Banner & al. 2000, figure 6)
!  1.6  the coefficient that transforms  SQRT(B) to Banner et al. (2000)'s epsilon
!  1/2  factor to correct overestimation of Banner et al. (2000)'s breaking probability due to zero-crossing analysis
REAL, PARAMETER :: BRKPBCOEF=28.16

!! Wave-turbulence interaction coefficient
REAL, PARAMETER :: SSDSC5  = 0.0

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE ABORT1                        !! TERMINATE PROCESSING.
   MODULE PROCEDURE ABORT1
END INTERFACE

INTERFACE AKI                           !! WAVE NUMBER FROM FREQUENCY AND DEPTH.
   MODULE PROCEDURE AKI
END INTERFACE

interface flush1
   module procedure flush1              !! enforce output
end interface

INTERFACE DIFDATE                       !! COMPUTE TIME DIFFERENCE.
   MODULE PROCEDURE DIFDATE
END INTERFACE

INTERFACE INCDATE                       !! UPDATE DATE TIME GROUP.
   MODULE PROCEDURE INCDATE
END INTERFACE

INTERFACE OPEN_FILE                     !! OPEN A FILE.
   MODULE  PROCEDURE OPEN_FILE
END INTERFACE

INTERFACE PRINT_ARRAY                   !! FORMATED OUTPUT OF AN ARRAY. 
   MODULE PROCEDURE PRINT_ARRAY_C
   MODULE PROCEDURE PRINT_ARRAY_CC
   MODULE PROCEDURE PRINT_ARRAY_R
   MODULE PROCEDURE PRINT_ARRAY_RC
   MODULE PROCEDURE PRINT_ARRAY_L
   MODULE PROCEDURE PRINT_ARRAY_I
   MODULE PROCEDURE PRINT_ARRAY_IC
END INTERFACE

INTERFACE PRINT_SPECTRUM                !! PRINTS A SPECTRUM. 
   MODULE  PROCEDURE PRINT_SPECTRUM
   MODULE  PROCEDURE PRINT_SPECTRUM_C
END INTERFACE

INTERFACE REDUCED_TO_REGULAR            !! INTERPOLATE REDUCED TO REGULAR GRID. 
   MODULE PROCEDURE REDUCED_TO_REGULAR
   MODULE PROCEDURE REDUCED_TO_REGULAR_C
END INTERFACE
PUBLIC REDUCED_TO_REGULAR

INTERFACE PRINT_GENERAL_MODULE          !! PRINT OUTPUT OF MODULE DATA.
   MODULE PROCEDURE PRINT_GENERAL_MODULE
END INTERFACE
PUBLIC PRINT_GENERAL_MODULE

INTERFACE SET_GENERAL_MODULE            !! TRANSFER OF DATA INTO MODULE.
   MODULE PROCEDURE SET_GENERAL_MODULE
END INTERFACE
PUBLIC SET_GENERAL_MODULE

INTERFACE SORT               !! SORT A DATA ARRAY BY AN INDEX ARRAY.
   MODULE PROCEDURE SORTI       !! INTEGER DATA ARRAY
END INTERFACE
PUBLIC SORT

INTERFACE SORTIN             !! FIND THE INDEX ARRAY WHICH SORTS THE DATA.
   MODULE PROCEDURE SORTINI     !! INTEGER DATA ARRAY
END INTERFACE
PUBLIC SORTIN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     B.  EXPLICIT INTERFACES.                                                 !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE
      
! ---------------------------------------------------------------------------- !

   SUBROUTINE INITMDL                    !! INITIALIZES THE WAM MODEL.
   END SUBROUTINE INITMDL

! ---------------------------------------------------------------------------- !

   SUBROUTINE PRINT_PREPROC_OUTPUT  !! PRINT PREPROC_OUTPUT.
   END SUBROUTINE PRINT_PREPROC_OUTPUT 
      
! ---------------------------------------------------------------------------- !

   SUBROUTINE READ_BOUNDARY_INPUT     !! READ BOUNDARY VALUE INPUT FILE.
   END SUBROUTINE READ_BOUNDARY_INPUT
   
! ---------------------------------------------------------------------------- !

   SUBROUTINE READ_WAM_USER              !! READ USER INPUT FOR WAM.
   END SUBROUTINE READ_WAM_USER
   
! ---------------------------------------------------------------------------- !

   SUBROUTINE WAMODEL                !! TIME INTEGRATION OF WAVE FIELDS.
   END SUBROUTINE WAMODEL
   
! ---------------------------------------------------------------------------- !

   SUBROUTINE WAVEMDL          !! SUPERVISES EXECUTION OF THE WAVE MODEL.
   END SUBROUTINE WAVEMDL
   
END INTERFACE 

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE ABORT1

! ---------------------------------------------------------------------------- !
!                                                                              !
!   ABORT1 - STOP PROCESSING                                                   !
!                                                                              !
!     H. GUNTHER     GKSS    SEPTEMBER 2000                                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!     CLOSE ALL INPUT AND OUTPUT UNITS AND STOP                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!      NONE                                                                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     LOCAL VARIABLES.
!     ----------------

LOGICAL :: DA
INTEGER :: I

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.0 CHECK FILE STATUS AND CLOSE FILE IF FILE IS ASSIGNED.                !
!         ------------------------------------------------------               !

WRITE (IU06,*) ' *****************************************************'
WRITE (IU06,*) ' *                                                   *'
WRITE (IU06,*) ' *      PROGRAM ABORTS  PROGRAM ABORTS               *'
WRITE (IU06,*) ' *      ==============================               *'
WRITE (IU06,*) ' *                                                   *'
WRITE (IU06,*) ' *         CHECK MESSAGE PRINTED ABOVE.              *'
WRITE (IU06,*) ' *                                                   *'
WRITE (IU06,*) ' *****************************************************'

DO I=1,99
   INQUIRE (UNIT=I, EXIST=DA)
   IF (DA) CLOSE (UNIT=I)
END DO
   call MPI_finalize (i)

STOP 1
END SUBROUTINE ABORT1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

REAL FUNCTION AKI (OM, BETA)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   AKI - FUNCTION TO COMPUTE WAVE NUMBER.                                     !
!                                                                              !
!     G. KOMEN, P. JANSSEN   KNMI        01/06/1986                            !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       WAVE NUMBER AS FUNCTION OF CIRCULAR FREQUENCY AND WATER DEPTH.         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NEWTONS METHOD TO SOLVE THE DISPERSION RELATION IN SHALLOW WATER.      !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

REAL, INTENT(IN) :: OM    !! CIRCULAR FREQUENCY.
REAL, INTENT(IN) :: BETA  !! WATER DEPTH.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: EBS = 0.0001  !! RELATIVE ERROR LIMIT OF NEWTON'S METHOD.

REAL :: AKP, BO, TH, STH

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. START WITH MAXIMUM FROM DEEP AND EXTREM SHALLOW WATER WAVE NUMBER.    !
!        ------------------------------------------------------------------    !

AKI   = MAX ( OM**2/(4.*G), OM/(2.*SQRT(G*BETA)) )

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. ITERATION LOOP.                                                       !
!        ---------------                                                       !

AKP = 10000.
DO WHILE (ABS(AKP-AKI).GT.EBS*AKI)
   BO = BETA*AKI
   IF (BO.GT.DKMAX) THEN
      AKI = OM**2/G
      EXIT
   ELSE
      AKP = AKI
      TH = G*AKI*TANH(BO)
      STH = SQRT(TH)
      AKI = AKI+(OM-STH)*STH*2./(TH/AKI+G*BO/COSH(BO)**2)
   END IF
END DO

END FUNCTION AKI

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine flush1 (iu)
   integer :: iu
   call flush (iu)
end subroutine flush1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE DIFDATE (CDATE1, CDATE2, ISHIFT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   DIFDATE - TO COMPUTE TIME DIFFERENCE.                                      !
!                                                                              !
!     H. GUNTHER   GKSS/ECMWF  NOVEMBER 1989                                   !
!     H. GUNTHER   GKSS   NOVEMBER 1999    FT90 AND CENTURY AND SECONDS.       !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE THE SECONDS BETWEEN THE INPUT DATES.                           !
!       DATES HAVE TO BE IN CONSECUTIVE YEARS.                                 !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!      REFERENCES.                                                             !
!      -----------                                                             !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!    INTERFACE VARIABLES.
!    --------------------

CHARACTER (LEN=14), INTENT(IN)  :: CDATE1 !! DATE TIME GROUP (YYYYMMDDHHMMSS).
CHARACTER (LEN=14), INTENT(IN)  :: CDATE2 !! DATE TIME GROUP (YYYYMMDDHHMMSS).
INTEGER,            INTENT(OUT) :: ISHIFT !! DIFFERENCE IN SECONDS 
                                          !! (CDATE2-CDATE1).

! ---------------------------------------------------------------------------- !
! 
!    LOCAL VARIABLES.
!    ----------------

INTEGER, SAVE      ::  MON(12) =(/31,28,31,30,31,30,31,31,30,31,30,31/)

INTEGER            ::  YEAR1, MONTH1, DAY1, HOUR1, MINUTE1, SECOND1
INTEGER            ::  YEAR2, MONTH2, DAY2, HOUR2, MINUTE2, SECOND2

INTEGER            :: M, MDAY
CHARACTER (LEN=14) :: CDT1, CDT2
LOGICAL            :: SWITCH

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.0 CHANGE DATE TIME GROUPS TO ENSURE THAT THE SECOND IS LARGER.         !
!         ------------------------------------------------------------         !

IF (CDATE1 .GT. CDATE2) THEN
   CDT1 = CDATE2
   CDT2 = CDATE1
   SWITCH = .TRUE.
ELSE IF (CDATE2 .GT. CDATE1) THEN
   CDT1 = CDATE1
   CDT2 = CDATE2
   SWITCH = .FALSE.
ELSE
   ISHIFT = 0
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.0 SPLITE DATE TIME GROUP INTO SECONDS, MINUTE, HOUR, DAY, MONTH, YEAR. !
!         -------------------------------------------------------------------- !

READ (CDT1,'(I4,5I2)') YEAR1, MONTH1, DAY1, HOUR1, MINUTE1, SECOND1
READ (CDT2,'(I4,5I2)') YEAR2, MONTH2, DAY2, HOUR2, MINUTE2, SECOND2


! ---------------------------------------------------------------------------- !
!
!     3.0 COMPUTE DIFFERENCE BETWEEN DAY, HOUR ,MINITE AND SECOND.
!         --------------------------------------------------------

ISHIFT = (((DAY2-DAY1)*24+HOUR2-HOUR1)*60+MINUTE2-MINUTE1)*60+SECOND2-SECOND1

! ---------------------------------------------------------------------------- !
! 
!   4.0 ADD DIFFERENCE FROM MONTH.
!       --------------------------

IF (YEAR2.GT.YEAR1) THEN

!   4.1 START AND END MONTH ARE IN DIFFERENT YEARS.

   DO M = MONTH1,12
      MDAY = MON(M)
      IF (M.EQ.2 .AND. MOD(YEAR1,4).EQ.0) THEN
         IF (MOD(YEAR1,400).EQ.0 .OR. MOD(YEAR1,100).NE.0) MDAY = 29
      END IF
      ISHIFT =ISHIFT + MDAY*24*3600
   END DO
   DO M = 1,MONTH2-1
      MDAY =MON(M)
      IF (M.EQ.2 .AND. MOD(YEAR2,4).EQ.0) THEN
         IF (MOD(YEAR2,400).EQ.0 .OR. MOD(YEAR2,100).NE.0) MDAY = 29
      END IF
      ISHIFT = ISHIFT + MDAY*24*3600
   END DO

ELSE

!   4.2 START AND END MONTH ARE IN THE SAME YEAR.

   DO M = MONTH1,MONTH2-1
      MDAY = MON(M)
      IF (M.EQ.2 .AND. MOD(YEAR1,4).EQ.0) THEN
         IF (MOD(YEAR1,400).EQ.0 .OR. MOD(YEAR1,100).NE.0) MDAY = 29
      END IF
      ISHIFT = ISHIFT + MDAY*24*3600
   END DO

END IF

! ---------------------------------------------------------------------------- !
! 
!   5.0 CHANGE SIGN OF DIFFERENCE IF DATES WERE EXCHANGED.
!       --------------------------------------------------

IF (SWITCH) ISHIFT = -ISHIFT

END SUBROUTINE DIFDATE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE INCDATE (CDATE, ISHIFT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   INCDATE - TO UPDATE DATE TIME GROUP                                        !
!                                                                              !
!     L. BERTOTTI, P.JANSSEN.                                                  !
!                                                                              !
!     H. GUNTHER   ECMWF  NOVEMBER 1989    NEGATIVE INCREMENTS.                !
!     H. GUNTHER   GKSS   NOVEMBER 2001    FT90 AND CENTURY AND SECONDS.       !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       UPDATING DATE TIME GROUP.                                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

CHARACTER (LEN=14), INTENT(INOUT) :: CDATE  !! DATE TIME GROUP (YYYYMMDDHHMMSS).
INTEGER,            INTENT(IN)    :: ISHIFT !! TIME INCREMENT IN SECONDS.

! ---------------------------------------------------------------------------- !
! 
!     LOCAL VARIABLES.
!     ----------------

INTEGER, SAVE ::  MON(12) =(/31,28,31,30,31,30,31,31,30,31,30,31/)

INTEGER ::  YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
INTEGER ::  DELT, MDAY

! ---------------------------------------------------------------------------- !
! 
!   1.0 RETURN IF TIME INCREMENT IS ZERO.
!       ---------------------------------

DELT = ISHIFT
IF (ABS(DELT).EQ.0) RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    2.0 SPLITE DATE TIME GROUP INTO SECONDS, MINUTE, HOUR, DAY, MONTH, YEAR. !
!         -------------------------------------------------------------------- !

READ (CDATE,'(I4,5I2)') YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.0 ADD AND CHECK SECONDS.                                               !
!         ----------------------                                               !
!                                                                              !
!   2.1 IF SECONDS ARE BETWEEN 0 AND 60 RETURN
!       --------------------------------------

SECOND = SECOND + DELT
IF (SECOND.GE.0. .AND. SECOND.LT.60.) THEN
   WRITE(CDATE,'(I4.4,5I2.2)') YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
   RETURN
END IF

!*  2.2 NEW MIMUTES AND SECONDS.
!       ------------------------

DELT = MODULO(SECOND,60)
MINUTE = MINUTE +(SECOND-DELT)/60
SECOND = DELT

! ---------------------------------------------------------------------------- !
! 
!   3.0 CHECK MINUTES.
!       --------------

IF (MINUTE.GE.60) THEN
   HOUR = HOUR + MINUTE/60        !! MINUTES > 59 ==> NEW HOURS.
ELSE IF (MINUTE.LT.0) THEN
   HOUR = HOUR + (MINUTE-59)/60   !! MINUTES < 0  ==> NEW HOURS.
ELSE
   WRITE (CDATE,'(I4.4,5I2.2)') YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
   RETURN                         !! ALL DONE  ==>  RETURN
END IF 
MINUTE = MODULO(MINUTE,60)        !! NEW MINUTES.

! ---------------------------------------------------------------------------- !
! 
!   4.0 CHECK HOURS.
!       ------------

IF (HOUR.GE.24) THEN
   DAY =  DAY + HOUR/24           !! HOURS > 23 ==> NEW DAYS.
ELSE IF (HOUR.LT.0) THEN
   DAY =  DAY + (HOUR-23)/24      !! HOURS < 0  ==> NEW DAYS.
ELSE
   WRITE (CDATE,'(I4.4,5I2.2)') YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
   RETURN                         !! ALL DONE  ==>  RETURN
END IF 
HOUR = MODULO(HOUR,24)            !! NEW HOURS.

! ---------------------------------------------------------------------------- !
! 
!   5.0 CHECK DAYS.
!       ------------
! 
!   5.1 IF DAYS ARE GREATER THAN DAYS OF MONTH. NEW DAY AND MONTH AND YEAR.
!       -------------------------------------------------------------------

MDAY = MON(MONTH)
IF (MONTH.EQ.2 .AND. MOD(YEAR,4).EQ.0) THEN
   IF (MOD(YEAR,400).EQ.0 .OR. MOD(YEAR,100).NE.0) MDAY = 29
END IF

DO WHILE (DAY > MDAY)
   DAY =  DAY - MDAY
   MONTH = MONTH+1
   IF (MONTH.GE.13) THEN
      YEAR = YEAR+1
      MONTH = MONTH-12
   END IF
   MDAY = MON(MONTH)
   IF (MONTH.EQ.2 .AND. MOD(YEAR,4).EQ.0) THEN
      IF (MOD(YEAR,400).EQ.0 .OR. MOD(YEAR,100).NE.0) MDAY = 29
   END IF
END DO

!   5.2 IF DAYS ARE LESS THAN 1. NEW DAY AND MONTH AND YEAR.
!       ----------------------------------------------------

DO WHILE ( DAY < 1)
   MONTH = MONTH-1
   IF (MONTH.EQ.0) THEN
      MONTH = 12
      YEAR = YEAR-1
   END IF
   MDAY = MON(MONTH)
   IF (MONTH.EQ.2 .AND. MOD(YEAR,4).EQ.0) THEN
      IF (MOD(YEAR,400).EQ.0 .OR. MOD(YEAR,100).NE.0) MDAY = 29
   END IF
   DAY = DAY + MDAY
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6.0 COMPOSE NEW DATE TIME GROUP.                                         !
!         ----------------------------                                         !

WRITE (CDATE,'(I4.4,5I2.2)') YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

END SUBROUTINE INCDATE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE OPEN_FILE (IU06, IUNIT, FILEID, CDATE, STAT, IFAIL, FORMA)

! ---------------------------------------------------------------------------- !
!
!   OPEN_FILE - OPEN A FILE.
!
!     HEINZ GUNTHER       ECMWF       OCTOBER 1989
!     HEINZ GUNTHER       ECMWF       OCTOBER 1990  MIGRATION TO YMP
!                                                   NEW FILE NAMES
!     H. GUNTHER    GKSS              NOVEMBER 1999  NEW DATES AND FT90.
!
!     PURPOSE.
!     --------
!
!         INCLUDE DATE AND TIME IN A FILE MANE AND ASSIGN THE FILE TO A UNIT. 
!
!     METHOD.
!     -------
!
!        THE FILENAME IS BUILT FROM THE FILEID FOLLOWED BY 'YYYYMMDDHHMMSS' 
!        WHERE YYYYMMDDHHMMSS IS THE TIME OF THE LAST FIELD STORED IN THE FILE.
!        THE FULL PATH NAME IS  USERID/PATH/RUNID/FILENAME
!
!     REFERENCES.
!     -----------
!
!         NONE
!
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

IMPLICIT NONE

INTEGER,            INTENT(IN)  :: IU06       !! PROTOCOL UNIT
INTEGER,            INTENT(IN)  :: IUNIT      !! FORTRAN UNIT FOR FILE.
CHARACTER (LEN=*),  INTENT(IN)  :: FILEID     !! FILE IDENTIFIER
CHARACTER (LEN=14), INTENT(IN)  :: CDATE      !! DATE (YYYYMMDDHHMMSS) OF FILE.
CHARACTER (LEN=*),  INTENT(IN)  :: STAT       !! FILE STATUS (NEW, OLD, UNKNOWN)
INTEGER,            INTENT(OUT) :: IFAIL      !! ERROR FLAG    = 0 NO ERROR
                                              !!               ­ 0 OPEN ERROR
CHARACTER (LEN=*),  INTENT(IN), OPTIONAL :: FORMA !! FORMATTED OR UNFORMATED

! ---------------------------------------------------------------------------- !
!
!     LOCAL VARIABLES.
!     ----------------

CHARACTER (LEN=255) :: FILENA
INTEGER             :: LEN, LEF
CHARACTER (LEN=11)  :: FORMAH*11

! ---------------------------------------------------------------------------- !
!
!     1. CONSTRUCT FULL FILE NAME.
!        -------------------------

FILENA = ' '
LEF   = LEN_TRIM(FILEID)
LEN = 1
IF (LEF.NE.0) THEN
   FILENA(LEN:LEN+LEF-1) = FILEID(1:LEF)
   LEN = LEN+LEF
END IF
FILENA(LEN:LEN+13) = CDATE
LEN = LEN+13
FORMAH = 'UNFORMATTED'
IF (PRESENT(FORMA)) FORMAH = FORMA

! ---------------------------------------------------------------------------- !
!
!     2. OPEN THE FILE.
!        --------------

IFAIL = 0
OPEN (UNIT=IUNIT, FILE=FILENA(1:LEN), FORM=FORMAH, STATUS=STAT, IOSTAT=IFAIL)
IF (IFAIL.NE.0) THEN
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) ' +                                           +'
   WRITE(IU06,*) ' +      WARNING ERROR IN --OPEN_FILE--       +'
   WRITE(IU06,*) ' +      ==============================       +'
   WRITE(IU06,*) ' +                                           +'
   WRITE(IU06,*) ' + COULD NOT OPEN FILE                       +'
   WRITE(IU06,*) ' + FILENAME                         : ', FILENA
   WRITE(IU06,*) ' + ERROR CODE FROM OPEN IS IOSTAT = : ', IFAIL
   WRITE(IU06,*) ' +                                           +'
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++'
ELSE   

   WRITE(IU06,*) ' SUB. OPEN_FILE: A FILE WAS CONNECTED TO UNIT =', IUNIT,     &
&                ' FILE NAME IS: ',FILENA(1:LEN)
END IF   

END SUBROUTINE OPEN_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_ARRAY_C (IUOUT, CDATE, TITL, ARRAY,                           &
                          AMOWEP, AMOSOP, AMOEAP, AMONOP, NG_R)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     PRINT_ARRAY - FORMATED OUTPUT OF AN ARRAY. (CHARACTER VERSION)           !
!                                                                              !
!     H. GUNTHER       ECMWF    NOVEMBER 1989                                  !
!     H. GUNTHER       GKSS     FEBRARY  2002   CHANGED TO FT90                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       FORMATED OUTPUT OF AN ARRAY.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       A TWO DIMENSIONAL ARRAY IS PRINTED WITH A MAXIMUM OF                   !
!       NPTS COLUMNS PER PAGE.                                                 !
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

INTEGER,            INTENT(IN)  :: IUOUT          !! OUTPUT UNIT.
CHARACTER (LEN=14), INTENT(IN)  :: CDATE          !! DATE (YYYYMMDDHHMMSS).
CHARACTER (LEN=*) , INTENT(IN)  :: TITL           !! HEADER TO BE PRINTED.
CHARACTER (LEN=1) , INTENT(IN)  :: ARRAY(:,:)     !! ARRAY TO BE PRINTED.
REAL,               INTENT(IN)  :: AMOWEP         !! WEST LONGITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMOSOP         !! SOUTH LATITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMOEAP         !! EAST LONGITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMONOP         !! NORTH LATITUDE (DEGREE).
INTEGER, OPTIONAL,  INTENT(IN)  :: NG_R(:)        !! NO. OF POINTS

! ---------------------------------------------------------------------------- !
!
!*    LOCAL VARIABLES.
!     ----------------


INTEGER    :: NGX                     !! FIRST DIMENSION  USED.
INTEGER    :: NGY                     !! SECOND DIMENSION USED.
INTEGER    :: I, NPAGE, IA, IE, L, LEN, LEN1
REAL       :: DLAMA

INTEGER, PARAMETER :: NPTS = 120

! ---------------------------------------------------------------------------- !
!
!*    1. COMPUTE LATITUDES AND LONGITUDES.
!         --------------------------------

NGX = SIZE(ARRAY,1)
NGY = SIZE(ARRAY,2)
IF (NGX.GT.1) THEN
   DLAMA = (AMOEAP-AMOWEP)/REAL(NGX-1)
ELSE
   DLAMA = 0.
END IF

! ---------------------------------------------------------------------------- !
!
!*    2. PRINT ARRAY.
!        ------------

NPAGE = (NGX+NPTS-1)/NPTS
LEN = LEN_TRIM(TITL)
LEN1 = LEN_TRIM(CDATE)

DO L = 1,NPAGE
   IA = (L-1)*NPTS
   IE = MIN(IA+NPTS,NGX)
   IA = IA + 1
   WRITE (IUOUT,'(''1'',4X,A,2X,A,5X,''PAGE '',I2,/)') CDATE(1:LEN1),          &
&                                                      TITL(1:LEN), L
   WRITE (IUOUT,'(2X,''LONGITUDE IS FROM '',F7.2,'' TO '',F7.2)')              &
&                        AMOWEP+REAL(IA-1)*DLAMA, AMOWEP+REAL(IE-1)*DLAMA
   WRITE (IUOUT,'(2X,''LATITUDE  IS FROM '',F7.2,'' TO '',F7.2)') AMONOP, AMOSOP
   WRITE (IUOUT,*) ' '
   WRITE (IUOUT,'(2X,130I1)') (MOD(I,10),I=IA,IE)
   DO I = NGY,1,-1
      IF (PRESENT(NG_R)) THEN
         IF (NG_R(I).GE.IE) THEN
            WRITE (IUOUT,'(1X,I1,130A1)') MOD(I,10),ARRAY(IA:IE,I)
         ELSE IF (NG_R(I).GE.IA) THEN
            WRITE (IUOUT,'(1X,I1,130A1)') MOD(I,10),ARRAY(IA:NG_R(I),I)
         END IF
      ELSE
         WRITE (IUOUT,'(1X,I1,130A1)') MOD(I,10),ARRAY(IA:IE,I)
      END IF
   END DO
   WRITE (IUOUT,'(2X,130I1)') (MOD(I,10),I=IA,IE)
END DO

END SUBROUTINE PRINT_ARRAY_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_ARRAY_CC (IUOUT, CDATE, TITL, ARRAY,                          &
                          AMOWEP, AMOSOP, AMOEAP, AMONOP, NG_R)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     PRINT_ARRAY - FORMATED OUTPUT OF AN ARRAY. (CHARACTER VERSION)           !
!                                                                              !
!     H. GUNTHER       ECMWF    NOVEMBER 1989                                  !
!     H. GUNTHER       GKSS     FEBRARY  2002   CHANGED TO FT90                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       FORMATED OUTPUT OF AN ARRAY.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       A TWO DIMENSIONAL ARRAY IS PRINTED WITH A MAXIMUM OF                   !
!       NPTS COLUMNS PER PAGE.                                                 !
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

INTEGER,            INTENT(IN)  :: IUOUT          !! OUTPUT UNIT.
CHARACTER (LEN=14), INTENT(IN)  :: CDATE          !! DATE (YYYYMMDDHHMMSS).
CHARACTER (LEN=*) , INTENT(IN)  :: TITL           !! HEADER TO BE PRINTED.
CHARACTER (LEN=1) , INTENT(IN)  :: ARRAY(:,:)     !! ARRAY TO BE PRINTED.
INTEGER,            INTENT(IN)  :: AMOWEP         !! WEST LONGITUDE (DEGREE).
INTEGER,            INTENT(IN)  :: AMOSOP         !! SOUTH LATITUDE (DEGREE).
INTEGER,            INTENT(IN)  :: AMOEAP         !! EAST LONGITUDE (DEGREE).
INTEGER,            INTENT(IN)  :: AMONOP         !! NORTH LATITUDE (DEGREE).
INTEGER, OPTIONAL,  INTENT(IN)  :: NG_R(:)        !! NO. OF POINTS

! ---------------------------------------------------------------------------- !
!
!*    LOCAL VARIABLES.
!     ----------------


INTEGER    :: NGX                     !! FIRST DIMENSION  USED.
INTEGER    :: NGY                     !! SECOND DIMENSION USED.
INTEGER    :: I, NPAGE, IA, IE, L, LEN, LEN1
INTEGER    :: DLAMA
character (len=len_coor) :: formtext1, formtext2

INTEGER, PARAMETER :: NPTS = 120

! ---------------------------------------------------------------------------- !
!
!*    1. COMPUTE LATITUDES AND LONGITUDES.
!         --------------------------------

NGX = SIZE(ARRAY,1)
NGY = SIZE(ARRAY,2)
IF (NGX.GT.1) THEN
   DLAMA = (AMOEAP-AMOWEP)/(NGX-1)
ELSE
   DLAMA = 0.
END IF

! ---------------------------------------------------------------------------- !
!
!*    2. PRINT ARRAY.
!        ------------

NPAGE = (NGX+NPTS-1)/NPTS
LEN = LEN_TRIM(TITL)
LEN1 = LEN_TRIM(CDATE)

DO L = 1,NPAGE
   IA = (L-1)*NPTS
   IE = MIN(IA+NPTS,NGX)
   IA = IA + 1
   WRITE (IUOUT,'(''1'',4X,A,2X,A,5X,''PAGE '',I2,/)') CDATE(1:LEN1),          &
&                                                      TITL(1:LEN), L
   formtext1 = write_coor_text (amowep+(ia-1)*dlama)
   formtext2 = write_coor_text (amowep+(ie-1)*dlama)
   WRITE (IUOUT,'(2X,''LONGITUDE IS FROM '',A,'' TO '',A)') formtext1,formtext2
   formtext1 = write_coor_text (amonop)
   formtext2 = write_coor_text (amosop)
   WRITE (IUOUT,'(2X,''LATITUDE  IS FROM '',A,'' TO '',A)') formtext1,formtext2
   WRITE (IUOUT,*) ' '
   WRITE (IUOUT,'(2X,130I1)') (MOD(I,10),I=IA,IE)
   DO I = NGY,1,-1
      IF (PRESENT(NG_R)) THEN
         IF (NG_R(I).GE.IE) THEN
            WRITE (IUOUT,'(1X,I1,130A1)') MOD(I,10),ARRAY(IA:IE,I)
         ELSE IF (NG_R(I).GE.IA) THEN
            WRITE (IUOUT,'(1X,I1,130A1)') MOD(I,10),ARRAY(IA:NG_R(I),I)
         END IF
      ELSE
         WRITE (IUOUT,'(1X,I1,130A1)') MOD(I,10),ARRAY(IA:IE,I)
      END IF
   END DO
   WRITE (IUOUT,'(2X,130I1)') (MOD(I,10),I=IA,IE)
END DO

END SUBROUTINE PRINT_ARRAY_CC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_ARRAY_L (IUOUT, CDATE, TITL, ARRAY,                           &
                          AMOWEP, AMOSOP, AMOEAP, AMONOP, NG_R)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     PRINT_ARRAY - FORMATED OUTPUT OF AN ARRAY. (LOGICAL VERSION)             !
!                                                                              !
!     H. GUNTHER       ECMWF    NOVEMBER 1989                                  !
!     H. GUNTHER       GKSS     FEBRARY  2002   CHANGED TO FT90                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       FORMATED OUTPUT OF AN ARRAY.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       A TWO DIMENSIONAL ARRAY IS PRINTED WITH A MAXIMUM OF                   !
!       NPTS COLUMNS PER PAGE.                                                 !
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

INTEGER,            INTENT(IN)  :: IUOUT          !! OUTPUT UNIT.
CHARACTER (LEN=14), INTENT(IN)  :: CDATE          !! DATE (YYYYMMDDHHMMSS).
CHARACTER (LEN=*) , INTENT(IN)  :: TITL           !! HEADER TO BE PRINTED.
LOGICAL,            INTENT(IN)  :: ARRAY(:,:)     !! ARRAY TO BE PRINTED.
REAL,               INTENT(IN)  :: AMOWEP         !! WEST LONGITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMOSOP         !! SOUTH LATITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMOEAP         !! EAST LONGITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMONOP         !! NORTH LATITUDE (DEGREE).
INTEGER, OPTIONAL,  INTENT(IN)  :: NG_R(:)        !! NO. OF POINTS

! ---------------------------------------------------------------------------- !
!
!*    LOCAL VARIABLES.
!     ----------------


INTEGER    :: NGX                     !! FIRST DIMENSION  USED.
INTEGER    :: NGY                     !! SECOND DIMENSION USED.
INTEGER    :: I, NPAGE, IA, IE, L, LEN, LEN1
REAL       :: DLAMA

INTEGER, PARAMETER :: NPTS = 120

! ---------------------------------------------------------------------------- !
!
!*    1. COMPUTE LATITUDES AND LONGITUDES.
!         --------------------------------

NGX = SIZE(ARRAY,1)
NGY = SIZE(ARRAY,2)
IF (NGX.GT.1) THEN
   DLAMA = (AMOEAP-AMOWEP)/REAL(NGX-1)
ELSE
   DLAMA = 0.
END IF

! ---------------------------------------------------------------------------- !
!
!*    2. PRINT ARRAY.
!        ------------

NPAGE = (NGX+NPTS-1)/NPTS
LEN = LEN_TRIM(TITL)
LEN1 = LEN_TRIM(CDATE)

DO L = 1,NPAGE
   IA = (L-1)*NPTS
   IE = MIN(IA+NPTS,NGX)
   IA = IA + 1
   WRITE (IUOUT,'(''1'',4X,A,2X,A,5X,''PAGE '',I2,/)') CDATE(1:LEN1),          &
&                                                      TITL(1:LEN), L
   WRITE (IUOUT,'(2X,''LONGITUDE IS FROM '',F7.2,'' TO '',F7.2)')              &
&                        AMOWEP+REAL(IA-1)*DLAMA, AMOWEP+REAL(IE-1)*DLAMA
   WRITE (IUOUT,'(2X,''LATITUDE  IS FROM '',F7.2,'' TO '',F7.2)') AMONOP, AMOSOP
   WRITE (IUOUT,*) ' '
   WRITE (IUOUT,'(2X,130I1)') (MOD(I,10),I=IA,IE)
   DO I = NGY,1,-1
      IF (PRESENT(NG_R)) THEN
         IF (NG_R(I).GE.IE) THEN
            WRITE (IUOUT,'(1X,I1,130L1)') MOD(I,10),ARRAY(IA:IE,I)
         ELSE IF (NG_R(I).GE.IA) THEN
            WRITE (IUOUT,'(1X,I1,130L1)') MOD(I,10),ARRAY(IA:NG_R(I),I)
         END IF
      ELSE
         WRITE (IUOUT,'(1X,I1,130L1)') MOD(I,10),ARRAY(IA:IE,I)
      END IF
   END DO
   WRITE (IUOUT,'(2X,130I1)') (MOD(I,10),I=IA,IE)
END DO

END SUBROUTINE PRINT_ARRAY_L

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_ARRAY_R (IUOUT, CDATE, TITL, ARRAY,                           &
                          AMOWEP, AMOSOP, AMOEAP, AMONOP, SCALE, ZMISS, NG_R)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     PRINT_ARRAY - FORMATTED OUTPUT OF AN ARRAY.                              !
!                                                                              !
!     H. GUNTHER       ECMWF    NOVEMBER 1989                                  !
!     H. GUNTHER       GKSS     FEBRARY  2002   CHANGED TO FT90                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       FORMATED OUTPUT OF AN ARRAY.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       A TWO DIMENSIONAL ARRAY IS PRINTED WITH A MAXIMUM OF                   !
!       NPTS COLUMNS PER PAGE.                                                 !
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

INTEGER,            INTENT(IN)  :: IUOUT          !! OUTPUT UNIT.
CHARACTER (LEN=14), INTENT(IN)  :: CDATE          !! DATE (YYYYMMDDHHMMSS).
CHARACTER (LEN=*) , INTENT(IN)  :: TITL           !! HEADER TO BE PRINTED.
REAL,               INTENT(IN)  :: ARRAY(:,:)     !! ARRAY TO BE PRINTED.
REAL,               INTENT(IN)  :: AMOWEP         !! WEST LONGITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMOSOP         !! SOUTH LATITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMOEAP         !! EAST LONGITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMONOP         !! NORTH LATITUDE (DEGREE).
REAL, OPTIONAL,     INTENT(IN)  :: SCALE          !! SCALING FACTOR.
REAL, OPTIONAL,     INTENT(IN)  :: ZMISS          !! MISSING VALUE.
INTEGER, OPTIONAL,  INTENT(IN)  :: NG_R(:)        !! NO. OF POINTS

! ---------------------------------------------------------------------------- !
!
!     LOCAL VARIABLES.
!     ----------------


INTEGER    :: NGX                     !! FIRST DIMENSION  USED.
INTEGER    :: NGY                     !! SECOND DIMENSION USED.
INTEGER    :: IARRAY(SIZE(ARRAY,1),SIZE(ARRAY,2)), ILON(SIZE(ARRAY,1))
REAL       :: YLAT(SIZE(ARRAY,2))
INTEGER    :: I, J, NPAGE, ISTART, IEND, NP, LEN, LEN1
REAL       :: DLAMA, DPHIA, SCALEH

INTEGER, PARAMETER :: NPTS = 30

! ---------------------------------------------------------------------------- !
!
!      1. COMPUTE LATITUDES AND LONGITUDES.
!         --------------------------------

NGX = SIZE(ARRAY,1)
NGY = SIZE(ARRAY,2)

IF (NGX.GT.1) THEN
   DLAMA = (AMOEAP-AMOWEP)/REAL(NGX-1)
ELSE
   DLAMA = 0.
END IF
DO I = 1,NGX
   ILON(I) = NINT(AMOWEP + (I-1)*DLAMA)
END DO
IF (NGX.GT.1) THEN
   DPHIA = (AMONOP-AMOSOP)/REAL(NGY-1)
ELSE
   DPHIA = 0.
END IF
DO J = 1,NGY
   YLAT(J) = AMOSOP + REAL(J-1)*DPHIA
END DO

! ---------------------------------------------------------------------------- !
!
!     2. SCALE DATA ARRAY.
!        -----------------

SCALEH = 1.
IF (PRESENT(SCALE)) THEN
   IF (SCALE.EQ.0.) THEN
      IF (PRESENT(ZMISS)) THEN
         SCALEH = MAXVAL(ABS(ARRAY),MASK=ARRAY.GT.ZMISS)
      ELSE
         SCALEH = MAXVAL(ABS(ARRAY)) 
      END IF
      IF (SCALEH.NE.0) THEN
         SCALEH = 100./(10.**NINT(LOG10(SCALEH)))
      ELSE
         SCALEH = 1.
      END IF
   ELSE
      SCALEH = SCALE
   END IF

   IF (PRESENT(ZMISS)) THEN
      WHERE (ARRAY.GT.ZMISS)
         IARRAY = NINT(SCALEH*ARRAY)
      ELSEWHERE
         IARRAY = NINT(ARRAY)
      END WHERE
   ELSE
      IARRAY = NINT(SCALEH*ARRAY)
   END IF
ELSE
   IARRAY = NINT(ARRAY)
END IF

! ---------------------------------------------------------------------------- !
!
!     3. PRINT ARRAY.
!        ------------

NPAGE = (NGX+NPTS-1)/NPTS
LEN  = LEN_TRIM(TITL)
LEN1 = LEN_TRIM(CDATE)

ISTART = -NPTS+1
IEND   = ISTART+NPTS-1
DO NP = 1,NPAGE
   IF (ABS(SCALEH).GT. 99999. .OR. ABS(SCALEH).LT. 0.01) THEN
      WRITE (IUOUT,'(''1'',4X,A,2X,2A,E10.3,5X,''PAGE '',I2,/)')               &
&                   CDATE(1:LEN1), TITL(1:LEN),' MULTIPLIED BY: ', SCALEH, NP
   ELSE
      WRITE (IUOUT,'(''1'',4X,A,2X,2A,F10.2,5X,''PAGE '',I2,/)')               &
&                   CDATE(1:LEN1), TITL(1:LEN),' MULTIPLIED BY: ', SCALEH, NP
   END IF

   ISTART = ISTART+NPTS
   IEND   = MIN(IEND+NPTS,NGX)
   WRITE (IUOUT,'(4X,30I4)') (I,I=ISTART,IEND)
    DO J = NGY,1,-1
      IF (PRESENT(NG_R)) THEN
         IF (NG_R(J).GE.IEND) THEN
            WRITE (IUOUT,'(I4,30I4)') J,IARRAY(ISTART:IEND,J)
         ELSE IF (NG_R(J).GE.ISTART) THEN
            WRITE (IUOUT,'(I4,30I4)') J,IARRAY(ISTART:NG_R(J),J)
         END IF
      ELSE
         WRITE (IUOUT,'(I4,30I4)') J,IARRAY(ISTART:IEND,J)
      END IF
   END DO
END DO

END SUBROUTINE PRINT_ARRAY_R

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_ARRAY_RC (IUOUT, CDATE, TITL, ARRAY,                          &
                          AMOWEP, AMOSOP, AMOEAP, AMONOP, SCALE, ZMISS, NG_R)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     PRINT_ARRAY - FORMATTED OUTPUT OF AN ARRAY.                              !
!                                                                              !
!     H. GUNTHER       ECMWF    NOVEMBER 1989                                  !
!     H. GUNTHER       GKSS     FEBRARY  2002   CHANGED TO FT90                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       FORMATED OUTPUT OF AN ARRAY.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       A TWO DIMENSIONAL ARRAY IS PRINTED WITH A MAXIMUM OF                   !
!       NPTS COLUMNS PER PAGE.                                                 !
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

INTEGER,            INTENT(IN)  :: IUOUT          !! OUTPUT UNIT.
CHARACTER (LEN=14), INTENT(IN)  :: CDATE          !! DATE (YYYYMMDDHHMMSS).
CHARACTER (LEN=*) , INTENT(IN)  :: TITL           !! HEADER TO BE PRINTED.
REAL,               INTENT(IN)  :: ARRAY(:,:)     !! ARRAY TO BE PRINTED.
INTEGER,            INTENT(IN)  :: AMOWEP         !! WEST LONGITUDE (M_SEC).
INTEGER,            INTENT(IN)  :: AMOSOP         !! SOUTH LATITUDE (M_SEC).
INTEGER,            INTENT(IN)  :: AMOEAP         !! EAST LONGITUDE (M_SEC).
INTEGER,            INTENT(IN)  :: AMONOP         !! NORTH LATITUDE (M_SEC).
REAL, OPTIONAL,     INTENT(IN)  :: SCALE          !! SCALING FACTOR.
REAL, OPTIONAL,     INTENT(IN)  :: ZMISS          !! MISSING VALUE.
INTEGER, OPTIONAL,  INTENT(IN)  :: NG_R(:)        !! NO. OF POINTS

! ---------------------------------------------------------------------------- !
!
!     LOCAL VARIABLES.
!     ----------------


INTEGER    :: NGX                     !! FIRST DIMENSION  USED.
INTEGER    :: NGY                     !! SECOND DIMENSION USED.
INTEGER    :: IARRAY(SIZE(ARRAY,1),SIZE(ARRAY,2))
INTEGER    :: I, J, NPAGE, ISTART, IEND, NP, LEN, LEN1
INTEGER    :: DLAMA, DPHIA
REAL       :: SCALEH
character (len=len_coor) :: formtext1, formtext2

INTEGER, PARAMETER :: NPTS = 30

! ---------------------------------------------------------------------------- !
!
!      1. COMPUTE LATITUDES AND LONGITUDES.
!         --------------------------------

NGX = SIZE(ARRAY,1)
NGY = SIZE(ARRAY,2)

IF (NGX.GT.1) THEN
   DLAMA = (AMOEAP-AMOWEP)/(NGX-1)
ELSE
   DLAMA = 0.
END IF
IF (NGX.GT.1) THEN
   DPHIA = (AMONOP-AMOSOP)/(NGY-1)
ELSE
   DPHIA = 0.
END IF

! ---------------------------------------------------------------------------- !
!
!     2. SCALE DATA ARRAY.
!        -----------------

SCALEH = 1.
IF (PRESENT(SCALE)) THEN
   IF (SCALE.EQ.0.) THEN
      IF (PRESENT(ZMISS)) THEN
         SCALEH = MAXVAL(ABS(ARRAY),MASK=ARRAY.GT.ZMISS)
      ELSE
         SCALEH = MAXVAL(ABS(ARRAY))
     END IF

      IF (SCALEH.NE.0) THEN
         SCALEH = 100./(10.**NINT(LOG10(SCALEH)))
      ELSE
         SCALEH = 1.
      END IF
   ELSE
      SCALEH = SCALE
   END IF

   IF (PRESENT(ZMISS)) THEN
      WHERE (ARRAY.GT.ZMISS)
         IARRAY = NINT(SCALEH*ARRAY)
      ELSEWHERE
         IARRAY = NINT(ARRAY)
      END WHERE
   ELSE
      IARRAY = NINT(SCALEH*ARRAY)
   END IF

ELSE
   IARRAY = NINT(ARRAY)
END IF

! ---------------------------------------------------------------------------- !
!
!     3. PRINT ARRAY.
!        ------------

NPAGE = (NGX+NPTS-1)/NPTS
LEN  = LEN_TRIM(TITL)
LEN1 = LEN_TRIM(CDATE)

ISTART = -NPTS+1
IEND   = ISTART+NPTS-1
DO NP = 1,NPAGE
   IF (ABS(SCALEH).GT. 99999. .OR. ABS(SCALEH).LT. 0.01) THEN
      WRITE (IUOUT,'(''1'',4X,A,2X,2A,E10.3,5X,''PAGE '',I2,/)')               &
&                   CDATE(1:LEN1), TITL(1:LEN),' MULTIPLIED BY: ', SCALEH, NP
   ELSE
      WRITE (IUOUT,'(''1'',4X,A,2X,2A,F10.2,5X,''PAGE '',I2,/)')               &
&                   CDATE(1:LEN1), TITL(1:LEN),' MULTIPLIED BY: ', SCALEH, NP
   END IF

   ISTART = ISTART+NPTS
   IEND   = MIN(IEND+NPTS,NGX)
   formtext1 = write_coor_text (amowep+(istart-1)*dlama)
   formtext2 = write_coor_text (amowep+(iend-1)*dlama)
   WRITE (IUOUT,'(2X,''LONGITUDE IS FROM '',A,'' TO '',A)') formtext1,formtext2
   formtext1 = write_coor_text (amonop)
   formtext2 = write_coor_text (amosop)
   WRITE (IUOUT,'(2X,''LATITUDE  IS FROM '',A,'' TO '',A)') formtext1,formtext2

   WRITE (IUOUT,'(4X,30I4)') (I,I=ISTART,IEND)
   DO J = NGY,1,-1
      IF (PRESENT(NG_R)) THEN
         IF (NG_R(J).GE.IEND) THEN
            WRITE (IUOUT,'(I4,30I4)') J,IARRAY(ISTART:IEND,J)
         ELSE IF (NG_R(J).GE.ISTART) THEN
            WRITE (IUOUT,'(I4,30I4)') J,IARRAY(ISTART:NG_R(J),J)
         END IF
      ELSE
         WRITE (IUOUT,'(I4,30I4)') J,IARRAY(ISTART:IEND,J)
      END IF
   END DO
END DO

END SUBROUTINE PRINT_ARRAY_RC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_ARRAY_I (IUOUT, CDATE, TITL, ARRAY,                           &
                          AMOWEP, AMOSOP, AMOEAP, AMONOP, NG_R)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     PRINT_ARRAY - FORMATTED OUTPUT OF AN ARRAY.                              !
!                                                                              !
!     H. GUNTHER       ECMWF    NOVEMBER 1989                                  !
!     H. GUNTHER       GKSS     FEBRARY  2002   CHANGED TO FT90                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       FORMATED OUTPUT OF AN ARRAY.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       A TWO DIMENSIONAL ARRAY IS PRINTED WITH A MAXIMUM OF                   !
!       NPTS COLUMNS PER PAGE.                                                 !
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

INTEGER,            INTENT(IN)  :: IUOUT          !! OUTPUT UNIT.
CHARACTER (LEN=14), INTENT(IN)  :: CDATE          !! DATE (YYYYMMDDHHMMSS).
CHARACTER (LEN=*) , INTENT(IN)  :: TITL           !! HEADER TO BE PRINTED.
INTEGER,            INTENT(IN)  :: ARRAY(:,:)     !! ARRAY TO BE PRINTED.
REAL,               INTENT(IN)  :: AMOWEP         !! WEST LONGITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMOSOP         !! SOUTH LATITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMOEAP         !! EAST LONGITUDE (DEGREE).
REAL,               INTENT(IN)  :: AMONOP         !! NORTH LATITUDE (DEGREE).
INTEGER, OPTIONAL,  INTENT(IN)  :: NG_R(:)        !! NO. OF POINTS

! ---------------------------------------------------------------------------- !
!
!     LOCAL VARIABLES.
!     ----------------


INTEGER    :: NGX                     !! FIRST DIMENSION  USED.
INTEGER    :: NGY                     !! SECOND DIMENSION USED.
INTEGER    :: I, J, NPAGE, ISTART, IEND, NP, LEN, LEN1
REAL       :: DLAMA

INTEGER, PARAMETER :: NPTS = 30

! ---------------------------------------------------------------------------- !
!
!      1. COMPUTE LATITUDES AND LONGITUDES.
!         --------------------------------

NGX = SIZE(ARRAY,1)
NGY = SIZE(ARRAY,2)

IF (NGX.GT.1) THEN
   DLAMA = (AMOEAP-AMOWEP)/REAL(NGX-1)
ELSE
   DLAMA = 0.
END IF

! ---------------------------------------------------------------------------- !
!
!     2. PRINT ARRAY.
!        ------------

NPAGE = (NGX+NPTS-1)/NPTS
LEN  = LEN_TRIM(TITL)
LEN1 = LEN_TRIM(CDATE)

ISTART = -NPTS+1
IEND   = ISTART+NPTS-1
DO NP = 1,NPAGE
   WRITE (IUOUT,'(''1'',4X,A,2X,A,5X,''PAGE '',I2,/)') CDATE(1:LEN1),          &
&                                  TITL(1:LEN), NP
   ISTART = ISTART+NPTS
   IEND   = MIN(IEND+NPTS,NGX)
   WRITE (IUOUT,'(2X,''LONGITUDE IS FROM '',F8.4,'' TO '',F8.4)') amowep+(istart-1)*dlama, amowep+(iend-1)*dlama
   WRITE (IUOUT,'(2X,''LATITUDE  IS FROM '',F8.4,'' TO '',F8.4)') AMONOP, AMOSOP
   WRITE (IUOUT,'(7X,''I='',30I4)') (I,I=ISTART,IEND)
   DO J = NGY,1,-1
      IF (PRESENT(NG_R)) THEN
         IF (NG_R(J).GE.IEND) THEN
            WRITE (IUOUT,'(I4,30I4)') MOD(I,10),ARRAY(ISTART:IEND,J)
         ELSE IF (NG_R(J).GE.ISTART) THEN
            WRITE (IUOUT,'(I4,30I4)') MOD(I,10),ARRAY(ISTART:NG_R(J),J)
         END IF
      ELSE
         WRITE (IUOUT,'(I4,30I4)') MOD(I,10),ARRAY(ISTART:IEND,J)
      END IF
   END DO
END DO

END SUBROUTINE PRINT_ARRAY_I

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_ARRAY_IC (IUOUT, CDATE, TITL, ARRAY,                           &
                          AMOWEP, AMOSOP, AMOEAP, AMONOP, NG_R)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     PRINT_ARRAY - FORMATTED OUTPUT OF AN ARRAY.                              !
!                                                                              !
!     H. GUNTHER       ECMWF    NOVEMBER 1989                                  !
!     H. GUNTHER       GKSS     FEBRARY  2002   CHANGED TO FT90                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       FORMATED OUTPUT OF AN ARRAY.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       A TWO DIMENSIONAL ARRAY IS PRINTED WITH A MAXIMUM OF                   !
!       NPTS COLUMNS PER PAGE.                                                 !
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

INTEGER,            INTENT(IN)  :: IUOUT          !! OUTPUT UNIT.
CHARACTER (LEN=14), INTENT(IN)  :: CDATE          !! DATE (YYYYMMDDHHMMSS).
CHARACTER (LEN=*) , INTENT(IN)  :: TITL           !! HEADER TO BE PRINTED.
INTEGER,            INTENT(IN)  :: ARRAY(:,:)     !! ARRAY TO BE PRINTED.
INTEGER,            INTENT(IN)  :: AMOWEP         !! WEST LONGITUDE (M_SEC).
INTEGER,            INTENT(IN)  :: AMOSOP         !! SOUTH LATITUDE (M_SEC).
INTEGER,            INTENT(IN)  :: AMOEAP         !! EAST LONGITUDE (M_SEC).
INTEGER,            INTENT(IN)  :: AMONOP         !! NORTH LATITUDE (M_SEC).
INTEGER, OPTIONAL,  INTENT(IN)  :: NG_R(:)        !! NO. OF POINTS

! ---------------------------------------------------------------------------- !
!
!     LOCAL VARIABLES.
!     ----------------


INTEGER    :: NGX                     !! FIRST DIMENSION  USED.
INTEGER    :: NGY                     !! SECOND DIMENSION USED.
INTEGER    :: I,J, NPAGE, ISTART, IEND, NP, LEN, LEN1
INTEGER    :: DLAMA
character (len=len_coor) :: ftext1, ftext2

INTEGER, PARAMETER :: NPTS = 30

! ---------------------------------------------------------------------------- !
!
!      1. COMPUTE LATITUDES AND LONGITUDES.
!         --------------------------------

NGX = SIZE(ARRAY,1)
NGY = SIZE(ARRAY,2)

IF (NGX.GT.1) THEN
   DLAMA = (AMOEAP-AMOWEP)/(NGX-1)
ELSE
   DLAMA = 0.
END IF


! ---------------------------------------------------------------------------- !
!
!     2. PRINT ARRAY.
!        ------------

NPAGE = (NGX+NPTS-1)/NPTS
LEN  = LEN_TRIM(TITL)
LEN1 = LEN_TRIM(CDATE)

ISTART = -NPTS+1
IEND   = ISTART+NPTS-1
DO NP = 1,NPAGE
   WRITE (IUOUT,'(''1'',4X,A,2X,A,5X,''PAGE '',I2,/)') CDATE(1:LEN1),          &
&                                  TITL(1:LEN), NP
   ISTART = ISTART+NPTS
   IEND   = MIN(IEND+NPTS,NGX)

   ftext1 = write_coor_text (amowep+(istart-1)*dlama)
   ftext2 = write_coor_text (amowep+(iend-1)*dlama)
   WRITE (IUOUT,'(2X,''LONGITUDE IS FROM '',A,'' TO '',A)') ftext1, ftext2

   ftext1 = write_coor_text (amonop)
   ftext2 = write_coor_text (amosop)
   WRITE (IUOUT,'(2X,''LATITUDE  IS FROM '',A,'' TO '',A)') ftext1, ftext2

   WRITE (IUOUT,'(4X,30I4)') (I,I=ISTART,IEND)
   DO J = NGY,1,-1
      IF (PRESENT(NG_R)) THEN
         IF (NG_R(J).GE.IEND) THEN
            WRITE (IUOUT,'(I4,30I4)') J,ARRAY(ISTART:IEND,J)
         ELSE IF (NG_R(J).GE.ISTART) THEN
            WRITE (IUOUT,'(I4,30I4)') J,ARRAY(ISTART:NG_R(J),J)
         END IF
      ELSE
         WRITE (IUOUT,'(I4,30I4)') J,ARRAY(ISTART:IEND,J)
      END IF
   END DO
END DO

END SUBROUTINE PRINT_ARRAY_IC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_SPECTRUM (IUNIT, CDATE, LONG, LAT, TITL, FR, THETA, SPEC,     &
&                          U10, UDIR, US, DEPTH, CURR_SP, CURR_DIR,            &
&                          HS, PPER, MPER, TM1, TM2, MDIR, SPRE)

! ---------------------------------------------------------------------------- !
! 
!   PRINT_SPECTRUM -  PRINTS A SPECTRUM.
! 
!         M. DE LAS HERAS  KNMI/PCM  FEBRUARY  1990
!
!     PURPOSE.
!     --------
!
!        PRINT A WAVE MODEL SPECTRUM.
!
!     METHOD.
!     -------
!
!       NONE.
!
!     REFERENCE.
!     ----------
!
!       NONE.
!
! ---------------------------------------------------------------------------- !
!
!     INTERFACE VARIABLES.
!     --------------------

INTEGER,            INTENT(IN) :: IUNIT     !! OUTPUT UNIT.
CHARACTER (LEN=14), INTENT(IN) :: CDATE     !! DATE OF SPECTRUM (YYYYMMDDHHMMSS).
REAL,               INTENT(IN) :: LONG      !! LONGITUDE OF SPECTRUM (DEGREE).
REAL,               INTENT(IN) :: LAT       !! LATITUDE OF SPECTRUM (DEGREE).
CHARACTER (LEN=40), INTENT(IN) :: TITL      !! TITLE.
REAL,               INTENT(IN) :: SPEC(:,:) !! SPECTRUM.
REAL,               INTENT(IN) :: FR(:)     !! FREQUENCY ARRAY IN HERTZ.
REAL,               INTENT(IN) :: THETA(:)  !! DIRECTION ARRAY IN RAD.
REAL,               INTENT(IN) :: U10       !! WIND SPEED U10 (METRES/SECOND).
REAL,               INTENT(IN) :: UDIR      !! WIND DIRECTION (DEGREES).
REAL,               INTENT(IN) :: US        !! FRICTION VELOCITY (METRES/SECOND).
REAL,               INTENT(IN) :: DEPTH     !! DEPTH (METRES).
REAL,               INTENT(IN) :: CURR_SP   !! CURRENT SPEED (METRES/SECONDS).
REAL,               INTENT(IN) :: CURR_DIR  !! CURRENT DIRECTION (DEGREE).
REAL,               INTENT(IN) :: HS        !! SIGNIFICANT WAVE HEIGHT (METRES).
REAL,               INTENT(IN) :: PPER      !! PEAK PERIOD (S).
REAL,               INTENT(IN) :: MPER      !! MEAN WAVE PERIOD (S).
REAL,               INTENT(IN) :: TM1       !! TM1 PERIOD (S).
REAL,               INTENT(IN) :: TM2       !! TM2 PERIOD (S).
REAL,               INTENT(IN) :: MDIR      !! MEAN WAVE DIRECTION (DEGREE).
REAL,               INTENT(IN) :: SPRE      !! DIRECTIONAL SPREAD (DEGREE).

! ---------------------------------------------------------------------------- !
! 
!     LOCAL VARIABLES.
!     ----------------

INTEGER, PARAMETER :: IPDIR = 12   !! NUMBER OF DIRECTIONS PRINTED PER LINE.

INTEGER            :: KL              !! NUMBER OF DIRECTIONS.
INTEGER            :: ML              !! NUMBER OF FREQUENCIES.
INTEGER            :: IPE, IP, M, LEN
REAL               :: DELTH
REAL               :: ANG(SIZE(THETA))     !! DIRECTIONS IN DEGREE.
REAL               :: ODSPEC(SIZE(FR))     !! 1-D SPECTRUM.  (M*M/HERTZ)
REAL               :: CDIR(SIZE(FR))       !! 1-D SPECTRUM.  (M*M/HERTZ)
REAL               :: SDIR(SIZE(FR))       !! 1-D SPECTRUM.  (M*M/HERTZ)
CHARACTER (LEN=60) :: FORM1, FORM2, FORM3  !! VARIABLE FORMATS.

! ---------------------------------------------------------------------------- !
! 
!     1. INITIALISE DIRECTIONS.
!        ----------------------

KL = SIZE(THETA)
ML = SIZE(FR)
DELTH = ZPI/REAL(KL)
ANG = THETA*DEG

! ---------------------------------------------------------------------------- !
!
!     2. COMPUTE 1-D SPECTRUM.
!        ---------------------

ODSPEC = SUM(SPEC, DIM=1)
ODSPEC = ODSPEC*DELTH
CDIR = MATMUL(COS(THETA),SPEC)
SDIR = MATMUL(SIN(THETA),SPEC)
CDIR = ATAN2(SDIR,CDIR)*DEG
WHERE (CDIR.LT.0.) CDIR = CDIR+360.
WHERE (CDIR.GT.359.95) CDIR = 0.

! ---------------------------------------------------------------------------- !
!
!*    3. PRINT SPECTRUM.
!        ---------------

WRITE(IUNIT,'(''1'',A40,'' DATE: '',A14,     &
&             ''   LONG.: '',F7.2,'' LAT.: '',F6.2)') TITL, CDATE, LONG, LAT
WRITE(IUNIT,'('' U10 = '',F5.2,''M/S  UDIR = '',F5.0,''DEG  USTAR = '',F5.2,   &
&             ''M/S'')') U10, UDIR, US
WRITE(IUNIT,'('' DEPTH = '',F6.1,''M    CURR_SPEED = '',F5.2,                  &
&             ''M/S  CURR_DIR = '',F5.0,''DEG'')') DEPTH, CURR_SP, CURR_DIR
WRITE(IUNIT,'(''  HS = '',F5.2,''M    PPER = '',F5.2,''S     MPER = '',F5.2,   &
&             ''S  TM1 = '',F5.2,''S  TM2 = '',F5.2,''S  MDIR = '',F5.0,       &
&             ''DEG  SPREAD = '',F5.0,''DEG'')')  HS, PPER, MPER, TM1, TM2,    &
&                           MDIR, SPRE

IPE = 0
IF (KL.GT.IPDIR) THEN
   FORM1 = '(1X,''DIR (DEG)'',T11,12F8.1,''  DIR (DEG)'')'
   WRITE (FORM1(21:22),'(I2)') IPDIR
   FORM2 = '(1X,F7.4,2X,12F8.3,F9.4)'
   WRITE (FORM2(13:14),'(I2)') IPDIR
   FORM3 = '(1X,''FREQ (HZ)'',T11, 96X,''  FREQ (HZ)'')'
   WRITE (FORM3(21:23),'(I3)') IPDIR*8
   DO IP = 1,KL-IPDIR,IPDIR
      IPE = IP+IPDIR-1
      WRITE(IUNIT,FORM1) ANG(IP:IPE)
      WRITE(IUNIT,FORM3)
      WRITE(IUNIT,FORM2) (FR(M),SPEC(IP:IPE,M),FR(M),M=1,ML)
      WRITE(IUNIT,FORM1) ANG(IP:IPE)
      WRITE(IUNIT,'(1X)')
   END DO
END IF

IP = IPE+1
IPE = KL
LEN = IPE-IP+1
FORM1 = '(1X,''DIR (DEG)'',T11,12F8.1,''  DIR (DEG)'')'
WRITE (FORM1(21:22),'(I2)') LEN
FORM2 = '(1X,F7.4,2X,12F8.3,F9.4,F10.3,F6.0)'
WRITE (FORM2(13:14),'(I2)') LEN
FORM3 = '(1X,''FREQ (HZ)'',T11, 96X,''  FREQ (HZ) 1-DSPEC M-DIR'')'
WRITE (FORM3(21:23),'(I3)') LEN*8
WRITE (IUNIT,FORM1) ANG(IP:IPE)
WRITE (IUNIT,FORM3)
WRITE (IUNIT,FORM2) (FR(M),SPEC(IP:IPE,M),FR(M),ODSPEC(M),CDIR(M),M=1,ML)
WRITE (IUNIT,FORM1) ANG(IP:IPE)

END SUBROUTINE PRINT_SPECTRUM

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_SPECTRUM_C (IUNIT, CDATE, LONG, LAT, TITL, FR, THETA, SPEC,   &
&                          U10, UDIR, US, DEPTH, CURR_SP, CURR_DIR,            &
&                          HS, PPER, MPER, TM1, TM2, MDIR, SPRE)

! ---------------------------------------------------------------------------- !
! 
!   PRINT_SPECTRUM -  PRINTS A SPECTRUM.
! 
!         M. DE LAS HERAS  KNMI/PCM  FEBRUARY  1990
!
!     PURPOSE.
!     --------
!
!        PRINT A WAVE MODEL SPECTRUM.
!
!     METHOD.
!     -------
!
!       NONE.
!
!     REFERENCE.
!     ----------
!
!       NONE.
!
! ---------------------------------------------------------------------------- !
!
!     INTERFACE VARIABLES.
!     --------------------

INTEGER,            INTENT(IN) :: IUNIT     !! OUTPUT UNIT.
CHARACTER (LEN=14), INTENT(IN) :: CDATE     !! DATE OF SPECTRUM (YYYYMMDDHHMMSS).
INTEGER,            INTENT(IN) :: LONG      !! LONGITUDE OF SPECTRUM (DEGREE).
INTEGER,            INTENT(IN) :: LAT       !! LATITUDE OF SPECTRUM (DEGREE).
CHARACTER (LEN=40), INTENT(IN) :: TITL      !! TITLE.
REAL,               INTENT(IN) :: SPEC(:,:) !! SPECTRUM.
REAL,               INTENT(IN) :: FR(:)     !! FREQUENCY ARRAY IN HERTZ.
REAL,               INTENT(IN) :: THETA(:)  !! DIRECTION ARRAY IN RAD.
REAL,               INTENT(IN) :: U10       !! WIND SPEED U10 (METRES/SECOND).
REAL,               INTENT(IN) :: UDIR      !! WIND DIRECTION (DEGREES).
REAL,               INTENT(IN) :: US        !! FRICTION VELOCITY (METRES/SECOND).
REAL,               INTENT(IN) :: DEPTH     !! DEPTH (METRES).
REAL,               INTENT(IN) :: CURR_SP   !! CURRENT SPEED (METRES/SECONDS).
REAL,               INTENT(IN) :: CURR_DIR  !! CURRENT DIRECTION (DEGREE).
REAL,               INTENT(IN) :: HS        !! SIGNIFICANT WAVE HEIGHT (METRES).
REAL,               INTENT(IN) :: PPER      !! PEAK PERIOD (S).
REAL,               INTENT(IN) :: MPER      !! MEAN WAVE PERIOD (S).
REAL,               INTENT(IN) :: TM1       !! TM1 PERIOD (S).
REAL,               INTENT(IN) :: TM2       !! TM2 PERIOD (S).
REAL,               INTENT(IN) :: MDIR      !! MEAN WAVE DIRECTION (DEGREE).
REAL,               INTENT(IN) :: SPRE      !! DIRECTIONAL SPREAD (DEGREE).

! ---------------------------------------------------------------------------- !
! 
!     LOCAL VARIABLES.
!     ----------------

INTEGER, PARAMETER :: IPDIR = 12   !! NUMBER OF DIRECTIONS PRINTED PER LINE.

INTEGER            :: KL              !! NUMBER OF DIRECTIONS.
INTEGER            :: ML              !! NUMBER OF FREQUENCIES.
INTEGER            :: IPE, IP, M, LEN
REAL               :: DELTH
REAL               :: ANG(SIZE(THETA))     !! DIRECTIONS IN DEGREE.
REAL               :: ODSPEC(SIZE(FR))     !! 1-D SPECTRUM.  (M*M/HERTZ)
REAL               :: CDIR(SIZE(FR))       !! 1-D SPECTRUM.  (M*M/HERTZ)
REAL               :: SDIR(SIZE(FR))       !! 1-D SPECTRUM.  (M*M/HERTZ)
CHARACTER (LEN=60) :: FORM1, FORM2, FORM3  !! VARIABLE FORMATS.
character (len=len_coor) :: ftext1, ftext2

! ---------------------------------------------------------------------------- !
! 
!     1. INITIALISE DIRECTIONS.
!        ----------------------

KL = SIZE(THETA)
ML = SIZE(FR)
DELTH = ZPI/REAL(KL)
ANG = THETA*DEG

! ---------------------------------------------------------------------------- !
!
!     2. COMPUTE 1-D SPECTRUM.
!        ---------------------

ODSPEC = SUM(SPEC, DIM=1)
ODSPEC = ODSPEC*DELTH
CDIR = MATMUL(COS(THETA),SPEC)
SDIR = MATMUL(SIN(THETA),SPEC)
CDIR = ATAN2(SDIR,CDIR)*DEG
WHERE (CDIR.LT.0.) CDIR = CDIR+360.
WHERE (CDIR.GT.359.95) CDIR = 0.

! ---------------------------------------------------------------------------- !
!
!*    3. PRINT SPECTRUM.
!        ---------------

ftext1 = write_coor_text (long)
ftext2 = write_coor_text (lat)
WRITE(IUNIT,'(''1'',A40,'' DATE: '',A14,     &
&             ''   LONG.: '',A,'' LAT.: '',A)') TITL, CDATE, ftext1, ftext2
WRITE(IUNIT,'('' U10 = '',F5.2,''M/S  UDIR = '',F5.0,''DEG  USTAR = '',F5.2,   &
&             ''M/S'')') U10, UDIR, US
WRITE(IUNIT,'('' DEPTH = '',F6.1,''M    CURR_SPEED = '',F5.2,                  &
&             ''M/S  CURR_DIR = '',F5.0,''DEG'')') DEPTH, CURR_SP, CURR_DIR
WRITE(IUNIT,'(''  HS = '',F5.2,''M    PPER = '',F5.2,''S     MPER = '',F5.2,   &
&             ''S  TM1 = '',F5.2,''S  TM2 = '',F5.2,''S  MDIR = '',F5.0,       &
&             ''DEG  SPREAD = '',F5.0,''DEG'')')  HS, PPER, MPER, TM1, TM2,    &
&                           MDIR, SPRE

IPE = 0
IF (KL.GT.IPDIR) THEN
   FORM1 = '(1X,''DIR (DEG)'',T11,12F8.1,''  DIR (DEG)'')'
   WRITE (FORM1(21:22),'(I2)') IPDIR
   FORM2 = '(1X,F7.4,2X,12F8.3,F9.4)'
   WRITE (FORM2(13:14),'(I2)') IPDIR
   FORM3 = '(1X,''FREQ (HZ)'',T11, 96X,''  FREQ (HZ)'')'
   WRITE (FORM3(21:23),'(I3)') IPDIR*8
   DO IP = 1,KL-IPDIR,IPDIR
      IPE = IP+IPDIR-1
      WRITE(IUNIT,FORM1) ANG(IP:IPE)
      WRITE(IUNIT,FORM3)
      WRITE(IUNIT,FORM2) (FR(M),SPEC(IP:IPE,M),FR(M),M=1,ML)
      WRITE(IUNIT,FORM1) ANG(IP:IPE)
      WRITE(IUNIT,'(1X)')
   END DO
END IF

IP = IPE+1
IPE = KL
LEN = IPE-IP+1
FORM1 = '(1X,''DIR (DEG)'',T11,12F8.1,''  DIR (DEG)'')'
WRITE (FORM1(21:22),'(I2)') LEN
FORM2 = '(1X,F7.4,2X,12F8.3,F9.4,F10.3,F6.0)'
WRITE (FORM2(13:14),'(I2)') LEN
FORM3 = '(1X,''FREQ (HZ)'',T11, 96X,''  FREQ (HZ) 1-DSPEC M-DIR'')'
WRITE (FORM3(21:23),'(I3)') LEN*8
WRITE (IUNIT,FORM1) ANG(IP:IPE)
WRITE (IUNIT,FORM3)
WRITE (IUNIT,FORM2) (FR(M),SPEC(IP:IPE,M),FR(M),ODSPEC(M),CDIR(M),M=1,ML)
WRITE (IUNIT,FORM1) ANG(IP:IPE)

END SUBROUTINE PRINT_SPECTRUM_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE REDUCED_TO_REGULAR_C (REDUCED_GRID, REGULAR_GRID,                   &
&                                NLON_RG, ZDELLO, XDELLO, ZMISS, DIR)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   REDUCED_TO_REGULAR -  INTERPOLATES FROM REDUCED TO REGULAR GRID.           !
!                                                                              !
!         H. GUENTHER    HZG   DECEMBER 2010                                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        INTERPOLATE PARAMETER FIELD FROM A REDUCED TO A REGULAR GRID.         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       LINEAR INTERPOLATION.                                                  !
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

REAL,    INTENT(IN)  :: REDUCED_GRID(:,:) !! INPUT FIELD ON REDUCED GRID.
REAL,    INTENT(OUT) :: REGULAR_GRID(:,:) !! OUTPUT FIELD ON REGULAR GRID.
INTEGER, INTENT(IN)  :: NLON_RG(:)        !! NO. OF LONGITUDES FOR EACH LAT.
INTEGER, INTENT(IN)  :: ZDELLO(:)         !! LONGITUDE INCREMENTS FOR EACH LAT.
INTEGER, INTENT(IN)  :: XDELLO            !! REGULAR LONGITUDE INCREMENT.
REAL,    INTENT(IN)  :: ZMISS             !! MISSING VALUE.
LOGICAL, INTENT(IN)  :: DIR               !! .TRUE. IF INPUT ARE DIRECTIONS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, DIMENSION(SIZE(ZDELLO))   :: ZDELLO_R
REAL    :: XDELLO_R

ZDELLO_R(:) = REAL(ZDELLO(:))
XDELLO_R    = REAL(XDELLO)

CALL REDUCED_TO_REGULAR (REDUCED_GRID, REGULAR_GRID,                           &
&                        NLON_RG, ZDELLO_R, XDELLO_R, ZMISS, DIR)

END SUBROUTINE REDUCED_TO_REGULAR_C

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE REDUCED_TO_REGULAR (REDUCED_GRID, REGULAR_GRID,                     &
&                              NLON_RG, ZDELLO, XDELLO, ZMISS, DIR)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   REDUCED_TO_REGULAR -  INTERPOLATES FROM REDUCED TO REGULAR GRID.           !
!                                                                              !
!         H. GUENTHER    HZG   DECEMBER 2010                                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        INTERPOLATE PARAMETER FIELD FROM A REDUCED TO A REGULAR GRID.         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       LINEAR INTERPOLATION.                                                  !
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

REAL,    INTENT(IN)  :: REDUCED_GRID(:,:) !! INPUT FIELD ON REDUCED GRID.
REAL,    INTENT(OUT) :: REGULAR_GRID(:,:) !! OUTPUT FIELD ON REGULAR GRID.
INTEGER, INTENT(IN)  :: NLON_RG(:)        !! NO. OF LONGITUDES FOR EACH LAT.
REAL,    INTENT(IN)  :: ZDELLO(:)         !! LONGITUDE INCREMENTS FOR EACH LAT.
REAL,    INTENT(IN)  :: XDELLO            !! REGULAR LONGITUDE INCREMENT.
REAL,    INTENT(IN)  :: ZMISS             !! MISSING VALUE.
LOGICAL, INTENT(IN)  :: DIR               !! .TRUE. IF INPUT ARE DIRECTIONS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: I, K, I_RED1, I_RED2
REAL    :: WH1, WH2
LOGICAL, SAVE :: FRSTIME = .TRUE.
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: I_RED_1, I_RED_2
REAL,    SAVE, ALLOCATABLE, DIMENSION(:,:) :: WH_1, WH_2
REAL,          ALLOCATABLE, DIMENSION(:,:) :: S_DIR, C_DIR

IF (FRSTIME) CALL INIT

IF (DIR) THEN
   ALLOCATE (S_DIR(SIZE(REDUCED_GRID,1),SIZE(REDUCED_GRID,2)))
   ALLOCATE (C_DIR(SIZE(REDUCED_GRID,1),SIZE(REDUCED_GRID,2)))
   DO K = 1,SIZE(REDUCED_GRID,2)
      DO I = 1,SIZE(REDUCED_GRID,1)
         IF (REDUCED_GRID(I,K).NE.ZMISS) THEN
            S_DIR(I,K) = SIN(REDUCED_GRID(I,K)*RAD)
            C_DIR(I,K) = COS(REDUCED_GRID(I,K)*RAD)
         ELSE
            S_DIR(I,K) = ZMISS
         END IF
      END DO
   END DO
   DO K = 1,SIZE(REGULAR_GRID,2)
      DO I = 1,SIZE(REGULAR_GRID,1)
         I_RED1 = I_RED_1(I,K)
         I_RED2 = I_RED_2(I,K)
         WH1 = WH_1(I,K)
         WH2 = WH_2(I,K)
         IF (S_DIR(I_RED1,K).EQ.ZMISS .AND. S_DIR(I_RED2,K).EQ.ZMISS) THEN
            REGULAR_GRID(I,K) = ZMISS
         ELSE IF (S_DIR(I_RED1,K).EQ.ZMISS) THEN
            IF (WH2.GE.0.5) THEN
               REGULAR_GRID(I,K) = REDUCED_GRID(I_RED2,K)
            ELSE
               REGULAR_GRID(I,K) = ZMISS
            END IF
         ELSE IF (S_DIR(I_RED2,K).EQ.ZMISS) THEN
            IF (WH1.GE.0.5) THEN
               REGULAR_GRID(I,K) = REDUCED_GRID(I_RED1,K)
            ELSE
               REGULAR_GRID(I,K) = ZMISS
            END IF
         ELSE
            REGULAR_GRID(I,K) = ATAN2(WH1*S_DIR(I_RED1,K)+WH2*S_DIR(I_RED2,K), &
&                               WH1*C_DIR(I_RED1,K)+WH2*C_DIR(I_RED2,K))*DEG
            IF (REGULAR_GRID(I,K).LT.0.) THEN
               REGULAR_GRID(I,K) = REGULAR_GRID(I,K)+360.
            END IF
            IF (REGULAR_GRID(I,K).GT.359.95) THEN
               REGULAR_GRID(I,K) = 0.
            END IF
         END IF
      END DO
   END DO
   DEALLOCATE (S_DIR, C_DIR)

ELSE

   DO K = 1,SIZE(REGULAR_GRID,2)
      DO I = 1,SIZE(REGULAR_GRID,1)
         I_RED1 = I_RED_1(I,K)
         I_RED2 = I_RED_2(I,K)
         WH1 = WH_1(I,K)
         WH2 = WH_2(I,K)
         IF (REDUCED_GRID(I_RED1,K).EQ.ZMISS .AND.                             &
&            REDUCED_GRID(I_RED2,K).EQ.ZMISS) THEN
            REGULAR_GRID(I,K) = ZMISS
         ELSE IF (REDUCED_GRID(I_RED1,K).EQ.ZMISS) THEN
            IF (WH2.GE.0.5) THEN
               REGULAR_GRID(I,K) = REDUCED_GRID(I_RED2,K)
            ELSE
               REGULAR_GRID(I,K) = ZMISS
            END IF
         ELSE IF (REDUCED_GRID(I_RED2,K).EQ.ZMISS) THEN
            IF (WH1.GE.0.5) THEN
               REGULAR_GRID(I,K) = REDUCED_GRID(I_RED1,K)
            ELSE
               REGULAR_GRID(I,K) = ZMISS
            END IF
         ELSE
            REGULAR_GRID(I,K) = WH1*REDUCED_GRID(I_RED1,K)                     &
&                             + WH2*REDUCED_GRID(I_RED2,K)
         END IF
      END DO
   END DO
END IF

CONTAINS

SUBROUTINE INIT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   INIT REDUCED_TO_REGULAR -  INTERPOLATES FROM REDUCED TO REGULAR GRID.      !
!                                                                              !
!         H. GUENTHER    HZG   DECEMBER 2010                                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        COMPUTE THE INTERPOLATION INDICES AND WEIGHTS TO BE USED BY           !
!        SUB. REDUCED TO A REGULAR GRID.                                       !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       LINEAR INTERPOLATION.                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCALE VARIABLES.                                                        !
!     -----------------                                                        !

   INTEGER :: N_RED
   REAL    :: X, D_RED
   LOGICAL :: IPER

   IPER = SIZE(REDUCED_GRID,1)*XDELLO.EQ. M_S_PER
   ALLOCATE (I_RED_1(SIZE(REGULAR_GRID,1),SIZE(REGULAR_GRID,2)))
   ALLOCATE (I_RED_2(SIZE(REGULAR_GRID,1),SIZE(REGULAR_GRID,2)))
   ALLOCATE (WH_1(SIZE(REGULAR_GRID,1),SIZE(REGULAR_GRID,2)))
   ALLOCATE (WH_2(SIZE(REGULAR_GRID,1),SIZE(REGULAR_GRID,2)))

   DO K = 1,SIZE(REGULAR_GRID,2)
      N_RED = NLON_RG(K)
      D_RED = XDELLO/ZDELLO(K)

      DO I = 1,SIZE(REGULAR_GRID,1)
         X = REAL(I-1)*D_RED
         I_RED1 = FLOOR(X)
         I_RED2 = CEILING(X)
         IF (I_RED1.EQ.I_RED2) THEN
            WH1 = 1.
            WH2 = 0.
         ELSE
            WH1 = REAL(I_RED1)-X+1.
            WH2 = 1.-WH1
         END IF
         IF (IPER) THEN
            I_RED1 = MOD(I_RED1+4*N_RED, N_RED)
            I_RED2 = MOD(I_RED2+4*N_RED, N_RED)
         ELSE
            IF (I_RED1.LT.0 .AND. WH1.LT. 0.1E-4) THEN
               I_RED1 = 0
               WH1 = 0.
               WH2 = 1.
            END IF
            IF (I_RED2.GE.N_RED .AND. WH2.LT. 0.1E-4) THEN
               I_RED2 = N_RED-1
               WH1 = 1.
               WH2 = 0.
            END IF
         END IF
         I_RED_1(I,K) = I_RED1+1
         I_RED_2(I,K) = I_RED2+1
         WH_1(I,K) = WH1
         WH_2(I,K) = WH2
      END DO
   END DO

   FRSTIME = .FALSE.

END SUBROUTINE INIT

END SUBROUTINE REDUCED_TO_REGULAR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_GENERAL_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PRINT_GENERAL_MODULE - PRINT STATUS OF WAM_GENERAL_MODULE.                 !
!                                                                              !
!     H.GUNTHER            HZG          APRIL 2015                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       MAKE A PRINTER OUTPUT OF THE DATA STORED IN WAM_GENERAL_MODULE.        !
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
!    1. WAM MODEL GLOBAL CONSTANTS.
!       ---------------------------

WRITE(IU06,'(/,'' -------------------------------------------------'')')
WRITE(IU06,*)'                WAM GENARAL MODULE DATA'
WRITE(IU06,'(  '' -------------------------------------------------'')')
WRITE(IU06,*) ' '
WRITE(IU06,*) ' ACCELERATION OF GRAVITY [M/S**2] .............. G = ', G
WRITE(IU06,*) ' PI ........................................... PI = ', PI
WRITE(IU06,*) ' 2.* PI ...................................... ZPI = ',  ZPI
WRITE(IU06,*) ' SQRT(PI) ................................ ZPISQRT = ', ZPISQRT
WRITE(IU06,*) ' COVERTION FROM RADIANS TO DEGREE ............ DEG = ', DEG
WRITE(IU06,*) ' COVERTION FROM DEGREE TO RADIANS ............ RAD = ', RAD
WRITE(IU06,*) ' EARTH CIRCUMFERENCE [M] .....................CIRC = ', CIRC
WRITE(IU06,*) ' EARTH RADIUS [M] .............................. R = ', R
WRITE(IU06,*) ' EULER CONSTANT .......................... GAMMA_E = ', GAMMA_E !! WAM-MAX

! ---------------------------------------------------------------------------- !
!
!    2. PARAMETERS FOR COUPLING.
!       ------------------------

WRITE(IU06,*) ' '
WRITE(IU06,*) ' AIR DENSITY ............................... ROAIR = ', ROAIR
WRITE(IU06,*) ' WATER DENSITY ........................... ROWATER = ', ROWATER
WRITE(IU06,*) ' ROAIR/ROWATER............................... XEPS = ',  XEPS
WRITE(IU06,*) ' 1./XEPS ................................. XINVEPS = ', XINVEPS
WRITE(IU06,*) ' PARAMETER FOR WIND INPUT ................ BETAMAX = ', BETAMAX
WRITE(IU06,*) ' SHIFTS GROWTH CURVE ........................ ZALP = ', ZALP
WRITE(IU06,*) ' MINIMUM CHARNOCK CONSTANT ................. ALPHA = ', ALPHA
WRITE(IU06,*) ' SHELTERING COEFFICIENT .............. TAUWSHELTER = ', TAUWSHELTER
WRITE(IU06,*) ' VON KARMAN CONSTANT ...................... XKAPPA = ', XKAPPA
WRITE(IU06,*) ' DEFAULT CHARNOCK VALUE FOR ICE AND DRY SEA  RCHAR = ', RCHAR

! ---------------------------------------------------------------------------- !
!
!    3. PARAMETERS FOR BFI.
!       -------------------

WRITE(IU06,*) ' '
WRITE(IU06,*) ' MAXIMUM VALUE OF DEPTH*WAVENUMBER ......... DKMAX = ', DKMAX
WRITE(IU06,*) ' MINIMUM VALUE ALLOWED FOR BFI SQUARED .... BF2MIN = ', BF2MIN
WRITE(IU06,*) ' MAXIMUM VALUE ALLOWED FOR BFI SQUARED .... BF2MAX = ', BF2MAX
WRITE(IU06,*) ' MINIMUM VALUE ALLOWED FOR KURTOSIS ........ C4MIN = ', C4MIN
WRITE(IU06,*) ' MAXIMUM VALUE ALLOWED FOR KURTOSIS ........ C4MAX = ', C4MAX

END SUBROUTINE PRINT_GENERAL_MODULE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_GENERAL_MODULE (B_MAX)

REAL, INTENT(IN) :: B_MAX

BETAMAX = B_MAX

END SUBROUTINE SET_GENERAL_MODULE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SORTI (iarray, index)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SORTI  - SORTS AN INTEGER ARRAY.
!                                                                              !
!     HEINZ GUNTHER    ECMWF/GKSS    JANUARY 1991
!                                                                              !
!    PURPOSE
!    -------
!                                                                              !
!        SORT AN INTEGER ARRAY BY AN INDEX ARRAY.
!                                                                              !
!     METHOD
!     ------
!                                                                              !
!       THE DATA IN THE INPUT ARRAY ARE SORTED AS DEFINED IN THE INPUT
!       INDEX ARRAY.
!       THIS ARRAY MAY BE GENEGATED BY SUB. SORTINI OR SORTINR.
!                                                                              !
!     EXTERNALS
!     ---------
!                                                                              !
!          NONE
!                                                                              !
!     REFERENCES
!     ----------
!                                                                              !
!          NONE
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     INTERFACE VARIABLES.
!     --------------------

INTEGER, DIMENSION(:), INTENT(INOUT) :: IARRAY  !! INPUT/ OUTPUT DATA ARRAY.
INTEGER, DIMENSION(:), INTENT(IN)    :: INDEX   !! INPUT INDEX ARRAY.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCALE VARIABLES.                                                        !
!     -----------------                                                        !

INTEGER, ALLOCATABLE, DIMENSION(:) :: IWORK

! ---------------------------------------------------------------------------- !
!                                                                              !
!    1. SORT DATA INTO WORK ARRAY.                                             !
!       --------------------------                                             !

ALLOCATE (IWORK(SIZE(IARRAY)))

IWORK(:) = IARRAY(INDEX(:))

! ---------------------------------------------------------------------------- !
!                                                                              !
!    2. COPY BACK TO DATA ARRAY.                                               !
!       ------------------------                                               !

IARRAY(:) = IWORK(:)

DEALLOCATE (IWORK)

END SUBROUTINE SORTI

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SORTINI (iarray, index)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     SORTINI* - SORTS AN INDEX ARRAY FORM AN INTEGER DATA ARRAY.
!                                                                              !
!     HEINZ GUNTHER    ECMWF/GKSS    JANUARY 1991
!                                                                              !
!     PURPOSE
!     -------
!                                                                              !
!        FIND THE INDEX ARRAY WHICH SORTS THE DATA.
!                                                                              !
!     METHOD
!     ------
!                                                                              !
!       THE INPUT ARRAY IS SCANNED AND AN INDEX ARRAY IS GENERATED,
!       WHICH CONTAINS THE INDEX OF THE INPUT ARRAY TO SORT THE DATA
!       IN INCREASING ORDER. THE INPUT ARRAY IS NOT CHANGED.
!                                                                              !
!     EXTERNALS
!     ---------
!                                                                              !
!          NONE
!                                                                              !
!     REFERENCES
!     ----------
!                                                                              !
!          NONE
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     INTERFACE VARIABLES.
!     --------------------

INTEGER, DIMENSION(:), INTENT(IN)   :: IARRAY  !! INPUT DATA TO BE SCANNED.
INTEGER, DIMENSION(:), INTENT(OUT)  :: INDEX   !! UTPUT ARRAY OF INDICES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCALE VARIABLES.                                                        !
!     -----------------                                                        !

INTEGER :: N, I, J, IH, IHM

! ---------------------------------------------------------------------------- !
!                                                                              !
!*   1. INITIAL INDEX ARRAY.
!       -------------------
N = SIZE(IARRAY)
DO I=1,N
   INDEX(I) = I
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!    2. CHECK FIRST TWO DATA.
!       ---------------------

IF (N.LT.2) RETURN
IF (IARRAY(2).LT.IARRAY(1)) THEN
   INDEX(2)=1
   INDEX(1)=2
ENDIF
IF (N.EQ.2) RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. DO REST OF DATA.
!       ----------------

OUTER: DO I=3,N
   IH = INDEX(I)
   IHM = INDEX(I-1)
   IF (IARRAY(IH).LT.IARRAY(IHM)) THEN
      INNER: DO J=I-1,2,-1
         INDEX(J+1) = INDEX(J)
         IF (IARRAY(INDEX(J-1)).LE.IARRAY(IH)) THEN
            INDEX(J) = IH
            CYCLE OUTER
         ENDIF
      END DO INNER
      INDEX(2) = INDEX(1)
      INDEX(1) = IH
   ENDIF
END DO OUTER

END SUBROUTINE SORTINI

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_GENERAL_MODULE

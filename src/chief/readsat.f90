SUBROUTINE READSAT (IU80, CDATE, RLAT, RLON, SWH, WS, EOFD)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READSAT  - READ MEASUREMENTS.                                              !
!                                                                              !
!     H. GUNTHER     GKSS/ECMWF       AUGUST 1991                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       READ ONE DATA RECORD FOR ASSIMILATION.                                 !
!                                                                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!         NONE.                                                                !
!                                                                              !
!      EXTERNALS.                                                              !
!      ----------                                                              !
!                                                                              !
!         NONE.                                                                !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE.
!     ----------

IMPLICIT NONE

INTEGER, PARAMETER :: KIND_D = 8

INTEGER,            INTENT(IN)  :: IU80  !! INPUT UNIT FOR MEASUREMENTS.
CHARACTER(LEN=14),  INTENT(OUT) :: CDATE !! DATE OF MEASUREMENT (YYYYMMDDHHMMSS).
REAL (KIND=KIND_D), INTENT(OUT) :: RLAT  !! LATITUDE OF MEASUREMENT (DEGREE).
REAL (KIND=KIND_D), INTENT(OUT) :: RLON  !! LONGITUDE OF MEASUREMENT (DEGREES).
REAL,               INTENT(OUT) :: SWH   !! WAVE HEIGHT.
REAL,               INTENT(OUT) :: WS    !! WIND SPEED.
LOGICAL,            INTENT(OUT) :: EOFD  !! END OF DATA. (.TRUE. , .FALSE.) 

! ---------------------------------------------------------------------------- !
! 
!     LOCAL VARIABLES.
!     ----------------

INTEGER  :: IOS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALIZE.                                                           !
!        -----------                                                           !

IOS = 0

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. READ THE MEASUREMENT.                                                 !
!        ---------------------                                                 !

READ (IU80,'(A14,4F10.4)',IOSTAT=IOS) CDATE, RLAT, RLON, SWH, WS

EOFD = (IOS.NE.0)

END SUBROUTINE READSAT

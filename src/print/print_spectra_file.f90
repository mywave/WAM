PROGRAM PRINT_SPECTRA_FILE
! ---------------------------------------------------------------------------- !
!                                                                              !
!     PRINT_SPECTRA_FILE - PRINTS SPECTRA FROM WAM-MODEL OUTPUT.               !
!                                                                              !
!      H. GUNTHER     GKSS/ECMWF  DECEMBER 1989                                !
!                     HZG         DECEMBER 2010      RE-ORGANISED              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        POSTPROCESSING OF WAM MODEL SPECTRA OUTPUT.                           !
!                                                                              !
!     INTERFACE.                                                               !
!     ----------                                                               !
!                                                                              !
!          IU01    INPUT UNIT OF SPECTRA FILE.                                 !
!          IU05    USER INPUT FILE.                                            !
!          IU06    PRINTER OUTPUT.                                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THIS PROGRAM READS THE WAM MODEL SPECTRA OUTPUTS AND EXTRACTS          !
!       SPECTRA AT SPECIFIED LOCATIONS AND TIMES.                              !
!       THE FILES ARE DYNAMICALLY ASSIGNED BY OPEN_FILE.                       !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!     -----------                                                              !
USE WAM_COORDINATE_MODULE

USE WAM_GENERAL_MODULE, ONLY:  &
&       INCDATE,               &  !! UPDATES A DATE/TIME GROUP.
&       PRINT_SPECTRUM,        &  !! PRINT A SPECTRUM.
&       OPEN_FILE                 !! OPEN A FILE.

USE WAM_PRINT_MODULE,   ONLY:  &
&       PRINT_SPECTRA_USER        !! PRINT A PROTOCOLL OF USER SETTINGS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY: IU05, FILE05, IU06, FILE06, ITEST
USE WAM_PRINT_MODULE, ONLY: CDATEA, CDATEE, IDELDO,                            &
&                           IU01, FILE01,  CDTFILE, IDFILE,                    &
&                           NOUT_S, PFLAG_S, TITL_S, NOUTT, COUTT,             &
&                           NOUTP, OUTLONG, OUTLAT, NAME, CFLAG_S,             &
&                           KL, ML, CO, FR, THETA,                             &
&                           SPEC_LAT, SPEC_LON, SPEC_DATE,                     &
&                           SPEC, U10, UDIR, US, DEPTH, CSPEED, CDIR,          &
&                           HS, PPER, MPER, TM1, TM2, MDIR, SPRE, TAUW,        &
&                           SPEC_SEA,                                          &
&                           HS_SEA, PPER_SEA, MPER_SEA, TM1_SEA, TM2_SEA,      &
&                           MDIR_SEA, SPRE_SEA,                                &
&                           SPEC_SWELL,                                        &
&                           HS_SWELL, PPER_SWELL, MPER_SWELL, TM1_SWELL,       &
&                           TM2_SWELL, MDIR_SWELL, SPRE_SWELL

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER     :: IFAIL,NFAIL     !! OPEN ERROR
LOGICAL     :: IEOF            !! END OF FILE ENCOUNTED IN SUB. READ_SPECTRUM

INTEGER            :: I
CHARACTER (LEN=40) :: HEADER

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITALISATION.                                                        !
!        --------------                                                        !

!     1.1 SET USER INPUT AND PROTOCOLL FILE NAMES.                             !

FILE05 = 'Spectra_User'
FILE06 = 'Spectra_Prot'

!     1.2  OPEN USER FILE AND READ USER INPUT.                                 !

OPEN (UNIT=IU06, FILE=FILE06, FORM='FORMATTED', STATUS="UNKNOWN")
CALL READ_SPECTRA_USER
CALL PRINT_SPECTRA_USER

!     1.3 FIRST AND LAST OUTPUT DATE.                                         !!

IF (NOUTT.GT.0) THEN
   CDATEE = COUTT(1)
   CDATEA = COUTT(1)
   DO I = 1,NOUTT
      IF (COUTT(I).LT.CDATEA) CDATEA = COUTT(I)
      IF (COUTT(I).GT.CDATEE) CDATEE = COUTT(I)
   END DO
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LOOP OVER INPUT FILES.                                                !
!        ----------------------                                                !
NFAIL=0
FILES: DO

!     2.1 FETCH FILE.                                                          !

   CALL OPEN_FILE (IU06, IU01, FILE01, CDTFILE, 'OLD', IFAIL)
!   IF (IFAIL.NE.0) STOP
   
! tbruns 26.01.2012
   IF (IFAIL.EQ.0) THEN

!     2.2  LOOP OVER OUTPUT TIMES.                                             !

   TIMES: DO

!     2.2.1 READ IN ONE OUTPUT TIME.                                           !

      CALL READ_SPECTRA_FILE (IU01, IEOF)

!     2.2.2 END OF FILE ENCOUNTED?                                             !

      IF (IEOF) EXIT TIMES   !! END OF FILE ENCOUNTERED
      IF (ITEST.GT.0) THEN
         WRITE (IU06,*) 'SUB. READ_SPECTRA_FILE DONE'
         WRITE (IU06,*) 'NEXT OUTPUT DATE, SPEC_DATE, SPEC_LAT, SPEC_LON: ',   &
&                        CDATEA, SPEC_DATE, SPEC_LAT, SPEC_LON
      END IF

!     2.2.3 OUTPUT TIME FOUND?                                                 !

      IF (SPEC_DATE.LT.CDATEA) CYCLE TIMES
      DO WHILE (SPEC_DATE.GT.CDATEA)
         CALL NEXT_OUTPUT_TIME
         IF (CDATEA.GT.CDATEE) EXIT FILES
         IF (SPEC_DATE.LT.CDATEA) CYCLE TIMES
      END DO

!     2.2.4 OUTPUT LOCATION?                                                   !

      LOCATION: DO I = 1,NOUTP
         IF (MOD(OUTLONG(I)-SPEC_LON+2*M_S_PER,+M_S_PER).EQ.0 .AND.           &
&           OUTLAT(I).EQ.SPEC_LAT) THEN
            IF (PFLAG_S(1) .AND. CFLAG_S(1)) THEN
               HEADER = TITL_S(1)(1:20)//NAME(I)
               CALL PRINT_SPECTRUM (IU06, SPEC_DATE, SPEC_LON, SPEC_LAT,       &
&                  HEADER, FR, THETA, SPEC, U10, UDIR, US, DEPTH, CSPEED, CDIR,&
&                  HS, PPER, MPER, TM1, TM2, MDIR, SPRE)
            END IF
            IF (PFLAG_S(2) .AND. CFLAG_S(2)) THEN
               HEADER = TITL_S(2)(1:20)//NAME(I)
               CALL PRINT_SPECTRUM (IU06, SPEC_DATE, SPEC_LON, SPEC_LAT,       &
&                  HEADER, FR, THETA, SPEC_SEA, U10, UDIR, US, DEPTH, CSPEED,  &
&                  CDIR, HS_SEA, PPER_SEA, MPER_SEA, TM1_SEA, TM2_SEA,         &
&                  MDIR_SEA, SPRE_SEA)
            END IF
            IF (PFLAG_S(3) .AND. CFLAG_S(3)) THEN
               HEADER = TITL_S(3)(1:20)//NAME(I)
               CALL PRINT_SPECTRUM (IU06, SPEC_DATE, SPEC_LON, SPEC_LAT,       &
&                  HEADER, FR, THETA, SPEC_SWELL, U10, UDIR, US, DEPTH, CSPEED,&
&                  CDIR, HS_SWELL, PPER_SWELL, MPER_SWELL, TM1_SWELL,          &
&                  TM2_SWELL, MDIR_SWELL, SPRE_SWELL)
            END IF
            CYCLE TIMES
         END IF
       END DO LOCATION

   END DO TIMES
  
   CLOSE (UNIT=IU01, STATUS='KEEP')   !! CLOSE OLD FILE
   ELSE
     NFAIL=NFAIL+1
     IF(NFAIL.gt.100) EXIT FILES
   ENDIF

   IF (CDATEA.EQ.CDATEE) EXIT FILES   !! ALL DONE?
   IF (IDFILE.GT.0) THEN
      CALL INCDATE (CDTFILE, IDFILE)  !! INCREMENT DATE FOR THE NEXT FILE.
   ELSE
      EXIT FILES
   END IF
END DO FILES

WRITE (*,*) ' PROGRAM PRINT_SPECTRA: ALL DONE'

! ---------------------------------------------------------------------------- !

CONTAINS

   SUBROUTINE NEXT_OUTPUT_TIME

   CHARACTER (LEN=14) :: IHH

      IF (NOUTT.EQ.0) THEN
         CALL INCDATE (CDATEA,IDELDO)
      ELSE
         IHH = '99999999999999'
         DO I=1,NOUTT
            IF (COUTT(I).GT.CDATEA .AND. COUTT(I).LT.IHH) IHH = COUTT(I)
         END DO
         CDATEA = IHH
      END IF

   END  SUBROUTINE NEXT_OUTPUT_TIME

END PROGRAM PRINT_SPECTRA_FILE

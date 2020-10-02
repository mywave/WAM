SUBROUTINE READ_SPECTRA_FILE (IUNIT, IEOF)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    READ_SPECTRA_FILE -  READS A SPECTRUM FROM WAVE MODEL OUTPUT FILE         !
!                                                                              !
!      M. DE LAS HERAS  KNMI/PCM  FEBRUARY  1990                               !
!      H. GUNTHER       HZG       DECEMBER 2010      RE-ORGANISED              !
!                                                                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       POST PROCESSING ROUTINE FOR WAVE MODEL.                                !
!       READ SPECTRA FROM WAMODEL OUTPUT FILE                                  !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       UNFORMATTED READ.                                                      !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!      ----------                                                              !

USE WAM_GENERAL_MODULE, ONLY:  &
&       ABORT1                  !! TERMINATES PROCESSING.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY: IU06

USE WAM_PRINT_MODULE, ONLY: KL, ML, CO, FR, THETA,                             &
&                           SPEC_LAT, SPEC_LON, SPEC_DATE,                     &
&                           NOUT_S, PFLAG_S,                                   &
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
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER,    INTENT(IN)  :: IUNIT    !! INPUT UNIT.
LOGICAL,    INTENT(OUT) :: IEOF     !! .TRUE. IF END OF FILE ENCOUNTED.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: I, IOS
REAL    :: XANG, XFRE, TH1, FR1
REAL, PARAMETER :: RAD = 3.1415927/180.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DATA HEADER FROM WAVE MODEL SPECTRA OUTPUT.                           !
!        -------------------------------------------                           !

IEOF = .FALSE.

READ (IUNIT,IOSTAT=IOS) SPEC_LON, SPEC_LAT, SPEC_DATE, XANG, XFRE, TH1, FR1, CO
IF (IOS.LT.0) THEN
   IEOF = .TRUE.
   RETURN
ELSE IF (IOS.GT.0) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_SPECTRA_FILE   *'
   WRITE(IU06,*) '*     =====================================   *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  READ ERROR IN FIRST HEADER RECORD          *'
   WRITE(IU06,*) '*  FROM  WAM SPECTRA OUTPUT FILE.             *'
   WRITE(IU06,*) '*  IOSTAT = ', IOS
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

READ (IUNIT, IOSTAT=IOS) PFLAG_S(1:NOUT_S)
IF (IOS.NE.0) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_SPECTRA_FILE   *'
   WRITE(IU06,*) '*     =====================================   *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  READ ERROR IN SECOND HEADER RECDRD         *'
   WRITE(IU06,*) '*  FROM  WAM SPECTRA OUTPUT FILE.             *'
   WRITE(IU06,*) '*  IOSTAT = ', IOS
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

KL = NINT(XANG)
ML = NINT(XFRE)
READ (IUNIT, IOSTAT=IOS) U10, UDIR, US, DEPTH, CSPEED, CDIR
IF (IOS.NE.0) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_SPECTRA_FILE   *'
   WRITE(IU06,*) '*     =====================================   *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  READ ERROR IN THIRD HEADER RECDRD         *'
   WRITE(IU06,*) '*  FROM  WAM SPECTRA OUTPUT FILE.             *'
   WRITE(IU06,*) '*  IOSTAT = ', IOS
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK ARRAYS.                                                         !
!        -------------                                                         !

IF (KL.LE.0 .OR. ML.LE.0) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_SPECTRA_FILE   *'
   WRITE(IU06,*) '*     =====================================   *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  NUMBER OF FREQUENCIES AND/OR DIRECTIONS    *'
   WRITE(IU06,*) '*  IS NOT POSITIVE.                           *'
   WRITE(IU06,*) '*  CHECK USER INPUT FOR FILE NAME             *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

IF (.NOT.ALLOCATED(FR)   ) THEN
   ALLOCATE (FR(1:ML))
   FR(1) = FR1              !! COMPUTE FREQUENCIES.
   DO I=2,ML
      FR(I) = CO*FR(I-1)
   END DO
END IF
IF (.NOT.ALLOCATED(THETA)) THEN
    ALLOCATE (THETA(1:KL))  !! COMPUTE DIRECTIONS.
    DO I = 1, KL
      THETA(I) = (TH1 + REAL(I-1)*360./XANG)*RAD
   END DO
END IF
IF (PFLAG_S(1).AND..NOT.ALLOCATED(SPEC)      ) ALLOCATE (SPEC(1:KL,1:ML))
IF (PFLAG_S(2).AND..NOT.ALLOCATED(SPEC_SEA)  ) ALLOCATE (SPEC_SEA(1:KL,1:ML))
IF (PFLAG_S(3).AND..NOT.ALLOCATED(SPEC_SWELL)) ALLOCATE (SPEC_SWELL(1:KL,1:ML))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. SPECTRA FROM WAVE MODEL SPECTRA OUTPUT.                               !
!        ---------------------------------------                               !

IF (PFLAG_S(1)) THEN
   READ(IUNIT, IOSTAT=IOS) HS, PPER, MPER, TM1, TM2, MDIR, SPRE
   READ(IUNIT, IOSTAT=IOS) SPEC
END IF
IF (PFLAG_S(2)) THEN
   READ(IUNIT, IOSTAT=IOS) HS_SEA, PPER_SEA, MPER_SEA, TM1_SEA, TM2_SEA,       &
&                          MDIR_SEA, SPRE_SEA
   READ(IUNIT, IOSTAT=IOS) SPEC_SEA
END IF
IF (PFLAG_S(3)) THEN
   READ(IUNIT, IOSTAT=IOS) HS_SWELL, PPER_SWELL, MPER_SWELL, TM1_SWELL,        &
&                          TM2_SWELL, MDIR_SWELL, SPRE_SWELL
   READ(IUNIT, IOSTAT=IOS) SPEC_SWELL
END IF
IF (IOS.NE.0) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_SPECTRA_FILE   *'
   WRITE(IU06,*) '*     =====================================   *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  READ ERROR IN SPECTRA DATA RECDRD          *'
   WRITE(IU06,*) '*  FROM  WAM SPECTRA OUTPUT FILE.             *'
   WRITE(IU06,*) '*  IOSTAT = ', IOS
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

END SUBROUTINE READ_SPECTRA_FILE

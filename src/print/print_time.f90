PROGRAM PRINT_TIME

! ---------------------------------------------------------------------------- !
!                                                                              !
!      PRINT_TIME - PRINT TIMESERIE FROM GRIDDED WAMODEL OUTPUT.               !
!                                                                              !
!      H. GUNTHER     GKSS/ECMWF  DECEMBER 1989                                !
!                     HZG         DECEMBER 2010      RE-ORGANISED              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!          POSTPROCESSING OF WAM MODEL INTEGRATED DATA.                        !
!          PRINT TIME SERIES FROM GRIDDED WAMODEL OUTPUT.                      !
!                                                                              !
!     INTERFACE.                                                               !
!     ----------                                                               !
!                                                                              !
!          IU01    INPUT UNIT OF INTEGRATED PARAMETER FILE.                    !
!          IU05    USER INPUT FILE.                                            !
!          IU06    PRINTER OUTPUT.                                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!          THIS PROGRAM READS THE WAM MODEL GRIDDED OUTPUT FILES AND           !
!          EXTRACTS TIMESERIES AT SPECIFIED LOCATIONS USING NEAREST NEIGHBOUR. !
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

USE WAM_COORDINATE_MODULE          !! COORDINATE TYPE AND PROCEDURES

USE WAM_GENERAL_MODULE, ONLY:  &
&       DIFDATE,               &  !! TIME DIFFERENCE.
&       INCDATE,               &  !! UPDATES A DATE/TIME GROUP.
&       OPEN_FILE                 !! OPEN A FILE.

USE WAM_PRINT_MODULE,   ONLY:  &
&       PRINT_TIME_USER,       &  !! PRINT A PROTOCOLL OF USER SETTINGS.
&       PREPARE_EXTRACTION

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY: IU05, FILE05, IU06, FILE06, ITEST
USE WAM_PRINT_MODULE, ONLY: CDATEA, CDATEE, IDELDO,                            &
&                           NOUTT, COUTT, NOUT_P,                              &
&                           NOUTP, OUTLONG, OUTLAT, NAME,                      &
&                           IU01, FILE01, CDTFILE, IDFILE,                     &
&                           PFLAG_P, CDTINTT, GRID

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

LOGICAL    :: IEOF            !! END OF FILE ENCOUNTED IN SUB. INGRID.
INTEGER    :: IFAIL           !! OPEN ERROR
INTEGER    :: IDELS           !! OUTPUT PERIOD.
INTEGER    :: NOU             !! NUMBER OF OUTPUT TIMESTEPS.
INTEGER    :: IS, IT, IP      !! COUNTER OF SITE, TIME, AND PARAMETER LOOPS.

REAL,     ALLOCATABLE :: SERIE(:,:,:) !! GATHERED TIMESERIES AT OUTPUT SITES.
INTEGER,  ALLOCATABLE :: I(:)         !! 1. GRID POINT INDICES FOR OUTPUT SITES.
INTEGER,  ALLOCATABLE :: K(:)         !! 2. GRID POINT INDICES FOR OUTPUT SITES.

REAL, PARAMETER :: ZMISS = -999.      !! INDICATOR OF MISSING VALUES.
character (len=13) :: ftext1, ftext2

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITALISATION.                                                        !
!        --------------                                                        !

!     1.1 SET USER INPUT AND PROTOCOLL FILE NAMES.                             !

FILE05 = 'Time_User'
FILE06 = 'Time_Prot'

!     1.2  OPEN USER FILE AND READ USER INPUT.                                 !

OPEN (UNIT=IU06, FILE=FILE06, FORM="FORMATTED", STATUS="UNKNOWN")
CALL READ_TIME_USER
CALL PRINT_TIME_USER

!     1.3 ALLOCATE OUTPUT ARRAYS.                                              !

CALL DIFDATE(CDATEA,CDATEE,IDELS)    !! OUTPUT PERIOD.
NOUTT = INT(IDELS/IDELDO)+1          !! NUMBER OF OUTPUT TIMESTEPS.

ALLOCATE (COUTT(NOUTT))
ALLOCATE (SERIE(NOUTT,NOUT_P,NOUTP))
SERIE = ZMISS
NOU = 0

ALLOCATE (I(NOUTP))
ALLOCATE (K(NOUTP))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LOOP OVER INPUT FILES.                                                !
!        ----------------------                                                !

FILES: DO

!     2.1 FETCH FILE.                                                          !

   CALL OPEN_FILE (IU06, IU01, FILE01, CDTFILE, 'OLD', IFAIL)
   IF (IFAIL.NE.0) EXIT FILES

!     2.2  LOOP OVER OUTPUT TIMES.                                             !

   TIMES: DO

!     2.2.1 READ IN WIND AND WAVE FIELDS.                                      !

      CALL READ_GRID_FILE (IU01, IEOF)

      IF (IEOF) EXIT TIMES     !! IF END OF FILE ENCOUNTED THEN EXIT TIME LOOP
      IF (ITEST.GT.0) THEN
         WRITE (IU06,*) 'SUB. READ_GRID_FILE DONE'
         WRITE (IU06,*) 'NEXT OUTPUT DATE, PARAMETER INPUT DATE: ',            &
&                        CDATEA, CDTINTT
      END IF

!     2.2.3 OUTPUT TIME FOUND?                                                 !

      IF (CDTINTT.LT.CDATEA) CYCLE TIMES
      DO WHILE (CDTINTT.GT.CDATEA)
         CALL INCDATE (CDATEA,IDELDO)
         IF (CDATEA.GT.CDATEE) EXIT FILES   !! ALL DONE?
         IF (CDTINTT.LT.CDATEA) CYCLE TIMES
      END DO

!     2.2.4. EXTRACT DATA. IF FIRST TIME COMPUTE GRID POINT INDICES.           !


      NOU = NOU+1
      IF (NOU.EQ.1) THEN
         CALL PREPARE_EXTRACTION (OUTLONG(1:NOUTP),OUTLAT(1:NOUTP), I, K)
         WRITE(IU06,*) ' '
         DO IS = 1,NOUTP
            IF (I(IS).EQ.-1) THEN
               ftext1 = write_coor_text (outlat(is))
               ftext2 = write_coor_text (outlong(is))
               WRITE(IU06,'('' SITE: '',A20,'' LAT. = '',A,'' LONG. = '' , &
&                            A, '' IS NOT IN GRID'')') NAME(IS),           &
&                            ftext1, ftext2
            END IF
         END DO
      END IF

      IF (NOU.GT.NOUTT) THEN
         NOU = NOUTT
         WRITE(IU06,'(''0'')')
         WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++'
         WRITE(IU06,*) ' +                                             +'
         WRITE(IU06,*) ' +         WARNING ERROR IN PRINT_TIME         +'
         WRITE(IU06,*) ' +         ===========================         +'
         WRITE(IU06,*) ' +                                             +'
         WRITE(IU06,*) ' + NUMBER OF OUTPUT TIMES EXCEEDS DIMENSION    +'
         WRITE(IU06,*) ' + DIMENSION IS  NOUTT = ', NOUTT
         WRITE(IU06,*) ' + LAST OUTPUT TIME IS ', COUTT(NOUTT)
         WRITE(IU06,*) ' +                                             +'
         WRITE(IU06,*) ' +      LATER TIMES ARE NOT PROCESSED          +'
         WRITE(IU06,*) ' +                                             +'
         WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++'
         EXIT FILES
      END IF

!     2.2.5. GATHER DATA FOR ALL SELECTED GRID POINTS IN OUTPUT ARRAYS.        !

      COUTT(NOU) = CDTINTT
      POINT: DO IS=1,NOUTP
         IF (I(IS).EQ.-1) CYCLE POINT
         DO IP = 1,NOUT_P
            IF (PFLAG_P(IP)) SERIE(NOU,IP,IS) = GRID(I(IS),K(IS),IP)
         END DO
      END DO POINT

!     2.2.6. NEXT OUTPUT TIME.                                                 !

      CALL INCDATE (CDATEA,IDELDO)      !! INCREMENT DATE FOR THE NEXT OUTPUT.
      IF (CDATEA.GT.CDATEE) EXIT FILES  !! ALL DONE?

   END DO TIMES

   CLOSE (UNIT=IU01, STATUS='KEEP')     !! CLOSE OLD FILE
   IF (IDFILE.GT.0) THEN
      CALL INCDATE (CDTFILE, IDFILE)    !! INCREMENT DATE FOR THE NEXT FILE.
   ELSE
      EXIT FILES
   END IF

END DO FILES

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. OUTPUT OF TIME SERIES.                                                !
!        ----------------------                                                !

SITE: DO IS = 1,NOUTP
   IF (I(IS).EQ.-1)  CYCLE SITE
   DO IT = 1,NOU
      IF (MOD(IT-1,60).EQ.0) THEN
         WRITE(IU06,'(''1'')')
         ftext1 = write_coor_text (outlat(is))
         ftext2 = write_coor_text (outlong(is))
         WRITE(IU06,'('' TIME SERIES OF INTEGRATED PARAMETERS AT '',A20,       &
&                     '' LAT. = '',A,'' LONG. = '',A)')  NAME(IS),             &
&                     ftext1, ftext2
         WRITE(IU06,*) ' '
         WRITE(IU06,*) '|----- DATE ----|-------- WIND ----------',            &
&                      '|-----| CURRENTS',                                     &
&                      '|--------------- WAVES -----------------',             &
&                      '|--------------- WINDSEA ---------------',             &
&                      '|--------------- SWELL ----------------|'

         WRITE(IU06,*) '                 U10  DIR.  US   CD  TAUW',            &
&                      ' DEPTH SPEED DIR.',                                    &
&                      '   HS    TP    TM   TM1   TM2  DIR. SPR.',             &
&                      '   HS    TP    TM   TM1   TM2  DIR. SPR.',             &
&                      '   HS    TP    TM   TM1   TM2  DIR. SPR.'
         WRITE(IU06,*) ' YYYYMMDDHHMMSS [M/S][DEG][M/S][1000] [%]',            &
&                      '  [M] [M/S][DEG]',                                     &
&                      '  [M]   [S]   [S]   [S]   [S] [DEG][DEG]',             &
&                      '  [M]   [S]   [S]   [S]   [S] [DEG][DEG]',             &
&                      '  [M]   [S]   [S]   [S]   [S] [DEG][DEG]'
      END IF
      WRITE(IU06,'(2X,A14,F6.1,F5.0,F5.2,F5.2,F5.2, F6.0, F5.2, F5.0,          &
&                 3(5F6.2,2F5.0))' )                                           &
&         COUTT(IT),                                                           &
&         SERIE(IT,1:3,IS), SERIE(IT,4,IS)*1000, SERIE(IT,16,IS),              &
&         SERIE(IT,6,IS), MAX(SERIE(IT,7,IS),-9.99), SERIE(IT,8,IS),           &
&         SERIE(IT,9:15,IS), SERIE(IT,17:23,IS), SERIE(IT,25:31,IS)

   END DO
END DO SITE

WRITE (*,*) ' PROGRAM PRINT_TIME: ALL DONE'

END PROGRAM PRINT_TIME

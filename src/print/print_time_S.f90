PROGRAM PRINT_TIME_S

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
&                           NOUTT, COUTT, NOUT_P, TITL_P, SCAL_P,              &
&                           NOUTP, OUTLONG, OUTLAT, NAME,                      &
&                           IU01, FILE01, CDTFILE, IDFILE,                     &
&                           CFLAG_P, PFLAG_P, CDTINTT, GRID

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
INTEGER    :: IOUT, II
INTEGER    :: IOUT_P(NOUT_P)

REAL,     ALLOCATABLE :: SERIE(:,:,:) !! GATHERED TIMESERIES AT OUTPUT SITES.
INTEGER,  ALLOCATABLE :: I(:)         !! 1. GRID POINT INDICES FOR OUTPUT SITES.
INTEGER,  ALLOCATABLE :: K(:)         !! 2. GRID POINT INDICES FOR OUTPUT SITES.

REAL, PARAMETER :: ZMISS = -999.      !! INDICATOR OF MISSING VALUES.
character (len=13) :: ftext1, ftext2

character (len=8) :: FORMA(1:NOUT_P)
character (len=8) :: FORMB(1:NOUT_P)
character (len=8) :: FORDI(1:NOUT_P)
character (len=600) :: PARAM_ZEILE
character (len=600) :: HEAD_ZEILE
character (len=600) :: DIMEN_ZEILE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITALISATION.                                                        !
!        --------------                                                        !

!     1.1 SET USER INPUT AND PROTOCOLL FILE NAMES.                             !

FILE05 = 'Time_User_S'
FILE06 = 'Time_Prot_S'

!     1.2  OPEN USER FILE AND READ USER INPUT.                                 !

OPEN (UNIT=IU06, FILE=FILE06, FORM="FORMATTED", STATUS="UNKNOWN")
CALL READ_TIME_USER_S
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

FORMB(1)  = '    WIND'
FORMB(2)  = '    WIND'
FORMB(3)  = '    WIND'
FORMB(4)  = '        '
FORMB(5)  = '        '
FORMB(6)  = '        '
FORMB(7)  = '    CURR'
FORMB(8)  = '    CURR'
FORMB(9)  = '        '
FORMB(10) = '        '
FORMB(11) = '        '
FORMB(12) = '        '
FORMB(13) = '        '
FORMB(14) = '        '
FORMB(15) = '        '
FORMB(16) = '        '
FORMB(17) = '     SEA'
FORMB(18) = '     SEA'
FORMB(19) = '     SEA'
FORMB(20) = '     SEA'
FORMB(21) = '     SEA'
FORMB(22) = '     SEA'
FORMB(23) = '     SEA'
FORMB(24) = '        '
FORMB(25) = '   SWELL'
FORMB(26) = '   SWELL'
FORMB(27) = '   SWELL'
FORMB(28) = '   SWELL'
FORMB(29) = '   SWELL'
FORMB(30) = '   SWELL'
FORMB(31) = '   SWELL'
FORMB(32) = '        '
FORMB(33) = '    GODA'
FORMB(34) = '        '
FORMB(35) = '        '
FORMB(36) = '        '
FORMB(37) = '        '
FORMB(38) = '   INTER'
FORMB(39) = '    PEAK'
FORMB(40) = '        '
FORMB(41) = '   SW 1 '
FORMB(42) = '   SW 1 '
FORMB(43) = '   SW 1 '
FORMB(44) = '   SW 2 '
FORMB(45) = '   SW 2 '
FORMB(46) = '   SW 2 '
FORMB(47) = '   SW 3 '
FORMB(48) = '   SW 3 '
FORMB(49) = '   SW 3 '
FORMB(50) = '        '
FORMB(51) = '   RADT '
FORMB(52) = '   RADT '
FORMB(53) = '   RADT '
FORMB(54) = '        '
FORMB(55) = '   WAVF '
FORMB(56) = '   WAVF '
FORMB(57) = '   STOK '
FORMB(58) = '   STOK '
FORMB(59) = '    ETO '
FORMB(60) = '    ETO '
FORMB(61) = '   MTOO '
FORMB(62) = '   MTOO '
FORMB(63) = '    ETO '
FORMB(64) = '   MTOO '
FORMB(65) = '   MTOO '
FORMB(66) = '        '
FORMB(67) = '  CRMAX '
FORMB(68) = '   HMAX '
FORMB(69) = '  CRMAX '
FORMB(70) = '   HMAX '

FORMA(1)  = '     U10'
FORMA(2)  = '    DIR.'
FORMA(3)  = '      US'
FORMA(4)  = '      CD'
FORMA(5)  = '    CHAR'
FORMA(6)  = '   DEPTH'
FORMA(7)  = '   SPEED'
FORMA(8)  = '    DIR.'
FORMA(9)  = '      HS'
FORMA(10) = '      TP'
FORMA(11) = '      TM'
FORMA(12) = '     TM1'
FORMA(13) = '     TM2'
FORMA(14) = '    DIR.'
FORMA(15) = '    SPR.'
FORMA(16) = '    TAUW'
FORMA(17) = '      HS'
FORMA(18) = '      TP'
FORMA(19) = '      TM'
FORMA(20) = '     TM1'
FORMA(21) = '     TM2'
FORMA(22) = '    DIR.'
FORMA(23) = '    SPR.'
FORMA(24) = '        '
FORMA(25) = '      HS'
FORMA(26) = '      TP'
FORMA(27) = '      TM'
FORMA(28) = '     TM1'
FORMA(29) = '     TM2'
FORMA(30) = '    DIR.'
FORMA(31) = '    SPR.'
FORMA(32) = '     ZO '
FORMA(33) = '      QP'
FORMA(34) = '    KURT'
FORMA(35) = '     BFI'
FORMA(36) = '    HMAX'
FORMA(37) = '    TMAX'
FORMA(38) = '      TP'
FORMA(39) = '    DIR.'
FORMA(40) = '    MSQS'
FORMA(41) = '      HS'
FORMA(42) = '     TM1'
FORMA(43) = '     DIR.'
FORMA(44) = '      HS'
FORMA(45) = '     TM1'
FORMA(46) = '    DIR.'
FORMA(47) = '      HS'
FORMA(48) = '     TM1'
FORMA(49) = '    DIR.'
FORMA(50) = '        '
FORMA(51) = '    SXX '
FORMA(52) = '    SYY '
FORMA(53) = '    SXY '
FORMA(54) = '        '
FORMA(55) = '    WFX '
FORMA(56) = '    WFY '
FORMA(57) = '    SDX '
FORMA(58) = '    SDY '
FORMA(59) = '     OC '
FORMA(60) = '     WA '
FORMA(61) = '      X '
FORMA(62) = '      Y '
FORMA(63) = '    BOT '
FORMA(64) = '   BOT_X'
FORMA(65) = '   BOT_Y'
FORMA(66) = '        '
FORMA(67) = ' FORRIS '
FORMA(68) = '  NAESS '
FORMA(69) = '  STQD  '
FORMA(70) = '  STQD  '

FORDI(1)  = '   [M/S]'
FORDI(2)  = '   [DEG]'
FORDI(3)  = '   [M/S]'
FORDI(4)  = '  [1000]'
FORDI(5)  = '  [1000]'
FORDI(6)  = '     [M]'
FORDI(7)  = '   [M/S]'
FORDI(8)  = '   [DEG]'
FORDI(9)  = '     [M]'
FORDI(10) = '     [S]'
FORDI(11) = '     [S]'
FORDI(12) = '     [S]'
FORDI(13) = '     [S]'
FORDI(14) = '   [DEG]'
FORDI(15) = '   [DEG]'
FORDI(16) = '     [%]'
FORDI(17) = '     [M]'
FORDI(18) = '     [S]'
FORDI(19) = '     [S]'
FORDI(20) = '     [S]'
FORDI(21) = '     [S]'
FORDI(22) = '   [DEG]'
FORDI(23) = '   [DEG]'
FORDI(24) = '        '
FORDI(25) = '     [M]'
FORDI(26) = '     [S]'
FORDI(27) = '     [S]'
FORDI(28) = '     [S]'
FORDI(29) = '     [S]'
FORDI(30) = '   [DEG]'
FORDI(31) = '   [DEG]'
FORDI(32) = '     [M]'
FORDI(33) = '        '
FORDI(34) = '   [100]'
FORDI(35) = '        '
FORDI(36) = '     [M]'
FORDI(37) = '     [S]'
FORDI(38) = '     [S]'
FORDI(39) = '   [DEG]'
FORDI(40) = '   [100]'
FORDI(41) = '     [M]'
FORDI(42) = '     [S]'
FORDI(43) = '   [DEG]'
FORDI(44) = '     [M]'
FORDI(45) = '     [S]'
FORDI(46) = '   [DEG]'
FORDI(47) = '     [M]'
FORDI(48) = '     [S]'
FORDI(49) = '   [DEG]'
FORDI(50) = '        '
FORDI(51) = '  KG/S/S'
FORDI(52) = '  KG/S/S'
FORDI(53) = '  KG/S/S'
FORDI(54) = '        '
FORDI(55) = '   N/M/M'
FORDI(56) = '   N/M/M'
FORDI(57) = '   [M/S]'
FORDI(58) = '   [M/S]'
FORDI(59) = 'KG/S/S/S'
FORDI(60) = 'KG/S/S/S'
FORDI(61) = 'KG/M/S/S'
FORDI(62) = 'KG/M/S/S'
FORDI(63) = 'KG/S/S/S'
FORDI(64) = 'KG/M/S/S'
FORDI(65) = 'KG/M/S/S'
FORDI(66) = '        '
FORDI(67) = '     [M]'
FORDI(68) = '     [M]'
FORDI(69) = '     [M]'
FORDI(70) = '     [M]'

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

         IOUT = 0
         DO IP = 1,NOUT_P
            IF (CFLAG_P(IP) .AND. PFLAG_P(IP)) THEN
               IOUT = IOUT+1
               IOUT_P(IOUT)=IP
            ELSE IF (CFLAG_P(IP) .AND. .NOT.PFLAG_P(IP)) THEN
                 WRITE(IU06, *) TITL_P(IP) , ' : REQUESTED PARAMETER IS NOT IN FILE'
            END IF
         END DO
         WRITE(IU06, *)
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

HEAD_ZEILE = ' '
PARAM_ZEILE = ' '
DIMEN_ZEILE = ' '
DIMEN_ZEILE(1:15) = ' YYYYMMDDHHMMSS'
HEAD_ZEILE(1:15)  = '       DATE    '

IS = 16
DO II = 1,IOUT
   HEAD_ZEILE(IS:IS+7)   = FORMB(IOUT_P(II))
   PARAM_ZEILE(IS:IS+7)  = FORMA(IOUT_P(II))
   DIMEN_ZEILE(IS:IS+7)  = FORDI(IOUT_P(II))
   IS = IS+8
END DO

SERIE(:,4:5,:) = SERIE(:,4:5,:)*1000.
SERIE(:,6,:)   = MIN(SERIE(:,6,:),999.)
SERIE(:,34,:) = SERIE(:,34,:)*100.
WHERE (SERIE(:,38,:).NE.0.)
   SERIE(:,38,:)   = 1./SERIE(:,38,:)
ELSEWHERE
   SERIE(:,38,:) = -1.
ENDWHERE
SERIE(:,40,:) = SERIE(:,40,:)*100.

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

         WRITE(IU06,*) HEAD_ZEILE
         WRITE(IU06,*) PARAM_ZEILE
         WRITE(IU06,*) DIMEN_ZEILE
      END IF

      WRITE(IU06,'(2X,A14,66F8.3)' ) COUTT(IT), (SERIE(IT,IOUT_P(II),IS),II = 1,IOUT)

   END DO
END DO SITE

WRITE (*,*) ' PROGRAM PRINT_TIME_S: ALL DONE'

END PROGRAM PRINT_TIME_S

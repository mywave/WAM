PROGRAM PRINT_GRID_FILE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     PRINT_GRID_FILE - PRINTS MAPS FROM GRIDDED WAMODEL OUTPUT.               !
!                                                                              !
!      H. GUNTHER     GKSS/ECMWF  DECEMBER 1989                                !
!                     HZG         DECEMBER 2010      RE-ORGANISED              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       POSTPROCESSING OF WAM MODEL INTEGRATED DATA.                           !
!       PRINTS MAPS FROM GRIDDED WAMODEL OUTPUT.                               !
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
!          PRINTS SELECTED FIELDS OF INTEGRATED DATA AT SPECIFIED TIMES.       !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!        NONE.                                                                 !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!      ----------                                                              !

USE WAM_GENERAL_MODULE, ONLY:  &
&       INCDATE,               &  !! UPDATES A DATE/TIME GROUP.
&       PRINT_ARRAY,           &  !! PRINT AN ARRAY.
&       OPEN_FILE,             &  !! !! OPEN A FILE.
&       REDUCED_TO_REGULAR        !! INTERPOLATES FROM REDUCED TO REGULAR GRID.                 

USE WAM_PRINT_MODULE,   ONLY:  &
&       PRINT_GRID_USER           !! PRINT A PROTOCOLL OF USER SETTINGS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY: IU05, FILE05, IU06, FILE06, ITEST
USE WAM_PRINT_MODULE, ONLY: CDATEA, CDATEE, IDELDO,                            &
&                           NOUTT, COUTT,  NOUT_P, CFLAG_P, TITL_P, SCAL_P,    &
&                           IU01, FILE01, CDTFILE, IDFILE,                     &
&                           NX, NY, AMOWEP, AMOSOP, AMOEAP, AMONOP,            &
&                           XDELLA, XDELLO, NLON_RG, ZDELLO,                   &
&                           PFLAG_P, CDTINTT, GRID, REGULAR

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

LOGICAL      :: IEOF             !! END OF FILE ENCOUNTED IN SUB. INGRID.
INTEGER      :: IFAIL            !! OPEN ERROR
INTEGER      :: I                !! LOOP COUNTER.
LOGICAL,SAVE :: FRSTIME = .TRUE.

LOGICAL,SAVE :: REDUCED_GRID = .FALSE.
REAL, ALLOCATABLE :: R_GRID(:,:)
REAL :: ZMISS = -995.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITALISATION.                                                        !
!        --------------                                                        !

!     1.1 SET USER INPUT AND PROTOCOLL FILE NAMES.                             !

FILE05 = 'Grid_User'
FILE06 = 'Grid_Prot'

!     1.2  OPEN USER FILE AND READ USER INPUT.                                 !

OPEN (UNIT=IU06, FILE=FILE06, FORM="FORMATTED", STATUS="UNKNOWN")
CALL READ_GRID_USER
CALL PRINT_GRID_USER

!     1.3 FIRST AND LAST OUTPUT DATE.                                          !

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
!     2. LOOP OVER OUTPUT FILES.                                               !
!        -----------------------                                               !

FILES: DO

!     2.1 FETCH FILE.                                                          !

   CALL OPEN_FILE (IU06, IU01, FILE01, CDTFILE, 'OLD', IFAIL)
   IF (IFAIL.NE.0) STOP

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


!     2.2.2 OUTPUT TIME FOUND?                                                 !

      IF (CDTINTT.LT.CDATEA) CYCLE TIMES
      DO WHILE (CDTINTT.GT.CDATEA)
         CALL NEXT_OUTPUT_TIME
         IF (CDATEA.GT.CDATEE) EXIT FILES   !! ALL DONE?
         IF (CDTINTT.LT.CDATEA) CYCLE TIMES
      END DO

!     2.2.3 DO OUTPUT OF REQUESTED FIELDS.                                     !

      WRITE (IU06,*) ' '
      IF (FRSTIME) THEN
         DO I = 1,NOUT_P
           IF (.NOT.PFLAG_P(I) .AND. CFLAG_P(I)) THEN
              WRITE(IU06,*) TITL_P(I), 'IS NOT STORED IN FILE'
            END IF
         END DO
	 XDELLO = (AMOEAP-AMOWEP) / REAL(NX-1)
	 REDUCED_GRID = MINVAL(NLON_RG).NE.NX
	 IF (REDUCED_GRID) THEN
            WRITE (IU06,*) ' INPUT DATA ARE ON A REDUCED GRID'
            IF (REGULAR) THEN
               WRITE (IU06,*) ' DATA ARE INTERPOLATED TO A REGULAR GRID'
            ELSE
               WRITE (IU06,*) '  DATA ARE NOT INTERPOLATED TO A REGULAR GRID'	    
               REDUCED_GRID = .FALSE.
	    END IF
         ELSE
            WRITE (IU06,*) ' INPUT DATA ARE ON A REGULAR GRID'
	 END IF
         REDUCED_GRID = REDUCED_GRID .AND. REGULAR
	 IF (REDUCED_GRID) ALLOCATE (R_GRID(NX,NY))
         FRSTIME = .FALSE.
         WRITE (IU06,*) ' '
      END IF
      DO I = 1,NOUT_P
         IF (PFLAG_P(I) .AND. CFLAG_P(I)) THEN
            IF (REDUCED_GRID) THEN
	           CALL REDUCED_TO_REGULAR (GRID(:,:,I), R_GRID,                   &
&                                       NLON_RG, ZDELLO, XDELLO)
               IF (I.EQ.5 ) GRID(:,:,I) = MIN(GRID(:,:,I),999.)
	           CALL PRINT_ARRAY (IU06, CDTINTT, TITL_P(I), R_GRID,             &
&                             AMOWEP, AMOSOP, AMOEAP, AMONOP, SCAL_P(I), ZMISS)
            ELSE 
               IF (I.EQ.5 ) GRID(:,:,I) = MIN(GRID(:,:,I),999.)
	           CALL PRINT_ARRAY (IU06, CDTINTT, TITL_P(I), GRID(:,:,I),        &
&                             AMOWEP, AMOSOP, AMOEAP, AMONOP, SCAL_P(I), ZMISS)
            END IF
         END IF
      END DO

!     2.2.7 NEXT OUTPUT TIME.                                                  !

      CALL NEXT_OUTPUT_TIME
      IF (CDATEA.GT.CDATEE) EXIT FILES  !! ALL DONE?

   END DO TIMES

   CLOSE (UNIT=IU01, STATUS='KEEP')     !! CLOSE OLD FILE
   IF (IDFILE.GT.0) THEN
      CALL INCDATE (CDTFILE, IDFILE)    !! INCREMENT DATE FOR THE NEXT FILE.
   ELSE
      EXIT FILES
   END IF
END DO FILES

WRITE (*,*) ' PROGRAM PRINT_GRID: ALL DONE'

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

END PROGRAM PRINT_GRID_FILE

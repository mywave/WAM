SUBROUTINE READ_GRID_FILE (IUNIT, IEOF)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      READ_GRID_FILE - READS WAVE MODEL OUTPUT FILE OF INTEGRATED DATA        !
!                                                                              !
!      H. GUNTHER     GKSS/ECMWF  DECEMBER 1989                                !
!                     HZG         DECEMBER 2010      RE-ORGANISED              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       POST PROCESSING ROUTINE FOR WAVE MODEL.                                !
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

USE WAM_PRINT_MODULE, ONLY:  &
&       ABORT1                  !! TERMINATES PROCESSING.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY: IU06
USE WAM_PRINT_MODULE, ONLY: NX, NY, NLON_RG, AMOWEP, AMOSOP, AMOEAP, AMONOP,   &
&                           CDTINTT, NOUT_P, PFLAG_P,                          &
&                           XDELLA, XDELLO, ZDELLO, GRID

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)  :: IUNIT !! INPUT UNIT FOR WAVE MODEL GRID OUTPUT FILE.
LOGICAL, INTENT(OUT) :: IEOF  !! END OF FILE INDICATOR.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: I, IOS
REAL    :: DNX, DNY
character (len=14) :: cstart

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DATA HEADER FROM WAVE MODEL GRIDDED OUTPUT.                           !
!        -------------------------------------------                           !

IEOF = .FALSE.

READ(IUNIT,IOSTAT=IOS) CDTINTT, DNX, DNY, AMOWEP, AMOSOP, AMOEAP, AMONOP, cstart
IF (IOS.LT.0) THEN
   IEOF = .TRUE.
   RETURN
ELSE IF (IOS.GT.0) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_GRID_FILE      *'
   WRITE(IU06,*) '*     ==================================      *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  READ ERROR IN FIRST HEADER RECORD          *'
   WRITE(IU06,*) '*  FROM GRIDDED WAM OUTPUT FILE.              *'
   WRITE(IU06,*) '*  IOSTAT = ', IOS
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

NX = NINT(DNX)
NY = NINT(DNY)
IF (.NOT.ALLOCATED(NLON_RG)) ALLOCATE(NLON_RG(1:NY))
IF (.NOT.ALLOCATED(ZDELLO) ) ALLOCATE(ZDELLO(1:NY) )

READ(IUNIT,IOSTAT=IOS) NLON_RG, ZDELLO
IF (IOS.LT.0) THEN
   IEOF = .TRUE.
   RETURN
ELSE IF (IOS.GT.0) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_GRID_FILE      *'
   WRITE(IU06,*) '*     ==================================      *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  READ ERROR IN SECOND HEADER RECDRD         *'
   WRITE(IU06,*) '*  FROM GRIDDED WAM OUTPUT FILE.              *'
   WRITE(IU06,*) '*  IOSTAT = ', IOS
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

READ(IUNIT, IOSTAT=IOS) PFLAG_P(1:NOUT_P)
IF (IOS.NE.0) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_GRID_FILE      *'
   WRITE(IU06,*) '*     ==================================      *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  READ ERROR IN THIRD HEADER RECDRD          *'
   WRITE(IU06,*) '*  FROM GRIDDED WAM OUTPUT FILE.              *'
   WRITE(IU06,*) '*  IOSTAT = ', IOS
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

IF (ANY(PFLAG_P(1:NOUT_P)) .AND. .NOT.ALLOCATED(GRID)) THEN
   ALLOCATE(GRID(NX,NY,NOUT_P))
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK ARRAYS.                                                         !
!        -------------                                                         !

IF (ALLOCATED(GRID) .AND. NX.NE.SIZE(GRID,1) .AND. NY.NE.SIZE(GRID,2)) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_GRID_FILE      *'
   WRITE(IU06,*) '*     ==================================      *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  GRID ARRAY DIMENSIONS ARE NOT CONSISTENT   *'
   WRITE(IU06,*) '*  DIMENSIONS READ FROM FILE ARE NX, NY = ', NY,NY
   WRITE(IU06,*) '*  DIMENSIONS OF ARRAY IS        NX, NY = ',                 &
&                                                     SIZE(GRID,1),SIZE(GRID,2)
   WRITE(IU06,*) '*  CHECK USER INPUT FOR FILE NAME             *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. DATA FROM WAVE MODEL GRIDDED OUTPUT.                                  !
!        ------------------------------------                                  !

DO I =1,NOUT_P
   IF (PFLAG_P(I)) THEN
      READ (IUNIT, IOSTAT=IOS) GRID(:,:,I)
      IF (IOS.NE.0) THEN
         WRITE(IU06,*) '***********************************************'
         WRITE(IU06,*) '*                                             *'
         WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_GRID_FILE      *'
         WRITE(IU06,*) '*     ==================================      *'
         WRITE(IU06,*) '*                                             *'
         WRITE(IU06,*) '*  READ ERROR IN GRID DATA RECORD             *'
         WRITE(IU06,*) '*  FROM GRIDDED WAM OUTPUT FILE FOR PARAMETER *'
         WRITE(IU06,*) '*       I = ', I
         WRITE(IU06,*) '*  IOSTAT = ', IOS
         WRITE(IU06,*) '*                                             *'
         WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
         WRITE(IU06,*) '***********************************************'
         CALL ABORT1
      END IF
   END IF
END DO

END SUBROUTINE READ_GRID_FILE

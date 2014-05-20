SUBROUTINE READ_RADIATION_FILE (IUNIT, IEOF)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      READ_RADIATION__FILE - READS RADIATION STRESS DATA OUTPUT FILE          !
!                                                                              !
!      H. GUNTHER     GKSS        DECEMBER 2005                                !
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
&                           CDTINTT, NOUT_R, PFLAG_R,                          &
&                           XDELLA, XDELLO, ZDELLO, GRID

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)  :: IUNIT !! INPUT UNIT RADIATION STRESS OUTPUT FILE.
LOGICAL, INTENT(OUT) :: IEOF  !! END OF FILE INDICATOR.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: I, IOS
REAL    :: DNX, DNY


! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DATA HEADER FROM WAVE MODEL RADIATION STRESS OUTPUT.                  !
!        ----------------------------------------------------                  !

IEOF = .FALSE.

READ(IUNIT,IOSTAT=IOS) CDTINTT, DNX, DNY, AMOWEP, AMOSOP, AMOEAP, AMONOP
IF (IOS.LT.0) THEN
   IEOF = .TRUE.
   RETURN
ELSE IF (IOS.GT.0) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_RAD_FILE       *'
   WRITE(IU06,*) '*     =================================       *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  READ ERROR IN FIRST HEADER RECORD          *'
   WRITE(IU06,*) '*  FROM RADIATION STRESS WAM OUTPUT FILE.     *'
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
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_RAD_FILE       *'
   WRITE(IU06,*) '*     =================================       *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  READ ERROR IN SECOND HEADER RECDRD         *'
   WRITE(IU06,*) '*  FROM RADIATION STRESS WAM OUTPUT FILE.     *'
   WRITE(IU06,*) '*  IOSTAT = ', IOS
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

READ(IUNIT, IOSTAT=IOS) PFLAG_R(1:NOUT_R)
IF (IOS.NE.0) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. READ_RAD_FILE       *'
   WRITE(IU06,*) '*     =================================       *'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  READ ERROR IN THIRD HEADER RECDRD          *'
   WRITE(IU06,*) '*  FROM RADIATION STRESS WAM OUTPUT FILE.     *'
   WRITE(IU06,*) '*  IOSTAT = ', IOS
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
   WRITE(IU06,*) '***********************************************'
   CALL ABORT1
END IF

IF (ANY(PFLAG_R(1:NOUT_R)).AND..NOT.ALLOCATED(GRID))  THEN
   ALLOCATE(GRID(NX,NY,NOUT_R))
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK ARRAYS.                                                         !
!        -------------                                                         !

IF (ALLOCATED(GRID) .AND. NX.NE.SIZE(GRID,1) .AND. NY.NE.SIZE(GRID,2)) THEN
   WRITE(IU06,*) '***********************************************'
   WRITE(IU06,*) '*                                             *'
   WRITE(IU06,*) '*  FATAL ERROR IN SUB. READ_RADIATION_FILE    *'
   WRITE(IU06,*) '*  =======================================    *'
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
!     3. DATA FROM WAVE MODEL RADIATION STRESS OUTPUT.                         !
!        ---------------------------------------------                         !

DO I =1,NOUT_R
   IF (PFLAG_R(I)) THEN
      READ (IUNIT, IOSTAT=IOS) GRID(:,:,I)
      IF (IOS.NE.0) THEN
         WRITE(IU06,*) '***********************************************'
         WRITE(IU06,*) '*                                             *'
         WRITE(IU06,*) '*  FATAL ERROR IN SUB. READ_RADIATION_FILE    *'
         WRITE(IU06,*) '*  =======================================    *'
         WRITE(IU06,*) '*                                             *'
         WRITE(IU06,*) '*  READ ERROR IN RADIATION DATA RECDRD        *'
         WRITE(IU06,*) '*  FROM RADIATION STRESS WAM OUTPUT FILE.     *'
         WRITE(IU06,*) '*       I = ', I
         WRITE(IU06,*) '*  IOSTAT = ', IOS
         WRITE(IU06,*) '*                                             *'
         WRITE(IU06,*) '*    PROGRAM ABORTS.   PROGRAM ABORTS.        *'
         WRITE(IU06,*) '***********************************************'
         CALL ABORT1
      ENDIF
   ENDIF
END DO

END SUBROUTINE READ_RADIATION_FILE

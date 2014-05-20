MODULE WAM_BOUNDARY_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE STORES THE BOUNDARY INPUT VALUES FOR A FINE GRID RUN.          !
!   THE VALUES WERE PRODUCED BY A PREVIOUS COARSE GRID RUN.                    !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1,                  &  !! TERMINATES PROCESSING.
&       DIFDATE,                 &  !! COMPUTE TIME DIFFERENCE.
&       OPEN_FILE,               &  !! OPENS A FILE.
&       INCDATE                     !! INCREMENT A DATE.

USE WAM_INTERFACE_MODULE, ONLY:  &
&       FEMEAN,                  &  !! COMPUTATION OF MEAN FREQUENCY.
&       INTSPEC,                 &  !! INTERPOLATE A SPECTRUM.
&       MEAN_DIRECTION,          &  !! COMPUTATION OF MEAN DIRECTION AND SPREAD.
&       TOTAL_ENERGY                !! COMPUTATION OF TOTAL ENERGY.
 
use wam_mpi_comp_module, only:   &  !! gather boundary values
&       mpi_gather_bound

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_FILE_MODULE,          ONLY: IU06, ITEST, FILE02, IU10, FILE10,         &
&                                   IU19, FILE19
USE WAM_OUTPUT_SET_UP_MODULE, ONLY: IDEL_OUT
USE WAM_TIMOPT_MODULE,        ONLY: CDATEE, CDTPRO, IDELPRO, IDELT, COLDSTART
USE WAM_FRE_DIR_MODULE,       ONLY: KL, ML, CO, FR, TH
USE WAM_GRID_MODULE,          ONLY: NX, NY, XDELLA, XDELLO,                    &
&                                   AMOWEP, AMOSOP, AMOEAP, AMONOP, IPER
USE WAM_MODEL_MODULE,         ONLY: FL3, DEPTH
USE WAM_NEST_MODULE,          ONLY: COARSE, FINE, N_NEST, n_code, MAX_NEST,    &
&                                   NBINP, NBOUNF,IJARF,IBFL,IBFR, BFW,        &
&                                   NBOUNC, BLNGC, BLATC, IJARC 

use wam_mpi_module,  only: irank, i_out_b_spec, nijs, nijl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. BOUNDARY VALUES OUTPUT FOR A FOLLOWING FINE GRID RUN.                 !
!        -----------------------------------------------------                 !

INTEGER            :: IDEL_B_OUT = -1   !! TIMESTEP TO SAVE BOUNDARY VALUES.
CHARACTER (LEN=14) :: CDT_B_OUT  = ' '  !! NEXT DATE TO SAVE BOUNDARY VALUES.
PUBLIC CDT_B_OUT

INTEGER            :: IDEL_BO_FILE = -1  !! TIMESTEP TO SAVE BOUNDARY FILES.
CHARACTER (LEN=14) :: CDT_BO_FILE  = ' ' !! NEXT DATE TO SAVE BOUNDARY FILES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. BOUNDARY VALUES FROM COARSE GRID.                                     !
!        ---------------------------------                                     !

INTEGER            :: IDEL_B_INP    !! TIMESTEP OF BOUNDARY INPUT VALUES.
INTEGER            :: IDEL_BI_FILE  !! TIMESTEP OF BOUNDARY FILE.
CHARACTER (LEN=14) :: CDT_BI_FILE   !! NEXT DATE TO FETCH BOUNDARY FILE.
REAL, ALLOCATABLE  :: XLAT(:)       !! LATITUDES OF INPUT SPECTRA.
REAL, ALLOCATABLE  :: XLON(:)       !! LONGITUDES OF INPUT SPECTRA.

!            FIRST TIME SPECTRA FROM COARSE GRID.                              !

CHARACTER (LEN=14) :: CDATE1 = ' '  !! DATE OF SPECTRA.
REAL, ALLOCATABLE  :: F1(:,:,:)     !! SPECTRA FROM COARSE GRID.
REAL, ALLOCATABLE  :: FMEAN1(:)     !! MEAN FREQUENCIES FROM COARSE GRID.
REAL, ALLOCATABLE  :: EMEAN1(:)     !! TOTAL ENERGIES FROM COARSE GRID.
REAL, ALLOCATABLE  :: THQ1(:)       !! MEAN DIRECTIONS FROM COARSE GRID (RAD).

!            SECOND TIME SPECTRA FROM COARSE GRID.                             !

CHARACTER (LEN=14) :: CDATE2 = ' '  !! DATE OF SPECTRA.
REAL, ALLOCATABLE  :: F2(:,:,:)     !! SPECTRA FROM COARSE GRID.
REAL, ALLOCATABLE  :: FMEAN2(:)     !! MEAN FREQUENCIES FROM COARSE GRID.
REAL, ALLOCATABLE  :: EMEAN2(:)     !! TOTAL ENERGIES FROM COARSE GRID.
REAL, ALLOCATABLE  :: THQ2(:)       !! MEAN DIRECTIONS FROM COARSE GRID (RAD).

!            TIME INTERPOLATED COARSE GRID SPECTRA.                            !

CHARACTER (LEN=14) :: CDATEI = ' '  !! DATE OF SPECTRA.
REAL, ALLOCATABLE  :: FI(:,:,:)     !! INTERPOLATED SPECTRUM.
REAL, ALLOCATABLE  :: FMEANI(:)     !! MEAN FREQUENCIES FROM COARSE GRID.
REAL, ALLOCATABLE  :: EMEANI(:)     !! TOTAL ENERGIES FROM COARSE GRID.
REAL, ALLOCATABLE  :: THQI(:)       !! MEAN DIRECTIONS FROM COARSE GRID (RAD).

PUBLIC CDT_BI_FILE, XLON, XLAT, IDEL_B_INP, IDEL_BI_FILE
PUBLIC CDATE2, EMEAN2, THQ2, FMEAN2, F2

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE BOUNDARY_INPUT                  !! BOUNDARY VALUE INPUT.
   MODULE PROCEDURE BOUNDARY_INPUT
END INTERFACE
PUBLIC BOUNDARY_INPUT

INTERFACE BOUNDARY_OUTPUT                 !! OUTPUT OF BOUNDARY VALUES.
   MODULE PROCEDURE BOUNDARY_OUTPUT
END INTERFACE
PUBLIC BOUNDARY_OUTPUT

INTERFACE PREPARE_BOUNDARY                !! PREPARES BOUNDARY MODULE.
   MODULE PROCEDURE PREPARE_BOUNDARY 
END INTERFACE
PUBLIC PREPARE_BOUNDARY

INTERFACE PRINT_BOUNDARY_STATUS           !! PRINTS BOUNDARY MODULE.
   MODULE PROCEDURE PRINT_BOUNDARY_STATUS 
END INTERFACE
PUBLIC PRINT_BOUNDARY_STATUS

INTERFACE SET_BOUNDARY_OUTPUT_TIMESTEPS   !! DEFINES BOUNDARY OUTPUT STEPS.
   MODULE PROCEDURE SET_BOUNDARY_OUTPUT_TIMESTEPS 
END INTERFACE
PUBLIC SET_BOUNDARY_OUTPUT_TIMESTEPS

INTERFACE
   SUBROUTINE READ_BOUNDARY_INPUT         !! READS ONE SET OF BOUNDARY SPECTRA.
   END SUBROUTINE READ_BOUNDARY_INPUT         
END INTERFACE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE SAVE_BOUNDARY_FILE         !! SAVES AND OPENS BOUNDARY OUTPUT FILES.
   MODULE PROCEDURE SAVE_BOUNDARY_FILE 
END INTERFACE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE BOUNDARY_INPUT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   BOUNDARY_INPUT - BOUNDARY VALUE INPUT INTO THE WAM MODEL.                  !
!                                                                              !
!     H. GUNTHER    GKSS/ECMWF   JANUARY 1991                                  !
!     H. GUNTHER    GKSS         JANUARY 2002    FT90                          !
!                                          - TIME INTERPOLATION OF SPECTRA.    !
!     E. MYKLEBUST               FEBRUARY 2005   MPI PARALLELIZATION           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO INTERPOLATE BOUNDARY SPECTRA IN TIME AND SPACE AND INSERT THE       !
!       BOUNDARY SPECTRA INTO THE WAM MODEL FIELD.                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE SUB. READS A COMPLETE SET OF BOUNDARY VALUES IF NECCESSARY.        !
!       IT PERFORMS THE TIME AND SPACE INTERPOLATION AND INSERTS THE SPECTA    !
!       INTO THE MODEL FIELD. INDICES AND WEIGHTS NECESSARY FOR THE SPACE      !
!       INTERPOLATION AND STORAGE ARE PRECOMPUTED IN PROG. PREPROC.            !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER      :: IJ, IDEL1L, IJF, IBCL, IBCR
REAL         :: DEL12, DEL1L, FMEAN, EMEAN, THQ

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. READ NEW INPUT IF REQUIRED.                                           !
!        ---------------------------                                           !

DO WHILE (CDATE2.LT.CDTPRO)
   CDATE1 = CDATE2
   FMEAN1 = FMEAN2
   EMEAN1 = EMEAN2
   THQ1   = THQ2
   F1 = F2
   CALL READ_BOUNDARY_INPUT
   IF (ITEST.GT.3)  WRITE (IU06,*)                                            &
&  '       SUB. BOUNDARY_INPUT: SECOND BOUNDARY VALUES READ CDATE2 = ', CDATE2
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TIME INTERPOLATION.                                                   !
!        -------------------                                                   !

IF (CDTPRO.EQ.CDATE1 ) THEN
   CDATEI = CDATE1
   FMEANI(1:) = FMEAN1
   EMEANI(1:) = EMEAN1
   THQI  (1:) = THQ1
   FI(:,:,1:) = F1
ELSE IF (CDTPRO.EQ.CDATE2 .OR. CDATE1.EQ.' ') THEN
   CDATEI = CDATE1
   FMEANI(1:) = FMEAN2
   EMEANI(1:) = EMEAN2
   THQI  (1:) = THQ2
   FI(:,:,1:) = F2
ELSE IF (CDTPRO.GT.CDATE1 .AND. CDTPRO.LT.CDATE2) THEN
   CDATEI = CDTPRO
   CALL DIFDATE (CDATE1, CDATEI, IDEL1L)
   DEL1L = REAL(IDEL1L)
   DEL12 = REAL(IDEL_B_INP)
   DO IJ = 1,NBINP
      CALL INTSPEC (DEL12, DEL1L,                                              &
&                   F1(:,:,IJ), FMEAN1(IJ), EMEAN1(IJ), THQ1(IJ),              &
&                   F2(:,:,IJ), FMEAN2(IJ), EMEAN2(IJ), THQ2(IJ),              &
&                   FI(:,:,IJ), FMEANI(IJ), EMEANI(IJ), THQI(IJ))
   END DO
   IF (ITEST.GT.3)  WRITE (IU06,*)                                             &
&  '       SUB. BOUNDARY_INPUT: TIME INTERPOLATION DONE     CDATEI = ', CDATEI
ELSE
   WRITE (IU06,*) '*******************************************'
   WRITE (IU06,*) '*                                         *'
   WRITE (IU06,*) '*   FATAL ERROR SUB. BOUNDARY_INPUT.      *'
   WRITE (IU06,*) '*   ================================      *'
   WRITE (IU06,*) '* DATES DO NOT MATCH.                     *'
   WRITE (IU06,*) '* DATE OF FIRST SPECTRA IS  CDATE1 =  ', CDATE1
   WRITE (IU06,*) '* MODEL DATE IS             CDTPRO =  ', CDTPRO
   WRITE (IU06,*) '* DATE OF SECOND SPECTRA IS CDATE2 =  ', CDATE2
   WRITE (IU06,*) '*                                         *'
   WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
   WRITE (IU06,*) '*                                         *'
   WRITE (IU06,*) '*******************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. SPACE INTERPOLATION.                                                  !
!        --------------------                                                  !

DEL12 = 1.
DO IJ = 1,NBOUNF
      IJF = IJARF(IJ)
      IF(IJF.GE.nijs .and. IJF.LE.nijl) THEN
         IBCL = IBFL(IJ)
         IBCR = IBFR(IJ)
         DEL1L = BFW(IJ)
         CALL INTSPEC (DEL12, DEL1L,                                          &
&                   FI(:,:,IBCL), FMEANI(IBCL), EMEANI(IBCL), THQI(IBCL),     &
&                   FI(:,:,IBCR), FMEANI(IBCR), EMEANI(IBCR), THQI(IBCR),     &
&                   FL3(IJF,:,:), FMEAN, EMEAN, THQ)
      END IF
END DO

END SUBROUTINE BOUNDARY_INPUT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE BOUNDARY_OUTPUT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   BOUNDARY_OUTPUT - OUTPUT OF THE COARSE GRID BOUNDARY VALUES.               !
!                                                                              !
!     R. PORTZ     MPI          JANUARY 1991                                   !
!     H. GUENTHER  GKSS         JANUARY 2002    FT90.                          !
!     A. Behrens   MSC/GKSS     December 2003   Message passing                !
!     E. Myklebust              November 2004   MPI parallelization            !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        WRITE THE BOUNDARY VALUE OUTPUT FILE.                                 !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       SEQUENCIAL UNFORMATTED WRITE TO UNIT.                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, itag, iz, i, iu
REAL    :: FBC(1:MAX_NEST,SIZE(FL3,2),SIZE(FL3,3),N_NEST)
REAL    :: THQC   (1:MAX_NEST)
REAL    :: EMEANC (1:MAX_NEST)
REAL    :: FMEANC (1:MAX_NEST)
logical, save :: first = .true.

save iz
if (first) then
   iz = 0
   first = .false.
endif
iz = iz+1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. GATHER BOUNDARY SPECTRA.                                              !
!        ------------------------                                              !

itag  = 850+iz
call mpi_gather_bound (i_out_b_spec, itag, fl3, fbc)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN PARAMETERS.                                              !
!        ------------------------                                              !

if (irank==i_out_b_spec) then

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. WRITE BOUNDARY SPECTRA (ascii or binary code)                         !
!        ---------------------------------------------                         !

   IU = IU19-1
   DO I = 1, N_NEST
      CALL TOTAL_ENERGY (FBC(1:NBOUNC(I),1:KL,1:ML,I), EMEANC(1:NBOUNC(I)))
      CALL FEMEAN (FBC(1:NBOUNC(I),1:KL,1:ML,I), EMEANC(1:NBOUNC(I)),          &
&                  FM=FMEANC(1:NBOUNC(I)))
      CALL MEAN_DIRECTION (FBC(1:NBOUNC(I),1:KL,1:ML,I), THQ=THQC(1:NBOUNC(I)))
      IU = IU+1
      if (n_code(i)==1) then
         do ij = 1,nbounc(i)
            write (iu,*) BLNGC(IJ,I), BLATC(IJ,I), CDTPRO, EMEANC(IJ),         &
&                        THQC(IJ), FMEANC(IJ)
            write (iu,*) FBC(IJ,1:KL,1:ML,I)
         enddo
      else
         DO IJ = 1,NBOUNC(I)
            WRITE(IU) BLNGC(IJ,I), BLATC(IJ,I), CDTPRO, EMEANC(IJ),            &
&                     THQC(IJ), FMEANC(IJ)
            WRITE(IU) FBC(IJ,1:KL,1:ML,I)
         END DO
      endif
   END DO
endif

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. NEXT DATE FOR OUTPUT.                                                 !
!        ---------------------                                                 !

CALL INCDATE (CDT_B_OUT, IDEL_B_OUT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. SAVE BOUNDARY FILE.                                                   !
!        -------------------                                                   !

IF (CDT_BO_FILE.EQ.CDTPRO .OR. CDATEE.EQ.CDTPRO) THEN
   CALL SAVE_BOUNDARY_FILE 
END IF

END SUBROUTINE BOUNDARY_OUTPUT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_BOUNDARY

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PREPARE_BOUNDARY - PREPARES THE BOUNDARY MODULE FOR COARSE AND FINE GRIDS. !
!                                                                              !
!     H. GUENTHER  GKSS         JANUARY 2010                                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        TO INITIAL THE BOUNDARY MODULE.                                       !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       - CHECK THE BOUNDAY OPTINS AND TIMESTEPS.                              !
!       - SETS DEFAULT VALUES IF NECCESSARY.                                   !
!       - OPENS FIRST OUTPUT FILE (COARSE GRID)                                !
!       - READS FIRST BOUNDARY VALUES (FINE GRID)                              !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CHECK NEST OPTION WITH PREPROC OUTPUT.                                !
!        --------------------------------------                                !

IF (COARSE .AND. N_NEST.LE.0) THEN
   COARSE = .FALSE.
   WRITE (IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) '+                                             +'
   WRITE (IU06,*) '+    WARNING ERROR SUB. PREPARE_BOUNDARY.     +'
   WRITE (IU06,*) '+    ====================================     +'
   WRITE (IU06,*) '+                                             +'
   WRITE (IU06,*) '+ NESTS ARE NOT DEFINED IN THE PREPROC FILE   +'
   WRITE (IU06,*) '+ BUT COARSE GRID RUN IS REQUESTED.           +'
   WRITE (IU06,*) '+                                             +'
   WRITE (IU06,*) '+  MODEL OPTION CHANGED TO COARSE = .FALSE.   +'
   WRITE (IU06,*) '+                                             +'
   WRITE (IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
END IF
IF (FINE .AND. NBINP.LE.0) THEN
   WRITE (IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) '+                                             +'
   WRITE (IU06,*) '+    WARNING ERROR SUB. PREPARE_BOUNDARY.     +'
   WRITE (IU06,*) '+    ====================================     +'
   WRITE (IU06,*) '+                                             +'
   WRITE (IU06,*) '+ NESTS ARE NOT DEFINED IN THE PREPROC FILE   +'
   WRITE (IU06,*) '+ BUT FINE GRID RUN IS REQUESTED.             +'
   WRITE (IU06,*) '+                                             +'
   WRITE (IU06,*) '+  MODEL OPTION CHANGED TO FINE = .FALSE.     +'
   WRITE (IU06,*) '+                                             +'
   WRITE (IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
   FINE = .FALSE.
END IF

IF (.NOT. COARSE .AND. .NOT. FINE) RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK COARSE GRID OUTPUT TIMESTEPS AND OPEN FIRST OUTPUT FILE.        !
!        --------------------------------------------------------------        !

IF (COARSE) THEN
 
   IF (IDEL_B_OUT   .LE. 0) IDEL_B_OUT   = IDELPRO 
   IF (IDEL_BO_FILE .LE. 0) IDEL_BO_FILE = MAX(IDEL_OUT,0)

   IF (MOD(IDEL_B_OUT,IDELPRO).NE.0) THEN
      WRITE (IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) '+                                             +'
      WRITE (IU06,*) '+    WARNING ERROR SUB. PREPARE_BOUNDARY.     +'
      WRITE (IU06,*) '+    ====================================     +'
      WRITE (IU06,*) '+                                             +'
      WRITE (IU06,*) ' + OUTPUT TIMESTEP FOR BOUDNDARY VALUES IS    +'
      WRITE (IU06,*) ' + NOT A MULTIPLE OF PROPAGATION TIMESTEP.    +'
      WRITE (IU06,*) ' + OUTPUT      TIMESTEP WAS: ',IDEL_B_OUT, ' SECONDS'
      WRITE (IU06,*) ' + PROPAGATION TIMESTEP  IS: ',IDELPRO, ' SECONDS'
      WRITE (IU06,*) ' +                                            +'
      WRITE (IU06,*) ' + TIMESTEP IS CHANGED TO NEAREST MULTIPLE OF +'
      WRITE (IU06,*) ' + THE PROPAGATION TIMESTEP.                  +'
      IDEL_B_OUT = MAX(NINT(REAL(IDEL_B_OUT)/REAL(IDELPRO)),1)*IDELPRO
      WRITE (IU06,*) ' + NEW OUTPUT TIMESTEP   IS: ',IDEL_B_OUT, ' SECONDS'
      WRITE (IU06,*) ' +                                            +'
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++'
   END IF

   IF (MOD(IDEL_BO_FILE,IDEL_B_OUT).NE.0) THEN
      WRITE (IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) '+                                             +'
      WRITE (IU06,*) '+    WARNING ERROR SUB. PREPARE_BOUNDARY.     +'
      WRITE (IU06,*) '+    ====================================     +'
      WRITE (IU06,*) '+                                             +'
      WRITE (IU06,*) ' + FILE SAVE TIMESTEP FOR BOUDNDARY VALUES IS +'
      WRITE (IU06,*) ' + NOT A MULTIPLE OF BOUDNDARY OUTPUT TIMESTEP.+'
      WRITE (IU06,*) ' + FILE SAVE   TIMESTEP WAS: ',IDEL_BO_FILE, ' SECONDS'
      WRITE (IU06,*) ' + OUTPUT TIMESTEP       IS: ',IDEL_B_OUT, ' SECONDS'
      WRITE (IU06,*) ' +                                            +'
      WRITE (IU06,*) ' + TIMESTEP IS CHANGED TO NEAREST MULTIPLE OF +'
      WRITE (IU06,*) ' + THE OUTPUT TIMESTEP.                       +'
      IDEL_BO_FILE = MAX(NINT(REAL(IDEL_BO_FILE)/REAL(IDEL_B_OUT)),1)*IDEL_B_OUT
      WRITE (IU06,*) ' + NEW FILE TIMESTEP     IS: ',IDEL_BO_FILE,' SECONDS'
      WRITE (IU06,*) ' +                                            +'
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++'
   END IF

!     TIME COUNTER FOR NEXT OUTPUT.                                         !

   CDT_B_OUT = ' '
   IF (IDEL_B_OUT .GT. 0.) THEN
      CDT_B_OUT = CDTPRO
      IF (.NOT.COLDSTART) CALL INCDATE (CDT_B_OUT, IDEL_B_OUT)
   END IF

!     TIME COUNTER FOR NEXT OUTPUT FILE AND OPEN FIRST OUTPUT FILE.        !

   CDT_BO_FILE = ' '
   IF (IDEL_BO_FILE.EQ.0.) THEN
      CDT_BO_FILE = CDATEE
   ELSE
      CDT_BO_FILE = CDTPRO
      IF (COLDSTART) CALL INCDATE (CDT_BO_FILE, -IDEL_BO_FILE)
   END IF
   CALL SAVE_BOUNDARY_FILE
   IF (COLDSTART) CALL BOUNDARY_OUTPUT
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ALLOCATE FINE GRID INPUT ARRAYS AND DO FIRST BOUNDAY INPUT.           !
!        -----------------------------------------------------------           !

IF (FINE) THEN
   IF (.NOT.ALLOCATED(XLON)  ) ALLOCATE (XLON(1:NBINP))
   IF (.NOT.ALLOCATED(XLAT)  ) ALLOCATE (XLAT(1:NBINP))
   IF (.NOT.ALLOCATED(F1)    ) ALLOCATE (F1(1:KL,1:ML,1:NBINP))
   IF (.NOT.ALLOCATED(FMEAN1)) ALLOCATE (FMEAN1(1:NBINP))
   IF (.NOT.ALLOCATED(EMEAN1)) ALLOCATE (EMEAN1(1:NBINP))
   IF (.NOT.ALLOCATED(THQ1)  ) ALLOCATE (THQ1(1:NBINP))
   IF (.NOT.ALLOCATED(F2)    ) ALLOCATE (F2(1:KL,1:ML,1:NBINP))
   IF (.NOT.ALLOCATED(FMEAN2)) ALLOCATE (FMEAN2(1:NBINP))
   IF (.NOT.ALLOCATED(EMEAN2)) ALLOCATE (EMEAN2(1:NBINP))
   IF (.NOT.ALLOCATED(THQ2)  ) ALLOCATE (THQ2(1:NBINP))
   IF (.NOT.ALLOCATED(FI)    ) ALLOCATE (FI(1:KL,1:ML,0:NBINP))
   IF (.NOT.ALLOCATED(FMEANI)) ALLOCATE (FMEANI(0:NBINP))
   IF (.NOT.ALLOCATED(EMEANI)) ALLOCATE (EMEANI(0:NBINP))
   IF (.NOT.ALLOCATED(THQI)  ) ALLOCATE (THQI(0:NBINP))
   FI(:,:,0) = 0.
   FMEANI(0) = 0.
   EMEANI(0) = 0.
   THQI(0)   = 0.
   
   IDEL_BI_FILE = 0
   CDT_BI_FILE = CDTPRO
   CDATE2 = '99991231235959'
   CALL READ_BOUNDARY_INPUT
   DO WHILE (CDATE2.LT.CDTPRO)
      CALL READ_BOUNDARY_INPUT
   END DO
   IF (ITEST.GT.3)  WRITE (IU06,*)                                            &
&  '       SUB. PREPARE_BOUNDARY: FIRST BOUNDARY VALUES READ CDATE2 = ', CDATE2

END IF

END SUBROUTINE PREPARE_BOUNDARY

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_BOUNDARY_STATUS

WRITE (IU06,*) ' '
WRITE (IU06,*) ' ------------------------------------------------- '
WRITE (IU06,*) '              BOUNDARY MODULE STATUS:'
WRITE (IU06,*) ' ------------------------------------------------- '

IF (COARSE) THEN
   WRITE (IU06,*) ' '
   WRITE (IU06,*) ' BOUNDARY OUTPUT IS DONE EVERY.....: ', IDEL_B_OUT,' SECONDS'
   WRITE (IU06,*) ' NEXT OUTPUT DATE IS...............: ', CDT_B_OUT
   WRITE (IU06,*) ' BOUNDARY FILE IS SAVED EVERY......: ', IDEL_BO_FILE,' SECONDS'
   WRITE (IU06,*) ' NEXT DATE TO SAVE BOUNDARY FILE IS: ', CDT_BO_FILE
ELSE
   WRITE(IU06,*) ' COARSE GRID BOUNDAY OUTPUT IS NOT PROCESSED.'
END IF

IF (FINE) THEN
   WRITE (IU06,*) ' '
   WRITE(IU06,*) ' BOUNDAY INPUT IS READ EVERY........: ', IDEL_B_INP,' SECONDS'
   WRITE(IU06,*) ' LAST INPUT DATE IS.................: ', CDATE2
   WRITE(IU06,*) ' NEXT DATE TO FETCH BOUNDARY FILE IS: ', CDT_BI_FILE
   WRITE(IU06,*) ' BOUNDARY FILE IS FETCHED EVERY.....: ', IDEL_BI_FILE,' SECONDS'
ELSE
   WRITE(IU06,*) ' FINE GRID BOUNDAY INPUT IS NOT PROCESSED.'
END IF

END SUBROUTINE PRINT_BOUNDARY_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_BOUNDARY_OUTPUT_TIMESTEPS (STEP_OUTPUT, STEP_FILE)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SET_BOUNDARY_OUTPUT_TIMESTEPS - TRANSFERS TIMESTEPS TO BOUNDARY MODULE.    !
!                                                                              !
!     H. GUENTHER  GKSS         JANUARY 2010                                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        TO DEFINEL THE COARSE GRID BOUNDARY OUTPUT TIMESTEPS.                 !
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
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, OPTIONAL, INTENT(IN) :: STEP_OUTPUT !! TIMESTEP BOUNDARY VALUES.
INTEGER, OPTIONAL, INTENT(IN) :: STEP_FILE   !! TIMESTEP BOUNDARY FILES.

! ---------------------------------------------------------------------------- !

IF (PRESENT(STEP_OUTPUT)) THEN
   IDEL_B_OUT  = STEP_OUTPUT 
ELSE
   IDEL_B_OUT  = -1
END IF
IF (PRESENT(STEP_FILE)) THEN
   IDEL_BO_FILE = STEP_FILE
ELSE
   IDEL_BO_FILE = -1
END IF

END SUBROUTINE SET_BOUNDARY_OUTPUT_TIMESTEPS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SAVE_BOUNDARY_FILE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SAVE_BOUNDARY_FILE - SAVES A BOUNDARY VALUE OUPUT FILE.                    !
!                                                                              !
!     H. GUNTHER    GKSS/ECMWF    OCTOBER 1989                                 !
!     P. JANSSEN    KNMI          OCTOBER 1990   YMP-MODIFICATION              !
!     H. GUNTHER    GKSS/ECMWF    OCTOBER 1990   NEW FILE NAMES.               !
!     H. GUNTHER    GKSS          NOVEMBER 1999  NEW DATES AND FT90.           !
!     E. Myklebust                November 2004  MPI parallelization           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO SAVE THE BOUNDARY VALUE OUTPUT FILE OF A COARSE GRID RUN.           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!        THE ASSIGNED FILE IS CLOSED AND A NEW ONE IS ASSIGNED TO THE UNIT,    !
!        IF THE MODEL DATE IS BEFORE THE END OF RUN DATE.                      !
!        THE NEW FILES ARE OPENED BY SUB. OPEN_FILE.                           !
!                                                                              !
!                                                                              !
!     REFERENCES.                                                              !
!      -----------                                                             !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLE.                                                          !
!     ---------------                                                          !

INTEGER :: IFAIL                               !! OPEN ERROR CODE
INTEGER :: I, IU, ix, len

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CLOSE OLD FILE.                                                       !
!     ------------------                                                       !

IU = IU19-1
DO I = 1, N_NEST
   IU = IU+1
   if (n_nest>1) then
      len = len_trim (file19)
      if (len<=3) then
         write (file19(2:3),'(i2.2)') i
      else
         ix = index (file19,'/',.true.)        !! find last occurence of '/'
         write (file19(ix+2:ix+3),'(i2.2)') i  !! new file name
      endif
   endif
   CLOSE (UNIT=IU, STATUS ="KEEP")             !! BOUNDARY VALUE FILE
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. OPEN A NEW FILE IF MODEL DATE IS BEFORE END DATE.                     !
!        -------------------------------------------------                     !

if (cdtpro<cdatee) then
   CALL INCDATE(CDT_BO_FILE, IDEL_BO_FILE)
   IU = IU19-1
    DO I = 1, N_NEST
       IU = IU+1
       if (n_nest>1) then
          len = len_trim (file19)
          if (len<=3) then
             write (file19(2:3),'(i2.2)') i
          else
             ix = index (file19,'/',.true.)        !! find last occurence of '/'
             write (file19(ix+2:ix+3),'(i2.2)') i  !! new file name
          endif
       endif
       if (n_code(i)==1) then
          call open_file (iu06, iu, file19, cdt_bo_file, 'unknown',            &
&                         ifail, 'formatted')
       else
          CALL OPEN_FILE (IU06, IU, FILE19, CDT_BO_FILE, 'UNKNOWN', IFAIL)
       endif
       IF (IFAIL.NE.0) CALL ABORT1

!     WRITE HEADER (ascii or binary code)

      if (irank==i_out_b_spec) then      
         if (n_code(i)==1) then
            write (iu,*) real(kl), real(ml), th(1), fr(1),                     &
&                        co, real(nbounc(i)), real(idel_b_out),                &
&                        real(idel_bo_file)
            write (iu,*) depth(ijarc(1:nbounc(i),i))
         else
            WRITE (IU) REAL(KL), REAL(ML), TH(1), FR(1), CO, REAL(NBOUNC(I)),  &
&                      REAL(IDEL_B_OUT), REAL(IDEL_BO_FILE)
            write (iu) depth(ijarc(1:nbounc(i),i))
         endif
      endif
   END DO
END IF

END SUBROUTINE SAVE_BOUNDARY_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_BOUNDARY_MODULE

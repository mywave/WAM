MODULE WAM_RADIATION_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS THE ROUTINES AND DATA FOR THE RADIATION STRESS.       !
!   WAVE  FORCE AND STOKES DRIFT.                                              !
!                                                                              !
! ---------------------------------------------------------------------------- !

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  METHODS FROM BASIC MODULES.                                          !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

USE WAM_GENERAL_MODULE,     ONLY:  &
&             ABORT1,              &  !! TERMINATES PROCESSING.
&             INCDATE,             &  !! UPDATES DATE/TIME GROUP.
&             OPEN_FILE,           &  !! OPENS A FILE
&             PRINT_ARRAY             !! PRINTS AN ARRAY

USE WAM_ICE_MODULE,         ONLY:  &
&             PUT_ICE                 !! PUTS ICE INDICATOR INTO DATA FILED.

USE WAM_TOPO_MODULE,        ONLY:  &
&             PUT_DRY                 !! PUTS DRY INDICATOR INTO DATA FILED.

USE WAM_INTERFACE_MODULE,   ONLY:  &
&             STOKES_DRIFT            !! COMPUTES STOKES DRIFT.

USE wam_mpi_comp_module,    ONLY:  &
&             MPI_GATHER_BLOCK,    &
&             MPI_EXCHNG_V

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B.  DATA FROM BASIC MODULES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST, IU27, FILE27
USE WAM_FRE_DIR_MODULE, ONLY: ML, KL, FR, SINTH, COSTH, DFIM, TFAK, TCGOND,    &
&                             DFIMOFR, GOM, C
USE WAM_GENERAL_MODULE, ONLY: G, ZPI, ROWATER
USE WAM_GRID_MODULE,    ONLY: NX, NY, NSEA, AMOWEP, AMOSOP, AMOEAP, AMONOP,    &
&                             L_S_MASK, NLON_RG, ZDELLO,                       &
&                             DELPHI, DELLAM, KLAT, KLON, KXLT
USE WAM_TIMOPT_MODULE,  ONLY: CDATEE, CDTPRO, IDELPRO, SPHERICAL_RUN,          &
&                             SHALLOW_RUN, COLDSTART
USE WAM_MODEL_MODULE,   ONLY: DEPTH, INDEP
USE WAM_OUTPUT_SET_UP_MODULE, ONLY: IDEL_OUT
use wam_propagation_module,   only: dco
USE WAM_NEST_MODULE,          ONLY: FINE, NBOUNF, IJARF

use wam_mpi_module,           only: irank, NINF, NSUP, nijs, nijl, i_out_rad

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C.  MODULE DATA.                                                         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

IMPLICIT NONE
include 'mpif.h'
PRIVATE

CHARACTER (LEN=14) , PARAMETER :: ZERO = ' '
INTEGER :: I

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. FILE AND PRINT OUTPUT.                                                !
!        ----------------------                                                !

INTEGER, PARAMETER :: N_OUT = 8   !! MAXIMUM NUMBER OF OUTPUT PARAMETER.

LOGICAL :: FFLAG(N_OUT) = (/(.FALSE.,I=1,N_OUT)/) !! FILE OUTPUT FLAG FOR EACH
                                                  !! PARAMETER.
LOGICAL :: PFLAG(N_OUT) = (/(.FALSE.,I=1,N_OUT)/) !! PRINT OUTPUT FLAG FOR EACH
                                                  !! PARAMETER.
LOGICAL :: CFLAG(N_OUT) = (/(.FALSE.,I=1,N_OUT)/) !! = FFLAG .OR. PFLAG

INTEGER            :: DEL_RAD_OUT = 0  !! OUTPUT TIME INCREMENT.
CHARACTER (LEN=14) :: CDTOUT = ZERO    !! NEXT DATE TO WRITE RADIATION STRESS.
PUBLIC CDTOUT

INTEGER            :: DELFIL = 0       !! FILE TIME INCREMENT.
CHARACTER (LEN=14) :: CDTFIL = ZERO    !! NEXT DATE TO DISPOSE FILE.

CHARACTER(LEN=60), DIMENSION(N_OUT) :: TITL = (/                     &
& ' RADIATION STRESS TENSOR SXX ( KG/S/S )                     ',    &   !!  1
& ' RADIATION STRESS TENSOR SYY ( KG/S/S )                     ',    &   !!  2
& ' RADIATION STRESS TENSOR SXY ( KG/S/S )                     ',    &   !!  3
& ' DUMMY                                                      ',    &   !!  4
& ' X-COMP. WAVE FORCE PER SURFACE UNIT ( N/M/M )              ',    &   !!  5
& ' Y-COMP. WAVE FORCE PER SURFACE UNIT ( N/M/M )              ',    &   !!  6
& ' X-COMP. STOKES DRIFT ( M/S )                               ',    &   !!  7
& ' Y-COMP. STOKES DRIFT ( M/S )                               '/)       !!  8

REAL, PARAMETER, DIMENSION(N_OUT) :: SCAL = (/                       &
&                      0.            ,    &   !!  1
&                      0.            ,    &   !!  2
&                      0.            ,    &   !!  3
&                      1.            ,    &   !!  4
&                      0.            ,    &   !!  5
&                      0.            ,    &   !!  6
&                      0.            ,    &   !!  7
&                      0.            /)       !!  8

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. DATA FIELDS.                                                          !
!        ------------                                                          !

CHARACTER (LEN=14) :: CDT_RAD = ZERO  !! DATE OF RADIATION STRESSES.

REAL, ALLOCATABLE  :: SXX(:)         !! RADIATION STRESS TENSOR XX
                                     !! AT ALL SEA POINTS.
REAL, ALLOCATABLE  :: SYY(:)         !! RADIATION STRESS TENSOR YY
                                     !! AT ALL SEA POINTS.
REAL, ALLOCATABLE  :: SXY(:)         !! RADIATION STRESS TENSOR XY = YX
                                     !! AT ALL SEA POINTS.
REAL, ALLOCATABLE  :: TAU_X(:)       !! WAVE FORCE X COMPONENT
                                     !! AT ALL SEA POINTS.
REAL, ALLOCATABLE  :: TAU_Y(:)       !! WAVE FORCE Y COMPONENT
                                     !! AT ALL SEA POINTS.
REAL, ALLOCATABLE  :: STOKES_X(:)    !! STOKES DRIFT X COMPONENT
                                     !! AT ALL SEA POINTS.
REAL, ALLOCATABLE  :: STOKES_Y(:)    !! STOKES DRIFT Y COMPONENT
                                     !! AT ALL SEA POINTS.
REAL, SAVE         :: ZMISS=-999999. !! MISSING VALUE (ICE OR DRY POINT)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE RADIATION_STRESS                 !! COMPUTES RATDIATION STRESS
   MODULE PROCEDURE RADIATION_STRESS
END INTERFACE
PUBLIC :: RADIATION_STRESS

INTERFACE SET_RADIATION_TIMES              !! SETS TIMES FOR COMPUTATION.
   MODULE PROCEDURE SET_RADIATION_TIMES
END INTERFACE
PUBLIC SET_RADIATION_TIMES

INTERFACE SET_RADIATION_FILE               !! SETS FILE FOR RADIATION OUTPUT.
   MODULE PROCEDURE SET_RADIATION_FILE
END INTERFACE
PUBLIC SET_RADIATION_FILE

INTERFACE PRINT_RADIATION_MODULE           !! PRINTS RADIATION SETTING.
   MODULE PROCEDURE PRINT_RADIATION_MODULE
END INTERFACE
PUBLIC PRINT_RADIATION_MODULE

INTERFACE PREPARE_RADIATION                !! PREPARES RADIATION MODULE.
   MODULE PROCEDURE PREPARE_RADIATION
END INTERFACE
PUBLIC PREPARE_RADIATION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTERFACE GRADIENT_RAD                  !! CALCULATES GRADIENTS.
   MODULE PROCEDURE GRADIENT_RAD
END INTERFACE

INTERFACE SAVE_RADIATION_FILE           !! SAVES AND OPENS OUTPUT FILE.
   MODULE PROCEDURE SAVE_RADIATION_FILE
END INTERFACE

INTERFACE RADIATION_OUTPUT              !! PRINTS GRIDDED STRESSES.
   MODULE PROCEDURE RADIATION_OUTPUT
END INTERFACE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE RADIATION_STRESS (FL3)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   RADIATION_STRESS - COMPUTAION OF RADIATION STRESS                          !
!                                                                              !
!     H. GUNTHER         GKSS/ECMWF         JUNE 1990                          !
!     C.SCHNEGGENBURGER  GKSS MODIFICATIONS FOR K-MODEL JUNE 1996              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE AND CONTROL THE RADIATION STRESS COMPUTATIONS.              !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !
!                                                                              !
!       GRADIENTS              - COMPUTES GRADIENTS.                           !
!       PRINT_RADIATION_GRID   - PRINTER OUTPUT.                               !
!       UPDATE_RADIATION_TIMES - NEXT RADIATION STRESS TIME.                   !
!       SAVE_RADIATION_FILE    - NEW OUTPUT FILE.                              !
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
!     INTERFACE VARIABLE.                                                      !
!     -------------------                                                      !

REAL,    INTENT(IN) :: FL3(nijs:nijl,KL,ML)  !! BLOCK OF SPECTRA.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, K, M
REAL    :: CGDCS
REAL, ALLOCATABLE, DIMENSION(:) :: F2
REAL, ALLOCATABLE, DIMENSION(:) :: CGDC  

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. RADIATION STRESS TENSOR.                                              !
!        --------------------------                                            !

CDT_RAD = CDTPRO

IF (ANY(CFLAG(1:6))) THEN
   ALLOCATE (SXX(NINF:NSUP))   
   ALLOCATE (SYY(NINF:NSUP))   
   ALLOCATE (SXY(NINF:NSUP))   
   ALLOCATE (F2(NIJS:NIJL))
   ALLOCATE (CGDC(NIJS:NIJL))

   F2 = 0.
   SXX = 0.
   SYY = 0.
   SXY = 0.
     
   IF(SHALLOW_RUN) THEN

      DO M = 1,ML
         CGDC(NIJS:NIJL) = TCGOND(INDEP(NIJS:NIJL),M)                          &
&                        * TFAK(INDEP(NIJS:NIJL),M) * DFIMOFR(M)/ZPI
                                             !! GROUP VELOCITY / PHASE VELOCITY
         DO K = 1,KL
            DO IJ = NIJS,NIJL
               SXX(IJ) = SXX(IJ) + SINTH(K)*SINTH(K)*CGDC(IJ)*FL3(IJ,K,M)
               SYY(IJ) = SYY(IJ) + COSTH(K)*COSTH(K)*CGDC(IJ)*FL3(IJ,K,M)
               SXY(IJ) = SXY(IJ) + SINTH(K)*COSTH(K)*CGDC(IJ)*FL3(IJ,K,M)
               F2 (IJ) = F2 (IJ) + (CGDC(IJ)-0.5*DFIM(M))*FL3(IJ,K,M)
            END DO
         END DO
      END DO
   ELSE
      DO M = 1,ML
         CGDCS = GOM(M) / C(M) * DFIM(M)   !! GROUP VELOCITY / PHASE VELOCITY
         DO K = 1,KL
            DO IJ = NIJS,NIJL
               SXX(IJ) = SXX(IJ) + SINTH(K)*SINTH(K)*CGDCS*FL3(IJ,K,M)
               SYY(IJ) = SYY(IJ) + COSTH(K)*COSTH(K)*CGDCS*FL3(IJ,K,M)
               SXY(IJ) = SXY(IJ) + SINTH(K)*COSTH(K)*CGDCS*FL3(IJ,K,M)
               F2 (IJ) = F2 (IJ) + (CGDCS -0.5*DFIM(M))*FL3(IJ,K,M)
            END DO
         END DO
      END DO
   END IF

   DO IJ = NIJS,NIJL
      SXX(IJ) = ROWATER*G*(SXX(IJ) + F2(IJ))
      SYY(IJ) = ROWATER*G*(SYY(IJ) + F2(IJ))
      SXY(IJ) = ROWATER*G*SXY(IJ)
   END DO

   CALL PUT_DRY(SXX(nijs:nijl), NIJS, NIJL, ZMISS)   !! FLAG DRY POINTS
   CALL PUT_DRY(SYY(nijs:nijl), NIJS, NIJL, ZMISS)
   CALL PUT_DRY(SXY(nijs:nijl), NIJS, NIJL, ZMISS)
   CALL PUT_ICE(SXX(nijs:nijl), NIJS, NIJL, ZMISS)   !! FLAG ICE POINTS
   CALL PUT_ICE(SYY(nijs:nijl), NIJS, NIJL, ZMISS)
   CALL PUT_ICE(SXY(nijs:nijl), NIJS, NIJL, ZMISS)
    
   IF (FINE) THEN              !! FLAG BOUNDARY INPUT POINTS (IF FINE GRID)
      DO I = 1, NBOUNF
         IF (IJARF(I).LT.NIJS .OR. IJARF(I).GT.NIJL) CYCLE
         SXX(IJARF(I)) = ZMISS
         SYY(IJARF(I)) = ZMISS
         SXY(IJARF(I)) = ZMISS
      END DO
   END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. WAVE FORCE VECTOR.                                                     !
!       ------------------                                                     !

   IF (ANY(CFLAG(5:6))) THEN

      ALLOCATE (TAU_X(nijs:nijl))   
      ALLOCATE (TAU_Y(nijs:nijl))   

      call mpi_exchng_V (SXX(ninf:nsup)) !! TOP AND BOTTOM FROM OTHER PROCESSES
      call mpi_exchng_V (SYY(ninf:nsup))
      call mpi_exchng_V (SXY(ninf:nsup))

      CALL GRADIENT_RAD (FIELD=SXX, D_LON=TAU_X)
      CALL GRADIENT_RAD (FIELD=SYY, D_LAT=TAU_Y)
      CALL GRADIENT_RAD (FIELD=SXY, D_LAT=CGDC,  D_LON=F2)

      WHERE (TAU_X.NE.ZMISS .AND. CGDC.NE.ZMISS) 
         TAU_X = -(TAU_X + CGDC )/ROWATER
      ELSEWHERE
         TAU_X = ZMISS
      END WHERE
      WHERE (TAU_Y.NE.ZMISS .AND. F2.NE.ZMISS) 
         TAU_Y = -(F2    + TAU_Y)/ROWATER
      ELSEWHERE
         TAU_Y = ZMISS
      END WHERE

      CALL PUT_DRY(TAU_X, NIJS, NIJL, ZMISS)
      CALL PUT_DRY(TAU_Y, NIJS, NIJL, ZMISS)
      CALL PUT_ICE(TAU_X, NIJS, NIJL, ZMISS)
      CALL PUT_ICE(TAU_Y, NIJS, NIJL, ZMISS)
   END IF

   IF (ALLOCATED (F2  )) DEALLOCATE (F2)
   IF (ALLOCATED (CGDC)) DEALLOCATE (CGDC)

END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTE STOKES DRIFT.                                                 !
!        ---------------------                                                 !

IF (ANY(CFLAG(7:8))) THEN
   ALLOCATE (STOKES_X(nijs:nijl))   
   ALLOCATE (STOKES_Y(nijs:nijl))   

   IF (SHALLOW_RUN) THEN
      CALL STOKES_DRIFT (FL3, STOKES_X(NIJS:NIJL), STOKES_Y(NIJS:NIJL),        &
&                        INDEP(NIJS:NIJL))
   ELSE
      CALL STOKES_DRIFT (FL3, STOKES_X(NIJS:NIJL), STOKES_Y(NIJS:NIJL))
   END IF

   CALL PUT_DRY(STOKES_X, NIJS, NIJL, ZMISS)   !! FLAG DRY POINTS
   CALL PUT_DRY(STOKES_Y, NIJS, NIJL, ZMISS)
   CALL PUT_ICE(STOKES_X, NIJS, NIJL, ZMISS)   !! FLAG ICE POINTS
   CALL PUT_ICE(STOKES_Y, NIJS, NIJL, ZMISS)
     
   IF (FINE) THEN              !! FLAG BOUNDARY INPUT POINTS (IF FINE GRID)
      DO I = 1, NBOUNF
         IF (IJARF(I).LT.NIJS .OR. IJARF(I).GT.NIJL) CYCLE
         STOKES_X(IJARF(I)) = ZMISS
         STOKES_Y(IJARF(I)) = ZMISS
      END DO
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. OUTPUT RADIATION STRESS.                                              !
!        ------------------------                                              !

CALL RADIATION_OUTPUT
IF (ITEST.GE.4) THEN
   WRITE(IU06,*) '       SUB. RADIATION_STRESS: RADIATION_OUTPUT DONE.'
END IF

IF (ALLOCATED(SXX     )) DEALLOCATE(SXX)
IF (ALLOCATED(SYY     )) DEALLOCATE(SYY)
if (allocated(sxy     )) deallocate(sxy)
IF (ALLOCATED(TAU_X   )) DEALLOCATE(TAU_X)
IF (ALLOCATED(TAU_Y   )) DEALLOCATE(TAU_Y)
IF (ALLOCATED(STOKES_X)) DEALLOCATE(STOKES_X)
IF (ALLOCATED(STOKES_Y)) DEALLOCATE(STOKES_Y)

END SUBROUTINE RADIATION_STRESS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_RADIATION_TIMES (OUT_DEL, PF, FF)

INTEGER,           INTENT(IN) :: OUT_DEL    !! TIME INCREMENT FOR
                                            !! RADIATION STRESS OUTPUT.

LOGICAL, OPTIONAL, INTENT(IN) :: PF(N_OUT)  !! .TRUE. IF PRINTER OUTPUT.
LOGICAL, OPTIONAL, INTENT(IN) :: FF(N_OUT)  !! .TRUE. IF FILE OUTPUT.

! ---------------------------------------------------------------------------- !

DEL_RAD_OUT = MAX(OUT_DEL, 0)

IF (PRESENT(PF)) PFLAG = PF
IF (PRESENT(FF)) FFLAG = FF

PFLAG(4) = .FALSE.
FFLAG(4) = .FALSE.
CFLAG = FFLAG .OR. PFLAG

END SUBROUTINE SET_RADIATION_TIMES

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SET_RADIATION_FILE (FILE_INC, NAME, UNIT)

INTEGER,   INTENT(IN), OPTIONAL :: FILE_INC     !! TIME INCREMENT TO SAVE
                                                !! OUTPUT FILE. IF .LE. 0. 
                                                !! FILE IS SAVED AT 
                                                !! THE END OF THE RUN
CHARACTER,  INTENT(IN), OPTIONAL :: NAME*(*)    !! RADIATION STRESS FILE NAME
INTEGER,    INTENT(IN), OPTIONAL :: UNIT        !! LOGICAL FILE UNIT NO.

! ---------------------------------------------------------------------------- !

DELFIL   = 0
IF (PRESENT(FILE_INC)) DELFIL = MAX(FILE_INC, 0)

FILE27 = 'RAD'
IF (PRESENT(NAME)) THEN
   IF (LEN_TRIM(NAME).GT.0) FILE27 = TRIM(NAME)
END IF

IU27 = 27
IF (PRESENT(UNIT)) THEN
   IF (UNIT.GT.0) IU27 = UNIT
END IF

END SUBROUTINE SET_RADIATION_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_RADIATION_MODULE

INTEGER :: I

WRITE (IU06,*) '  '
WRITE (IU06,*) ' ------------------------------------------------'
WRITE (IU06,*) '    MODEL SPECIFICATIONS FOR RADIATION STRESS:'
WRITE (IU06,*) ' ------------------------------------------------'
WRITE (IU06,*) '  '

IF (DEL_RAD_OUT.LE.0.) THEN
   WRITE(IU06,*) ' OUTPUT OF RADIATION STRESS IS NOT REQUESTED '
   RETURN
END IF

WRITE(IU06,*) ' TO PRINTER AND/OR FILE : ', TRIM(FILE27),'YYYYMMDDHHMMSS'
WRITE(IU06,*) ' OUTPUT OF RADIATION STRESS EVERY ..: ', DEL_RAD_OUT, ' [S]'
WRITE(IU06,*) ' THE NEXT OUTPUT DATE IS ...........: ', CDTOUT
IF (DELFIL.GT.0.) THEN
   WRITE(IU06,*) ' A NEW FILE WILL BE USED EVERY .....: ', DELFIL, ' [S]'
   WRITE(IU06,*) ' THE PRESENT FILE DATE .............: ', CDTFIL
ELSE
   WRITE(IU06,*) ' THE FILE DATE IS AT END OF RUN.....: ', CDTFIL
END IF

WRITE(IU06,*) '                                                           ',   &
&                                                        'PRINTER     UNIT'
DO I=1,N_OUT
   WRITE(IU06,*) TITL(I),'....', PFLAG( I),'......', FFLAG( I)
END DO

WRITE(IU06,*) '  '
IF (CDT_RAD.NE.ZERO) THEN
   WRITE(IU06,*) ' DATE OF RADIATION STRESSES STORED IN MODULE IS: ', CDT_RAD
ELSE
   WRITE(IU06,*) ' RADIATION STRESSES ARE NOT STORED IN MODULE'
END IF
WRITE(IU06,*) '  '

END SUBROUTINE PRINT_RADIATION_MODULE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_RADIATION (FL3)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLE.                                                      !
!     -------------------                                                      !

REAL,    INTENT(IN) :: FL3(:,:,:)  !! BLOCK OF SPECTRA.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. OUTPUT REQUESTED?                                                     !
!        ------------------                                                    !

IF (.NOT.ANY(CFLAG)) THEN
   CDT_RAD = ' '
   DEL_RAD_OUT = -1
   PFLAG = .FALSE.
   FFLAG = .FALSE.
   CFLAG = .FALSE.
   RETURN
ELSE
   IF (DEL_RAD_OUT.LE. 0) DEL_RAD_OUT = IDELPRO     
   IF (DELFIL     .LE. 0) DELFIL      = MAX(IDEL_OUT,0)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. CHECK INCREMENTS.                                                     !
!        -----------------                                                     !

IF (((DEL_RAD_OUT/IDELPRO)*IDELPRO).NE.DEL_RAD_OUT) THEN
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                                       +'
   WRITE(IU06,*) '+     WARNING ERROR IN SUB. PREPARE_RADIATION           +'
   WRITE(IU06,*) '+     =======================================           +'
   WRITE(IU06,*) '+                                                       +'
   WRITE(IU06,*) '+  OUTPUT INCREMENT OF RADIATION STRESS IS NOT A        +'
   WRITE(IU06,*) '+  MULTIPLE OF PROPAGATION INCREMENT.                   +'
   WRITE(IU06,*) '+  OUTPUT INCREMENT CHANGED.                            +'
   WRITE(IU06,*) '+  OLD OUTPUT INCREMENT WAS : ' ,DEL_RAD_OUT
   DEL_RAD_OUT = MAX(((DEL_RAD_OUT/IDELPRO)*IDELPRO), IDELPRO)
   WRITE(IU06,*) '+  NEW OUTPUT INCREMENT IS  : ' ,DEL_RAD_OUT
   WRITE(IU06,*) '+                                                       +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF   

IF (((DELFIL/DEL_RAD_OUT)*DEL_RAD_OUT).NE.DELFIL) THEN
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) '+                                                       +'
   WRITE(IU06,*) '+     WARNING ERROR IN SUB. PREPARE_RADIATION           +'
   WRITE(IU06,*) '+     =======================================           +'
   WRITE(IU06,*) '+                                                       +'
   WRITE(IU06,*) '+  FILE INCREMENT OF RADIATION OUTPUT IS NOT A          +'
   WRITE(IU06,*) '+  MULTIPLE OF OUTPUT INCREMENT.                        +'
   WRITE(IU06,*) '+  FILE INCREMENT CHANGED.                              +'
   WRITE(IU06,*) '+  OLD FILE INCREMENT WAS : ' ,DELFIL
   DELFIL = MAX(((DELFIL/DEL_RAD_OUT)*DEL_RAD_OUT), DEL_RAD_OUT)
   WRITE(IU06,*) '+  NEW FILE INCREMENT IS  : ' ,DELFIL
   WRITE(IU06,*) '+                                                       +'
   WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. TIME COUNTER FOR NEXT OUTPUT.                                         !
!        -----------------------------                                         !

CDTOUT = ' '
IF (DEL_RAD_OUT .GT. 0.) THEN
   CDTOUT = CDTPRO
   IF (.NOT.COLDSTART) CALL INCDATE (CDTOUT, DEL_RAD_OUT)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. TIME COUNTER FOR NEXT OUTPUT FILES AND OPEN FIRST OUTPUT FILE.        !
!        --------------------------------------------------------------        !

CDTFIL = ' '
IF (ANY(FFLAG)) THEN
   IF (DELFIL.EQ.0.) THEN
      CDTFIL = CDATEE
   ELSE
      CDTFIL = CDTPRO
      IF (COLDSTART)  CALL INCDATE (CDTFIL, -DELFIL)
   END IF
   CALL SAVE_RADIATION_FILE
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. COMPUTE FIRST RADIATION STRESS.                                       !
!        -------------------------------                                       !

IF (CDTPRO.EQ.CDTOUT) THEN
   CALL RADIATION_STRESS (FL3)
   IF (ITEST.GE.2) THEN
      WRITE (IU06,*) '   SUB. PREPARE_RADIATION: RADIATION_STRESS DONE '
   END IF
END IF

END SUBROUTINE PREPARE_RADIATION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE GRADIENT_RAD (FIELD, D_LAT, D_LON)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   GRADIENT_RAD - CALCULATES GRADIENT.                                        !
!                                                                              !
!     K.P. HUBBERT              AUGUST   1988                                  !
!     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.           !
!     H. GUNTHER    HZG         DECEMBER 2010  MISSING VALUES.                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO CALCULATE THE GRADIENT VECTOR FIELD FROM INPUT FIELD.               !
!                                                                              !
!     METHOD.                                                                  !
!     ------                                                                   !
!                                                                              !
!       CENTRAL DIFFERENCING. IF A VALUE IS MISSING TRY FIRST ORDER DIFFENCES  !
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

REAL,           INTENT(IN)  :: FIELD(NINF:NSUP)  !! INPUT FIELD.
REAL, OPTIONAL, INTENT(OUT) :: D_LAT(nijs:nijl)  !! LATITUDE  DERIVATIVE.
REAL, OPTIONAL, INTENT(OUT) :: D_LON(nijs:nijl)  !! LONGITUDE DERIVATIVE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, IJ1, IJ2

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NORTH SOUTH GRADIENTS.                                                !
!        ----------------------                                                !

IF (PRESENT(D_LAT)) THEN
   D_LAT = ZMISS

   DO IJ = nijs,nijl
      IF (FIELD(IJ).EQ.ZMISS) CYCLE
      IJ1 = KLAT(IJ,1)
      IJ2 = KLAT(IJ,2)
      IF (IJ1.LE.ninf-1 .AND. IJ2.LE.ninf-1) CYCLE

      IF (IJ1.GT.ninf-1 .AND. IJ2.GT.ninf-1) THEN
         IF (FIELD(IJ2).NE.ZMISS .AND. FIELD(IJ1).NE.ZMISS) THEN
            D_LAT(IJ) = (FIELD(IJ2)-FIELD(IJ1))/(2.*DELPHI)
      ELSE IF (FIELD(IJ1).NE.ZMISS .AND. FIELD(IJ2).EQ.ZMISS) THEN
            D_LAT(IJ) = (FIELD(IJ)-FIELD(IJ1))/DELPHI
         ELSE IF (FIELD(IJ1).EQ.ZMISS .AND. FIELD(IJ2).NE.ZMISS) THEN
            D_LAT(IJ) = (FIELD(IJ2)-FIELD(IJ))/DELPHI
         END IF

      ELSE IF (IJ1.GT.ninf-1) THEN
         IF (FIELD(IJ1).NE.ZMISS) THEN
            D_LAT(IJ) = (FIELD(IJ)-FIELD(IJ1))/DELPHI
         END IF

      ELSE IF (IJ2.GT.ninf-1) THEN
         IF (FIELD(IJ2).NE.ZMISS) THEN
            D_LAT(IJ) = (FIELD(IJ2)-FIELD(IJ))/DELPHI
         END IF
      END IF

   END DO
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. EAST WEST GRADIENTS.                                                  !
!        --------------------                                                  !

IF (PRESENT(D_LON)) THEN
   D_LON = ZMISS

   DO IJ = nijs,nijl
      IF (FIELD(IJ).EQ.ZMISS) CYCLE
      IJ1 = KLON(IJ,1)
      IJ2 = KLON(IJ,2)
      IF (IJ1.LE.ninf-1 .AND. IJ2.LE.ninf-1) CYCLE
 
      IF (IJ1.GT.ninf-1 .AND. IJ2.GT.ninf-1) THEN
         IF (FIELD(IJ2).NE.ZMISS .AND. FIELD(IJ1).NE.ZMISS) THEN
            D_LON(IJ) = (FIELD(IJ2)-FIELD(IJ1))/(2.*DELLAM(KXLT(IJ)))
         ELSE IF (FIELD(IJ1).NE.ZMISS .AND. FIELD(IJ2).EQ.ZMISS) THEN
            D_LON(IJ) = (FIELD(IJ)-FIELD(IJ1))/DELLAM(KXLT(IJ))
         ELSE IF (FIELD(IJ1).EQ.ZMISS .AND. FIELD(IJ2).NE.ZMISS) THEN
            D_LON(IJ) = (FIELD(IJ2)-FIELD(IJ))/DELLAM(KXLT(IJ))
         END IF

      ELSE IF (IJ1.GT.ninf-1) THEN
         IF (FIELD(IJ1).NE.ZMISS) THEN
            D_LON(IJ) = (FIELD(IJ)-FIELD(IJ1))/DELLAM(KXLT(IJ))
         END IF

      ELSE IF (IJ2.GT.ninf-1) THEN
         IF (FIELD(IJ2).NE.ZMISS) THEN
            D_LON(IJ) = (FIELD(IJ2)-FIELD(IJ))/DELLAM(KXLT(IJ))
         END IF
      END IF
      IF (SPHERICAL_RUN.and.d_lon(ij).ne.zmiss) D_LON(IJ) = D_LON(IJ)*DCO(IJ)

   END DO
END IF

END SUBROUTINE GRADIENT_RAD

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SAVE_RADIATION_FILE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SAVE_RADIATION_FILE - TO SAVES RADIATION OUTPUT FILES.                     !
!                                                                              !
!     HEINZ GUNTHER       GKSS  JANUARY   2005                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO CLOSE ASSIGNED RADIATION OUTPUT FILE AND OPEN NEW A FILE            !
!       IF NECCESSARY.                                                         !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !
!                                                                              !
!       ABORT1    - TERMINATES PROCESSING.                                     !
!       OPEN_FILE - OPEN A NEW FILE.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        ASSIGNED FILES ARE CLOSED AND NEW FILES ARE OPENED IF THE             !
!        MODEL RUN IS NOT FINISHED.                                            !
!                                                                              !
!        THE NEW FILE DATE IS INCEREMENTED.                                    !
!        THE FILE NAME CONVENTION IS EXPLAINED IN SUB OPENFIL.                 !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTEGER           :: IFAIL
CHARACTER*7, SAVE :: STAT='UNKNOWN'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. SAVE THE OLD RADIATION FILE.                                          !
!        ----------------------------                                          !

CLOSE(UNIT=IU27,STATUS='KEEP')

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. OPEN A NEW RADIATION FILE.                                            !
!        --------------------------                                            !

IF (CDTPRO.LT.CDATEE) THEN
   CALL INCDATE (CDTFIL, DELFIL)
   CALL OPEN_FILE (IU06, IU27, FILE27, CDTFIL, STAT, IFAIL)
   IF (IFAIL.NE.0) CALL ABORT1
END IF

END SUBROUTINE SAVE_RADIATION_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE RADIATION_OUTPUT

integer :: ierr
REAL, ALLOCATABLE, DIMENSION(:,:) :: GRID
REAL, ALLOCATABLE, DIMENSION(:)   :: BLOCK_TOTAL

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. WRITE OUTPUT.                                                         !
!        -------------                                                         !

IF (irank==i_out_rad) then
   ALLOCATE (GRID(NX,NY))
   ALLOCATE (BLOCK_TOTAL(1:nsea))
   IF (ANY(FFLAG)) THEN
      WRITE (IU27) CDT_RAD, REAL(NX), REAL(NY), AMOWEP, AMOSOP, AMOEAP, AMONOP
      WRITE (IU27) NLON_RG, ZDELLO
      WRITE (IU27) FFLAG
   END IF
END IF

IF (PFLAG(1) .OR. FFLAG(1)) THEN
   if (irank==i_out_rad) then
      CALL mpi_gather_block(i_out_rad, SXX(NIJS:NIJL), BLOCK_TOTAL)
   else
      CALL mpi_gather_block(i_out_rad, SXX(NIJS:NIJL))
   end if
   call mpi_barrier (mpi_comm_world, ierr)

   if (irank==i_out_rad) then
      GRID = UNPACK (BLOCK_TOTAL, L_S_MASK, ZMISS)
      IF (FFLAG(1)) WRITE (IU27) GRID
      IF (PFLAG(1)) CALL PRINT_ARRAY (IU06, CDT_RAD, TITL(1), GRID,          &
&                        AMOWEP, AMOSOP, AMOEAP, AMONOP, SCAL(1), ZMISS)
   END IF
END IF

IF (PFLAG(2) .OR. FFLAG(2)) THEN
   if (irank==i_out_rad) then
      CALL mpi_gather_block(i_out_rad, SYY(NIJS:NIJL), BLOCK_TOTAL)
   else
      CALL mpi_gather_block(i_out_rad, SYY(NIJS:NIJL))
   end if
   call mpi_barrier (mpi_comm_world, ierr)

   if (irank==i_out_rad) then
      GRID = UNPACK (BLOCK_TOTAL, L_S_MASK, ZMISS)
      IF (FFLAG(2)) WRITE (IU27) GRID
      IF (PFLAG(2)) CALL PRINT_ARRAY (IU06, CDT_RAD, TITL(2), GRID,          &
&                        AMOWEP, AMOSOP, AMOEAP, AMONOP, SCAL(2), ZMISS)
   END IF
END IF

IF (PFLAG(3) .OR. FFLAG(3)) THEN
   if (irank==i_out_rad) then
      CALL mpi_gather_block(i_out_rad, SXY(NIJS:NIJL), BLOCK_TOTAL)
   else
      CALL mpi_gather_block(i_out_rad, SXY(NIJS:NIJL))
   end if
   call mpi_barrier (mpi_comm_world, ierr)

   if (irank==i_out_rad) then
      GRID = UNPACK (BLOCK_TOTAL, L_S_MASK, ZMISS)
      IF (FFLAG(3)) WRITE (IU27) GRID
      IF (PFLAG(3)) CALL PRINT_ARRAY (IU06, CDT_RAD, TITL(3), GRID,          &
&                        AMOWEP, AMOSOP, AMOEAP, AMONOP, SCAL(3), ZMISS)
   END IF
END IF

IF (PFLAG(5) .OR. FFLAG(5)) THEN
   if (irank==i_out_rad) then
      CALL mpi_gather_block(i_out_rad, TAU_X(NIJS:NIJL), BLOCK_TOTAL)
   else
      CALL mpi_gather_block(i_out_rad, TAU_X(NIJS:NIJL))
   end if
   call mpi_barrier (mpi_comm_world, ierr)

   if (irank==i_out_rad) then
      GRID = UNPACK (BLOCK_TOTAL, L_S_MASK, ZMISS)
      IF (FFLAG(5)) WRITE (IU27) GRID
      IF (PFLAG(5)) CALL PRINT_ARRAY (IU06, CDT_RAD, TITL(5), GRID,          &
&                        AMOWEP, AMOSOP, AMOEAP, AMONOP, SCAL(5), ZMISS)
   END IF
END IF
   
IF (PFLAG(6) .OR. FFLAG(6)) THEN
   if (irank==i_out_rad) then
      CALL mpi_gather_block(i_out_rad, TAU_Y(NIJS:NIJL), BLOCK_TOTAL)
   else
      CALL mpi_gather_block(i_out_rad, TAU_Y(NIJS:NIJL))
   end if
   call mpi_barrier (mpi_comm_world, ierr)

   if (irank==i_out_rad) then
      GRID = UNPACK (BLOCK_TOTAL, L_S_MASK, ZMISS)
      IF (FFLAG(6)) WRITE (IU27) GRID
      IF (PFLAG(6)) CALL PRINT_ARRAY (IU06, CDT_RAD, TITL(6), GRID,          &
&                        AMOWEP, AMOSOP, AMOEAP, AMONOP, SCAL(6), ZMISS)
   END IF
END IF
   
IF (PFLAG(7) .OR. FFLAG(7)) THEN
   if (irank==i_out_rad) then
      CALL mpi_gather_block(i_out_rad, STOKES_X(NIJS:NIJL), BLOCK_TOTAL)
   else
      CALL mpi_gather_block(i_out_rad, STOKES_X(NIJS:NIJL))
   end if
   call mpi_barrier (mpi_comm_world, ierr)

   if (irank==i_out_rad) then
      GRID = UNPACK (BLOCK_TOTAL, L_S_MASK, ZMISS)
      IF (FFLAG(7)) WRITE (IU27) GRID
      IF (PFLAG(7)) CALL PRINT_ARRAY (IU06, CDT_RAD, TITL(7), GRID,          &
&                        AMOWEP, AMOSOP, AMOEAP, AMONOP, SCAL(7), ZMISS)
   END IF
END IF

IF (PFLAG(8) .OR. FFLAG(8)) THEN
   if (irank==i_out_rad) then
      CALL mpi_gather_block(i_out_rad, STOKES_Y(NIJS:NIJL), BLOCK_TOTAL)
   else
      CALL mpi_gather_block(i_out_rad, STOKES_Y(NIJS:NIJL))
   end if
   call mpi_barrier (mpi_comm_world, ierr)

   if (irank==i_out_rad) then
      GRID = UNPACK (BLOCK_TOTAL, L_S_MASK, ZMISS)
      IF (FFLAG(8)) WRITE (IU27) GRID
      IF (PFLAG(8)) CALL PRINT_ARRAY (IU06, CDT_RAD, TITL(8), GRID,          &
&                        AMOWEP, AMOSOP, AMOEAP, AMONOP, SCAL(8), ZMISS)
   END IF
END IF

IF (ALLOCATED(BLOCK_TOTAL)) DEALLOCATE(BLOCK_TOTAL)
IF (ALLOCATED(GRID)) DEALLOCATE(GRID)


! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. UP-DATE OUTPUT TIME.                                                  !
!        --------------------                                                  !

CALL INCDATE (CDTOUT, DEL_RAD_OUT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. SAVE OUTPUT FILE.                                                      !
!        ----------------                                                      !

IF (CDTFIL.EQ.CDTPRO) CALL SAVE_RADIATION_FILE

END SUBROUTINE RADIATION_OUTPUT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_RADIATION_MODULE

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

USE WAM_ICE_MODULE,         ONLY:  &
&             PUT_ICE                 !! PUTS ICE INDICATOR INTO DATA FILED.

USE WAM_TOPO_MODULE,        ONLY:  &
&             PUT_DRY                 !! PUTS DRY INDICATOR INTO DATA FILED.

USE wam_mpi_comp_module,    ONLY:  &
&             MPI_EXCHNG

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B.  DATA FROM BASIC MODULES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

USE WAM_FILE_MODULE,    ONLY: IU06, ITEST
USE WAM_FRE_DIR_MODULE, ONLY: ML, KL, FR, SINTH, COSTH, DFIM, DFIMOFR, GOM, C
USE WAM_GENERAL_MODULE, ONLY: G, ZPI, ROWATER
USE WAM_GRID_MODULE,    ONLY: NX, NY, NSEA, AMOWEP, AMOSOP, AMOEAP, AMONOP,    &
&                             L_S_MASK, NLON_RG, ZDELLO, DELPHI, DELLAM,       &
&                             KLAT, KLON, KXLT, WLAT, REDUCED_GRID
USE WAM_TIMOPT_MODULE,  ONLY: CDATEE, CDTPRO, IDELPRO, SPHERICAL_RUN,          &
&                             SHALLOW_RUN, COLDSTART, LCFLX
USE WAM_MODEL_MODULE,   ONLY: DEPTH, INDEP
use wam_propagation_module,   only: dco
USE WAM_NEST_MODULE,          ONLY: FINE, NBOUNF, IJARF
USE WAM_TABLES_MODULE,        ONLY: TFAK, TCGOND

use wam_mpi_module,           only: irank, NINF, NSUP, nijs, nijl
USE WAM_OUTPUT_SET_UP_MODULE, ONLY: CFLAG_P, ZMISS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C.  MODULE DATA.                                                         !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

IMPLICIT NONE
include 'mpif.h'
PRIVATE

INTEGER :: I

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DATA FIELDS.                                                          !
!        ------------                                                          !

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

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE RADIATION_STRESS                 !! COMPUTES RATDIATION STRESS
   MODULE PROCEDURE RADIATION_STRESS
END INTERFACE
PUBLIC :: RADIATION_STRESS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

INTERFACE GRADIENT_RAD                  !! CALCULATES GRADIENTS.
   MODULE PROCEDURE GRADIENT_RAD
END INTERFACE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE RADIATION_STRESS (FL3, BLOCK)

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
!       TO COMPUTE AND CONTROL THE RADIATION STRESS                            !
!       AND WAVE FORCE COMPUTATIONS.                                           !
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

REAL,    INTENT(IN)  :: FL3(NIJS:NIJL,KL,ML)  !! BLOCK OF SPECTRA.
REAL,    INTENT(OUT) :: BLOCK(NIJS:NIJL,51:56)  !! OUTPUT PARAMETERS.

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

   IF (CFLAG_P(51)) BLOCK(NIJS:NIJL,51) = SXX(NIJS:NIJL)
   IF (CFLAG_P(52)) BLOCK(NIJS:NIJL,52) = SYY(NIJS:NIJL)
   IF (CFLAG_P(53)) BLOCK(NIJS:NIJL,53) = SXY(NIJS:NIJL)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. WAVE FORCE VECTOR.                                                     !
!       ------------------                                                     !

   IF (ANY(CFLAG_P(55:56))) THEN

      ALLOCATE (TAU_X(nijs:nijl))   
      ALLOCATE (TAU_Y(nijs:nijl))   

      call mpi_exchng (SXX(ninf:nsup)) !! TOP AND BOTTOM FROM OTHER PROCESSES
      call mpi_exchng (SYY(ninf:nsup))
      call mpi_exchng (SXY(ninf:nsup))

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

      IF (CFLAG_P(55)) BLOCK(NIJS:NIJL,55) = TAU_X(NIJS:NIJL)
      IF (CFLAG_P(56)) BLOCK(NIJS:NIJL,56) = TAU_Y(NIJS:NIJL)
   END IF

   IF (ALLOCATED (F2  )) DEALLOCATE (F2)
   IF (ALLOCATED (CGDC)) DEALLOCATE (CGDC)



! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. OUTPUT RADIATION STRESS.                                              !
!        ------------------------                                              !

IF (ITEST.GE.4) THEN
   WRITE(IU06,*) '       SUB. RADIATION_STRESS DONE.'
END IF

IF (ALLOCATED(SXX     )) DEALLOCATE(SXX)
IF (ALLOCATED(SYY     )) DEALLOCATE(SYY)
if (allocated(sxy     )) deallocate(sxy)
IF (ALLOCATED(TAU_X   )) DEALLOCATE(TAU_X)
IF (ALLOCATED(TAU_Y   )) DEALLOCATE(TAU_Y)

END SUBROUTINE RADIATION_STRESS

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
REAL, OPTIONAL, INTENT(OUT) :: D_LAT(NIJS:NIJL)  !! LATITUDE  DERIVATIVE.
REAL, OPTIONAL, INTENT(OUT) :: D_LON(NIJS:NIJL)  !! LONGITUDE DERIVATIVE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: NLAND, IJ, IP, IM, IP2, IM2, KX
REAL    :: D_LAT_P, D_LAT_M

NLAND = NINF-1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NORTH SOUTH GRADIENTS.                                                !
!        ----------------------                                                !

IF (PRESENT(D_LAT)) THEN

   IF (REDUCED_GRID) THEN

      DO IJ = NIJS,NIJL
         IF (FIELD(IJ).EQ.ZMISS) CYCLE
         IP  = KLAT(IJ,2,1)
         IM  = KLAT(IJ,1,1)
         IP2 = KLAT(IJ,2,2)
         IM2 = KLAT(IJ,1,2)
         IF (IM.LE.NLAND .AND. IM2.LE.NLAND .AND.  &
&            IP.LE.NLAND .AND. IP2.LE.NLAND) CYCLE

         IF (IP.GT.NLAND .AND. IP2.GT.NLAND) THEN
            IF (FIELD(IP).NE.ZMISS .AND. FIELD(IP2).NE.ZMISS) THEN
               D_LAT_P = WLAT(IJ,2)*FIELD(IP)+(1.-WLAT(IJ,2))*FIELD(IP2)
            ELSE IF (FIELD(IP).NE.ZMISS) THEN
               D_LAT_P = FIELD(IP)
            ELSE IF (FIELD(IP2).NE.ZMISS) THEN
               D_LAT_P = FIELD(IP2)
            ELSE
               D_LAT_P = ZMISS
            END IF
         ELSE IF (IP.GT.NLAND) THEN
            D_LAT_P = FIELD(IP)
         ELSE IF (IP2.GT.NLAND) THEN
            D_LAT_P = FIELD(IP2)
         ELSE
           D_LAT_P = ZMISS
         END IF

         IF (IM.GT.NLAND .AND. IM2.GT.NLAND) THEN
            IF (FIELD(IM).NE.ZMISS .AND. FIELD(IM2).NE.ZMISS) THEN
               D_LAT_M = WLAT(IJ,1)*FIELD(IM)+(1.-WLAT(IJ,1))*FIELD(IM2)
            ELSE IF (FIELD(IM).NE.ZMISS) THEN
               D_LAT_M = FIELD(IM)
            ELSE IF (FIELD(IM2).NE.ZMISS) THEN
               D_LAT_M = FIELD(IM2)
            ELSE
               D_LAT_M = ZMISS
            END IF
         ELSE IF (IM.GT.NLAND) THEN
            D_LAT_M = FIELD(IM)
         ELSE IF (IP2.GT.NLAND) THEN
            D_LAT_M = FIELD(IM2)
         ELSE
           D_LAT_M = ZMISS
         END IF

         IF (D_LAT_P.NE.ZMISS .AND. D_LAT_P.NE.ZMISS) THEN
            D_LAT(IJ) = (D_LAT_P - D_LAT_M)/(2.*DELPHI)
         ELSE IF (D_LAT_M.NE.ZMISS) THEN
            D_LAT(IJ) = (FIELD(IJ)-D_LAT_M)/DELPHI
         ELSE IF (D_LAT_P.NE.ZMISS) THEN
            D_LAT(IJ) = (D_LAT_P-FIELD(IJ))/DELPHI
         ELSE
            D_LAT(IJ) = ZMISS
         END IF
      END DO

   ELSE

      DO IJ = NIJS,NIJL
         D_LAT(IJ) = ZMISS
         IF (FIELD(IJ).EQ.ZMISS) CYCLE
         IM = KLAT(IJ,1,1)
         IP = KLAT(IJ,2,1)
         IF (IM.LE.NLAND .AND. IP.LE.NLAND) CYCLE

         IF (IM.GT.NLAND .AND. IP.GT.NLAND) THEN
            IF (FIELD(IP).NE.ZMISS .AND. FIELD(IM).NE.ZMISS) THEN
               D_LAT(IJ) = (FIELD(IP)-FIELD(IM))/(2.*DELPHI)
            ELSE IF (FIELD(IM).NE.ZMISS .AND. FIELD(IP).EQ.ZMISS) THEN
               D_LAT(IJ) = (FIELD(IJ)-FIELD(IM))/DELPHI
            ELSE IF (FIELD(IM).EQ.ZMISS .AND. FIELD(IP).NE.ZMISS) THEN
               D_LAT(IJ) = (FIELD(IP)-FIELD(IJ))/DELPHI
            END IF

         ELSE IF (IM.GT.NLAND) THEN
            IF (FIELD(IM).NE.ZMISS) THEN
               D_LAT(IJ) = (FIELD(IJ)-FIELD(IM))/DELPHI
            END IF

         ELSE IF (IP.GT.NLAND) THEN
            IF (FIELD(IP).NE.ZMISS) THEN
               D_LAT(IJ) = (FIELD(IP)-FIELD(IJ))/DELPHI
            END IF
         END IF
      END DO

   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. EAST WEST GRADIENTS.                                                  !
!        --------------------                                                  !

IF (PRESENT(D_LON)) THEN

   DO IJ = NIJS,NIJL
      D_LON(IJ) = ZMISS
      IF (FIELD(IJ).EQ.ZMISS) CYCLE
      IM = KLON(IJ,1)
      IP = KLON(IJ,2)
      KX  = KXLT(IJ)
      IF (IM.LE.NLAND .AND. IP.LE.NLAND) CYCLE
 
      IF (IM.GT.NLAND .AND. IP.GT.NLAND) THEN
         IF (FIELD(IP).NE.ZMISS .AND. FIELD(IM).NE.ZMISS) THEN
            D_LON(IJ) = (FIELD(IP)-FIELD(IM))/(2.*DELLAM(KX))
         ELSE IF (FIELD(IM).NE.ZMISS) THEN
            D_LON(IJ) = (FIELD(IJ)-FIELD(IM))/DELLAM(KX)
         ELSE IF (FIELD(IP).NE.ZMISS) THEN
            D_LON(IJ) = (FIELD(IP)-FIELD(IJ))/DELLAM(KX)
         END IF

      ELSE IF (IM.GT.NLAND) THEN
         IF (FIELD(IM).NE.ZMISS) THEN
            D_LON(IJ) = (FIELD(IJ)-FIELD(IM))/DELLAM(KX)
         END IF

      ELSE IF (IP.GT.NLAND) THEN
         IF (FIELD(IP).NE.ZMISS) THEN
            D_LON(IJ) = (FIELD(IP)-FIELD(IJ))/DELLAM(KX)
         END IF
      END IF
      IF (SPHERICAL_RUN.AND.D_LON(IJ).NE.ZMISS) D_LON(IJ) = D_LON(IJ)*DCO(IJ)

   END DO
END IF

END SUBROUTINE GRADIENT_RAD

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_RADIATION_MODULE

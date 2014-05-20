MODULE WAM_PROPAGATION_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL VARIABLES, CONSTANTS AND PROCEDURES FOR THE       !
!   PROPAGATION.                                                               !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE,  ONLY: &
&       ABORT1,                &    !! TERMINATES PROCESSING.
&       flush1                      !! enforces output

USE WAM_TOPO_MODULE,  ONLY: &
&       PUT_DRY                     !! .

use wam_mpi_comp_module, only: &
&       mpi_exchng,            &
&       mpi_gather_cfl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE, ONLY: PI, CIRC, ZPI, R, DEG
USE WAM_FRE_DIR_MODULE, ONLY: ML, KL, FR, TH, GOM, DELTH, DELTR, COSTH, SINTH, &
&                             TCGOND, TFAK, TSIHKD, NDEPTH, DEPTHA, DEPTHD
USE WAM_GRID_MODULE,    ONLY: NSEA, IXLG, KXLT, KLAT, KLON,                    &
&                             DELPHI, DELLAM, SINPH, COSPH, ONE_POINT 
USE WAM_TIMOPT_MODULE,  ONLY: IDELPRO, IDELT,                                  &
&                             SPHERICAL_RUN, SHALLOW_RUN,                      &
&                             REFRACTION_D_RUN, REFRACTION_C_RUN
USE WAM_FILE_MODULE,    ONLY: IU06, ITEST
USE WAM_MODEL_MODULE,   ONLY: DEPTH, INDEP, U, V

use wam_mpi_module,     only: petotal, irank, nijs, nijl, ninf, nsup

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
include 'mpif.h'
PRIVATE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NUMBER OF ADVECTION ITERATIONS PER CALL OF WAMODEL.                   !
!        ---------------------------------------------------                   !

INTEGER, PUBLIC, SAVE :: NADV       = 0
REAL,    SAVE :: NEWIDELPRO =-1. !! TIMESTEP WAM PROPAGATION IN SECONDS.
INTEGER, SAVE :: COUNTER    =-1  !! NO. OF PROPAGATION SUB-TIMESTEP IN SECONDS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LATITUDE LONGITUDE PROPAGATION.                                       !
!        -------------------------------                                       !

REAL, ALLOCATABLE :: CGOND(:,:)  !! SHALLOW WATER GROUP VELOCITIES.
REAL, PUBLIC, ALLOCATABLE :: DCO(:)    !! 1./ COS(LATITUDE).
REAL, PUBLIC, ALLOCATABLE :: DPSN(:,:) !! COS(LATITUDE SOUTH)/COS(LATITUDE) 
                                       !! second index = 1.
                                       !! COS(LATITUDE NORTH)/COS(LATITUDE) 
                                       !! second index = 2.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. DEPTH AND CURRENT PART OF THETA DOT.                                  !
!        ------------------------------------                                  !

REAL, ALLOCATABLE :: THDD(:,:,:)  !! THETA DOT.
REAL, ALLOCATABLE :: SIDC(:,:,:)  !! SIGMA DOT.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. FREQUENCY, DIRECTION NEIGHTBOURS.                                     !
!        ------------------------------------                                  !

INTEGER, PUBLIC, ALLOCATABLE :: MP(:)    !! FREQUENCY INDEX +1.
INTEGER, PUBLIC, ALLOCATABLE :: MM(:)    !! FREQUENCY INDEX -1.
INTEGER, PUBLIC, ALLOCATABLE :: KP(:)    !! DIRECTION INDEX +1.
INTEGER, PUBLIC, ALLOCATABLE :: KM(:)    !! DIRECTION INDEX -1.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. CFL INFORMATION.                                                      !
!        ----------------                                                      !

REAL, PARAMETER   :: XLIMIT = 0.999  !! MAXIMUM ALLOWED CFL NUMBER
REAL, ALLOCATABLE :: CFLMAX(:)       !! MAXIMUM CFL NUMBER on each process.
INTEGER :: MAXPOINT(3)               !! POINT, DIRECTION AND FREQUENCY INDEX
                                     !! OF MAXIMUM CFL NUMBER on each process.
REAL    :: total_cfl = 0.            !! MAXIMUM CFL NUMBER frm all processes.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE GRADIENT                    !! CALCULATES GRADIENTS.
   MODULE PROCEDURE GRADIENT
END INTERFACE
PUBLIC GRADIENT

INTERFACE PREPARE_PROPAGATION         !! PREPARES PROPAGATION.
   MODULE PROCEDURE PREPARE_PROPAGATION
END INTERFACE
PUBLIC PREPARE_PROPAGATION

INTERFACE PRINT_PROPAGATION_STATUS    !! PRINTS PROPAGATION STATUS.
   MODULE PROCEDURE PRINT_PROPAGATION_STATUS
END INTERFACE
PUBLIC PRINT_PROPAGATION_STATUS

INTERFACE PROPAGS                     !! PROPAGATION.
   MODULE PROCEDURE PROPAGS
END INTERFACE
PUBLIC PROPAGS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE CHECK_CFL                   !! CFL CHECK.
   MODULE PROCEDURE CHECK_CFL
END INTERFACE
PRIVATE CHECK_CFL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE GRADIENT (FIELD, D_LAT, D_LON)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   GRADIENT - CALCULATES GRADIENT.                                            !
!                                                                              !
!     K.P. HUBBERT              AUGUST   1988                                  !
!     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO CALCULATE THE GRADIENT VECTOR FIELD FROM INPUT FIELD.               !
!                                                                              !
!     METHOD.                                                                  !
!     ------                                                                   !
!                                                                              !
!       CENTRAL DIFFERENCING.                                                  !
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

REAL,    INTENT(IN)  :: FIELD(1:NSEA)  !! INPUT FIELD.
REAL,    INTENT(OUT) :: D_LAT(nijs:nijl)  !! LATITUDE  DERIVATIVE.
REAL,    INTENT(OUT) :: D_LON(nijs:nijl)  !! LONGITUDE DERIVATIVE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NORTH SOUTH GRADIENTS.                                                !
!        ----------------------                                                !

D_LAT(nijs:nijl) = 0.0
WHERE (KLAT(nijs:nijl,1).GT.NINF-1 .AND. KLAT(nijs:nijl,2).GT.NINF-1)          &
&    D_LAT(nijs:nijl) = (FIELD(KLAT(nijs:nijl,2))-FIELD(KLAT(nijs:nijl,1)))    &
&                     / (2.*DELPHI)
CALL PUT_DRY (D_LAT, nijs, nijl, 0.)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. EAST WEST GRADIENTS.                                                  !
!        --------------------                                                  !

D_LON(nijs:nijl) = 0.0
WHERE (KLON(nijs:nijl,2).GT.NINF-1 .AND. KLON(nijs:nijl,1).GT.NINF-1)          &
&     D_LON(nijs:nijl) = (FIELD(KLON(nijs:nijl,2))-FIELD( KLON(nijs:nijl,1)))  &
&                      / (2.*DELLAM(KXLT(nijs:nijl)))
CALL PUT_DRY ( D_LON, nijs, nijl, 0.)

END SUBROUTINE GRADIENT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_PROPAGATION

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PREPARE_PROPAGATION - PREPARES DOT TERMS FROM DEPTH AND CURRENT GRADIENT.  !
!                                                                              !
!     H. GUNTHER   GKSS/ECMWF   17/02/91                                       !
!     H. GUNTHER   GKSS         DECEMBER 2001    FT90                          !
!     E. MYKLEBUST              NOVEMBER 2005    MPI PARALLELIZATION
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF REFRACTION DOT TERMS FOR PROPAGATION.                   !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE DEPTH AND CURRENT GRADIENTS ARE COMPUTED.                          !
!       DEPENDING OF THE MODEL OPTIONS THE DEPTH AND CURRENT REFRACTION PART   !
!       FOR THETA DOT AND THE COMPLETE SIGMA DOT TERM ARE COMPUTED.            !
!       FOR A MULTIBLOCK VERSION ALL TERMS ARE WRITTEN TO MASS STORAGE (IU15). !
!       FOR A ONE BLOCK MODEL THE SIGMA DOT TERM IS WRITTEN ONLY THE OTHER     !
!       TERMS ARE STORED IN WAM_MODEL_MODULE.                                  !
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

INTEGER :: M, K
REAL    :: SD, CD, SS, SC, CC

INTEGER, ALLOCATABLE, DIMENSION(:) :: INDEP_G
REAL,    ALLOCATABLE, DIMENSION(:) :: TEMP, THDC
REAL,    ALLOCATABLE, DIMENSION(:) :: DDPHI, DDLAM
REAL,    ALLOCATABLE, DIMENSION(:) :: DUPHI, DULAM, DVPHI, DVLAM

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. SHALLOW WATER TABLE INDICES (IF SHALLOW WATER).                       !
!        ----------------------------                                          !

IF (SHALLOW_RUN) THEN
   IF (ALLOCATED(INDEP)) DEALLOCATE (INDEP)
   ALLOCATE (INDEP(ninf:nsup))
   INDEP(ninf:nsup) = NINT(LOG(DEPTH(ninf:nsup)/DEPTHA)/LOG(DEPTHD)+1.)
   INDEP(ninf:nsup) = MAX(INDEP(ninf:nsup),1)
   INDEP(ninf:nsup) = MIN(INDEP(ninf:nsup),NDEPTH)
ELSE
   !Indep must be allocated to be able to select subarray in call to IMPLSCH
   IF (ALLOCATED(INDEP)) DEALLOCATE (INDEP)
   ALLOCATE (INDEP(ninf:nsup))
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. GROUP VELOCITIES (IF SHALLOW WATER).                                  !
!         -----------------------------------                                  !

IF (SHALLOW_RUN) THEN
   IF (ALLOCATED(CGOND)) DEALLOCATE(CGOND)      
   ALLOCATE(CGOND(ninf-1:nsup,1:ML))
   CGOND(ninf-1,:) = TCGOND(NDEPTH,:)
   CGOND(ninf:nsup,:) = TCGOND(INDEP(ninf:nsup),:)
   DO M = 1,ML
      CALL PUT_DRY (CGOND(ninf:nsup,M), ninf,nsup, TCGOND(NDEPTH,M))
   END DO
END IF

IF (ONE_POINT) RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CHECK REFRACTION AND IF REFRECTION IS OFF DO CFL CHECK AND RETURN.    !
!        ------------------------------------------------------------------    !

IF (REFRACTION_C_RUN .AND. (.NOT.ALLOCATED(U) .OR. .NOT.ALLOCATED(V))) THEN
   WRITE(IU06,*) '*****************************************************'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*      FATAL ERROR IN SUB.PREPARE_PROPAGATION       *'
   WRITE(IU06,*) '*      ======================================       *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '* A RUN WITH CURRENT REFRACTION IS REQUESTED,       *'
   WRITE(IU06,*) '* BUT CURRENT FIELDS ARE NOT IN WAM_MODEL_MODULE.   *'
   WRITE(IU06,*) '* IS A CURRENT FILE NAME DEFINED IN THE USER INPUT? *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*         PROGRAM ABORTS   PROGRAM ABORTS           *'
   WRITE(IU06,*) '*                                                   *'
   WRITE(IU06,*) '*****************************************************'
   CALL ABORT1
END IF

IF (.NOT.REFRACTION_D_RUN .AND. .NOT.REFRACTION_C_RUN) THEN
   CALL CHECK_CFL
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTE DEPTH AND CURRENT GRADIENTS.                                  !
!        ------------------------------------                                  !

IF (SHALLOW_RUN) THEN
   ALLOCATE (DDPHI(nijs:nijl), DDLAM(nijs:nijl))
   CALL GRADIENT (DEPTH, DDPHI, DDLAM)
   IF (SPHERICAL_RUN) DDLAM = DDLAM*DCO(nijs:nijl)
END IF
IF (REFRACTION_C_RUN) THEN
   ALLOCATE (DUPHI(nijs:nijl), DULAM(nijs:nijl))
   ALLOCATE (DVPHI(nijs:nijl), DVLAM(nijs:nijl))
   CALL GRADIENT (U, DUPHI, DULAM)
   IF (SPHERICAL_RUN) DULAM = DULAM*DCO(nijs:nijl)
   CALL GRADIENT (V, DVPHI, DVLAM)
   IF (SPHERICAL_RUN) DVLAM = DVLAM*DCO(nijs:nijl)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. COMPUTE DEPTH INDEX TO SMOOTH GRADIENTS.                              !
!        ----------------------------------------                              !

IF (SHALLOW_RUN) THEN
   ALLOCATE (TEMP(nijs:nijl))
   ALLOCATE (INDEP_G(NIJS:NIJL))

   TEMP(NIJS:NIJL) = 4.*DEPTH(NIJS:NIJL)
   INDEP_G(NIJS:NIJL) = 4
   WHERE (KLAT(NIJS:NIJL,1).GT.NINF-1 .AND. KLAT(NIJS:NIJL,2).GT.NINF-1)
      TEMP(NIJS:NIJL) = TEMP(NIJS:NIJL)+DEPTH(KLAT(NIJS:NIJL,2))               &
&                                      +DEPTH(KLAT(NIJS:NIJL,1))
      INDEP_G(NIJS:NIJL) = INDEP_G(NIJS:NIJL) + 2
   END WHERE

   WHERE (KLON(NIJS:NIJL,1).GT.NINF-1 .AND. KLON(NIJS:NIJL,2).GT.NINF-1)
      TEMP(NIJS:NIJL) = TEMP(NIJS:NIJL)+DEPTH(KLON(NIJS:NIJL,2))               &
&                                      +DEPTH(KLON(NIJS:NIJL,1))
      INDEP_G(NIJS:NIJL) = INDEP_G(NIJS:NIJL) + 2
   END WHERE

   TEMP(NIJS:NIJL) = TEMP(NIJS:NIJL)/REAL(INDEP_G(NIJS:NIJL))
   INDEP_G(NIJS:NIJL) = NINT(LOG(TEMP(NIJS:NIJL)/DEPTHA)/LOG(DEPTHD)+1.)
   INDEP_G(NIJS:NIJL) = MAX(INDEP_G(NIJS:NIJL),1)
   INDEP_G(NIJS:NIJL) = MIN(INDEP_G(NIJS:NIJL),NDEPTH)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. DEPTH OF THETA DOT.                                                   !
!        -------------------                                                   !

IF (SHALLOW_RUN) THEN
   IF (.NOT.ALLOCATED(THDD) ) ALLOCATE (THDD(nijs:nijl,1:KL,1:ML))
   DO  K = 1,KL
      TEMP(nijs:nijl) = (SINTH(K)+SINTH(KP(K)))*DDPHI(nijs:nijl)               &
&                     - (COSTH(K)+COSTH(KP(K)))*DDLAM(nijs:nijl)
       DO M = 1,ML
           THDD(nijs:nijl,K,M) = TEMP(nijs:nijl) * TSIHKD(INDEP_G(nijs:nijl),M)
       END DO
   END DO

END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. DEPTH PART OF SIGMA DOT.                                              !
!        ------------------------                                              !

IF (SHALLOW_RUN .AND. REFRACTION_C_RUN) THEN
   TEMP(nijs:nijl) = V(nijs:nijl)*DDPHI(nijs:nijl)+U(nijs:nijl)*DDLAM(nijs:nijl)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     8. CURRENT PART OF SIGMA DOT AND THETA DOT.                              !
!        -----------------------------------------                             !


IF (REFRACTION_C_RUN) THEN

!     8.1 CURRENT PART OF THETA DOT.                                           !
!         --------------------------                                           !

   IF (SHALLOW_RUN) THEN
      IF (.NOT.ALLOCATED(THDC) ) ALLOCATE (THDC(nijs:nijl))
      IF (.NOT.ALLOCATED(SIDC) ) ALLOCATE (SIDC(nijs:nijl,1:KL,1:ML))
   ELSE
      IF (.NOT.ALLOCATED(SIDC) ) ALLOCATE (SIDC(nijs:nijl,1:KL,ML:ML))
      IF (.NOT.ALLOCATED(THDD) ) ALLOCATE (THDD(nijs:nijl,1:KL,ML:ML))
   END IF

   DO  K = 1,KL
      SS  = SINTH(K)**2 + SINTH(KP(K))**2
      SC  = SINTH(K)*COSTH(K) + SINTH(KP(K))*COSTH(KP(K))
      CC  = COSTH(K)**2 +COSTH(KP(K))**2

      IF (SHALLOW_RUN) THEN
         THDC(nijs:nijl) = SS*DUPHI(nijs:nijl) + SC*DVPHI(nijs:nijl)          &
&                        - SC*DULAM(nijs:nijl) - CC*DVLAM(nijs:nijl)
         DO M =1,ML
            THDD(nijs:nijl,K,M) = THDD(nijs:nijl,K,M) + THDC(nijs:nijl)
         END DO
      ELSE
          THDD(nijs:nijl,K,ML) = SS*DUPHI(nijs:nijl) + SC*DVPHI(nijs:nijl)     &
&                       - SC*DULAM(nijs:nijl) - CC*DVLAM(nijs:nijl)
      END IF
   END DO

!     8.2 CURRENT PART OF SIGMA DOT.                                           !
!         --------------------------                                           !

   DO  K = 1,KL
      SD = SINTH(K)
      CD = COSTH(K)
      SS  = SD**2
      SC  = SD*CD
      CC  = CD**2
      SIDC(nijs:nijl,K,ML) = - SC*DUPHI(nijs:nijl) - CC*DVPHI(nijs:nijl)       &
&                            - SS*DULAM(nijs:nijl) - SC*DVLAM(nijs:nijl)
   END DO

!    8.3 ADD SHALLOW WATER TO CURRENT PART OF SIGMA DOT.                       !
!        -----------------------------------------------                       !

   IF (SHALLOW_RUN) THEN
      DO M = 1,ML
         DO  K = 1,KL
            SIDC(nijs:nijl,K,M) = (SIDC(nijs:nijl,K,ML)                        &
&                   * TCGOND(INDEP(nijs:nijl),M)                               &
&                   + TEMP(nijs:nijl) * TSIHKD(INDEP_G(nijs:nijl),M))          &
&                   * TFAK(INDEP_G(nijs:nijl),M)
         END DO
      END DO
   END IF
END IF

IF (ALLOCATED(INDEP_G) ) DEALLOCATE (INDEP_G)
IF (ALLOCATED(THDC) ) DEALLOCATE (THDC)

IF (ALLOCATED(TEMP) ) DEALLOCATE(TEMP) 
IF (ALLOCATED(DDPHI)) DEALLOCATE(DDPHI) 
IF (ALLOCATED(DDLAM)) DEALLOCATE(DDLAM) 
IF (ALLOCATED(DUPHI)) DEALLOCATE(DUPHI) 
IF (ALLOCATED(DULAM)) DEALLOCATE(DULAM) 
IF (ALLOCATED(DVPHI)) DEALLOCATE(DVPHI) 
IF (ALLOCATED(DVLAM)) DEALLOCATE(DVLAM) 

! ---------------------------------------------------------------------------- !
!                                                                              !
!     9. CHECK CFL.                                                            !
!        ----------                                                            !

CALL CHECK_CFL

END SUBROUTINE PREPARE_PROPAGATION

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PRINT_PROPAGATION_STATUS

WRITE(IU06,*) '  '
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '              PROPAGATION MODULE STATUS:'
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '  '
IF (NADV.LE.0) THEN
   WRITE(IU06,*) '              PROPAGATION MODULE NOT INITIALIZED'
ELSE
   IF (ONE_POINT) THEN
      WRITE(IU06,*) ' '
      WRITE(IU06,*) '  ONE POINT MODEL: PROPAGATION IS NOT INITIALIZED'
      WRITE(IU06,*) '  NUMBER OF SOURCE TIME STEPS IN ONE CALL ',              &
&                                     'OF SUB WAVEMDL IS ' , 'NADV = ', NADV
   ELSE
      WRITE(IU06,*) '  NUMBER OF PROPAGATION STEPS IN ONE CALL ',              &
&                                     'OF SUB WAVEMDL IS ' , 'NADV = ', NADV
      WRITE(IU06,*) ' '
      WRITE(IU06,*) '  MAX. CFL NUMBER IS total_cfl = ', total_cfl
      WRITE(IU06,*) '  '
      IF (COUNTER.GT.1) THEN
          WRITE(IU06,*) ' PROPAGATION IS DIVIDED INTO ',COUNTER,' STEPS'    
          WRITE(IU06,*) ' INTERNAL PROPAGATION STEP IS ',NEWIDELPRO,' SECONDS'
          WRITE(IU06,*) '  '
      END IF
   END IF
END IF

END SUBROUTINE PRINT_PROPAGATION_STATUS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PROPAGS (F3)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PROPAGS - COMPUTATION OF A PROPAGATION TIME STEP.                          !
!                                                                              !
!     S.D. HASSELMANN.                                                         !
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER                               !
!                                                                              !
!     MODIFIED BY   H. GUNTHER   01/06/90    -   LAND POINTS ARE TAKEN         !
!                             OUT OF BLOCKS AND REFRACTION INTEGRATION         !
!                             CORRECTED FOR N-S AND S-N PROPAGATION.           !
!                                                                              !
!     K.P. HUBBERT                /07/89    -   DEPTH AND CURRENT              !
!     S. HASSELMANN   MPIFM       /04/90        REFRACTION SHALLOW             !
!                                                                              !
!     H. GUNTHER   GKSS/ECMWF   17/01/91    -   MODIFIED FOR CYCLE_4           !
!                                                                              !
!     H. GUENTHER   GKSS   FEBRUARY 2002    -   FT 90                          !
!                                                                              !
!     A. Behrens  MSC/GKSS  January 2004    -   MPI parallelization            !
!                                                                              !
!     E. Myklebust         February 2005    -   Optimization                   !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF A PROPAGATION TIME STEP.                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       FIRST ORDER FLUX SCHEME.                                               !
!                                                                              !
!     INTERNAL SUBROUTINES FOR PROPAGATION.                                    !
!     -------------------------------------                                    !
!                                                                              !
!         SPHERICAL GRID.                                                      !
!                                                                              !
!       *P_SPHER_DEEP*         DEEP WATER WITHOUT CURRENT REFRACTION.          !
!       *P_SPHER_SHALLOW*      SHALLOW WATER WITHOUT CURRENT REFRACTION.       !
!       *P_SPHER_DEEP_CURR*    DEEP WATER WITH CURRENT REFRACTION.             !
!       *P_SPHER_SHALLOW_CURR* SHALLOW WATER WITH DEPTH AND CURRENT REFRACTION.!
!                                                                              !
!         CARTESIAN GRID.                                                      !
!                                                                              !
!       *P_CART_DEEP*          DEEP WATER WITHOUT CURRENT REFRACTION.          !
!       *P_CART_SHALLOW*       SHALLOW WATER WITHOUT CURRENT REFRACTION.       !
!       *P_CART_DEEP_CUR*      DEEP WATER WITH CURRENT REFRACTION.             !
!       *P_CART_SHALLOW_CUR*   SHALLOW WATER WITH DEPTH AND CURRENT REFRACTION.!
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

IMPLICIT NONE

REAL,    INTENT(INOUT)   :: F3(NIJS:NIJL,KL,ML)     !! SPECTRUM AT TIME T+DELT.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: I, IJ, IS, IT, IC, ID, K, M, ierr

INTEGER :: KMK,KPK,MMM,MPM !! KM(K), ... NEEDED FOR COMPILER TO VECTORIZE LOOPS

REAL   :: DELPH0, DELTH0, DELTHR, DELFR0
REAL   :: SD, CD, SDA, CDA
REAL   :: DFP, DFM, CGS, CGC, SM, SP

REAL   :: DTC     !! WEIGHT OF CENTER.
REAL   :: DPN     !! WEIGHT OF NORTH POINT.
REAL   :: DPS     !! WEIGHT OF SOUTH POINT.
REAL   :: DLE     !! WEIGHT OF EAST  POINT.
REAL   :: DLW     !! WEIGHT OF WEST POINT.
REAL   :: DTP     !! WEIGHT OF DIRECTION +1.
REAL   :: DTM     !! WEIGHT OF DIRECTION -1.
REAL   :: DOP     !! WEIGHT OF FREQUENCY +1.
REAL   :: DOM     !! WEIGHT OF FREQUENCY -1.

REAL   :: DELLA0(NIJS:NIJL)  !! DELTA TIME / LONGITUTE INCREMENTS. (S/M).

REAL, ALLOCATABLE, DIMENSION(:) :: DRGM !! DIR. PART OF THETA DOT FROM GRID.
REAL, ALLOCATABLE, DIMENSION(:) :: WOK1 !! WORK ARRAY.
REAL, ALLOCATABLE, DIMENSION(:) :: WOK2 !! WORK ARRAY.

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: f1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. INITIAL.                                                              !
!        --------                                                              !

allocate (f1(ninf-1:nsup,1:kl,1:ml))
   
IF (SPHERICAL_RUN)  ALLOCATE (DRGM(NIJS:NIJL))
IF (REFRACTION_C_RUN) THEN
   ALLOCATE (WOK1(ninf-1:nsup))
   ALLOCATE (WOK2(ninf-1:nsup))
END IF

DO I = 1, COUNTER

   F1(nijs:nijl,:,:) = F3(:,:,:)           !! COPY INPUT SPECTRA
   call mpi_exchng (f1(ninf:nsup,:,:))
   call mpi_barrier (mpi_comm_world, ierr)
   F1(ninf-1,:,:) = 0.                     !! SPECTRUM AT LAND TO ZERO.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.0 SELECT CASES.                                                        !
!         -------------                                                        !

   IF (SPHERICAL_RUN) THEN

!     1.1 PROPAGATION ON SPHERICAL GRID.                                       !
!         ------------------------------                                       !

      IF (.NOT.REFRACTION_C_RUN) THEN
         IF (SHALLOW_RUN) THEN
            CALL P_SPHER_SHALLOW  !! SHALLOW WATER WITHOUT CURRENT REFRACTION.
         ELSE
            CALL P_SPHER_DEEP     !! DEEP WATER WITHOUT CURRENT REFRACTION.
         END IF
      ELSE
         IF (SHALLOW_RUN) THEN
            CALL P_SPHER_SHALLOW_CURR !! SHALLOW WATER WITH DEPTH AND CURRENT REF.
         ELSE
            CALL P_SPHER_DEEP_CURR    !! DEEP WATER WITH CURRENT REFRACTION.
         END IF
      END IF

   ELSE

!     1.2 PROPAGATION ON CARTESIAN GRID.                                       !
!         ------------------------------                                       !

      IF (.NOT.REFRACTION_C_RUN) THEN
         IF (SHALLOW_RUN) THEN
            CALL P_CART_SHALLOW      !! SHALLOW WATER WITHOUT CURRENT REFRACTION.
         ELSE
            CALL P_CART_DEEP         !! DEEP WATER WITHOUT CURRENT REFRACTION.
         END IF
      ELSE
         IF (SHALLOW_RUN) THEN
            CALL P_CART_SHALLOW_CUR  !! SHALLOW WATER WITH DEPTH AND CURRENT REF.
         ELSE
            CALL P_CART_DEEP_CUR     !! DEEP WATER WITH CURRENT REFRACTION.
         END IF
      END IF
   END IF
END DO

DEALLOCATE (f1)
IF (SPHERICAL_RUN)    DEALLOCATE (DRGM)
IF (REFRACTION_C_RUN) DEALLOCATE (WOK1, WOK2)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.0 RETURN.                                                              !
!         -------                                                              !

RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.0 INTERNAL SUBROUTINES.                                                !
!         ---------------------                                                !

CONTAINS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.1 PROPAGATION FOR CARTESIAN GRID WITHOUT CURRENT REFRACTION (DEEP).    !
!         -----------------------------------------------------------------    !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE P_CART_DEEP

   DELLA0(nijs:nijl) = NEWIDELPRO/DELLAM(KXLT(nijs:nijl))
   DELPH0            = NEWIDELPRO/DELPHI

   FRE: DO M = 1,ML                   !! LOOP OVER FREQUENCIES.
      DIR: DO K = 1,KL
         SD = SINTH(K)
         CD = COSTH(K)

         IF (SD.LT.0) THEN                  !! INDEX FOR ADJOINING LONGITUDE.
            IS = 2
         ELSE
            IS = 1
         END IF
         IF (CD.LT.0) THEN                  !! INDEX FOR ADJOINING LATITUDE.
            IC = 2
         ELSE
            IC = 1
         END IF

         DPN = ABS(CD*DELPH0*GOM(M))
         SD  = ABS(SD*GOM(M))
         DO IJ = nijs,nijl
            DLE = SD*DELLA0(IJ)
            DTC = 1. - DLE - DPN

            F3(IJ,K,M) = DTC*F1(IJ,K,M)                                      &
&                      + DPN*F1(KLAT(IJ,IC),K,M)                             &
&                      + DLE*F1(KLON(IJ,IS),K,M)
         END DO
      END DO DIR
   END DO FRE

END SUBROUTINE P_CART_DEEP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.2 PROPAGATION FOR CARTESIAN GRID WITHOUT CURRENT REFRACTION (SHALLOW). !
!         -------------------------------------------------------------------- !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE P_CART_SHALLOW

   DELLA0(nijs:nijl) = 0.5*NEWIDELPRO/DELLAM(KXLT(nijs:nijl))
   DELPH0 = 0.5*NEWIDELPRO/DELPHI
   DELTHR = 0.5*NEWIDELPRO/DELTH

   FRE: DO M = 1,ML
   DIR: DO K = 1,KL
      SD = SINTH(K)
      CD = COSTH(K)

      IF (SD.LT.0) THEN                  !! INDEX FOR ADJOINING LONGITUDE.
         IS = 2
         IT = 1
      ELSE
         IS = 1
         IT = 2
      END IF
      IF (CD.LT.0) THEN                  !! INDEX FOR ADJOINING LATITUDE.
         IC = 2
         ID = 1
      ELSE
         IC = 1
         ID = 2
      END IF

      SD = ABS(SD)
      CD = ABS(CD*DELPH0)

      IF (REFRACTION_D_RUN) THEN
         DO IJ = nijs,nijl
            DLE = SD*(CGOND(KLON(IJ,IS),M) + CGOND(IJ,M))  * DELLA0(IJ)
            DPN = CD*(CGOND(KLAT(IJ,IC),M) + CGOND(IJ,M))

            DTC = 1. - SD*(CGOND(KLON(IJ,IT),M) + CGOND(IJ,M)) * DELLA0(IJ)    &
&                    - CD*(CGOND(KLAT(IJ,ID),M) + CGOND(IJ,M))
            KPK = KP(K)
            KMK = KM(K)
            DTP = THDD(IJ,K  ,M)*DELTHR
            DTM = THDD(IJ,KMK,M)*DELTHR

            DTC = DTC - MAX(0.,DTP) + MIN (0.,DTM)
            DTP = - MIN (0.,DTP)
            DTM =   MAX (0.,DTM)

            F3(IJ,K,M) = DTC*F1(IJ,K,M)                                        &
&                      + DPN*F1(KLAT(IJ,IC),K  ,M) + DLE*F1(KLON(IJ,IS),K  ,M) &
&                      + DTP*F1(IJ         ,KPK,M) + DTM*F1(IJ         ,KMK,M)

         END DO
      ELSE
         DO IJ = nijs,nijl
            DLE = SD*(CGOND(KLON(IJ,IS),M) + CGOND(IJ,M))  * DELLA0(IJ)
            DPN = CD*(CGOND(KLAT(IJ,IC),M) + CGOND(IJ,M))

            DTC = 1. - SD*(CGOND(KLON(IJ,IT),M) + CGOND(IJ,M)) * DELLA0(IJ)    &
&                    - CD*(CGOND(KLAT(IJ,ID),M) + CGOND(IJ,M))

            F3(IJ,K,M) = DTC*F1(IJ,K,M)                                        &
&                      + DPN*F1(KLAT(IJ,IC),K  ,M) + DLE*F1(KLON(IJ,IS),K  ,M)
         END DO
      END IF

   END DO DIR
   END DO FRE

END SUBROUTINE P_CART_SHALLOW

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.3 PROPAGATION FOR CARTESIAN GRID WITH CURRENT REFRACTION (DEEP).       !
!         --------------------------------------------------------------       !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE P_CART_DEEP_CUR

   DELLA0(nijs:nijl) = 0.5*NEWIDELPRO/DELLAM(KXLT(nijs:nijl))
   DELPH0 = 0.5*NEWIDELPRO/DELPHI
   DELTHR = 0.5*NEWIDELPRO/DELTH
   DELFR0 = 0.5*NEWIDELPRO*2.1/0.2

   FRE: DO M = 1,ML
      MMM = MM(M)
      MPM = MP(M)
      DIR: DO K = 1,KL
         KMK = KM(K)
         KPK = KP(K)

         SD = SINTH(K)
         CD = COSTH(K)

         CGS = GOM(M)*SD                !! GROUP VELOCITIES.
         CGC = GOM(M)*CD

         WOK1(ninf-1)    = CGS
         WOK1(ninf:nsup) = U(ninf:nsup) + CGS 
         WOK2(ninf-1)    = CGC*DELPH0
         WOK2(ninf:nsup) = (V(ninf:nsup) + CGC)*DELPH0

         DO IJ=NIJS,NIJL
            DLW = (WOK1(IJ) + WOK1(KLON(IJ,1))) * DELLA0(IJ)
            DLE = (WOK1(IJ) + WOK1(KLON(IJ,2))) * DELLA0(IJ)
            DTC = 1. - MAX(0.,DLE) + MIN(0.,DLW)
            DLE = -MIN(0.,DLE)
            DLW =  MAX(0.,DLW)

            DPS = WOK2(IJ) + WOK2(KLAT(IJ,1))
            DPN = WOK2(IJ) + WOK2(KLAT(IJ,2))
            DTC = DTC - MAX(0.,DPN) + MIN(0.,DPS)
            DPN = -MIN(0.,DPN)
            DPS =  MAX(0.,DPS)

            DTP = THDD(IJ,K  ,ML)*DELTHR
            DTM = THDD(IJ,KMK,ML)*DELTHR
            DTC = DTC - MAX(0.,DTP) + MIN(0.,DTM)
            DTP = -MIN(0.,DTP)
            DTM =  MAX(0.,DTM)

            DOM = SIDC(IJ,K,ML) * DELFR0
            DTC = DTC - ABS(DOM)
            DOP = -MIN(0.,DOM) / 1.1
            DOM =  MAX(0.,DOM) * 1.1

            F3(IJ,K,M) = DTC * F1(IJ,K,M)                &
&                      + DPN * F1(KLAT(IJ,2),K,M)        &
&                      + DPS * F1(KLAT(IJ,1),K,M)        &
&                      + DLE * F1(KLON(IJ,2),K,M)        &
&                      + DLW * F1(KLON(IJ,1),K,M)        &
&                      + DTP * F1(IJ,KPK,M  )            &
&                      + DTM * F1(IJ,KMK,M  )            &
&                      + DOP * F1(IJ,K  ,MPM)            &
&                      + DOM * F1(IJ,K  ,MMM)
         END DO
      END DO DIR
   END DO FRE

END SUBROUTINE P_CART_DEEP_CUR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.4 PROPAGATION FOR CARTESIAN GRID WITH CURRENT REFRACTION (SHALLOW).    !
!         -----------------------------------------------------------------    !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE P_CART_SHALLOW_CUR

   DELLA0(nijs:nijl) = 0.5*NEWIDELPRO/DELLAM(KXLT(nijs:nijl))
   DELPH0 = 0.5*NEWIDELPRO/DELPHI
   DELTHR = 0.5*NEWIDELPRO/DELTH
   DELFR0 = 0.5*NEWIDELPRO/(0.1*ZPI)

   FRE: DO M = 1,ML
      MMM = MM(M)
      MPM = MP(M)
      DFP = DELFR0/FR(M)
      DFM = DELFR0/FR(MMM)
      DIR: DO K = 1,KL
         SD = SINTH(K)
         CD = COSTH(K)

         KMK = KM(K)
         KPK = KP(K)

         WOK1(ninf-1)    = SD*CGOND(ninf-1,M)
         WOK2(ninf-1)    = CD*CGOND(ninf-1,M)*DELPH0
         DO IJ = ninf,nsup
            WOK1(IJ) =  U(IJ) + SD*CGOND(IJ,M)
            WOK2(IJ) = (V(IJ) + CD*CGOND(IJ,M)) * DELPH0
         END DO

         DO IJ = nijs,nijl
            DLW = (WOK1(IJ) + WOK1(KLON(IJ,1))) * DELLA0(IJ)
            DLE = (WOK1(IJ) + WOK1(KLON(IJ,2))) * DELLA0(IJ)
            DTC = 1. - MAX(0.,DLE) + MIN(0.,DLW)
            DLE = -MIN(0.,DLE)
            DLW =  MAX(0.,DLW)

            DPS = WOK2(IJ) + WOK2(KLAT(IJ,1))
            DPN = WOK2(IJ) + WOK2(KLAT(IJ,2))
            DTC = DTC - MAX(0.,DPN) + MIN(0.,DPS)
            DPN = -MIN(0.,DPN)
            DPS =  MAX(0.,DPS)

            DTP = THDD(IJ,K  ,M)*DELTHR
            DTM = THDD(IJ,KMK,M)*DELTHR
            DTC = DTC - MAX(0.,DTP) + MIN(0.,DTM)
            DTP = -MIN(0.,DTP)
            DTM =  MAX(0.,DTM)

            DOP = (SIDC(IJ,K,M) + SIDC(IJ,K,MPM))*DFP
            DOM = (SIDC(IJ,K,M) + SIDC(IJ,K,MMM))*DFM
            DTC = DTC - MAX(0.,DOP) + MIN(0.,DOM)
            DOP = -MIN(0.,DOP)/1.1
            DOM =  MAX(0.,DOM)*1.1

            F3(IJ,K,M) = DTC * F1(IJ,K,M )               &
&                      + DPN * F1(KLAT(IJ,2),K,M)        &
&                      + DPS * F1(KLAT(IJ,1),K,M)        &
&                      + DLE * F1(KLON(IJ,2),K,M)        &
&                      + DLW * F1(KLON(IJ,1),K,M)        &
&                      + DTP * F1(IJ,KPK,M  )            &
&                      + DTM * F1(IJ,KMK,M  )            &
&                      + DOP * F1(IJ,K  ,MPM)            &
&                      + DOM * F1(IJ,K  ,MMM)
         END DO
      END DO DIR
   END DO FRE

END SUBROUTINE P_CART_SHALLOW_CUR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.5 PROPAGATION FOR SPHERICAL GRID WITHOUT CURRENT REFRACTION (DEEP).    !
!         -----------------------------------------------------------------    !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE P_SPHER_DEEP

   DELTH0 = 0.5*NEWIDELPRO/DELTR
   DELPH0 = 0.5*NEWIDELPRO/DELPHI

   DO IJ = nijs,nijl
      DELLA0(IJ) = NEWIDELPRO*DCO(IJ)/DELLAM(KXLT(IJ))
      DRGM(IJ)   = SINPH(KXLT(IJ))*DCO(IJ)
   END DO

   FRE: DO M = 1,ML
      DIR: DO K = 1,KL
         SD = SINTH(K)
         CD = COSTH(K)
         SDA = ABS(SD)
         CDA = ABS(CD*DELPH0)

         KMK = KM(K)
         KPK = KP(K)

         IF (SD.LT.0) THEN
            IS = 2
         ELSE
            IS = 1
         END IF
         IF (CD.LT.0) THEN
            IC = 2
            ID = 1
         ELSE
            IC = 1
            ID = 2
         END IF

         SP  = DELTH0*(SINTH(K)+SINTH(KPK))         !! GRID REFRACTION WEIGHTS.
         SM  = DELTH0*(SINTH(K)+SINTH(KMK))

         DO IJ = nijs,nijl
            DTC = CDA*(DPSN(IJ,ID) + 1.)
            DPN = CDA*(DPSN(IJ,IC) + 1.)  !! LATITUDE WEIGHTS.

            DLE = SDA*DELLA0(IJ)  !! LONGITUDE WEIGHTS.

            DTP = DRGM(IJ)*SP
            DTM = DRGM(IJ)*SM
            DTC = DTC + DLE + MAX(0. , DTP) - MIN(0. , DTM)
            DTP = -MIN(0. , DTP)
            DTM =  MAX(0. , DTM)

            F3(IJ,K,M) = (1. - GOM(M)*DTC) * F1(IJ,K,M)      &
&                      + GOM(M)*(  DPN * F1(KLAT(IJ,IC),K,M) &
&                                + DLE * F1(KLON(IJ,IS),K,M) &
&                                + DTP * F1(IJ,KPK,M)        &
&                                + DTM * F1(IJ,KMK,M) )
         END DO
      END DO DIR
   END DO FRE

END SUBROUTINE P_SPHER_DEEP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.6 PROPAGATION FOR SPHERICAL GRID WITHOUT CURRENT REFRACTION (SHALLOW). !
!         -------------------------------------------------------------------- !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE P_SPHER_SHALLOW

   DELTH0 = 0.5*NEWIDELPRO/DELTR
   DELPH0 = 0.5*NEWIDELPRO/DELPHI
   DELTHR = 0.5*NEWIDELPRO/DELTH

   DO IJ = nijs,nijl
      DELLA0(IJ) = 0.5*NEWIDELPRO*DCO(IJ)/DELLAM(KXLT(IJ))
      DRGM(IJ)   = SINPH(KXLT(IJ))*DCO(IJ)
   END DO

   FRE: DO M = 1,ML
      DIR: DO K = 1,KL
         SD = SINTH(K)
         CD = COSTH(K)
         SDA = ABS(SD)
         CDA = ABS(CD*DELPH0)

         KMK=KM(K)
         KPK=KP(K)
         SP  = DELTH0*(SINTH(K)+SINTH(KPK))   !! PRE_COMPUTE GRID REFRACTION.
         SM  = DELTH0*(SINTH(K)+SINTH(KMK))

         IF (SD.LT.0) THEN                    !! INDEX FOR ADJOINING POINTS.
            IS = 2
            IT = 1
         ELSE
            IS = 1
            IT = 2
         END IF
         IF (CD.LT.0) THEN
            IC = 2
            ID = 1
         ELSE
            IC = 1
            ID = 2
         END IF

         IF (REFRACTION_D_RUN) THEN                 !! DEPTH REFRACTION
            DO IJ = nijs,nijl
               DTC = 1. - CDA*(CGOND(KLAT(IJ,ID),M)*DPSN(IJ,ID) + CGOND(IJ,M))
               DPN =      CDA*(CGOND(KLAT(IJ,IC),M)*DPSN(IJ,IC) + CGOND(IJ,M))

               DTC = DTC - SDA*(CGOND(KLON(IJ,IT),M) + CGOND(IJ,M))*DELLA0(IJ)
               DLE =       SDA*(CGOND(KLON(IJ,IS),M) + CGOND(IJ,M))*DELLA0(IJ)

               DTP = DRGM(IJ)*CGOND(IJ,M)           !! REFRACTION WEIGHTS
               DTM = DTP*SM + THDD(IJ,KMK,M)*DELTHR 
               DTP = DTP*SP + THDD(IJ,K  ,M)*DELTHR 
               DTC = DTC - MAX(0. , DTP) + MIN(0. , DTM)
               DTP = -MIN(0. , DTP)
               DTM =  MAX(0. , DTM)

               F3(IJ,K,M) = DTC * F1(IJ,K,M)                &
&                         + DPN * F1(KLAT(IJ,IC),K,M)       &
&                         + DLE * F1(KLON(IJ,IS),K,M)       &
&                         + DTP * F1(IJ,KPK,M)              &
&                         + DTM * F1(IJ,KMK,M)

            END DO

         ELSE

            DO IJ = nijs,nijl
               DTC = 1. - CDA*(CGOND(KLAT(IJ,ID),M)*DPSN(IJ,ID) + CGOND(IJ,M))
               DPN =      CDA*(CGOND(KLAT(IJ,IC),M)*DPSN(IJ,IC) + CGOND(IJ,M))

               DTC = DTC - SDA*(CGOND(KLON(IJ,IT),M) + CGOND(IJ,M))*DELLA0(IJ)
               DLE =       SDA*(CGOND(KLON(IJ,IS),M) + CGOND(IJ,M))*DELLA0(IJ)

               DTP = DRGM(IJ)*CGOND(IJ,M) !! REFRACTION WEIGHTS
               DTM = DTP*SM
               DTP = DTP*SP

               DTC = DTC - MAX(0. , DTP) + MIN(0. , DTM)
               DTP = -MIN(0. , DTP)
               DTM =  MAX(0. , DTM)

               F3(IJ,K,M) = DTC * F1(IJ,K,M)                &
&                         + DPN * F1(KLAT(IJ,IC),K,M)       &
&                         + DLE * F1(KLON(IJ,IS),K,M)       &
&                         + DTP * F1(IJ,KPK,M)              &
&                         + DTM * F1(IJ,KMK,M)

            END DO

         END IF
      END DO DIR
   END DO FRE

END SUBROUTINE P_SPHER_SHALLOW

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.7 PROPAGATION FOR SPHERICAL GRID WITH CURRENT REFRACTION (DEEP).       !
!         --------------------------------------------------------------       !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE P_SPHER_DEEP_CURR

   DELPH0 = 0.5*NEWIDELPRO/DELPHI
   DELTH0 = 0.5*NEWIDELPRO/DELTR
   DELTHR = 0.5*NEWIDELPRO/DELTH
   DELFR0 = 0.5*NEWIDELPRO*2.1/0.2

   DO IJ = nijs,nijl
      DELLA0(IJ) = 0.5*NEWIDELPRO*DCO(IJ)/DELLAM(KXLT(IJ))
      DRGM(IJ)   = SINPH(KXLT(IJ))*DCO(IJ)
   END DO

   FRE: DO M = 1,ML
      MPM = MP(M)
      MMM = MM(M)
      DIR: DO K = 1,KL
         SD = SINTH(K)
         CD = COSTH(K)
         KPK = KP(K)
         KMK = KM(K)
         SP = DELTH0*(SINTH(K)+SINTH(KPK))*GOM(M)     !! PRE-COMPUTE GRID REFRACTION.
         SM = DELTH0*(SINTH(K)+SINTH(KMK))*GOM(M)

         CGS = GOM(M)*SD                       !! GROUP VELOCITIES.
         CGC = GOM(M)*CD

         WOK1(ninf-1)    = CGS
         WOK1(ninf:nsup) = U(ninf:nsup) + CGS
         WOK2(ninf-1)    = CGC*DELPH0
         WOK2(ninf:nsup) = (V(ninf:nsup) + CGC) * DELPH0

         DO IJ = nijs,nijl
            DLW = (WOK1(IJ) + WOK1(KLON(IJ,1))) * DELLA0(IJ)
            DLE = (WOK1(IJ) + WOK1(KLON(IJ,2))) * DELLA0(IJ)
            DTC = 1. - MAX(0.,DLE) + MIN(0.,DLW)
            DLE = -MIN(0.,DLE)
            DLW =  MAX(0.,DLW)

            DPS = WOK2(IJ)+WOK2(KLAT(IJ,1))*DPSN(IJ,1)
            DPN = WOK2(IJ)+WOK2(KLAT(IJ,2))*DPSN(IJ,2)
            DTC = DTC - MAX(0.,DPN) + MIN(0.,DPS)
            DPN = -MIN(0.,DPN)
            DPS =  MAX(0.,DPS)

            DTP = SP*DRGM(IJ) + THDD(IJ,K  ,ML)*DELTHR
            DTM = SM*DRGM(IJ) + THDD(IJ,KMK,ML)*DELTHR
            DTC = DTC - MAX(0.,DTP) + MIN(0.,DTM)
            DTP = -MIN(0.,DTP)
            DTM =  MAX(0.,DTM)

            DOM =  SIDC(IJ,K,ML) * DELFR0
            DTC =  DTC - ABS(DOM)
            DOP = -MIN(0.,DOM)/1.1
            DOM =  MAX(0.,DOM)*1.1

            F3(IJ,K,M) = DTC * F1(IJ,K,M)                &
&                      + DPN * F1(KLAT(IJ,2),K,M)        &
&                      + DPS * F1(KLAT(IJ,1),K,M)        &
&                      + DLE * F1(KLON(IJ,2),K,M)        &
&                      + DLW * F1(KLON(IJ,1),K,M)        &
&                      + DTP * F1(IJ,KPK,M  )            &
&                      + DTM * F1(IJ,KMK,M  )            &
&                      + DOP * F1(IJ,K  ,MPM)            &
&                      + DOM * F1(IJ,K  ,MMM)

         END DO
      END DO DIR
   END DO FRE

END SUBROUTINE P_SPHER_DEEP_CURR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.8 PROPAGATION FOR SPHERICAL GRID WITH CURRENT REFRACTION (SHALLOW).    !
!         -----------------------------------------------------------------    !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE P_SPHER_SHALLOW_CURR

   DELPH0 = 0.5*NEWIDELPRO/DELPHI
   DELTH0 = 0.5*NEWIDELPRO/DELTR
   DELTHR = 0.5*NEWIDELPRO/DELTH
   DELFR0 = 0.5*NEWIDELPRO/(0.1*ZPI)

   DO IJ = nijs,nijl
      DELLA0(IJ) = 0.5*NEWIDELPRO*DCO(IJ)/DELLAM(KXLT(IJ))
      DRGM(IJ)   = SINPH(KXLT(IJ))*DCO(IJ)
   END DO

   FRE: DO M = 1,ML
      MPM = MP(M)
      MMM = MM(M)
      DFP = DELFR0/FR(M)
      DFM = DELFR0/FR(MMM)
      DIR: DO K = 1,KL
         SD = SINTH(K)
         CD = COSTH(K)
         KPK = KP(K)
         KMK = KM(K)
         SP = DELTH0*(SINTH(K)+SINTH(KPK))      !! GRID REFRACTION.
         SM = DELTH0*(SINTH(K)+SINTH(KMK))

         WOK1(ninf-1)    = SD*CGOND(ninf-1,M)
         WOK1(ninf:nsup) = U(ninf:nsup)+SD*CGOND(ninf:nsup,M)
         WOK2(ninf-1)    = CD*CGOND(ninf-1,M)*DELPH0
         WOK2(ninf:nsup) = (V(ninf:nsup)+CD*CGOND(ninf:nsup,M))*DELPH0

         DO IJ = nijs,nijl

            DLW = ( WOK1(IJ) + WOK1(KLON(IJ,1)) ) * DELLA0(IJ)
            DLE = ( WOK1(IJ) + WOK1(KLON(IJ,2)) ) * DELLA0(IJ)
            DTC = 1. - MAX(0.,DLE) + MIN(0.,DLW)
            DLE = -MIN(0.,DLE)
            DLW =  MAX(0.,DLW)

            DPS = WOK2(IJ)+WOK2(KLAT(IJ,1))*DPSN(IJ,1)
            DPN = WOK2(IJ)+WOK2(KLAT(IJ,2))*DPSN(IJ,2)
            DTC = DTC - MAX(0.,DPN) + MIN(0.,DPS)
            DPN = -MIN(0.,DPN)
            DPS =  MAX(0.,DPS)

            DTP = SP*DRGM(IJ)*CGOND(IJ,M) + THDD(IJ,K  ,M)*DELTHR
            DTM = SM*DRGM(IJ)*CGOND(IJ,M) + THDD(IJ,KMK,M)*DELTHR
            DTC = DTC - MAX(0.,DTP) + MIN(0.,DTM)
            DTP = -MIN(0.,DTP)
            DTM =  MAX(0.,DTM)

            DOP = (SIDC(IJ,K,M) + SIDC(IJ,K,MPM))*DFP
            DOM = (SIDC(IJ,K,M) + SIDC(IJ,K,MMM))*DFM
            DTC = DTC - MAX(0.,DOP) + MIN(0.,DOM)
            DOP = -MIN(0.,DOP)/1.1
            DOM =  MAX(0.,DOM)*1.1

            F3(IJ,K,M) = DTC * F1(IJ,K,M )               &
&                      + DPN * F1(KLAT(IJ,2),K,M)        &
&                      + DPS * F1(KLAT(IJ,1),K,M)        &
&                      + DLE * F1(KLON(IJ,2),K,M)        &
&                      + DLW * F1(KLON(IJ,1),K,M)        &
&                      + DTP * F1(IJ,KPK,M  )            &
&                      + DTM * F1(IJ,KMK,M  )            &
&                      + DOP * F1(IJ,K  ,MPM)            &
&                      + DOM * F1(IJ,K  ,MMM)
         END DO
      END DO DIR
   END DO FRE

END SUBROUTINE P_SPHER_SHALLOW_CURR

! ---------------------------------------------------------------------------- !

END SUBROUTINE PROPAGS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CHECK_CFL

! ---------------------------------------------------------------------------- !
!                                                                              !
!     CHECK_CFL  - CFL CHECK.                                                  !
!                                                                              !
!     H. GUNTHER   GKSS   FEBRUARY 2002                                        !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO CHECK THE PROPAGATION TIME STEP.                                    !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       FOR EACH GRID POINT, FREQUENCY AND DIRECTION THE RELATIVE LOSS OF      !
!       ENERGY PER TIME STEP IS COMPUTED AND CHECKED TO BE LESS THAN 1.        !
!                                                                              !
!     INTERNAL SUBROUTINES.                                                    !
!     ---------------------                                                    !
!                                                                              !
!         SPHERICAL GRID.                                                      !
!                                                                              !
!       *C_SPHER_DEEP*         DEEP WATER WITHOUT CURRENT REFRACTION.          !
!       *C_SPHER_SHALLOW*      SHALLOW WATER WITHOUT CURRENT REFRACTION.       !
!       *C_SPHER_DEEP_CURR*    DEEP WATER WITH CURRENT REFRACTION.             !
!       *C_SPHER_SHALLOW_CURR* SHALLOW WATER WITH DEPTH AND CURRENT REFRACTION.!
!                                                                              !
!         CARTESIAN GRID.                                                      !
!                                                                              !
!       *C_CART_DEEP*          DEEP WATER WITHOUT CURRENT REFRACTION.          !
!       *C_CART_SHALLOW*       SHALLOW WATER WITHOUT CURRENT REFRACTION.       !
!       *C_CART_DEEP_CUR*      DEEP WATER WITH CURRENT REFRACTION.             !
!       *C_CART_SHALLOW_CUR*   SHALLOW WATER WITH DEPTH AND CURRENT REFRACTION.!
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


! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, K, M, IS, IC
INTEGER :: KMK,KPK,MMM,MPM !! KM(K), ... NEEDED FOR COMPILER TO VECTORIZE LOOPS

REAL    :: DELPH0, DELTH0, DELTHR 
REAL    :: DELFR0, DFP, DFM, CGS, CGC, SM, SP
REAL    :: SD, CD, SDA, CDA
REAL    :: DTC, DTP, DTM

REAL, ALLOCATABLE, DIMENSION(:) :: DPH, DLA
REAL, ALLOCATABLE, DIMENSION(:) :: DELLA0 
REAL, ALLOCATABLE, DIMENSION(:) :: DRGM

REAL, ALLOCATABLE, DIMENSION(:,:,:)    :: CF  !! CFL NUMBERS.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. INITIAL.                                                              !
!        --------                                                              !

ALLOCATE (DELLA0(nijs:nijl)) 
IF (SPHERICAL_RUN) ALLOCATE (DRGM  (nijs:nijl))
ALLOCATE (CF    (nijs:nijl,1:KL,1:ML)) 

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.0 SELECT CASE.                                                         !
!         -------------                                                        !

IF (SPHERICAL_RUN) THEN

!     1.1 CFL CHECK FOR SPHERICAL GRID.                                        !
!         -----------------------------                                        !

   IF (.NOT.REFRACTION_C_RUN) THEN
      IF (SHALLOW_RUN) THEN
         CALL C_SPHER_SHALLOW  !! SHALLOW WATER WITHOUT CURRENT REFRACTION.
      ELSE
         ALLOCATE(DLA   (nijs:nijl))
         CALL C_SPHER_DEEP     !! DEEP WATER WITHOUT CURRENT REFRACTION.
         DEALLOCATE (DLA)
      END IF
   ELSE
      ALLOCATE(DPH   (ninf-1:nsup))
      ALLOCATE(DLA   (ninf-1:nsup))
      IF (SHALLOW_RUN) THEN
         CALL C_SPHER_SHALLOW_CURR !! SHALLOW WATER WITH DEPTH AND CURRENT REF.
      ELSE
         CALL C_SPHER_DEEP_CURR    !! DEEP WATER WITH CURRENT REFRACTION.
      END IF
      DEALLOCATE (DPH, DLA)
   END IF

ELSE

!     1.2 CFL CHECK FOR CARTESIAN GRID.                                        !
!         -----------------------------                                        !

   IF (.NOT.REFRACTION_C_RUN) THEN
      IF (SHALLOW_RUN) THEN
         CALL C_CART_SHALLOW      !! SHALLOW WATER WITHOUT CURRENT REFRACTION.
      ELSE
         CALL C_CART_DEEP         !! DEEP WATER WITHOUT CURRENT REFRACTION.
      END IF

   ELSE
      ALLOCATE(DPH   (ninf-1:nsup))
      ALLOCATE(DLA   (ninf-1:nsup))
      IF (SHALLOW_RUN) THEN
         CALL C_CART_SHALLOW_CUR  !! SHALLOW WATER WITH DEPTH AND CURRENT REF.
      ELSE
         CALL C_CART_DEEP_CUR     !! DEEP WATER WITH CURRENT REFRACTION.
      END IF
      DEALLOCATE (DPH, DLA)
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.0 EVALUATE CFL.                                                        !
!         -------------                                                        !

ALLOCATE (CFLMAX(1:petotal))

CFLMAX(irank) = MAXVAL(CF(nijs:nijl,1:KL,1:ML))
MAXPOINT = MAXLOC(CF(nijs:nijl,1:KL,1:ML))
if (petotal/=1) then
   call mpi_gather_cfl (CFLMAX)
   total_cfl = maxval(CFLMAX)
else
   total_cfl = CFLMAX(irank)
end if

IF (total_cfl .GT. XLIMIT) THEN
   WRITE(IU06,*)
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) ' +                                                   +'
   WRITE(IU06,*) ' +           WARNING ERROR IN SUB. CFLCHECK          +'
   WRITE(IU06,*) ' +           ==============================          +'
   WRITE(IU06,*) ' +                                                   +'
   WRITE(IU06,*) ' +            VIOLATIONS OF CFL-CRITERION            +'
   WRITE(IU06,*) ' +                                                   +'
   WRITE(IU06,*) ' + Maxium clf at all processes = ', total_cfl
   WRITE(IU06,*) ' +                                                   +'
   WRITE(IU06,*) ' + Process Number = ', irank
   WRITE(IU06,*) ' +    MAXIMUM CFL = ', CFLMAX(irank)
   WRITE(IU06,*) ' +       AT POINT = ', MAXPOINT(1)
   WRITE(IU06,*) ' +           KLON = ', IXLG(MAXPOINT(1))
   WRITE(IU06,*) ' +           KLAT = ', KXLT(MAXPOINT(1))
   WRITE(IU06,*) ' +      DIRECTION = ', TH(MAXPOINT(2))*DEG
   WRITE(IU06,*) ' +      FREQUENCY = ', FR(MAXPOINT(3))
   WRITE(IU06,*) ' +                                                   +'
   WRITE(IU06,*) ' +                 PROGRAM CONTINUES                 +'
   WRITE(IU06,*) ' +      WITH MODIFIED PROPAGATION TIME STEP          +'

   COUNTER = CEILING (total_cfl)
   NEWIDELPRO = REAL(IDELPRO)/REAL(COUNTER)
   total_cfl = total_cfl/COUNTER

   WRITE(IU06,*) ' +  PROPAGATION TIME STEP IS DIVIDED BY ',COUNTER 
   WRITE(IU06,*) ' +  NEW PROPAGATION TIME STEP IS        ',NEWIDELPRO
   WRITE(IU06,*) ' +                                                   +'
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++++++'

ELSE
   NEWIDELPRO = REAL(IDELPRO)
   COUNTER = 1
END IF

DEALLOCATE (DELLA0) 
IF (SPHERICAL_RUN) DEALLOCATE (DRGM)
DEALLOCATE (CF    ) 
DEALLOCATE (CFLMAX) 

RETURN

! ---------------------------------------------------------------------------- !

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.0 INTERNAL SUBROUTINES.                                                !
!         ---------------------                                                !

CONTAINS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.1 CFL CHECK FOR CARTESIAN GRID WITHOUT CURRENT REFRACTION (DEEP).      !
!         ---------------------------------------------------------------      !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE C_CART_DEEP

   DELLA0(nijs:nijl) = REAL(IDELPRO)/DELLAM(KXLT(nijs:nijl))
   DELPH0 = REAL(IDELPRO)/DELPHI

   FRE: DO M = 1,ML
      DIR: DO K = 1,KL
         CF(nijs:nijl,K,M) = ( ABS(SINTH(K)*DELLA0(nijs:nijl))                 &
&                            + ABS(COSTH(K)*DELPH0))*GOM(M)
      END DO DIR
   END DO FRE

END SUBROUTINE C_CART_DEEP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.2 CFL CHECK FOR CARTESIAN GRID WITHOUT CURRENT REFRACTION (SHALLOW).   !
!         ------------------------------------------------------------------   !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE C_CART_SHALLOW

   DELLA0(nijs:nijl) = 0.5*REAL(IDELPRO)/DELLAM(KXLT(nijs:nijl))
   DELPH0 = 0.5*REAL(IDELPRO)/DELPHI
   DELTHR = 0.5*REAL(IDELPRO)/DELTH

   DIR: DO K = 1,KL
      KPK = KP(K)
      KMK = KM(K)
      SD = SINTH(K)
      CD = COSTH(K)
      IF (SD.LT.0) THEN                  !! INDEX FOR ADJOINING LONGITUDE.
         IS = 1
      ELSE
         IS = 2
      END IF
      IF (CD.LT.0) THEN                  !! INDEX FOR ADJOINING LATITUDE.
         IC = 1
      ELSE
         IC = 2
      END IF

      SD = ABS(SD)
      CD = ABS(CD*DELPH0)

      FRE: DO M = 1,ML
         CF(nijs:nijl,K,M) = SD*( CGOND(KLON(nijs:nijl,IS),M)                 &
&                               + CGOND(nijs:nijl,M))*DELLA0(nijs:nijl)       &
&                          + CD*( CGOND(KLAT(nijs:nijl,IC),M)                 &
&                               + CGOND(nijs:nijl,M))
         IF (REFRACTION_D_RUN) THEN
            CF(nijs:nijl,K,M) = CF(nijs:nijl,K,M)                             &
&                             + DELTHR  * ( MAX(0.,THDD(nijs:nijl,K  ,M))     &
&                                          -MIN(0.,THDD(nijs:nijl,KMK,M)))
         END IF
      END DO FRE
   END DO DIR

END SUBROUTINE C_CART_SHALLOW

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.3 CFL CHECK FOR CARTESIAN GRID WITH CURRENT REFRACTION (DEEP).         !
!         ------------------------------------------------------------         !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE C_CART_DEEP_CUR

   DELLA0(nijs:nijl) = 0.5*REAL(IDELPRO)/DELLAM(KXLT(nijs:nijl))
   DELPH0 = 0.5*REAL(IDELPRO)/DELPHI
   DELTHR = 0.5*REAL(IDELPRO)/DELTH
   DELFR0 = 0.5*REAL(IDELPRO)*2.1/0.2

   DIR: DO K = 1,KL
      KPK = KP(K)
      KMK = KM(K)
      SD = SINTH(K)
      CD = COSTH(K)

      FRE: DO M = 1,ML

         CGS = GOM(M)*SD               !! GROUP VELOCITIES.
         CGC = GOM(M)*CD

         DLA(ninf-1) = CGS
         DPH(ninf-1) = CGC*DELPH0
         DLA(ninf:nsup) = (U(ninf:nsup) + CGS)
         DPH(ninf:nsup) = (V(ninf:nsup) + CGC)*DELPH0

         DO IJ = nijs,nijl
            CF(IJ,K,M) = (  MAX(0. , DLA(IJ) + DLA(KLON(IJ,2)))               &
&                         - MIN(0. , DLA(IJ) + DLA(KLON(IJ,1)))) * DELLA0(IJ) &
&                      + MAX(0. , DPH(IJ) + DPH(KLAT(IJ,2)))                  &
&                      - MIN(0. , DPH(IJ) + DPH(KLAT(IJ,1)))                  &
&                      + (  MAX(0. , THDD(IJ,K  ,ML))                         &
&                         - MIN(0. , THDD(IJ,KMK,ML)) )*DELTHR                &
&                      + ABS(SIDC(IJ,K,ML)*DELFR0)

         END DO
      END DO FRE
   END DO DIR

END SUBROUTINE C_CART_DEEP_CUR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.4 CFL CHECK FOR CARTESIAN GRID WITH CURRENT REFRACTION (SHALLOW).      !
!         ---------------------------------------------------------------      !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE C_CART_SHALLOW_CUR

   DELLA0(nijs:nijl) = 0.5*REAL(IDELPRO)/DELLAM(KXLT(nijs:nijl))
   DELPH0 = 0.5*REAL(IDELPRO)/DELPHI
   DELTHR = 0.5*REAL(IDELPRO)/DELTH
   DELFR0 = 0.5*REAL(IDELPRO)/(0.1*ZPI)

   DIR: DO K = 1,KL
      KPK = KP(K)
      KMK = KM(K)
      SD = SINTH(K)
      CD = COSTH(K)

      FRE: DO M = 1,ML
         MMM = MM(M)
         MPM = MP(M)
         DFP = DELFR0/FR(M)
         DFM = DELFR0/FR(MMM)

         DLA(ninf-1) = SD*CGOND(ninf-1,M)
         DPH(ninf-1) = CD*CGOND(ninf-1,M)*DELPH0
         DLA(ninf:nsup) = (U(ninf:nsup) + SD*CGOND(ninf:nsup,M))
         DPH(ninf:nsup) = (V(ninf:nsup) + CD*CGOND(ninf:nsup,M)) * DELPH0

         DO IJ = nijs,nijl
            CF(IJ,K,M) = (   MAX(0. , DLA(IJ) + DLA(KLON(IJ,2)))              &
&                          - MIN(0. , DLA(IJ) + DLA(KLON(IJ,1))))*DELLA0(IJ)  &
&                      + MAX(0. , DPH(IJ) + DPH(KLAT(IJ,2)))                  &
&                      - MIN(0. , DPH(IJ) + DPH(KLAT(IJ,1)))                  &
&                      + (   MAX(0. , THDD(IJ,K  ,M))                         &
&                          - MIN(0. , THDD(IJ,KMK,M)))*DELTHR                 &
&                      + MAX(0. , (SIDC(IJ,K,M) + SIDC(IJ,K,MPM))*DFP)        &
&                      - MIN(0. , (SIDC(IJ,K,M) + SIDC(IJ,K,MMM))*DFM)
         END DO
      END DO FRE
   END DO DIR

END SUBROUTINE C_CART_SHALLOW_CUR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.5 CFL CHECK FOR SPHERICAL GRID WITHOUT CURRENT REFRACTION (DEEP).      !
!         ---------------------------------------------------------------      !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE C_SPHER_DEEP

   DELTH0 = 0.5*REAL(IDELPRO)/DELTR
   DELPH0 = 0.5*REAL(IDELPRO)/DELPHI

   DO IJ = nijs,nijl 
      DELLA0(IJ) = REAL(IDELPRO)*DCO(IJ)/DELLAM(KXLT(IJ))
      DRGM(IJ)   = SINPH(KXLT(IJ))*DCO(IJ)
   END DO

   DIR: DO K = 1,KL
      KPK = KP(K)
      KMK = KM(K)
      SP  = DELTH0*(SINTH(K)+SINTH(KPK))
      SM  = DELTH0*(SINTH(K)+SINTH(KMK))
      SDA = ABS(SINTH(K))
      CDA = ABS(COSTH(K))*DELPH0
      IF (COSTH(K).GT.0.) THEN
         IC = 2         
      ELSE
         IC = 1
      END IF

      DLA(nijs:nijl) = MAX(0. , DRGM(nijs:nijl)*SP)                            &
&                    - MIN(0. , DRGM(nijs:nijl)*SM)                            &
&                    + CDA*(DPSN(nijs:nijl,IC) + 1.)                           &
&                    + SDA*DELLA0(nijs:nijl)

      FRE: DO M = 1,ML
         CF(nijs:nijl,K,M) = DLA(nijs:nijl)*GOM(M)
      END DO FRE
   END DO DIR

END SUBROUTINE C_SPHER_DEEP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.6 CFL CHECK FOR SPHERICAL GRID WITHOUT CURRENT REFRACTION (SHALLOW).   !
!         ------------------------------------------------------------------   !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE C_SPHER_SHALLOW

   DELTH0 = 0.5*REAL(IDELPRO)/DELTR
   DELPH0 = 0.5*REAL(IDELPRO)/DELPHI
   DELTHR = 0.5*REAL(IDELPRO)/DELTH

   DO IJ = nijs,nijl 
      DELLA0(IJ) = 0.5*REAL(IDELPRO)*DCO(IJ)/DELLAM(KXLT(IJ))
      DRGM(IJ)   = SINPH(KXLT(IJ))*DCO(IJ)
   END DO

   DIR: DO K = 1,KL
      KPK = KP(K)
      KMK = KM(K)
      SD = SINTH(K)
      CD = COSTH(K)
      SDA = ABS(SD)
      CDA = ABS(CD)*DELPH0

!     COMPUTE GRID REFRACTION.                                                 !

      SP  = DELTH0*(SINTH(K)+SINTH(KPK))
      SM  = DELTH0*(SINTH(K)+SINTH(KMK))

      IF (SD.LT.0) THEN          !! INDEX FOR ADJOINING POINTS.
         IS = 1
      ELSE
         IS = 2
      END IF
      IF (CD.LT.0.) THEN
         IC = 1
      ELSE
         IC = 2
      END IF

      IF (REFRACTION_D_RUN) THEN  !! DEPTH REFRACTION
         FRE1: DO M = 1,ML
            DO IJ = nijs,nijl 
               DTC = CDA * (CGOND(KLAT(IJ,IC),M)*DPSN(IJ,IC) + CGOND(IJ,M))    &
                   + SDA * (CGOND(KLON(IJ,IS),M) + CGOND(IJ,M))*DELLA0(IJ)
               DTP = DRGM(IJ)*SP*CGOND(IJ,M) + THDD(IJ,K  ,M)*DELTHR
               DTM = DRGM(IJ)*SM*CGOND(IJ,M) + THDD(IJ,KMK,M)*DELTHR

               CF(IJ,K,M) = DTC + MAX (0.,DTP) - MIN (0.,DTM)
            END DO
         END DO FRE1

      ELSE                          !! WITHOUT DEPTH REFRACTION

         FRE2: DO M = 1,ML
            DO IJ = nijs,nijl 
               DTC = CDA * (CGOND(KLAT(IJ,IC),M)*DPSN(IJ,IC) + CGOND(IJ,M))    &
                   + SDA * (CGOND(KLON(IJ,IS),M) + CGOND(IJ,M))*DELLA0(IJ)

               DTP = DRGM(IJ)*SP*CGOND(IJ,M)   !! REFRACTION WEIGHTS.
               DTM = DRGM(IJ)*SM*CGOND(IJ,M)

               CF(IJ,K,M) = DTC + MAX(0.,DTP) - MIN(0.,DTM)
            END DO
         END DO FRE2
      END IF

   END DO DIR

END SUBROUTINE C_SPHER_SHALLOW

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.7 CFL CHECK FOR SPHERICAL GRID WITH CURRENT REFRACTION (DEEP).         !
!         ------------------------------------------------------------         !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE C_SPHER_DEEP_CURR

   DELPH0 = 0.5*REAL(IDELPRO)/DELPHI
   DELTH0 = 0.5*REAL(IDELPRO)/DELTR
   DELTHR = 0.5*REAL(IDELPRO)/DELTH
   DELFR0 = 0.5*REAL(IDELPRO)*2.1/0.2

   DO IJ = nijs,nijl
      DELLA0(IJ) = 0.5*REAL(IDELPRO)*DCO(IJ)/DELLAM(KXLT(IJ))
      DRGM(IJ) = SINPH(KXLT(IJ))*DCO(IJ)
   END DO

   DIR: DO K = 1,KL
      KPK = KP(K)
      KMK = KM(K)
      SD = SINTH(K)
      CD = COSTH(K)

      SP = DELTH0*(SINTH(K)+SINTH(KPK))      !! GRID REFRACTION.
      SM = DELTH0*(SINTH(K)+SINTH(KMK))

      FRE: DO M = 1,ML

         CGS = GOM(M)*SD                       !! GROUP VELOCITIES.
         CGC = GOM(M)*CD

         DLA(ninf-1) = CGS
         DPH(ninf-1) = CGC*DELPH0
         DLA(ninf:nsup) = (U(ninf:nsup) + CGS)
         DPH(ninf:nsup) = (V(ninf:nsup) + CGC) * DELPH0

         DO IJ = nijs,nijl
            DTP = DRGM(IJ)*SP*GOM(M) + DELTHR*THDD(IJ,K  ,ML)   !! REFRACTION WEIGHTS.
            DTM = DRGM(IJ)*SM*GOM(M) + DELTHR*THDD(IJ,KMK,ML)
            CF(IJ,K,M) = MAX(0. , DLA(IJ) + DLA(KLON(IJ,2)))*DELLA0(IJ)      &
&                      - MIN(0. , DLA(IJ) + DLA(KLON(IJ,1)))*DELLA0(IJ)      &
&                      + MAX(0. , DPH(IJ) + DPH(KLAT(IJ,2))*DPSN(IJ,2))      &
&                      - MIN(0. , DPH(IJ) + DPH(KLAT(IJ,1))*DPSN(IJ,1))      &
&                      + MAX(0. , DTP)                                       &
&                      - MIN(0. , DTM)                                       &
&                      + ABS(SIDC(IJ,K,ML) * DELFR0)

         END DO
      END DO FRE
   END DO DIR

END SUBROUTINE C_SPHER_DEEP_CURR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.8 CFL CHECK FOR SPHERICAL GRID WITH CURRENT REFRACTION (SHALLOW).      !
!         ---------------------------------------------------------------      !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE C_SPHER_SHALLOW_CURR

   DELPH0 = 0.5*REAL(IDELPRO)/DELPHI
   DELTH0 = 0.5*REAL(IDELPRO)/DELTR
   DELTHR = 0.5*REAL(IDELPRO)/DELTH
   DELFR0 = 0.5*REAL(IDELPRO)/(0.1*ZPI)

   DO IJ = nijs,nijl
      DELLA0(IJ) = 0.5*REAL(IDELPRO)*DCO(IJ)/DELLAM(KXLT(IJ))
      DRGM(IJ) = SINPH(KXLT(IJ))*DCO(IJ)
   END DO

   DIR: DO K = 1,KL
      KPK = KP(K)
      KMK = KM(K)
      SD = SINTH(K)
      CD = COSTH(K)

      SP = DELTH0*(SINTH(K)+SINTH(KPK))    !! GRID REFRACTION.
      SM = DELTH0*(SINTH(K)+SINTH(KMK))

      FRE: DO M = 1,ML
         MMM = MM(M)
         MPM = MP(M)
         DFP = DELFR0/FR(M)
         DFM = DELFR0/FR(MMM)

         DLA(ninf-1) = SD*CGOND(ninf-1,M)
         DPH(ninf-1) = CD*CGOND(ninf-1,M)*DELPH0
         DLA(ninf:nsup) = (U(ninf:nsup)+SD*CGOND(ninf:nsup,M))
         DPH(ninf:nsup) = (V(ninf:nsup)+CD*CGOND(ninf:nsup,M)) * DELPH0

         DO IJ = nijs,nijl
            DTP = DRGM(IJ)*SP*CGOND(IJ,M) + THDD(IJ,K  ,M)*DELTHR 
            DTM = DRGM(IJ)*SM*CGOND(IJ,M) + THDD(IJ,KMK,M)*DELTHR

            CF(IJ,K,M) = MAX(0. , DLA(IJ) + DLA(KLON(IJ,2)))*DELLA0(IJ)        &
&                      - MIN(0. , DLA(IJ) + DLA(KLON(IJ,1)))*DELLA0(IJ)        &
&                      + MAX(0. , DPH(IJ) + DPH(KLAT(IJ,2))*DPSN(IJ,2))        &
&                      - MIN(0. , DPH(IJ) + DPH(KLAT(IJ,1))*DPSN(IJ,1))        &
&                      + MAX(0. , DTP)                                         &
&                      - MIN(0. , DTM)                                         &
&                      + MAX(0. , (SIDC(IJ,K,M) + SIDC(IJ,K,MPM))*DFP)         &
&                      - MIN(0. , (SIDC(IJ,K,M) + SIDC(IJ,K,MMM))*DFM)
         END DO
      END DO FRE
   END DO DIR

END SUBROUTINE C_SPHER_SHALLOW_CURR

! ---------------------------------------------------------------------------- !

END SUBROUTINE CHECK_CFL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_PROPAGATION_MODULE

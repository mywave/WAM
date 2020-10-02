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
&                             MPM, KPM, JXO,JYO
USE WAM_GRID_MODULE,    ONLY: NSEA, KFROMIJ, KLAT, KLON, WLAT,                 &
&                             DELPHI, DELLAM, SINPH, COSPH, ONE_POINT,         &
&                             REDUCED_GRID, OBSLAT, OBSLON
USE WAM_TIMOPT_MODULE,  ONLY: IDELPRO, IDELT,                                  &
&                             SPHERICAL_RUN, SHALLOW_RUN,                      &
&                             REFRACTION_D_RUN, REFRACTION_C_RUN, L_OBSTRUCTION
USE WAM_FILE_MODULE,    ONLY: IU06, ITEST
USE WAM_MODEL_MODULE,   ONLY: DEPTH, INDEP, U, V
USE WAM_TABLES_MODULE,  ONLY: TCGOND, TFAK, TSIHKD, NDEPTH, DEPTHA, DEPTHD
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

REAL,         ALLOCATABLE :: CGOND(:,:)  !! SHALLOW WATER GROUP VELOCITIES.
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
!     4. WEIGHTS IN ADVECTION SCHEME.                                          !
!        ----------------------------                                          !

REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: COEF_LAT !! N-S DIRECTIONS.
!                            COEF_LAT(IJ,K,M,IC,ICL)
!                            IJ: GRID POINT
!                            K : DIRECTION
!                            M : FREQUENCY
!                            IC : NORTH SOUTH INDEX
!                            ICL : 1 FOR THE CLOSEST GRID POINT
!                                  2 FOR THE SECOND CLOSEST GRID POINT
REAL, ALLOCATABLE, DIMENSION(:,:,:,:)   :: COEF_LON !! E-W DIRECTIONS.
!                            COEF_LON(IJ,K,M,IC)
!                            IJ: GRID POINT
!                            K : DIRECTION
!                            M : FREQUENCY
!                            IC : EAST WEST INDEX
REAL, ALLOCATABLE, DIMENSION(:,:,:,:)   :: COEF_THETA !! REFRACTION TERM
!                            COEF_THETA(IJ,K,M,IC)
!                            IJ: GRID POINT
!                            K : DIRECTION
!                            M : FREQUENCY
!                            IC:  K+-1 INDEX
REAL, ALLOCATABLE, DIMENSION(:,:,:,:)   :: COEF_SIGD  !! FREQUENCY SHIFTING TERM
!                            COEF_SIGD(IJ,K,M,IC)
!                            IJ: GRID POINT
!                            K : DIRECTION
!                            M : FREQUENCY
!                            IC:  M+-1 INDEX
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: COEF_SUM  !! SUM OF WEIGHTS AT CENTRE
!                            COEF_SUM(IJ,K,M)
!                            IJ: GRID POINT
!                            K : DIRECTION
!                            M : FREQUENCY

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. CFL INFORMATION.                                                      !
!        ----------------                                                      !

REAL, PARAMETER   :: XLIMIT = 0.999  !! MAXIMUM ALLOWED CFL NUMBER
REAL, ALLOCATABLE :: CFLMAX(:)       !! MAXIMUM CFL NUMBER on each process.
REAL    :: total_cfl = 0.            !! MAXIMUM CFL NUMBER for all processes.

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

INTERFACE PREPARE_PROPAGATION_COEF    !! PROPAGATION SCHEME WEIGHTS.
   MODULE PROCEDURE PREPARE_PROPAGATION_COEF
END INTERFACE
PRIVATE PREPARE_PROPAGATION_COEF

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

REAL,    INTENT(IN)  :: FIELD(NINF:NSUP)  !! INPUT FIELD.
REAL,    INTENT(OUT) :: D_LAT(nijs:nijl)  !! LATITUDE  DERIVATIVE.
REAL,    INTENT(OUT) :: D_LON(nijs:nijl)  !! LONGITUDE DERIVATIVE.

INTEGER :: NLAND, IJ, IP, IM, IP2, IM2, KX

NLAND = NINF-1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NORTH SOUTH GRADIENTS.                                                !
!        ----------------------                                                !

IF (REDUCED_GRID) THEN
   DO IJ = NIJS,NIJL
      IP  = KLAT(IJ,2,1)
      IM  = KLAT(IJ,1,1)
      IP2 = KLAT(IJ,2,2)
      IM2 = KLAT(IJ,1,2)
      IF (IP.NE.NLAND .AND. IM.NE.NLAND .AND.                                  &
&         IP2.NE.NLAND .AND. IM2.NE.NLAND) THEN
         D_LAT(IJ) = (WLAT(IJ,2)*FIELD(IP)+(1.-WLAT(IJ,2))*FIELD(IP2)          &
&                - WLAT(IJ,1)*FIELD(IM)-(1.-WLAT(IJ,1))*FIELD(IM2))/(2.*DELPHI)
      ELSEIF (IP.NE.NLAND .AND. IM.NE.NLAND ) THEN
         D_LAT(IJ) = (FIELD(IP)-FIELD(IM))/(2.*DELPHI)
      ELSEIF (IP2.NE.NLAND .AND. IM2.NE.NLAND ) THEN
         D_LAT(IJ) = (FIELD(IP2)-FIELD(IM2))/(2.*DELPHI)
      ELSE
         D_LAT(IJ) = 0.0
      ENDIF
   ENDDO

ELSE

   DO IJ = NIJS,NIJL
      IP  = KLAT(IJ,2,1)
      IM  = KLAT(IJ,1,1)
      IF (IP.NE.NLAND .AND. IM.NE.NLAND) THEN
         D_LAT(IJ) = (FIELD(IP) - FIELD(IM))/(2.*DELPHI)
      ELSE
         D_LAT(IJ) = 0.0
      ENDIF
   ENDDO
ENDIF
CALL PUT_DRY (D_LAT, nijs, nijl, 0.)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. EAST WEST GRADIENTS.                                                  !
!        --------------------                                                  !

DO IJ = NIJS,NIJL
   IP = KLON(IJ,2)
   IM = KLON(IJ,1)
   KX  = KFROMIJ(IJ)
   IF (IP.NE.NLAND .AND. IM.NE.NLAND) THEN
      D_LON(IJ) = (FIELD(IP)-FIELD(IM))/(2.*DELLAM(KX))
   ELSE
      D_LON(IJ) = 0.0
   ENDIF
ENDDO

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

INTEGER :: NLAND, IJ, IP, IM, IP2, IM2, IZ
REAL    :: D_MEAN

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
   CALL PREPARE_PROPAGATION_COEF
   IF (ITEST.GE.3) THEN
      WRITE(IU06,*) '    SUB. PREPARE_PROPAGATION: PREPARE_PROPAGATION_COEF DONE'
   END IF
   CALL CHECK_CFL
   IF (ITEST.GE.3) WRITE(IU06,*) '    SUB. PREPARE_PROPAGATION: CHECK_CFL DONE'
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
   CALL GRADIENT (V, DVPHI, DVLAM)
   IF (SPHERICAL_RUN) THEN
      DULAM = DULAM*DCO(nijs:nijl)
      DVLAM = DVLAM*DCO(nijs:nijl)
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. COMPUTE MEAN DEPTH INDEX TO SMOOTH GADIENTS.                          !
!        --------------------------------------------                          !

IF (SHALLOW_RUN) THEN
   ALLOCATE (INDEP_G(NIJS:NIJL))
   NLAND = NINF-1

IF (REDUCED_GRID) THEN

   DO IJ = NIJS,NIJL
      D_MEAN = 4.*DEPTH(IJ)
      IZ = 4
      IP  = KLAT(IJ,2,1)
      IM  = KLAT(IJ,1,1)
      IP2 = KLAT(IJ,2,2)
      IM2 = KLAT(IJ,1,2)
      IF (IP.NE.NLAND .AND. IM.NE.NLAND .AND.                                 &
&         IP2.NE.NLAND .AND. IM2.NE.NLAND) THEN
         D_MEAN = D_MEAN + WLAT(IJ,2)*DEPTH(IP)+(1.-WLAT(IJ,2))*DEPTH(IP2)    &
&               + WLAT(IJ,1)*DEPTH(IM)+(1.-WLAT(IJ,1))*DEPTH(IM2)
         IZ = IZ + 2
      ELSEIF (IP.NE.NLAND .AND. IM.NE.NLAND ) THEN
         D_MEAN = D_MEAN + DEPTH(IP)+DEPTH(IM)
         IZ = IZ + 2
      ELSEIF (IP2.NE.NLAND .AND. IM2.NE.NLAND ) THEN
         D_MEAN = D_MEAN + DEPTH(IP2)+DEPTH(IM2)
         IZ = IZ + 2
      ENDIF

      IP = KLON(IJ,2)
      IM = KLON(IJ,1)
      IF (IP.NE.NLAND .AND. IM.NE.NLAND) THEN
         D_MEAN = D_MEAN + DEPTH(IP)+DEPTH(IM)
         IZ = IZ + 2
      ENDIF

      D_MEAN = D_MEAN/REAL(IZ)
      IZ = NINT(LOG(D_MEAN/DEPTHA)/LOG(DEPTHD)+1.)
      IZ = MAX(IZ, 1)
      INDEP_G(IJ) = MIN(IZ, NDEPTH)
   ENDDO
ELSE
   DO IJ = NIJS,NIJL
      D_MEAN = 4.*DEPTH(IJ)
      IZ = 4
      IP  = KLAT(IJ,2,1)
      IM  = KLAT(IJ,1,1)
      IF (IP.NE.NLAND .AND. IM.NE.NLAND) THEN
         D_MEAN = D_MEAN + DEPTH(IP)+DEPTH(IM)
         IZ = IZ + 2
      ENDIF

      IP = KLON(IJ,2)
      IM = KLON(IJ,1)
      IF (IP.NE.NLAND .AND. IM.NE.NLAND) THEN
         D_MEAN = D_MEAN + DEPTH(IP)+DEPTH(IM)
         IZ = IZ + 2
      ENDIF

      D_MEAN = D_MEAN/REAL(IZ)
      IZ = NINT(LOG(D_MEAN/DEPTHA)/LOG(DEPTHD)+1.)
      IZ = MAX(IZ, 1)
      INDEP_G(IJ) = MIN(IZ, NDEPTH)
   ENDDO
END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. DEPTH OF THETA DOT.                                                   !
!        -------------------                                                   !

IF (SHALLOW_RUN) THEN
   ALLOCATE (TEMP(nijs:nijl))
   IF (.NOT.ALLOCATED(THDD) ) ALLOCATE (THDD(nijs:nijl,1:KL,1:ML))
   DO K = 1,KL
      TEMP(nijs:nijl) = SINTH(K)*DDPHI(nijs:nijl) - COSTH(K)*DDLAM(nijs:nijl)
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
      SS  = SINTH(K)**2
      SC  = SINTH(K)*COSTH(K)
      CC  = COSTH(K)**2

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
&                   * CGOND(nijs:nijl,M)                                       &
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
!     9. PREPARE PROPAGATION COEF.                                             !
!        -------------------------                                             !

CALL PREPARE_PROPAGATION_COEF
IF (ITEST.GE.3) THEN
   WRITE(IU06,*) '    SUB. PREPARE_PROPAGATION: PREPARE_PROPAGATION_COEF DONE'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!    10. CHECK CFL.                                                            !
!        ----------                                                            !

CALL CHECK_CFL
IF (ITEST.GE.3) WRITE(IU06,*) '    SUB. PREPARE_PROPAGATION: CHECK_CFL DONE'

IF (ALLOCATED(THDD) ) DEALLOCATE (THDD)
IF (ALLOCATED(SIDC) ) DEALLOCATE (SIDC)

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
!       *PROP_CURR*         WITH CURRENT REFRACTION.                           !
!       *PROP_NO_CURR*      WITHOUT CURRENT.                                   !
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

INTEGER :: I,ierr

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: f1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. INITIAL.                                                              !
!        --------                                                              !

allocate (f1(ninf-1:nsup,1:kl,1:ml))

DO I = 1, COUNTER

   F1(nijs:nijl,:,:) = F3(:,:,:)           !! COPY INPUT SPECTRA
     
   call mpi_exchng (f1(ninf:nsup,:,:))
   call mpi_barrier (mpi_comm_world, ierr)

   F1(ninf-1,:,:) = 0.                     !! SPECTRUM AT LAND TO ZERO.
    
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.0 SELECT CASES.                                                        !
!         -------------                                                        !

   IF (REFRACTION_C_RUN) THEN
      CALL PROP_CURR            !! WITH CURRENT REFRACTION.
   ELSE
      CALL PROP_NO_CURR         !! WITHOUT CURRENT REFRACTION.
   END IF

END DO

DEALLOCATE (f1)

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
!     3.1 PROPAGATION WITHOUT CURRENT REFRACTION.                              !
!         ---------------------------------------                              !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE PROP_NO_CURR

INTEGER :: IJ, K, M, KP, KM, IS, IC

DIR: DO K = 1,KL
   KP = KPM(K, 1)
   KM = KPM(K,-1)
   IS = JXO(K,1)         !! INDEX FOR ADJOINING LONGITUDE.
   IC = JYO(K,1)         !! INDEX FOR ADJOINING LATITUDE.

   FRE: DO M = 1,ML
      DO IJ = NIJS,NIJL
         F3(IJ,K,M) = COEF_SUM(IJ,K,M)     * F1(IJ,K,M)                        &
&                   + COEF_LAT(IJ,K,M,1,1) * F1(KLAT(IJ,IC,1),K,M)             &
&                   + COEF_LON(IJ,K,M,1)   * F1(KLON(IJ,IS)  ,K,M)

      END DO

      IF (REFRACTION_D_RUN .OR. SPHERICAL_RUN) THEN
         DO IJ = nijs,nijl
            F3(IJ,K,M) = F3(IJ,K,M)                                            &
&                      + COEF_THETA(IJ,K,M,1) *F1(IJ          ,KP ,M)          &
&                      + COEF_THETA(IJ,K,M,2) *F1(IJ          ,KM ,M)
         END DO
      END IF

      IF (REDUCED_GRID) THEN
         DO IJ = NIJS,NIJL
            F3(IJ,K,M) = F3(IJ,K,M)                                            &
&                      + COEF_LAT(IJ,K,M,1,2) * F1(KLAT(IJ,IC,2),K,M)
         END DO
      END IF
   END DO FRE
END DO DIR

END SUBROUTINE PROP_NO_CURR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.2 PROPAGATION WITH CURRENT REFRACTION.                                 !
!         ------------------------------------                                 !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE PROP_CURR

INTEGER :: IJ, K, M, KP, KM, MM, MP

FRE: DO M = 1,ML
   MM = MPM(M,-1)
   MP = MPM(M, 1)
   DIR: DO K = 1,KL
      KP = KPM(K, 1)
      KM = KPM(K,-1)
      DO IJ = NIJS,NIJL
         F3(IJ,K,M) = COEF_SUM(IJ,K,M)     * F1(IJ,K,M )                      &
&                   + COEF_LAT(IJ,K,M,2,1) * F1(KLAT(IJ,2,1),K,M)             &
&                   + COEF_LAT(IJ,K,M,1,1) * F1(KLAT(IJ,1,1),K,M)             &
&                   + COEF_LON(IJ,K,M,2)   * F1(KLON(IJ,2),K,M)               &
&                   + COEF_LON(IJ,K,M,1)   * F1(KLON(IJ,1),K,M)               &
&                   + COEF_THETA(IJ,K,M,1) * F1(IJ,KP,M  )                    &
&                   + COEF_THETA(IJ,K,M,2) * F1(IJ,KM,M  )                    &
&                   + COEF_SIGD(IJ,K,M,1)  * F1(IJ,K  ,MP)                    &
&                   + COEF_SIGD(IJ,K,M,2)  * F1(IJ,K  ,MM)
      END DO
      IF (REDUCED_GRID) THEN
         DO IJ = NIJS,NIJL
            F3(IJ,K,M) = F3(IJ,K,M )                                          &
&                      + COEF_LAT(IJ,K,M,2,2) * F1(KLAT(IJ,2,2),K,M)          &
&                      + COEF_LAT(IJ,K,M,1,2) * F1(KLAT(IJ,1,2),K,M)
         END DO
      END IF
   END DO DIR
END DO FRE

END SUBROUTINE PROP_CURR

! ---------------------------------------------------------------------------- !

END SUBROUTINE PROPAGS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_PROPAGATION_COEF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     PREPARE_PROPAGATION_COEF - COMPUTATION OF THE PROPAGATION SCHEME WEIGHTS.!
!                                                                              !
!*    PURPOSE.
!     --------
!                                                                              !
!       PRE-COMPUTATION OF THE UPSTREAM WEIGHT
!       USED IN THE PROPAGATION FOR A GIVEN TIME STEP.
!                                                                              !
!**   INTERFACE.
!     ----------
!                                                                              !
!     METHOD.
!     -------
!                                                                              !
!     REFERENCE.
!     ----------
!                                                                              !
!       NONE.
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.0 SELECT CASES.                                                        !
!         -------------                                                        !

   IF (SPHERICAL_RUN) THEN

!     1.1 PROPAGATION ON SPHERICAL GRID.                                       !
!         ------------------------------                                       !

      IF (.NOT.REFRACTION_C_RUN) THEN
         IF (SHALLOW_RUN) THEN
            CALL PREP_SPHER_SHALLOW  !! SHALLOW WATER WITHOUT CURRENT REFRACTION.
         ELSE
            CALL PREP_SPHER_DEEP     !! DEEP WATER WITHOUT CURRENT REFRACTION.
         END IF
      ELSE
         IF (SHALLOW_RUN) THEN
            IF (REDUCED_GRID) THEN
               CALL PREP_SPHER_SHALLOW_CURR !! SHALLOW WATER WITH DEPTH AND CURRENT REF.
            ELSE
               CALL PREP_SPHER_SHALLOW_CURR_REG !! SHALLOW WATER WITH DEPTH AND CURRENT REF.
            END IF
         ELSE
            IF (REDUCED_GRID) THEN
               CALL PREP_SPHER_DEEP_CURR    !! DEEP WATER WITH CURRENT REFRACTION.
            ELSE
               CALL PREP_SPHER_DEEP_CURR_REG    !! DEEP WATER WITH CURRENT REFRACTION.
            END IF
         END IF
      END IF

   ELSE

!     1.2 PROPAGATION ON CARTESIAN GRID.                                       !
!         ------------------------------                                       !

      IF (.NOT.REFRACTION_C_RUN) THEN
         IF (SHALLOW_RUN) THEN
            CALL PREP_CART_SHALLOW      !! SHALLOW WATER WITHOUT CURRENT REFRACTION.
         ELSE
            CALL PREP_CART_DEEP         !! DEEP WATER WITHOUT CURRENT REFRACTION.
         END IF
      ELSE
         IF (SHALLOW_RUN) THEN
            CALL PREP_CART_SHALLOW_CUR  !! SHALLOW WATER WITH DEPTH AND CURRENT REF.
         ELSE
            CALL PREP_CART_DEEP_CUR     !! DEEP WATER WITH CURRENT REFRACTION.
         END IF
      END IF
   END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.0 RETURN.                                                              !
!         -------                                                              !

RETURN

CONTAINS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.1 PROPAGATION FOR CARTESIAN GRID WITHOUT CURRENT REFRACTION (DEEP).    !
!         -----------------------------------------------------------------    !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE PREP_CART_DEEP

INTEGER :: IJ, K, M
REAL    :: DELPRO, DELPH0, SD, CD
REAL    :: DPN, DLE
REAL, ALLOCATABLE, DIMENSION(:) :: DELLA0

ALLOCATE(DELLA0(NIJS:NIJL))

DELPRO = REAL(IDELPRO)
DELPH0 = DELPRO/DELPHI

DO IJ = NIJS,NIJL
   DELLA0(IJ) = DELPRO/DELLAM(KFROMIJ(IJ))
END DO

IF (.NOT.ALLOCATED(COEF_LAT))   ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:1,1:1))
IF (.NOT.ALLOCATED(COEF_LON))   ALLOCATE(COEF_LON(NIJS:NIJL,1:KL,1:ML,1:1))
IF (.NOT.ALLOCATED(COEF_SUM))   ALLOCATE(COEF_SUM(NIJS:NIJL,1:KL,1:ML))

FRE: DO M = 1,ML                   !! LOOP OVER FREQUENCIES.
   DIR: DO K = 1,KL
      SD = SINTH(K)
      CD = COSTH(K)

      DPN = ABS(CD*DELPH0*GOM(M))
      SD  = ABS(SD*GOM(M))
      DO IJ = nijs,nijl
         DLE = SD*DELLA0(IJ)
         COEF_SUM(IJ,K,M)     = 1. - DLE - DPN
         COEF_LON(IJ,K,M,1)   = DLE
         COEF_LAT(IJ,K,M,1,1) = DPN
      END DO
   END DO DIR
END DO FRE
DEALLOCATE(DELLA0)

END SUBROUTINE PREP_CART_DEEP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.2 PROPAGATION FOR CARTESIAN GRID WITHOUT CURRENT REFRACTION (SHALLOW). !
!         -------------------------------------------------------------------- !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE PREP_CART_SHALLOW

INTEGER :: IJ, K, M, KP, KM, IS, IC, IT, ID
REAL    :: DELPRO, DELPH0, DELTHR, SD, CD
REAL    :: DTP, DTM
REAL, ALLOCATABLE, DIMENSION(:) :: DELLA0

ALLOCATE(DELLA0(NIJS:NIJL))

DELPRO = REAL(IDELPRO)
DELPH0 = 0.5*DELPRO/DELPHI
DELTHR = 0.5*DELPRO/DELTH

DO IJ = NIJS,NIJL
   DELLA0(IJ) = 0.5*DELPRO/DELLAM(KFROMIJ(IJ))
END DO

IF (.NOT.ALLOCATED(COEF_LAT))   ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:1,1:1))
IF (.NOT.ALLOCATED(COEF_LON))   ALLOCATE(COEF_LON(NIJS:NIJL,1:KL,ML,1:1))
IF (.NOT.ALLOCATED(COEF_SUM))   ALLOCATE(COEF_SUM(NIJS:NIJL,1:KL,1:ML))

IF (REFRACTION_D_RUN) THEN
   IF (.NOT.ALLOCATED(COEF_THETA)) ALLOCATE(COEF_THETA(NIJS:NIJL,1:KL,1:ML,1:2))
END IF

FRE: DO M = 1,ML
   DIR: DO K = 1,KL
      SD = SINTH(K)
      CD = COSTH(K)

      IS = JXO(K,1)         !! INDEX FOR ADJOINING LONGITUDE.
      IT = JXO(K,2)
      IC = JYO(K,1)         !! INDEX FOR ADJOINING LATITUDE.
      ID = JYO(K,2)

      SD = ABS(SD)
      CD = ABS(CD*DELPH0)

      DO IJ = nijs,nijl
         COEF_LON(IJ,K,M,1)   = SD*(CGOND(KLON(IJ,IS),M) + CGOND(IJ,M)) * DELLA0(IJ)
         COEF_LAT(IJ,K,M,1,1) = CD*(CGOND(KLAT(IJ,IC,1),M) + CGOND(IJ,M))
         COEF_SUM(IJ,K,M)     = 1.                                                      &
&                             - SD*(CGOND(KLON(IJ,IT),M) + CGOND(IJ,M)) * DELLA0(IJ)    &
&                             - CD*(CGOND(KLAT(IJ,ID,1),M) + CGOND(IJ,M))
      END DO

      IF (REFRACTION_D_RUN) THEN
         KP = KPM(K, 1)
         KM = KPM(K,-1)
         DO IJ = nijs,nijl
            DTP = (THDD(IJ,K ,M)+THDD(IJ,KP ,M))*DELTHR
            DTM = (THDD(IJ,K ,M)+THDD(IJ,KM ,M))*DELTHR

            COEF_SUM(IJ,K,M) = COEF_SUM(IJ,K,M) - MAX(0.,DTP) + MIN (0.,DTM)
            COEF_THETA(IJ,K,M,1) =  - MIN (0.,DTP)
            COEF_THETA(IJ,K,M,2) =    MAX (0.,DTM)
         END DO
      END IF

   END DO DIR
END DO FRE

IF (L_OBSTRUCTION) THEN
   DO M = 1,ML
      DO K = 1,KL
         IS = JXO(K,1)         !! INDEX FOR ADJOINING LONGITUDE.
         IC = JYO(K,1)         !! INDEX FOR ADJOINING LATITUDE.
         DO IJ = NIJS,NIJL
            COEF_LAT(IJ,K,M,1,1) = COEF_LAT(IJ,K,M,1,1)*OBSLAT(IJ,IC,M)
            COEF_LON(IJ,K,M,1)   = COEF_LON(IJ,K,M,1) * OBSLON(IJ,IS,M)
         END DO
      END DO
   END DO
END IF

DEALLOCATE(DELLA0)

END SUBROUTINE PREP_CART_SHALLOW

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.3 PROPAGATION FOR CARTESIAN GRID WITH CURRENT REFRACTION (DEEP).       !
!         --------------------------------------------------------------       !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE PREP_CART_DEEP_CUR

INTEGER :: IJ, K, M, KP, KM, MP, MM
REAL    :: DELPRO, DELPH0, DELTHR, DELFR0, SD, CD, CGS, CGC
REAL    :: DTC, DPN, DPS, DLE, DLW, DTP, DTM, DOM
REAL, ALLOCATABLE, DIMENSION(:) :: DELLA0
REAL, ALLOCATABLE, DIMENSION(:) :: WOK1
REAL, ALLOCATABLE, DIMENSION(:) :: WOK2

ALLOCATE(DELLA0(NIJS:NIJL))

DELPRO = REAL(IDELPRO)
DELPH0 = 0.5*DELPRO/DELPHI
DELTHR = 0.5*DELPRO/DELTH
DELFR0 = 0.5*DELPRO*2.1/0.2

DO IJ = NIJS,NIJL
   DELLA0(IJ) = 0.5*DELPRO/DELLAM(KFROMIJ(IJ))
END DO

IF (.NOT.ALLOCATED(COEF_LAT))   ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:2,1:1))
IF (.NOT.ALLOCATED(COEF_LON))   ALLOCATE(COEF_LON(NIJS:NIJL,1:KL,ML,1:2))
IF (.NOT.ALLOCATED(COEF_THETA)) ALLOCATE(COEF_THETA(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SIGD))  ALLOCATE(COEF_SIGD(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SUM))   ALLOCATE(COEF_SUM(NIJS:NIJL,1:KL,1:ML))

ALLOCATE(WOK1(ninf-1:nsup))
ALLOCATE(WOK2(ninf-1:nsup))

   FRE: DO M = 1,ML
      MM = MPM(M,-1)
      MP = MPM(M, 1)
      DIR: DO K = 1,KL
         KP = KPM(K, 1)
         KM = KPM(K,-1)

         SD = SINTH(K)
         CD = COSTH(K)

         CGS = GOM(M)*SD                !! GROUP VELOCITIES.
         CGC = GOM(M)*CD

         WOK1(ninf-1:nsup) = U(ninf-1:nsup) + CGS
         WOK2(ninf-1:nsup) = (V(ninf-1:nsup) + CGC)*DELPH0

         DO IJ=NIJS,NIJL
            DLW = (WOK1(IJ) + WOK1(KLON(IJ,1))) * DELLA0(IJ)
            DLE = (WOK1(IJ) + WOK1(KLON(IJ,2))) * DELLA0(IJ)
            DTC = 1. - MAX(0.,DLE) + MIN(0.,DLW)
            COEF_LON(IJ,K,M,2) = -MIN(0.,DLE)
            COEF_LON(IJ,K,M,1) =  MAX(0.,DLW)

            DPS = WOK2(IJ) + WOK2(KLAT(IJ,1,1))
            DPN = WOK2(IJ) + WOK2(KLAT(IJ,2,1))
            DTC = DTC - MAX(0.,DPN) + MIN(0.,DPS)
            COEF_LAT(IJ,K,M,2,1) = -MIN(0.,DPN)
            COEF_LAT(IJ,K,M,1,1) =  MAX(0.,DPS)

            DTP = (THDD(IJ,K  ,ML)+THDD(IJ,KP ,ML))*DELTHR
            DTM = (THDD(IJ,K  ,ML)+THDD(IJ,KM ,ML))*DELTHR
            DTC = DTC - MAX(0.,DTP) + MIN(0.,DTM)
            COEF_THETA(IJ,K,M,1) = -MIN(0.,DTP)
            COEF_THETA(IJ,K,M,2) =  MAX(0.,DTM)

            DOM = SIDC(IJ,K,ML) * DELFR0
            COEF_SUM(IJ,K,M)    = DTC - ABS(DOM)
            COEF_SIGD(IJ,K,M,1) = -MIN(0.,DOM) / 1.1
            COEF_SIGD(IJ,K,M,2) =  MAX(0.,DOM) * 1.1

         END DO
      END DO DIR
   END DO FRE

DEALLOCATE(DELLA0)
DEALLOCATE(WOK1)
DEALLOCATE(WOK2)

END SUBROUTINE PREP_CART_DEEP_CUR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.4 PROPAGATION FOR CARTESIAN GRID WITH CURRENT REFRACTION (SHALLOW).    !
!         -----------------------------------------------------------------    !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE PREP_CART_SHALLOW_CUR

INTEGER :: IJ, K, M, KP, KM, MP, MM
REAL    :: DELPRO, DELPH0, DELTHR, DELFR0, SD, CD
REAL    :: DTC, DPN, DPS, DLE, DLW, DTP, DTM, DOP, DOM, DFP, DFM
REAL, ALLOCATABLE, DIMENSION(:) :: DELLA0
REAL, ALLOCATABLE, DIMENSION(:) :: WOK1
REAL, ALLOCATABLE, DIMENSION(:) :: WOK2

ALLOCATE(DELLA0(NIJS:NIJL))

DELPRO = REAL(IDELPRO)
DELPH0 = 0.5*DELPRO/DELPHI
DELTHR = 0.5*DELPRO/DELTH
DELFR0 = 0.5*DELPRO/(0.1*ZPI)

DO IJ = NIJS,NIJL
   DELLA0(IJ) = 0.5*DELPRO/DELLAM(KFROMIJ(IJ))
END DO

IF (.NOT.ALLOCATED(COEF_LAT))   ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:2,1:1))
IF (.NOT.ALLOCATED(COEF_LON))   ALLOCATE(COEF_LON(NIJS:NIJL,1:KL,ML,1:2))
IF (.NOT.ALLOCATED(COEF_THETA)) ALLOCATE(COEF_THETA(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SIGD))  ALLOCATE(COEF_SIGD(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SUM))   ALLOCATE(COEF_SUM(NIJS:NIJL,1:KL,1:ML))

ALLOCATE(WOK1(ninf-1:nsup))
ALLOCATE(WOK2(ninf-1:nsup))

   FRE: DO M = 1,ML
      MM = MPM(M,-1)
      MP = MPM(M, 1)
      DFP = DELFR0/FR(M)
      DFM = DELFR0/FR(MM)
      DIR: DO K = 1,KL
         SD = SINTH(K)
         CD = COSTH(K)

         KP = KPM(K, 1)
         KM = KPM(K,-1)

         DO IJ = ninf-1,nsup
            WOK1(IJ) =  U(IJ) + SD*CGOND(IJ,M)
            WOK2(IJ) = (V(IJ) + CD*CGOND(IJ,M)) * DELPH0
         END DO

         DO IJ = nijs,nijl
            DLW = (WOK1(IJ) + WOK1(KLON(IJ,1))) * DELLA0(IJ)
            DLE = (WOK1(IJ) + WOK1(KLON(IJ,2))) * DELLA0(IJ)
            DTC = 1. - MAX(0.,DLE) + MIN(0.,DLW)
            COEF_LON(IJ,K,M,1) =  MAX(0.,DLW)
            COEF_LON(IJ,K,M,2) = -MIN(0.,DLE)

            DPS = WOK2(IJ) + WOK2(KLAT(IJ,1,1))
            DPN = WOK2(IJ) + WOK2(KLAT(IJ,2,1))
            DTC = DTC - MAX(0.,DPN) + MIN(0.,DPS)
            COEF_LAT(IJ,K,M,1,1) =  MAX(0.,DPS)
            COEF_LAT(IJ,K,M,2,1) = -MIN(0.,DPN)

            DTP = (THDD(IJ,K  ,M)+THDD(IJ,KP ,M))*DELTHR
            DTM = (THDD(IJ,K  ,M)+THDD(IJ,KM ,M))*DELTHR
            DTC = DTC - MAX(0.,DTP) + MIN(0.,DTM)
            COEF_THETA(IJ,K,M,1) = -MIN(0.,DTP)
            COEF_THETA(IJ,K,M,2) =  MAX(0.,DTM)

            DOP = (SIDC(IJ,K,M) + SIDC(IJ,K,MP))*DFP
            DOM = (SIDC(IJ,K,M) + SIDC(IJ,K,MM))*DFM
            COEF_SIGD(IJ,K,M,1) = -MIN(0.,DOP) / 1.1
            COEF_SIGD(IJ,K,M,2) =  MAX(0.,DOM) * 1.1

            COEF_SUM(IJ,K,M)    = DTC - MAX(0.,DOP) + MIN(0.,DOM)

         END DO
      END DO DIR
   END DO FRE

IF (L_OBSTRUCTION) THEN
   DO M = 1,ML
      DO K = 1,KL
         DO IJ = NIJS,NIJL
            COEF_LAT(IJ,K,M,1,1) = COEF_LAT(IJ,K,M,1,1)*OBSLAT(IJ,1,M)
            COEF_LAT(IJ,K,M,2,1) = COEF_LAT(IJ,K,M,2,1)*OBSLAT(IJ,2,M)
            COEF_LON(IJ,K,M,1)   = COEF_LON(IJ,K,M,1)  *OBSLON(IJ,1,M)
            COEF_LON(IJ,K,M,2)   = COEF_LON(IJ,K,M,2)  *OBSLON(IJ,2,M)
         END DO
      END DO
   END DO
END IF

DEALLOCATE(DELLA0)
DEALLOCATE(WOK1)
DEALLOCATE(WOK2)

END SUBROUTINE PREP_CART_SHALLOW_CUR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.5 PROPAGATION FOR SPHERICAL GRID WITHOUT CURRENT REFRACTION (DEEP).    !
!         -----------------------------------------------------------------    !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE PREP_SPHER_DEEP

INTEGER :: IJ, K, M, KP, KM, IC, ID
REAL    :: DELPRO, DELTH0, DELPH0, SD, CD, SDA, CDA, SP, SM
REAL    :: DTC, DPN, DLE, DTP, DTM

REAL, ALLOCATABLE, DIMENSION(:) :: DELLA0
REAL, ALLOCATABLE, DIMENSION(:) :: DRGM
REAL, ALLOCATABLE, DIMENSION(:,:) :: WLATM1

DELPRO = REAL(IDELPRO)
DELTH0 = 0.5*DELPRO/DELTR
DELPH0 = 0.5*DELPRO/DELPHI

ALLOCATE(DELLA0(NIJS:NIJL))
ALLOCATE(DRGM(NIJS:NIJL))

DO IJ = NIJS,NIJL
   DELLA0(IJ) = DELPRO*DCO(IJ)/DELLAM(KFROMIJ(IJ))
   DRGM(IJ)   = SINPH(KFROMIJ(IJ))*DCO(IJ)
END DO

IF (REDUCED_GRID) THEN
   IF (.NOT.ALLOCATED(COEF_LAT))   ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:1,1:2))
   ALLOCATE(WLATM1(NIJS:NIJL,2))
   WLATM1(NIJS:NIJL,:) = 1.- WLAT(NIJS:NIJL,:)
ELSE
   IF (.NOT.ALLOCATED(COEF_LAT))   ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:1,1:1))
ENDIF
IF (.NOT.ALLOCATED(COEF_LON))   ALLOCATE(COEF_LON(NIJS:NIJL,1:KL,ML,1:1))
IF (.NOT.ALLOCATED(COEF_THETA)) ALLOCATE(COEF_THETA(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SUM))   ALLOCATE(COEF_SUM(NIJS:NIJL,1:KL,1:ML))

DIR: DO K = 1,KL
   SD = SINTH(K)
   CD = COSTH(K)
   SDA = ABS(SD)
   CDA = ABS(CD*DELPH0)

   KP = KPM(K, 1)
   KM = KPM(K,-1)

   IC = JYO(K,1)         !! INDEX FOR ADJOINING LATITUDE.
   ID = JYO(K,2)

   SP  = DELTH0*(SINTH(K)+SINTH(KP))         !! GRID REFRACTION WEIGHTS.
   SM  = DELTH0*(SINTH(K)+SINTH(KM))

   FRE: DO M = 1,ML

      DO IJ = NIJS,NIJL
         DPN = CDA*(DPSN(IJ,IC) + 1.)  !! LATITUDE WEIGHTS.
         DLE = SDA*DELLA0(IJ)          !! LONGITUDE WEIGHTS.
         DTC = CDA*(DPSN(IJ,ID) + 1.) + DLE

         DTP = DRGM(IJ)*SP
         DTM = DRGM(IJ)*SM
         DTC = DTC + MAX(0. , DTP) - MIN(0. , DTM)
         COEF_LAT(IJ,K,M,1,1)  =  GOM(M)*DPN
         COEF_LON(IJ,K,M,1)    =  GOM(M)*DLE
         COEF_THETA(IJ,K,M,1)  = -GOM(M)*MIN(0. , DTP)
         COEF_THETA(IJ,K,M,2)  =  GOM(M)*MAX(0. , DTM)

         COEF_SUM(IJ,K,M)      = 1. - GOM(M)*DTC

      END DO
      IF (REDUCED_GRID) THEN
         DO IJ = NIJS,NIJL
            COEF_LAT(IJ,K,M,1,2) = COEF_LAT(IJ,K,M,1,1)*WLATM1(IJ,IC)
            COEF_LAT(IJ,K,M,1,1) = COEF_LAT(IJ,K,M,1,1)*WLAT(IJ,IC)
         END DO
      END IF
   END DO FRE
END DO DIR

DEALLOCATE(DELLA0)
DEALLOCATE(DRGM)
IF (ALLOCATED(WLATM1)) DEALLOCATE(WLATM1)

END SUBROUTINE PREP_SPHER_DEEP

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.6 PROPAGATION FOR SPHERICAL GRID WITHOUT CURRENT REFRACTION (SHALLOW). !
!         -------------------------------------------------------------------- !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE PREP_SPHER_SHALLOW

INTEGER :: IJ, K, M, KP, KM, IS, IT, IC, ID
REAL    :: DELPRO, DELTH0, DELPH0, DELTHR, SD, CD, SDA, CDA, SP, SM
REAL    :: DTC, DTP, DTM

REAL, ALLOCATABLE, DIMENSION(:)     :: DELLA0
REAL, ALLOCATABLE, DIMENSION(:)     :: DRGM
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: CGKLAT
REAL, ALLOCATABLE, DIMENSION(:,:)   :: WLATM1

DELPRO = REAL(IDELPRO)
DELTH0 = 0.5*DELPRO/DELTR
DELPH0 = 0.5*DELPRO/DELPHI
DELTHR = 0.5*DELPRO/DELTH

ALLOCATE(DELLA0(NIJS:NIJL))
ALLOCATE(DRGM(NIJS:NIJL))

DO IJ = NIJS,NIJL
   DELLA0(IJ) = 0.5*DELPRO*DCO(IJ)/DELLAM(KFROMIJ(IJ))
   DRGM(IJ)   = SINPH(KFROMIJ(IJ))*DCO(IJ)
END DO

ALLOCATE(CGKLAT(NIJS:NIJL,ML,2))
IF (REDUCED_GRID) THEN
   IF (.NOT.ALLOCATED(COEF_LAT)) ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:1,1:2))
   ALLOCATE(WLATM1(NIJS:NIJL,2))
   WLATM1(NIJS:NIJL,:) = 1.- WLAT(NIJS:NIJL,:)

   DO IC = 1,2
      DO M = 1,ML
         DO IJ = NIJS,NIJL
            CGKLAT(IJ,M,IC) = CGOND(IJ,M) + DPSN(IJ,IC)*                       &
&                            (WLAT(IJ,IC)  *CGOND(KLAT(IJ,IC,1),M) +           &
&                             WLATM1(IJ,IC)*CGOND(KLAT(IJ,IC,2),M))
         ENDDO
      ENDDO
   ENDDO
ELSE
   IF (.NOT.ALLOCATED(COEF_LAT)) ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:1,1:1))
   DO IC = 1,2
      DO M = 1,ML
         DO IJ = NIJS,NIJL
            CGKLAT(IJ,M,IC) = CGOND(IJ,M) + DPSN(IJ,IC)*CGOND(KLAT(IJ,IC,1),M)
         END DO
      END DO
   END DO
END IF

IF (.NOT.ALLOCATED(COEF_LON))   ALLOCATE(COEF_LON(NIJS:NIJL,1:KL,ML,1:1))
IF (.NOT.ALLOCATED(COEF_THETA)) ALLOCATE(COEF_THETA(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SUM))   ALLOCATE(COEF_SUM(NIJS:NIJL,1:KL,1:ML))


DIR: DO K = 1,KL
   SD = SINTH(K)
   CD = COSTH(K)
   SDA = ABS(SD)
   CDA = ABS(CD*DELPH0)

   KP = KPM(K, 1)
   KM = KPM(K,-1)
   SP  = DELTH0*(SINTH(K)+SINTH(KP))   !! PRE_COMPUTE GRID REFRACTION.
   SM  = DELTH0*(SINTH(K)+SINTH(KM))

   IS = JXO(K,1)         !! INDEX FOR ADJOINING LONGITUDE.
   IT = JXO(K,2)
   IC = JYO(K,1)         !! INDEX FOR ADJOINING LATITUDE.
   ID = JYO(K,2)

   IF (REFRACTION_D_RUN) THEN                 !! DEPTH REFRACTION
      FRE1: DO M = 1,ML
         DO IJ = NIJS,NIJL
            DTC = 1. - CDA * CGKLAT(IJ,M,ID)
            COEF_LAT(IJ,K,M,1,1) = CDA * CGKLAT(IJ,M,IC)

            DTC = DTC - SDA*(CGOND(KLON(IJ,IT),M) + CGOND(IJ,M)) * DELLA0(IJ)
            COEF_LON(IJ,K,M,1) = SDA*(CGOND(KLON(IJ,IS),M) + CGOND(IJ,M))       &
&                             *DELLA0(IJ)

            DTP = DRGM(IJ)*CGOND(IJ,M)           !! REFRACTION WEIGHTS
            DTM = DTP*SM + (THDD(IJ,K ,M)+THDD(IJ,KM ,M))*DELTHR
            DTP = DTP*SP + (THDD(IJ,K ,M)+THDD(IJ,KP ,M))*DELTHR
            DTC = DTC - MAX(0. , DTP) + MIN(0. , DTM)
            COEF_THETA(IJ,K,M,1) = -MIN(0. , DTP)
            COEF_THETA(IJ,K,M,2) =  MAX(0. , DTM)

            COEF_SUM(IJ,K,M)     = DTC

         END DO
      END DO FRE1

   ELSE

      FRE2: DO M = 1,ML
         DO IJ = NIJS,NIJL
            DTC = 1. - CDA * CGKLAT(IJ,M,ID)
            COEF_LAT(IJ,K,M,1,1) = CDA * CGKLAT(IJ,M,IC)

            DTC = DTC - SDA*(CGOND(KLON(IJ,IT),M) + CGOND(IJ,M)) *DELLA0(IJ)
            COEF_LON(IJ,K,M,1) = SDA*(CGOND(KLON(IJ,IS),M) + CGOND(IJ,M))     &
&                             *DELLA0(IJ)

            DTP = DRGM(IJ)*CGOND(IJ,M) !! REFRACTION WEIGHTS
            DTM = DTP*SM
            DTP = DTP*SP

            DTC = DTC - MAX(0. , DTP) + MIN(0. , DTM)
            COEF_THETA(IJ,K,M,1) = -MIN(0. , DTP)
            COEF_THETA(IJ,K,M,2) =  MAX(0. , DTM)

            COEF_SUM(IJ,K,M)      = DTC

         END DO
      END DO FRE2
   END IF
END DO DIR

IF (L_OBSTRUCTION) THEN
   DO M = 1,ML
      DO K = 1,KL
         IS = JXO(K,1)         !! INDEX FOR ADJOINING LONGITUDE.
         IC = JYO(K,1)         !! INDEX FOR ADJOINING LATITUDE.
         DO IJ = NIJS,NIJL
            COEF_LAT(IJ,K,M,1,1) = COEF_LAT(IJ,K,M,1,1)*OBSLAT(IJ,IC,M)
            COEF_LON(IJ,K,M,1)   = COEF_LON(IJ,K,M,1) * OBSLON(IJ,IS,M)
         END DO
      END DO
   END DO
END IF

IF (REDUCED_GRID) THEN
   DO M = 1,ML
      DO K = 1,KL
         IC = JYO(K,1)         !! INDEX FOR ADJOINING LATITUDE.
         DO IJ = NIJS,NIJL
            COEF_LAT(IJ,K,M,1,2)  = COEF_LAT(IJ,K,M,1,1) * WLATM1(IJ,IC)
            COEF_LAT(IJ,K,M,1,1)  = COEF_LAT(IJ,K,M,1,1) * WLAT(IJ,IC)
         END DO
      END DO
   END DO
END IF

DEALLOCATE(CGKLAT)
DEALLOCATE(DELLA0)
DEALLOCATE(DRGM)
IF (ALLOCATED(WLATM1))DEALLOCATE(WLATM1)

END SUBROUTINE PREP_SPHER_SHALLOW

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.7 PROPAGATION FOR SPHERICAL GRID WITH CURRENT REFRACTION (DEEP).       !
!         --------------------------------------------------------------       !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE PREP_SPHER_DEEP_CURR

INTEGER :: IJ, K, M, KP, KM, MM, MP
REAL    :: DELPRO, DELTH0, DELPH0, DELTHR, DELFR0, CGS, CGC, SD, CD, SP, SM
REAL    :: DTC, DPN, DPN2, DPS, DPS2, DLE, DLW, DTP, DTM, DOM

REAL, ALLOCATABLE, DIMENSION(:)     :: DELLA0
REAL, ALLOCATABLE, DIMENSION(:)     :: DRGM
REAL, ALLOCATABLE, DIMENSION(:)     :: WOK1
REAL, ALLOCATABLE, DIMENSION(:)     :: WOK2
REAL, ALLOCATABLE, DIMENSION(:,:)   :: WLATM1

ALLOCATE(WLATM1(NIJS:NIJL,2))
WLATM1(NIJS:NIJL,:) = 1.- WLAT(NIJS:NIJL,:)

DELPRO = REAL(IDELPRO)
DELPH0 = 0.5*DELPRO/DELPHI
DELTH0 = 0.5*DELPRO/DELTR
DELTHR = 0.5*DELPRO/DELTH
DELFR0 = 0.5*DELPRO*2.1/0.2

ALLOCATE(DELLA0(NIJS:NIJL))
ALLOCATE(DRGM(NIJS:NIJL))
DO IJ = NIJS,NIJL
   DELLA0(IJ) = 0.5*DELPRO*DCO(IJ)/DELLAM(KFROMIJ(IJ))
   DRGM(IJ)   = SINPH(KFROMIJ(IJ))*DCO(IJ)
END DO

ALLOCATE(WOK1(NINF-1:NSUP))
ALLOCATE(WOK2(NINF-1:NSUP))

IF (.NOT.ALLOCATED(COEF_LAT))   ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:2,1:2))
IF (.NOT.ALLOCATED(COEF_LON))   ALLOCATE(COEF_LON(NIJS:NIJL,1:KL,ML,1:2))
IF (.NOT.ALLOCATED(COEF_THETA)) ALLOCATE(COEF_THETA(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SIGD))  ALLOCATE(COEF_SIGD(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SUM))   ALLOCATE(COEF_SUM(NIJS:NIJL,1:KL,1:ML))

FRE: DO M = 1,ML
   MM = MPM(M,-1)
   MP = MPM(M, 1)
   DIR: DO K = 1,KL
      SD = SINTH(K)
      CD = COSTH(K)
      KP = KPM(K, 1)
      KM = KPM(K,-1)
      SP = DELTH0*(SINTH(K)+SINTH(KP))*GOM(M)     !! PRE-COMPUTE GRID REFRACTION.
      SM = DELTH0*(SINTH(K)+SINTH(KM))*GOM(M)

      CGS = GOM(M)*SD                       !! GROUP VELOCITIES.
      CGC = GOM(M)*CD

      WOK1(ninf-1:nsup) = U(ninf-1:nsup) + CGS
      WOK2(ninf-1:nsup) = (V(ninf-1:nsup) + CGC) * DELPH0

      DO IJ = NIJS,NIJL
         DLW = (WOK1(IJ) + WOK1(KLON(IJ,1))) * DELLA0(IJ)
         DLE = (WOK1(IJ) + WOK1(KLON(IJ,2))) * DELLA0(IJ)
         DTC = 1. - MAX(0.,DLE) + MIN(0.,DLW)
         COEF_LON(IJ,K,M,2) = -MIN(0.,DLE)
         COEF_LON(IJ,K,M,1) =  MAX(0.,DLW)

         DPS  = (WOK2(IJ)+WOK2(KLAT(IJ,1,1))*DPSN(IJ,1))*WLAT(IJ,1)
         DPS2 = (WOK2(IJ)+WOK2(KLAT(IJ,1,2))*DPSN(IJ,1))*WLATM1(IJ,1)
         DPN  = (WOK2(IJ)+WOK2(KLAT(IJ,2,1))*DPSN(IJ,2))*WLAT(IJ,2)
         DPN2 = (WOK2(IJ)+WOK2(KLAT(IJ,2,2))*DPSN(IJ,2))*WLATM1(IJ,2)
         DTC = DTC - MAX(0.,DPN) - MAX(0.,DPN2) + MIN(0.,DPS) + MIN(0.,DPS2)
         COEF_LAT(IJ,K,M,2,1) = -MIN(0.,DPN)
         COEF_LAT(IJ,K,M,2,2) = -MIN(0.,DPN2)
         COEF_LAT(IJ,K,M,1,1) =  MAX(0.,DPS)
         COEF_LAT(IJ,K,M,1,2) =  MAX(0.,DPS2)

         DTP = SP*DRGM(IJ) + (THDD(IJ,K ,ML)+THDD(IJ,KP,ML))*DELTHR
         DTM = SM*DRGM(IJ) + (THDD(IJ,K ,ML)+THDD(IJ,KM,ML))*DELTHR
         DTC = DTC - MAX(0.,DTP) + MIN(0.,DTM)
         COEF_THETA(IJ,K,M,1) = -MIN(0.,DTP)
         COEF_THETA(IJ,K,M,2) =  MAX(0.,DTM)

         DOM =  SIDC(IJ,K,ML) * DELFR0
         DTC =  DTC - ABS(DOM)
         COEF_SIGD(IJ,K,M,1) = -MIN(0.,DOM)/1.1
         COEF_SIGD(IJ,K,M,2) =  MAX(0.,DOM)*1.1

         COEF_SUM(IJ,K,M)      = DTC
      END DO
   END DO DIR
END DO FRE

DEALLOCATE(WLATM1)
DEALLOCATE(WOK1)
DEALLOCATE(WOK2)
DEALLOCATE(DELLA0)
DEALLOCATE(DRGM)

END SUBROUTINE PREP_SPHER_DEEP_CURR

! ---------------------------------------------------------------------------- !

SUBROUTINE PREP_SPHER_DEEP_CURR_REG

INTEGER :: IJ, K, M, KP, KM, MM, MP
REAL    :: DELPRO, DELTH0, DELPH0, DELTHR, DELFR0, CGS, CGC, SD, CD, SP, SM
REAL    :: DTC, DPN, DPS, DLE, DLW, DTP, DTM, DOM

REAL, ALLOCATABLE, DIMENSION(:)     :: DELLA0
REAL, ALLOCATABLE, DIMENSION(:)     :: DRGM
REAL, ALLOCATABLE, DIMENSION(:)     :: WOK1
REAL, ALLOCATABLE, DIMENSION(:)     :: WOK2

DELPRO = REAL(IDELPRO)
DELPH0 = 0.5*DELPRO/DELPHI
DELTH0 = 0.5*DELPRO/DELTR
DELTHR = 0.5*DELPRO/DELTH
DELFR0 = 0.5*DELPRO*2.1/0.2

ALLOCATE(DELLA0(NIJS:NIJL))
ALLOCATE(DRGM(NIJS:NIJL))
DO IJ = NIJS,NIJL
   DELLA0(IJ) = 0.5*DELPRO*DCO(IJ)/DELLAM(KFROMIJ(IJ))
   DRGM(IJ)   = SINPH(KFROMIJ(IJ))*DCO(IJ)
END DO

ALLOCATE(WOK1(NINF-1:NSUP))
ALLOCATE(WOK2(NINF-1:NSUP))

IF (.NOT.ALLOCATED(COEF_LAT))   ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:2,1:1))
IF (.NOT.ALLOCATED(COEF_LON))   ALLOCATE(COEF_LON(NIJS:NIJL,1:KL,ML,1:2))
IF (.NOT.ALLOCATED(COEF_THETA)) ALLOCATE(COEF_THETA(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SIGD))  ALLOCATE(COEF_SIGD(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SUM))   ALLOCATE(COEF_SUM(NIJS:NIJL,1:KL,1:ML))

FRE: DO M = 1,ML
   MM = MPM(M,-1)
   MP = MPM(M, 1)
   DIR: DO K = 1,KL
      SD = SINTH(K)
      CD = COSTH(K)
      KP = KPM(K, 1)
      KM = KPM(K,-1)
      SP = DELTH0*(SINTH(K)+SINTH(KP))*GOM(M)     !! PRE-COMPUTE GRID REFRACTION.
      SM = DELTH0*(SINTH(K)+SINTH(KM))*GOM(M)

      CGS = GOM(M)*SD                       !! GROUP VELOCITIES.
      CGC = GOM(M)*CD

      WOK1(ninf-1:nsup) = U(ninf-1:nsup) + CGS
      WOK2(ninf-1:nsup) = (V(ninf-1:nsup) + CGC) * DELPH0

      DO IJ = NIJS,NIJL
         DLW = (WOK1(IJ) + WOK1(KLON(IJ,1))) * DELLA0(IJ)
         DLE = (WOK1(IJ) + WOK1(KLON(IJ,2))) * DELLA0(IJ)
         DTC = 1. - MAX(0.,DLE) + MIN(0.,DLW)
         COEF_LON(IJ,K,M,1) =  MAX(0.,DLW)
         COEF_LON(IJ,K,M,2) = -MIN(0.,DLE)

         DPS  = (WOK2(IJ)+WOK2(KLAT(IJ,1,1))*DPSN(IJ,1))
         DPN  = (WOK2(IJ)+WOK2(KLAT(IJ,2,1))*DPSN(IJ,2))
         DTC = DTC - MAX(0.,DPN) + MIN(0.,DPS)
         COEF_LAT(IJ,K,M,1,1) =  MAX(0.,DPS)
         COEF_LAT(IJ,K,M,2,1) = -MIN(0.,DPN)

         DTP = SP*DRGM(IJ) + (THDD(IJ,K  ,ML)+THDD(IJ,KP,ML))*DELTHR
         DTM = SM*DRGM(IJ) + (THDD(IJ,K  ,ML)+THDD(IJ,KM,ML))*DELTHR
         DTC = DTC - MAX(0.,DTP) + MIN(0.,DTM)
         COEF_THETA(IJ,K,M,1) = -MIN(0.,DTP)
         COEF_THETA(IJ,K,M,2) =  MAX(0.,DTM)

         DOM =  SIDC(IJ,K,ML) * DELFR0
         DTC =  DTC - ABS(DOM)
         COEF_SIGD(IJ,K,M,1) = -MIN(0.,DOM)/1.1
         COEF_SIGD(IJ,K,M,2) =  MAX(0.,DOM)*1.1

         COEF_SUM(IJ,K,M) = DTC

      END DO
   END DO DIR
END DO FRE

DEALLOCATE(WOK1)
DEALLOCATE(WOK2)
DEALLOCATE(DELLA0)
DEALLOCATE(DRGM)

END SUBROUTINE PREP_SPHER_DEEP_CURR_REG

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.8 PROPAGATION FOR SPHERICAL GRID WITH CURRENT REFRACTION (SHALLOW).    !
!         -----------------------------------------------------------------    !
!                                                                              !
! ---------------------------------------------------------------------------- !

SUBROUTINE PREP_SPHER_SHALLOW_CURR

INTEGER :: IJ, K, M, KP, KM, MM, MP, IC
REAL    :: DELPRO, DELTH0, DELPH0, DELTHR, DELFR0, SD, CD, SP, SM
REAL    :: DTC, DPN, DPN2, DPS, DPS2, DLE, DLW, DTP, DTM, DOP, DOM, DFP, DFM

REAL, ALLOCATABLE, DIMENSION(:)   :: DELLA0
REAL, ALLOCATABLE, DIMENSION(:)   :: DRGM
REAL, ALLOCATABLE, DIMENSION(:)   :: WOK1
REAL, ALLOCATABLE, DIMENSION(:)   :: WOK2
REAL, ALLOCATABLE, DIMENSION(:,:) :: WLATM1

ALLOCATE(WLATM1(NIJS:NIJL,2))
WLATM1(NIJS:NIJL,:) = 1.- WLAT(NIJS:NIJL,:)

DELPRO = REAL(IDELPRO)
DELPH0 = 0.5*DELPRO/DELPHI
DELTH0 = 0.5*DELPRO/DELTR
DELTHR = 0.5*DELPRO/DELTH
DELFR0 = 0.5*DELPRO/(0.1*ZPI)

ALLOCATE(DELLA0(NIJS:NIJL))
ALLOCATE(DRGM(NIJS:NIJL))
DO IJ = NIJS,NIJL
   DELLA0(IJ) = 0.5*DELPRO*DCO(IJ)/DELLAM(KFROMIJ(IJ))
   DRGM(IJ)   = SINPH(KFROMIJ(IJ))*DCO(IJ)
END DO

ALLOCATE(WOK1(ninf-1:nsup))
ALLOCATE(WOK2(ninf-1:nsup))

IF (.NOT.ALLOCATED(COEF_LAT))   ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:2,1:2))
IF (.NOT.ALLOCATED(COEF_LON))   ALLOCATE(COEF_LON(NIJS:NIJL,1:KL,ML,1:2))
IF (.NOT.ALLOCATED(COEF_THETA)) ALLOCATE(COEF_THETA(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SIGD))  ALLOCATE(COEF_SIGD(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SUM))   ALLOCATE(COEF_SUM(NIJS:NIJL,1:KL,1:ML))

FRE: DO M = 1,ML
   MM = MPM(M,-1)
   MP = MPM(M, 1)
   DFP = DELFR0/FR(M)
   DFM = DELFR0/FR(MM)
   DIR: DO K = 1,KL
      SD = SINTH(K)
      CD = COSTH(K)
      KP = KPM(K, 1)
      KM = KPM(K,-1)
      SP = DELTH0*(SINTH(K)+SINTH(KP))      !! GRID REFRACTION.
      SM = DELTH0*(SINTH(K)+SINTH(KM))

      WOK1(ninf-1:nsup) = U(ninf-1:nsup)+SD*CGOND(ninf-1:nsup,M)
      WOK2(ninf-1:nsup) = (V(ninf-1:nsup)+CD*CGOND(ninf-1:nsup,M))*DELPH0

      DO IJ = NIJS,NIJL

         DLW = ( WOK1(IJ) + WOK1(KLON(IJ,1)) ) * DELLA0(IJ)
         DLE = ( WOK1(IJ) + WOK1(KLON(IJ,2)) ) * DELLA0(IJ)
         DTC = 1. - MAX(0.,DLE) + MIN(0.,DLW)
         COEF_LON(IJ,K,M,1) =  MAX(0.,DLW)
         COEF_LON(IJ,K,M,2) = -MIN(0.,DLE)

         ic = 1
         DPS  = wok2(IJ) + DPSN(IJ,IC)*                                    &
&                           (WLAT(IJ,IC)  *wok2(KLAT(IJ,IC,1)) +           &
&                            WLATM1(IJ,IC)*wok2(KLAT(IJ,IC,2)))
         DPS2 = DPS*WLATM1(IJ,1)
         DPS  = DPS*WLAT(IJ,1)

         ic = 2
         DPN = wok2(IJ) + DPSN(IJ,IC)*                                     &
&                           (WLAT(IJ,IC)  *wok2(KLAT(IJ,IC,1)) +           &
&                            WLATM1(IJ,IC)*wok2(KLAT(IJ,IC,2)))
         DPN2 = DPN*WLATM1(IJ,2)
         DPN  = DPN*WLAT(IJ,2)

         DTC = DTC - MAX(0.,DPN) - MAX(0.,DPN2) + MIN(0.,DPS) + MIN(0.,DPS2)
         COEF_LAT(IJ,K,M,1,1) =  MAX(0.,DPS)
         COEF_LAT(IJ,K,M,1,2) =  MAX(0.,DPS2)
         COEF_LAT(IJ,K,M,2,1) = -MIN(0.,DPN)
         COEF_LAT(IJ,K,M,2,2) = -MIN(0.,DPN2)

         DTP = SP*DRGM(IJ)*CGOND(IJ,M) + (THDD(IJ,K ,M)+THDD(IJ,KP,M))*DELTHR
         DTM = SM*DRGM(IJ)*CGOND(IJ,M) + (THDD(IJ,K ,M)+THDD(IJ,KM,M))*DELTHR
         DTC = DTC - MAX(0.,DTP) + MIN(0.,DTM)
         COEF_THETA(IJ,K,M,1) = -MIN(0.,DTP)
         COEF_THETA(IJ,K,M,2) =  MAX(0.,DTM)

         DOP = (SIDC(IJ,K,M) + SIDC(IJ,K,MP))*DFP
         DOM = (SIDC(IJ,K,M) + SIDC(IJ,K,MM))*DFM
         DTC = DTC - MAX(0.,DOP) + MIN(0.,DOM)
         COEF_SIGD(IJ,K,M,1) = -MIN(0.,DOP)/1.1
         COEF_SIGD(IJ,K,M,2) =  MAX(0.,DOM)*1.1

         COEF_SUM(IJ,K,M) = DTC
      END DO
   END DO DIR
END DO FRE

IF (L_OBSTRUCTION) THEN
   DO M = 1,ML
      DO K = 1,KL
         DO IJ = NIJS,NIJL
            COEF_LAT(IJ,K,M,1,1) = COEF_LAT(IJ,K,M,1,1)*OBSLAT(IJ,1,M)
            COEF_LAT(IJ,K,M,1,2) = COEF_LAT(IJ,K,M,1,2)*OBSLAT(IJ,1,M)
            COEF_LAT(IJ,K,M,2,1) = COEF_LAT(IJ,K,M,2,1)*OBSLAT(IJ,2,M)
            COEF_LAT(IJ,K,M,2,2) = COEF_LAT(IJ,K,M,2,2)*OBSLAT(IJ,2,M)
            COEF_LON(IJ,K,M,1)   = COEF_LON(IJ,K,M,1)  *OBSLON(IJ,1,M)
            COEF_LON(IJ,K,M,2)   = COEF_LON(IJ,K,M,2)  *OBSLON(IJ,2,M)
         END DO
      END DO
   END DO
END IF

DEALLOCATE(WLATM1)
DEALLOCATE(DELLA0)
DEALLOCATE(DRGM)
DEALLOCATE(WOK1)
DEALLOCATE(WOK2)

END SUBROUTINE PREP_SPHER_SHALLOW_CURR

! ---------------------------------------------------------------------------- !

SUBROUTINE PREP_SPHER_SHALLOW_CURR_REG

INTEGER :: IJ, K, M, KP, KM, MM, MP, IC
REAL    :: DELPRO, DELTH0, DELPH0, DELTHR, DELFR0, SD, CD, SP, SM
REAL    :: DTC, DPN, DPS, DLE, DLW, DTP, DTM, DOP, DOM, DFP, DFM

REAL, ALLOCATABLE, DIMENSION(:)   :: DELLA0
REAL, ALLOCATABLE, DIMENSION(:)   :: DRGM
REAL, ALLOCATABLE, DIMENSION(:)   :: WOK1
REAL, ALLOCATABLE, DIMENSION(:)   :: WOK2

DELPRO = REAL(IDELPRO)
DELPH0 = 0.5*DELPRO/DELPHI
DELTH0 = 0.5*DELPRO/DELTR
DELTHR = 0.5*DELPRO/DELTH
DELFR0 = 0.5*DELPRO/(0.1*ZPI)

ALLOCATE(DELLA0(NIJS:NIJL))
ALLOCATE(DRGM(NIJS:NIJL))
DO IJ = NIJS,NIJL
   DELLA0(IJ) = 0.5*DELPRO*DCO(IJ)/DELLAM(KFROMIJ(IJ))
   DRGM(IJ)   = SINPH(KFROMIJ(IJ))*DCO(IJ)
END DO

ALLOCATE(WOK1(ninf-1:nsup))
ALLOCATE(WOK2(ninf-1:nsup))

IF (.NOT.ALLOCATED(COEF_LAT))   ALLOCATE(COEF_LAT(NIJS:NIJL,1:KL,1:ML,1:2,1:1))
IF (.NOT.ALLOCATED(COEF_LON))   ALLOCATE(COEF_LON(NIJS:NIJL,1:KL,ML,1:2))
IF (.NOT.ALLOCATED(COEF_THETA)) ALLOCATE(COEF_THETA(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SIGD))  ALLOCATE(COEF_SIGD(NIJS:NIJL,1:KL,1:ML,1:2))
IF (.NOT.ALLOCATED(COEF_SUM))   ALLOCATE(COEF_SUM(NIJS:NIJL,1:KL,1:ML))

FRE: DO M = 1,ML
   MM = MPM(M,-1)
   MP = MPM(M, 1)
   DFP = DELFR0/FR(M)
   DFM = DELFR0/FR(MM)

   DIR: DO K = 1,KL
      SD = SINTH(K)
      CD = COSTH(K)
      KP = KPM(K, 1)
      KM = KPM(K,-1)
      SP = DELTH0*(SINTH(K)+SINTH(KP))      !! GRID REFRACTION.
      SM = DELTH0*(SINTH(K)+SINTH(KM))

      WOK1(ninf-1:nsup) = U(ninf-1:nsup)+SD*CGOND(ninf-1:nsup,M)
      WOK2(ninf-1:nsup) = (V(ninf-1:nsup)+CD*CGOND(ninf-1:nsup,M))*DELPH0

      DO IJ = NIJS,NIJL

         DLW = ( WOK1(IJ) + WOK1(KLON(IJ,1)) ) * DELLA0(IJ)
         DLE = ( WOK1(IJ) + WOK1(KLON(IJ,2)) ) * DELLA0(IJ)
         DTC = 1. - MAX(0.,DLE) + MIN(0.,DLW)
         COEF_LON(IJ,K,M,1) =  MAX(0.,DLW)
         COEF_LON(IJ,K,M,2) = -MIN(0.,DLE)

         ic = 1
         DPS = wok2(IJ) + DPSN(IJ,IC)*wok2(KLAT(IJ,IC,1))
         ic = 2
         DPN = wok2(IJ) + DPSN(IJ,IC)*wok2(KLAT(IJ,IC,1))

         DTC = DTC - MAX(0.,DPN) + MIN(0.,DPS)
         COEF_LAT(IJ,K,M,1,1)  =  MAX(0.,DPS)
         COEF_LAT(IJ,K,M,2,1)  = -MIN(0.,DPN)

         DTP = SP*DRGM(IJ)*CGOND(IJ,M) + (THDD(IJ,K ,M)+THDD(IJ,KP,M))*DELTHR
         DTM = SM*DRGM(IJ)*CGOND(IJ,M) + (THDD(IJ,K ,M)+THDD(IJ,KM,M))*DELTHR
         DTC = DTC - MAX(0.,DTP) + MIN(0.,DTM)
         COEF_THETA(IJ,K,M,1) = -MIN(0.,DTP)
         COEF_THETA(IJ,K,M,2) =  MAX(0.,DTM)

         DOP = (SIDC(IJ,K,M) + SIDC(IJ,K,MP))*DFP
         DOM = (SIDC(IJ,K,M) + SIDC(IJ,K,MM))*DFM
         DTC = DTC - MAX(0.,DOP) + MIN(0.,DOM)
         COEF_SIGD(IJ,K,M,1) = -MIN(0.,DOP)/1.1
         COEF_SIGD(IJ,K,M,2) =  MAX(0.,DOM)*1.1

         COEF_SUM(IJ,K,M) = DTC

      END DO
   END DO DIR
END DO FRE

IF (L_OBSTRUCTION) THEN
   DO M = 1,ML
      DO K = 1,KL
         DO IJ = NIJS,NIJL
            COEF_LAT(IJ,K,M,1,1) = COEF_LAT(IJ,K,M,1,1)*OBSLAT(IJ,1,M)
            COEF_LAT(IJ,K,M,2,1) = COEF_LAT(IJ,K,M,2,1)*OBSLAT(IJ,2,M)
            COEF_LON(IJ,K,M,1)   = COEF_LON(IJ,K,M,1)  *OBSLON(IJ,1,M)
            COEF_LON(IJ,K,M,2)   = COEF_LON(IJ,K,M,2)  *OBSLON(IJ,2,M)
         END DO
      END DO
   END DO
END IF

DEALLOCATE(DELLA0)
DEALLOCATE(DRGM)
DEALLOCATE(WOK1)
DEALLOCATE(WOK2)

END SUBROUTINE PREP_SPHER_SHALLOW_CURR_REG

END SUBROUTINE PREPARE_PROPAGATION_COEF

! ---------------------------------------------------------------------------- !

SUBROUTINE CHECK_CFL

! ---------------------------------------------------------------------------- !
!                                                                              !
!     CHECK_CFL  - CFL CHECK.                                                  !
!                                                                              !
!     H. GUNTHER   GKSS   FEBRUARY 2002                                        !
!     H. GUNTHER   HZG    MAY 2017                                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO CHECK THE PROPAGATION WEIGHTS AND MODIFY PROPAGATION TIME STEP,     !
!       IF NECESSARY                                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       FOR EACH GRID POINT, FREQUENCY AND DIRECTION THE RELATIVE LOSS OF      !
!       ENERGY PER TIME STEP IS COMPUTED AND CHECKED TO BE LESS THAN 1.        !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL    :: COEF_LAT_1_1_MAX, COEF_LAT_1_1_MIN
REAL    :: COEF_LAT_1_2_MAX, COEF_LAT_1_2_MIN
REAL    :: COEF_LON_1_1_MAX, COEF_LON_1_1_MIN
REAL    :: COEF_LON_1_2_MAX, COEF_LON_1_2_MIN
REAL    :: COEF_THETA_1_MAX, COEF_THETA_1_MIN
REAL    :: COEF_THETA_2_MAX, COEF_THETA_2_MIN
REAL    :: COEF_SIGD_1_MAX,  COEF_SIGD_1_MIN
REAL    :: COEF_SIGD_2_MAX,  COEF_SIGD_2_MIN
REAL    :: COEF_SUM_MAX,     COEF_SUM_MIN
REAL    :: COEF_MAX,         COEF_MIN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.0 MAXIMA AND MINIMA OF ALL WEIGTHS.                                    !
!         ---------------------------------                                    !

COEF_LAT_1_1_MAX = MAXVAL(COEF_LAT(:,:,:,:,1))
COEF_LAT_1_1_MIN = MINVAL(COEF_LAT(:,:,:,:,1))
IF (REDUCED_GRID) THEN
  COEF_LAT_1_2_MAX = MAXVAL(COEF_LAT(:,:,:,:,2))
  COEF_LAT_1_2_MIN = MINVAL(COEF_LAT(:,:,:,:,2))
ELSE
   COEF_LAT_1_2_MAX = 0.1
   COEF_LAT_1_2_MIN = 0.1
END IF

COEF_LON_1_1_MAX = MAXVAL(COEF_LON(:,:,:,1))
COEF_LON_1_1_MIN = MINVAL(COEF_LON(:,:,:,1))
IF (REFRACTION_C_RUN) THEN
   COEF_LON_1_2_MAX = MAXVAL(COEF_LON(:,:,:,2))
   COEF_LON_1_2_MIN = MINVAL(COEF_LON(:,:,:,2))
ELSE
   COEF_LON_1_2_MAX = 0.1
   COEF_LON_1_2_MIN = 0.1
END IF
IF (ALLOCATED(COEF_THETA)) THEN
   COEF_THETA_1_MAX = MAXVAL(COEF_THETA(:,:,:,1))
   COEF_THETA_1_MIN = MINVAL(COEF_THETA(:,:,:,1))
   COEF_THETA_2_MAX = MAXVAL(COEF_THETA(:,:,:,2))
   COEF_THETA_2_MIN = MINVAL(COEF_THETA(:,:,:,2))
ELSE
   COEF_THETA_1_MAX = 0.1
   COEF_THETA_1_MIN = 0.1
   COEF_THETA_2_MAX = 0.1
   COEF_THETA_2_MIN = 0.1
END IF
IF (ALLOCATED(COEF_SIGD)) THEN
   COEF_SIGD_1_MAX = MAXVAL(COEF_SIGD(:,:,:,1))
   COEF_SIGD_1_MIN = MINVAL(COEF_SIGD(:,:,:,1))
   COEF_SIGD_2_MAX = MAXVAL(COEF_SIGD(:,:,:,2))
   COEF_SIGD_2_MIN = MINVAL(COEF_SIGD(:,:,:,2))
ELSE
   COEF_SIGD_1_MAX = 0.1
   COEF_SIGD_1_MIN = 0.1
   COEF_SIGD_2_MAX = 0.1
   COEF_SIGD_2_MIN = 0.1
END IF

COEF_SUM_MAX     = 1.-MINVAL(COEF_SUM(:,:,:))
COEF_SUM_MIN     = 1.-MAXVAL(COEF_SUM(:,:,:))

COEF_MAX = MAX(COEF_LAT_1_1_MAX, COEF_LAT_1_2_MAX,                            &
&              COEF_LON_1_1_MAX, COEF_LON_1_2_MAX,                            &
&              COEF_THETA_1_MAX, COEF_THETA_2_MAX,                            &
&              COEF_SIGD_1_MAX,  COEF_SIGD_2_MAX, COEF_SUM_MAX)
COEF_MIN = MIN(COEF_LAT_1_1_MIN, COEF_LAT_1_2_MIN,                            &
&              COEF_LON_1_2_MIN, COEF_LON_1_2_MIN,                            &
&              COEF_THETA_1_MIN, COEF_THETA_2_MIN,                            &
&              COEF_SIGD_1_MIN,  COEF_SIGD_2_MIN, COEF_SUM_MIN)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.0 CHECK FOR NEGATIVE WEIGHTS.                                          !
!         ---------------------------                                          !

ALLOCATE (CFLMAX(1:petotal))

CFLMAX(irank) = COEF_MIN
if (petotal/=1) then
   call mpi_gather_cfl (CFLMAX)
   total_cfl = minval(CFLMAX)
else
   total_cfl = CFLMAX(irank)
end if

IF (total_cfl.LT.0.) THEN
   WRITE(IU06,*) ' **********************************************************'
   WRITE(IU06,*) ' *                                                        *'
   WRITE(IU06,*) ' *           FATAL ERROR IN SUB. CFLCHECK                 *'
   WRITE(IU06,*) ' *           ================================             *'
   WRITE(IU06,*) ' *                                                        *'
   WRITE(IU06,*) ' *            VIOLATIONS OF CFL-CRITERION                 *'
   WRITE(IU06,*) ' *                                                        *'
   WRITE(IU06,*) ' * negative propagation coefficent detected.              *'
   WRITE(IU06,*) ' *                                                        *'
   WRITE(IU06,*) ' * TOTAl Minimum over all processors = ', total_cfl
   WRITE(IU06,*) ' * coefficents on Process Number     = ', irank
   WRITE(IU06,*) ' * COEF_LAT_1_1_MIN = ', COEF_LAT_1_1_MIN
   IF (REDUCED_GRID) THEN
      WRITE(IU06,*) ' * COEF_LAT_1_2_MIN = ', COEF_LAT_1_2_MIN
   END IF
   WRITE(IU06,*) ' * COEF_LON_1_1_MIN = ', COEF_LON_1_1_MIN
   IF (REDUCED_GRID) THEN
      WRITE(IU06,*) ' * COEF_LON_1_2_MIN = ', COEF_LON_1_2_MIN
   END IF
   IF (ALLOCATED(COEF_THETA)) THEN
      WRITE(IU06,*) ' * COEF_THETA_1_MIN = ', COEF_THETA_1_MIN
      WRITE(IU06,*) ' * COEF_THETA_2_MIN = ', COEF_THETA_2_MIN
   END IF
   IF (ALLOCATED(COEF_SIGD)) THEN
      WRITE(IU06,*) ' * COEF_SIGD_1_MIN = ', COEF_SIGD_1_MIN
      WRITE(IU06,*) ' * COEF_SIGD_2_MIN = ', COEF_SIGD_2_MIN
   END IF
   WRITE(IU06,*) ' * COEF_SUM_MIN     = ', COEF_SUM_MIN
   WRITE(IU06,*) ' *                                                        *'
   WRITE(IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS               *'
   WRITE(IU06,*) ' *                                                        *'
   WRITE(IU06,*) ' **********************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.0 EVALUATE CFL.                                                        !
!         -------------                                                        !

CFLMAX(irank) = COEF_MAX
if (petotal/=1) then
   call mpi_gather_cfl (CFLMAX)
   total_cfl = maxval(CFLMAX)
else
   total_cfl = CFLMAX(irank)
end if

IF (total_cfl.GT.XLIMIT) THEN
   WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE(IU06,*) ' +                                                   +'
   WRITE(IU06,*) ' +           WARNING ERROR IN SUB. CFLCHECK          +'
   WRITE(IU06,*) ' +           ==============================          +'
   WRITE(IU06,*) ' +                                                   +'
   WRITE(IU06,*) ' +            VIOLATIONS OF CFL-CRITERION            +'
   WRITE(IU06,*) ' +                                                   +'
   WRITE(IU06,*) ' + TOTAl Maximum over all processors = ', total_cfl
   WRITE(IU06,*) ' + coefficents on process number     = ', irank
   WRITE(IU06,*) ' + COEF_LAT_1_1_MAX = ', COEF_LAT_1_1_MAX
   IF (REDUCED_GRID) THEN
      WRITE(IU06,*) ' + COEF_LAT_1_2_MAX = ', COEF_LAT_1_2_MAX
   END IF
   WRITE(IU06,*) ' + COEF_LON_1_1_MAX = ', COEF_LON_1_1_MAX
   IF (REDUCED_GRID) THEN
      WRITE(IU06,*) ' * COEF_LON_1_2_MAX = ', COEF_LON_1_2_MAX
   END IF
   IF (ALLOCATED(COEF_THETA)) THEN
      WRITE(IU06,*) ' + COEF_THETA_1_MAX = ', COEF_THETA_1_MAX
      WRITE(IU06,*) ' + COEF_THETA_2_MAX = ', COEF_THETA_2_MAX
   END IF
   IF (ALLOCATED(COEF_SIGD)) THEN
      WRITE(IU06,*) ' * COEF_SIGD_1_MAX = ', COEF_SIGD_1_MAX
      WRITE(IU06,*) ' * COEF_SIGD_2_MAX = ', COEF_SIGD_2_MAX
   END IF
   WRITE(IU06,*) ' + COEF_SUM_MAX     = ', COEF_SUM_MAX
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
   COUNTER = 1
END IF

IF (COUNTER.NE.1) THEN
   COEF_LAT(:,:,:,:,:) = COEF_LAT(:,:,:,:,:) /REAL(COUNTER)
   COEF_LON(:,:,:,:)   = COEF_LON(:,:,:,:)   /REAL(COUNTER)
   IF (ALLOCATED(COEF_THETA)) COEF_THETA(:,:,:,:) = COEF_THETA(:,:,:,:) /REAL(COUNTER)
   IF (ALLOCATED(COEF_SIGD)) THEN
       COEF_SIGD(:,:,:,:) = COEF_SIGD(:,:,:,:) /REAL(COUNTER)
   ENDIF
   COEF_SUM(:,:,:)     = 1. - (1.-COEF_SUM(:,:,:))/REAL(COUNTER)
END IF

DEALLOCATE (CFLMAX)

END SUBROUTINE CHECK_CFL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_PROPAGATION_MODULE

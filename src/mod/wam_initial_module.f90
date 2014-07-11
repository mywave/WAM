MODULE WAM_INITIAL_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS ALL DATA WHICH ARE USED BY THE PRESET PROGRAM.        !
!   ALL PROCEDURES ARE INCLUDED TO COMPUTE THE WAM MODEL COLDSTART FILE,       !
!   TO COMPUTE SAVE AND CONNECT RESTART FILES TO THE WAM MODEL.                !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE,    ONLY:  &
&       ABORT1,                   &  !! TERMINATES PROCESSING.
&       INCDATE                      !! UPDATES DATE/TIME GROUP.

USE WAM_RESTART_MODULE,    ONLY:  &
&       CONNECT_RESTART              !! CONNECT RESTART FILE TO THE WAM MODEL.

USE WAM_COLDSTART_MODULE,  ONLY:  &
&       PREPARE_COLDSTART            !! PREPARES COLDSTART FIELDS FOR WAMODEL.

USE WAM_ICE_MODULE,        ONLY:  &
&       PUT_ICE,                  &  !! PUTS ICE INDICATOR INTO DATA FILED.
&       GET_ICE                      !! GETS A NEW ICE FIELD.

USE WAM_TOPO_MODULE,       ONLY:  &
&       WAM_TOPO,                 &  !! READ AND PREPARE DEPTHS.
&       PUT_DRY,                  &  !! PUTS DRY INDICATOR INTO DATA FILED.
&       FIND_DRY_POINTS              !! FINDS DRY POINTS.

USE WAM_CURRENT_MODULE,    ONLY:  &
&       WAM_CURRENT                  !! READ AND PREPARE CURRENTS.

use wam_special_module,    only:  &
&       chready                      !! wait for wind/ice files 
    
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE, ONLY: DEG

USE WAM_FRE_DIR_MODULE, ONLY: KL, ML, FR, CO, TH, DELTH, DELTR, COSTH, SINTH,  &
&                             GOM, C, INV_LOG_CO,                              &
&                             DF, DF_FR2, DF_FR,                               &
&                             DFIM, DFIMOFR, DFIM_FR2, DFIM_FR, FR5, FRM5,     &
&                             FMIN, MO_TAIL, MM1_TAIL, MP1_TAIL, MP2_TAIL,     &
&                             NDEPTH, TCGOND, TSIHKD, TFAK, TFAC_ST

USE WAM_GRID_MODULE,    ONLY: HEADER, NX, NY, NSEA, NLON_RG,                   &
&                             XDELLA, DELLAM, XDELLO, ZDELLO, DELPHI,          &
&                             AMOWEP, AMOSOP, AMOEAP, AMONOP, IPER,            &
&                             SINPH, COSPH, DEPTH_B, KLAT, KLON, IXLG, KXLT,   &
&                             L_S_MASK, ONE_POINT, REDUCED_GRID

USE WAM_NEST_MODULE,    ONLY: N_NEST, MAX_NEST, N_NAME, n_code,                &
                              NBOUNC, IJARC, BLATC, BLNGC, DLAMAC, DPHIAC,     &
&                             N_SOUTH, N_NORTH, N_EAST, N_WEST, N_ZDEL,        &
&                             NBINP, NBOUNF, C_NAME, BLNGF, BLATF, IJARF,      &
&                             IBFL, IBFR, BFW

USE WAM_MODEL_MODULE,   ONLY: FL3, U10, UDIR, TAUW, USTAR, Z0, DEPTH, INDEP,   &
&                             U, V

USE WAM_TIMOPT_MODULE,  ONLY: CDATEA, CDATEE, CDTPRO, IDELPRO, IDELT, IDEL_WAM,&
&                             COLDSTART, SPHERICAL_RUN, SHALLOW_RUN,           &
&                             REFRACTION_C_RUN,                                &
&                             CDATEWO, ifcst,                                  &
&                             CDTA, TOPO_RUN, CD_TOPO_NEW,                     &
&                             CDCA, CURRENT_RUN, CD_CURR_NEW 

USE WAM_FILE_MODULE,    ONLY: FILE03, IU06, ITEST, IU07, FILE07, FILE08, FILE09

USE WAM_WIND_MODULE,    ONLY: IDELWI, IDELWO

USE WAM_TOPO_MODULE,    ONLY: IDELTI, IDELTO

USE WAM_CURRENT_MODULE, ONLY: IDELCI, IDELCO

USE WAM_PROPAGATION_MODULE, ONLY: NADV, DCO, DPSN, MP, MM, KP, KM

USE WAM_ICE_MODULE,         ONLY: ICE_RUN

use wam_mpi_module,     only: nijs, nijl, ninf, nsup
use wam_special_module, only: readyf

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
PRIVATE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE PREPARE_START               !! PREPARES START FIELDS FOR WAMODEL.
   MODULE PROCEDURE PREPARE_START
END INTERFACE
PUBLIC PREPARE_START

INTERFACE READ_PREPROC_FILE           !! READS PREPROC OUTPUT FILE.
   MODULE PROCEDURE READ_PREPROC_FILE
END INTERFACE
PUBLIC READ_PREPROC_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE PREPARE_FIRST_DEPTH             !! PREPARES START FIELDS FOR WAMODEL.
   MODULE PROCEDURE PREPARE_FIRST_DEPTH
END INTERFACE
PRIVATE PREPARE_FIRST_DEPTH

INTERFACE PREPARE_FIRST_CURRENT             !! PREPARES START FIELDS FOR WAMODEL.
   MODULE PROCEDURE PREPARE_FIRST_CURRENT
END INTERFACE
PRIVATE PREPARE_FIRST_CURRENT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_START

INTEGER :: M, K, LEN
LOGICAL :: ERROR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. ALLOCATE ARRAYS.                                                      !
!        ----------------                                                      !

IF (ALLOCATED(U10   )) DEALLOCATE(U10  )
IF (ALLOCATED(UDIR  )) DEALLOCATE(UDIR )
IF (ALLOCATED(USTAR )) DEALLOCATE(USTAR)
IF (ALLOCATED(Z0    )) DEALLOCATE(Z0   )
IF (ALLOCATED(TAUW  )) DEALLOCATE(TAUW )
IF (ALLOCATED(DEPTH )) DEALLOCATE(DEPTH)
IF (ALLOCATED(INDEP )) DEALLOCATE(INDEP)
IF (ALLOCATED(FL3   )) DEALLOCATE(FL3  )
IF (ALLOCATED(U     )) DEALLOCATE(U    )
IF (ALLOCATED(V     )) DEALLOCATE(V    )
     
ALLOCATE (FL3(nijs:nijl, 1:KL,1:ML))
     
ALLOCATE (U10(nijs:nijl))
allocate (udir(nijs:nijl))
allocate (ustar(nijs:nijl))
allocate (z0(nijs:nijl))
allocate (tauw(nijs:nijl))

ALLOCATE(DEPTH(1:NSEA))
ALLOCATE(INDEP(1:NSEA))
ALLOCATE(U    (1:NSEA))
ALLOCATE(V    (1:NSEA))

USTAR = 0.
TAUW  = 0.
Z0    = 0.
U     = 0.
V     = 0.
DEPTH = 999.0

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. NEIGHBOUR INDEX FOR FREQUENCY AND DIRECTION.                          !
!        --------------------------------------------                          !

IF (.NOT.ALLOCATED(MP)) THEN
   ALLOCATE(MP(1:ML))
   MP = (/(MIN(M,ML),M=2,ML+1)/)
END IF
IF (.NOT.ALLOCATED(MM)) THEN
   ALLOCATE(MM(1:ML))
   MM = (/(MAX(M,1),M=0,ML-1)/)
END IF
IF (.NOT.ALLOCATED(KP)) THEN
   ALLOCATE(KP(1:KL))
   KP = CSHIFT((/(K,K=1,KL)/),1)
END IF
IF (.NOT.ALLOCATED(KM)) THEN
   ALLOCATE(KM(1:KL))
   KM = CSHIFT(KP,-2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COSINE OF LATITUDE FACTORS (IF SPHERICAL GRID).                       !
!        -----------------------------------------------                       !

IF (SPHERICAL_RUN .AND. .NOT.ONE_POINT) THEN
   IF (.NOT.ALLOCATED(DCO))  ALLOCATE(DCO(NINF:NSUP))
   IF (.NOT.ALLOCATED(DPSN)) ALLOCATE(DPSN(nijs:nijl,2))

   DCO(ninf:nsup) = 1./COSPH(KXLT(ninf:nsup))           !! COSINE OF LATITUDE.
   DPSN = 1.
   WHERE (KLAT(nijs:nijl,1).NE.ninf-1) DPSN(nijs:nijl,1) =                     &
&         DCO(nijs:nijl)/DCO(KLAT(nijs:nijl,1))         !! COS PHI FACTOR SOUTH.
   WHERE (KLAT(nijs:nijl,2).NE.ninf-1) DPSN(nijs:nijl,2) =                     &
&         DCO(nijs:nijl)/DCO(KLAT(nijs:nijl,2))         !! COS PHI FACTOR NORTH.
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. GENERATE START FIELDS OR READ RESTART FILE.                           !
!        -------------------------------------------                           !

IF (COLDSTART) THEN
   CALL PREPARE_COLDSTART
   IF (ITEST.GE.2) WRITE(IU06,*) '    SUB. PREPARE_START: PREPARE_COLDSTART DONE'
ELSE
   CALL CONNECT_RESTART
   IF (ITEST.GE.2) WRITE(IU06,*) '    SUB. PREPARE_START: CONNECT_RESTART DONE'
END IF

! ---------------------------------------------------------------------------- !
!
!     5. ICE INFORMATION.
!        ----------------

ICE_RUN = .FALSE.
LEN = LEN_TRIM(FILE03)
IF (LEN.GT.0) THEN
   INQUIRE (FILE=FILE03(1:LEN), EXIST=ICE_RUN)
   IF (.NOT.ICE_RUN) THEN
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +        WARNING ERROR SUB. PREPARE_START.         +'
      WRITE (IU06,*) ' +        =================================         +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' + THE ICE INPUT FILE DEFINED IN THE USER INPUT     +'
      WRITE (IU06,*) ' + DOES NOT EXIST.                                  +'
      WRITE (IU06,*) ' + FILENAME IS FIL03 = ', TRIM(FILE03)
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +               MODEL CONTINUES                    +'
      WRITE (IU06,*) ' +                 WITHOUT ICE                      +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
   END IF
END IF

if (readyf) then
   call chready (ifcst)  !! wait for first ice/wind-file
   write (iu06,*) ' +++ WAM waits for wind ready files'
else
   write (iu06,*) ' +++ WAM runs without wind ready files'
endif
   
IF (ICE_RUN) THEN
   CALL GET_ICE
   IF (ITEST.GE.2) WRITE(IU06,*) '    SUB. PREPARE_START: GET_ICE DONE'

   CALL PUT_ICE (FL3, 0.)
   IF (ITEST.GE.2) WRITE(IU06,*) '    SUB. PREPARE_START: ICE INSERTED'
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. PREPARE FIRST DEPTH FIELD AND DATE FOR NEXT DEPTH FIELD.              !
!        --------------------------------------------------------              !

CALL PREPARE_FIRST_DEPTH
IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '    SUB. PREPARE_START: PREPARE_FIRST_DEPTH DONE '
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. PREPARE FIRST CURRENT FIELD AND DATE FOR NEXT CURRENT FIELD.          !
!        ------------------------------------------------------------          !

CALL PREPARE_FIRST_CURRENT
IF (ITEST.GE.2) THEN
   WRITE (IU06,*) '    SUB. PREPARE_START: PREPARE_FIRST_CURRENT DONE '
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     8. NUMBER OF PROPAGATION TIME STEPS PER WAM MODEL CALL.                  !
!        ----------------------------------------------------                  !

IF (ONE_POINT) THEN
   IDELPRO = IDELT
END IF

NADV = MAX(IDELWI, IDELCI, IDELTI, IDELPRO, IDELT)

ERROR = (MOD(NADV,IDELPRO).NE.0) .OR. (MOD(NADV,IDELT ).NE.0)
IF (IDELWI.GT.0)   ERROR = ERROR .OR. (MOD(NADV,IDELWI).NE.0)
IF (IDELTI.GT.0)   ERROR = ERROR .OR. (MOD(NADV,IDELTI).NE.0)
IF (IDELCI.GT.0)   ERROR = ERROR .OR. (MOD(NADV,IDELCI).NE.0)


IF (ERROR) THEN
   WRITE(IU06,*) '*************************************************'
   WRITE(IU06,*) '*                                               *'
   WRITE(IU06,*) '*     FATAL ERROR IN SUB. PREPARE_START         *'
   WRITE(IU06,*) '*     =================================         *'
   WRITE(IU06,*) '*                                               *'
   WRITE(IU06,*) '* TIMESTEPS DO NOT SYNCRONIZE AT THERE MAXIMUM. *'
   WRITE(IU06,*) '* WIND INPUT TIMESTEP      : ', IDELWI
   WRITE(IU06,*) '* DEPTH INPUT TIMESTEP     : ', IDELTI
   WRITE(IU06,*) '* CURRENT INPUT TIMESTEP   : ', IDELCI
   WRITE(IU06,*) '* SOURCE FUNCTION TIMESTEP : ', IDELT
   WRITE(IU06,*) '* PROPAGATION TIMESTEP     : ', IDELPRO
   WRITE(IU06,*) '*                                               *'
   WRITE(IU06,*) '*        PROGRAM ABORTS  PROGRAM ABORTS         *'
   WRITE(IU06,*) '*                                               *'
   WRITE(IU06,*) '*************************************************'
   CALL ABORT1
END IF

NADV = MAX(NADV/IDELPRO,1)
IDEL_WAM = NADV*IDELPRO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     9. INITIALIZE DATE FOR NEXT WIND FIELD.                                  !
!        ------------------------------------                                  !

CDATEWO = CDTPRO
IF (IDELT.GT.IDELWO) THEN
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +        WARNING ERROR SUB. PREPARE_START.         +'
   WRITE (IU06,*) ' +        =================================         +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + WIND OUTPUT TIMESTEP     : ', IDELWO
   WRITE (IU06,*) ' +    IS LESS THAN                                  +'
   WRITE (IU06,*) ' + SOURCE FUNCTION TIMESTEP : ', IDELT
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +            MODEL CONTINUES WITH                  +'
   WRITE (IU06,*) ' + WIND OUTPUT TIMESTEP = SOURECE FUNCTION TIMESTEP.+'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
   IDELWO = IDELT
END IF
IF (IDELT.LT.IDELWO) CALL INCDATE(CDATEWO,IDELWO/2)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    10. INITIALIZE DATE FOR NEXT DEPTH FIELD.                                 !
!        -------------------------------------                                 !

IF (IDELTI.LE.0 .OR. .NOT.TOPO_RUN) THEN
   CD_TOPO_NEW = '99991231235900'
   TOPO_RUN = .FALSE.
   IDELTI = 0
   IDELTO = 0
ELSE
   CD_TOPO_NEW = CDTPRO
   IF (IDELPRO.GT.IDELTO) THEN
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +        WARNING ERROR SUB. PREPARE_START.         +'
      WRITE (IU06,*) ' +        =================================         +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' + TOPO OUTPUT TIMESTEP  : ', IDELTO
      WRITE (IU06,*) ' +    IS LESS THAN                                  +'
      WRITE (IU06,*) ' + PROPAGATION TIMESTEP  : ', IDELPRO
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +            MODEL CONTINUES WITH                  +'
      WRITE (IU06,*) ' +   TOPO OUTPUT TIMESTEP = PROPAGATION TIMESTEP.   +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      IDELTO = IDELPRO
   END IF
   IF (IDELPRO.LT.IDELTO) CALL INCDATE(CD_TOPO_NEW,IDELTO/2)
 END IF
IF (ITEST.GE.3) THEN
   WRITE (IU06,*) '        NEXT DEPTH DATE IS  CD_TOTO_NEW = ', CD_TOPO_NEW
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!    11. INITIALIZE DATE FOR NEXT CURRENT FIELD.                               !
!        ---------------------------------------                               !

IF (IDELCI.LE.0 .OR. .NOT.CURRENT_RUN) THEN
   CD_CURR_NEW = '99991231235900'
   CURRENT_RUN = .FALSE.
   IDELCI = 0
   IDELCO = 0
ELSE
   CD_CURR_NEW = CDTPRO
   IF (IDELPRO.GT.IDELCO) THEN
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +        WARNING ERROR SUB. PREPARE_START.         +'
      WRITE (IU06,*) ' +        =================================         +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' + CURRENT OUTPUT TIMESTEP  : ', IDELCO
      WRITE (IU06,*) ' +    IS LESS THAN                                  +'
      WRITE (IU06,*) ' + PROPAGATION TIMESTEP     : ', IDELPRO
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +            MODEL CONTINUES WITH                  +'
      WRITE (IU06,*) ' + CURRENT OUTPUT TIMESTEP = PROPAGATION TIMESTEP.  +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      IDELCO = IDELPRO
   END IF
   IF (IDELPRO.LT.IDELCO) CALL INCDATE(CD_CURR_NEW ,IDELCO/2)
END IF
IF (ITEST.GE.3) THEN
   WRITE (IU06,*) '        NEXT CURRENT DATE IS CD_CURR_NEW = ', CD_CURR_NEW
END IF

END SUBROUTINE PREPARE_START

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE READ_PREPROC_FILE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READ_PREPROC_FILE -  READ OUTPUT FILE FROM PREPROC.                        !
!                                                                              !
!     H. GUNTHER      GKSS/ECMWF      MAY 1990                                 !
!     H. GUNTHER      GKSS        OCTOBER 2000  FT90                           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       INPUT OF PREPROC OUTPUT.                                               !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       UNFORMATTED READ FROM FILE07.                                          !
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

INTEGER  :: IOS, LEN, I

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. OPEN GRID_INFO FILE FROM PREPROC OUTPUT.                              !
!        ----------------------------------------                              !

IOS = 0
LEN = LEN_TRIM(FILE07)
OPEN (UNIT=IU07, FILE=FILE07(1:LEN), FORM='UNFORMATTED', STATUS='OLD',         &
&                                                                 IOSTAT=IOS)
IF (IOS.NE.0) THEN
   WRITE (IU06,*) ' ****************************************************'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' *     FATAL ERROR IN SUB. READ_PREPROC_FILE        *'
   WRITE (IU06,*) ' *     =====================================        *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' * PREPROC OUTPUT FILE COULD NOT BE OPENED          *'
   WRITE (IU06,*) ' *    ERROR CODE IS IOSTAT = ', IOS
   WRITE (IU06,*) ' *    FILE NAME IS  FILE07 = ', FILE07(1:LEN)
   WRITE (IU06,*) ' *    UNIT IS         IU07 = ', IU07
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' *         PROGRAM ABORTS  PROGRAM ABORTS           *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' ****************************************************'
   CALL ABORT1
END IF
READ (IU07) HEADER

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. READ COARSE GRID BOUNDARY OUTPUT INFORMATION.                         !
!        ---------------------------------------------                         !

READ(IU07) N_NEST, MAX_NEST
IF (.NOT.ALLOCATED(NBOUNC )) ALLOCATE (NBOUNC(N_NEST))
IF (.NOT.ALLOCATED(N_NAME )) ALLOCATE (N_NAME(N_NEST))
if (.not.allocated(n_code )) allocate (n_code(n_nest))
IF (.NOT.ALLOCATED(IJARC  )) ALLOCATE (IJARC(MAX_NEST,N_NEST))
IF (.NOT.ALLOCATED(N_SOUTH)) ALLOCATE (N_SOUTH(N_NEST))
IF (.NOT.ALLOCATED(N_NORTH)) ALLOCATE (N_NORTH(N_NEST))
IF (.NOT.ALLOCATED(N_EAST )) ALLOCATE (N_EAST(N_NEST))
IF (.NOT.ALLOCATED(N_WEST )) ALLOCATE (N_WEST(N_NEST))
IF (.NOT.ALLOCATED(BLNGC  )) ALLOCATE (BLNGC(MAX_NEST,N_NEST))
IF (.NOT.ALLOCATED(BLATC  )) ALLOCATE (BLATC(MAX_NEST,N_NEST))
IF (.NOT.ALLOCATED(N_ZDEL )) ALLOCATE (N_ZDEL(MAX_NEST,N_NEST))
DO I=1,N_NEST
   READ(IU07) NBOUNC(I), N_NAME(I), n_code(i)
   IF (NBOUNC(I).GT.0) THEN
      READ(IU07) IJARC(1:NBOUNC(I),I)
      READ(IU07) XDELLO, XDELLA, N_SOUTH(I), N_NORTH(I), N_EAST(I), N_WEST(I), &
&                BLNGC(1:NBOUNC(I),I), BLATC(1:NBOUNC(I),I), N_ZDEL(1:NBOUNC(I),I)
   END IF
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. READ FINE GRID BOUNDARY INPUT INFORMATION.                            !
!        ------------------------------------------                            !

READ (UNIT=IU07) NBOUNF, NBINP, C_NAME

IF (NBOUNF.GT.0) THEN
   IF (.NOT.ALLOCATED(BLNGF)) ALLOCATE (BLNGF(NBOUNF))
   IF (.NOT.ALLOCATED(BLATF)) ALLOCATE (BLATF(NBOUNF))
   IF (.NOT.ALLOCATED(IJARF)) ALLOCATE (IJARF(NBOUNF))
   IF (.NOT.ALLOCATED(IBFL )) ALLOCATE (IBFL(NBOUNF))
   IF (.NOT.ALLOCATED(IBFR )) ALLOCATE (IBFR(NBOUNF))
   IF (.NOT.ALLOCATED(BFW  )) ALLOCATE (BFW(NBOUNF))
   READ (IU07) BLNGF(1:NBOUNF), BLATF(1:NBOUNF), IJARF(1:NBOUNF),              &
&             IBFL(1:NBOUNF), IBFR(1:NBOUNF), BFW(1:NBOUNF)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. READ FREQUENCY DIRECTION GRID.                                        !
!        ------------------------------                                        !

READ (IU07) ML, KL
IF ( .NOT.ALLOCATED(FR)      ) ALLOCATE( FR(ML)      )
IF ( .NOT.ALLOCATED(DFIM)    ) ALLOCATE( DFIM(ML)    )
IF ( .NOT.ALLOCATED(GOM)     ) ALLOCATE( GOM(ML)     )
IF ( .NOT.ALLOCATED(C)       ) ALLOCATE( C(ML)       )
IF ( .NOT.ALLOCATED(TH)      ) ALLOCATE( TH(KL)      )
IF ( .NOT.ALLOCATED(COSTH)   ) ALLOCATE( COSTH(KL)   )
IF ( .NOT.ALLOCATED(SINTH)   ) ALLOCATE( SINTH(KL)   )
IF ( .NOT.ALLOCATED(DF)      ) ALLOCATE( DF(ML)      )
IF ( .NOT.ALLOCATED(DF_FR)   ) ALLOCATE( DF_FR(ML)   )
IF ( .NOT.ALLOCATED(DF_FR2)  ) ALLOCATE( DF_FR2(ML)  )
IF ( .NOT.ALLOCATED(DFIM)    ) ALLOCATE( DFIM(ML)    )
IF ( .NOT.ALLOCATED(DFIMOFR) ) ALLOCATE( DFIMOFR(ML) )
IF ( .NOT.ALLOCATED(DFIM_FR) ) ALLOCATE( DFIM_FR(ML) )
IF ( .NOT.ALLOCATED(DFIM_FR2)) ALLOCATE( DFIM_FR2(ML))
IF ( .NOT.ALLOCATED(FR5)     ) ALLOCATE( FR5(ML)     )
IF ( .NOT.ALLOCATED(FRM5)    ) ALLOCATE( FRM5(ML)    )
READ (IU07)  FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH, INV_LOG_CO,     &
&            DF, DF_FR, DF_FR2, DFIM, DFIMOFR, DFIM_FR, DFIM_FR2, FR5, FRM5,   &
&            FMIN, MO_TAIL, MM1_TAIL, MP1_TAIL, MP2_TAIL

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. READ GRID INFORMATION.                                                !
!        ----------------------                                                !

READ (IU07) NX, NY, NSEA, IPER, ONE_POINT, REDUCED_GRID

IF ( .NOT.ALLOCATED(L_S_MASK)) ALLOCATE( L_S_MASK(1:NX,1:NY) )
IF ( .NOT.ALLOCATED(NLON_RG) ) ALLOCATE( NLON_RG(1:NY) )
IF ( .NOT.ALLOCATED(DELLAM)  ) ALLOCATE( DELLAM(1:NY)  )
IF ( .NOT.ALLOCATED(ZDELLO)  ) ALLOCATE( ZDELLO(1:NY) )
IF ( .NOT.ALLOCATED(IXLG)    ) ALLOCATE( IXLG(1:NSEA) )
IF ( .NOT.ALLOCATED(KXLT)    ) ALLOCATE( KXLT(1:NSEA) )
IF ( .NOT.ALLOCATED(SINPH)   ) ALLOCATE( SINPH(1:NY)  )
IF ( .NOT.ALLOCATED(COSPH)   ) ALLOCATE( COSPH(1:NY)  )
IF ( .NOT.ALLOCATED(KLAT)    ) ALLOCATE( KLAT(1:NSEA,1:2) )
IF ( .NOT.ALLOCATED(KLON)    ) ALLOCATE( KLON(1:NSEA,1:2) )
IF ( .NOT.ALLOCATED(DEPTH_B) ) ALLOCATE( DEPTH_B(1:NSEA) )

READ (IU07) NLON_RG
READ (IU07) DELPHI, DELLAM, SINPH, COSPH, AMOWEP, AMOSOP, AMOEAP, AMONOP,      &
&           XDELLA, XDELLO, ZDELLO
READ (IU07) IXLG, KXLT, L_S_MASK
READ (IU07) KLAT, KLON, DEPTH_B

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. READ SHALLOW WATER TABLES.                                            !
!        --------------------------                                            !

IF ( .NOT.ALLOCATED(TCGOND) ) ALLOCATE( TCGOND(NDEPTH,ML) )
IF ( .NOT.ALLOCATED(TFAK)   ) ALLOCATE( TFAK(NDEPTH,ML)   )
IF ( .NOT.ALLOCATED(TSIHKD) ) ALLOCATE( TSIHKD(NDEPTH,ML) )
IF ( .NOT.ALLOCATED(TFAC_ST)) ALLOCATE( TFAC_ST(NDEPTH,ML))

READ (IU07) TCGOND, TFAK, TSIHKD, TFAC_ST

! ---------------------------------------------------------------------------- !
!                                                                              !
!     9. CLOSE FILE AND RETURN.                                                !
!        ----------------------                                                !

CLOSE (UNIT=IU07, STATUS='KEEP')

END SUBROUTINE READ_PREPROC_FILE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_FIRST_DEPTH

INTEGER  :: LEN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. DEEP WATER RUN: RETURN.                                               !
!        -----------------------                                               !

IF (.NOT.SHALLOW_RUN) THEN
   TOPO_RUN = .FALSE.
   CDTA = ' '
   RETURN 
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. SHALLOW WATER RUN: CHECK TOPO FILE.                                   !
!        -----------------------------------                                   !

LEN = LEN_TRIM(FILE08)
TOPO_RUN = .FALSE.
IF (LEN.GT.0) INQUIRE (FILE=FILE08(1:LEN), EXIST=TOPO_RUN)
IF (LEN.GT.0 .AND. .NOT.TOPO_RUN) THEN
   WRITE (IU06,*) ' ****************************************************'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' *        FATAL ERROR SUB. PREPARE_FIRST_DEPTH.     *'
   WRITE (IU06,*) ' *        =====================================     *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' * THE TOPO INPUT FILE DEFINED IN THE USER INPUT    *'
   WRITE (IU06,*) ' * DOES NOT EXIST.                                  *'
   WRITE (IU06,*) ' *    FILE NAME IS     FILE08 = ', TRIM(FILE08)
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' * CHANGE FILE NAME TO BLANK FOR BASIC DEPTH OR     *'
   WRITE (IU06,*) ' * CORRECT FILE NAME IN THE USER INPUT.             *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' *       PROGRAM ABORTS  PROGRAM ABORTS             *'
   WRITE (IU06,*) ' *                                                  *'
   WRITE (IU06,*) ' ****************************************************'
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. SHALLOW WATER RUN: A TOPO FILE DOES NOT EXIST.                        !
!        ----------------------------------------------                        !

IF (.NOT.TOPO_RUN) THEN
   IF (CDTA.EQ.' ') THEN

!     3.1 DEPTH IS NOT DEFINED: USED BASIC DEPTH FROM PREPROC.                        !

      DEPTH = DEPTH_B
      IF (ITEST.GE.3) THEN
         WRITE (IU06,*) '    SUB.PREPARE_FIRST_DEPTH:'
         WRITE (IU06,*) '        BASIC DEPTH FROM PREPROC'
         WRITE (IU06,*) '        DEPTH IS STATIONARY'
      END IF
   ELSE

!     3.2 DEPTH IS DEFINED FROM RESTART: USED DEPTH FROM RESTART.                        !

      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +      WARNING ERROR SUB. PREPARE_FIRST_DEPTH.     +'
      WRITE (IU06,*) ' +      =======================================     +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' + A RUN WITH STATIONARY DEPTH WAS REQUESTED.       +'
      WRITE (IU06,*) ' + (A TOPO FILE IS NOT DEFINED IN THE USER INPUT)   +'
      WRITE (IU06,*) ' + A DEPTH FIELD WAS FOUND IN THE RESTART FILE.     +'
      WRITE (IU06,*) ' + THIS DEPTH FIELD IS NOT THE BASIC DEPTH FIELD.   +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +               MODEL CONTINUES                    +'
      WRITE (IU06,*) ' +   USING DEPTH FROM RESTART FOR THE FULL RUN.     +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'

      IF (ITEST.GE.3) THEN
         WRITE (IU06,*) '    SUB.PREPARE_FIRST_DEPTH:'
         WRITE (IU06,*) '        DEPTH FROM RESTART FILE'
         WRITE (IU06,*) '        DEPTH IS STATIONARY'
      END IF
   END IF

   CALL FIND_DRY_POINTS
   IF (ITEST.GE.3) THEN
      WRITE (IU06,*) '    SUB. PREPARE_FIRST_DEPTH: FIND_DRY_POINTS DONE'
   END IF
   CALL PUT_DRY (FL3, 0.)
   RETURN

END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. SHALLOW WATER RUN: A TOPO FILE EXISTS.                                !
!        --------------------------------------                                !

IF (TOPO_RUN) THEN
   IF (COLDSTART) THEN

!     4.1 COLD START: PREPARE FIRST DEPTH FIELD.                               !

      CDTA = CDATEA
      CALL WAM_TOPO (DEPTH, CDTA)
      IF (ITEST.GE.3) THEN
         WRITE (IU06,*) '    SUB.PREPARE_FIRST_DEPTH: WAM_TOPO DONE'
         WRITE (IU06,*) '        FIRST DEPTH FIELD PROCESSED'
         WRITE (IU06,*) '        DATE IS ................... CDTA = ', CDTA
      END IF
   ELSE

!     4.2 HOT START.                                                            !

      IF (CDTA.EQ.' ') THEN

!     4.2.1 FIRST DEPTH FIELD IS NOT IN RESTART: PREPARE FIRST DEPTH FIELD.     !

         WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
         WRITE (IU06,*) ' +                                                  +'
         WRITE (IU06,*) ' +      WARNING ERROR SUB. PREPARE_FIRST_DEPTH.     +'
         WRITE (IU06,*) ' +      =======================================     +'
         WRITE (IU06,*) ' +                                                  +'
         WRITE (IU06,*) ' + A HOT START FOR A MODEL RUN WITH DEPTH DIFFERENT +'
         WRITE (IU06,*) ' + FROM THE BASIC DEPTH IN PREPROC IS REQUESTED.    +'
         WRITE (IU06,*) ' + BUT THE FIRST DEPTH FIELD DOES NOT EXIST IN      +'
         WRITE (IU06,*) ' + THE RESTART FILE.                                +'
         WRITE (IU06,*) ' + (PREVIOUS RUN WAS WITH BASIC DEPTH)              +'
         WRITE (IU06,*) ' +                                                  +'
         WRITE (IU06,*) ' +                                                  +'
         WRITE (IU06,*) ' +               MODEL CONTINUES                    +'
         WRITE (IU06,*) ' +     USING DEPTH FIELD FROM TOPO INPUT FILE       +'
         WRITE (IU06,*) ' +                                                  +'
         WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'

         CDTA = CDATEA
         CALL WAM_TOPO (DEPTH, CDTA)
         IF (ITEST.GE.3) THEN
            WRITE (IU06,*) '    SUB.PREPARE_FIRST_DEPTH: WAM_TOPO DONE'
            WRITE (IU06,*) '        FIRST DEPTH FIELD PROCESSED'
            WRITE (IU06,*) '        DATE IS ................... CDTA = ', CDTA
          END IF

      ELSE IF (CDTA.NE.CDATEA) THEN

!     4.2.2 FIRST DEPTH FIELD IS IN RESTART DATES DO NOT MATCH.                 !
 
         WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
         WRITE (IU06,*) ' +                                                  +'
         WRITE (IU06,*) ' +      WARNING ERROR SUB. PREPARE_FIRST_DEPTH.     +'
         WRITE (IU06,*) ' +      =======================================     +'
         WRITE (IU06,*) ' +                                                  +'
         WRITE (IU06,*) ' + A HOT START FOR A MODEL RUN WITH DEPTH DIFFERENT +'
         WRITE (IU06,*) ' + FROM THE BASIC DEPTH IN PREPROC IS REQUESTED.    +'
         WRITE (IU06,*) ' + BUT THE DATE OF THE FIRST DEPTH FIELD IN THE     +'
         WRITE (IU06,*) ' + RESTART FILE IS NOT THE START DATE.              +'
         WRITE (IU06,*) ' + (PREVIOUS RUN WAS WITH STATIONARY DEPTH)         +'
         WRITE (IU06,*) ' + DEPTH DATE FROM RESTART IS CDTA = ', CDTA
         WRITE (IU06,*) ' +                                                  +'
         WRITE (IU06,*) ' +               MODEL CONTINUES                    +'
         WRITE (IU06,*) ' +          CHANGING THE DEPTH DATE                 +'
         WRITE (IU06,*) ' +                                                  +'
         WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
         CDTA = CDATEA
         IF (ITEST.GE.3) THEN
            WRITE (IU06,*) '    SUB.PREPARE_FIRST_DEPTH:'
            WRITE (IU06,*) '        DEPTH FROM RESTART FILE'
            WRITE (IU06,*) '        DATE IS ................... CDTA = ', CDTA
         END IF
      END IF
   END IF      

END IF

CALL FIND_DRY_POINTS
IF (ITEST.GE.3) THEN
   WRITE (IU06,*) '    SUB. PREPARE_FIRST_DEPTH: FIND_DRY_POINTS DONE'
END IF
CALL PUT_DRY (FL3, 0.)

END SUBROUTINE PREPARE_FIRST_DEPTH

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PREPARE_FIRST_CURRENT

INTEGER  :: LEN

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. RUN WITHOUT CURRENT REFRACTION: RETURN.                               !
!        ---------------------------------------                               !

IF (.NOT. REFRACTION_C_RUN) THEN
   CURRENT_RUN = .FALSE.
   CDCA = ' '
   RETURN 
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. RUN WITH CURRENT REFRACTION:CHECK CURRENT FILE.                       !
!        -----------------------------------------------                       !

LEN = LEN_TRIM(FILE09)
CURRENT_RUN = .FALSE.
IF (LEN.GT.0) INQUIRE (FILE=FILE09(1:LEN), EXIST=CURRENT_RUN)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COLD START.                                                           !
!        -----------                                                           !
      
IF (COLDSTART) THEN

   IF (.NOT. CURRENT_RUN) THEN

!     3.1 COLD START AND FILE DOES NOT EXIST: ABORT.                           !

      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *      FATAL ERROR SUB. PREPARE_FIRST_CURRENT.     *'
      WRITE (IU06,*) ' *      =======================================     *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * A CURRENT INPUT FILE IS NOT DEFINED IN THE USER  *'
      WRITE (IU06,*) ' * INPUT OR DOES NOT EXIST. BUT A RUN WITH CURRENT  *'
      WRITE (IU06,*) ' * REFRACTION IS REQUESTED.                         *'
      WRITE (IU06,*) ' *    FILE NAME IS     FILE09 = ', TRIM(FILE09)
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * SWITCH CURRENT REFRACTION OF OR                  *'
      WRITE (IU06,*) ' * CORRECT FILE NAME IN THE USER INPUT.             *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *       PROGRAM ABORTS  PROGRAM ABORTS             *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1
   ELSE

!     3.2 COLD START AND FILE EXISTS: DO FIRST CURRENT FIELD.                  !

      CDCA = CDATEA
      CALL WAM_CURRENT (U, V, CDCA)
      IF (ITEST.GE.3) THEN
         WRITE (IU06,*) '  '
         WRITE (IU06,*) '    SUB. PREPARE_FIRST_CURRENT: WAM_CURRENT DONE'
         WRITE (IU06,*) '        FIRST CURRENT FIELD PROCESSED'
         WRITE (IU06,*) '        DATE IS ................... CDCA = ', CDCA
      END IF
      RETURN
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. HOT START.                                                            !
!        -----------                                                           !

IF (.NOT. CURRENT_RUN) THEN

!     4.1 HOT START AND FILE DOES NOT EXIST.                                   !

   IF (CDCA.EQ.' ') THEN

!     4.1.1 HOT START AND CURRENTS NOT IN RESTART: ABORT.                      !
  
      WRITE (IU06,*) ' ****************************************************'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *      FATAL ERROR SUB. PREPARE_FIRST_CURRENT.     *'
      WRITE (IU06,*) ' *      =======================================     *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' * A HOT START WITH CURRENT REFRACTION IS REQUESTED.*'
      WRITE (IU06,*) ' * A CURRENT INPUT FILE IS NOT DEFINED IN THE USER  *'
      WRITE (IU06,*) ' * INPUT OR DOES NOT EXIST AND A CURRENT FIELD DOES *'
      WRITE (IU06,*) ' * NOT EXIST IN THE RESTART FILE.                   *'
      WRITE (IU06,*) ' *    FILE NAME IS     FILE09 = ', TRIM(FILE09)
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' *       PROGRAM ABORTS  PROGRAM ABORTS             *'
      WRITE (IU06,*) ' *                                                  *'
      WRITE (IU06,*) ' ****************************************************'
      CALL ABORT1   
   ELSE 

!     4.1.2 HOT START AND CURRENTS ARE IN RESTART: CONTINUE.                     !

      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +     WARNING ERROR SUB. PREPARE_FIRST_CURRENT.    +'
      WRITE (IU06,*) ' +     =========================================    +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' + A HOT START WITH CURRENT REFRACTION IS REQUESTED.+'
      WRITE (IU06,*) ' + A CURRENT INPUT FILE IS NOT DEFINED IN THE USER  +'
      WRITE (IU06,*) ' + INPUT OR DOES NOT EXIST.                         +'
      WRITE (IU06,*) ' + A CURRENT FIELD WAS FOUND IN THE RESTART FILE.   +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' +               MODEL CONTINUES                    +'
      WRITE (IU06,*) ' +     USING CURRENT FIELD FROM RESTART FILE.       +'
      WRITE (IU06,*) ' +           CURRENTS ARE STATIONARY                +'
      WRITE (IU06,*) ' +                                                  +'
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      IDELCI = 0
      IF (ITEST.GE.3) THEN
         WRITE (IU06,*) '  '
         WRITE (IU06,*) '    SUB. PREPARE_FIRST_CURRENT:'
         WRITE (IU06,*) '        FIRST CURRENT FIELD FROM RESTART'
         WRITE (IU06,*) '        DATE IS ................... CDCA = ', CDCA
      END IF
      RETURN
   END IF
END IF

!     4.2 HOT START AND FILE EXISTS.                                       !

IF (CDCA.EQ.' ') THEN

!     4.1.1 HOT START, CURRENTS ARE NOT IN RESTART: FIRST FIELD FROM FILE. !

   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +     WARNING ERROR SUB. PREPARE_FIRST_CURRENT.    +'
   WRITE (IU06,*) ' +     =========================================    +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + A HOT START WITH CURRENT REFRACTION IS REQUESTED.+'
   WRITE (IU06,*) ' + THE FIRST CURRENT FIELD DOES NOT EXIST IN THE    +'
   WRITE (IU06,*) ' + RESTART FILE.                                    +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +               MODEL CONTINUES                    +'
   WRITE (IU06,*) ' + USING FIRST CURRENT FIELD FROM CURRENT INPUT FILE+'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF   

!     4.1.1 HOT START, CURRENTS ARE IN RESTART: CHECK FILE DATE. !

IF (CDCA.NE.' ' .AND. CDCA.NE.CDATEA) THEN
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +     WARNING ERROR SUB. PREPARE_FIRST_CURRENT.    +'
   WRITE (IU06,*) ' +     =========================================    +'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' + A HOT START WITH CURRENT REFRACTION IS REQUESTED.+'
   WRITE (IU06,*) ' + THE DATE OF CURRENT FIELD IN THE RESTART FILE    +'
   WRITE (IU06,*) ' + IS NOT THE START DATE, BECAUSE THE PREVIOUS RUN  +'
   WRITE (IU06,*) ' + WAS WITH STATIONARY CURRENTS.                    +'
   WRITE (IU06,*) ' + CURRENT DATE FROM RESTART IS CDCA = ', CDCA
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' +               MODEL CONTINUES                    +'
   WRITE (IU06,*) ' + USING FIRST CURRENT FIELD FROM CURRENT INPUT FILE+'
   WRITE (IU06,*) ' +                                                  +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
END IF

CDCA = CDATEA
CALL WAM_CURRENT (U, V, CDCA)
IF (ITEST.GE.3) THEN
   WRITE (IU06,*) '  '
   WRITE (IU06,*) '    SUB. PREPARE_FIRST_CURRENT: WAM_CURRENT DONE'
   WRITE (IU06,*) '        FIRST CURRENT FIELD PROCESSED'
   WRITE (IU06,*) '        DATE IS ................... CDCA = ', CDCA
END IF

END SUBROUTINE PREPARE_FIRST_CURRENT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_INITIAL_MODULE

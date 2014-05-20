MODULE WAM_ASSI_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS THE O-I DATA ASSIMILATION FOR THE WAM MODEL.          !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_COORDINATE_MODULE           !! COORDINATE PROCEDURES

USE WAM_GENERAL_MODULE,   ONLY:  &
&                ABORT1,         &
&                INCDATE,        &
&                DIFDATE,        &
&                PRINT_ARRAY,    &  !! PRINTER OUTPUT OF AN ARRAY.
&                OPEN_FILE

USE WAM_INTERFACE_MODULE, ONLY:  &
&       TOTAL_ENERGY,            &  !! COMPUTATION OF TOTAL ENERGY.
&       COS2_SPR                    !! COSINE**2 SPREADING FUNCTION

USE WAM_OUTPUT_SET_UP_MODULE, ONLY:  &
&       SAVE_OUTPUT_FILES           !! CLOSES AND OPENS OUTPUT FILES.

USE WAM_OUTPUT_MODULE, ONLY:     &
&       MODEL_OUTPUT_CONTROL        !! CONTROLS MODEL OUTPUT.

USE WAM_ICE_MODULE,       ONLY:  & 
&       PUT_ICE

USE WAM_TOPO_MODULE,      ONLY:  & 
&       PUT_DRY
  
use wam_mpi_comp_module,  only:  &
&       mpi_gather_block,        &  !! gather blocked data
&       mpi_gather_oifl             !! broadcast grid field onto all processes

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

USE WAM_GENERAL_MODULE,  ONLY: G, PI, ZPI, RAD, DEG
USE WAM_FILE_MODULE,     ONLY: IU06, ITEST
USE WAM_TIMOPT_MODULE,   ONLY: CDATEA, CDATEE, IDELPRO, CDTPRO, SPHERICAL_RUN
USE WAM_GRID_MODULE,     ONLY: NX, NY, NSEA, AMOWEP, AMOSOP, AMOEAP, AMONOP,   &
&                              XDELLA, ZDELLO, NLON_RG, IPER,                  &
&                              L_S_MASK, SINPH, COSPH
USE WAM_FRE_DIR_MODULE,  ONLY: FR, DFIM, DELTH, TH
USE WAM_MODEL_MODULE,    ONLY: U10, UDIR, USTAR
USE WAM_ICE_MODULE,      ONLY: ICE_RUN 
USE WAM_TOPO_MODULE,     ONLY: N_DRY 

USE WAM_OUTPUT_SET_UP_MODULE, ONLY: CDTINTT, CDTSPT, CDT_OUT, FFLAG20, FFLAG25
USE WAM_ASSI_SET_UP_MODULE,   ONLY: &
&                              IASSI, CDATAA, CDATAE, CDTASS, IDELASS,         &
&                              DIST, LMAX, SIGOBS, SIGMOD,                     &
&                              LLON, LLAT, NDIM2,                              &
&                              IU80, FILE80,                                   &
&                              FG_OUTPUT, IU30, FILE30,IU35, FILE35

use wam_mpi_module,      only: irank, nijs, nijl

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE
include 'mpif.h'

INTEGER, PRIVATE :: LENT
			
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     D.  PUBLIC INTERFACES.                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE WAMASSI                        !! SUPERVISES WAVE DATA ASSIMILATION.
   MODULE PROCEDURE WAMASSI
END INTERFACE
PUBLIC WAMASSI

INTERFACE
   SUBROUTINE READSAT (IU80, CDATE, RLAT, RLON, SWH, WS, EOFD) !! READ SAT INPUT.
   INTEGER, PARAMETER :: KIND_D = 8
   INTEGER,           INTENT(IN)  :: IU80  !! INPUT UNIT FOR MEASUREMENTS.
   CHARACTER(LEN=14), INTENT(OUT) :: CDATE !! DATE OF MEASUREMENT (YYYYMMDDHHMMSS).
   REAL (KIND=KIND_D),INTENT(OUT) :: RLAT  !! LATITUDE OF MEASUREMENT (DEGREE). 
   REAL (KIND=KIND_D),INTENT(OUT) :: RLON  !! LONGITUDE OF MEASUREMENT (DEGREES).
   REAL,              INTENT(OUT) :: SWH   !! WAVE HEIGHT.
   REAL,              INTENT(OUT) :: WS    !! WIND SPEED. 
   LOGICAL,           INTENT(OUT) :: EOFD  !! END OF DATA. (.TRUE. , .FALSE.)
   END SUBROUTINE READSAT
END INTERFACE
PUBLIC  READSAT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     E.  PRIVATE INTERFACES.                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE GRDATA            !! MEASUREMENTS TO MODEL GRID.
   MODULE PROCEDURE GRDATA
END INTERFACE
PRIVATE GRDATA

INTERFACE OIFIELD            !! OPTIMUM INTERPOLATION.
   MODULE PROCEDURE OIFIELD
END INTERFACE
PRIVATE OIFIELD

INTERFACE ANALYSE            !! COMPUTE FIELD AT (I,J) BY OPTIMUM INTERPOLATION.
   MODULE PROCEDURE ANALYSE
END INTERFACE
PRIVATE ANALYSE

INTERFACE FDUR               !! WINDSEA DURATION.
   MODULE PROCEDURE FDUR
END INTERFACE
PRIVATE FDUR

INTERFACE FUSTAR             !! FRICTION VELOCITY AND WINDSEA MEAN PERIOD.
   MODULE PROCEDURE FUSTAR
END INTERFACE
PRIVATE FUSTAR

INTERFACE FWSEA              !! PROVIDE THE WIND-SEA PART OF THE SPECTRUM.
   MODULE PROCEDURE FWSEA
END INTERFACE
PRIVATE  FWSEA

INTERFACE F4SPEC             !! PRODUCE A WINDSEA SPECTRUM FOR ASSIMILATION.
   MODULE PROCEDURE F4SPEC
END INTERFACE
PRIVATE F4SPEC

INTERFACE UPDATE             !! CONSISTENT UPDATE OF WAVE AND WIND FIELD.
   MODULE PROCEDURE UPDATE
END INTERFACE
PUBLIC UPDATE                                 !

INTERFACE UPSPEC             !! MODIFY THE SPECTRUM BY STRETCHING AND SCALING.
   MODULE PROCEDURE UPSPEC
END INTERFACE
PRIVATE UPSPEC

INTERFACE WSMFEN             !! MEAN FREQUENCY AND ENERGY OF WINDSEA SPECTRUM.
   MODULE PROCEDURE WSMFEN
END INTERFACE
PRIVATE WSMFEN

INTERFACE EN             !! DIMENSIONLESS ENERGY FROM DIMENSIONLESS FREQUENCY.
   MODULE PROCEDURE EN
END INTERFACE
PRIVATE EN 

INTERFACE YNU            !! DIMENSIONLESS FREQUENCY FROM DIMENSIONLESS ENERGY.
   MODULE PROCEDURE YNU
END INTERFACE
PRIVATE YNU

INTERFACE SYMINV       !! INVERT THE LOWER TRIANGLE OF A SYMMETRIC MATRIX.
   MODULE PROCEDURE SYMINV
END INTERFACE
PRIVATE SYMINV

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PUBLIC MODULE PROCEDURES.                                             !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WAMASSI (fl3)

! ---------------------------------------------------------------------------- !
!                                                                              !
!        WAMASSI - SUPERVISES EXECUTION OF MAIN MODULES                        !
!                  OF THE WAVE DATA ASSIMILATION.                              !
!                                                                              !
!      P. LIONELLO     ECMWF     APRIL 1990                                    !
!      J. BIDLOT       ECMWF     FEBRARY 1997 MODULE and ATMOSPHERIC COUPLING  !
!      J. BIDLOT       ECMWF     MARCH 1997  ADD SAVSTRESS AND SAVSPEC         !
!      A. Behrens      GKSS      September 2010 MPI                            !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!         TO ANALYSE THE WAVE FIELD SUBSTITUTING THE FIRST GUESS               !
!         SPECTRA WITH ANALYSED SPECTRA AND FIRST GUESS VALUES                 !
!         WITH ANALYSED VALUES IN THE GRID                                     !
!                                                                              !
!     INTERFACE.                                                               !
!     ----------                                                               !
!                                                                              !
!      CALL WAMASSI                                                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !
!                                                                              !
!        ANALYSE  - UTILITY OF OIFIELD . IT PRODUCES THE OPTIMALLY             !
!                   INTERPOLATED VALUE AT A SINGLE GRIDPOINT.                  !
!        FDUR     - COMPUTES WINDSEA DURATION.                                 !
!        FEMEAN   - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.          !
!        FUSTAR   - COMPUTES ANALYSED USTAR ACCORDING TO ANALYSED              !
!                     WINDSEA ENERGY.                                          !
!        FWSEA    - DETERMINES IF WINDSEA IS PRESENT IN THE SPECTRUM,,         !
!                     COMPUTING EVENTUALLY WINDSEA ENERGY AND MEAN FREQ.       !
!        GRDATA   - TRANSFERS MEASUREMENTS TO THE MODEL GRID.                  !
!        INCDATE  - UPDATE DATE TIME GROUP.                                    !
!        OIFIELD  -  MERGES MODEL AND MEASUREMENTS BY O.I.                     !
!        TOTAL_ENERGY   - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.      !
!        UPDATE   - ANALYSES THE WAM MODEL SPECTRA.                            !
!        UPSPEC   - PRODUCES THE UPDATED SPECTRUM.                             !
!        WSMFEN   - UTILITY OF FWSEA , COMPUTING  WINDSEA ENERGY AND           !
!                   DURATION IN EACH GRID POINT.                               !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!          NONE                                                                !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------

real, dimension (:,:,:), intent(inout) :: fl3
     
! ---------------------------------------------------------------------------- !
!                                                                              !
!      LOCAL VARIABLES.
!      ----------------

REAL, PARAMETER :: ZMISS =  -999.       !! MISSING VALUE.

REAL :: ETOT(1:nsea)                    !! FIRST GUESS WAVE ENERGIES
REAL :: ETOIB(1:nsea)                   !! BLOCK OF UPDATED WAVE ENERGIES.  
REAL :: USTAR_OLD(nijs:nijl)            !! BLOCK OF FIRST GUESS  USTARS. 
REAL :: USOIB(1:nsea)                   !! UPDATED USTARS. 
real :: xustar(1:nsea)                  !! complete ustar
real :: xu10(1:nsea)                    !! complete u10
REAL :: CDG(NX,NY)                      !! DRAG COEF.
real :: xetot(nijs:nijl)

REAL, ALLOCATABLE, DIMENSION(:,:) :: USOI    !! USTAR FIELD FROM O.I.
REAL, ALLOCATABLE, DIMENSION(:,:) :: USME    !! USTAR FIELD FROM MEASUREMENTS.
REAL, ALLOCATABLE, DIMENSION(:,:) :: WHOI    !! WAVE HEIGHT FIELD FROM O.I.
REAL, ALLOCATABLE, DIMENSION(:,:) :: WHME    !! WAVE HEIGHT FIELD FROM MEASURMENTS 
REAL, ALLOCATABLE, DIMENSION(:,:) :: WHGTTG  !! FIRST GUESS WAVE HEIGHT. 
REAL, ALLOCATABLE, DIMENSION(:,:) :: USTARG  !! FIRST GUESS USTAR.

character (len=100) :: titl                  !! WORK STRING FOR PRINTING
integer :: isend, irecv

! ----------------------------------------------------------------------       !
   
WRITE (IU06,*) '  '
WRITE (IU06,*) ' START OF DATA ASSIMILATION: DATE IS CDTPRO: ', CDTPRO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. MODEL OUTPUT OF FIRST GUESS TO DISK AND/OR PRINTER.                   !
!        ---------------------------------------------------                   !

IF (FG_OUTPUT) THEN
   IF (CDTINTT.EQ.CDTPRO .OR. CDTSPT.EQ.CDTPRO) THEN
      CALL MODEL_OUTPUT_CONTROL (FL3, IU30, IU35)
      if (cdt_out==cdtpro .and. cdtpro<cdatae) then
         CALL SAVE_OUTPUT_FILES (IU30, FILE30, IU35, FILE35)
      END IF
      IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. WAMASSI: MODEL_OUTPUT_CONTROL DONE'
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN PARAMETERS.                                              !
!        ------------------------                                              !

ALLOCATE (WHGTTG(NX,NY), USTARG(NX,NY))

CALL TOTAL_ENERGY (FL3, ETOT(nijs:nijl))
    
xetot = etot(nijs:nijl)
   
isend = 1
irecv = 1
    
call mpi_gather_block (irecv, xetot, etot)   !! gather fields on irecv
call mpi_gather_block (irecv, ustar, xustar)
call mpi_gather_block (irecv, u10, xu10)
    
if (irank==irecv) then
   WHGTTG = UNPACK (4.*SQRT(ETOT), L_S_MASK, ZMISS)   !! MAKE HS GRIDDED FIELD.
   USTARG = UNPACK (xUSTAR, L_S_MASK, ZMISS)          !! MAKE USTAR GRIDDED FIELD.
   CDG = UNPACK((xUSTAR**2 +0.0001)/(xU10**2+0.01), L_S_MASK, ZMISS)
endif                           
    
call mpi_gather_oifl (isend, whgttg)         !! send fields to all processes
call mpi_gather_oifl (isend, ustarg)
call mpi_gather_oifl (isend, cdg)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. MEASUREMENT ARE TRANSFERRED TO THE MODEL GRID.                        !
!        ----------------------------------------------                        !

ALLOCATE (WHME(NX,NY), USME(NX,NY))
ALLOCATE (WHOI(NX,NY), USOI(NX,NY))

CALL GRDATA (WHME, USME, CDG)

IF (ITEST.GE.4) THEN
   TITL = ' GRDATA :    WAVE HEIGHT ( METRES )'
   CALL PRINT_ARRAY (IU06, CDTPRO, TITL, WHME,                                 &
&                                AMOWEP, AMOSOP, AMOEAP, AMONOP, 10.)
   TITL = ' GRDATA USTAR      ( METRES/SECOND )'
   CALL PRINT_ARRAY (IU06, CDTPRO, TITL, USME,                                 &
&                                AMOWEP, AMOSOP, AMOEAP, AMONOP, 100.)
END IF
IF (ITEST.GE.3) WRITE(IU06,*) '      SUB. WAMASSI: MEASUREMENTS ---> GRID'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. MEASUREMENTS AND MODEL ARE MERGED BY OPTIMUM INTERPOLATION.           !
!        -----------------------------------------------------------           !

CALL OIFIELD (USME, USTARG, USOI)

IF (ITEST.GE.3)  THEN
   WRITE(IU06,*) '      SUB. WAMASSI: SUB. OIFIELD DONE FOR WINDS'
END IF

CALL OIFIELD (WHME, WHGTTG, WHOI)
IF (ITEST.GE.3) THEN
   WRITE(IU06,*) '      SUB. WAMASSI: SUB. OIFIELD DONE FOR WAVE HEIGHTS'
END IF

IF (ITEST.GE.4) THEN
   TITL = ' OIFIELD:    WAVE HEIGHT ( METRES )'
   CALL PRINT_ARRAY (IU06, CDTPRO, TITL, WHOI,                                 &
&                                AMOWEP, AMOSOP, AMOEAP, AMONOP, 10.)
   TITL = ' OIFIELD: USTAR    ( METRES/SECOND )'
   CALL PRINT_ARRAY (IU06, CDTPRO, TITL, USOI,                                 &
&                                AMOWEP, AMOSOP, AMOEAP, AMONOP, 100.)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. ANALYSING THE SPECTRA.                                                !
!        ----------------------                                                !

ETOIB = PACK (WHOI, L_S_MASK)    !! BLOCK UPDATED WAVE HEIGHTS
WHERE (ETOIB.GT.0.) 
   ETOIB = ETOIB**2/16.          !! UPDATED WAVE ENERGIES
ELSEWHERE
   ETOIB = -99.
ENDWHERE

USTAR_OLD = USTAR                !! SAVE FIRST GUESS WINDS
USOIB = PACK (USOI, L_S_MASK)    !! BLOCK UPDATED WINDS

CALL UPDATE (FL3, ETOT(nijs:nijl), ETOIB(nijs:nijl), USOIB(nijs:nijl),         &
&            USTAR, UDIR)

IF (ITEST.GE.3) WRITE(IU06,*) '      SUB. WAMASSI: SUB. UPDATE DONE'

!           THE WIND FIELD CURRENTLY USED IS UPDATED.                          !

U10 = USTAR /(USTAR_OLD+0.00000001)* U10

IF (ITEST.GE.3) WRITE(IU06,*) '      SUB. WAMASSI: WINDS REPLACED '

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. UPDATED ASSIMILATION TIME.                                            !
!        ---------------------------                                           !

CALL INCDATE (CDTASS, IDELASS)

IF (CDTASS.GT.CDATAE) THEN
   IASSI = 0
   WRITE(IU06,*) '   END OF ASSIMILATION FOR THIS MODEL RUN'
   IF (FFLAG20) CLOSE (UNIT=IU30, STATUS ="KEEP")
   IF (FFLAG25) CLOSE (UNIT=IU35, STATUS ="KEEP")
ELSE
   WRITE(IU06,*) '   END OF ASSIMILATION NEXT ASSIMILATION AT ', CDTASS
END IF
WRITE (IU06,*) ' '

! ---------------------------------------------------------------------------- !
!                                                                              !
!     7. FREE UP SOME MEMORY.                                                  !
!        ----------------------                                                !
     
DEALLOCATE (USOI, USME, WHOI, WHME)

END SUBROUTINE WAMASSI

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     F. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE GRDATA (WHME, USME, cdg)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     GRDATA - TRANSFERS MEASUREMENTS TO MODEL GRID,                           !
!              DISCARDING UNRELAIABLE DATA.                                    !
!                                                                              !
!     P.LIONELLO     ECMWF       APRIL 1990                                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       PREPARE THE DATA FOR OPTIMAL INTERPOLATION.                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        A BOX OF SIZE DELLA*DELLO IS CENTERED IN EACH GRID POINT.             !
!        THE AVERAGE VALUE OF THE SATELLITE MEASUREMENTS INSIDE THE            !
!        BOX IS TAKEN TO UPDATE THE WAM MODEL. THE VALUE IS DISCARDED          !
!        IF THE ROOT MEAN SQUARE ERROR IS TOO BIG.                             !
!        MOREOVER THE VALUE IS DISCARDED IF THERE ARE LESS THEN FOUR           !
!        MEASUREMENTS INSIDE THE BOX.                                          !
!                                                                              !
!      EXTERNALS.                                                              !
!      ----------                                                              !
!                                                                              !
!         READSAT  - READ A MEASUREMENT.                                       !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.

REAL, INTENT(OUT) :: WHME(NX,NY)  !! FIELD OF MEASURED WAVE HEIGHTS.
REAL, INTENT(OUT) :: USME(NX,NY)  !! FIELD OF MEASURED FRICTION VELOCITIES.
real, intent(in)  :: cdg(nx,ny)   !! drag coefficient

! ---------------------------------------------------------------------------- !
! 
!      LOCAL VARIABLES.

REAL, PARAMETER :: ZMISS =  -999.    !! MISSING VALUE
INTEGER         :: IFAIL
LOGICAL         :: EOFD

CHARACTER*14    :: CBEGINDT  !! BEGINNING DATE OF THE TIME WINDOW
                             !! DURING WHICH DATA ARE USED FOR ANALYSIS.
CHARACTER*14    :: CENDDT    !! END DATE OF THE TIME WINDOW
                             !! DURING WHICH DATA ARE SAVE FOR ANALYSIS.
CHARACTER*14    :: CDATE     !! DATE OF DATA RECORD
CHARACTER*14    :: CDATEO    !! DATE OF PREVIOUS DATA RECORD

REAL    :: WHSE(NX,NY)    !! WAVE HEIGHT ROOT MEAN SQUARE ERROR.
REAL    :: USSE(NX,NY)    !! USTAR ROOT MEAN SQUARE ERROR.
INTEGER :: NUMBWH(NX,NY)  !! NUMBER OF HS MEASUREMENTS INSIDE A BOX.
INTEGER :: NUMBUS(NX,NY)  !! NUMBER OF US MEASUREMENTS INSIDE A BOX.
INTEGER :: NR(NX,NY)      !! 

INTEGER :: ICOUNT, ICOUNTT, ICOUNTD, ICOUNTB, ICOUNTQCF
INTEGER :: IOLD, JOLD, J, I, NUWH, NUUS
INTEGER :: NN, NWH, NUS, ISHIFT
INTEGER :: M_RLAT, M_RLON
REAL (KIND=KIND_D) :: RLAT, RLON
REAL               :: USO, US, WSO, WS, SWHO, SWH
REAL               :: WHSE2, USSE2, RFU, RFWM, XX

LOGICAL :: ID_MAP(NX,NY)     !! MAP OF ICE, DRY AND LAND POINTS
REAL    :: ID_MAP_AC(NSEA)   !! MAP OF ICE AND DRY POINTS

!----------------------------------------------------------------------------- !
!                                                                              !
!                                                                              !
!     1. FIX TIME WINDOW FOR DATA INPUT.                                       !
!        --------------------------------                                      !

CBEGINDT = CDTPRO
CALL INCDATE (CBEGINDT,-IDELASS/2)
CENDDT = CBEGINDT
CALL INCDATE (CENDDT,IDELASS)

WRITE(IU06,*) '    SUB. GRDATA : SELECTING MEASUREMENTS '
WRITE(IU06,*) '                  PERIOD FROM ', CBEGINDT,' TO ', CENDDT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INITIALIZE ARRAYS, COUNTERS.                                          !
!        ----------------------------                                          !

ICOUNT    = 0  !! NUMBER OF DATA IN BOXES AROUND SEA POINTS
ICOUNTT   = 0  !! NUMBER OF DATA IN SELECTED PERIOD
ICOUNTD   = 0  !! NUMBER OF RECORDS READ FROM DATA FILE
ICOUNTB   = 0  !! NUMBER OF DATA BEFORE SELECTED PERIOD
ICOUNTQCF = 0  !! NUMBER OF GOOD DATA ABOVE LAND, ICE, DRY

WHME = 0.
USME = 0.
WHSE = 0.
USSE = 0.
NUMBWH = 0
NUMBUS = 0
NR = 0

IOLD = 0
JOLD = 0
CDATEO = ' '

!     MAKE A MAP OF ICE, DRY AND LANDPOINTS.

IF (ICE_RUN .OR. N_DRY.GT.0) THEN
   ID_MAP_AC = 0.
   IF (ICE_RUN) CALL PUT_ICE (ID_MAP_AC, 1.)
   IF (N_DRY.GT.0) CALL PUT_DRY (ID_MAP_AC, 1.)
   ID_MAP = UNPACK(ID_MAP_AC, L_S_MASK, 1.) .EQ. 1.
ELSE
   ID_MAP = .NOT. L_S_MASK
END IF

!      OPEN INPUT FILE.

IFAIL = 0
CALL OPEN_FILE (IU06, IU80, FILE80, CDTPRO, 'OLD', IFAIL, FORMA='FORMATTED')
IF (IFAIL.NE.0) THEN
   CALL ABORT1
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. READ THE MEASUREMENT AND CUMULATE THE DATA ON GRID POINTS.            !
!        ----------------------------------------------------------            !

DATA_LOOP: DO

   CALL READSAT (IU80, CDATE, RLAT, RLON, SWH, WS, EOFD)
   M_RLAT = DEG_TO_M_SEC (RLAT) 
   M_RLON = DEG_TO_M_SEC (RLON)

   IF (EOFD) EXIT DATA_LOOP

   ICOUNTD = ICOUNTD+1
   IF (CDATE .LT. CBEGINDT) THEN       !! CHECK DATE.
      ICOUNTB = ICOUNTB + 1
      CYCLE DATA_LOOP
   END IF

   IF (CDATE.GT.CENDDT) EXIT DATA_LOOP

!       DATE IS INSIDE THE TIME WINDOW.                                        !

   ICOUNTT = ICOUNTT+1                                                    
   J = NINT(REAL(M_RLAT-AMOSOP)/REAL(XDELLA)) + 1 

   IF ((J.GT.NY) .OR. (J.LT.1)) CYCLE DATA_LOOP
   
   
   I = MOD(M_RLON-AMOWEP+2*M_S_PER,M_S_PER)                                        
   I = NINT(REAL(I)/REAL(ZDELLO(J))) + 1                                                
   IF (IPER .AND. I.EQ.NLON_RG(J)+1) I = 1                                                  
   IF ((I.GT.NLON_RG(J)) .OR. (I.LT.1)) CYCLE DATA_LOOP

!     THE INDICES ARE INSIDE THE GRID THE MEASUREMENT IS CUMULATED             !

!     IF THE GRID POINT IS SET TO BE ICE OR DRY OR                             !
!     MEASURENTS ARE VERY CLOSE TO LAND AND PASSED THE QC AS                   !
!     BEING ON SEA BUT REGARDED AS BEING ON LAND BY GRDATA.                    !
!     THESE ARE ELIMINATED QUIETLY HERE.                                       !

   IF (ID_MAP(I,J)) THEN
      ICOUNTQCF = ICOUNTQCF + 1
      CYCLE DATA_LOOP
   END IF

   ICOUNT = ICOUNT+1

   IF (I.NE.IOLD .OR. J.NE.JOLD) THEN

!     NEW GRID BOX                                                             !
!     DATA ARE IGNORED IF TIME STEP IS LESS THEN 3S AND                        !
!     WAVE HEIGHT JUMPS MORE THEN 3M.                                          !

      IF (CDATEO .NE. ' ') THEN
         CALL DIFDATE (CDATEO, CDATE, ISHIFT)
      ELSE
         ISHIFT = 999
      END IF   
      IF (ISHIFT.LT.3) THEN 
         IF (SWH-SWHO.GT.2.0) THEN
            NR(I,J) = NR(I,J) + 1
            CYCLE DATA_LOOP
         END IF
      END IF
      IOLD = I
      JOLD = J
   ELSE

!     OLD GRID BOX                                                             !
!     SUDDEN POSITVE STEPS ARE REFUSED.                                        !
!                                                                              !
      IF (SWH-SWHO.GT.1.0) THEN
         NR(I,J) = NR(I,J) + 1
         CYCLE DATA_LOOP
       END IF
   END IF

   IF (SWH.GE.0.) THEN
      WHME(I,J) = WHME(I,J) + SWH
      WHSE(I,J) = WHSE(I,J) + SWH**2
      NUMBWH(I,J) = NUMBWH(I,J) + 1
      IF (WS.GE.0.) THEN
          US = SQRT(CDG(I,J))*WS
          USME(I,J) = USME(I,J) + US
          USSE(I,J) = USSE(I,J) + US**2
          NUMBUS(I,J) = NUMBUS(I,J) + 1
      END IF
   END IF

!     A PREVIOUS SUDDEN NEGATIVE STEP IS COMPENSATED.                          !

   IF (SWHO-SWH.GT.2.0) THEN
      NR(I,J) = NR(I,J) + 1
      WHME(I,J) = WHME(I,J) - SWHO
      WHSE (I,J) = WHSE(I,J) - SWHO**2
      NUMBWH(I,J) = NUMBWH(I,J) - 1
      IF (WSO.GE.0.) THEN
         USO = SQRT(CDG(I,J))*WSO
         USME(I,J) = USME(I,J) - USO
         USSE (I,J) = USSE(I,J) - USO**2
         NUMBUS(I,J) = NUMBUS(I,J) - 1
      END IF
   END IF
   SWHO = SWH
   WSO = WS
   CDATEO = CDATE

END DO DATA_LOOP

!     CLOSE THE DATA FILE                                                      !

CLOSE (IU80)

NUWH = SUM(NUMBWH)
NUUS = SUM(NUMBUS)

IF (ITEST.GE.4 ) THEN
   WRITE (IU06,*) '    SUB. GRDATA: ALL DATA READ AND UNIT CLOSED'
END IF
WRITE (IU06,*) '    SUB: GRDATA: INPUT STATISTICS:'
WRITE (IU06,*) '       NUMBER OF RECORDS READ FROM DATA FILE IS ....', ICOUNTD
WRITE (IU06,*) '       NUMBER OF DATA IN SELECTED PERIOD IS ........', ICOUNTT
WRITE (IU06,*) '       NUMBER OF DATA BEFORE SELECTED PERIOD IS ....', ICOUNTB
WRITE (IU06,*) '       NUMBER OF DATA IN BOXES AROUND SEA POINTS IS ', ICOUNT
WRITE (IU06,*) '       NUMBER OF GOOD DATA ABOVE LAND, ICE, DRY ... ', ICOUNTQCF
WRITE (IU06,*) '       NUMBER OF DATA ACCEPTED BY "GRDATA"  IS FOR WH :',      &
&                                                   NUWH,' AND FOR US: ', NUUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. DETERMINE AVERAGE VALUE AND ERROR OF THE MEASUREMENTS IN THE          !
!        BOX AROUND GRIDPOINT I,J.                                             !
!        ------------------------------------------------------------          !
!                                                                              !

DO J = 1,NY
   DO I = 1,NLON_RG(J)
      NN = NUMBWH(I,J)
      IF (NN.GE.4) THEN
         WHME(I,J) = WHME(I,J)/REAL(NN)                   !! MEAN VALUE.
         WHSE2 =  (WHSE(I,J)-WHME(I,J)**2*REAL(NN)) / (REAL(NN)-1.)
         WHSE(I,J) = SQRT( MAX(WHSE2,0.))                 !! STANDARD DEVIATION.
      END IF
   END DO
END DO

DO J = 1,NY
   DO I = 1,NLON_RG(J)
      NN = NUMBUS(I,J)
      IF (NN.GE.4) THEN
         USME(I,J) = USME(I,J)/REAL(NN)                  !! MEAN VALUE.
         USSE2 =  (USSE(I,J)-USME(I,J)**2*REAL(NN)) / (REAL(NN)-1.)
         USSE(I,J) = SQRT( MAX(USSE2,0.))                !! STANDARD DEVIATION.
      END IF
   END DO
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. FLAG GRID POINTS, IF THERE ARE FEW OR STRONGLY SCATTERED DATA.        !
!        --------------------------------------------------------------        !

NWH = 0
DO J = 1,NY
   DO I = 1,NLON_RG(J)
      NN = NUMBWH(I,J)
      IF (NN.EQ.0) THEN                       !! IF THERE ARE NO MEASUREMENTS
         WHME(I,J) = -3.
      ELSE IF (NN.GT.0 .AND. NN.LT.4) THEN    !! IF THERE ARE FEW MEASUREMENTS.
         WHME(I,J) = -2.
      ELSE IF (NN.GE.4) THEN                  !! IF THE VARIANCE IS TOO LARGE
                                              !! OR SWH IS TOO SMALL  
					      !! OR THERE ARE TOO MANY SPIKES.
         XX = REAL(NR(I,J))/REAL(NN)
         RFWM = MAX(0.5, 0.25*WHME(I,J))
         IF (WHSE(I,J).GT.RFWM .OR. WHME(I,J).LT.0.50 .OR. XX.GT.0.1) THEN
            WHME(I,J) = -WHME(I,J)
         ELSE
            NWH = NWH+1
         END IF
      END IF
   END DO
END DO

NUS = 0
DO J = 1,NY
   DO I = 1,NLON_RG(J)
      NN = NUMBUS(I,J)
      IF (NN.EQ.0) THEN                        !! IF THERE ARE NO MEASUREMENTS 
         USME(I,J) = -3.
      ELSE IF (NN.GT.0 .AND. NN.LT.4) THEN     !! IF THERE ARE FEW MEASUREMENTS.
         USME(I,J) = -2.
      ELSE IF (NN.GE.4) THEN                   !! IF THE VARIANCE IS TOO LARGE
                                               !! OR THERE IS NO CORRESPONDING SWH.
         IF (USME(I,J).GT.0.0001) RFU = USSE(I,J)/USME(I,J)
         IF (RFU.GE.0.5 .OR. WHME(I,J).LE.0.) THEN
            USME(I,J) = -USME(I,J)
         ELSE
            NUS = NUS + 1
         END IF
      END IF
   END DO
END DO

WRITE (IU06,*) '    SUB. GRDATA: MEASUREMENTS ARE SELECTED:'
WRITE (IU06,*) '                 WAVE HEIGHT AT NWH = ', NWH,' GRIDPOINTS'
WRITE (IU06,*) '                 WIND SPEED  AT NUS = ', NUS,' GRIDPOINTS'

END SUBROUTINE GRDATA

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE OIFIELD (XME, XMO, XOI)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     OIFIELD - OPTIMUM INTERPOLATION.                                         !
!                                                                              !
!     P.LIONELLO     ECMWF       APRIL 1990                                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO PRODUCE A MAP OF THE FIELD X, MERGING MEASUREMENT                   !
!       AND MODEL , BY OPTIMUM INTERPOLATION. THE ARRAY XOI                    !
!       AT THE END OF THE SUBROUTINE CONTAINS THE VALUES                       !
!       TO BE USED TO ANALYSE THE SPECTRA, HAVING NEGATIVE RETURN              !
!       CODES WHERE O.I. PRODUCED NO RESULTS.                                  !
!       THE FIRST GUESS FIELD XMO IS NOT MODIFIED                              !
!       IN THIS SUBROUTINE.                                                    !
!                                                                              !
!      EXTERNALS.                                                              !
!      ----------                                                              !
!                                                                              !
!        SYMINV - INVERT LOWER TRIANGLE OF A SYMMETRIC POSITIVE DEFINITE MATRIX!
!        ANALYSE - COMPUTE FIELD AT (I,J) BY OPTIMUM INTERPOLATION.            !
!                                                                              !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.

REAL,    INTENT(IN)  :: XME(NX,NY)   !! FIELD FROM MEASUREMENTS.
REAL,    INTENT(IN)  :: XMO(NX,NY)   !! FIELD FROM MODEL (FIRST GUESS).
REAL,    INTENT(OUT) :: XOI(NX,NY)   !! FIELD FROM O.I.

! ---------------------------------------------------------------------------- !
! 
!     LOCAL VARIABLES

INTEGER :: IMEAS(NDIM2), JMEAS(NDIM2)
REAL    :: D, P(NDIM2,NDIM2), V(NDIM2)

INTEGER :: J, I, NOBS, J1, J2, I1, I2, IX, JX, NEWI

REAL :: D0BS, COND, EPS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALISE SEARCH DISTANCE AND CORRELATION LENGTH.                    !
!     -----------------------------------------------------                    !

!       THE MEASUREMENTS ARE STORED IN THE O.I. FIELD.                         !
!       NEGATIVE VALUES WILL REMAIN WHERE O.I. WILL PRODUCE NO RESULTS         !

XOI = XME

!       DETERMINE OBSERVATION CORRELATION MATRIX  D                            !

D = (SIGOBS/SIGMOD)**2
V(1:NDIM2) = 0.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LOOP OVER GRID POINTS.                                                !
!        ----------------------                                                !

BLOCK: DO J = 1,NY
   POINT: DO I = 1,NLON_RG(J)

      IF (XMO(I,J).LT.0.01) CYCLE POINT
      NOBS = 0
      IMEAS = 0
      JMEAS = 0

!      FIND NUMBER AND POSITION OF OBSERVATIONS IN A BOX                       !
!      OF SIZE 2*LMAX+1 AROUND ANALYSIS POINT I,J.                             !

      DO JX = MAX(J - LLAT, 1), MIN(J + LLAT, NY)
         NEWI = NINT((I-1)*REAL(ZDELLO(J))/REAL(ZDELLO(JX)))+1
         DO IX = MAX(1,NEWI-LLON(JX)), MIN(NLON_RG(JX),NEWI+LLON(JX))
            IF (XME(IX,JX) .LE. 0.  ) CYCLE
            IF (XMO(IX,JX) .LT. 0.01) CYCLE
            NOBS = NOBS+1
            IMEAS(NOBS) = IX
            JMEAS(NOBS) = JX
         END DO

         IF (IPER) THEN
            DO IX = 1, NEWI+LLON(JX)-NLON_RG(JX)
               IF (XME(IX,JX).LE.0.  ) CYCLE
               IF (XMO(IX,JX).LT.0.01) CYCLE
               NOBS = NOBS+1
               IMEAS(NOBS) = IX
               JMEAS(NOBS) = JX
            END DO

            DO IX = NEWI-LLON(JX)+NLON_RG(JX), NLON_RG(JX)
               IF (XME(IX,JX).LE.0.  ) CYCLE
               IF (XMO(IX,JX).LT.0.01) CYCLE
               NOBS = NOBS+1
               IMEAS(NOBS) = IX
               JMEAS(NOBS) = JX
            END DO
         END IF
      END DO

      IF (NOBS.EQ.0) CYCLE POINT

!     DETERMINE MODEL CORRELATION MATRIX AT OBSERVATION POINTS.                !


      DO IX = 1,NOBS
         DO JX = IX+1,NOBS
            I1 = IMEAS(IX)
            I2 = IMEAS(JX)
            J1 = JMEAS(IX)
            J2 = JMEAS(JX)

            IF (SPHERICAL_RUN) THEN
               D0BS = REAL((I1-1)*ZDELLO(J1)-(I2-1)*ZDELLO(J2))
               D0BS = D0BS / REAL(M_DEGREE)*RAD
               D0BS = COS(D0BS)*COSPH(J1)*COSPH(J2) + SINPH(J1)*SINPH(J2)
               D0BS = ACOS(D0BS)
            ELSE
               D0BS = (REAL((I1-1)*ZDELLO(J1)-(I2-1)*ZDELLO(J2))/REAL(M_DEGREE))**2      &
&                   + (REAL((J1-J2)*XDELLA)/REAL(M_DEGREE))**2

               D0BS =  SQRT(D0BS)*RAD
            END IF
            P(IX,JX) = EXP(-D0BS/DIST)
            P(JX,IX) = P(IX,JX)
         END DO
         P(IX,IX) = 1.
      END DO

!   DETERMINE M=P+D                                                            !

      DO IX = 1,NOBS
         P(IX,IX) = P(IX,IX)+D
      END DO

!   INVERSE MATRIX M                                                           !

      COND = 0.
      EPS = 0.
      CALL SYMINV (P(1:NOBS,1:NOBS), NOBS, NOBS, COND, V)

!   PRODUCE FIELD AT POINT I,J ACCORDING TO OPTIMUM INTERPOLATION              !

      CALL ANALYSE (P(1:NOBS,1:NOBS), NOBS, XME, XMO, XOI(I,J),                &
&                   IMEAS(1:NOBS), JMEAS(1:NOBS), I, J, DIST)

   END DO POINT
END DO BLOCK

END SUBROUTINE OIFIELD

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE ANALYSE (XM, NOBS, XME, XMO, XOI, IMEAS, JMEAS, I, J, DIST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      ANALYSE - COMPUTE FIELD AT (I,J) BY OPTIMUM INTERPOLATION.              !
!                                                                              !
!     PIERO LIONELLO      ECMWF     JUNE 1990                                  !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE THE FIELD IN THE POINT (I,J) BY OPTIMUM INTERPOLATION.         !
!       ANALYSIS IS DONE FOR ONE GRID POINT.                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.

REAL,    INTENT(INOUT) :: XM(NOBS,NOBS) !! INVERSE OF (P+M) : P=MODEL CORRELATION
                                        !! MATRIX , M=MEASUREMENTS CORRELATION MATRIX
INTEGER, INTENT(IN)    :: NOBS          !! NUMBER OF OBSERVATION IN THE
                                        !! BOX AROUND POINT (I,J)
REAL,    INTENT(IN)    :: XME(NX,NY)    !! FIELD FROM MEASUREMENT
REAL,    INTENT(IN)    :: XMO(NX,NY)    !! FIELD FROM MODEL
REAL,    INTENT(OUT)   :: XOI           !! RESULT OF O.I.
INTEGER, INTENT(IN)    :: IMEAS(NOBS)   !! X-INDICES OF MEASUREMENTS
INTEGER, INTENT(IN)    :: JMEAS(NOBS)   !! Y-INDICES OF MEASUREMENTS
INTEGER, INTENT(IN)    :: I             !! POINT WHERE O.I.FIELD IS PRODUCED
INTEGER, INTENT(IN)    :: J             !! POINT WHERE O.I.FIELD IS PRODUCED
REAL,    INTENT(IN)    :: DIST          !! RADIUS OF INFLUENCE IN RAD

! ---------------------------------------------------------------------------- !
! 
!      LOCAL VARIABLES

INTEGER         :: IOBS, I1, I2, J1, J2, JOBS
REAL            :: D0BS, W

! ---------------------------------------------------------------------------- !

DO I1 = 1,NOBS-1
   DO J1 = I1+1,NOBS
      XM(I1,J1) = XM(J1,I1)
   END DO
END DO

XOI = XMO(I,J)
I1 = I
J1 = J 
DO IOBS = 1,NOBS
   I2 = IMEAS(IOBS)
   J2 = JMEAS(IOBS)

   IF ((I1.EQ.I2) .AND. (J1.EQ.J2)) THEN
      W = 1.
   ELSE
      IF (SPHERICAL_RUN) THEN
         D0BS = REAL((I1-1)*ZDELLO(J1)-(I2-1)*ZDELLO(J2))
         D0BS = D0BS/REAL(M_DEGREE)*RAD              
         D0BS = COS(D0BS)*COSPH(J1)*COSPH(J2) + SINPH(J1)*SINPH(J2)
         D0BS = ACOS(D0BS)
      ELSE
         D0BS = (REAL((I1-1)*ZDELLO(J1)-(I2-1)*ZDELLO(J2))/REAL(M_DEGREE))**2            &
&                   + (REAL((J1-J2)*XDELLA)/REAL(M_DEGREE))**2
         D0BS = SQRT(D0BS)*RAD
      END IF
      W       = EXP(-D0BS/DIST)
   END IF

   DO JOBS = 1,NOBS
      I2 = IMEAS(JOBS)
      J2 = JMEAS(JOBS)
      XOI = XOI + W*XM(IOBS,JOBS)*(XME(I2,J2)-XMO(I2,J2))
   END DO
END DO

END SUBROUTINE ANALYSE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FDUR (EWFG, T, USMO)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      FDUR* - EVALUATE THE WINDSEA DURATION FROM THE MODEL GROTH CURVE.       !
!                                                                              !
!     P.LIONELLO     ECMWF       APRIL 1990                                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       EVALUATE THE WINDSEA DURATION FROM THE MODEL GROWTH CURVE.             !
!                                                                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        THE DURATION IS DERIVED BY THE GROWTH CURVE                           !
!                ESTAR=EO*TANH(A*TSTAR**B).                                    !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.

REAL,    INTENT(IN)  :: EWFG(:)  !! WINDSEA ENERGY (FIRST GUESS). 
REAL,    INTENT(OUT) :: T(:)     !! WINDSEA DURATION (OUTPUT).
REAL,    INTENT(IN)  :: USMO(:)  !! FIRST GUESS USTAR.

! ---------------------------------------------------------------------------- !
! 
!      LOCAL VARIABLES

REAL, PARAMETER :: EO = 955.          !! PARAMETER OF THE ENERGY GROWTH CURVE.
REAL, PARAMETER :: A = 6.02*10.**(-5) !! PARAMETER OF THE ENERGY GROWTH CURVE.
REAL, PARAMETER :: B = 0.695          !! PARAMETER OF THE ENERGY GROWTH CURVE.

INTEGER :: IJ
REAL    :: TSTAR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALIZE DURATION.                                                  !
!        -------------------                                                   !

T = -9.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LOOP OVER POINTS.                                                     !
!        -----------------                                                     !

DO IJ = 1,SIZE(EWFG)
   IF (EWFG(IJ).LE.0.01) CYCLE

!        DETERMINE THE EFFECTIVE DURATION OF THE WINDSEA ACCORDING TO          !
!        THE WAM MODEL GROWTH CURVE, IF THERE WINDSEA IS PRESENT.              !

   TSTAR = EWFG(IJ)*G**2/USMO(IJ)**4
   TSTAR = MIN (TSTAR/EO, .999999)
   TSTAR = (1.+TSTAR)/(1.-TSTAR)
   TSTAR = (ALOG(TSTAR)/(2.*A))**(1./B)
   T(IJ) = USMO(IJ)/G*TSTAR
END DO

END SUBROUTINE FDUR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FUSTAR (USMO, USA, USOI, EWFG, EWOI, EWA, T, FMWA)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    FUSTAR - ESTIMATE THE ANALYSED FRICTION VELOCITY FROM THE  DURATION       !
!             AND THE "MEASURED" WINDSEA ENERGY.                               !
!             ESTIMATE THE ANALYSED WINDSEA MEAN FREQUENCY                     !
!             FROM THE "MEASURED" WINDSEA ENERGY                               !
!                                                                              !
!     P.LIONELLO     ECMWF       FEBRUARY 1989                                 !
!         (FROM MODIFICATION OF A CODE BY P.JANSSEN                            !
!           AND P.LIONELLO - SUMMER '87 - )                                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        ESTIMATE THE ANALYSED FRICTION VELOCITY FROM THE URATION AND          !
!        THE "MEASURED" WINDSEA ENERGY. ESTIMATE THE ANALYSED WINDSEA          !
!        MEAN FREQUENCY FROM THE "MEASURED" WINDSEA ENERGY.                    !
!                                                                              !
!      METHOD.                                                                 !
!      -------                                                                 !
!                                                                              !
!       THE DURATION OF THE WIND SEA IS USED TO COMPUTE THE FRICTION           !
!       VELOCITY BY NEWTONS METHOD FROM THE MODEL DURATION CURVE AND           !
!       THE MEASURED ENERGY.  (FIVE ITERATIONS ARE DONE)                       !
!       THIS VALUE IS USED TO UPDATE THE MODEL IF IT IS                        !
!       IN REASONABLE AGREEMENT WITH THE MEASUREMENT FRICTION VELOCITY.        !
!       IF THERE IS NO AGREEMENT THE MEASURED WIND SPEED IS USED               !
!       TO COMPUTE THE WIND SEA.                                               !
!                                                                              !
!      EXTERNALS.                                                              !
!      ----------                                                              !
!                                                                              !
!        YNU     - DIMENSIONLESS FREQUENCY AS FUNCTION OF DIMENSIONLESS ENERGY.!
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.

REAL,    INTENT(IN)  :: USMO(:) !! FIRST GUESS FRICTION VELOCITY.
REAL,    INTENT(OUT) :: USA(:)  !! ANALISED FRICTION VELOCITY.
REAL,    INTENT(IN)  :: USOI(:) !! MEASURED FRICTION VELOCITY.
REAL,    INTENT(IN)  :: EWFG(:) !! WINDSEA ENERGY (FIRST GUESS).
REAL,    INTENT(IN)  :: EWOI(:) !! ESTIMATED WIND-SEA ENERGY FROM MEASUREMENTS.
REAL,    INTENT(OUT) :: EWA(:)  !! ANALISED WIND-SEA ENERGY.
REAL,    INTENT(IN)  :: T(:)    !! DURATION.
REAL,    INTENT(OUT) :: FMWA(:) !! ANALYSED WINDSEA MEAN FREQUENCY.

! ---------------------------------------------------------------------------- !
!  
!      LOCAL VARIABLES

REAL, PARAMETER :: EO = 955.          !! PARAMETERS OF THE ENERGY GROWTH CURVE 
REAL, PARAMETER :: A = 6.02*10.**(-5) !! ESTAR = EO*TANH(A*TSTAR**B).
REAL, PARAMETER :: B = 0.695

INTEGER :: IJ, ITER, NNW
REAL    :: UMODEL, USTANAL, EWSR, TX, TSTAR, XX, TT, RKSI, SIGNR, SIGNM2, SIGNM1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALIZE OUTPUT FIELDS.                                             !
!        -------------------------                                             !

NNW = 0
USA  =   -99.
EWA  =  -999.
FMWA = -9999.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LOOP OVER GRID POINTS.                                                !
!        ----------------------                                                !

DO IJ = 1,SIZE(EWFG)
   IF (EWFG(IJ).LE.0.01) CYCLE

!     2.1  DETERMINE FRICTION VELOCITY FROM DURATION GROWTH CURVE              !
!          IF WINDSEA IS PRESENT.                                              !
!           ------------------------------------------------------             !

   UMODEL = USMO(IJ)
   USTANAL = UMODEL
   EWSR = EWOI(IJ)
   TX = T(IJ)
   DO ITER = 1,5
      TSTAR = G*TX/USTANAL
      XX = A*TSTAR**B
      IF (XX.GE.2.) THEN
         USTANAL = (EWSR*G**2/EO)**(0.25)
      ELSE
         TT = TANH(XX)
         RKSI = EWSR*G**2/USTANAL**4/EO
         USTANAL = USTANAL*(1.-(RKSI-TT)/(B*XX*(1.-TT**2)-4.*TT))
      END IF
   END DO

!     2.2  IF SIGNS OF (USA-USMO) AND (USOI-USMO) DIFFER, THEN                 !
!          THE MODEL IS WRONG ABOUT THE CORRECT RATIO BETWEEN                  !
!          WINDSEA AND SWELL. THEREFORE NO ASSIMILATION IS DONE.               !
!          -----------------------------------------------------               !

   IF (USOI(IJ).GT.0.) THEN
      SIGNR  = SIGN(1.,(USTANAL-UMODEL))
      SIGNM2 = SIGN(1.,(0.8*USOI(IJ)-UMODEL))
      SIGNM1 = SIGN(1.,(1.2*USOI(IJ)-UMODEL))
      IF (SIGNR.NE.SIGNM1.AND.SIGNR.NE.SIGNM2) THEN
         NNW = NNW+1
      ELSE
         USA(IJ) = USTANAL
         EWA(IJ) = EWSR
         FMWA(IJ) = YNU(EWSR*G**2/USTANAL**4)*G/USTANAL
      END IF
   ELSE
      USA(IJ) = USTANAL
      EWA(IJ) = EWSR
      FMWA(IJ) = YNU(EWSR*G**2/USTANAL**4)*G/USTANAL
   END IF
END DO

IF (ITEST.GE.4) THEN
   WRITE(IU06,*) '          SUB. FUSTAR: WINDSEA / USTAR ASSIMILATION ',       &
&                                                 'FAILED AT ', NNW, ' POINTS'
END IF

END SUBROUTINE FUSTAR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FWSEA (F, EMEAN, EWFG, THEW, FMWFG, ETOI, USMO, THMO)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      FWSEA - ANALYSE THE MODEL SPECTRA TO PROVIDE THE WIND-SEA PART          !
!              OF THE SPECTRUM, ITS ENERGY, ITS MEAN FREQUENCY.                !
!                                                                              !
!      P.LIONELLO     ECMWF         FEBRUARY 1989                              !
!        (MODIFICATION OF A CODE BY P.JANSSEN                                  !
!         AND P.LIONELLO - SUMMER '87 - )                                      !
!      J.BIDLOT       ECMWF         JANUARY 1997 : CORRECTION FOR INITIAL      !
!                                                  WAVE DIRECTION NOT BEING 0  !
!                                                                              !
!      METHOD.                                                                 !
!      -------                                                                 !
!                                                                              !
!      A PEAK IN THE SPECTRUM IS SOUGHT AROUND THE WIND DIRECTION.             !
!      IF A PEAK IS FOUND ,ITS FREQUENCY IS COMPARED TO THE P.M.               !
!      FREQUENCY TO ESTABLISH IF IT IS WINDSEA.                                !
!      THE WHOLE MODEL SPECTRUM IN THE REGION AROUND THE PEAK IS ASSUMED       !
!      TO BE WINDSEA. IN THE REMAINIG PART OF THE SPECTRUM THE MINIMUM         !
!      OF THE MODEL SPECTRUM AND THE CORRESPONDING JONSWAP SPECTRUM IS         !
!      TAKEN AS WINDSEA                                                        !
!      WIND SEA ENERGU AND MEAN FREQUENCY ARE COMPUTED BY THE EXTERNAL         !
!      WSMFEN                                                                  !
!                                                                              !
!     F        :     MODEL SPECTRUM (FIRST GUESS)                              !
!     FSEA     :     WIND-SEA PART OF THE MODEL SPECTRUM  (FIRST GUESS)        !
!     EWFG     :     FIRST GUESS WIND-SEA ENERGY                               !
!     THEW     :     DIRECTION OF THE WIND-SEA PEAK (FIRST GUESS)              !
!     THMO     :     WIND DIRECTION                                            !
!     USMO     :     USTAR (FIRST GUESS)                                       !
!     FMWFG    :     WINDSEA MEAN FREQUENCY (FIRST GUESS)                      !
!     ETOI     :     TOTAL ENERGY FROM O.I. (ANALYSIS)                         !
!                                                                              !
!      EXTERNALS.                                                              !
!      ----------                                                              !
!                                                                              !
!        F4SPEC  - WINDSEA SPECTRUM FOR ASSIMILATION.                          !
!        WSMFEN  - MEAN FREQUENCY AND ENERGY OF THE WINDSEA SPECTRUM           !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.

REAL,    INTENT(IN)  :: F(:,:,:) !! MODEL SPECTRUM (FIRST GUESS).
REAL,    INTENT(IN)  :: EMEAN(:) !! FIRST GUESS ENERGY.
REAL,    INTENT(OUT) :: EWFG(:)  !! FIRST GUESS WIND-SEA ENERGY.
REAL,    INTENT(OUT) :: THEW(:)  !! DIRECTION OF THE WIND-SEA PEAK (FIRST GUESS)
REAL,    INTENT(IN)  :: THMO(:)  !! WIND DIRECTION  .
REAL,    INTENT(IN)  :: USMO(:)  !! USTAR (FIRST GUESS).
REAL,    INTENT(OUT) :: FMWFG(:) !! WINDSEA MEAN FREQUENCY (FIRST GUESS).
REAL,    INTENT(IN)  :: ETOI(:)  !! TOTAL ENERGY FROM O.I. (ANALYSIS).


! ---------------------------------------------------------------------------- !
! 
!      LOCAL VARIABLES

INTEGER :: KL, ML
INTEGER :: NWP, IJ, IDW, MP, KP, II, K, M, I, KPP, KP1, KM1, JP

REAL :: GAM, SAO, X
REAL :: FR1, THEWI, FSMAX, FMODXX, FMOD, FMOD3, FMOD1, X1, X2, X3, X4, X5, X6
REAL :: FPWW, XNU, THE, DDIR, DTHE, SPINT, ALPHA, GAMMA, SA, FSH1, XLL

REAL :: FSEA(SIZE(F,2), SIZE(F,3))   !! WIND-SEA PART  (FIRST GUESS).
REAL :: FJONS(SIZE(F,2), SIZE(F,3))


REAL, PARAMETER :: CGAM = 0.00608
REAL, PARAMETER :: ESTARM = 955., A = 6.02*10.**(-5), B = 0.695
REAL, PARAMETER :: EGAM = 1.439, ESAO = 2.338
REAL, PARAMETER :: XL11 = 0.0953101

!     CGAM       : P.M. FREQUENCY IN THE WAM MODEL                             !
!     ESTARM,A,B : PARAMETERS OF THE ENERGY GROWTH CURVE                       !
!     EGAM,ESAO  : PARAMETERS IN THE WINDSEA SPECTRUM SHAPE                    !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INLINE FUNCTIONS.                                                        !
!     -----------------                                                        !
!                                                                              !
!     OVERSHOOT AS FUNCTION OF DIMENSIONLESS FREQUENCY.                        !

GAM(X) = MAX (1.,1+2.70*(1-(CGAM/X)**EGAM))

!     PEAK EXTENSION AS FUNCTION OF DIMENSIONLESS ENERGY.                      !

SAO(X) = MAX (.02,.02+.16*(1-(CGAM/X)**ESAO))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. PARAMETERS AND INITIALIZATION.                                        !
!        ------------------------------                                        !

KL = SIZE(F,2)
ML = SIZE(F,3)

FR1 = FR(1)
NWP = 0
THEW  = 999999999.
EWFG  = -999.
FMWFG = -9999.

! ---------------------------------------------------------------------------- !
!                                                                              !
!                                                                              !
!     2. LOOP ON THE POINTS INSIDE THE BLOCK.                                  !
!        ------------------------------------                                  !

POINT: DO IJ = 1,SIZE(F,1)

!     2.1 SKIP THE LAND POINTS.                                                !
!         ---------------------                                                !

   IF (EMEAN(IJ).LE.0.01.OR.ETOI(IJ).LE.0.01) CYCLE POINT

!*    2.2 INITIALIZE SPECTRUM BY ZERO.                                         !
!         ----------------------------                                         !

   FSEA(1:KL,1:ML) = 0.

!     2.3 SEARCH A PEAK OF THE SPECTRUM AROUND THE WIND DIRECTION.             !
!         --------------------------------------------------------             !

   THEWI = THMO(IJ)
   IDW = NINT(MOD(THEWI-TH(1),ZPI)/DELTH)+1
   FSMAX = 0.
   MP = 1
   KP = 1
   DO II = IDW-1,IDW+1
      K = II
      IF (K.GT.KL) K = K-KL
      IF (K.LE.0)    K = K+KL
      DO M = 1,ML
         FMODXX = F(IJ,K,M)
         IF (FMODXX.GT.FSMAX) THEN
            FSMAX = FMODXX
            MP = M
            KP = K
         END IF
      END DO
   END DO

!     2.4 IF A PEAK HAS BEEN FOUND THE PEAK FREQUENCY IS EVALUATED.            !
!         ---------------------------------------------------------            !

   IF (MP.LE.1 .OR. MP.GE.ML) CYCLE POINT

      FMOD=0.
      FMOD3=0.
      FMOD1=0.

!     2.4.1 THE SPECTRUM IS INTEGRATED AROUND THE PEAK TO OBTAIN               !
!           A BETTER ESTIMATE OF THE PEAK IN THE 1-DIM SPECTRUM.               !
!           ----------------------------------------------------               !

      DO I = -1,1
         KPP=KP+I
         IF (KPP.EQ.KL+1) KPP = 1
         IF (KPP.EQ.0)      KPP = KL
         FMOD  = F(IJ,KPP,MP  )+FMOD
         FMOD3 = F(IJ,KPP,MP+1)+FMOD3
         FMOD1 = F(IJ,KPP,MP-1)+FMOD1
      END DO

!     2.4.2 CHECK IF THERE IS A PEAK IN THE 1-DIM SPECTRUM.                    !
!           IF NOT THEN THE 2D SPECTRUM IS USED TO COMPUTE                     !
!           THE PEAK FREQUENCY.                                                !
!           -----------------------------------------------                    !

      IF (FMOD.LE.FMOD3 .OR. FMOD.LE.FMOD1) THEN
         FMOD  = F(IJ,KP,MP  )
         FMOD3 = F(IJ,KP,MP+1)
         FMOD1 = F(IJ,KP,MP-1)
      END IF
      X1 = FR(MP)**2-FR(MP+1)**2
      X2 = FR(MP)-FR(MP+1)
      X3 = FMOD-FMOD3
      X4 = FR(MP-1)**2-FR(MP+1)**2
      X5 = FR(MP-1)-FR(MP+1)
      X6 = FMOD1-FMOD3
      FPWW = 0.5*(X3*X4-X1*X6)/(X3*X5-X2*X6)
      XNU = USMO(IJ)*FPWW/G

!     2.4.3 EVALUATE PEAK DIRECTION                                            !
!           IF PEAK FREQUENCY IS GE PM FREQUENCY.                              !
!           -------------------------------------                              !

      IF (XNU.LT.CGAM) CYCLE POINT
         THE = TH(KP)
         KP1 = KP+1
         KM1 = KP-1
         IF (KP1.GT.KL) KP1=KP1-KL
         IF (KM1.LT.1) KM1=KM1+KL
         FMOD  = F(IJ,KP ,MP)
         FMOD3 = F(IJ,KP1,MP)
         FMOD1 = F(IJ,KM1,MP)
         DDIR = DELTH*(FMOD3-FMOD1)/(2.*FMOD-FMOD1-FMOD3)
         THE = THE+DDIR/2.

         DTHE = MOD(ABS(THE-THEWI),ZPI)
         IF (DTHE.GT.PI) DTHE = ZPI-DTHE
         IF (DTHE.GT.DELTH) CYCLE POINT  
     
         NWP = NWP+1
!                                                                              !
!     2.4.4 BUILD THE CORRESPONDING WIND-SEA SPECTRUM                          !
!           A PARAMETRIZATION OF THE MODEL SPECTRUM HAVING A "TOBA"            !
!           F**(-4) TAIL IS BUILT .                                            !

         SPINT = 0.
!                                                                              !
!             TOBA CONSTANT ACCORDING TO THE WAM MODEL                         !
!                                                                              !
         ALPHA = .145
         GAMMA = GAM(XNU)
         SA = SAO(XNU)
         CALL F4SPEC (ALPHA, FPWW, GAMMA, SA, THE, USMO(IJ), SPINT, FJONS)
         FSEA(1:KL,1:ML) = MIN(F(IJ,1:KL,1:ML),FJONS(1:KL,1:ML))

!     2.4.5 AROUND THE WIND SEA PEAK THE WHOLE MODEL SPECTRUM IS ASSUMED       !
!           TO BE WINDSEA.                                                     !
!           ------------------------------------------------------------       !

         IDW = NINT(MOD(THE-TH(1),ZPI)/DELTH)+1
         DO II = IDW-2,IDW+2
            K = MOD(II,KL)
            IF (K.LE.0)  K=K+KL
            FSH1 = COS(THE-TH(K))
            IF (FSH1.GT.0.001) THEN
               XLL = 0.7*FPWW/(FSH1*FR1)
               JP = NINT(1.+ALOG(XLL)/XL11)
               IF (JP.LE.ML) THEN
                  JP = MAX (JP,1)
                  FSEA(K,JP:ML) = F(IJ,K,JP:ML)
               END IF
            END IF
         END DO

!     2.4.6 THE WIND DIRECTION IS STORED.                                      !
!           -----------------------------                                      !

         THEW(IJ) = THE
!                                                                              !
!     2.4.7 COMPUTATION OF THE WINDSEA ENERGY WINDSEA MEAN FREQUENCY.          !
!           ---------------------------------------------------------          !

         CALL WSMFEN (FSEA, EWFG(IJ), FMWFG(IJ), USMO(IJ))

END DO POINT

IF (ITEST.GE.4) THEN
   WRITE(IU06,*) '          SUB. FWSEA: WIND SEA FOUND AT ',NWP, ' POINTS'
END IF

END SUBROUTINE FWSEA

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE F4SPEC (ALFA, FM, GAMMA, SA, THETAQ, USTAR, SPINT, SPT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      F4SPEC* - PRODUCE A WINDSEA SPECTRUM FOR ASSIMILATION.                  !
!                                                                              !
!      P. LIONELLO    ECMWF        FEBRUARY 1990                               !
!      (DERIVED FROM JSPEC OF PETER JANSSEN  1987.)                            !
!                                                                              !
!      METHOD.                                                                 !
!      -------                                                                 !
!                                                                              !
!      THE ONE DIMENSIONAL PARAMETRIC SPECTRUM IS EVALUATED AND                !
!      DISTRIBUTED AROUND THE MEAN SEA DIRECTION.                              !
!      THE TOTAL ENERGY IS COMPUTED.                                           !
!                                                                              !
!      EXTERNALS.                                                              !
!      ----------                                                              !
!                                                                              !
!         *COS2_SPR*  -  SUB TO COMPUTE COS**2 SPREADING FUNCTION.             !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!     INTERFACE VARIABLES.

IMPLICIT NONE

REAL, INTENT(IN)  :: ALFA        !! TOBAS' CONSTANT.
REAL, INTENT(IN)  :: FM          !! PEAK FREQUENCY.
REAL, INTENT(IN)  :: GAMMA       !! OVERSHOOT PARAMETER.
REAL, INTENT(IN)  :: SA          !! WIDTH PARAMETER.
REAL, INTENT(IN)  :: THETAQ      !! MEAN WINDSEA DIRECTION.
REAL, INTENT(IN)  :: USTAR       !! FRICTION VELOCITY.
REAL, INTENT(OUT) :: SPINT       !! TOTAL ENERGY.
REAL, INTENT(OUT) :: SPT(:,:)  !! 2-D SPECTRUM.

! ---------------------------------------------------------------------------- !
!
!      LOCAL VARIABLES

INTEGER :: ML, KL
INTEGER :: M, K, ML2, KL2
REAL    :: FAK, FRH, EARG, FJON, FMPF, FJONH
REAL    :: DELFR(SIZE(SPT,2))    !! FREQUENCY INCREMENTS.
REAL    :: ET(SIZE(SPT,2))       !! 1-D TOBA SPECTRUM.
REAL    :: ST(SIZE(SPT,1))       !! SPREADING FUNCTION. 

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE FREQUENCY INTERVALLS.                                         !
!        -----------------------------                                         !

KL = SIZE(SPT,1)
ML = SIZE(SPT,2)
DELFR(1:ML) = FR(1:ML)*0.1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE TOBA SPECTRUM.                                                !
!        ----------------------                                                !

FAK = ALFA*G*USTAR/(ZPI**3)
DO  M = 1,ML
   FRH = FR(M)
   EARG = (.5*((FRH-FM)/(SA*FM))**2)
   FJON = GAMMA**EXP(-EARG)
   FMPF = 4./6.*(FM/FRH)**6
   FJONH = EXP(-FMPF)
   ET(M) = FAK*FJONH*FJON/FRH**4
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTATION OF SPREADING FUNCTION.                                    !
!        ----------------------------------                                    !

CALL COS2_SPR (TH, THETAQ, ST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTATION OF 2-D SPECTRUM.                                          !
!        ---------------------------                                           !

DO M = 1,ML
   DO K = 1,KL
      SPT(K,M) = ET(M) * ST(K)
   END DO
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. COMPUTE TOTAL ENERGY.                                                 !
!        --------------------                                                  !
!                                                                              !
!     5.1 INTEGRATE OVER DIRECTIONS.                                           !
!         --------------------------                                           !

KL2 = KL-2
SPINT = 0.
DO M = 1,ML
   ET(M) = SPT(1,M)
   DO K = 2,KL2,2
      ET(M) = ET(M)+4.*SPT(K,M)+2.*SPT(K+1,M)
   END DO
   ET(M) = (ET(M)+4.*SPT(KL,M)+SPT(1,M))*DELTH/3.
END DO

!*    5.2 INTEGRATE 1-D SPECTRUM OVER FREQUENCIES.                             !
!         ----------------------------------------                             !

SPINT = ET(1)*DELFR(1)
ML2 = ML-2
DO M = 2,ML2,2
   SPINT = SPINT+(4.*ET(M)*DELFR(M)+2.*ET(M+1)*DELFR(M+1))
END DO

IF (2*(ML/2).EQ.ML) THEN
   SPINT = (SPINT+4.*ET(ML)*DELFR(ML))/3.
ELSE
   SPINT = (SPINT+4.*ET(ML-1)*DELFR(ML-1) + ET(ML)*DELFR(ML))/3.
ENDIF

END SUBROUTINE F4SPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE UPDATE (F, EMEAN, ETOI, USOI, USMO, THMO)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      UPDATE - ANALYSE THE WAVE SPECTRUM, PRODUCING A CONSISTENT              !
!               UPDATE OF WAVE AND WIND FIELD.                                 !
!                                                                              !
!     PIERO LIONELLO      ECMWF     JUNE 1990                                  !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO ANALYSE THE WAVE SPECTRUM, PRODUCING A CONSISTENT                   !
!       UPDATE OF WAVE AND WIND FIELD.                                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!      EXTERNALS.                                                              !
!      ----------                                                              !
!                                                                              !
!          FDUR    -  FRICTION VELOCITY AND WINDSEA MEAN PERIOD.               !
!          FUSTAR  -  WINDSEA DURATION.                                        !
!          FWSEA   -  PROVIDE THE WIND-SEA PART OF THE SPECTRUM.               !
!          F4SPEC  -  RODUCE A WINDSEA SPECTRUM FOR ASSIMILATION.              !
!          WSMFEN  -  MEAN FREQUENCY AND ENERGY OF WINDSEA SPECTRUM.           !
!          UPSPEC  -  MODIFY THE SPECTRUM BY STRETCHING AND SCALING.           !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.

REAL,    INTENT(INOUT) :: F(:,:,:) !! SPECTRUM (FIRST GUESS IN INPUT ;
                                   !!           ANALYSIS IN OUTPUT ).
REAL,    INTENT(IN)    :: EMEAN(:) !! WAVE ENERGY (FROM FIRST GUESS).
REAL,    INTENT(inout) :: ETOI(:)  !! WAVE ENERGY (FROM O.I.).
REAL,    INTENT(IN)    :: USOI(:)  !! USTAR ( FROM O.I.).
REAL,    INTENT(INOUT) :: USMO(:)  !! USTAR ( FIRST GUESS IN INPUT ;
                                   !!         ANALYSIS IN OUTPUT ).
REAL,    INTENT(IN)    :: THMO(:)  !! WIND DIRECTION.

! ---------------------------------------------------------------------------- !
! 
!      LOCAL VARIABLES

REAL :: EWFG(SIZE(F,1))  !! WINDSEA ENERGY (FIRST GUESS).
REAL :: EWOI(SIZE(F,1))  !! WINDSEA ENERGY (MEASUREMENT).
REAL :: EWA(SIZE(F,1))   !! WINDSEA ENERGY (ANALYSIS).
REAL :: FMWFG(SIZE(F,1)) !! MEAN FREQUENCY OF THE WINDSEA (FIRST GUESS).
REAL :: FMWA(SIZE(F,1))  !! MEAN FREQUENCY OF THE WINDSEA (ANALYSIS).
REAL :: THEW(SIZE(F,1))  !! DIRECTION OF THE WINDSEA PEAK (FIRST GUESS).
REAL :: TDUR(SIZE(F,1))  !! WINDSEA DURATION.
REAL :: USA(SIZE(F,1))   !! ESTIMATE OF USTAR FROM THE FIRST GUESS ENERGY.
                         !! AND DURATION (I.E. ANALYSED USTAR).

! ---------------------------------------------------------------------------- !
! 
!     1. FIND THE WINDSEA.                                                     !
!        -----------------                                                     !

CALL FWSEA (F, EMEAN, EWFG, THEW, FMWFG, ETOI, USMO, THMO)

! ---------------------------------------------------------------------------- !
! 
!     2. CORRECT OVERESTIMATE OF FIRST GUESS WINDSEA ENERGY.                   !
!        ---------------------------------------------------                   !

EWFG = MIN(EWFG,EMEAN)

IF (ITEST.GE.4) WRITE (IU06,*) '      SUB. UPDATE: ',                          &
&                              ' WINDSEA ENERGY FIXED IN SUB. FWSEA'

! ---------------------------------------------------------------------------- !
! 
!     3. FIND THE DURATION.                                                    !
!        ------------------                                                    !

CALL FDUR (EWFG, TDUR, USMO)

IF (ITEST.GE.4) WRITE (IU06,*) '      SUB. UPDATE: ',                          &
&                              ' WINDSEA DURATION FIXED IN SUB. FDUR'

! ---------------------------------------------------------------------------- !
! 
!     4. THE RATIO WINDSEA/SWELL IS ASSUMED CORRECT.                           !
!        -------------------------------------------                           !

WHERE (EMEAN.GT.0.) EWOI = EWFG/EMEAN*ETOI

! ---------------------------------------------------------------------------- !
! 
!     5. FIND THE NEW USTAR AND THE NEW WINDSEA MEAN FREQUENCY.                !
!        ------------------------------------------------------                !

CALL FUSTAR (USMO, USA, USOI, EWFG, EWOI, EWA, TDUR, FMWA)

WHERE (USA.GT.0.) USMO = USA

IF (ITEST.GE.4) WRITE (IU06,*) '      SUB. UPDATE: ',                          &
&                              ' USTAR COMPUTED IN SUB. FUSTAR'

! ---------------------------------------------------------------------------- !
! 
!     6. THE ASSIMILATION IS ADJUSTED ON THE DOMINANT PART OF THE SPECTRUM.    !
!        ------------------------------------------------------------------    !
!                                                                              !
!        THE FOLLOWING RETURN CODE ARE POSSIBLE :                              !
!          THERE IS NO WINDSEA ->                                              !
!                 EWFG  =  -999.                                               !
!                 FMWFG = -9999.                                               !
!          THE RATIO WINDSEA/SWELL IS NOT CONSISTENT WITH USOI ->              !
!                 EWA   =  -999.                                               !
!                 FMWA  = -9999.                                               !

WHERE (EWA.GT.0. .AND. EWA.LT.0.5*ETOI)
      EWFG  =  -999.
      FMWFG = -9999.
ENDWHERE

!     IF THERE WERE WINDSEA BUT THE RATIO WINDSEA SWELL                        !
!     WAS WRONG, DO NOT ASSIMILATE.                                            !

WHERE (EWFG.GT.0. .AND. EWA.LT.0.) ETOI = -999.

! ---------------------------------------------------------------------------- !
! 
!     7. UPDATE THE SPECTRA.                                                   !
!        -------------------                                                   !

CALL UPSPEC (F, EMEAN, ETOI, FMWFG, FMWA)

IF (ITEST.GE.4) WRITE (IU06,*) '      SUB. UPDATE: ',                          &
&                              ' NEW SPECTRA COMPUTED IN SUB. UPSPEC'

END SUBROUTINE UPDATE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE UPSPEC (F, ETFG, ETA, FMWFG, FMWA)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     UPSPEC - MODIFY THE SPECTRUM BY STRETCHING AND SCALING.                  !
!                                                                              !
!     P.LIONELLO      ECMWF         APRIL 1990                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!         TO MODIFY THE SPECTRUM BY STRETCHING AND SCALING.                    !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!       A NEW SPECTRUM IN THE FORM                                             !
!          F   (IPOINT,F,K) =  A F   (IPOINT,BF,K)                             !
!           NEW                   OLD                                          !
!       IS BUILD.                                                              !
!                                                                              !
!       IF THERE IS MAINLY SWELL THE SPECTRUM IS UPDATED USING THE             !
!       AVERAGE STEEPNESS CRITERIUM. TO CONSERVE EXACTLY THE                   !
!       STEEPNESS AND CHANGE THE ENERGY THE CONSTANT A AND B                   !
!       MUST BE GIVEN BY                                                       !
!                                                                              !
!                A = (ETA/ETFG)**1.25  B = (ETA/ETFG)**.25                     !
!                                                                              !
!       A SMALL CORRECTION , WHICH IS SUGGESTED BY THE MODEL                   !
!       DECAY CURVE ,IS ACTUALLY INTRODUCED ACCORDING TO                       !
!                                                                              !
!            DELTA = 1 - .006 * ( HNEW - HOLD )                                !
!                                                                              !
!                A = DELTA * (ETA/ETFG)**1.25                                  !
!                B = DELTA * (ETA/ETFG)**.25                                   !
!                                                                              !
!       IF THERE IS MAINLY WINDSEA, THE STRETCHING CONSTANT B IS               !
!       COMPUTED TO PRODUCE IN THE ANALYSED  SPECTRUM THE ANALYSED             !
!       MEAN FREQUENCY DERIVED BY THE MODEL GROWTH CURVE                       !
!       ( AS CARRIED OUT IN THE SUBROUTINE FUSTAR ). IN THIS CASE              !
!                                                                              !
!                  B = FMWFG/FMWA                                              !
!                  A = (ETA/ETFG)*B                                            !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.

REAL,    INTENT(INOUT) :: F(:,:,:) !! INPUT: THE OLD SWELL SPECTRUM.
                                   !! OUTPUT: THE ANALISED SPECTRUM.
REAL,    INTENT(IN)    :: ETFG(:)  !! FIRST GUESS TOTAL ENERGY.
REAL,    INTENT(IN)    :: ETA(:)   !! ANALISED TOTAL ENERGY.
REAL,    INTENT(IN)    :: FMWFG(:) !! MEAN FREQUENCY OF THE WINDSEA 
                                   !! FROM THEFIRST GUESS SPECTRUM.
REAL,    INTENT(IN)    :: FMWA(:)  !! ANALYSED MEAN FREQUENCY OF WINDSEA.

! ---------------------------------------------------------------------------- !
!
!      LOCAL VARIABLES

REAL, PARAMETER :: XL11 = 0.0953101    !! ALOG(1.1)

INTEGER :: KL, ML
INTEGER :: IJ, M, M1, M2, K
REAL    :: FR1, FNEW, HNEW, HOLD, FOLD, XDELTA, XR, XB, FU, DFRE, F1, F2, DE
REAL    :: FTEMP(SIZE(F,2),SIZE(F,3))

! ---------------------------------------------------------------------------- !

KL = SIZE(F,2)
ML = SIZE(F,3)
FR1 = FR(1)

! ---------------------------------------------------------------------------- !
! 
!     1. LOOP OVER GRID POINTS.                                                !
!        ---------------------                                                 !

POINT: DO IJ = 1,SIZE(F,1)

!     1.1 SKIP LAND POINTS AND POINTS WHERE THERE ARE NO RELIABLE DATA.        !
!         -------------------------------------------------------------        !

   IF (ETFG(IJ).LE.0.001.OR.ETA(IJ).LE.0.001) CYCLE POINT
   FNEW = FMWA(IJ)
   HNEW = 4.*SQRT(ETA(IJ))
   HOLD = 4.*SQRT(ETFG(IJ))
   FOLD = FMWFG(IJ)

!     1.2 COMPUTE SCALING AND STRETCHING FACTORS.                              !
!         ---------------------------------------                              !

   IF (FOLD.LE.0.) THEN

!     1.2.1 THE SPECTRUM IS MAINLY SWELL (THE WINDSEA MEAN FREQUENCY           !
!           OF THE FIRST GUESS SPECTRUM IS NEGATIVE).                          !
!           --------------------------------------------------------           !

      XDELTA = 1.-0.006*(HNEW-HOLD)
      XR = XDELTA*(HNEW/HOLD)**2.5
      XB = XDELTA* SQRT(HNEW/HOLD)
   ELSE
!                                                                              !
!     1.2.2 THE SPECTRUM IS MAINLY WINDSEA.                                    !
!           -------------------------------                                    !

      XR = (HNEW/HOLD)**2*FOLD/FNEW
      XB = FOLD/FNEW
   END IF

!     1.3 LOOP OVER NEW FREQUENCIES.                                           !
!         --------------------------                                           !

   DO M = 1,ML
      FU = FR(M)*XB
      M1 = IFIX(LOG(FU/FR1)/XL11)+1
      IF (M1.GE.ML.OR.M1.LT.1) THEN

!     1.3.1 LOOP OVER DIRECTIONS FOR FREQUENCIES GETTING ENERGY FORM           !
!           FREQUENCIES OUT OF RANGE.                                          !
!           ---------------------------------------------------------          !

         FTEMP(1:KL,M) = 0.

      ELSE
!                                                                              !
!     1.3.2 LOOP OVER DIRECTIONS FOR FREQUENCIES GETTING ENERGY.               !
!           ----------------------------------------------------               !
      
         M2 = M1+1
         DFRE = FR(M1)*.1
         DO K=1,KL
            F1 = F(IJ,K,M1)
            F2 = F(IJ,K,M2)
            DE = (F2-F1)/DFRE * (FU-FR(M1))
            FTEMP(K,M) = XR*MAX(0.,F1+DE)
         END DO

      END IF

   END DO  !! BRANCH BACK FOR NEXT FREQUENCY.

!     1.4 THE UPDATED SPECTRUM IS STORED.                                      !
!         -------------------------------                                      !

   F(IJ,1:KL,1:ML) = FTEMP(1:KL,1:ML)

END DO POINT         !! BRANCH BACK FOR NEXT GRIDPOINT.

END SUBROUTINE UPSPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WSMFEN (FSEA, EW, FM, USTT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      WSMFEN* - COMPUTES MEAN FREQUENCY AND ENERGY OF THE WINDSEA             !
!                SPECTRUM.                                                     !
!                                                                              !
!     P.LIONELLO     ECMWF       APRIL 1990                                    !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTES MEAN FREQUENCY AND ENERGY OF THE WINDSEA SPECTRUM.            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       ENERGY AND MEAN FREQUENCY ARE FIRST COMPUTED BY INTEGRATION            !
!       OF THE WINDSEA PART OF THE SPECTRUM:                                   !
!       THIS IMPLIES AN UNDERESTIMATION BOTH OF ENERGY AND MEAN                !
!       FREQUENCY. THE UNDERESTIMATED VALUES ARE USED IN THE GROWTH            !
!       CURVE PRODUCING OVERSTIMATES OF BOTH ENERGY AND MEAN FREQUENCY.        !
!       THE AVERAGE OF THE TWO ESTIMATES IS TAKEN TO PROVIDE A BEST            !
!       ESTIMATE.                                                              !
!                                                                              !
!     EXTERNAL.                                                                !
!     ---------                                                                !
!                                                                              !
!      EN  - DIMENSIONLESS ENERGY AS FUNCTION OF DIMENSIONLESS MEAN FREQUENCY. !
!      YNU - DIMENSIONLESS FREQUENCY AS FUNCTION OF DIMENSIONLESS MEAN ENERGY. !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.

IMPLICIT NONE

REAL, INTENT(IN)  :: FSEA(:,:)      !! WINDSEA PART OF THE WAVE SPECTRUM.
REAL, INTENT(OUT) :: EW             !! WINDSEA ENERGY.
REAL, INTENT(OUT) :: FM             !! MEAN FREQUENCY OF THE WINDSEA.
REAL, INTENT(IN)  :: USTT           !! FRICTION VELOCITY.

! ---------------------------------------------------------------------------- !
! 
!      LOCAL VARIABLES

REAL, PARAMETER :: EPSMIN=0.1E-32

INTEGER :: ML, KL
INTEGER :: M, K
REAL    :: DELT25, TEMP, SPFB, XTEMP, SPINTDI, FREQDI

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATING THE WINDSEA PART OF THE SPECTRUM.                         !
!        (UNDERESTIMATE IS IMPLIED)                                            !
!        ---------------------------------------------                         !


KL = SIZE(FSEA,1)
ML = SIZE(FSEA,2)

!     1.1 COMPUTATION OF THE WIND-SEA ENERGY.                                  !
!         -----------------------------------                                  !

EW = EPSMIN
DO M = 1,ML
   TEMP = 0.
   DO K = 1,KL
      TEMP = TEMP+FSEA(K,M)
   END DO
   EW = EW+TEMP*DFIM(M)
END DO

DELT25 = FR(ML)/3.*DELTH
EW = EW+DELT25*TEMP

!     1.2 COMPUTATION OF THE MEAN FREQUENCY.                                   !
!         ----------------------------------                                   !

FM = EPSMIN
DO M = 1,ML
   SPFB = 0.
   DO K = 1,KL
      SPFB = SPFB+FSEA(K,M)
   END DO
   SPFB = SPFB*DELTH
   FM = FM+SPFB*FR(M)**2.*(1.1-1./1.1)/2.
END DO

FM = FM+SPFB*FR(ML)**2./3.
FM = FM/EW

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. ESTIMATES FROM MODEL RELATIONS.                                       !
!        -------------------------------                                       !
!                                                                              !
!     2.1 ENERGY IS DERIVED FROM UNDERESTIMATED MEAN FREQUENCY.                !
!         -----------------------------------------------------                !

XTEMP = FM*USTT/G
IF (XTEMP.LE.0.) WRITE(*,*) FM,USTT,EW
SPINTDI = USTT**4/G**2 * EN(FM*USTT/G)

!     2.2 MEAN FREQUENCY IS DERIVED FROM THE UNDERESTIMATED ENERGY.            !
!         ---------------------------------------------------------            !

FREQDI = G/USTT * YNU(EW*G**2/USTT**4)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. FINAL ESTIMATE.                                                       !
!        ---------------                                                       !
!                                                                              !
!     3.1 AVERAGING THE ENERGY ESTIMATES.                                      !
!         -------------------------------                                      !

EW = (EW+SPINTDI)*0.5

!     3.2 AVERAGING THE MEAN FREQUENCY ESTIMATES.                              !
!         ---------------------------------------                              !

FM = (FM+FREQDI)*0.5

END SUBROUTINE WSMFEN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

FUNCTION EN (X) RESULT(Z)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     EN - DIMENSIONLESS ENERGY AS FUNCTION OF DIMENSIONLESS FREQUENCY.        !
!                                                                              !
!     H. GUNTHER     GKSS      FEBRUARY 2004                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL, INTENT(IN) :: X    !! DIMENSIONLESS FREQUENCY.
REAL             :: Z    !! DIMENSIONLESS ENERGY

REAL, PARAMETER :: AF = 0.000168  !! PARAMETER OF THE ENERGY - MEAN FREQUENCY RELATION.
REAL, PARAMETER :: BF = -3.27     !! PARAMETER OF THE ENERGY - MEAN FREQUENCY RELATION.

! ---------------------------------------------------------------------------- !

Z = AF*X**BF               !! ENERGY AS FUNCTION OF MEAN FREQUENCY.

END FUNCTION EN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

FUNCTION YNU (X) RESULT(Z)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     YNU - DIMENSIONLESS FREQUENCY AS FUNCTION OF DIMENSIONLESS ENERGY.       !
!                                                                              !
!     H. GUNTHER     GKSS      FEBRUARY 2004                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL, INTENT(IN) :: X    !! DIMENSIONLESS ENERGY
REAL             :: Z    !! DIMENSIONLESS FREQUENCY.

! ---------------------------------------------------------------------------- !
!                                                                              !
!      LOCAL VARIABLES                                                         !
!      ---------------                                                         !

REAL, PARAMETER :: AF = 0.000168  !! PARAMETER OF THE ENERGY - MEAN FREQUENCY RELATION.
REAL, PARAMETER :: BF = -3.27     !! PARAMETER OF THE ENERGY - MEAN FREQUENCY RELATION.

! ---------------------------------------------------------------------------- !

Z = (X/AF)**(1./BF)       !! MEAN FREQUENCY AS FUNCTION OF ENERGY.

END FUNCTION YNU

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SYMINV (A, NDIM, N, COND, V)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INVERT IN PLACE THE LOWER TRIANGLE OF A (I.E. A(I,J) I.GE.J)             !
!     A IS A SYMMETRIC POSITIVE DEFINITE MATRIX                                !
!     THE UPPER TRIANGLE OF A (IE. A(I,J) I.LT.J) IS NOT USED OR ALTERED       !
!                                                                              !
!     THIS VERSION IS OPTIMIZED FOR THE CDC FTN COMPILER.                      !
!     IT REQUIRES THE FUNCTION 'NUMARG' FROM  'ECLIB'                          !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)    :: NDIM
INTEGER, INTENT(IN)    :: N
REAL,    INTENT(OUT)   :: COND
REAL,    INTENT(INOUT) :: A(NDIM,NDIM)
REAL,    INTENT(OUT)   :: V(NDIM)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      LOCAL VARIABLES                                                         !
!      ---------------                                                         !

INTEGER :: IMX, J, I, K, NDIA
REAL    :: EPS, ZMX, X, S

! ---------------------------------------------------------------------------- !

IF(N.LT.1) RETURN
EPS = 0.

NDIA = NDIM+1
IMX = ISAMAX (N, A(1,1), NDIA)
ZMX = A(IMX,IMX)*FLOAT(N*N)

J = 1
IF (A(1,1).LE.EPS) THEN
   COND = -FLOAT(J)
   RETURN
END IF

SELECT CASE (N)

   CASE(1)
      V(1) = 1./A(1,1)
      A(1,1) = V(1)

   CASE (2)
      V(1) = 1./A(1,1)     
      X = A(2,1)*V(1)
      A(2,2) = A(2,2)-X*A(2,1)
      A(2,1) = X
      IF (A(2,2).LE.EPS) THEN
         COND = -FLOAT(J)
         RETURN
      END IF
      V(2) = 1./A(2,2)

      DO J = 1,2
         A(J,J) = V(J)
      END DO

      A(2,1) = -A(2,1)
      S = A(1,1)
      X = A(2,2)*A(2,1)
      S = S + A(2,2)*X
      A(2,1) = X
      A(1,1) = S

   CASE DEFAULT

!**** 1. ROOT FREE CHOLESKY DECOMPOSITION  A = L D L(TRANSPOSE)                !

      V(1) = 1./A(1,1)     
      X = A(2,1)*V(1)
      A(2,2) = A(2,2)-X*A(2,1)
      A(2,1) = X
      IF (A(2,2).LE.EPS) THEN
         COND = -FLOAT(J)
         RETURN
      END IF
      V(2) = 1./A(2,2)

      DO I = 3,N
         DO J = 3,I
            S = A(I,J-1)
            DO K = 3,J
               S = S-A(I,K-2)*A(J-1,K-2)
            END DO
            A(I,J-1) = S
         END DO

         S = A(I,I)
         DO J = 2,I
            X = A(I,J-1)*V(J-1)
            S = S - X*A(I,J-1)
            A(I,J-1) = X
         END DO
         A(I,I) = S

         IF (A(I,I).LE.EPS) THEN
            COND = -FLOAT(J)
            RETURN
         END IF
         V(I) = 1./A(I,I)
      END DO

!**** 2.  COPY INVERSE OF D WHICH HAS ALREADY BEEN CALCULATED.                 !

      DO J = 1,N
         A(J,J) = V(J)
      END DO

!**** 3.  INVERSION OF L                                                       !

      A(2,1) = -A(2,1)
      DO I = 2,N - 1
         DO J = 2,I
            A(I+1,J-1) = -A(I+1,J-1) - DOT_PRODUCT(A(I+1,J:I), A(J:I,J-1))
         END DO
         A(I+1,I) = -A(I+1,I)
      END DO

!**** 4.  INV A = INV L(TRANSPOSE) * INV D * INV L                             !

      DO J = 2,N
         S = A(J-1,J-1)
         DO I = J,N
            X = A(I,I)*A(I,J-1)
            S = S + A(I,J-1)*X
            A(I,J-1) = X
         END DO
         A(J-1,J-1) = S

         DO I = J,N - 1
            A(I,J-1) = A(I,J-1) + DOT_PRODUCT(A(I+1:N,I), A(I+1:N,J-1))
         END DO
      END DO

END SELECT

IMX = ISAMAX (N, A(1,1), NDIA)
COND = 1./ABS(A(IMX,IMX)*ZMX)

RETURN

CONTAINS

! ---------------------------------------------------------------------------- !

   FUNCTION ISAMAX (N, A, M) RESULT(IS)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     FIND THE LARGEST ABSOLUTE ELEMENT OF A , SPACED M WORDS APART            !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

   INTEGER, INTENT(IN) :: N       !! 
   INTEGER, INTENT(IN) :: M       !!
   REAL,    INTENT(IN) :: A(*)    !!
   INTEGER             :: IS

! ---------------------------------------------------------------------------- !
!                                                                              !
!      LOCAL VARIABLES                                                         !
!      ---------------                                                         !

   INTEGER :: I, INDEX

! ---------------------------------------------------------------------------- !

   IS = 1
   IF (N.LE.1) RETURN

   INDEX = 1 + M
   DO I = 2,N
      IF (ABS(A(INDEX)).GE.ABS(A(IS))) IS = I
      INDEX = INDEX + M
   END DO

   END FUNCTION ISAMAX

END SUBROUTINE SYMINV

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_ASSI_MODULE

SUBROUTINE WAMODEL

! ---------------------------------------------------------------------------- !
!                                                                              !
!    WAMODEL  - 3-G WAM MODEL - TIME INTEGRATION OF WAVE FIELDS.               !
!                                                                              !
!     S.D. HASSELMANN  MPI       1.12.85                                       !
!                                                                              !
!     G. KOMEN         KNMI         6.86  MODIFIED FOR SHALLOW WATER           !
!     P. JANSSEN                          ASPECTS.                             !
!                                                                              !
!     S.D. HASSELMANN  MPI       15.2.87  MODIFIED FOR CYBER 205.              !
!                                                                              !
!     P. LIONELLO      ISDGM      6.3.87  MODIFIED TO OUTPUT SWELL.            !
!                                                                              !
!     S.D. HASSELMANN  MPI        1.6.87  ALL VERSIONS COMBINED INTO           !
!                                         ONE MODEL. DEEP AND SHALLOW          !
!                                         WATER , CRAY AND CYBER 205           !
!                                         VERSION.                             !
!                                                                              !
!     CYCLE_2 MODICIFATIONS:                                                   !
!     ----------------------                                                   !
!                                                                              !
!     L. ZAMBRESKY     GKSS        10.87  OPTIMIZED FOR CRAY, CYBER 205        !
!     H. GUNTHER                                                               !
!                                                                              !
!     A. SPEIDEL       MPI          4.88  VARIABLE DIMENSIONS, INTERNAL        !
!                                         CHECKS (CFL-CRITERION).              !
!                                                                              !
!     A. SPEIDEL       MPI         11.88  CHANGES FOR CRAY-2.                  !
!                                                                              !
!     K. HUBBERT       POL          6.89  DEPTH AND CURRENT REFRACTION.        !
!                                         PRECALCULATION OF TERMS IN           !
!                                         *PROPDOT*.                           !
!                                         SOLVE WAVE ACTION EQUATION           !
!                                         FOR CURRENT REFRACTION.              !
!                                                                              !
!     CYCLE_3 MODICIFATIONS:                                                   !
!     ----------------------                                                   !
!                                                                              !
!     R. PORTZ , S.D. HASSELMANN   MPI          1990                           !
!                                                                              !
!      - RESTRUCTURE MODEL TO CALL THE ACTUAL INTEGRATION IN TIME              !
!        AS A SUBROUTINE: WAMODEL. A SHELL PROGRAM "WAMSHELL" READS            !
!        OUTPUT FROM PREPROC AND COMPUTES THE WIND ARRAYS FOR THE              !
!        INTEGRATION PERIOD FROM PREWIND, WHICH HAS BEEN INCORPORATED          !
!        AS A SUBROUTINE.                                                      !
!      - ALL INTERMEDIATE AND RESTART I/O IS DONE IN THE SUBROUTINE            !
!        WAMODEL AND INPREST.                                                  !
!      - THE COMMON BLOCK IN THE PREPROCESSOR AND MODEL ARE MADE               !
!        COMPATIBLE.                                                           !
!      - THE COMPUTATION OF SEVERAL PARAMETERS HAS BEEN TRANSFERRED            !
!        FROM THE MODEL TO PREPROC.                                            !
!      - DEPTH AND CURRENT REFRACTION HAS BEEN INCORPORATED INTO THE           !
!        MODEL.                                                                !
!      - OPEN BOUNDARIES ARE INCORPORATED IN THE MODEL.                        !
!      - SEVERAL MINOR ERRORS HAVE BEEN REMOVED.                               !
!      - THE BUFFERED I/O FOR THE CYBER 205 HAS BEEN CHANGED INTO A            !
!        BINARY READ AND WRITE.                                                !
!                                                                              !
!     CYCLE_4 MODICIFATIONS:                                                   !
!     ----------------------                                                   !
!                                                                              !
!     L. ZAMBRESKY   GKSS/ECMWF   6.89  ECMWF SUB VERSION                      !
!                                       BASED ON CYCLE_2.                      !
!                                                                              !
!     H. GUNTHER     GKSS/ECMWF 10.89  ECMWF SUB VERSION REORGANIZED.          !
!                                      - COMMON BLOCK STRUCTURE.               !
!                                      - BLOCKING STRUCTURE.                   !
!                                      - TIME COUNTING.                        !
!                                      - GRIDDED OUTPUT FIELDS.                !
!                                      - HEADERS ADDED TO OUTPUT FILES.        !
!                                      - ERRORS IN PROPAGATION CORRECTED       !
!                                                                              !
!     P.A.E.M. JANSSEN KNMI      1990  COUPLED MODEL.                          !
!                                                                              !
!     H. GUNTHER     GKSS/ECMWF  8.91  LOGARITHMIC DEPTH TABLES.               !
!                                      MPI CYCLE_3 AND ECMWF VERSIONS          !
!                                      COMBINED INTO CYCLE_4.                  !
!                                                                              !
!     H. GUNTHER       GKSS    DECEMBER 2001                                   !
!                                    - DKRZ AND ECMWF VERSIONS INTEGRATED.     !
!                                    - FT 90 VERSION.                          !
!                                    - DYNAMICAL ARRAYS.                       !
!                                    - ICE FIELDS INCORPORATED                 !
!                                    - INTEGRATED PARAMETER EXTENTED BY        !
!                                         TM1, TM2 PERIODS                     !
!                                         MEAN DIRECTIONAL SPREAD              !
!                                    - OUTPUT OF INTEGRATED PARAMETER          !
!                                      FOR TOTAL, SEA AND SWELL IN ONE FILE.   !
!                                    - OUTPUT OF TOTAL, SEA AND SWELL SPECTRA  !
!                                      IN ONE FILE.                            !
!                                    - SCRATCH FILES FOR WIND STORAGE REMOVED. !
!                                    - DATE TIME GROUPS EXTENTED BY CENTURAY   !
!                                      AND SECONDS.                            !
!                                    - PROPAAGATION TIMESTEP CAN BE LESS THAN  !
!                                      SOURCE FUNCTION TIME STEP.              !
!                                    - TIME INTERPOLATION OF BOUNDARY SPECTRA. !
!                                    - SPECIFICATION OF SPECTRA OUTPUT SITES   !
!                                      MOVED FROM PREPROC TO CHIEF.            !
!                                                                              !
!     R. LALBEHARRY     MSC    MARCH    2003                                   !
!                                    - Modified to save restart file once only.!
!                                                                              !
!     A. Behrens   MSC/GKSS    December 2003   Message passing with MPI        !
!                                                                              !   
!     E. Myklebust             NOVEMBER 2004 MPI parallelization               !      
!                                                                              !
!     H. GUNTHER       GKSS    JANUARY 2010  CYCLE 4.5.3                       !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF THE 2-D FREQUENCY-DIRECTION WAVE SPECTRUM AT ALL        !
!       GRID POINTS FOR A GIVEN INITIAL SPECTRUM AND FORCING SURFACE WINDS.    !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       GRID POINTS ARE LAT - LONG,VECTORIZATION IS ACHIEVED BY RUNNING        !
!       THROUGH THE GRID POINTS IN AN INNER LOOP ORGANIZED AS 1-D ARRAY.       !
!                                                                              !
!       ALL COMPONENTS OF THE SPECTRUM ARE COMPUTED PROGNOSTICALLY FROM        !
!       THE SPECTRAL TRANSPORT EQUATION UP TO A VARIABLE CUT-OFF               !
!       FREQUENCY. BEYOND THE PROGNOSTIC CUTOFF A DIAGNOSTIC F**-5 TAIL        !
!       IS ATTACHED CONTINUOUSLY FOR EACH DIRECTION,                           !
!                                                                              !
!       SOURCE FUNCTIONS ARE TAKEN FROM KOMEN ET AL(1984)                      !
!                                                                              !
!       THE NONLINEAR TRANSFER IS PARAMETERIZED BY THE DISCRETE INTER-         !
!       ACTION APPROXIMATION OF HASSELMANN ET AL (1985B)                       !
!                                                                              !
!       THE SOURCE FUNCTION AND THE ADVECTION TERM ARE INTEGRATED ON TWO       !
!       DIFFERENT TIME STEP LEVELS AND WITH DIFFERENT METHODS.                 !
!                                                                              !
!       THE SOURCE FUNCTIONS ARE INTEGRATED IMPLICITLY ACCORDING TO            !
!       HASSELMANN AND HASSELMANN (1985A),-THE RELEVANT FUNCTIONAL             !
!       DERIVATIVES OF THE INDIVIDUAL SOURCE FUNCTIONS REQUIRED FOR THE        !
!       SOLUTION OF THE IMPLICIT EQUATION ARE COMPUTED WITHIN THE SOURCE       !
!       FUNCTION SUBS,  THE TIME STEP IS TYPICALLY 20 MIN,                     !
!                                                                              !
!       THE ADVECTION IS INTEGRATED BY A FIRST ORDER UPWIND SCHEME,ALSO        !
!       ACCORDING TO HASSELMANN AND HASSELMANN (1985A), THE ADVECTIVE          !
!       TIMESTEP IS DEPENDENT ON THE FREQUENCY AND SPATIAL GRID IN             !
!       ACCORDANCE WITH CFL,                                                   !
!                                                                              !
!       WINDS ARE TAKEN EVERY WIND TIME STEP. IF THE WIND TIME STEP IS GREATER !
!       THAN THE SOURCE TERM TIME STEP THE INTERMEDIATE STEPS  ARE INTEGRATED  !
!       WITH CONSTANT WINDS.                                                   !
!       WIND TIME STEP, PROPAGATION TIME STEP AND SOURCE TERM TIME STEP        !
!       SHOULD HAVE INTEGER RATIOS, THEY ARE GIVEN IN SECONDS.                 !
!                                                                              !
!       ZERO ENERGY INFLUX IS ASSUMED AT COAST LINES. OPEN BOUNDARIES          !
!       ARE INCORPORATED IN THE MODEL, IF IT RUNS AS A NESTED GRID.            !
!                                                                              !
!       SEA POINTS ARE COUNTED ALONG LINES OF LATITUDES FROM LEFT COAST        !
!       TO RIGHT COAST WORKING FROM SOUTH TO NORTH.                            !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       SNYDER, R.L., F.W. DOBSON, J.A. ELLIOT, AND R.B. LONG:                 !
!          ARRAY MEASUREMENTS OF ATMOSPHERIC PRESSURE FLUCTUATIONS             !
!          ABOVE SURFACE GRAVITY WAVES. J.FLUID MECH. 102, 1-59 ,1981.         !
!       G. KOMEN, S. HASSELMANN, K. HASSELMANN:                                !
!          ON THE EXISTENCE OF A FULLY DEVELOPED WIND SEA SPECTRUM.            !
!          JPO,1984.                                                           !
!       S. HASSELMANN, K. HASSELMANN, J.H. ALLENDER, T.P. BARNETT:             !
!          IMPROVED METHODS OF COMPUTING AND PARAMETERIZING THE                !
!          NONLINEAR ENERGY TRANSFER IN A GRAVITY WAVE SPECTRUM.               !
!          JPO, 1985.                                                          !
!       S. HASSELMANN, K. HASSELMANN: A GLOBAL WAVE MODEL,                     !
!          WAM REPORT,JUNE,30/1985.                                            !
!       P. JANSSEN, G. KOMEN: A SHALLOW WATER EXTENSION OF THE                 !
!          3-G WAM-MODEL. WAM REPORT 1985.                                     !
!       THE WAMDI GROUP: THE WAM MODEL - A THIRD GENERATION OCEAN              !
!          WAVE PREDICTION MODEL. JPO, VOL. 18, NO. 12, 1988.                  !
!       P.A.E.M JANSSEN: JPO, 1989 AND 1991.                                   !
!       K. HASSELMANN: TRANSPORT EQUATION OF FINITE DEPTH SURFACE              !
!          WAVE SPECTRUM IN TIME DPENDANT CURRENT AND DEPTH FIELD USING        !
!          NONCANONICAL SPACIAL (SPHERICAL) AND WAVE NUMBER (FRQUENCY-         !
!          DIRECTION) COORDINATES. WAM REPROT 1988.                            !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!      ----------                                                              !

USE WAM_BOUNDARY_MODULE,      ONLY: &
&       BOUNDARY_INPUT,             & !! INPUT OF BOUNDARY VALUES.
&       BOUNDARY_OUTPUT               !! OUTPUT OF BOUNDARY VALUES.

USE WAM_CURRENT_MODULE,       ONLY: &
&       GET_CURRENT                   !! GETS A NEW CURRENT FIELD.

USE WAM_GENERAL_MODULE,       ONLY: &
&       INCDATE,                    & !! UPDATE DATE TIME GROUP.
&       flush1                        !! enforce output.

USE WAM_ICE_MODULE,           ONLY: &
&       PUT_ICE,                    & !! PUTS ICE INDICATOR INTO DATA FIELD.
&       GET_ICE                       !! GETS A NEW ICE FIELD.

USE WAM_OUTPUT_SET_UP_MODULE, ONLY: &
&       SAVE_OUTPUT_FILES,          & !! CLOSES AND OPENS OUTPUT FILES.
&       UPDATE_OUTPUT_TIME            !! UPDATES OUTPUT TIMES.

USE WAM_OUTPUT_MODULE,        ONLY: &
&       MODEL_OUTPUT_CONTROL          !! CONTROLS MODEL OUTPUT.

USE WAM_PROPAGATION_MODULE,   ONLY: &
&       PROPAGS,                    & !! PROPAGATION SCHEME.
&       PREPARE_PROPAGATION

USE WAM_RADIATION_MODULE,     ONLY: &
&       RADIATION_STRESS              !! COMPUTE RQADIATION STRESS.

USE WAM_RESTART_MODULE,       ONLY: &
&       SAVE_RESTART_FILE             !! SAVES RESTART FILES.

USE WAM_SOURCE_MODULE,        ONLY: &
&       IMPLSCH                       !! INTEGRATION OF SOURCE FUNCTION.

USE WAM_TOPO_MODULE,          ONLY: &
&       PUT_DRY,                    & !! PUTS DRY INDICATOR INTO DATA FILED.
&       GET_TOPO                      !! GETS A NEW DEPTH FIELD.

USE WAM_WIND_MODULE,          ONLY: &
&       GET_WIND                      !! GETS A NEW WIND FIELD.

use wam_assi_module,          only: & 
&       wamassi                       !! performs data assimilation

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_BOUNDARY_MODULE,      ONLY: CDT_B_OUT

USE WAM_FILE_MODULE,          ONLY: IU06, ITEST, IU20, IU25, FILE20, FILE25

USE WAM_ICE_MODULE,           ONLY: ICE_RUN, CD_ICE_NEW

USE WAM_MODEL_MODULE,         ONLY: FL3, U10, UDIR, USTAR, TAUW, Z0,           &
&                                   DEPTH, INDEP

USE WAM_NEST_MODULE,          ONLY: COARSE, FINE

USE WAM_OUTPUT_SET_UP_MODULE, ONLY: CDTINTT, CDTSPT, CDT_OUT, IDEL_OUT,        &
&                                   fflag20, fflag25

USE WAM_PROPAGATION_MODULE,   ONLY: NADV

USE WAM_RADIATION_MODULE,     ONLY: CDTOUT

USE WAM_RESTART_MODULE,       ONLY: CDT_RES

USE WAM_TIMOPT_MODULE,        ONLY: CDATEE, CDTPRO, CDTSOU, IDELPRO, IDELT,    &
&                                   SHALLOW_RUN, CDATEWO,                      &
&                                   CDTA, TOPO_RUN, CD_TOPO_NEW,               &
&                                   CDCA, CURRENT_RUN, CD_CURR_NEW, cdtstop

use wam_grid_module,          only: one_point
use wam_mpi_module,           only: nijs, nijl
use wam_assi_set_up_module,   only: iassi, cdtass

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER             :: KADV
CHARACTER (LEN=14)  :: CDTSOE
LOGICAL             :: NEW_DEPTH_OR_CURR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. ADVECTION TIME LOOP.                                                  !
!        --------------------                                                  !

PROP: DO KADV = 1,NADV

   IF (ITEST.GE.2) THEN
      WRITE(IU06,*) '   SUB. WAMODEL: START OF PROPAGATION AT        ',        &
&                   'CDTPRO = ',CDTPRO
   END IF

!     1.1 UPDATE TIMES.                                                        !
!         -------------                                                        !

   CALL INCDATE (CDTPRO,IDELPRO)     !! UPDATE END DATE OF PROPAGATION.

!     1.2 NEW DEPTH AND/OR CURRENT DATA.                                       !
!         ------------------------------                                       !

   NEW_DEPTH_OR_CURR = .FALSE.
   IF (TOPO_RUN) THEN
      IF (CDTPRO.GE.CD_TOPO_NEW) THEN
         CALL GET_TOPO (CD_TOPO_NEW)
         NEW_DEPTH_OR_CURR = .TRUE.
         IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '   SUB. WAMODEL: NEW DEPTH FIELD CDTA = ',CDTA
         END IF
      END IF
   END IF
   IF (CURRENT_RUN) THEN
      IF (CDTPRO.GE.CD_CURR_NEW) THEN
         CALL GET_CURRENT (CD_CURR_NEW)
         NEW_DEPTH_OR_CURR = .TRUE.
         IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '   SUB. WAMODEL: NEW CURRENT FIELD CDCA = ',CDCA
         END IF
      END IF
   END IF
   IF (NEW_DEPTH_OR_CURR) CALL PREPARE_PROPAGATION

!     1.3 COMPUTE OF PROPAGATION.                                              !
!         -----------------------                                              !
   
    IF (.NOT. ONE_POINT) THEN
       CALL PROPAGS (FL3)
    END IF

!     1.4 SOURCE INTEGRATION AND LOOP.                                         !
!         -----------------------------                                        !

   CDTSOE = CDTSOU                  !! END DATE OF SOURCE INTEGRATION.
   CALL INCDATE (CDTSOE,IDELT)

   PHYSICS: DO WHILE (CDTSOE.LE.CDTPRO)
      IF (ITEST.GE.2) THEN
         WRITE(IU06,*) '   SUB. WAMODEL: START OF SOURCE INTEGRATION AT ',     &
&                      'CDTSOU = ',CDTSOU
      END IF

      IF (CDTSOE.GE.CDATEWO) THEN         !! NEW WINDS IF NEEDED
         CALL GET_WIND (CDATEWO)
      END IF
!AB
      call implsch (fl3, u10, udir, tauw, ustar, z0,                          &
&                                         depth(nijs:nijl), indep(nijs:nijl))

      CDTSOU = CDTSOE                     !! UPDATE SOURCE TIME.
      CALL INCDATE (CDTSOE,IDELT)

   END DO PHYSICS

   IF (ITEST.GE.2) then
      WRITE(IU06,*)  '   SUB. WAMODEL: SOURCE INTEGRATION FINISHED'
      call flush1 (iu06)
   endif

!     1.5 INPUT OF BOUNDARY VALUES.                                            !
!         -------------------------                                            !

   IF (FINE) THEN
      CALL BOUNDARY_INPUT
      IF (ITEST.GE.2) THEN
         WRITE(IU06,*) '   SUB. WAMODEL: BOUNDARY VALUES (FINE GRID) INSERTED'
      END IF
   END IF

!     1.6 SET SPECTRA AT ICE POINTS TO ZERO.                                   !
!         ----------------------------------                                   !

   IF (ICE_RUN) THEN
      IF (CDTPRO.GT.CD_ICE_NEW) THEN
         CALL GET_ICE                 !! NEW ICE DATA
         IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '   SUB. WAMODEL: NEW ICE FIELD '
         END IF
      END IF
      CALL PUT_ICE (FL3, 0.)
      IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. WAMODEL: ICE INSERTED'
   END IF

!     1.7 SET SPECTRA AT DRY POINTS TO ZERO.                                   !
!         ----------------------------------                                   !

   IF (SHALLOW_RUN) THEN
      CALL PUT_DRY (FL3, 0.)
      IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. WAMODEL: DRY INSERTED'
   END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3.0  DATA ASSIMILATION.                                                  !
!          ------------------                                                  !

   if (iassi==1) then
      if (cdtass==cdtpro) then
         call wamassi (fl3)
         if (itest>=2) write (iu06,*) '   SUB. WAMODEL: ASSIMILATION DONE'
      endif
   endif   
     
!     1.7 OUTPUT OF BOUNDARY POINTS.                                           !
!        --------------------------                                            !

   IF (COARSE) THEN
      IF (CDT_B_OUT.EQ.CDTPRO) THEN
         CALL BOUNDARY_OUTPUT
         IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '   SUB. WAMODEL: BOUNDARY OUTPUT (COURSE GRID) DONE'
         END IF
      END IF
   END IF

!     1.8 MODEL OUTPUT TO DISK AND/OR PRINTED.                                 !
!         ------------------------------------                                 !

   IF (CDTINTT.EQ.CDTPRO .OR. CDTSPT.EQ.CDTPRO) THEN
      CALL MODEL_OUTPUT_CONTROL (fl3, iu20, iu25)

      IF (CDT_OUT.EQ.CDTPRO) THEN
         CALL SAVE_OUTPUT_FILES (IU20, FILE20, IU25, FILE25)
         CALL INCDATE (CDT_OUT, IDEL_OUT)
      END IF
      CALL UPDATE_OUTPUT_TIME                     !! UPDATE OUTPUT TIMES.
      IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. WAMODEL: MODEL_OUTPUT_CONTROL DONE'
      call flush1 (iu06)
   END IF

!     1.9 SAVE RECOVERY FILES WHEN TIME REACHES THE SAVE DATE.                 !
!         ----------------------------------------------------                 !

   if (cdt_res==cdtpro.and.cdt_res<=cdtstop) call save_restart_file

!     1.10 RADIATION STRESS.                                                   !
!          -----------------                                                   !

   IF (CDTPRO.EQ.CDTOUT) THEN
      CALL RADIATION_STRESS (fl3)
      IF (ITEST.GE.2) THEN
         WRITE (IU06,*) '   SUB. WAMODEL: RADIATION_STRESS DONE '
      END IF
   END IF

!     1.11 PRINT TIME.                                                         !
!          -----------                                                         !

   WRITE (IU06,'(/,3X,''!!!!!!!!!!!!!! WAVE FIELDS INTEGRATED  DATE IS: '',    &
&                A14,''  !!!!!!!!!!!!!! '')') CDTPRO

   call flush1 (iu06)
END DO PROP

! ---------------------------------------------------------------------------- !

END SUBROUTINE WAMODEL

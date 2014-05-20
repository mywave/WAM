SUBROUTINE PRINT_WAM_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PRINT_WAM_STATUS - PRINT STATUS OF ALL MODULES USED IN WAM.                !
!                                                                              !
!     H. GUNTHER     GKSS/ECMWF     NOVEMBER 1989                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO PROVIDE AN OVERVIEW OF THE STATUS AND DATA STORED IN THE MODULES OF !
!       THE WAM MODEL.                                                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       DONE BY CALL TO THE PRIND SUBROUTINES INCULDED IN THE DIFFERENT        !
!       MODULES.                                                               !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!     -----------                                                              !

USE WAM_BOUNDARY_MODULE,      ONLY: PRINT_BOUNDARY_STATUS
USE WAM_COLDSTART_MODULE,     ONLY: PRINT_COLDSTART_STATUS
USE WAM_CURRENT_MODULE,       ONLY: PRINT_CURRENT_STATUS
USE WAM_FILE_MODULE,          ONLY: PRINT_FILE_STATUS
USE WAM_FRE_DIR_MODULE,       ONLY: PRINT_FRE_DIR_STATUS
USE WAM_GRID_MODULE,          ONLY: PRINT_GRID_STATUS
USE WAM_ICE_MODULE,           ONLY: PRINT_ICE_STATUS
USE WAM_NEST_MODULE,          ONLY: PRINT_NEST_STATUS
USE WAM_OUTPUT_SET_UP_MODULE, ONLY: PRINT_OUTPUT_STATUS
USE WAM_PROPAGATION_MODULE,   ONLY: PRINT_PROPAGATION_STATUS
USE WAM_RADIATION_MODULE,     ONLY: PRINT_RADIATION_MODULE
USE WAM_RESTART_MODULE,       ONLY: PRINT_RESTART_STATUS
USE WAM_SOURCE_MODULE,        ONLY: PRINT_SOURCE_STATUS
USE WAM_TIMOPT_MODULE,        ONLY: PRINT_TIMOPT_STATUS
USE WAM_TOPO_MODULE,          ONLY: PRINT_TOPO_STATUS
USE WAM_WIND_MODULE,          ONLY: PRINT_WIND_STATUS
use wam_assi_set_up_module,   only: print_assimilation_status

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. PRINT STATUS.                                                         !
!        -------------                                                         !

CALL PRINT_TIMOPT_STATUS
CALL PRINT_COLDSTART_STATUS
CALL PRINT_FRE_DIR_STATUS
CALL PRINT_GRID_STATUS
CALL PRINT_WIND_STATUS
CALL PRINT_TOPO_STATUS
CALL PRINT_CURRENT_STATUS
CALL PRINT_ICE_STATUS
CALL PRINT_NEST_STATUS
CALL PRINT_BOUNDARY_STATUS
CALL PRINT_OUTPUT_STATUS
CALL PRINT_FILE_STATUS
CALL PRINT_RESTART_STATUS
CALL PRINT_PROPAGATION_STATUS
CALL PRINT_SOURCE_STATUS
CALL PRINT_RADIATION_MODULE
call print_assimilation_status

END SUBROUTINE PRINT_WAM_STATUS

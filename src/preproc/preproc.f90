PROGRAM PREPROC

! ---------------------------------------------------------------------------- !
!                                                                              !
!   PROGRAM PREPROC - PREPARE DATA (BUT NOT WINDS) FOR INPUT TO WAM WAVE MODEL.!
!                                                                              !
!     SUSANNE HASSELMANN  MPI     JUNE 1986.                                   !
!                                                                              !
!     ANNEGRET SPEIDEL    MPI  OCTOBER 1988. MODFIED FOR CYCLE_2.              !
!                                                                              !
!     K. HUBBERT          POL     JUNE 1989  DEPTH AND CURRENT                 !
!                                            REFRACTION.                       !
!                                                                              !
!     H. GUNTHER   ECMWF/GKSS    APRIL 1990  LAND POINTS ARE REMOVED           !
!                                            FROM BLOCKS AND THE CODE          !
!                                            HAS BEEN RESTRUCTURED.            !
!                                                                              !
!     R. PORTZ     MPI         JANUARY 1991  NESTED GRID OPTION.               !
!                                                                              !
!     H. GUNTHER   ECMWF/GKSS    APRIL 1991  CYCLE_4 MODIFICATIONS.            !
!                                            MULTI-PART REMOVED.               !
!                                            NEW SOURCE FUNCTIONS.             !
!                                            LOG. DEPTH TABLE.                 !
!     H. GUNTHER  GKSS         JANUARY 2002  FT90.                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO ARRANGE A GRID FOR THE WAM WAVE MODEL AND COMPUTE ALL FIXED MODEL   !
!       PARAMETERS WHICH ARE STORED IN WAM_CONST_MODULE.                       !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       A REPRESENTATIVE TOPOGRAPHIC DATA SET ON LAT-LONG COORDINATES          !
!       CONTAINING THE MODEL SQUARE BOX REGION IS READ IN. THE MODEL REGION IS !
!       EXTRACTED AND INTERPOLATED ONTO GIVEN LAT-LONG GRID INCREMENTS.        !
!       THE PROGRAM CHECKS FOR A PERIODIC LATITUDE GRID. IF THE GRID IS NOT    !
!       PERIODIC A CLOSED BASIN IS ASSUMED.                                    !
!       THE PROGRAM DOES NOT DISTINGUISH BETWEEN DEEP AND SHALLOW WATER.       !
!                                                                              !
!       -BLOCK STRUCTURE :                                                     !
!        GRID POINTS ARE COLLECTED INTO A 1-DIMENSIONAL ARRAY, GRID POINTS     !
!        (ONLY SEAPOINTS) ARE COUNTED ALONG LATITUDES FROM WEST TO EAST        !
!        WORKING FROM SOUTH TO NORTH.                                          !
!                                                                              !
!       -NESTED GRIDS: THE GRID GENERATED CAN BE A                             !
!         - COARSE GRID WHICH MEANS OUTPUT OF BOUNDARY SPECTRA                 !
!                       FOR A FOLLOW UP FINE GRID RUN.                         !
!         - FINE   GRID WHICH MEANS INPUT OF  BOUNDARY SPECTRA                 !
!                       FROM AN EARLIER COARSE GRID RUN.                       !
!         - COARSE AND A FINE GRID                                             !
!                                                                              !
!       - REFRACTION: IF A CURRENT FILE IS PROVIDED TO THE PROGRAM             !
!         THE CURRENT FIELD IS READ, INTERPOLATED TO THE MODEL GRID.           !
!                                                                              !
!     INTERFACE.                                                               !
!     ----------                                                               !
!                                                                              !
!     INPUT, OUTPUT FILES AND UNITS.                                           !
!     ------------------------------                                           !
!                                                                              !
!       ALL UNITS AND FILE NAMES ARE STORED IN WAM_FILE_MODULE.                !
!       ALL FILES ARE DYNAMICALLY ASSIGNED TO A UNIT.                          !
!                                                                              !
!       THE PROGRAM ASSIGNS TO UNIT IU05 = 5 THE INPUT FILE "Preproc_User" AND !
!       TO UNIT IU06 = 6 THE OUTPUT FILE "Preproc_Prot".                       !
!       THESE FILE NAMES ARE DEFINED IN SECTION 1 OF THIS PROGRAM              !
!                                                                              !
!       ALL OTHER FILE NAMES CAN BE REDEFINED OR HAVE TO BE DEFINED IN THE     !
!       INPUT FILE "Preproc_User"                                              !
!                                                                              !
!     MODULES.                                                                 !
!     --------                                                                 !
!                                                                              !
!       WAM_FILE_MODULE         - WAM FILES AND UNITS                          !
!       WAM_CONST_MODULE        - WAM MODEL CONSTANTS                          !
!       WAM_INTERFACE_MODULE    - TO PROVIDE SUBROUTINES TO PREPROC            !
!       PREPROC_MODULE          - VARIABLES AND CONSTANTS FOR PREPROC          !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !

USE WAM_FILE_MODULE,      ONLY:  &
&       SET_USER_FILE,           &  !!  SET USER FILE
&       SET_PROTOCOL_FILE           !!  SET PROTOCOL FILE

USE PREPROC_MODULE, ONLY:        &
&       READ_PREPROC_USER,       &  !! READ USER INPUT FILE.
&       READ_TOPOGRAPHY,         &  !! READ TOPOGRAPHY INPUT FILE.
&       PREPARE_CONST,           &  !! PREPARE WAM CONST MODULE.
&       PRINT_PREPROC_STATUS,    &  !! PRINT WAM CONST MODULE
&       WRITE_PREPROC_FILE          !! WRITES PREPROC OUTPUT FILE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY: IU06, FILE06, ITEST

IMPLICIT NONE
include "mpif.h"
integer ierror

call mpi_init (ierror)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CHANGE DEFAULT FILES TO PREPROC VALUES.                               !
!        ---------------------------------------                               !

CALL SET_USER_FILE    ('Preproc_User')    !! SET USER FILE
CALL SET_PROTOCOL_FILE('Preproc_Prot')    !! SET PROTOCOL FILE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. OPEN PROTOCOLL FILE.                                                  !
!        --------------------                                                  !

OPEN (UNIT=IU06, FILE=FILE06, FORM='FORMATTED', STATUS='UNKNOWN')

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. READ INPUT AND PRINT PROTOCOL.                                        !
!        ------------------------------                                        !

CALL READ_PREPROC_USER
IF (ITEST.GT.0) WRITE (IU06,*) ' SUB READ_INPUT_PREPROC DONE'

CALL READ_TOPOGRAPHY
IF (ITEST.GT.0) WRITE (IU06,*) ' SUB READ_TOPOGRAPHY DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. PREPARE WAM CONST MODULE.                                             !
!        -------------------------                                             !

CALL PREPARE_CONST
IF (ITEST.GT.0) WRITE (IU06,*) ' SUB PREPARE_CONST DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. WRITE THE COMPUTED INFORMATION TO FILE.                               !
!        ---------------------------------------                               !

CALL WRITE_PREPROC_FILE
IF (ITEST.GT.0) WRITE (IU06,*) ' SUB WRITE_PREPROC_FILE DONE'

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. PRINT THE COMPUTED INFORMATION.                                       !
!        -------------------------------                                       !

WRITE(IU06,*) '  '
WRITE(IU06,*) ' ------------------------------------------------- '
WRITE(IU06,*) '              MODULE STATUS AT END OF PREPROC:'
WRITE(IU06,*) ' ------------------------------------------------- '
CALL PRINT_PREPROC_STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6. END OF PREPROC.                                                       !
!        ---------------                                                       !

WRITE (IU06,*) ' '
WRITE (IU06,*) ' PROGRAM PREPROC: ALL DONE'

call mpi_finalize (ierror)
END PROGRAM PREPROC

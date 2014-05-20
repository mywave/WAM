MODULE WAM_MODEL_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS: 
!   THE TIME STATUS OF INTEGRATION, INTEGRATION OPTIONS, AND
!   THE ACTUAL MODEL WAVE SPECTRA, WIND, DEPTH AND CURRENT INFORMATION. 
!
! ---------------------------------------------------------------------------- !

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. MODEL SPECTRA.                                                        !
!        --------------                                                        !

REAL,    ALLOCATABLE :: FL3 (:,:,:) !! SPECTRA.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. WIND FIELDS USED IN THE WAM-MODEL.                                    !
!        ----------------------------------                                    !

REAL,    ALLOCATABLE :: U10  (:)   !! WIND SPEED [M/S].
REAL,    ALLOCATABLE :: UDIR (:)   !! WIND DIRECTION IN [RAD].
                                   !! OCEANOGRAPHIC NOTATION (POINTING ANGLE OF
                                   !! WIND VECTOR,  CLOCKWISE FROM NORTH).
REAL,    ALLOCATABLE :: USTAR (:)  !! FRICTION VELOCITY [M/S].
REAL,    ALLOCATABLE :: Z0    (:)  !! ROUGHNESS LENGTH [M]. 
REAL,    ALLOCATABLE :: TAUW  (:)  !! WAVE STRESS IN (M/S)**2 

REAL,    ALLOCATABLE :: DEPTH (:)  !! WATER DEPTH IN (M)
INTEGER, ALLOCATABLE :: INDEP(:)   !! DEPTH INDEX.

REAL,    ALLOCATABLE :: U (:)      !! WEST-EAST COMP. OF CURRENTS IN (M/S)
REAL,    ALLOCATABLE :: V (:)      !! SOUTH-NORTH COMP. OF CURRENTS IN (M/S)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_MODEL_MODULE

MODULE WAM_FLUX_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS: 
!   ENERGY AND MOMENTUM FLUXES IN THE WAM-MODEL.
!
! ---------------------------------------------------------------------------- !

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. ENERGY AND MOMENTUM FLUXES IN THE WAM-MODEL.                          !
!        --------------------------------------------                          !

REAL,    ALLOCATABLE :: PHIOC(:)    !! ENERGY FLUX TO OCEAN.
REAL,    ALLOCATABLE :: PHIAW(:)    !! ENERGY FLUX FROM WIND TO WAVES.
REAL,    ALLOCATABLE :: TAUOC_X(:)  !! MOMENTUM FLUX INTO OCEAN.
REAL,    ALLOCATABLE :: TAUOC_Y(:)  !! MOMENTUM FLUX INTO OCEAN.
REAL,    ALLOCATABLE :: PHIBOT(:)   !! BOTTOM ENERGY FLUX TO OCEAN.
REAL,    ALLOCATABLE :: TAUBOT_X(:) !! BOTTOM MOMENTUM FLUX INTO OCEAN.
REAL,    ALLOCATABLE :: TAUBOT_Y(:) !! BOTTOM MOMENTUM FLUX INTO OCEAN.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_FLUX_MODULE

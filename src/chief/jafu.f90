! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

   INTEGER FUNCTION JAFU (CL, J, IAN)

   ! ------------------------------------------------------------------------- !
   !                                                                           !
   !   JAFU - FUNCTION TO COMPUTE THE INDEX ARRAY FOR THE ANGLES OF THE        !
   !          INTERACTING WAVENUMBERS.                                         !
   !                                                                           !
   !     S. HASSELMANN        MPIFM        01/12/1985.                         !
   !                                                                           !
   !     PURPOSE.                                                              !
   !     --------                                                              !
   !                                                                           !
   !       INDICES DEFINING BINS IN FREQUENCY AND DIRECTION PLANE INTO         !
   !       WHICH NONLINEAR ENERGY TRANSFER INCREMENTS ARE STORED. NEEDED       !
   !       FOR COMPUTATION OF THE NONLINEAR ENERGY TRANSFER.                   !
   !                                                                           !
   !     METHOD.                                                               !
   !     -------                                                               !
   !                                                                           !
   !       SEE REFERENCE.                                                      !
   !                                                                           !
   !     REFERENCE.                                                            !
   !     ----------                                                            !
   !                                                                           !
   !        S. HASSELMANN AND K. HASSELMANN,JPO, 1985 B.                       !
   !                                                                           !
   ! ------------------------------------------------------------------------- !

   REAL,    INTENT(IN) :: CL    !! WEIGHTS.
   INTEGER, INTENT(IN) :: J     !! INDEX IN ANGULAR ARRAY.
   INTEGER, INTENT(IN) :: IAN   !! NUMBER OF ANGLES IN ARRAY.

   JAFU = J + INT(CL)
   IF (JAFU.LE.0)   JAFU = JAFU+IAN
   IF (JAFU.GT.IAN) JAFU = JAFU-IAN

   END FUNCTION JAFU

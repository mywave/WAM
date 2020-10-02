MODULE WAM_OUTPUT_PARAMETER_MODULE

! ---------------------------------------------------------------------------- !
!                                                                              !
!   THIS MODULE CONTAINS: OUTPUT PARAMTER NAMES AND SCALES FOR PRINTING.       !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     A.  EXTERNALS.                                                           !
!                                                                              !
! ---------------------------------------------------------------------------- !

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     B. VARIABLES FROM OTHER MODULES.                                         !
!                                                                              !
! ---------------------------------------------------------------------------- !

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     C. MODULE VARIABLES.                                                     !
!                                                                              !
! ---------------------------------------------------------------------------- !

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. NUMBER OF INTEGRATED PARAMETER.                                       !
!        -------------------------------                                       !

INTEGER, PARAMETER :: NOUT_P = 70

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TITLE FOR OUTPUT PARAMETER AND SPECTRA.                               !
!        ---------------------------------------                               !

CHARACTER(LEN=60), DIMENSION(NOUT_P) :: TITL_P = (/                  &
& ' WIND SPEED U10 ( METRES/SECOND )                           ',    &   !!  1
& ' WIND DIRECTION ( DEGREE FROM NORTH TO )                    ',    &   !!  2
& ' FRICTION VELOCITY ( METRES/SECOND )                        ',    &   !!  3
& ' DRAG COEFFICIENT ( PROMILLE )                              ',    &   !!  4
& ' CHARNOCK PARAMETER                                         ',    &   !!  5
& ' WATER DEPTH (METRES) (DEEPER THAN 999M ARE PRINTED AS 999) ',    &   !!  6
& ' CURRENT SPEED ( METRES/SECOND )                            ',    &   !!  7
& ' CURRENT DIRECTION ( DEGREE FROM NORTH TO )                 ',    &   !!  8
& ' SIGNIFICANT WAVE HEIGHT ( METRES )                         ',    &   !!  9
& ' WAVE PEAK PERIOD ( SECONDS )                               ',    &   !! 10
& ' WAVE MEAN PERIOD (SECONDS )                                ',    &   !! 11
& ' WAVE TM1 PERIOD ( SECONDS )                                ',    &   !! 12
& ' WAVE TM2 PERIOD ( SECONDS )                                ',    &   !! 13
& ' WAVE DIRECTION ( DEGREE FROM NORTH TO )                    ',    &   !! 14
& ' DIRECTIONAL SPREAD ( DEGREES )                             ',    &   !! 15
& ' NORMALISED WAVE STRESS ( % )                               ',    &   !! 16
& ' SEA SIGNIFICANT WAVE HEIGHT ( METRES )                     ',    &   !! 17
& ' SEA PEAK PERIOD ( SECONDS )                                ',    &   !! 18
& ' SEA MEAN PERIOD ( SECONDS )                                ',    &   !! 19
& ' SEA TM1 PERIOD ( SECONDS )                                 ',    &   !! 20
& ' SEA TM2 PERIOD (  SECONDS )                                ',    &   !! 21
& ' SEA DIRECTION ( DEGREE FROM NORTH TO )                     ',    &   !! 22
& ' SEA DIRECTIONAL SPREAD ( DEGREES )                         ',    &   !! 23
& ' DUMMY                                                      ',    &   !! 24
& ' SWELL SIGNIFICANT WAVE HEIGHT ( METRES )                   ',    &   !! 25
& ' SWELL PEAK PERIOD ( SECONDS )                              ',    &   !! 26
& ' SWELL MEAN PERIOD ( SECONDS )                              ',    &   !! 27
& ' SWELL TM1 PERIOD ( SECONDS )                               ',    &   !! 28
& ' SWELL TM2 PERIOD ( SECONDS )                               ',    &   !! 29
& ' SWELL DIRECTION ( DEGREE FROM NORTH TO )                   ',    &   !! 30
& ' SWELL DIRECTIONAL SPREAD ( DEGREES )                       ',    &   !! 31
& ' ROUGHNESS LENGTH Z0 ( METRES )                             ',    &   !! 32
& ' GODA PEAKEDNESS PARAMETER                                  ',    &   !! 33
& ' KURTOSIS                                                   ',    &   !! 34
& ' BENJAMIN-FEIR INDEX                                        ',    &   !! 35
& ' NORMALIZED MAXIMUM WAVE HEIGHT                             ',    &   !! 36
& ' MAXIMUM WAVE PERIOD ( SECONDS )                            ',    &   !! 37
& ' PEAK FREQUENCY (INTERPOLATED) ( HZ )                       ',    &   !! 38
& ' PEAK DIRECTION ( DEGREE FROM NORTH TO )                    ',    &   !! 39
& ' MEAN SQUARE SLOPE                                          ',    &   !! 40
& ' FIRST SWELL SIGNIFICANT WAVE HEIGHT ( METRES )             ',    &   !! 41
& ' FIRST SWELL TM1 PERIOD ( SECONDS )                         ',    &   !! 42
& ' FIRST SWELL DIRECTION ( DEGREE FROM NORTH TO )             ',    &   !! 43
& ' SECOND SWELL SIGNIFICANT WAVE HEIGHT ( METRES )            ',    &   !! 44
& ' SECOND SWELL TM1 PERIOD ( SECONDS )                        ',    &   !! 45
& ' SECOND SWELL DIRECTION ( DEGREE FROM NORTH TO )            ',    &   !! 46
& ' THIRD SWELL SIGNIFICANT WAVE HEIGHT ( METRES )             ',    &   !! 47
& ' THIRD SWELL TM1 PERIOD ( SECONDS )                         ',    &   !! 48
& ' THIRD SWELL DIRECTION ( DEGREE FROM NORTH TO )             ',    &   !! 49
& ' DUMMY                                                      ',    &   !! 50
& ' RADIATION STRESS TENSOR SXX ( KG/S/S )                     ',    &   !! 51
& ' RADIATION STRESS TENSOR SYY ( KG/S/S )                     ',    &   !! 52
& ' RADIATION STRESS TENSOR SXY ( KG/S/S )                     ',    &   !! 53
& ' DUMMY                                                      ',    &   !! 54
& ' X-COMP. WAVE FORCE PER SURFACE UNIT ( N/M/M )              ',    &   !! 55
& ' Y-COMP. WAVE FORCE PER SURFACE UNIT ( N/M/M )              ',    &   !! 56
& ' X-COMP. STOKES DRIFT ( M/S )                               ',    &   !! 57
& ' Y-COMP. STOKES DRIFT ( M/S )                               ',    &   !! 58
& ' ENERGY FLUX TO OCEAN ( KG/S/S/S )                          ',    &   !! 59
& ' TOTAL ENERGY FLUX FROM WIND TO WAVES ( KG/S/S/S )          ',    &   !! 60
& ' X-COMP.  MOMENTUM FLUX INTO OCEAN ( KG/M/S/S )             ',    &   !! 61
& ' Y-COMP.  MOMENTUM FLUX INTO OCEAN ( KG/M/S/S )             ',    &   !! 62
& ' ENERGY FLUX FROM WAVES TO BOTTOM ( KG/S/S/S )              ',    &   !! 63
& ' X-COMP.  MOMENTUM FLUX FROM WAVES TO BOTTOM ( KG/M/S/S )   ',    &   !! 64
& ' Y-COMP.  MOMENTUM FLUX FROM WAVES TO BOTTOM ( KG/M/S/S )   ',    &   !! 65
& ' DUMMY                                                      ',    &   !! 66
& ' CREST MAX (TIME, FORRISTALL)                               ',    &   !! 67  !! WAM-MAX
& ' HMAX (TIME, NAESS)                                         ',    &   !! 68  !! WAM-MAX
& ' MAXIMUM CREST H.- SPACE-TIME (STQD)                        ',    &   !! 69  !! WAM-MAX
& ' MAXIMUM WAVE H.- SPACE-TIME (STQD)                         '/)       !! 70  !! WAM-MAX

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. SCALING FACTORS FOR OUTPUT PARAMETER.                                 !
!        -------------------------------------                                 !

REAL, PARAMETER, DIMENSION(NOUT_P) :: SCAL_P = (/                              &
&                     10.            ,    &   !!  1
&                      1.            ,    &   !!  2
&                    100.            ,    &   !!  3
&                  10000.            ,    &   !!  4
&                  10000.            ,    &   !!  5
&                      1.            ,    &   !!  6
&                    100.            ,    &   !!  7
&                      1.            ,    &   !!  8
&                     10.            ,    &   !!  9
&                     10.            ,    &   !! 10
&                     10.            ,    &   !! 11
&                     10.            ,    &   !! 12
&                     10.            ,    &   !! 13
&                      1.            ,    &   !! 14
&                      1.            ,    &   !! 15
&                    100.            ,    &   !! 16
&                     10.            ,    &   !! 17
&                     10.            ,    &   !! 18
&                     10.            ,    &   !! 19
&                     10.            ,    &   !! 20
&                     10.            ,    &   !! 21
&                      1.            ,    &   !! 22
&                      1.            ,    &   !! 23
&                      1.            ,    &   !! 24
&                     10.            ,    &   !! 25
&                     10.            ,    &   !! 26
&                     10.            ,    &   !! 27
&                     10.            ,    &   !! 28
&                     10.            ,    &   !! 29
&                      1.            ,    &   !! 30
&                      1.            ,    &   !! 31
&                  10000.            ,    &   !! 32
&                     10.            ,    &   !! 33
&                    100.            ,    &   !! 34
&                     10.            ,    &   !! 35
&                     10.            ,    &   !! 36
&                     10.            ,    &   !! 37
&                   1000.            ,    &   !! 38
&                      1.            ,    &   !! 39
&                   1000.            ,    &   !! 40
&                     10.            ,    &   !! 41
&                     10.            ,    &   !! 42
&                      1.            ,    &   !! 43
&                     10.            ,    &   !! 44
&                     10.            ,    &   !! 45
&                      1.            ,    &   !! 46
&                     10.            ,    &   !! 47
&                     10.            ,    &   !! 48
&                      1.            ,    &   !! 49
&                      1.            ,    &   !! 50
&                      0.            ,    &   !! 51
&                      0.            ,    &   !! 52
&                      0.            ,    &   !! 53
&                      1.            ,    &   !! 54
&                      0.            ,    &   !! 55
&                      0.            ,    &   !! 56
&                    100.            ,    &   !! 57
&                    100.            ,    &   !! 58
&                    100.            ,    &   !! 59
&                    100.            ,    &   !! 60
&                    100.            ,    &   !! 61
&                    100.            ,    &   !! 62
&                    100.            ,    &   !! 63
&                    100.            ,    &   !! 64
&                    100.            ,    &   !! 65
&                    100.            ,    &   !! 66
&                     10.            ,    &   !! 67  !! WAM-MAX
&                     10.            ,    &   !! 68  !! WAM-MAX
&                     10.            ,    &   !! 69  !! WAM-MAX
&                     10.            /)       !! 70  !! WAM-MAX

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. FLAG FOR DIRECTION PARAMTERS.                                         !
!        -----------------------------                                         !

LOGICAL, DIMENSION(NOUT_P) :: DIR_FLAG = (/                          &
& .FALSE.   ,    &   !!  1
& .TRUE.    ,    &   !!  2
& .FALSE.   ,    &   !!  3
& .FALSE.   ,    &   !!  4
& .FALSE.   ,    &   !!  5
& .FALSE.   ,    &   !!  6
& .FALSE.   ,    &   !!  7
& .TRUE.    ,    &   !!  8
& .FALSE.   ,    &   !!  9
& .FALSE.   ,    &   !! 10
& .FALSE.   ,    &   !! 11
& .FALSE.   ,    &   !! 12
& .FALSE.   ,    &   !! 13
& .TRUE.    ,    &   !! 14
& .FALSE.   ,    &   !! 15
& .FALSE.   ,    &   !! 16
& .FALSE.   ,    &   !! 17
& .FALSE.   ,    &   !! 18
& .FALSE.   ,    &   !! 19
& .FALSE.   ,    &   !! 20
& .FALSE.   ,    &   !! 21
& .TRUE.    ,    &   !! 22
& .FALSE.   ,    &   !! 23
& .FALSE.   ,    &   !! 24
& .FALSE.   ,    &   !! 25
& .FALSE.   ,    &   !! 26
& .FALSE.   ,    &   !! 27
& .FALSE.   ,    &   !! 28
& .FALSE.   ,    &   !! 29
& .TRUE.    ,    &   !! 30
& .FALSE.   ,    &   !! 31
& .FALSE.   ,    &   !! 32
& .FALSE.   ,    &   !! 33
& .FALSE.   ,    &   !! 34
& .FALSE.   ,    &   !! 35
& .FALSE.   ,    &   !! 36
& .FALSE.   ,    &   !! 37
& .FALSE.   ,    &   !! 38
& .TRUE.    ,    &   !! 39
& .FALSE.   ,    &   !! 40
& .FALSE.   ,    &   !! 41
& .FALSE.   ,    &   !! 42
& .TRUE.    ,    &   !! 43
& .FALSE.   ,    &   !! 44
& .FALSE.   ,    &   !! 45
& .TRUE.    ,    &   !! 46
& .FALSE.   ,    &   !! 47
& .FALSE.   ,    &   !! 48
& .TRUE.    ,    &   !! 49
& .FALSE.   ,    &   !! 50
& .FALSE.   ,    &   !! 51
& .FALSE.   ,    &   !! 52
& .FALSE.   ,    &   !! 53
& .FALSE.   ,    &   !! 54
& .FALSE.   ,    &   !! 55
& .FALSE.   ,    &   !! 56
& .FALSE.   ,    &   !! 57
& .FALSE.   ,    &   !! 58
& .FALSE.   ,    &   !! 59
& .FALSE.   ,    &   !! 60
& .FALSE.   ,    &   !! 61
& .FALSE.   ,    &   !! 62
& .FALSE.   ,    &   !! 63
& .FALSE.   ,    &   !! 64
& .FALSE.   ,    &   !! 65
& .FALSE.   ,    &   !! 66
& .FALSE.   ,    &   !! 67
& .FALSE.   ,    &   !! 68
& .FALSE.   ,    &   !! 69
& .FALSE.   /)       !! 70

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. NUMBER OF OUTPUT SPECTRA.                                             !
!        -------------------------                                             !

INTEGER, PARAMETER :: NOUT_S = 4

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. TITLE OF OUTPUT SPECTRA.                                              !
!        -----------------------                                               !

CHARACTER(LEN=60), DIMENSION(NOUT_S) :: TITL_S = (/                  &
& ' SPECTRUM                                                   ',    &   !!  1
& ' SEA SPECTRUM                                               ',    &   !!  2
& ' SWELL SPECTRUM                                             ',    &   !!  3
& ' DUMMY                                                      '/)       !!  4

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_OUTPUT_PARAMETER_MODULE

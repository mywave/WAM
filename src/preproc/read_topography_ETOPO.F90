SUBROUTINE READ_TOPOGRAPHY

! ---------------------------------------------------------------------------- !
!                                                                              !
!    READ_TOPOGRAPHY - ETOPO_INPUT.                                            !
!                                                                              !
!        ETOPO_INPUT - READS THE FILES                                         !
!                      etopo2_2006apr.dat AND 'reference_levels'               !
!                      AND MERGES THE DATA                                     !
!                                                                              !
!                                                                              !
!      J. BIDLOT       ECMWF     JUNE 2007                                     !
!      H. GUENTHER     HZG       JANUARY 2015                                  !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO READ THE FILES  etopo2_2006apr.dat AND 'reference_levels' ,         !
!       TO MERGE THE DATA AND TRANSFER THE DATA TO THE PREPROC MODULE.         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       FILE08, WHICH IS DEFINED IN THE USER INPUT, IS ASSIGNED TO IU08.       !
!       THE FILE MUST STORE:                                                   !
!           1. RECORD: THE TOPOGRAPHIC DATA HEADER.                            !
!           FOLLOWING RECORDS: THE DEPTH DATA MATRIX.                          !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY SUB. SET_TOPOGRAPHY !
!       TO THE WAM_GRID_MODULE:                                                !
!         INTEGER    :: N_LON      !! NUMBER OF LONGITUDES.                    !
!         INTEGER    :: N_LAT      !! NUMBER OF LATITUDES.                     !
!         REAL*8     :: D_LAT      !! LATITUDE INCREMENT.                      !
!         REAL*8     :: D_LON      !! LONGITUDE INCREMENT.                     !
!         REAL*8     :: SOUTH      !! SOUTH LATITUDE.                          !
!         REAL*8     :: NORTH      !! NORTH LATITUDE.                          !
!         REAL*8     :: WEST       !! WEST LONGITUDE.                          !
!         REAL*8     :: EAST       !! EAST LONGITUDE.                          !
!         REAL       :: D_MAP(:,:) !! WATER DEPTH [M].                         !
!                                                                              !
!      ALL INCREMENTS, LATITUDES AND LONGITUDES MUST BE REAL*8 IN DEGREES      !
!      OR  INTEGER IN M_SEC, OR A CHARCTER STRING                              !
!                                                                              !
!      THE TOPOGRAPHY MUST BE ON A REGULAR LATITUDE-LONGITUDE GRID ARRANGED    !
!      FROM  WEST TO EAST AND FROM SOUTH TO NORTH, WHICH IS                    !
!      THE DEPTH ARRAY "D_MAP(I,K)" MUST BE ORDERED AS                         !
!                 (    1,    1 ) <==> SOUTH WEST                               !
!                 (N_LON,    1 ) <==> SOUTH EAST                               !
!                 (    1, N_LAT) <==> NORTH WEST                               !
!                 (N_LON, N_LAT) <==> NORTH EAST                               !
!       POSITIVE VALUES ARE SEA DEPTHS AND NEGATIVE VALUES ARE LAND.           !
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

USE WAM_COORDINATE_MODULE           !! COORDINATE PROCEDURES

USE WAM_GENERAL_MODULE,   ONLY:  &
&       ABORT1                      !! TERMINATES PROCESSING.

USE PREPROC_MODULE,       ONLY:  &
&       SET_TOPOGRAPHY              !! TRANSFERS DEPTH DATA TO MODULE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,   ONLY: IU06, IU08, FILE08

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!      LOCAL VARIABLES.
!      ----------------

INTEGER, PARAMETER :: ILON=10801  !! no. of longitudes in etopo2_2006apr.dat
INTEGER, PARAMETER :: ILAT=5400   !! no. of latitudes  in etopo2_2006apr.dat

INTEGER, ALLOCATABLE :: IDEPTH(:,:)  !! Depth data of etopo2_2006apr.dat.
REAL (KIND=KIND_D)   ::  AALON(ILON)       !! longitudes in etopo2_2006apr.dat
REAL (KIND=KIND_D)   ::  AALAT(ILAT)       !! latitudes  in etopo2_2006apr.dat
INTEGER ::  WEST
INTEGER ::  EAST
INTEGER ::  NORTH
INTEGER ::  SOUTH
INTEGER ::  D_LAT
INTEGER ::  D_LON

INTEGER, PARAMETER :: NREF=300       !! Max No. of entries in reference_levels
INTEGER            :: NREFERENCE     !! No. of entries in reference_levels
INTEGER            :: LEVEL(NREF)    !! New reference levels
INTEGER            :: NDEPTH(NREF)   !! new depths
REAL               :: XINF(NREF)     !! West longitude of areas.
REAL               :: XSUP(NREF)     !! East longitude of areas.
REAL               :: YINF(NREF)     !! South Latitude of areas.
REAL               :: YSUP(NREF)     !! North Latitude of areas.
CHARACTER(LEN=72)  :: LOCATION(NREF) !! Area names.
REAL, ALLOCATABLE :: RDEPTH(:,:)  !! Depth data of etopo2_2006apr.dat.

INTEGER :: I, J, IJ, NTOT, IC, ICOUNT, IR, IOS
INTEGER :: IDUM(15)   !! Work array for etopo2_2006apr.da input.
REAL    :: XI, YJ
REAL (KIND=KIND_D) :: X
CHARACTER(LEN=LEN_COOR) :: C
CHARACTER(LEN=20) :: FILE_REF

! ---------------------------------------------------------------------------- !
!
!     1. COMPUTE THE LATITUDE AND LONGITUDE POINTS.
!        ------------------------------------------

X = -180.
WEST  = DEG_TO_M_SEC (X)
X =  180.
EAST  = DEG_TO_M_SEC (X)
X = 90.
NORTH = DEG_TO_M_SEC (X)
X = -90.
SOUTH = DEG_TO_M_SEC (X)
C = '+000:02:00.00'
D_LON = READ_COOR_TEXT (C)
D_LAT = READ_COOR_TEXT (C)
SOUTH = SOUTH+D_LAT

DO I = 1,ILON
   AALON(I) = -180.0 + (FLOAT(I)-1.0)*2.0/60.0
END DO
DO J = 1,ILAT
   AALAT(J) = 90.0 - (FLOAT(J) -1.0)*2.0/60.0
END DO

! ---------------------------------------------------------------------------- !
!
!     2. READ INPUT DATA etopo2_2006apr.dat
!        ----------------------------------

ALLOCATE (IDEPTH(ILON,ILAT))

OPEN (IU08,FILE=FILE08,STATUS='OLD',IOSTAT=IOS)
IF (IOS.NE.0) THEN
   WRITE (IU06,*) ' ********************************************************'
   WRITE (IU06,*) ' *                                                      *'
   WRITE (IU06,*) ' *          FATAL  ERROR IN SUB. ETOPO_INPUT.           *'
   WRITE (IU06,*) ' *          =================================           *'
   WRITE (IU06,*) ' *                                                      *'
   WRITE (IU06,*) ' * ETOPO DEPTH DATA FILE DOES NOT EXIST.                *'
   WRITE (IU06,*) ' + FILE NAME IS: ', FILE08
   WRITE (IU06,*) ' *                                                      *'
   WRITE (IU06,*) ' *           PROGRAM ABORTS  PROGRAM ABORTS             *'
   WRITE (IU06,*) ' *                                                      *'
   WRITE (IU06,*) ' ********************************************************'
   CALL ABORT1
END IF

NTOT = ILON*ILAT
DO IJ = 1, NTOT, 15
   READ (IU08,'(15I6)',END=101) (IDUM(IC),IC=1,15)
   DO IC = 1,15
      ICOUNT = IJ+IC-1
      J = ICOUNT/ILON+1
      I = MOD(ICOUNT,ILON)
      IF (I.EQ.0) THEN
         I = ILON
         J = J-1
      ENDIF
      IDEPTH(I,J) = IDUM(IC)
   END DO
END DO
101   CONTINUE

CLOSE (UNIT=IU08)

! ---------------------------------------------------------------------------- !
!
!     3. READ REFERENCE LEVEL DEFINITION FILE (DEFINING REGIONS WHERE
!        MEAN SEA LEVEL SHOULD NOT BE THE REFERENCE LEVEL)
!        -------------------------------------------------------------

FILE_REF = 'reference_levels'

OPEN (IU08,FILE=FILE_REF,STATUS='OLD',IOSTAT=IOS)

IF (IOS.NE.0) THEN
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   WRITE (IU06,*) ' +                                                      +'
   WRITE (IU06,*) ' +           WARNING ERROR SUB. ETOPO_INPUT.            +'
   WRITE (IU06,*) ' +           ===============================            +'
   WRITE (IU06,*) ' +                                                      +'
   WRITE (IU06,*) ' +     REFERENCE LEVEL DEFINITION FILE DOES NOT EXIST.  +'
   WRITE (IU06,*) ' +     FILE NAME IS: ', FILE_REF
   WRITE (IU06,*) ' +                                                      +'
   WRITE (IU06,*) ' +                 MODEL CONTINUES                      +'
   WRITE (IU06,*) ' +              WITHOUT CORRECTIONS.                    +'
   WRITE (IU06,*) ' +                                                      +'
   WRITE (IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   NREFERENCE = 0
ELSE
   DO IR=1,NREF
      READ (IU08,*,END=1000,ERR=1000) XINF(IR), YINF(IR), XSUP(IR), YSUP(IR), &
&                               LEVEL(IR), NDEPTH(IR), LOCATION(IR)
   ENDDO
1000  NREFERENCE = IR-1
   CLOSE (UNIT=IU08)
END IF
WRITE(6,*) 'READ ',NREFERENCE,' NEW REFERENCE LEVELS'

! ---------------------------------------------------------------------------- !
!
!     4. merge referenc data into etotopo data.
!        --------------------------------------

DO IR = 1,NREFERENCE
   WRITE (*,*) XINF(IR), YINF(IR), XSUP(IR), YSUP(IR),                         &
&              LEVEL(IR), NDEPTH(IR), LOCATION(IR)
   DO J=1,ILAT
      YJ=AALAT(J)
      IF (YJ.GE.YINF(IR) .AND. YJ.LE.YSUP(IR)) THEN
         DO I=1,ILON
            XI=AALON(I)
            IF (XI.GE.XINF(IR) .AND. XI.LE.XSUP(IR)) THEN
               IF (IDEPTH(I,J).LE.LEVEL(IR)) THEN
                  IF (NDEPTH(IR).NE.0) THEN
                     IDEPTH(I,J) = NDEPTH(IR)
                  ELSE
                     IDEPTH(I,J) = IDEPTH(I,J)-LEVEL(IR)
                  ENDIF
               ELSE
                  IDEPTH(I,J) = IDEPTH(I,J)-LEVEL(IR)
               ENDIF
            ENDIF
         ENDDO
      ENDIF
   ENDDO
ENDDO

! ---------------------------------------------------------------------------- !
!
!     5. SOUTH OF 64S ALL NON DEEP POINTS ARE SET TO LAND TO AVOID
!        THE PROBLEM WITH PERMANENT ICE SHEET.
!        ----------------------------------------------------------

DO J = 1,ILAT
   IF (AALAT(J).LE.-64.0) THEN
      DO I=1,ILON
         IF (IDEPTH(I,J).GE.-250) IDEPTH(I,J)=1
      ENDDO
   ENDIF
ENDDO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. TRANSFER TO MODEULE.                                                  !
!        --------------------                                                  !

ALLOCATE (RDEPTH(ILON-1,ILAT))
RDEPTH(1:ILON-1,1:ILAT) = -REAL(IDEPTH(1:ILON-1,ILAT:1:-1))
DEALLOCATE (IDEPTH)
EAST = EAST-D_LON
CALL SET_TOPOGRAPHY (ILON-1, ILAT, D_LON, D_LAT, SOUTH, NORTH, WEST, EAST,     &
&                    RDEPTH)

DEALLOCATE (RDEPTH)

END SUBROUTINE READ_TOPOGRAPHY

MODULE UTIL_TNORGWD_MOD

USE YOMNORGWD, ONLY : TNORGWD

INTERFACE SAVE
MODULE PROCEDURE SAVE_TNORGWD
END INTERFACE

INTERFACE LOAD
MODULE PROCEDURE LOAD_TNORGWD
END INTERFACE

INTERFACE COPY
MODULE PROCEDURE COPY_TNORGWD
END INTERFACE

INTERFACE SIZE
MODULE PROCEDURE SIZE_TNORGWD
END INTERFACE


CONTAINS

SUBROUTINE SAVE_TNORGWD (KLUN, YD)

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (TNORGWD), INTENT (IN) :: YD

WRITE (KLUN) YD%NORGWD_SCHEME
WRITE (KLUN) YD%NORGWD_PRMAX
WRITE (KLUN) YD%NORGWD_DZ
WRITE (KLUN) YD%NORGWD_PTROPO
WRITE (KLUN) YD%NORGWD_NTROPO
WRITE (KLUN) YD%NORGWD_RUWMAX
WRITE (KLUN) YD%NORGWD_SAT
WRITE (KLUN) YD%NORGWD_RDISS
WRITE (KLUN) YD%NORGWD_DELTAT
WRITE (KLUN) YD%NORGWD_KMIN
WRITE (KLUN) YD%NORGWD_KMAX
WRITE (KLUN) YD%NORGWD_CMIN
WRITE (KLUN) YD%NORGWD_CMAX
WRITE (KLUN) YD%NORGWD_PLAUNCH
WRITE (KLUN) YD%NORGWD_NLAUNCH
WRITE (KLUN) YD%NORGWD_PNOVERDIF
WRITE (KLUN) YD%NORGWD_NNOVERDIF
WRITE (KLUN) YD%NORGWD_DZFRON
WRITE (KLUN) YD%NORGWD_GFRON
WRITE (KLUN) YD%NORGWD_GB
END SUBROUTINE
SUBROUTINE LOAD_TNORGWD (KLUN, YD)

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (TNORGWD), INTENT (OUT) :: YD

READ (KLUN) YD%NORGWD_SCHEME
READ (KLUN) YD%NORGWD_PRMAX
READ (KLUN) YD%NORGWD_DZ
READ (KLUN) YD%NORGWD_PTROPO
READ (KLUN) YD%NORGWD_NTROPO
READ (KLUN) YD%NORGWD_RUWMAX
READ (KLUN) YD%NORGWD_SAT
READ (KLUN) YD%NORGWD_RDISS
READ (KLUN) YD%NORGWD_DELTAT
READ (KLUN) YD%NORGWD_KMIN
READ (KLUN) YD%NORGWD_KMAX
READ (KLUN) YD%NORGWD_CMIN
READ (KLUN) YD%NORGWD_CMAX
READ (KLUN) YD%NORGWD_PLAUNCH
READ (KLUN) YD%NORGWD_NLAUNCH
READ (KLUN) YD%NORGWD_PNOVERDIF
READ (KLUN) YD%NORGWD_NNOVERDIF
READ (KLUN) YD%NORGWD_DZFRON
READ (KLUN) YD%NORGWD_GFRON
READ (KLUN) YD%NORGWD_GB
END SUBROUTINE

SUBROUTINE COPY_TNORGWD (YD, LDCREATED)

IMPLICIT NONE
TYPE (TNORGWD), INTENT (IN) :: YD
LOGICAL, OPTIONAL, INTENT (IN) :: LDCREATED
LOGICAL :: LLCREATED

LLCREATED = .FALSE.
IF (PRESENT (LDCREATED)) THEN
  LLCREATED = LDCREATED
ENDIF
IF (.NOT. LLCREATED) THEN
  !$acc enter data create (YD)
  !$acc update device (YD)
ENDIF




















END SUBROUTINE
INTEGER*8 FUNCTION SIZE_TNORGWD (YD, CDPATH, LDPRINT) RESULT (KSIZE)

IMPLICIT NONE
TYPE (TNORGWD),     INTENT (IN) :: YD
CHARACTER(LEN=*), INTENT (IN) :: CDPATH
LOGICAL,          INTENT (IN) :: LDPRINT
INTEGER*8 :: ISIZE, JSIZE

KSIZE = 0
ISIZE = KIND (YD%NORGWD_SCHEME) * LEN (YD%NORGWD_SCHEME)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_SCHEME'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_PRMAX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_PRMAX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_DZ)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_DZ'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_PTROPO)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_PTROPO'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_NTROPO)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_NTROPO'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_RUWMAX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_RUWMAX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_SAT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_SAT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_RDISS)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_RDISS'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_DELTAT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_DELTAT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_KMIN)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_KMIN'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_KMAX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_KMAX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_CMIN)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_CMIN'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_CMAX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_CMAX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_PLAUNCH)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_PLAUNCH'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_NLAUNCH)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_NLAUNCH'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_PNOVERDIF)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_PNOVERDIF'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_NNOVERDIF)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_NNOVERDIF'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_DZFRON)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_DZFRON'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_GFRON)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_GFRON'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NORGWD_GB)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NORGWD_GB'
ENDIF
KSIZE = KSIZE + ISIZE
END FUNCTION

END MODULE

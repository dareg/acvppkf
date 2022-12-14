MODULE UTIL_MODEL_PHYSICS_MF_TYPE_MOD

USE MODEL_PHYSICS_MF_MOD, ONLY : MODEL_PHYSICS_MF_TYPE

INTERFACE SAVE
MODULE PROCEDURE SAVE_MODEL_PHYSICS_MF_TYPE
END INTERFACE

INTERFACE LOAD
MODULE PROCEDURE LOAD_MODEL_PHYSICS_MF_TYPE
END INTERFACE

INTERFACE COPY
MODULE PROCEDURE COPY_MODEL_PHYSICS_MF_TYPE
END INTERFACE

INTERFACE SIZE
MODULE PROCEDURE SIZE_MODEL_PHYSICS_MF_TYPE
END INTERFACE


CONTAINS

SUBROUTINE SAVE_MODEL_PHYSICS_MF_TYPE (KLUN, YD)
USE UTIL_SL_STRUCT_MOD
USE UTIL_TARPHY_MOD
USE UTIL_TCVMNH_MOD
USE UTIL_TLOUIS_MOD
USE UTIL_TMSE_MOD
USE UTIL_TNORGWD_MOD
USE UTIL_TPARAR_MOD
USE UTIL_TPHY_MOD
USE UTIL_TPHY0_MOD
USE UTIL_TPHY1_MOD
USE UTIL_TPHY2_MOD
USE UTIL_TPHY3_MOD
USE UTIL_TPHYDS_MOD
USE UTIL_TSIMPHL_MOD
USE UTIL_TTOPH_MOD
USE UTIL_TVDOZ_MOD
IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (MODEL_PHYSICS_MF_TYPE), INTENT (IN) :: YD

CALL SAVE (KLUN, YD%YRPHY)
CALL SAVE (KLUN, YD%YRPHY0)
CALL SAVE (KLUN, YD%YRPHY1)
CALL SAVE (KLUN, YD%YRPHY2)
CALL SAVE (KLUN, YD%YRPHY3)
CALL SAVE (KLUN, YD%YRPHYDS)
CALL SAVE (KLUN, YD%YRCVMNH)
CALL SAVE (KLUN, YD%YRTOPH)
CALL SAVE (KLUN, YD%YRVDOZ)
CALL SAVE (KLUN, YD%YRSIMPHL)
CALL SAVE (KLUN, YD%YRARPHY)
CALL SAVE (KLUN, YD%YRPARAR)
CALL SAVE (KLUN, YD%YRMSE)
CALL SAVE (KLUN, YD%YRLOUIS)
CALL SAVE (KLUN, YD%YRNORGWD)
CALL SAVE (KLUN, YD%YRGR)
END SUBROUTINE
SUBROUTINE LOAD_MODEL_PHYSICS_MF_TYPE (KLUN, YD)
USE UTIL_SL_STRUCT_MOD
USE UTIL_TARPHY_MOD
USE UTIL_TCVMNH_MOD
USE UTIL_TLOUIS_MOD
USE UTIL_TMSE_MOD
USE UTIL_TNORGWD_MOD
USE UTIL_TPARAR_MOD
USE UTIL_TPHY_MOD
USE UTIL_TPHY0_MOD
USE UTIL_TPHY1_MOD
USE UTIL_TPHY2_MOD
USE UTIL_TPHY3_MOD
USE UTIL_TPHYDS_MOD
USE UTIL_TSIMPHL_MOD
USE UTIL_TTOPH_MOD
USE UTIL_TVDOZ_MOD

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (MODEL_PHYSICS_MF_TYPE), INTENT (OUT) :: YD

CALL LOAD (KLUN, YD%YRPHY)
CALL LOAD (KLUN, YD%YRPHY0)
CALL LOAD (KLUN, YD%YRPHY1)
CALL LOAD (KLUN, YD%YRPHY2)
CALL LOAD (KLUN, YD%YRPHY3)
CALL LOAD (KLUN, YD%YRPHYDS)
CALL LOAD (KLUN, YD%YRCVMNH)
CALL LOAD (KLUN, YD%YRTOPH)
CALL LOAD (KLUN, YD%YRVDOZ)
CALL LOAD (KLUN, YD%YRSIMPHL)
CALL LOAD (KLUN, YD%YRARPHY)
CALL LOAD (KLUN, YD%YRPARAR)
CALL LOAD (KLUN, YD%YRMSE)
CALL LOAD (KLUN, YD%YRLOUIS)
CALL LOAD (KLUN, YD%YRNORGWD)
CALL LOAD (KLUN, YD%YRGR)
END SUBROUTINE

SUBROUTINE COPY_MODEL_PHYSICS_MF_TYPE (YD, LDCREATED)
USE UTIL_SL_STRUCT_MOD
USE UTIL_TARPHY_MOD
USE UTIL_TCVMNH_MOD
USE UTIL_TLOUIS_MOD
USE UTIL_TMSE_MOD
USE UTIL_TNORGWD_MOD
USE UTIL_TPARAR_MOD
USE UTIL_TPHY_MOD
USE UTIL_TPHY0_MOD
USE UTIL_TPHY1_MOD
USE UTIL_TPHY2_MOD
USE UTIL_TPHY3_MOD
USE UTIL_TPHYDS_MOD
USE UTIL_TSIMPHL_MOD
USE UTIL_TTOPH_MOD
USE UTIL_TVDOZ_MOD
IMPLICIT NONE
TYPE (MODEL_PHYSICS_MF_TYPE), INTENT (IN) :: YD
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
CALL COPY (YD%YRPHY, LDCREATED=.TRUE.)

CALL COPY (YD%YRPHY0, LDCREATED=.TRUE.)

CALL COPY (YD%YRPHY1, LDCREATED=.TRUE.)

CALL COPY (YD%YRPHY2, LDCREATED=.TRUE.)

CALL COPY (YD%YRPHY3, LDCREATED=.TRUE.)

CALL COPY (YD%YRPHYDS, LDCREATED=.TRUE.)

CALL COPY (YD%YRCVMNH, LDCREATED=.TRUE.)

CALL COPY (YD%YRTOPH, LDCREATED=.TRUE.)

CALL COPY (YD%YRVDOZ, LDCREATED=.TRUE.)

CALL COPY (YD%YRSIMPHL, LDCREATED=.TRUE.)

CALL COPY (YD%YRARPHY, LDCREATED=.TRUE.)

CALL COPY (YD%YRPARAR, LDCREATED=.TRUE.)

CALL COPY (YD%YRMSE, LDCREATED=.TRUE.)

CALL COPY (YD%YRLOUIS, LDCREATED=.TRUE.)

CALL COPY (YD%YRNORGWD, LDCREATED=.TRUE.)

CALL COPY (YD%YRGR, LDCREATED=.TRUE.)

END SUBROUTINE
INTEGER*8 FUNCTION SIZE_MODEL_PHYSICS_MF_TYPE (YD, CDPATH, LDPRINT) RESULT (KSIZE)
USE UTIL_SL_STRUCT_MOD
USE UTIL_TARPHY_MOD
USE UTIL_TCVMNH_MOD
USE UTIL_TLOUIS_MOD
USE UTIL_TMSE_MOD
USE UTIL_TNORGWD_MOD
USE UTIL_TPARAR_MOD
USE UTIL_TPHY_MOD
USE UTIL_TPHY0_MOD
USE UTIL_TPHY1_MOD
USE UTIL_TPHY2_MOD
USE UTIL_TPHY3_MOD
USE UTIL_TPHYDS_MOD
USE UTIL_TSIMPHL_MOD
USE UTIL_TTOPH_MOD
USE UTIL_TVDOZ_MOD
IMPLICIT NONE
TYPE (MODEL_PHYSICS_MF_TYPE),     INTENT (IN) :: YD
CHARACTER(LEN=*), INTENT (IN) :: CDPATH
LOGICAL,          INTENT (IN) :: LDPRINT
INTEGER*8 :: ISIZE, JSIZE

KSIZE = 0
JSIZE = 0
ISIZE = SIZE (YD%YRPHY, CDPATH//'%YRPHY', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRPHY'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRPHY0, CDPATH//'%YRPHY0', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRPHY0'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRPHY1, CDPATH//'%YRPHY1', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRPHY1'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRPHY2, CDPATH//'%YRPHY2', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRPHY2'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRPHY3, CDPATH//'%YRPHY3', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRPHY3'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRPHYDS, CDPATH//'%YRPHYDS', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRPHYDS'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRCVMNH, CDPATH//'%YRCVMNH', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRCVMNH'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRTOPH, CDPATH//'%YRTOPH', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRTOPH'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRVDOZ, CDPATH//'%YRVDOZ', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRVDOZ'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRSIMPHL, CDPATH//'%YRSIMPHL', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRSIMPHL'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRARPHY, CDPATH//'%YRARPHY', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRARPHY'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRPARAR, CDPATH//'%YRPARAR', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRPARAR'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRMSE, CDPATH//'%YRMSE', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRMSE'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRLOUIS, CDPATH//'%YRLOUIS', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRLOUIS'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRNORGWD, CDPATH//'%YRNORGWD', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRNORGWD'
ENDIF
JSIZE = 0
ISIZE = SIZE (YD%YRGR, CDPATH//'%YRGR', .FALSE.)
JSIZE = JSIZE + ISIZE
KSIZE = KSIZE + ISIZE
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') JSIZE
  WRITE (*, *) TRIM (CDPATH)//'%YRGR'
ENDIF
END FUNCTION

END MODULE

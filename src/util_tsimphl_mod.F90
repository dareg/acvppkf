MODULE UTIL_TSIMPHL_MOD

USE YOMSIMPHL, ONLY : TSIMPHL

INTERFACE SAVE
MODULE PROCEDURE SAVE_TSIMPHL
END INTERFACE

INTERFACE LOAD
MODULE PROCEDURE LOAD_TSIMPHL
END INTERFACE

INTERFACE COPY
MODULE PROCEDURE COPY_TSIMPHL
END INTERFACE

INTERFACE SIZE
MODULE PROCEDURE SIZE_TSIMPHL
END INTERFACE


CONTAINS

SUBROUTINE SAVE_TSIMPHL (KLUN, YD)

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (TSIMPHL), INTENT (IN) :: YD

WRITE (KLUN) YD%LSIMPH
WRITE (KLUN) YD%LTRAJPS
WRITE (KLUN) YD%LTRAJPST
WRITE (KLUN) YD%LSMOOTHD
WRITE (KLUN) YD%LSMOOTHA
WRITE (KLUN) YD%LSMOOTHB
WRITE (KLUN) YD%LCVRASP
WRITE (KLUN) YD%LGWDSP
WRITE (KLUN) YD%LRAYSP
WRITE (KLUN) YD%LSTRASP
WRITE (KLUN) YD%LVDIFSP
WRITE (KLUN) YD%LVDIFSPNL
WRITE (KLUN) YD%LRRMESSP
WRITE (KLUN) YD%LCLOUDS
WRITE (KLUN) YD%LGWDSPNL
WRITE (KLUN) YD%LSTRASPN
WRITE (KLUN) YD%LPROCLDTL
WRITE (KLUN) YD%LMELTTL
WRITE (KLUN) YD%LMELTNL
WRITE (KLUN) YD%LMICROTL
WRITE (KLUN) YD%LTRAJRAIN
WRITE (KLUN) YD%LTRAJCOND
WRITE (KLUN) YD%LNEBCVPPKF
WRITE (KLUN) YD%LCOLLECTL
WRITE (KLUN) YD%LEVAPTL
WRITE (KLUN) YD%LSMOOTHEVP
WRITE (KLUN) YD%LIGELREPRO
WRITE (KLUN) YD%LCVRASBM
WRITE (KLUN) YD%LCONSENTH
WRITE (KLUN) YD%LAPPROXCONV
WRITE (KLUN) YD%RHCRIT1S
WRITE (KLUN) YD%RHCRIT2S
WRITE (KLUN) YD%TADJ
WRITE (KLUN) YD%RMINEVP
WRITE (KLUN) YD%DELTAH
WRITE (KLUN) YD%RMODULQCPROG
END SUBROUTINE
SUBROUTINE LOAD_TSIMPHL (KLUN, YD)

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (TSIMPHL), INTENT (OUT) :: YD

READ (KLUN) YD%LSIMPH
READ (KLUN) YD%LTRAJPS
READ (KLUN) YD%LTRAJPST
READ (KLUN) YD%LSMOOTHD
READ (KLUN) YD%LSMOOTHA
READ (KLUN) YD%LSMOOTHB
READ (KLUN) YD%LCVRASP
READ (KLUN) YD%LGWDSP
READ (KLUN) YD%LRAYSP
READ (KLUN) YD%LSTRASP
READ (KLUN) YD%LVDIFSP
READ (KLUN) YD%LVDIFSPNL
READ (KLUN) YD%LRRMESSP
READ (KLUN) YD%LCLOUDS
READ (KLUN) YD%LGWDSPNL
READ (KLUN) YD%LSTRASPN
READ (KLUN) YD%LPROCLDTL
READ (KLUN) YD%LMELTTL
READ (KLUN) YD%LMELTNL
READ (KLUN) YD%LMICROTL
READ (KLUN) YD%LTRAJRAIN
READ (KLUN) YD%LTRAJCOND
READ (KLUN) YD%LNEBCVPPKF
READ (KLUN) YD%LCOLLECTL
READ (KLUN) YD%LEVAPTL
READ (KLUN) YD%LSMOOTHEVP
READ (KLUN) YD%LIGELREPRO
READ (KLUN) YD%LCVRASBM
READ (KLUN) YD%LCONSENTH
READ (KLUN) YD%LAPPROXCONV
READ (KLUN) YD%RHCRIT1S
READ (KLUN) YD%RHCRIT2S
READ (KLUN) YD%TADJ
READ (KLUN) YD%RMINEVP
READ (KLUN) YD%DELTAH
READ (KLUN) YD%RMODULQCPROG
END SUBROUTINE

SUBROUTINE COPY_TSIMPHL (YD, LDCREATED)

IMPLICIT NONE
TYPE (TSIMPHL), INTENT (IN) :: YD
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
INTEGER*8 FUNCTION SIZE_TSIMPHL (YD, CDPATH, LDPRINT) RESULT (KSIZE)

IMPLICIT NONE
TYPE (TSIMPHL),     INTENT (IN) :: YD
CHARACTER(LEN=*), INTENT (IN) :: CDPATH
LOGICAL,          INTENT (IN) :: LDPRINT
INTEGER*8 :: ISIZE, JSIZE

KSIZE = 0
ISIZE = KIND (YD%LSIMPH)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LSIMPH'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LTRAJPS)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LTRAJPS'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LTRAJPST)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LTRAJPST'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LSMOOTHD)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LSMOOTHD'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LSMOOTHA)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LSMOOTHA'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LSMOOTHB)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LSMOOTHB'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LCVRASP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LCVRASP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LGWDSP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LGWDSP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LRAYSP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LRAYSP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LSTRASP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LSTRASP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LVDIFSP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LVDIFSP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LVDIFSPNL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LVDIFSPNL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LRRMESSP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LRRMESSP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LCLOUDS)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LCLOUDS'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LGWDSPNL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LGWDSPNL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LSTRASPN)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LSTRASPN'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LPROCLDTL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LPROCLDTL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LMELTTL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LMELTTL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LMELTNL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LMELTNL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LMICROTL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LMICROTL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LTRAJRAIN)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LTRAJRAIN'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LTRAJCOND)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LTRAJCOND'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LNEBCVPPKF)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LNEBCVPPKF'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LCOLLECTL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LCOLLECTL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LEVAPTL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LEVAPTL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LSMOOTHEVP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LSMOOTHEVP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LIGELREPRO)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LIGELREPRO'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LCVRASBM)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LCVRASBM'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LCONSENTH)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LCONSENTH'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LAPPROXCONV)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LAPPROXCONV'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RHCRIT1S)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RHCRIT1S'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RHCRIT2S)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RHCRIT2S'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%TADJ)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%TADJ'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RMINEVP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RMINEVP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%DELTAH)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%DELTAH'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RMODULQCPROG)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RMODULQCPROG'
ENDIF
KSIZE = KSIZE + ISIZE
END FUNCTION

END MODULE

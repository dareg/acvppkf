MODULE UTIL_TPHY1_MOD

USE YOMPHY1, ONLY : TPHY1

INTERFACE SAVE
MODULE PROCEDURE SAVE_TPHY1
END INTERFACE

INTERFACE LOAD
MODULE PROCEDURE LOAD_TPHY1
END INTERFACE

INTERFACE COPY
MODULE PROCEDURE COPY_TPHY1
END INTERFACE

INTERFACE SIZE
MODULE PROCEDURE SIZE_TPHY1
END INTERFACE


CONTAINS

SUBROUTINE SAVE_TPHY1 (KLUN, YD)

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (TPHY1), INTENT (IN) :: YD

WRITE (KLUN) YD%GF3
WRITE (KLUN) YD%GF4
WRITE (KLUN) YD%TREF4
WRITE (KLUN) YD%RCTVEG
WRITE (KLUN) YD%RGL
WRITE (KLUN) YD%SODELX
WRITE (KLUN) YD%GCZ0H
WRITE (KLUN) YD%ALBGLA
WRITE (KLUN) YD%ALBMAX
WRITE (KLUN) YD%ALBMER
WRITE (KLUN) YD%ALBMED
WRITE (KLUN) YD%ALBMIN
WRITE (KLUN) YD%ALCRIN
WRITE (KLUN) YD%ALRCN1
WRITE (KLUN) YD%ALRCN2
WRITE (KLUN) YD%EA
WRITE (KLUN) YD%EC2REF
WRITE (KLUN) YD%EMCRIN
WRITE (KLUN) YD%EMMGLA
WRITE (KLUN) YD%EMMMER
WRITE (KLUN) YD%EWFC
WRITE (KLUN) YD%EWWILT
WRITE (KLUN) YD%GA
WRITE (KLUN) YD%GC1
WRITE (KLUN) YD%GC1S1
WRITE (KLUN) YD%GC1S2
WRITE (KLUN) YD%GC1S3
WRITE (KLUN) YD%GC1S4
WRITE (KLUN) YD%GC1Y1
WRITE (KLUN) YD%GTSVAP
WRITE (KLUN) YD%GVEGMX
WRITE (KLUN) YD%GLAIMX
WRITE (KLUN) YD%GNEIMX
WRITE (KLUN) YD%GWPIMX
WRITE (KLUN) YD%GCGEL
WRITE (KLUN) YD%GC2
WRITE (KLUN) YD%GC2REF
WRITE (KLUN) YD%GC3
WRITE (KLUN) YD%GC31
WRITE (KLUN) YD%GC32
WRITE (KLUN) YD%GCONV
WRITE (KLUN) YD%GF1
WRITE (KLUN) YD%GWFC
WRITE (KLUN) YD%GWLEX
WRITE (KLUN) YD%GWLMX
WRITE (KLUN) YD%GWWILT
WRITE (KLUN) YD%G1B
WRITE (KLUN) YD%G1CGSAT
WRITE (KLUN) YD%G1C1SAT
WRITE (KLUN) YD%G1P
WRITE (KLUN) YD%G1WSAT
WRITE (KLUN) YD%G2B
WRITE (KLUN) YD%G2CGSAT
WRITE (KLUN) YD%G2C1SAT
WRITE (KLUN) YD%G2P
WRITE (KLUN) YD%G2WSAT
WRITE (KLUN) YD%G3CGSAT
WRITE (KLUN) YD%GSNC1
WRITE (KLUN) YD%GSNC2
WRITE (KLUN) YD%HSOL
WRITE (KLUN) YD%HSOLIWR
WRITE (KLUN) YD%HSOLIT0
WRITE (KLUN) YD%OMTPRO
WRITE (KLUN) YD%OMWPRO
WRITE (KLUN) YD%RC1MAX
WRITE (KLUN) YD%RCTGLA
WRITE (KLUN) YD%RCGMAX
WRITE (KLUN) YD%RD1
WRITE (KLUN) YD%RD2GLA
WRITE (KLUN) YD%RD2MER
WRITE (KLUN) YD%RHOMAX
WRITE (KLUN) YD%RHOMIN
WRITE (KLUN) YD%RSMAX
WRITE (KLUN) YD%RTINER
WRITE (KLUN) YD%RZ0GLA
WRITE (KLUN) YD%RZ0MER
WRITE (KLUN) YD%RZHZ0G
WRITE (KLUN) YD%RZHZ0M
WRITE (KLUN) YD%TMERGL
WRITE (KLUN) YD%TOEXP
WRITE (KLUN) YD%TOLIN
WRITE (KLUN) YD%WCRIN
WRITE (KLUN) YD%WCRINC
WRITE (KLUN) YD%WCRING
WRITE (KLUN) YD%WNEW
WRITE (KLUN) YD%WPMX
WRITE (KLUN) YD%WSMX
WRITE (KLUN) YD%XCRINR
WRITE (KLUN) YD%XCRINV
WRITE (KLUN) YD%LALBMERCLIM
WRITE (KLUN) YD%LIMC
WRITE (KLUN) YD%LIMW
WRITE (KLUN) YD%LC1VAP
WRITE (KLUN) YD%LCLS_HS
WRITE (KLUN) YD%NTVGLA
WRITE (KLUN) YD%NTVMER
WRITE (KLUN) YD%GCGELS
WRITE (KLUN) YD%GVEGMXS
WRITE (KLUN) YD%GLAIMXS
WRITE (KLUN) YD%GNEIMXS
WRITE (KLUN) YD%ALB1
WRITE (KLUN) YD%ALB2
WRITE (KLUN) YD%RLAIMX
WRITE (KLUN) YD%RLAI
WRITE (KLUN) YD%ACLS_HS
WRITE (KLUN) YD%NCHSP
END SUBROUTINE
SUBROUTINE LOAD_TPHY1 (KLUN, YD)

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (TPHY1), INTENT (OUT) :: YD

READ (KLUN) YD%GF3
READ (KLUN) YD%GF4
READ (KLUN) YD%TREF4
READ (KLUN) YD%RCTVEG
READ (KLUN) YD%RGL
READ (KLUN) YD%SODELX
READ (KLUN) YD%GCZ0H
READ (KLUN) YD%ALBGLA
READ (KLUN) YD%ALBMAX
READ (KLUN) YD%ALBMER
READ (KLUN) YD%ALBMED
READ (KLUN) YD%ALBMIN
READ (KLUN) YD%ALCRIN
READ (KLUN) YD%ALRCN1
READ (KLUN) YD%ALRCN2
READ (KLUN) YD%EA
READ (KLUN) YD%EC2REF
READ (KLUN) YD%EMCRIN
READ (KLUN) YD%EMMGLA
READ (KLUN) YD%EMMMER
READ (KLUN) YD%EWFC
READ (KLUN) YD%EWWILT
READ (KLUN) YD%GA
READ (KLUN) YD%GC1
READ (KLUN) YD%GC1S1
READ (KLUN) YD%GC1S2
READ (KLUN) YD%GC1S3
READ (KLUN) YD%GC1S4
READ (KLUN) YD%GC1Y1
READ (KLUN) YD%GTSVAP
READ (KLUN) YD%GVEGMX
READ (KLUN) YD%GLAIMX
READ (KLUN) YD%GNEIMX
READ (KLUN) YD%GWPIMX
READ (KLUN) YD%GCGEL
READ (KLUN) YD%GC2
READ (KLUN) YD%GC2REF
READ (KLUN) YD%GC3
READ (KLUN) YD%GC31
READ (KLUN) YD%GC32
READ (KLUN) YD%GCONV
READ (KLUN) YD%GF1
READ (KLUN) YD%GWFC
READ (KLUN) YD%GWLEX
READ (KLUN) YD%GWLMX
READ (KLUN) YD%GWWILT
READ (KLUN) YD%G1B
READ (KLUN) YD%G1CGSAT
READ (KLUN) YD%G1C1SAT
READ (KLUN) YD%G1P
READ (KLUN) YD%G1WSAT
READ (KLUN) YD%G2B
READ (KLUN) YD%G2CGSAT
READ (KLUN) YD%G2C1SAT
READ (KLUN) YD%G2P
READ (KLUN) YD%G2WSAT
READ (KLUN) YD%G3CGSAT
READ (KLUN) YD%GSNC1
READ (KLUN) YD%GSNC2
READ (KLUN) YD%HSOL
READ (KLUN) YD%HSOLIWR
READ (KLUN) YD%HSOLIT0
READ (KLUN) YD%OMTPRO
READ (KLUN) YD%OMWPRO
READ (KLUN) YD%RC1MAX
READ (KLUN) YD%RCTGLA
READ (KLUN) YD%RCGMAX
READ (KLUN) YD%RD1
READ (KLUN) YD%RD2GLA
READ (KLUN) YD%RD2MER
READ (KLUN) YD%RHOMAX
READ (KLUN) YD%RHOMIN
READ (KLUN) YD%RSMAX
READ (KLUN) YD%RTINER
READ (KLUN) YD%RZ0GLA
READ (KLUN) YD%RZ0MER
READ (KLUN) YD%RZHZ0G
READ (KLUN) YD%RZHZ0M
READ (KLUN) YD%TMERGL
READ (KLUN) YD%TOEXP
READ (KLUN) YD%TOLIN
READ (KLUN) YD%WCRIN
READ (KLUN) YD%WCRINC
READ (KLUN) YD%WCRING
READ (KLUN) YD%WNEW
READ (KLUN) YD%WPMX
READ (KLUN) YD%WSMX
READ (KLUN) YD%XCRINR
READ (KLUN) YD%XCRINV
READ (KLUN) YD%LALBMERCLIM
READ (KLUN) YD%LIMC
READ (KLUN) YD%LIMW
READ (KLUN) YD%LC1VAP
READ (KLUN) YD%LCLS_HS
READ (KLUN) YD%NTVGLA
READ (KLUN) YD%NTVMER
READ (KLUN) YD%GCGELS
READ (KLUN) YD%GVEGMXS
READ (KLUN) YD%GLAIMXS
READ (KLUN) YD%GNEIMXS
READ (KLUN) YD%ALB1
READ (KLUN) YD%ALB2
READ (KLUN) YD%RLAIMX
READ (KLUN) YD%RLAI
READ (KLUN) YD%ACLS_HS
READ (KLUN) YD%NCHSP
END SUBROUTINE

SUBROUTINE COPY_TPHY1 (YD, LDCREATED)

IMPLICIT NONE
TYPE (TPHY1), INTENT (IN) :: YD
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
INTEGER*8 FUNCTION SIZE_TPHY1 (YD, CDPATH, LDPRINT) RESULT (KSIZE)

IMPLICIT NONE
TYPE (TPHY1),     INTENT (IN) :: YD
CHARACTER(LEN=*), INTENT (IN) :: CDPATH
LOGICAL,          INTENT (IN) :: LDPRINT
INTEGER*8 :: ISIZE, JSIZE

KSIZE = 0
ISIZE = KIND (YD%GF3) * SIZE (YD%GF3)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GF3'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GF4) * SIZE (YD%GF4)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GF4'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%TREF4) * SIZE (YD%TREF4)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%TREF4'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RCTVEG) * SIZE (YD%RCTVEG)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RCTVEG'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RGL) * SIZE (YD%RGL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RGL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%SODELX) * SIZE (YD%SODELX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%SODELX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GCZ0H) * SIZE (YD%GCZ0H)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GCZ0H'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ALBGLA)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ALBGLA'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ALBMAX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ALBMAX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ALBMER)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ALBMER'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ALBMED)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ALBMED'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ALBMIN)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ALBMIN'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ALCRIN)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ALCRIN'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ALRCN1)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ALRCN1'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ALRCN2)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ALRCN2'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%EA)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%EA'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%EC2REF)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%EC2REF'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%EMCRIN)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%EMCRIN'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%EMMGLA)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%EMMGLA'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%EMMMER)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%EMMMER'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%EWFC)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%EWFC'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%EWWILT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%EWWILT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GA)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GA'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC1)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC1'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC1S1)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC1S1'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC1S2)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC1S2'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC1S3)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC1S3'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC1S4)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC1S4'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC1Y1)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC1Y1'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GTSVAP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GTSVAP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GVEGMX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GVEGMX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GLAIMX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GLAIMX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GNEIMX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GNEIMX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GWPIMX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GWPIMX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GCGEL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GCGEL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC2)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC2'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC2REF)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC2REF'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC3)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC3'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC31)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC31'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GC32)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GC32'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GCONV)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GCONV'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GF1)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GF1'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GWFC)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GWFC'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GWLEX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GWLEX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GWLMX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GWLMX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GWWILT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GWWILT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G1B)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G1B'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G1CGSAT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G1CGSAT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G1C1SAT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G1C1SAT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G1P)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G1P'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G1WSAT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G1WSAT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G2B)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G2B'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G2CGSAT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G2CGSAT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G2C1SAT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G2C1SAT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G2P)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G2P'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G2WSAT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G2WSAT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%G3CGSAT)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%G3CGSAT'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GSNC1)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GSNC1'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GSNC2)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GSNC2'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%HSOL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%HSOL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%HSOLIWR)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%HSOLIWR'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%HSOLIT0)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%HSOLIT0'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%OMTPRO)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%OMTPRO'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%OMWPRO)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%OMWPRO'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RC1MAX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RC1MAX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RCTGLA)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RCTGLA'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RCGMAX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RCGMAX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RD1)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RD1'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RD2GLA)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RD2GLA'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RD2MER)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RD2MER'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RHOMAX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RHOMAX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RHOMIN)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RHOMIN'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RSMAX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RSMAX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RTINER)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RTINER'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RZ0GLA)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RZ0GLA'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RZ0MER)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RZ0MER'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RZHZ0G)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RZHZ0G'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RZHZ0M)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RZHZ0M'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%TMERGL)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%TMERGL'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%TOEXP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%TOEXP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%TOLIN)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%TOLIN'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%WCRIN)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%WCRIN'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%WCRINC)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%WCRINC'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%WCRING)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%WCRING'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%WNEW)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%WNEW'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%WPMX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%WPMX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%WSMX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%WSMX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%XCRINR)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%XCRINR'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%XCRINV)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%XCRINV'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LALBMERCLIM)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LALBMERCLIM'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LIMC)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LIMC'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LIMW)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LIMW'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LC1VAP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LC1VAP'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%LCLS_HS)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%LCLS_HS'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NTVGLA)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NTVGLA'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NTVMER)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NTVMER'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GCGELS)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GCGELS'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GVEGMXS)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GVEGMXS'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GLAIMXS)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GLAIMXS'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%GNEIMXS)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%GNEIMXS'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ALB1)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ALB1'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ALB2)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ALB2'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RLAIMX)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RLAIMX'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%RLAI)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%RLAI'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%ACLS_HS)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%ACLS_HS'
ENDIF
KSIZE = KSIZE + ISIZE
ISIZE = KIND (YD%NCHSP)
IF (LDPRINT) THEN
  WRITE (*, '(I10," ")', ADVANCE='NO') ISIZE
  WRITE (*, *) TRIM (CDPATH)//'%NCHSP'
ENDIF
KSIZE = KSIZE + ISIZE
END FUNCTION

END MODULE

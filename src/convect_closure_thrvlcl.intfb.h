INTERFACE

SUBROUTINE CONVECT_CLOSURE_THRVLCL( CVPEXT, CST, D,          &
PPRES, PTH, PRV, PZ, OWORK1,        &
PTHLCL, PRVLCL, PZLCL, PTLCL, PTELCL,&
KLCL, KDPL, KPBL )
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_CST, ONLY : CST_T
USE MODD_CONVPAREXT, ONLY : CONVPAREXT
USE MODD_CST, ONLY: CST_T
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
TYPE(CONVPAREXT),           INTENT(IN) :: CVPEXT
TYPE(CST_T),                INTENT(IN) :: CST
TYPE(DIMPHYEX_T),           INTENT(IN) :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PTH
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PRV
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PPRES
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PZ
INTEGER, DIMENSION(D%NIT),   INTENT(IN) :: KDPL
INTEGER, DIMENSION(D%NIT),   INTENT(IN) :: KPBL
LOGICAL, DIMENSION(D%NIT),   INTENT(IN) :: OWORK1
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PTHLCL
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PRVLCL
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PZLCL
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PTLCL
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PTELCL
INTEGER, DIMENSION(D%NIT),  INTENT(OUT):: KLCL
END SUBROUTINE CONVECT_CLOSURE_THRVLCL

END INTERFACE

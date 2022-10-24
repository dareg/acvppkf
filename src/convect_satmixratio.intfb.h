INTERFACE

SUBROUTINE CONVECT_SATMIXRATIO(CST, D, PPRES, PT, PEW, PLV, PLS, PCPH)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_CST, ONLY : CST_T
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
TYPE(CST_T),            INTENT(IN) :: CST
TYPE(DIMPHYEX_T),       INTENT(IN) :: D
REAL, DIMENSION(D%NIT),  INTENT(IN) :: PPRES
REAL, DIMENSION(D%NIT),  INTENT(IN) :: PT
REAL, DIMENSION(D%NIT),  INTENT(OUT):: PEW
REAL, DIMENSION(D%NIT),  INTENT(OUT):: PLV
REAL, DIMENSION(D%NIT),  INTENT(OUT):: PLS
REAL, DIMENSION(D%NIT),  INTENT(OUT):: PCPH
END SUBROUTINE CONVECT_SATMIXRATIO

END INTERFACE

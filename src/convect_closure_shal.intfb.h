INTERFACE

SUBROUTINE CONVECT_CLOSURE_SHAL( CVP_SHAL, CVPEXT, CST, D,         &
PPRES, PDPRES, PZ, PLMASS,           &
PTHL, PTH, PRW, PRC, PRI, OTRIG1,           &
PTHC, PRWC, PRCC, PRIC, PWSUB,              &
KLCL, KDPL, KPBL, KCTL,                     &
PUMF, PUER, PUDR, PUTHL, PURW,              &
PURC, PURI, PCAPE, PTIMEC, KFTSTEPS         )
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_CST, ONLY : CST_T
USE MODD_CONVPAR_SHAL, ONLY : CONVPAR_SHAL
USE MODD_CONVPAREXT, ONLY : CONVPAREXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
TYPE(CONVPAR_SHAL),        INTENT(IN) :: CVP_SHAL
TYPE(CONVPAREXT),          INTENT(IN) :: CVPEXT
TYPE(CST_T),               INTENT(IN) :: CST
TYPE(DIMPHYEX_T),          INTENT(IN) :: D
INTEGER, DIMENSION(D%NIT),  INTENT(IN) :: KLCL
INTEGER, DIMENSION(D%NIT),  INTENT(IN) :: KCTL
INTEGER, DIMENSION(D%NIT),  INTENT(IN) :: KDPL
INTEGER, DIMENSION(D%NIT),  INTENT(IN) :: KPBL
REAL, DIMENSION(D%NIT),  INTENT(INOUT) :: PTIMEC
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PTHL
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PTH
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PRW
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PRC
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PRI
LOGICAL, DIMENSION(D%NIT),  INTENT(IN) :: OTRIG1
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PPRES
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PDPRES
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PLMASS
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PZ
REAL, DIMENSION(D%NIT),     INTENT(IN)  :: PCAPE
INTEGER,                INTENT(OUT)   :: KFTSTEPS
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT):: PUMF
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT):: PUER
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT):: PUDR
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PUTHL
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PURW
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PURC
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PURI
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PTHC
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PRWC
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PRCC
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PRIC
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PWSUB
END SUBROUTINE CONVECT_CLOSURE_SHAL

END INTERFACE
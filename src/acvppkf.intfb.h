INTERFACE
SUBROUTINE ACVPPKF( YDCST, YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
 & PAPRSF, PAPHIF, PDELP, PR, PT, PQ,&
 & PQL, PQI, PU, PV, PVERVEL, PCP, PTKE,&
 & PDIFCQ, PDIFCS, PFCCQL, PFCCQN, PPRODTH,&
 & KNLAB, PQCPP, PNEBPP,&
 & KNND) 
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE YOMCST , ONLY : TCST
TYPE (TCST), INTENT (IN) :: YDCST
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
INTEGER(KIND=JPIM) ,INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM) ,INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM) ,INTENT(IN) :: KLON
INTEGER(KIND=JPIM) ,INTENT(IN) :: KTDIA
INTEGER(KIND=JPIM) ,INTENT(IN) :: KLEV
REAL(KIND=JPRB) ,INTENT(IN) :: PAPRSF (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PAPHIF (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PDELP (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PR (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PT (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQ (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQL (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQI (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PU (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PV (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PVERVEL(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PCP (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PTKE (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PDIFCQ (KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PDIFCS (KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCCQL (KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCCQN (KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PPRODTH(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PQCPP (KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PNEBPP (KLON,KLEV)
INTEGER(KIND=JPIM) ,INTENT(OUT) :: KNLAB (KLON,KLEV)
INTEGER(KIND=JPIM) ,INTENT(OUT) :: KNND (KLON)
END SUBROUTINE ACVPPKF
END INTERFACE
!OPTIONS XOPT(NOEVAL)
!-----------------------------------------------------------------
SUBROUTINE WRAP_ACVPPKF(NSTEPMAX, NPROMA, LAST_BLK_NPROMA, NFLEV, NBLKS)

!-----------------------------------------------------------------

!  Authors  : E. Bazile and P. Bechtold  (CNRM/GMAP et L.A.)

!-----------------------------------------------------------------

!  Modified : 
!  05/2002    phased with CONVECTION call for IFS/ECMWF 
!            (routine cucalln.F90 calling both Tiedtke convection scheme
!             and present scheme)
!             ouput of present scheme (updraft QL and QV) provides also
!             necessary parameters for Tiedtke prognostic cloud scheme
!  03/2002  P. Marquet.  new  ZFHMLTS, ZFHEVPP in CPFHPRS (for Lopez)
!  03/2002  P. Marquet.  new  LKFDEEP, LKFSHAL
!  03/2002  P. Marquet.  "call deep_convection" 
!                      > "call convection" (deep + shallow)
!  09/2006  E. Bazile : Appel de la routine de shallow convection d'AROME
!                       uniquement
!  04/2008  E. Bazile : calcul du terme de production thermique PPROTH
!  10/2008  Y. Bouteloup & F. Bouyssel : Correction of bugs in initialization
!  07/2009  E. Bazile : TKE en entree de KFB et W fct de W_conv
!  K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!  10/2009  F. Bouyssel : Limitation on maximal TKE value
!  02/2010  E. Bazile : Correction for W without TKE scheme.
!  04/2010  F. Bouyssel : Bug correction on KNLAB computation
!  09/2010  O. Spaniel : Bug correction in expression SQRT(MIN)
!  04/2011  F. Bouyssel : Correction of a jlon loop (kidia,kfdia)
!  12/2012  E. Bazile   : Modif of W_turb and qc and cc fct of mass flux.
!     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.

!  Peter.Bechtold@ecmwf.int

! Sequence de  routines :
! aplpar > acvppkf > convection_shal

! iv)  Momentum transport:
!      Option LLUVTRANS: c'est possible d'utiliser maintenant
!           mais pas encore bien teste. Donc par defaut mettre
!           LLUVTRANS=.FALSE.

! vi)   Traceurs passifs - chimie: 
!      Cette partie est utilisee uniquement dans MOCAGE et dans MESONH
!      Si on ne veut pas de transport de traceurs (ex. Ozone,CO) dans 
!      ARPEGE/ECMWF IFS, mettre tout siplement OCHTRANS=FALSE et KCH1=0 
!      (nombre de traceurs). PCH1 (traceur) et PCH1TEN (tendance 
!      convective du traceur) ont alors les dimensions
!      (KLON,KLEV,KCH1=0) qui ne prennent pas de place.

!-----------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMCST  , ONLY :  TCST
USE YOMCT3               , ONLY : NSTEP
USE UTIL_TCST_MOD
USE UTIL_MODEL_PHYSICS_MF_TYPE_MOD
!-----------------------------------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: NSTEPMAX, NPROMA, LAST_BLK_NPROMA, NFLEV, NBLKS

!TYPE (TCST), INTENT (IN) :: YDCST
!TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
!INTEGER(KIND=JPIM) ,INTENT(IN)    :: KIDIA 
!INTEGER(KIND=JPIM) ,INTENT(IN)    :: KFDIA 
!INTEGER(KIND=JPIM) ,INTENT(IN)    :: KLON
!INTEGER(KIND=JPIM) ,INTENT(IN)    :: KTDIA 
!INTEGER(KIND=JPIM) ,INTENT(IN)    :: KLEV
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PAPRSF (KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PAPHIF (KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PDELP  (KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PR     (KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PT     (KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PQ     (KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PQL    (KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PQI    (KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PU     (KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PV     (KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PCP    (KLON,KLEV)
!REAL(KIND=JPRB)    ,INTENT(IN)    :: PTKE   (KLON,KLEV)
!REAL(KIND=JPRB)    ,INTENT(INOUT) :: PDIFCQ (KLON,0:KLEV) 
!REAL(KIND=JPRB)    ,INTENT(INOUT) :: PDIFCS (KLON,0:KLEV) 
!REAL(KIND=JPRB)    ,INTENT(INOUT) :: PFCCQL (KLON,0:KLEV) 
!REAL(KIND=JPRB)    ,INTENT(INOUT) :: PFCCQN (KLON,0:KLEV) 
!REAL(KIND=JPRB)    ,INTENT(INOUT) :: PPRODTH(KLON,0:KLEV)
!REAL(KIND=JPRB)    ,INTENT(INOUT) :: PQCPP  (KLON,KLEV)
!REAL(KIND=JPRB)    ,INTENT(INOUT) :: PNEBPP (KLON,KLEV)
!INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KNLAB  (KLON,KLEV) 
!INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KNND   (KLON) 

INTEGER(KIND=JPIM) :: KIDIA 
INTEGER(KIND=JPIM) :: KFDIA 
INTEGER(KIND=JPIM) :: KLON
INTEGER(KIND=JPIM) :: KTDIA 
INTEGER(KIND=JPIM) :: KLEV
REAL(KIND=JPRB)    :: PAPRSF (KLON,KLEV) 
REAL(KIND=JPRB)    :: PAPHIF (KLON,KLEV) 
REAL(KIND=JPRB)    :: PDELP  (KLON,KLEV) 
REAL(KIND=JPRB)    :: PR     (KLON,KLEV) 
REAL(KIND=JPRB)    :: PT     (KLON,KLEV) 
REAL(KIND=JPRB)    :: PQ     (KLON,KLEV) 
REAL(KIND=JPRB)    :: PQL    (KLON,KLEV) 
REAL(KIND=JPRB)    :: PQI    (KLON,KLEV) 
REAL(KIND=JPRB)    :: PU     (KLON,KLEV) 
REAL(KIND=JPRB)    :: PV     (KLON,KLEV) 
REAL(KIND=JPRB)    :: PVERVEL(KLON,KLEV) 
REAL(KIND=JPRB)    :: PCP    (KLON,KLEV)
REAL(KIND=JPRB)    :: PTKE   (KLON,KLEV)
REAL(KIND=JPRB)    :: PDIFCQ (KLON,0:KLEV) 
REAL(KIND=JPRB)    :: PDIFCS (KLON,0:KLEV) 
REAL(KIND=JPRB)    :: PFCCQL (KLON,0:KLEV) 
REAL(KIND=JPRB)    :: PFCCQN (KLON,0:KLEV) 
REAL(KIND=JPRB)    :: PPRODTH(KLON,0:KLEV)
REAL(KIND=JPRB)    :: PQCPP  (KLON,KLEV)
REAL(KIND=JPRB)    :: PNEBPP (KLON,KLEV)
INTEGER(KIND=JPIM) :: KNLAB  (KLON,KLEV) 
INTEGER(KIND=JPIM) :: KNND   (KLON) 

INTEGER :: FH
INTEGER, PARAMETER :: CHECK = 123456789

#include "acvppkf.intfb.h"

CALL ACVPPKF(YDCST,YDML_PHY_MF, KIDIA, KFDIA, KLON, KTDIA, KLEV,  &
& PAPRSF, PAPHIF, PDELP,  &
& PR, PT, PQ, PQL, PQI, PU,                 &
& PV, PVERVEL, PCP, PTKE, &
& PDIFCQ, PDIFCS, PFCCQL, PFCCQN, PPRODTH, KNLAB, PQCPP, PNEBPP,                       &
& KNND)


END SUBROUTINE WRAP_ACVPPKF

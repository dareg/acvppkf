!OPTIONS XOPT(NOEVAL)
!-----------------------------------------------------------------
SUBROUTINE ACVPPKF( YDCST, YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV, &
 !-----------------------------------------------------------------
 ! - INPUT  2D .
 & PAPRSF, PAPHIF, PDELP, PR, PT, PQ, &
 & PQL, PQI, PU, PV, PVERVEL, PCP, PTKE, &
 ! - OUTPUT 2D .
 & PDIFCQ, PDIFCS, PFCCQL, PFCCQN, PPRODTH, &
 & KNLAB, PQCPP, PNEBPP,&
 ! - OUTPUT 1D .
 & KNND)

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
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST  , ONLY :  TCST
USE YOMLSFORC, ONLY : LMUSCLFA,NMUSCLFA
USE MODD_CONVPAR_SHAL, ONLY : LLSMOOTH, XA25, XATPERT, XAW, XBTPERT, XBW, XCDEPTH, XCDEPTH_D, XCRAD, XDTPERT, XENTR, &
& XNHGAM, XSTABC, XSTABT, XTFRZ1, XTFRZ2, XWTRIG, XZLCL, XZPBL
!-----------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KLON
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KTDIA 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KLEV
REAL(KIND=JPRB)    ,INTENT(IN)    :: PAPRSF (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PAPHIF (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PDELP  (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PR     (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PT     (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQ     (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQL    (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQI    (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PU     (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PV     (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PVERVEL(KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(IN)    :: PCP    (KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(IN)    :: PTKE   (KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PDIFCQ (KLON,0:KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PDIFCS (KLON,0:KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PFCCQL (KLON,0:KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PFCCQN (KLON,0:KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PPRODTH(KLON,0:KLEV)
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PQCPP  (KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PNEBPP (KLON,KLEV)
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KNLAB  (KLON,KLEV) 
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KNND   (KLON) 

!-----------------------------------------------------------------

LOGICAL :: LLREFRESH_ALL, LLDOWN, LLUVTRANS, LLOCHTRANS, LLCONDWT

REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZDTCONV, ZVMD, ZWMD, ZSMD, ZTDCP, ZEPS, ZDQCDT, ZDTLDT

INTEGER(KIND=JPIM) :: JLON, JLEV, I_KBDIA, IKICE
INTEGER, PARAMETER :: I_KCH1 = 0

INTEGER(KIND=JPIM) :: I_KCOUNT(KLON)
INTEGER(KIND=JPIM) :: I_KCLTOP(KLON)
INTEGER(KIND=JPIM) :: I_KCLBAS(KLON)

REAL(KIND=JPRB) :: ZCAPE(KLON)

REAL(KIND=JPRB) :: ZW     (KLON,KLEV)
REAL(KIND=JPRB) :: ZDTDT  (KLON,KLEV)
REAL(KIND=JPRB) :: ZDQVDT (KLON,KLEV)
REAL(KIND=JPRB) :: ZDQLDT (KLON,KLEV)
REAL(KIND=JPRB) :: ZDQIDT (KLON,KLEV)
REAL(KIND=JPRB) :: ZDUDT  (KLON,KLEV)
REAL(KIND=JPRB) :: ZDVDT  (KLON,KLEV)
REAL(KIND=JPRB) :: ZUMF   (KLON,KLEV)
REAL(KIND=JPRB) :: ZUQV   (KLON,KLEV)
REAL(KIND=JPRB) :: ZUQL   (KLON,KLEV)
REAL(KIND=JPRB) :: ZDPSG  (KLON,KLEV)
REAL(KIND=JPRB) :: ZLV    (KLON,KLEV)
REAL(KIND=JPRB) :: ZLS    (KLON,KLEV)
REAL(KIND=JPRB) :: ZQC    (KLON,KLEV)
REAL(KIND=JPRB) :: ZBETA  (KLON,KLEV)
REAL(KIND=JPRB) :: ZAPHIF (KLON,KLEV)

REAL(KIND=JPRB) :: ZTHETA (KLON,0:KLEV+1)
REAL(KIND=JPRB) :: ZRHO   (KLON,0:KLEV+1)

REAL(KIND=JPRB) :: ZCH1   (KLON,KLEV,0)
REAL(KIND=JPRB) :: ZCH1TEN(KLON,KLEV,0)
REAL(KIND=JPRB) :: ZTKECLS(KLON)
REAL(KIND=JPRB) :: ZUMFMAX(KLON)

INTEGER  :: JI, JK, JKP, JN  ! loop index
! Local arrays (upside/down) necessary for change of ECMWF arrays to convection arrays
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZT     ! grid scale T at time t  (K)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZRV    ! grid scale water vapor  (kg/kg)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZRC    ! grid scale r_c mixing ratio (kg/kg)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZRI    ! grid scale r_i mixing ratio (kg/kg)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZU     ! grid scale horiz. wind u (m/s)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZV     ! grid scale horiz. wind v (m/s)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZW     ! grid scale vertical velocity (m/s)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZPABS  ! grid scale pressure (Pa)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZZZ    ! height of model layer (m)

REAL , DIMENSION(KLON,KLEV) :: SHAL_ZTTEN  ! convective temperat. tendency (K/s)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZRVTEN ! convective r_v tendency (1/s)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZRCTEN ! convective r_c tendency (1/s)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZRITEN ! convective r_i tendency (1/s)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZUTEN  ! convective u tendency (m/s^2)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZVTEN  ! convective m tendency (m/s^2)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZUMF   ! updraft mass flux   (kg/s m2)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZURV   ! water vapor in updrafts (kg/kg)
REAL , DIMENSION(KLON,KLEV) :: SHAL_ZURCI  ! total condensate in updrafts (kg/kg)
INTEGER,  DIMENSION(KLON)   :: SHAL_ICLTOP ! cloud top level (number of model level)
INTEGER,  DIMENSION(KLON)   :: SHAL_ICLBAS ! cloud base level(number of model level)
REAL , DIMENSION(KLON,KLEV,I_KCH1):: SHAL_ZCH1     ! grid scale chemical species
REAL , DIMENSION(KLON,KLEV,I_KCH1):: SHAL_ZCH1TEN  ! chemical convective tendency

! special for shallow convection
REAL , DIMENSION(KLON,KLEV,I_KCH1) :: SHAL_ZCH1TENS
INTEGER,  DIMENSION(KLON)   :: SHAL_ICLBASS, SHAL_ICLTOPS
!-----------------------------------------------------------------

#include "fcttrm.func.h"
#include "wrscmr.intfb.h"

!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACVPPKF',0,ZHOOK_HANDLE)
ASSOCIATE(RKFBNBX=>YDML_PHY_MF%YRPHY0%RKFBNBX, RKFBTAU=>YDML_PHY_MF%YRPHY0%RKFBTAU, RQLCR=>YDML_PHY_MF%YRPHY0%RQLCR, &
 & RPRTH=>YDML_PHY_MF%YRPHY0%RPRTH, ECTMIN=>YDML_PHY_MF%YRPHY0%ECTMIN, AECLS4=>YDML_PHY_MF%YRPHY0%AECLS4, &
 & TSPHY=>YDML_PHY_MF%YRPHY2%TSPHY, &
 & LSMOOTH=>YDML_PHY_MF%YRCVMNH%LSMOOTH, &
 & OTADJS=>YDML_PHY_MF%YRCVMNH%OTADJS, &
 & LSETTADJ=>YDML_PHY_MF%YRCVMNH%LSETTADJ, &
 & RATM=>YDCST%RATM, RCPD=>YDCST%RCPD, RCPV=>YDCST%RCPV, RCS=>YDCST%RCS, RCW=>YDCST%RCW, &
 & RG=>YDCST%RG, RKAPPA=>YDCST%RKAPPA, RLSZER=>YDCST%RLSZER, RLVZER=>YDCST%RLVZER, &
 & LECT=>YDML_PHY_MF%YRPHY%LECT, LCVDD=>YDML_PHY_MF%YRPHY%LCVDD)
!-----------------------------------------------------------------

ZVMD=RCPV-RCPD
ZWMD=RCW-RCPD
ZSMD=RCS-RCPD

ZUMFMAX(:)=1.E-12_JPRB
ZW(:,:)=0.0_JPRB
ZTKECLS(:)=0.0_JPRB
IF (LECT) THEN
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      ZW(JLON,JLEV) = SQRT(MIN(3.0_JPRB,MAX(ECTMIN,PTKE(JLON,JLEV)))/AECLS4)
    ENDDO
  ENDDO
  DO JLON=KIDIA,KFDIA
     ZTKECLS(JLON)=PTKE(JLON,KLEV) 
  ENDDO
ENDIF
DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    ZDPSG (JLON,JLEV) = PDELP(JLON,JLEV)/RG
    ZAPHIF(JLON,JLEV) = PAPHIF(JLON,JLEV)/RG
  ENDDO
ENDDO

DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    ZDTDT (JLON,JLEV) = 0.0_JPRB
    ZDQVDT(JLON,JLEV) = 0.0_JPRB
    ZDQLDT(JLON,JLEV) = 0.0_JPRB
    ZDQIDT(JLON,JLEV) = 0.0_JPRB
    ZDUDT (JLON,JLEV) = 0.0_JPRB
    ZDVDT (JLON,JLEV) = 0.0_JPRB
    ZUMF  (JLON,JLEV) = 0.0_JPRB
    ZUQV  (JLON,JLEV) = 0.0_JPRB
    ZUQL  (JLON,JLEV) = 0.0_JPRB
    ZLV   (JLON,JLEV) = FOLH(PT(JLON,JLEV),0.0_JPRB)
    ZLS   (JLON,JLEV) = FOLH(PT(JLON,JLEV),1.0_JPRB)
    ZQC   (JLON,JLEV) = PQL(JLON,JLEV)+PQI(JLON,JLEV)
  ENDDO
ENDDO

DO JLON=1,KLON
  I_KCOUNT(JLON) = 0
  ZCAPE   (JLON) = 0.0_JPRB
ENDDO

I_KBDIA=1
IKICE=1
!I_KCH1=0
LLREFRESH_ALL=.TRUE.
LLDOWN=LCVDD
LLUVTRANS=.FALSE. ! not yet well tested but possible to use
LLOCHTRANS=.FALSE.
ZDTCONV=TSPHY

I_KCLTOP(:)  = 1 ! set default value when no convection
I_KCLBAS(:)  = 1 ! can be changed  depending on user
SHAL_ICLTOP(:)  = 1
SHAL_ICLBAS(:)  = 1
SHAL_ICLTOPS(:) = 1
SHAL_ICLBASS(:) = 1

!*       2.   Flip arrays upside-down as  first vertical level in convection is 1
!             --------------------------------------------------------------------

DO JK = 1, KLEV
  JKP = KLEV - JK + 1
  DO JI = KIDIA, KFDIA
    SHAL_ZPABS(JI,JKP) = PAPRSF(JI,JK)
    SHAL_ZZZ(JI,JKP)   = ZAPHIF(JI,JK)
    SHAL_ZT(JI,JKP)    = PT(JI,JK)
    SHAL_ZRV(JI,JKP)   = PQ(JI,JK) / ( 1.0 - PQ(JI,JK) ) ! transform specific humidity
    SHAL_ZRC(JI,JKP)   = PQL(JI,JK) / ( 1.0 - PQL(JI,JK) ) ! in mixing ratio
    SHAL_ZRI(JI,JKP)   = PQI(JI,JK) / ( 1.0 - PQI(JI,JK) )
    SHAL_ZU(JI,JKP)    = PU(JI,JK)
    SHAL_ZV(JI,JKP)    = PV(JI,JK)
    SHAL_ZW(JI,JKP)    = ZW(JI,JK)
  ENDDO
ENDDO
IF ( LLOCHTRANS ) THEN
  DO JK = 1, KLEV
    JKP = KLEV - JK + 1
    DO JN = 1, I_KCH1
      DO JI = KIDIA, KFDIA
        SHAL_ZCH1(JI,JKP,JN) = ZCH1(JI,JK,JN)
      ENDDO
    ENDDO
  ENDDO
ENDIF

  I_KCOUNT(:)     =0
  SHAL_ZTTEN(:,:)    =0.0
  SHAL_ZRVTEN(:,:)   =0.0
  SHAL_ZRCTEN(:,:)   =0.0
  SHAL_ZRITEN(:,:)   =0.0
  SHAL_ZUTEN(:,:)    =0.0
  SHAL_ZVTEN(:,:)    =0.0
  SHAL_ZUMF(:,:)     =0.0
  SHAL_ZURV(:,:)     =0.0
  SHAL_ZURCI(:,:)    =0.0
  SHAL_ZCH1TEN(:,:,:)=0.0
  ZCAPE(:)      =0.0

!*       4.b  Call shallow convection routine
!             -------------------------------

  CALL INI_CONVPAR
  XA25=YDML_PHY_MF%YRCVMNH%XA25
  XCRAD=YDML_PHY_MF%YRCVMNH%XCRAD
  XCDEPTH=YDML_PHY_MF%YRCVMNH%XCDEPTH
  XCDEPTH_D=YDML_PHY_MF%YRCVMNH%XCDEPTH_D
  XDTPERT=YDML_PHY_MF%YRCVMNH%XDTPERT
  XATPERT=YDML_PHY_MF%YRCVMNH%XATPERT
  XBTPERT=YDML_PHY_MF%YRCVMNH%XBTPERT
  XENTR=YDML_PHY_MF%YRCVMNH%XENTR
  XZLCL=YDML_PHY_MF%YRCVMNH%XZLCL
  XZPBL=YDML_PHY_MF%YRCVMNH%XZPBL
  XWTRIG=YDML_PHY_MF%YRCVMNH%XWTRIG
  XNHGAM=YDML_PHY_MF%YRCVMNH%XNHGAM
  XTFRZ1=YDML_PHY_MF%YRCVMNH%XTFRZ1
  XTFRZ2=YDML_PHY_MF%YRCVMNH%XTFRZ2
  XSTABT=YDML_PHY_MF%YRCVMNH%XSTABT
  XSTABC=YDML_PHY_MF%YRCVMNH%XSTABC
  XAW=YDML_PHY_MF%YRCVMNH%XAW
  XBW=YDML_PHY_MF%YRCVMNH%XBW
  LLSMOOTH=LSMOOTH

  CALL SHALLOW_CONVECTION( KLON, KLEV, KIDIA, KFDIA, I_KBDIA, KTDIA,        &
   & IKICE, LSETTADJ, OTADJS,         &
   & SHAL_ZPABS, SHAL_ZZZ,ZTKECLS,                         &
   & SHAL_ZT, SHAL_ZRV, SHAL_ZRC, SHAL_ZRI, SHAL_ZW,                      &
   & SHAL_ZTTEN, SHAL_ZRVTEN, SHAL_ZRCTEN, SHAL_ZRITEN,          &
   & SHAL_ICLTOPS, SHAL_ICLBASS, SHAL_ZUMF, &
   & LLOCHTRANS, I_KCH1, SHAL_ZCH1, SHAL_ZCH1TENS )

DO JI = KIDIA, KFDIA
  SHAL_ICLTOP(JI)   = MAX(SHAL_ICLTOP(JI), SHAL_ICLTOPS(JI))
  SHAL_ICLBAS(JI)   = MAX(SHAL_ICLBAS(JI), SHAL_ICLBASS(JI))
ENDDO

!*       6.  Reflip arrays to ECMWF/ARPEGE vertical structure
!            change mixing ratios to sepcific humidity

DO JK = 1, KLEV
  JKP = KLEV - JK + 1
  DO JI = KIDIA, KFDIA
    ZDTDT(JI,JK)  = SHAL_ZTTEN(JI,JKP)
   ! don't transform back to specific hum, does not conserve integrals
    ZDQVDT(JI,JK) = SHAL_ZRVTEN(JI,JKP) ! / ( 1.0 + ZRV(JI,JKP) ) ** 2
    ZDQLDT(JI,JK) = SHAL_ZRCTEN(JI,JKP) ! / ( 1.0 + ZRC(JI,JKP) ) ** 2
    ZDQIDT(JI,JK) = SHAL_ZRITEN(JI,JKP) ! / ( 1.0 + ZRI(JI,JKP) ) ** 2
    ZDUDT(JI,JK)  = SHAL_ZUTEN(JI,JKP)
    ZDVDT(JI,JK)  = SHAL_ZVTEN(JI,JKP)
    ZUMF(JI,JK)   = SHAL_ZUMF(JI,JKP)
    ZUQV(JI,JK)   = SHAL_ZURV(JI,JKP) / ( 1.0 + SHAL_ZURV(JI,JKP) )
    ZUQL(JI,JK)  = SHAL_ZURCI(JI,JKP)/ ( 1.0 + SHAL_ZURCI(JI,JKP) )
  ENDDO
ENDDO

DO JI = KIDIA, KFDIA
  JK = SHAL_ICLTOP(JI)
  I_KCLTOP(JI) = KLEV - JK + 1
  JK = SHAL_ICLBAS(JI)
  I_KCLBAS(JI) = KLEV - JK + 1
  IF ( SHAL_ICLTOP(JI) == 1 ) I_KCLTOP(JI) = 1
  IF ( SHAL_ICLBAS(JI) == 1 ) I_KCLBAS(JI) = 1
ENDDO

IF ( LLOCHTRANS ) THEN
  DO JK = 1, KLEV
    JKP = KLEV - JK + 1
    DO JN = 1, I_KCH1
      DO JI = KIDIA, KFDIA
        ZCH1TEN(JI,JK,JN) = SHAL_ZCH1TEN(JI,JKP,JN)
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZMF_shal',ZUMF,KLON,KLEV)

! Calcul de la production thermique pour la TKE
PPRODTH(:,:)=0.0_JPRB
IF (RPRTH > 0._JPRB) THEN
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      ZBETA (JLON,JLEV) = (RATM/PAPRSF(JLON,JLEV))**(RKAPPA)
      ZTHETA(JLON,JLEV) = PT(JLON,JLEV)*ZBETA(JLON,JLEV)
      ZRHO  (JLON,JLEV) = PAPRSF(JLON,JLEV)/(PR(JLON,JLEV)*PT(JLON,JLEV))
    ENDDO
  ENDDO
  DO JLON=KIDIA,KFDIA
    ZTHETA(JLON,0)      = ZTHETA(JLON,1)
    ZRHO  (JLON,0)      = ZRHO  (JLON,1)
    ZTHETA(JLON,KLEV+1) = ZTHETA(JLON,KLEV)
    ZRHO  (JLON,KLEV+1) = ZRHO  (JLON,KLEV)
  ENDDO

  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      ZTDCP=ZVMD*ZDQVDT(JLON,JLEV)+ZWMD*ZDQLDT(JLON,JLEV)+ZSMD*ZDQIDT(JLON,JLEV)
      ZDQCDT=ZDQLDT(JLON,JLEV)+ZDQIDT(JLON,JLEV)
      ZDTLDT=ZBETA (JLON,JLEV)&
        &  * ( ZDTDT(JLON,JLEV) + ZLV(JLON,JLEV)/PCP(JLON,JLEV)&
        &    * ( ZQC(JLON,JLEV)*ZTDCP/PCP(JLON,JLEV)-ZDQCDT ) )
      PPRODTH(JLON,JLEV)=PPRODTH(JLON,JLEV-1)-ZDPSG(JLON,JLEV)*ZDTLDT
    ENDDO
  ENDDO

  DO JLEV=0,KLEV
    DO JLON=KIDIA,KFDIA
      PPRODTH(JLON,JLEV)=PPRODTH(JLON,JLEV)*RG*4._JPRB&
        & / ( ZRHO  (JLON,JLEV) + ZRHO  (JLON,JLEV+1) )&
        & / ( ZTHETA(JLON,JLEV) + ZTHETA(JLON,JLEV+1) )
    ENDDO
  ENDDO
ENDIF ! Fin du calcul de la production thermique pour la TKE

DO JLON=KIDIA,KFDIA
  DO JLEV=1,KLEV
    KNLAB(JLON,JLEV)=1-MAX(0,MIN(1,(I_KCLTOP(JLON)-JLEV)*(I_KCLBAS(JLON)-JLEV)))
  ENDDO
  KNND(JLON)=MIN(1,I_KCLTOP(JLON)-1)
ENDDO

! Calcul de la nebulosite et de l'eau condensee
ZEPS=1.E-12_JPRB
IF (RKFBTAU > 0._JPRB) THEN
  DO JLEV=1,KLEV
!DEC$ IVDEP
    DO JLON=KIDIA,KFDIA
      ZDQCDT=MAX(0.0_JPRB,ZDQLDT(JLON,JLEV)+ZDQIDT(JLON,JLEV))
      PQCPP (JLON,JLEV)=RKFBTAU*ZDQCDT*FLOAT(KNLAB(JLON,JLEV))
      PNEBPP(JLON,JLEV)=MAX(ZEPS,MIN(RKFBNBX,PQCPP(JLON,JLEV)/RQLCR))
      ZUMFMAX(JLON)=MAX(ZUMFMAX(JLON),ZUMF(JLON,JLEV))
    ENDDO
  ENDDO
  IF (.NOT.LSMOOTH) THEN
    DO JLEV=1,KLEV
!DEC$ IVDEP
      DO JLON=KIDIA,KFDIA
        PQCPP (JLON,JLEV)=PQCPP(JLON,JLEV)&
         & *MIN(1.0_JPRB,ZUMF(JLON,JLEV)/ZUMFMAX(JLON))
        PNEBPP(JLON,JLEV)=MAX(ZEPS,MIN(RKFBNBX,PQCPP(JLON,JLEV)/RQLCR))
      ENDDO
    ENDDO
  ENDIF
ENDIF ! Fin du calcul de la nebulosite et de l'eau condensee

LLCONDWT=.FALSE.
IF (.NOT.LLCONDWT) THEN
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      ZDQVDT(JLON,JLEV)=ZDQVDT(JLON,JLEV)+ZDQLDT(JLON,JLEV)+ZDQIDT(JLON,JLEV)
      ZDTDT(JLON,JLEV)=ZDTDT(JLON,JLEV) - (ZLV(JLON,JLEV)*ZDQLDT(JLON,JLEV)&
        & + ZLS(JLON,JLEV)*ZDQIDT(JLON,JLEV)) / PCP(JLON,JLEV)
      ZDQLDT(JLON,JLEV)=0.0_JPRB
      ZDQIDT(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      PFCCQL(JLON,JLEV)=PFCCQL(JLON,JLEV-1)+ZDPSG(JLON,JLEV)*ZDQLDT(JLON,JLEV)
      PFCCQN(JLON,JLEV)=PFCCQN(JLON,JLEV-1)+ZDPSG(JLON,JLEV)*ZDQIDT(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF

DO JLEV=1,KLEV
!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
    PDIFCQ(JLON,JLEV)=PDIFCQ(JLON,JLEV-1)-ZDPSG(JLON,JLEV)*ZDQVDT(JLON,JLEV)
    ZTDCP=ZVMD*ZDQVDT(JLON,JLEV)+ZWMD*ZDQLDT(JLON,JLEV)+ZSMD*ZDQIDT(JLON,JLEV)
    PDIFCS(JLON,JLEV)=PDIFCS(JLON,JLEV-1)-ZDPSG(JLON,JLEV)&
     & * (ZDTDT(JLON,JLEV)*(PCP(JLON,JLEV)+TSPHY*ZTDCP)+PT(JLON,JLEV)*ZTDCP)
  ENDDO
ENDDO

DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    PDIFCQ(JLON,JLEV)=PDIFCQ(JLON,JLEV)-PFCCQL(JLON,JLEV)-PFCCQN(JLON,JLEV)
    PDIFCS(JLON,JLEV)=PDIFCS(JLON,JLEV)&
     & + RLVZER*PFCCQL(JLON,JLEV) + RLSZER*PFCCQN(JLON,JLEV)
  ENDDO
ENDDO

!-----------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACVPPKF',1,ZHOOK_HANDLE)
END SUBROUTINE ACVPPKF

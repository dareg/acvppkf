!OPTIONS XOPT(NOEVAL)
!-----------------------------------------------------------------
SUBROUTINE ACVPPKF(YDCST,YDML_PHY_MF,CST,D,NSV,CONVPAR,KTDIA, &
 !-----------------------------------------------------------------
 ! - INPUT  2D .
 & PAPRSF, PAPHIF, PDELP, PR, PT, PQ, &
 & PQL, PQI, PU, PV, PCP, PTKE, &
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
USE MODD_CONVPAR, ONLY : CONVPAR_T
USE MODD_CONVPAR_SHAL, ONLY : CONVPAR_SHAL
USE MODD_CST, ONLY: CST_T
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_NSV, ONLY: NSV_T
!-----------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST)         ,INTENT (IN)   :: YDCST
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
TYPE(CST_T)        ,INTENT(IN)    :: CST
TYPE(DIMPHYEX_T)   ,INTENT(IN)    :: D
TYPE(NSV_T)        ,INTENT(IN)    :: NSV
TYPE(CONVPAR_T)    ,INTENT(IN)    :: CONVPAR
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KTDIA 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PAPRSF (D%NIT,D%NKT) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PAPHIF (D%NIT,D%NKT) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PDELP  (D%NIT,D%NKT) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PR     (D%NIT,D%NKT) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PT     (D%NIT,D%NKT) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQ     (D%NIT,D%NKT) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQL    (D%NIT,D%NKT) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQI    (D%NIT,D%NKT) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PU     (D%NIT,D%NKT) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PV     (D%NIT,D%NKT) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PCP    (D%NIT,D%NKT)
REAL(KIND=JPRB)    ,INTENT(IN)    :: PTKE   (D%NIT,D%NKT)
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PDIFCQ (D%NIT,0:D%NKT) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PDIFCS (D%NIT,0:D%NKT) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PFCCQL (D%NIT,0:D%NKT) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PFCCQN (D%NIT,0:D%NKT) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PPRODTH(D%NIT,0:D%NKT)
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PQCPP  (D%NIT,D%NKT)
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PNEBPP (D%NIT,D%NKT)
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KNLAB  (D%NIT,D%NKT) 
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KNND   (D%NIT) 

!-----------------------------------------------------------------

LOGICAL :: LLREFRESH_ALL, LLDOWN, LLUVTRANS, LLOCHTRANS, LLCONDWT

REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZDTCONV, ZVMD, ZWMD, ZSMD, ZTDCP, ZEPS, ZDQCDT, ZDTLDT

INTEGER(KIND=JPIM) :: JLON, JLEV, I_KBDIA, IKICE
INTEGER, PARAMETER :: I_KCH1 = 0

INTEGER(KIND=JPIM) :: I_KCOUNT(D%NIT)
INTEGER(KIND=JPIM) :: I_KCLTOP(D%NIT)
INTEGER(KIND=JPIM) :: I_KCLBAS(D%NIT)

REAL(KIND=JPRB) :: ZCAPE(D%NIT)

REAL(KIND=JPRB) :: ZW     (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZDTDT  (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZDQVDT (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZDQLDT (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZDQIDT (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZDUDT  (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZDVDT  (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZUMF   (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZUQV   (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZUQL   (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZDPSG  (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZLV    (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZLS    (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZQC    (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZBETA  (D%NIT,D%NKT)
REAL(KIND=JPRB) :: ZAPHIF (D%NIT,D%NKT)

REAL(KIND=JPRB) :: ZTHETA (D%NIT,0:D%NKT+1)
REAL(KIND=JPRB) :: ZRHO   (D%NIT,0:D%NKT+1)

REAL(KIND=JPRB) :: ZCH1   (D%NIT,D%NKT,I_KCH1)
REAL(KIND=JPRB) :: ZCH1TEN(D%NIT,D%NKT,I_KCH1)
REAL(KIND=JPRB) :: ZTKECLS(D%NIT)
REAL(KIND=JPRB) :: ZUMFMAX(D%NIT)

INTEGER  :: JI, JK, JKP, JN  ! loop index
! Local arrays (upside/down) necessary for change of ECMWF arrays to convection arrays
REAL , DIMENSION(D%NIT,D%NKT) :: ZT     ! grid scale T at time t  (K)
REAL , DIMENSION(D%NIT,D%NKT) :: ZRV    ! grid scale water vapor  (kg/kg)
REAL , DIMENSION(D%NIT,D%NKT) :: ZRC    ! grid scale r_c mixing ratio (kg/kg)
REAL , DIMENSION(D%NIT,D%NKT) :: ZRI    ! grid scale r_i mixing ratio (kg/kg)
REAL , DIMENSION(D%NIT,D%NKT) :: ZU     ! grid scale horiz. wind u (m/s)
REAL , DIMENSION(D%NIT,D%NKT) :: ZV     ! grid scale horiz. wind v (m/s)
REAL , DIMENSION(D%NIT,D%NKT) :: ZZW    ! grid scale vertical velocity (m/s)
REAL , DIMENSION(D%NIT,D%NKT) :: ZPABS  ! grid scale pressure (Pa)
REAL , DIMENSION(D%NIT,D%NKT) :: ZZZ    ! height of model layer (m)

REAL , DIMENSION(D%NIT,D%NKT) :: ZTTEN  ! convective temperat. tendency (K/s)
REAL , DIMENSION(D%NIT,D%NKT) :: ZRVTEN ! convective r_v tendency (1/s)
REAL , DIMENSION(D%NIT,D%NKT) :: ZRCTEN ! convective r_c tendency (1/s)
REAL , DIMENSION(D%NIT,D%NKT) :: ZRITEN ! convective r_i tendency (1/s)
REAL , DIMENSION(D%NIT,D%NKT) :: ZUTEN  ! convective u tendency (m/s^2)
REAL , DIMENSION(D%NIT,D%NKT) :: ZVTEN  ! convective m tendency (m/s^2)
REAL , DIMENSION(D%NIT,D%NKT) :: ZZUMF  ! updraft mass flux   (kg/s m2)
REAL , DIMENSION(D%NIT,D%NKT) :: ZURV   ! water vapor in updrafts (kg/kg)
REAL , DIMENSION(D%NIT,D%NKT) :: ZURCI  ! total condensate in updrafts (kg/kg)
INTEGER,  DIMENSION(D%NIT)   :: ICLTOP ! cloud top level (number of model level)
INTEGER,  DIMENSION(D%NIT)   :: ICLBAS ! cloud base level(number of model level)
REAL , DIMENSION(D%NIT,D%NKT,I_KCH1):: SHAL_ZCH1     ! grid scale chemical species
REAL , DIMENSION(D%NIT,D%NKT,I_KCH1):: SHAL_ZCH1TEN  ! chemical convective tendency
! special for shallow convection
REAL , DIMENSION(D%NIT,D%NKT,I_KCH1) :: SHAL_ZCH1TENS
INTEGER,  DIMENSION(D%NIT)   :: ICLBASS, ICLTOPS
!
TYPE(CONVPAR_SHAL) :: CVP_SHAL
!-----------------------------------------------------------------

#include "fcttrm.func.h"
#include "wrscmr.intfb.h"
#include "shallow_convection.h"

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
  DO JLEV=1,D%NKT
    DO JLON=D%NIB,D%NIE
      ZW(JLON,JLEV) = SQRT(MIN(3.0_JPRB,MAX(ECTMIN,PTKE(JLON,JLEV)))/AECLS4)
    ENDDO
  ENDDO
  DO JLON=D%NIB,D%NIE
     ZTKECLS(JLON)=PTKE(JLON,D%NKT) 
  ENDDO
ENDIF
DO JLEV=1,D%NKT
  DO JLON=D%NIB,D%NIE
    ZDPSG (JLON,JLEV) = PDELP(JLON,JLEV)/RG
    ZAPHIF(JLON,JLEV) = PAPHIF(JLON,JLEV)/RG
  ENDDO
ENDDO

DO JLEV=1,D%NKT
  DO JLON=D%NIB,D%NIE
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

DO JLON=1,D%NIE
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
ICLTOP(:)  = 1
ICLBAS(:)  = 1
ICLTOPS(:) = 1
ICLBASS(:) = 1

!*       2.   Flip arrays upside-down as  first vertical level in convection is 1
!             --------------------------------------------------------------------

DO JK = 1, D%NKT
  JKP = D%NKT - JK + 1
  DO JI = D%NIB, D%NIE
    ZPABS(JI,JKP) = PAPRSF(JI,JK)
    ZZZ(JI,JKP)   = ZAPHIF(JI,JK)
    ZT(JI,JKP)    = PT(JI,JK)
    ZRV(JI,JKP)   = PQ(JI,JK) / ( 1.0 - PQ(JI,JK) ) ! transform specific humidity
    ZRC(JI,JKP)   = PQL(JI,JK) / ( 1.0 - PQL(JI,JK) ) ! in mixing ratio
    ZRI(JI,JKP)   = PQI(JI,JK) / ( 1.0 - PQI(JI,JK) )
    ZU(JI,JKP)    = PU(JI,JK)
    ZV(JI,JKP)    = PV(JI,JK)
    ZZW(JI,JKP)   = ZW(JI,JK)
  ENDDO
ENDDO
IF ( LLOCHTRANS ) THEN
  DO JK = 1, D%NKT
    JKP = D%NKT - JK + 1
    DO JN = 1, I_KCH1
      DO JI = D%NIB, D%NIE
        SHAL_ZCH1(JI,JKP,JN) = ZCH1(JI,JK,JN)
      ENDDO
    ENDDO
  ENDDO
ENDIF

  I_KCOUNT(:)   =0
  ZTTEN(:,:)    =0.0
  ZRVTEN(:,:)   =0.0
  ZRCTEN(:,:)   =0.0
  ZRITEN(:,:)   =0.0
  ZUTEN(:,:)    =0.0
  ZVTEN(:,:)    =0.0
  ZZUMF(:,:)    =0.0
  ZURV(:,:)     =0.0
  ZURCI(:,:)    =0.0
  SHAL_ZCH1TEN(:,:,:)=0.0
  ZCAPE(:)      =0.0

!*       4.b  Call shallow convection routine
!             -------------------------------

  CVP_SHAL%XA25=YDML_PHY_MF%YRCVMNH%XA25
  CVP_SHAL%XCRAD=YDML_PHY_MF%YRCVMNH%XCRAD
  CVP_SHAL%XCDEPTH=YDML_PHY_MF%YRCVMNH%XCDEPTH
  CVP_SHAL%XCDEPTH_D=YDML_PHY_MF%YRCVMNH%XCDEPTH_D
  CVP_SHAL%XDTPERT=YDML_PHY_MF%YRCVMNH%XDTPERT
  CVP_SHAL%XATPERT=YDML_PHY_MF%YRCVMNH%XATPERT
  CVP_SHAL%XBTPERT=YDML_PHY_MF%YRCVMNH%XBTPERT
  CVP_SHAL%XENTR=YDML_PHY_MF%YRCVMNH%XENTR
  CVP_SHAL%XZLCL=YDML_PHY_MF%YRCVMNH%XZLCL
  CVP_SHAL%XZPBL=YDML_PHY_MF%YRCVMNH%XZPBL
  CVP_SHAL%XWTRIG=YDML_PHY_MF%YRCVMNH%XWTRIG
  CVP_SHAL%XNHGAM=YDML_PHY_MF%YRCVMNH%XNHGAM
  CVP_SHAL%XTFRZ1=YDML_PHY_MF%YRCVMNH%XTFRZ1
  CVP_SHAL%XTFRZ2=YDML_PHY_MF%YRCVMNH%XTFRZ2
  CVP_SHAL%XSTABT=YDML_PHY_MF%YRCVMNH%XSTABT
  CVP_SHAL%XSTABC=YDML_PHY_MF%YRCVMNH%XSTABC
  CVP_SHAL%XAW=YDML_PHY_MF%YRCVMNH%XAW
  CVP_SHAL%XBW=YDML_PHY_MF%YRCVMNH%XBW
  CVP_SHAL%LLSMOOTH=LSMOOTH

  CALL SHALLOW_CONVECTION(CVP_SHAL, CST, D, NSV, CONVPAR, I_KBDIA,     &
                          KTDIA, IKICE, LSETTADJ, OTADJS, ZPABS, ZZZ,  &
                          ZTKECLS, ZT, ZRV, ZRC, ZRI, ZZW, ZTTEN,      &
                          ZRVTEN, ZRCTEN, ZRITEN, ICLTOPS, ICLBASS,    &
                          ZZUMF, LLOCHTRANS, I_KCH1, SHAL_ZCH1,        &
                          SHAL_ZCH1TENS)

DO JI = D%NIB, D%NIE
  ICLTOP(JI)   = MAX(ICLTOP(JI), ICLTOPS(JI))
  ICLBAS(JI)   = MAX(ICLBAS(JI), ICLBASS(JI))
ENDDO

!*       6.  Reflip arrays to ECMWF/ARPEGE vertical structure
!            change mixing ratios to sepcific humidity

DO JK = 1, D%NKT
  JKP = D%NKT - JK + 1
  DO JI = D%NIB, D%NIE
    ZDTDT(JI,JK)  = ZTTEN(JI,JKP)
   ! don't transform back to specific hum, does not conserve integrals
    ZDQVDT(JI,JK) = ZRVTEN(JI,JKP) ! / ( 1.0 + ZRV(JI,JKP) ) ** 2
    ZDQLDT(JI,JK) = ZRCTEN(JI,JKP) ! / ( 1.0 + ZRC(JI,JKP) ) ** 2
    ZDQIDT(JI,JK) = ZRITEN(JI,JKP) ! / ( 1.0 + ZRI(JI,JKP) ) ** 2
    ZDUDT(JI,JK)  = ZUTEN(JI,JKP)
    ZDVDT(JI,JK)  = ZVTEN(JI,JKP)
    ZUMF(JI,JK)   = ZZUMF(JI,JKP)
    ZUQV(JI,JK)   = ZURV(JI,JKP) / ( 1.0 + ZURV(JI,JKP) )
    ZUQL(JI,JK)  = ZURCI(JI,JKP)/ ( 1.0 + ZURCI(JI,JKP) )
  ENDDO
ENDDO

DO JI = D%NIB, D%NIE
  JK = ICLTOP(JI)
  I_KCLTOP(JI) = D%NKT - JK + 1
  JK = ICLBAS(JI)
  I_KCLBAS(JI) = D%NKT - JK + 1
  IF ( ICLTOP(JI) == 1 ) I_KCLTOP(JI) = 1
  IF ( ICLBAS(JI) == 1 ) I_KCLBAS(JI) = 1
ENDDO

IF ( LLOCHTRANS ) THEN
  DO JK = 1, D%NKT
    JKP = D%NKT - JK + 1
    DO JN = 1, I_KCH1
      DO JI = D%NIB, D%NIE
        ZCH1TEN(JI,JK,JN) = SHAL_ZCH1TEN(JI,JKP,JN)
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZMF_shal',ZUMF,D%NIT,D%NKT)

! Calcul de la production thermique pour la TKE
PPRODTH(:,:)=0.0_JPRB
IF (RPRTH > 0._JPRB) THEN
  DO JLEV=1,D%NKT
    DO JLON=D%NIB,D%NIE
      ZBETA (JLON,JLEV) = (RATM/PAPRSF(JLON,JLEV))**(RKAPPA)
      ZTHETA(JLON,JLEV) = PT(JLON,JLEV)*ZBETA(JLON,JLEV)
      ZRHO  (JLON,JLEV) = PAPRSF(JLON,JLEV)/(PR(JLON,JLEV)*PT(JLON,JLEV))
    ENDDO
  ENDDO
  DO JLON=D%NIB,D%NIE
    ZTHETA(JLON,0)      = ZTHETA(JLON,1)
    ZRHO  (JLON,0)      = ZRHO  (JLON,1)
    ZTHETA(JLON,D%NKT+1) = ZTHETA(JLON,D%NKT)
    ZRHO  (JLON,D%NKT+1) = ZRHO  (JLON,D%NKT)
  ENDDO

  DO JLEV=1,D%NKT
    DO JLON=D%NIB,D%NIE
      ZTDCP=ZVMD*ZDQVDT(JLON,JLEV)+ZWMD*ZDQLDT(JLON,JLEV)+ZSMD*ZDQIDT(JLON,JLEV)
      ZDQCDT=ZDQLDT(JLON,JLEV)+ZDQIDT(JLON,JLEV)
      ZDTLDT=ZBETA (JLON,JLEV)&
        &  * ( ZDTDT(JLON,JLEV) + ZLV(JLON,JLEV)/PCP(JLON,JLEV)&
        &    * ( ZQC(JLON,JLEV)*ZTDCP/PCP(JLON,JLEV)-ZDQCDT ) )
      PPRODTH(JLON,JLEV)=PPRODTH(JLON,JLEV-1)-ZDPSG(JLON,JLEV)*ZDTLDT
    ENDDO
  ENDDO

  DO JLEV=0,D%NKT
    DO JLON=D%NIB,D%NIE
      PPRODTH(JLON,JLEV)=PPRODTH(JLON,JLEV)*RG*4._JPRB&
        & / ( ZRHO  (JLON,JLEV) + ZRHO  (JLON,JLEV+1) )&
        & / ( ZTHETA(JLON,JLEV) + ZTHETA(JLON,JLEV+1) )
    ENDDO
  ENDDO
ENDIF ! Fin du calcul de la production thermique pour la TKE

DO JLON=D%NIB,D%NIE
  DO JLEV=1,D%NKT
    KNLAB(JLON,JLEV)=1-MAX(0,MIN(1,(I_KCLTOP(JLON)-JLEV)*(I_KCLBAS(JLON)-JLEV)))
  ENDDO
  KNND(JLON)=MIN(1,I_KCLTOP(JLON)-1)
ENDDO

! Calcul de la nebulosite et de l'eau condensee
ZEPS=1.E-12_JPRB
IF (RKFBTAU > 0._JPRB) THEN
  DO JLEV=1,D%NKT
!DEC$ IVDEP
    DO JLON=D%NIB,D%NIE
      ZDQCDT=MAX(0.0_JPRB,ZDQLDT(JLON,JLEV)+ZDQIDT(JLON,JLEV))
      PQCPP (JLON,JLEV)=RKFBTAU*ZDQCDT*FLOAT(KNLAB(JLON,JLEV))
      PNEBPP(JLON,JLEV)=MAX(ZEPS,MIN(RKFBNBX,PQCPP(JLON,JLEV)/RQLCR))
      ZUMFMAX(JLON)=MAX(ZUMFMAX(JLON),ZUMF(JLON,JLEV))
    ENDDO
  ENDDO
  IF (.NOT.LSMOOTH) THEN
    DO JLEV=1,D%NKT
!DEC$ IVDEP
      DO JLON=D%NIB,D%NIE
        PQCPP (JLON,JLEV)=PQCPP(JLON,JLEV)&
         & *MIN(1.0_JPRB,ZUMF(JLON,JLEV)/ZUMFMAX(JLON))
        PNEBPP(JLON,JLEV)=MAX(ZEPS,MIN(RKFBNBX,PQCPP(JLON,JLEV)/RQLCR))
      ENDDO
    ENDDO
  ENDIF
ENDIF ! Fin du calcul de la nebulosite et de l'eau condensee

LLCONDWT=.FALSE.
IF (.NOT.LLCONDWT) THEN
  DO JLEV=1,D%NKT
    DO JLON=D%NIB,D%NIE
      ZDQVDT(JLON,JLEV)=ZDQVDT(JLON,JLEV)+ZDQLDT(JLON,JLEV)+ZDQIDT(JLON,JLEV)
      ZDTDT(JLON,JLEV)=ZDTDT(JLON,JLEV) - (ZLV(JLON,JLEV)*ZDQLDT(JLON,JLEV)&
        & + ZLS(JLON,JLEV)*ZDQIDT(JLON,JLEV)) / PCP(JLON,JLEV)
      ZDQLDT(JLON,JLEV)=0.0_JPRB
      ZDQIDT(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,D%NKT
    DO JLON=D%NIB,D%NIE
      PFCCQL(JLON,JLEV)=PFCCQL(JLON,JLEV-1)+ZDPSG(JLON,JLEV)*ZDQLDT(JLON,JLEV)
      PFCCQN(JLON,JLEV)=PFCCQN(JLON,JLEV-1)+ZDPSG(JLON,JLEV)*ZDQIDT(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF

DO JLEV=1,D%NKT
!DEC$ IVDEP
  DO JLON=D%NIB,D%NIE
    PDIFCQ(JLON,JLEV)=PDIFCQ(JLON,JLEV-1)-ZDPSG(JLON,JLEV)*ZDQVDT(JLON,JLEV)
    ZTDCP=ZVMD*ZDQVDT(JLON,JLEV)+ZWMD*ZDQLDT(JLON,JLEV)+ZSMD*ZDQIDT(JLON,JLEV)
    PDIFCS(JLON,JLEV)=PDIFCS(JLON,JLEV-1)-ZDPSG(JLON,JLEV)&
     & * (ZDTDT(JLON,JLEV)*(PCP(JLON,JLEV)+TSPHY*ZTDCP)+PT(JLON,JLEV)*ZTDCP)
  ENDDO
ENDDO

DO JLEV=1,D%NKT
  DO JLON=D%NIB,D%NIE
    PDIFCQ(JLON,JLEV)=PDIFCQ(JLON,JLEV)-PFCCQL(JLON,JLEV)-PFCCQN(JLON,JLEV)
    PDIFCS(JLON,JLEV)=PDIFCS(JLON,JLEV)&
     & + RLVZER*PFCCQL(JLON,JLEV) + RLSZER*PFCCQN(JLON,JLEV)
  ENDDO
ENDDO

!-----------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACVPPKF',1,ZHOOK_HANDLE)
END SUBROUTINE ACVPPKF

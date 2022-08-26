INTERFACE
SUBROUTINE CONVECTION_SHAL( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,             &
 & PDTCONV, LD_OREFRESH_ALL, LD_ODOWN, KICE,   &
 & LD_OSETTADJ, PTADJS, LSMOOTH,                  &
 & PA25, PCRAD, PCDEPTH, PCDEPTH_D, PDTPERT, PATPERT, PBTPERT, PENTR, &
 & PZLCL, PZPBL, PWTRIG, PNHGAM, PTFRZ1, PTFRZ2, PSTABT, PSTABC, PAW, PBW, &
 & PPABS, PZZ, PTKECLS,                           &
 & PT, PRV, PRC, PRI, PU, PV, PW,                      &
 & KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,              &
 & PUMF, PCAPE, KCLTOP, KCLBAS, &
 & PURV, PURCI,                                        &
 & LD_OUVTRANS, PUTEN, PVTEN,                             &
 & LD_OCHTRANS, KCH1, PCH1, PCH1TEN )

USE PARKIND1  ,ONLY : JPIM     ,JPRB                                                                   
IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON ! horizontal dimension
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV ! vertical dimension
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA ! value of the first point in x
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA ! value of the last point in x
INTEGER(KIND=JPIM),INTENT(IN)    :: KBDIA ! vertical  computations start at
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA ! vertical computations can be
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTCONV ! Interval of time between two
LOGICAL           ,INTENT(IN)    :: LD_OREFRESH_ALL ! refresh or not all 
LOGICAL           ,INTENT(IN)    :: LD_ODOWN ! take or not convective
INTEGER(KIND=JPIM),INTENT(IN)    :: KICE ! flag for ice ( 1 = yes, 
LOGICAL           ,INTENT(IN)    :: LD_OSETTADJ ! logical to set convective
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTADJS ! user defined shal. adjustment time (s)
LOGICAL           ,INTENT(IN)    :: LSMOOTH ! Supposed to be necessary ...
REAL (KIND=JPRB)  ,INTENT(IN)    :: PA25, PCRAD, PCDEPTH, PCDEPTH_D, PDTPERT, &
& PATPERT, PBTPERT,PENTR, &
& PZLCL, PZPBL, PWTRIG, PNHGAM, PTFRZ1, PTFRZ2, PSTABT, PSTABC, PAW , PBW
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPABS(KLON,KLEV) ! grid scale pressure (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZZ(KLON,KLEV) ! geopotential (m2/s2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTKECLS(KLON) ! TKE in CLS (m2/s2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) ! grid scale T at time t  (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRV(KLON,KLEV) ! grid scale water vapor  (kg/kg)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRC(KLON,KLEV) ! grid scale r_c (kg/kg)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRI(KLON,KLEV) ! grid scale r_i (kg/kg)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KLEV) ! grid scale horiz. wind u (m/s) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON,KLEV) ! grid scale horiz. wind v (m/s)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PW(KLON,KLEV) ! grid scale vertical velocity (m/s)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCOUNT(KLON) ! convective counter(recompute
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTTEN(KLON,KLEV) ! convective temperat. tendency (K/s)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRVTEN(KLON,KLEV) ! convective r_v tendency (1/s)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRCTEN(KLON,KLEV) ! convective r_c tendency (1/s)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRITEN(KLON,KLEV) ! convective r_i tendency (1/s)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUMF(KLON,KLEV) ! updraft mass flux   (kg/s m2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KLON) ! CAPE (J/kg)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCLTOP(KLON) ! cloud top level (number of model level)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCLBAS(KLON) ! cloud base level(number of model level)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PURV(KLON,KLEV) ! water vapor in updraft (kg/kg)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PURCI(KLON,KLEV) ! total condensate in updraft (kg/kg)
LOGICAL           ,INTENT(IN)    :: LD_OUVTRANS ! flag to compute convective
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUTEN(KLON,KLEV) ! convecctive u tendency (m/s^2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVTEN(KLON,KLEV) ! convecctive v tendency (m/s^2)
LOGICAL           ,INTENT(IN)    :: LD_OCHTRANS ! flag to compute convective
INTEGER(KIND=JPIM),INTENT(IN)    :: KCH1 ! number of species
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCH1(KLON,KLEV,KCH1) ! grid scale chemical species
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCH1TEN(KLON,KLEV,KCH1) ! chemical convective tendency

END SUBROUTINE CONVECTION_SHAL
END INTERFACE


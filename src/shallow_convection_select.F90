SUBROUTINE SHALLOW_CONVECTION_SELECT(CVP_SHAL, CVPEXT, CST, D, NSV, CONVPAR,    &
                                     ICONV, KICE, OSETTADJ, PTADJS,    &
                                     PPABST, PZZ, PTT, PRVT, PRCT,     &
                                     PRIT, PTTEN, PRVTEN, PRCTEN,      &
                                     PRITEN, KCLTOP, KCLBAS, PUMF,     &
                                     OCH1CONV, KCH1, PCH1, PCH1TEN,    &
                                     PRDOCP, PTHT, PSTHV, PSTHES,      &
                                     ISDPL, ISPBL, ISLCL, PSTHLCL,     &
                                     PSTLCL, PSRVLCL, PSWLCL, PSZLCL,  &
                                     PSTHVELCL, GTRIG1)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_CONVPAR, ONLY: CONVPAR_T
USE MODD_CONVPAR_SHAL, ONLY: CONVPAR_SHAL
USE MODD_CONVPAREXT, ONLY: CONVPAREXT
USE MODD_CST, ONLY: CST_T
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_NSV, ONLY: NSV_T

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
TYPE(CONVPAR_SHAL),              INTENT(IN)   :: CVP_SHAL
TYPE(CONVPAREXT),                INTENT(IN)   :: CVPEXT
TYPE(CST_T),                     INTENT(IN)   :: CST
TYPE(DIMPHYEX_T),                INTENT(IN)   :: D
TYPE(NSV_T),                     INTENT(IN)   :: NSV
TYPE(CONVPAR_T),                 INTENT(IN)   :: CONVPAR
INTEGER,                         INTENT(IN)   :: ICONV    ! number of convective columns 
INTEGER,                         INTENT(IN)   :: KICE     ! flag for ice ( 1 = yes,
                                                          !                0 = no ice )
LOGICAL,                         INTENT(IN)   :: OSETTADJ ! logical to set convective
                                                          ! adjustment time by user
REAL,                            INTENT(IN)   :: PTADJS   ! user defined adjustment time
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(IN)   :: PTT      ! grid scale temperature at t
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(IN)   :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(IN)   :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(IN)   :: PRIT     ! grid scale r_i "
                                                        ! velocity (m/s)
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(IN)   :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(IN)   :: PZZ      ! height of model layer (m)
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(INOUT):: PTTEN  ! convective temperature
                                                      ! tendency (K/s)
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
INTEGER, DIMENSION(D%NIT),       INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(D%NIT),       INTENT(INOUT):: KCLBAS ! cloud base level
                                                        ! they are given a value of
                                                        ! 0 if no convection
REAL, DIMENSION(D%NIT,D%NKT),    INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
!
LOGICAL,                         INTENT(IN)   :: OCH1CONV ! include tracer transport
INTEGER,                         INTENT(IN)   :: KCH1     ! number of species
REAL, DIMENSION(D%NIT,D%NKT,KCH1), INTENT(IN)   :: PCH1! grid scale chemical species
REAL, DIMENSION(D%NIT,D%NKT,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)

REAL,                             INTENT(IN)   :: PRDOCP   ! R_d/C_p
REAL, DIMENSION(D%NIT,D%NKT),     INTENT(IN)   :: PTHT, PSTHV, PSTHES  ! grid scale theta, theta_v
INTEGER, DIMENSION(D%NIT)  ,      INTENT(IN)   :: ISDPL   ! index for parcel departure level
INTEGER, DIMENSION(D%NIT)  ,      INTENT(IN)   :: ISPBL   ! index for source layer top
INTEGER, DIMENSION(D%NIT)  ,      INTENT(IN)   :: ISLCL   ! index for lifting condensation level
REAL, DIMENSION(D%NIT)     ,      INTENT(IN)   :: PSTHLCL ! updraft theta at LCL/L
REAL, DIMENSION(D%NIT)     ,      INTENT(IN)   :: PSTLCL  ! updraft temp. at LCL
REAL, DIMENSION(D%NIT)     ,      INTENT(IN)   :: PSRVLCL ! updraft rv at LCL
REAL, DIMENSION(D%NIT)     ,      INTENT(IN)   :: PSWLCL  ! updraft w at LCL
REAL, DIMENSION(D%NIT)     ,      INTENT(IN)   :: PSZLCL  ! LCL height
REAL, DIMENSION(D%NIT)     ,      INTENT(IN)   :: PSTHVELCL! envir. theta_v at LCL
LOGICAL, DIMENSION(D%NIT)  ,      INTENT(INOUT):: GTRIG1  ! logical mask for convection
!
!
!*       0.2   Declarations of local fixed memory variables :
!
INTEGER  :: JI, JL                  ! horizontal loop index
INTEGER  :: JN                      ! number of tracers
INTEGER  :: JK, JKM                 ! vertical loop index
INTEGER  :: IKB, IKE                ! vertical loop bounds
!
INTEGER, DIMENSION(ICONV)    :: ISORT
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER, DIMENSION(ICONV)    :: IDPL    ! index for parcel departure level
INTEGER, DIMENSION(ICONV)    :: IPBL    ! index for source layer top
INTEGER, DIMENSION(ICONV)    :: ILCL    ! index for lifting condensation level
INTEGER, DIMENSION(ICONV)    :: ICTL    ! index for cloud top level
INTEGER, DIMENSION(D%NIT)    :: IMINCTL ! min between index for cloud top level
                                        ! and lifting condensation level
!
! grid scale variables
REAL, DIMENSION(ICONV,D%NKT)  :: ZZ      ! height of model layer (m)
REAL, DIMENSION(ICONV,D%NKT)  :: ZPRES   ! grid scale pressure
REAL, DIMENSION(ICONV,D%NKT)  :: ZTT     ! temperature
REAL, DIMENSION(ICONV,D%NKT)  :: ZTH     ! grid scale theta
REAL, DIMENSION(ICONV,D%NKT)  :: ZTHV    ! grid scale theta_v
REAL, DIMENSION(ICONV,D%NKT)  :: ZTHES   ! grid scale saturated theta_e
REAL, DIMENSION(ICONV,D%NKT)  :: ZRV     ! grid scale water vapor (kg/kg)
REAL, DIMENSION(ICONV,D%NKT)  :: ZRC     ! grid scale cloud water (kg/kg)
REAL, DIMENSION(ICONV,D%NKT)  :: ZRI     ! grid scale cloud ice (kg/kg)
!
! updraft variables
REAL, DIMENSION(ICONV,D%NKT) :: ZUMF    ! updraft mass flux (kg/s)
REAL, DIMENSION(ICONV)       :: ZTHLCL  ! updraft theta at LCL
REAL, DIMENSION(ICONV)       :: ZTLCL   ! updraft temp. at LCL
REAL, DIMENSION(ICONV)       :: ZRVLCL  ! updraft rv at LCL
REAL, DIMENSION(ICONV)       :: ZWLCL   ! updraft w at LCL
REAL, DIMENSION(ICONV)       :: ZZLCL   ! LCL height
REAL, DIMENSION(ICONV)       :: ZTHVELCL! envir. theta_v at LCL
!
REAL, DIMENSION(ICONV,D%NKT)  :: ZTHC    ! conv. adj. grid scale theta
REAL, DIMENSION(ICONV,D%NKT)  :: ZRVC    ! conv. adj. grid scale r_w
REAL, DIMENSION(ICONV,D%NKT)  :: ZRCC    ! conv. adj. grid scale r_c
REAL, DIMENSION(ICONV,D%NKT)  :: ZRIC    ! conv. adj. grid scale r_i
!
! Chemical Tracers:
REAL, DIMENSION(D%NIT,D%NKT,KCH1)  :: ZPCH1TEN
!
TYPE(DIMPHYEX_T) :: ZD
!
!-------------------------------------------------------------------------------
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "shallow_convection_compute.h"

IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_SELECT',0,ZHOOK_HANDLE)

IKB = 1 + CVPEXT%JCVEXB
IKE = D%NKT - CVPEXT%JCVEXT
! Gather grid scale and updraft base variables in arrays using mask GTRIG
JL=1
DO JI=D%NIB,D%NIE
  IF(GTRIG1(JI))THEN
    ISORT(JL) = JI
    JL=JL+1
  ENDIF
ENDDO

DO JK = IKB, IKE
DO JI = D%NIB, ICONV
    GTRIG1(JI)    = GTRIG1(ISORT(JI))
    ZZ(JI,JK)     = PZZ(ISORT(JI),JK)
    ZPRES(JI,JK)  = PPABST(ISORT(JI),JK)
    ZTT(JI,JK)    = PTT(ISORT(JI),JK)
    ZTH(JI,JK)    = PTHT(ISORT(JI),JK)
    ZTHES(JI,JK)  = PSTHES(ISORT(JI),JK)
    ZRV(JI,JK)    = PRVT(ISORT(JI),JK)
    ZRC(JI,JK)    = PRCT(ISORT(JI),JK)
    ZRI(JI,JK)    = PRIT(ISORT(JI),JK)
    ZTHV(JI,JK)   = PSTHV(ISORT(JI),JK)
END DO
END DO

DO JI = D%NIB,ICONV
    IDPL(JI)      = ISDPL(ISORT(JI))
    IPBL(JI)      = ISPBL(ISORT(JI))
    ILCL(JI)      = ISLCL(ISORT(JI))
    ZTHLCL(JI)    = PSTHLCL(ISORT(JI))
    ZTLCL(JI)     = PSTLCL(ISORT(JI))
    ZRVLCL(JI)    = PSRVLCL(ISORT(JI))
    ZWLCL(JI)     = PSWLCL(ISORT(JI))
    ZZLCL(JI)     = PSZLCL(ISORT(JI))
    ZTHVELCL(JI)  = PSTHVELCL(ISORT(JI))
END DO

ZD=D
ZD%NIT=ICONV
ZD%NIE=ICONV
CALL SHALLOW_CONVECTION_COMPUTE(CVP_SHAL, CVPEXT, CST, ZD, NSV, CONVPAR, KICE,  &
                                OSETTADJ, PTADJS, ZPRES, ZZ, ZTT, ZRV, &
                                ZRC, ZRI, OCH1CONV, KCH1, PCH1, PRDOCP,&
                                ZTH, ZTHV, ZTHES, IDPL, IPBL, ILCL,    &
                                ZTHLCL, ZTLCL, ZRVLCL, ZWLCL, ZZLCL,   &
                                ZTHVELCL, GTRIG1, ZUMF, ZTHC, ZRVC,    &
                                ZRCC, ZRIC, ICTL, IMINCTL, ZPCH1TEN)
DO JK = IKB, IKE
DO JI = D%NIB, ICONV
    PTTEN(ISORT(JI),JK)   = ZTHC(JI,JK)
    PRVTEN(ISORT(JI),JK)  = ZRVC(JI,JK)
    PRCTEN(ISORT(JI),JK)  = ZRCC(JI,JK)
    PRITEN(ISORT(JI),JK)  = ZRIC(JI,JK)
END DO
END DO

DO JI = D%NIB, ICONV
    KCLTOP(ISORT(JI)) = ICTL(JI)
    KCLBAS(ISORT(JI)) = IMINCTL(JI)
END DO

IF ( OCH1CONV ) THEN
  JKM = IKE
  DO JN = 1, KCH1
    DO JK = IKB, IKE
      DO JI = D%NIB, ICONV
          PCH1TEN(ISORT(JI),JK,JN) = ZPCH1TEN(JI,JK,JN)
      END DO
    END DO
  END DO
END IF

DO JK = IKB, IKE
DO JI = D%NIB, ICONV
    PUMF(ISORT(JI),JK) = ZUMF(JI,JK)
END DO
END DO

IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_SELECT',1,ZHOOK_HANDLE)
END SUBROUTINE SHALLOW_CONVECTION_SELECT

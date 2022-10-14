SUBROUTINE SHALLOW_CONVECTION_SELECT( KLON, ICONV, KLEV, KIDIA, KFDIA, KICE, OSETTADJ,&
                                      PTADJS, PPABST, PZZ, PTT, PRVT,   &
                                      PRCT, PRIT, PTTEN, PRVTEN, PRCTEN,&
                                      PRITEN, KCLTOP, KCLBAS, PUMF,     &
                                      OCH1CONV, KCH1, PCH1, PCH1TEN,    &
                                      IKB, IKE, IFTSTEPS, ZRDOCP, ZTHT, &
                                      ZSTHV, ZSTHES, ISDPL, ISPBL,      &
                                      ISLCL, ZSTHLCL, ZSTLCL, ZSRVLCL,  &
                                      ZSWLCL, ZSZLCL, ZSTHVELCL, GTRIG1)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                         INTENT(IN)   :: KLON     ! horizontal dimension
INTEGER,                         INTENT(IN)   :: ICONV    ! number of convective columns 
INTEGER,                         INTENT(IN)   :: KLEV     ! vertical dimension
INTEGER,                         INTENT(IN)   :: KIDIA    ! value of the first point in x
INTEGER,                         INTENT(IN)   :: KFDIA    ! value of the last point in x
INTEGER,                         INTENT(IN)   :: KICE     ! flag for ice ( 1 = yes,
                                                          !                0 = no ice )
LOGICAL,                         INTENT(IN)   :: OSETTADJ ! logical to set convective
                                                          ! adjustment time by user
REAL,                            INTENT(IN)   :: PTADJS   ! user defined adjustment time
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)   :: PTT      ! grid scale temperature at t
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)   :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)   :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)   :: PRIT     ! grid scale r_i "
                                                          ! velocity (m/s)
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)   :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)   :: PZZ      ! height of model layer (m)
REAL, DIMENSION(KLON,KLEV),      INTENT(INOUT):: PTTEN  ! convective temperature
                                                        ! tendency (K/s)
REAL, DIMENSION(KLON,KLEV),      INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(KLON,KLEV),      INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(KLON,KLEV),      INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
INTEGER, DIMENSION(KLON),        INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(KLON),        INTENT(INOUT):: KCLBAS ! cloud base level
                                                        ! they are given a value of
                                                        ! 0 if no convection
REAL, DIMENSION(KLON,KLEV),      INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
!
LOGICAL,                         INTENT(IN)   :: OCH1CONV ! include tracer transport
INTEGER,                         INTENT(IN)   :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN)   :: PCH1! grid scale chemical species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)

INTEGER, INTENT(IN)                           :: IKB, IKE ! vertical loop bounds
INTEGER, INTENT(INOUT)                        :: IFTSTEPS ! only used for chemical tracers
REAL   , INTENT(IN)                           :: ZRDOCP   ! R_d/C_p
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)   :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta, theta_v
INTEGER, DIMENSION(KLON)  ,      INTENT(IN)   :: ISDPL   ! index for parcel departure level
INTEGER, DIMENSION(KLON)  ,      INTENT(IN)   :: ISPBL   ! index for source layer top
INTEGER, DIMENSION(KLON)  ,      INTENT(IN)   :: ISLCL   ! index for lifting condensation level
REAL, DIMENSION(KLON)     ,      INTENT(IN)   :: ZSTHLCL ! updraft theta at LCL/L
REAL, DIMENSION(KLON)     ,      INTENT(IN)   :: ZSTLCL  ! updraft temp. at LCL
REAL, DIMENSION(KLON)     ,      INTENT(IN)   :: ZSRVLCL ! updraft rv at LCL
REAL, DIMENSION(KLON)     ,      INTENT(IN)   :: ZSWLCL  ! updraft w at LCL
REAL, DIMENSION(KLON)     ,      INTENT(IN)   :: ZSZLCL  ! LCL height
REAL, DIMENSION(KLON)     ,      INTENT(IN)   :: ZSTHVELCL! envir. theta_v at LCL
LOGICAL, DIMENSION(KLON)  ,      INTENT(INOUT):: GTRIG1  ! logical mask for convection
!
!
!*       0.2   Declarations of local fixed memory variables :
!
INTEGER  :: JI, JL                  ! horizontal loop index
INTEGER  :: JN                      ! number of tracers
INTEGER  :: JK, JKM                 ! vertical loop index
!
LOGICAL, DIMENSION(KLON)           :: GTRIG  ! 2D logical mask for trigger test
INTEGER, DIMENSION(ICONV)          :: ISORT
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER, DIMENSION(ICONV)    :: IDPL    ! index for parcel departure level
INTEGER, DIMENSION(ICONV)    :: IPBL    ! index for source layer top
INTEGER, DIMENSION(ICONV)    :: ILCL    ! index for lifting condensation level
INTEGER, DIMENSION(ICONV)    :: ICTL    ! index for cloud top level
INTEGER, DIMENSION(KLON)     :: IMINCTL ! min between index for cloud top level
                                        ! and lifting condensation level
!
! grid scale variables
REAL, DIMENSION(ICONV,KLEV)  :: ZZ      ! height of model layer (m)
REAL, DIMENSION(ICONV,KLEV)  :: ZPRES   ! grid scale pressure
REAL, DIMENSION(ICONV,KLEV)  :: ZTT     ! temperature
REAL, DIMENSION(ICONV,KLEV)  :: ZTH     ! grid scale theta
REAL, DIMENSION(ICONV,KLEV)  :: ZTHV    ! grid scale theta_v
REAL, DIMENSION(ICONV,KLEV)  :: ZTHES   ! grid scale saturated theta_e
REAL, DIMENSION(ICONV,KLEV)  :: ZRV     ! grid scale water vapor (kg/kg)
REAL, DIMENSION(ICONV,KLEV)  :: ZRC     ! grid scale cloud water (kg/kg)
REAL, DIMENSION(ICONV,KLEV)  :: ZRI     ! grid scale cloud ice (kg/kg)
!
! updraft variables
REAL, DIMENSION(ICONV,KLEV)  :: ZUMF    ! updraft mass flux (kg/s)
REAL, DIMENSION(ICONV)       :: ZTHLCL  ! updraft theta at LCL
REAL, DIMENSION(ICONV)       :: ZTLCL   ! updraft temp. at LCL
REAL, DIMENSION(ICONV)       :: ZRVLCL  ! updraft rv at LCL
REAL, DIMENSION(ICONV)       :: ZWLCL   ! updraft w at LCL
REAL, DIMENSION(ICONV)       :: ZZLCL   ! LCL height
REAL, DIMENSION(ICONV)       :: ZTHVELCL! envir. theta_v at LCL
!
! closure variables
REAL, DIMENSION(ICONV)       :: ZTIMEC  ! advective time period
!
REAL, DIMENSION(ICONV,KLEV)  :: ZTHC    ! conv. adj. grid scale theta
REAL, DIMENSION(ICONV,KLEV)  :: ZRVC    ! conv. adj. grid scale r_w
REAL, DIMENSION(ICONV,KLEV)  :: ZRCC    ! conv. adj. grid scale r_c
REAL, DIMENSION(ICONV,KLEV)  :: ZRIC    ! conv. adj. grid scale r_i
!
! Chemical Tracers:
REAL, DIMENSION(ICONV,KLEV,KCH1) :: ZCH1    ! grid scale chemical specy (kg/kg)
REAL, DIMENSION(ICONV,KLEV,KCH1) :: ZCH1C   ! conv. adjust. chemical specy 1
REAL, DIMENSION(KLON,KLEV,KCH1)  :: ZPCH1TEN
!
!-------------------------------------------------------------------------------
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_SELECT',0,ZHOOK_HANDLE)

! Gather grid scale and updraft base variables in arrays using mask GTRIG
JL=1
DO JI=KIDIA,KFDIA
  IF(GTRIG1(JI))THEN
    ISORT(JL) = JI
    JL=JL+1
  ENDIF
ENDDO
GTRIG1 = .TRUE.

DO JK = IKB, IKE
DO JI = KIDIA, ICONV
    ZZ(JI,JK)     = PZZ(ISORT(JI),JK)
    ZPRES(JI,JK)  = PPABST(ISORT(JI),JK)
    ZTT(JI,JK)    = PTT(ISORT(JI),JK)
    ZTH(JI,JK)    = ZTHT(ISORT(JI),JK)
    ZTHES(JI,JK)  = ZSTHES(ISORT(JI),JK)
    ZRV(JI,JK)    = PRVT(ISORT(JI),JK)
    ZRC(JI,JK)    = PRCT(ISORT(JI),JK)
    ZRI(JI,JK)    = PRIT(ISORT(JI),JK)
    ZTHV(JI,JK)   = ZSTHV(ISORT(JI),JK)
END DO
END DO

DO JI = KIDIA,ICONV
    IDPL(JI)      = ISDPL(ISORT(JI))
    IPBL(JI)      = ISPBL(ISORT(JI))
    ILCL(JI)      = ISLCL(ISORT(JI))
    ZTHLCL(JI)    = ZSTHLCL(ISORT(JI))
    ZTLCL(JI)     = ZSTLCL(ISORT(JI))
    ZRVLCL(JI)    = ZSRVLCL(ISORT(JI))
    ZWLCL(JI)     = ZSWLCL(ISORT(JI))
    ZZLCL(JI)     = ZSZLCL(ISORT(JI))
    ZTHVELCL(JI)  = ZSTHVELCL(ISORT(JI))
END DO

CALL SHALLOW_CONVECTION_COMPUTE(ICONV, KLEV, KIDIA, ICONV, KICE,       &
                                OSETTADJ, PTADJS, ZPRES, ZZ, ZTT, ZRV, &
                                ZRC, ZRI, OCH1CONV, KCH1, PCH1, IKB,   &
                                IKE, IFTSTEPS, ZRDOCP, ZTH, ZTHV,      &
                                ZTHES, IDPL, IPBL, ILCL, ZTHLCL, ZTLCL,&
                                ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL, GTRIG1,&
                                ZTIMEC, ZCH1, ZCH1C, ZUMF, ZTHC, ZRVC, &
                                ZRCC, ZRIC, ICTL, IMINCTL, ZPCH1TEN)
DO JK = IKB, IKE
DO JI = KIDIA, ICONV
    PTTEN(ISORT(JI),JK)   = ZTHC(JI,JK)
    PRVTEN(ISORT(JI),JK)  = ZRVC(JI,JK)
    PRCTEN(ISORT(JI),JK)  = ZRCC(JI,JK)
    PRITEN(ISORT(JI),JK)  = ZRIC(JI,JK)
END DO
END DO

DO JI = KIDIA, ICONV
    KCLTOP(ISORT(JI)) = ICTL(JI)
    KCLBAS(ISORT(JI)) = IMINCTL(JI)
END DO

IF ( OCH1CONV ) THEN
  JKM = IKE
  DO JN = 1, KCH1
    DO JK = IKB, IKE
      DO JI = KIDIA, ICONV
          PCH1TEN(ISORT(JI),JK,JN) = ZPCH1TEN(JI,JK,JN)
      END DO
    END DO
  END DO
END IF

DO JK = IKB, IKE
DO JI = KIDIA, ICONV
    PUMF(ISORT(JI),JK) = ZUMF(JI,JK)
END DO
END DO

IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_SELECT',1,ZHOOK_HANDLE)
END SUBROUTINE SHALLOW_CONVECTION_SELECT

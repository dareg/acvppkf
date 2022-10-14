SUBROUTINE SHALLOW_CONVECTION_ALL( KLON, KLEV, KIDIA, KFDIA, KICE, OSETTADJ, PTADJS,  &
                                   PPABST, PZZ, PTT, PRVT, PRCT, PRIT,  &
                                   PTTEN, PRVTEN, PRCTEN, PRITEN,       &
                                   KCLTOP, KCLBAS, PUMF, OCH1CONV, KCH1,&
                                   PCH1, PCH1TEN, IKB, IKE, IFTSTEPS,   &
                                   ZRDOCP, ZTHT, ZSTHV, ZSTHES, ISDPL,  &
                                   ISPBL, ISLCL, ZSTHLCL, ZSTLCL,       &
                                   ZSRVLCL, ZSWLCL, ZSZLCL, ZSTHVELCL,  &
                                   GTRIG1)

USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                         INTENT(IN)   :: KLON     ! horizontal dimension
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
LOGICAL, DIMENSION(KLON)  ,      INTENT(IN)   :: GTRIG1  ! logical mask for convection
!
!
REAL, DIMENSION(KLON)              :: ZWORK2 ! work array
INTEGER  :: JI                      ! horizontal loop index
INTEGER  :: JN                      ! number of tracers
INTEGER  :: JK, JKM                 ! vertical loop index
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER, DIMENSION(KLON)  :: ICTL    ! index for cloud top level
INTEGER, DIMENSION(KLON)  :: IMINCTL ! min between index for cloud top level
                                     ! and lifting condensation level
!
! updraft variables
REAL, DIMENSION(KLON,KLEV)  :: ZUMF    ! updraft mass flux (kg/s)
!
! closure variables
REAL, DIMENSION(KLON)       :: ZTIMEC  ! advective time period
!
REAL, DIMENSION(KLON,KLEV)  :: ZTHC    ! conv. adj. grid scale theta
REAL, DIMENSION(KLON,KLEV)  :: ZRVC    ! conv. adj. grid scale r_w
REAL, DIMENSION(KLON,KLEV)  :: ZRCC    ! conv. adj. grid scale r_c
REAL, DIMENSION(KLON,KLEV)  :: ZRIC    ! conv. adj. grid scale r_i
!
! Chemical Tracers:
REAL, DIMENSION(KLON,KLEV,KCH1):: ZCH1    ! grid scale chemical specy (kg/kg)
REAL, DIMENSION(KLON,KLEV,KCH1):: ZCH1C   ! conv. adjust. chemical specy 1
REAL, DIMENSION(KLON,KLEV,KCH1):: ZPCH1TEN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_ALL',0,ZHOOK_HANDLE)

CALL SHALLOW_CONVECTION_COMPUTE(KLON, KLEV, KIDIA, KFDIA, KICE,        &
                                OSETTADJ, PTADJS, PPABST, PZZ, PTT,    &
                                PRVT, PRCT, PRIT, OCH1CONV, KCH1, PCH1,&
                                IKB, IKE, IFTSTEPS, ZRDOCP, ZTHT,      &
                                ZSTHV, ZSTHES, ISDPL, ISPBL, ISLCL,    &
                                ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL,      &
                                ZSZLCL, ZSTHVELCL, GTRIG1, ZTIMEC,     &
                                ZCH1, ZCH1C, ZUMF, ZTHC, ZRVC, ZRCC,   &
                                ZRIC, ICTL, IMINCTL, ZPCH1TEN)
DO JK = IKB, IKE
DO JI = KIDIA,KFDIA
  IF(GTRIG1(JI))THEN
    PTTEN(JI,JK)   = ZTHC(JI,JK)
    PRVTEN(JI,JK)  = ZRVC(JI,JK)
    PRCTEN(JI,JK)  = ZRCC(JI,JK)
    PRITEN(JI,JK)  = ZRIC(JI,JK)
  ENDIF
END DO
END DO

DO JI = KIDIA,KFDIA
  IF(GTRIG1(JI))THEN
    KCLTOP(JI) = ICTL(JI)
    KCLBAS(JI) = IMINCTL(JI)
  ENDIF
END DO

IF ( OCH1CONV ) THEN
  JKM = IKE
  DO JN = 1, KCH1
    DO JK = IKB, IKE
      DO JI = KIDIA,KFDIA
        IF(GTRIG1(JI))THEN
          PCH1TEN(JI,JK,JN) = ZPCH1TEN(JI,JK,JN)
        ENDIF
      END DO
    END DO
  END DO
END IF

ZWORK2(:) = 1.
DO JK = IKB, IKE
DO JI = KIDIA,KFDIA
  IF(GTRIG1(JI))THEN
    IF ( KCLTOP(JI) <= IKB+1 ) ZWORK2(JI) = 0.
    PUMF(JI,JK) = ZUMF(JI,JK) * ZWORK2(JI)
  ENDIF
END DO
END DO
IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_ALL',1,ZHOOK_HANDLE)
END SUBROUTINE SHALLOW_CONVECTION_ALL

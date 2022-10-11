!     ######spl
    SUBROUTINE SHALLOW_CONVECTION( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
                                   KICE, OSETTADJ, PTADJS,               &
                                   PPABST, PZZ, PTKECLS,                          &
                                   PTT, PRVT, PRCT, PRIT, PWT,                    &
                                   PTTEN, PRVTEN, PRCTEN, PRITEN,                 &
                                   KCLTOP, KCLBAS, PUMF,                          &
                                   OCH1CONV, KCH1, PCH1, PCH1TEN                  )
    USE PARKIND1, ONLY : JPRB
    USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!   ###############################################################################
!
!!**** Monitor routine to compute all convective tendencies by calls
!!     of several subroutines.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      tendencies. The routine first prepares all necessary grid-scale
!!      variables. The final convective tendencies are then computed by
!!      calls of different subroutines.
!!
!!
!!**  METHOD
!!    ------
!!      We start by selecting convective columns in the model domain through
!!      the call of routine TRIGGER_FUNCT. Then, we allocate memory for the
!!      convection updraft and downdraft variables and gather the grid scale
!!      variables in convective arrays.
!!      The updraft and downdraft computations are done level by level starting
!!      at the  bottom and top of the domain, respectively.
!!      All computations are done on MNH thermodynamic levels. The depth
!!      of the current model layer k is defined by DP(k)=P(k-1)-P(k)
!!
!!
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_TRIGGER_SHAL
!!    CONVECT_SATMIXRATIO
!!    CONVECT_UPDRAFT_SHAL
!!        CONVECT_CONDENS
!!        CONVECT_MIXING_FUNCT
!!    CONVECT_CLOSURE_SHAL
!!        CONVECT_CLOSURE_THRVLCL
!!        CONVECT_CLOSURE_ADJUST_SHAL
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                   ! gravity constant
!!          XPI                  ! number Pi
!!          XP00                 ! reference pressure
!!          XRD, XRV             ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV           ! specific heat for dry air and water vapor
!!          XRHOLW               ! density of liquid water
!!          XALPW, XBETAW, XGAMW ! constants for water saturation pressure
!!          XTT                  ! triple point temperature
!!          XLVTT, XLSTT         ! vaporization, sublimation heat constant
!!          XCL, XCI             ! specific heat for liquid water and ice
!!
!!      Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT       ! extra levels on the vertical boundaries
!!
!!      Module MODD_CONVPAR
!!          XA25                 ! reference grid area
!!          XCRAD                ! cloud radius
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Bechtold, 1997 : Meso-NH scientific  documentation (31 pp)
!!      Fritsch and Chappell, 1980, J. Atmos. Sci., Vol. 37, 1722-1761.
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol. 47, 2784-2801.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol. 24, 165-170.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Peter Bechtold 15/11/96 replace theta_il by enthalpy
!!         "        10/12/98 changes for ARPEGE
!!         "        01/01/02 Apply conservation correction
!!   F Bouyssel     05/11/08 Modifications for reproductibility
!!   E. Bazile      20/07/09 Input of TKECLS.
!!   F. Bouyssel    08/11/13 Modifications for reproductibility
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : XALPW, XBETAW, XCPD, XGAMW, XP00, XRD, XRV
USE MODD_CONVPAREXT, ONLY : JCVEXB, JCVEXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                    INTENT(IN) :: KIDIA    ! value of the first point in x
INTEGER,                    INTENT(IN) :: KFDIA    ! value of the last point in x
INTEGER,                    INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                  ! KBDIA that is at least 1
INTEGER,                    INTENT(IN) :: KTDIA    ! vertical computations can be
                                                   ! limited to KLEV + 1 - KTDIA
                                                   ! default=1
                                                   ! scheme
INTEGER,                    INTENT(IN) :: KICE     ! flag for ice ( 1 = yes,
                                                   !                0 = no ice )
LOGICAL,                    INTENT(IN) :: OSETTADJ ! logical to set convective
                                                   ! adjustment time by user
REAL,                       INTENT(IN) :: PTADJS   ! user defined adjustment time
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical
                                                   ! velocity (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m)
REAL, DIMENSION(KLON),      INTENT(IN) :: PTKECLS  ! TKE in the CLS  (m2/s2)
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperature
                                                   ! tendency (K/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLBAS ! cloud base level
                                                   ! they are given a value of
                                                   ! 0 if no convection
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
!
LOGICAL,                    INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                    INTENT(IN) :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN) :: PCH1! grid scale chemical species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)
!
!
!*       0.2   Declarations of local fixed memory variables :
!
INTEGER  :: IKB, IKE                ! vertical loop bounds
INTEGER  :: JI                      ! horizontal loop index
INTEGER  :: JK                      ! vertical loop index
INTEGER  :: IFTSTEPS                ! only used for chemical tracers
INTEGER  :: ICONV
REAL     :: ZEPS, ZEPSA             ! R_d / R_v, R_v / R_d
REAL     :: ZRDOCP                  ! R_d/C_p
!
REAL, DIMENSION(KLON,KLEV)         :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta, theta_v
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER, DIMENSION(KLON)  :: ISDPL   ! index for parcel departure level
INTEGER, DIMENSION(KLON)  :: ISPBL   ! index for source layer top
INTEGER, DIMENSION(KLON)  :: ISLCL   ! index for lifting condensation level
REAL, DIMENSION(KLON)     :: ZSTHLCL ! updraft theta at LCL/L
REAL, DIMENSION(KLON)     :: ZSTLCL  ! updraft temp. at LCL
REAL, DIMENSION(KLON)     :: ZSRVLCL ! updraft rv at LCL
REAL, DIMENSION(KLON)     :: ZSWLCL  ! updraft w at LCL
REAL, DIMENSION(KLON)     :: ZSZLCL  ! LCL height
REAL, DIMENSION(KLON)     :: ZSTHVELCL! envir. theta_v at LCL
!
LOGICAL, DIMENSION(KLON)    :: GTRIG1  ! logical mask for convection
REAL                        :: ZES     ! saturation vapor mixng ratio
!
!-------------------------------------------------------------------------------
!
!
!*       0.3    Compute loop bounds
!               -------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION',0,ZHOOK_HANDLE)
JCVEXB = MAX( 0, KBDIA - 1 )
IKB = 1 + JCVEXB
JCVEXT = MAX( 0, KTDIA - 1)
IKE = KLEV - JCVEXT
!
!*       0.7    Reset convective tendencies to zero if convective
!               counter becomes negative
!               -------------------------------------------------
!
PTTEN(:,:)  = 0.
PRVTEN(:,:) = 0.
PRCTEN(:,:) = 0.
PRITEN(:,:) = 0.
PUMF(:,:)   = 0.
KCLTOP(:)  = 0
KCLBAS(:)  = 0

IF ( OCH1CONV ) THEN
  PCH1TEN(:,:,:) = 0.
END IF
!
!
!*       1.     Initialize  local variables
!               ----------------------------
!
ZEPS   = XRD / XRV
ZEPSA  = XRV / XRD
ZRDOCP = XRD / XCPD
!
!-------------------------------------------------------------------------------
!
!*       1.1    Set up grid scale theta, theta_v, theta_es
!               ------------------------------------------
!
ZTHT(:,:) = 300.
ZSTHV(:,:)= 300.
ZSTHES(:,:)= 400.
DO JK = IKB, IKE
DO JI = KIDIA, KFDIA
  IF ( PPABST(JI,JK) > 40.E2 ) THEN
    ZTHT(JI,JK)  = PTT(JI,JK) * ( XP00 / PPABST(JI,JK) ) ** ZRDOCP
    ZSTHV(JI,JK) = ZTHT(JI,JK) * ( 1. + ZEPSA * PRVT(JI,JK) ) /              &
                   ( 1. + PRVT(JI,JK) + PRCT(JI,JK) + PRIT(JI,JK) )
!
        ! use conservative Bolton (1980) formula for theta_e
        ! it is used to compute CAPE for undilute parcel ascent
        ! For economical reasons we do not use routine CONVECT_SATMIXRATIO here
!
    ZES = EXP( XALPW - XBETAW / PTT(JI,JK) - XGAMW * LOG( PTT(JI,JK) ) )
    ZES = MIN( 1., ZEPS * ZES / ( PPABST(JI,JK) - ZES ) )
    ZSTHES(JI,JK) = PTT(JI,JK) * ( ZTHT(JI,JK) / PTT(JI,JK) ) **             &
              ( 1. - 0.28 * ZES ) * EXP( ( 3374.6525 / PTT(JI,JK) - 2.5403 ) &
                                        * ZES * ( 1. + 0.81 * ZES ) )
  END IF
END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*       2.     Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------
!
!*       2.3    Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------
!
ISLCL(:) = MAX( IKB, 2 )   ! initialize DPL PBL and LCL
ISDPL(:) = IKB
ISPBL(:) = IKB
!
CALL CONVECT_TRIGGER_SHAL(  KLON, KLEV, KIDIA, KFDIA,                    &
                            PPABST, ZTHT, ZSTHV, ZSTHES,                 &
                            PRVT, PWT, PZZ, PTKECLS,             &
                            ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSZLCL, &
                            ZSTHVELCL, ISLCL, ISDPL, ISPBL, GTRIG1)
ICONV = COUNT(GTRIG1(:))
IF(ICONV==0)THEN
  ! Do nothing if there are no selected columns
ELSE IF (ICONV < KLON/2) THEN
  CALL SHALLOW_CONVECTION_SELECT( KLON, ICONV, KLEV, KIDIA, KFDIA, KICE, OSETTADJ,&
                                  PTADJS, PPABST, PZZ, PTT, PRVT,   &
                                  PRCT, PRIT, PTTEN, PRVTEN, PRCTEN,&
                                  PRITEN, KCLTOP, KCLBAS, PUMF,     &
                                  OCH1CONV, KCH1, PCH1, PCH1TEN,    &
                                  IKB, IKE, IFTSTEPS, ZRDOCP, ZTHT, &
                                  ZSTHV, ZSTHES, ISDPL, ISPBL,      &
                                  ISLCL, ZSTHLCL, ZSTLCL, ZSRVLCL,  &
                                  ZSWLCL, ZSZLCL, ZSTHVELCL, GTRIG1)
ELSE
  CALL SHALLOW_CONVECTION_ALL( KLON, KLEV, KIDIA, KFDIA, KICE, OSETTADJ, PTADJS,  &
                               PPABST, PZZ, PTT, PRVT, PRCT, PRIT,  &
                               PTTEN, PRVTEN, PRCTEN, PRITEN,       &
                               KCLTOP, KCLBAS, PUMF, OCH1CONV, KCH1,&
                               PCH1, PCH1TEN, IKB, IKE, IFTSTEPS,   &
                               ZRDOCP, ZTHT, ZSTHV, ZSTHES, ISDPL,  &
                               ISPBL, ISLCL, ZSTHLCL, ZSTLCL,       &
                               ZSRVLCL, ZSWLCL, ZSZLCL, ZSTHVELCL,  &
                               GTRIG1)
ENDIF
IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION',1,ZHOOK_HANDLE)
CONTAINS
!
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
LOGICAL, DIMENSION(KLON)  ,      INTENT(INOUT):: GTRIG1  ! logical mask for convection
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
                                ZRIC, ICTL)
DO JK = IKB, IKE
DO JI = KIDIA,KFDIA
  IF(GTRIG1(JI) .EQV. .TRUE.)THEN
    PTTEN(JI,JK)   = ZTHC(JI,JK)
    PRVTEN(JI,JK)  = ZRVC(JI,JK)
    PRCTEN(JI,JK)  = ZRCC(JI,JK)
    PRITEN(JI,JK)  = ZRIC(JI,JK)
  ENDIF
END DO
END DO

DO JI = KIDIA,KFDIA
  IF(GTRIG1(JI) .EQV. .TRUE.)THEN
    KCLTOP(JI) = ICTL(JI)
    KCLBAS(JI) = MIN(ISLCL(JI), ICTL(JI))
  ENDIF
END DO

IF ( OCH1CONV ) THEN
  JKM = IKE
  DO JN = 1, KCH1
    DO JK = IKB, IKE
      DO JI = KIDIA,KFDIA
        IF(GTRIG1(JI) .EQV. .TRUE.)THEN
          PCH1TEN(JI,JK,JN) = (ZCH1C(JI,JK,JN)-ZCH1(JI,JK,JN) ) / ZTIMEC(JI)
        ENDIF
      END DO
    END DO
  END DO
END IF

ZWORK2(:) = 1.
DO JK = IKB, IKE
DO JI = KIDIA,KFDIA
  IF(GTRIG1(JI) .EQV. .TRUE.)THEN
    IF ( KCLTOP(JI) <= IKB+1 ) ZWORK2(JI) = 0.
    PUMF(JI,JK) = ZUMF(JI,JK) * ZWORK2(JI)
  ENDIF
END DO
END DO
IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_ALL',1,ZHOOK_HANDLE)
END SUBROUTINE SHALLOW_CONVECTION_ALL
!
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
INTEGER, DIMENSION(KLON)           :: ISORT
REAL, DIMENSION(KLON)              :: ZWORK2 ! work array
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER, DIMENSION(ICONV)    :: IDPL    ! index for parcel departure level
INTEGER, DIMENSION(ICONV)    :: IPBL    ! index for source layer top
INTEGER, DIMENSION(ICONV)    :: ILCL    ! index for lifting condensation level
INTEGER, DIMENSION(ICONV)    :: ICTL    ! index for cloud top level
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
!
!-------------------------------------------------------------------------------
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_SELECT',0,ZHOOK_HANDLE)

! Gather grid scale and updraft base variables in arrays using mask GTRIG
GTRIG(:) = GTRIG1(:)
GTRIG1 = .TRUE.
JL=1
DO JI=KIDIA,KFDIA
  IF(GTRIG(JI))THEN
    ISORT(JI) = JL
    JL=JL+1
  ENDIF
ENDDO

DO JK = IKB, IKE
DO JI = KIDIA, KFDIA
  IF(GTRIG(JI))THEN
    ZZ(ISORT(JI),JK)     = PZZ(JI,JK)
    ZPRES(ISORT(JI),JK)  = PPABST(JI,JK)
    ZTT(ISORT(JI),JK)    = PTT(JI,JK)
    ZTH(ISORT(JI),JK)    = ZTHT(JI,JK)
    ZTHES(ISORT(JI),JK)  = ZSTHES(JI,JK)
    ZRV(ISORT(JI),JK)    = MAX( 0., PRVT(JI,JK) )
    ZRC(ISORT(JI),JK)    = MAX( 0., PRCT(JI,JK) )
    ZRI(ISORT(JI),JK)    = MAX( 0., PRIT(JI,JK) )
    ZTHV(ISORT(JI),JK)   = ZSTHV(JI,JK)
  ENDIF
END DO
END DO

DO JI = KIDIA,KFDIA
  IF(GTRIG(JI))THEN
    IDPL(ISORT(JI))      = ISDPL(JI)
    IPBL(ISORT(JI))      = ISPBL(JI)
    ILCL(ISORT(JI))      = ISLCL(JI)
    ZTHLCL(ISORT(JI))    = ZSTHLCL(JI)
    ZTLCL(ISORT(JI))     = ZSTLCL(JI)
    ZRVLCL(ISORT(JI))    = ZSRVLCL(JI)
    ZWLCL(ISORT(JI))     = ZSWLCL(JI)
    ZZLCL(ISORT(JI))     = ZSZLCL(JI)
    ZTHVELCL(ISORT(JI))  = ZSTHVELCL(JI)
  ENDIF
END DO

CALL SHALLOW_CONVECTION_COMPUTE(ICONV, KLEV, KIDIA, ICONV, KICE,       &
                                OSETTADJ, PTADJS, ZPRES, ZZ, ZTT, ZRV, &
                                ZRC, ZRI, OCH1CONV, KCH1, PCH1, IKB,   &
                                IKE, IFTSTEPS, ZRDOCP, ZTH, ZTHV,      &
                                ZTHES, IDPL, IPBL, ILCL, ZTHLCL, ZTLCL,&
                                ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL, GTRIG1,&
                                ZTIMEC, ZCH1, ZCH1C, ZUMF, ZTHC, ZRVC, &
                                ZRCC, ZRIC, ICTL)
DO JK = IKB, IKE
DO JI = KIDIA, KFDIA
  IF(GTRIG(JI))THEN
    PTTEN(JI,JK)   = ZTHC(ISORT(JI),JK)
    PRVTEN(JI,JK)  = ZRVC(ISORT(JI),JK)
    PRCTEN(JI,JK)  = ZRCC(ISORT(JI),JK)
    PRITEN(JI,JK)  = ZRIC(ISORT(JI),JK)
  ENDIF
END DO
END DO

ILCL(:) = MIN( ILCL(:), ICTL(:) )
DO JI = KIDIA, KFDIA
  IF(GTRIG(JI))THEN
    KCLTOP(JI) = ICTL(ISORT(JI))
    KCLBAS(JI) = ILCL(ISORT(JI))
  ENDIF
END DO

IF ( OCH1CONV ) THEN
  JKM = IKE
  DO JN = 1, KCH1
    DO JK = IKB, IKE
      DO JI = KIDIA, KFDIA
        IF(GTRIG(JI))THEN
          PCH1TEN(JI,JK,JN) = (ZCH1C(ISORT(JI),JK,JN)-ZCH1(ISORT(JI),JK,JN) ) / ZTIMEC(ISORT(JI))
        ENDIF
      END DO
    END DO
  END DO
END IF

ZWORK2(:) = 1.
DO JK = IKB, IKE
DO JI = KIDIA, KFDIA
  IF(GTRIG(JI))THEN
    IF ( KCLTOP(JI) <= IKB+1 ) ZWORK2(JI) = 0.
    PUMF(JI,JK) = ZUMF(ISORT(JI),JK) * ZWORK2(JI)
  ENDIF
END DO
END DO

IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_SELECT',1,ZHOOK_HANDLE)
END SUBROUTINE SHALLOW_CONVECTION_SELECT
!
SUBROUTINE SHALLOW_CONVECTION_COMPUTE( KLON, KLEV, KIDIA, KFDIA, KICE, OSETTADJ, PTADJS,  &
                                   PPABST, PZZ, PTT, PRVT, PRCT, PRIT,  &
                                   OCH1CONV, KCH1,&
                                   PCH1, IKB, IKE, IFTSTEPS,   &
                                   ZRDOCP, ZTHT, ZSTHV, ZSTHES, ISDPL,  &
                                   ISPBL, ISLCL, ZSTHLCL, ZSTLCL,       &
                                   ZSRVLCL, ZSWLCL, ZSZLCL, ZSTHVELCL,  &
                                   GTRIG1, ZTIMEC, ZCH1, ZCH1C, ZUMF,   &
                                   ZTHC, ZRVC, ZRCC, ZRIC, ICTL &
                                   )

USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_CST, ONLY : XCI, XCL, XCPD, XCPV, XG, XLSTT, XLVTT, XP00, XTT
USE MODD_NSV, ONLY : NSV_LGBEG,NSV_LGEND
USE MODD_CONVPAR_SHAL, ONLY : LLSMOOTH, XA25, XCTIME_SHAL

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                         INTENT(IN)  :: KLON     ! horizontal dimension
INTEGER,                         INTENT(IN)  :: KLEV     ! vertical dimension
INTEGER,                         INTENT(IN)  :: KIDIA    ! value of the first point in x
INTEGER,                         INTENT(IN)  :: KFDIA    ! value of the last point in x
INTEGER,                         INTENT(IN)  :: KICE     ! flag for ice ( 1 = yes,
                                                         !                0 = no ice )
LOGICAL,                         INTENT(IN)  :: OSETTADJ ! logical to set convective
                                                         ! adjustment time by user
REAL,                            INTENT(IN)  :: PTADJS   ! user defined adjustment time
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)  :: PTT      ! grid scale temperature at t
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)  :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)  :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)  :: PRIT     ! grid scale r_i "
                                                         ! velocity (m/s)
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)  :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)  :: PZZ      ! height of model layer (m)
                                                       ! tendency (K/s)
                                                       ! they are given a value of
                                                       ! 0 if no convection
!
LOGICAL,                         INTENT(IN)  :: OCH1CONV ! include tracer transport
INTEGER,                         INTENT(IN)  :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN)  :: PCH1     ! grid scale chemical species

INTEGER, INTENT(IN)                          :: IKB, IKE ! vertical loop bounds
INTEGER, INTENT(INOUT)                       :: IFTSTEPS ! only used for chemical tracers
REAL   , INTENT(IN)                          :: ZRDOCP   ! R_d/C_p
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)  :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta, theta_v
INTEGER, DIMENSION(KLON)  ,      INTENT(IN)  :: ISDPL   ! index for parcel departure level
INTEGER, DIMENSION(KLON)  ,      INTENT(IN)  :: ISPBL   ! index for source layer top
INTEGER, DIMENSION(KLON)  ,      INTENT(IN)  :: ISLCL   ! index for lifting condensation level
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: ZSTHLCL ! updraft theta at LCL/L
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: ZSTLCL  ! updraft temp. at LCL
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: ZSRVLCL ! updraft rv at LCL
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: ZSWLCL  ! updraft w at LCL
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: ZSZLCL  ! LCL height
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: ZSTHVELCL! envir. theta_v at LCL
LOGICAL, DIMENSION(KLON)  ,      INTENT(IN)  :: GTRIG1  ! logical mask for convection
REAL, DIMENSION(KLON),           INTENT(OUT) :: ZTIMEC  ! advective time period
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(OUT) :: ZCH1    ! grid scale chemical specy (kg/kg)
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(OUT) :: ZCH1C   ! conv. adjust. chemical specy 1
REAL, DIMENSION(KLON,KLEV),      INTENT(OUT) :: ZUMF    ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV),      INTENT(OUT) :: ZTHC    ! conv. adj. grid scale theta
REAL, DIMENSION(KLON,KLEV),      INTENT(OUT) :: ZRVC    ! conv. adj. grid scale r_w
REAL, DIMENSION(KLON,KLEV),      INTENT(OUT) :: ZRCC    ! conv. adj. grid scale r_c
REAL, DIMENSION(KLON,KLEV),      INTENT(OUT) :: ZRIC    ! conv. adj. grid scale r_i
INTEGER, DIMENSION(KLON),        INTENT(OUT) :: ICTL    ! index for cloud top level
!
!
REAL, DIMENSION(KLON)              :: ZWORK2, ZWORK2B ! work array
REAL                               :: ZW1     ! work variable
INTEGER  :: JI                      ! horizontal loop index
INTEGER  :: JN                      ! number of tracers
INTEGER  :: JK, JKM, JKP            ! vertical loop index
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER, DIMENSION(KLON)  :: IETL    ! index for zero buoyancy level
INTEGER, DIMENSION(KLON)  :: ILFS    ! index for level of free sink
!
! grid scale variables
REAL, DIMENSION(KLON,KLEV)  :: ZDPRES  ! pressure difference between
                                              ! bottom and top of layer (Pa)
REAL, DIMENSION(KLON,KLEV)  :: ZTHL    ! grid scale enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV)  :: ZRW     ! grid scale total water (kg/kg)
!
! updraft variables
REAL, DIMENSION(KLON,KLEV)  :: ZUER    ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV)  :: ZUDR    ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV)  :: ZUTHL   ! updraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV)  :: ZUTHV   ! updraft theta_v (K)
REAL, DIMENSION(KLON,KLEV)  :: ZURW    ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV)  :: ZURC    ! updraft cloud water (kg/kg)
REAL, DIMENSION(KLON,KLEV)  :: ZURI    ! updraft cloud ice   (kg/kg)
REAL, DIMENSION(KLON)       :: ZCAPE   ! available potent. energy
!
! downdraft variables
REAL, DIMENSION(KLON,KLEV)  :: ZDMF    ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV)  :: ZDER    ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV)  :: ZDDR    ! downdraft detrainment (kg/s)
!
! closure variables
REAL, DIMENSION(KLON,KLEV)  :: ZLMASS  ! mass of model layer (kg)
!
REAL, DIMENSION(KLON,KLEV)  :: ZWSUB   ! envir. compensating subsidence (Pa/s)
!
LOGICAL, DIMENSION(KLON)    :: GTRIG2  ! logical mask for convection
REAL, DIMENSION(KLON)       :: ZCPH    ! specific heat C_ph
REAL, DIMENSION(KLON)       :: ZLV, ZLS! latent heat of vaporis., sublim.
!
! Chemical Tracers:
REAL, DIMENSION(KLON,KCH1)     :: ZWORK3  ! conv. adjust. chemical specy 1
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_COMPUTE',0,ZHOOK_HANDLE)
ZDPRES = 0.0
ZTHL  = 0.0
ZRW   = 0.0
!
!*           3.2    Compute pressure difference
!                   ---------------------------------------------------
!
ZDPRES(:,IKB) = 0.
DO JK = IKB + 1, IKE
  ZDPRES(:,JK)  = PPABST(:,JK-1) - PPABST(:,JK)
END DO
!
!*           3.3   Compute environm. enthalpy and total water = r_v + r_i + r_c
!                  ----------------------------------------------------------
!
DO JK = IKB, IKE, 1
  ZRW(:,JK)  = MAX(0., PRVT(:,JK)) + MAX(0., PRCT(:,JK)) + MAX(0., PRIT(:,JK))
  ZCPH(:)    = XCPD + XCPV * ZRW(:,JK)
  ZLV(:)     = XLVTT + ( XCPV - XCL ) * ( PTT(:,JK) - XTT ) ! compute L_v
  ZLS(:)     = XLSTT + ( XCPV - XCI ) * ( PTT(:,JK) - XTT ) ! compute L_i
  ZTHL(:,JK) = ZCPH(:) * PTT(:,JK) + ( 1. + ZRW(:,JK) ) * XG * PZZ(:,JK) &
               - ZLV(:) * MAX(0., PRCT(:,JK)) - ZLS(:) * MAX(0., PRIT(:,JK))
END DO
!
!-------------------------------------------------------------------------------
!
!*           4.     Compute updraft properties
!                   ----------------------------
!
!*           4.1    Set mass flux at LCL ( here a unit mass flux with w = 1 m/s )
!                   -------------------------------------------------------------
!
CALL CONVECT_UPDRAFT_SHAL( KLON, KLEV,                                     &
                           KICE, PPABST, ZDPRES, PZZ, ZTHL, ZSTHV, ZSTHES, ZRW, &
                           ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSZLCL, ZSTHVELCL,   &
                           XA25 * 1.E-3, GTRIG2, ISLCL, ISDPL, ISPBL,                &
                           ZUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZURW,            &
                           ZURC, ZURI, ZCAPE, ICTL, IETL, GTRIG1                    )

ZDMF(:,:) = 0.
ZDER(:,:) = 0.
ZDDR(:,:) = 0.
ILFS(:)   = IKB
DO JK = IKB, IKE
  ZLMASS(:,JK)  = XA25 * ZDPRES(:,JK) / XG  ! mass of model layer
END DO
ZLMASS(:,IKB) = ZLMASS(:,IKB+1)
!
!-------------------------------------------------------------------------------
!
!*           5.     Compute downdraft properties
!                   ----------------------------
!
  ZTIMEC(:) = XCTIME_SHAL
  IF ( OSETTADJ ) ZTIMEC(:) = PTADJS
!
!*           7.     Determine adjusted environmental values assuming
!                   that all available buoyant energy must be removed
!                   within an advective time step ZTIMEC.
!                   ---------------------------------------------------
!
  CALL CONVECT_CLOSURE_SHAL( KLON, KLEV,                         &
                             PPABST, ZDPRES, PZZ, XA25, ZLMASS,    &
                             ZTHL, ZTHT, ZRW, PRCT, PRIT, GTRIG2,    &
                             ZTHC, ZRVC, ZRCC, ZRIC, ZWSUB,       &
                             ISLCL, ISDPL, ISPBL, ICTL,              &
                             ZUMF, ZUER, ZUDR, ZUTHL, ZURW,       &
                             ZURC, ZURI, ZCAPE, ZTIMEC, IFTSTEPS  )
!
!-------------------------------------------------------------------------------
!
!*           8.     Determine the final grid-scale (environmental) convective
!                   tendencies and set convective counter
!                   --------------------------------------------------------
!
!
!*           8.1    Grid scale tendencies
!                   ---------------------
!
          ! in order to save memory, the tendencies are temporarily stored
          ! in the tables for the adjusted grid-scale values
!
DO JK = IKB, IKE
   ZTHC(:,JK) = ( ZTHC(:,JK) - ZTHT(:,JK) ) / ZTIMEC(:)             &
     * ( PPABST(:,JK) / XP00 ) ** ZRDOCP ! change theta in temperature
   ZRVC(:,JK) = ( ZRVC(:,JK) - ZRW(:,JK) + MAX(0., PRCT(:,JK)) + MAX(0., PRIT(:,JK)) ) &
                                        / ZTIMEC(:)

   ZRCC(:,JK) = ( ZRCC(:,JK) - MAX(0., PRCT(:,JK)) ) / ZTIMEC(:)
   ZRIC(:,JK) = ( ZRIC(:,JK) - MAX(0., PRIT(:,JK)) ) / ZTIMEC(:)
END DO
!
!
!*           8.2    Apply conservation correction
!                   -----------------------------
!
          ! adjustment at cloud top to smooth possible discontinuous profiles at PBL inversions
          ! (+ - - tendencies for moisture )
!
!
IF (LLSMOOTH) THEN
  DO JI = KIDIA,KFDIA
     JK = ICTL(JI)
     JKM= MAX(2,ICTL(JI)-1)
     JKP= MAX(2,ICTL(JI)-2)
     ZRVC(JI,JKM) = ZRVC(JI,JKM) + .5 * ZRVC(JI,JK)
     ZRCC(JI,JKM) = ZRCC(JI,JKM) + .5 * ZRCC(JI,JK)
     ZRIC(JI,JKM) = ZRIC(JI,JKM) + .5 * ZRIC(JI,JK)
     ZTHC(JI,JKM) = ZTHC(JI,JKM) + .5 * ZTHC(JI,JK)
     ZRVC(JI,JKP) = ZRVC(JI,JKP) + .3 * ZRVC(JI,JK)
     ZRCC(JI,JKP) = ZRCC(JI,JKP) + .3 * ZRCC(JI,JK)
     ZRIC(JI,JKP) = ZRIC(JI,JKP) + .3 * ZRIC(JI,JK)
     ZTHC(JI,JKP) = ZTHC(JI,JKP) + .3 * ZTHC(JI,JK)
     ZRVC(JI,JK)  = .2 * ZRVC(JI,JK)
     ZRCC(JI,JK)  = .2 * ZRCC(JI,JK)
     ZRIC(JI,JK)  = .2 * ZRIC(JI,JK)
     ZTHC(JI,JK)  = .2 * ZTHC(JI,JK)
  END DO
ENDIF
!
!
          ! Compute vertical integrals - Fluxes
!
JKM = IKE
ZWORK2(:) = 0.
ZWORK2B(:) = 0.
DO JK = IKB+1, JKM
  JKP = JK + 1
  DO JI = KIDIA,KFDIA
    IF ( JK <= ICTL(JI) ) THEN
    ZW1 =  ZRVC(JI,JK) + ZRCC(JI,JK) + ZRIC(JI,JK)
    ZWORK2(JI) = ZWORK2(JI) +  ZW1 *          & ! moisture
                                .5 * (PPABST(JI,JK-1) - PPABST(JI,JKP)) / XG
    ZW1 = ( XCPD + XCPV * ZRW(JI,JK) )* ZTHC(JI,JK)   - &
          ( XLVTT + ( XCPV - XCL ) * ( PTT(JI,JK) - XTT ) ) * ZRCC(JI,JK) - &
          ( XLSTT + ( XCPV - XCL ) * ( PTT(JI,JK) - XTT ) ) * ZRIC(JI,JK)
    ZWORK2B(JI) = ZWORK2B(JI) + ZW1 *         & ! energy
                                .5 * (PPABST(JI,JK-1) - PPABST(JI,JKP)) / XG
    END IF
  END DO
END DO
!
          ! Budget error (integral must be zero)
!
DO JI = KIDIA,KFDIA
  IF ( ICTL(JI) > IKB+1 ) THEN
    JKP = ICTL(JI)
    ZW1 = XG / ( PPABST(JI,IKB) - PPABST(JI,JKP) - &
              .5 * (ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
    ZWORK2(JI) =  ZWORK2(JI) * ZW1
    ZWORK2B(JI) = ZWORK2B(JI)* ZW1
  END IF
END DO
!
          ! Apply uniform correction
!
DO JK = JKM, IKB+1, -1
DO JI = KIDIA,KFDIA
  IF ( ICTL(JI) > IKB+1 .AND. JK <= ICTL(JI) ) THEN
    ZRVC(JI,JK) = ZRVC(JI,JK) - ZWORK2(JI)                                ! moisture
    ZTHC(JI,JK) = ZTHC(JI,JK) - ZWORK2B(JI) /  XCPD                       ! enthalpy
  END IF
END DO
END DO
!
!*           8.7    Compute convective tendencies for Tracers
!                   ------------------------------------------
!
IF ( OCH1CONV ) THEN
  DO JK = IKB, IKE
  DO JI = KIDIA,KFDIA
    IF(GTRIG1(JI) .EQV. .TRUE.)THEN
      ZCH1(JI,JK,:) = PCH1(JI,JK,:)
    ENDIF
  END DO
  END DO
  CALL CONVECT_CHEM_TRANSPORT( KLON, KLEV, KCH1, ZCH1, ZCH1C,          &
                               ISDPL, ISPBL, ISLCL, ICTL, ILFS, ILFS,      &
                               ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,      &
                               ZTIMEC, XA25, ZDMF(:,1), ZLMASS, ZWSUB, &
                               IFTSTEPS )
!
!
!*           8.8    Apply conservation correction
!                   -----------------------------
!
          ! Compute vertical integrals
!
  JKM = IKE
  DO JN = 1, KCH1
    IF(JN < NSV_LGBEG .OR. JN>NSV_LGEND-1) THEN ! no correction for xy lagrangian variables
      ZWORK3(:,JN) = 0.
      ZWORK2(:)    = 0.
      DO JK = IKB+1, JKM
        JKP = JK + 1
        DO JI = KIDIA,KFDIA
          ZW1 = .5 * (PPABST(JI,JK-1) - PPABST(JI,JKP))
          ZWORK3(JI,JN) = ZWORK3(JI,JN) + (ZCH1C(JI,JK,JN)-ZCH1(JI,JK,JN)) * ZW1
          ZWORK2(JI)    = ZWORK2(JI)    + ABS(ZCH1C(JI,JK,JN)) * ZW1
        END DO
      END DO
!
! Apply concentration weighted correction
!
      DO JK = JKM, IKB+1, -1
        DO JI = KIDIA,KFDIA
          IF ( ICTL(JI) > IKB+1 .AND. JK <= ICTL(JI) ) THEN
            ZCH1C(JI,JK,JN) = ZCH1C(JI,JK,JN) -   &
                              ZWORK3(JI,JN)*ABS(ZCH1C(JI,JK,JN))/MAX(1.E-30,ZWORK2(JI))
          END IF
        END DO
      END DO
    END IF
  END DO
END IF

DO JK = IKB, IKE
  ZUMF(:,JK)  = ZUMF(:,JK) / XA25 ! Mass flux per unit area
END DO
IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_COMPUTE',1,ZHOOK_HANDLE)
END SUBROUTINE SHALLOW_CONVECTION_COMPUTE
END SUBROUTINE SHALLOW_CONVECTION

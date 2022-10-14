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
END SUBROUTINE SHALLOW_CONVECTION

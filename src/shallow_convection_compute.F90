SUBROUTINE SHALLOW_CONVECTION_COMPUTE(CVP_SHAL, CVPEXT, KLON, KLEV, KIDIA, KFDIA, KICE, OSETTADJ, PTADJS,  &
                                   PPABST, PZZ, PTT, PRVT, PRCT, PRIT,  &
                                   OCH1CONV, KCH1,&
                                   PCH1, IKB, IKE, IFTSTEPS,   &
                                   PRDOCP, PTHT, PSTHV, PSTHES, ISDPL,  &
                                   ISPBL, ISLCL, PSTHLCL, PSTLCL,       &
                                   PSRVLCL, PSWLCL, PSZLCL, PSTHVELCL,  &
                                   GTRIG1, PUMF,   &
                                   PTHC, PRVC, PRCC, PRIC, ICTL, IMINCTL, &
                                   PPCH1TEN)

USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_CST, ONLY : XCI, XCL, XCPD, XCPV, XG, XLSTT, XLVTT, XP00, XTT
USE MODD_NSV, ONLY : NSV_LGBEG,NSV_LGEND
USE MODD_CONVPAR_SHAL, ONLY : CONVPAR_SHAL
USE MODD_CONVPAREXT, ONLY: CONVPAREXT

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
TYPE(CONVPAR_SHAL),              INTENT(IN)  :: CVP_SHAL
TYPE(CONVPAREXT),                INTENT(IN)  :: CVPEXT
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
REAL   , INTENT(IN)                          :: PRDOCP   ! R_d/C_p
REAL, DIMENSION(KLON,KLEV),      INTENT(IN)  :: PTHT, PSTHV, PSTHES  ! grid scale theta, theta_v
INTEGER, DIMENSION(KLON)  ,      INTENT(IN)  :: ISDPL   ! index for parcel departure level
INTEGER, DIMENSION(KLON)  ,      INTENT(IN)  :: ISPBL   ! index for source layer top
INTEGER, DIMENSION(KLON)  ,      INTENT(IN)  :: ISLCL   ! index for lifting condensation level
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: PSTHLCL ! updraft theta at LCL/L
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: PSTLCL  ! updraft temp. at LCL
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: PSRVLCL ! updraft rv at LCL
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: PSWLCL  ! updraft w at LCL
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: PSZLCL  ! LCL height
REAL, DIMENSION(KLON)     ,      INTENT(IN)  :: PSTHVELCL! envir. theta_v at LCL
LOGICAL, DIMENSION(KLON)  ,      INTENT(IN)  :: GTRIG1  ! logical mask for convection
REAL, DIMENSION(KLON,KLEV),      INTENT(OUT) :: PUMF    ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV),      INTENT(OUT) :: PTHC    ! conv. adj. grid scale theta
REAL, DIMENSION(KLON,KLEV),      INTENT(OUT) :: PRVC    ! conv. adj. grid scale r_w
REAL, DIMENSION(KLON,KLEV),      INTENT(OUT) :: PRCC    ! conv. adj. grid scale r_c
REAL, DIMENSION(KLON,KLEV),      INTENT(OUT) :: PRIC    ! conv. adj. grid scale r_i
INTEGER, DIMENSION(KLON),        INTENT(OUT) :: ICTL    ! index for cloud top level
INTEGER, DIMENSION(KLON),        INTENT(OUT) :: IMINCTL ! min between index for cloud top level
                                                        ! and lifting condensation level
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(OUT) :: PPCH1TEN
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
REAL, DIMENSION(KLON)       :: ZTIMEC  ! advective time period
!
REAL, DIMENSION(KLON,KLEV)  :: ZWSUB   ! envir. compensating subsidence (Pa/s)
!
LOGICAL, DIMENSION(KLON)    :: GTRIG2  ! logical mask for convection
REAL, DIMENSION(KLON)       :: ZCPH    ! specific heat C_ph
REAL, DIMENSION(KLON)       :: ZLV, ZLS! latent heat of vaporis., sublim.
!
! Chemical Tracers:
REAL, DIMENSION(KLON,KCH1)     :: ZWORK3  ! conv. adjust. chemical specy 1
REAL, DIMENSION(KLON,KLEV,KCH1):: ZCH1    ! grid scale chemical specy (kg/kg)
REAL, DIMENSION(KLON,KLEV,KCH1):: ZCH1C   ! conv. adjust. chemical specy 1


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
CALL CONVECT_UPDRAFT_SHAL( CVP_SHAL, CVPEXT, KLON, KLEV, KIDIA, KFDIA,                  &
                           KICE, PPABST, ZDPRES, PZZ, ZTHL, PSTHV, PSTHES, ZRW, &
                           PSTHLCL, PSTLCL, PSRVLCL, PSWLCL, PSZLCL, PSTHVELCL,   &
                           CVP_SHAL%XA25 * 1.E-3, GTRIG2, ISLCL, ISDPL, ISPBL,                &
                           PUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZURW,            &
                           ZURC, ZURI, ZCAPE, ICTL, IETL, GTRIG1                    )

ZDMF(:,:) = 0.
ZDER(:,:) = 0.
ZDDR(:,:) = 0.
ILFS(:)   = IKB
DO JK = IKB, IKE
  ZLMASS(:,JK)  = CVP_SHAL%XA25 * ZDPRES(:,JK) / XG  ! mass of model layer
END DO
ZLMASS(:,IKB) = ZLMASS(:,IKB+1)
!
!-------------------------------------------------------------------------------
!
!*           5.     Compute downdraft properties
!                   ----------------------------
!
  ZTIMEC(:) = CVP_SHAL%XCTIME_SHAL
  IF ( OSETTADJ ) ZTIMEC(:) = PTADJS
!
!*           7.     Determine adjusted environmental values assuming
!                   that all available buoyant energy must be removed
!                   within an advective time step ZTIMEC.
!                   ---------------------------------------------------
!
  CALL CONVECT_CLOSURE_SHAL( CVP_SHAL, CVPEXT, KLON, KLEV, KIDIA, KFDIA, &
                             PPABST, ZDPRES, PZZ, ZLMASS,    &
                             ZTHL, PTHT, ZRW, PRCT, PRIT, GTRIG2,    &
                             PTHC, PRVC, PRCC, PRIC, ZWSUB,       &
                             ISLCL, ISDPL, ISPBL, ICTL,              &
                             PUMF, ZUER, ZUDR, ZUTHL, ZURW,       &
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
  DO JI = KIDIA,KFDIA
   PTHC(JI,JK) = ( PTHC(JI,JK) - PTHT(JI,JK) ) / ZTIMEC(JI)             &
     * ( PPABST(JI,JK) / XP00 ) ** PRDOCP ! change theta in temperature
   PRVC(JI,JK) = ( PRVC(JI,JK) - ZRW(JI,JK) + MAX(0., PRCT(JI,JK)) + MAX(0., PRIT(JI,JK)) ) &
                                        / ZTIMEC(JI)

   PRCC(JI,JK) = ( PRCC(JI,JK) - MAX(0., PRCT(JI,JK)) ) / ZTIMEC(JI)
   PRIC(JI,JK) = ( PRIC(JI,JK) - MAX(0., PRIT(JI,JK)) ) / ZTIMEC(JI)
   ENDDO
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
IF (CVP_SHAL%LLSMOOTH) THEN
  DO JI = KIDIA,KFDIA
     JK = ICTL(JI)
     JKM= MAX(2,ICTL(JI)-1)
     JKP= MAX(2,ICTL(JI)-2)
     PRVC(JI,JKM) = PRVC(JI,JKM) + .5 * PRVC(JI,JK)
     PRCC(JI,JKM) = PRCC(JI,JKM) + .5 * PRCC(JI,JK)
     PRIC(JI,JKM) = PRIC(JI,JKM) + .5 * PRIC(JI,JK)
     PTHC(JI,JKM) = PTHC(JI,JKM) + .5 * PTHC(JI,JK)
     PRVC(JI,JKP) = PRVC(JI,JKP) + .3 * PRVC(JI,JK)
     PRCC(JI,JKP) = PRCC(JI,JKP) + .3 * PRCC(JI,JK)
     PRIC(JI,JKP) = PRIC(JI,JKP) + .3 * PRIC(JI,JK)
     PTHC(JI,JKP) = PTHC(JI,JKP) + .3 * PTHC(JI,JK)
     PRVC(JI,JK)  = .2 * PRVC(JI,JK)
     PRCC(JI,JK)  = .2 * PRCC(JI,JK)
     PRIC(JI,JK)  = .2 * PRIC(JI,JK)
     PTHC(JI,JK)  = .2 * PTHC(JI,JK)
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
    ZW1 =  PRVC(JI,JK) + PRCC(JI,JK) + PRIC(JI,JK)
    ZWORK2(JI) = ZWORK2(JI) +  ZW1 *          & ! moisture
                                .5 * (PPABST(JI,JK-1) - PPABST(JI,JKP)) / XG
    ZW1 = ( XCPD + XCPV * ZRW(JI,JK) )* PTHC(JI,JK)   - &
          ( XLVTT + ( XCPV - XCL ) * ( PTT(JI,JK) - XTT ) ) * PRCC(JI,JK) - &
          ( XLSTT + ( XCPV - XCL ) * ( PTT(JI,JK) - XTT ) ) * PRIC(JI,JK)
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
    PRVC(JI,JK) = PRVC(JI,JK) - ZWORK2(JI)                                ! moisture
    PTHC(JI,JK) = PTHC(JI,JK) - ZWORK2B(JI) /  XCPD                       ! enthalpy
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
  CALL CONVECT_CHEM_TRANSPORT( CVPEXT, KLON, KLEV, KIDIA, KFDIA, KCH1, ZCH1, ZCH1C,&
                               ISDPL, ISPBL, ISLCL, ICTL, ILFS, ILFS,      &
                               PUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,      &
                               ZTIMEC, CVP_SHAL%XA25, ZDMF(:,1), ZLMASS, ZWSUB, &
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
          PPCH1TEN(JI,JK,JN) = (ZCH1C(JI,JK,JN)-ZCH1(JI,JK,JN) ) / ZTIMEC(JI)
          IF(.NOT. GTRIG1(JI)) PPCH1TEN(JI,JK,JN) = 0.
        END DO
      END DO
    END IF
  END DO
END IF

DO JK = IKB, IKE
  DO JI=KIDIA, KFDIA
    IF (ICTL(JI) <= IKB+1) THEN
      PUMF(JI,JK) = 0
    ELSE
      PUMF(JI,JK)  = PUMF(JI,JK) / CVP_SHAL%XA25 ! MASS FLUX PER UNIT AREA
    ENDIF
  ENDDO
END DO

DO JI=KIDIA, KFDIA
  IMINCTL(JI) = MIN(ISLCL(JI), ICTL(JI))
ENDDO

DO JK = IKB, IKE
DO JI = KIDIA,KFDIA
  IF(.NOT. GTRIG1(JI))THEN
    PTHC(JI, JK) = 0.
    PRVC(JI, JK) = 0.
    PRCC(JI, JK) = 0.
    PRIC(JI, JK) = 0.
  ENDIF
ENDDO
ENDDO
DO JI = KIDIA,KFDIA
  IF(.NOT. GTRIG1(JI))THEN
    ICTL(JI) = 0.
    IMINCTL(JI) = 0.
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_COMPUTE',1,ZHOOK_HANDLE)
END SUBROUTINE SHALLOW_CONVECTION_COMPUTE

!     ######spl
      SUBROUTINE CONVECT_CLOSURE_THRVLCL( CVPEXT, CST, KLON, KLEV, KIDIA, KFDIA,          &
                                          PPRES, PTH, PRV, PZ, OWORK1,        &
                                         PTHLCL, PRVLCL, PZLCL, PTLCL, PTELCL,&
                                          KLCL, KDPL, KPBL )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #########################################################################
!
!!**** Determine thermodynamic properties at new LCL
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the thermodynamic
!!      properties at the new lifting condensation level LCL
!!   
!!
!!
!!**  METHOD
!!    ------
!!    see CONVECT_TRIGGER_FUNCT
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! Reference pressure
!!          XRD, XRV           ! Gaz  constants for dry air and water vapor
!!          XCPD               ! Cpd (dry air)
!!          XTT                ! triple point temperature
!!          XBETAW, XGAMW      ! constants for vapor saturation pressure
!!
!!      Module MODD_CONVPAR
!!          XA25               ! reference grid area
!!          XZLCL              ! lowest allowed pressure difference between
!!                             ! surface and LCL
!!          XZPBL              ! minimum mixed layer depth to sustain convection
!!          XWTRIG             ! constant in vertical velocity trigger
!!
!!      Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine TRIGGER_FUNCT)
!!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!!      F Bouyssel  08/11/13 Modifications for reproductibility
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : CST_T
USE MODD_CONVPAREXT, ONLY : CONVPAREXT
USE MODD_CST, ONLY: CST_T
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CONVPAREXT),           INTENT(IN) :: CVPEXT
TYPE(CST_T),                INTENT(IN) :: CST
INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
INTEGER,                    INTENT(IN) :: KIDIA ! value of the first point in x
INTEGER,                    INTENT(IN) :: KFDIA ! value of the last point in x
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! theta
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRV   ! vapor mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of grid point (m)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KDPL  ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KPBL  ! " vert. index of source layer top
LOGICAL, DIMENSION(KLON),   INTENT(IN) :: OWORK1! logical mask 
!
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHLCL ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PRVLCL ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PZLCL  ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTLCL  ! temperature at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTELCL ! environm. temp. at LCL (K)
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KLCL   ! contains vert. index of LCL
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK, JKM, JKMIN, JKMAX      ! vertical loop index
INTEGER :: JI                         ! horizontal loop index 
INTEGER :: IKB, IKE              ! horizontal + vertical loop bounds
REAL    :: ZEPS           ! R_d / R_v
REAL    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(KLON) :: ZPLCL    ! pressure at LCL
REAL, DIMENSION(KLON) :: ZTMIX    ! mixed layer temperature
REAL, DIMENSION(KLON) :: ZEVMIX   ! mixed layer water vapor pressure 
REAL, DIMENSION(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
REAL, DIMENSION(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
REAL, DIMENSION(KLON) :: ZDP      ! pressure between LCL and model layer
REAL, DIMENSION(KLON) :: ZWORK1, ZWORK2     ! work arrays
!
!
!-------------------------------------------------------------------------------
!
!*       0.3    Compute array bounds
!               --------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CONVECT_CLOSURE_THRVLCL',0,ZHOOK_HANDLE)
IKB = 1 + CVPEXT%JCVEXB 
IKE = KLEV - CVPEXT%JCVEXT 
!
!
!*       1.     Initialize local variables
!               --------------------------
!
ZEPS      = CST%XRD / CST%XRV
ZCPORD    = CST%XCPD / CST%XRD
ZRDOCP    = CST%XRD / CST%XCPD
!
ZDPTHMIX(KIDIA:KFDIA) = 0.
ZPRESMIX(KIDIA:KFDIA) = 0.
PTHLCL(KIDIA:KFDIA)   = 300.
PTLCL(KIDIA:KFDIA)    = 300.
PTELCL(KIDIA:KFDIA)   = 300.
PRVLCL(KIDIA:KFDIA)   = 0.
PZLCL(KIDIA:KFDIA)    = PZ(KIDIA:KFDIA,IKB)
ZTMIX(KIDIA:KFDIA)    = 230.
ZPLCL(KIDIA:KFDIA)    = 1.E4 
KLCL(KIDIA:KFDIA)     = IKB + 1
!
!
!*       2.     Construct a mixed layer as in TRIGGER_FUNCT
!               -------------------------------------------
!
     JKMAX=IKE
     JKMIN=IKB
     DO JK = IKB + 1, JKMAX
        JKM = JK + 1
        DO JI = KIDIA, KFDIA
        IF ( JK >= KDPL(JI) .AND. JK <= KPBL(JI) ) THEN
!           
            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            PTHLCL(JI)   = PTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            PRVLCL(JI)   = PRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)
!
        END IF
        END DO
     END DO
!
!
DO JI=KIDIA,KFDIA
  IF ( OWORK1(JI) ) THEN
!
        ZPRESMIX(JI) = ZPRESMIX(JI) / ZDPTHMIX(JI)
        PTHLCL(JI)   = PTHLCL(JI)   / ZDPTHMIX(JI)
        PRVLCL(JI)   = PRVLCL(JI)   / ZDPTHMIX(JI)
!
!*       3.1    Use an empirical direct solution ( Bolton formula )
!               to determine temperature and pressure at LCL.
!               NotaJI the adiabatic saturation temperature is not
!                     equal to the dewpoint temperature
!               --------------------------------------------------
!
!
        ZTMIX(JI)  = PTHLCL(JI) * ( ZPRESMIX(JI) / CST%XP00 ) ** ZRDOCP
        ZEVMIX(JI) = PRVLCL(JI) * ZPRESMIX(JI) / ( PRVLCL(JI) + ZEPS )
        ZEVMIX(JI) = MAX( 1.E-8, ZEVMIX(JI) )
        ZWORK1(JI) = ALOG( ZEVMIX(JI) / 613.3 )
              ! dewpoint temperature
        ZWORK1(JI) = ( 4780.8 - 32.19 * ZWORK1(JI) ) / ( 17.502 - ZWORK1(JI) ) 
              ! adiabatic saturation temperature
        PTLCL(JI)  = ZWORK1(JI) - ( .212 + 1.571E-3 * ( ZWORK1(JI) - CST%XTT )      &
                  - 4.36E-4 * ( ZTMIX(JI) - CST%XTT ) ) * ( ZTMIX(JI) - ZWORK1(JI) )
        PTLCL(JI)  = MIN( PTLCL(JI), ZTMIX(JI) )
        ZPLCL(JI)  = CST%XP00 * ( PTLCL(JI) / PTHLCL(JI) ) ** ZCPORD
!
  END IF
ENDDO
!
     ZPLCL(KIDIA:KFDIA) = MIN( 2.E5, MAX( 10., ZPLCL(KIDIA:KFDIA) ) ) ! bound to avoid overflow
!
!
!*       3.2    Correct PTLCL in order to be completely consistent
!               with MNH saturation formula
!               --------------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( CST, KLON, KIDIA, KFDIA, ZPLCL, PTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
     DO JI=KIDIA,KFDIA
       IF( OWORK1(JI) ) THEN
        ZWORK2(JI) = ZWORK1(JI) / PTLCL(JI) * ( CST%XBETAW / PTLCL(JI) - CST%XGAMW ) ! dr_sat/dT
        ZWORK2(JI) = ( ZWORK1(JI) - PRVLCL(JI) ) /                              &
                        ( 1. + ZLV(JI) / ZCPH(JI) * ZWORK2(JI) ) 
        PTLCL(JI)  = PTLCL(JI) - ZLV(JI) / ZCPH(JI) * ZWORK2(JI)
       END IF
     ENDDO
!
!
!*       3.3    If PRVLCL is oversaturated set humidity and temperature
!               to saturation values.
!               -------------------------------------------------------
!
    CALL CONVECT_SATMIXRATIO( CST, KLON, KIDIA, KFDIA, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
    DO JI=KIDIA,KFDIA
    IF( OWORK1(JI) .AND. PRVLCL(JI) > ZWORK1(JI) ) THEN
        ZWORK2(JI) = ZWORK1(JI) / ZTMIX(JI) * ( CST%XBETAW / ZTMIX(JI) - CST%XGAMW ) ! dr_sat/dT
        ZWORK2(JI) = ( ZWORK1(JI) - PRVLCL(JI) ) /                              &
                        ( 1. + ZLV(JI) / ZCPH(JI) * ZWORK2(JI) )
        PTLCL(JI)  = ZTMIX(JI) + ZLV(JI) / ZCPH(JI) * ZWORK2(JI)
        PRVLCL(JI) = PRVLCL(JI) - ZWORK2(JI)
        ZPLCL(JI)  = ZPRESMIX(JI)
        PTHLCL(JI) = PTLCL(JI) * ( CST%XP00 / ZPLCL(JI) ) ** ZRDOCP
      END IF
    ENDDO
!
!
!*        4.1   Determine  vertical loop index at the LCL 
!               -----------------------------------------
!
     DO JK = JKMIN, IKE - 1
        DO JI = KIDIA, KFDIA
        IF ( ZPLCL(JI) <= PPRES(JI,JK) .AND. OWORK1(JI) ) THEN
            KLCL(JI)  = JK + 1
            PZLCL(JI) = PZ(JI,JK+1)
        END IF
        END DO
     END DO
!
!
!*        4.2   Estimate height and environmental temperature at LCL
!               ----------------------------------------------------
!
    DO JI = KIDIA, KFDIA
        JK   = KLCL(JI)
        JKM  = JK - 1
        ZDP(JI)     = ALOG( ZPLCL(JI) / PPRES(JI,JKM) ) /                     &
                      ALOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZWORK1(JI)  = PTH(JI,JK)  * ( PPRES(JI,JK) / CST%XP00 ) ** ZRDOCP
        ZWORK2(JI)  = PTH(JI,JKM) * ( PPRES(JI,JKM) / CST%XP00 ) ** ZRDOCP
        ZWORK1(JI)  = ZWORK2(JI) + ( ZWORK1(JI) - ZWORK2(JI) ) * ZDP(JI) 
           ! we compute the precise value of the LCL
           ! The precise height is between the levels KLCL and KLCL-1.
        ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
    END DO
    DO JI=KIDIA,KFDIA
    IF( OWORK1(JI) ) THEN
       PTELCL(JI) = ZWORK1(JI)
       PZLCL(JI)  = ZWORK2(JI)
    END IF
    ENDDO
!        
!
!
IF (LHOOK) CALL DR_HOOK('CONVECT_CLOSURE_THRVLCL',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECT_CLOSURE_THRVLCL

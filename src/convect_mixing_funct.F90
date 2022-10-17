!     ######spl
      SUBROUTINE CONVECT_MIXING_FUNCT( KLON, KIDIA, KFDIA,  &
                                       PMIXC, KMF, PER, PDR )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #######################################################
!
!!**** Determine the area under the distribution function
!!     KMF = 1 : gaussian  KMF = 2 : triangular distribution function
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the entrainment and
!!      detrainment rate by evaluating the are under the distribution
!!      function. The integration interval is limited by the critical
!!      mixed fraction PMIXC
!!
!!
!!
!!**  METHOD
!!    ------
!!      Use handbook of mathemat. functions by Abramowitz and Stegun, 1968
!!
!!
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine MIXING_FUNCT)
!!      Abramovitz and Stegun (1968), handbook of math. functions
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,               INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,               INTENT(IN) :: KIDIA  ! value of the first point in x
INTEGER,               INTENT(IN) :: KFDIA  ! value of the last point in x
INTEGER,               INTENT(IN) :: KMF    ! switch for dist. function
REAL, DIMENSION(KLON), INTENT(IN) :: PMIXC  ! critical mixed fraction
!
REAL, DIMENSION(KLON), INTENT(OUT):: PER    ! normalized entrainment rate
REAL, DIMENSION(KLON), INTENT(OUT):: PDR    ! normalized detrainment rate
!
!*       0.2   Declarations of local variables :
!
REAL    :: ZSIGMA = 0.166666667                   ! standard deviation
REAL    :: ZFE    = 4.931813949                   ! integral normalization
REAL    :: ZSQRTP = 2.506628,  ZP  = 0.33267      ! constants
REAL    :: ZA1    = 0.4361836, ZA2 =-0.1201676    ! constants
REAL    :: ZA3    = 0.9372980, ZT1 = 0.500498     ! constants
REAL    :: ZE45   = 0.01111                       ! constant
!
REAL, DIMENSION(KLON) :: ZX, ZY, ZW1, ZW2         ! work variables
REAL    :: ZW11
INTEGER :: JI
!
!
!-------------------------------------------------------------------------------
!
!       1.     Use gaussian function for KMF=1
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CONVECT_MIXING_FUNCT',0,ZHOOK_HANDLE)
IF( KMF == 1 ) THEN
    ! ZX(:)  = ( PMIXC(:) - 0.5 ) / ZSIGMA
      ZX(KIDIA:KFDIA)  = 6. * PMIXC(KIDIA:KFDIA) - 3.
      ZW1(KIDIA:KFDIA) = 1. / ( 1.+ ZP * ABS ( ZX(KIDIA:KFDIA) ) )
      ZY(KIDIA:KFDIA)  = EXP( -0.5 * ZX(KIDIA:KFDIA) * ZX(KIDIA:KFDIA) )
      ZW2(KIDIA:KFDIA) = ZA1 * ZW1(KIDIA:KFDIA) + ZA2 * ZW1(KIDIA:KFDIA) * ZW1(KIDIA:KFDIA) +                   &
               ZA3 * ZW1(KIDIA:KFDIA) * ZW1(KIDIA:KFDIA) * ZW1(KIDIA:KFDIA)
      ZW11   = ZA1 * ZT1 + ZA2 * ZT1 * ZT1 + ZA3 * ZT1 * ZT1 * ZT1
ENDIF
!
DO JI=KIDIA, KFDIA
  IF ( KMF == 1 .AND. ZX(JI) >= 0. ) THEN
          PER(JI) = ZSIGMA * ( 0.5 * ( ZSQRTP - ZE45 * ZW11                 &
                   - ZY(JI) * ZW2(JI) ) + ZSIGMA * ( ZE45 - ZY(JI) ) )        &
                   - 0.5 * ZE45 * PMIXC(JI) * PMIXC(JI)
          PDR(JI) = ZSIGMA*( 0.5 * ( ZY(JI) * ZW2(JI) - ZE45 * ZW11   )       &
                   + ZSIGMA * ( ZE45 - ZY(JI) ) )                           &
                   - ZE45 * ( 0.5 + 0.5 * PMIXC(JI) * PMIXC(JI) - PMIXC(JI) )
  END IF
ENDDO
DO JI=KIDIA, KFDIA
IF ( KMF == 1 .AND. ZX(JI) < 0. ) THEN
        PER(JI) = ZSIGMA*( 0.5 * ( ZY(JI) * ZW2(JI) - ZE45 * ZW11   )       &
                 + ZSIGMA * ( ZE45 - ZY(JI) ) )                           &
                 - 0.5 * ZE45 * PMIXC(JI) * PMIXC(JI)
        PDR(JI) = ZSIGMA * ( 0.5 * ( ZSQRTP - ZE45 * ZW11 - ZY(JI)         &
                 * ZW2(JI) ) + ZSIGMA * ( ZE45 - ZY(JI) ) )                &
                 - ZE45 * ( 0.5 + 0.5 * PMIXC(JI) * PMIXC(JI) - PMIXC(JI) )
  END IF
ENDDO

!
      PER(KIDIA:KFDIA) = PER(KIDIA:KFDIA) * ZFE
      PDR(KIDIA:KFDIA) = PDR(KIDIA:KFDIA) * ZFE
!
!
!       2.     Use triangular function KMF=2
!              -------------------------------
!
!     not yet released
!
!
IF (LHOOK) CALL DR_HOOK('CONVECT_MIXING_FUNCT',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECT_MIXING_FUNCT

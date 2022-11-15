!     ######spl
      MODULE MODD_CONVPAR
!     ###################
!
!!****  *MODD_CONVPAR* - Declaration of convection constants
!!
!!    PURPOSE
!!    -------
!      The purpose of this declarative module is to declare  the
!      constants in the deep convection parameterization.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_CONVPAR)
!!
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  15/11/96
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
TYPE CONVPAR_T
REAL :: XA25        ! 25 km x 25 km reference grid area
!
REAL :: XCRAD       ! cloud radius
REAL :: XCDEPTH     ! minimum necessary cloud depth
REAL :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)
!
REAL :: XZLCL       ! maximum allowed allowed height
                    ! difference between departure level and surface
REAL :: XZPBL       ! minimum mixed layer depth to sustain convection
REAL :: XWTRIG      ! constant in vertical velocity trigger
!
!
REAL :: XNHGAM      ! accounts for non-hydrost. pressure
                    ! in buoyancy term of w equation
                    ! = 2 / (1+gamma)
REAL :: XTFRZ1      ! begin of freezing interval
REAL :: XTFRZ2      ! end of freezing interval
!
REAL :: XRHDBC      ! relative humidity below cloud in downdraft
!
REAL :: XRCONV      ! constant in precipitation conversion
REAL :: XSTABT      ! factor to assure stability in  fractional time
                    ! integration, routine CONVECT_CLOSURE
REAL :: XSTABC      ! factor to assure stability in CAPE adjustment,
                    !  routine CONVECT_CLOSURE
REAL :: XUSRDPTH    ! pressure thickness used to compute updraft
                    ! moisture supply rate for downdraft
REAL :: XMELDPTH    ! layer (Pa) through which precipitation melt is
                    ! allowed below  melting level
REAL :: XUVDP       ! constant for pressure perturb in momentum transport
END TYPE CONVPAR_T

!Keep global variables for parts of the code not ported to the type yet
REAL :: XA25        ! 25 km x 25 km reference grid area
!
REAL :: XCRAD       ! cloud radius
REAL :: XCDEPTH     ! minimum necessary cloud depth
REAL :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)
!
REAL :: XZLCL       ! maximum allowed allowed height
                    ! difference between departure level and surface
REAL :: XZPBL       ! minimum mixed layer depth to sustain convection
REAL :: XWTRIG      ! constant in vertical velocity trigger
!
!
REAL :: XNHGAM      ! accounts for non-hydrost. pressure
                    ! in buoyancy term of w equation
                    ! = 2 / (1+gamma)
REAL :: XTFRZ1      ! begin of freezing interval
REAL :: XTFRZ2      ! end of freezing interval
!
REAL :: XRHDBC      ! relative humidity below cloud in downdraft
!
REAL :: XRCONV      ! constant in precipitation conversion
REAL :: XSTABT      ! factor to assure stability in  fractional time
                    ! integration, routine CONVECT_CLOSURE
REAL :: XSTABC      ! factor to assure stability in CAPE adjustment,
                    !  routine CONVECT_CLOSURE
REAL :: XUSRDPTH    ! pressure thickness used to compute updraft
                    ! moisture supply rate for downdraft
REAL :: XMELDPTH    ! layer (Pa) through which precipitation melt is
                    ! allowed below  melting level
REAL :: XUVDP       ! constant for pressure perturb in momentum transport
!
END MODULE MODD_CONVPAR

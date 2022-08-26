!*   ------------------------------------------------------------------
!
!     FONCTION DE LA PARTITION EAU/GLACE
!     ICE/LIQUID WATER PARTITION FUNCTION
!
REAL(KIND=JPRB) :: PTARG1,P_RDTFAC, FONICE
!
FONICE ( PTARG1, P_RDTFAC ) = 1.0_JPRB - EXP ( - (YDCST%RTT-MIN(YDCST%RTT,PTARG1))**2 &
  & * (1.0_JPRB/(2.0_JPRB*(YDCST%RDT*P_RDTFAC)**2)) )
!*
!     -----------------------------------------------------------------

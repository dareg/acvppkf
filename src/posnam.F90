SUBROUTINE POSNAM(KULNAM,CDNAML)

!**** *POSNAM* - position namelist file for reading

!     Purpose.
!     --------
!     To position namelist file at correct place for reading
!     namelist CDNAML. Replaces use of Cray specific ability
!     to skip to the correct namelist.

!**   Interface.
!     ----------
!        *CALL* *POSNAM*(..)

!        Explicit arguments :     KULNAM - file unit number (input)
!        --------------------     CDNAML - namelist name    (input)

!        Implicit arguments :     None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   None
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 93-06-22
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      R. El Khatib 04-08-10 Apply norms + proper abort if namelist is missing
!      P. Marguinaud   Proxy to POSNAME
!     --------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULNAM 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDNAML 


INTEGER(KIND=JPIM) :: ISTAT

#include "posname.intfb.h"


CALL POSNAME (KULNAM, CDNAML, ISTAT)

SELECT CASE (ISTAT)
  CASE (0)
  CASE (1)
  CASE DEFAULT
END SELECT

END SUBROUTINE POSNAM


SUBROUTINE POSNAME(KULNAM,CDNAML,KSTAT)

!**** *POSNAME* - position namelist file for reading; return error code
!                 if namelist is not found

!     Purpose.
!     --------
!     To position namelist file at correct place for reading
!     namelist CDNAML. Replaces use of Cray specific ability
!     to skip to the correct namelist.

!**   Interface.
!     ----------
!        *CALL* *POSNAME*(..)

!        Explicit arguments :     KULNAM - file unit number (input)
!        --------------------     CDNAML - namelist name    (input)
!                                 KSTAT  - non-zero if namelist not found
!                                          1 = namelist not found

!      P.Marguinaud 22-Nov-2010
!     --------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULNAM 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDNAML 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSTAT



CHARACTER (LEN = 40) ::  CLINE
CHARACTER (LEN =  1) ::  CLTEST

INTEGER(KIND=JPIM) :: ILEN, IND1, ISTATUS, ISCAN

!      -----------------------------------------------------------

!*       1.    POSITION FILE
!              -------------


KSTAT = 0

CLINE='                                        '
REWIND(KULNAM)
ILEN=LEN(CDNAML)
ISTATUS=0
ISCAN=0
DO WHILE (ISTATUS==0 .AND. ISCAN==0)
  READ(KULNAM,'(A)',IOSTAT=ISTATUS) CLINE
  SELECT CASE (ISTATUS)
  CASE (:-1)
    KSTAT=1
    ISCAN=-1
  CASE (0)
    IF (INDEX(CLINE(1:10),'&') == 0) THEN
      ISCAN=0
    ELSE
      IND1=INDEX(CLINE,'&'//CDNAML)
      IF (IND1 == 0) THEN
        ISCAN=0
      ELSE
        CLTEST=CLINE(IND1+ILEN+1:IND1+ILEN+1)
        IF (   (LGE(CLTEST,'0').AND.LLE(CLTEST,'9')) &
         & .OR.(LGE(CLTEST,'A').AND.LLE(CLTEST,'Z')) ) THEN
          ISCAN=0
        ELSE
          ISCAN=1
        ENDIF
      ENDIF
    ENDIF
  CASE (1:)
  END SELECT
ENDDO
BACKSPACE(KULNAM)

!     ------------------------------------------------------------------

END SUBROUTINE POSNAME


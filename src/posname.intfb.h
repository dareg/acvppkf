INTERFACE
SUBROUTINE POSNAME(KULNAM,CDNAML,KSTAT)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KULNAM
CHARACTER(LEN=*) ,INTENT(IN) :: CDNAML
INTEGER(KIND=JPIM),INTENT(OUT) :: KSTAT
END SUBROUTINE POSNAME
END INTERFACE

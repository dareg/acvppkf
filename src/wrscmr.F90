SUBROUTINE WRSCMR(KUNIT,CDNOM,PIN,KLON,KLEN)


!  Ecriture d'une colonne, à l'intérieur du modèle 1D (SCUM/MUSC)

! KLEN peut valoir NFLEV ou NFLEV+1 

!---------------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
!USE YOMLUN_IFSAUX  , ONLY : NULOUT
!---------------------------------------------------------------------------
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KUNIT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIN(KLON,KLEN) 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDNOM
!---------------------------------------------------------------------------

!REAL(KIND=JPRB) :: ZOUT(KLEN)
!INTEGER(KIND=JPIM) :: JLEN
!---------------------------------------------------------------------------

!DO JLEN=1,KLEN
!   ZOUT(JLEN)=PIN(1,JLEN)
!ENDDO
!
!CALL LFAECRR(KUNIT,CDNOM,ZOUT,KLEN)
!---------------------------------------------------------------------------
END SUBROUTINE WRSCMR

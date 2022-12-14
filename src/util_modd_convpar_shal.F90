    MODULE UTIL_MODD_CONVPAR_SHAL
    USE MODD_CONVPAR_SHAL
    INTEGER :: CHECK_REF = 123456789
    CONTAINS

    SUBROUTINE SAVE_MODD_CONVPAR_SHAL(FH)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FH
    WRITE(FH) CHECK_REF

    WRITE(FH) XA25
    WRITE(FH) XCRAD
    WRITE(FH) XCTIME_SHAL
    WRITE(FH) XCDEPTH
    WRITE(FH) XCDEPTH_D
    WRITE(FH) XDTPERT
    WRITE(FH) XATPERT
    WRITE(FH) XBTPERT
    WRITE(FH) XENTR
    WRITE(FH) XZLCL
    WRITE(FH) XZPBL
    WRITE(FH) XWTRIG
    WRITE(FH) XNHGAM
    WRITE(FH) XTFRZ1
    WRITE(FH) XTFRZ2
    WRITE(FH) XSTABT
    WRITE(FH) XSTABC
    WRITE(FH) XAW
    WRITE(FH) XBW
    WRITE(FH) LLSMOOTH
    WRITE(FH) CHECK_REF
    END SUBROUTINE SAVE_MODD_CONVPAR_SHAL

    SUBROUTINE LOAD_MODD_CONVPAR_SHAL(FH)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FH
    INTEGER :: CHECK
    READ(FH) CHECK
    IF(CHECK /= CHECK_REF)WRITE(*,*)__LINE__,__FILE__,"WRONG CONTROL NUMBER",CHECK

    READ(FH) XA25
    READ(FH) XCRAD
    READ(FH) XCTIME_SHAL
    READ(FH) XCDEPTH
    READ(FH) XCDEPTH_D
    READ(FH) XDTPERT
    READ(FH) XATPERT
    READ(FH) XBTPERT
    READ(FH) XENTR
    READ(FH) XZLCL
    READ(FH) XZPBL
    READ(FH) XWTRIG
    READ(FH) XNHGAM
    READ(FH) XTFRZ1
    READ(FH) XTFRZ2
    READ(FH) XSTABT
    READ(FH) XSTABC
    READ(FH) XAW
    READ(FH) XBW
    READ(FH) LLSMOOTH
    READ(FH) CHECK
    IF(CHECK /= CHECK_REF)WRITE(*,*)__LINE__,__FILE__,"WRONG CONTROL NUMBER",CHECK
    END SUBROUTINE LOAD_MODD_CONVPAR_SHAL
END MODULE UTIL_MODD_CONVPAR_SHAL
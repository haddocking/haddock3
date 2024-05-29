C
C ==========================================================================
C
      SUBROUTINE DPSTOR(X)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'machvar.inc'
      DOUBLE PRECISION  X
C
C FORCES THE ARGUMENT VALUE X TO BE STORED IN MEMORY LOCATION V.
C Source code taken from fftpack41/tfftpk.f
C
      DPTMPV=X
      RETURN
      END
C
C ==========================================================================
C
      DOUBLE PRECISION FUNCTION DPTRUNC(X)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'machvar.inc'
      DOUBLE PRECISION  X
C
C DPTRUNC IS A PORTABLE FORTRAN FUNCTION WHICH TRUNCATES A VALUE TO THE
C MACHINE DOUBLE PRECISION WORD SIZE, REGARDLESS OF WHETHER LONGER
C PRECISION INTERNAL REGISTERS ARE USED FOR FLOATING POINT ARITHMETIC IN
C COMPUTING THE VALUE INPUT TO DPTRUNC.  THE METHOD USED IS TO FORCE A
C STORE INTO MEMORY BY USING A COMMON BLOCK IN ANOTHER SUBROUTINE.
C Source code taken from fftpack41/tfftpk.f
C
      CALL DPSTOR(X)
      DPTRUNC=DPTMPV
      RETURN
      END
C
C ==========================================================================
C
      SUBROUTINE SETFPEPS
C
      IMPLICIT NONE
C I/O
      INCLUDE 'numbers.inc'
      INCLUDE 'machvar.inc'
C
C Determine the machine epsilon i.e. the smallest FPEPS such that both
C 1+FPEPS and 1-FPEPS are different from 1.
C
C local
      DOUBLE PRECISION ONEDP
      DOUBLE COMPLEX DBCOMP
C
C external
C NONE
C
C begin
C
C ONEDP doesn't really need to be defined - EPSILON() only needs to know
C the type
      ONEDP=(1.00)
      FPEPS = EPSILON(ONEDP)
      DBCOMP = DCMPLX(ZERO,ZERO)
      CALL DECLAR('FP_EPSILON', 'DP', ' ', DBCOMP, FPEPS)
      RETURN
      END
C
C ==========================================================================
C

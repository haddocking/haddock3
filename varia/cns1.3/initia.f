      SUBROUTINE ATMINI(ISTART,ISTOP)
C
C Subroutine sets coordinate arrays and harmonic constraint array
C to default values for atom numbers between ISTART and ISTOP.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'heap.inc'
      INTEGER ISTART, ISTOP
C local
      INTEGER I
      DOUBLE PRECISION    ANUM
      PARAMETER (ANUM=9999.0D0)
C begin
      DO I=ISTART,ISTOP
      X(I)=ANUM
      Y(I)=ANUM
      Z(I)=ANUM
      WMAIN(I)=0.0D0
      QMAIN(I)=1.0D0
      XCOMP(I)=ANUM
      YCOMP(I)=ANUM
      ZCOMP(I)=ANUM
      WCOMP(I)=0.0D0
      QCOMP(I)=1.0D0
      KCNSTR(I)=0.0D0
      KZCNSTR(I)=0.0D0
      RMSD(I)=0.0D0
      REFX(I)=ANUM
      REFY(I)=ANUM
      REFZ(I)=ANUM
      FBETA(I)=0.0D0
      DX(I)=0.0D0
      DY(I)=0.0D0
      DZ(I)=0.0D0
      XV(I)=0.0D0
      YV(I)=0.0D0
      ZV(I)=0.0D0
      IMOVE(I)=0
      END DO
      CALL ATMIN2(HEAP(IINTER(1)),ISTART,ISTOP)
      RETURN
      END
C
      SUBROUTINE ATMIN2(INTERE,ISTART,ISTOP)
C
      IMPLICIT NONE
C I/O
      INTEGER INTERE(*), ISTART, ISTOP
C local
      INTEGER I
C begin
      DO I=ISTART,ISTOP
      INTERE(I)=0
      ENDDO
C
      RETURN
      END

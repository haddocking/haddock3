C ===============
      SUBROUTINE ECOLLAPSE(EC)
C 
C Front end for ECOLLAPSE2
C
C John Kuszewski 11/12/98
C ===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'collapse.inc'
      INCLUDE 'heap.inc'
C i/o
      DOUBLE PRECISION EC
C local vbls
C begin
      CALL ECOLLAPSE2(EC, HEAP(CLPSSELPTR))
      RETURN
      END
C ===============
      SUBROUTINE ECOLLAPSE2(EC, SELECTEDATOMS)
C
C Energy term for a simple collapse-inducing force
C
C Just forces the radius of gyration to be equal to
C the given value
C
C Perhaps I'll someday change this to handle mass-weighting
C and selections of atoms instead of all atoms
C
C JJK 6/14/96
C ===============
      IMPLICIT NONE
C i/o
      DOUBLE PRECISION EC
      INTEGER SELECTEDATOMS(*)
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'collapse.inc'
      INCLUDE 'deriv.inc'
C local vbls 
      DOUBLE PRECISION CENX, CENY, CENZ, SUM, N, RGYR
      DOUBLE PRECISION DEDXI, DEDYI, DEDZI
      INTEGER COUNT
C begin
C
C first, calculate the centroid of the structure
C
      N = NSELECTED
      CENX = ZERO
      CENY = ZERO
      CENZ = ZERO
      DO COUNT = 1, NSELECTED
           CENX = CENX + X(SELECTEDATOMS(COUNT))
           CENY = CENY + Y(SELECTEDATOMS(COUNT))
           CENZ = CENZ + Z(SELECTEDATOMS(COUNT))
      END DO
      CENX = CENX / N
      CENY = CENY / N
      CENZ = CENZ / N
C
C Now calculate radius of gyration
C
C Rg = RMS (xi - xcm)
C    = Sqrt(Sum(ri - rcm)^2)) / Sqrt(N)
C
C where N is the number of items in the sum
C
      SUM = ZERO
      DO COUNT = 1, NSELECTED
           SUM = SUM + (X(SELECTEDATOMS(COUNT)) - CENX)**2 + 
     &                 (Y(SELECTEDATOMS(COUNT)) - CENY)**2 + 
     &                 (Z(SELECTEDATOMS(COUNT)) - CENZ)**2
      END DO
      RGYR = SQRT(SUM) / SQRT(N)
      EC = KCOLLAPSE * (RGYR - RTARGET)**2
C
C Now calculate forces
C
      DO COUNT = 1, NSELECTED
           DEDXI = 2 * KCOLLAPSE * (RGYR - RTARGET) * 
     &          (X(SELECTEDATOMS(COUNT)) - CENX) /
     &          (SQRT(N) * SQRT(SUM))
           DEDYI = 2 * KCOLLAPSE * (RGYR - RTARGET) * 
     &          (Y(SELECTEDATOMS(COUNT)) - CENY) /
     &          (SQRT(N) * SQRT(SUM))
           DEDZI = 2 * KCOLLAPSE * (RGYR - RTARGET) * 
     &          (Z(SELECTEDATOMS(COUNT)) - CENZ) /
     &          (SQRT(N) * SQRT(SUM))
           DX(SELECTEDATOMS(COUNT)) = DX(SELECTEDATOMS(COUNT)) + DEDXI
           DY(SELECTEDATOMS(COUNT)) = DY(SELECTEDATOMS(COUNT)) + DEDYI
           DZ(SELECTEDATOMS(COUNT)) = DZ(SELECTEDATOMS(COUNT)) + DEDZI
      END DO
      RETURN
      END
C ================
      SUBROUTINE READCOLLAPSE
C
C Reads the input for collapse-energy stuff.
C
C JJK 6/14/96
C ================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'collapse.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
C local vbls
C begin
C
      CALL PUSEND('COLLAPSE>')
      DO WHILE (.NOT.DONE)
           CALL NEXTWD('COLLAPSE>')
           CALL MISCOM('COLLAPSE>',USED)
           IF (.NOT.USED) THEN
C
                IF (WD(1:4).EQ.'HELP') THEN
C-DOCUMENTATION-SOURCE-BEGIN
                WRITE(DUNIT,'(8X,A)')
     &'COLLapse <collapse-statement> END ',
     &'<COLLapse-statement>:== ',
     &'   FORCe <real> | sets the force constant (default 1.0)',
     &'   TARGet <real> | sets the target radius of gyration',
     &'      | (default 1.0)',
     &'   SELEction <sel> | selects the atoms which are to be',
     &'      | used to determine Rgyr and which get the forces',
     &'      | (default: no atoms selected) '
C-DOCUMENTATION-SOURCE-END
C
C set the force constant
C
                ELSE IF (WD(1:4).EQ.'FORC') THEN
                     CALL NEXTF('force constant =', KCOLLAPSE)
C
C set the target radius of gyration
C
                ELSE IF (WD(1:4).EQ.'TARG') THEN
                     CALL NEXTF('target radius of gyration', RTARGET)
C
C get the selection of atoms to include 
C
                ELSE IF (WD(1:4).EQ.'SELE') THEN
                     IF (CLPSSELPTR.EQ.0) THEN
                          CLPSSELPTR=ALLHP(INTEG4(NATOM))
                     END IF
                     CALL READCOLLAPSE2(HEAP(CLPSSELPTR))
C     
C check for END statement
C
                ELSE
                     CALL CHKEND('COUPLINGS>', DONE)
                END IF
           END IF
      END DO
      DONE = .FALSE.
      RETURN
      END
C ==============
      SUBROUTINE READCOLLAPSE2 (SEL)
C
C Reads the selection of atoms to include in the collapse term
C 
C John Kuszewski 11/12/98
C ==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'collapse.inc'
C i/o
      INTEGER SEL(*)
C local vbls
C begin
      CALL SELCTA(SEL, NSELECTED, X, Y, Z,.TRUE.)
      CALL MAKIND(SEL, NATOM, NSELECTED) 
      RETURN
      END
C ==============
      SUBROUTINE COLLAPSEINIT
C
C initializes stuff for the collapse energy term
C
C JJK 6/14/96
C ==============
      IMPLICIT NONE
C include files
      INCLUDE 'collapse.inc'
C begin
      KCOLLAPSE = 1.0D0
      RTARGET = 1.0D0
      NSELECTED = 0
      CLPSSELPTR = 0
      RETURN 
      END

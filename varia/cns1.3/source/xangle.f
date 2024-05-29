C     ------------------------------------------------------------------
      SUBROUTINE EXANGLES (EANGLES, WHICH)
C     ------------------------------------------------------------------
C     Calls EXANGLES2, which does the actual energy calculation
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "heap.inc"
      INCLUDE  "xangle.inc"

      DOUBLE PRECISION EANGLES
      CHARACTER*7 WHICH

      CALL EXANGLES2(EANGLES,HEAP(XANGLESJPTR),HEAP(XANGLESKPTR),
     &                       HEAP(XANGLESJLST),HEAP(XANGLESKLST),
     &                       HEAP(XANGLESOBS1PTR),HEAP(XANGLESOBS2PTR),
     &                       HEAP(XANGLESERRPTR),
     &                       HEAP(CALCXANGLESPTR),WHICH)

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE EXANGLES2 (EANGLES,ATOMJPTR,ATOMKPTR,ATOMJLST,ATOMKLST,
     &                      XANGLESOBS1,XANGLESOBS2,XANGLESERR,
     &                      XANGLESCALC,WHICH)
C     ------------------------------------------------------------------
C     Calculates RDC-ANGLES energies
C
C     Energies are of the form:
C        E = K1*(SIN(ANG)**2)
C     where:
C        K1 = Force constant, 
C        ANG = The angle between the target vector and the actual vector
C
C     Note that: 
C        Atom J = Detected nucleus (e.g. H)
C        Atom K = Coupled nucleus (e.g. N)
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE


      INCLUDE  "cns.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "coord.inc"
      INCLUDE  "numbers.inc"
      INCLUDE  "deriv.inc"
      INCLUDE  "xangle.inc"
      INCLUDE  "consta.inc"

      INTEGER ATOMJPTR(*),ATOMKPTR(*)
      INTEGER ATOMJLST(*),ATOMKLST(*)
      DOUBLE PRECISION XANGLESOBS1(*),XANGLESOBS2(*),XANGLESERR(*)
      DOUBLE PRECISION XANGLESCALC(*),EANGLES
      CHARACTER*7 WHICH
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      INTEGER COUNT,CLASS,MM,NN,II,I,J,K,L,NA,M,N
      PARAMETER (NA=20)
      INTEGER DMM,DNN,OUTFLAG(NA),MCOUNT,NCOUNT,PCOUNT
      DOUBLE PRECISION XJ,XK,YJ,YK,ZJ,ZK,CALCXANGLES,E,
     &                 DFXJ(NA),DFYJ(NA),DFZJ(NA),
     &                 DFXK(NA),DFYK(NA),DFZK(NA),DF1,
     &                 O1,O2,O3,O4,O5,O6,O7,O8,O9,O10,O11,O12,
     &                 PS,DJK,DJK2,DJK3
      DOUBLE PRECISION OBSXANGLES1,OBSXANGLES2,ERRXANGLES,K1,
     &                 A,B,
     &                 DT,DP,DELTA,DELTAXANGLES,
     &                 DD
C     ------------------------------------------------------------------
C     Zero out partial energy
C     ------------------------------------------------------------------
      EANGLES = ZERO

      CLASS = 1
      K1=XANGLESFORCES(1)

      COUNT=0

      DO WHILE(COUNT.LT.ANGLESNUM)
         COUNT = COUNT+1
C     ------------------------------------------------------------------
C     Reset individual E to zero
C     ------------------------------------------------------------------
         E=0

         DO WHILE ((XANGLESASSNDX(CLASS).LT.COUNT).OR.
     &             (XANGLESASSNDX(CLASS).EQ.0))
                    CLASS=CLASS+1
         END DO

         IF (XANGLESASSNDX(CLASS).EQ.XANGLESASSNDX(CLASS-1)) THEN
             COUNT=COUNT-1
         ENDIF

         K1=XANGLESFORCES(CLASS)
C     ------------------------------------------------------------------
C     Note there should only be one atom for J and K
C     ------------------------------------------------------------------
         J=ATOMJPTR(COUNT)+1
         K=ATOMKPTR(COUNT)+1

         XJ=X(ATOMJLST(J))
         XK=X(ATOMKLST(K))
         YJ=Y(ATOMJLST(J))
         YK=Y(ATOMKLST(K))
         ZJ=Z(ATOMJLST(J))
         ZK=Z(ATOMKLST(K))

         OBSXANGLES1=XANGLESOBS1(COUNT)
         OBSXANGLES2=XANGLESOBS2(COUNT)
         ERRXANGLES=XANGLESERR(COUNT)
C     ------------------------------------------------------------------
C     Initialize calculated RDC-angles and counter
C     ------------------------------------------------------------------
         CALCXANGLES=0
         II=0
         MCOUNT=0
         NCOUNT=0
C     ------------------------------------------------------------------
C     Check for correct permutations of paired atoms to get the nucleus-
C     nucleus vector. This depends solely on the order of the assigned
C     atoms. OUTFLAG=1 indicates that the permutation is not allowed.
C     ------------------------------------------------------------------
         DMM=ATOMJPTR(COUNT+1)-ATOMJPTR(COUNT)
         DNN=ATOMKPTR(COUNT+1)-ATOMKPTR(COUNT)
         DO MM=ATOMJPTR(COUNT)+1,ATOMJPTR(COUNT+1)
            MCOUNT=MCOUNT+1
            NCOUNT=0
            DO NN=ATOMKPTR(COUNT)+1,ATOMKPTR(COUNT+1)
               NCOUNT=NCOUNT+1
               II=II+1
               IF ((DMM.GT.1).AND.(DNN.GT.1)) THEN
                    IF (MCOUNT.EQ.NCOUNT) THEN
                        OUTFLAG(II)=0
                    ELSE
                        OUTFLAG(II)=1
                    ENDIF
               ELSE
                    OUTFLAG(II)=0
               END IF
            END DO
         END DO

         II=0
         PCOUNT=0
         DO MM=ATOMJPTR(COUNT)+1,ATOMJPTR(COUNT+1)
            DO NN=ATOMKPTR(COUNT)+1,ATOMKPTR(COUNT+1)
               II=II+1
               IF (OUTFLAG(II).NE.1) THEN
                   PCOUNT=PCOUNT+1
C     ------------------------------------------------------------------
C     Calculate the function and its derivatives
C     ------------------------------------------------------------------
                   OBSXANGLES1=OBSXANGLES1*PI/180.
                   IF (OBSXANGLES1.LT.0.) OBSXANGLES1=OBSXANGLES1+2.*PI
                   OBSXANGLES2=OBSXANGLES2*PI/180.
                   IF (OBSXANGLES2.LT.0.) OBSXANGLES2=OBSXANGLES2+2.*PI

                   O1=SIN(OBSXANGLES1)*COS(OBSXANGLES2)
                   O2=SIN(OBSXANGLES1)*SIN(OBSXANGLES2)
                   O3=COS(OBSXANGLES1)
                   O4=XJ-XK
                   O5=YJ-YK
                   O6=ZJ-ZK
                   O7=O4**2
                   O8=O5**2
                   O9=O6**2
                   DJK2=O7+O8+O9
                   DJK=SQRT(DJK2)
                   DJK3=DJK**3
                   PS=O1*O4+O2*O5+O3*O6
                   O10=O1*DJK-PS*O4/DJK
                   O11=O2*DJK-PS*O5/DJK
                   O12=O3*DJK-PS*O6/DJK

                   DFXJ(II)=-2.*PS*O10/DJK3
                   DFYJ(II)=-2.*PS*O11/DJK3
                   DFZJ(II)=-2.*PS*O12/DJK3
                   DFXK(II)=-DFXJ(II)
                   DFYK(II)=-DFYJ(II)
                   DFZK(II)=-DFZJ(II)
C     ------------------------------------------------------------------
C     Calculate the cosine square of the angle
C     ------------------------------------------------------------------
                   CALCXANGLES=(PS/DJK)**2
               END IF
            END DO
         END DO

         IF (WHICH.EQ.'ANALYZE') THEN
             XANGLESCALC(COUNT)=CALCXANGLES
         END IF
C     ------------------------------------------------------------------
C     The deviation is the sine square of the angle
C     ------------------------------------------------------------------
         DELTA=(1.-CALCXANGLES)
C     ------------------------------------------------------------------
C     Adjust the deviation based on the error range (expressed as the
C     accepted deviation of the cosine square from 1)
C     ------------------------------------------------------------------
         IF ((DELTA.GT.0.000).AND.(DELTA.GT.ERRXANGLES)) THEN
              DELTAXANGLES=DELTA-ERRXANGLES
         ELSE
              DELTAXANGLES=0.0
         END IF

         DD=DELTAXANGLES

         II=0
         PCOUNT=0
         DO MM=ATOMJPTR(COUNT)+1,ATOMJPTR(COUNT+1)
            DO NN=ATOMKPTR(COUNT)+1,ATOMKPTR(COUNT+1)
               II=II+1
               IF (OUTFLAG(II).NE.1) THEN
                   PCOUNT=PCOUNT+1
                   IF (DD.EQ.0.0) THEN
                       E=E+0.0
                       DF1=0.0
                   ELSE
                       E=E+K1*DELTAXANGLES
                       DF1=K1
                   END IF
               END IF
            END DO
         END DO
C     ------------------------------------------------------------------
C     Accumulate energy
C     ------------------------------------------------------------------
         EANGLES=EANGLES+E
C     ------------------------------------------------------------------
C     Now update forces if in energy/force mode
C     ------------------------------------------------------------------
         IF (WHICH.NE.'ANALYZE') THEN
             II=0
             DO MM=ATOMJPTR(COUNT)+1,ATOMJPTR(COUNT+1)
                DO NN=ATOMKPTR(COUNT)+1,ATOMKPTR(COUNT+1)
                   II=II+1
                   IF (OUTFLAG(II).NE.1) THEN
                       DX(ATOMJLST(J))=DX(ATOMJLST(J))+DF1*DFXJ(II)
                       DY(ATOMJLST(J))=DY(ATOMJLST(J))+DF1*DFYJ(II)
                       DZ(ATOMJLST(J))=DZ(ATOMJLST(J))+DF1*DFZJ(II)
                       DX(ATOMKLST(K))=DX(ATOMKLST(K))+DF1*DFXK(II)
                       DY(ATOMKLST(K))=DY(ATOMKLST(K))+DF1*DFYK(II)
                       DZ(ATOMKLST(K))=DZ(ATOMKLST(K))+DF1*DFZK(II)
                   END IF
                END DO
             END DO
         END IF
      END DO

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE READXANGLES
C     ------------------------------------------------------------------
C     Reads in RDC-ANGLES information
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "comand.inc"
      INCLUDE  "xangle.inc"
      INCLUDE  "funct.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "heap.inc"
      INCLUDE  "numbers.inc"

      INTEGER COUNT,SPTR,OLDCLASS,OLDMAXXANGLES
      DOUBLE PRECISION K1,CUTOFF
      CHARACTER*4 THENAME

      SPTR=ALLHP(INTEG4(NATOM))
      CALL PUSEND('XANGLE>')
862   CONTINUE
      CALL NEXTWD('XANGLE>')
      CALL MISCOM('XANGLE>',USED)
      IF (.NOT.USED) THEN
C     ------------------------------------------------------------------
C     Documentation
C     ------------------------------------------------------------------
          IF (WD(1:4).EQ.'HELP') THEN
              WRITE(DUNIT,'(10X,A)')
     &              ' XANGLE {<ANGLE-STATEMENT>} END ',
     &              ' <ANGLE-STATEMENT>:== ',
     &              ' ASSIGN <SEL> <SEL> <REAL> <REAL> <REAL>',
     &              ' * Restraint: Detected_nucleus Coupled_nucleus',
     &              ' Theta Phi Err *',
     &              ' CLASsification <NAME>',
     &              ' * Starts a new class *',
     &              ' FORCeconstant <REAL>',
     &              ' * Force constant for the current class *',
     &              ' NREStraints <INTEGER>',
     &              ' * Number of slots for restraints to allocate *',
     &              ' PRINt THREshold <REAL>',
     &              ' * Prints violations larger than the THREshold *',
     &              ' RESEt',
     &              ' * Erases the restraint table, but keeps NRES *'
C     ------------------------------------------------------------------
C     Get class name. Determine if it's an already-defined class.
C     Insert a new class if it's not.
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'CLAS') THEN
              OLDCLASS=CURCLASS
              CALL NEXTA4('Class name =',THENAME)
              MODE=NEW
              DO COUNT=1,NCLASSES
                 IF (XANGLESCLASSNAMES(COUNT).EQ.THENAME) THEN
                     MODE=UPDATE
                     CURCLASS=COUNT
                  END IF
              END DO
              IF (MODE.EQ.NEW) THEN
C     ------------------------------------------------------------------
C     Make sure you can't add more than the maximum number of classes
C     ------------------------------------------------------------------
                  IF (OLDCLASS.EQ.MAXXANGLESCLASSES) THEN
                      CALL DSPERR('XANGLE','Too many classes.')
                      CALL DSPERR('XANGLE',
     &                     'Increase MAXXANGLESCLASSES and recompile.')
                      CALL WRNDIE(-5, 'READXANGLES',
     &                     'Too many RDC-angles classes.')
                   END IF
                   NCLASSES=NCLASSES+1
                   CURCLASS=NCLASSES
                   XANGLESCLASSNAMES(CURCLASS)=THENAME
C     ------------------------------------------------------------------
C     If this isn't the first class, close off the old class
C     ------------------------------------------------------------------
                   IF (NCLASSES.GT.1) THEN
                       XANGLESASSNDX(OLDCLASS)=ANGLESNUM
                   END IF
              END IF
C     ------------------------------------------------------------------
C     Set force constant for current class
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'FORC') THEN
C     ------------------------------------------------------------------
C     Start a default class if there isn't one defined
C     ------------------------------------------------------------------
              IF (CURCLASS.EQ.0) THEN
                  NCLASSES = 1
                  CURCLASS = 1
              END IF
              CALL NEXTF('Force constant =',K1)
              WRITE(DUNIT,'(A,A,A,F8.3)') 
     &              'Setting force constant for class ', 
     &              XANGLESCLASSNAMES(CURCLASS),' to ',K1
                    XANGLESFORCES(CURCLASS)=K1
C     ------------------------------------------------------------------
C     Reset RDC-angles database
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'RESE') THEN
               CALL XANGLESDEFAULTS
               CALL ALLOCXANGLES(0,MAXXANGLES)
C     ------------------------------------------------------------------
C     Change number of assignment slots
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'NRES') THEN
               OLDMAXXANGLES=MAXXANGLES
               CALL NEXTI('Number of slots =',MAXXANGLES)
               CALL ALLOCXANGLES(OLDMAXXANGLES,MAXXANGLES)
               WRITE(DUNIT,'(A,I8,A)')
     &            'XANGLE: Allocating space for',MAXXANGLES,
     &            'number of restraints.'
C     ------------------------------------------------------------------
C     Read in an assignment
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'ASSI') THEN
C     ------------------------------------------------------------------
C     Make sure you can't add more assignments than you have slots for
C     ------------------------------------------------------------------
               IF (XANGLESMX.EQ.MAXXANGLES) THEN
                   CALL DSPERR('XANGLE','Too many assignments.')
                   CALL DSPERR('XANGLE',
     &                         'Increasing NREStraints by 100.')
                   OLDMAXXANGLES=MAXXANGLES
                   MAXXANGLES=MAXXANGLES+100
                   CALL ALLOCXANGLES(OLDMAXXANGLES,MAXXANGLES)
               END IF
C     ------------------------------------------------------------------
C     If there isn't a class specified, start a default class
C     ------------------------------------------------------------------
               IF (CURCLASS.EQ.0) THEN
                   NCLASSES=1
                   CURCLASS=1
               END IF
               CALL READXANGLES2(HEAP(XANGLESJPTR),HEAP(XANGLESKPTR),
     &                           HEAP(XANGLESJLST),HEAP(XANGLESKLST),
     &                           HEAP(XANGLESOBS1PTR),
     &                           HEAP(XANGLESOBS2PTR),
     &                           HEAP(XANGLESERRPTR),HEAP(SPTR))
C     ------------------------------------------------------------------
C     Print violations
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'PRIN') THEN
               CALL NEXTWD('PRINT>')
               IF (WD(1:4).NE.'THRE') THEN
                   CALL DSPERR('XANGLE',
     &                  'PRINt expects THREshold parameter.')
               ELSE
                   CALL NEXTF('THREshold =',CUTOFF)
                   IF (CUTOFF.LT.ZERO) THEN
                       CALL DSPERR('XANGLE',
     &                             'Cutoff must be positive.')
                       CUTOFF = ABS(CUTOFF)
                   END IF
                   CALL NEXTA4('ALL OR CLASs>',THENAME)
                   IF (THENAME(1:3).EQ.'ALL') THEN
                       DO COUNT=1,NCLASSES
                          PRINTCLASS(COUNT)=.TRUE.
                       END DO
                   ELSE IF (THENAME(1:4).EQ.'CLAS') THEN
                       CALL NEXTA4('CLASS NAME =',THENAME)
                       DO COUNT=1,NCLASSES
                          IF (XANGLESCLASSNAMES(COUNT).EQ.THENAME) THEN
                              PRINTCLASS(COUNT)=.TRUE.
                          ELSE
                              PRINTCLASS(COUNT)=.FALSE.
                          END IF
                       END DO
                    ELSE
                       DO COUNT=1,NCLASSES
                          PRINTCLASS(COUNT)=.TRUE.
                       END DO
                    END IF
                    CALL PRINTXANGLES(CUTOFF,HEAP(CALCXANGLESPTR),
     &                                HEAP(XANGLESOBS1PTR),
     &                                HEAP(XANGLESOBS2PTR),
     &                                HEAP(XANGLESERRPTR),
     &                                HEAP(XANGLESJPTR),
     &                                HEAP(XANGLESKPTR),
     &                                HEAP(XANGLESJLST),
     &                                HEAP(XANGLESKLST))
               END IF
C     ------------------------------------------------------------------
C     Check for END statement
C     ------------------------------------------------------------------
          ELSE
              CALL CHKEND('XANGLE>', DONE)
          END IF
      END IF
      IF (.NOT.(DONE)) GOTO 862
      DONE=.FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE ALLOCXANGLES (OLDSIZE, NEWSIZE)
C     ------------------------------------------------------------------
C     Resets RDC-angles arrays to hold size entries
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "funct.inc"
      INCLUDE  "xangle.inc"

      INTEGER OLDSIZE,NEWSIZE

      IF (OLDSIZE.NE.0) THEN
           CALL FREHP(XANGLESJPTR,INTEG4(OLDSIZE))
           CALL FREHP(XANGLESKPTR,INTEG4(OLDSIZE))

           CALL FREHP(XANGLESJLST,INTEG4(OLDSIZE))
           CALL FREHP(XANGLESKLST,INTEG4(OLDSIZE))

           CALL FREHP(XANGLESOBS1PTR,IREAL8(OLDSIZE))
           CALL FREHP(XANGLESOBS2PTR,IREAL8(OLDSIZE))
           CALL FREHP(XANGLESERRPTR,IREAL8(OLDSIZE))
           CALL FREHP(CALCXANGLESPTR,IREAL8(OLDSIZE))
      END IF

      XANGLESJPTR=ALLHP(INTEG4(NEWSIZE))
      XANGLESKPTR=ALLHP(INTEG4(NEWSIZE))

      XANGLESJLST=ALLHP(INTEG4(NEWSIZE))
      XANGLESKLST=ALLHP(INTEG4(NEWSIZE))

      XANGLESOBS1PTR=ALLHP(IREAL8(NEWSIZE))
      XANGLESOBS2PTR=ALLHP(IREAL8(NEWSIZE))
      XANGLESERRPTR=ALLHP(IREAL8(NEWSIZE))
      CALCXANGLESPTR=ALLHP(IREAL8(NEWSIZE))

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE XANGLESDEFAULTS
C     ------------------------------------------------------------------
C     Sets up defaults
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "xangle.inc"

      INTEGER COUNT

      MODE=NEW
      MAXXANGLES=200
      XANGLESMX=200
      ANGLESNUM=0
      NCLASSES=0
      CURCLASS=0
      DO COUNT=1,MAXXANGLESCLASSES
           XANGLESCLASSNAMES(COUNT)='DEFAULT'
           XANGLESASSNDX(COUNT)=0
           XANGLESFORCES(COUNT)=5.0
      END DO

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE READXANGLES2(ATOMJPTR,ATOMKPTR,ATOMJLST,ATOMKLST,
     &                        XANGLESOBS1,XANGLESOBS2,XANGLESERR,SEL)
C     ------------------------------------------------------------------
C     Reads actual RDC-angles assignments into arrays
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "coord.inc"
      INCLUDE  "xangle.inc"

      INTEGER ATOMJPTR(*),ATOMKPTR(*),ATOMJLST(*),ATOMKLST(*),SEL(*)
      INTEGER JTMP,KTMP,LTMP,II
      DOUBLE PRECISION XANGLESOBS1(*),XANGLESOBS2(*),XANGLESERR(*)
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      INTEGER NSEL,INSERTPOS,COUNT,CURSTOP,OTHERSTOP,NFLAG
      INTEGER I
      DOUBLE PRECISION XANGLESO1,XANGLESO2,XANGLESE
C     ------------------------------------------------------------------
C     If we're in UPDATE mode, make a space for the new line
C     ------------------------------------------------------------------
      NFLAG = 0
      IF (MODE.EQ.UPDATE) THEN
          DO COUNT=ANGLESNUM+1,XANGLESASSNDX(CURCLASS)+1,-1
             ATOMJPTR(COUNT)=ATOMJPTR(COUNT-1)
             ATOMKPTR(COUNT)=ATOMKPTR(COUNT-1)
             XANGLESOBS1(COUNT)=XANGLESOBS1(COUNT-1)
             XANGLESOBS2(COUNT)=XANGLESOBS2(COUNT-1)
             XANGLESERR(COUNT)=XANGLESERR(COUNT-1)
           END DO
           CURSTOP=XANGLESASSNDX(CURCLASS)
           DO COUNT=1,NCLASSES
                OTHERSTOP=XANGLESASSNDX(COUNT)
                IF (OTHERSTOP.GT.CURSTOP) THEN
                     XANGLESASSNDX(COUNT)=OTHERSTOP+1
                END IF
           END DO
           XANGLESASSNDX(CURCLASS)=CURSTOP+1
           INSERTPOS=CURSTOP
           ANGLESNUM=ANGLESNUM + 1
      ELSE
           ANGLESNUM=ANGLESNUM + 1
           INSERTPOS=ANGLESNUM
           XANGLESASSNDX(CURCLASS)=INSERTPOS
      END IF
C     ------------------------------------------------------------------
C     Reading in the atom selection in the constraint table
C     ------------------------------------------------------------------
      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN 
          CALL DSPERR('XANGLE',
     &                'More than 1 atom in SEL for atom J. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XANGLE','Error with atom selection')
          NFLAG=1
      END IF
      CALL MAKIND(SEL, NATOM, NSEL) 
      IF (INSERTPOS.EQ.1) ATOMJPTR(INSERTPOS)=0
      II=ATOMJPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XANGLESMX) XANGLESMX=II
      ATOMJLST(II)=SEL(1)
      ATOMJPTR(INSERTPOS+1)=II

      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN 
          CALL DSPERR('XANGLE',
     &                'More than 1 atom in SEL for atom K. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XANGLE','Error with atom selection')
          NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL) 
      IF (INSERTPOS.EQ.1) ATOMKPTR(INSERTPOS)=0
      II=ATOMKPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XANGLESMX) XANGLESMX=II
      ATOMKLST(II)=SEL(1)
      ATOMKPTR(INSERTPOS+1)=II
C     ------------------------------------------------------------------
C     Reading in the observed Theta and Phi
C     ------------------------------------------------------------------
      CALL NEXTF('Observed THETA =',XANGLESO1)
      CALL NEXTF('Observed PHI =',XANGLESO2)
      CALL NEXTF('Error on sine square =',XANGLESE)

      XANGLESOBS1(INSERTPOS)=XANGLESO1
      XANGLESOBS2(INSERTPOS)=XANGLESO2
      XANGLESERR(INSERTPOS)=XANGLESE
C     ------------------------------------------------------------------
C     Check for error in atom selection. If there is one then reset the
C     counter for restraint
C     ------------------------------------------------------------------
      IF (NFLAG.EQ.1) THEN
          ANGLESNUM=ANGLESNUM-1
      END IF

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE XANGLESINIT
C     ------------------------------------------------------------------
C     Initializes RDC-angles stuff
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "xangle.inc"

      CALL XANGLESDEFAULTS
      CALL ALLOCXANGLES(0,MAXXANGLES)

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE PRINTXANGLES (CUTOFF,XANGLESCALC,XANGLESOBS1,
     &                         XANGLESOBS2,XANGLESERR,
     &                         ATOMJPTR,ATOMKPTR,ATOMJLST,ATOMKLST)
C     ------------------------------------------------------------------
C     Prints RDC-angles with deviation of cosine square from 1 larger
C     than cutoff, calculates RMSD and puts it into $RESULT
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "xangle.inc"
      INCLUDE  "comand.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "numbers.inc"

      DOUBLE PRECISION CUTOFF,XANGLESCALC(*),XANGLESOBS1(*),
     &                        XANGLESOBS2(*),XANGLESERR(*)
      INTEGER ATOMJLST(*),ATOMKLST(*)
      INTEGER ATOMJPTR(*),ATOMKPTR(*)
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      DOUBLE PRECISION CALCXANGLES,OBSXANGLES1,OBSXANGLES2,
     &                 DELTAXANGLES,DELTA,DP
      INTEGER COUNT,CLASS,J,K,NUM,II
      DOUBLE PRECISION RMS,VIOLS,ERRXANGLES,XANGLESENERGY
      DOUBLE COMPLEX DUMMY2
      LOGICAL PRINTTHISCLASS

      RMS=ZERO
      VIOLS=ZERO
      NUM=0
C     ------------------------------------------------------------------
C     Make sure that the array of calculated RDC-angles is up to date
C     ------------------------------------------------------------------
      CALL EXANGLES(XANGLESENERGY,'ANALYZE')
      WRITE (PUNIT,'(A)') 'The following RDC-angles have'
      WRITE (PUNIT,'(A)') 'deviation of cosine square from 1'
      WRITE (PUNIT,'(A)') 'larger than or equal to the cutoff:'
C     ------------------------------------------------------------------
C     Write out first class heading
C     ------------------------------------------------------------------
      CLASS=1
      PRINTTHISCLASS=PRINTCLASS(CLASS)
      IF (PRINTTHISCLASS) THEN
           WRITE (PUNIT,'(A,A)') 'Class ',XANGLESCLASSNAMES(1)
      END IF
C     ------------------------------------------------------------------
C     For every RDC-angles entry...
C     ------------------------------------------------------------------
      COUNT=0
      DO WHILE(COUNT.LT.ANGLESNUM)
          COUNT=COUNT+1
C     ------------------------------------------------------------------
C     Is this the start of a new class?
C     ------------------------------------------------------------------
          IF (XANGLESASSNDX(CLASS).LT.COUNT) THEN
              CLASS=CLASS+1
              PRINTTHISCLASS=PRINTCLASS(CLASS)
              IF (XANGLESASSNDX(CLASS).EQ.XANGLESASSNDX(CLASS-1)) THEN
                  COUNT=COUNT-1
              END IF
              IF (PRINTTHISCLASS) THEN
                  WRITE(PUNIT,'(A,A)') 'Class ',XANGLESCLASSNAMES(CLASS)
              END IF
          END IF
C     ------------------------------------------------------------------
C     If this assignment is in a class to be printed
C     and make sure there is an entry for that class
C     ------------------------------------------------------------------
          IF((PRINTTHISCLASS).AND.
     &       (XANGLESASSNDX(CLASS).NE.XANGLESASSNDX(CLASS-1))) THEN
C     ------------------------------------------------------------------
C     Always update RMSD
C     ------------------------------------------------------------------
              CALCXANGLES=XANGLESCALC(COUNT)
              OBSXANGLES1=XANGLESOBS1(COUNT)
              OBSXANGLES2=XANGLESOBS2(COUNT)
              ERRXANGLES=XANGLESERR(COUNT)
              DP=(1.-CALCXANGLES)

              IF ((DP.GT.0.000).AND.(DP.GT.ERRXANGLES)) THEN
                   DELTAXANGLES=DP-ERRXANGLES
              ELSE
                   DELTAXANGLES=0.0
              END IF

              RMS=RMS+DELTAXANGLES**2
              NUM=NUM+1
C     ------------------------------------------------------------------
C     Print out deviations larger than cutoff
C     and update number of violations
C     ------------------------------------------------------------------
              IF (ABS(DELTAXANGLES).GE.CUTOFF) THEN
                  J=ATOMJLST(ATOMJPTR(COUNT)+1)
                  K=ATOMKLST(ATOMKPTR(COUNT)+1)
                  WRITE(PUNIT,'(A,A)') '==============================',
     &                                 '=============================='
                  WRITE(PUNIT,'(A)') ' Set-J-atoms'
                  DO II=ATOMJPTR(COUNT)+1,ATOMJPTR(COUNT+1)
                     J=ATOMJLST(II)
                     WRITE(PUNIT,'(9X,4(1X,A))') SEGID(J),RESID(J),
     &                                           RES(J),TYPE(J)
                  END DO
                  WRITE(PUNIT,'(A)') ' Set-K-atoms'
                  DO II=ATOMKPTR(COUNT)+1,ATOMKPTR(COUNT+1)
                     K=ATOMKLST(II)
                     WRITE(PUNIT,'(9X,4(1X,A))') SEGID(K),RESID(K),
     &                                           RES(K),TYPE(K)
                  END DO
                  WRITE(PUNIT,'(2(2X,A,1X,F8.3))')
     &                  'Calc: ',CALCXANGLES,'Obs: ',1.
                  WRITE(PUNIT,'(2X,A,1X,F8.3,2X,A,1X,F8.3)')
     &                  'Error: ',ERRXANGLES,'Delta: ',DELTAXANGLES
                  VIOLS=VIOLS+ONE
              END IF
          END IF
      END DO

      IF (NUM.GT.0) THEN
          RMS=SQRT(RMS/NUM)
      ELSE
          RMS=0.0
      ENDIF

      CALL DECLAR('RESULT','DP',' ',DUMMY2,RMS)
      CALL DECLAR('VIOLATIONS','DP',' ',DUMMY2,VIOLS)

      RETURN
      END

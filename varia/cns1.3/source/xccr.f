C     ------------------------------------------------------------------
      SUBROUTINE EXCCR (ECCR, WHICH)
C     ------------------------------------------------------------------
C     Calls EXCCR2, which does the actual energy calculation
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "xccr.inc"
      INCLUDE  "heap.inc"

      DOUBLE PRECISION ECCR
      CHARACTER*7 WHICH

      CALL EXCCR2(ECCR,HEAP(XCCRIPTR),HEAP(XCCRJPTR),HEAP(XCCRKPTR),
     &                 HEAP(XCCRILST),HEAP(XCCRJLST),HEAP(XCCRKLST),
     &                 HEAP(XCCROBSPTR),HEAP(XCCRERRPTR),
     &                 HEAP(CALCXCCRPTR),WHICH)

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE EXCCR2 (ECCR,ATOMIPTR,ATOMJPTR,ATOMKPTR,
     &                        ATOMILST,ATOMJLST,ATOMKLST,
     &                        XCCROBS,XCCRERR,
     &                        XCCRCALC,WHICH)
C     ------------------------------------------------------------------
C     Calculates cross correlation rate energies
C
C     Energies are of the form:
C        E = K*(DELTACCR**2)
C     where: 
C        K = FORCE CONSTANT 
C        DELTACCR = CALCULATED CCR - OBSERVED CCR
C     and cross correlation rate function is defined as:
C        CCR = KM*(3*COS(THETA)^2-1)/DIST^3
C
C     Note that:
C        Atom I = Metal 
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
      INCLUDE  "xccr.inc"
      INCLUDE  "consta.inc"

      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*)
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*)
      DOUBLE PRECISION XCCROBS(*),XCCRERR(*)
      DOUBLE PRECISION XCCRCALC(*),ECCR
      CHARACTER*7 WHICH
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      INTEGER COUNT,CLASS,MM,NN,II,I,J,K,L,NA,M,N
      PARAMETER (NA=20)
      INTEGER DMM,DNN,OUTFLAG(NA),MCOUNT,NCOUNT,PCOUNT
      DOUBLE PRECISION XI,XJ,XK,YI,YJ,YK,ZI,ZJ,ZK,CALCXCCR,E,
     &                 DFXI(NA),DFYI(NA),DFZI(NA),DFXJ(NA),DFYJ(NA),
     &                 DFZJ(NA),DFXK(NA),DFYK(NA),DFZK(NA),
     &                 DJK2,DJI2,DJI5,DJI,DF1,PS,DJK,DJI3,
     &                 O1,O2,O3,O4,O5,O6,O7,O8,O9,O10,O11,O12,O13,O14,
     &                 O15,O16,O17,O18,O19,O20,O21,O22,O23,O24,O25,O26,
     &                 O27,O28,O29,O30,O31,O32,O33,O34,O35,O36,O37,O38,
     *                 O39
      DOUBLE PRECISION OBSXCCR,ERRXCCR,K1,COEF1,A,B,DELTA,
     &                 DP,DELTAXCCR,DD
C     ------------------------------------------------------------------
C     Zero out partial energy
C     ------------------------------------------------------------------
      ECCR = ZERO
  
      CLASS = 1
      K1=XCCRFORCES(1)
      COEF1=XCCRCOEF1(1)

      COUNT=0
      DO WHILE(COUNT.LT.CCRNUM)
         COUNT=COUNT+1
C     ------------------------------------------------------------------
C     Reset individual E to zero
C     ------------------------------------------------------------------
         E=0

         DO WHILE ((XCCRASSNDX(CLASS).LT.COUNT).OR.
     &             (XCCRASSNDX(CLASS).EQ.0))
            CLASS=CLASS+1
         END DO

         IF (XCCRASSNDX(CLASS).EQ.XCCRASSNDX(CLASS-1)) THEN
             COUNT=COUNT-1
         ENDIF

         K1=XCCRFORCES(CLASS)
         COEF1=XCCRCOEF1(CLASS)
C     ------------------------------------------------------------------
C     Note there should only be one atom for I, J, and K
C     ------------------------------------------------------------------
         I=ATOMIPTR(COUNT)+1
         J=ATOMJPTR(COUNT)+1
         K=ATOMKPTR(COUNT)+1

         XI=X(ATOMILST(I))
         XJ=X(ATOMJLST(J))
         XK=X(ATOMKLST(K))
         YI=Y(ATOMILST(I))
         YJ=Y(ATOMJLST(J))
         YK=Y(ATOMKLST(K))
         ZI=Z(ATOMILST(I))
         ZJ=Z(ATOMJLST(J))
         ZK=Z(ATOMKLST(K))

         OBSXCCR=XCCROBS(COUNT)
         ERRXCCR=XCCRERR(COUNT)
C     ------------------------------------------------------------------
C     Initialize calculated cross correlation rates and counter
C     ------------------------------------------------------------------
         CALCXCCR=0
         II=0
         MCOUNT=0
         NCOUNT=0
C     ------------------------------------------------------------------
C     Check for correct permutations of paired atoms to get the metal-
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
                   O1=XJ-XK
                   O2=XJ-XI
                   O3=YJ-YK
                   O4=YJ-YI
                   O5=ZJ-ZK
                   O6=ZJ-ZI
                   O7=O2**2
                   O8=O4**2
                   O9=O6**2
                   O10=O1**2
                   O11=O3**2
                   O12=O5**2
                   PS=O1*O2+O3*O4+O5*O6
                   DJK2=O10+O11+O12
                   DJI2=O7+O8+O9
                   DJI5=SQRT(DJI2**5)
                   DJI3=SQRT(DJI2**3)
                   DJI=SQRT(DJI2)
                   O13=2.*XJ-XI-XK
                   O14=6.*PS*O13
                   O15=2.*O1*DJI2
                   O16=2.*O2*DJK2
                   O17=DJK2*DJI5
                   O18=3.*PS**2-DJK2*DJI2
                   O19=2.*O1*DJI5
                   O20=5.*O2*DJK2*DJI3
                   O21=O17**2
                   O22=6.*PS*(-O2)
                   O23=6.*PS*(-O1)
                   O24=2.*YJ-YI-YK
                   O25=6.*PS*O24
                   O26=2.*O3*DJI2
                   O27=2.*O4*DJK2
                   O28=2.*O3*DJI5
                   O29=5.*O4*DJK2*DJI3
                   O30=6.*PS*(-O4)
                   O31=6.*PS*(-O3)
                   O32=2.*ZJ-ZI-ZK
                   O33=6.*PS*O32
                   O34=2.*O5*DJI2
                   O35=2.*O6*DJK2
                   O36=2.*O5*DJI5
                   O37=5.*O6*DJK2*DJI3
                   O38=6.*PS*(-O6)
                   O39=6.*PS*(-O5)
                   DJK=SQRT(DJK2)
 
                   DFXI(II)=((O23+O16)*O17-O18*O20)/O21
                   DFYI(II)=((O31+O27)*O17-O18*O29)/O21
                   DFZI(II)=((O39+O35)*O17-O18*O37)/O21
                   DFXJ(II)=((O14-O15-O16)*O17-O18*(O19+O20))/O21
                   DFYJ(II)=((O25-O26-O27)*O17-O18*(O28+O29))/O21
                   DFZJ(II)=((O33-O34-O35)*O17-O18*(O36+O37))/O21
                   DFXK(II)=((O22+O15)*O17+O18*O19)/O21
                   DFYK(II)=((O30+O26)*O17+O18*O28)/O21
                   DFZK(II)=((O38+O34)*O17+O18*O36)/O21
C     ------------------------------------------------------------------
C     Coefficient
C     ------------------------------------------------------------------
                   A=COEF1
C     ------------------------------------------------------------------
C     Calculate the cross correlation rate
C     ------------------------------------------------------------------
                   CALCXCCR=A*((3.*(PS/(DJK*DJI))**2-1.)/DJI3)
               END IF
            END DO
         END DO

         IF (WHICH.EQ.'ANALYZE') THEN
             XCCRCALC(COUNT)=CALCXCCR
         END IF

         DELTA=(CALCXCCR-OBSXCCR)
C     ------------------------------------------------------------------
C     Adjust the deviation based on the error range
C     ------------------------------------------------------------------
         IF ((DELTA.LT.0.000).AND.(ABS(DELTA).GT.ERRXCCR)) THEN
              DELTAXCCR=DELTA+ERRXCCR
         ELSE IF ((DELTA.GT.0.000).AND.(DELTA.GT.ERRXCCR)) THEN
              DELTAXCCR=DELTA-ERRXCCR
         ELSE
              DELTAXCCR=0.0
         END IF

         DD=DELTAXCCR

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
C     ------------------------------------------------------------------
C     Depending on CCRFLAG, the weight may be proportional to DIST^3
C     ------------------------------------------------------------------
                       IF (CCRFLG) THEN
                           E=E+K1*(DELTAXCCR**2)*DJI3
                           DF1=2*K1*DELTAXCCR*DJI3
                       ELSE
                          E=E+K1*(DELTAXCCR**2)
                          DF1=2*K1*DELTAXCCR
                       END IF
                   END IF
               END IF
            END DO
         END DO
C     ------------------------------------------------------------------
C     Accumulate energy
C     ------------------------------------------------------------------
         ECCR=ECCR+E
C     ------------------------------------------------------------------
C     Now update forces if in energy/force mode
C     ------------------------------------------------------------------
         IF (WHICH.NE.'ANALYZE') THEN
             II=0
             DO MM=ATOMJPTR(COUNT)+1,ATOMJPTR(COUNT+1)
                DO NN=ATOMKPTR(COUNT)+1,ATOMKPTR(COUNT+1)
                   II=II+1
                   IF (OUTFLAG(II).NE.1) THEN
                       DX(ATOMILST(I))=DX(ATOMILST(I))+DF1*DFXI(II)
                       DY(ATOMILST(I))=DY(ATOMILST(I))+DF1*DFYI(II)
                       DZ(ATOMILST(I))=DZ(ATOMILST(I))+DF1*DFZI(II)
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
      SUBROUTINE READXCCR
C     ------------------------------------------------------------------
C     Reads in cross correlation rate information
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "comand.inc"
      INCLUDE  "xccr.inc"
      INCLUDE  "funct.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "heap.inc"
      INCLUDE  "numbers.inc"

      INTEGER COUNT,SPTR,OLDCLASS,OLDMAXXCCR,FNCLASS
      DOUBLE PRECISION K1,CUTOFF,COEF1
      CHARACTER*4 THENAME
      CHARACTER*132 NOMEFILE
      INTEGER WEISWI
      
      NOMEFILE='XCCR.OUT'

      SPTR=ALLHP(INTEG4(NATOM))
      CALL PUSEND('XCCR>')
862   CONTINUE
      CALL NEXTWD('XCCR>')
      CALL MISCOM('XCCR>',USED)
      IF (.NOT.USED) THEN
C     ------------------------------------------------------------------
C     Documentation
C     ------------------------------------------------------------------
          IF (WD(1:4).EQ.'HELP') THEN
              WRITE(DUNIT,'(10X,A)')
     &              ' XCCR {<CCR-STATEMENT>} END ',
     &              ' <CCR-STATEMENT>:== ',
     &              ' ASSIGN <SEL> <SEL> <SEL> <REAL> <REAL>',
     &              ' * Restraint: Metal Coupled_nucleus',
     &              ' Detected_nucleus CCR Err *',
     &              ' CLASsification <NAME>',
     &              ' * Starts a new class *',
     &              ' WEIP <1|0>',
     &              ' * Switch weight proportional to R^3 (1=Y,0=N) *',
     &              ' COEFficient <REAL>',
     &              ' * Proportionality constant *',
     &              ' FORCeconstant <REAL>',
     &              ' * Force constant for the current class *',
     &              ' NREStraints <INTEGER>',
     &              ' * Number of slots for restraints to allocate *',
     &              ' PRINt THREshold <REAL>',
     &              ' * Prints violations larger than the THREshold *',
     &              ' RESEt',
     &              ' * Erases the restraint table, but keeps NRES *',
     &              ' FRUN',
     &              ' * Runs FANTACROSS *'
C     ------------------------------------------------------------------
C     About FANTACROSS
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'FRUN') THEN
              CALL NEXTI('FANTACROSS on class number =',FNCLASS)
              IF (FNCLASS.GT.NCLASSES.OR.FNCLASS.EQ.0) THEN
                  PRINT*,'%FRUN-ERR: This class does not exist...'
              ELSE
                  CALL FANTACCR(HEAP(XCCRIPTR),HEAP(XCCRJPTR),
     &                          HEAP(XCCRKPTR),HEAP(XCCRILST),
     &                          HEAP(XCCRJLST),HEAP(XCCRKLST),
     &                          HEAP(XCCROBSPTR),HEAP(XCCRERRPTR),
     &                          FNCLASS,NOMEFILE)
              END IF
C     ------------------------------------------------------------------
C     Get class name. Determine if it's an already-defined class.
C     Insert a new class if it's not.
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'CLAS') THEN
              OLDCLASS=CURCLASS
              CALL NEXTA4('Class name =',THENAME)
              MODE=NEW
              DO COUNT=1,NCLASSES
                 IF (XCCRCLASSNAMES(COUNT).EQ.THENAME) THEN
                     MODE=UPDATE
                     CURCLASS=COUNT
                 END IF
              END DO
              IF (MODE.EQ.NEW) THEN
C     ------------------------------------------------------------------
C     Make sure you can't add more than the maximum number of classes
C     ------------------------------------------------------------------
                  IF (OLDCLASS.EQ.MAXXCCRCLASSES) THEN
                      CALL DSPERR('XCCR','Too many classes.')
                      CALL DSPERR('XCCR',
     &                     'Increase MAXXCCRCLASSES and recompile.')
                      CALL WRNDIE(-5, 'READXCCR',
     &                      'Too many CCR classes.')
                  END IF
                  NCLASSES=NCLASSES+1
                  CURCLASS=NCLASSES
                  XCCRCLASSNAMES(CURCLASS)=THENAME
C     ------------------------------------------------------------------
C     If this isn't the first class, close off the old class
C     ------------------------------------------------------------------
                  IF (NCLASSES.GT.1) THEN
                      XCCRASSNDX(OLDCLASS)=CCRNUM
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
     &              'Setting force const for class ', 
     &              XCCRCLASSNAMES(CURCLASS),' to ',K1
                    XCCRFORCES(CURCLASS)=K1
C     ------------------------------------------------------------------
C     Set weight proportional to DIST^3
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'WEIP') THEN
              CALL NEXTI('Weight proportional to R^3 (Y=1,N=0):',WEISWI)
              IF (WEISWI.EQ.1) THEN
                  CCRFLG=.TRUE.
                  WRITE (DUNIT,'(A)') 'CCR weight proportional to R^3.'
              ELSE IF (WEISWI.EQ.0) THEN
                  CCRFLG=.FALSE.
                  WRITE (DUNIT,'(A)') 'CCR weight not dependent on R.'
              ELSE
                  WRITE (DUNIT,'(A)') 'Unknown switch. Default weight.'
              END IF
C     ------------------------------------------------------------------
C     Set coefficient constant for current class
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'COEF') THEN
              CALL NEXTF('CCR COEFficient =',COEF1)
C     ------------------------------------------------------------------
C     Start a default class if there isn't one already defined
C     ------------------------------------------------------------------
              IF (CURCLASS.EQ.0) THEN
                  NCLASSES = 1
                  CURCLASS = 1
              END IF
              WRITE(DUNIT,'(A,A,A,F8.3)')
     &              'Setting coefficient for class ',
     &               XCCRCLASSNAMES(CURCLASS),'to',COEF1
                     XCCRCOEF1(CURCLASS)=COEF1
C     ------------------------------------------------------------------
C     Reset cross correlation rates database
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'RESE') THEN
              CALL XCCRDEFAULTS
              CALL ALLOCXCCR(0,MAXXCCR)
C     ------------------------------------------------------------------
C     Change number of assignment slots
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'NRES') THEN
               OLDMAXXCCR=MAXXCCR
               CALL NEXTI('Number of slots =',MAXXCCR)
               CALL ALLOCXCCR(OLDMAXXCCR,MAXXCCR)
               WRITE(DUNIT,'(A,I8,A)')
     &            'XCCR: Allocating space for',MAXXCCR,
     &            'number of restraints.'
C     ------------------------------------------------------------------
C     Read in an assignment
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'ASSI') THEN
C     ------------------------------------------------------------------
C     Make sure you can't add more assignments than you have slots for
C     ------------------------------------------------------------------
              IF (XCCRMX.EQ.MAXXCCR) THEN
                  CALL DSPERR('XCCR','Too many assignments.')
                  CALL DSPERR('XCCR',
     &                        'Increasing NREStraints by 100.')
                  OLDMAXXCCR=MAXXCCR
                  MAXXCCR=MAXXCCR+100
                  CALL ALLOCXCCR(OLDMAXXCCR,MAXXCCR)
              END IF
C     ------------------------------------------------------------------
C     If there isn't a class specified, start a default class
C     ------------------------------------------------------------------
              IF (CURCLASS.EQ.0) THEN
                  NCLASSES=1
                  CURCLASS=1
              END IF
              CALL READXCCR2(HEAP(XCCRIPTR),HEAP(XCCRJPTR), 
     &                       HEAP(XCCRKPTR),HEAP(XCCRILST),
     &                       HEAP(XCCRJLST),HEAP(XCCRKLST),
     &                       HEAP(XCCROBSPTR),HEAP(XCCRERRPTR),
     &                       HEAP(SPTR))
C     ------------------------------------------------------------------
C     Print violations
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'PRIN') THEN
              CALL NEXTWD('PRINT>')
              IF (WD(1:4).NE.'THRE') THEN
                  CALL DSPERR('XCCR',
     &                 'PRINt expects THREshold parameter.')
              ELSE
                  CALL NEXTF('THREshold =',CUTOFF)
                  IF (CUTOFF.LT.ZERO) THEN
                      CALL DSPERR('XCCR',
     &                            'Cutoff must be positive.')
                      CUTOFF = ABS(CUTOFF)
                  END IF
                  CALL NEXTA4('ALL OR CLASs>',THENAME)
                  IF (THENAME(1:3).EQ.'ALL') THEN
                      DO COUNT=1,NCLASSES
                         PRINTCLASS(COUNT)=.TRUE.
                      END DO
                  ELSE IF (THENAME(1:4).EQ.'CLAS') THEN
                      CALL NEXTA4('Class name =',THENAME)
                      DO COUNT=1,NCLASSES
                         IF (XCCRCLASSNAMES(COUNT).EQ.THENAME) THEN
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
                  CALL PRINTXCCR(CUTOFF,HEAP(CALCXCCRPTR),
     &                           HEAP(XCCROBSPTR),HEAP(XCCRERRPTR),
     &                           HEAP(XCCRIPTR),HEAP(XCCRJPTR),
     &                           HEAP(XCCRKPTR),HEAP(XCCRILST),
     &                           HEAP(XCCRJLST),HEAP(XCCRKLST))
              END IF
C     ------------------------------------------------------------------
C     Check for END statement
C     ------------------------------------------------------------------
          ELSE
              CALL CHKEND('XCCR>',DONE)
          END IF
      END IF
      IF (.NOT.(DONE)) GOTO 862
      DONE=.FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE ALLOCXCCR (OLDSIZE,NEWSIZE)
C     ------------------------------------------------------------------
C     Resets cross correlation rate arrays to hold size entries
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE "funct.inc"
      INCLUDE "xccr.inc"

      INTEGER OLDSIZE,NEWSIZE

      IF (OLDSIZE.NE.0) THEN
           CALL FREHP(XCCRIPTR,INTEG4(OLDSIZE))
           CALL FREHP(XCCRJPTR,INTEG4(OLDSIZE))
           CALL FREHP(XCCRKPTR,INTEG4(OLDSIZE))

           CALL FREHP(XCCRILST,INTEG4(OLDSIZE))
           CALL FREHP(XCCRJLST,INTEG4(OLDSIZE))
           CALL FREHP(XCCRKLST,INTEG4(OLDSIZE))

           CALL FREHP(XCCROBSPTR,IREAL8(OLDSIZE))
           CALL FREHP(XCCRERRPTR,IREAL8(OLDSIZE))
           CALL FREHP(CALCXCCRPTR,IREAL8(OLDSIZE))
      END IF
      XCCRIPTR=ALLHP(INTEG4(NEWSIZE))
      XCCRJPTR=ALLHP(INTEG4(NEWSIZE))
      XCCRKPTR=ALLHP(INTEG4(NEWSIZE))

      XCCRILST=ALLHP(INTEG4(NEWSIZE))
      XCCRJLST=ALLHP(INTEG4(NEWSIZE))
      XCCRKLST=ALLHP(INTEG4(NEWSIZE))

      XCCROBSPTR=ALLHP(IREAL8(NEWSIZE))
      XCCRERRPTR=ALLHP(IREAL8(NEWSIZE))
      CALCXCCRPTR=ALLHP(IREAL8(NEWSIZE))

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE XCCRDEFAULTS
C     ------------------------------------------------------------------
C     Sets up defaults
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE "xccr.inc"

      INTEGER COUNT

      MODE=NEW
      MAXXCCR=200
      XCCRMX=200
      CCRNUM=0
      NCLASSES=0
      CURCLASS=0
      CCRFLG=.FALSE.
      DO COUNT=1,MAXXCCRCLASSES
         XCCRCLASSNAMES(COUNT)='DEFAULT'
         XCCRASSNDX(COUNT)=0
         XCCRFORCES(COUNT)=5.0
         XCCRCOEF1(COUNT)=1.0
      END DO

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE READXCCR2(ATOMIPTR,ATOMJPTR,ATOMKPTR,
     &                     ATOMILST,ATOMJLST,ATOMKLST,
     &                     XCCROBS,XCCRERR,SEL)
C     ------------------------------------------------------------------
C     Reads actual cross correlation rate assignments into arrays
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE "cns.inc"
      INCLUDE "mtf.inc"
      INCLUDE "coord.inc"
      INCLUDE "xccr.inc"

      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*)
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*)
      INTEGER SEL(*)
      INTEGER ITMP,JTMP,KTMP,II
      DOUBLE PRECISION XCCROBS(*),XCCRERR(*)
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      INTEGER NSEL,INSERTPOS,COUNT,CURSTOP,OTHERSTOP,NFLAG
      INTEGER I
      DOUBLE PRECISION XCCRO,XCCRE
C     ------------------------------------------------------------------
C     If we're in UPDATE mode, make a space for the new line
C     ------------------------------------------------------------------
      NFLAG = 0
      IF (MODE.EQ.UPDATE) THEN
          DO COUNT=CCRNUM+1,XCCRASSNDX(CURCLASS)+1,-1
                ATOMIPTR(COUNT)=ATOMIPTR(COUNT-1)
                ATOMJPTR(COUNT)=ATOMJPTR(COUNT-1)
                ATOMKPTR(COUNT)=ATOMKPTR(COUNT-1)
                XCCROBS(COUNT)=XCCROBS(COUNT-1)
                XCCRERR(COUNT)=XCCRERR(COUNT-1)
          END DO
          CURSTOP=XCCRASSNDX(CURCLASS)
          DO COUNT=1,NCLASSES
             OTHERSTOP=XCCRASSNDX(COUNT)
             IF (OTHERSTOP.GT.CURSTOP) THEN
                 XCCRASSNDX(COUNT)=OTHERSTOP+1
             END IF
          END DO
          XCCRASSNDX(CURCLASS)=CURSTOP+1
          INSERTPOS=CURSTOP
          CCRNUM=CCRNUM+1
      ELSE
           CCRNUM=CCRNUM+1
           INSERTPOS=CCRNUM
           XCCRASSNDX(CURCLASS)=INSERTPOS
      END IF
C     ------------------------------------------------------------------
C     Reading in the atom selection in the restraint table
C     ------------------------------------------------------------------
      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN 
          CALL DSPERR('XCCR',
     &                'More than 1 atom in SEL for atom I. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XCCR','Error with atom selection')
          NFLAG = 1
      END IF
      CALL MAKIND(SEL, NATOM, NSEL) 
      IF (INSERTPOS.EQ.1) ATOMIPTR(INSERTPOS)=0
      II=ATOMIPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XCCRMX) XCCRMX=II
      ATOMILST(II) = SEL(1)
      ATOMIPTR(INSERTPOS+1)=II

      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN 
          CALL DSPERR('XCCR',
     &                'More than 1 atom in SEL for atom J. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XCCR','Error with atom selection')
          NFLAG = 1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL) 
      IF (INSERTPOS.EQ.1) ATOMJPTR(INSERTPOS)=0
      II=ATOMJPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XCCRMX) XCCRMX=II
      ATOMJLST(II)=SEL(1)
      ATOMJPTR(INSERTPOS+1)=II

      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN 
          CALL DSPERR('XCCR',
     &                'More than 1 atom in SEL for atom K. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XCCR','Error with atom selection')
          NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL) 
      IF (INSERTPOS.EQ.1) ATOMKPTR(INSERTPOS)=0
      II=ATOMKPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XCCRMX) XCCRMX=II
      ATOMKLST(II)=SEL(1)
      ATOMKPTR(INSERTPOS+1)=II
C     ------------------------------------------------------------------
C     Reading in the observed cross correlation rate
C     ------------------------------------------------------------------
      CALL NEXTF('Observed CCR =',XCCRO)
      CALL NEXTF('Error in CCR =',XCCRE)

      XCCROBS(INSERTPOS)=XCCRO
      XCCRERR(INSERTPOS)=XCCRE
C     ------------------------------------------------------------------
C     Check for error in atom selection. If there is one then reset the
C     counter for restraint
C     ------------------------------------------------------------------
      IF (NFLAG.EQ.1) THEN
          CCRNUM=CCRNUM-1
      END IF

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE XCCRINIT
C     ------------------------------------------------------------------
C     Initializes cross correlation rate stuff
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE "xccr.inc"

      CALL XCCRDEFAULTS
      CALL ALLOCXCCR(0,MAXXCCR)

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE PRINTXCCR (CUTOFF,XCCRCALC,XCCROBS,XCCRERR,
     &                      ATOMIPTR,ATOMJPTR,ATOMKPTR,
     &                      ATOMILST,ATOMJLST,ATOMKLST)
C     ------------------------------------------------------------------
C     Prints cross correlation rates with dveiation larger than cutoff,
C     calculates RMSD and puts it into $RESULT
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "xccr.inc"
      INCLUDE  "comand.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "numbers.inc"

      DOUBLE PRECISION CUTOFF,XCCRCALC(*),XCCROBS(*),XCCRERR(*)
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*)
      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*)
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      DOUBLE PRECISION CALCXCCR,OBSXCCR,DELTAXCCR,DELTA,DP
      INTEGER COUNT,CLASS,I,J,K,NUM,II
      DOUBLE PRECISION RMS,VIOLS,ERRXCCR,XCCRENERGY
      DOUBLE COMPLEX DUMMY2
      LOGICAL PRINTTHISCLASS

      RMS=ZERO
      VIOLS=ZERO
      NUM=0
C     ------------------------------------------------------------------
C     Make sure that the array of calculated CCR is up to date
C     ------------------------------------------------------------------
      CALL EXCCR(XCCRENERGY,'ANALYZE')
      WRITE (PUNIT,'(A)') 'The following cross correlation rates have'
      WRITE (PUNIT,'(A)') 'deviation larger than or'
      WRITE (PUNIT,'(A)') 'equal to the cutoff:'
C     ------------------------------------------------------------------
C     Write out first class heading
C     ------------------------------------------------------------------
      CLASS=1
      PRINTTHISCLASS=PRINTCLASS(CLASS)
      IF (PRINTTHISCLASS) THEN
           WRITE (PUNIT,'(A,A)') 'Class ',XCCRCLASSNAMES(1)
      END IF
C     ------------------------------------------------------------------
C     For every cross correlation rate entry...
C     ------------------------------------------------------------------
      COUNT = 0
      DO WHILE(COUNT.LT.CCRNUM)
         COUNT=COUNT+1
C     ------------------------------------------------------------------
C     Is this the start of a new class?
C     ------------------------------------------------------------------
         IF (XCCRASSNDX(CLASS).LT.COUNT) THEN
             CLASS=CLASS+1
             PRINTTHISCLASS=PRINTCLASS(CLASS)
             IF (XCCRASSNDX(CLASS).EQ.XCCRASSNDX(CLASS-1)) THEN
                 COUNT=COUNT - 1
             END IF
             IF (PRINTTHISCLASS) THEN
                  WRITE (PUNIT,'(A,A)') 'Class ',XCCRCLASSNAMES(CLASS)
             END IF
         END IF
C     ------------------------------------------------------------------
C     If this assignment is in a class to be printed
C     and make sure there is an entry for that class
C     ------------------------------------------------------------------
         IF ((PRINTTHISCLASS).AND.
     &       (XCCRASSNDX(CLASS).NE.XCCRASSNDX(CLASS-1))) THEN
C     ------------------------------------------------------------------
C     Always update RMSD
C     ------------------------------------------------------------------
              CALCXCCR=XCCRCALC(COUNT)
              OBSXCCR=XCCROBS(COUNT)
              ERRXCCR=XCCRERR(COUNT)
              DP=(CALCXCCR-OBSXCCR)

              IF ((DP.LT.0.000).AND.(ABS(DP).GT.ERRXCCR)) THEN
                   DELTAXCCR=DP+ERRXCCR
              ELSE IF ((DP.GT.0.000).AND.(DP.GT.ERRXCCR)) THEN
                   DELTAXCCR=DP-ERRXCCR
              ELSE
                   DELTAXCCR=0.0
              END IF

              RMS=RMS+DELTAXCCR**2
              NUM=NUM+1
C     ------------------------------------------------------------------
C     Print out deviations larger than cutoff
C     and update number of violations
C     ------------------------------------------------------------------
              IF (ABS(DELTAXCCR).GE.CUTOFF) THEN
                  I=ATOMILST(ATOMIPTR(COUNT)+1)
                  J=ATOMJLST(ATOMJPTR(COUNT)+1)
                  K=ATOMKLST(ATOMKPTR(COUNT)+1)
                  WRITE(PUNIT,'(A,A)') '==============================',
     &                                 '=============================='
                  WRITE(PUNIT,'(A)') ' Set-I-atoms'
                  DO II=ATOMIPTR(COUNT)+1,ATOMIPTR(COUNT+1)
                     I=ATOMILST(II)
                     WRITE(PUNIT,'(9X,4(1X,A))') SEGID(I),RESID(I),
     &                                           RES(I),TYPE(I)
                  END DO
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
     &                  'Calc: ',CALCXCCR,'Obs: ',OBSXCCR
                  WRITE(PUNIT,'(2X,A,1X,F8.3,2X,A,1X,F8.3)')
     &                  'Error: ',ERRXCCR,'Delta: ',DELTAXCCR
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

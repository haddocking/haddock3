C     ------------------------------------------------------------------
      SUBROUTINE FANTACCR (ATOMIPTR,ATOMJPTR,ATOMKPTR, 
     &                     ATOMILST,ATOMJLST,ATOMKLST, 
     &                     XCCROBS,XCCRERR,NCLASS,NOMEFILE)
C     ------------------------------------------------------------------
C     Fits cross correlation rate coefficient on a given structure
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
      INCLUDE  "fantaxplor.inc"
      INCLUDE  'heap.inc'

      DOUBLE PRECISION AK
      LOGICAL CSWITCH,CSAVE,CERR
      COMMON /CRUN/AK 
      COMMON /CPARAM/CSWITCH,CSAVE,CERR
      INTEGER CVIOLS,IA
      COMMON /CVIOL/CVIOLS
      CHARACTER NOMEFILE*132
      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*) 
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*)

      INTEGER COUNT,MM,NN,I,J,K,L,NA,M,N,CLASS,NORES
      INTEGER NCLASS,NNCO,NNCC
      DOUBLE PRECISION XCCROBS(*),XCCRERR(*)
      DOUBLE PRECISION OBSXCCR,ERRXCCR
      DOUBLE PRECISION XI,XJ,XK,YI,YJ,YK,ZI,ZJ,ZK
      DOUBLE PRECISION XFA(1000),YFA(1000),ZFA(1000),XMFA(1000),
     &                 YMFA(1000),ZMFA(1000),ERRFA(1000),OBSFA(1000),
     &                 FENPCS(20,2000),WPROT(1000),
     &                 XFAN(1000),YFAN(1000),ZFAN(1000)
      DOUBLE PRECISION O1,O2,O3,O4,O5,O6,O7,O8,O9,O10,O11,O12
      DOUBLE PRECISION PS,DJK2,DJI2,DJK(1000),DJI(1000),ANGFACT(1000)
      DOUBLE COMPLEX DUMMY2
C     ------------------------------------------------------------------
C     CLASS must reamin this during the cycle of variables exchange.
C     The class is defined by NCLASS.
C     ------------------------------------------------------------------
      CLASS=1
      NNCO=0
      NNCC=0
      COUNT=0

      PRINT*,'This is FANTACROSS.'
      PRINT*,'CCRNUM is:',CCRNUM

      DO WHILE(COUNT.LT.CCRNUM)
         COUNT=COUNT+1
         DO WHILE ((XCCRASSNDX(CLASS).LT.COUNT).OR.
     &             (XCCRASSNDX(CLASS).EQ.0))
            CLASS=CLASS+1
         END DO
         IF (XCCRASSNDX(CLASS).EQ.XCCRASSNDX(CLASS-1)) THEN
             COUNT=COUNT-1
         ENDIF

         IF (NCLASS.EQ.CLASS) THEN 
             NNCO=NNCO+1
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
             OBSFA(NNCO)=XCCROBS(COUNT)
             OBSXCCR=XCCROBS(COUNT)
             ERRFA(NNCO)=XCCRERR(COUNT)
             ERRXCCR=XCCRERR(COUNT)
         END IF

         DO MM=ATOMJPTR(COUNT)+1,ATOMJPTR(COUNT+1)
            DO NN=ATOMKPTR(COUNT)+1,ATOMKPTR(COUNT+1)
               IF (NCLASS.EQ.CLASS) THEN
                   NNCC=NNCC+1
                   XMFA(NNCC)=X(ATOMILST(MM))
                   XFA(NNCC)=X(ATOMJLST(NN))
                   XFAN(NNCC)=X(ATOMKLST(NN))
                   YMFA(NNCC)=Y(ATOMILST(MM))
                   YFA(NNCC)=Y(ATOMJLST(NN))
                   YFAN(NNCC)=Y(ATOMKLST(NN))
                   ZMFA(NNCC)=Z(ATOMILST(MM))
                   ZFA(NNCC)=Z(ATOMJLST(NN))
                   ZFAN(NNCC)=Z(ATOMKLST(NN))
                   NORES=ATOMJLST(NN)
                   NRESIDFA(NNCC)=RESID(NORES)
                   NRESFA(NNCC)=RES(NORES)
                   NTYPEFA(NNCC)=TYPE(NORES)
               END IF
            END DO
         END DO
      END DO      

      IF (NNCC.EQ.NNCO) THEN
          DO IA=1,CCRNUM
             O1=XFA(IA)-XFAN(IA)
             O2=XFA(IA)-XMFA(IA)
             O3=YFA(IA)-YFAN(IA)
             O4=YFA(IA)-YMFA(IA)
             O5=ZFA(IA)-ZFAN(IA)
             O6=ZFA(IA)-ZMFA(IA)
             O7=O2**2
             O8=O4**2
             O9=O6**2
             O10=O1**2
             O11=O3**2
             O12=O5**2
             PS=O1*O2+O3*O4+O5*O6
             DJK2=O10+O11+O12
             DJK(IA)=SQRT(DJK2)
             DJI2=O7+O8+O9
             DJI(IA)=SQRT(DJI2)
             ANGFACT(IA)=(3.*(PS/(DJK(IA)*DJI(IA)))**2-1.)
           END DO 
           CALL FANTALINCCR(DJI,DJK,ANGFACT,NNCO,OBSFA,ERRFA,NOMEFILE,
     &                      WPROT,0)
      ELSE
         PRINT*,'Error in atom count...'
      END IF

      RETURN
      END

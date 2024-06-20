C     ------------------------------------------------------------------
      SUBROUTINE FANTADIPO(ATOMIPTR,ATOMJPTR,ATOMKPTR,ATOMLPTR,ATOMMPTR,
     &                     ATOMNPTR,ATOMILST,ATOMJLST,ATOMKLST,ATOMLLST,
     &                     ATOMMLST,ATOMNLST,XRDCOBS,XRDCERR,
     &                     NCLASS,NOMEFILE)
C     ------------------------------------------------------------------
C     Fits residual dipolar couplings tensor on a given structure
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "coord.inc"
      INCLUDE  "numbers.inc"
      INCLUDE  "deriv.inc"
      INCLUDE  "xdipo_rdc.inc"
      INCLUDE  "consta.inc"
      INCLUDE  "fantaxplor.inc"
      INCLUDE  'heap.inc'

      CHARACTER NOMEFILE*132
      DOUBLE PRECISION APARA,APERP
      COMMON /FRUN/APERP,APARA
      INTEGER FDSWITCH,FDSAVE,FDERR
      COMMON /FDPARAM/FDSWITCH,FDSAVE,FDERR
      INTEGER FDVIOLS
      COMMON /DVIOL/FDVIOLS
      DOUBLE PRECISION FENERGTOT
      COMMON /FENERG/FENERGTOT
      DOUBLE PRECISION FDENT(20,2000),FDPER(20,2000),FDPAR(20,2000),
     &                 FDIPOE(20,2000) 
      INTEGER DCONTSAVE(20)
      COMMON /FDMEDIA/FDENT,FDPER,FDPAR,FDIPOE,DCONTSAVE
      DOUBLE PRECISION DIPOE
      COMMON /DIPOENERGY/ DIPOE
      DOUBLE COMPLEX DUMMY2
      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*)
      INTEGER ATOMLPTR(*),ATOMMPTR(*),ATOMNPTR(*)
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*)
      INTEGER ATOMLLST(*),ATOMMLST(*),ATOMNLST(*)
      INTEGER COUNT,MM,NN,I,J,K,L,NA,M,N,CLASS,NORES
      INTEGER NCLASS,NNCO,NNCC
      DOUBLE PRECISION XRDCOBS(*),XRDCERR(*)
      DOUBLE PRECISION XI,XJ,XK,XL,XM,XN,YI,YJ,YK,YL,YM,YN,
     &                 ZI,ZJ,ZK,ZL,ZM,ZN
      DOUBLE PRECISION XFA(1000),YFA(1000),ZFA(1000),XMFA(1000),
     &                 YMFA(1000),ZMFA(1000),ERRFA(1000),OBSFA(1000),
     &                 WPROT(1000)
C     ------------------------------------------------------------------
C     CLASS must remain this during the cycle of variables exchange.
C     The class is defined by NCLASS.
C     ------------------------------------------------------------------
      CLASS=1
      NNCO=0
      NNCC=0
      COUNT=0

      PRINT*,'This is FANTALIN for RDC.'
      PRINT*,'RDCNUM is:',RDCNUM

      DO WHILE(COUNT.LT.RDCNUM)
         COUNT=COUNT+1
         DO WHILE ((XRDCASSNDX(CLASS).LT.COUNT).OR.
     &             (XRDCASSNDX(CLASS).EQ.0))
            CLASS=CLASS+1
         END DO
         IF (XRDCASSNDX(CLASS).EQ.XRDCASSNDX(CLASS-1)) THEN
             COUNT=COUNT-1
         ENDIF

         IF (NCLASS.EQ.CLASS) THEN 
             NNCO=NNCO+1
             I=ATOMIPTR(COUNT)+1
             J=ATOMJPTR(COUNT)+1
             K=ATOMKPTR(COUNT)+1
             L=ATOMLPTR(COUNT)+1
             XI=X(ATOMILST(I))
             XJ=X(ATOMJLST(J))
             XK=X(ATOMKLST(K))
             XL=X(ATOMLLST(L))
             YI=Y(ATOMILST(I))
             YJ=Y(ATOMJLST(J))
             YK=Y(ATOMKLST(K))
             YL=Y(ATOMLLST(L))
             ZI=Z(ATOMILST(I))
             ZJ=Z(ATOMJLST(J))
             ZK=Z(ATOMKLST(K))
             ZL=Z(ATOMLLST(L))
             OBSFA(NNCO)=XRDCOBS(COUNT)
             ERRFA(NNCO)=XRDCERR(COUNT)
         END IF

         DO MM=ATOMMPTR(COUNT)+1,ATOMMPTR(COUNT+1)
            DO NN=ATOMNPTR(COUNT)+1,ATOMNPTR(COUNT+1)
               IF (NCLASS.EQ.CLASS) THEN
                   NNCC=NNCC+1
                   XMFA(NNCC)=X(ATOMMLST(MM))
                   XFA(NNCC)=X(ATOMNLST(NN))
                   YMFA(NNCC)=Y(ATOMMLST(MM))
                   YFA(NNCC)=Y(ATOMNLST(NN))
                   ZMFA(NNCC)=Z(ATOMMLST(MM))
                   ZFA(NNCC)=Z(ATOMNLST(NN))
                   NORES=ATOMNLST(NN)
                   NRESIDFA(NNCC)=RESID(NORES)
                   NRESFA(NNCC)=RES(NORES)
                   NTYPEFA(NNCC)=TYPE(NORES)
               END IF
            END DO
         END DO
      END DO      

      IF (NNCC.EQ.NNCO) THEN
          CALL FANTALINEXE(XFA,YFA,ZFA,XMFA,YMFA,ZMFA,NNCO,OBSFA,ERRFA,
     &                     0,NOMEFILE,WPROT,0,DIPFLG)
          PRINT*,'FDERR is=',FDERR
          IF (FDERR.EQ.1) THEN
              CALL FANTAERRORE(XFA,YFA,ZFA,XMFA,YMFA,ZMFA,NNCO,OBSFA,
     &                         ERRFA,0,DIPFLG)
          END IF
      ELSE
          PRINT*,'Error in atom count...'
      END IF

      PRINT*,'FDSWITCH is:',FDSWITCH
      IF (FDSWITCH.EQ.1) THEN
          XRDCCOEF1(NCLASS)=APARA
          XRDCCOEF2(NCLASS)=APERP
          PRINT*,'For class:',NCLASS
          PRINT*,'New value for A1 =',APARA
          PRINT*,'New value for A2 =',APERP
      END IF

      IF (FDSAVE.EQ.1) THEN
          DCONTSAVE(NCLASS)=DCONTSAVE(NCLASS)+1
          FDENT(NCLASS,DCONTSAVE(NCLASS))=FENERGTOT
          FDIPOE(NCLASS,DCONTSAVE(NCLASS))=DIPOE
          FDPAR(NCLASS,DCONTSAVE(NCLASS))=APARA
          FDPER(NCLASS,DCONTSAVE(NCLASS))=APERP
          PRINT*,'XDIPO_RDC energy =',DIPOE
          PRINT*,'Total energy =',FENERGTOT
          PRINT*,'DCONTSAVE is:',DCONTSAVE(NCLASS)
      END IF

      CALL DECLAR('DCHIAX','DP',' ', DUMMY2,APARA)
      CALL DECLAR('DCHIRH','DP',' ', DUMMY2,APERP)

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE FANTADIPOMEDIA (FMA,PERC)
C     ------------------------------------------------------------------
C     Averages tensor parameters on a given sub-ensemble of structures
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION VETE(2,2000),VETT(3,2000)
      DOUBLE PRECISION FDENT(20,2000),FDPER(20,2000),FDPAR(20,2000),
     &                 FDIPOE(20,2000)
      COMMON /FDMEDIA/ FDENT,FDPER,FDPAR,FDIPOE,DCONTSAVE
      INTEGER DCONTSAVE(20)
      INTEGER FMA,COUNT,CB,NUM,COUNT1,CONFPRE
      DOUBLE PRECISION PERC
      DOUBLE PRECISION M1PAR,M1PER,DEVPAR
      DOUBLE COMPLEX DUMMY2
C     ------------------------------------------------------------------
C     Through FMA the class on which the average must be done is passed.
C     If it is zero, all the counters are reset.
C     ------------------------------------------------------------------
      IF (FMA.EQ.0) THEN
          DO COUNT=1,19
             DCONTSAVE(COUNT)=0
          END DO
          PRINT*,'Mean vectors reset.'
      ELSE
C     ------------------------------------------------------------------
C     The best are selected
C     ------------------------------------------------------------------
          DO COUNT=1,DCONTSAVE(FMA)
             VETE(1,COUNT)=0
          END DO

          DO COUNT=1,DCONTSAVE(FMA)
            NUM=DCONTSAVE(FMA)
            DO WHILE (VETE(1,NUM).EQ.1)
               NUM=NUM-1
            END DO
            DO COUNT1=1,DCONTSAVE(FMA)-1
               IF ((FDENT(FMA,COUNT1).LE.FDENT(FMA,NUM)).AND.
     &             (VETE(1,COUNT1).EQ.0.0)) THEN
                    NUM=COUNT1
               END IF
            END DO
            PRINT*,'NUM is:',NUM
            VETE(1,NUM)=1
            VETE(2,COUNT)=FDENT(FMA,NUM)
            VETT(1,COUNT)=FDPAR(FMA,NUM)
            VETT(2,COUNT)=FDPER(FMA,NUM)
         END DO

         CONFPRE=NINT((DBLE(DCONTSAVE(FMA))*(PERC/100.)))
         IF (CONFPRE.EQ.0) THEN
             CONFPRE=1
         END IF
         PRINT*,'CONFPRE is:',CONFPRE
         PRINT*,'PERC is:',PERC
         PRINT*,'DBLE(CONTSAVE(FMA)) is:',DBLE(DCONTSAVE(FMA))

         DO COUNT=1,DCONTSAVE(FMA)
            PRINT*,'Number',COUNT,'Energy',VETE(1,COUNT),VETE(2,COUNT)
            PRINT*,'A1',VETT(1,COUNT),' A2',VETT(2,COUNT)
         END DO

         CB=0
         M1PAR=0
         M1PER=0
         DO COUNT=1,CONFPRE
            M1PAR=M1PAR+VETT(1,COUNT)
            M1PER=M1PER+VETT(2,COUNT)
            CB=CB+1
            PRINT*,'A1 =',VETT(1,COUNT),'A2 =',VETT(2,COUNT)
         END DO
         M1PAR=M1PAR/CB
         M1PER=M1PER/CB

         CB=0
         DEVPAR=0
         DO COUNT=1,CONFPRE
            DEVPAR=DEVPAR+(VETT(1,COUNT)-M1PAR)**2
         END DO
         DEVPAR=DEVPAR/CONFPRE
         PRINT*,'Averaged new A1 =',M1PAR
         PRINT*,'Averaged new A2 =',M1PER
         PRINT*,'Standard deviation for new A1 =',SQRT(DEVPAR)

         CALL DECLAR('DCHIAX','DP',' ',DUMMY2,M1PAR)
         CALL DECLAR('DCHIRH','DP',' ',DUMMY2,M1PER)

         END IF

         RETURN
         END

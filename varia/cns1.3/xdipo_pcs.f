C     ------------------------------------------------------------------
      SUBROUTINE FANTALIN (ATOMIPTR,ATOMJPTR,ATOMKPTR,ATOMLPTR,ATOMMPTR,
     &                     ATOMNPTR,ATOMILST,ATOMJLST,ATOMKLST,ATOMLLST,
     &                     ATOMMLST,ATOMNLST,XPCSOBS,XPCSERR,NCLASS,
     &                     NOMEFILE)
C     ------------------------------------------------------------------
C     Fits pseudocontact shift tensor on a given structure
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "coord.inc"
      INCLUDE  "numbers.inc"
      INCLUDE  "deriv.inc"
      INCLUDE  "xdipo_pcs.inc"
      INCLUDE  "consta.inc"
      INCLUDE  "fantaxplor.inc"
      INCLUDE  'heap.inc'

      DOUBLE PRECISION APARA,APERP
      COMMON /FRUN/APERP,APARA
      INTEGER FSWITCH,FSAVE,FERR
      COMMON /FPARAM/FSWITCH,FSAVE,FERR
      INTEGER FVIOLS
      COMMON /VIOL/FVIOLS
      DOUBLE PRECISION FENERGTOT
      COMMON /FENERG/ FENERGTOT
      DOUBLE PRECISION XPCSE
      COMMON /XPCSENERGY/ XPCSE
      DOUBLE PRECISION FENT(20,2000),FPER(20,2000),FPAR(20,2000),
     &                 FXPCSE(20,2000)
      INTEGER CONTSAVE(20)
      COMMON /FMEDIA/FENT,FPER,FPAR,FXPCSE,CONTSAVE
      CHARACTER NOMEFILE*132
      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*),ATOMLPTR(*)
      INTEGER ATOMMPTR(*),ATOMNPTR(*)
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*),ATOMLLST(*)
      INTEGER ATOMMLST(*),ATOMNLST(*)
      INTEGER COUNT,MM,NN,I,J,K,L,NA,M,N,CLASS,NORES
      INTEGER NCLASS,NNCO,NNCC
      DOUBLE PRECISION XPCSOBS(*),XPCSERR(*)
      DOUBLE PRECISION OBSXPCS,ERRXPCS
      DOUBLE PRECISION XI,XJ,XK,XL,XM,XN,YI,YJ,YK,YL,YM,YN,
     &                 ZI,ZJ,ZK,ZL,ZM,ZN
      DOUBLE PRECISION XFA(1000),YFA(1000),ZFA(1000),
     &                 XMFA(1000),YMFA(1000),ZMFA(1000),
     &                 ERRFA(1000),OBSFA(1000),
     &                 WPROT(1000)
      DOUBLE COMPLEX DUMMY2
C     ------------------------------------------------------------------
C     CLASS must remain this during the cycle of variables exchange.
C     The class is defined by NCLASS.
C     ------------------------------------------------------------------
      CLASS=1
      NNCO=0
      NNCC=0
      COUNT=0

      PRINT*,'This is FANTALIN for PCS.'
      PRINT*,'PCSNUM is:',PCSNUM

      DO WHILE (COUNT.LT.PCSNUM)
         COUNT=COUNT+1
         DO WHILE ((XPCSASSNDX(CLASS).LT.COUNT).OR.
     &            (XPCSASSNDX(CLASS).EQ.0))
            CLASS=CLASS+1
         END DO
         IF (XPCSASSNDX(CLASS).EQ.XPCSASSNDX(CLASS-1)) THEN
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
             OBSFA(NNCO)=XPCSOBS(COUNT)
             OBSXPCS=XPCSOBS(COUNT)
             ERRFA(NNCO)=XPCSERR(COUNT)
             ERRXPCS=XPCSERR(COUNT)
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
     &                     1,NOMEFILE,WPROT,0,PCSFLG)
          PRINT*,'FERR is:',FERR
          IF (FERR.EQ.1) THEN
              CALL FANTAERRORE(XFA,YFA,ZFA,XMFA,YMFA,ZMFA,NNCO,OBSFA,
     &                         ERRFA,1,PCSFLG)
          END IF
      ELSE
          PRINT*,'Error in atom count...'
      END IF

      PRINT *,'FSWITCH is:',FSWITCH
      IF (FSWITCH.EQ.1) THEN
          XPCSCOEF1(NCLASS)=APARA
          XPCSCOEF2(NCLASS)=APERP
          PRINT*,'For class:',NCLASS
          PRINT*,'New value for A1 =',APARA
          PRINT*,'New value for A2 =',APERP
      END IF

      IF (FSAVE.EQ.1) THEN
          CONTSAVE(NCLASS)=CONTSAVE(NCLASS)+1
          FENT(NCLASS,CONTSAVE(NCLASS))=FENERGTOT
          FXPCSE(NCLASS,CONTSAVE(NCLASS))=XPCSE
          FPAR(NCLASS,CONTSAVE(NCLASS))=APARA
          FPER(NCLASS,CONTSAVE(NCLASS))=APERP
          PRINT*,'Value of A1 =',FPAR(NCLASS,CONTSAVE(NCLASS))
          PRINT*,'Value of A2 =',FPER(NCLASS,CONTSAVE(NCLASS))
          PRINT*,'Total energy =',FENERGTOT
          PRINT*,'XDIPO_PCS energy =',XPCSE
          PRINT*,'CONTSAVE is:',CONTSAVE(NCLASS)
      END IF

      CALL DECLAR('CHIAX','DP',' ',DUMMY2,APARA)
      CALL DECLAR('CHIRH','DP',' ',DUMMY2,APERP)

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE FANTAMEDIA (FMA,PERC)
C     ------------------------------------------------------------------
C     Averages tensor parameters on a given sub-ensemble of structures
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION VETE(2,2000),VETT(3,2000)
      DOUBLE PRECISION FENT(20,2000),FPER(20,2000),FPAR(20,2000),
     &                 FXPCSE(20,2000)
      COMMON /FMEDIA/FENT,FPER,FPAR,FXPCSE,CONTSAVE
      INTEGER CONTSAVE(20)
      INTEGER FMA,COUNT,CB,COUNT1,NUM,CONFPRE
      DOUBLE PRECISION PERC
      DOUBLE PRECISION M1PAR,M1PER,DEVPAR
      DOUBLE COMPLEX DUMMY2
C     ------------------------------------------------------------------
C     Through FMA the class on which the average must be done is passed.
C     If it is zero, all the counters are reset.
C     ------------------------------------------------------------------
      IF (FMA.EQ.0) THEN
          DO COUNT=1,19
             CONTSAVE(COUNT)=0
          END DO
          PRINT*,'Mean vectors reset.'
      ELSE
C     ------------------------------------------------------------------
C     The best are selected
C     ------------------------------------------------------------------
          DO COUNT=1,CONTSAVE(FMA)
             VETE(1,COUNT)=0
          END DO

          DO COUNT=1,CONTSAVE(FMA)
             NUM=CONTSAVE(FMA)
             DO WHILE (VETE(1,NUM).EQ.1)
                NUM=NUM-1
             END DO
             DO COUNT1=1,CONTSAVE(FMA)-1
                IF ((FENT(FMA,COUNT1).LE.FENT(FMA,NUM)).AND.
     &              (VETE(1,COUNT1).EQ.0.0)) THEN
                     NUM=COUNT1
                END IF
             END DO
             PRINT*,'NUM is:',NUM
             VETE(1,NUM)=1
             VETE(2,COUNT)=FENT(FMA,NUM)
             VETT(1,COUNT)=FPAR(FMA,NUM)
             VETT(2,COUNT)=FPER(FMA,NUM)
          END DO

          CONFPRE=NINT((DBLE(CONTSAVE(FMA))*(PERC/100.)))
          IF (CONFPRE.EQ.0) THEN
              CONFPRE=1
          END IF
          PRINT*,'CONFPRE is:',CONFPRE
          PRINT*,'PERC is:',PERC
          PRINT*,'DBLE(CONTSAVE(FMA)) is:',DBLE(CONTSAVE(FMA))

          DO COUNT=1,CONTSAVE(FMA)
             PRINT*,'NUMBER',COUNT,'ENERGY',VETE(1,COUNT),VETE(2,COUNT)
             PRINT*,'A1',VETT(1,COUNT),' A2',VETT(2,COUNT)
          END  DO

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

          CALL DECLAR('CHIAX','DP',' ',DUMMY2,M1PAR)
          CALL DECLAR('CHIRH','DP',' ',DUMMY2,M1PER)

      END IF

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE FANTAERRORE (XPX,XPY,XPZ,XPMX,XPMY,XPMZ,NNAT,XOBS,
     &                        XTOLPROT,METHOD,LOCFLG)
C     ------------------------------------------------------------------
C     Calculates error on tensor parameters through Monte Carlo method
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION APARA,APERP
      COMMON /FRUN/ APERP,APARA
      DOUBLE PRECISION ERRPARA(10000),ERRRH(10000)
      DOUBLE PRECISION FEPERC
      INTEGER FENUMC
      COMMON /FERRORE/FEPERC,FENUMC
      DOUBLE PRECISION FDEPERC
      INTEGER FDENUMC
      COMMON /FDERRORE/FDEPERC,FDENUMC
      DOUBLE PRECISION XPX(1000),XPY(1000),XPZ(1000),
     &                 XOBS(1000),XTOLPROT(1000)
      CHARACTER NOMEFILE*132
      DOUBLE PRECISION  XPMX(1000),XPMY(1000),XPMZ(1000),WPROT(1000)
      INTEGER NNAT,METHOD,I,A,J,IA,AA,SEL
      DOUBLE PRECISION RAND,MEDERRPARA,DEVERRPARA,DEVERRRH,MEDERRRH
      LOGICAL LOCFLG

      IF (METHOD.EQ.0) THEN
          FENUMC=FDENUMC
          FEPERC=FDEPERC
      END IF
      NOMEFILE='MONTECARLO.OUT'
      PRINT*,'FEPERC is:',FEPERC

      DO A=1,FENUMC
         J=0
         DO I=1,NNAT
            WPROT(I)=1
         END DO
         AA=NINT(FEPERC*NNAT)
         PRINT*,'AA is:',AA
         IA=0
         DO WHILE (IA.LT.AA)
            SEL=NINT(RAND()*DBLE(NNAT))
            PRINT*,'SEL IS:',SEL
            IF (WPROT(SEL).EQ.1) THEN
                WPROT(SEL)=0
                IA=IA+1
            END IF
         END DO
         CALL FANTALINEXE(XPX,XPY,XPZ,XPMX,XPMY,XPMZ,NNAT,XOBS,XTOLPROT,
     &                    METHOD,NOMEFILE,WPROT,1,LOCFLG)
         ERRPARA(A)=APARA
         ERRRH(A)=APERP
         PRINT*,'A1 =',APARA
         PRINT*,'A2 =',APERP
      END DO

      MEDERRPARA=0
      MEDERRRH=0
      DO A=1,FENUMC
         MEDERRPARA=ERRPARA(A)+MEDERRPARA
         MEDERRRH=ERRRH(A)+MEDERRRH
      END DO
      MEDERRPARA=MEDERRPARA/DBLE(FENUMC)
      MEDERRRH=MEDERRRH/DBLE(FENUMC)

      DEVERRRH=0
      DEVERRPARA=0
      DO A=1,FENUMC
         DEVERRPARA=DEVERRPARA+(ERRPARA(A)-MEDERRPARA)**2
         DEVERRRH=DEVERRRH+(ERRRH(A)-MEDERRRH)**2
      END DO
      DEVERRPARA=DEVERRPARA/DBLE(FENUMC)
      DEVERRRH=DEVERRRH/DBLE(FENUMC)

      PRINT*,'Standard deviation for A1 =',SQRT(DEVERRPARA)
      PRINT*,'Standard deviation for A2 =',SQRT(DEVERRRH)

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE EXPCS (EDI,WHICH)
C     ------------------------------------------------------------------
C     Calls EXPCS2, which does the actual energy calculation
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE "cns.inc"
      INCLUDE "xdipo_pcs.inc"
      INCLUDE "heap.inc"

      DOUBLE PRECISION EDI
      CHARACTER*7 WHICH

      CALL EXPCS2 (EDI,HEAP(XPCSIPTR),HEAP(XPCSJPTR),HEAP(XPCSKPTR),
     &                 HEAP(XPCSLPTR),HEAP(XPCSMPTR),HEAP(XPCSNPTR),
     &                 HEAP(XPCSILST),HEAP(XPCSJLST),HEAP(XPCSKLST),
     &                 HEAP(XPCSLLST),HEAP(XPCSMLST),HEAP(XPCSNLST),
     &                 HEAP(XPCSOBSPTR),HEAP(XPCSERRPTR),
     &                 HEAP(CALCXPCSPTR),WHICH)
      
      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE EXPCS2 (EDI,ATOMIPTR,ATOMJPTR,ATOMKPTR,
     &                       ATOMLPTR,ATOMMPTR,ATOMNPTR,
     &                       ATOMILST,ATOMJLST,ATOMKLST,
     &                       ATOMLLST,ATOMMLST,ATOMNLST,
     &                       XPCSOBS,XPCSERR,XPCSCALC,WHICH)
C     ------------------------------------------------------------------
C     Calculates pseudocontact shift energies
C
C     Energies are of the form:
C        E = K*(DELTAPCS**2)
C     where:
C        K = FORCE CONSTANT
C        DELTAPCS = CALCULATED PCS - OBSERVED PCS
C     and pseudocontact shift function is defined as:
C        PCS=A1*(3*COS(THETA)^2-1)+(3/2)*A2*(SIN(THETA)^2*COS(2*PHI))/R3
C
C     WHICH is a flag that switches between energy/force calculations
C     and PCS calculations for violations
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "coord.inc"
      INCLUDE  "numbers.inc"
      INCLUDE  "deriv.inc"
      INCLUDE  "xdipo_pcs.inc"
      INCLUDE  "consta.inc"

      DOUBLE PRECISION XPCSE
      COMMON /XPCSENERGY/XPCSE
      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*)
      INTEGER ATOMLPTR(*),ATOMMPTR(*),ATOMNPTR(*)
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*)
      INTEGER ATOMLLST(*),ATOMMLST(*),ATOMNLST(*)
      DOUBLE PRECISION XPCSOBS(*),XPCSERR(*)
      DOUBLE PRECISION XPCSCALC(*),EDI
      CHARACTER*7 WHICH
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      INTEGER COUNT,CLASS,MM,NN,II,I,J,K,L,NA,M,N
      PARAMETER (NA=20)
      INTEGER DMM,DNN,OUTFLAG(NA),MCOUNT,NCOUNT,PCOUNT
      DOUBLE PRECISION XI,XJ,XK,XL,XM,XN,YI,YJ,YK,YL,YM,YN,
     &                 ZI,ZJ,ZK,ZL,ZM,ZN,
     &                 E,DF1(NA),DF2(NA),DF3(NA),
     &                 PHI(NA),DPXI(NA),DPYI(NA),DPZI(NA),
     &                 DPXJ(NA),DPYJ(NA),DPZJ(NA),DPXK(NA),
     &                 DPYK(NA),DPZK(NA),DPXL(NA),DPYL(NA),
     &                 DPZL(NA),DPXM(NA),DPYM(NA),DPZM(NA),
     &                 DPXN(NA),DPYN(NA),DPZN(NA),
     &                 THETA(NA),DTXI(NA),DTYI(NA),DTZI(NA),
     &                 DTXJ(NA),DTYJ(NA),DTZJ(NA),DTXK(NA),
     &                 DTYK(NA),DTZK(NA),DTXL(NA),DTYL(NA),
     &                 DTZL(NA),DTXM(NA),DTYM(NA),DTZM(NA),
     &                 DTXN(NA),DTYN(NA),DTZN(NA),DRXM(NA),
     &                 DRYM(NA),DRZM(NA),DRXN(NA),DRYN(NA),
     &                 DRZN(NA),DIST,
     &                 O1,O2,O3,O4,O5,O6,O7,O8,O9,O10,
     &                 O11,O12,O13,O14,O15,O16,O17,O18,O19,
     &                 O20,O21,O22,O23,O24,O25,O26,O27,O28,
     &                 O29,O30,O31,O32,O33,O34,O35,O36,O37,
     &                 O38,O39,O40,O41,O42,O43,O44,O45,O46,
     &                 O47,O48,O49,O50,O51,O52,O53,O54,O55,
     &                 O56,O57,O58,O59,O60,O61,O62,O63,O64,
     &                 O65,O66,O67,O68,O69,O70,O71,O72,O73,
     &                 O74,O75,O76,O77,O78,O79,O80,O81,O82
      DOUBLE PRECISION OBSXPCS,ERRXPCS,K1,
     &                 COEF1,COEF2,A,B,DELTA,
     &                 DT,DP,DELTAXPCS,
     &                 DD,DR,CALCXPCS
C     ------------------------------------------------------------------
C     Zero out partial energy
C     ------------------------------------------------------------------
      EDI=ZERO

      CLASS=1
      K1=XPCSFORCES(1)
      COEF1=XPCSCOEF1(1)
      COEF2=XPCSCOEF2(1)

      COUNT=0
      DO WHILE (COUNT.LT.PCSNUM)
         COUNT=COUNT+1
C     ------------------------------------------------------------------
C     Reset individual E to zero
C     ------------------------------------------------------------------
         E=0

         DO WHILE ((XPCSASSNDX(CLASS).LT.COUNT).OR.
     &             (XPCSASSNDX(CLASS).EQ.0))
            CLASS=CLASS+1
         END DO

         IF (XPCSASSNDX(CLASS).EQ.XPCSASSNDX(CLASS-1)) THEN
             COUNT=COUNT-1
         END IF

         K1=XPCSFORCES(CLASS)
         COEF1=XPCSCOEF1(CLASS)
         COEF2=XPCSCOEF2(CLASS)
C     ------------------------------------------------------------------
C     Note there should only be one atom for I, J, K, and L
C     ------------------------------------------------------------------
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

         OBSXPCS=XPCSOBS(COUNT)
         ERRXPCS=XPCSERR(COUNT)
C     ------------------------------------------------------------------
C     Initialize calculated pseudocontact shift and counter
C     ------------------------------------------------------------------
         CALCXPCS=0
         II=0
         MCOUNT=0
         NCOUNT=0
C     ------------------------------------------------------------------
C     Check for correct permutations of paired atoms to get the metal-
C     nucleus vector. This depends solely on the order of the assigned 
C     atoms. OUTFLAG=1 indicates that the permutation is not allowed.
C     ------------------------------------------------------------------
         DMM=ATOMMPTR(COUNT+1)-ATOMMPTR(COUNT)
         DNN=ATOMNPTR(COUNT+1)-ATOMNPTR(COUNT)
         DO MM=ATOMMPTR(COUNT)+1,ATOMMPTR(COUNT+1)
            MCOUNT=MCOUNT+1
            NCOUNT=0
            DO NN=ATOMNPTR(COUNT)+1,ATOMNPTR(COUNT+1)
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
         DO MM=ATOMMPTR(COUNT)+1,ATOMMPTR(COUNT+1)
            DO NN=ATOMNPTR(COUNT)+1,ATOMNPTR(COUNT+1)
               II=II+1
               IF (OUTFLAG(II).NE.1) THEN
                  PCOUNT=PCOUNT+1
                  XM=X(ATOMMLST(MM))
                  XN=X(ATOMNLST(NN))
                  YM=Y(ATOMMLST(MM))
                  YN=Y(ATOMNLST(NN))
                  ZM=Z(ATOMMLST(MM))
                  ZN=Z(ATOMNLST(NN))
C     ------------------------------------------------------------------
C     Calculate THETA, PHI and derivatives with respect to X, Y, Z of 
C     atoms I, J, K, L, M, N. This was done with Mathematica by Nico!
C     ------------------------------------------------------------------
                  O1=-XI+XJ
                  O2=O1**2
                  O3=-YI+YJ
                  O4=O3**2
                  O5=-ZI+ZJ
                  O6=O5**2
                  O7=O2+O4+O6
                  O8=SQRT(O7)
                  O9=1/O8
                  O10=-XM+XN
                  O11=O1*O10
                  O12=-YM+YN
                  O13=O12*O3
                  O14=-ZM+ZN
                  O15=O14*O5
                  O16=O11+O13+O15
                  O17=O10**2
                  O18=O12**2
                  O19=O14**2
                  O20=O17+O18+O19
                  O21=SQRT(O20)
                  O22=1/O21
                  O23=1/O7
                  O24=O16**2
                  O25=1/O20
                  O26=O23*O24*O25
                  O27=SQRT(1.D0-O26)
                  O28=1/O27
                  O29=XM-XN
                  O30=O7**(3.D0/2.D0)
                  O31=1/O30
                  O32=O1*O16*O22*O31
                  O33=YM-YN
                  O34=O16*O22*O3*O31
                  O35=ZM-ZN
                  O36=O16*O22*O31*O5
                  O37=O20**(3.D0/2.D0)
                  O38=1/O37
                  O39=O10*O16*O38*O9
                  O40=O12*O16*O38*O9
                  O41=O14*O16*O38*O9
                  O42=-XI+XK
                  O43=O42**2
                  O44=-YI+YK
                  O45=O44**2
                  O46=-ZI+ZK
                  O47=O46**2
                  O48=O43+O45+O47
                  O49=SQRT(O48)
                  O50=-XI+XL
                  O51=O50**2
                  O52=-YI+YL
                  O53=O52**2
                  O54=-ZI+ZL
                  O55=O54**2
                  O56=O51+O53+O55
                  O57=SQRT(O56)
                  O58=1/O57
                  O59=O10*O42
                  O60=O12*O44
                  O61=O14*O46
                  O62=O59+O60+O61
                  O63=1/O62
                  O64=O10*O50
                  O65=O12*O52
                  O66=O14*O54
                  O67=O64+O65+O66
                  O68=O62**2
                  O69=1/O68
                  O70=O56**(3.D0/2.D0)
                  O71=1/O70
                  O72=O49*O50*O63*O67*O71
                  O73=1/O49
                  O74=O42*O58*O63*O67*O73
                  O75=1/O56
                  O76=O67**2
                  O77=O48*O69*O75*O76
                  O78=1/(1.D0+O77)
                  O79=O49*O52*O63*O67*O71
                  O80=O44*O58*O63*O67*O73
                  O81=O49*O54*O63*O67*O71
                  O82=O46*O58*O63*O67*O73

                  THETA(II)=ACOS(O16*O22*O9)

                  DTXI(II)=-(O28*(O32+O22*O29*O9))
                  DTYI(II)=-(O28*(O34+O22*O33*O9))
                  DTZI(II)=-(O28*(O36+O22*O35*O9))
                  DTXJ(II)=-(O28*(-O32+O10*O22*O9))
                  DTYJ(II)=-(O28*(-O34+O12*O22*O9))
                  DTZJ(II)=-(O28*(-O36+O14*O22*O9))
                  DTXK(II)=0
                  DTYK(II)=0
                  DTZK(II)=0
                  DTXL(II)=0
                  DTYL(II)=0
                  DTZL(II)=0
                  DTXM(II)=-(O28*(O39+O22*O9*(XI-XJ)))
                  DTYM(II)=-(O28*(O40+O22*O9*(YI-YJ)))
                  DTZM(II)=-(O28*(O41+O22*O9*(ZI-ZJ)))
                  DTXN(II)=-(O28*(-O39+O1*O22*O9))
                  DTYN(II)=-(O28*(-O40+O22*O3*O9))
                  DTZN(II)=-(O28*(-O41+O22*O5*O9))

                  PHI(II)=ATAN(O49*O58*O63*O67)

                  DPXI(II)=O78*(O29*O49*O58*O63-O29*O49*O58*O67*O69+
     &                          O72-O74)
                  DPYI(II)=O78*(O33*O49*O58*O63-O33*O49*O58*O67*O69+
     &                          O79-O80)
                  DPZI(II)=O78*(O35*O49*O58*O63-O35*O49*O58*O67*O69+
     &                          O81-O82)
                  DPXJ(II)=0
                  DPYJ(II)=0
                  DPZJ(II)=0
                  DPXK(II)=(-(O10*O49*O58*O67*O69)+O74)*O78
                  DPYK(II)=O78*(-(O12*O49*O58*O67*O69)+O80)
                  DPZK(II)=O78*(-(O14*O49*O58*O67*O69)+O82)
                  DPXL(II)=(O10*O49*O58*O63-O72)*O78
                  DPYL(II)=O78*(O12*O49*O58*O63-O79)
                  DPZL(II)=O78*(O14*O49*O58*O63-O81)
                  DPXM(II)=O78*(-(O49*O58*O67*O69*(XI-XK))+
     &                            O49*O58*O63*(XI-XL))
                  DPYM(II)=O78*(-(O49*O58*O67*O69*(YI-YK))+
     &                            O49*O58*O63*(YI-YL))
                  DPZM(II)=O78*(-(O49*O58*O67*O69*(ZI-ZK))+
     &                            O49*O58*O63*(ZI-ZL))
                  DPXN(II)=(O49*O50*O58*O63-O42*O49*O58*O67*O69)*O78
                  DPYN(II)=(O49*O52*O58*O63-O44*O49*O58*O67*O69)*O78
                  DPZN(II)=(O49*O54*O58*O63-O46*O49*O58*O67*O69)*O78
C     ------------------------------------------------------------------
C     Calculate the distance between the two atoms M, N (metal-nucleus)
C     and calculate the derivatives as well
C     ------------------------------------------------------------------
                  O1=-XM+XN
                  O2=O1**2
                  O3=-YM+YN
                  O4=O3**2
                  O5=-ZM+ZN
                  O6=O5**2
                  O7=O2+O4+O6
                  O8=SQRT(O7)
                  O9=1/O8
                  O10=O1*O9
                  O11=O3*O9
                  O12=O5*O9
                  DIST=O8
                  DRXM(II)=-O10
                  DRYM(II)=-O11
                  DRZM(II)=-O12
                  DRXN(II)=O10
                  DRYN(II)=O11
                  DRZN(II)=O12
C     ------------------------------------------------------------------
C     Calculate the pseudocontact shift
C     ------------------------------------------------------------------
                  A=COEF1
                  B=COEF2
                  CALCXPCS=(A*(3.*COS(THETA(II))*COS(THETA(II))-1.)+
     &                      B*(3./2.)*(SIN(THETA(II))*SIN(THETA(II))*
     &                      COS(2.*PHI(II))))/(DIST**3.0)
               END IF
            END DO
         END DO

         IF (WHICH.EQ.'ANALYZE') THEN
            XPCSCALC(COUNT)=CALCXPCS
         END IF

         DELTA=(CALCXPCS-OBSXPCS)
C     ------------------------------------------------------------------
C     Adjust the deviation based on the error range
C     ------------------------------------------------------------------
         IF ((DELTA.LT.0.000).AND.(ABS(DELTA).GT.ERRXPCS)) THEN
              DELTAXPCS=DELTA+ERRXPCS
         ELSE IF ((DELTA.GT.0.000).AND.(DELTA.GT.ERRXPCS)) THEN
              DELTAXPCS=DELTA-ERRXPCS
         ELSE
              DELTAXPCS=0.0
         END IF
         DD=DELTAXPCS

         II = 0
         PCOUNT = 0
         DO MM=ATOMMPTR(COUNT)+1,ATOMMPTR(COUNT+1)
            DO NN=ATOMNPTR(COUNT)+1,ATOMNPTR(COUNT+1)
               II = II+1
               IF (OUTFLAG(II).NE.1) THEN
                   PCOUNT=PCOUNT + 1
C     ------------------------------------------------------------------
C     Taking care of derivative
C     ------------------------------------------------------------------
                   DT=(1.0/(DIST**3))*(A*(-6.)*SIN(THETA(II))*
     &                   COS(THETA(II))+B*3.*SIN(THETA(II))*
     &                   COS(THETA(II))*COS(2.*PHI(II)))
                   DP=(1.0/DIST**3)*(B*(-3.)*SIN(THETA(II))*
     &                   SIN(THETA(II))*SIN(2.*PHI(II)))
                   DR=(-3.0/(DIST**4))*(A*(3.*COS(THETA(II))*
     &                   COS(THETA(II))-1.)+B*(3./2.)*
     &                   (SIN(THETA(II))*SIN(THETA(II))*
     &                   COS(2.*PHI(II))))
                  IF (DD.EQ.0.0) THEN
                      E=E+0.0
                      DF1(II)=0.0
                      DF2(II)=0.0
                      DF3(II)=0.0
                  ELSE
                      E=E+K1*(DELTAXPCS**2)
                      DF1(II)=2*K1*DELTAXPCS*DT
                      DF2(II)=2*K1*DELTAXPCS*DP
                      DF3(II)=2*K1*DELTAXPCS*DR
                  END IF
               END IF
            END DO
         END DO
C     ------------------------------------------------------------------
C     Accumulate energy
C     ------------------------------------------------------------------
         EDI=EDI+E
C     ------------------------------------------------------------------
C     Now update forces if in energy/force mode
C     ------------------------------------------------------------------
         IF (WHICH.NE.'ANALYZE') THEN
             II = 0
             DO MM=ATOMMPTR(COUNT)+1,ATOMMPTR(COUNT+1)
                DO NN=ATOMNPTR(COUNT)+1,ATOMNPTR(COUNT+1)
                   II=II+1
                   IF (OUTFLAG(II).NE.1) THEN
                       DX(ATOMILST(I))=DX(ATOMILST(I))+
     &                                 DF1(II)*DTXI(II)+DF2(II)*DPXI(II)
                       DY(ATOMILST(I))=DY(ATOMILST(I))+
     &                                 DF1(II)*DTYI(II)+DF2(II)*DPYI(II)
                       DZ(ATOMILST(I))=DZ(ATOMILST(I))+
     &                                 DF1(II)*DTZI(II)+DF2(II)*DPZI(II)
                       DX(ATOMJLST(J))=DX(ATOMJLST(J))+
     &                                 DF1(II)*DTXJ(II)+DF2(II)*DPXJ(II)
                       DY(ATOMJLST(J))=DY(ATOMJLST(J))+
     &                                 DF1(II)*DTYJ(II)+DF2(II)*DPYJ(II)
                       DZ(ATOMJLST(J))=DZ(ATOMJLST(J))+
     &                                 DF1(II)*DTZJ(II)+DF2(II)*DPZJ(II)
                       DX(ATOMKLST(K))=DX(ATOMKLST(K))+
     &                                 DF1(II)*DTXK(II)+DF2(II)*DPXK(II)
                       DY(ATOMKLST(K))=DY(ATOMKLST(K))+
     &                                 DF1(II)*DTYK(II)+DF2(II)*DPYK(II)
                       DZ(ATOMKLST(K))=DZ(ATOMKLST(K))+
     &                                 DF1(II)*DTZK(II)+DF2(II)*DPZK(II)
                       DX(ATOMLLST(L))=DX(ATOMLLST(L))+
     &                                 DF1(II)*DTXL(II)+DF2(II)*DPXL(II)
                       DY(ATOMLLST(L))=DY(ATOMLLST(L))+
     &                                 DF1(II)*DTYL(II)+DF2(II)*DPYL(II)
                       DZ(ATOMLLST(L))=DZ(ATOMLLST(L))+
     &                                 DF1(II)*DTZL(II)+DF2(II)*DPZL(II)
                       DX(ATOMMLST(MM))=DX(ATOMMLST(MM))+
     &                                  DF1(II)*DTXM(II)+DF2(II)*
     &                                  DPXM(II)+DF3(II)*DRXM(II)
                       DY(ATOMMLST(MM))=DY(ATOMMLST(MM))+
     &                                  DF1(II)*DTYM(II)+DF2(II)*
     &                                  DPYM(II)+DF3(II)*DRYM(II)
                       DZ(ATOMMLST(MM))=DZ(ATOMMLST(MM))+
     &                                  DF1(II)*DTZM(II)+DF2(II)*
     &                                  DPZM(II)+DF3(II)*DRZM(II)
                       DX(ATOMNLST(NN))=DX(ATOMNLST(NN))+
     &                                  DF1(II)*DTXN(II)+DF2(II)*
     &                                  DPXN(II)+DF3(II)*DRXN(II)
                       DY(ATOMNLST(NN))=DY(ATOMNLST(NN))+
     &                                  DF1(II)*DTYN(II)+DF2(II)*
     &                                  DPYN(II)+DF3(II)*DRYN(II)
                       DZ(ATOMNLST(NN))=DZ(ATOMNLST(NN))+
     &                                  DF1(II)*DTZN(II)+DF2(II)*
     &                                  DPZN(II)+DF3(II)*DRZN(II)

                   END IF
                END DO
             END DO
         END IF
      END DO

      XPCSE=EDI

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE READXPCS
C     ------------------------------------------------------------------
C     Reads in pseudocontact shift information
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE "cns.inc"
      INCLUDE "comand.inc"
      INCLUDE "xdipo_pcs.inc"
      INCLUDE "funct.inc"
      INCLUDE "mtf.inc"
      INCLUDE "heap.inc"
      INCLUDE "numbers.inc"

      INTEGER FSWITCH,FSAVE,FERR
      COMMON /FPARAM/FSWITCH,FSAVE,FERR
      INTEGER FVIOLS
      COMMON /VIOL/FVIOLS
      DOUBLE PRECISION FEPERC
      INTEGER FENUMC
      COMMON /FERRORE/FEPERC,FENUMC

      INTEGER COUNT,SPTR,OLDCLASS,OLDMAXXPCS,FNCLASS
      DOUBLE PRECISION K1,CUTOFF,COEF1,COEF2
      DOUBLE PRECISION FMEDPERC
      CHARACTER*4 THENAME
      CHARACTER*132 NOMEFILE
      INTEGER TOLSWI

      NOMEFILE='XDIPO_PCS.OUT'

      SPTR=ALLHP(INTEG4(NATOM))
      CALL PUSEND('XDIPO_PCS>')
862   CONTINUE
      CALL NEXTWD('XDIPO_PCS>')
      CALL MISCOM('XDIPO_PCS>',USED)
      IF (.NOT.USED) THEN
C     ------------------------------------------------------------------
C     Documentation
C     ------------------------------------------------------------------
          IF (WD(1:4).EQ.'HELP') THEN
              WRITE(DUNIT,'(10X,A)')
     &              ' XPCS {<PCS-STATEMENT>} END ',
     &              ' <PCS-STATEMENT>:== ',
     &              ' ASSIgn <SEL> <SEL> <SEL> <SEL> <SEL>',
     &              ' <REAL> <REAL>',
     &              ' * Restraint: Metal Z X Y Nucleus PCS Err *',
     &              ' CLASsification <NAME>',
     &              ' * Starts a new class *',
     &              ' TOLL <1|0>',
     &              ' * Switch tolerance in FANTALIN (1=Yes, 0=No)',
     &              ' COEFficient <REAL> <REAL>',
     &              ' * Tensor parameters: A1 A2 *',
     &              ' FORCeconstant <REAL>',
     &              ' * Force constant for the current class *',
     &              ' NREStraints <INTEGER>',
     &              ' * Number of slots for restraints to allocate *',
     &              ' PRINt THREshold <REAL>',
     &              ' * Prints violations larger than the threshold *',
     &              ' RESEt',
     &              ' * Erases the restraint table, but keeps NRES *',
     &              ' SAVE',
     &              ' * Filename to save PCS values *',
     &              ' FMED',
     &              ' * Calculates average tensor parameters *',
     &              ' ERRON',
     &              ' * Switches Monte Carlo error evaluation ON *',
     &              ' ERROFF',
     &              ' * Switches Monte Carlo error evaluation OFF *',
     &              ' FON',
     &              ' * Switches tensor auto-update mode ON *',
     &              ' FOFF',
     &              ' * Switches tensor auto-update mode OFF *',
     &              ' SON',
     &              ' * Switches saving mode ON *',
     &              ' SOFF',
     &              ' * Switches saving mode OFF *',
     &              ' FRUN',
     &              ' * Runs FANTALIN *'
C     ------------------------------------------------------------------
C     About FANTALIN
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'TOLL') THEN
              CALL NEXTI('Tolerance setting in FANTALIN (1/0) =',TOLSWI)
              IF (TOLSWI.EQ.1) THEN
                  PCSFLG=.TRUE.
                  WRITE (DUNIT,'(A)') 'Tolerance ON in FANTALIN.'
              ELSE IF (TOLSWI.EQ.0) THEN
                  PCSFLG=.FALSE.
                  WRITE (DUNIT,'(A)') 'Tolerance OFF in FANTALIN.'
              ELSE
                  WRITE (DUNIT,'(A)') 'Unknown switch. Using default.'
              END IF
          ELSE IF (WD(1:4).EQ.'SAVE') THEN
              CALL NEXTFI('FILENAME =',NOMEFILE)
              PRINT*,'Name of the file to save PCS values :',NOMEFILE
          ELSE IF (WD(1:3).EQ.'FME') THEN
              CALL NEXTF('Percentage for average (0/100) =',FMEDPERC)
              CALL NEXTI('Average on class number =',FNCLASS)
              IF (FNCLASS.GT.NCLASSES) THEN
                  PRINT*,'%FRUN-ERR: This class does not exist...'
              ELSE
                  CALL FANTAMEDIA(FNCLASS,FMEDPERC)
              END IF
          ELSE IF (WD(1:5).EQ.'ERRON') THEN
              CALL NEXTI('Number of cycles =',FENUMC)
              CALL NEXTF('Percentage to discard (0/100) =',FEPERC)
              FEPERC=FEPERC/100.0
              FERR=1
              PRINT*,'Monte Carlo error evaluation ON'
          ELSE IF (WD(1:6).EQ.'ERROFF') THEN
              FERR=0
              PRINT*,'Monte Carlo error evaluation OFF'
          ELSE IF (WD(1:3).EQ.'FON') THEN
              FSWITCH=1
              PRINT*,'Tensor auto-update mode ON'
          ELSE IF (WD(1:4).EQ.'FOFF') THEN
              FSWITCH=0
              PRINT*,'Tensor auto-update mode OFF'
          ELSE IF (WD(1:3).EQ.'SON') THEN
              FSAVE=1
              PRINT*,'Saving mode ON'
          ELSE IF (WD(1:4).EQ.'SOFF') THEN
              FSAVE=0
              PRINT*,'Saving mode OFF'
          ELSE IF (WD(1:4).EQ.'FRUN') THEN
              CALL NEXTI('FANTALIN on class number =',FNCLASS)
              IF (FNCLASS.GT.NCLASSES.OR.FNCLASS.EQ.0) THEN
                  PRINT*,'%FRUN-ERR: This class does not exist...'
             ELSE
                  CALL FANTALIN(HEAP(XPCSIPTR),HEAP(XPCSJPTR),
     &                          HEAP(XPCSKPTR),HEAP(XPCSLPTR),
     &                          HEAP(XPCSMPTR),HEAP(XPCSNPTR),
     &                          HEAP(XPCSILST),HEAP(XPCSJLST),
     &                          HEAP(XPCSKLST),HEAP(XPCSLLST),
     &                          HEAP(XPCSMLST),HEAP(XPCSNLST),
     &                          HEAP(XPCSOBSPTR),HEAP(XPCSERRPTR),
     &                          FNCLASS,NOMEFILE )
             END  IF
C     ------------------------------------------------------------------
C     Get class name. Determine if it's an already-defined class.
C     Insert a new class if it's not.
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'CLAS') THEN
              OLDCLASS=CURCLASS
              CALL NEXTA4('Class name =',THENAME)
              MODE=NEW
              DO COUNT=1,NCLASSES
                 IF (XPCSCLASSNAMES(COUNT).EQ.THENAME) THEN
                     MODE=UPDATE
                     CURCLASS=COUNT
                 END IF
              END DO
              IF (MODE.EQ.NEW) THEN
C     ------------------------------------------------------------------
C     Make sure you can't add more than the maximum number of classes
C     ------------------------------------------------------------------
                  IF (OLDCLASS.EQ.MAXXPCSCLASSES) THEN
                      CALL DSPERR('XDIPO_PCS','Too many classes.')
                      CALL DSPERR('XDIPO_PCS',
     &                     'Increase MAXXPCSCLASSES and recompile.')
                      CALL WRNDIE(-5,'READXPCS',
     &                     'Too many anisotropy classes.')
                  END IF
                  NCLASSES=NCLASSES+1
                  CURCLASS=NCLASSES
                  XPCSCLASSNAMES(CURCLASS)=THENAME

C     ------------------------------------------------------------------
C     If this isn't the first class, close off the old class
C     ------------------------------------------------------------------
                  IF (NCLASSES.GT.1) THEN
                      XPCSASSNDX(OLDCLASS)=PCSNUM
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
                  NCLASSES=1
                  CURCLASS=1
              END IF
              CALL NEXTF('Force constant =',K1)
              WRITE(DUNIT,'(A,A,A,F8.3)')
     &              'Setting force constant for class ',
     &              XPCSCLASSNAMES(CURCLASS),' to ',K1
                    XPCSFORCES(CURCLASS)=K1
C     ------------------------------------------------------------------
C     Set coefficient constant for current class
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'COEF') THEN
              CALL NEXTF('Coefficient A1 =',COEF1)
              CALL NEXTF('Coefficient A2 =',COEF2)
C     ------------------------------------------------------------------
C     Start a default class if there isn't one defined
C     ------------------------------------------------------------------
              IF (CURCLASS.EQ.0) THEN
                  NCLASSES=1
                  CURCLASS=1
              END IF
              WRITE(DUNIT,'(A,A,A,F8.3,F8.3)')
     &              'Setting coefficients for class ',
     &              XPCSCLASSNAMES(CURCLASS),' to ',COEF1,COEF2
                    XPCSCOEF1(CURCLASS)=COEF1
                    XPCSCOEF2(CURCLASS)=COEF2
C     ------------------------------------------------------------------
C     Reset pseudocontact shifts database
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'RESE') THEN
              CALL XPCSDEFAULTS
              CALL ALLOCXPCS(0,MAXXPCS)
C     ------------------------------------------------------------------
C     Change number of assignment slots
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'NRES') THEN
               OLDMAXXPCS=MAXXPCS
               CALL NEXTI('Number of slots =',MAXXPCS)
               CALL ALLOCXPCS(OLDMAXXPCS,MAXXPCS)
               WRITE(DUNIT,'(A,I8,A)')
     &               'XDIPO_PCS: Allocating space for',MAXXPCS,
     &               'number of restraints.'
C     ------------------------------------------------------------------
C     Read in an assignment
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'ASSI') THEN
C     ------------------------------------------------------------------
C     Make sure you can't add more assignments than you have slots for
C     ------------------------------------------------------------------
              IF (XPCSMX.EQ.MAXXPCS) THEN
                  CALL DSPERR('XDIPO_PCS','Too many assignments.')
                  CALL DSPERR('XDIPO_PCS',
     &                        'Increasing NREStraints by 100.')
                  OLDMAXXPCS=MAXXPCS
                  MAXXPCS=MAXXPCS+100
                  CALL ALLOCXPCS(OLDMAXXPCS,MAXXPCS)
              END IF
C     ------------------------------------------------------------------
C     If there isn't a class specified, start a default class
C     ------------------------------------------------------------------
              IF (CURCLASS.EQ.0) THEN
                  NCLASSES=1
                  CURCLASS=1
              END IF
              CALL READXPCS2(HEAP(XPCSIPTR),HEAP(XPCSJPTR),
     &                       HEAP(XPCSKPTR),HEAP(XPCSLPTR),
     &                       HEAP(XPCSMPTR),HEAP(XPCSNPTR),
     &                       HEAP(XPCSILST),HEAP(XPCSJLST),
     &                       HEAP(XPCSKLST),HEAP(XPCSLLST),
     &                       HEAP(XPCSMLST),HEAP(XPCSNLST),
     &                       HEAP(XPCSOBSPTR),HEAP(XPCSERRPTR),
     &                       HEAP(SPTR))
C     ------------------------------------------------------------------
C     Print violations
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'PRIN') THEN
              CALL NEXTWD('PRINT>')
              IF (WD(1:4).NE.'THRE') THEN
                  CALL DSPERR('XDIPO_PCS',
     &                 'Print expects THREshold parameter.')
              ELSE
                  CALL NEXTF('THREshold =',CUTOFF)
                  IF (CUTOFF.LT.ZERO) THEN
                      CALL DSPERR('XDIPO_PCS',
     &                     'Cutoff must be positive.')
                      CUTOFF=ABS(CUTOFF)
                  END IF
                  CALL NEXTA4('ALL OR CLASs>',THENAME)
                  IF (THENAME(1:3).EQ.'ALL') THEN
                      DO COUNT=1,NCLASSES
                         PRINTCLASS(COUNT)=.TRUE.
                      END DO
                  ELSE IF (THENAME(1:4).EQ.'CLAS') THEN
                      CALL NEXTA4('CLASS NAME =',THENAME)
                      DO COUNT=1,NCLASSES
                         IF (XPCSCLASSNAMES(COUNT).EQ.THENAME) THEN
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
                  CALL PRINTXPCS(CUTOFF,HEAP(CALCXPCSPTR),
     &                           HEAP(XPCSOBSPTR),HEAP(XPCSERRPTR),
     &                           HEAP(XPCSIPTR),HEAP(XPCSJPTR),
     &                           HEAP(XPCSKPTR),HEAP(XPCSLPTR),
     &                           HEAP(XPCSMPTR),HEAP(XPCSNPTR),
     &                           HEAP(XPCSILST),HEAP(XPCSJLST),
     &                           HEAP(XPCSKLST),HEAP(XPCSLLST),
     &                           HEAP(XPCSMLST),HEAP(XPCSNLST))
              END IF
C     ------------------------------------------------------------------
C     Check for END statement
C     ------------------------------------------------------------------
          ELSE
              CALL CHKEND('XDIPO_PCS>',DONE)
          END IF
      END IF
      IF (.NOT.DONE) GOTO 862
      DONE=.FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))
      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE ALLOCXPCS (OLDSIZE,NEWSIZE)
C     ------------------------------------------------------------------
C     Reset pseudocontact shift arrays to hold size entries
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "funct.inc"
      INCLUDE  "xdipo_pcs.inc"

      INTEGER OLDSIZE,NEWSIZE

      IF (OLDSIZE.NE.0) THEN
          CALL FREHP(XPCSIPTR,INTEG4(OLDSIZE))
          CALL FREHP(XPCSJPTR,INTEG4(OLDSIZE))
          CALL FREHP(XPCSKPTR,INTEG4(OLDSIZE))
          CALL FREHP(XPCSLPTR,INTEG4(OLDSIZE))
          CALL FREHP(XPCSMPTR,INTEG4(OLDSIZE))
          CALL FREHP(XPCSNPTR,INTEG4(OLDSIZE))

          CALL FREHP(XPCSILST,INTEG4(OLDSIZE))
          CALL FREHP(XPCSJLST,INTEG4(OLDSIZE))
          CALL FREHP(XPCSKLST,INTEG4(OLDSIZE))
          CALL FREHP(XPCSLLST,INTEG4(OLDSIZE))
          CALL FREHP(XPCSMLST,INTEG4(OLDSIZE))
          CALL FREHP(XPCSNLST,INTEG4(OLDSIZE))

          CALL FREHP(XPCSOBSPTR,IREAL8(OLDSIZE))
          CALL FREHP(XPCSERRPTR,IREAL8(OLDSIZE))
          CALL FREHP(CALCXPCSPTR,IREAL8(OLDSIZE))
      END IF
      XPCSIPTR=ALLHP(INTEG4(NEWSIZE))
      XPCSJPTR=ALLHP(INTEG4(NEWSIZE))
      XPCSKPTR=ALLHP(INTEG4(NEWSIZE))
      XPCSLPTR=ALLHP(INTEG4(NEWSIZE))
      XPCSMPTR=ALLHP(INTEG4(NEWSIZE))
      XPCSNPTR=ALLHP(INTEG4(NEWSIZE))

      XPCSILST=ALLHP(INTEG4(NEWSIZE))
      XPCSJLST=ALLHP(INTEG4(NEWSIZE))
      XPCSKLST=ALLHP(INTEG4(NEWSIZE))
      XPCSLLST=ALLHP(INTEG4(NEWSIZE))
      XPCSMLST=ALLHP(INTEG4(NEWSIZE))
      XPCSNLST=ALLHP(INTEG4(NEWSIZE))

      XPCSOBSPTR=ALLHP(IREAL8(NEWSIZE))
      XPCSERRPTR=ALLHP(IREAL8(NEWSIZE))
      CALCXPCSPTR=ALLHP(IREAL8(NEWSIZE))

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE XPCSDEFAULTS
C     ------------------------------------------------------------------
C     Sets up defaults
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "xdipo_pcs.inc"

      INTEGER COUNT

      MODE=NEW
      MAXXPCS=200
      XPCSMX=200
      PCSNUM=0
      NCLASSES=0
      CURCLASS=0
      PCSFLG=.FALSE.
      DO COUNT=1,MAXXPCSCLASSES
         XPCSCLASSNAMES(COUNT)='DEFAULT'
         XPCSASSNDX(COUNT)=0
         XPCSFORCES(COUNT)=5.0
         XPCSCOEF1(COUNT)=1.0
         XPCSCOEF2(COUNT)=1.0
      END DO

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE READXPCS2(ATOMIPTR,ATOMJPTR,ATOMKPTR,ATOMLPTR,ATOMMPTR,
     &                     ATOMNPTR,ATOMILST,ATOMJLST,ATOMKLST,ATOMLLST,
     &                     ATOMMLST,ATOMNLST,XPCSOBS,XPCSERR,SEL)
C     ------------------------------------------------------------------
C     Reads actual pseudocontact shift assignments into arrays
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "coord.inc"
      INCLUDE  "xdipo_pcs.inc"

      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*),ATOMLPTR(*)
      INTEGER ATOMMPTR(*),ATOMNPTR(*),ATOMILST(*)
      INTEGER ATOMJLST(*),ATOMKLST(*),ATOMLLST(*),ATOMMLST(*)
      INTEGER ATOMNLST(*),SEL(*)
      INTEGER ITMP,JTMP,KTMP,LTMP,II
      DOUBLE PRECISION XPCSOBS(*),XPCSERR(*)
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      INTEGER NSEL,INSERTPOS,COUNT,CURSTOP,OTHERSTOP,NFLAG
      INTEGER I
      DOUBLE PRECISION XPCSO,XPCSE
C     ------------------------------------------------------------------
C     If we're in update mode, make a space for the new line
C     ------------------------------------------------------------------
      NFLAG=0
      IF (MODE.EQ.UPDATE) THEN
          DO COUNT=PCSNUM+1,XPCSASSNDX(CURCLASS)+1,-1
            ATOMIPTR(COUNT)=ATOMIPTR(COUNT-1)
            ATOMJPTR(COUNT)=ATOMJPTR(COUNT-1)
            ATOMKPTR(COUNT)=ATOMKPTR(COUNT-1)
            ATOMLPTR(COUNT)=ATOMLPTR(COUNT-1)
            ATOMMPTR(COUNT)=ATOMMPTR(COUNT-1)
            ATOMNPTR(COUNT)=ATOMNPTR(COUNT-1)
            XPCSOBS(COUNT)=XPCSOBS(COUNT-1)
            XPCSERR(COUNT)=XPCSERR(COUNT-1)
          END DO
          CURSTOP=XPCSASSNDX(CURCLASS)
          DO COUNT=1,NCLASSES
             OTHERSTOP=XPCSASSNDX(COUNT)
             IF (OTHERSTOP.GT.CURSTOP) THEN
                 XPCSASSNDX(COUNT)=OTHERSTOP+1
             END IF
          END DO
          XPCSASSNDX(CURCLASS)=CURSTOP+1
          INSERTPOS=CURSTOP
          PCSNUM=PCSNUM+1
      ELSE
          PCSNUM=PCSNUM+1
          INSERTPOS=PCSNUM
          XPCSASSNDX(CURCLASS)=INSERTPOS
      END IF
C     ------------------------------------------------------------------
C     Reading in the atom selection in the restraint table 
C     ------------------------------------------------------------------
      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_PCS',
     &                'More than 1 atom in SEL for atom I. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
         CALL DSPERR('XDIPO_PCS','Error with atom selection')
         NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL)
      IF (INSERTPOS.EQ.1) ATOMIPTR(INSERTPOS)=0
      II=ATOMIPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XPCSMX) XPCSMX=II
      ATOMILST(II)=SEL(1)
      ATOMIPTR(INSERTPOS+1)=II
C     ------------------------------------------------------------------
C     Both the atoms I and M are the metal center
C     ------------------------------------------------------------------
      IF (INSERTPOS.EQ.1) ATOMMPTR(INSERTPOS)=0
      II=ATOMMPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XPCSMX) XPCSMX=II
      ATOMMLST(II)=SEL(1)
      ATOMMPTR(INSERTPOS+1)=II

      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_PCS',
     &                'More than 1 atom in SEL for atom J. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
         CALL DSPERR('XDIPO_PCS','Error with atom selection')
         NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL)
      IF (INSERTPOS.EQ.1) ATOMJPTR(INSERTPOS)=0
      II=ATOMJPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XPCSMX) XPCSMX=II
      ATOMJLST(II)=SEL(1)
      ATOMJPTR(INSERTPOS+1)=II

      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_PCS',
     &                'More than 1 atom in SEL for atom K. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XDIPO_PCS','Error with atom selection')
          NFLAG=1
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (INSERTPOS.EQ.1) ATOMKPTR(INSERTPOS)=0
      II=ATOMKPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XPCSMX) XPCSMX=II
      ATOMKLST(II)=SEL(1)
      ATOMKPTR(INSERTPOS+1)=II

      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_PCS',
     &                'More than 1 atom in SEL for atom L. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XDIPO_PCS','Error with atom selection')
          NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL)
      IF (INSERTPOS.EQ.1) ATOMLPTR(INSERTPOS)=0
      II=ATOMLPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XPCSMX) XPCSMX=II
      ATOMLLST(II)=SEL(1)
      ATOMLPTR(INSERTPOS+1)=II
C     ------------------------------------------------------------------
C     Next selection defines the nucleus whose PCS was measured
C     ------------------------------------------------------------------
      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_PCS',
     &                'More than 1 atom in SEL for atom N. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XDIPO_PCS','Error with atom selection')
          NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL)
      IF (INSERTPOS.EQ.1) ATOMNPTR(INSERTPOS)=0
      II=ATOMNPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XPCSMX) XPCSMX=II
      ATOMNLST(II)=SEL(1)
      ATOMNPTR(INSERTPOS+1)=II
C     ------------------------------------------------------------------
C     Reading in the observed pseudocontact shift
C     ------------------------------------------------------------------
      CALL NEXTF('Observed PCS =', XPCSO)
      CALL NEXTF('Error in PCS =', XPCSE)

      XPCSOBS(INSERTPOS)=XPCSO
      XPCSERR(INSERTPOS)=XPCSE
C     ------------------------------------------------------------------
C     Check for error atom selection. If there is one then reset the
C     counter for restraint
C     ------------------------------------------------------------------
      IF (NFLAG.EQ.1) THEN
          PCSNUM=PCSNUM-1
      END IF

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE XPCSINIT
C     ------------------------------------------------------------------
C     Initializes pseudocontact shift stuff
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "xdipo_pcs.inc"

      CALL XPCSDEFAULTS
      CALL ALLOCXPCS(0,MAXXPCS)

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE PRINTXPCS (CUTOFF,XPCSCALC,XPCSOBS,XPCSERR,
     &                      ATOMIPTR,ATOMJPTR,ATOMKPTR,
     &                      ATOMLPTR,ATOMMPTR,ATOMNPTR,
     &                      ATOMILST,ATOMJLST,ATOMKLST,
     &                      ATOMLLST,ATOMMLST,ATOMNLST)
C     ------------------------------------------------------------------
C     Prints pseudocontact shifts with deviation larger than cutoff,
C     calculates RMSD and puts it into $RESULT
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "xdipo_pcs.inc"
      INCLUDE  "comand.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "numbers.inc"

      INTEGER FVIOLS
      COMMON  /VIOL/FVIOLS
      DOUBLE PRECISION CUTOFF,XPCSCALC(*),XPCSOBS(*),XPCSERR(*)
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*)
      INTEGER ATOMLLST(*),ATOMMLST(*),ATOMNLST(*)
      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*)
      INTEGER ATOMLPTR(*),ATOMMPTR(*),ATOMNPTR(*)
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      DOUBLE PRECISION CALCXPCS,OBSXPCS,DELTAXPCS,DELTA,DP
      INTEGER COUNT,CLASS,I,J,K,L,M,N,NUM,II
      DOUBLE PRECISION RMS,VIOLS,ERRXPCS,XPCSENERGY
      DOUBLE COMPLEX DUMMY2
      LOGICAL PRINTTHISCLASS

      RMS=ZERO
      VIOLS=ZERO
      NUM=0
C     ------------------------------------------------------------------
C     Make sure that the array of calculated PCS is up to date
C     ------------------------------------------------------------------
      CALL EXPCS(XPCSENERGY,'ANALYZE')
      WRITE (PUNIT,'(A)') 'The following pseudocontact shifts have'
      WRITE (PUNIT,'(A)') 'deviation larger than or'
      WRITE (PUNIT,'(A)') 'equal to the cutoff:'
C     ------------------------------------------------------------------
C     Write out first class heading
C     ------------------------------------------------------------------
      CLASS=1
      PRINTTHISCLASS=PRINTCLASS(CLASS)
      IF (PRINTTHISCLASS) THEN
          WRITE (PUNIT,'(A,A)') 'Class ',XPCSCLASSNAMES(1)
      END IF
C     ------------------------------------------------------------------
C     For every pseudocontact shift entry...
C     ------------------------------------------------------------------
      COUNT = 0
      DO WHILE (COUNT.LT.PCSNUM)
         COUNT=COUNT+1
C     ------------------------------------------------------------------
C     Is this the start of a new class?
C     ------------------------------------------------------------------
         IF (XPCSASSNDX(CLASS).LT.COUNT) THEN
             CLASS=CLASS+1
            PRINTTHISCLASS=PRINTCLASS(CLASS)
            IF (XPCSASSNDX(CLASS).EQ.XPCSASSNDX(CLASS-1)) THEN
                COUNT=COUNT-1
            END IF
            IF (PRINTTHISCLASS) THEN
                WRITE (PUNIT,'(A,A)') 'Class ',XPCSCLASSNAMES(CLASS)
            END IF
         END IF
C     ------------------------------------------------------------------
C     If this assignment is in a class to be printed
C     and make sure there is an entry for that class
C     ------------------------------------------------------------------
         IF ((PRINTTHISCLASS).AND.
     &       (XPCSASSNDX(CLASS).NE.XPCSASSNDX(CLASS-1))) THEN
C     ------------------------------------------------------------------
C     Always update RMSD
C     ------------------------------------------------------------------
            CALCXPCS=XPCSCALC(COUNT)
            OBSXPCS=XPCSOBS(COUNT)
            ERRXPCS=XPCSERR(COUNT)
            DP=(CALCXPCS-OBSXPCS)

            IF ((DP.LT.0.000).AND.(ABS(DP).GT.ERRXPCS)) THEN
                 DELTAXPCS=DP+ERRXPCS
            ELSE IF ((DP.GT.0.000).AND.(DP.GT.ERRXPCS)) THEN
                 DELTAXPCS=DP-ERRXPCS
            ELSE
                 DELTAXPCS=0.0
            END IF

            RMS=RMS+DELTAXPCS**2
            NUM=NUM+1
C     ------------------------------------------------------------------
C     Print out deviations larger than cutoff
C     and update number of violations
C     ------------------------------------------------------------------
            IF (ABS(DELTAXPCS).GE.CUTOFF) THEN
                I=ATOMILST(ATOMIPTR(COUNT)+1)
                J=ATOMJLST(ATOMJPTR(COUNT)+1)
                K=ATOMKLST(ATOMKPTR(COUNT)+1)
                L=ATOMLLST(ATOMLPTR(COUNT)+1)
                WRITE(PUNIT,'(A,A)') '===============================',
     &                               '==============================='
                WRITE(PUNIT,'(A)') ' Set-M-atoms'
                DO II=ATOMMPTR(COUNT)+1,ATOMMPTR(COUNT+1)
                   M=ATOMMLST(II)
                   WRITE(PUNIT,'(9X,4(1X,A))') SEGID(M),RESID(M),RES(M),
     &                                         TYPE(M)
                END DO
                WRITE(PUNIT,'(A)') ' Set-N-atoms'
                DO II=ATOMNPTR(COUNT)+1,ATOMNPTR(COUNT+1)
                   N=ATOMNLST(II)
                   WRITE(PUNIT,'(9X,4(1X,A))') SEGID(N),RESID(N),RES(N),
     &                                         TYPE(N)
                END DO

                WRITE(PUNIT,'(2(2X,A,1X,F8.3))')
     &                'Calc: ',CALCXPCS,'Obs: ',OBSXPCS
                WRITE(PUNIT,'(2X,A,1X,F8.3,2X,A,1X,F8.3)')
     &                'Error: ',ERRXPCS,'Delta: ',DELTAXPCS
                VIOLS = VIOLS + ONE
            END IF
         END IF
      END DO

      IF (NUM.GT.0) THEN
         RMS=SQRT(RMS/NUM)
      ELSE
         RMS=0.0
      END IF

      FVIOLS=VIOLS
      PRINT *,'FVIOLS =',FVIOLS

      CALL DECLAR('RESULT','DP',' ',DUMMY2,RMS)
      CALL DECLAR('VIOLATIONS','DP',' ',DUMMY2,VIOLS)

      RETURN
      END

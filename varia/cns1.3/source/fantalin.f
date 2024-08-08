C     ------------------------------------------------------------------
      SUBROUTINE FANTALINEXE(XPX,XPY,XPZ,XPMX,XPMY,XPMZ,NNAT,XOBS,
     &                       XTOLPROT,METHOD,NOMEFILE1,ERWPROT,ERRORE,
     &                       LOCFLG)
C     ------------------------------------------------------------------
C     Handles the stuff for the calculation of tensor parameters
C
C     METHOD switches between PCS(1) and RDC(0)
C     CX,CY,CZ(STRUCTURE,ATOMNUMBER) keep coordinates of nuclei
C     FX,FY,FZ(STRUCTURE,ATOMNUMBER) keep coordinates of metal
C     NAT is a vector with the number of atoms
C     OBS is a vector with the observed pseudocontact shifts
C     TOLPROT is a vector with the errors on PCS
C     NOMEFILE1 is the name of the file to save deviations
C     ERWPROT is a vector with the weights
C     ERRORE switches the USE(1) of weights
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'sup2.inc'
      INCLUDE 'fantaxplor.inc'

      INTEGER NNAT,METHOD
      DOUBLE PRECISION XPX(1000),XPY(1000),XPZ(1000),
     &                 XOBS(1000),XTOLPROT(1000),ERWPROT(1000)
      DOUBLE PRECISION XPMX(1000),XPMY(1000),XPMZ(1000)
      CHARACTER NOMEFILE1*132 
      INTEGER ERRORE
      LOGICAL LOCFLG

      PRINT*,'                    '
      PRINT*,'===---FANTALIN---==='
C     ------------------------------------------------------------------
C     IHP is the number of PCS
C     ------------------------------------------------------------------
      MDIPO=METHOD
      PRINT*,'MDIPO is:',MDIPO
      PRINT*,'ERRORE flag is:',ERRORE
      NAT=NNAT
      IHP=NNAT
      NHP=NNAT
      PRINT*,'IHP is:',IHP
C     ------------------------------------------------------------------
C     NFE is the number of paramagnetic centers
C     ------------------------------------------------------------------
      NFE=1
      NIA=0
C     ------------------------------------------------------------------
C     Reference system not active
C     ------------------------------------------------------------------
      NSYSTEM=0
C     ------------------------------------------------------------------
C     Calculation is done on a single structure
C     ------------------------------------------------------------------
      NSTR=1
C     ------------------------------------------------------------------
C     Errors are set to 0
C     ------------------------------------------------------------------
      PERC=0
C     ------------------------------------------------------------------
C     Grid setting
C     ------------------------------------------------------------------
      NGRID=1
C     ------------------------------------------------------------------
C     Weights and multiplicity are set to 1
C     ------------------------------------------------------------------
      DO NIA=1,NAT
         MLPROT(NIA)=1
         IF(ERRORE.EQ.0) THEN
            WPROT(NIA)=1
         ELSE
            WPROT(NIA)=ERWPROT(NIA)
         END IF
      END DO

      OPEN(322,FILE='FANTALIN.LOG')
C     ------------------------------------------------------------------
C     Vectors are passed from Xplor-NIH core
C     ------------------------------------------------------------------
      DO NIA=1,NAT
         OBS(NIA)=XOBS(NIA)
C     ------------------------------------------------------------------
C     The use of tolerance in fitting depends on this flag
C     ------------------------------------------------------------------
         IF (LOCFLG) THEN
             TOLPROT(NIA)=XTOLPROT(NIA)
         ELSE
             TOLPROT(NIA)=0
         END IF
         CX(NSTR,NIA)=XPX(NIA)
         CY(NSTR,NIA)=XPY(NIA)
         CZ(NSTR,NIA)=XPZ(NIA)
         WRITE (322,123) NIA,OBS(NIA),CX(NSTR,NIA),CY(NSTR,NIA),
     &                   CZ(NSTR,NIA),TOLPROT(NIA),WPROT(NIA)
         WRITE (322,55) NRESIDFA(NIA),NRESFA(NIA),NTYPEFA(NIA)
123      FORMAT (2X,I3,2X,F8.3,2X,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3,F8.3)
 55      FORMAT (1X,'RESID',1X,A4,1X,'RES',1X,A4,1X,'TYPE',1X,A4)
C     ------------------------------------------------------------------
C     Only one paramagnetic center is passed
C     ------------------------------------------------------------------
         IF (METHOD.EQ.1) THEN
             FX(NSTR,1)=XPMX(NIA)
             FY(NSTR,1)=XPMY(NIA)
             FZ(NSTR,1)=XPMZ(NIA)
         ELSE
             FX(NSTR,NIA)=XPMX(NIA)
             FY(NSTR,NIA)=XPMY(NIA)
             FZ(NSTR,NIA)=XPMZ(NIA)
         END IF
         WRITE (322,124) FX(NSTR,1),FY(NSTR,1),FZ(NSTR,1)
124      FORMAT (2X,F8.3,2X,F8.3,2X,F8.3)
      END DO
      CLOSE (322)
      FILENAME3=NOMEFILE1
      PRINT*,'NOMEFILE1 is:',NOMEFILE1
      OPEN(744,FILE=FILENAME3)

      CALL SIMPLEXRUN()
      CALL SCAMB(CHIPARA,CHIPERP)

      RETURN    
      END  
C     ------------------------------------------------------------------
      SUBROUTINE SCAMB(A,B)
C     ------------------------------------------------------------------
C     The value of A2 is always negative
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      DOUBLE PRECISION A,B,APARA,APERP
      COMMON /FRUN/APERP,APARA

      APARA=A
      APERP=-B

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE SIMPLEXRUN()
C     ------------------------------------------------------------------
C     Calculates the susceptibility tensor of paramagnetic system  
C     which mimimize the experimental-to-calculated squared shift  
C     difference
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'sup2.inc'

      DIMENSION SIMP(MP,NP),EXTR(NP),VAL(MP)  
      DIMENSION ZRR(3,3)  

      OLDRESID=1.E+9  
      RESID=2.E+9  
      IVIOLATION=2000000  
      IOLDVIO=1000000  
C     ------------------------------------------------------------------
C     Initial conditions grid
C     ------------------------------------------------------------------
      DELTA=3.14159/DBLE(NGRID+1)  
      DO 1000 L=1,NGRID  
         DO 999 NK=1,NFE  
            SIMP(1,(NK-1)*MPAR+1)=1000.*(1.-2*RAND())  
            SIMP(1,(NK-1)*MPAR+2)=1000.*(1.-2*RAND())  
            SIMP(1,(NK-1)*MPAR+3)=1000.*(1.-2*RAND())  
            SIMP(1,(NK-1)*MPAR+4)=1000.*(1.-2*RAND())  
            SIMP(1,(NK-1)*MPAR+5)=1000.*(1.-2*RAND())  
999      CONTINUE  
         DO 10 I=1,MP  
            DO 10 J=1,NP  
               SIMP(I,J)=SIMP(1,J)  
               IF (I.EQ.1) THEN  
                   SIMP(1,J)=SIMP(1,J)*0.8  
               ENDIF  
               IF (I.EQ.J+1) THEN  
                   SIMP(I,J)=SIMP(I,J)/0.8*1.2  
               ENDIF  
10       CONTINUE  
         DO 11 I=1,MP  
            DO 12 J=1,NP  
12             EXTR(J)=SIMP(I,J)  
               VAL(I)=CASH(EXTR)  
11       CONTINUE  
         CALL SIMPCALC(SIMP,VAL)  
C     ------------------------------------------------------------------
C     Criterion of minimization of optimal solution
C     ------------------------------------------------------------------
         IF (OLDRESID.GE.RESID) THEN  
             OLDRESID=RESID  
             IOLDVIO=IVIOLATION  
             DO 13 I=1,NFE  
                OPTPHI(I)=PHI(I)  
                OPTTETA(I)=TETA(I)  
                OPTOMEGA(I)=OMEGA(I)  
                OPTA1(I)=A1DIP(I)  
                OPTA2(I)=A2DIP(I)  
13           CONTINUE  
         ENDIF  
1000  CONTINUE

      DO 998 NK=1,NFE  
         EXTR((NK-1)*MPAR+1)=OPTPHI(NK)  
         EXTR((NK-1)*MPAR+2)=OPTTETA(NK)  
         EXTR((NK-1)*MPAR+3)=OPTOMEGA(NK)  
         EXTR((NK-1)*MPAR+4)=OPTA1(NK)  
         EXTR((NK-1)*MPAR+5)=OPTA2(NK)  
998   CONTINUE 

      PRINT*,' '  
      WRITE(*,FMT='(1X,A)') ' *****    BEST SOLUTION    *****'  
      PRINT*,' '  
      WRITE(*,20) IOLDVIO,OLDRESID  
      WRITE(*,21) (OPTPHI(I),OPTTETA(I),OPTOMEGA(I),I=1,NFE)  
      WRITE(*,22) (OPTA1(I),OPTA2(I),I=1,NFE)  

      CALL DIAG(ZRR,OPTPHI(1),OPTTETA(1),OPTOMEGA(1),OPTA1(1),  
     &          OPTA2(1),CHIPARA,CHIPERP) 
      
      WRITE(*,*) 'A1 = ',CHIPARA,', ',CHIPARA*12*3.1416/1E4,'E-32'
      WRITE(*,*) 'A2 = ',-CHIPERP,', ',-CHIPERP*12*3.1416/1E4,'E-32'
      WRITE(744,*) 'A1 = ',CHIPARA,', ',CHIPARA*12*3.1416/1E4,'E-32'
      WRITE(744,*) 'A2 = ',-CHIPERP,', ',-CHIPERP*12*3.1416/1E4,'E-32'

      X4=ZRR(1,1)  
      Y4=ZRR(1,2)  
      Z4=ZRR(1,3)  
      X4=ZRR(1,1)  
      Y4=ZRR(2,1)  
      Z4=ZRR(3,1)  
      X4=ZRR(2,1)  
      Y4=ZRR(2,2)  
      Z4=ZRR(2,3)  
      X4=ZRR(1,2)  
      Y4=ZRR(2,2)  
      Z4=ZRR(3,2)  
      X4=ZRR(3,1)  
      Y4=ZRR(3,2)  
      Z4=ZRR(3,3)  
      X4=ZRR(1,3)  
      Y4=ZRR(2,3)  
      Z4=ZRR(3,3)  
      DO 777 I=1,NFE  
         P=OPTPHI(I)  
         T=OPTTETA(I)  
         O=OPTOMEGA(I)  
         ASXX(I,1)=COS(P)*COS(O)  
         ASXX(I,2)=SIN(P)*COS(O)  
         ASXX(I,3)=SIN(O)  
         ASYY(I,1)=-COS(T)*SIN(P)-SIN(O)*COS(P)*SIN(T)  
         ASYY(I,2)=COS(T)*COS(P)-SIN(O)*SIN(P)*SIN(T)  
         ASYY(I,3)=SIN(T)*COS(O)  
         ASZZ(I,1)=SIN(T)*SIN(P)-SIN(O)*COS(P)*COS(T)  
         ASZZ(I,2)=-SIN(T)*COS(P)-SIN(O)*SIN(P)*COS(T)  
         ASZZ(I,3)=COS(T)*COS(O)  
777   CONTINUE  

      CALL PRRES  
  
20    FORMAT (1X,' VIOLATIONS ',I6,1X,' RESIDUALS ',F10.3)  
21    FORMAT (1X,'C1= ',F15.3,1X,'C2= ',F15.3,1X,'C3= ',F15.3)  
22    FORMAT (1X,'C4= ',F15.3,1X,'C5= ',F15.3)  
  
      ERR=CASH(EXTR)  
      CALL DISPSHIFT

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE SIMPCALC(P,Y)
C     ------------------------------------------------------------------
C     Simplex calculation
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'sup2.inc'
      INCLUDE 'fantaxplor.inc'  

      PARAMETER (AL=1.,BE=.5,GA=2.)  
      DIMENSION P(MP,NP),Y(MP),PR(MP),PT(MP),PE(MP),  
     &          PC(MP),PK(MP),PBAR(MP),OLDP(NP),EXTR(MP)  
C     ------------------------------------------------------------------
C     Algorithm initialization
C     ------------------------------------------------------------------
      ITER=0  
      RNDOM=0  
1     PERM=1  
      DO 1000 WHILE (PERM.GT.0)  
         PERM=0  
         DO 1100 I=1,NP  
            IF (Y(I).GT.Y(I+1)) THEN  
                APP=Y(I+1)  
                Y(I+1)=Y(I)  
                Y(I)=APP  
                DO 1200 J=1,NP  
                   APP=P(I+1,J)  
                   P(I+1,J)=P(I,J)  
                   P(I,J)=APP 
1200            CONTINUE  
                PERM=PERM+1  
            ENDIF  
1100     CONTINUE    
1000  CONTINUE  

      AVERAGE=0  
      ASTDEV=0  
      DO 100 I=1,MP  
100   AVERAGE=AVERAGE+Y(I)  
      AVERAGE=AVERAGE/MP  
      DO 101 I=1,MP  
101   ASTDEV=ASTDEV+(Y(I)-AVERAGE)**2.  
      ERR=DSQRT(ASTDEV/MP)  
      RESID=Y(1)    
      IF (ITER.GT.39999.OR.ERR.LT.1.D-9) THEN  
          WRITE(*,200) ITER,IVIOLATION,RESID  
          WRITE (*,201)   
     &          (P(1,(I-1)*MPAR+1),  
     &           P(1,(I-1)*MPAR+2),  
     &           P(1,(I-1)*MPAR+3),  
     &           P(1,(I-1)*MPAR+4),P(1,(I-1)*MPAR+5),I=1,NFE)  
200       FORMAT (1X,'ITERS ',I6,1X,' VIOLATIONS ',I6,  
     &            1X,'RESIDUALS= ',F10.3)  
201       FORMAT (1X,'C1= ',F9.3,1X,'C2= ',F9.3,1X,'C3= ',F9.3,  
     &            1X,'C4= ',F9.3,1X,'C5= ',F9.3)  
          DO 997 NK=1,NFE  
             PHI(NK)=P(1,(NK-1)*MPAR+1)  
             TETA(NK)=P(1,(NK-1)*MPAR+2)  
             OMEGA(NK)=P(1,(NK-1)*MPAR+3)  
             A1DIP(NK)=P(1,(NK-1)*MPAR+4)  
             A2DIP(NK)=P(1,(NK-1)*MPAR+5)  
997       CONTINUE  
          RETURN  
      ENDIF  

      ITER=ITER+1  
      DO 12 J=1,NP  
12    PBAR(J)=0  
      DO 13 I=1,NP  
         DO 13 J=1,NP  
13    PBAR(J)=PBAR(J)+P(I,J)  
      DO 14 J=1,NP  
         PBAR(J)=PBAR(J)/NP  
         PR(J)=(1.+AL)*PBAR(J)-AL*P(MP,J)  
14    CONTINUE  
      YR=CASH(PR)  
      DO 15 J=1,NP  
         PK(J)=PR(J)  
15    CONTINUE  
      YK=YR  
      IF (YR.LT.Y(NP)) THEN  
          IF (YR.LT.Y(1)) THEN  
              DO 16 J=1,NP  
                 PE(J)=GA*PR(J)+(1.-GA)*PBAR(J)  
16            CONTINUE  
              YE=CASH(PE)  
              ITER=ITER+1  
              IF (YE.LT.Y(1)) THEN  
                  DO 17 J=1,NP  
                     PK(J)=PE(J)  
17                CONTINUE  
                  YK=YE 
              ENDIF  
          ENDIF  
      ELSE   
          DO 18 J=1,NP  
             PT(J)=P(MP,J)  
18        CONTINUE  
          YT=Y(MP)  
          IF (YR.LT.YT) THEN  
              DO 19 J=1,NP  
                 PT(J)=PR(J)  
19            CONTINUE  
              YT=YR  
          ENDIF  
          DO 20 J=1,NP  
             PC(J)=BE*PT(J)+(1.-BE)*PBAR(J)  
20        CONTINUE  
          YC=CASH(PC)  
          ITER=ITER+1  
          IF (YC.LT.Y(NP)) THEN  
              DO 21 J=1,NP  
                 PK(J)=PC(J)  
21            CONTINUE  
              YK=YC  
          ELSE  
              DO 22 I=1,NP  
                 OLDP(I)=P(2,I)  
22            CONTINUE  
              DO 24 I=2,NP   
                 DO 23 J=1,NP  
                    P(I,J)=.5*(P(1,J)+P(I,J))  
23               CONTINUE  
                 PR(I)=P(I,I)    
                 Y(I)=CASH(PR)  
24            CONTINUE  
              ITER=ITER+NP-1  
              DO 25 I=1,NP  
                 OLDP(I)=ABS(OLDP(I)-P(2,I))  
25            CONTINUE  
              OLDPOPPA=0   
              DO 251 I=1,NP-1  
                 OLDPOPPA=OLDPOPPA+OLDP(I)  
251           CONTINUE  
              IF (ABS(OLDPOPPA).LT.1.D-8.AND.RNDOM.LT.10) THEN  
                  RNDOM=RNDOM+1  
                  WRITE(*,FMT='(BN,A,F16.6)') 'RANDOM  ',Y(1)  
                  DO 2001 IL=1,MP  
                     DO 2002 IM=1,NP  
                        RMR=1+RAND()  
                        P(IL,IM)=P(IL,IM)*RMR  
2002                 CONTINUE  
2001              CONTINUE     
                  DO 2011 IL=1,MP  
                     DO 2012 JM=1,NP  
2012                 EXTR(JM)=P(IL,JM)  
                     Y(IL)=CASH(EXTR)  
2011              CONTINUE  
                  GOTO 1   
              ENDIF   
              IF (ABS(OLDPOPPA).LT.1.D-9) THEN  
                  ITER=50000  
              ENDIF  
              DO 26 J=1,NP  
                 PK(J)=0.5*(P(1,J)+P(MP,J))  
26            CONTINUE  
              YK=CASH(PK)  
              ITER=ITER+1  
          ENDIF  
      ENDIF  
      DO 27 J=1,NP  
         P(MP,J)=PK(J)  
27    CONTINUE  
      Y(MP)=YK   
      GOTO 1  
      END  
C     ------------------------------------------------------------------
      FUNCTION CASH(VETT)  
C     ------------------------------------------------------------------
C     Calculation of deviations and violations
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'sup2.inc'  

      DIMENSION VETT(NP)
      INTEGER M1  

      IVIOLATION=0  
      TMP1=0  
      I=1  
      DO WHILE (I.LE.NHP*NSTR)  
         SHIFT(I)=0.0  
         I=I+1  
      ENDDO  

      DO 1 N=1,NSTR  
         DO 2 M=1,NFE  
            P=VETT((M-1)*MPAR+1)  
            T=VETT((M-1)*MPAR+2)  
            O=VETT((M-1)*MPAR+3)  
            A1D=VETT((M-1)*MPAR+4)  
            A2D=VETT((M-1)*MPAR+5)  
            TMP2=0.0  
            I=1  
            IHP=(N-1)*NHP 
            DO 10 WHILE (I.LE.NHP)  
               NPROTML=MLPROT(IHP+I)
               IF (MDIPO.EQ.0) THEN
                   M1=I
               ELSE
                   M1=M
               END IF
               IF (NPROTML.GT.1) THEN  
                   COSTMULT=1/DBLE(NPROTML)  
                   APP=0.0  
                   DO 11 L=1,NPROTML
                      IF (PRIMO.EQ.0) THEN  
                          X4=CX(N,I)-FX(N,M1)
                          Y4=CY(N,I)-FY(N,M1)  
                          Z4=CZ(N,I)-FZ(N,M1)  
                          R(N,I)=SQRT(X4**2+Y4**2+Z4**2)  
                          TETA1(N,I)=ACOS(Z4/R(N,I))  
                          PHI1(N,I)=ACOS(X4/(R(N,I)*SIN(TETA1(N,I))))  
                          IF (Y4.LT.0) PHI1(N,I)=2*3.1416-PHI1(N,I)  
                      ENDIF  
                      G1=(1.-3*COS(TETA1(N,I))**2)/R(N,I)**3  
                      G2=(SIN(TETA1(N,I))**2*COS(2*PHI1(N,I)))/
     &                    R(N,I)**3  
                      G3=2*(SIN(TETA1(N,I))**2*SIN(2*PHI1(N,I)))/
     &                    R(N,I)**3  
                      G4=2*(SIN(2*TETA1(N,I))*COS(PHI1(N,I)))/
     &                    R(N,I)**3  
                      G5=2*(SIN(2*TETA1(N,I))*SIN(PHI1(N,I)))/
     &                    R(N,I)**3  
                      APP=APP+3*(P*G1+T*G2+O*G3+A1D*G4+A2D*G5)/2.  
                      I=I+1  
11                 CONTINUE  
                   DO 12 L=1,NPROTML  
                      SHIFT(IHP+I-L)=SHIFT(IHP+I-L)+APP*COSTMULT  
12                 CONTINUE  
               ELSE  
                   IF (PRIMO.EQ.0) THEN  
                       X4=CX(N,I)-FX(N,M1)  
                       Y4=CY(N,I)-FY(N,M1)  
                       Z4=CZ(N,I)-FZ(N,M1)  
                       R(N,I)=SQRT(X4**2+Y4**2+Z4**2)  
                      TETA1(N,I)=ACOS(Z4/R(N,I))  
                      PHI1(N,I)=ACOS(X4/(R(N,I)*SIN(TETA1(N,I))))  
                      IF (Y4.LT.0) PHI1(N,I)=2*3.1416-PHI1(N,I)  
                   ENDIF  
                   G1=(1.-3*COS(TETA1(N,I))**2)/R(N,I)**3  
                   G2=(SIN(TETA1(N,I))**2*COS(2*PHI1(N,I)))/R(N,I)**3  
                   G3=2*(SIN(TETA1(N,I))**2*SIN(2*PHI1(N,I)))/R(N,I)**3
                   G4=2*(SIN(2*TETA1(N,I))*COS(PHI1(N,I)))/R(N,I)**3  
                   G5=2*(SIN(2*TETA1(N,I))*SIN(PHI1(N,I)))/R(N,I)**3  
                   SHIFT(IHP+I)=SHIFT(IHP+I)+
     &                          3*(P*G1+T*G2+O*G3+A1D*G4+A2D*G5)/2.   
                   I=I+1  
               ENDIF  
10          CONTINUE  
2        CONTINUE  
         DO 3 I=1,NHP  
            TMP2=ABS(SHIFT(IHP+I)-OBS(IHP+I))-TOLPROT(IHP+I)  
            IF (TMP2.GT.0.0) THEN  
                IVIOLATION=IVIOLATION+1  
                TMP1=TMP1+TMP2**2*WPROT(IHP+I)/DBLE(MLPROT(IHP+I))  
            ENDIF  
3        CONTINUE   
1     CONTINUE  
      IF (PRIMO.EQ.0) PRIMO=1  
  
      CASH=TMP1  

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE DISPSHIFT 
C     ------------------------------------------------------------------
C     Display of violations
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'sup2.inc'
      INCLUDE 'fantaxplor.inc'  

      CHARACTER FILESTR*21,SS*40  
      DIMENSION VV(MAXSTR,MAXOS),VMAX(MAXOS),VMIN(MAXOS)  
CCC      DIMENSION NNMAX(MAXOS),NNMIN(MAXOS)  

      RM=0  
      TMP1=0  
      SMED=0  
      STDEV=0  
      DO 2 I=1,NSTR  
         DO 2 J=1,NHP  
            TMP2=ABS(OBS(J+(I-1)*NHP)-SHIFT(J+(I-1)*NHP))-
     &           TOLPROT(J+(I-1)*NHP)  
            IF (TMP2.GT.0.0) THEN  
                IF (SHIFT(J+(I-1)*NHP).GT.OBS(J+(I-1)*NHP)) THEN  
                    VV(I,J)=TMP2  
                ELSE  
                    VV(I,J)=-TMP2  
                ENDIF  
            ELSE  
                VV(I,J)=0.  
            ENDIF  
2     CONTINUE  
      DO 44 I=1,NHP  
         VMAX(I)=0.  
         VMIN(I)=1000.  
44    CONTINUE  
      DO 4 J=1,NHP  
         DO 4 I=1,NSTR  
            IF (ABS(VV(I,J)).GT.ABS(VMAX(J))) THEN   
                VMAX(J)=VV(I,J)  
            ENDIF  
            IF ((ABS(VV(I,J)).LT.ABS(VMIN(J))).AND.(VV(I,J).NE.0.)) THEN
                 VMIN(J)=VV(I,J)   
            ENDIF  
4     CONTINUE     
      DO 7 J=1,NHP  
         DO 8 I=1,NSTR  
              IF (ABS(VV(I,J)).GT.0.) THEN  
                  SS(I:I+1)="*"  
              ELSE  
                  SS(I:I+1)=" "  
              ENDIF  
8        CONTINUE  
         IF (VMIN(J).EQ.1000.) THEN  
             VMIN(J)=0.  
         ENDIF  
  
7     CONTINUE      
  
      WRITE(744,'(1X)')   
      WRITE(744,'(1X)')   
      WRITE(744,'(1X)')   
      WRITE(744,'(1X,A)')
     &      '    PROTONS              OBS.    CALC.    ERR.^2'

      FILESTR='USELESS.FILE'
      NIA=0
      REWIND(1)  
      DO 3 J=1,NHP*NSTR  
         TMP2=ABS(OBS(J)-SHIFT(J))-TOLPROT(J)
         NIA=NIA+1
         IF (TMP2.GT.0.0) THEN  
             IVIOLATION=IVIOLATION+1  
             RM=RM+1  
             SMED=SMED+TMP2**2*WPROT(J)/MLPROT(J)  
             TMP1=TMP2**2*WPROT(J)/MLPROT(J)  
             WRITE (744,'(1X,A4,1X,A4,1X,A4,1X,F7.3,1X,F7.3,1X,F10.5)')
     &              NRESIDFA(NIA),NRESFA(NIA),NTYPEFA(NIA),
     &              OBS(J),SHIFT(J),TMP1  
15           FORMAT(1X,A,1X,F7.3,1X,F7.3,1X,F7.3)  
         ELSE  
             WRITE (744,'(1X,A4,1X,A4,1X,A4,1X,F7.3,1X,F7.3)')  
     &              NRESIDFA(NIA),NRESFA(NIA),NTYPEFA(NIA),
     &              OBS(J),SHIFT(J)  
         ENDIF  
3     CONTINUE  
      SMED=SMED/RM  

      DO 41 J=1,NHP*NSTR  
         TMP2=ABS(OBS(J)-SHIFT(J))-TOLPROT(J)  
         IF (TMP2.GT.0.0) THEN  
             STDEV=STDEV+(TMP2**2-MED)**2*WPROT(J)/MLPROT(J)  
         ENDIF  
41    CONTINUE  
      STDEV=STDEV/RM  

      PRINT*,' '  
      WRITE(744,'(A)') ' '   
      WRITE(*,11) 'Mean value of error =',SMED  
11    FORMAT(1X,A,1X,F10.6)  
      WRITE(744,12) 'Mean value of error =',SMED  
12    FORMAT(1X,A,1X,F10.6)  
      WRITE(*,13) 'Standard deviation =',STDEV  
13    FORMAT(1X,A,1X,F10.6)      
      WRITE(744,14) 'Standard deviation =',STDEV  
14    FORMAT(1X,A,1X,F10.6)           
      CLOSE(744)  

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE PRRES  
C     ------------------------------------------------------------------
C     SAXX, SAYY, SAZZ : Direction cosines of rotation matrices
C     SAXSYS, SAYSYS, SAZSYS : Direction cosines of reference systems
C 
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'sup2.inc'
      INCLUDE 'fantaxplor.inc'  

      DIMENSION SAXX(3),SAYY(3),SAZZ(3),  
     &          APPX(3),APPY(3),  
     &          SAXSYS(3),SAYSYS(3),SAZSYS(3)  
      CHARACTER NAMATOM*4,NAMRESS*3  
      PI=3.1415927  

99    I=0      
      DO 10 N=1,NFE  
         A1=OPTA1(N)  
         A2=OPTA2(N)  
         IF (NSYSTEM.GT.1) THEN  
             DO 1 I=1,3  
                SAXSYS(I)=AXSYS(N,I)  
                SAYSYS(I)=AYSYS(N,I)  
                SAZSYS(I)=AZSYS(N,I)  
1            CONTINUE  
         ELSE  
             DO 2 I=1,3  
                SAXSYS(I)=AXSYS(1,I)  
                SAYSYS(I)=AYSYS(1,I)  
                SAZSYS(I)=AZSYS(1,I)  
2            CONTINUE  
         ENDIF  
         DO 3 I=1,3  
            SAXX(I)=ASXX(N,I)  
            SAYY(I)=ASYY(N,I)  
            SAZZ(I)=ASZZ(N,I)  
3        CONTINUE  
         SDZX=SAXX(2)*SAYY(3)-SAXX(3)*SAYY(2)  
         IF (SDZX/SAZZ(1).GT.0.0) THEN  
C     ------------------------------------------------------------------
C     Right-handed
C     ------------------------------------------------------------------
             PRINT*,' '  
         ELSE  
C     ------------------------------------------------------------------
C     Left-handed
C     ------------------------------------------------------------------
             PRINT *,' '  
             DO 4 I=1,3  
                APPX(I)=SAXX(I)  
                APPY(I)=SAYY(I)  
4            CONTINUE  
             DO 5 I=1,3  
                SAXX(I)=APPY(I)  
                SAYY(I)=APPX(I)  
5            CONTINUE  
             A2=-A2  
         ENDIF  
         PRINT*,  
     &   '***********************************************************'
         WRITE(744,'(A)') ' '  
         WRITE(744,'(A)')  
     &   '**********************************************************'  
      NUMBERATOM=NUMBERATOM+1  
      NAMATOM=' CEN'  
      NAMRESS='PMC'  
      NORES=NORES+N  
10    CONTINUE  
22    FORMAT (1X,'A1= ',F9.3,2X,A1,F7.3,A,2X,'A2= ',F9.3,2X,A1,F7.3,A)
23    FORMAT (1X,'ARCCOS X^XT = ',F8.3,4X,'ARCCOS Y^YT = ',F8.3,  
     &        4X,'ARCCOS Z^ZT = ',F8.3)  
  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE PANG(N,NR)  
C     ------------------------------------------------------------------
C     Angle permutations
C  
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'sup2.inc'
      INCLUDE 'fantaxplor.inc'  
      PI=3.1415927   

      IF (N.EQ.0) THEN  
C     ------------------------------------------------------------------
C     No permutation  
C     ------------------------------------------------------------------
          IF (ASXX(NR,3).EQ.0.0) THEN  
              O=0.0  
              P=ASIN(ASXX(NR,2))  
              T=ASIN(DSQRT(ASZZ(NR,1)**2+ASZZ(NR,2)**2))  
          ELSE IF (ASXX(NR,3).EQ.1.0) THEN  
              O=PI/2   
              P=ASIN(DSQRT(ASYY(NR,1)**2+ASZZ(NR,1)**2))  
              T=ASIN(DSQRT(ASZZ(NR,1)**2+ASZZ(NR,2)**2))  
          ELSE IF (ASXX(NR,3).EQ.-1.0) THEN  
              O=-PI/2  
              P=ASIN(DSQRT(ASYY(NR,1)**2+ASZZ(NR,1)**2))  
              T=ASIN(DSQRT(ASZZ(NR,1)**2+ASZZ(NR,2)**2))  
          ELSE  
              O=ASIN(ASXX(NR,3))  
              P=ASIN(DSQRT((ASYY(NR,1)**2+ASZZ(NR,1)**2-ASXX(NR,3)**2)/
     &          (1.-ASXX(NR,3)**2)))  
              T=ASIN(DSQRT((ASZZ(NR,1)**2+ASZZ(NR,2)**2-ASXX(NR,3)**2)/
     &          (1.-ASXX(NR,3)**2)))  
          ENDIF  
      ELSE IF (N.EQ.1) THEN  
C     ------------------------------------------------------------------
C     Permutation 1
C     ------------------------------------------------------------------
          IF (ASZZ(NR,3).EQ.0.0) THEN  
              O=0.0  
              P=ASIN(ASZZ(NR,2))  
              T=ASIN(DSQRT(ASYY(NR,1)**2+ASYY(NR,2)**2))  
        ELSE IF (ASZZ(NR,3).EQ.1.0) THEN  
              O=PI/2  
              P=ASIN(DSQRT(ASXX(NR,1)**2+ASYY(NR,1)**2))  
              T=ASIN(DSQRT(ASYY(NR,1)**2+ASYY(NR,2)**2))  
        ELSE IF (ASZZ(NR,3).EQ.-1.0) THEN  
              O=-PI/2  
              P=ASIN(DSQRT(ASXX(NR,1)**2+ASYY(NR,1)**2))  
              T=ASIN(DSQRT(ASYY(NR,1)**2+ASYY(NR,2)**2))  
        ELSE  
              O=ASIN(ASZZ(NR,3))  
              P=DASIN(DSQRT((ASXX(NR,1)**2+ASYY(NR,1)**2-ASZZ(NR,3)**2)/
     &          (1.-ASZZ(NR,3)**2)))  
              T=ASIN(DSQRT((ASYY(NR,1)**2+ASYY(NR,2)**2-ASZZ(NR,3)**2)/
     &          (1.-ASZZ(NR,3)**2)))  
        ENDIF  
      ELSE  
C     ------------------------------------------------------------------
C     Permutation 2
C     ------------------------------------------------------------------
        IF (ASYY(NR,3).EQ.0.0) THEN  
            O=0.0  
            P=ASIN(ASYY(NR,2))  
            T=ASIN(DSQRT(ASXX(NR,1)**2+ASXX(NR,2)**2))  
        ELSE IF (ASYY(NR,3).EQ.1.0) THEN  
            O=PI/2  
            P=ASIN(DSQRT(ASZZ(NR,1)**2+ASXX(NR,1)**2))  
            T=ASIN(DSQRT(ASXX(NR,1)**2+ASXX(NR,2)**2))  
        ELSE IF (ASYY(NR,3).EQ.-1.0) THEN  
            O=-PI/2  
            P=ASIN(DSQRT(ASZZ(NR,1)**2+ASXX(NR,1)**2))  
            T=ASIN(DSQRT(ASXX(NR,1)**2+ASXX(NR,2)**2))  
        ELSE  
            O=ASIN(ASYY(NR,3))  
            P=ASIN(DSQRT((ASZZ(NR,1)**2+ASXX(NR,1)**2-ASYY(NR,3)**2)/  
     &        (1.-ASYY(NR,3)**2)))  
            T=ASIN(DSQRT((ASXX(NR,1)**2+ASXX(NR,2)**2-ASYY(NR,3)**2)/  
     &        (1.-ASYY(NR,3)**2)))  
        ENDIF  
      ENDIF  

      SXX=COS(P)*COS(O)  
      SXY=SIN(P)*COS(O)  
      SXZ=SIN(O)  
      SYX=-COS(T)*SIN(P)-SIN(O)*COS(P)*SIN(T)  
      SYY=COS(T)*COS(P)-SIN(O)*SIN(P)*SIN(T)  
      SYZ=SIN(T)*COS(O)  
      SZX=SIN(T)*SIN(P)-SIN(O)*COS(P)*COS(T)  
      SZY=-SIN(T)*COS(P)-SIN(O)*SIN(P)*COS(T)  
      SZZ=COS(T)*COS(O)  
      VZX=SXY*SYZ-SXZ*SYY  
      VZY=SXZ*SYX-SXX*SYZ  
      VZZ=SXX*SYY-SXY*SYX  
      CDIR=SZX*VZX+SZY*VZY+SZZ*VZZ  
      IF (CDIR.LT.0.0) THEN  
          P=-P  
          PRINT*,'Changing PHI' 
      ENDIF  
      SXX=COS(P)*COS(O)  
      SXY=SIN(P)*COS(O)  
      SXZ=SIN(O)  
      SYX=-COS(T)*SIN(P)-SIN(O)*COS(P)*SIN(T)  
      SYY=COS(T)*COS(P)-SIN(O)*SIN(P)*SIN(T)  
      SYZ=SIN(T)*COS(O)  
      SZX=SIN(T)*SIN(P)-SIN(O)*COS(P)*COS(T)  
      SZY=-SIN(T)*COS(P)-SIN(O)*SIN(P)*COS(T)  
      SZZ=COS(T)*COS(O)  
      VZX=SXY*SYZ-SXZ*SYY  
      VZY=SXZ*SYX-SXX*SYZ  
      VZZ=SXX*SYY-SXY*SYX  
      CDIR=SZX*VZX+SZY*VZY+SZZ*VZZ  
      IF (CDIR.LT.0.0) THEN  
          T=-T  
          PRINT*,'Changing THETA'
      ENDIF  
      WRITE (*,21) P,T,O  
21    FORMAT(1X,'PHI= ',F7.3,4X,'THETA= ',F7.3,4X,'OMEGA= ',F7.3)  
  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE DIAG(ZRR,R1,R2,R3,R4,R5,CHIPARA,CHIPERP)  
C     ------------------------------------------------------------------
C     Handles matrix diagonalization
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)  
      DIMENSION WR(3)  
      DIMENSION ARR(3,3),ARI(3,3)  
      DIMENSION ZRR(3,3),ZRI(3,3)  
      DIMENSION WK1(3),WK2(3),WK3(3)  
C     ------------------------------------------------------------------
C     Energy matrix for hyperfine coupling  
C     ------------------------------------------------------------------
      ARR(1,1)=0  
      ARR(2,2)=-R2  
      ARR(3,3)=(ARR(2,2)-3*R1)/2  
      ARR(1,2)=R3  
      ARR(1,3)=R4  
      ARR(2,3)=R5  
      ARR(2,1)=ARR(1,2)  
      ARR(3,2)=ARR(2,3)  
      ARR(3,1)=ARR(1,3)  
      NMX=3  
      MBRANC=3  
C     ------------------------------------------------------------------
C     Energy matrix diagonalization
C     ------------------------------------------------------------------
      DO 45 I=1,NMX  
         DO 45 J=1,NMX  
            ARI(I,J)=0  
45    CONTINUE  

      CALL F02AXF(ARR,MBRANC,ARI,MBRANC,NMX,WR,ZRR,MBRANC,ZRI,MBRANC,
     &            WK1,WK2,WK3,0)  
  
      CHIPARA=WR(3)-(WR(2)+WR(1))/2  
      CHIPERP=WR(2)-WR(1)
      IF (DABS(CHIPERP).GT.2./3.*DABS(CHIPARA)) THEN
          CHIPARA=DABS(WR(1)-(WR(3)+WR(2))/2)
          CHIPERP=DABS(WR(3)-WR(2))
      END IF
      WRITE(6,*)'Eigenvalues:',WR(1),WR(2),WR(3)  

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE F01BCF(N,TOL,Z,IZ,W,IW,D,E,C,S)  
C     ------------------------------------------------------------------
C     MARK 3 RELEASE. NAG COPYRIGHT 1972.  
C     MARK 4 REVISED.  
C     MARK 4.5 REVISED.
C     MARK 5C REVISED.  
C     MARK 6B REVISED IER-125 IER-127 (Jul 1978).
C     MARK 11 REVISED VECTORISATION (Jan 1984).  
C     MARK 11.5 (F77) REVISED (Sept 1985).
C  
C     Reduces a complex Hermitian matrix to real tridiagonal form from
C     which the eigenvalues and eigenvectors can be found using
C     subroutine F02AYF. The Hermitian matrix A=A(1) is reduced to the
C     tridiagonal matrix A(N-1) by (N-2) unitary transformations.
C     The Householder reduction itself does not give a real tridiagonal
C     matrix, because the off-diagonal elements are complex. They are
C     subsequently made real by a diagonal transformation.
C
C     April 1st. 1972  
C     ------------------------------------------------------------------
C     Scalar arguments
C     ------------------------------------------------------------------
      DOUBLE PRECISION  TOL  
      INTEGER           IW, IZ, N  
C     ------------------------------------------------------------------
C     Array arguments
C     ------------------------------------------------------------------
      DOUBLE PRECISION  C(N), D(N), E(N), S(N), W(IW,N), Z(IZ,N)  
C     ------------------------------------------------------------------
C     Local scalars
C     ------------------------------------------------------------------
      DOUBLE PRECISION  CO, F, FI, FR, G, GI, GR, H, HH, R, SI  
      INTEGER           I, II, J, K, L  
C     ------------------------------------------------------------------
C     External subroutines  
C     ------------------------------------------------------------------
      EXTERNAL          F01BCY, F01BCZ  
C     ------------------------------------------------------------------
C     Intrinsic functions 
C     ------------------------------------------------------------------
      INTRINSIC         ABS, SQRT  
C     ------------------------------------------------------------------
C     Executable statements 
C     ------------------------------------------------------------------
      DO 20 I = 1, N  
         D(I) = Z(N,I)  
         E(I) = -W(N,I)  
20    CONTINUE  
      IF (N.EQ.1) GO TO 540  
      DO 360 II = 2, N  
         I = N - II + 2  
         L = I - 2  
         G = 0.0D0  
         FR = D(I-1)  
         FI = E(I-1)  
         IF (L.EQ.0) GO TO 60  
         DO 40 K = 1, L  
            G = G + D(K)*D(K) + E(K)*E(K)  
40       CONTINUE  
60       H = G + FR*FR + FI*FI  
C     ------------------------------------------------------------------
C     L is now I-1  
C     ------------------------------------------------------------------
         L = L + 1  
         IF (ABS(FR)+ABS(FI).NE.0.0D0) GO TO 80  
         R = 0.0D0  
         CO = 1.0D0  
         C(I) = 1.0D0  
         SI = 0.0D0  
         S(I) = 0.0D0  
         GO TO 140  
80       IF (ABS(FR).LT.ABS(FI)) GO TO 100  
         R = ABS(FR)*SQRT(1.0D0+(FI/FR)**2)  
         GO TO 120  
100      R = ABS(FI)*SQRT(1.0D0+(FR/FI)**2)  
120      SI = FI/R  
         S(I) = -SI  
         CO = FR/R  
         C(I) = CO  
140      IF (G.LE.TOL) GO TO 280  
         G = -SQRT(H)  
         E(I) = G  
C     ------------------------------------------------------------------
C     E(I) has its final real value  
C     ------------------------------------------------------------------
         H = H - R*G  
C     ------------------------------------------------------------------
C        S*S + SR  
C     ------------------------------------------------------------------
         D(I-1) = (R-G)*CO  
         E(I-1) = (R-G)*SI  
         DO 160 J = 1, L  
            Z(J,I) = D(J)  
            W(J,I) = E(J)  
160      CONTINUE  

         CALL F01BCZ(Z,IZ,W,IW,L,D,E,C,S)  
C     ------------------------------------------------------------------
C     Form P  
C     ------------------------------------------------------------------
         DO 180 J = 1, L  
            C(J) = C(J)/H  
            S(J) = S(J)/H  
180      CONTINUE  
         FR = 0.0D0  
         DO 200 J = 1, L  
            FR = FR + C(J)*D(J) + S(J)*E(J)  
200      CONTINUE  
C     ------------------------------------------------------------------
C     Form K  
C     ------------------------------------------------------------------
         HH = FR/(H+H)  
C     ------------------------------------------------------------------
C     Form Q  
C     ------------------------------------------------------------------
         DO 220 J = 1, L  
            C(J) = C(J) - HH*D(J)  
            S(J) = S(J) - HH*E(J)  
220      CONTINUE  
C     ------------------------------------------------------------------
C     Now form reduced A  
C     ------------------------------------------------------------------
         DO 260 J = 1, L  
            FR = D(J)  
            FI = E(J)  
            GR = C(J)  
            GI = S(J)  
            DO 240 K = J, L  
               Z(K,J) = (((Z(K,J)-GR*D(K))-GI*E(K))-FR*C(K)) - FI*S(K)  
               W(K,J) = (((W(K,J)-GR*E(K))+GI*D(K))-FR*S(K)) + FI*C(K)  
240         CONTINUE  
            D(J) = Z(L,J)  
            Z(I,J) = 0.0D0  
            E(J) = -W(L,J)  
            W(I,J) = 0.0D0  
            W(J,J) = 0.0D0  
260      CONTINUE  
         GO TO 340  
280      E(I) = R  
         H = 0.0D0  
         DO 300 J = 1, L  
            Z(J,I) = D(J)  
            W(J,I) = E(J)  
300      CONTINUE  
         DO 320 J = 1, L  
            Z(I,J) = 0.0D0  
            D(J) = Z(I-1,J)  
            W(I,J) = 0.0D0  
            E(J) = -W(I-1,J)  
320      CONTINUE  
340      D(I) = H  
360   CONTINUE  
C     ------------------------------------------------------------------
C     We now form the product of the Householder matrices,
C     overwriting on Z and W
C     ------------------------------------------------------------------
      DO 500 I = 2, N  
         L = I - 1  
         Z(N,L) = Z(L,L)  
         Z(L,L) = 1.0D0  
         W(N,L) = E(L)  
         W(L,L) = 0.0D0  
         H = D(I)  
         IF (H.EQ.0.0D0) GO TO 460  
         DO 380 K = 1, L  
            D(K) = 0.0D0  
            E(K) = 0.0D0  
380      CONTINUE  

         CALL F01BCY(Z,IZ,W,IW,L,L,Z(1,I),W(1,I),D,E)  

         DO 400 K = 1, L  
            D(K) = D(K)/H  
            E(K) = -E(K)/H  
400      CONTINUE  
         DO 440 J = 1, L  
            DO 420 K = 1, L  
               Z(K,J) = Z(K,J) - Z(K,I)*D(J) + W(K,I)*E(J)  
               W(K,J) = W(K,J) - Z(K,I)*E(J) - W(K,I)*D(J)  
420         CONTINUE  
440      CONTINUE  
460      DO 480 J = 1, L  
            Z(J,I) = 0.0D0  
            W(J,I) = 0.0D0  
480      CONTINUE  
500   CONTINUE  
      W(N,N) = E(N)  
      DO 520 I = 1, N  
         D(I) = Z(N,I)  
         Z(N,I) = 0.0D0  
         E(I) = W(N,I)  
         W(N,I) = 0.0D0  
520   CONTINUE  
540   Z(N,N) = 1.0D0  
      W(N,N) = 0.0D0  
      E(1) = 0.0D0  
C     ------------------------------------------------------------------
C     Now we multiply by the COSTHETA + I SINTHETA column factors
C     ------------------------------------------------------------------
      CO = 1.0D0  
      SI = 0.0D0  
      IF (N.EQ.1) RETURN  
      DO 580 I = 2, N  
         F = CO*C(I) - SI*S(I)  
         SI = CO*S(I) + SI*C(I)  
         CO = F  
         DO 560 J = 1, N  
            F = Z(J,I)*CO - W(J,I)*SI  
            W(J,I) = Z(J,I)*SI + W(J,I)*CO  
            Z(J,I) = F  
560      CONTINUE  
580   CONTINUE  

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE F01BCY(AR,IAR,AI,IAI,M,N,BR,BI,CR,CI)  
C     ------------------------------------------------------------------
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.  
C     MARK 11.5 (F77) REVISED (SEPT 1985).  
C  
C     Computes C=C+(A**H)*B (COMPLEX) where:
C     A is rectangular M by N  
C     C must be distinct from B  
C     ------------------------------------------------------------------
C     Scalar arguments  
C     ------------------------------------------------------------------
      INTEGER           IAI, IAR, M, N  
C     ------------------------------------------------------------------
C     Array arguments  
C     ------------------------------------------------------------------
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), BI(M), BR(M), CI(N), CR(N)
C     ------------------------------------------------------------------  
C     Local scalars  
C     ------------------------------------------------------------------
      DOUBLE PRECISION  XI, XR  
      INTEGER           I, J  
C     ------------------------------------------------------------------
C     Executable statements  
C     ------------------------------------------------------------------
      DO 40 I = 1, N  
         XR = CR(I)  
         XI = CI(I)  
         DO 20 J = 1, M  
            XR = XR + AR(J,I)*BR(J) + AI(J,I)*BI(J)  
            XI = XI + AR(J,I)*BI(J) - AI(J,I)*BR(J)  
20       CONTINUE  
         CR(I) = XR  
         CI(I) = XI  
40    CONTINUE  

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE F01BCZ(AR,IAR,AI,IAI,N,BR,BI,CR,CI)  
C     ------------------------------------------------------------------
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.  
C     MARK 11.5 (F77) REVISED (SEPT 1985).  
C  
C     Computes C=A*B (COMPLEX) where:
C     A is a Hermitian N-by-N matrix whose lower triangle is stored in A
C     C must be distinct from B  
C     ------------------------------------------------------------------
C     Scalar arguments  
C     ------------------------------------------------------------------
      INTEGER           IAI, IAR, N  
C     ------------------------------------------------------------------
C     Array arguments  
C     ------------------------------------------------------------------
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), BI(N), BR(N), CI(N), CR(N)
C     ------------------------------------------------------------------  
C     Local scalars  
C     ------------------------------------------------------------------
      DOUBLE PRECISION  YI, YR  
      INTEGER           I, IP1, J, NM1  
C     ------------------------------------------------------------------
C     Executable statements  
C     ------------------------------------------------------------------
      DO 20 I = 1, N  
         CR(I) = 0.0D0  
         CI(I) = 0.0D0  
20    CONTINUE  
      IF (N.EQ.1) GO TO 100  
      NM1 = N - 1  
      DO 80 I = 1, NM1  
         DO 40 J = I, N  
            CR(J) = CR(J) + AR(J,I)*BR(I) - AI(J,I)*BI(I)  
            CI(J) = CI(J) + AR(J,I)*BI(I) + AI(J,I)*BR(I)  
40       CONTINUE  
         YR = CR(I)  
         YI = CI(I)  
         IP1 = I + 1  
         DO 60 J = IP1, N  
            YR = YR + AR(J,I)*BR(J) + AI(J,I)*BI(J)  
            YI = YI + AR(J,I)*BI(J) - AI(J,I)*BR(J)  
60       CONTINUE  
         CR(I) = YR  
         CI(I) = YI  
80    CONTINUE  
100   CR(N) = CR(N) + AR(N,N)*BR(N) - AI(N,N)*BI(N)  
      CI(N) = CI(N) + AR(N,N)*BI(N) + AI(N,N)*BR(N)  

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE F02AXF(AR,IAR,AI,IAI,N,WR,VR,IVR,VI,IVI,WK1,WK2,WK3,  
     &                  IFAIL)  
C     ------------------------------------------------------------------
C     MARK 3 RELEASE. NAG COPYRIGHT 1972.  
C     MARK 4.5 REVISED.
C     MARK 9 REVISED IER-327 (Sep 1981).  
C     MARK 11.5 (F77) REVISED (Sep 1985).
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (Apr 1988).  
C  
C     Eigenvalues and eigenvectors of a complex Hermitian matrix  
C     1st April 1972  
C     ------------------------------------------------------------------
C     Parameters  
C     ------------------------------------------------------------------
      CHARACTER*6       SRNAME  
      PARAMETER         (SRNAME='F02AXF')  
C     ------------------------------------------------------------------
C     Scalar arguments   
C     ------------------------------------------------------------------
      INTEGER           IAI, IAR, IFAIL, IVI, IVR, N  
C     ------------------------------------------------------------------
C     Array arguments  
C     ------------------------------------------------------------------
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), VI(IVI,N), VR(IVR,N),  
     &                  WK1(N), WK2(N), WK3(N), WR(N)  
C     ------------------------------------------------------------------
C     Local scalars  
C     ------------------------------------------------------------------
      DOUBLE PRECISION  A, B, MAX, SQ, SUM, XXXX  
      INTEGER           I, ISAVE, J  
C     ------------------------------------------------------------------
C     Locla arrays   
C     ------------------------------------------------------------------
      CHARACTER*1       P01REC(1)  
C     ------------------------------------------------------------------
C     External functions  
C     ------------------------------------------------------------------
      DOUBLE PRECISION  X02AJF, X02AKF  
      INTEGER           P01ABF  
      EXTERNAL          X02AJF, X02AKF, P01ABF  
C     ------------------------------------------------------------------
C     External subroutines   
C     ------------------------------------------------------------------
      EXTERNAL          F01BCF, F02AYF  
C     ------------------------------------------------------------------
C     Intrinsic functions   
C     ------------------------------------------------------------------
      INTRINSIC         SQRT  
C     ------------------------------------------------------------------
C     Executable statements   
C     ------------------------------------------------------------------
      ISAVE = IFAIL  
      DO 40 I = 1, N  
         IF (AI(I,I).NE.0.0D0) GO TO 140  
         DO 20 J = 1, I  
            VR(I,J) = AR(I,J)  
            VI(I,J) = AI(I,J)  
20       CONTINUE  
40    CONTINUE  

      CALL F01BCF(N,X02AKF()/X02AJF(),VR,IVR,VI,IVI,WR,WK1,WK2,WK3)  
      CALL F02AYF(N,X02AJF(),WR,WK1,VR,IVR,VI,IVI,IFAIL)  

      IF (IFAIL.EQ.0) GO TO 60  
      IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)  
      RETURN  

60    DO 120 I = 1, N  
         SUM = 0.0D0  
         MAX = 0.0D0  
         DO 80 J = 1, N  
            SQ = VR(J,I)*VR(J,I) + VI(J,I)*VI(J,I)  
            SUM = SUM + SQ  
            IF (SQ.LE.MAX) GO TO 80  
            MAX = SQ  
            A = VR(J,I)  
            B = VI(J,I)  
80       CONTINUE  
         IF (SUM.EQ.0.0D0) GO TO 120  
         SUM = 1.0D0/SQRT(SUM*MAX)  
         DO 100 J = 1, N  
            SQ = SUM*(VR(J,I)*A+VI(J,I)*B)  
            VI(J,I) = SUM*(VI(J,I)*A-VR(J,I)*B)  
            VR(J,I) = SQ  
100      CONTINUE  
120   CONTINUE  
      RETURN  
140   IFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)  

      RETURN  
      END
C     ------------------------------------------------------------------
      SUBROUTINE F02AYF(N,EPS,D,E,Z,IZ,W,IW,IFAIL)  
C     ------------------------------------------------------------------
C     MARK 3 RELEASE. NAG COPYRIGHT 1972.  
C     MARK 4 REVISED.  
C     MARK 4.5 REVISED. 
C     MARK 9 REVISED IER-326 (SEP 1981).  
C     MARK 11.5 (F77) REVISED (SEP 1985).  
C  
C     This subroutine finds the eigenvalues and eigenvectors of a  
C     Hermitian matrix, which has been reduced to a real tridiagonal
C     matrix, T, given with its diagonal elements in the array D(N) and
C     its sub-diagonal elements in the last N-1 stores of the array
C     E(N), using QL transformations. The eigenvalues are overwritten
C     on the diagonal elements in the array D in ascending order. The
C     real and imaginary parts of the eigenvectors are formed in the 
C     arrays Z,W(N,N) respectively, overwriting the accumulated
C     transformations as supplied by the subroutine F01BCF. The
C     subroutine will fail if all eigenvalues take more than 30*N
C     iterations.
C     1st April 1972  
C     ------------------------------------------------------------------
C     Parameters   
C     ------------------------------------------------------------------
      CHARACTER*6       SRNAME  
      PARAMETER         (SRNAME='F02AYF')  
C     ------------------------------------------------------------------
C     Scalar arguments  
C     ------------------------------------------------------------------
      DOUBLE PRECISION  EPS  
      INTEGER           IFAIL, IW, IZ, N  
C     ------------------------------------------------------------------
C     Array arguments   
C     ------------------------------------------------------------------
      DOUBLE PRECISION  D(N), E(N), W(IW,N), Z(IZ,N)  
C     ------------------------------------------------------------------
C     Local scalars   
C     ------------------------------------------------------------------
      DOUBLE PRECISION  B, C, F, G, H, P, R, S  
      INTEGER           I, I1, II, ISAVE, J, K, L, M, M1  
C     ------------------------------------------------------------------
C     Local arrays   
C     -----------------------------------------------------------------
      CHARACTER*1       P01REC(1)  
C     ------------------------------------------------------------------
C     External functions   
C     ------------------------------------------------------------------
      INTEGER           P01ABF  
      EXTERNAL          P01ABF  
C     ------------------------------------------------------------------
C     Intrinsic functions   
C     ------------------------------------------------------------------
      INTRINSIC         ABS, SQRT  
C     ------------------------------------------------------------------
C     Executable statements  
C     ------------------------------------------------------------------
      ISAVE = IFAIL  
      IF (N.EQ.1) GO TO 40  
      DO 20 I = 2, N  
         E(I-1) = E(I)  
20    CONTINUE  
40    E(N) = 0.0D0  
      B = 0.0D0  
      F = 0.0D0  
      J = 30*N  
      DO 300 L = 1, N  
         H = EPS*(ABS(D(L))+ABS(E(L)))  
         IF (B.LT.H) B = H  
C     ------------------------------------------------------------------
C     Look for small sub-diagonal elements
C     ------------------------------------------------------------------
         DO 60 M = L, N  
            IF (ABS(E(M)).LE.B) GO TO 80  
60       CONTINUE  
80       IF (M.EQ.L) GO TO 280  
100      IF (J.LE.0) GO TO 400  
         J = J - 1  
C     ------------------------------------------------------------------
C     Form shift
C     ------------------------------------------------------------------
         G = D(L)  
         H = D(L+1) - G  
         IF (ABS(H).GE.ABS(E(L))) GO TO 120  
         P = H*0.5D0/E(L)  
         R = SQRT(P*P+1.0D0)  
         H = P + R  
         IF (P.LT.0.0D0) H = P - R  
         D(L) = E(L)/H  
         GO TO 140  
120      P = 2.0D0*E(L)/H  
         R = SQRT(P*P+1.0D0)  
         D(L) = E(L)*P/(1.0D0+R)  
140      H = G - D(L)  
         I1 = L + 1  
         IF (I1.GT.N) GO TO 180  
         DO 160 I = I1, N  
            D(I) = D(I) - H  
160      CONTINUE  
180      F = F + H  
C     ------------------------------------------------------------------
C     QL transformation  
C     ------------------------------------------------------------------
         P = D(M)  
         C = 1.0D0  
         S = 0.0D0  
         M1 = M - 1  
         DO 260 II = L, M1  
            I = M1 - II + L  
            G = C*E(I)  
            H = C*P  
            IF (ABS(P).LT.ABS(E(I))) GO TO 200  
            C = E(I)/P  
            R = SQRT(C*C+1.0D0)  
            E(I+1) = S*P*R  
            S = C/R  
            C = 1.0D0/R  
            GO TO 220  
200         C = P/E(I)  
            R = SQRT(C*C+1.0D0)  
            E(I+1) = S*E(I)*R  
            S = 1.0D0/R  
            C = C/R  
220         P = C*D(I) - S*G  
            D(I+1) = H + S*(C*G+S*D(I))  
C     ------------------------------------------------------------------
C     Form vector  
C     ------------------------------------------------------------------
            DO 240 K = 1, N  
               H = Z(K,I+1)  
               Z(K,I+1) = S*Z(K,I) + C*H  
               Z(K,I) = C*Z(K,I) - S*H  
               H = W(K,I+1)  
               W(K,I+1) = S*W(K,I) + C*H  
               W(K,I) = C*W(K,I) - S*H  
240         CONTINUE  
260      CONTINUE  
         E(L) = S*P  
         D(L) = C*P  
         IF (ABS(E(L)).GT.B) GO TO 100  
280      D(L) = D(L) + F  
300   CONTINUE  
C     ------------------------------------------------------------------
C     Sort eigenvalues and eigenvectors  
C     ------------------------------------------------------------------
      DO 380 I = 1, N  
         K = I  
         P = D(I)  
         I1 = I + 1  
         IF (I1.GT.N) GO TO 340  
         DO 320 J = I1, N  
            IF (D(J).GE.P) GO TO 320  
            K = J  
            P = D(J)  
320      CONTINUE  
340      IF (K.EQ.I) GO TO 380  
         D(K) = D(I)  
         D(I) = P  
         DO 360 J = 1, N  
            P = Z(J,I)  
            Z(J,I) = Z(J,K)  
            Z(J,K) = P  
            P = W(J,I)  
            W(J,I) = W(J,K)  
            W(J,K) = P  
360      CONTINUE  
380   CONTINUE  
      RETURN  
400   IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)  

      RETURN  
      END  
C     ------------------------------------------------------------------
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)  
C     ------------------------------------------------------------------
C     MARK 11.5 (F77) RELEASE. NAG COPYRIGHT 1986.  
C     MARK 13 REVISED IER-621 (APR 1988).  
C     MARK 13B REVISED IER-668 (AUG 1988).  
C  
C     P01ABF is the error-handling routine for the NAG library.  
C
C     P01ABF either returns the value of IERROR through the routine  
C     name (soft failure), or terminates execution of the program  
C     (hard failure). Diagnostic messages may be output.  
C  
C     If IERROR = 0 (successful exit from the calling routine),  
C     the value 0 is returned through the routine name, and no  
C     message is ouptut.
C  
C     If IERROR is non-zero (abnormal exit from the calling routine),  
C     the action taken depends on the value of IFAIL.  
C  
C     IFAIL =  1: Soft failure, silent exit (i.e. no messages are  
C                 output)  
C     IFAIL = -1: Soft failure, noisy exit (i.e. messages are output)  
C     IFAIL =-13: Soft failure, noisy exit but standard messages from  
C                 P01ABF are suppressed  
C     IFAIL =  0: Hard failure, noisy exit  
C  
C     For compatibility with certain routines included before MARK 12  
C     P01ABF also allows an alternative specification of IFAIL in which  
C     it is regarded as a decimal integer with least significant digits  
C     CBA. Then:
C  
C     A = 0: Hard failure  A = 1: Soft failure  
C     B = 0: Silent exit   B = 1: Noisy exit  
C  
C     Except that hrad failure now always implies a noisy exit.  
C  
C     S.Hammarling, M.P.Hooper and J.J.Du Croz, NAG Central Office.  
C     ------------------------------------------------------------------
C     Scalar arguments  
C     ------------------------------------------------------------------
      INTEGER                 IERROR, IFAIL, NREC  
      CHARACTER*(*)           SRNAME  
C     ------------------------------------------------------------------
C     Array arguments   
C     ------------------------------------------------------------------
      CHARACTER*(*)           REC(*)  
C     ------------------------------------------------------------------
C     Local scalars   
C     ------------------------------------------------------------------
      INTEGER                 I, NERR  
      CHARACTER*72            MESS  
C     ------------------------------------------------------------------
C     External subroutines   
C     ------------------------------------------------------------------
      EXTERNAL                P01ABZ, X04AAF, X04BAF  
C     ------------------------------------------------------------------
C     Intrinsic functions   
C     ------------------------------------------------------------------
      INTRINSIC               ABS, MOD  
C     ------------------------------------------------------------------
C     Executable statements   
C     ------------------------------------------------------------------
      IF (IERROR.NE.0) THEN  
C     ------------------------------------------------------------------
C     Abnormal exit from calling routine  
C     ------------------------------------------------------------------
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.  
     &       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN  
C     ------------------------------------------------------------------
C     Noisy exit  
C     ------------------------------------------------------------------
            CALL X04AAF(0,NERR)  
            DO 20 I = 1, NREC  
               CALL X04BAF(NERR,REC(I))  
20          CONTINUE  
            IF (IFAIL.NE.-13) THEN  
               WRITE (MESS,FMT=99999) SRNAME, IERROR  
               CALL X04BAF(NERR,MESS)  
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN  
C     ------------------------------------------------------------------
C     Hard failure  
C     ------------------------------------------------------------------
                  CALL X04BAF(NERR,  
     &                     ' ** NAG HARD FAILURE - EXECUTION TERMINATED'  
     &                        )  
                  CALL P01ABZ  
               ELSE  
C     ------------------------------------------------------------------
C     Soft failure  
C     ------------------------------------------------------------------
                  CALL X04BAF(NERR,  
     &                        ' ** NAG SOFT FAILURE - CONTROL RETURNED')  
               END IF  
            END IF  
         END IF  
      END IF  
      P01ABF = IERROR  
      RETURN  

99999 FORMAT (' ** ABNORMAL EXIT FROM NAG LIBRARY ROUTINE ',A,': IFAIL',  
     &  ' =',I6)  

      END  
C     ------------------------------------------------------------------
      SUBROUTINE P01ABZ  
C     ------------------------------------------------------------------
C     MARK 11.5 (F77) RELEASE. NAG COPYRIGHT 1986.  
C  
C     Terminates execution when a hard failure occurs.  
C  
C     ******************** IMPLEMENTATION NOTE ********************  
C     THE FOLLOWING STOP STATEMENT MAY BE REPLACED BY A CALL TO AN  
C     IMPLEMENTATION-DEPENDENT ROUTINE TO DISPLAY A MESSAGE AND/OR  
C     TO ABORT THE PROGRAM.  
C     *************************************************************  
C     ------------------------------------------------------------------
C     Executable statements
C     ------------------------------------------------------------------
      STOP  
      END  
C     ------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION X02AJF()  
C     ------------------------------------------------------------------
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.  
C  
C     Returns  (1/2)*B**(1-P)  if ROUNDS is .TRUE.  
C     Returns  B**(1-P)  otherwise  
C     ------------------------------------------------------------------
C     Executable statements  
C     ------------------------------------------------------------------
      X02AJF =  0.11102230246251568000D-015  

      RETURN  
      END  
C     ------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION X02AKF()  
C     ------------------------------------------------------------------
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.  
C  
C     Returns  B**(EMIN-1)  (The smallest positive model number)  
C     ------------------------------------------------------------------
C     Executable statements  
C     ------------------------------------------------------------------
      X02AKF =  0.22250738585072014000D-307  

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE X04AAF(I,NERR)  
C     ------------------------------------------------------------------
C     MARK 7 RELEASE. NAG COPYRIGHT 1978. 
C     MARK 7C REVISED IER-190 (MAY 1979).
C     MARK 11.5 (F77) REVISED (SEPT 1985).  
C
C     If I = 0, sets NERR to current error message unit number  
C     (stored in NERR1).  
C     If I = 1, changes current error message unit number to  
C     value specified by NERR.  
C     ------------------------------------------------------------------
C     Scalar arguments  
C     ------------------------------------------------------------------
      INTEGER           I, NERR  
C     ------------------------------------------------------------------
C     Local scalars   
C     ------------------------------------------------------------------
      INTEGER           NERR1  
C     ------------------------------------------------------------------
C     Save statement   
C     ------------------------------------------------------------------
      SAVE              NERR1  
C     ------------------------------------------------------------------
C     Data statements   
C     ------------------------------------------------------------------
      DATA              NERR1/0/  
C     ------------------------------------------------------------------
C     Executable statements   
C     ------------------------------------------------------------------
      IF (I.EQ.0) NERR = NERR1  
      IF (I.EQ.1) NERR1 = NERR  

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE X04BAF(NOUT,REC)  
C     ------------------------------------------------------------------
C     MARK 11.5 (F77) RELEASE. NAG COPYRIGHT 1986.  
C  
C     X04BAF writes teh contents of REC to the unit dfeined by NOUT.  
C  
C     Trailing blanks are not output, except that if REC is entirely  
C     blank, a single blank character is output.  
C     If NOUT.LT.0, i.e. if NOUT is not a valid FORTRAN unit identifier,  
C     then no output occurs.  
C     ------------------------------------------------------------------
C     Scalar arguments  
C     ------------------------------------------------------------------
      INTEGER           NOUT  
      CHARACTER*(*)     REC  
C     ------------------------------------------------------------------
C     Local scalars   
C     ------------------------------------------------------------------
      INTEGER           I  
C     ------------------------------------------------------------------
C     Intrinsic functions   
C     ------------------------------------------------------------------
      INTRINSIC         LEN  
C     ------------------------------------------------------------------
C     Executable statements   
C     ------------------------------------------------------------------
      IF (NOUT.GE.0) THEN  
C     ------------------------------------------------------------------
C     Remove trailing blanks  
C     ------------------------------------------------------------------
         DO 20 I = LEN(REC), 2, -1  
            IF (REC(I:I).NE.' ') GO TO 40  
20       CONTINUE  
C     ------------------------------------------------------------------
C     Write record to external file  
C     ------------------------------------------------------------------
40       WRITE (NOUT,FMT=99999) REC(1:I)  
      END IF  
      RETURN  

99999 FORMAT (A)
        
      END

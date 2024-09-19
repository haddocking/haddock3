C     ------------------------------------------------------------------
      SUBROUTINE FANTALINCCR(DJI,DJK,ANGFACT,NNAT,XOBS,XTOLPROT,
     &                       NOMEFILE1,ERWPROT,ERRORE)
C     ------------------------------------------------------------------
C     Handles the stuff for the calculation of CCR constant coefficient
C
C     CX,CY,CZ(STRUCTURE,ATOMNUMBER) keep coordinates of nuclei
C     FX,FY,FZ(STRUCTURE,ATOMNUMBER) keep coordinates of metal
C     NAT is a vector with the number of atoms
C     OBS is a vector with the observed cross correlation rates
C     TOLPROT is a vector with the errors on CCR
C     NOMEFILE1 is the name of the file to save deviations
C     ERWPROT is a vector with the weights
C     ERRORE switches the use(1) of weights
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'supccr.inc'
      INCLUDE 'fantaxplor.inc'

      INTEGER NNAT,METHOD
      DOUBLE PRECISION AK,A
      COMMON /CRUN/AK
      DOUBLE PRECISION DJI(1000),DJK(1000),ANGFACT(1000),
     &                 XOBS(1000),XTOLPROT(1000),ERWPROT(1000)
      CHARACTER NOMEFILE1*132 
      INTEGER ERRORE

      DO I=1,1000
         AK1(I)=ANGFACT(I)
         RK1(I)=DJI(I)
         RK2(I)=DJK(I)
      END DO

      PRINT*,'                      '
      PRINT*,'===---FANTACROSS---==='
C     ------------------------------------------------------------------
C     IHP is the number of CCR
C     ------------------------------------------------------------------
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

      OPEN(322,FILE='FANTACROSS.LOG')
C     ------------------------------------------------------------------
C     Vectors are passed from Xplor-NIH core
C     ------------------------------------------------------------------
      DO NIA=1,NAT
         OBS(NIA)=XOBS(NIA)
C     ------------------------------------------------------------------
C     No tolerance in fitting
C     ------------------------------------------------------------------
         TOLPROT(NIA)=0
      END DO
      CLOSE (322)
      FILENAME3=NOMEFILE1
      PRINT*,'NOMEFILE1 is:',NOMEFILE1
      OPEN(744,FILE=FILENAME3)

      CALL CSIMPLEXRUN()

      RETURN    
      END  
C     ------------------------------------------------------------------
      SUBROUTINE CSIMPLEXRUN()  
C     ------------------------------------------------------------------
C     Calculates the CCR constant coefficient which minimizes the 
C     experimental-to-calculated squared difference
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'supccr.inc'

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
               VAL(I)=CCASH(EXTR)  
11       CONTINUE  
         CALL CSIMPCALC(SIMP,VAL)  
C     ------------------------------------------------------------------
C     Criterion of minimization of optimal solution
C     ------------------------------------------------------------------
         IF (OLDRESID.GE.RESID) THEN  
             OLDRESID=RESID  
             IOLDVIO=IVIOLATION  
             DO 13 I=1,NFE  
                OPTPHI(I)=PHI(I)  
13           CONTINUE  
         ENDIF  
1000  CONTINUE  

      DO 998 NK=1,NFE  
         EXTR((NK-1)*MPAR+1)=OPTPHI(NK)  
998   CONTINUE 

      PRINT*,' '  
      WRITE(*,FMT='(1X,A)') ' *****    BEST SOLUTION    *****'  
      PRINT*,' '  
      WRITE(*,20) IOLDVIO,OLDRESID  
      WRITE(*,21) (OPTPHI(I),I=1,NFE)  
      WRITE(*,*) 'Value of K =',OPTPHI(1)
      WRITE(744,*) 'Value of K =',OPTPHI(1)

20    FORMAT (1X,' VIOLATIONS ',I6,1X,' RESIDUALS ',F10.3)  
21    FORMAT (1X,'K = ',F15.3)  
  
      ERR=CCASH(EXTR)  
      CALL CDISPSHIFT

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE CSIMPCALC(P,Y)  
C     ------------------------------------------------------------------
C     Simplex calculation
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'supccr.inc'
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
          WRITE (*,200) ITER,IVIOLATION,RESID  
          WRITE (*,201) (P(1,(I-1)*MPAR+5),I=1,NFE)  
200       FORMAT (1X,'Iters ',I6,1X,' Violations ',I6,  
     &            1X,'Residuals= ',F10.3)  
201       FORMAT (1X,'K = ',F9.3)
          DO 997 NK=1,NFE  
             PHI(NK)=P(1,(NK-1)*MPAR+1)  
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
      YR=CCASH(PR)  
      DO 15 J=1,NP  
         PK(J)=PR(J)  
15    CONTINUE  
      YK=YR  
      IF (YR.LT.Y(NP)) THEN  
          IF (YR.LT.Y(1)) THEN  
              DO 16 J=1,NP  
                 PE(J)=GA*PR(J)+(1.-GA)*PBAR(J)  
16            CONTINUE  
              YE=CCASH(PE)  
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
          YC=CCASH(PC)  
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
                 Y(I)=CCASH(PR)  
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
                     Y(IL)=CCASH(EXTR)  
2011              CONTINUE  
                  GOTO 1   
              ENDIF   
              IF (ABS(OLDPOPPA).LT.1.D-9) THEN  
                  ITER=50000  
              ENDIF  
              DO 26 J=1,NP  
                 PK(J)=0.5*(P(1,J)+P(MP,J))  
26            CONTINUE  
              YK=CCASH(PK)  
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
      FUNCTION CCASH(VETT)  
C     ------------------------------------------------------------------
C     Calculation of deviations and violations
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'supccr.inc'  

      DIMENSION VETT(NP)
      INTEGER M1  

      IVIOLATION=0  
      TMP1=0  
      I=1  
      DO WHILE (I.LE.NHP*NSTR)  
         CSHIFT(I)=0.0  
         I=I+1  
      ENDDO  

      DO 1 N=1,NSTR  
         DO 2 M=1,NFE  
            P=VETT((M-1)*MPAR+1)  
            TMP2=0.0  
            I=1  
            IHP=(N-1)*NHP 
            DO 10 WHILE (I.LE.NHP)  
               NPROTML=MLPROT(IHP+I)
               CSHIFT(IHP+I)=CSHIFT(IHP+I)+AK1(IHP+I)/(RK2(IHP+I)**3*
     &                       RK1(IHP+I)**3)
               I=I+1  
10          CONTINUE  
2        CONTINUE  
         DO 3 I=1,NHP  
            TMP2=ABS(CSHIFT(IHP+I)-OBS(IHP+I))-TOLPROT(IHP+I)  
            IF (TMP2.GT.0.0) THEN  
                IVIOLATION=IVIOLATION+1  
                TMP1=TMP1+TMP2**2*WPROT(IHP+I)/DBLE(MLPROT(IHP+I))  
            ENDIF  
3        CONTINUE   
1     CONTINUE  
  
      CCASH=TMP1  

      RETURN  
      END  
C     ------------------------------------------------------------------
      SUBROUTINE CDISPSHIFT 
C     ------------------------------------------------------------------
C     Display of violations
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INCLUDE 'supccr.inc'
      INCLUDE 'fantaxplor.inc'  

      CHARACTER FILESTR*21,SS*40  
      DIMENSION VV(MAXSTR,MAXOS),VMAX(MAXOS),VMIN(MAXOS)  

      RM=0  
      TMP1=0  
      SMED=0  
      STDEV=0  
      DO 2 I=1,NSTR  
         DO 2 J=1,NHP  
            TMP2=ABS(OBS(J+(I-1)*NHP)-CSHIFT(J+(I-1)*NHP))-
     &           TOLPROT(J+(I-1)*NHP)  
            IF (TMP2.GT.0.0) THEN  
                IF (CSHIFT(J+(I-1)*NHP).GT.OBS(J+(I-1)*NHP)) THEN  
                    VV(I,J)=TMP2  
                ELSE  
                    VV(I,J)=-TMP2  
                ENDIF  
            ELSE  
                VV(I,J)=0.  
            ENDIF  
2     CONTINUE  
      DO 44 I=1, NHP  
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
         TMP2=ABS(OBS(J)-CSHIFT(J))-TOLPROT(J)
         NIA=NIA+1
         IF (TMP2.GT.0.0) THEN  
             IVIOLATION=IVIOLATION+1  
             RM=RM+1  
             SMED=SMED+TMP2**2*WPROT(J)/MLPROT(J)  
             TMP1=TMP2**2*WPROT(J)/MLPROT(J)  
             WRITE (744,'(1X,A4,1X,A4,1X,A4,1X,F7.3,1X,F7.3,1X,F10.5)')
     &              NRESIDFA(NIA),NRESFA(NIA),NTYPEFA(NIA),
     &              OBS(J),CSHIFT(J),TMP1
15           FORMAT(1X,A,1X,F7.3,1X,F7.3,1X,F7.3)  
         ELSE  
             WRITE (744,'(1X,A4,1X,A4,1X,A4,1X,F7.3,1X,F7.3)')  
     &              NRESIDFA(NIA),NRESFA(NIA),NTYPEFA(NIA),
     &              OBS(J),CSHIFT(J)
        ENDIF  
3     CONTINUE  
      SMED=SMED/RM  

      DO 41 J=1,NHP*NSTR  
         TMP2=ABS(OBS(J)-CSHIFT(J))-TOLPROT(J)  
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

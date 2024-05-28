C     ------------------------------------------------------------------
      SUBROUTINE EXRDC (EDI,WHICH)
C     ------------------------------------------------------------------
C     Calls EXRDC2, which does the actual energy calculation
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "xdipo_rdc.inc"
      INCLUDE  "heap.inc"

      DOUBLE PRECISION EDI
      CHARACTER*7 WHICH

      CALL EXRDC2(EDI,HEAP(XRDCIPTR),HEAP(XRDCJPTR),HEAP(XRDCKPTR),
     &                HEAP(XRDCLPTR),HEAP(XRDCMPTR),HEAP(XRDCNPTR),
     &                HEAP(XRDCILST),HEAP(XRDCJLST),HEAP(XRDCKLST),
     &                HEAP(XRDCLLST),HEAP(XRDCMLST),HEAP(XRDCNLST),
     &                HEAP(XRDCOBSPTR),HEAP(XRDCERRPTR),
     &                HEAP(CALCXRDCPTR),WHICH)

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE EXRDC2 (EDI,ATOMIPTR,ATOMJPTR,ATOMKPTR,ATOMLPTR, 
     &                       ATOMMPTR,ATOMNPTR,ATOMILST,ATOMJLST,
     &                       ATOMKLST,ATOMLLST,ATOMMLST,ATOMNLST,
     &                       XRDCOBS,XRDCERR,XRDCCALC,WHICH)
C     ------------------------------------------------------------------
C     Calculates residual dipolar coupling energies
C
C     Energies are of the form:
C        E=K*(DELTARDC**2)
C     where :
C        K=FORCE CONSTANT
C        DELTARDC=CALCULATED RDC-OBSERVED RDC
C     and residual dipolar coupling function is defined as:
C        RDC=A1*(3*COS(THETA)^2-1)+(3/2)*A2*(SIN(THETA)^2*COS(2*PHI))
C
C     WHICH is a flag that switches between energy/force calculations
C     and RDC calculations for violations
C
C     Note that:
C        Atom I=Origin of the tensor
C        Atom J=Z direction
C        Atom K=X direction
C        Atom L=Y direction
C        Atom M=Coupled nucleus (e.g. N)
C        Atom N=Detected nucleus (e.g. H)
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

      DOUBLE PRECISION DIPOE
      COMMON /DIPOENERGY/DIPOE
      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*)
      INTEGER ATOMLPTR(*),ATOMMPTR(*),ATOMNPTR(*)
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*)
      INTEGER ATOMLLST(*),ATOMMLST(*),ATOMNLST(*)
      DOUBLE PRECISION XRDCOBS(*),XRDCERR(*)
      DOUBLE PRECISION XRDCCALC(*),EDI
      CHARACTER*7 WHICH
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      INTEGER COUNT,CLASS,MM,NN,II,I,J,K,L,NA,M,N,SWAP
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
     &                 O1,O2,O3,O4 ,O5 ,O6 ,O7 ,O8,O9,O10,
     &                 O11,O12,O13,O14,O15,O16,O17,O18,O19,
     &                 O20,O21,O22,O23,O24,O25,O26,O27,O28,
     &                 O29,O30,O31,O32,O33,O34,O35,O36,O37,
     &                 O38,O39,O40,O41,O42,O43,O44,O45,O46,
     &                 O47,O48,O49,O50,O51,O52,O53,O54,O55,
     &                 O56,O57,O58,O59,O60,O61,O62,O63,O64,
     &                 O65,O66,O67,O68,O69,O70,O71,O72,O73,
     &                 O74,O75,O76,O77,O78,O79,O80,O81,O82
      DOUBLE PRECISION OBSXRDC,ERRXRDC,K1,
     &                 COEF1,COEF2,A,B,DELTA,
     &                 DT,DP,DELTAXRDC,
     &                 DD,DR,CALCXRDC
C     ------------------------------------------------------------------
C     Zero out partial energy
C     ------------------------------------------------------------------
      EDI=ZERO

      CLASS=1
      K1=XRDCFORCES(1)
      COEF1=XRDCCOEF1(1)
      COEF2=XRDCCOEF2(1)

      COUNT=0

      DO WHILE (COUNT.LT.RDCNUM)
         COUNT=COUNT+1
C     ------------------------------------------------------------------
C     Reset individual E to zero
C     ------------------------------------------------------------------
         E=0

         DO WHILE ((XRDCASSNDX(CLASS).LT.COUNT).OR.
     &             (XRDCASSNDX(CLASS).EQ.0))
            CLASS=CLASS+1
         END DO

         IF (XRDCASSNDX(CLASS).EQ.XRDCASSNDX(CLASS-1)) THEN
             COUNT=COUNT-1
         ENDIF

         K1=XRDCFORCES(CLASS)
         COEF1=XRDCCOEF1(CLASS)
         COEF2=XRDCCOEF2(CLASS)
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

         OBSXRDC=XRDCOBS(COUNT)
         ERRXRDC=XRDCERR(COUNT)
C     ------------------------------------------------------------------
C     Initialzie calculated residual dipolar coupling and counter
C     ------------------------------------------------------------------
         CALCXRDC=0
         II=0
         MCOUNT=0
         NCOUNT=0
C     ------------------------------------------------------------------
C     Check for correct permutations of paired atoms to get the nucleus-
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
                   DTYN(II)= -(O28*(-O40+O22*O3*O9))
                   DTZN(II)=-(O28*(-O41+O22*O5*O9))
  
                   PHI(II)=ATAN(O49*O58*O63*O67)

                   DPXI(II)=O78*(O29*O49*O58*O63-O29*O49*O58*O67*O69+
     &                           O72-O74)
                   DPYI(II)=O78*(O33*O49*O58*O63-O33*O49*O58*O67*O69+
     &                           O79-O80)
                   DPZI(II)=O78*(O35*O49*O58*O63-O35*O49*O58*O67*O69+
     &                           O81-O82)
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
     &                             O49*O58*O63*(XI-XL))
                   DPYM(II)=O78*(-(O49*O58*O67*O69*(YI-YK))+
     &                             O49*O58*O63*(YI-YL))
                   DPZM(II)=O78*(-(O49*O58*O67*O69*(ZI-ZK))+
     &                             O49*O58*O63*(ZI-ZL))
                   DPXN(II)=(O49*O50*O58*O63-O42*O49*O58*O67*O69)*O78
                   DPYN(II)=(O49*O52*O58*O63-O44*O49*O58*O67*O69)*O78
                   DPZN(II)=(O49*O54*O58*O63-O46*O49*O58*O67*O69)*O78
C     ------------------------------------------------------------------
C     Calculate the distance between the two atoms M, N (coupled nuclei)
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
C     Calculate the residual dipolar coupling
C     ------------------------------------------------------------------
                   A=COEF1
                   B=COEF2
                   CALCXRDC=A*(3.*COS(THETA(II))*COS(THETA(II))-1.)+
     &                      B*(3./2.)*(SIN(THETA(II))*SIN(THETA(II))*
     &                      COS(2.*PHI(II)))
               END IF
            END DO
         END DO

         IF (WHICH.EQ.'ANALYZE') THEN
             XRDCCALC(COUNT)=CALCXRDC
         END IF

         DELTA=CALCXRDC-OBSXRDC
C     ------------------------------------------------------------------
C     Adjust the deviation based on the error range
C     ------------------------------------------------------------------
         IF ((DELTA.LT.0.000).AND.(ABS(DELTA).GT.ERRXRDC)) THEN
              DELTAXRDC=DELTA+ERRXRDC
         ELSE IF ((DELTA.GT.0.000).AND.(DELTA.GT.ERRXRDC)) THEN
              DELTAXRDC=DELTA-ERRXRDC
         ELSE
              DELTAXRDC=0.0
         END IF

         DD=DELTAXRDC

         II=0
         PCOUNT=0
         DO MM=ATOMMPTR(COUNT)+1,ATOMMPTR(COUNT+1)
            DO NN=ATOMNPTR(COUNT)+1,ATOMNPTR(COUNT+1)
               II=II+1
               IF (OUTFLAG(II).NE.1) THEN
                   PCOUNT=PCOUNT+1
C     ------------------------------------------------------------------
C     Taking care of derivative
C     ------------------------------------------------------------------
                   DT=(A*(-6.)*SIN(THETA(II))*
     &                 COS(THETA(II))+B*3.*SIN(THETA(II))*
     &                 COS(THETA(II))*COS(2.*PHI(II)))
                   DP=(B*(-3.)*SIN(THETA(II))*
     &                 SIN(THETA(II))*SIN(2.*PHI(II)))
                   DR=0.0
                   IF (DD.EQ.0.0) THEN
                       E=E+0.0
                       DF1(II)=0.0
                       DF2(II)=0.0
                       DF3(II)=0.0
                   ELSE
                       E=E+K1*(DELTAXRDC**2)
                       DF1(II)=2*K1*DELTAXRDC*DT
                       DF2(II)=2*K1*DELTAXRDC*DP
                       DF3(II)=2*K1*DELTAXRDC*DR
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
             II=0
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

      DIPOE=EDI

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE READXRDC
C     ------------------------------------------------------------------
C     Reads in residual dipolar coupling information
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "comand.inc"
      INCLUDE  "xdipo_rdc.inc"
      INCLUDE  "funct.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "heap.inc"
      INCLUDE  "numbers.inc"

      INTEGER FDSWITCH,FDSAVE,FDERR
      COMMON  /FDPARAM/FDSWITCH,FDSAVE,FDERR
      DOUBLE PRECISION FDEPERC
      INTEGER FDENUMC
      COMMON /FDERRORE/FDEPERC,FDENUMC

      INTEGER COUNT,SPTR,OLDCLASS,OLDMAXXRDC,FNCLASS,FNUMC
      DOUBLE PRECISION K1,K2,CUTOFF,COEF1,COEF2
      DOUBLE PRECISION FMEDPERC,FPERC
      CHARACTER*4 THENAME
      CHARACTER*132 NOMEFILE
      INTEGER TOLSWI

      NOMEFILE='XDIPO_RDC.OUT'

      SPTR=ALLHP(INTEG4(NATOM))
      CALL PUSEND('XDIPO_RDC>')
862   CONTINUE
      CALL NEXTWD('XDIPO_RDC>')
      CALL MISCOM('XDIPO_RDC>',USED)
      IF (.NOT.USED) THEN
C     ------------------------------------------------------------------
C     Documentation
C     ------------------------------------------------------------------
          IF (WD(1:4).EQ.'HELP') THEN
              WRITE(DUNIT,'(10X,A)')
     &              ' XRDC {<RDC-STATEMENT>} END ',
     &              ' <RDC-STATEMENT>:== ',
     &              ' ASSIgn <SEL> <SEL> <SEL> <SEL> <SEL> <SEL>',
     &              ' <REAL> <REAL>',
     &              ' * Restraint: Metal Z X Y Coupled_nucleus',
     &              ' Detected_nucleus RDC Err *',
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
     &              ' * Prints violations larger than the THREshold *',
     &              ' RESEt',
     &              ' * Erases the restraint table, but keeps NRES *',
     &              ' SAVE',
     &              ' * Filename to save RDC values *',
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
                  DIPFLG=.TRUE.
                  WRITE (DUNIT,'(A)') 'Tolerance ON in FANTALIN.'
              ELSE IF (TOLSWI.EQ.0) THEN
                  DIPFLG=.FALSE.
                  WRITE (DUNIT,'(A)') 'Tolerance OFF in FANTALIN.'
              ELSE
                  WRITE (DUNIT,'(A)') 'Unknown switch. Using default.'
              END IF
          ELSE IF (WD(1:4).EQ.'SAVE') THEN
              CALL NEXTFI('Filename =',NOMEFILE)
              PRINT*,'Name of the file to save RDC values :',NOMEFILE
          ELSE IF (WD(1:3).EQ.'FME') THEN
              CALL NEXTF('Percentage for average (0/100) =',FMEDPERC)
              CALL NEXTI('Average on class number =',FNCLASS)
              IF (FNCLASS.GT.NCLASSES) THEN
                  PRINT*,'%FRUN-ERR: This class does not exist...'
              ELSE
                  CALL FANTADIPOMEDIA(FNCLASS,FMEDPERC)
              END IF
          ELSE IF (WD(1:5).EQ.'ERRON') THEN
              CALL NEXTI('Number of cycles =',FDENUMC)
              CALL NEXTF('Percentage to discard (0/100) =',FDEPERC)
              FDEPERC=FDEPERC/100.0
              FDERR=1
              PRINT*,'Monte Carlo error evaluation ON'
          ELSE IF (WD(1:6).EQ.'ERROFF') THEN
              FDERR=0
              PRINT*,'Monte Carlo error evaluation OFF'
          ELSE IF (WD(1:3).EQ.'FON') THEN
              FDSWITCH=1
              PRINT*,'Tensor auto-update mode ON'
          ELSE IF (WD(1:4).EQ.'FOFF') THEN
              FDSWITCH=0
              PRINT*,'Tensor auto-update mode OFF'
          ELSE IF (WD(1:3).EQ.'SON') THEN
              FDSAVE=1
              PRINT*,'Saving mode ON'
          ELSE IF (WD(1:4).EQ.'SOFF') THEN
              FDSAVE=0
              PRINT*,'Saving mode OFF'
          ELSE IF (WD(1:4).EQ.'FRUN') THEN
              CALL NEXTI('FANTALIN on class number =',FNCLASS)
              IF (FNCLASS.GT.NCLASSES.OR.FNCLASS.EQ.0) THEN
                  PRINT*,'%FRUN-ERR: This class does not exist...'
             ELSE
                  CALL FANTADIPO(HEAP(XRDCIPTR),HEAP(XRDCJPTR),
     &                          HEAP(XRDCKPTR),HEAP(XRDCLPTR),
     &                          HEAP(XRDCMPTR),HEAP(XRDCNPTR),
     &                          HEAP(XRDCILST),HEAP(XRDCJLST),
     &                          HEAP(XRDCKLST),HEAP(XRDCLLST),
     &                          HEAP(XRDCMLST),HEAP(XRDCNLST),
     &                          HEAP(XRDCOBSPTR),HEAP(XRDCERRPTR),
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
                  IF (XRDCCLASSNAMES(COUNT).EQ.THENAME) THEN
                      MODE=UPDATE
                      CURCLASS=COUNT
                  END IF
               END DO
               IF (MODE.EQ.NEW) THEN
C     ------------------------------------------------------------------
C     Make sure you can't add more than the maximum number of classes
C     ------------------------------------------------------------------
                   IF (OLDCLASS.EQ.MAXXRDCCLASSES) THEN
                       CALL DSPERR('XDIPO_RDC','Too many classes.')
                       CALL DSPERR('XDIPO_RDC',
     &                      'Increase MAXXRDCCLASSES and recompile.')
                       CALL WRNDIE(-5, 'READXRDC',
     &                      'Too many anisotropy classes.')
                   END IF
                   NCLASSES=NCLASSES+1
                   CURCLASS=NCLASSES
                   XRDCCLASSNAMES(CURCLASS)=THENAME
C     ------------------------------------------------------------------
C     If this isn't the first class, close off the old class
C     ------------------------------------------------------------------
                   IF (NCLASSES.GT.1) THEN
                       XRDCASSNDX(OLDCLASS)=RDCNUM
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
     &              XRDCCLASSNAMES(CURCLASS),' to ',K1
                    XRDCFORCES(CURCLASS)=K1
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
     &              XRDCCLASSNAMES(CURCLASS),' to ',COEF1,COEF2
                    XRDCCOEF1(CURCLASS)=COEF1
                    XRDCCOEF2(CURCLASS)=COEF2
C     ------------------------------------------------------------------
C     Reset residual dipolar couplings database
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'RESE') THEN
              CALL XRDCDEFAULTS
              CALL ALLOCXRDC(0,MAXXRDC)
C     ------------------------------------------------------------------
C     Change number of assignment slots
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'NRES') THEN
              OLDMAXXRDC=MAXXRDC
              CALL NEXTI('Number of slots =',MAXXRDC)
              CALL ALLOCXRDC(OLDMAXXRDC,MAXXRDC)
              WRITE(DUNIT,'(A,I8,A)')
     &             'XDIPO_RDC: Allocating space for',MAXXRDC,
     &             'number of restraints.'
C     ------------------------------------------------------------------
C     Read in an assignment
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'ASSI') THEN
C     ------------------------------------------------------------------
C     Make sure you can't add more assignments than you have slots for
C     ------------------------------------------------------------------
              IF (XRDCMX.EQ.MAXXRDC) THEN
                  CALL DSPERR('XDIPO_RDC','Too many assignments.')
                  CALL DSPERR('XDIPO_RDC',
     &                        'Increasing NREStraints by 100.')
                  OLDMAXXRDC=MAXXRDC
                  MAXXRDC=MAXXRDC+100
                  CALL ALLOCXRDC(OLDMAXXRDC,MAXXRDC)
              END IF
C     ------------------------------------------------------------------
C     If there isn't a class specified, start a default class
C     ------------------------------------------------------------------
              IF (CURCLASS.EQ.0) THEN
                  NCLASSES=1
                  CURCLASS=1
              END IF
              CALL READXRDC2(HEAP(XRDCIPTR),HEAP(XRDCJPTR), 
     &                       HEAP(XRDCKPTR),HEAP(XRDCLPTR),
     &                       HEAP(XRDCMPTR),HEAP(XRDCNPTR),
     &                       HEAP(XRDCILST),HEAP(XRDCJLST),
     &                       HEAP(XRDCKLST),HEAP(XRDCLLST),
     &                       HEAP(XRDCMLST),HEAP(XRDCNLST),
     &                       HEAP(XRDCOBSPTR),HEAP(XRDCERRPTR),
     &                       HEAP(SPTR))
C     ------------------------------------------------------------------
C     Print violations
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'PRIN') THEN
              CALL NEXTWD('PRINT>')
              IF (WD(1:4).NE.'THRE') THEN
                  CALL DSPERR('XDIPO_RDC',
     &                 'PRINt expects THREshold parameter.')
              ELSE
                  CALL NEXTF('THREshold =',CUTOFF)
                  IF (CUTOFF.LT.ZERO) THEN
                      CALL DSPERR('XDIPO_RDC',
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
                         IF (XRDCCLASSNAMES(COUNT).EQ.THENAME) THEN
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
                  CALL PRINTXRDC(CUTOFF,HEAP(CALCXRDCPTR),
     &                           HEAP(XRDCOBSPTR),HEAP(XRDCERRPTR),
     &                           HEAP(XRDCIPTR),HEAP(XRDCJPTR),
     &                           HEAP(XRDCKPTR),HEAP(XRDCLPTR),
     &                           HEAP(XRDCMPTR),HEAP(XRDCNPTR),
     &                           HEAP(XRDCILST),HEAP(XRDCJLST),
     &                           HEAP(XRDCKLST),HEAP(XRDCLLST),
     &                           HEAP(XRDCMLST),HEAP(XRDCNLST))
                END IF
C     ------------------------------------------------------------------
C     Check for END statement
C     ------------------------------------------------------------------
          ELSE
              CALL CHKEND('XDIPO_RDC>',DONE)
          END IF
      END IF
      IF (.NOT.DONE) GOTO 862
      DONE=.FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))
      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE ALLOCXRDC (OLDSIZE,NEWSIZE)
C     ------------------------------------------------------------------
C     Reset residual dipolar coupling arrays to hold size entries

C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE  "funct.inc"
      INCLUDE  "xdipo_rdc.inc"

      INTEGER OLDSIZE, NEWSIZE

      IF (OLDSIZE.NE.0) THEN
          CALL FREHP(XRDCIPTR,INTEG4(OLDSIZE))
          CALL FREHP(XRDCJPTR,INTEG4(OLDSIZE))
          CALL FREHP(XRDCKPTR,INTEG4(OLDSIZE))
          CALL FREHP(XRDCLPTR,INTEG4(OLDSIZE))
          CALL FREHP(XRDCMPTR,INTEG4(OLDSIZE))
          CALL FREHP(XRDCNPTR,INTEG4(OLDSIZE))

          CALL FREHP(XRDCILST,INTEG4(OLDSIZE))
          CALL FREHP(XRDCJLST,INTEG4(OLDSIZE))
          CALL FREHP(XRDCKLST,INTEG4(OLDSIZE))
          CALL FREHP(XRDCLLST,INTEG4(OLDSIZE))
          CALL FREHP(XRDCMLST,INTEG4(OLDSIZE))
          CALL FREHP(XRDCNLST,INTEG4(OLDSIZE))

          CALL FREHP(XRDCOBSPTR,IREAL8(OLDSIZE))
          CALL FREHP(XRDCERRPTR,IREAL8(OLDSIZE))
          CALL FREHP(CALCXRDCPTR,IREAL8(OLDSIZE))
      END IF
      XRDCIPTR=ALLHP(INTEG4(NEWSIZE))
      XRDCJPTR=ALLHP(INTEG4(NEWSIZE))
      XRDCKPTR=ALLHP(INTEG4(NEWSIZE))
      XRDCLPTR=ALLHP(INTEG4(NEWSIZE))
      XRDCMPTR=ALLHP(INTEG4(NEWSIZE))
      XRDCNPTR=ALLHP(INTEG4(NEWSIZE))

      XRDCILST=ALLHP(INTEG4(NEWSIZE))
      XRDCJLST=ALLHP(INTEG4(NEWSIZE))
      XRDCKLST=ALLHP(INTEG4(NEWSIZE))
      XRDCLLST=ALLHP(INTEG4(NEWSIZE))
      XRDCMLST=ALLHP(INTEG4(NEWSIZE))
      XRDCNLST=ALLHP(INTEG4(NEWSIZE))

      XRDCOBSPTR=ALLHP(IREAL8(NEWSIZE))
      XRDCERRPTR=ALLHP(IREAL8(NEWSIZE))
      CALCXRDCPTR=ALLHP(IREAL8(NEWSIZE))

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE XRDCDEFAULTS
C     ------------------------------------------------------------------
C     Sets up defaults
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "xdipo_rdc.inc"

      INTEGER COUNT

      MODE=NEW
      MAXXRDC=200
      XRDCMX=200
      RDCNUM=0
      NCLASSES=0
      CURCLASS=0
      DIPFLG=.FALSE.
      DO COUNT=1, MAXXRDCCLASSES
           XRDCCLASSNAMES(COUNT)='DEFAULT'
           XRDCASSNDX(COUNT)=0
           XRDCFORCES(COUNT)=5.0
           XRDCCOEF1(COUNT)=1.0
           XRDCCOEF2(COUNT)=1.0
      END DO

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE READXRDC2(ATOMIPTR,ATOMJPTR,ATOMKPTR,ATOMLPTR,ATOMMPTR,
     &                     ATOMNPTR,ATOMILST,ATOMJLST,ATOMKLST,ATOMLLST,
     &                     ATOMMLST,ATOMNLST,XRDCOBS,XRDCERR,SEL)
C     ------------------------------------------------------------------
C     Reads actual residual dipolar coupling assignments into arrays
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "coord.inc"
      INCLUDE  "xdipo_rdc.inc"

      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*)
      INTEGER ATOMLPTR(*),ATOMMPTR(*),ATOMNPTR(*)
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*)
      INTEGER ATOMLLST(*),ATOMMLST(*),ATOMNLST(*),SEL(*)
      INTEGER ITMP,JTMP,KTMP,LTMP,II
      DOUBLE PRECISION XRDCOBS(*),XRDCERR(*)
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      INTEGER NSEL,INSERTPOS,COUNT,CURSTOP,OTHERSTOP,NFLAG
      INTEGER I
      DOUBLE PRECISION XRDCO,XRDCE
C     ------------------------------------------------------------------
C     If we're in UPDATE mode, make a space for the new line
C     ------------------------------------------------------------------
      NFLAG=0
      IF (MODE.EQ.UPDATE) THEN
          DO COUNT=RDCNUM+1,XRDCASSNDX(CURCLASS)+1, -1
             ATOMIPTR(COUNT)=ATOMIPTR(COUNT-1)
             ATOMJPTR(COUNT)=ATOMJPTR(COUNT-1)
             ATOMKPTR(COUNT)=ATOMKPTR(COUNT-1)
             ATOMLPTR(COUNT)=ATOMLPTR(COUNT-1)
             ATOMMPTR(COUNT)=ATOMMPTR(COUNT-1)
             ATOMNPTR(COUNT)=ATOMNPTR(COUNT-1)
             XRDCOBS(COUNT)=XRDCOBS(COUNT-1)
             XRDCERR(COUNT)=XRDCERR(COUNT-1)
           END DO
           CURSTOP=XRDCASSNDX(CURCLASS)
           DO COUNT=1,NCLASSES
                OTHERSTOP=XRDCASSNDX(COUNT)
                IF (OTHERSTOP.GT.CURSTOP) THEN
                    XRDCASSNDX(COUNT)=OTHERSTOP+1
                END IF
           END DO
           XRDCASSNDX(CURCLASS)=CURSTOP+1
           INSERTPOS=CURSTOP
           RDCNUM=RDCNUM+1
      ELSE
           RDCNUM=RDCNUM+1
           INSERTPOS=RDCNUM
           XRDCASSNDX(CURCLASS)=INSERTPOS
      END IF
C     ------------------------------------------------------------------
C     Reading in the atom selection in the restraint table
C     ------------------------------------------------------------------
      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_RDC',
     &                'More than 1 atom in SEL for atom I. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
         CALL DSPERR('XDIPO_RDC','Error with atom selection')
         NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL)
      IF (INSERTPOS.EQ.1) ATOMIPTR(INSERTPOS)=0
      II=ATOMIPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XRDCMX) XRDCMX=II
      ATOMILST(II)=SEL(1)
      ATOMIPTR(INSERTPOS+1)=II

      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_RDC',
     &                'More than 1 atom in SEL for atom J. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
         CALL DSPERR('XDIPO_RDC','Error with atom selection')
         NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL)
      IF (INSERTPOS.EQ.1) ATOMJPTR(INSERTPOS)=0
      II=ATOMJPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XRDCMX) XRDCMX=II
      ATOMJLST(II)=SEL(1)
      ATOMJPTR(INSERTPOS+1)=II

      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_RDC',
     &                'More than 1 atom in SEL for atom K. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XDIPO_RDC','Error with atom selection')
          NFLAG=1
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (INSERTPOS.EQ.1) ATOMKPTR(INSERTPOS)=0
      II=ATOMKPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XRDCMX) XRDCMX=II
      ATOMKLST(II)=SEL(1)
      ATOMKPTR(INSERTPOS+1)=II

      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_RDC',
     &                'More than 1 atom in SEL for atom L. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XDIPO_RDC','Error with atom selection')
          NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL)
      IF (INSERTPOS.EQ.1) ATOMLPTR(INSERTPOS)=0
      II=ATOMLPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XRDCMX) XRDCMX=II
      ATOMLLST(II)=SEL(1)
      ATOMLPTR(INSERTPOS+1)=II
C     ------------------------------------------------------------------
C     Next two selections define the vector of the coupled nuclei. Note
C     that the order of atoms picked is important.
C     ------------------------------------------------------------------
      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_RDC',
     &                'More than 1 atom in SEL for atom M. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XDIPO_RDC','Error with atom selection')
          NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL)
      IF (INSERTPOS.EQ.1) ATOMMPTR(INSERTPOS)=0
      II=ATOMMPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XRDCMX) XRDCMX=II
      ATOMMLST(II)=SEL(1)
      ATOMMPTR(INSERTPOS+1)=II

      CALL SELCTA(SEL,NSEL,X,Y,Z,.TRUE.)
      IF (NSEL.GT.1) THEN
          CALL DSPERR('XDIPO_RDC',
     &                'More than 1 atom in SEL for atom N. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
          CALL DSPERR('XDIPO_RDC','Error with atom selection')
          NFLAG=1
      END IF
      CALL MAKIND(SEL,NATOM,NSEL)
      IF (INSERTPOS.EQ.1) ATOMNPTR(INSERTPOS)=0
      II=ATOMNPTR(INSERTPOS)
      II=II+1
      IF (II.GT.XRDCMX) XRDCMX=II
      ATOMNLST(II)=SEL(1)
      ATOMNPTR(INSERTPOS+1)=II
C     ------------------------------------------------------------------
C     Reading in the observed residual dipolar coupling
C     ------------------------------------------------------------------
      CALL NEXTF('Observed RDC =',XRDCO)
      CALL NEXTF('Error in RDC =',XRDCE)

      XRDCOBS(INSERTPOS)=XRDCO
      XRDCERR(INSERTPOS)=XRDCE
C     ------------------------------------------------------------------
C     Check for error atom selection. If there is one then reset the
C     counter for restraint
C     ------------------------------------------------------------------
      IF (NFLAG.EQ.1) THEN
          RDCNUM=RDCNUM-1
      END IF

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE XRDCINIT
C     ------------------------------------------------------------------
C     Initializes residual dipolar coupling stuff
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "xdipo_rdc.inc"

      CALL XRDCDEFAULTS
      CALL ALLOCXRDC(0,MAXXRDC)

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE PRINTXRDC (CUTOFF,XRDCCALC,XRDCOBS,XRDCERR,
     &                      ATOMIPTR,ATOMJPTR,ATOMKPTR,
     &                      ATOMLPTR,ATOMMPTR,ATOMNPTR,
     &                      ATOMILST,ATOMJLST,ATOMKLST,
     &                      ATOMLLST,ATOMMLST,ATOMNLST)
C     ------------------------------------------------------------------
C     Prints residual dipolar couplings with deviation larger than cutoff,
C     calculates RMSD and puts it into $RESULT
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "xdipo_rdc.inc"
      INCLUDE  "comand.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "numbers.inc"

      DOUBLE PRECISION CUTOFF,XRDCCALC(*),XRDCOBS(*),XRDCERR(*)
      INTEGER ATOMILST(*),ATOMJLST(*),ATOMKLST(*)
      INTEGER ATOMLLST(*),ATOMMLST(*),ATOMNLST(*)
      INTEGER ATOMIPTR(*),ATOMJPTR(*),ATOMKPTR(*)
      INTEGER ATOMLPTR(*),ATOMMPTR(*),ATOMNPTR(*)
C     ------------------------------------------------------------------
C     Local variables
C     ------------------------------------------------------------------
      DOUBLE PRECISION CALCXRDC,OBSXRDC,DELTAXRDC,DELTA,DP
      INTEGER COUNT,CLASS,I,J,K,L,M,N,NUM,II
      DOUBLE PRECISION RMS,VIOLS,ERRXRDC,XRDCENERGY
      DOUBLE COMPLEX DUMMY2
      LOGICAL PRINTTHISCLASS

      RMS=ZERO
      VIOLS=ZERO
      NUM=0
C     ------------------------------------------------------------------
C     Make sure that the array of calculated RDC is up to date
C     ------------------------------------------------------------------
      CALL EXRDC(XRDCENERGY,'ANALYZE')
      WRITE (PUNIT,'(A)') 'The following RDC''s have'
      WRITE (PUNIT,'(A)') 'deviation larger than or'
      WRITE (PUNIT,'(A)') 'equal to the cutoff:'
C     ------------------------------------------------------------------
C     Write out first class heading
C     ------------------------------------------------------------------
      CLASS=1
      PRINTTHISCLASS=PRINTCLASS(CLASS)
      IF (PRINTTHISCLASS) THEN
          WRITE (PUNIT,'(A,A)') 'Class ',XRDCCLASSNAMES(1)
      END IF
C     ------------------------------------------------------------------
C     For every residual dipolar couplig entryY...
C     ------------------------------------------------------------------
      COUNT=0
      DO WHILE (COUNT.LT.RDCNUM)
         COUNT=COUNT+1
C     ------------------------------------------------------------------
C     Is this the start of a new class?
C     ------------------------------------------------------------------
         IF (XRDCASSNDX(CLASS).LT.COUNT) THEN
             CLASS=CLASS+1
             PRINTTHISCLASS=PRINTCLASS(CLASS)
             IF (XRDCASSNDX(CLASS).EQ.XRDCASSNDX(CLASS-1)) THEN
                 COUNT=COUNT-1
             END IF
             IF (PRINTTHISCLASS) THEN
                 WRITE (PUNIT,'(A,A)') 'Class ',XRDCCLASSNAMES(CLASS)
             END IF
         END IF
C     ------------------------------------------------------------------
C     If this assignment is in a class to be printed
C     and make sure there is an entry for that class
C     ------------------------------------------------------------------
         IF ((PRINTTHISCLASS).AND.
     &       (XRDCASSNDX(CLASS).NE.XRDCASSNDX(CLASS-1))) THEN
C     ------------------------------------------------------------------
C     Always update RMSD
C     ------------------------------------------------------------------
              CALCXRDC=XRDCCALC(COUNT)
              OBSXRDC=XRDCOBS(COUNT)
              ERRXRDC=XRDCERR(COUNT)
              DP=(CALCXRDC-OBSXRDC)

              IF ((DP.LT.0.000).AND.(ABS(DP).GT.ERRXRDC)) THEN
                   DELTAXRDC=DP+ERRXRDC
              ELSE IF ((DP.GT.0.000).AND.(DP.GT.ERRXRDC)) THEN
                   DELTAXRDC=DP-ERRXRDC
              ELSE
                   DELTAXRDC=0.0
              END IF

              RMS=RMS+DELTAXRDC**2
              NUM=NUM+1
C     ------------------------------------------------------------------
C     Print out deviations larger than cutoff
C     and update number of violations
C     ------------------------------------------------------------------
              IF (ABS(DELTAXRDC).GE.CUTOFF) THEN
                  I=ATOMILST(ATOMIPTR(COUNT)+1)
                  J=ATOMJLST(ATOMJPTR(COUNT)+1)
                  K=ATOMKLST(ATOMKPTR(COUNT)+1)
                  L=ATOMLLST(ATOMLPTR(COUNT)+1)
                  WRITE(PUNIT,'(A,A)') '==============================',
     &                                 '=============================='
                  WRITE(PUNIT,'(A)') 'Set-M-atoms'
                  DO II=ATOMMPTR(COUNT)+1,ATOMMPTR(COUNT+1)
                     M=ATOMMLST(II)
                     WRITE(PUNIT,'(9X,4(1X,A))') SEGID(M),RESID(M),
     &                                           RES(M),TYPE(M)
                  END DO
                  WRITE(PUNIT,'(A)') 'Set-N-atoms'
                  DO II=ATOMNPTR(COUNT)+1,ATOMNPTR(COUNT+1)
                     N=ATOMNLST(II)
                     WRITE(PUNIT,'(9X,4(1X,A))') SEGID(N),RESID(N),
     &                                           RES(N),TYPE(N)
                  END DO

                  WRITE(PUNIT,'(2(2X,A,1X,F8.3))')
     &                  'Calc: ',CALCXRDC,'Obs: ',OBSXRDC
                  WRITE(PUNIT,'(2X,A,1X,F8.3,2X,A,1X,F8.3)')
     &                  'Error: ',ERRXRDC,'Delta: ',DELTAXRDC
                  VIOLS=VIOLS+ONE
              END IF
         END IF
      END DO

      IF (NUM.GT.0) THEN
          RMS=SQRT(RMS / NUM)
      ELSE
          RMS=0.0
      ENDIF

      CALL DECLAR('RESULT','DP',' ',DUMMY2,RMS)
      CALL DECLAR('VIOLATIONS','DP',' ',DUMMY2,VIOLS)

      RETURN
      END

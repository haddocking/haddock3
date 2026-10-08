      SUBROUTINE NOESET
C
C Sets up NOE restraint force field for dynamics and minimization.
C Also performs elementary analysis to determine possible NOE's for
C a given structure.
C
C Author: Axel Brunger
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
C local
      INTEGER ISLCT, JSLCT
C begin
      ISLCT=ALLHP(INTEG4(NATOM))
      JSLCT=ALLHP(INTEG4(NATOM))
      CALL NOESE2(HEAP(ISLCT),HEAP(JSLCT))
      CALL FREHP(JSLCT,INTEG4(NATOM))
      CALL FREHP(ISLCT,INTEG4(NATOM))
      RETURN
      END
C====================================================================
      SUBROUTINE NOESE2(ISLCT,JSLCT)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INTEGER ISLCT(*), JSLCT(*)
C local
      INTEGER I, II, NISLCT, NJSLCT, ICL1, ICL2, NN, ITEMP, IEMODE
      CHARACTER*4 CLASS, CLASS2, STEMP, SRESET, DENMODE
      DOUBLE PRECISION SCALE, RTEMP, NOEDST, THRESH, GAMMA, KAPPA
C AB added to print the number of NOE restraints
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
      LOGICAL MATCH
C parameter
      DOUBLE PRECISION ZERO, ONE, SIX
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, SIX=6.0D0)
C begin
C
C allocate space
      IF (NOEMAX.EQ.0) CALL NOEHP(200)
C all defaults are set in NOERES!
      CALL PUSEND('NOE>')
      DO WHILE (.NOT. DONE)
      CALL NEXTWD('NOE>')
      CALL MISCOM('NOE>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-noe')
C
C====================================================================
      ELSE IF (WD(1:4).EQ.'NRES') THEN
      NN=NOEMAX
      CALL NEXTI('NREStraints=',NN)
      IF (NN.NE.NOEMAX) THEN
      WRITE(6,'(A,I7,A)')
     & ' NOE: allocating space for ',NN,' restraints.'
      CALL NOEHP(NN)
      END IF
C====================================================================
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      CALL NOERES
C====================================================================
      ELSE IF (WD(1:4).EQ.'TEMP') THEN
      CALL NEXTF('TEMPerature=',NOETEM)
C====================================================================
      ELSE IF (WD(1:4).EQ.'CEIL') THEN
      CALL NEXTF('CEILing=',NOECEI)
C====================================================================
      ELSE IF (WD(1:4).EQ.'CLAS') THEN
      IF (NOECCN.EQ.1.AND.NOECNM(1).EQ.'NIL'.AND.NOENUM.EQ.0) THEN
      ELSE IF (NOECCN+1.GE.NOECMX) THEN
      CALL WRNDIE(-5,'NOESET',
     & 'exceeded NOECMX parameter --> recompile program')
      ELSE
      NOECCN=NOECCN+1
      END IF
      CALL NEXTA4('CLASsification=',NOECNM(NOECCN))
      DO I=1,NOECCN-1
      CALL EQSTWC(NOECNM(I),4,NOECNM(NOECCN),4,1,1,MATCH)
      IF (MATCH) THEN
      CALL WRNDIE(+5,'NOESET',
     & 'duplicate class name specified. ')
      END IF
      END DO
C
C set the default values for the class specifications
      NOESCA(NOECCN)=1.0D0
      NOEAVE(NOECCN)=NOEM6
      NOEPOT(NOECCN)=NOEBIH
      NOEEXP(NOECCN)=2
      NOECNS(NOECCN)=1.0D0
      NOEOFF(NOECCN)=0.0D0
      NOEMOF(NOECCN)=0.0D0
      NOESWI(NOECCN)=10.0D0
      NOEASY(NOECCN)=0.0D0
      NOESOE(NOECCN)=2
      NOEMSW(NOECCN)=10.0D0
      NOEMAS(NOECCN)=0.0D0
      NOEMSO(NOECCN)=2
      NMONO(NOECCN)=1
      NAVEXP(NOECCN)=SIX
      OREXP(NOECCN)=SIX
C
C reset time and running-average stuff
      DAVEXP(NOECCN)=3
      RAVEXP(NOECCN)=3
      MEMTAU(NOECCN)=1000.
      TMODE(NOECCN)='OFF'
      RMODE(NOECCN)='OFF'
      TFMODE(NOECCN)='CONS'
C
C initialize high-dim and 3-d stuff
C   noecor is number of noes in this class for highdim only
C   noebar sets hight of barrier in highdim pot.

      NOECOR(NOECCN)=2
      NOEBAR(NOECCN)=0.1
      IF (NOECCN.EQ.1.AND.NOENUM.GT.0) THEN
      WRITE(6,'(A)')
     &' %NOE-ERR: NIL-classified ASSIgn entries preceed this statement'
      END IF
C====================================================================
      ELSE IF (WD(1:4).EQ.'ASSI') THEN
      CALL NOEASS(NOENUM,NOEORN,NOEMAX,NOECCN,NOESMX,HEAP(HPNORR),
     &     HEAP(HPNIPR),
     &     HEAP(HPNILS),HEAP(HPNJPR),HEAP(HPNJLS),HEAP(HPNCND),
     &     HEAP(HPNRAV),HEAP(HPNRRV),HEAP(HPNNSP),HEAP(HPNDIS),
     &     HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNWGH),HEAP(HPNCV),
     &     HEAP(HPNVIO),
     &     NPEAKI,HEAP(HPNPID),HEAP(HPNPP1),HEAP(HPNPP2),
     &     HEAP(HPNVOL),NATOM,ISLCT,NISLCT,
     &     JSLCT,NJSLCT,X,Y,Z,HEAP(HPDENINIT))
C===================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(A,I6,A,I6,A,/,A,F12.3,A,I8)')
     & ' NOE: total number of restraints:',NOENUM,
     & ' partitioned into ',NOECCN,' classes',
     & ' NOE: ceiling=',NOECEI,' current allocation=',NOEMAX
      IF (NOEICV.GT.0) THEN
      WRITE(6,'(A,/,A,I6)')
     & ' NOE: data are partitioned into working set and test set.',
     & ' NOE: test set number=',NOEICV
      END IF
      DBPREC=NOENUM
      CALL DECLAR( 'NUMNOE', 'DP', ' ', DBCOMP, DBPREC )
C===================================================================
      ELSE IF (WD(1:4).EQ.'SCAL') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('SCALE=',SCALE)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOESCA(I)=SCALE
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'AVEX') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
CCC modification ATB 4/27/08
      CALL NEXTF('AVEXponent=',RTEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      NAVEXP(I)=RTEMP
      END IF
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'OREX') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
CCC modification ATB 4/27/08
      CALL NEXTF('OREXponent=',RTEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      OREXP(I)=RTEMP
      END IF
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'AVER') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTA4('AVERaging=',STEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      IF (STEMP.EQ.'R-6') THEN
      WRITE(6,'(A)')
     & 'R-6 distance-averaging.'
      NOEAVE(I)=NOEM6
      ELSE IF (STEMP.EQ.'R-3') THEN
      NOEAVE(I)=NOEM3
      ELSE IF (STEMP.EQ.'SUM') THEN
      NOEAVE(I)=NOESUM
      ELSE IF (STEMP.EQ.'CENT') THEN
      NOEAVE(I)=NOECEN
      ELSE
      CALL DSPERR('NOE','unknown averaging option')
      END IF
      END IF
      END DO
      IF (STEMP.EQ.'CENT'.AND.EMODE.EQ.'ON') THEN
      WRITE(6,'(A)')
     &'center averaging not compatible with ensemble-averaging',
     &'==> ensemble averaging turned OFF'
      EMODE='OFF'
      ENDIF
C====================================================================
      ELSE IF (WD(1:4).EQ.'MONO') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTI('MONOmers=',NN)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NMONO(I)=NN
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'TAVE') THEN
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
C defaults
C
C parsing
      CALL PUSEND('TAVErage>')
      DO WHILE (.NOT. DONE)
      CALL NEXTWD('TAVErage>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-noe-taverage')
C
      ELSE IF (WD(1:3).EQ.'TAU') THEN
      CALL NEXTI('TAU [steps]=',NN)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      MEMTAU(I)=NN
      END IF
      END DO
C
      ELSE IF (WD(1:3).EQ.'EXP') THEN
      CALL NEXTI('EXPOnent=',NN)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      DAVEXP(I)=NN
      END IF
      END DO
C
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      CALL NEXTWD('RESEt=')


      IF (WD(1:4).EQ.'CURR') THEN
      SRESET='CURR'
      ELSEIF (WD(1:4).EQ.'CONS') THEN
      SRESET='CONS'
      ELSE
      CALL DSPERR('taverage','unknown reset option.')
      END IF
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      CALL RESETN(I,HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNDIS),SRESET)
      END IF
      END DO
C
      ELSE IF (WD(1:2).EQ.'ON') THEN
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      TMODE(I)='ON'
      END IF
      END DO
C
      ELSE IF (WD(1:3).EQ.'OFF') THEN
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      TMODE(I)='OFF'
      END IF
      END DO
C
      ELSE IF (WD(1:4).EQ.'STAT') THEN
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      TMODE(I)='STAT'
      END IF
      END DO
C
      ELSE IF (WD(1:4).EQ.'FORC') THEN
      CALL NEXTWD('FORCe=')
      IF (WD(1:4).EQ.'CONS') THEN
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      TFMODE(I)='CONS'
      END IF
      END DO
      ELSEIF (WD(1:4).EQ.'NONC') THEN
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      TFMODE(I)='NONC'
      END IF
      END DO
      ELSE
      CALL DSPERR('taverage','unknown force option.')
      END IF
C
      ELSE
      CALL CHKEND('TAVErage>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
C====================================================================
      ELSE IF (WD(1:4).EQ.'RAVE') THEN
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
C defaults
C
C parsing
      CALL PUSEND('RAVErage>')
      DO WHILE (.NOT. DONE)
      CALL NEXTWD('RAVErage>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-noe-raverage')
C
      ELSE IF (WD(1:3).EQ.'EXP') THEN
      CALL NEXTI('EXPOnent=',NN)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      RAVEXP(I)=NN
      END IF
      END DO
C
      ELSE IF (WD(1:2).EQ.'ON') THEN
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      RMODE(I)='ON'
      END IF
      END DO
C
      ELSE IF (WD(1:3).EQ.'OFF') THEN
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      RMODE(I)='OFF'
      END IF
      END DO
C
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      CALL NEXTWD('RESEt=')
      IF (WD(1:4).EQ.'CURR') THEN
      SRESET='CURR'
      ELSEIF (WD(1:4).EQ.'CONS') THEN
      SRESET='CONS'
      ELSE
      CALL DSPERR('raverage','unknown reset option.')
      END IF
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      CALL RESETN(I,HEAP(HPNCND),HEAP(HPNRRV),HEAP(HPNDIS),SRESET)
      CALL RESETC(I,HEAP(HPNNSP),HEAP(HPNCND))
      END IF
      END DO
C
      ELSE
      CALL CHKEND('RAVErage>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
C====================================================================
      ELSE IF (WD(1:4).EQ.'POTE') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTA4('POTEntial=',STEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) THEN
      IF (STEMP.EQ.'BIHA') THEN
      NOEPOT(I)=NOEBIH
      ELSE IF (STEMP.EQ.'SQUA') THEN
      NOEPOT(I)=NOESQW
      ELSE IF (STEMP.EQ.'SOFT') THEN
      NOEPOT(I)=NOESOF
      ELSE IF (STEMP.EQ.'SYMM') THEN
      NOEPOT(I)=NOESYM
      ELSE IF (STEMP.EQ.'HIGH') THEN
      NOEPOT(I)=NOEHDI
      ELSE IF (STEMP.EQ.'LOGN') THEN
      NOEPOT(I)=NOELOG
      ELSE IF (STEMP.EQ.'3DPO') THEN
      NOEPOT(I)=NOE3DP
      ELSE
      CALL DSPERR('NOE','unknown potential option')
      END IF
      END IF
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'SQCO') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('SQUAre-constant=',RTEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOECNS(I)=RTEMP
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'SQEX') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTI('SQEXponent=',II)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOEEXP(I)=II
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'SQOF') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('SQOFfset=',RTEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOEOFF(I)=RTEMP
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'SQMO') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('SQMOFfset=',RTEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOEMOF(I)=RTEMP
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'SOEX') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTI('SOEXponent=',II)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOESOE(I)=II
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'RSWI') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('RSWItch=',RTEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOESWI(I)=RTEMP
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'ASYM') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('ASYMptote=',RTEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOEASY(I)=RTEMP
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'MSOE') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTI('SOEXponent=',II)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOEMSO(I)=II
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'MRSW') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('RSWItch=',RTEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOEMSW(I)=RTEMP
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'MASY') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('ASYMptote=',RTEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOEMAS(I)=RTEMP
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'NCOU') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTI('NCOUnt=',II)
      IF (II.GT.20) THEN
      CALL WRNDIE(-5,'NOE>','maximum value NCOUnt is 20. ')
      END IF
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOECOR(I)=II
      END DO
C====================================================================
      ELSE IF (WD(1:4).EQ.'BHIG') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTF('BHIGht=',RTEMP)
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) NOEBAR(I)=RTEMP
      END DO

C====================================================================
      ELSE IF (WD(1:4).EQ.'DIST') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      CALL NEXTA4('class-name=',CLASS2)
      CALL NEXTF('switch-value=',NOEDST)
      ICL1=0
      ICL2=0
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) ICL1=I
      CALL EQSTWC(NOECNM(I),4,CLASS2,4,1,1,MATCH)
      IF (MATCH) ICL2=I
      END DO
      IF ((ICL1*ICL2).GT.0) THEN
      CALL DISNOE(ICL1,ICL2,NOEDST,HEAP(HPNCND))

      END IF
C====================================================================
      ELSE IF (WD(1:4).EQ.'COUN') THEN
C
C get the class name:
      CALL NEXTA4('class-name=',CLASS)
      ICL1=0
      DO I=1,NOECCN
      CALL EQSTWC(NOECNM(I),4,CLASS,4,1,1,MATCH)
      IF (MATCH) ICL1=I
      END DO
      IF (ICL1.GT.0) THEN
      CALL COUNTV(ICL1,HEAP(HPNCND))
      END IF
C====================================================================
      ELSE IF (WD(1:4).EQ.'PRIN') THEN
      CALL NEXTWD('PRINt>')
      IF (WD(1:4).NE.'THRE') THEN
      WRITE(6,'(A)') ' %NOE-ERR: print expects THREshold parameter'
      ELSE
      CALL NEXTF('THREshold=',THRESH)
      CALL NOEPRI(THRESH,0)
      IF (NOEICV.GT.0) CALL NOEPRI(THRESH,1)
      END IF
C====================================================================
CCC modification ATB 5/02/08
      ELSE IF (WD(1:4).EQ.'OUTD') THEN
      CALL NOEOUT(0)
      IF (NOEICV.GT.0) CALL NOEOUT(1)
C====================================================================
CCC modification ATB 4/25/08
      ELSE IF (WD(1:3).EQ.'DEN') THEN
      CALL NEXTWD('DEN>')
      IF (WD(1:4).NE.'INIT'.AND.WD(1:4).NE.'UPDA') THEN
      WRITE(6,'(A)') ' %NOE-ERR: DEN expects INIT or UPDAte'
      ELSE
      DENMODE=WD(1:4)
      GAMMA=0.0D0
      KAPPA=0.0D0
      IF (DENMODE.EQ.'UPDA') THEN
      CALL NEXTWD('DEN>')
      IF (WD(1:4).NE.'GAMM') THEN
      WRITE(6,'(A)') ' %NOE-ERR: DEN UPDAte expects GAMMa parameter'
      ELSE
      CALL NEXTF('GAMMa=',GAMMA)
      CALL NEXTWD('DEN>')
      IF (WD(1:4).NE.'KAPP') THEN
      WRITE(6,'(A)') ' %NOE-ERR: DEN UPDAte expects KAPPa parameter'
      ELSE
      CALL NEXTF('KAPPa=',KAPPA)
      END IF
      END IF
      END IF 
      CALL NOEDEN(DENMODE,GAMMA,KAPPA,0)
      IF (NOEICV.GT.0) CALL NOEDEN(DENMODE,GAMMA,KAPPA,1)
      END IF
C====================================================================
      ELSE IF (WD(1:4).EQ.'ANAL') THEN
      CALL NEXTA4('ANALyse=',RANA)
C====================================================================
      ELSE IF (WD(1:2).EQ.'CV') THEN
      CALL NEXTI('CV excluded partition number:',NOEICV)
C====================================================================
      ELSE IF (WD(1:4).EQ.'PART') THEN
      CALL NEXTI('number of PARTitions:',ITEMP)
      ITEMP=MAX(0,ITEMP)
      CALL NOECVS(NOENUM,HEAP(HPNCV),ITEMP)
      IF (ITEMP.EQ.0) THEN
      NOEICV=0
      END IF
C====================================================================
      ELSE IF (WD(1:4).EQ.'PRED') THEN
      CALL NOEPRD(ISLCT,JSLCT)
C====================================================================
      ELSE IF (WD(1:4).EQ.'ENSE') THEN
C
C parsing
      CALL PUSEND('ENSEmble>')
      DO WHILE (.NOT. DONE)
      CALL NEXTWD('ENSEmble>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-noe-ensemble')
C
      ELSE IF (WD(1:2).EQ.'ON') THEN
      EMODE='ON'
      WRITE(6,'(A,A3)') 'ensemble-averaging turned ',EMODE
      IEMODE = 0
      DO I=1,NOECCN
      IF (NOEAVE(I).EQ.NOECEN) THEN
      IEMODE = 1
      END IF
      END DO
      IF (IEMODE.EQ.1) THEN
      WRITE(6,'(A)')
     & ' center averaging not compatible with ensemble-averaging',
     & ' ==> ensemble-averaging turned OFF'
      EMODE='OFF'
      END IF
C
      ELSE IF (WD(1:3).EQ.'OFF') THEN
      EMODE='OFF'
      DIENS=' '
      IDIMER=0
      WRITE(6,'(A,A3)') 'ensemble-averaging turned ',EMODE
C
      ELSE IF (WD(1:4).EQ.'DIME') THEN
      DIENS='DIME'
      IDIMER=IDIMER+1
      IF (IDIMER.GT.MAXMOD) THEN
      WRITE(6,'(A)')
     & ' maximum number of allowed interactions for ',
     & ' ensemble-averaged restraints exceeded'
      CALL WRNDIE(-1,'ENSE',
     & ' exceeded MAXMOD parameter --> recompile program')
      END IF
      CALL NEXTA4('segid #1=',ENSMOD(IDIMER,1))
      CALL NEXTA4('segid #2=',ENSMOD(IDIMER,2))
      WRITE(6,'(A,A4,A,A4)')
     & ' DIMEr NOE restraints allowed between segid''s ',
     & ENSMOD(IDIMER,1),' and ',ENSMOD(IDIMER,2)
C
      ELSE IF (WD(1:4).EQ.'OCCU') THEN
      CALL NEXTLO('OCCUpancy=', OCCU)
      IF (OCCU) THEN
      WRITE(6,'(A)') 'occupancy refinement turned ON'
      ELSE
      WRITE(6,'(A)') 'occupancy refinement turned OFF'
      END IF
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      IDIMER = 0
      DIENS = ' '
C
      ELSE
      CALL CHKEND('ENSEmble>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      ELSE
      CALL CHKEND('NOE>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
      RETURN
      END
C====================================================================
      SUBROUTINE NOEASS(NOENUM,NOEORN,NOEMAX,NOECCN,NOESMX,NOEORR,
     &           NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,NOERAV,NOERRV,
     &           NOESTP,NOEDIS,NOELOW,NOEHIG,NOEWGH,NOECV,NOEVIO,
     &           NPEAKI,
     &           NOEPID,NOEPP1,NOEPP2,NOEVOL,
     &           NATOM,ISLCT,NISLCT,
     &           JSLCT,NJSLCT,X,Y,Z,NOEDENINIT)
C
C Subroutine parses an NOE ASSIgn statement and
C puts the information into the NOE restraints list
C
C Author: Axel T. Brunger
C         extension to handle "or" restraints by Michael Nilges
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'timer.inc'
C I/O
      INTEGER NOENUM, NOEORN, NOEMAX, NOECCN, NOESMX, NOEORR(*)
      INTEGER NOEIPR(*), NOEILS(*), NOEJPR(*), NOEJLS(*), NOECND(*)
      DOUBLE PRECISION NOERAV(*), NOERRV(*), NOEDIS(*)
      DOUBLE PRECISION NOELOW(*), NOEHIG(*), NOEWGH(*)
      INTEGER NOECV(*), NOEVIO(*), NOESTP(*), NPEAKI, NOEPID(*)
      DOUBLE PRECISION NOEPP1(*),NOEPP2(*),NOEVOL(*)
      INTEGER NATOM, ISLCT(*), NISLCT, JSLCT(*), NJSLCT
      DOUBLE PRECISION X(*), Y(*), Z(*), NOEDENINIT(*)
C local
      INTEGER I, II, JJ, J, TISLCT, TJSLCT
      LOGICAL SUCCES,OK
C begin
      NPEAKI=NPEAKI+1
      SUCCES=.FALSE.
      TISLCT=0
      TJSLCT=0
C
      IF (NOEORN.GE.NOEMAX-1) THEN
      WRITE(6,'(A)')
     & ' %NOE-ERR: allocation for NOE-restraints exceeded',
     & '              increase NRES parameter and run again'
      CALL WRNDIE(-1,'NOE>',
     & 'exceeded allocation for NOE-restraints')
      ELSE
      NOENUM=NOENUM+1
      END IF
C
      IF(NOENUM.EQ.1) NOEORR(NOENUM)=0
C
C parse fixed format portion of NOE restraint
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(JSLCT,NJSLCT,X,Y,Z,.TRUE.)
      CALL NEXTF('Distance=',NOEDIS(NOENUM))
      CALL NEXTF('lower-average-distance=',NOELOW(NOENUM))
      CALL NEXTF('higher-average-distance=',NOEHIG(NOENUM))
C
C initialize initial DEN distance
      NOEDENINIT(NOENUM)=ZERO

C
C initialize time-average and running average stuff
C
C set the time-averaged distance and
C running-average term to the NOE constraint distance
      NOERAV(NOENUM)=NOEDIS(NOENUM)
      NOERRV(NOENUM)=NOEDIS(NOENUM)
C
C Number of times the running-average term has been calculated
      NOESTP(NOENUM)=1
C
C assign the index into the class list
      NOECND(NOENUM)=NOECCN
C include the distance in the working set
C Default value for NOECV is -1 (needed to allow partitioning of
C single classes only) AB April 2004
      NOECV(NOENUM)=-1
C initialize number of violations
      NOEVIO(NOENUM)=0
C initialize peak id
      NOEPID(NOENUM)=NPEAKI
C initialize frequencies
      NOEPP1(NOENUM)=-999.0D0
      NOEPP2(NOENUM)=-999.0D0
C initialize volume
      NOEVOL(NOENUM)=ZERO
C initialize Iweight
      NOEWGH(NOENUM)=ONE
C
C set up the I and J lists
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1.AND.JSLCT(I).EQ.1) THEN
      WRITE(6,'(A)') ' %NOESET-ERR: i and j set have to be disjoint'
      NISLCT=0
      NJSLCT=0
      END IF
      END DO
C
      IF (NISLCT.EQ.0.OR.NJSLCT.EQ.0) THEN
C
C modified 12/02/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' %NOESET-INFO: One of the two selections has no atoms.'
      END IF
      ELSE
      SUCCES=.TRUE.
      NOEORN=NOEORN+1
C
C update maximum number of protons in each group
      NOESMX=MAX(NOESMX,2*NISLCT)
      NOESMX=MAX(NOESMX,2*NJSLCT)
C
      CALL MAKIND(ISLCT,NATOM,NISLCT)
      CALL MAKIND(JSLCT,NATOM,NJSLCT)
C
C now setup the I-list
      IF (NOENUM.EQ.1) NOEIPR(NOENUM)=0
      II=NOEIPR(NOEORN)
      DO I=1,NISLCT
      IF (II.GE.NOEMAX-1) THEN
      WRITE(6,'(A)')
     & ' %NOE-ERR: allocation for NOE-restraints exceeded',
     & '              increase NRES parameter and run again'
      CALL WRNDIE(-1,'NOE',
     & 'exceeded allocation for NOE-restraints')
      ELSE
      II=II+1
      END IF
      NOEILS(II)=ISLCT(I)
      END DO
      NOEIPR(NOEORN+1)=II
C
C now setup the J-list
      IF (NOENUM.EQ.1) NOEJPR(NOENUM)=0
      JJ=NOEJPR(NOEORN)
      DO J=1,NJSLCT
      IF (JJ.GE.NOEMAX-1) THEN
      WRITE(6,'(A)')
     & ' %NOE-ERR: allocation for NOE-restraints exceeded',
     & '              increase NRES parameter and run again'
      CALL WRNDIE(-1,'NOE',
     & 'exceeded allocation for NOE-restraints')
      ELSE
      JJ=JJ+1
      END IF
      NOEJLS(JJ)=JSLCT(J)
      END DO
      NOEJPR(NOEORN+1)=JJ
      END IF
C
C default parameters are parsed, check for additional keywords
      OK=.FALSE.
      DO WHILE (.NOT. OK)
      CALL NEXTWD('NOE>')
      IF (WD(1:4).EQ.'PEAK') THEN
      CALL NEXTI('PEAKid=',NPEAKI)
      ELSE IF (WD(1:4).EQ.'PPM1') THEN
      CALL NEXTF('PPM(F1)=',NOEPP1(NOENUM))
      ELSE IF (WD(1:4).EQ.'PPM2') THEN
      CALL NEXTF('PPM(F2)=',NOEPP2(NOENUM))
      ELSE IF (WD(1:4).EQ.'VOLU') THEN
      CALL NEXTF('VOLUME=',NOEVOL(NOENUM))
      ELSE IF (WD(1:4).EQ.'WEIG') THEN
      CALL NEXTF('WEIGHT=',NOEWGH(NOENUM))
      ELSE IF (WD(1:4).EQ.'CV  ') THEN
      CALL NEXTI('CV=',NOECV(NOENUM))
      ELSE IF (WD(1:4).EQ.'OR  ') THEN
C
C add a restraint to the NOEORR array
      IF ((NOEORN.GE.NOEMAX-1)) THEN
      WRITE(6,'(A)')
     & ' %NOE-ERR: allocation for NOE-restraints exceeded',
     & '              increase NRES parameter and run again'
      CALL WRNDIE(-1,'NOE>',
     & 'exceeded allocation for NOE-restraints')
      ELSE
      NOEORN=NOEORN+1
      END IF
      TISLCT=TISLCT+NISLCT
      TJSLCT=TJSLCT+NJSLCT
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(JSLCT,NJSLCT,X,Y,Z,.TRUE.)
C
C the rest is very similar to the normal parsing of selections above.
C set up the I and J lists
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1.AND.JSLCT(I).EQ.1) THEN
      WRITE(6,'(A)') ' %NOESET-ERR: i and j set have to be disjoint'
      NISLCT=0
      NJSLCT=0
      END IF
      END DO
C
      IF (NISLCT.EQ.0.OR.NJSLCT.EQ.0) THEN
C
C modified 12/02/08, ATB
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     &  ' %NOESET-INFO: One of the two selections has no atoms.'
      END IF
      NOEORN=NOEORN-1
      ELSE
      SUCCES=.TRUE.
C
C update maximum number of protons in each group
      TISLCT=TISLCT+NISLCT
      TJSLCT=TJSLCT+NJSLCT
      NOESMX=MAX(NOESMX,2*TISLCT)
      NOESMX=MAX(NOESMX,2*TJSLCT)
C
      CALL MAKIND(ISLCT,NATOM,NISLCT)
      CALL MAKIND(JSLCT,NATOM,NJSLCT)
C
C now setup the I-list
      II=NOEIPR(NOEORN)
      DO I=1,NISLCT
      IF (II.GE.NOEMAX-1) THEN
      WRITE(6,'(A)')
     & ' %NOE-ERR: allocation for NOE-restraints exceeded',
     & '              increase NRES parameter and run again'
      CALL WRNDIE(-1,'NOE',
     & 'exceeded allocation for NOE-restraints')
      ELSE
      II=II+1
      END IF
      NOEILS(II)=ISLCT(I)
      END DO
      NOEIPR(NOEORN+1)=II
C
C now setup the J-list
      JJ=NOEJPR(NOEORN)
      DO J=1,NJSLCT
      IF (JJ.GE.NOEMAX-1) THEN
      WRITE(6,'(A)')
     & ' %NOE-ERR: allocation for NOE-restraints exceeded',
     & '              increase NRES parameter and run again'
      CALL WRNDIE(-1,'NOE',
     & 'exceeded allocation for NOE-restraints')
      ELSE
      JJ=JJ+1
      END IF
      NOEJLS(JJ)=JSLCT(J)
      END DO
      NOEJPR(NOEORN+1)=JJ
      END IF

      ELSE
      CALL SAVEWD
      OK=.TRUE.
      END IF
      END DO
C
      IF (.NOT.SUCCES) THEN
C
C modified, ATB 12/02/08
CCC      WRITE (6,'(A,I5,2F10.3)') ' %NOE-ERR: problem at ',
CCC     &    NPEAKI, NOEPP1(NOENUM), NOEPP2(NOENUM)
      NOENUM=NOENUM-1
C we do not need to reduce NOEORN since it is only advanced if SUCCES
      ELSE
      NOEPID(NOENUM)=NPEAKI
      NOEORR(NOENUM+1)=NOEORN
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE NOEPRD(ISLCT,JSLCT)
C
C Parsing routine for NOE PREDict
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'consta.inc'
      INTEGER ISLCT(*), JSLCT(*)
C local
      EXTERNAL EXCHR
      INTEGER MBOUNDS, I, UNIT
      PARAMETER (MBOUNDS=10)
      DOUBLE PRECISION CUTON, CUTOFF, BOUNDS(MBOUNDS), SCRAMBL
      DOUBLE PRECISION BSYMM, BUNSYM
      INTEGER NISLCT, NJSLCT, NBOUNDS
      LOGICAL UNIQUE, BYPASS
      CHARACTER*4 SMODE
C parameter
      DOUBLE PRECISION ZERO, TEN
      PARAMETER (ZERO=0.0D0, TEN=10.0D0)
C begin
C defaults
      SMODE='TOPO'
C topology mode
      CUTON=ZERO
      CUTOFF=TEN
      CALL FILL4(ISLCT,NATOM,1)
      CALL FILL4(JSLCT,NATOM,1)
C assignment mode
      NBOUNDS=1
      BOUNDS(1)=ZERO
      SCRAMBL=ZERO
      BSYMM=ZERO
      BUNSYM=ZERO
      UNIQUE=.FALSE.
      BYPASS=.FALSE.
      OFILE=' '
C
C parsing
      CALL PUSEND('PREDict>')
      DO WHILE (.NOT. DONE)
      CALL NEXTWD('PREDict>')
      CALL MISCOM('PREDict>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-noe-predict')
C
      ELSE IF (WD(1:4).EQ.'TO  ') THEN
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      ELSE IF (WD(1:4).EQ.'FROM') THEN
      CALL SELCTA(JSLCT,NJSLCT,X,Y,Z,.TRUE.)
      ELSE IF (WD(1:5).EQ.'CUTON') THEN
      CALL NEXTF('CUTON=',CUTON)
      ELSE IF (WD(1:5).EQ.'CUTOF') THEN
      CALL NEXTF('CUTOFf=',CUTOFF)
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      CALL NEXTA4('MODE=',SMODE)
      ELSE IF (WD(1:5).EQ.'BSYMM') THEN
      CALL NEXTF('BSYMMetry=',BSYMM)
      ELSE IF (WD(1:6).EQ.'BUNSYM') THEN
      CALL NEXTF('BUNSYMetry=',BUNSYM)
      ELSE IF (WD(1:4).EQ.'BOUN') THEN
      IF (NBOUNDS.GE.MBOUNDS) THEN
      CALL WRNDIE(-5,'NOEPRD>',
     &    'maximum number of bounds (MBOUNDS) exceeded.')
      ELSE
      NBOUNDS=NBOUNDS+1
      END IF
      CALL NEXTF('BOUNdary=',BOUNDS(NBOUNDS))
C
C the 0 bound was already specified as default for NBOUNDS=1
      IF (ABS(BOUNDS(NBOUNDS)).LT.RSMALL) THEN
      NBOUNDS=NBOUNDS-1
      END IF
      ELSE IF (WD(1:4).EQ.'SCRA') THEN
      CALL NEXTF('SCRAmble=',SCRAMBL)
      ELSE IF (WD(1:4).EQ.'UNIQ') THEN
      CALL NEXTLO('UNIQue=',UNIQUE)
      ELSE IF (WD(1:4).EQ.'BYPA') THEN
      CALL NEXTLO('BYPAss=',BYPASS)
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE
      CALL CHKEND('PREDict>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (SMODE.EQ.'TOPO') THEN
C
C topology mode
      IF (NOENUM.GT.0) THEN
      WRITE(6,'(A)')
     & ' NOEPRD-warning: NOE restraints database erased.'
      END IF
C
C allocate sufficient space
      CALL NOEHP(MAX(NISLCT,NJSLCT))
      CALL NOEPRE(ISLCT,JSLCT,CUTON,CUTOFF,HEAP(HPNORR),
     &     HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),HEAP(HPNJLS),
     &     HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),HEAP(HPNNSP),
     &     HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNWGH),
     &     HEAP(HPNCV))
C
C reset the NOE database
      CALL NOERES
C
      ELSE
C use existing assignments for the predictions
C
      IF (OFILE.NE.' ') THEN
      CALL ASSFIL(OFILE,UNIT,'WRITE','FORMATTED',ERROR)
      ELSE
      UNIT=0
      END IF
C
C sort the bounds
      CALL SORT(NBOUNDS,EXCHR,ORDERR,BOUNDS,1,0,0,0,0,0,0)
      WRITE(6,'(A)') ' The following bounds will be used:'
      WRITE(6,'(A,10F8.3)') ' ',(BOUNDS(I),I=1,NBOUNDS)
      WRITE(6,'(A,F8.3,A,L1,A,L1,A,F8.3,A,F8.3)')
     & ' SCRAmbling=',SCRAMBL,'  UNIQue=',UNIQUE,'  BYPAss=',BYPASS,
     & ' BSYMMetry_error=',BSYMM,' BUNSYMetry_error=',BUNSYM
C
      CALL NOEPD2(HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),
     &     HEAP(HPNJPR),HEAP(HPNJLS),HEAP(HPNCND),
     &     HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),
     &     HEAP(HPNCV),
     &     NBOUNDS,BOUNDS,SCRAMBL,BSYMM,BUNSYM,UNIQUE,BYPASS,
     &     UNIT)
      END IF
C
      RETURN
      END
C====================================================================
      SUBROUTINE NOEPD2(NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,
     &           NOECND,NOEDIS,NOELOW,NOEHIG,
     &           NOECV,NBOUNDS,BOUNDS,SCRAMBL,BSYMM,BUNSYM,
     &           UNIQUE,BYPASS,PUNIT)
C
C Routine "predicts" NOEs based on current cartesian coordinates.
C Distance ranges are declared according to the boundaries in
C array BOUNDS.  It is assumed that BOUNDS is sorted.  The
C boundaries define distances ranges by
C lower-bound=0
C upper-bound=BOUNDS (i)
C central-distance=BOUNDS(i)/2.
C
C Assigned distance restraints (NOEs) are pointed to the distance range
C that satisfies lower-bound < R < upper-bound where R is
C computed with the specified averaging option.
C
C If scramble is specified, distance restraints are randomly
C assigned to one range or the other if they are within SCRAMble
C around a boundary.
C
C If unique is specified (i.e., true) then only
C the distance restraints are "touched" with unique
C assignments.
C
C If BYPAss is specified then the NOEs are untouched.
C
C Only NOEs in the working set will be predicted.
C
C If the percentage NOE distance error (BSYMM,BUNSYM) is greater
C than zero then boundaries will be calculated for the
C "predicted NOEs.  There are two distances ranges, BUNSYM and
C BSYMM, for unsymmetrical boundaries and symmetrical boundaries.
C For BUNSYM
C lower boundary= (distance * BUNSYM) or
C                 [(distance - (distance * BUNSYM))>1.6]
C upper boundary= (distance * BUNSYM)
C For BSYMM
C lower boundary= (distance * BUNSYM)
C upper boundary= (distance * BUNSYM)
C
C Author: Axel T. Brunger, Gregory L. Warren
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'consta.inc'
      INTEGER NOEORR(*), NOEIPR(*), NOEILS(*), NOEJPR(*)
      INTEGER NOEJLS(*), NOECND(*)
      DOUBLE PRECISION NOEDIS(*), NOELOW(*), NOEHIG(*)
      INTEGER NOECV(*), NBOUNDS
      DOUBLE PRECISION BOUNDS(*), SCRAMBL
      DOUBLE PRECISION BSYMM, BUNSYM
      LOGICAL UNIQUE, BYPASS
      INTEGER PUNIT
C local
      INTEGER N, I, K, J, CLASS, NB , KLAST
      DOUBLE PRECISION OLDDIS, OLDLOW, OLDHIG, DISTANC, RNUM, SMALLST
      DOUBLE PRECISION NEWLB, NEWUB, EN, DIFF, LOWLIM
      LOGICAL QFIRST, QUNCHG
C parameter
      DOUBLE PRECISION ZERO, TWO, HALF
      PARAMETER (ZERO=0.0D0, TWO=2.0D0, HALF=0.5D0, LOWLIM=1.6D0)
C begin
      QUNCHG=.TRUE.
C
C loop over all classes
      DO CLASS=1,NOECCN
      IF (NOEPOT(CLASS).EQ.NOESYM.OR.
     &    NOEPOT(CLASS).EQ.NOEHDI.OR.
     &    NOEPOT(CLASS).EQ.NOE3DP) THEN
      CALL WRNDIE(-5,'NOEPD2',
     & 'PREDict not possible with POTEntial=SYMM, HIGH, or 3DPO')
      ELSE
C
C loop over all NOE's belonging to this class
      DO N=1,NOENUM
      IF (NOECV(N).NE.NOEICV.AND.NOECND(N).EQ.CLASS) THEN
C
C store old values
      OLDDIS=NOEDIS(N)
      OLDLOW=NOELOW(N)
      OLDHIG=NOEHIG(N)
C
C get distance (is returned in PCDATA(PCGEOM))
      CALL ENOE(EN,'ANAL',N)
      DISTANC=PCDATA(PCGEOM)
C
C check if UNIQUE or BYPASS is specified
      IF ((.NOT.BYPASS).AND.(.NOT.UNIQUE.OR.
     &    (NOEORR(N+1)-NOEORR(N).EQ.1.AND.
     &     NOEIPR(N+1)-(NOEIPR(N)+1).EQ.0 .AND.
     &     NOEJPR(N+1)-(NOEJPR(N)+1).EQ.0))) THEN
C
C generate the distance ranges using the percent NOE
C distance error value
      IF ((BSYMM.GT.ZERO).OR.(BUNSYM.GT.ZERO)) THEN
C test for low limit of distance
        IF (DISTANC.LT.LOWLIM) THEN
          DISTANC=LOWLIM
        END IF
C Function used generate symmetrical (BSYMM)
C  or unsymmetrical (BUNSYM) bounds
        IF (BSYMM.GT.ZERO) THEN
          NOEDIS(N)=DISTANC
          NOEHIG(N)=(DISTANC*BSYMM)
          NOELOW(N)=(DISTANC*BSYMM)
        ELSE
        NOEDIS(N)=DISTANC
        NOEHIG(N)=DISTANC*BUNSYM
        NOELOW(N)=DISTANC-(DISTANC*BUNSYM)
        IF (NOELOW(N).LT.LOWLIM) THEN
            NOELOW(N)=DISTANC-LOWLIM
        ELSE
            NOELOW(N)=DISTANC*BUNSYM
        END IF
        END IF
        QUNCHG=.FALSE.
      ELSE
C or find nearest boundary
      NB=1
      SMALLST=R4BIG
      DO I=1,NBOUNDS
      DIFF=DISTANC-BOUNDS(I)
      IF (ABS(DIFF).LT.ABS(SMALLST)) THEN
      SMALLST=DIFF
      NB=I
      END IF
      END DO
C
C scramble boundary assignment if required
      IF (ABS(SMALLST).LT.SCRAMBL) THEN
      CALL GGUBFS(RNUM)
      IF (RNUM.GE.HALF) THEN
      NB=MIN(NBOUNDS,NB+1)
      END IF
      ELSE
C
C no scrambling -- assign to appropriate range
C
C distance range is 0 < r < BOUNDS(NB+1) if DIFF is greater than
C BOUNDS(NB).
      IF (SMALLST.GE.ZERO) THEN
      NB=MIN(NBOUNDS,NB+1)
      END IF
      END IF
C
C generate the distance bounds for the restraint
      NEWLB=ZERO
      NEWUB=BOUNDS(NB)
      NOEDIS(N)=(NEWUB-NEWLB)/TWO
      NOELOW(N)=NOEDIS(N)-NEWLB
      NOEHIG(N)=NEWUB-NOEDIS(N)
C
      QUNCHG=.FALSE.
C
      END IF
C
C distance restraint is unchanged
      QUNCHG=.TRUE.
      END IF
C
      IF (PUNIT.GT.0) THEN
      WRITE(PUNIT,'(A)') ' '
      QFIRST=.TRUE.
      KLAST=NOEILS(NOEIPR(N+1))
      DO I=NOEIPR(N)+1,NOEIPR(N+1)
      K=NOEILS(I)
      IF (QFIRST) THEN
      QFIRST=.FALSE.
      IF (K.NE.KLAST) THEN
      WRITE(PUNIT,'(7A)') ' ASSIgn ((segid "',SEGID(K),
     & '" and resid ',RESID(K),' and name ',TYPE(K),')'
      ELSE
      WRITE(PUNIT,'(7A)') ' ASSIgn ((segid "',SEGID(K),
     & '" and resid ',RESID(K),' and name ',TYPE(K),'))'
      END IF
      ELSE
      IF (K.NE.KLAST) THEN
      WRITE(PUNIT,'(7A)') '      OR (segid "',SEGID(K),
     & '" and resid ',RESID(K),' and name ',TYPE(K),')'
      ELSE
      WRITE(PUNIT,'(7A)') '      OR (segid "',SEGID(K),
     & '" and resid ',RESID(K),' and name ',TYPE(K),'))'
      END IF
      END IF
      END DO
      QFIRST=.TRUE.
      KLAST=NOEJLS(NOEJPR(N+1))
      DO J=NOEJPR(N)+1,NOEJPR(N+1)
      K=NOEJLS(J)
      IF (QFIRST) THEN
      QFIRST=.FALSE.
      IF (K.NE.KLAST) THEN
      WRITE(PUNIT,'(7A)') '        ((segid "',SEGID(K),
     & '" and resid ',RESID(K),' and name ',TYPE(K),')'
      ELSE
      WRITE(PUNIT,'(7A)') '        ((segid "',SEGID(K),
     & '" and resid ',RESID(K),' and name ',TYPE(K),'))'
      END IF
      ELSE
      IF (K.NE.KLAST) THEN
      WRITE(PUNIT,'(7A)') '      OR (segid "',SEGID(K),
     & '" and resid ',RESID(K),' and name ',TYPE(K),')'
      ELSE
      WRITE(PUNIT,'(7A)') '      OR (segid "',SEGID(K),
     & '" and resid ',RESID(K),' and name ',TYPE(K),'))'
      END IF
      END IF
      END DO
      IF (QUNCHG) THEN
      WRITE(PUNIT,'(A,3F10.4,A,F10.4,A)')
     & '    ',NOEDIS(N),NOELOW(N),NOEHIG(N),
     & ' { unchanged } { R=',DISTANC,'}'
      ELSE
      WRITE(PUNIT,'(A,3F10.4,A,3F10.4,A)')
     & '    ',NOEDIS(N),NOELOW(N),NOEHIG(N),
     & ' { previous: ',OLDDIS,OLDLOW,OLDHIG,' }'
      END IF
      END IF
C
      END IF
      END DO
      END IF
      END DO
C
      RETURN
      END
C====================================================================
CMNTODO: implement noeorr array, for simplest case
      SUBROUTINE NOEPRE(ISLCT,JSLCT,CUTON,CUTOFF,
     &           NOEORR,NOEIPR,NOEILS,
     &           NOEJPR,NOEJLS,NOECND,NOERAV,NOERRV,NOESTP,NOEDIS,
     &           NOELOW,NOEHIG,NOEWGH,NOECV)
C
C predicts NOE distances
C Author: Axel T. Brunger
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'numbers.inc'
      INTEGER ISLCT(*), JSLCT(*)
      DOUBLE PRECISION CUTON, CUTOFF
      INTEGER NOEORR(*),NOEIPR(*), NOEILS(*), NOEJPR(*)
      INTEGER NOEJLS(*), NOECND(*)
      DOUBLE PRECISION NOERAV(*), NOERRV(*), NOEDIS(*)
      DOUBLE PRECISION NOELOW(*), NOEHIG(*), NOEWGH(*)
      INTEGER NOESTP(*), NOECV(*)
C local
      INTEGER IG, JG, IAT, JAT, IGSLCT, JGSLCT, KGSLCT, LGSLCT, K
      DOUBLE PRECISION EN
      INTEGER DDX, DDY, DDZ
C
C begin
C
C quick-and-dirty fix for uninitialized NOEWGH (RWGK/PDA 07/99)
      DO K=1,NOEMAX
        NOEWGH(K)=ONE
      END DO
C
C loop over group-group pairs
      DO IG=1,NGRP
      IGSLCT=0
      KGSLCT=0
      DO IAT=IGPBS(IG)+1,IGPBS(IG+1)
      IGSLCT=IGSLCT+ISLCT(IAT)
      KGSLCT=KGSLCT+JSLCT(IAT)
      END DO
      IF (IGSLCT.GT.0) THEN
      DO JG=1,NGRP
      JGSLCT=0
      LGSLCT=0
      DO JAT=IGPBS(JG)+1,IGPBS(JG+1)
      JGSLCT=JGSLCT+JSLCT(JAT)
      LGSLCT=LGSLCT+ISLCT(JAT)
      END DO
      IF (JGSLCT.GT.0) THEN
C
C avoid duplications for group pairs
      IF (.NOT.(KGSLCT.GT.0.AND.LGSLCT.GT.0.AND.IG.GE.JG)) THEN
C
C now fill a dummy NOE table
      NOENUM=1
      NOEORR(1)=0
      NOEORR(2)=1
      NOEIPR(1)=0
      NOEIPR(2)=0
      DO IAT=IGPBS(IG)+1,IGPBS(IG+1)
      IF (ISLCT(IAT).GT.0) THEN
      NOEIPR(2)=NOEIPR(2)+1
      NOEILS(NOEIPR(2))=IAT
      END IF
      END DO
      NOEJPR(1)=0
      NOEJPR(2)=0
      DO JAT=IGPBS(JG)+1,IGPBS(JG+1)
      IF (JSLCT(JAT).GT.0) THEN
      NOEJPR(2)=NOEJPR(2)+1
      NOEJLS(NOEJPR(2))=JAT
      END IF
      END DO
      NOELOW(1)=1.0D0
      NOEHIG(1)=1.0D0
      NOECND(1)=1
      NOESCA(1)=1.0D0
      NOEAVE(1)=NOEM6
      NOEPOT(1)=NOEBIH
      NOEDIS(1)=1.0D0
C
C symmetry stuff
      NMONO(1)=1
C
C reset time and running-average stuff
      NOESTP(1)=0
      DAVEXP(1)=3
      RAVEXP(1)=3
      MEMTAU(1)=1000.
      TMODE(1)='OFF'
      RMODE(1)='OFF'
      TFMODE(1)='CONS'
C
C update maximum number of protons in each group
      NOESMX=MAX(NOEIPR(2),NOESMX)
      NOESMX=MAX(NOEJPR(2),NOESMX)
      DDX=ALLHP(IREAL8(NOESMX*NOESMX))
      DDY=ALLHP(IREAL8(NOESMX*NOESMX))
      DDZ=ALLHP(IREAL8(NOESMX*NOESMX))
      CALL ENOE2(EN,'ANAL',1,HEAP(DDX),HEAP(DDY),HEAP(DDZ),
     &     NOEORR,NOEIPR,NOEILS,
     &     NOEJPR,NOEJLS,NOECND,
     &     NOERAV,NOERRV,NOESTP,NOEDIS,
     &     NOELOW,NOEHIG,NOEWGH,NOECV)
      CALL FREHP(DDX,IREAL8(NOESMX*NOESMX))
      CALL FREHP(DDY,IREAL8(NOESMX*NOESMX))
      CALL FREHP(DDZ,IREAL8(NOESMX*NOESMX))
      IF (PCDATA(PCGEOM).GT.CUTON.AND.PCDATA(PCGEOM).LT.CUTOFF) THEN
      WRITE(PUNIT,'(A)') ' ========================================'
      WRITE(PUNIT,'(A)') ' set-i-atoms'
      DO IAT=NOEIPR(1)+1,NOEIPR(2)
      K=NOEILS(IAT)
      WRITE(PUNIT,'(9X,4(1X,A))') SEGID(K),RESID(K),RES(K),TYPE(K)
      END DO
      WRITE(PUNIT,'(A)') ' set-j-atoms'
      DO JAT=NOEJPR(1)+1,NOEJPR(2)
      K=NOEJLS(JAT)
      WRITE(PUNIT,'(9X,4(1X,A))') SEGID(K),RESID(K),RES(K),TYPE(K)
      END DO
      WRITE(PUNIT,'(A,F8.3)')
     &  ' ( ( <R^-6> )^-1/6 )        =        ',PCDATA(PCGEOM)
      END IF
C
      END IF
      END IF
      END DO
      END IF
      END DO
C
      RETURN
      END
C====================================================================
      SUBROUTINE NOEPRI(THRESH,ITEST)
C
C print a list of current NOE restraints
C Front-end for NOEPR2
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'noe.inc'
      DOUBLE PRECISION THRESH
      INTEGER ITEST
C begin
      CALL NOEPR2(THRESH,ITEST,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),
     &           HEAP(HPNVIO),HEAP(HPNCV), HEAP(HPNPID))

      RETURN
      END
C====================================================================
      SUBROUTINE NOEPR2(THRESH,ITEST,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,
     &           NOEJLS,NOECND,NOERAV,NOERRV,
     &           NOEDIS,NOELOW,NOEHIG,
     &           NOEVIO,NOECV,NOEPID)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'numbers.inc'
      DOUBLE PRECISION THRESH
      INTEGER ITEST
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOECV(*), NOEPID(*)
C
C local
      DOUBLE PRECISION EN, NOERMS, NOERM2(NOECMX)
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
      INTEGER NOENU2(NOECMX)
      INTEGER NOEVIT
      INTEGER N, K, I, J, CLASS, NOENUML, NORR
      CHARACTER*6 SSAVE
      CHARACTER*11 SSPOT
      LOGICAL QHEAD
C begin
      NOERMS=ZERO
      NOEVIT=0
      NOENUML=0
      IF (NOENUM.GT.0) THEN
C
      IF (NOEICV.GT.0) THEN
      IF (ITEST.EQ.0) THEN
      WRITE(PUNIT,'(A)')

     & ' $$$$$$$$$$$$$$$$$$ working set $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      ELSE
      WRITE(PUNIT,'(A,I5,A)')

     & ' $$$$$$$$$$$$$$$$$$$$ test set (TEST=',NOEICV,
     & ')  $$$$$$$$$$$$$$$$$$$$$$'
      END IF
      END IF
C
C loop over all classes

      DO CLASS=1,NOECCN
      NOENU2(CLASS)=0
      NOERM2(CLASS)=ZERO
      NOEVIO(CLASS)=0
      IF (NOEPOT(CLASS).EQ.NOESYM.OR.
     &    NOEPOT(CLASS).EQ.NOEHDI.OR.
     &    NOEPOT(CLASS).EQ.NOE3DP) THEN
      CALL WRNDIE(-5,'NOEPRI',
     & 'PRINt not possible with POTEntial set to SYMM, HIGH, or 3DPO')
      ELSE
C
C loop over all NOE's belonging to this class
      QHEAD=.TRUE.
      DO N=1,NOENUM
      IF ((ITEST.EQ.0.AND.NOECV(N).NE.NOEICV).OR.
     &    (ITEST.EQ.1.AND.NOECV(N).EQ.NOEICV)) THEN
      IF (NOECND(N).EQ.CLASS) THEN
      NOENUML=NOENUML+1
C
C get energy and other results
      CALL ENOE(EN,'ANAL',N)
      NOERMS=NOERMS+PCDATA(PCDEVI)**2
      NOERM2(CLASS)=NOERM2(CLASS)+PCDATA(PCDEVI)**2
      NOENU2(CLASS)=NOENU2(CLASS)+1
      IF (ABS(PCDATA(PCDEVI)).GT.THRESH) THEN
      NOEVIO(CLASS)=NOEVIO(CLASS)+1
      NOEVIT=NOEVIT+1
      END IF
C
C if deviation greater than threshold then print the NOE
      IF (ABS(PCDATA(PCDEVI)).GT.THRESH) THEN
C
C check whether we've to print the header for this class
      IF (QHEAD) THEN
      QHEAD=.FALSE.
      WRITE(PUNIT,'(/2A/,3A)')
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',
     & '+++++++++++++++++++',
     & ' +++++++++++++++++++++++++++++++ CLASS ',NOECNM(CLASS),
     & ' +++++++++++++++++++++++++++++++++++'
      IF (NOEAVE(NOECND(N)).EQ.NOEM6) THEN
      SSAVE='R-6   '
      ELSE IF (NOEAVE(NOECND(N)).EQ.NOEM3) THEN
      SSAVE='R-3   '
      ELSE IF (NOEAVE(NOECND(N)).EQ.NOESUM) THEN
      SSAVE='sum'
      ELSE IF (NOEAVE(NOECND(N)).EQ.NOECEN) THEN
      SSAVE='center'
      END IF
      IF (NOEPOT(NOECND(N)).EQ.NOEBIH) THEN
      SSPOT='biharmonic '
      ELSE IF (NOEPOT(NOECND(N)).EQ.NOESQW) THEN
      SSPOT='square-well'
      ELSE IF (NOEPOT(NOECND(N)).EQ.NOESOF) THEN
      SSPOT='soft-square'
      ELSE IF (NOEPOT(NOECND(N)).EQ.NOELOG) THEN
      SSPOT='log-normal'
      END IF
      WRITE(PUNIT,'(A,F8.3,4A)')
     & ' for this class: SCALe=',NOESCA(NOECND(N)),
     & ' AVERage=',SSAVE,' POTEntial=',SSPOT
      IF (NOEPOT(NOECND(N)).EQ.NOEBIH) THEN
      WRITE(PUNIT,'(A,F8.3)')
     & '                 TEMPerature=',NOETEM
      ELSE IF (NOEPOT(NOECND(N)).EQ.NOESQW) THEN
      WRITE(PUNIT,'(A,F8.3,A,I5,A,2F8.3)')
     & '                 SQCOnstant=',NOECNS(NOECND(N)),
     & ' SQEXponent=',NOEEXP(NOECND(N)),
     & ' SQOFfsets(+/-)=',NOEOFF(NOECND(N)), NOEMOF(NOECND(N))
      ELSE IF (NOEPOT(NOECND(N)).EQ.NOESOF) THEN
      WRITE(PUNIT,'(A,F8.3,A,I5,A,2F8.3/A,I5,A,F8.3,A,F8.3)')
     & '                 SQCOnstant=',NOECNS(NOECND(N)),
     & ' SQEXponent=',NOEEXP(NOECND(N)),
     & ' SQOFfsets(+/-)=',NOEOFF(NOECND(N)), NOEMOF(NOECND(N)),
     & '                 SOEXponent=',NOESOE(NOECND(N)),
     & ' RSWItch=',NOESWI(NOECND(N)),
     & ' ASYMptote=',NOEASY(NOECND(N))
      ELSE IF (NOEPOT(NOECND(N)).EQ.NOELOG) THEN
      WRITE(PUNIT,'(A,F8.3)')
     &        '          SQCOnstant=',NOECNS(NOECND(N))
      END IF
      WRITE(PUNIT,'(A)') ' '
      END IF
C
      WRITE(PUNIT,'(A,I5,A)')
     & ' ========== restraint ',NOEPID(N),' =========='
      DO NORR=NOEORR(N)+1, NOEORR(N+1)
      WRITE(PUNIT,'(A)') ' set-i-atoms'
      DO I=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      K=NOEILS(I)
      WRITE(PUNIT,'(9X,4(1X,A))') SEGID(K),RESID(K),RES(K),TYPE(K)
      END DO
      WRITE(PUNIT,'(A)') ' set-j-atoms'
      DO J=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      K=NOEJLS(J)
      WRITE(PUNIT,'(9X,4(1X,A))') SEGID(K),RESID(K),RES(K),TYPE(K)
      END DO
      END DO
      END IF
C
C
      IF (ABS(PCDATA(PCDEVI)).GT.THRESH) THEN
C
      WRITE(PUNIT,'(A,F8.3,A,F5.2,A,F5.2,A,F5.2,2A,F8.3,A,F8.3)')
     &  ' R<average>=',PCDATA(PCGEOM),
     &  ' NOE=',NOEDIS(N),' (-',NOELOW(N),'/+',
     &    NOEHIG(N),')',
     &  ' Delta=',PCDATA(PCDEVI),
     & '  E(NOE)=',PCDATA(PCENER)
      IF (RMODE(NOECND(N)).NE.'OFF') THEN
      IF (RANA.EQ.'RAVE') THEN
      WRITE(PUNIT,'(A,F8.3,A)')
     &  ' running-averaged <R>=',NOERRV(N),
     &  ' (used in energy calculation)'
      ELSE
      WRITE(PUNIT,'(A,F8.3)')
     &  ' running-averaged <R>=',NOERRV(N)
      END IF
      END IF
      IF (TMODE(NOECND(N)).NE.'OFF') THEN
      IF (RANA.EQ.'TAVE') THEN
      WRITE(PUNIT,'(A,F8.3,A)')
     &  ' time-averaged <R>=',NOERAV(N),
     &  ' (used in energy calculation)'
      ELSE
      WRITE(PUNIT,'(A,F8.3)')
     &  ' time-averaged <R>=',NOERAV(N)
      END IF
      END IF
      END IF
      END IF
      END IF
      END DO
      END IF
      END DO
C
      WRITE(PUNIT,'(/A,F8.3,A,F4.1,A,I6,A,I6,A)')
     & ' NOEPRI: RMS diff. =',SQRT(NOERMS/NOENUML),
     & ',  #(violat.>',THRESH,')=',NOEVIT,' of ',NOENUML,' NOEs'
      IF (ITEST.EQ.0) THEN
      DBPREC=NOEVIT
      CALL DECLAR( 'VIOLATIONS', 'DP', ' ', DBCOMP, DBPREC )
      DBPREC =SQRT(NOERMS/NOENUML)
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
      CALL DECLAR( 'RMS', 'DP', ' ', DBCOMP, DBPREC )
      DBPREC=NOENUML
      CALL DECLAR( 'NUMBER', 'DP', ' ', DBCOMP, DBPREC )
      ELSE
      DBPREC=NOEVIT
      CALL DECLAR( 'TEST_VIOLATIONS', 'DP', ' ', DBCOMP, DBPREC )
      DBPREC =SQRT(NOERMS/NOENUML)
      CALL DECLAR( 'TEST_RMS', 'DP', ' ', DBCOMP, DBPREC )
      END IF
      DO I=1,NOECCN
      IF (NOENU2(I).GT.0) THEN
      WRITE(PUNIT,'(3A,F8.3,A,F4.1,A,I6,A,I6,A)')
     & ' NOEPRI: RMS diff. class ',
     & NOECNM(I),' =',SQRT(NOERM2(I)/NOENU2(I)),
     & ',  #(viol.>',THRESH,
     & ')=',NOEVIO(I),' of ',NOENU2(I),' NOEs'
      END IF
      END DO
      END IF
      RETURN
      END
C====================================================================
      SUBROUTINE NOEOUT(ITEST)
C
C print a list of current NOE restraints
C as assign statements
C Front-end for NOEPR3
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'noe.inc'
      INTEGER ITEST
C begin
      CALL NOEPR3(ITEST,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),
     &           HEAP(HPNVIO),HEAP(HPNCV), HEAP(HPNPID))

      RETURN
      END
C====================================================================
      SUBROUTINE NOEPR3(ITEST,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,
     &           NOEJLS,NOECND,NOERAV,NOERRV,
     &           NOEDIS,NOELOW,NOEHIG,
     &           NOEVIO,NOECV,NOEPID)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'numbers.inc'
      INTEGER ITEST
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOECV(*), NOEPID(*)
C
C local
      INTEGER N, K, I, J, CLASS, NORR
      INTEGER I1, I2
C begin
      IF (NOENUM.GT.0) THEN
C
C loop over all classes
      DO CLASS=1,NOECCN
      IF (NOEPOT(CLASS).EQ.NOESYM.OR.
     &    NOEPOT(CLASS).EQ.NOEHDI.OR.
     &    NOEPOT(CLASS).EQ.NOE3DP) THEN
      CALL WRNDIE(-5,'NOEPRI',
     & 'DENOut not possible with POTEntial set to SYMM, HIGH, or 3DPO')
      ELSE
C
C loop over all NOE's belonging to this class
      DO N=1,NOENUM
      IF ((ITEST.EQ.0.AND.NOECV(N).NE.NOEICV).OR.
     &    (ITEST.EQ.1.AND.NOECV(N).EQ.NOEICV)) THEN
      IF (NOECND(N).EQ.CLASS) THEN
C
C     WE ONLY TAKE THE FIRST FROM EACH SET
      NORR=NOEORR(N)+1
C
C     DO I=NOEIPR(NORR)+1,NOEIPR(NORR+1)
C
C     WE ONLY TAKE THE FIRST FROM EACH SET
      I=NOEIPR(NORR)+1
      K=NOEILS(I)
      I1=K
C
C     WE ONLY TAKE THE FIRST FROM EACH SET
      J=NOEJPR(NORR)+1
      K=NOEJLS(J)
      I2=K
C
CCC      WRITE(PUNIT,'(A,I6,A,I6,A,F6.3,A,F6.3,A,F6.3)') 
CCC     & 'ASSIGn (KNOWn AND ID ',
CCC     & I1, ' ) (KNOWn AND ID ', I2, ' ) ', NOEDIS(N), 
CCC     & ' ', NOELOW(N), ' ', NOEHIG(N)
C
C modified, ATB, 12/01/08
      WRITE(PUNIT,'(13A,F7.4,A,F6.3,A,F6.3)')
     &  ' ASSIgn ( KNOWn and ATOM "', SEGID(I1),
     &          '" "',RESID(I1),'" "',
     &  TYPE(I1),'" ) ( KNOWn and ATOM "',SEGID(I2),
     &          '" "',RESID(I2),'" "',
     &  TYPE(I2),'" ) ',NOEDIS(N), 
     & ' ', NOELOW(N), ' ', NOEHIG(N)
C
      END IF
      END IF
      END DO
      END IF
      END DO
      END IF
      RETURN
      END
C====================================================================
      SUBROUTINE NOEDEN(DENMODE,GAMMA,KAPPA,ITEST)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'noe.inc'
      CHARACTER*4 DENMODE
      DOUBLE PRECISION GAMMA, KAPPA
      INTEGER ITEST
C begin
      CALL NOEDEN2(DENMODE,GAMMA,KAPPA,ITEST,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),
     &           HEAP(HPNVIO),HEAP(HPNCV), HEAP(HPNPID),
     &           HEAP(HPDENINIT))

      RETURN
      END
C====================================================================
      SUBROUTINE NOEDEN2(DENMODE,GAMMA,KAPPA,ITEST,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,
     &           NOEJLS,NOECND,NOERAV,NOERRV,
     &           NOEDIS,NOELOW,NOEHIG,
     &           NOEVIO,NOECV,NOEPID,NOEDENINIT)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'numbers.inc'
      CHARACTER*4 DENMODE
      DOUBLE PRECISION GAMMA,KAPPA
      INTEGER ITEST
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOECV(*), NOEPID(*)
C initial den distance
      DOUBLE PRECISION NOEDENINIT(*)
C
C local
      DOUBLE PRECISION EN, NOERMS, NOERM2(NOECMX), THRESH
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
      INTEGER NOENU2(NOECMX)
      INTEGER NOEVIT
      INTEGER N, K, I, J, CLASS, NOENUML, NORR
      CHARACTER*6 SSAVE
      CHARACTER*11 SSPOT
      LOGICAL QHEAD
C begin
      THRESH=-0.1
      NOERMS=ZERO
      NOEVIT=0
      NOENUML=0
      IF (NOENUM.GT.0) THEN
C
C
C loop over all classes
      DO CLASS=1,NOECCN
      NOENU2(CLASS)=0
      NOERM2(CLASS)=ZERO
      NOEVIO(CLASS)=0
      IF (NOEPOT(CLASS).EQ.NOESYM.OR.
     &    NOEPOT(CLASS).EQ.NOEHDI.OR.
     &    NOEPOT(CLASS).EQ.NOE3DP) THEN
      CALL WRNDIE(-5,'NOEPRI',
     & 'DEN not possible with POTEntial set to SYMM, HIGH, or 3DPO')
      ELSE
C
C loop over all NOE's belonging to this class
      QHEAD=.TRUE.
      DO N=1,NOENUM
      IF ((ITEST.EQ.0.AND.NOECV(N).NE.NOEICV).OR.
     &    (ITEST.EQ.1.AND.NOECV(N).EQ.NOEICV)) THEN
      IF (NOECND(N).EQ.CLASS) THEN
      NOENUML=NOENUML+1
C
C get NOE energy, current distance (in PCDATA(PCGEOM), and other results
      CALL ENOE(EN,'ANAL',N)
      NOERMS=NOERMS+PCDATA(PCDEVI)**2
      NOERM2(CLASS)=NOERM2(CLASS)+PCDATA(PCDEVI)**2
      NOENU2(CLASS)=NOENU2(CLASS)+1
      IF (ABS(PCDATA(PCDEVI)).GT.THRESH) THEN
      NOEVIO(CLASS)=NOEVIO(CLASS)+1
      NOEVIT=NOEVIT+1
      END IF
C
      IF (DENMODE.EQ.'INIT') THEN
C
C initialization
      NOEDENINIT(N)=PCDATA(PCGEOM)
CCCCxxx      NOEDIS(N)=PCDATA(PCGEOM)
      ELSE
C
C update step
      NOEDIS(N)=NOEDIS(N)+KAPPA*(GAMMA*(PCDATA(PCGEOM)-NOEDIS(N))+ 
     &          (ONE-GAMMA)*(NOEDENINIT(N)-NOEDIS(N)) )       
      END IF
C
      END IF      
      END IF
      END DO
      END IF
      END DO
C
      END IF
      RETURN
      END
C====================================================================
      SUBROUTINE ENOE(EN,ANALYS,NA)
C
C Routine computes force field for NOE restraints
C
C Front-end routine for ENOE2
C
C Author: Axel Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'mtf.inc'
      DOUBLE PRECISION EN
      CHARACTER*4 ANALYS
      INTEGER NA
C local
      INTEGER DDX, DDY, DDZ
C begin
      DDX=ALLHP(IREAL8(NOESMX*NOESMX))
      DDY=ALLHP(IREAL8(NOESMX*NOESMX))
      DDZ=ALLHP(IREAL8(NOESMX*NOESMX))
      CALL ENOE2(EN,ANALYS,NA,HEAP(DDX),HEAP(DDY),HEAP(DDZ),
     &     HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),
     &     HEAP(HPNJPR),HEAP(HPNJLS),HEAP(HPNCND),
     &     HEAP(HPNRAV),HEAP(HPNRRV),HEAP(HPNNSP),HEAP(HPNDIS),
     &     HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNWGH),HEAP(HPNCV))
      CALL FREHP(DDX,IREAL8(NOESMX*NOESMX))
      CALL FREHP(DDY,IREAL8(NOESMX*NOESMX))
      CALL FREHP(DDZ,IREAL8(NOESMX*NOESMX))
      RETURN
      END
C====================================================================
      SUBROUTINE ENOE2(EN,ANALYS,NA,DDX,DDY,DDZ,
     &           NOEORR,NOEIPR,
     &           NOEILS,NOEJPR,NOEJLS,NOECND,NOERAV,NOERRV,NOESTP,
     &           NOEDIS,NOELOW,NOEHIG,NOEWGH,NOECV)
C
C Modified: Alexandre Bonvin
C           inclusion of constraint interactions for multimodel
C           refinement
C Modified: Michael Nilges
C           to allow ambiguous restraints using "OR"
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      DOUBLE PRECISION EN
      CHARACTER*4 ANALYS
      INTEGER NA
      DOUBLE PRECISION DDX(*), DDY(*), DDZ(*)
      INTEGER NOEORR(*),NOEIPR(*), NOEILS(*), NOEJPR(*)
      INTEGER NOEJLS(*), NOECND(*)
      DOUBLE PRECISION NOERAV(*), NOERRV(*), NOEDIS(*)
      DOUBLE PRECISION NOELOW(*), NOEHIG(*), NOEWGH(*)
      INTEGER NOESTP(*), NOECV(*)
C local
      DOUBLE PRECISION XI, YI, ZI, XJ, YJ, ZJ, WENS
      DOUBLE PRECISION RN, RNN, RN2, XIJ, YIJ, ZIJ, SIJ, RIJ
      DOUBLE PRECISION DF, DDF, C
      DOUBLE PRECISION EL, DEVI, LW
      DOUBLE PRECISION AAA, BBB, NEXP, NSOE, SAVEXP
      DOUBLE PRECISION LOGDEVI
      INTEGER I, J, N, NSTART, NSTOP, II, JJ, LI, LJ
      INTEGER IPRI, NPRIME, ICLASS
      INTEGER L1,L2,INORR,NNORR,NORR
      DOUBLE PRECISION XIJP(2), YIJP(2), ZIJP(2), RIJP(2)
      LOGICAL QAVER, FIRST
      DOUBLE PRECISION  VRN(20)
      INTEGER HN
C begin
C
      EN=ZERO
C
      IF (ANALYS.EQ.'ANAL') THEN
      NSTART=NA
      NSTOP=NA
      ELSE
      NSTART=1
      NSTOP=NOENUM
      END IF
C
      FIRST=.TRUE.
      HN=0
C
      DO N=NSTART,NSTOP
C
      IF (((NOECV(N).NE.NOEICV).OR.(ANALYS.EQ.'ANAL'))
     &    .AND.((NOEPOT(NOECND(N)).EQ.NOEHDI).OR.
     &         (NOEPOT(NOECND(N)).EQ.NOE3DP))) THEN
      CALL MULTID(N,HN,EN,RN,NOEORR,
     &           NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,NOEDIS,
     &           NOELOW,NOEHIG,NOEWGH,VRN)
C
      ELSEIF (((NOECV(N).NE.NOEICV).OR.(ANALYS.EQ.'ANAL').OR.
     &    (ANALYS.EQ.'ASSI')).AND.(NOEPOT(NOECND(N)).EQ.NOESYM)) THEN
C
C symmetry potential
C essentially, the whole code is duplicated to allow the computation
C of the NOE symmetry term. Two subsequent NOE restraints are constrained
C to have a specified difference equal.
C The input format of the NOE list is identical to that of
C the usual NOE terms, but the interpretation of the terms is different
C first restraint:   r => signed distance difference (a positive
C distance
C   difference means that the first distance is bigger than the second)
C   d- => -error; d+ => +error.
C second restraint:  r, d- and d+ are not used.
C (MN)
C
C check if this is the first or the second time round
      IF (FIRST) THEN
      FIRST=.FALSE.
C
C check if the following or previous is also NOESYM
      IF ((NOEPOT(NOECND(N+1)).NE.NOESYM)
     &     .AND.(ANALYS.NE.'ANAL')) THEN
      CALL WRNDIE(-5,'NOE',
     & 'error in SYMMETRY potential, check NOE table')
C
      ELSE
      ICLASS=NOECND(N)
C
C get the two distances
      DO IPRI=1,2
C
C calculate the indices properly with the new NOEORR stuff.
C nprime needs to be replaced by noeorr(nprime)+1
      NPRIME=N+IPRI-1
      II=NOEIPR(NOEORR(NPRIME)+1)+1-NOEIPR(NOEORR(NPRIME)+1+1)
      JJ=NOEJPR(NOEORR(NPRIME)+1)+1-NOEJPR(NOEORR(NPRIME)+1+1)
      QAVER=(II.NE.0.OR.JJ.NE.0)
C if more than one pair of atoms in constraint
      IF (QAVER) THEN
      LI=0
      XI=ZERO
      YI=ZERO
      ZI=ZERO
C ...loop over all pairs of atoms belonging to restraint N
      DO II=NOEIPR(NOEORR(NPRIME)+1)+1,NOEIPR(NOEORR(NPRIME)+1+1)
      I=NOEILS(II)
      XI=XI+X(I)
      YI=YI+Y(I)
      ZI=ZI+Z(I)
      LI=LI+1
      END DO
      XI=XI/LI
      YI=YI/LI
      ZI=ZI/LI
      LJ=0
      XJ=ZERO
      YJ=ZERO
      ZJ=ZERO
      DO JJ=NOEJPR(NOEORR(NPRIME)+1)+1,NOEJPR(NOEORR(NPRIME)+1+1)
      J=NOEJLS(JJ)
      XJ=XJ+X(J)
      YJ=YJ+Y(J)
      ZJ=ZJ+Z(J)
      LJ=LJ+1
      END DO
      XJ=XJ/LJ
      YJ=YJ/LJ
      ZJ=ZJ/LJ
      XIJP(IPRI)=XI-XJ
      YIJP(IPRI)=YI-YJ
      ZIJP(IPRI)=ZI-ZJ
      RIJP(IPRI)=MAX(RSMALL,
     &           SQRT(XIJP(IPRI)**2+YIJP(IPRI)**2+ZIJP(IPRI)**2))
C  if only one pair of atoms in constraint
      ELSE
      I=NOEILS(NOEIPR(NOEORR(NPRIME)+1)+1)
      J=NOEJLS(NOEJPR(NOEORR(NPRIME)+1)+1)
      XIJP(IPRI)=X(I)-X(J)
      YIJP(IPRI)=Y(I)-Y(J)
      ZIJP(IPRI)=Z(I)-Z(J)
      RIJP(IPRI)=XIJP(IPRI)**2+YIJP(IPRI)**2+ZIJP(IPRI)**2
      RIJP(IPRI)=MAX(RSMALL, SQRT(RIJP(IPRI)))
      ENDIF
C
      END DO
C
      RN=RIJP(1)-RIJP(2)
C
C do square well potential with soft asymptote
C
      C=MIN(NOEWGH(N)*NOESCA(ICLASS)*NOECNS(ICLASS),NOECEI)
C
      DEVI=MAX(ZERO,RN-NOEDIS(N)-
     &            MAX(ZERO,NOEHIG(N)-NOEOFF(ICLASS)))+
     &           MIN(ZERO,RN-NOEDIS(N)+
     &            MAX(ZERO,NOELOW(N)-NOEMOF(ICLASS)))
      IF ((DEVI.LE.NOESWI(ICLASS))
     &          .AND.(DEVI.GE.-NOESWI(ICLASS))) THEN
      EL=C*DEVI**NOEEXP(ICLASS)
      DF=NOEEXP(ICLASS)*C*DEVI**(NOEEXP(ICLASS)-1)
      ELSE
      AAA=(FLOAT(NOEEXP(ICLASS))
     &           /FLOAT(NOESOE(ICLASS))+1.0D0)
     &           *NOESWI(ICLASS)**(NOEEXP(ICLASS))
     &           -NOEASY(ICLASS)
     &           *(1.0D0+1.0D0/FLOAT(NOESOE(ICLASS)))
     &           *NOESWI(ICLASS)
      BBB=-FLOAT(NOEEXP(ICLASS))
     &           /FLOAT(NOESOE(ICLASS))
     &           *NOESWI(ICLASS)
     &           **(NOEEXP(ICLASS)+NOESOE(ICLASS))
     &           +NOEASY(ICLASS)/FLOAT(NOESOE(ICLASS))
     &            *NOESWI(ICLASS)**(NOESOE(ICLASS)+1)
      EL=C*(AAA+BBB/ABS(DEVI)**(NOESOE(ICLASS))
     &          +NOEASY(ICLASS)*ABS(DEVI))
      DF=C*(-NOESOE(ICLASS)*BBB/ABS(DEVI)**(NOESOE(ICLASS)+1)
     &          +NOEASY(ICLASS))*SIGN(1.0D0, DEVI)
      END IF
C
C accumulate energy
      EN=EN+EL
C
C accumulate derivatives
      DO IPRI=1,2
C
      NPRIME=N+IPRI-1
      LI=NOEIPR(NOEORR(NPRIME)+1+1)-NOEIPR(NOEORR(NPRIME)+1)
      LJ=NOEJPR(NOEORR(NPRIME)+1+1)-NOEJPR(NOEORR(NPRIME)+1)
      QAVER=(LI.GT.1.OR.LJ.GT.1)
C if more than one pair of atoms in constraint
      IF (QAVER) THEN
C
      DO II=NOEIPR(NOEORR(NPRIME)+1)+1,NOEIPR(NOEORR(NPRIME)+1+1)
      I=NOEILS(II)
      DX(I)=DX(I)+DF*XIJP(IPRI)/(RIJP(IPRI)*LI)
      DY(I)=DY(I)+DF*YIJP(IPRI)/(RIJP(IPRI)*LI)
      DZ(I)=DZ(I)+DF*ZIJP(IPRI)/(RIJP(IPRI)*LI)
      END DO
      DO JJ=NOEJPR(NOEORR(NPRIME)+1)+1,NOEJPR(NOEORR(NPRIME)+1+1)
      J=NOEJLS(JJ)
      DX(J)=DX(J)-DF*XIJP(IPRI)/(RIJP(IPRI)*LJ)
      DY(J)=DY(J)-DF*YIJP(IPRI)/(RIJP(IPRI)*LJ)
      DZ(J)=DZ(J)-DF*ZIJP(IPRI)/(RIJP(IPRI)*LJ)
      END DO
C if only one pair of atoms in constraint
      ELSE
      I=NOEILS(NOEIPR(NOEORR(NPRIME)+1)+1)
      J=NOEJLS(NOEJPR(NOEORR(NPRIME)+1)+1)
      DX(I)=DX(I)+DF*XIJP(IPRI)/RIJP(IPRI)
      DY(I)=DY(I)+DF*YIJP(IPRI)/RIJP(IPRI)
      DZ(I)=DZ(I)+DF*ZIJP(IPRI)/RIJP(IPRI)
      DX(J)=DX(J)-DF*XIJP(IPRI)/RIJP(IPRI)
      DY(J)=DY(J)-DF*YIJP(IPRI)/RIJP(IPRI)
      DZ(J)=DZ(J)-DF*ZIJP(IPRI)/RIJP(IPRI)
      ENDIF
C the following takes care of a -1 in the chain rule
      DF=-DF
      END DO
C
      END IF
      ELSE
      FIRST=.TRUE.
      END IF
C
      ELSEIF (((NOECV(N).NE.NOEICV).OR.(ANALYS.EQ.'ANAL').OR.
     &    (ANALYS.EQ.'ASSI')).AND.(NOEPOT(NOECND(N)).NE.NOESYM)) THEN
C
C mn_new we add an outer loop that allows calculation of multiple
C restraints connected by "or" statements
C
C is more than one pair of atoms in constraint ?
C noeipr(noeorr(N)+1) replaces noeipr(N)
      QAVER=.TRUE.
      NNORR=NOEORR(N+1)-NOEORR(N)
      IF (NNORR.EQ.1) THEN
         NORR=NOEORR(N)+1
         II=NOEIPR(NORR)+1-NOEIPR(NORR+1)
         JJ=NOEJPR(NORR)+1-NOEJPR(NORR+1)
         QAVER=(II.NE.0.OR.JJ.NE.0)
      END IF
C
      IF (.NOT.QAVER) THEN
         I=NOEILS(NOEIPR(NOEORR(N)+1)+1)
         J=NOEJLS(NOEJPR(NOEORR(N)+1)+1)
         RN=SQRT(MAX(((X(I)-X(J))*(X(I)-X(J)) +
     &                (Y(I)-Y(J))*(Y(I)-Y(J)) +
     &                (Z(I)-Z(J))*(Z(I)-Z(J))),RSMALL))
         DF=ONE
C we have more than one atom and/or more than one possible
      ELSE
         L1=0
         RN=ZERO
         INORR=0
         DO NORR=NOEORR(N)+1, NOEORR(N+1)
            INORR=INORR+1
C
            IF (NOEAVE(NOECND(N)).EQ.NOEM6
     &           .OR.(NOEAVE(NOECND(N)).EQ.NOEM3)
     &           .OR.(NOEAVE(NOECND(N)).EQ.NOESUM)) THEN
C do r-6/r-3 averaging
C ++++++++++++++++
C R = ( <R-6> ) -1/6
         IF (NOEAVE(NOECND(N)).EQ.NOEM3) THEN
          SAVEXP=THREE
         ELSEIF (NOEAVE(NOECND(N)).EQ.NOEM6) THEN
          SAVEXP=SIX
         ELSE
          SAVEXP=NAVEXP(NOECND(N))
         END IF
C
         L2=0
         LW=ZERO
         RN2=ZERO
C
C loop over all pairs of atoms belonging to restraint N
         DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
          DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
           I=NOEILS(II)
           J=NOEJLS(JJ)
C
C ensemble averaging:
C setup weight for each selections ij
C with value 0 if i and j does not belong to identical
C segids otherwise 1 or occupancy if OCCU = .true.
           IF (EMODE .EQ.'ON') THEN
              CALL ENSBLE(WENS,I,J)
           ELSE
            WENS = ONE
           END IF
C
           IF (WENS.GE.ZERO) THEN
            XIJ=X(I)-X(J)
            YIJ=Y(I)-Y(J)
            ZIJ=Z(I)-Z(J)
            SIJ=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
            RIJ=ONE/(SIJ**(SAVEXP/TWO))
            RN2=RN2+RIJ*WENS
C local derivative accumulation
            LW=LW+WENS
            L1=L1+1
            L2=L2+1
            RIJ=RIJ/SIJ
            DDX(L1)=RIJ*XIJ
            DDY(L1)=RIJ*YIJ
            DDZ(L1)=RIJ*ZIJ
           END IF
C
          END DO
         END DO
C
C for NOESUM we do not average (divide by L)
C if we have a multimeric protein, we have to divide by the
C number of molecules to get the scaling right.
C the exponent is (-OREXP(NOECND(N))-1+SAVEXP+1)
C if savexp=orexp, the derivative here is one as it should
C if NNORR=1, the use of orexp is avoided, and the distance is calculated
C the standard way.
         IF (NNORR.EQ.1) THEN
          IF (NOEAVE(NOECND(N)).EQ.NOESUM) THEN
           RN2=MAX((NMONO(NOECND(N))/RN2)**(ONE/SAVEXP),RSMALL)
           DDF=ONE/NMONO(NOECND(N))
          ELSE
           RN2=MAX((LW/RN2)**(ONE/SAVEXP),RSMALL)
           DDF=ONE/LW
          END IF
         ELSE
          IF (NOEAVE(NOECND(N)).EQ.NOESUM) THEN
           IF (SAVEXP.EQ.OREXP(NOECND(N))) THEN
            RN2=MAX(RN2/NMONO(NOECND(N)),RSMALL)
            DDF=ONE/NMONO(NOECND(N))
           ELSE
            RN2=MAX((NMONO(NOECND(N))/RN2)**(ONE/SAVEXP),RSMALL)
            DDF=RN2**(SAVEXP-OREXP(NOECND(N)))/NMONO(NOECND(N))
           END IF
          ELSE
           RN2=MAX((LW/RN2)**(ONE/SAVEXP),RSMALL)
           DDF=RN2**(SAVEXP-OREXP(NOECND(N)))/LW
          END IF
         END IF
C
C derivatives
         L1=L1-L2
         DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
          DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
           I=NOEILS(II)
           J=NOEJLS(JJ)
C
C ensemble averaging:
C setup weight for each selections ij
C with value 0 if i and j does not belong to identical
C segids otherwise 1 or occupancy if OCCU = .true.
           IF (EMODE .EQ.'ON') THEN
              CALL ENSBLE(WENS,I,J)
           ELSE
            WENS = ONE
           END IF
C
           IF (WENS.GE.ZERO) THEN
             L1=L1+1
             DDX(L1)=DDX(L1)*DDF
             DDY(L1)=DDY(L1)*DDF
             DDZ(L1)=DDZ(L1)*DDF
           END IF
C
          END DO
         END DO
C
        ELSE IF (NOEAVE(NOECND(N)).EQ.NOECEN) THEN
C do center averaging
C +++++++++++++++++++
C R= (r1   -    r2      )
C       center    center
C
         LI=0
         XI=ZERO
         YI=ZERO
         ZI=ZERO
C ...loop over all pairs of atoms belonging to restraint N
         DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
          I=NOEILS(II)
          XI=XI+X(I)
          YI=YI+Y(I)
          ZI=ZI+Z(I)
          LI=LI+1
         END DO
         XI=XI/LI
         YI=YI/LI
         ZI=ZI/LI
         LJ=0
         XJ=ZERO
         YJ=ZERO
         ZJ=ZERO
         DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
          J=NOEJLS(JJ)
          XJ=XJ+X(J)
          YJ=YJ+Y(J)
          ZJ=ZJ+Z(J)
          LJ=LJ+1
         END DO
         XJ=XJ/LJ
         YJ=YJ/LJ
         ZJ=ZJ/LJ
         XIJ=XI-XJ
         YIJ=YI-YJ
         ZIJ=ZI-ZJ
         RN2=MAX(SQRT(XIJ**2+YIJ**2+ZIJ**2),RSMALL)
C chain rule
         IF (NNORR.GT.1) THEN
          DDF=RN2**(-OREXP(NOECND(N))-1)
         ELSE
          DDF=ONE
         END IF
         DDX(INORR)=XIJ/(RN2*LI)*DDF
         DDY(INORR)=YIJ/(RN2*LI)*DDF
         DDZ(INORR)=ZIJ/(RN2*LI)*DDF
         DDX(NNORR+INORR)=XIJ/(RN2*LJ)*DDF
         DDY(NNORR+INORR)=YIJ/(RN2*LJ)*DDF
         DDZ(NNORR+INORR)=ZIJ/(RN2*LJ)*DDF
C
        ELSE
         WRITE(6,'(A)') ' %ENOE-ERR: averaging - check code'
        END IF
C
        IF (NNORR.GT.1) THEN
         IF (NOEAVE(NOECND(N)).EQ.NOESUM) THEN
          IF (SAVEXP.EQ.OREXP(NOECND(N))) THEN
           RN=RN+RN2
          ELSE
           RN=RN+RN2**(-OREXP(NOECND(N)))
          END IF
         ELSE
          RN=RN+RN2**(-OREXP(NOECND(N)))
         END IF
        ELSE
         RN=RN2
         IF (NOEAVE(NOECND(N)).EQ.NOECEN) THEN
          DF=ONE
         ELSE
          DF=RN**(SAVEXP+1)
         END IF
        END IF

C this is the end of the "OR" loop
       END DO
C
       IF ((NNORR.GT.1)) THEN
        RN2=MAX((ONE/RN)**(ONE/OREXP(NOECND(N))),RSMALL)
        DF=RN2/RN
       END IF
C
       RN=RN2
      END IF
C
C running-average distance
C ++++++++++++++++++++++
      IF (RMODE(NOECND(N)).NE.'OFF'.AND.ANALYS.NE.'ANAL') THEN
       NOESTP(N)=NOESTP(N)+1
       NOERRV(N)=(((NOESTP(N)-1)*(NOERRV(N)**(-RAVEXP(NOECND(N))))+
     &          RN**(-RAVEXP(NOECND(N)))) /
     &          NOESTP(N))**(-ONE/RAVEXP(NOECND(N)))
      END IF
C
C time-averaged distance
C ++++++++++++++++++++++
C Modified by SIOD, 17/11/92
      RNN=RN
      IF (TMODE(NOECND(N)).NE.'OFF'.AND.ANALYS.NE.'ANAL') THEN
        RN=(RN**(-DAVEXP(NOECND(N)))/MEMTAU(NOECND(N))
     &   +(EXP(-ONE/MEMTAU(NOECND(N))))*NOERAV(N)**(-DAVEXP(NOECND(N)))
     &   )**(-ONE/DAVEXP(NOECND(N)))
       IF (TFMODE(NOECND(N)).EQ.'CONS') THEN
        DF=DF*(RN/RNN)**FOUR/MEMTAU(NOECND(N))
       END IF
C
       IF (TMODE(NOECND(N)).NE.'STAT') THEN
        NOERAV(N)=RN
       END IF
      END IF
C
C time-averaging analyse mode only
C ++++++++++++++++++++++++++++++++
C A. Bonvin 17/05/95
      IF (ANALYS.EQ.'ANAL') THEN
       IF (RMODE(NOECND(N)).NE.'OFF'.AND.RANA.EQ.'RAVE') THEN
        RN=NOERRV(N)
        RNN=RN
       END IF
      IF (TMODE(NOECND(N)).NE.'OFF'.AND.RANA.EQ.'TAVE') THEN
       RN=NOERAV(N)
       RNN=RN
       END IF
      END IF
C
C
      IF ((NOECV(N).NE.NOEICV).OR.(ANALYS.EQ.'ANAL')) THEN
C potential function and derivative calculation only for working set
C Mod: if ANALYS mode test set is also included for cross-validation
C A.Bonvin 17/05/95
C
       IF (NOEPOT(NOECND(N)).EQ.NOEBIH) THEN
C do biharmonic potential
C +++++++++++++++++++++++
C
C                      -1         2            c(low)   R < d
C    E= S K T (2 c **2 ) ( R - d )     c   =
C          b      ij                    ij     c(high)  R > d
C
        IF (RN.GT.NOEDIS(N)) THEN
         C=NOEHIG(N)
        ELSE
         C=NOELOW(N)
        END IF
C
C ...convert the constant
        C=ABS(NOEWGH(N)*(NOESCA(NOECND(N))
     &      *KBOLTZ*NOETEM*0.5D0)/(C*C))
        C=MIN(C,NOECEI)
C
        DEVI=RN-NOEDIS(N)
        EL=C*DEVI*DEVI
        DF=DF*TWO*C*DEVI
C
       ELSE IF (NOEPOT(NOECND(N)).EQ.NOESQW) THEN
C do square well potential
C ++++++++++++++++++++++++
C
C                   exp           R-(d+dplus)  d+dplus <  R
C E= const ( delta )    , delta = 0            d-dminus < R <  d+dplus
C                                 (d-dminus)-R   R < d-dminus
C
        DEVI=MAX(ZERO,
     &     RN-(NOEDIS(N)+MAX(ZERO,NOEHIG(N)-NOEOFF(NOECND(N)))))+
     &     MIN(ZERO,RN
     &     -(NOEDIS(N)-MAX(ZERO,NOELOW(N)-NOEMOF(NOECND(N)))))
        C=MIN(NOEWGH(N)*NOESCA(NOECND(N))*NOECNS(NOECND(N)),NOECEI)
        EL=C*DEVI**NOEEXP(NOECND(N))
        DF=DF*NOEEXP(NOECND(N))*C*DEVI**(NOEEXP(NOECND(N))-1)
C
      ELSEIF (NOEPOT(NOECND(N)).EQ.NOELOG) THEN
C do lognormal potential
C +++++++++++++++++++++++
C                    2
C E= const log((d/R))
C
C we calculate the deviation and the deviation of logarithms but keep only
C the deviation of logarithms
        ICLASS=NOECND(N)
        C=MIN(NOEWGH(N)*NOESCA(ICLASS)*NOECNS(ICLASS),NOECEI)
        DEVI=RN-NOEDIS(N)
        LOGDEVI=LOG(RN/NOEDIS(N))
        DF=DF*C*TWO*LOGDEVI/RN
        EL=C*LOGDEVI*LOGDEVI
        DEVI=LOGDEVI
C
       ELSE IF (NOEPOT(NOECND(N)).EQ.NOESOF) THEN
C do square well potential with soft asymptote
C ++++++++++++++++++++++++++++++++++++++++++++
C Modifications: Michael Nilges
C                  softexp
C E= a + b / delta         + c delta,
C                     delta = R-(d+dplus)+rswitch, R>d+dplus+ rswitch
C                   exp          ( R-(d+dplus) , d+dplus < R < d+dplus+rswitch
C E= const ( delta )    , delta= ( 0           , d-dminus < R < d+dplus
C                                ( (d-dminus)-R, R < d-dminus
C
C
        ICLASS=NOECND(N)
        C=MIN(NOEWGH(N)*NOESCA(ICLASS)*NOECNS(ICLASS),NOECEI)
        DEVI=MAX(ZERO,RN
     &       -(NOEDIS(N)+MAX(ZERO,NOEHIG(N)-NOEOFF(ICLASS))))+
     &       MIN(ZERO,RN
     &       -(NOEDIS(N)-MAX(ZERO,NOELOW(N)-NOEMOF(ICLASS))))
C
        IF (DEVI.LT.-NOEMSW(ICLASS)) THEN
         NEXP=FLOAT(NOEEXP(ICLASS))
         NSOE=FLOAT(NOEMSO(ICLASS))
         AAA=(NEXP/NSOE+ONE)*NOEMSW(ICLASS)**NEXP
     &      +NOEMAS(ICLASS)*(ONE+ONE/NSOE)*NOEMSW(ICLASS)
C MN sign corrected for nsoe <> 1 (12-02-99)
         BBB=+NEXP/NSOE*NOEMSW(ICLASS)**(NEXP+NSOE)
     &      +NOEMAS(ICLASS)/NSOE*NOEMSW(ICLASS)**(NSOE+ONE)
         EL=C*(AAA-BBB/ABS(DEVI)**NSOE-NOEMAS(ICLASS)*ABS(DEVI))
         DF=DF*C*(NSOE*BBB/ABS(DEVI)**(NSOE+1)-NOEMAS(ICLASS))
     &     *SIGN(1.0D0, DEVI)
C
        ELSE IF (DEVI.LE.NOESWI(ICLASS)) THEN
         EL=C*DEVI**NOEEXP(ICLASS)
         DF=DF*NOEEXP(ICLASS)*C*DEVI**(NOEEXP(ICLASS)-1)
C
        ELSE
         NEXP=NOEEXP(ICLASS)
         NSOE=NOESOE(ICLASS)
         AAA=(NEXP/NSOE+ONE)*NOESWI(ICLASS)**NEXP
     &      -NOEASY(ICLASS)*(ONE+ONE/NSOE)*NOESWI(ICLASS)
         BBB=-NEXP/NSOE*NOESWI(ICLASS)**(NEXP+NSOE)
     &      +NOEASY(ICLASS)/NSOE*NOESWI(ICLASS)**(NSOE+1)
         EL=C*(AAA+BBB/DEVI**NSOE+NOEASY(ICLASS)*DEVI)
         DF=DF*C*(-NSOE*BBB/DEVI**(NSOE+1)+NOEASY(ICLASS))
        END IF
C
       ELSE
        WRITE(6,'(A)') ' %ENOE-ERR: potential - check code'
       END IF

C
C accumulate energy
       EN=EN+EL
C
C accumulate derivatives
C
C mn_new we add an outer loop that allows calculation of multiple
C restraints connected by "or" statements
C here: calculation of derivatives
       IF (.NOT.QAVER) THEN
        I=NOEILS(NOEIPR(NOEORR(N)+1)+1)
        J=NOEJLS(NOEJPR(NOEORR(N)+1)+1)
        DF=DF/RN
        DX(I)= DX(I)+DF*(X(I)-X(J))
        DX(J)= DX(J)-DF*(X(I)-X(J))
        DY(I)= DY(I)+DF*(Y(I)-Y(J))
        DY(J)= DY(J)-DF*(Y(I)-Y(J))
        DZ(I)= DZ(I)+DF*(Z(I)-Z(J))
        DZ(J)= DZ(J)-DF*(Z(I)-Z(J))
C
       ELSE
C
        L1=0
        INORR=0
        NNORR=NOEORR(N+1)-NOEORR(N)
        DO NORR=NOEORR(N)+1, NOEORR(N+1)
         INORR=INORR+1
C
         IF (NOEAVE(NOECND(N)).EQ.NOEM6
     &    .OR.(NOEAVE(NOECND(N)).EQ.NOEM3)
     &    .OR.(NOEAVE(NOECND(N)).EQ.NOESUM)) THEN
C do r-6/r-3 averaging
C ++++++++++++++++
C
          DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
           DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
            I=NOEILS(II)
            J=NOEJLS(JJ)
C
C ensemble averaging:
C setup weight for each selections ij
C with value 0 if i and j does not belong to have identical
C segids otherwise 1
            IF (EMODE.EQ.'ON') THEN
               CALL ENSBLE(WENS,I,J)
            ELSE
             WENS = ONE
            END IF
C
            IF (WENS.GE.ZERO) THEN
             L1=L1+1
             DX(I)=DX(I)+DF*DDX(L1)
             DX(J)=DX(J)-DF*DDX(L1)
             DY(I)=DY(I)+DF*DDY(L1)
             DY(J)=DY(J)-DF*DDY(L1)
             DZ(I)=DZ(I)+DF*DDZ(L1)
             DZ(J)=DZ(J)-DF*DDZ(L1)
            END IF
           END DO
          END DO
C
         ELSE IF (NOEAVE(NOECND(N)).EQ.NOECEN) THEN
C do center averaging
C +++++++++++++++++++
C
          DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
           I=NOEILS(II)
           DX(I)=DX(I)+DF*DDX(INORR)
           DY(I)=DY(I)+DF*DDY(INORR)
           DZ(I)=DZ(I)+DF*DDZ(INORR)
          END DO
          DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
           J=NOEJLS(JJ)
           DX(J)=DX(J)-DF*DDX(NNORR+INORR)
           DY(J)=DY(J)-DF*DDY(NNORR+INORR)
           DZ(J)=DZ(J)-DF*DDZ(NNORR+INORR)
          END DO
C
         ELSE
          WRITE(6,'(A)') ' %ENOE-ERR: averaging - check code'
         END IF
C
C end of OR loop for force calculation
        END DO
       END IF
C
      END IF
      END IF
      END DO
C
      IF (ANALYS.EQ.'ANAL') THEN
      PCDATA(PCGEOM)=RNN
      PCDATA(PCENER)=EN
      PCDATA(PCEQUI)=NOEDIS(NA)
      PCDATA(PCCONS)=C
      PCDATA(PCDEVI)=-DEVI
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE MULTID(N,HN,EN,RN,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOEDIS,NOELOW,NOEHIG,
     &           NOEWGH,VRN)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      DOUBLE PRECISION EN
      INTEGER NOEORR(*), NOEIPR(*), NOEILS(*), NOEJPR(*)
      INTEGER NOEJLS(*), NOECND(*)
      DOUBLE PRECISION NOEDIS(*), NOELOW(*), NOEHIG(*),NOEWGH(*)
      DOUBLE PRECISION VRN(20)
C local
      DOUBLE PRECISION RN
      DOUBLE PRECISION EL
      DOUBLE PRECISION DEVI
      DOUBLE PRECISION ELOW, EHIG, HE
      DOUBLE PRECISION ZZ(20),ZZZ(20),ZZZZ(20),DZZ(20),DZZZ(20)
      DOUBLE PRECISION D1ZZZ(20), D2ZZZ(20)
      DOUBLE PRECISION LIM1,LIM2,E1,E2,FAKT1,FAKT2,HV1,HV2,ELIM1
      DOUBLE PRECISION HA,HB,HC,HDZ,HDZZ,HDZZZ,HZ,HZZ,HZZZ,HZD,HZZD
      DOUBLE PRECISION VXI(20),VYI(20),VZI(20),VLJ(20),VXJ(20),VYJ(20)
      DOUBLE PRECISION VZJ(20),VXIJ(20),VYIJ(20),VZIJ(20)
      DOUBLE PRECISION VDF(20)
      DOUBLE PRECISION VSIJ(20),VRIJ(20),VDDX(20,20)
      DOUBLE PRECISION VDDY(20,20),VDDZ(20,20)
      INTEGER VLI(20)
      INTEGER I, J, N, L, II, JJ
      INTEGER HN, K, BV
C parameter
      DOUBLE PRECISION SIX, SEVEN, ZERO, ONE, TWO, THREE, TWELVE
      DOUBLE PRECISION FOUR, TEN, P15
      PARAMETER (SIX=6.0D0, SEVEN=7.0D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (THREE=3.0D0, ZERO=0.0D0, TWELVE=12.0D0, FOUR=4.0D0)
      PARAMETER (TEN=10.D0, P15=1.5D0)
C begin
      IF ((NOEORR(N+1)-NOEORR(N)).GT.1) THEN
      WRITE(6,'(A)') '%ENOE-ERR: OR not allowed in MULTID'
C
      ELSEIF (NOEPOT(NOECND(N)).EQ.NOEHDI) THEN
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C do the highdimensional potential for stereospecific assignment
C
      HN = HN + 1
C HN counts number of read in noes in higdim pot.; if HN=1 new
C potential is started to
C be constructed with the  next noecor noes. if done then HN set to 0.
      E1=ZERO
      E2=ZERO
      IF (HN.EQ.1) THEN
C only at the beginning of a new class of noes
C
      DO L=1, NOECOR(NOECND(N)),2
C loop over the whole noe class, pairwise
C
      DO K=L,L+1
C loop for each noe pair
C
      IF (NOEAVE(NOECND(N)).NE.NOECEN) THEN
      WRITE(6,'(A)') '%ENOE-ERR: averaging - check code'
      END IF
C do center averaging
C +++++++++++++++++++
C R= (r1   -    r2      )
C       center    center
C
      VLI(K)=0
      VXI(K)=ZERO
      VYI(K)=ZERO
      VZI(K)=ZERO
C ...loop over all pairs of atoms belonging to restraint N
CSGI-specific optimization
C*$*SCALAR OPTIMIZE(0)
      DO II=NOEIPR(N+K-1)+1,NOEIPR(N+K)
      I=NOEILS(II)
      VXI(K)=VXI(K)+X(I)
      VYI(K)=VYI(K)+Y(I)
      VZI(K)=VZI(K)+Z(I)
      VLI(K)=VLI(K)+1
      END DO
CSGI-specific optimization
C*$*SCALAR OPTIMIZE(3)
      VXI(K)=VXI(K)/VLI(K)
      VYI(K)=VYI(K)/VLI(K)
      VZI(K)=VZI(K)/VLI(K)
      VLJ(K)=0
      VXJ(K)=ZERO
      VYJ(K)=ZERO
      VZJ(K)=ZERO
      DO JJ=NOEJPR(N+K-1)+1,NOEJPR(N+K)
      J=NOEJLS(JJ)
      VXJ(K)=VXJ(K)+X(J)
      VYJ(K)=VYJ(K)+Y(J)
      VZJ(K)=VZJ(K)+Z(J)
      VLJ(K)=VLJ(K)+1
      END DO
      VXJ(K)=VXJ(K)/VLJ(K)
      VYJ(K)=VYJ(K)/VLJ(K)
      VZJ(K)=VZJ(K)/VLJ(K)
      VXIJ(K)=VXI(K)-VXJ(K)
      VYIJ(K)=VYI(K)-VYJ(K)
      VZIJ(K)=VZI(K)-VZJ(K)
      VRN(K)=SQRT(VXIJ(K)**2+VYIJ(K)**2+VZIJ(K)**2)
      END DO
C***************************************
      LIM1=0.1D0
      LIM2=0.2D0
      FAKT1=ONE
      FAKT2=ONE
C LIM1 is distance from where on the square function is set linear
C LIM2 distance between the two potential minima from where on hight
C is set.
C FAKT1 variable for barrier calculation
C FAKT2 variable for setting the gradient
C****************************************
      ZZ(L)=ZERO
      ZZZ(L)=ZERO
      ZZZZ(L)=ZERO
      ZZ(L+1)=ZERO
      ZZZ(L+1)=ZERO
      ZZZZ(L+1)=ZERO
      DZZ(L)=ZERO
      DZZ(L+1)=ZERO
      DZZZ(L)=ZERO
      DZZZ(L+1)=ZERO
      D1ZZZ(L)=ZERO
      D1ZZZ(L+1)=ZERO
      D2ZZZ(L)=ZERO
      D2ZZZ(L+1)=ZERO
      DEVI =ZERO
C calculating FAKT1
      IF (NOELOW(N+L-1).LT.NOELOW(N+L)) THEN
      ELIM1= NOELOW(N+L) - NOELOW(N+L-1) + LIM1
      ELSE
      ELIM1= NOELOW(N+L-1)-NOELOW(N+L) + LIM1
      END IF
C
C ATB
      BV=0
C
      IF ((NOEHIG(N+L-1).GT.ZERO).AND.(NOEHIG(N+L).GT.ZERO)) THEN
      IF(NOEHIG(N+L-1)*SQRT(TWO).LT.LIM2) THEN
      BV=0
      ELSE
      BV=1
      HV1=NOEDIS(L+N)-NOELOW(L+N-1)
      IF (HV1.GT.ZERO) THEN
      HV1=HV1/TWO
      ELSE
      HV1=ZERO
      END IF
      HV2=SQRT(TWO)*NOEHIG(N+L-1)/TWO
      END IF
      ELSE IF ((NOEHIG(N+L-1).LT.ZERO).AND.(NOEHIG(N+L).LT.ZERO)) THEN
      IF (ABS(NOEHIG(N+L))*SQRT(TWO).LT.LIM2) THEN
      BV=0
      ELSE
      BV=1
      HV1=NOEDIS(L+N-1)-NOELOW(L+N)
      IF (HV1.GT.ZERO) THEN
      HV1=HV1/TWO
      ELSE
      HV1=ZERO
      END IF
      HV2=SQRT(TWO)*ABS(NOEHIG(N+L))/TWO
      END IF
      END IF
C
      IF (BV.EQ.1) THEN
      IF (HV1.GT.LIM1) THEN
      HZ=LIM1**2+2*LIM1*(HV1-LIM1)
      ELSE
      HZ=HV1**2
      END IF
      IF (HV2.GT.LIM1) THEN
      HZZ=LIM1**2 + 2*LIM1*(HV2-LIM1)
      ELSE
      HZZ=HV2**2
      END IF
      FAKT1=2*SQRT(NOEBAR(NOECND(N)))/NOECOR(NOECND(N))/(2*HZ+HZZ)
      ELSE
      FAKT1=ONE
      END IF
C calculate FAKT2
      IF (NOELOW(L+N).GT.NOELOW(L+N-1)) THEN
      HDZ=NOELOW(N+L) + LIM1 - NOELOW(N+L-1)
      HDZZ = LIM1
      IF (NOEHIG(N+L-1).GT.ZERO) THEN
      HDZZZ=NOEHIG(N+L-1)/SQRT(TWO)
      ELSE
      HDZZZ=ZERO
      END IF
C
      ELSE IF (NOELOW(L+N-1).GE.NOELOW(N+L)) THEN
      HDZ=NOELOW(N+L-1)+LIM1-NOELOW(N+L)
      HDZZ=LIM1
      IF (NOEHIG(N+L).LT.ZERO) THEN
      HDZZZ=ABS(NOEHIG(N+L))/SQRT(TWO)
      ELSE
      HDZZZ=ZERO
      END IF
      END IF
C
      HZZ=LIM1**2 *FAKT1
      HZZD=FAKT1*TWO*LIM1
      IF (HDZ.GT.LIM1) THEN
      HZ=(LIM1**2 + TWO*LIM1*(HDZ-LIM1))*FAKT1
      HZD=FAKT1*TWO*LIM1
      ELSE
      HZ=HDZ**2 *FAKT1
      HZD=TWO*FAKT1*HDZ
      END IF
      IF (HDZZZ.GT.LIM1) THEN
      HZZZ=(LIM1**2 + TWO*LIM1*(HDZZZ-LIM1))*FAKT1
      ELSE
      HZZZ=FAKT1*HDZZZ**2
      END IF
      HA=HZZ*HZZD
      HB= HZZ*HZD + HZ*HZZD + HZZZ*HZZD
      HC= HZ*HZD + HZD*HZZZ - ELIM1
      FAKT2=(-HB + SQRT(HB**2 - FOUR*HA*HC))/(TWO*HA)
C
      IF (BV.EQ.0) THEN
      FAKT2=ONE
      FAKT1=(2 * LIM1 *
     &          (LIM1**2 +           2*LIM1*(HDZ-LIM1) + HDZZZ**2 +
     &           LIM1**2 + LIM1**2 + 2*LIM1*(HDZ-LIM1) + HDZZZ**2 +
     &           LIM1**2))
      FAKT1=SQRT(ELIM1/FAKT1)
      END IF
C
      IF ((BV.NE.0).AND.(FAKT2.LE.ZERO)) THEN
      FAKT2=0.0001D0
      END IF
C
C calculate noe(i1)x
      IF (VRN(L).GT.NOELOW(N+L-1)) THEN
      DEVI=VRN(L)-NOELOW(N+L-1)
      IF (DEVI.LT.LIM1) THEN
      ZZ(L)=(DEVI)**2
      DZZ(L)=TWO*DEVI
      ELSE
      ZZ(L)=LIM1**2  + (DEVI-LIM1)*2*LIM1
      DZZ(L)=TWO*LIM1
      END IF
      ELSE
      IF (VRN(L).LT.NOEDIS(N+L-1)) THEN
      DEVI=NOEDIS(N+L-1)-VRN(L)
      IF (DEVI.LT.LIM1) THEN
      ZZ(L)=(DEVI)**2
      DZZ(L)=-TWO*DEVI
      ELSE
      ZZ(L)= LIM1**2 + (DEVI-LIM1)*2*LIM1
      DZZ(L)=-TWO*LIM1
      END IF
      END IF
      END IF
C calculate noe(i1)y
      IF (VRN(L+1).GT.NOELOW(N+L)) THEN
      DEVI=VRN(L+1)-NOELOW(N+L)
      IF (DEVI.LT.LIM1) THEN
      ZZZ(L)=(DEVI)**2
      DZZZ(L)=TWO*DEVI
      ELSE
      ZZZ(L)=LIM1**2  + (DEVI-LIM1)*2*LIM1
      DZZZ(L)=TWO*LIM1
      END IF
      ELSE
      IF (VRN(L+1).LT.NOEDIS(N+L)) THEN
      DEVI=NOEDIS(N+L)-VRN(L+1)
      IF (DEVI.LT.LIM1) THEN
      ZZZ(L)=(DEVI)**2
      DZZZ(L)=-TWO*DEVI
      ELSE
      ZZZ(L)= LIM1**2 + (DEVI-LIM1)*2*LIM1
      DZZZ(L)=-TWO*LIM1
      END IF
      END IF
      END IF
C delta noe(i1)
      IF (VRN(L+1).LT.VRN(L)+NOEHIG(N+L-1)) THEN
      DEVI=(-VRN(L+1)+VRN(L)+NOEHIG(N+L-1))/SQRT(TWO)
      IF (DEVI.LT.LIM1) THEN
      ZZZZ(L)=DEVI**2
      D1ZZZ(L)=TWO*DEVI/SQRT(TWO)
      D2ZZZ(L)=-TWO*DEVI/SQRT(TWO)
      ELSE
      ZZZZ(L)=LIM1**2 + (DEVI-LIM1)*2*LIM1
      D1ZZZ(L)=TWO*LIM1/SQRT(TWO)
      D2ZZZ(L)=-TWO*LIM1/SQRT(TWO)
      END IF
      ELSE IF (VRN(L+1).GT.VRN(L)+NOEHIG(N+L)) THEN
      DEVI=(VRN(L+1)-VRN(L)-NOEHIG(N+L))/SQRT(TWO)
      IF (DEVI.LT.LIM1) THEN
      ZZZZ(L)=DEVI**2
      D1ZZZ(L)=-TWO*DEVI/SQRT(TWO)
      D2ZZZ(L)=TWO*DEVI/SQRT(TWO)
      ELSE
      ZZZZ(L)=LIM1**2+(DEVI-LIM1)*2*LIM1
      D1ZZZ(L)=-TWO*LIM1/SQRT(TWO)
      D2ZZZ(L)=TWO*LIM1/SQRT(TWO)
      END IF
      END IF
C
      IF (VRN(L).GT.NOELOW(N+L)) THEN
      DEVI=VRN(L)-NOELOW(N+L)
      IF (DEVI.LT.LIM1) THEN
      ZZ(L+1)=DEVI**2
      DZZ(L+1)=TWO*DEVI
      ELSE
      ZZ(L+1)=LIM1**2  + (DEVI-LIM1)*2*LIM1
      DZZ(L+1)=TWO*LIM1
      END IF
      ELSE
      IF (VRN(L).LT.NOEDIS(N+L)) THEN
      DEVI=NOEDIS(N+L)-VRN(L)
      IF (DEVI.LT.LIM1) THEN
      ZZ(L+1)=DEVI**2
      DZZ(L+1)=-TWO*DEVI
      ELSE
      ZZ(L+1)= LIM1**2 + (DEVI-LIM1)*2*LIM1
      DZZ(L+1)=-TWO*LIM1
      END IF
      END IF
      END IF
C
      IF (VRN(L+1).GT.NOELOW(N+L-1)) THEN
      DEVI=VRN(L+1)-NOELOW(N+L-1)
      IF (DEVI.LT.LIM1) THEN
      ZZZ(L+1)=DEVI**2
      DZZZ(L+1)=TWO*DEVI
      ELSE
      ZZZ(L+1)=LIM1**2  + (DEVI-LIM1)*2*LIM1
      DZZZ(L+1)=TWO*LIM1
      END IF
      ELSE
      IF (VRN(L+1).LT.NOEDIS(N+L-1)) THEN
      DEVI=NOEDIS(N+L-1)-VRN(L+1)
      IF (DEVI.LT.LIM1) THEN
      ZZZ(L+1)=DEVI**2
      DZZZ(L+1)=-TWO*DEVI
      ELSE
      ZZZ(L+1)= LIM1**2 + (DEVI-LIM1)*2*LIM1
      DZZZ(L+1)=-TWO*LIM1
      END IF
      END IF
      END IF
C
      IF (VRN(L+1).LT.VRN(L)-NOEHIG(N+L)) THEN
      DEVI=(-VRN(L+1)+VRN(L)-NOEHIG(N+L))/SQRT(TWO)
      IF (DEVI.LT.LIM1) THEN
      ZZZZ(L+1)=DEVI**2
      D1ZZZ(L+1)=2*DEVI/SQRT(TWO)
      D2ZZZ(L+1)=-2*DEVI/SQRT(TWO)
      ELSE
      ZZZZ(L+1)=LIM1**2 + (DEVI-LIM1)*2*LIM1
      D1ZZZ(L+1)=2*LIM1/SQRT(TWO)
      D2ZZZ(L+1)=-2*LIM1/SQRT(TWO)
      END IF
      ELSE IF (VRN(L+1).GT.VRN(L)-NOEHIG(N+L-1)) THEN
      DEVI=(VRN(L+1)-VRN(L)+NOEHIG(N+L-1))/SQRT(TWO)
      IF (DEVI.LT.LIM1) THEN
      ZZZZ(L+1)=DEVI**2
      D1ZZZ(L+1)=-2*DEVI/SQRT(TWO)
      D2ZZZ(L+1)=2*DEVI/SQRT(TWO)
      ELSE
      ZZZZ(L+1)=LIM1**2+(DEVI-LIM1)*2*LIM1
      D1ZZZ(L+1)=-2*LIM1/SQRT(TWO)
      D2ZZZ(L+1)=2*LIM1/SQRT(TWO)
      END IF
      END IF
C
      DZZ(L)=DZZ(L)*FAKT1
      DZZ(L+1)=DZZ(L+1)*FAKT1
      DZZZ(L)=DZZZ(L)*FAKT1
      DZZZ(L+1)=DZZZ(L+1)*FAKT1
      D1ZZZ(L)=D1ZZZ(L)*FAKT1
      D2ZZZ(L)=D2ZZZ(L)*FAKT1
      D1ZZZ(L+1)=D1ZZZ(L+1)*FAKT1
      D2ZZZ(L+1)=D2ZZZ(L+1)*FAKT1
C multiply FAKT2
      IF ((VRN(L).LT.NOEDIS(L+N-1)).AND.(VRN(L).LT.NOEDIS(L+N))) THEN
      IF (NOEDIS(L+N).LE.NOEDIS(L+N-1)) THEN
      ZZ(L+1)=FAKT2*ZZ(L+1)
      DZZ(L+1)=FAKT2*DZZ(L+1)
      ELSE IF (NOEDIS(L+N-1).LE.NOEDIS(L+N)) THEN
      ZZ(L) = FAKT2*ZZ(L)
      DZZ(L)=FAKT2*DZZ(L)
      END IF
      ELSE IF ((VRN(L).GT.NOELOW(L+N))
     &  .AND.(VRN(L).GT.NOELOW(L+N-1))) THEN
      IF (NOELOW(L+N).GE.NOELOW(L+N-1)) THEN
      ZZ(L+1)=FAKT2*ZZ(L+1)
      DZZ(L+1)=FAKT2*DZZ(L+1)
      ELSE IF (NOELOW(L+N-1).GE.NOELOW(L+N)) THEN
      ZZ(L)=FAKT2*ZZ(L)
      DZZ(L)=FAKT2*DZZ(L)
      END IF
      END IF
C
      IF ((VRN(L+1).LT.NOEDIS(L+N-1))
     &        .AND.(VRN(L+1).LT.NOEDIS(L+N))) THEN
      IF (NOEDIS(L+N).LE.NOEDIS(L+N-1)) THEN
      ZZZ(L)=FAKT2*ZZZ(L)
      DZZZ(L)=DZZZ(L)*FAKT2
      ELSE IF (NOEDIS(L+N-1).LE.NOEDIS(L+N)) THEN
      ZZZ(L+1) = FAKT2*ZZZ(L+1)
      DZZZ(L+1)= DZZZ(L+1)*FAKT2
      END IF
      ELSE IF ((VRN(L+1).GT.NOELOW(L+N))
     &        .AND.(VRN(L+1).GT.NOELOW(L+N-1))) THEN
      IF (NOELOW(L+N).GE.NOELOW(L+N-1)) THEN
      ZZZ(L)=FAKT2*ZZZ(L)
      DZZZ(L)=FAKT2*DZZZ(L)
      ELSE IF (NOELOW(L+N-1).GE.NOELOW(L+N)) THEN
      ZZZ(L+1)=FAKT2*ZZZ(L+1)
      DZZZ(L+1)=FAKT2*DZZZ(L+1)
      END IF
      END IF
C
      IF (((VRN(L+1).GT.ABS(NOEHIG(L+N))+VRN(L)).AND.
     &     (VRN(L+1).GT.VRN(L)+ABS(NOEHIG(L+N-1)))).OR.
     &     ((VRN(L).GT.ABS(NOEHIG(L+N))+VRN(L+1)).AND.
     &     (VRN(L).GT.ABS(NOEHIG(L+N-1))+VRN(L+1)))) THEN
      IF (ZZZZ(L+1).LE.ZZZZ(L)) THEN
      ZZZZ(L+1)=ZZZZ(L+1)*FAKT2
      D1ZZZ(L+1)=D1ZZZ(L+1)*FAKT2
      D2ZZZ(L+1)=D2ZZZ(L+1)*FAKT2
      ELSE
      ZZZZ(L)=ZZZZ(L)*FAKT2
      D1ZZZ(L)=D1ZZZ(L)*FAKT2
      D2ZZZ(L)=D2ZZZ(L)*FAKT2
      END IF
      END IF
      ZZZZ(L)=ZZZZ(L)*P15
      ZZZZ(L+1)=ZZZZ(L+1)*P15
      D1ZZZ(L)=D1ZZZ(L)*P15
      D2ZZZ(L)=D2ZZZ(L)*P15
      D1ZZZ(L+1)=D1ZZZ(L+1)*P15
      D2ZZZ(L+1)=D2ZZZ(L+1)*P15
      E1=E1 + FAKT1*(ZZ(L)+ZZZ(L) + ZZZZ(L))
      E2=E2 + FAKT1*(ZZ(L+1)+ZZZ(L+1) + ZZZZ(L+1))
      END DO
      DO L=1,NOECOR(NOECND(N)),2
      VDF(L)=NOEWGH(N)*NOESCA(NOECND(N))*(E2*(DZZ(L)+D1ZZZ(L))+
     &                              E1*(DZZ(L+1)+D1ZZZ(L+1)))
      VDF(L+1)=NOEWGH(N)*NOESCA(NOECND(N))*(E2*(DZZZ(L) + D2ZZZ(L))
     &                            + E1*(DZZZ(L+1) + D2ZZZ(L+1)))
      END DO
C*************************************
      EL = E1*E2* NOEWGH(N)*NOESCA(NOECND(N))
C accumulate energy
      EN=EN+EL
      DO K=1,NOECOR(NOECND(N))
C
C accumulate derivatives
C
C do center averaging
C +++++++++++++++++++
      L=0
      DO II=NOEIPR(N+K-1)+1,NOEIPR(N+K)
      I=NOEILS(II)
      DX(I)=DX(I)+VDF(K)*VXIJ(K)/(VRN(K)*VLI(K))
      DY(I)=DY(I)+VDF(K)*VYIJ(K)/(VRN(K)*VLI(K))
      DZ(I)=DZ(I)+VDF(K)*VZIJ(K)/(VRN(K)*VLI(K))
      END DO
      DO JJ=NOEJPR(N+K-1)+1,NOEJPR(N+K)
      J=NOEJLS(JJ)
      DX(J)=DX(J)-VDF(K)*VXIJ(K)/(VRN(K)*VLJ(K))
      DY(J)=DY(J)-VDF(K)*VYIJ(K)/(VRN(K)*VLJ(K))
      DZ(J)=DZ(J)-VDF(K)*VZIJ(K)/(VRN(K)*VLJ(K))
      END DO
      END DO
C
      ELSE IF (HN.EQ.NOECOR(NOECND(N))) THEN
      HN=0
      END IF
C +++++++++++++++++++++++++++++++++++++++++++++++++++++
C end highdimensional potential
C
      ELSE IF (NOEPOT(NOECND(N)).EQ.NOE3DP) THEN
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  do 3dpotential
C
      HN = HN + 1
      E1=ZERO
      E2=ZERO
      IF (HN.EQ.1) THEN
C only at the beginning of a new class of noes
C
      DO K=1,2
C loop for each noe pair
C
      IF (NOEAVE(NOECND(N)).EQ.NOEM6) THEN
C do r-6 averaging
C ++++++++++++++++
C      WRITE(PUNIT,'(A)') 'DO R-6 AVERAGING '
C R = ( <R-6> ) -1/6
      VLI(K)=0
      VRN(K)=ZERO
C loop over all pairs of atoms belonging to restraint N
      DO II=NOEIPR(N+K-1)+1,NOEIPR(N+K)
      DO JJ=NOEJPR(N+K-1)+1,NOEJPR(N+K)
      I=NOEILS(II)
      J=NOEJLS(JJ)
      VXIJ(K)=X(I)-X(J)
      VYIJ(K)=Y(I)-Y(J)
      VZIJ(K)=Z(I)-Z(J)
      VSIJ(K)=VXIJ(K)**2+VYIJ(K)**2+VZIJ(K)**2
      VRIJ(K)=ONE/(VSIJ(K)**(THREE))
      VRN(K)=VRIJ(K)+VRN(K)
C local derivative accumulation
      VLI(K)=VLI(K)+1
      VRIJ(K)=VRIJ(K)/VSIJ(K)
      VDDX(VLI(K),K)=VRIJ(K)*VXIJ(K)
      VDDY(VLI(K),K)=VRIJ(K)*VYIJ(K)
      VDDZ(VLI(K),K)=VRIJ(K)*VZIJ(K)
      END DO
      END DO
C
      VRN(K)=(VLI(K)/VRN(K))**(ONE/SIX)
C
      ELSE IF (NOEAVE(NOECND(N)).EQ.NOECEN) THEN

C do center averaging
C +++++++++++++++++++
C R= (r1   -    r2      )
C       center    center
C
      VLI(K)=0
      VXI(K)=ZERO
      VYI(K)=ZERO
      VZI(K)=ZERO
C ...loop over all pairs of atoms belonging to restraint N
CSGI-specific optimization
C*$*SCALAR OPTIMIZE(0)
      DO II=NOEIPR(N+K-1)+1,NOEIPR(N+K)
      I=NOEILS(II)
      VXI(K)=VXI(K)+X(I)
      VYI(K)=VYI(K)+Y(I)
      VZI(K)=VZI(K)+Z(I)
      VLI(K)=VLI(K)+1
      END DO
CSGI-specific optimization
C*$*SCALAR OPTIMIZE(3)
      VXI(K)=VXI(K)/VLI(K)
      VYI(K)=VYI(K)/VLI(K)
      VZI(K)=VZI(K)/VLI(K)
      VLJ(K)=0
      VXJ(K)=ZERO
      VYJ(K)=ZERO
      VZJ(K)=ZERO
      DO JJ=NOEJPR(N+K-1)+1,NOEJPR(N+K)
      J=NOEJLS(JJ)
      VXJ(K)=VXJ(K)+X(J)
      VYJ(K)=VYJ(K)+Y(J)
      VZJ(K)=VZJ(K)+Z(J)
      VLJ(K)=VLJ(K)+1
      END DO
      VXJ(K)=VXJ(K)/VLJ(K)
      VYJ(K)=VYJ(K)/VLJ(K)
      VZJ(K)=VZJ(K)/VLJ(K)
      VXIJ(K)=VXI(K)-VXJ(K)
      VYIJ(K)=VYI(K)-VYJ(K)
      VZIJ(K)=VZI(K)-VZJ(K)
      VRN(K)=SQRT(VXIJ(K)**2+VYIJ(K)**2+VZIJ(K)**2)
      END IF
      END DO
C
      IF (NOEDIS(N).LE.NOELOW(N)) NOELOW(N)=NOEDIS(N)/TEN
      ELOW=(NOESWI(NOECND(N))*(NOEDIS(N)-NOELOW(N)))**(-ONE/TWELVE)
      EHIG=(NOESWI(NOECND(N))*(NOEDIS(N)+NOEHIG(N)))**(-ONE/TWELVE)
      E2=(VRN(1)**(-SIX)*VRN(2)**(-SIX))**(-ONE/TWELVE)
      IF (E2.GT.ELOW) THEN
      HE=E2-ELOW
      ELSE IF (E2.LT.EHIG) THEN
      HE=E2-EHIG
      ELSE
      HE=0
      END IF
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++
      VDF(1)=NOEWGH(N)*NOESCA(NOECND(N))*(HE)*E2**13*VRN(2)**(-SIX)*
     &             VRN(1)**(-SEVEN)
      VDF(2)=NOEWGH(N)*NOESCA(NOECND(N))*(HE)*E2**13*VRN(2)**(-SEVEN)*
     &             VRN(1)**(-SIX)
      EL = NOEWGH(N)*NOESCA(NOECND(N))*(HE)**2
C
C accumulate energy
      EN=EN+EL
      DO K=1,2
C
C accumulate derivatives
C
      IF (NOEAVE(NOECND(N)).EQ.NOEM6) THEN
C do r-6 averaging
C ++++++++++++++++
      VDF(K)=(VDF(K)*VRN(K)**(SEVEN))/VLI(K)
      VLI(K)=0
      DO II=NOEIPR(N+K-1)+1,NOEIPR(N+K)
      DO JJ=NOEJPR(N+K-1)+1,NOEJPR(N+K)
      I=NOEILS(II)
      J=NOEJLS(JJ)
      VLI(K)=VLI(K)+1
      DX(I)=DX(I)+VDF(K)*VDDX(VLI(K),K)
      DX(J)=DX(J)-VDF(K)*VDDX(VLI(K),K)
      DY(I)=DY(I)+VDF(K)*VDDY(VLI(K),K)
      DY(J)=DY(J)-VDF(K)*VDDY(VLI(K),K)
      DZ(I)=DZ(I)+VDF(K)*VDDZ(VLI(K),K)
      DZ(J)=DZ(J)-VDF(K)*VDDZ(VLI(K),K)
      END DO
      END DO
C
      ELSE IF (NOEAVE(NOECND(N)).EQ.NOECEN) THEN
C do center averaging
C +++++++++++++++++++
      L=0
      DO II=NOEIPR(N+K-1)+1,NOEIPR(N+K)
      I=NOEILS(II)
      DX(I)=DX(I)+VDF(K)*VXIJ(K)/(VRN(K)*VLI(K))
      DY(I)=DY(I)+VDF(K)*VYIJ(K)/(VRN(K)*VLI(K))
      DZ(I)=DZ(I)+VDF(K)*VZIJ(K)/(VRN(K)*VLI(K))
      END DO
      DO JJ=NOEJPR(N+K-1)+1,NOEJPR(N+K)
      J=NOEJLS(JJ)
      DX(J)=DX(J)-VDF(K)*VXIJ(K)/(VRN(K)*VLJ(K))
      DY(J)=DY(J)-VDF(K)*VYIJ(K)/(VRN(K)*VLJ(K))
      DZ(J)=DZ(J)-VDF(K)*VZIJ(K)/(VRN(K)*VLJ(K))
      END DO
      END IF
      END DO
C
      RN=VRN(1)
      ELSE IF (HN.EQ.2.0) THEN
      RN=VRN(2)
      HN=0
      END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C end 3dpot
C
      END IF
C
      RETURN
      END
C====================================================================
      SUBROUTINE DISNOE(ICL1,ICL2,NOEDST,NOECND)
C
C
C subroutine reclassifies noe restraints between two classes:
C If the actual deviation is larger than NOEDST, it goes into class 2,
C otherwise it goes into class 1.
C
C Author: Michael Nilges
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'comand.inc'
      INTEGER ICL1, ICL2
      DOUBLE PRECISION NOEDST
      INTEGER NOECND(*)
C local
      DOUBLE PRECISION EN
      INTEGER N, NEXCL, N1, N2
C begin
      NEXCL=0
      N1=0
      N2=0
      IF (NOENUM.GT.0) THEN
      DO N=1,NOENUM
      CALL ENOE(EN,'ANAL',N)
      IF (NOECND(N).EQ.ICL1) THEN
      IF (ABS(PCDATA(PCDEVI)).GE.NOEDST) THEN
      NOECND(N)=ICL2
      NEXCL=NEXCL+1
      N2=N2+1
      ELSE
      N1=N1+1
      END IF
      ELSE
      IF (NOECND(N).EQ.ICL2) THEN
      IF (ABS(PCDATA(PCDEVI)).LT.NOEDST) THEN
      NOECND(N)=ICL1
      NEXCL=NEXCL+1
      N1=N1+1
      ELSE
      N2=N2+1
      END IF
      END IF
      END IF
      END DO
      WRITE(6,'(A,I5,A)') ' DISNOE:',NEXCL,' restraints reclassified'
      WRITE(6,'(A,I5,A,A4)')
     &     '         ',N1,' restraints in class ',NOECNM(ICL1)
      WRITE(6,'(A,I5,A,A4)')
     &     '         ',N2,' restraints in class ',NOECNM(ICL2)
      END IF
      RETURN
      END
C====================================================================
      SUBROUTINE RESETN(ICL1,NOECND,NOERAV,NOEDIS,SRESET)
C
C
C This subroutine resets current values of either the NOERAV or NOERRV
C arrays; these contain the current (exponentially weighted)
C time averaged-distances, and the running-total distances respectively.
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'comand.inc'
      INTEGER ICL1
      INTEGER NOECND(*)
      DOUBLE PRECISION NOERAV(*), NOEDIS(*)
      CHARACTER*4 SRESET
C local
      INTEGER N
      DOUBLE PRECISION EN
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      IF (NOENUM.GT.0) THEN
      IF (SRESET.EQ.'CURR') THEN
      DO N=1,NOENUM
      IF (NOECND(N).EQ.ICL1) THEN
      CALL ENOE(EN,'ANAL',N)
      NOERAV(N)=PCDATA(PCGEOM)
      END IF
      END DO
      ELSEIF (SRESET.EQ.'CONS') THEN
      DO N=1,NOENUM
      IF (NOECND(N).EQ.ICL1) THEN
      NOERAV(N)=NOEDIS(N)
      END IF
      END DO
      END IF
      END IF
      RETURN
      END
C====================================================================
      SUBROUTINE RESETC(ICL1,NOESTP,NOECND)
C
C
C This subroutine resets current values of the NOESTP array.
C This array contains the current number of time/energy steps used
C to calculate the running-average.
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'comand.inc'
      INTEGER ICL1, NOESTP(*), NOECND(*)
C local
      INTEGER N
C begin
      DO N=1,NOENUM
      IF (NOECND(N).EQ.ICL1) THEN
      NOESTP(N)=0
      END IF
      END DO
      RETURN
      END
C====================================================================
      SUBROUTINE COUNTV(ICL1,NOECND)
C
C
C subroutine counts noe violations in a class
C and stores the smallest positive violation in $result.
C if there are no violations, $result is set to 0, if there
C are no restraints, $result is set to -1.
C
C Author: Michael Nilges
C
      IMPLICIT NONE
C I/O

      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'symbol.inc'
      INTEGER ICL1
      INTEGER NOECND(*)
C local
      DOUBLE PRECISION EN, ZERO, MINVIO
      INTEGER N, NVIOL, MONE, NCLASS
      PARAMETER (ZERO=0.0D0,MONE=-1.0D0)
      DOUBLE COMPLEX DBCOMP
C begin
      NVIOL=0
      NCLASS=0
      MINVIO=9999.0D0
      IF (NOENUM.GT.0) THEN
      DO N=1,NOENUM
      IF (NOECND(N).EQ.ICL1) THEN
      CALL ENOE(EN,'ANAL',N)
      NCLASS=NCLASS+1
      IF (ABS(PCDATA(PCDEVI)).GT.ZERO) THEN
      NVIOL=NVIOL+1
      MINVIO=MIN(MINVIO,ABS(PCDATA(PCDEVI)))
      END IF
      END IF
      END DO
      END IF
      IF(NVIOL.EQ.0)MINVIO=ZERO
      IF(NCLASS.EQ.0)MINVIO=MONE
      WRITE(6,'(A,I5,A,A4)')
     &     ' COUNTV:',NVIOL,' violations in class ',NOECNM(ICL1)
C
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, MINVIO )
      CALL DECLAR( 'MIN_VIOLATIONS', 'DP', ' ', DBCOMP, MINVIO )
      IF (NVIOL.GT.0) THEN
      WRITE(6,'(A,F8.3)') '        smallest violation=',MINVIO
      END IF
C
      RETURN
      END
C====================================================================
      SUBROUTINE NOEINI
C
C subroutine initializes the NOE list
C

C Author: Axel Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
C begin
      NOETEM=300.0D0
      NOECEI=30.0D0
C
C set heap pointers of dynamic data structure
C to zero
      HPNORR=0
      HPNIPR=0
      HPNILS=0
      HPNJPR=0
      HPNJLS=0
      HPNCND=0
      HPNRAV=0
      HPNRRV=0
      HPNNSP=0
      HPNDIS=0
      HPDENINIT=0
      HPNLOW=0
      HPNHIG=0
      HPNWGH=0
      HPNCV=0
      HPNVIO=0
      HPNPID=0
      HPNPP1=0
      HPNPP2=0
      HPNVOL=0
C
C restraint analysis
      HPNHGL=0
      HPNMAT=0
      HPNMA2=0
C
      HPNCAL=0
      HPNPSE=0
      HPNDIA=0
C
C reset other variables
      CALL NOERES
      RETURN
      END
C====================================================================
      SUBROUTINE NOERES
C
C subroutine resets the NOE list
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'consta.inc'
C begin
      CALL NOEHP(0)
C
      NOESMX=0
      NOECCN=1
      NOECNM(1)='NIL'
      NOESCA(1)=1.0D0
      NOEAVE(1)=NOEM6
      NOEPOT(1)=NOEBIH
      NOEEXP(1)=2
      NOECNS(1)=20.0D0
      NOEOFF(1)=0.0D0
      NOEMOF(1)=0.0D0
      NOESOE(1)=2
      NOEASY(1)=0.0D0
      NOESWI(1)=10.0D0
      NOEMSO(1)=2
      NOEMAS(1)=0.0D0
      NOEMSW(1)=10.0D0
      NOEICV=0
      VIOSTR = 0
      VIOEXC = 0.0D0
      VIOTHR = 0.0D0
      MATDIM = 0
      AMBLEV = 1.0D0
      AMBMAX = 1
      AMBMIN = 1
      AMBSTR = 0
      AMBCUT = R4BIG
      AVEMOD = 'AVER'
      AMBMOD = 'ALL '
      AMBFRM = 'OR-R'
      CALCNT=0
      CALDIS=6.0D0
      NMONO(1)=1
      NAVEXP(1)=6.0D0
      OREXP(1)=6.0D0
C
C reset time and running-average stuff
      DAVEXP(1)=3
      RAVEXP(1)=3
      MEMTAU(1)=1000.

      TMODE(1)='OFF'
      RMODE(1)='OFF'
      TFMODE(1)='CONS'
      RANA='CURR'
C
C reset ensemble-average stuff
      IDIMER=0
      DIENS=' '
      EMODE='OFF'
      OCCU=.FALSE.
C
      RETURN
      END
C====================================================================
      SUBROUTINE NOEHP(NOENEW)
C
C routine allocates space for NOE restraints list on the HEAP
C
C Author: Axel T. Brunger
C =======================
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'noe.inc'
C mtf.inc necessary for ARIA
      INCLUDE 'mtf.inc'
      INTEGER NOENEW
      LOGICAL HEAPDEBUG
C begin
C
C check whether space was already allocated
      HEAPDEBUG=.FALSE.
      IF (HEAPDEBUG) THEN
      WRITE (6,'(A,I10,I10)') 'HPNORR', HPNORR, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNIPR', HPNIPR, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNILS', HPNILS, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNJPR', HPNJPR, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNJLS', HPNJLS, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNCND', HPNCND, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNRAV', HPNRAV, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNRRV', HPNRRV, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNNSP', HPNNSP, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNDIS', HPNDIS, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPDENINIT', HPDENINIT, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNLOW', HPNLOW, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNHIG', HPNHIG, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNWGH', HPNWGH, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNCV', HPNCV, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNVIO', HPNVIO, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNPID', HPNPID, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNPP1', HPNPP1, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNPP2', HPNPP2, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNVOL', HPNVOL, NOEMAX
      WRITE (6,'(A,I10,I10)') 'HPNCAL', HPNCAL, NATOM
      WRITE (6,'(A,I10,I10)') 'HPNDIA', HPNDIA, NATOM
      WRITE (6,'(A,I10,I10)') 'HPNPSE', HPNPSE, NATOM
      WRITE (6,'(A,I10,I10)') 'HPNHGL', HPNHGL, NATOM
      WRITE (6,'(A,I10,I10)') 'HPNMA2', HPNMA2, MATDIM
      WRITE (6,'(A,I10,I10)') 'HPNMAT', HPNMAT, MATDIM
      END IF

      IF (HPNORR.NE.0) CALL FREHP(HPNORR,INTEG4(NOEMAX))
      IF (HPNIPR.NE.0) CALL FREHP(HPNIPR,INTEG4(NOEMAX))
      IF (HPNILS.NE.0) CALL FREHP(HPNILS,INTEG4(NOEMAX))
      IF (HPNJPR.NE.0) CALL FREHP(HPNJPR,INTEG4(NOEMAX))
      IF (HPNJLS.NE.0) CALL FREHP(HPNJLS,INTEG4(NOEMAX))
      IF (HPNCND.NE.0) CALL FREHP(HPNCND,INTEG4(NOEMAX))
      IF (HPNRAV.NE.0) CALL FREHP(HPNRAV,IREAL8(NOEMAX))
      IF (HPNRRV.NE.0) CALL FREHP(HPNRRV,IREAL8(NOEMAX))
      IF (HPNNSP.NE.0) CALL FREHP(HPNNSP,INTEG4(NOEMAX))
      IF (HPNDIS.NE.0) CALL FREHP(HPNDIS,IREAL8(NOEMAX))
      IF (HPDENINIT.NE.0) CALL FREHP(HPDENINIT,IREAL8(NOEMAX))
      IF (HPNLOW.NE.0) CALL FREHP(HPNLOW,IREAL8(NOEMAX))
      IF (HPNHIG.NE.0) CALL FREHP(HPNHIG,IREAL8(NOEMAX))
      IF (HPNWGH.NE.0) CALL FREHP(HPNWGH,IREAL8(NOEMAX))
      IF (HPNCV.NE.0)  CALL FREHP(HPNCV, INTEG4(NOEMAX))
      IF (HPNVIO.NE.0) CALL FREHP(HPNVIO,INTEG4(NOEMAX))
      IF (HPNPID.NE.0) CALL FREHP(HPNPID,INTEG4(NOEMAX))
      IF (HPNPP1.NE.0) CALL FREHP(HPNPP1,IREAL8(NOEMAX))
      IF (HPNPP2.NE.0) CALL FREHP(HPNPP2,IREAL8(NOEMAX))
      IF (HPNVOL.NE.0) CALL FREHP(HPNVOL,IREAL8(NOEMAX))
C free arrays for ARIA
      IF (HPNCAL.NE.0) CALL FREHP(HPNCAL,INTEG4(NATOM))
      IF (HPNDIA.NE.0) CALL FREHP(HPNDIA,IREAL8(NATOM))
      IF (HPNPSE.NE.0) CALL FREHP(HPNPSE,IREAL8(NATOM))
      IF (HPNHGL.NE.0) CALL FREHP(HPNHGL,INTEG4(NATOM))
      IF (HPNMA2.NE.0) CALL FREHP(HPNMA2,INTEG4(MATDIM))
      IF (HPNMAT.NE.0) CALL FREHP(HPNMAT,IREAL8(MATDIM))
C
      NOENUM=0
      NOEORN=0
      NOEMAX=0
      NPEAKI=0
C
      HPNORR=0
      HPNIPR=0
      HPNILS=0
      HPNJPR=0
      HPNJLS=0
      HPNCND=0
      HPNRAV=0
      HPNRRV=0
      HPNNSP=0
      HPNDIS=0
      HPDENINIT=0
      HPNLOW=0
      HPNHIG=0
      HPNWGH=0
      HPNCV=0
      HPNVIO=0
      HPNPID=0
      HPNPP1=0
      HPNPP2=0
      HPNVOL=0
      HPNCAL=0
      HPNPSE=0
      HPNDIA=0
      HPNHGL=0
      HPNMAT=0
      HPNMA2=0
C
C now allocate new space
      IF (NOENEW.GT.0) THEN
      NOEMAX=NOENEW
      HPNVOL=ALLHP(IREAL8(NOENEW))
      HPNPP2=ALLHP(IREAL8(NOENEW))
      HPNPP1=ALLHP(IREAL8(NOENEW))
      HPNPID=ALLHP(INTEG4(NOENEW))
      HPNVIO=ALLHP(INTEG4(NOENEW))
      HPNCV=ALLHP(INTEG4(NOENEW))
      HPNWGH=ALLHP(IREAL8(NOENEW))
      HPNHIG=ALLHP(IREAL8(NOENEW))
      HPNLOW=ALLHP(IREAL8(NOENEW))
      HPNRAV=ALLHP(IREAL8(NOENEW))
      HPNRRV=ALLHP(IREAL8(NOENEW))
      HPNNSP=ALLHP(INTEG4(NOENEW))
      HPNDIS=ALLHP(IREAL8(NOENEW))
      HPDENINIT=ALLHP(IREAL8(NOENEW))
      HPNCND=ALLHP(INTEG4(NOENEW))
      HPNJLS=ALLHP(INTEG4(NOENEW))
      HPNJPR=ALLHP(INTEG4(NOENEW))
      HPNILS=ALLHP(INTEG4(NOENEW))
      HPNIPR=ALLHP(INTEG4(NOENEW))
      HPNORR=ALLHP(INTEG4(NOENEW))
      END IF
      RETURN
      END
C====================================================================
      SUBROUTINE NOECVS(NOENUM,NOECV,PART)
C
C Routine partitions NOE data into PART sets.
C NOECV will contain integer numbers between 1 and PART.
C Default value for NOECV is -1 (needed to allow partitioning of 
C single classes only) AB April 2004
C
C Author: Axel T. Brunger
C
C
C     IMPLICIT NONE
C I/O
      INTEGER NOENUM, NOECV(*)
      INTEGER PART
C local
      INTEGER I, P, NP, NRETRY, NTRYTOT
      DOUBLE PRECISION RNUM
      NTRYTOT = 0
C begin
  100 CONTINUE
      NRETRY = 0
      IF (PART.GT.0 .AND. NOENUM .GT.0) THEN
      DO I=1,NOENUM
      CALL GGUBFS(RNUM)
      NOECV(I)=MAX(1,MIN(PART,INT(RNUM*PART)+1))
      END DO
C
      IF (PART .EQ. NOENUM) THEN
      DO I=1,NOENUM
      NOECV(I)=I
      END DO
      END IF
C
      DO P=1,PART
      NP=0
      DO I=1,NOENUM
      IF (NOECV(I).EQ.P) THEN
      NP=NP+1
      END IF
      END DO
      IF (NP .EQ. 0) THEN
        NRETRY = 1
      ENDIF
      WRITE(6,'(A,I3,A,I5,A)') ' For set ',P,
     & ' there are ',NP,' distance restraints.'
      END DO
      ELSE
      WRITE(6,'(A)')
     & ' Data are not partitioned or partitioning removed.'
      DO I=1,NOENUM
      NOECV(I)=-1
      END DO
      END IF
      NTRYTOT = NTRYTOT + NRETRY
      IF (NRETRY .GT. 0 .AND. NTRYTOT .LE. 10) THEN
      WRITE(6,'(A)')
     & ' Test set with 0 constraints! New trial...'
      GOTO 100
      ELSE IF (NTRYTOT .GT. 20) THEN
      CALL WRNDIE(-1,'NOECVS',
     & 'Unable to correctly partition the NOE data within 20 trials')
      ENDIF
C
      RETURN
      END
C====================================================================
      SUBROUTINE ENSBLE(WENS,I,J)
C
C  set weighting factors for ensemble-averaged restraints
C  to keep only intra model distances
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'coord.inc'
      DOUBLE PRECISION WENS
      INTEGER I, J, L
C
      IF (SEGID(I).NE.SEGID(J)) THEN
      WENS = -1.0
      IF (DIENS.EQ.'DIME') THEN
      DO L = 1, IDIMER
      IF ((SEGID(I).EQ.ENSMOD(L,1).AND.
     &     SEGID(J).EQ.ENSMOD(L,2)).OR.
     &    (SEGID(I).EQ.ENSMOD(L,2).AND.
     &     SEGID(J).EQ.ENSMOD(L,1))) THEN
      WENS = MIN(QMAIN(I),QMAIN(J))
      END IF
      END DO
      END IF
      ELSE
      WENS = MIN(QMAIN(I),QMAIN(J))
      END IF
      IF (.NOT. OCCU .AND. WENS .GE. 0.0) WENS = 1.0
C      WRITE(6,'(A,A5,A,A5,F5.2)')
C     & 'SEGID I=',SEGID(I),' SEGID J=',SEGID(J),WENS
      RETURN
      END
C


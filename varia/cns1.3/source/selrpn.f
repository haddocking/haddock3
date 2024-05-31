      SUBROUTINE SELCTA(FLAGS,NSELCT,X,Y,Z,QCOOR)
C
C Recursive atom selection routine.
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C =======================
C
C Parameter MAXDEP specifies the maximum depth of an expression
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INTEGER FLAGS(*)
      INTEGER NSELCT
      DOUBLE PRECISION X(*), Y(*), Z(*)
      LOGICAL QCOOR
C local
      INTEGER RETADR, QSTACK
      INTEGER RETMAX, MAXDEP, STKSIZ
      PARAMETER (RETMAX=100, MAXDEP=10)
C begin
      STKSIZ=MAX(MAXDEP*NATOM,MAXDEP)
      RETADR=ALLHP(INTEG4(RETMAX))
      QSTACK=ALLHP(ILOGIC(STKSIZ))
      CALL SELCT2(RETMAX,HEAP(RETADR),FLAGS,NSELCT,QCOOR,X,Y,Z,
     &            STKSIZ,HEAP(QSTACK))
      CALL FREHP(QSTACK,ILOGIC(STKSIZ))
      CALL FREHP(RETADR,INTEG4(RETMAX))
      END
C
      SUBROUTINE SELCT2(RETMAX,RETADR,FLAGS,NSELCT,QCOOR,
     &           XC,YC,ZC,STKSIZ,QSTACK)
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'ncs.inc'
      INTEGER RETMAX, RETADR(*), FLAGS(*), NSELCT
      LOGICAL QCOOR
      DOUBLE PRECISION XC(*), YC(*), ZC(*)
      INTEGER STKSIZ
      LOGICAL QSTACK(*)
C local
      INTEGER QSTKLN, QSTART, QINDX1, QINDX2, QLIFT
      INTEGER ISTK
      LOGICAL ELEMNT, OK, ERR, QKNOWN, XRERR
      LOGICAL QRANGE, QABS, QWCARD, QWCARD1, QWCARD2, QWCARD3, EXIST
      LOGICAL CLOOP
      INTEGER IAT, IATS1
      INTEGER IAT2, NUMB, NUMB1, NUMB2
      INTEGER I, INDX1, IGRP, IBT, JBT
      INTEGER J, IDLEN, K, IATS2, IIAT
      INTEGER LEVEL, ISYM, INSYM
      DOUBLE PRECISION RMAX, RMAX2, DIST2, DXL, DYL, DZL, R(3)
      DOUBLE PRECISION XREF, YREF, ZREF
C MODIFICATION: selection token lengths doubled for wildcards with
C backslash escapes, and for input of numeric resSeq > 9999.
      CHARACTER*8 SEG1,  RES1,  RESN1,  RESID1,  SEG2,  RESN2
      INTEGER LSEG1, LRESN1, LRESD1, LSEG2, LRESN2
      CHARACTER*8 RESID2, TYPE1, TYPE2
      INTEGER LRESD2,LTYPE1,LTYPE2
      INTEGER LRESD
      CHARACTER*2 ALT1, ALT2
      INTEGER LALT1, LALT2
      CHARACTER*1 ALPH, ALPH1, ALPH2
      CHARACTER*1 SBRA, SKET
      PARAMETER (SBRA='(',SKET=')')
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
      DOUBLE PRECISION XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
      DOUBLE PRECISION XRRF, YRRF, ZRRF, XRRFF, YRRFF, ZRRFF
C pointer
      INTEGER XF, YF, ZF, XFF, YFF, ZFF, XFFF, YFFF, ZFFF
C parameter
      DOUBLE PRECISION ZERO, ONE, BIG
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, BIG=10000.0D0)
C
C begin
C
C initialize working stack
      QSTART=1
      QSTKLN=NATOM
      QLIFT=0
C
      LEVEL=1
C
      CALL SENEXT
      IF (WD(1:4).EQ.'=   ') CALL SENEXT
      IF (WD(1:2).EQ.'? ') THEN
C
C-begin-show-the-current-set
C
C This procedure prints a listing of the atoms that are currently
C set.
      WRITE(6,7001)
7001  FORMAT(' The following atoms are currently set:')
C
      SEG1='    '
      RES1='    '
      RESN1='    '
      J=0
      DO IAT=1,NATOM
      IF (SEG1.NE.SEGID(IAT).OR.RES1.NE.RESID(IAT).OR.RESN1.NE.RES(IAT))
     1 THEN
      SEG1=SEGID(IAT)
      RES1=RESID(IAT)
      RESN1=RES(IAT)
      IF (J.GT.0) WRITE(6,'(2A)') ' ',WD(1:WDLEN)
      WDLEN=0
      CALL ADDST(WD,WDMAX,WDLEN,SEGID(IAT),4)
      CALL ADDST(WD,WDMAX,WDLEN,' ',1)
      CALL ADDST(WD,WDMAX,WDLEN,RESID(IAT),4)
      CALL ADDST(WD,WDMAX,WDLEN,' ',1)
      CALL ADDST(WD,WDMAX,WDLEN,RES(IAT),4)
      IDLEN = WDLEN
      J=0
      END IF
      IF (FLAGS(IAT).EQ.1) THEN
      IF (WDLEN+1+4 .GT. 78) THEN
      WRITE(6,'(2A)') ' ',WD(1:WDLEN)
      J = 0
      WDLEN = 0
      DO K=1,IDLEN
      CALL ADDST(WD,WDMAX,WDLEN,' ',1)
      END DO
      END IF
      CALL ADDST(WD,WDMAX,WDLEN,' ',1)
      CALL ADDST(WD,WDMAX,WDLEN,TYPE(IAT),4)
      J=J+1
      END IF
      END DO
      IF (J.GT.0) WRITE(6,'(2A)') ' ',WD(1:WDLEN)
C-end-show-the-current-set
C
C jump to exit
      GOTO 99999
C
      END IF
C
C check that expression begins with a "("
      IF (.NOT.WD(1:1).EQ.SBRA) THEN
      CALL DSPERR('SELRPN','selection has to begin with a "("')
      CALL WRNDIE(-5,'SELRPN','parsing error')
C jump to exit
      GOTO 99999
      END IF
C
      CALL SENEXT
C ** invoke procedure expression **
      RETADR(LEVEL)=1
      GOTO 1110
1111  CONTINUE
C ** return label                **
C
C
      IF (.NOT.WD(1:1).EQ.SKET) THEN
      CALL DSPERR('SELRPN','selection has to end with a ")"')
      CALL SAVEWD
      END IF
      DONE=.FALSE.
C
C ** return      **
      GOTO 9999
C ** return      **
C
C
C---- BEGIN PROCEDURE EXPRESSION --------------------------------------
1110  CONTINUE
      CALL SEDOWN(LEVEL,RETMAX)
C
C ** invoke procedure term    **
      RETADR(LEVEL)=3
      GOTO 2220
2223  CONTINUE
C ** return label             **
C
2233  IF (WD(1:2).NE.'OR') GOTO 1122
      CALL SENEXT
C
C ** invoke procedure term    **
      RETADR(LEVEL)=4
      GOTO 2220
2224  CONTINUE
C ** return label             **
C
C-begin-process-or-operator
      QINDX1=QSTART
      QINDX2=QSTART-QSTKLN
      DO IAT=1,NATOM
      QSTACK(QINDX2)=QSTACK(QINDX1).OR.QSTACK(QINDX2)
      QINDX1=QINDX1+1
      QINDX2=QINDX2+1
      END DO
      CALL SEBACK(QLIFT,QSTART,QSTKLN,ISTK)
C-end-process-or-operator
C
      GOTO 2233
1122  CONTINUE
C
C
      CALL SEUP(LEVEL,RETMAX)
C
C     return to address:
      IF (RETADR(LEVEL).EQ.1) GOTO 1111
      IF (RETADR(LEVEL).EQ.2) GOTO 1112
      WRITE(6,'(A)') 
     & ' %SELRPN-ERR: Return address unknown. Fatal coding error'
      CALL DIE
C---- END PROCEDURE EXPRESSION ---------------------------------------
C
C
C
C---- BEGIN PROCEDURE TERM -------------------------------------------
2220  CONTINUE
      CALL SEDOWN(LEVEL,RETMAX)
C
C ** invoke procedure factor   **
      RETADR(LEVEL)=5
      GOTO 3330
3335  CONTINUE
C ** return label              **
C
4455  IF (WD(1:3).NE.'AND') GOTO 5566
      CALL SENEXT
C
C ** invoke procedure factor   **
      RETADR(LEVEL)=6
      GOTO 3330
3336  CONTINUE
C ** return label              **
C
C-begin-process-and-operator
      QINDX1=QSTART
      QINDX2=QSTART-QSTKLN
      DO IAT=1,NATOM
      QSTACK(QINDX2)=QSTACK(QINDX1).AND.QSTACK(QINDX2)
      QINDX1=QINDX1+1
      QINDX2=QINDX2+1
      END DO
      CALL SEBACK(QLIFT,QSTART,QSTKLN,ISTK)
C-end-process-and-operator
C
      GOTO 4455
5566  CONTINUE
C
      CALL SEUP(LEVEL,RETMAX)
C
C     return to address:
      IF (RETADR(LEVEL).EQ.3) GOTO 2223
      IF (RETADR(LEVEL).EQ.4) GOTO 2224
      WRITE(6,'(A)') 
     & ' %SELRPN-ERR: Return address unknown. Fatal coding error.'
      CALL DIE
C---- END PROCEDURE TERM ----------------------------------------------
C
C
C
C
C
C---- BEGIN PROCEDURE FACTOR ------------------------------------------
3330  CONTINUE
      CALL SEDOWN(LEVEL,RETMAX)
C=NOT===================================================================
      IF (WD(1:3).NE.'NOT') GOTO 10
C
      CALL SENEXT
C
C ** invoke procedure factor  **
      RETADR(LEVEL)=7
      GOTO 3330
3337  CONTINUE
C ** return label             **
C
C-begin-process-not-operator
      QINDX1=QSTART
      DO IAT=1,NATOM
      QSTACK(QINDX1)=.NOT.QSTACK(QINDX1)
      QINDX1=QINDX1+1
      END DO
C-end-process-not-operator
C
      GOTO 77777
C=BYRE===================================================================
10    IF (WD(1:4).NE.'BYRE') GOTO 20
C
      CALL SENEXT
C
C ** invoke procedure factor  **
      RETADR(LEVEL)=8
      GOTO 3330
3338  CONTINUE
C ** return label             **
C
C-begin-process-byres-operator
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO INDX1=QSTART,QSTART+QSTKLN-1
      QSTACK(INDX1)=QSTACK(INDX1-QSTKLN)
      END DO
      DO INDX1=QSTART,QSTART+QSTKLN-1
      QSTACK(INDX1-QSTKLN)=.FALSE.
      END DO
      DO IAT2=1,NATOM
      IF (QSTACK(QSTART+IAT2-1)) THEN
      SEG1=SEGID(IAT2)
      RES1=RESID(IAT2)
      RESN1=RES(IAT2)
      DO IAT=1,NATOM
      IF (SEGID(IAT).EQ.SEG1) THEN
      IF (RESID(IAT).EQ.RES1.AND.RES(IAT).EQ.RESN1) THEN
      QSTACK(QSTART-QSTKLN+IAT-1)=.TRUE.
      END IF
      END IF
      END DO
      END IF
      END DO
      CALL SEBACK(QLIFT,QSTART,QSTKLN,ISTK)
C-end-process-byres-operator
C
      GOTO 77777
C=BYGR==============================================================
20    IF (WD(1:4).NE.'BYGR') GOTO 30
C
      CALL SENEXT
C
C ** invoke procedure factor  **
      RETADR(LEVEL)=9
      GOTO 3330
3339  CONTINUE
C ** return label             **
C
C-begin-process-bygroup-operator
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO INDX1=QSTART,QSTART+QSTKLN-1
      QSTACK(INDX1)=QSTACK(INDX1-QSTKLN)
      END DO
      DO INDX1=QSTART,QSTART+QSTKLN-1
      QSTACK(INDX1-QSTKLN)=.FALSE.
      END DO
      DO IGRP=1,NGRP
      IAT2=IGPBS(IGRP)
      ELEMNT=.FALSE.
      CLOOP = .TRUE.
      DO WHILE (CLOOP)
      IAT2=IAT2+1
      IF (QSTACK(QSTART+IAT2-1)) THEN
      DO IAT=IGPBS(IGRP)+1,IGPBS(IGRP+1)
      QSTACK(QSTART-QSTKLN+IAT-1)=.TRUE.
      END DO
      ELEMNT=.TRUE.
      END IF
      IF (ELEMNT.OR.IAT2.GE.IGPBS(IGRP+1)) CLOOP = .FALSE.
      END DO
      END DO
      CALL SEBACK(QLIFT,QSTART,QSTKLN,ISTK)
C-end-process-bygroup-operator
C
      GOTO 77777
C=BOUN=========================================================
30    IF (WD(1:4).NE.'BOND') GOTO 40
C
      CALL SENEXT
C
C ** invoke procedure factor  **
      RETADR(LEVEL)=10
      GOTO 3330
3340  CONTINUE
C ** return label             **
C
C BONDedto operator added by Michael Nilges, EMBL
C-begin-process-bondedto-operator
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO INDX1=QSTART,QSTART+QSTKLN-1
      QSTACK(INDX1)=QSTACK(INDX1-QSTKLN)
      END DO
      DO INDX1=QSTART,QSTART+QSTKLN-1
      QSTACK(INDX1-QSTKLN)=.FALSE.
      END DO
      DO I=1,NBOND
      IBT=IB(I)
      JBT=JB(I)
      IF (IBT.GT.0 .AND. JBT.GT.0) THEN
      IF(QSTACK(QSTART+JBT-1)) THEN
      QSTACK(QSTART-QSTKLN+IBT-1)=.TRUE.
      END IF
      IF(QSTACK(QSTART+IBT-1)) THEN
      QSTACK(QSTART-QSTKLN+JBT-1)=.TRUE.
      END IF
      END IF
      END DO
      CALL SEBACK(QLIFT,QSTART,QSTKLN,ISTK)
C-end-process-boundto-operator
C
      GOTO 77777
C=(============================================================
40    IF (WD(1:1).NE.SBRA) GOTO 50
C
      CALL SENEXT
C
C ** invoke procedure expression **
      RETADR(LEVEL)=2
      GOTO 1110
1112  CONTINUE
C ** return label                **
      IF (.NOT.WD(1:1).EQ.SKET) THEN
      WRITE(6,'(A)') ' %SELRPN-ERR: unbalanced parentheses.'
      CALL WRNDIE(-5,'SELRPN','parsing error')
C exit selection subroutine
      GOTO 99999
      END IF
      CALL SENEXT
      GOTO 77777
C==============================================================
C==========
C ALL token
C==========
50    IF (WD(1:3).EQ.'ALL') THEN
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      ELEMNT=.TRUE.
      QSTACK(QSTART+IAT-1)=ELEMNT
      END DO
      CALL SENEXT
C==========
C NONE token
C==========
      ELSE IF (WD(1:4).EQ.'NONE') THEN
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      ELEMNT=.FALSE.
      QSTACK(QSTART+IAT-1)=ELEMNT
      END DO
      CALL SENEXT
C
C===========
C PREV token
C===========
      ELSE IF (WD(1:4).EQ.'PREV') THEN
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      ELEMNT=(FLAGS(IAT).EQ.1)
      QSTACK(QSTART+IAT-1)=ELEMNT
      END DO
      CALL SENEXT
C
C=============
C RECALL token
C=============
      ELSE IF (WD(1:4).EQ.'RECA'.OR.WD(1:4).EQ.'STOR') THEN
C fill default (FALSE)
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=.FALSE.
      END DO
C
      IF (WD(1:7).EQ.'RECALL1'.OR.WD(1:6).EQ.'STORE1') THEN
      J=1
      ELSE IF (WD(1:7).EQ.'RECALL2'.OR.WD(1:6).EQ.'STORE2') THEN
      J=2
      ELSE IF (WD(1:7).EQ.'RECALL3'.OR.WD(1:6).EQ.'STORE3') THEN
      J=3
      ELSE IF (WD(1:7).EQ.'RECALL4'.OR.WD(1:6).EQ.'STORE4') THEN
      J=4
      ELSE IF (WD(1:7).EQ.'RECALL5'.OR.WD(1:6).EQ.'STORE5') THEN
      J=5
      ELSE IF (WD(1:7).EQ.'RECALL6'.OR.WD(1:6).EQ.'STORE6') THEN
      J=6
      ELSE IF (WD(1:7).EQ.'RECALL7'.OR.WD(1:6).EQ.'STORE7') THEN
      J=7
      ELSE IF (WD(1:7).EQ.'RECALL8'.OR.WD(1:6).EQ.'STORE8') THEN
      J=8
      ELSE IF (WD(1:7).EQ.'RECALL9'.OR.WD(1:6).EQ.'STORE9') THEN
      J=9
      ELSE
      CALL NEXTI('RECALL=',J)
      END IF
      IF (.NOT. (J.GE.1.AND.J.LE.MAXSTO)) THEN
      CALL WRNDIE(-5,'SELRPN','Store number out of range')
      ELSE IF (PTRSTO(J).EQ.0) THEN
      WRITE(6,'(A)') ' SELRPN: Store empty. No atoms selected.'
      ELSE IF (LENSTO(J).NE.NATOM) THEN
      CALL WRNDIE(-5,'SELRPN','dim. mismatch')
      ELSE
      CALL SELRCL(QSTACK(QSTART),HEAP(PTRSTO(J)),NATOM)
      END IF
      CALL SENEXT
C
C=========
C ID token
C=========
      ELSE IF (WD(1:4).EQ.'BYNU'.OR.WD(1:2).EQ.'ID') THEN
      CALL SENEXT
      IATS1=DECODI(WD,WDLEN,OK)
      ERR=.FALSE.
      IF (.NOT.OK.OR.IATS1.LE.0.OR.IATS1.GT.NATOM) THEN
      WRITE(6,'(2A)') ' %SELRPN-BYNU-ERR: no valid number ',WD(1:WDLEN)
      CALL WRNDIE(-5,'SELRPN','no valid number')
      ERR=.TRUE.
      END IF
      IATS2=IATS1
      CALL SENEXT
      QRANGE=(WD(1:1).EQ.':')
      IF (QRANGE) THEN
      CALL SENEXT
      IATS2=DECODI(WD,WDLEN,OK)
      IF (.NOT.OK.OR.IATS2.LE.0.OR.IATS2.GT.NATOM) THEN
      WRITE(6,'(2A)') ' %SELRPN-BYNU-ERR: no valid number ',WD(1:WDLEN)
      CALL WRNDIE(-5,'SELRPN','no valid number')
      ERR=.TRUE.
      END IF
      CALL SENEXT
      END IF
      IF (.NOT.ERR) THEN
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=.FALSE.
      END DO
      DO IAT=IATS1,IATS2
      QSTACK(QSTART+IAT-1)=.TRUE.
      END DO
      END IF
C
C============
C SEGID token
C============
      ELSE IF (WD(1:4).EQ.'SEGI') THEN
      CALL NEXTST('segment-id=',WD)
C
      QWCARD=.NOT.QQUOT
C
      CALL COPYST(SEG1,LEN(SEG1),LSEG1,WD,WDLEN)
      CALL SENEXT
      QRANGE=(WD(1:1).EQ.':')
      IF (QRANGE) THEN
      CALL NEXTST('segment-id=',WD)
      CALL COPYST(SEG2,LEN(SEG2),LSEG2,WD,WDLEN)
      CALL SENEXT
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=
     &    (SEG1(1:LSEG1).LE.SEGID(IAT).AND.SEGID(IAT).LE.SEG2(1:LSEG2))
      END DO
      ELSE
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
C
      IF (QWCARD) THEN
      CALL EQSTWC(SEGID,4,SEG1,LSEG1,QSTART,NATOM,QSTACK)
      ELSE
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=SEG1(1:LSEG1).EQ.SEGID(IAT)
      END DO
      END IF
C
      END IF
C
C============
C ALTID token
C============
      ELSE IF (WD(1:4).EQ.'ALTI') THEN
      CALL NEXTST('alternate-id=',WD)
C
      QWCARD=.NOT.QQUOT
C
      CALL COPYST(ALT1,LEN(ALT1),LALT1,WD,WDLEN)
      CALL SENEXT
      QRANGE=(WD(1:1).EQ.':')
      IF (QRANGE) THEN
      CALL NEXTST('alternate-id=',WD)
      CALL COPYST(ALT2,LEN(ALT1),LALT2,WD,WDLEN)
      CALL SENEXT
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=
     &    (ALT1(1:LALT1).LE.ALTID(IAT).AND.ALTID(IAT).LE.ALT2(1:LALT2))
      END DO
      ELSE
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
C
      IF (QWCARD) THEN
      CALL EQSTWC(ALTID,1,ALT1,LALT1,QSTART,NATOM,QSTACK)
      ELSE
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=ALT1(1:LALT1).EQ.ALTID(IAT)
      END DO
      END IF
C
      END IF
C
C===========
C RESN token
C===========
      ELSE IF (WD(1:4).EQ.'RESN') THEN
      CALL NEXTST('residue-name=',WD)
C
      QWCARD=.NOT.QQUOT
C
      CALL COPYST(RESN1,LEN(RESN1),LRESN1,WD,WDLEN)
      CALL SENEXT
      QRANGE=(WD(1:1).EQ.':')
      IF (QRANGE) THEN
      CALL NEXTST('residue-name=',WD)
      CALL COPYST(RESN2,LEN(RESN2),LRESN2,WD,WDLEN)
      CALL SENEXT
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=
     &   (RESN1(1:LRESN1).LE.RES(IAT).AND.RES(IAT).LE.RESN2(1:LRESN2))
      END DO
      ELSE
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
C
      IF (QWCARD) THEN
      CALL EQSTWC(RES,4,RESN1,LRESN1,QSTART,NATOM,QSTACK)
      ELSE
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=RESN1(1:LRESN1).EQ.RES(IAT)
      END DO
      END IF
C
      END IF
C
C============
C RESID token
C============
      ELSE IF (WD(1:4).EQ.'RESI') THEN
      CALL NEXTST('residue-id=',WD)
C
      QWCARD=.NOT.QQUOT
C
      CALL COPYST(RESID1,LEN(RESID1),LRESD1,WD,WDLEN)
      CALL SENEXT
      QRANGE=(WD(1:1).EQ.':')
      IF (QRANGE) THEN
      CALL SPLITI(RESID1,NUMB1,ALPH1,LRESD1,OK)
C
      IF (.NOT.OK) THEN
      WRITE(6,'(2A)') ' %SELRPN-ERR: invalid RESID ',RESID1
      CALL WRNDIE(-5,'SELRPN','invalid residue specification')
      ERR=.TRUE.
      END IF
C
      CALL NEXTST('residue-id=',WD)
      CALL COPYST(RESID2,LEN(RESID2),LRESD2,WD,WDLEN)
      CALL SPLITI(RESID2,NUMB2,ALPH2,LRESD2,OK)
C
      IF (.NOT.OK) THEN
      WRITE(6,'(2A)') ' %SELRPN-ERR: invalid RESID ',RESID2
      CALL WRNDIE(-5,'SELRPN','invalid residue specification')
      ERR=.TRUE.
      END IF
C
      CALL SENEXT
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      LRESD=4
      CALL SPLITI(RESID(IAT),NUMB,ALPH,LRESD,OK)
C
      IF (.NOT.OK) THEN
      WRITE(6,'(3A,I9)') 
     & ' %SELRPN-ERR: invalid RESId ',RESID(IAT),' for atom ',IAT
      CALL WRNDIE(-5,'SELRPN','invalid RESID.')
      ERR=.TRUE.
      END IF
C
      QSTACK(QSTART+IAT-1)=
     &      (NUMB1.LT.NUMB.OR.
     &      (NUMB1.EQ.NUMB.AND.ALPH1.LE.ALPH))
     &       .AND.
     &      (NUMB.LT.NUMB2.OR.
     &      (NUMB.EQ.NUMB2.AND.ALPH.LE.ALPH2))
      END DO
      ELSE
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
C
      IF (QWCARD) THEN
      CALL EQSTWC(RESID,4,RESID1,LRESD1,QSTART,NATOM,QSTACK)
      ELSE
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=RESID1(1:LRESD1).EQ.RESID(IAT)
      END DO
      END IF
C
      END IF
C
C===========
C NAME token
C===========
      ELSE IF (WD(1:4).EQ.'TYPE'.OR.WD(1:4).EQ.'NAME') THEN
      CALL NEXTST('atom-name=',WD)
C
      QWCARD=.NOT.QQUOT
C
      CALL COPYST(TYPE1,LEN(TYPE1),LTYPE1,WD,WDLEN)
      CALL SENEXT
      QRANGE=(WD(1:1).EQ.':')
      IF (QRANGE) THEN
      CALL NEXTST('atom-name=',WD)
      CALL COPYST(TYPE2,LEN(TYPE2),LTYPE2,WD,WDLEN)
      CALL SENEXT
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=
     &  TYPE1(1:LTYPE1).LE.TYPE(IAT).AND.TYPE(IAT).LE.TYPE2(1:LTYPE2)
      END DO
      ELSE
C
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
C
      IF (QWCARD) THEN
      CALL EQSTWC(TYPE,4,TYPE1,LTYPE1,QSTART,NATOM,QSTACK)
      ELSE
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=TYPE1(1:LTYPE1).EQ.TYPE(IAT)
      END DO
      END IF
C
      END IF
C
C===========
C CHEM token
C===========
      ELSE IF (WD(1:4).EQ.'CHEM') THEN
      CALL NEXTST('chemical-type=',WD)
C
      QWCARD=.NOT.QQUOT
C
      CALL COPYST(TYPE1,LEN(TYPE1),LTYPE1,WD,WDLEN)
      CALL SENEXT
      QRANGE=(WD(1:1).EQ.':')
      IF (QRANGE) THEN
      CALL NEXTST('chemical-type=',WD)
      CALL COPYST(TYPE2,LEN(TYPE2),LTYPE2,WD,WDLEN)
      CALL SENEXT
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=
     &   TYPE1(1:LTYPE1).LE.IAC(IAT).AND.IAC(IAT).LE.TYPE2(1:LTYPE2)
      END DO
      ELSE
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
C
      IF (QWCARD) THEN
      CALL EQSTWC(IAC,4,TYPE1,LTYPE1,QSTART,NATOM,QSTACK)
      ELSE
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=TYPE1(1:LTYPE1).EQ.IAC(IAT)
      END DO
      END IF
C
      END IF
C
C===========
C ATOM token
C===========
      ELSE IF (WD(1:4).EQ.'ATOM') THEN
      CALL NEXTST('segment-id=',WD)
C
      QWCARD1=.NOT.QQUOT
C
      CALL COPYST(SEG1,LEN(SEG1),LSEG1,WD,WDLEN)
      CALL NEXTST('residue-id=',WD)
C
      QWCARD2=.NOT.QQUOT
C
      CALL COPYST(RES1,LEN(RES1),LRESD1,WD,WDLEN)
      CALL NEXTST('atom-name=',WD)
C
      QWCARD3=.NOT.QQUOT
C
      CALL COPYST(WDT,WDTMAX,WDTLEN,WD,WDLEN)
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
C
      IF (QWCARD1) THEN
      CALL EQSTWC(SEGID,4,SEG1,LSEG1,QSTART,NATOM,QSTACK)
      ELSE
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=SEG1(1:LSEG1).EQ.SEGID(IAT)
      END DO
      END IF
C
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
C
      IF (QWCARD2) THEN
      CALL EQSTWC(RESID,4,RES1,LRESD1,QSTART,NATOM,QSTACK)
      ELSE
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=RES1(1:LRESD1).EQ.RESID(IAT)
      END DO
      END IF
C
      QINDX1=QSTART
      QINDX2=QSTART-QSTKLN
      DO IAT=1,NATOM
      QSTACK(QINDX2)=QSTACK(QINDX1).AND.QSTACK(QINDX2)
      QINDX1=QINDX1+1
      QINDX2=QINDX2+1
      END DO
      CALL SEBACK(QLIFT,QSTART,QSTKLN,ISTK)
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
C
      IF (QWCARD3) THEN
      CALL EQSTWC(TYPE,4,WDT,WDTLEN,QSTART,NATOM,QSTACK)
      ELSE
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=WDT(1:WDTLEN).EQ.TYPE(IAT)
      END DO
      END IF
C
      QINDX1=QSTART
      QINDX2=QSTART-QSTKLN
      DO IAT=1,NATOM
      QSTACK(QINDX2)=QSTACK(QINDX1).AND.QSTACK(QINDX2)
      QINDX1=QINDX1+1
      QINDX2=QINDX2+1
      END DO
      CALL SEBACK(QLIFT,QSTART,QSTKLN,ISTK)
      CALL SENEXT
C
C===========
C POIN token
C===========
      ELSE IF (WD(1:4).EQ.'POIN') THEN
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      CALL SELTVF('POINt-vector=',R)
      XREF=R(1)
      YREF=R(2)
      ZREF=R(3)
      CALL SENEXT
      IF (WD(1:3).EQ.'CUT') THEN
      CALL NEXTF('CUT=',RMAX)
      CALL SENEXT
      ELSE
      RMAX=8.0D0
      END IF
      IF (.NOT.QCOOR) THEN
      CALL WRNDIE(-5,'SELRPN','no coordinates present')
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=.FALSE.
      END DO
      ELSE
      RMAX2=RMAX*RMAX
      DO IAT=1,NATOM
      DXL=XREF-XC(IAT)
      DYL=YREF-YC(IAT)
      DZL=ZREF-ZC(IAT)
      DIST2=DXL*DXL+DYL*DYL+DZL*DZL
      ELEMNT=(DIST2.LE.RMAX2)
      QSTACK(QSTART+IAT-1)=ELEMNT
      END DO
      END IF
C
C===========
C ATTR token
C===========
      ELSE IF (WD(1:4).EQ.'ATTR') THEN
      CALL NEXTQL('ATTRibute=')
      IF (WD(1:3).EQ.'ABS') THEN
      QABS=.TRUE.
      CALL NEXTQL('ATTRibute=')
      ELSE
      QABS=.FALSE.
      END IF
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      IF (WD(1:4).EQ.'X   ') THEN
      CALL SELATT(X,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'Y   ') THEN
      CALL SELATT(Y,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'Z   ') THEN
      CALL SELATT(Z,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'WMAI') THEN
      CALL SELATT(WMAIN,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'QMAI') THEN
      CALL SELATT(QMAIN,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'B   ') THEN
      CALL SELATT(WMAIN,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'Q   ') THEN
      CALL SELATT(QMAIN,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'XCOM') THEN
      CALL SELATT(XCOMP,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'YCOM') THEN
      CALL SELATT(YCOMP,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'ZCOM') THEN
      CALL SELATT(ZCOMP,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'WCOM') THEN
      CALL SELATT(WCOMP,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'BCOM') THEN
      CALL SELATT(WCOMP,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'QCOM') THEN
      CALL SELATT(QCOMP,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'REFX') THEN
      CALL SELATT(REFX,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'REFY') THEN
      CALL SELATT(REFY,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'REFZ') THEN
      CALL SELATT(REFZ,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'RMSD') THEN
      CALL SELATT(RMSD,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'HARM') THEN
      CALL SELATT(KCNSTR,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'ZHAR') THEN
      CALL SELATT(KZCNSTR,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'MASS') THEN
      CALL SELATT(AMASS,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'CHAR') THEN
      CALL SELATT(CG,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:6).EQ.'STORE1') THEN
      CALL SELATT(HEAP(PTRSTO(1)),QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:6).EQ.'STORE2') THEN
      CALL SELATT(HEAP(PTRSTO(2)),QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:6).EQ.'STORE3') THEN
      CALL SELATT(HEAP(PTRSTO(3)),QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:6).EQ.'STORE4') THEN
      CALL SELATT(HEAP(PTRSTO(4)),QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:6).EQ.'STORE5') THEN
      CALL SELATT(HEAP(PTRSTO(5)),QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:6).EQ.'STORE6') THEN
      CALL SELATT(HEAP(PTRSTO(6)),QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:6).EQ.'STORE7') THEN
      CALL SELATT(HEAP(PTRSTO(7)),QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:6).EQ.'STORE8') THEN
      CALL SELATT(HEAP(PTRSTO(8)),QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:6).EQ.'STORE9') THEN
      CALL SELATT(HEAP(PTRSTO(9)),QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'DX  ') THEN
      CALL SELATT(DX,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'DY  ') THEN
      CALL SELATT(DY,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'DZ  ') THEN
      CALL SELATT(DZ,QSTART,NATOM,QSTACK,QABS)
      ELSE IF (WD(1:4).EQ.'VX  ') THEN
      DO I=1,NATOM
      XV(I)=XV(I)/TIMFAC
      END DO
      CALL SELATT(XV,QSTART,NATOM,QSTACK,QABS)
      DO I=1,NATOM
      XV(I)=XV(I)*TIMFAC
      END DO
      ELSE IF (WD(1:4).EQ.'VY  ') THEN
      DO I=1,NATOM
      YV(I)=YV(I)/TIMFAC
      END DO
      CALL SELATT(YV,QSTART,NATOM,QSTACK,QABS)
      DO I=1,NATOM
      YV(I)=YV(I)*TIMFAC
      END DO
      ELSE IF (WD(1:4).EQ.'VZ  ') THEN
      DO I=1,NATOM
      ZV(I)=ZV(I)/TIMFAC
      END DO
      CALL SELATT(ZV,QSTART,NATOM,QSTACK,QABS)
      DO I=1,NATOM
      ZV(I)=ZV(I)*TIMFAC
      END DO
      ELSE IF (WD(1:4).EQ.'FBET') THEN
      CALL SELATT(FBETA,QSTART,NATOM,QSTACK,QABS)
C
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_A1') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('A1',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_A2') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('A2',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_A3') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('A3',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_A4') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('A4',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_B1') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('B1',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_B2') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('B2',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_B3') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('B3',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_B4') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('B4',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_C') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('C',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_FP') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('FP',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE IF (WD(1:WDLEN).EQ.'SCATTER_FDP') THEN
      XF=ALLHP(IREAL8(NATOM))
      CALL XSCAFIL('FDP',NATOM,HEAP(XF),XRERR)
      CALL SELATT(HEAP(XF),QSTART,NATOM,QSTACK,QABS)
      CALL FREHP(XF,IREAL8(NATOM))
C
      ELSE
      CALL DSPERR('SELRPN',' unknown attribute')
      END IF
      CALL SENEXT
C
C============
C KNOWn token
C============
      ELSE IF (WD(1:4).EQ.'INIT'.OR.WD(1:4).EQ.'KNOW') THEN
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      ELEMNT=INITIA(IAT,X,Y,Z)
      QSTACK(QSTART+IAT-1)=ELEMNT
      END DO
      CALL SENEXT
C
C===========
C HYDR token
C===========
      ELSE IF (WD(1:4).EQ.'HYDR') THEN
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO IAT=1,NATOM
      ELEMNT=HYDROG(IAT)
      QSTACK(QSTART+IAT-1)=ELEMNT
      END DO
      CALL SENEXT
C
C===========
C TAG token
C===========
      ELSE IF (WD(1:4).EQ.'TAG ') THEN
C
C selects exactly one atom from each "entity".  An
C entity is defined as a set of atoms with equal segids,
C resids, and residue names.
C
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      SEG1='    '
      RES1='    '
      RESN1='    '
      DO IAT=1,NATOM
      IF (SEG1.NE.SEGID(IAT).OR.RES1.NE.RESID(IAT).OR.
     &      RESN1.NE.RES(IAT)) THEN
      SEG1=SEGID(IAT)
      RES1=RESID(IAT)
      RESN1=RES(IAT)
C
C check if it has been previously already defined
      EXIST=.FALSE.
      DO IIAT=1,IAT-1
      IF (QSTACK(QSTART+IIAT-1)) THEN
      IF (SEG1.EQ.SEGID(IIAT).AND.RES1.EQ.RESID(IIAT).AND.
     &      RESN1.EQ.RES(IIAT)) THEN
      EXIST=.TRUE.
      END IF
      END IF
      END DO
C
C
      IF (.NOT.EXIST) THEN
      QSTACK(QSTART+IAT-1)=.TRUE.
      ELSE
      QSTACK(QSTART+IAT-1)=.FALSE.
      END IF
      ELSE
      QSTACK(QSTART+IAT-1)=.FALSE.
      END IF
      END DO
      CALL SENEXT
C=====================================
      ELSE IF (WD(1:4).EQ.'FBOX') THEN
      CALL NEXTF('xmin=',XMIN)
      CALL NEXTF('xmax=',XMAX)
      CALL NEXTF('ymin=',YMIN)
      CALL NEXTF('ymax=',YMAX)
      CALL NEXTF('zmin=',ZMIN)
      CALL NEXTF('zmax=',ZMAX)
      XMIN=MAX(ZERO,XMIN)
      XMAX=MIN(ONE,XMAX)
      YMIN=MAX(ZERO,YMIN)
      YMAX=MIN(ONE,YMAX)
      ZMIN=MAX(ZERO,ZMIN)
      ZMAX=MIN(ONE,ZMAX)
      CALL SENEXT
C
      IF (.NOT.QCOOR) THEN
      CALL WRNDIE(-5,'SELRPN','no coordinates present')
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=.FALSE.
      END DO
      ELSE
      DO IAT=1,NATOM
C
      QSTACK(QSTART+IAT-1)=.FALSE.
C compute fractional coordinates
      CALL XPSYMO(1,XC(IAT),YC(IAT),ZC(IAT),1,XRRFF,YRRFF,ZRRFF)
C
C apply mod operation
      XRRFF=MOD(XRRFF+BIG,ONE)
      YRRFF=MOD(YRRFF+BIG,ONE)
      ZRRFF=MOD(ZRRFF+BIG,ONE)
      IF          (XRRFF.GE.XMIN.AND.XRRFF.LE.XMAX
     &        .AND.YRRFF.GE.YMIN.AND.YRRFF.LE.YMAX
     &        .AND.ZRRFF.GE.ZMIN.AND.ZRRFF.LE.ZMAX) THEN
      QSTACK(QSTART+IAT-1)=.TRUE.
      END IF
      END DO
      END IF
C=====================================
      ELSE IF (WD(1:5).EQ.'SFBOX') THEN
      CALL NEXTF('xmin=',XMIN)
      CALL NEXTF('xmax=',XMAX)
      CALL NEXTF('ymin=',YMIN)
      CALL NEXTF('ymax=',YMAX)
      CALL NEXTF('zmin=',ZMIN)
      CALL NEXTF('zmax=',ZMAX)
      XMIN=MAX(ZERO,XMIN)
      XMAX=MIN(ONE,XMAX)
      YMIN=MAX(ZERO,YMIN)
      YMAX=MIN(ONE,YMAX)
      ZMIN=MAX(ZERO,ZMIN)
      ZMAX=MIN(ONE,ZMAX)
      CALL SENEXT
C
      IF (.NOT.QCOOR) THEN
      CALL WRNDIE(-5,'SELRPN','no coordinates present')
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=.FALSE.
      END DO
      ELSE
      DO IAT=1,NATOM
C
      QSTACK(QSTART+IAT-1)=.FALSE.
C
C loop over non-crystallographic symmetry operators
      DO INSYM=1,XNNSYM
C
C compute non-crystallographic symmetry related atom (orthogonal coords)
      CALL NCSXYZ(1,XC(IAT),YC(IAT),ZC(IAT),INSYM,XRRF,YRRF,ZRRF)
C
C loop over crystallographic symmetry operators
      DO ISYM=1,XRNSYM
C
C compute symmetry related molecule (in fractional coordinates)
      CALL XPSYMO(1,XRRF,YRRF,ZRRF,ISYM,XRRFF,YRRFF,ZRRFF)
C
C apply mod operation
      XRRFF=MOD(XRRFF+BIG,ONE)
      YRRFF=MOD(YRRFF+BIG,ONE)
      ZRRFF=MOD(ZRRFF+BIG,ONE)
      IF          (XRRFF.GE.XMIN.AND.XRRFF.LE.XMAX
     &        .AND.YRRFF.GE.YMIN.AND.YRRFF.LE.YMAX
     &        .AND.ZRRFF.GE.ZMIN.AND.ZRRFF.LE.ZMAX) THEN
      QSTACK(QSTART+IAT-1)=.TRUE.
      END IF
      END DO
      END DO
      END DO
      END IF
C==================================
      ELSE
      WRITE(6,'(2A)') ' %SELRPN-ERR: Unrecognized word ',WD(1:WDLEN)
      CALL WRNDIE(-5,'SELRPN','parsing error')
C exit selection subroutine
      GOTO 99999
      END IF
C-end-parse-token
C
77777 CONTINUE
C=====================================
      IF (WD(1:4).EQ.'AROU') THEN
C
      IF (.NOT.QCOOR) THEN
      CALL WRNDIE(-5,'SELRPN','no coordinates present')
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=.FALSE.
      END DO
      CALL SENEXT
      ELSE
      CALL NEXTF('AROUnd=',RMAX)
      RMAX2=RMAX*RMAX
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO INDX1=QSTART,QSTART+QSTKLN-1
      QSTACK(INDX1)=QSTACK(INDX1-QSTKLN)
      END DO
      DO INDX1=QSTART,QSTART+QSTKLN-1
      QSTACK(INDX1-QSTKLN)=.FALSE.
      END DO
C
      QINDX2=QSTART
      QKNOWN=.TRUE.
      DO IAT2=1,NATOM
      IF (QSTACK(QINDX2)) THEN
      IF (.NOT. INITIA(IAT2,X,Y,Z)) THEN
      QKNOWN=.FALSE.
      ELSE
      QINDX1=QSTART-QSTKLN-1
C
C================
C-VECTORIZED-CODE
C================
      DO IAT=1,NATOM
      DXL=XC(IAT)-XC(IAT2)
      DYL=YC(IAT)-YC(IAT2)
      DZL=ZC(IAT)-ZC(IAT2)
      DIST2=DXL*DXL+DYL*DYL+DZL*DZL
      QINDX1=QINDX1+1
      QSTACK(QINDX1)=QSTACK(QINDX1).OR.(DIST2.LE.RMAX2)
      END DO
      END IF
      END IF
      QINDX2=QINDX2+1
      END DO
      END IF
      CALL SEBACK(QLIFT,QSTART,QSTKLN,ISTK)
      IF (.NOT.QKNOWN) THEN
      WRITE(6,'(A)')
     & ' %SELRPN-AROUnd-ERR: some atom coordinates unknown'
      END IF
      CALL SENEXT
C=====================================
      ELSE IF (WD(1:4).EQ.'SARO') THEN
C
      IF (.NOT.QCOOR) THEN
      CALL WRNDIE(-5,'SELRPN','no coordinates present')
      DO IAT=1,NATOM
      QSTACK(QSTART+IAT-1)=.FALSE.
      END DO
      CALL SENEXT
      ELSE
      CALL NEXTF('SAROund=',RMAX)
      RMAX2=RMAX*RMAX
      CALL SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
      DO INDX1=QSTART,QSTART+QSTKLN-1
      QSTACK(INDX1)=QSTACK(INDX1-QSTKLN)
      END DO
      DO INDX1=QSTART,QSTART+QSTKLN-1
      QSTACK(INDX1-QSTKLN)=.FALSE.
      END DO
C
C allocate temporary space
      XF=ALLHP(IREAL8(NATOM))
      YF=ALLHP(IREAL8(NATOM))
      ZF=ALLHP(IREAL8(NATOM))
      XFF=ALLHP(IREAL8(NATOM))
      YFF=ALLHP(IREAL8(NATOM))
      ZFF=ALLHP(IREAL8(NATOM))
      XFFF=ALLHP(IREAL8(NATOM))
      YFFF=ALLHP(IREAL8(NATOM))
      ZFFF=ALLHP(IREAL8(NATOM))
C
C loop over non-crystallographic symmetry operators
      DO INSYM=1,XNNSYM
C
C compute non-crystallographic symmetry related molecule (orthogonal coords)
      CALL NCSXYZ(NATOM,X,Y,Z,INSYM,HEAP(XF),HEAP(YF),HEAP(ZF))
C
C loop over crystallographic symmetry operators
      DO ISYM=1,XRNSYM
C
C compute symmetry related molecule (in fractional coordinates)
      CALL XPSYMO(NATOM,HEAP(XF),HEAP(YF),HEAP(ZF),ISYM,
     &            HEAP(XFF),HEAP(YFF),HEAP(ZFF))
C
      QINDX2=QSTART
      QKNOWN=.TRUE.
      DO IAT2=1,NATOM
      IF (QSTACK(QINDX2)) THEN
      IF (.NOT. INITIA(IAT2,X,Y,Z)) THEN
      QKNOWN=.FALSE.
      ELSE
      QINDX1=QSTART-QSTKLN-1
C
C compute the minimum image distance
      CALL XPIMAG(X(IAT2),Y(IAT2),Z(IAT2),
     &            1,NATOM,HEAP(XFF),HEAP(YFF),HEAP(ZFF),
     &            HEAP(XFFF),HEAP(YFFF),HEAP(ZFFF),.FALSE.,0)
      CALL SELMMM(NATOM,HEAP(XFFF),HEAP(YFFF),HEAP(ZFFF),QSTACK,
     &            QINDX1,RMAX2)
      END IF
      END IF
      QINDX2=QINDX2+1
      END DO
C
      END DO
      END DO
C
C deallocate temporary space
      CALL FREHP(XF,IREAL8(NATOM))
      CALL FREHP(YF,IREAL8(NATOM))
      CALL FREHP(ZF,IREAL8(NATOM))
      CALL FREHP(XFF,IREAL8(NATOM))
      CALL FREHP(YFF,IREAL8(NATOM))
      CALL FREHP(ZFF,IREAL8(NATOM))
      CALL FREHP(XFFF,IREAL8(NATOM))
      CALL FREHP(YFFF,IREAL8(NATOM))
      CALL FREHP(ZFFF,IREAL8(NATOM))
C
      END IF
      CALL SEBACK(QLIFT,QSTART,QSTKLN,ISTK)
      IF (.NOT.QKNOWN) THEN
      WRITE(6,'(A)')
     & ' %SELRPN-SAROund-ERR: some atom coordinates unknown'
      END IF
C
      CALL SENEXT
C=====================================
      END IF
C
      CALL SEUP(LEVEL,RETMAX)
C     return to address:
      IF (RETADR(LEVEL).EQ.5) GOTO 3335
      IF (RETADR(LEVEL).EQ.6) GOTO 3336
      IF (RETADR(LEVEL).EQ.7) GOTO 3337
      IF (RETADR(LEVEL).EQ.8) GOTO 3338
      IF (RETADR(LEVEL).EQ.9) GOTO 3339
      IF (RETADR(LEVEL).EQ.10) GOTO 3340
      WRITE(6,'(A)') 
     &' %SELRPN-ERR: Return address unknown. Fatal coding error.'
      CALL DIE
C---- END PROCEDURE FACTOR ----------------------------------------------
C
C
C
C
9999  CONTINUE
      QINDX1=QSTART
      DO IAT=1,NATOM
      IF (QSTACK(QINDX1)) THEN
      FLAGS(IAT)=1
      ELSE
      FLAGS(IAT)=0
      END IF
      QINDX1=QINDX1+1
      END DO
C
C error exit label
99999 CONTINUE
      NSELCT=0
      DO IAT=1,NATOM
      IF (FLAGS(IAT).EQ.1) NSELCT=NSELCT+1
      END DO
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I7,A,I7)')
     & ' SELRPN:',NSELCT,' atoms have been selected out of',NATOM
C
      END IF
      DPVAL = NSELCT
      CALL DECLAR( 'SELECT', 'DP', ' ', DCVAL, DPVAL )
C
      RETURN
      END
C
      SUBROUTINE SEDOWN(LEVEL,RETMAX)
C
C go down one level in SELRPN
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER LEVEL, RETMAX
C begin
      IF (LEVEL.GE.RETMAX) THEN
      WRITE(6,'(A)') ' %SELRPN-ERR: Return address overflow (RETMAX)'
      ELSE
      LEVEL=LEVEL+1
      END IF
      RETURN
      END
C===================================================================
      SUBROUTINE SEUP(LEVEL,RETMAX)
C go up one level in SELRPN
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER LEVEL, RETMAX
C begin
      LEVEL=LEVEL-1
      IF (LEVEL.LT.0) THEN
      WRITE(6,'(A)') ' %SELRPN-ERR: Return address underflow. Check'
      END IF
      RETURN
      END
C===================================================================
      SUBROUTINE SENEXT
C
C get next word  ( lets be able to type HELP at any place !! )
C Author: Axel T. Brunger
C ======================
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C begin
      CALL NEXTWD('SELRPN>')
C
      DO WHILE (WD(1:4).EQ.'HELP')
C
      CALL CNSHELP('atom-selection')
C
      CALL NEXTWD('SELRPN>')
      END DO
      RETURN
      END
C=====================================================================
      SUBROUTINE SEPUSH(QLIFT,QSTART,QSTKLN,ISTK,STKSIZ)
C
C push stack in SELRPN
C Author: Axel T. Brunger
C ======================
      IMPLICIT NONE
      INTEGER QLIFT, QSTART, QSTKLN, ISTK, STKSIZ
C begin
      QSTART=QSTART+QSTKLN
      IF (QSTART.GT.STKSIZ) THEN
      CALL WRNDIE(-5,'SELRPN>',
     & 'max. depth of selection expression exceeded')
      ELSE
      QLIFT=QLIFT+1
      ISTK=QLIFT
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE SEBACK(QLIFT,QSTART,QSTKLN,ISTK)
C
C back stack in SELRPN
C Author: Axel T. Brunger
C ======================
      IMPLICIT NONE
      INTEGER QLIFT, QSTART, QSTKLN, ISTK
C begin
      QLIFT=QLIFT-1
      QSTART=QSTART-QSTKLN
      IF (QSTART.LE.0) THEN
      CALL WRNDIE(-5,'SELRPN>',
     & 'fatal stack underflow, check code')
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE SELATT(ARRAY,QSTART,NATOM,QSTACK,QABS)
C
C routine selects atoms based on comparisons with values in ARRAY
C Author: Axel T. Brunger
C ======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      DOUBLE PRECISION ARRAY(*)
      INTEGER QSTART, NATOM
      LOGICAL QSTACK(*), QABS
C local
      INTEGER IROPR, IAT
      DOUBLE PRECISION TOL, CUTOFF, DATAPT
      CHARACTER*2 RELAT
C begin
      IROPR=1
      CALL NEXTQL('operator=')
      RELAT=WD
C
C check if we have a composite relat. operator, such as ">=" or "<=" or "=="
      IF (WD(1:1).EQ.'<'.OR.WD(1:1).EQ.'>') THEN
      CALL NEXTQL('operator=')
      IF (WD(1:1).EQ.'=') THEN
      RELAT(2:2)='='
      ELSE
      CALL SAVEWD
      END IF
      END IF
C
      IF (RELAT.EQ.'LT'.OR.RELAT.EQ.'<') THEN
      IROPR=1
      ELSE IF (RELAT.EQ.'LE'.OR.RELAT.EQ.'<=') THEN
      IROPR=2
      ELSE IF (RELAT.EQ.'GT'.OR.RELAT.EQ.'>') THEN
      IROPR=3
      ELSE IF (RELAT.EQ.'GE'.OR.RELAT.EQ.'>=') THEN
      IROPR=4
      ELSE IF (RELAT.EQ.'EQ'.OR.RELAT.EQ.'=') THEN
      IROPR=5
      ELSE IF (RELAT.EQ.'NE'.OR.RELAT.EQ.'#') THEN
      IROPR=6
      ELSE
      CALL DSPERR('SELRPN','unknown relational operator')
      END IF
      CALL NEXTF('comparison=',CUTOFF)
      TOL=MAX(RSMALL,ABS(RSMALL*CUTOFF))
      DO IAT=1,NATOM
      DATAPT=ARRAY(IAT)
      IF (QABS.AND.DATAPT.LT.0) DATAPT=-DATAPT
      IF ((1).EQ.(IROPR)) THEN
      QSTACK(QSTART+IAT-1)= ((DATAPT-CUTOFF).LT.-TOL)
      ELSE IF ((2).EQ.(IROPR)) THEN
      QSTACK(QSTART+IAT-1)= ((DATAPT-CUTOFF).LE.TOL)
      ELSE IF ((3).EQ.(IROPR)) THEN
      QSTACK(QSTART+IAT-1)= ((DATAPT-CUTOFF).GT.TOL)
      ELSE IF ((4).EQ.(IROPR)) THEN
      QSTACK(QSTART+IAT-1)= ((DATAPT-CUTOFF).GE.-TOL)
      ELSE IF ((5).EQ.(IROPR)) THEN
      QSTACK(QSTART+IAT-1)= (ABS(DATAPT-CUTOFF).LT.TOL)
      ELSE IF ((6).EQ.(IROPR)) THEN
      QSTACK(QSTART+IAT-1)= (ABS(DATAPT-CUTOFF).GE.TOL)
      END IF
      END DO
      RETURN
      END
C
      SUBROUTINE SELRCL(QSTACK,STORE,NATOM)
C
C Subroutine fills QSTACK according to the values in STORE
C QSTACK=(STORE.NE.0)
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output

      INCLUDE 'consta.inc'
      LOGICAL QSTACK(*)
      DOUBLE PRECISION STORE(*)
      INTEGER NATOM
C local
      INTEGER I
C begin
      DO I=1,NATOM
      QSTACK(I)=(ABS(STORE(I)).GT.RSMALL)
      END DO
      RETURN
      END
C====================================================================
      SUBROUTINE SELTVF(PROMPT,VECTOR)
C
C parses a 3-dimensional vector within selection
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C =======================
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      DOUBLE PRECISION VECTOR(*)
      CHARACTER*(*) PROMPT
C local
      INTEGER I
      LOGICAL OK
      LOGICAL CLOOP
      CHARACTER*1 SBRA, SKET
      PARAMETER (SBRA='(',SKET=')')
C begin
      ERROR=.FALSE.
C
      CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A,3G12.5,A)') PROMPT,'(',(VECTOR(I),I=1,3),')'
      ELSE
      DO I=1,3
      VECTOR(I)=0.0D0
      END DO
C
      IF (WD(1:1).EQ.SBRA) THEN
      I=0
      CALL NEXTWD(PROMPT)
      CLOOP = .TRUE.
      DO WHILE (CLOOP)
C parse the x,y,z form of <vector>
      IF (I.GE.3) THEN
      ERROR=.TRUE.
      WRITE(6,'(A)') ' %SELTVF-ERR: dimension of vector greater 3'
      ELSE
      I=I+1
      VECTOR(I)=DECODF(WD,WDLEN,OK)
      ERROR=.NOT.OK
      IF (ERROR) THEN
      WRITE(6,'(2A)') ' %SELTVF-ERR: numerical value expected : ',
     1               WD(1:WDLEN)
      END IF
      END IF
      CALL NEXTWD(PROMPT)
      IF (WD(1:1).EQ.SKET.OR.ERROR) CLOOP = .FALSE.
      END DO
C
      ELSE
      WRITE(6,'(A)') ' %SELTVF-ERR: vector starts with a "("'
      END IF
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE SELMMM(NATOM,XFF,YFF,ZFF,QSTACK,QINDX1,RMAX2)
C
C Routine checks if distances are within specified cutoff RMAX2
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER NATOM
      DOUBLE PRECISION XFF(*), YFF(*), ZFF(*)
      LOGICAL QSTACK(*)
      INTEGER QINDX1
      DOUBLE PRECISION RMAX2
C local
      INTEGER IAT
      DOUBLE PRECISION DIST2
C begin
      DO IAT=1,NATOM
      DIST2=XFF(IAT)**2+YFF(IAT)**2+ZFF(IAT)**2
      QINDX1=QINDX1+1
      QSTACK(QINDX1)=QSTACK(QINDX1).OR.(DIST2.LE.RMAX2)
      END DO
      RETURN
      END
C=====================================================================
      SUBROUTINE XSCAFIL(OP,NATOM,ARRAY,XRERR)
C
C Routine stores scatter factor OP into atomic array ARRAY.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
C I/O
      CHARACTER*(*) OP
      INTEGER NATOM
      DOUBLE PRECISION ARRAY(*)
      LOGICAL XRERR
C local
      INTEGER INDEX
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
      IF (HPATOF.EQ.0.OR.XRNATF.EQ.0) THEN
      CALL DSPERR ( 'PUTSCA',
     & 'SCATTer database undefined.' )
      XRERR = .TRUE.
      DO INDEX=1,NATOM
      ARRAY(INDEX)=ZERO
      END DO
C
      ELSE
      CALL XSCAFIL2(OP,NATOM,ARRAY,XRERR,XRNATF,HEAP(HPATOF),
     &            HEAP(HPINDF),XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP)
      END IF
      RETURN
      END
C===========================================================
      SUBROUTINE XSCAFIL2(OP,NATOM,ARRAY,XRERR,XRNATF,XRATOF,
     &                  XRINDF,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP)
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'vector.inc'
C I/O
      CHARACTER*(*) OP
      INTEGER NATOM
      DOUBLE PRECISION ARRAY(*)
      LOGICAL XRERR
      INTEGER XRNATF, XRATOF(*), XRINDF(*)
      INTEGER XRSM
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(XRSM), XRFDP(XRSM)
C local
      INTEGER IAT, IINDEX, INDEX
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C loop over all atoms
      DO INDEX=1,NATOM
C we need to lookup INDEX
      IINDEX=0
      DO IAT=1,XRNATF
      IF (XRATOF(IAT).EQ.INDEX) IINDEX=IAT
      END DO
C
      IF ( IINDEX.EQ.0 ) THEN
      ARRAY(INDEX)=ZERO
      ELSEIF (OP.EQ.'A1') THEN
      ARRAY(INDEX)=XRSA(XRINDF(IINDEX),1)
      ELSEIF (OP.EQ.'A2') THEN
      ARRAY(INDEX)=XRSA(XRINDF(IINDEX),2)
      ELSEIF (OP.EQ.'A3') THEN
      ARRAY(INDEX)=XRSA(XRINDF(IINDEX),3)
      ELSEIF (OP.EQ.'A4') THEN
      ARRAY(INDEX)=XRSA(XRINDF(IINDEX),4)
      ELSEIF (OP.EQ.'B1') THEN
      ARRAY(INDEX)=XRSB(XRINDF(IINDEX),1)
      ELSEIF (OP.EQ.'B2') THEN
      ARRAY(INDEX)=XRSB(XRINDF(IINDEX),2)
      ELSEIF (OP.EQ.'B3') THEN
      ARRAY(INDEX)=XRSB(XRINDF(IINDEX),3)
      ELSEIF (OP.EQ.'B4') THEN
      ARRAY(INDEX)=XRSB(XRINDF(IINDEX),4)
      ELSEIF (OP.EQ.'C') THEN
      ARRAY(INDEX)=XRSC(XRINDF(IINDEX))
      ELSEIF (OP.EQ.'FP') THEN
      ARRAY(INDEX)=XRFP(XRINDF(IINDEX))
      ELSEIF (OP.EQ.'FDP') THEN
      ARRAY(INDEX)=XRFDP(XRINDF(IINDEX))
      END IF
C
      END DO
C
      RETURN
      END
C

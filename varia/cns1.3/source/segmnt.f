      SUBROUTINE SEGMNT
C
C Segment generation and residue sequence interpreter
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'cstack.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
C local
C parameter
C SEGDIM is the maximum number of atoms involved during a
C single patching process, NADD is the maximum number
C of added atoms during a single patching process.
C MAX is the maximum number of link statements in all link files combined. 
      INTEGER SEGDIM, NADD, MAX
      PARAMETER (SEGDIM=100, NADD=50, MAX=6000)
C pointer
      INTEGER SET, ISET
      INTEGER LINK, LINKH, LINKT, FIRST
      INTEGER FIRSTT, LAST, LASTH
      INTEGER SLINKH, SLINKT, SFIRST, SLASTH
C begin
      LINK=CALLST(ICHAR4(MAX))
      LINKH=CALLST(ICHAR4(MAX))
      LINKT=CALLST(ICHAR4(MAX))
      FIRST=CALLST(ICHAR4(MAX))
      FIRSTT=CALLST(ICHAR4(MAX))
      LAST=CALLST(ICHAR4(MAX))
      LASTH=CALLST(ICHAR4(MAX))
      SLINKH=CALLST(ICHAR4(MAX))
      SLINKT=CALLST(ICHAR4(MAX))
      SFIRST=CALLST(ICHAR4(MAX))
      SLASTH=CALLST(ICHAR4(MAX))
      SET=CALLST(ICHAR4(SEGDIM+NADD))
      ISET=ALLHP(INTEG4(SEGDIM+NADD))
      CALL SEGMN2(NADD,SEGDIM,CSTACK(SET),HEAP(ISET),MAX,
     &           CSTACK(LINK),CSTACK(LINKH),CSTACK(LINKT),
     &           CSTACK(FIRST),CSTACK(FIRSTT),CSTACK(LAST),
     &           CSTACK(LASTH),CSTACK(SLINKH),
     &           CSTACK(SLINKT),CSTACK(SFIRST),CSTACK(SLASTH))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(MAX))
      CALL CFREST(ICHAR4(SEGDIM+NADD))
      CALL FREHP(ISET,INTEG4(SEGDIM+NADD))
C
      RETURN
      END
C
      SUBROUTINE SEGMN2(NADD,DIM,SET,ISET,MAX,
     &                  LINK, LINKH, LINKT, FIRST,
     &                  FIRSTT, LAST, LASTH, 
     &                  SLINKH, SLINKT, SFIRST, SLASTH)
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'ctitla.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      INTEGER NADD, DIM
      CHARACTER*4 SET(*)
      INTEGER ISET(*)
      INTEGER MAX
      CHARACTER*4 LINK(MAX), LINKH(MAX), LINKT(MAX), FIRST(MAX)
      CHARACTER*4 FIRSTT(MAX), LAST(MAX), LASTH(MAX)
      CHARACTER*4 SLINKH(MAX), SLINKT(MAX), SFIRST(MAX), SLASTH(MAX)
C local
      INTEGER NUMBER, NRES, NRSOLD
      INTEGER I, J, DUMMY, NSET
      INTEGER ILINK, IFIRST, ILAST
      LOGICAL FOUND, QFIRST, ECHOLD, SIDCHK, QSEP, QCHAIN, QLINK
      CHARACTER*4 RID, SID, REN, RENOLD, RIDOLD, PATNAM
      CHARACTER*4 SIDX, SIDOLD, RIDX, RENX, IUPX
      CHARACTER*5 A
      CHARACTER*1 CHAIN
      DOUBLE PRECISION XIN, YIN, ZIN, WIN
      INTEGER NATOM1, IJ
      LOGICAL QRESTART
      CHARACTER*100 ERRMSG
      INTEGER ERRMSG_LEN, CURRSTRM
C begin
C
C defaults
      SID='    '
      NRES=0
      SIDCHK=.FALSE.
C
C parsing
      CALL PUSEND('SEGMENT>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('SEGMENT>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-segment')
C
      ELSE IF (WD(1:4).EQ.'NAME') THEN
      CALL NEXTA4('SEGMent-name=',SID)
      SIDCHK=.TRUE.
      ELSE IF (WD(1:4).EQ.'MOLE') THEN
      REN='    '
      NUMBER=1
      CALL PUSEND('MOLECULE>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('MOLECULE>')
      IF (WD(1:4).EQ.'NAME') THEN
      CALL NEXTA4('MOLEcule=',REN)
      ELSE IF (WD(1:4).EQ.'NUMB') THEN
      CALL NEXTI('NUMBer=',NUMBER)
      ELSE
      CALL CHKEND('MOLECULE>',DONE)
      END IF
      END DO
      DONE=.FALSE.
      DO I=1,NUMBER
      NRES=NRES+1
C
C convert residue number into hybrid_36 format ATB 2/15/10
      WDT=' '
      CALL HY36ENCODE(4,NRES,WDT,ERRMSG,ERRMSG_LEN)
      IF (ERRMSG_LEN.NE.0) THEN
      WRITE(6,'(A)') ' SEGMNT-err: specified',
     & ' number of molecules exceeds RESId field'
      WRITE(6,'(A,A)') ' SEGMNT-err', ERRMSG(1:ERRMSG_LEN) 
      CALL WRNDIE(-5,'SEGMNT','Too many molecules specified.')
      END IF
      WDTLEN=4
      CALL TRIML(WDT,WDTLEN)
CC
CCCC      CALL ENCODI(NRES,WDT,WDTMAX,WDTLEN)
      CALL COPYST(RID,4,DUMMY,WDT,WDTLEN)
C
C now check whether residue is new at all
      FOUND=.FALSE.
      J=0
      DO WHILE (J.LT.NATOM.AND..NOT.FOUND)
      J=J+1
      FOUND=(SID.EQ.SEGID(J).AND.RID.EQ.RESID(J).AND.REN.EQ.RES(J))
      END DO
      IF (FOUND) THEN
      WRITE(6,'(7A)') ' %SEGMNT-ERR: attempt to enter duplicate ',
     1'residue ',SID,' ',RID,' ',REN
      ELSE
      CALL GENRES(REN,RID,SID)
      END IF
      END DO
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(A,I5,6A)')
     & ' SEGMNT: ',NUMBER,' ',REN,' molecules have ',
     & 'been added to segment "',SID,'"'
      END IF
C
      ELSE IF (WD(1:4).EQ.'CHAI') THEN
C
C defaults
      ILINK=0
      IFIRST=0
      ILAST=0
      QSEP=.FALSE.
      QCHAIN=.FALSE.
      QLINK=.TRUE.
C
C parsing
      CALL PUSEND('CHAIN>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('CHAIN>')
      CALL MISCOM('CHAIN>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'LINK') THEN
      IF (ILINK.GE.MAX) THEN
      WRITE(6,'(A)') ' %SEGMENT-ERR: max. number of LINK exceeded'
      ELSE
      ILINK=ILINK+1
      END IF
      CALL NEXTA4('LINK-patch=',LINK(ILINK))
      LINKH(ILINK)='*   '
      LINKT(ILINK)='*   '
      SLINKH(ILINK)='-   '
      SLINKT(ILINK)='+   '
      CALL PUSEND('LINK>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('LINK>')
      IF (WD(1:4).EQ.'HEAD') THEN
      CALL NEXTA4('HEAD-reference=',SLINKH(ILINK))
      CALL NEXTA4('HEAD-residues=',LINKH(ILINK))
      ELSE IF (WD(1:4).EQ.'TAIL') THEN
      CALL NEXTA4('TAIL-reference=',SLINKT(ILINK))
      CALL NEXTA4('TAIL-residues=',LINKT(ILINK))
      ELSE
      CALL CHKEND('LINK>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      ELSE IF (WD(1:4).EQ.'FIRS') THEN
      IF (IFIRST.GE.MAX) THEN
      WRITE(6,'(A)') ' %SEGMENT-ERR: max. number of FIRST exceeded'
      ELSE
      IFIRST=IFIRST+1
      END IF
      CALL NEXTA4('FIRSt-patch=',FIRST(IFIRST))
      SFIRST(IFIRST)='+   '
      FIRSTT(IFIRST)='*   '
      CALL PUSEND('FIRST>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('FIRST>')
      IF (WD(1:4).EQ.'TAIL') THEN
      CALL NEXTA4('TAIL-reference=',SFIRST(IFIRST))
      CALL NEXTA4('TAIL-residues=',FIRSTT(IFIRST))
      ELSE
      CALL CHKEND('FIRST>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      ELSE IF (WD(1:4).EQ.'LAST') THEN
      IF (ILAST.GE.MAX) THEN
      WRITE(6,'(A)') ' %SEGMENT-ERR: max. number of LAST exceeded'
      ELSE
      ILAST=ILAST+1
      END IF
      CALL NEXTA4('LAST-patch=',LAST(ILAST))
      LASTH(ILAST)='*   '
      SLASTH(ILAST)='-   '
      CALL PUSEND('LAST>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('LAST>')
      IF (WD(1:4).EQ.'HEAD') THEN
      CALL NEXTA4('HEAD-reference=',SLASTH(ILAST))
      CALL NEXTA4('HEAD-residues=',LASTH(ILAST))
      ELSE
      CALL CHKEND('LAST>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      ELSE IF (WD(1:4).EQ.'SEQU') THEN
      QFIRST=.TRUE.
      NRSOLD=NRES
      CALL PUSEND('SEQUENCE>')
      CALL NEXTWD('SEQUENCE>')
      IF (WD(1:3).EQ.'END') CALL CHKEND('SEQUENCE>',DONE)
      DO WHILE (.NOT.DONE)
      CALL COPYST(REN,4,DUMMY,WD,WDLEN)
      NRES=NRES+1
C
C convert residue number into hybrid_36 format ATB 2/15/10
      WDT=' '
      CALL HY36ENCODE(4,NRES,WDT,ERRMSG,ERRMSG_LEN)
      IF (ERRMSG_LEN.NE.0) THEN
      WRITE(6,'(A)') ' SEGMNT-err: specified',
     & ' number of residues exceeds RESId field'
      WRITE(6,'(A,A)') ' SEGMNT-err', ERRMSG(1:ERRMSG_LEN) 
      CALL WRNDIE(-5,'SEGMNT','Too many residues specified.')
      END IF
      WDTLEN=4
      CALL TRIML(WDT,WDTLEN)
CC
CCCC      CALL ENCODI(NRES,WDT,WDTMAX,WDTLEN)
      CALL COPYST(RID,4,DUMMY,WDT,WDTLEN)
C
C now check whether residue is new at all
      FOUND=.FALSE.
      J=0
      DO WHILE (J.LT.NATOM.AND..NOT.FOUND)
      J=J+1
      FOUND=(SID.EQ.SEGID(J).AND.RID.EQ.RESID(J).AND.REN.EQ.RES(J))
      END DO
      IF (FOUND) THEN
      WRITE(6,'(7A)') ' %SEGMNT-ERR: attempt to enter duplicate ',
     1'residue ',SID,' ',RID,' ',REN
      NRES=NRES-1
      ELSE
      IF (QFIRST) THEN
      CALL SEGFIR(REN,RID,SID,FIRST,FIRSTT,SFIRST,PATNAM,SET,
     &            NATOM1,IFIRST,NADD,ISET,NSET,DIM)
      QFIRST=.FALSE.
      ELSE
      CALL SEGLNK(NATOM1,ILINK,NSET,DIM,ISET,NADD,REN,RID,
     &      SID,LINK,LINKH,LINKT,SLINKH,SLINKT,RENOLD,RIDOLD,
     &      SET,PATNAM,.TRUE.)
      END IF
      RENOLD=REN
      RIDOLD=RID
      END IF
      CALL NEXTWD('SEQUence-element (terminate with END) =')
      IF (WD(1:3).EQ.'END') CALL CHKEND('SEQUENCE>',DONE)
      END DO
      DONE=.FALSE.
      IF (NRES.GT.NRSOLD) THEN
      CALL SEGLAS(NSET,NADD,DIM,ISET,ILAST,NATOM1,REN,RID,LAST,
     &           LASTH,SLASTH,SET,PATNAM)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I5,3A)') ' SEGMNT: ',NRES-NRSOLD,
     1 ' residues were inserted into segment "',SID,'"'
      END IF
      END IF
C
      ELSE IF (WD(1:4).EQ.'SEPA') THEN
      CALL NEXTLO('SEPArate-by-segid=',QSEP)
C
      ELSE IF (WD(1:4).EQ.'CONV') THEN
      CALL NEXTLO('CONVert-chainid-to-segid=',QCHAIN)
C
      ELSE IF (WD(1:4).EQ.'COOR') THEN
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' SEGMNT-info: sequence read from coordinate file'
      END IF
C
C initialize various counters
      NRSOLD=NRES
      SIDOLD='    '
      QFIRST=.TRUE.
      RID='    '
      REN='    '
      ECHOLD=.FALSE.
C 
      CURRSTRM=0
C
C parsing
      CALL PUSEND('COOR>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('COOR>')
CCCC
C check for missing end statement, ATB 2/21/10
      IF (CURRSTRM.EQ.0) THEN
      CURRSTRM=NSTRM
      END IF
      IF (CURRSTRM.NE.NSTRM.AND.WD(1:WDLEN).NE.'END') THEN
C
      WRITE(6,'(A)') 
     & ' SEGMNT-err: End of coordinate file without END statement.'
      WRITE(6,'(A)') 
     & '             Please insert END statement and re-run.'
      CALL WRNDIE(-5,'SEGMNT','Fatal error in reading coordinates.')
      END IF
CCCC
      IF (ECHOLD) QECHO=.TRUE.
      CALL MISCOM('COOR>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'ATOM'.OR.WD(1:4).EQ.'HETA') THEN
C
C this section reads the Brookhaven atom coordinates
      READ(COMLYN,'(12X,A4,1X,A4,A1,A5,3X,3F8.3,6X,F6.2,6X,A4)',
     &     ERR=8888) IUPX,RENX,CHAIN,A,XIN,YIN,ZIN,WIN,SIDX
C
C special handling of PDB RESSEQ and ICODE, ATB 02/15/10
      IF (A(5:5).EQ.' ') then
         RIDX=A(1:4)
      ELSE
         IF (A(1:1).NE.' ') then
           WRITE(6,'(A,A,A)')     
     & ' %CREAD-ERR: residue ID and insertion character ',
     & A,' exceed 4 characters.'           
           CALL WRNDIE(-5,'SEGMNT',
     &      'Unsupported PDB fields.')
         END IF
         RIDX=A(2:5) 
      END IF
C
C make RIDX left-justified
      IJ=4
      CALL TRIML(RIDX,IJ)
      CALL TRIMM(RIDX,IJ)
CCC      IJ=5
CCC      CALL TRIML(A,IJ)
CCC      CALL TRIMM(A,IJ)
CCC      IF (IJ.GT.4) THEN
CCC      WRITE(6,'(A,A,A)')
CCC     & ' %CREAD-ERR: residue ID and insertion character ',
CCC     & A,' exceed 4 characters.'
CCC      END IF
CCC      RIDX=A(1:4)
C
C make RENX left-justified, ATB 12/12/09
      I=4
      CALL TRIML(RENX,I)
C
      IF ( QCHAIN ) THEN
         IF (CHAIN.NE.' ') THEN
           SIDX=CHAIN
         END IF
      END IF
C
      IF (QSEP) THEN
      IF (SIDX.NE.SIDOLD) THEN
C
C make last patch if more than one residue in current chain
      IF (NRES.GT.NRSOLD) THEN
      CALL SEGLAS(NSET,NADD,DIM,ISET,ILAST,NATOM1,REN,RID,LAST,
     &           LASTH,SLASTH,SET,PATNAM)
      IF (WRNLEV.GE.5.AND.NRES-NRSOLD.GT.1) THEN
      WRITE(6,'(A,I5,3A)') 
     & ' SEGMNT-info: auto chain termination due to different segid.',
     & NRES-NRSOLD,
     & ' residues were inserted into segid "',SIDOLD,'"'
      END IF
      END IF
C
C initialize various counters
      NRSOLD=NRES
      QFIRST=.TRUE.
      RID='    '
      REN='    '
      ECHOLD=.FALSE.
      SIDOLD=SIDX
C
C modification, ATB 12/01/09: automatic chain termination if no link found
C to next residue
      ELSEIF (RID.NE.RIDX.OR.REN.NE.RENX) THEN
C test if a link exists between these two residues.  If not -> restart
      FOUND=.FALSE.
      I=0
      DO WHILE (.NOT. (FOUND.OR.I.GE.ILINK))
      I=I+1
      FOUND=(RENX.EQ.LINKT(I).AND.RENOLD.EQ.LINKH(I))
      END DO
      IF (.NOT.FOUND) THEN
C
C make last patch if more than one residue in current chain
      IF (NRES.GT.NRSOLD) THEN
      CALL SEGLAS(NSET,NADD,DIM,ISET,ILAST,NATOM1,REN,RID,LAST,
     &           LASTH,SLASTH,SET,PATNAM)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I5,3A)') 
     & ' SEGMNT-info: auto chain termination due to unavailable link. ',
     &  NRES-NRSOLD,
     & ' residues were inserted into segid "',SIDOLD,'"'
      END IF
      END IF
C
C initialize various counters
      NRSOLD=NRES
      QFIRST=.TRUE.
      RID='    '
      REN='    '
      ECHOLD=.FALSE.
      SIDOLD=SIDX
      END IF   
      END IF
      END IF
C
      CURSOR=COMLEN
      CALL SEGATM(NATOM1,ILINK,NSET,DIM,ISET,NADD,NRES,
     &           IFIRST,REN,RID,SID,LINK,LINKH,LINKT,SLINKH,SLINKT,
     &           FIRST,SFIRST,FIRSTT,RENOLD,RIDOLD,SET,PATNAM,RIDX,
     &           SIDX,RENX,QFIRST,SIDCHK,QLINK)
      IF (.NOT.QLINK) QLINK=.TRUE.
C
      ECHOLD=QECHO
      QECHO=.FALSE.
C
      ELSE IF (WD(1:3).EQ.'TER') THEN
C
C terminate chain and start new chain
C
C go to end of line
      CURSOR=COMLEN
C
C make last patch
      IF (NRES.GT.NRSOLD) THEN
      CALL SEGLAS(NSET,NADD,DIM,ISET,ILAST,NATOM1,REN,RID,LAST,
     &           LASTH,SLASTH,SET,PATNAM)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I5,3A)') 
     & ' SEGMNT-info: chain termination due to TER keyword. ',
     & NRES-NRSOLD,
     & ' residues were inserted into segid "',SID,'"'
      END IF
      END IF
C
C initialize various counters
      NRSOLD=NRES
      QFIRST=.TRUE.
      RID='    '
      REN='    '
      ECHOLD=.FALSE.
C
      ELSE IF (WD(1:5).EQ.'BREAK') THEN
C
C do not make link, but do not start new chain
C
C go to end of line
      CURSOR=COMLEN
C
C initialize various flags
      QFIRST=.FALSE.
      QLINK=.FALSE.
C
      ELSE IF (WD(1:5).EQ.'CRYST') THEN
C
C go to end of line
      CURSOR=COMLEN
      ELSE IF (WD(1:6).EQ.'ANISOU') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'ORIG') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'SCALE') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'HEADER') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'COMPND') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'SOURCE') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'KEYWDS') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'EXPDTA') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'DBREF') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'AUTHOR') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'REVDAT') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'JRNL') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'SEQADV') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'SEQRES') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'TITLE') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'MODRES') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'FTNOTE') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'HET   ') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'HETNAM') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'HETSYN') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'CAVEAT') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'FORMUL') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'HELIX') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'SHEET') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'TURN') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'SSBOND') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'SPLIT') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'SITE') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:4).EQ.'LINK') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:5).EQ.'MTRIX') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'CISPEP') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'CONECT') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'MASTER') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE IF (WD(1:6).EQ.'SPRSDE') THEN
C
C go to end of line
      CURSOR=COMLEN
C
      ELSE
C
      CALL CHKEND('COOR>',DONE)
      END IF
      END IF
      END DO
      GOTO 7777
8888  ERROR=.TRUE.
7777  CONTINUE
      IF (ERROR) WRITE(6,'(A)') ' %SEGMNT-ERR: ERROR during read'
      DONE=.FALSE.
C
C make last patch
      IF (NRES.GT.NRSOLD) THEN
      CALL SEGLAS(NSET,NADD,DIM,ISET,ILAST,NATOM1,REN,RID,LAST,
     &           LASTH,SLASTH,SET,PATNAM)
      IF (WRNLEV.GE.5) THEN     
      WRITE(6,'(A,I5,3A)') 
     & ' SEGMNT-info: chain termination due to END keyword. ',
     & NRES-NRSOLD,
     & ' residues were inserted into segid "',SID,'"'
      END IF
      END IF
C
      ELSE
      CALL CHKEND('CHAIN>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      ELSE
      CALL CHKEND('SEGMENT>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
C finally, clear the lists
      CALL SCRATC
      CALL SHOW
      RETURN
      END
C
      SUBROUTINE SEGFIR(REN,RID,SID,FIRST,FIRSTT,SFIRST,PATNAM,SET,
     &            NATOM1,IFIRST,NADD,ISET,NSET,DIM)
C
C generates residue REN and apply first patch (see SEGMNT above)
C Author: Axel Brunger
C I/O
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      CHARACTER*4 REN, RID, SID, FIRST(*), FIRSTT(*), SFIRST(*)
      CHARACTER*4 PATNAM, SET(*)
      INTEGER NATOM1, IFIRST, NADD, ISET(*), NSET, DIM
C local
      INTEGER I, J
      LOGICAL FOUND
C begin
      NATOM1=NATOM
      CALL GENRES(REN,RID,SID)
C
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(8A)') ' SEGMNT: residue ',REN,' ',RID,' has been ',
     & 'added to segid "',SID,'"'
      END IF
C
      IF (IFIRST.GT.0) THEN
C
C lookup patch FIRST residue
      FOUND=.FALSE.
      I=0
      DO WHILE (.NOT. (FOUND.OR.I.GE.IFIRST))
      I=I+1
      FOUND=(REN.EQ.FIRSTT(I))
      END DO
C disabled wildcard matching, ATB 12/01/09
CCC      IF (.NOT.FOUND) THEN
CCC      I=0
CCC      DO WHILE (.NOT. (FOUND.OR.I.GE.IFIRST))
CCC      I=I+1
CCC      CALL EQSTWC(REN,4,FIRSTT(I),4,1,1,FOUND)
CCC      END DO
CCC      END IF
      IF (.NOT.FOUND) THEN
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' %Segmnt-info: No matching FIRST patch found'
      END IF
      ELSE
      PATNAM=FIRST(I)
C
C generate SET, ISET for PATCH
      NSET=0
      DO J=NATOM1+1,NATOM
      IF (NSET.GE.DIM) THEN
      WRITE(6,'(A)') ' %SEGMNT-ERR: DIM exceeded'
      ELSE
      NSET=NSET+1
      SET(NSET)=SFIRST(I)
      ISET(NSET)=J
      END IF
      END DO
      CALL PATCH(PATNAM,NATOM+NADD,NSET,SET,ISET,.FALSE.)
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(6A)') ' SEGMNT: ',PATNAM,' has been patched to ',
     & REN,' ',RID
      END IF
      END IF
      END IF
      RETURN
      END
C
      SUBROUTINE SEGLNK(NATOM1,ILINK,NSET,DIM,ISET,NADD,REN,RID,
     &      SID,LINK,LINKH,LINKT,SLINKH,SLINKT,RENOLD,RIDOLD,
     &      SET,PATNAM,QLINK)
C
C generate residue ren and link it to chain head
C (see SEGMNT above)
C Author: Axel Brunger
C I/O
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INTEGER NATOM1, ILINK, NSET, DIM, ISET(*), NADD
      CHARACTER*4 REN, RID, SID, LINK(*), LINKH(*), LINKT(*)
      CHARACTER*4 SLINKH(*), SLINKT(*)
      CHARACTER*4 RENOLD, RIDOLD, SET(*), PATNAM
      LOGICAL QLINK
C local
      INTEGER I, J, NATOM2
      LOGICAL FOUND, COND
C begin
C
C
C NATOM1+1 keeps track of the first atom in the previously inserted 
C residue (NATOM1 is the end of next to last residue) -> store in NATOM2
      NATOM2=NATOM1
C
C NATOM1 is now the last atom in the previously inserted residue
      NATOM1=NATOM
C
C generate new residue REN
      CALL GENRES(REN,RID,SID)
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(9A)') ' SEGMNT: residue ',REN,' ',RID,' has been ',
     & 'added to segid "',SID,'"'
      END IF
C
      IF (ILINK.GT.0.AND.QLINK) THEN
C
C lookup patch link element between last and current residue
      FOUND=.FALSE.
      I=0
      DO WHILE (.NOT. (FOUND.OR.I.GE.ILINK))
      I=I+1
      FOUND=(REN.EQ.LINKT(I).AND.RENOLD.EQ.LINKH(I))
      END DO
CCC disabled wildcard matching, ATB 12/01/92
CCC      IF (.NOT.FOUND) THEN
CCC      I=0
CCC      DO WHILE (.NOT. (FOUND.OR.I.GE.ILINK))
CCC      I=I+1
CCC      CALL EQSTWC(REN,4,LINKT(I),4,1,1,COND)
CCC      CALL EQSTWC(RENOLD,4,LINKH(I),4,1,1,FOUND)
CCC      FOUND=FOUND.AND.COND
CCC      END DO
CCC      END IF
      IF (.NOT. FOUND) THEN
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' %Segmnt-info: no matching LINK found'
      END IF
      ELSE
      PATNAM=LINK(I)
C
C generate SET, ISET list for PATCH
      NSET=0
      DO J=NATOM2+1,NATOM1
      IF (NSET.GE.DIM) THEN
      WRITE(6,'(A)') ' %SEGMNT-ERR: DIM exceeded'
      ELSE
      NSET=NSET+1
      SET(NSET)=SLINKH(I)
      ISET(NSET)=J
      END IF
      END DO
      DO J=NATOM1+1,NATOM
      IF (NSET.GE.DIM) THEN
      WRITE(6,'(A)') ' %SEGMNT-ERR: DIM exceeded'
      ELSE
      NSET=NSET+1
      SET(NSET)=SLINKT(I)
      ISET(NSET)=J
      END IF
      END DO
      CALL PATCH(PATNAM,NATOM+NADD,NSET,SET,ISET,.FALSE.)
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(10A)') ' SEGMNT: ',PATNAM,' has been patched to ',
     & RENOLD,' ',RIDOLD,' and ',REN,' ',RID
      END IF
      END IF
      END IF
C
C Redetermine NATOM1 (end of next-to-last residue) in case it
C was changed by patch.
      NATOM1=NATOM-1
      DO WHILE( RESID(NATOM1).EQ.RESID(NATOM) .AND.
     &          SEGID(NATOM1).EQ.SEGID(NATOM) .AND. NATOM1.GT.1)
      NATOM1=NATOM1-1
      END DO     
      RETURN
      END
C
      SUBROUTINE SEGLAS(NSET,NADD,DIM,ISET,ILAST,NATOM1,REN,RID,LAST,
     &           LASTH,SLASTH,SET,PATNAM)
C
C lookup appropriate last patch for residue REN
C Author: Axel Brunger
C I/O
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INTEGER NSET, NADD, DIM, ISET(*), ILAST, NATOM1
      CHARACTER*4 REN, RID, LAST(*), LASTH(*), SLASTH(*), SET(*)
      CHARACTER*4 PATNAM
C local
      LOGICAL FOUND
      INTEGER I, J, L
C begin
      IF (ILAST.GT.0) THEN
      FOUND=.FALSE.
      I=0
      DO WHILE (.NOT. (FOUND.OR.I.GE.ILAST))
      I=I+1
      FOUND=(REN.EQ.LASTH(I))
      END DO
C disabled wildcard matching, ATB, 12/01/09
CCC      IF (.NOT.FOUND) THEN
CCC      I=0
CCC      DO WHILE (.NOT. (FOUND.OR.I.GE.ILAST))
CCC      I=I+1
CCC      CALL EQSTWC(REN,4,LASTH(I),4,1,1,FOUND)
CCC      END DO
CCC      END IF
      IF (.NOT.FOUND) THEN
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' %Segmnt-info: no matching LAST patch found'
      END IF
      ELSE
      PATNAM=LAST(I)
C
C generate SET, ISET list for PATCH
      NSET=0
      L=NATOM
      DO J=NATOM1+1,NATOM
      IF (NSET.GE.DIM) THEN
      WRITE(6,'(A)') ' %SEGMNT-ERR: DIM exceeded'
      ELSE
      NSET=NSET+1
      SET(NSET)=SLASTH(I)
      ISET(NSET)=J
      END IF
      END DO
      CALL PATCH(PATNAM,NATOM+NADD,NSET,SET,ISET,.FALSE.)
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(6A)') ' SEGMNT: ',PATNAM,' has been patched to ',
     & REN,' ',RID
      END IF
      END IF
      END IF
      RETURN
      END
C
      SUBROUTINE SEGATM(NATOM1,ILINK,NSET,DIM,ISET,NADD,NRES,
     &           IFIRST,REN,RID,SID,LINK,LINKH,LINKT,SLINKH,SLINKT,
     &           FIRST,SFIRST,FIRSTT,RENOLD,RIDOLD,SET,PATNAM,RIDX,
     &           SIDX,RENX,QFIRST,SIDCHK,QLINK)
C
C processes current atom line for COOR option
C (see routine SEGMNT above)
C Author: Axel Brunger
C  I/O
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INTEGER NATOM1, ILINK, NSET, DIM, ISET(*), NADD, NRES
      INTEGER IFIRST
      CHARACTER*4 REN, RID, SID, LINK(*), LINKH(*), LINKT(*)
      CHARACTER*4 SLINKH(*), SLINKT(*), FIRST(*), SFIRST(*)
      CHARACTER*4 FIRSTT(*), RENOLD, RIDOLD, SET(*), PATNAM
      CHARACTER*4 RIDX, SIDX, RENX
      LOGICAL QFIRST, SIDCHK, QLINK
C local
      INTEGER I
      LOGICAL FOUND
C begin
C
C now we check whether segment id's match and whether residue is new
      IF (SIDCHK) THEN
      IF (SID.NE.SIDX) THEN
      WRITE(6,'(A)')
     & ' %SEGMNT-ERR: segid in coordinate file does not match NAME.'
      CALL WRNDIE(-5,'SEGMNT',' check NAME specification.')
      END IF
      END IF
C
      IF (RID.NE.RIDX.OR.REN.NE.RENX) THEN
C
C now check whether residue is new at all
      FOUND=.FALSE.
      I=0
      DO WHILE (I.LT.NATOM.AND..NOT.FOUND)
      I=I+1
      FOUND=(SIDX.EQ.SEGID(I).AND.RIDX.EQ.RESID(I).AND.RENX.EQ.RES(I))
      END DO
      IF (FOUND) THEN
      WRITE(6,'(7A)') ' %SEGMNT-ERR: attempt to enter duplicate ',
     &'residue ',SIDX,' ',RIDX,' ',RENX
      CALL WRNDIE(-5,'SEGMNT',
     & ' check PDB file for residue duplications.')
      ELSE
      SID=SIDX
      RID=RIDX
      REN=RENX
      NRES=NRES+1
      IF (QFIRST) THEN
      CALL SEGFIR(REN,RID,SID,FIRST,FIRSTT,SFIRST,PATNAM,SET,
     &            NATOM1,IFIRST,NADD,ISET,NSET,DIM)
      QFIRST=.FALSE.
      ELSE
      CALL SEGLNK(NATOM1,ILINK,NSET,DIM,ISET,NADD,REN,RID,
     &            SID,LINK,LINKH,LINKT,SLINKH,SLINKT,RENOLD,RIDOLD,
     &            SET,PATNAM,QLINK)
      END IF
      RIDOLD=RID
      RENOLD=REN
      END IF
      END IF
      RETURN
      END

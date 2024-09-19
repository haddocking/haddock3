      SUBROUTINE GENRES(RESF,RESIDF,SEGIDF)
C
C adds a residue RESF to the molecular structure database
C with residue (segment) identifiers RESIDF (SEGIDF).
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'rtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'ener.inc'
      CHARACTER*4 RESF, RESIDF, SEGIDF
C local
      INTEGER I, IRES, J, IGRLST, NGRPFS, NATOMF
      INTEGER II, JJ, KK, LL, NL
      INTEGER NATPT, IPT
      LOGICAL FOUND
      INTEGER MARK
      PARAMETER (MARK=-9999)
C begin
C
C find the index of RESF in the topology file
      J=0
      FOUND=.FALSE.
      DO WHILE (J.LT.NRTRS.AND..NOT.FOUND)
      J=J+1
      FOUND=(AA(J).EQ.RESF.AND.RTRTYP(J).EQ.0)
      END DO
      IF (.NOT.FOUND) THEN
      WRITE(6,'(4A)') ' %GENRES-ERR: residue ',RESF,
     & ' not found in topology files. Add topology',
     & ' definition.'
      CALL WRNDIE(-5,'GENRES',
     & 'missing topology definition for residue')
      ELSE
      IRES=J
C
C insert atoms and atom specifications
      IGRLST=0
      NATOMF=NATOM+1
      NGRPFS=NGRP+1
      NATPT=0
      IF (IRES.GT.1) NATPT=NIC(1,IRES-1)
      DO I=NATPT+1,NIC(1,IRES)
      IF (NATOM.GE.MAXA) THEN
      CALL WRNDIE(-5,'GENRES',
     & 'exceeded MAXA parameter --> recompile program')
      ELSE
      NATOM=NATOM+1
      END IF
C
C temporary setup of IBLO ( so we don't have to search for the
C atom number I when setting up the exclusion list below )
      IBLO(NATOM)=I
C
C get chemical type, mass, iupac, ...
      IAC(NATOM)=MAC(I)
      AMASS(NATOM)=ARMASS(I)
      TYPE(NATOM)=FTP(I)(1:4)
      SEGID(NATOM)=SEGIDF
      RESID(NATOM)=RESIDF
      RES(NATOM)=RESF
      CG(NATOM)=CHG(I)
      LOOKUP(NATOM)=0
      QLOOKU(NATOM)=.TRUE.
      IMOVE(NATOM)=0
      CALL ATMINI(NATOM,NATOM)
C
C group list
      IF (IGRLST.NE.GRPR(I).OR.NGRP.EQ.0) THEN
      IGRLST=GRPR(I)
      IF (NGRP.GE.MAXGRP) THEN
      CALL WRNDIE(-5,'GENRES',
     & 'exceeded MAXGRP parameter --> recompile program')
      ELSE
      NGRP=NGRP+1
      END IF
      IGPBS(NGRP)=NATOM-1
      END IF
      END DO
      IGPBS(NGRP+1)=NATOM
C
C put in bonds
      IF (IRES.EQ.1) THEN
      I=1
      ELSE
      I=NIC(2,IRES-1)+1
      END IF
      DO IPT=I,NIC(2,IRES)
      II=SRCHC4(MIB(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      JJ=SRCHC4(MJB(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      IF (II.NE.MARK.AND.JJ.NE.MARK) THEN
      IF (NBOND.GE.MAXB) THEN
      CALL WRNDIE(-5,'GENRES',
     & 'exceeded MAXB parameter --> recompile program')
      ELSE
      NBOND=NBOND+1
      END IF
      IB(NBOND)=II
      JB(NBOND)=JJ
      QACBC(NBOND)=.TRUE.
      QACBB(NBOND)=.TRUE.
      ELSE
      WRITE(6,'(6A)') ' %GENRES-ERR: atom(s) in bond ',
     & MIB(IPT),' ',
     & MJB(IPT),' are undefined in topology residue ',RESF
      END IF
      END DO
C
C put in angles
      IF (IRES.EQ.1) THEN
      I=1
      ELSE
      I=NIC(3,IRES-1)+1
      END IF
      DO IPT=I,NIC(3,IRES)
      II=SRCHC4(MIT(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      JJ=SRCHC4(MJT(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      KK=SRCHC4(MKT(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      IF (II.NE.MARK.AND.JJ.NE.MARK.AND.KK.NE.MARK) THEN
      IF (NTHETA.GE.MAXT) THEN
      CALL WRNDIE(-5,'GENRES',
     & 'exceeded MAXT parameter --> recompile program')
      ELSE
      NTHETA=NTHETA+1
      END IF
      IT(NTHETA)=II
      JT(NTHETA)=JJ
      KT(NTHETA)=KK
      QACTC(NTHETA)=.TRUE.
      QACTB(NTHETA)=.TRUE.
      QACTUC(NTHETA)=.TRUE.
      QACTUB(NTHETA)=.TRUE.
      ELSE
      WRITE(6,'(8A)') ' %GENRES-ERR: atom(s) in angle ',
     & MIT(IPT),' ',MJT(IPT),' ',MKT(IPT),
     & ' are undefined in topology residue ',RESF
      END IF
      END DO
C
C put in dihedrals
      IF (IRES.EQ.1) THEN
      I=1
      ELSE
      I=NIC(4,IRES-1)+1
      END IF
      DO IPT=I,NIC(4,IRES)
      II=SRCHC4(MIP(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      JJ=SRCHC4(MJP(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      KK=SRCHC4(MKP(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      LL=SRCHC4(MLP(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      IF (II.NE.MARK.AND.JJ.NE.MARK.AND.KK.NE.MARK.AND.LL.NE.MARK) THEN
      IF (NPHI.GE.MAXP) THEN
      CALL WRNDIE(-5,'GENRES',
     & 'exceeded MAXP parameter --> recompile program')
      ELSE
      NPHI=NPHI+1
      END IF
      IP(NPHI)=II
      JP(NPHI)=JJ
      KP(NPHI)=KK
      LP(NPHI)=LL
      QACPB(NPHI)=.TRUE.
      QACPC(NPHI)=.TRUE.
      QACPD(NPHI)=.TRUE.
      ELSE
      WRITE(6,'(11A)') ' %GENRES-ERR: atom(s) in dihedral ',
     & MIP(IPT),' ',MJP(IPT),' ',MKP(IPT),' ',MLP(IPT),
     & ' are undefined in ','topology residue ',RESF
      END IF
      END DO
C
C put in improper dihedrals
      IF (IRES.EQ.1) THEN
      I=1
      ELSE
      I=NIC(5,IRES-1)+1
      END IF
      DO IPT=I,NIC(5,IRES)
      II=SRCHC4(MIM(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      JJ=SRCHC4(MJM(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      KK=SRCHC4(MKM(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      LL=SRCHC4(MLM(IPT)(1:4),TYPE,NATOMF,NATOM,MARK)
      IF (II.NE.MARK.AND.JJ.NE.MARK.AND.KK.NE.MARK.AND.LL.NE.MARK) THEN
      IF (NIMPHI.GE.MAXIMP) THEN
      CALL WRNDIE(-5,'GENRES',
     & 'exceeded MAXIMP parameter --> recompile program')
      ELSE
      NIMPHI=NIMPHI+1
      END IF
      IM(NIMPHI)=II
      JM(NIMPHI)=JJ
      KM(NIMPHI)=KK
      LM(NIMPHI)=LL
      QACIB(NIMPHI)=.TRUE.
      QACIC(NIMPHI)=.TRUE.
      QACID(NIMPHI)=.TRUE.
      ELSE
      WRITE(6,'(11A)') ' %GENRES-ERR: atom(s) in improper ',
     & MIM(IPT),' ',MJM(IPT),' ',MKM(IPT),' ',MLM(IPT),
     & ' are undefined in ','topology residue ',RESF
      END IF
      END DO
C
C put in nonbonding exclusions
      DO I=NATOMF,NATOM
      NL=IBLO(I)
      II=0
      IF (NL.GT.1) II=MXN(NL-1)
      JJ=MXN(NL)-II
      DO KK=1,JJ
      II=II+1
      LL=SRCHC4(MNB(II)(1:4),TYPE,NATOMF,NATOM,MARK)
      IF (LL.NE.MARK) THEN
      IF (NNB.GE.MAXNB) THEN
      CALL WRNDIE(-5,'GENRES',
     & 'exceeded MAXNB parameter --> recompile program')
      ELSE
      NNB=NNB+1
      END IF
      INB(NNB)=LL
      ELSE
      WRITE(6,'(4A)') ' %GENRES-ERR: atom in exclusion ',
     & MNB(II),' is undefined in topology residue ',RESF
      END IF
      END DO
      IBLO(I)=NNB
      END DO
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE PATSET
C
C Patch command parser
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'cstack.inc'
      INCLUDE 'heap.inc'
C local
      INTEGER ISLCT, SET, NADD, DIM
      PARAMETER (NADD=50)
C begin
C
C estimate the final number of atoms
C NADD is the maximum number of atoms to be added when invoking
C this routine once.
      DIM=MIN(MAXA,NATOM+NADD)
C
C dynamic
      ISLCT=ALLHP(INTEG4(DIM))
      SET=CALLST(ICHAR4(DIM))
      CALL PATSE2(DIM,HEAP(ISLCT),CSTACK(SET))
      CALL FREHP(ISLCT,INTEG4(DIM))
      CALL CFREST(ICHAR4(DIM))
      RETURN
      END
C
      SUBROUTINE PATSE2(DIM,ISLCT,SET)
C
      IMPLICIT NONE
C input/ouput
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INTEGER DIM, ISLCT(*)
      CHARACTER*4 SET(*)
C local
      CHARACTER*4 PATNAM, SETNAM, STEMP
      LOGICAL QSORT
      INTEGER I, NSLCT, NSET
C begin
      QSORT=.FALSE.
      DO I=1,NATOM
      SET(I)='    '
      END DO
      CALL NEXTA4('patch-residue=',PATNAM)
      CALL PUSEND('PATCH>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('PATCH>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-patch')
C
      ELSE IF (WD(1:4).EQ.'REFE') THEN
      CALL NEXTA4('reference=',SETNAM)
      CALL FILL4(ISLCT,NATOM,0)
      CALL SELCTA(ISLCT,NSLCT,X,Y,Z,.TRUE.)
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1) THEN
      IF (SET(I).NE.'    ') ERROR=.TRUE.
      SET(I)=SETNAM
      END IF
      END DO
      IF (ERROR) THEN
      WRITE(6,'(A)') ' %PATSET-ERR: referenced sets have to be disjoint'
      END IF
      ELSE IF (WD(1:4).EQ.'SORT') THEN
      CALL NEXTLO('SORT=',QSORT)
      ELSE
      CALL CHKEND('PATCH>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
C translate SET
      NSET=0
      DO I=1,NATOM
      IF (SET(I).NE.'    ') THEN
      NSET=NSET+1
      STEMP=SET(I)
      SET(NSET)=STEMP
      ISLCT(NSET)=I
      END IF
      END DO
C
      CALL PATCH(PATNAM,DIM,NSET,SET,ISLCT,QSORT)
C
C finally, clear the lists
      CALL SCRATC
      CALL SHOW
      RETURN
      END
C
      SUBROUTINE PATCH(PATNAM,DIM,NSET,SET,ISET,LSORT)
C
C General patching routine, patching data read from RTF.
C Deletions and addition of atoms, bonds, etc,
C modify atom specifications and group boundaries.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'rtf.inc'
      CHARACTER*4 PATNAM
      INTEGER DIM, NSET
      CHARACTER*4 SET(*)
      INTEGER ISET(*)
      LOGICAL LSORT
C local
      INTEGER GRPLST, MAP, INVMAP, DATOM, LINB, LIBLO, ATMLST
      INTEGER XNNB
C dynamic
      XNNB=NNB
      GRPLST=ALLHP(INTEG4(DIM))
      ATMLST=ALLHP(INTEG4(DIM))
      MAP=ALLHP(INTEG4(DIM))
      INVMAP=ALLHP(INTEG4(DIM))
      DATOM=ALLHP(INTEG4(DIM))
      LINB=ALLHP(INTEG4(XNNB))
      LIBLO=ALLHP(INTEG4(DIM))
      CALL PATCH2(PATNAM,DIM,NSET,SET,ISET,LSORT,
     &            HEAP(ATMLST),HEAP(GRPLST),
     &            HEAP(MAP),HEAP(INVMAP),HEAP(DATOM),
     &            HEAP(LINB),HEAP(LIBLO))
      CALL FREHP(LIBLO,INTEG4(DIM))
      CALL FREHP(LINB,INTEG4(XNNB))
      CALL FREHP(DATOM,INTEG4(DIM))
      CALL FREHP(INVMAP,INTEG4(DIM))
      CALL FREHP(MAP,INTEG4(DIM))
      CALL FREHP(ATMLST,INTEG4(DIM))
      CALL FREHP(GRPLST,INTEG4(DIM))
      RETURN
      END
C======================================================================
      SUBROUTINE PATCH2(PATNAM,DIM,NSET,SET,ISET,LSORT,
     &                 GRPLST,ATMLST,MAP,INVMAP,DATOM,LINB,
     &                 LIBLO)
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'rtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'cstack.inc'
      CHARACTER*4 PATNAM
      INTEGER DIM, NSET
      CHARACTER*4 SET(*)
      INTEGER ISET(*)
      LOGICAL LSORT
      INTEGER GRPLST(*), ATMLST(*)
      INTEGER MAP(*), INVMAP(*), DATOM(*)
      INTEGER LINB(*), LIBLO(*)
C local
      EXTERNAL  EXCH
      INTEGER   OFFSET, PATGRP
      INTEGER   GROUP, I, II, IPATC, IPT, K, INDEX
      INTEGER   J, JJ, KK, L, LL, NL, NDATOM, IDELTA, ALEN
      LOGICAL   FOUND, ERROR, BYPASS, QERROR
      CHARACTER*4 STEMP
      CHARACTER*8 ACTION
C parameter
      INTEGER   MARK, GPOFF
      PARAMETER (MARK=-9999, GPOFF=100)
      DOUBLE PRECISION ANUM
      PARAMETER (ANUM=9998.0D0)
C
C begin
      QERROR=.FALSE.
C
C get index of patch residue in topology file
      J=0
      FOUND=.FALSE.
      DO WHILE (J.LT.NRTRS.AND..NOT.FOUND)
      J=J+1
      FOUND=(AA(J).EQ.PATNAM.AND.RTRTYP(J).EQ.1)
      END DO
      IF (.NOT.FOUND) THEN
      WRITE(6,'(3A)')
     &' %PATCH-ERR: patch ',PATNAM,' not found in topology file '
      ELSE
      IPATC=J
C
      NDATOM=0
C
C Fill group list GRPLST with values corresponding to
C current molecular topology and store current order of atoms in ATMLST
C note: we assume that a single patch does not contain more
C than GPOFF sub-groups
      DO I=1,NGRP
      DO J=IGPBS(I)+1,IGPBS(I+1)
      GRPLST(J)=I*GPOFF
      ATMLST(J)=J
      END DO
      END DO
C
C Read in atoms of patch residue PATNAM
C Prepare atom mapping (add atoms, modify GRPLST, make list
C of atoms to be deleted)
      PATGRP=-1
      BYPASS=.TRUE.
C
      IF (IPATC.GT.1) THEN
      II=NIC(1,IPATC-1)+1
      ELSE
      II=1
      END IF
      DO I=II,NIC(1,IPATC)
      NL=PATOM(FTP(I),NSET,SET,ISET,TYPE,MARK)
      IF (DELAT(I).NE.'ADD ' .AND. NL.EQ.MARK) THEN
      CALL PATCHERR(QERROR,PATNAM,NSET,SET,ISET,DELAT(I),ACTION,ALEN)
      WRITE(6,'(5(1X,A))') '%PATCH-ERR: to be',ACTION(1:ALEN),'atom',
     &    FTP(I),'not found in molecular structure.'
      END IF
      IF (DELAT(I).EQ.'DELE') THEN
C
C Atom is to be deleted
      IF (NL.NE.MARK) THEN
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(7A)')
     &' PATCH: atom ',SEGID(NL),' ',RESID(NL),' ',TYPE(NL),
     &' will be deleted'
      END IF
      BYPASS=.FALSE.
      NDATOM=NDATOM+1
      DATOM(NDATOM)=NL
      END IF
C
C
      ELSE IF (DELAT(I).EQ.'MODI') THEN
C
C change atom specifications according to patch residue.
      IF (NL.NE.MARK) THEN
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(7A)')
     &' PATCH: atom ',SEGID(NL),' ',RESID(NL),' ',TYPE(NL),
     &' specifications will be modified'
      END IF
C
C dont change chemical type or charge if not specified
      IF (.NOT.MAC(I).EQ.'    ') THEN
      LOOKUP(NL)=0
      QLOOKU(NL)=.TRUE.
      IAC(NL)=MAC(I)
      AMASS(NL)=ARMASS(I)
      END IF
      IF (ABS(CHG(I)).LT.ANUM) CG(NL)=CHG(I)
C
C The following group re-assignment is bypassed if atom I in PRES
C has group number 0
      IF (GRPR(I).GT.0) THEN
      ATMLST(NL)=I
      BYPASS=.FALSE.
      IF (GRPR(I).NE.PATGRP) THEN
      PATGRP=GRPR(I)
      GROUP=PATGRP+GRPLST(NL)
      IF (PATGRP.GE.GPOFF) THEN
      CALL WRNDIE(-5,'PATCH',
     & 'exceeded GPOFF (routine PATCH2) parameter')
      END IF
      END IF
      GRPLST(NL)=GROUP
      END IF
      END IF
C
      ELSE IF (DELAT(I).EQ.'ADD ') THEN
C
C This atom should be added.
      IF (NATOM+1.GT.DIM) THEN
      CALL WRNDIE(-5,'PATCH',
     & 'exceeded NADD (routine PATCH) or MAXA parameter')
      ELSE
C
C check that atom won't be duplicated
      IF (NL.NE.MARK) THEN
      CALL PATCHERR(QERROR,PATNAM,NSET,SET,ISET,DELAT(I),ACTION,ALEN)
      WRITE(6,'(2A)') ' %PATCH-ERR: attempt to add duplicate atom =',
     1 FTP(I)
      CALL PATCHERRATM(FTP(I),NL)
      ELSE
C
C find appropriate SEGID, RESID for atom to be added
      FOUND=.FALSE.
      J=0
      DO WHILE (.NOT. (FOUND.OR.J.GE.NSET))
      J=J+1
      IF (SET(J).EQ.'NIL ') THEN
      FOUND=.TRUE.
      ELSE
      FOUND=(SET(J)(1:1).EQ.FTP(I)(1:1))
      END IF
      END DO
      IF (.NOT. FOUND) THEN
      CALL PATCHERR(QERROR,PATNAM,NSET,SET,ISET,DELAT(I),ACTION,ALEN)
      WRITE(6,'(3A)') ' %PATCH-ERR: no refe=',FTP(I)(1:1),
     &   ' set found for atom to be added'
      ELSE
C
C set up new atom
      BYPASS=.FALSE.
      NATOM=NATOM+1
      IAC(NATOM)=MAC(I)
      AMASS(NATOM)=ARMASS(I)
      CG(NATOM)=CHG(I)
      IMOVE(NATOM)=0
      LOOKUP(NATOM)=0
      QLOOKU(NATOM)=.TRUE.
      STEMP=RES(ISET(J))
      RES(NATOM)=STEMP
      STEMP=SEGID(ISET(J))
      SEGID(NATOM)=STEMP
      STEMP=RESID(ISET(J))
      RESID(NATOM)=STEMP
      IF (SET(J).EQ.'NIL ') THEN
      TYPE(NATOM)=FTP(I)(1:4)
      ELSE
      TYPE(NATOM)=FTP(I)(2:5)
      END IF
      NSET=NSET+1
      STEMP=SET(J)
      SET(NSET)=STEMP
      ISET(NSET)=NATOM
C
C first, all explicit nonbonded are set to zero
      IBLO(NATOM)=NNB
C
C now, the new atom has to be to assigned to a group number
C the case with GRPR(I)=0 is special
      IF (GRPR(I).EQ.0) THEN
      GRPLST(NATOM)=GRPLST(ISET(J))
      ATMLST(NATOM)=ATMLST(ISET(J))
      ELSE
      ATMLST(NATOM)=I
      IF (GRPR(I).NE.PATGRP) THEN
      PATGRP=GRPR(I)
      GROUP=PATGRP+GRPLST(ISET(J))
      IF (PATGRP.GE.GPOFF) THEN
      CALL WRNDIE(-5,'PATCH',
     & 'exceeded GPOFF (routine PATCH2) parameter')
      END IF
      END IF
      GRPLST(NATOM)=GROUP
      END IF
C
C Now reset coordinates, etc. for added atom:
      CALL ATMINI(NATOM,NATOM)
C
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(7A)') ' PATCH: atom ',SEGID(NATOM),' ',
     1 RESID(NATOM),' ',TYPE(NATOM),' added'
      END IF
      END IF
      END IF
      END IF
      END IF
      END DO
C
      IF (.NOT.BYPASS) THEN
      CALL SORT(NDATOM,EXCH,ORDER,DATOM,1,0,0,0,0,0,0)
C
C Now define a map as well as an inverse map to map the atom
C numbers.
C First, sort GRPLST, using quicksort (SORT4P conserves the order of
C equal elements, this is essential at this point!)
      CALL SORTP(NATOM,INVMAP,ORDER5,GRPLST,ATMLST,0,0,0,0,0,2)
      DO I=1,NATOM
      MAP(INVMAP(I))=I
      END DO
C
C account for deletions, reduce NATOM by NDATOM,
C deleted atoms are mapped to MARK
      OFFSET=1
      DO I=1,NATOM
      IF (OFFSET.LE.NDATOM) THEN
      IF (DATOM(OFFSET).EQ.I) THEN
      MAP(I)=MARK
      OFFSET=OFFSET+1
      END IF
      END IF
      END DO
      NATOM=NATOM-NDATOM
      OFFSET=0
      DO I=1,NATOM
      DO WHILE (MAP(INVMAP(I+OFFSET)).EQ.MARK)
      OFFSET=OFFSET+1
      END DO
      INVMAP(I)=INVMAP(I+OFFSET)
      END DO
      DO I=1,NATOM
      MAP(INVMAP(I))=I
      END DO
C
C
C now begin complete mapping of molecular topology, COORD, COORDC, CNST
C note that NATOM is already reduced by the number of deleted atoms.
      CALL MAPIC(MAP,INVMAP,GRPLST,LINB,LIBLO,MARK)
C
C finally map and compress the ISET array
      CALL AVALUE(MAP,ISET,NSET,MARK)
      II=NSET
      NSET=0
      DO I=1,II
      IF (ISET(I).NE.MARK) THEN
      NSET=NSET+1
      STEMP=SET(I)
      SET(NSET)=STEMP
      ISET(NSET)=ISET(I)
      END IF
      END DO
      END IF
C
C
C add in or delete internal coordinates using final atom sequence
C
C deletion of bonds, angles, dihedrals and improper dihedrals
C takes account of the symmetry operations:
C bonds:     (a,b)<-->(b,a)
C angles:    (a,b,c)<-->(c,b,a)
C dihedrals: (a,b,c,d)<-->(d,c,b,a)
C impropers: (a,b,c,d)<-->(d,c,b,a)
C
C put in bonds
C
      IF (IPATC.GT.1) THEN
      I=NIC(2,IPATC-1)+1
      ELSE
      I=1
      END IF
      DO IPT=I,NIC(2,IPATC)
      II=PATOM(MIB(IPT),NSET,SET,ISET,TYPE,MARK)
      JJ=PATOM(MJB(IPT),NSET,SET,ISET,TYPE,MARK)
      ERROR=.FALSE.
      IF (II.EQ.MARK.OR.JJ.EQ.MARK) THEN
      ERROR=.TRUE.
      ELSE IF (DELBD(IPT).EQ.'ADD ') THEN
      IF (NBOND.GE.MAXB) THEN
      CALL WRNDIE(-5,'PATCH',
     & 'exceeded MAXB parameter --> recompile program')
      ELSE
      NBOND=NBOND+1
      IB(NBOND)=II
      JB(NBOND)=JJ
      QACBB(NBOND)=.TRUE.
      QACBC(NBOND)=.TRUE.
      END IF
      ELSE IF (DELBD(IPT).EQ.'DELE') THEN
      INDEX=FIND52(IB,JB,0,0,0,II,JJ,0,0,0,NBOND,2,MARK)
      IF (INDEX.EQ.MARK) THEN
      INDEX=FIND52(IB,JB,0,0,0,JJ,II,0,0,0,NBOND,2,MARK)
      END IF
      IF (INDEX.NE.MARK) THEN
      IB(INDEX)=MARK
      ELSE
      ERROR=.TRUE.
      END IF
      END IF
      IF (ERROR) THEN
      CALL PATCHERR(QERROR,PATNAM,NSET,SET,ISET,DELBD(IPT),ACTION,ALEN)
      WRITE(6,'(6(1X,A))') '%PATCH-ERR: to be',ACTION(1:ALEN),'bond',
     &  MIB(IPT),MJB(IPT),'not found in molecular structure.'
      CALL PATCHERRATM(MIB(IPT),II)
      CALL PATCHERRATM(MJB(IPT),JJ)
      END IF
      END DO
C
C put in angles
      IF (IPATC.GT.1) THEN
      I=NIC(3,IPATC-1)+1
      ELSE
      I=1
      END IF
      DO IPT=I,NIC(3,IPATC)
      II=PATOM(MIT(IPT),NSET,SET,ISET,TYPE,MARK)
      JJ=PATOM(MJT(IPT),NSET,SET,ISET,TYPE,MARK)
      KK=PATOM(MKT(IPT),NSET,SET,ISET,TYPE,MARK)
      ERROR=.FALSE.
      IF (II.EQ.MARK.OR.JJ.EQ.MARK.OR.KK.EQ.MARK) THEN
      ERROR=.TRUE.
      ELSE IF (DELAN(IPT).EQ.'ADD ') THEN
      IF (NTHETA.GE.MAXT) THEN
      CALL WRNDIE(-5,'PATCH',
     & 'exceeded MAXT parameter --> recompile program')
      ELSE
      NTHETA=NTHETA+1
      IT(NTHETA)=II
      JT(NTHETA)=JJ
      KT(NTHETA)=KK
      QACTB(NTHETA)=.TRUE.
      QACTC(NTHETA)=.TRUE.
      QACTUB(NTHETA)=.TRUE.
      QACTUC(NTHETA)=.TRUE.
      END IF
      ELSE IF (DELAN(IPT).EQ.'DELE') THEN
      INDEX=FIND52(IT,JT,KT,0,0,II,JJ,KK,0,0,NTHETA,3,MARK)
      IF (INDEX.EQ.MARK) THEN
      INDEX=FIND52(IT,JT,KT,0,0,KK,JJ,II,0,0,NTHETA,3,MARK)
      END IF
      IF (INDEX.NE.MARK) THEN
      IT(INDEX)=MARK
      ELSE
      ERROR=.TRUE.
      END IF
      END IF
      IF (ERROR) THEN
      CALL PATCHERR(QERROR,PATNAM,NSET,SET,ISET,DELAN(IPT),ACTION,ALEN)
      WRITE(6,'(7(1X,A))') '%PATCH-ERR: to be',ACTION(1:ALEN),'angle',
     & MIT(IPT),MJT(IPT),MKT(IPT),'not found in molecular structure.'
      CALL PATCHERRATM(MIT(IPT),II)
      CALL PATCHERRATM(MJT(IPT),JJ)
      CALL PATCHERRATM(MKT(IPT),KK)
      END IF
      END DO
C
C     put in dihedrals
C
      IF (IPATC.GT.1) THEN
      I=NIC(4,IPATC-1)+1
      ELSE
      I=1
      END IF
      DO IPT=I,NIC(4,IPATC)
      ERROR=.FALSE.
      II=PATOM(MIP(IPT),NSET,SET,ISET,TYPE,MARK)
      JJ=PATOM(MJP(IPT),NSET,SET,ISET,TYPE,MARK)
      KK=PATOM(MKP(IPT),NSET,SET,ISET,TYPE,MARK)
      LL=PATOM(MLP(IPT),NSET,SET,ISET,TYPE,MARK)
      IF (II.EQ.MARK.OR.JJ.EQ.MARK.OR.KK.EQ.MARK.OR.LL.EQ.MARK) THEN
      ERROR=.TRUE.
      ELSE IF (DELPT(IPT).EQ.'ADD ') THEN
      IF (NPHI.GE.MAXP) THEN
      CALL WRNDIE(-5,'PATCH',
     & 'exceeded MAXP parameter --> recompile program')
      ELSE
      NPHI=NPHI+1
      END IF
      IP(NPHI)=II
      JP(NPHI)=JJ
      KP(NPHI)=KK
      LP(NPHI)=LL
      QACPB(NPHI)=.TRUE.
      QACPC(NPHI)=.TRUE.
      QACPD(NPHI)=.TRUE.
      ELSE IF (DELPT(IPT).EQ.'DELE') THEN
      INDEX=FIND52(IP,JP,KP,LP,0,II,JJ,KK,LL,0,NPHI,4,MARK)
      IF (INDEX.EQ.MARK) THEN
      INDEX=FIND52(IP,JP,KP,LP,0,LL,KK,JJ,II,0,NPHI,4,MARK)
      END IF
      IF (INDEX.NE.MARK) THEN
      IP(INDEX)=MARK
      ELSE
      ERROR=.TRUE.
      END IF
      END IF
      IF (ERROR) THEN
      CALL PATCHERR(QERROR,PATNAM,NSET,SET,ISET,DELPT(IPT),ACTION,ALEN)
      WRITE(6,'(8(1X,A))')'%PATCH-ERR: to be',ACTION(1:ALEN),'dihedral',
     & MIP(IPT),MJP(IPT),MKP(IPT),MLP(IPT),
     & 'not found in molecular structure.'
      CALL PATCHERRATM(MIP(IPT),II)
      CALL PATCHERRATM(MJP(IPT),JJ)
      CALL PATCHERRATM(MKP(IPT),KK)
      CALL PATCHERRATM(MLP(IPT),LL)
      END IF
      END DO
C
C put in improper dihedrals
      IF (IPATC.GT.1) THEN
      I=NIC(5,IPATC-1)+1
      ELSE
      I=1
      END IF
      DO IPT=I,NIC(5,IPATC)
      ERROR=.FALSE.
      II=PATOM(MIM(IPT),NSET,SET,ISET,TYPE,MARK)
      JJ=PATOM(MJM(IPT),NSET,SET,ISET,TYPE,MARK)
      KK=PATOM(MKM(IPT),NSET,SET,ISET,TYPE,MARK)
      LL=PATOM(MLM(IPT),NSET,SET,ISET,TYPE,MARK)
      IF (II.EQ.MARK.OR.JJ.EQ.MARK.OR.KK.EQ.MARK.OR.LL.EQ.MARK) THEN
      ERROR=.TRUE.
      ELSE IF (DELMT(IPT).EQ.'ADD ') THEN
      IF (NIMPHI.GE.MAXIMP) THEN
      CALL WRNDIE(-5,'PATCH',
     & 'exceeded MAXIMP parameter --> recompile program')
      ELSE
      NIMPHI=NIMPHI+1
      END IF
      IM(NIMPHI)=II
      JM(NIMPHI)=JJ
      KM(NIMPHI)=KK
      LM(NIMPHI)=LL
      QACIB(NIMPHI)=.TRUE.
      QACIC(NIMPHI)=.TRUE.
      QACID(NIMPHI)=.TRUE.
      ELSE IF (DELMT(IPT).EQ.'DELE') THEN
      INDEX=FIND52(IM,JM,KM,LM,0,II,JJ,KK,LL,0,NIMPHI,4,MARK)
      IF (INDEX.EQ.MARK) THEN
      INDEX=FIND52(IM,JM,KM,LM,0,LL,KK,JJ,II,0,NIMPHI,4,MARK)
      END IF
      IF (INDEX.NE.MARK) THEN
      IM(INDEX)=MARK
      ELSE
      ERROR=.TRUE.
      END IF
      ELSE
      ERROR=.TRUE.
      END IF
      IF (ERROR) THEN
      CALL PATCHERR(QERROR,PATNAM,NSET,SET,ISET,DELMT(IPT),ACTION,ALEN)
      WRITE(6,'(8(1X,A))')'%PATCH-ERR: to be',ACTION(1:ALEN),'improper',
     & MIM(IPT),MJM(IPT),MKM(IPT),MLM(IPT),
     & 'not found in molecular structure.'
      CALL PATCHERRATM(MIM(IPT),II)
      CALL PATCHERRATM(MJM(IPT),JJ)
      CALL PATCHERRATM(MKM(IPT),KK)
      CALL PATCHERRATM(MLM(IPT),LL)
      END IF
      END DO
C
C put in explicit nonbonding exclusions
      IF (IPATC.GT.1) THEN
      II=NIC(1,IPATC-1)
      ELSE
      II=0
      END IF
      DO J=II+1,NIC(1,IPATC)
C
      IF (J.GT.1) THEN
      L=MXN(J-1)
      ELSE
      L=0
      END IF
      IF (MXN(J).GT.L) THEN
C
C for this atom there are some exclusions ...
      JJ=PATOM(FTP(J),NSET,SET,ISET,TYPE,MARK)
      IF (JJ.EQ.MARK) THEN
      CALL PATCHERR(QERROR,PATNAM,NSET,SET,ISET,'DELE',ACTION,ALEN)
      WRITE(6,'(3A)')
     & ' %PATCH-ERR: to be deleted exclusion list atom ',FTP(J),
     & ' not found in molecular structure.'
      ELSE
C
C first, delete all previous exclusions for this atom
      IF (JJ.GT.1) THEN
      IDELTA=(IBLO(JJ)-IBLO(JJ-1))
      ELSE
      IDELTA=IBLO(JJ)
      END IF
      IF (IDELTA.GT.0) THEN
      NNB=NNB-IDELTA
      DO IPT=JJ,NATOM
      IBLO(IPT)=IBLO(IPT)-IDELTA
      END DO
      DO IPT=IBLO(JJ)+1,NNB
      INB(IPT)=INB(IPT+IDELTA)
      END DO
      END IF
C
C now, look-up all new pairs
      DO K=L+1,MXN(J)
      KK=PATOM(MNB(K),NSET,SET,ISET,TYPE,MARK)
      IF (KK.EQ.MARK) THEN
      CALL PATCHERR(QERROR,PATNAM,NSET,SET,ISET,'DELE',ACTION,ALEN)
      WRITE(6,'(3A)')
     & ' %PATCH-ERR: to be deleted exclusion list atom ',MNB(K),
     & ' not found in molecular structure.'
      ELSE
      IF (NNB.GE.MAXNB) THEN
      CALL WRNDIE(-5,'PATCH',
     & 'exceeded MAXNB parameter --> recompile program')
      ELSE
      NNB=NNB+1
      END IF
C
C finally, we have to insert the pair JJ, KK into the exclusion list
      DO IPT=JJ,NATOM
      IBLO(IPT)=IBLO(IPT)+1
      END DO
      DO IPT=NNB,IBLO(JJ)+1,-1
      INB(IPT)=INB(IPT-1)
      END DO
      INB(IBLO(JJ))=KK
      END IF
      END DO
      END IF
      END IF
C
      END DO
      IF (QERROR) WRITE(6,'(1X,78("%"),/)')
C
C remove all internal coordinates containing MARK
      CALL CMPRIC(MARK,LSORT)
      END IF
C
      RETURN
C
      END
C======================================================================
      SUBROUTINE PATCHERR(QERROR,PATNAM,NSET,SET,ISET,ID,ACTION,ALEN)
C
C Print a patch error header, if this is the first error.
C Also set the string ACTION(1:ALEN) to the action based on 
C the 4-character ID action identifier.
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      LOGICAL QERROR
      INTEGER NSET,ISET(NSET), ALEN
      CHARACTER*(*) PATNAM,SET(NSET),ID,ACTION
C local
      INTEGER I, J, K, IAT
C begin
      IF (ID.EQ.'DELE') THEN
      ACTION='deleted'
      ALEN=7
      ELSE IF (ID.EQ.'ADD ') THEN
      ACTION='added'
      ALEN=5
      ELSE
      ACTION='modified'
      ALEN=8
      END IF
      IF (.NOT.QERROR) THEN
      QERROR=.TRUE.
      WRITE(6,'(1X,78("%"))')
      WRITE(6,'(1X,2A)') '%PATCH-ERR: PATCH=',PATNAM
      I=0
      DO WHILE(I.LT.NSET)
      I=I+1
      J=I+1
      DO WHILE(SET(I).EQ.SET(J).AND.J.LT.NSET)
      J=J+1
      END DO
C     Write first and last atoms slected for each reference set,
C     or all selected atoms when message=all. Note that nothing will
C     be printed for a reference set that was not defined by the
C     patch invocation.
      IF (J-1.GT.I) THEN
      IAT = ISET(I)
      WRITE(6,'(1X,4A,3(1X,A),A)') '%PATCH-ERR: REFErence=',SET(I),
     &  '   first atom=(',SEGID(IAT),RESID(IAT),RES(IAT),TYPE(IAT),')'
      IF (WRNLEV.GT.5) THEN
      DO K = I+1, J-2
      IAT = ISET(K)
      WRITE(6,'(1X,A,29X,2A,3(1X,A),A)') '%PATCH-ERR:',
     &  '(',SEGID(IAT),RESID(IAT),RES(IAT),TYPE(IAT),')'
      END DO
      END IF
      IAT = ISET(J-1)
      WRITE(6,'(1X,A,19X,2A,3(1X,A),A)') '%PATCH-ERR:',
     &  'last atom=(',SEGID(IAT),RESID(IAT),RES(IAT),TYPE(IAT),')'
      ELSE
      IAT = ISET(I)
      WRITE(6,'(1X,4A,3(1X,A),A)') '%PATCH-ERR: REFErence=',SET(I),
     &  '   atom=(',SEGID(IAT),RESID(IAT),RES(IAT),TYPE(IAT),')'
      END IF
      I=J
      END DO
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE PATCHERRATM(NAME,IAT)
C
C Print error information about a patch atom selection
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      CHARACTER*(*) NAME
      INTEGER IAT
C begin
      IF (WRNLEV.GT.5) THEN
      IF (IAT.LE.0 .OR. IAT.GT.NATOM) THEN
      WRITE(6,'(3A)')
     &    ' %PATCH-ERR:      atom ',NAME,' =  not found'
      ELSE
      WRITE(6,'(4A,3(1X,A),A)')
     &    ' %PATCH-ERR:      atom ',NAME,' = (',
     &    SEGID(IAT),RESID(IAT),RES(IAT),TYPE(IAT),')'
      END IF
      END IF
      END
C======================================================================
      SUBROUTINE DELTIC
C
C Routine deletes atoms and all references to atoms
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
C local
      INTEGER FLAGS, MAP, INVMAP, GRPLST, LINB, LIBLO
      INTEGER XNATOM, XNNB
C dynamic
      XNATOM=NATOM
      XNNB=NNB
      FLAGS=ALLHP(INTEG4(NATOM))
      MAP=ALLHP(INTEG4(NATOM))
      INVMAP=ALLHP(INTEG4(NATOM))
      GRPLST=ALLHP(INTEG4(NATOM))
      LINB=ALLHP(INTEG4(NNB))
      LIBLO=ALLHP(INTEG4(NATOM))
      CALL DELTI2(HEAP(FLAGS),HEAP(MAP),HEAP(INVMAP),
     &            HEAP(GRPLST),HEAP(LINB),HEAP(LIBLO))
      CALL FREHP(LIBLO,INTEG4(XNATOM))
      CALL FREHP(LINB,INTEG4(XNNB))
      CALL FREHP(GRPLST,INTEG4(XNATOM))
      CALL FREHP(INVMAP,INTEG4(XNATOM))
      CALL FREHP(MAP,INTEG4(XNATOM))
      CALL FREHP(FLAGS,INTEG4(XNATOM))
      RETURN
      END
C
      SUBROUTINE DELTI2(FLAGS,MAP,INVMAP,GRPLST,LINB,LIBLO)
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'comand.inc'
      INTEGER FLAGS(*)
      INTEGER MAP(*), INVMAP(*)
      INTEGER GRPLST(*), LINB(*), LIBLO(*)
C local
      LOGICAL LSORT
      INTEGER I, J, II, OFFSET, NSLCT
C Temporary mark for deleted atoms:
      INTEGER MARK
      PARAMETER (MARK=-9999)
C
C begin
C
C defaults
      LSORT=.FALSE.
      CALL FILL4(FLAGS,NATOM,0)
      NSLCT=0
C
C parsing
      CALL PUSEND('DELETE>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('DELETE>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-delete')
C
      ELSE IF (WD(1:4).EQ.'SORT') THEN
      CALL NEXTLO('SORT=',LSORT)
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(FLAGS,NSLCT,X,Y,Z,.TRUE.)
      ELSE
      CALL CHKEND('DELETE>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      IF (NSLCT.GT.0) THEN
C
C construct MAP and INVMAP
      II=NATOM
      DO I=1,NATOM
      IF (FLAGS(I).EQ.1) THEN
      II=II-1
      MAP(I)=MARK
      END IF
      END DO
      NATOM=II
C
      OFFSET=0
      DO I=1,NATOM
      DO WHILE (FLAGS(I+OFFSET).EQ.1)
      OFFSET=OFFSET+1
      END DO
      INVMAP(I)=I+OFFSET
      END DO
      DO I=1,NATOM
      MAP(INVMAP(I))=I
      END DO
C
C fill temporary group list
      DO I=1,NGRP
      DO J=IGPBS(I)+1,IGPBS(I+1)
      GRPLST(J)=I
      END DO
      END DO
C
C now do mapping and compression of molecular topology, etc.
      CALL MAPIC(MAP,INVMAP,GRPLST,LINB,LIBLO,MARK)
      CALL CMPRIC(MARK,LSORT)
C
      END IF
C
C finally, scratch the lists
      CALL SCRATC
      CALL SHOW
      RETURN
      END
C======================================================================
      INTEGER FUNCTION PATOM(ATOM,NSET,SET,ISET,TYPE,MARK)
C
C tries to find ATOM=<reference><iupac> in SET, TYPE.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      CHARACTER*5 ATOM
      INTEGER NSET
      CHARACTER*4 SET(*)
      INTEGER ISET(*)
      CHARACTER*4 TYPE(*)
      INTEGER MARK
C local
      LOGICAL FOUND
      INTEGER I
C begin
      FOUND=.FALSE.
      I=0
      DO WHILE (.NOT. (FOUND.OR.I.GE.NSET))
      I=I+1
      IF (SET(I).EQ.'NIL ') THEN
      FOUND=(ATOM(1:4).EQ.TYPE(ISET(I)))
      ELSE
      FOUND=(SET(I)(1:1).EQ.ATOM(1:1))
      IF (FOUND) FOUND=(ATOM(2:5).EQ.TYPE(ISET(I)))
      END IF
      END DO
      IF (FOUND) THEN
      PATOM=ISET(I)
      ELSE
      PATOM=MARK
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE MAPIC(MAP,INVMAP,GRPLST,LINB,LIBLO,MARK)
C
C Subroutine maps molecular topology, coordinates, and
C constraints according to
C a MAP and its inverse INV_MAP. GRP_LST contains a group number for
C each atom and is used to redefine the group boundaries. See PATIC2
C for use of this routine.
C
C note: this routine may be modified when the molecular topology,
C coordinates, or constraints blocks are modified.
C after calling this routine a call to cmpric is necessary to
C delete marked atoms in all lists.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'cstack.inc'
      INCLUDE 'param.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'ener.inc'
      INTEGER MAP(*), INVMAP(*), GRPLST(*), LINB(*), LIBLO(*), MARK
C local
      INTEGER WORK, CWORK, LWORK
      INTEGER GROUP, I, II, J
C
C begin
      LWORK=MAX(NATOM,NBOND,NTHETA,NPHI,NIMPHI)
      WORK=ALLHP(IREAL8(LWORK))
      CWORK=CALLST(ICHAR4(NATOM))
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,1500)
1500  FORMAT(' MAPIC: Atom numbers being modified')
      END IF
C
C map all atom properties
      CALL AINDC4(INVMAP,IAC,NATOM,CSTACK(CWORK))
      CALL AINDX4(INVMAP,LOOKUP,NATOM,HEAP(WORK))
      CALL AINDL4(INVMAP,QLOOKU,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,AMASS,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,CG,NATOM,HEAP(WORK))
      CALL AINDC4(INVMAP,TYPE,NATOM,CSTACK(CWORK))
      CALL AINDC4(INVMAP,RES,NATOM,CSTACK(CWORK))
      CALL AINDC4(INVMAP,RESID,NATOM,CSTACK(CWORK))
      CALL AINDC4(INVMAP,SEGID,NATOM,CSTACK(CWORK))
      CALL AINDX4(INVMAP,IMOVE,NATOM,HEAP(WORK))
C
      CALL AINDR8(INVMAP,X,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,Y,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,Z,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,WMAIN,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,QMAIN,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,XCOMP,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,YCOMP,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,ZCOMP,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,WCOMP,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,QCOMP,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,DX,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,DY,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,DZ,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,XV,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,YV,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,ZV,NATOM,HEAP(WORK))
C
      CALL AINDR8(INVMAP,RMSD,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,KCNSTR,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,KZCNSTR,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,REFX,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,REFY,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,REFZ,NATOM,HEAP(WORK))
      CALL AINDR8(INVMAP,FBETA,NATOM,HEAP(WORK))
C
C delete stores
      DO I=1,MAXSTO
      IF (PTRSTO(I).NE.0) THEN
      CALL FREHP(PTRSTO(I),IREAL8(LENSTO(I)))
      PTRSTO(I)=0
      LENSTO(I)=0
      END IF
      END DO
C
C map atom lists
      CALL AVALUE(MAP,IB,NBOND,MARK)
      CALL AVALUE(MAP,JB,NBOND,MARK)
C
      CALL AVALUE(MAP,IT,NTHETA,MARK)
      CALL AVALUE(MAP,JT,NTHETA,MARK)
      CALL AVALUE(MAP,KT,NTHETA,MARK)
C
      CALL AVALUE(MAP,IP,NPHI,MARK)
      CALL AVALUE(MAP,JP,NPHI,MARK)
      CALL AVALUE(MAP,KP,NPHI,MARK)
      CALL AVALUE(MAP,LP,NPHI,MARK)
C
      CALL AVALUE(MAP,IM,NIMPHI,MARK)
      CALL AVALUE(MAP,JM,NIMPHI,MARK)
      CALL AVALUE(MAP,KM,NIMPHI,MARK)
      CALL AVALUE(MAP,LM,NIMPHI,MARK)
C
      CALL AINDX4(INVMAP,GRPLST,NATOM,HEAP(WORK))
C
C re-define INB, IBLO explicit non-bond-exclusion list.
      CALL AVALUE(MAP,INB,NNB,MARK)
C
      NNB=0
      DO I=1,NATOM
      IF (INVMAP(I).GT.1) THEN
      II=IBLO(INVMAP(I)-1)+1
      ELSE
      II=1
      END IF
      DO J=II,IBLO(INVMAP(I))
      IF (INB(J).NE.MARK) THEN
      NNB=NNB+1
      LINB(NNB)=INB(J)
      END IF
      END DO
      LIBLO(I)=NNB
      END DO
C
      DO I=1,NATOM
      IBLO(I)=LIBLO(I)
      END DO
      DO I=1,NNB
      INB(I)=LINB(I)
      END DO
C
C now re-define groups:
      GROUP=-1
      NGRP=0
      DO I=1,NATOM
      IF (GROUP.NE.GRPLST(I)) THEN
      GROUP=GRPLST(I)
      IF (NGRP.GE.MAXGRP) THEN
      CALL WRNDIE(-5,'MAPIC',
     & 'exceeded MAXGRP parameter --> recompile program')
      ELSE
      NGRP=NGRP+1
      IGPBS(NGRP)=I-1
      END IF
      END IF
      END DO
      IGPBS(NGRP+1)=NATOM
C
      CALL FREHP(WORK,IREAL8(LWORK))
      CALL CFREST(ICHAR4(NATOM))
      RETURN
      END
C
      SUBROUTINE CMPRIC(MARK,LSORT)
C
C     Routine deletes marked bonds, angles, ...
C     compresses molecular structure database.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INTEGER   MARK
      LOGICAL   LSORT
C local
      EXTERNAL  EXCH5
      INTEGER   I, II
      LOGICAL   CONDIT
C begin
C     compress bond lists,...
C
      II=0
      DO I=1,NBOND
      CONDIT=(IB(I).EQ.MARK.OR.JB(I).EQ.MARK)
      IF (.NOT.CONDIT) THEN
      II=II+1
      IB(II)=IB(I)
      JB(II)=JB(I)
      ACBB(II)=ACBB(I)
      ACBC(II)=ACBC(I)
      QACBB(II)=QACBB(I)
      QACBC(II)=QACBC(I)
      END IF
      END DO
      NBOND=II
C
      II=0
      DO I=1,NTHETA
      CONDIT=(IT(I).EQ.MARK.OR.JT(I).EQ.MARK.OR.KT(I).EQ.MARK)
      IF (.NOT. CONDIT) THEN
      II=II+1
      IT(II)=IT(I)
      JT(II)=JT(I)
      KT(II)=KT(I)
      ACTB(II)=ACTB(I)
      ACTC(II)=ACTC(I)
      ACTUB(II)=ACTUB(I)
      ACTUC(II)=ACTUC(I)
      QACTB(II)=QACTB(I)
      QACTC(II)=QACTC(I)
      QACTUB(II)=QACTUB(I)
      QACTUC(II)=QACTUC(I)
      END IF
      END DO
      NTHETA=II
C
      II=0
      DO I=1,NPHI
      CONDIT=(IP(I).EQ.MARK.OR.JP(I).EQ.MARK.OR.KP(I).EQ.MARK
     1        .OR.LP(I).EQ.MARK)
      IF (.NOT.CONDIT) THEN
      II=II+1
      IP(II)=IP(I)
      JP(II)=JP(I)
      KP(II)=KP(I)
      LP(II)=LP(I)
      ACPB(II)=ACPB(I)
      ACPC(II)=ACPC(I)
      ACPD(II)=ACPD(I)
      QACPB(II)=QACPB(I)
      QACPC(II)=QACPC(I)
      QACPD(II)=QACPD(I)
      END IF
      END DO
      NPHI=II
C
      II=0
      DO I=1,NIMPHI
      CONDIT=(IM(I).EQ.MARK.OR.JM(I).EQ.MARK.OR.KM(I).EQ.MARK
     1        .OR.LM(I).EQ.MARK)
      IF (.NOT.CONDIT) THEN
      II=II+1
      IM(II)=IM(I)
      JM(II)=JM(I)
      KM(II)=KM(I)
      LM(II)=LM(I)
      ACIB(II)=ACIB(I)
      ACIC(II)=ACIC(I)
      ACID(II)=ACID(I)
      QACIB(II)=QACIB(I)
      QACIC(II)=QACIC(I)
      QACID(II)=QACID(I)
      END IF
      END DO
      NIMPHI=II
C
C     finally sort all lists if LSORT is true
C
      IF (LSORT) THEN
      CALL SORT(NBOND,EXCH5,ORDER5,IB,JB,0,0,0,0,0,2)
      CALL SORT(NTHETA,EXCH5,ORDER5,IT,JT,KT,0,0,0,0,3)
      CALL SORT(NPHI,EXCH5,ORDER5,IP,JP,KP,LP,0,0,0,4)
      CALL SORT(NIMPHI,EXCH5,ORDER5,IM,JM,KM,LM,0,0,0,4)
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE DUPLIC
C
C parses duplicate statement
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
C local
      INTEGER ONATOM
C pointer
      INTEGER FLAGS, FLAGS2, INDEX, MAP, GRPLST
C begin
      ONATOM=NATOM
      FLAGS=ALLHP(INTEG4(ONATOM))
      FLAGS2=ALLHP(INTEG4(ONATOM+1))
      INDEX=ALLHP(INTEG4(ONATOM))
      MAP=ALLHP(INTEG4(ONATOM+1))
      GRPLST=ALLHP(INTEG4(ONATOM))
      CALL DUPLI2(HEAP(FLAGS),HEAP(FLAGS2),
     &            HEAP(INDEX),HEAP(MAP),HEAP(GRPLST))
      CALL FREHP(GRPLST,INTEG4(ONATOM))
      CALL FREHP(MAP,INTEG4(ONATOM+1))
      CALL FREHP(INDEX,INTEG4(ONATOM))
      CALL FREHP(FLAGS2,INTEG4(ONATOM+1))
      CALL FREHP(FLAGS,INTEG4(ONATOM))
      RETURN
      END
C======================================================================
      SUBROUTINE DUPLI2(FLAGS,FLAGS2,INDEX,MAP,GRPLST)
C
C See routine DUPLIC above.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'mtf.inc'
      INTEGER FLAGS(*), FLAGS2(0:*), INDEX(*), MAP(0:*), GRPLST(*)
C local
      INTEGER NFLAGS
      CHARACTER*4 NSEGID, NRESID
C begin
C
C defaults
      CALL FILL4(FLAGS,NATOM,0)
      NFLAGS=0
      NSEGID='    '
      NRESID='    '
C
C parsing
      CALL PUSEND('DUPLicate>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('DUPLicate>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-duplicate')
C
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
      ELSE IF (WD(1:4).EQ.'SEGI') THEN
      CALL NEXTA4('SEGId=',NSEGID)
      ELSE IF (WD(1:4).EQ.'RESI') THEN
      CALL NEXTA4('RESIdue=',NRESID)
      ELSE
      CALL CHKEND('DUPLicate>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      IF (NFLAGS.GT.0) THEN
      IF (NSEGID.EQ.'    '.AND.NRESID.EQ.'    ') THEN
      CALL WRNDIE(-5,'DUPLI2',
     & 'Non-blank SEGId or RESIdue (or both) have to be specified')
      ELSE
C
      CALL DUPLI3(NFLAGS,FLAGS,FLAGS2,INDEX,MAP,GRPLST,NSEGID,NRESID)
C
C scratch the lists
      CALL SCRATC
      CALL SHOW
      END IF
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE DUPLI3(NFLAGS,FLAGS,FLAGS2,
     &                  INDEX,MAP,GRPLST,NSEGID,NRESID)
C
C Routine duplicates the specified atoms.  It copies all atomic
C properties, such as X,Y,Z, MASS, etc.  It copies those bonds,
C angles, dihedrals for which at least one atom has been
C selected.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'funct.inc'
      INTEGER NFLAGS, FLAGS(*), FLAGS2(0:*), INDEX(*)
      INTEGER MAP(0:*), GRPLST(*)
      CHARACTER*4 NSEGID, NRESID
C local
      INTEGER IAT, II, IG, J, OGROUP, ONATOM, NN
C begin
C
C first make an index list of all selected atoms.
      FLAGS2(0)=0
      DO IAT=1,NATOM
      INDEX(IAT)=FLAGS(IAT)
      FLAGS2(IAT)=FLAGS(IAT)
      END DO
      CALL MAKIND(INDEX,NATOM,NFLAGS)
C
C initialize the atom number map
      MAP(0)=0
      DO IAT=1,NATOM
      MAP(IAT)=IAT
      END DO
C
C make a group list of the current atoms
      DO IG=1,NGRP
      DO IAT=IGPBS(IG)+1,IGPBS(IG+1)
      GRPLST(IAT)=IG
      END DO
      END DO
      ONATOM=NATOM
C
C duplicate atomic properties
      DO IAT=1,NFLAGS
      IF (NATOM.GE.MAXA) THEN
      CALL WRNDIE(-5,'DUPLI3',
     & 'exceeded MAXA parameter --> recompile program')
      ELSE
      NATOM=NATOM+1
      END IF
      IAC(NATOM)=IAC(INDEX(IAT))
      LOOKUP(NATOM)=LOOKUP(INDEX(IAT))
      QLOOKU(NATOM)=QLOOKU(INDEX(IAT))
      AMASS(NATOM)=AMASS(INDEX(IAT))
      CG(NATOM)=CG(INDEX(IAT))
      TYPE(NATOM)=TYPE(INDEX(IAT))
      RES(NATOM)=RES(INDEX(IAT))
      IF (NRESID.NE.'    ') THEN
      RESID(NATOM)=NRESID
      ELSE
      RESID(NATOM)=RESID(INDEX(IAT))
      END IF
      IF (NSEGID.NE.'    ') THEN
      SEGID(NATOM)=NSEGID
      ELSE
      SEGID(NATOM)=SEGID(INDEX(IAT))
      END IF
      IMOVE(NATOM)=IMOVE(INDEX(IAT))
      X(NATOM)=X(INDEX(IAT))
      Y(NATOM)=Y(INDEX(IAT))
      Z(NATOM)=Z(INDEX(IAT))
      WMAIN(NATOM)=WMAIN(INDEX(IAT))
      QMAIN(NATOM)=QMAIN(INDEX(IAT))
      XCOMP(NATOM)=XCOMP(INDEX(IAT))
      YCOMP(NATOM)=YCOMP(INDEX(IAT))
      ZCOMP(NATOM)=ZCOMP(INDEX(IAT))
      WCOMP(NATOM)=WCOMP(INDEX(IAT))
      QCOMP(NATOM)=QCOMP(INDEX(IAT))
      DX(NATOM)=DX(INDEX(IAT))
      DY(NATOM)=DY(INDEX(IAT))
      DZ(NATOM)=DZ(INDEX(IAT))
      XV(NATOM)=XV(INDEX(IAT))
      YV(NATOM)=YV(INDEX(IAT))
      ZV(NATOM)=ZV(INDEX(IAT))
      RMSD(NATOM)=RMSD(INDEX(IAT))
      KCNSTR(NATOM)=KCNSTR(INDEX(IAT))
      KZCNSTR(NATOM)=KZCNSTR(INDEX(IAT))
      REFX(NATOM)=REFX(INDEX(IAT))
      REFY(NATOM)=REFY(INDEX(IAT))
      REFZ(NATOM)=REFZ(INDEX(IAT))
      FBETA(NATOM)=FBETA(INDEX(IAT))
C
C update the map
      MAP(INDEX(IAT))=NATOM
      END DO
C
C delete VECTOR stores
      DO II=1,MAXSTO
      IF (PTRSTO(II).NE.0) THEN
      CALL FREHP(PTRSTO(II),IREAL8(LENSTO(II)))
      PTRSTO(II)=0
      LENSTO(II)=0
      END IF
      END DO
C
C duplicate bonds
      NN=NBOND
      DO II=1,NN
      IF (FLAGS2(IB(II)).EQ.1.OR.FLAGS2(JB(II)).EQ.1) THEN
      IF (NBOND.GE.MAXB) THEN
      CALL WRNDIE(-5,'DUPLI3',
     & 'exceeded MAXB parameter --> recompile program')
      ELSE
      NBOND=NBOND+1
      END IF
      IB(NBOND)=MAP(IB(II))
      JB(NBOND)=MAP(JB(II))
      ACBB(NBOND)=ACBB(II)
      ACBC(NBOND)=ACBC(II)
      QACBB(NBOND)=QACBB(II)
      QACBC(NBOND)=QACBC(II)
      END IF
      END DO
C
C duplicate angles
      NN=NTHETA
      DO II=1,NN
      IF (FLAGS2(IT(II)).EQ.1.OR.FLAGS2(JT(II)).EQ.1.OR.
     &      FLAGS2(KT(II)).EQ.1) THEN
      IF (NTHETA.GE.MAXT) THEN
      CALL WRNDIE(-5,'DUPLI3',
     & 'exceeded MAXT parameter --> recompile program')
      ELSE
      NTHETA=NTHETA+1
      END IF
      IT(NTHETA)=MAP(IT(II))
      JT(NTHETA)=MAP(JT(II))
      KT(NTHETA)=MAP(KT(II))
      ACTB(NTHETA)=ACTB(II)
      ACTC(NTHETA)=ACTC(II)
      ACTUB(NTHETA)=ACTUB(II)
      ACTUC(NTHETA)=ACTUC(II)
      QACTB(NTHETA)=QACTB(II)
      QACTC(NTHETA)=QACTC(II)
      QACTUB(NTHETA)=QACTUB(II)
      QACTUC(NTHETA)=QACTUC(II)
      END IF
      END DO
C
C duplicate dihedrals
      NN=NPHI
      DO II=1,NN
      IF (FLAGS2(IP(II)).EQ.1.OR.FLAGS2(JP(II)).EQ.1.OR.
     &      FLAGS2(KP(II)).EQ.1.OR.FLAGS2(LP(II)).EQ.1) THEN
      IF (NPHI.GE.MAXP) THEN
      CALL WRNDIE(-5,'DUPLI3',
     & 'exceeded MAXP parameter --> recompile program')
      ELSE
      NPHI=NPHI+1
      END IF
      IP(NPHI)=MAP(IP(II))
      JP(NPHI)=MAP(JP(II))
      KP(NPHI)=MAP(KP(II))
      LP(NPHI)=MAP(LP(II))
      ACPB(NPHI)=ACPB(II)
      ACPC(NPHI)=ACPC(II)
      ACPD(NPHI)=ACPD(II)
      QACPB(NPHI)=QACPB(II)
      QACPC(NPHI)=QACPC(II)
      QACPD(NPHI)=QACPD(II)
      END IF
      END DO
C
C duplicate impropers
      NN=NIMPHI
      DO II=1,NN
      IF (FLAGS2(IM(II)).EQ.1.OR.FLAGS2(JM(II)).EQ.1.OR.
     &      FLAGS2(KM(II)).EQ.1.OR.FLAGS2(LM(II)).EQ.1) THEN
      IF (NIMPHI.GE.MAXIMP) THEN
      CALL WRNDIE(-5,'DUPLI3',
     & 'exceeded MAXIMP parameter --> recompile program')
      ELSE
      NIMPHI=NIMPHI+1
      END IF
      IM(NIMPHI)=MAP(IM(II))
      JM(NIMPHI)=MAP(JM(II))
      KM(NIMPHI)=MAP(KM(II))
      LM(NIMPHI)=MAP(LM(II))
      ACIB(NIMPHI)=ACIB(II)
      ACIC(NIMPHI)=ACIC(II)
      ACID(NIMPHI)=ACID(II)
      QACIB(NIMPHI)=QACIB(II)
      QACIC(NIMPHI)=QACIC(II)
      QACID(NIMPHI)=QACID(II)
      END IF
      END DO
C
C duplicate explicit nonbonded exclusions
      DO IAT=1,NFLAGS
      IF (INDEX(IAT).GT.1) THEN
      II=IBLO(INDEX(IAT)-1)+1
      ELSE
      II=1
      END IF
      DO J=II,IBLO(INDEX(IAT))
      IF (NNB.GE.MAXNB) THEN
      CALL WRNDIE(-5,'DUPLI3',
     & 'exceeded MAXNB parameter --> recompile program')
      ELSE
      NNB=NNB+1
      END IF
      INB(NNB)=MAP(INB(J))
      END DO
      IBLO(MAP(INDEX(IAT)))=NNB
      END DO
C
C duplicate groups
C
C make a group list of the current atoms
      OGROUP=NGRP
      DO IAT=1,NFLAGS
      IF (GRPLST(INDEX(IAT)).NE.OGROUP) THEN
C
C make a new group
      IF (NGRP.GE.MAXGRP) THEN
      CALL WRNDIE(-5,'DUPLI3',
     & 'exceeded MAXGRP parameter --> recompile program')
      ELSE
      NGRP=NGRP+1
      END IF
      OGROUP=GRPLST(INDEX(IAT))
      IGPBS(NGRP)=MAP(INDEX(IAT))-1
      IF (IGPBS(NGRP).NE.ONATOM+IAT-1) THEN
      CALL WRNDIE(-5,'DUPL3','Programming error')
      END IF
      END IF
      END DO
      IGPBS(NGRP+1)=NATOM
C
      RETURN
      END

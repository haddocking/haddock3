      SUBROUTINE RSTRAN
C
C restraints parser.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C local
      DOUBLE PRECISION TMP
C pointer
      INTEGER IFLAGS
C parsing
      CALL NEXTWD('RESTraints>')
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-restraints')
C
      ELSEIF (WD(1:4).EQ.'PLAN') THEN
      CALL PLNPAR
      ELSE IF (WD(1:4).EQ.'DIHE') THEN
      IFLAGS=ALLHP(INTEG4(NATOM))
      CALL RSTDIH(HEAP(IFLAGS))
      CALL FREHP(IFLAGS,INTEG4(NATOM))
      ELSE IF (WD(1:4).EQ.'HARM') THEN
C
C turn on the harmonic energy term
      QENER(SSHARM)=.TRUE.
C
      CALL PUSEND('HARMonic>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('HARMonic>')
      CALL MISCOM('HARMonic>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-restraints-harmonic')
C
      ELSE IF (WD(1:1).EQ.'?') THEN
      WRITE (6,'(A,I3)')' HARM: exponent=',KCEXPN
      IF(QPLANAR) THEN
      WRITE (6,'(A)')
     &' HARM: planar restraints defined for atoms with neg. harmonic'
      WRITE (6,'(A)')'       energy constant,'
      WRITE (6,'(A,3F8.3)')'       plane defined by the normal vector',
     &                     PNORMAL
      END IF
C====
      ELSE IF (WD(1:4).EQ.'NORM') THEN
      CALL NEXTVF('NORMal-vector=', PNORMAL)
      TMP=SQRT(PNORMAL(1)**2+PNORMAL(2)**2+PNORMAL(3)**2)
      IF (TMP.LE.RSMALL) THEN
      QPLANAR=.FALSE.
      ELSE
      PNORMAL(1)=PNORMAL(1)/TMP
      PNORMAL(2)=PNORMAL(2)/TMP
      PNORMAL(3)=PNORMAL(3)/TMP
      QPLANAR=.TRUE.
      END IF
C====
      ELSE IF (WD(1:4).EQ.'EXPO') THEN
      CALL NEXTI('EXPOnent=',KCEXPN)
C====
      ELSE
      CALL CHKEND('HARMonic>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C===============================================
C==== A. Bonvin Z-position restraining within a slice
C===============================================
      ELSE IF (WD(1:4).EQ.'ZHAR') THEN
C
C turn on the harmonic energy term
      QENER(SSZHAR)=.TRUE.
C
      CALL PUSEND('ZHARmonic>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('ZHARmonic>')
      CALL MISCOM('ZHARmonic>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-restraints-zharmonic')
C
      ELSE IF (WD(1:1).EQ.'?') THEN
      WRITE (6,'(A,I3)')' ZHAR: exponent=',KZCEXP
C====
      ELSE IF (WD(1:4).EQ.'ZMAX') THEN
      CALL NEXTF('ZMAX value=', ZHARMAX)
      ELSE IF (WD(1:4).EQ.'ZMIN') THEN
      CALL NEXTF('ZMIN value=', ZHARMIN)
C====
      ELSE IF (WD(1:4).EQ.'EXPO') THEN
      CALL NEXTI('EXPOnent=',KZCEXP)
C====
      ELSE
      CALL CHKEND('ZHARmonic>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C===============================================
C===============================================
      ELSE
      CALL DSPERR('RSTRAN','Unkown Restraints Option.')
      END IF
C
      CALL EFLAGSYMBOLS
C
      RETURN
      END
C======================================================================
      SUBROUTINE FIX
C
C fixed atom setup parser.
C
C For syntax see main parsing loop "HELP"
C
C Authors: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C pointer
      INTEGER IFLAGS
C begin
      IFLAGS=ALLHP(INTEG4(NATOM))
      CALL FIX2(HEAP(IFLAGS))
      CALL FREHP(IFLAGS,INTEG4(NATOM))
      RETURN
      END
C======================================================================
      SUBROUTINE FIX2(IFLAGS)
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'update.inc'
      INCLUDE 'heap.inc'
      INTEGER IFLAGS(*)
C local
      INTEGER ISELCT, I, NFIXED
C parameter
      DOUBLE PRECISION ZERO, RAD, ONE
      PARAMETER (RAD=PI/180.0D0, ZERO=0.0D0, ONE=1.0D0)
C begin
C
C parsing
      CALL PUSEND('FIX>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('FIX>')
      CALL MISCOM('FIX>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-fix')
C
      ELSE IF (WD(1:1).EQ.'?') THEN
      NFIXED=0
      DO I=1,NATOM
      IF (IMOVE(I).EQ.1) THEN
      NFIXED=NFIXED+1
      END IF
      END DO
      WRITE (6,'(/A,I5,A)')' FIX : ',NFIXED,' atoms fixed'
C
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      DO I=1,NATOM
      IFLAGS(I)=MAX(0,IMOVE(I))
      END DO
      CALL SELCTA(IFLAGS,ISELCT,X,Y,Z,.TRUE.)
      NFIXED=ISELCT
      IF (ISELCT.GE.0) THEN
      DO I=1,NATOM
      IMOVE(I)=IFLAGS(I)
      END DO
      END IF
C
      ELSE
      CALL CHKEND('FIX>',DONE)
      END IF
C
      END IF
      END DO
      DONE=.FALSE.
C
      RETURN
      END
C====================================================================
      SUBROUTINE RSTDIH(IFLAGS)
C
C setup dihedral angle restraints.
C
C Author: Axel T. Brunger
C modified for cross-validation by Alexandre Bonvin 12/6/95
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'mtf.inc'
      INTEGER IFLAGS(*)
C local
      INTEGER NEWASS,ITEMP
C parameter
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C begin
C
C allocate space if not already allocated
      IF (MCMPHI.EQ.0) THEN
      CALL CDIHP(400)
      END IF
C parsing
      CALL PUSEND('DIHEDRAL>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('DIHEDRAL>')
      CALL MISCOM('DIHEDRAL>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-restraints-dihedral')
C
      ELSE IF (WD(1:4).EQ.'NASS') THEN
      NEWASS=MCMPHI
      CALL NEXTI('NASSign=',NEWASS)
      IF (NEWASS.NE.MCMPHI) THEN
      WRITE(6,'(A,I7,A)')
     & ' RSTDIH: allocating space for ',NEWASS,' assignments.'
      CALL CDIHP(NEWASS)
      END IF
C====================================================================
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      CALL CDIRES
      QENER(SSCDIH)=.FALSE.
C====================================================================
      ELSE IF (WD(1:4).EQ.'SCAL') THEN
      CALL NEXTF('SCALe=',CNSCA)
C====================================================================
      ELSE IF (WD(1:1).EQ.'?') THEN
      CALL CDIPRI(PUNIT,CNPHI,CNSCA,HEAP(ICIP),HEAP(ICJP),
     &            HEAP(ICKP),HEAP(ICLP),HEAP(ICICP),HEAP(ICCPD),
     &            HEAP(ICCPE),HEAP(ICCPB),HEAP(ICCPC),HEAP(ICCPO),
     &            -ONE,HEAP(HPDCV),0,DIHICV)
      IF (DIHICV.GT.0) THEN
      CALL CDIPRI(PUNIT,CNPHI,CNSCA,HEAP(ICIP),HEAP(ICJP),
     &            HEAP(ICKP),HEAP(ICLP),HEAP(ICICP),HEAP(ICCPD),
     &            HEAP(ICCPE),HEAP(ICCPB),HEAP(ICCPC),HEAP(ICCPO),
     &            -ONE,HEAP(HPDCV),1,DIHICV)
      END IF
C====================================================================
      ELSE IF (WD(1:4).EQ.'ASSI') THEN
      QENER(SSCDIH)=.TRUE.
      CALL CDIASS(MCMPHI,CNPHI,HEAP(ICIP),HEAP(ICJP),HEAP(ICKP),
     &            HEAP(ICLP),HEAP(ICICP),HEAP(ICCPD),HEAP(ICCPE),
     &            HEAP(ICCPB),HEAP(ICCPC),HEAP(ICCPO),HEAP(HPDCV),
     &            IFLAGS,NATOM,X,Y,Z)
C====================================================================
      ELSE IF (WD(1:2).EQ.'CV') THEN
      CALL NEXTI('CV excluded partition number:',DIHICV)
C A. Bonvin
C store working set
      IF (DIHICV.GT.0) THEN
      CALL STWCV(CNPHI,DIHICV,HEAP(ICIP),HEAP(ICJP),HEAP(ICKP),
     &           HEAP(ICLP),HEAP(ICICP),HEAP(ICCPD),HEAP(ICCPE),
     &           HEAP(ICCPB),HEAP(ICCPC),HEAP(ICCPO),HEAP(HPDCV),
     &           CWNPHI,HEAP(IWCIP),HEAP(IWCJP),HEAP(IWCKP),
     &           HEAP(IWCLP),HEAP(IWCCPD),HEAP(IWCCPE),
     &           HEAP(IWCCPB),HEAP(IWCCPC),HEAP(IWCCPO))
      END IF
C====================================================================
      ELSE IF (WD(1:4).EQ.'PART') THEN
      CALL NEXTI('number of PARTitions:',ITEMP)
      ITEMP=MAX(0,ITEMP)
      CALL DIHCVS(CNPHI,HEAP(HPDCV),ITEMP)
      IF (ITEMP.EQ.0) THEN
      DIHICV=0
      END IF
C====================================================================
      ELSE
      CALL CHKEND('DIHEDRAL>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      CALL EFLAGSYMBOLS
C
      RETURN
      END
C===================================================================
      SUBROUTINE CDIINI
C
C Routine initializes the CONStraints DIHEdral list
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
C begin
C initialize heap pointers of dynamic data structure
      ICIP=0
      ICJP=0
      ICKP=0
      ICLP=0
      ICICP=0
      ICCPD=0
      ICCPE=0
      ICCPB=0
      ICCPC=0
      ICCPO=0
      IWCIP=0
      IWCJP=0
      IWCKP=0
      IWCLP=0
      IWCCPD=0
      IWCCPE=0
      IWCCPB=0
      IWCCPC=0
      IWCCPO=0
      HPDCV=0
      CALL CDIRES
      RETURN
      END
C====================================================================
      SUBROUTINE CDIRES
C
C Routine resets the CONStraints DIHEdral list
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
C parameter
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C begin
      CNPHI=0
      CWNPHI=0
      CNSCA=ONE
      DIHICV=0
      RETURN
      END
C====================================================================
      SUBROUTINE CDIHP(NEW)
C
C Routine allocates space for CONStraints DIHEdral list
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INTEGER NEW
C begin
      IF (ICIP.NE.0) CALL FREHP(ICIP,INTEG4(MCMPHI))
      IF (ICJP.NE.0) CALL FREHP(ICJP,INTEG4(MCMPHI))
      IF (ICKP.NE.0) CALL FREHP(ICKP,INTEG4(MCMPHI))
      IF (ICLP.NE.0) CALL FREHP(ICLP,INTEG4(MCMPHI))
      IF (ICICP.NE.0) CALL FREHP(ICICP,INTEG4(LCICP))
      IF (ICCPD.NE.0) CALL FREHP(ICCPD,INTEG4(MCMPHI))
      IF (ICCPE.NE.0) CALL FREHP(ICCPE,INTEG4(MCMPHI))
      IF (ICCPB.NE.0) CALL FREHP(ICCPB,IREAL8(MCMPHI))
      IF (ICCPC.NE.0) CALL FREHP(ICCPC,IREAL8(MCMPHI))
      IF (ICCPO.NE.0) CALL FREHP(ICCPO,IREAL8(MCMPHI))
      IF (HPDCV.NE.0) CALL FREHP(HPDCV,INTEG4(MCMPHI))
C
      IF (IWCIP.NE.0) CALL FREHP(IWCIP,INTEG4(MCMPHI))
      IF (IWCJP.NE.0) CALL FREHP(IWCJP,INTEG4(MCMPHI))
      IF (IWCKP.NE.0) CALL FREHP(IWCKP,INTEG4(MCMPHI))
      IF (IWCLP.NE.0) CALL FREHP(IWCLP,INTEG4(MCMPHI))
      IF (IWCCPD.NE.0) CALL FREHP(IWCCPD,INTEG4(MCMPHI))
      IF (IWCCPE.NE.0) CALL FREHP(IWCCPE,INTEG4(MCMPHI))
      IF (IWCCPB.NE.0) CALL FREHP(IWCCPB,IREAL8(MCMPHI))
      IF (IWCCPC.NE.0) CALL FREHP(IWCCPC,IREAL8(MCMPHI))
      IF (IWCCPO.NE.0) CALL FREHP(IWCCPO,IREAL8(MCMPHI))
C
      MCMPHI=0
      CNPHI=0
      CWNPHI=0
      ICIP=0
      ICJP=0
      ICKP=0
      ICLP=0
      ICICP=0
      ICCPD=0
      ICCPE=0
      ICCPB=0
      ICCPC=0
      ICCPO=0
      LCICP=0
      IWCIP=0
      IWCJP=0
      IWCKP=0
      IWCLP=0
      IWCCPD=0
      IWCCPE=0
      IWCCPB=0
      IWCCPC=0
      IWCCPO=0
      HPDCV=0
C
C allocate new space
      IF (NEW.GT.0) THEN
      MCMPHI=NEW
      ICIP=ALLHP(INTEG4(MCMPHI))
      ICJP=ALLHP(INTEG4(MCMPHI))
      ICKP=ALLHP(INTEG4(MCMPHI))
      ICLP=ALLHP(INTEG4(MCMPHI))
      LCICP=NATOM
      ICICP=ALLHP(INTEG4(LCICP))
      ICCPD=ALLHP(INTEG4(MCMPHI))
      ICCPE=ALLHP(INTEG4(MCMPHI))
      ICCPB=ALLHP(IREAL8(MCMPHI))
      ICCPC=ALLHP(IREAL8(MCMPHI))
      ICCPO=ALLHP(IREAL8(MCMPHI))
C
      IWCIP=ALLHP(INTEG4(MCMPHI))
      IWCJP=ALLHP(INTEG4(MCMPHI))
      IWCKP=ALLHP(INTEG4(MCMPHI))
      IWCLP=ALLHP(INTEG4(MCMPHI))
      IWCCPD=ALLHP(INTEG4(MCMPHI))
      IWCCPE=ALLHP(INTEG4(MCMPHI))
      IWCCPB=ALLHP(IREAL8(MCMPHI))
      IWCCPC=ALLHP(IREAL8(MCMPHI))
      IWCCPO=ALLHP(IREAL8(MCMPHI))
      HPDCV=ALLHP(INTEG4(MCMPHI))
      END IF
      RETURN
      END
C====================================================================
      SUBROUTINE CDIPRI(PUNIT,CNPHI,CNSCA,CIP,CJP,CKP,CLP,
     &                  CICP,CCPD,CCPE,CCPB,CCPC,CCPO,THRESH,
     &                  DIHCV,ITEST,DIHICV)
C
C Routine prints the CONStraints DIHEdral list
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INTEGER PUNIT, CNPHI, ITEST, DIHICV
      DOUBLE PRECISION CNSCA
      INTEGER CIP(*), CJP(*), CKP(*), CLP(*), CICP(*), CCPD(*)
      INTEGER CCPE(*), DIHCV(*)
      DOUBLE PRECISION CCPB(*), CCPC(*), CCPO(*), THRESH
C local
      INTEGER I, NRMS, NVIOL
      DOUBLE PRECISION EDUMMY, RMS
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
C parameter
      DOUBLE PRECISION ZERO, ONE, RAD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, RAD=PI/180.0D0)
C begin
      RMS=ZERO
      NRMS=0
      NVIOL=0
C
      WRITE(6,'(A,I6)') ' Total number of dihedral angle restraints=',
     &        CNPHI
      WRITE(6,'(A,F10.4)') '  overall scale =',CNSCA
      IF (DIHICV.GT.0) THEN
      IF (ITEST.EQ.0) THEN
      WRITE(6,'(A)')
     & ' $$$$$$$$$$$$$$$$$$ working set $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      ELSE
      WRITE(6,'(A,I5,A)')
     & ' $$$$$$$$$$$$$$$$$$$$ test set (TEST=',DIHICV,
     & ')  $$$$$$$$$$$$$$$$$$$$$$'
      END IF
      END IF
      DO I=1,CNPHI
      IF ((ITEST.EQ.0.AND.DIHCV(I).NE.DIHICV).OR.
     &    (ITEST.EQ.1.AND.DIHCV(I).EQ.DIHICV)) THEN
      CALL ETOR(EDUMMY,EDUMMY,CIP(I),CJP(I),CKP(I),CLP(I),1,1,
     &     CCPC(I),
     &     CCPD(I),CCPB(I),'GENERAL',CCPO(I),CCPE(I),'ANAL',
     &     CNSCA,ONE)
      IF (ABS(PCDATA(PCDEVI)).GE.THRESH) THEN
      WRITE(PUNIT,'(A)') ' ========================================'
      WRITE(PUNIT,'(4(1X,A))')
     &  SEGID(CIP(I)),RESID(CIP(I)),RES(CIP(I)),TYPE(CIP(I))
      WRITE(PUNIT,'(4(1X,A))')
     &  SEGID(CJP(I)),RESID(CJP(I)),RES(CJP(I)),TYPE(CJP(I))
      WRITE(PUNIT,'(4(1X,A))')
     &  SEGID(CKP(I)),RESID(CKP(I)),RES(CKP(I)),TYPE(CKP(I))
      WRITE(PUNIT,'(4(1X,A))')
     &  SEGID(CLP(I)),RESID(CLP(I)),RES(CLP(I)),TYPE(CLP(I))
      WRITE(PUNIT,'(5(A,F9.3))')
     &  ' Dihedral=',PCDATA(PCGEOM),'  Energy=',PCDATA(PCENER),
     &  ' C=',PCDATA(PCCONS),' Equil=',PCDATA(PCEQUI), ' Delta=',
     &  PCDATA(PCDEVI)
      WRITE(PUNIT,'(A,F8.3,A,I3)')
     $  ' Range=',CCPO(I)/RAD,' Exponent=',CCPE(I)
      NVIOL=NVIOL+1
      END IF
      NRMS=NRMS+1
      RMS=RMS+PCDATA(PCDEVI)**2
      END IF
      END DO
C
      WRITE(6,'(A,I5)')
     &     ' Number of dihedral angle restraints=', NRMS
      IF (THRESH.GT.ZERO) THEN
      WRITE(6,'(A,F8.3,A,I5)')
     &     ' Number of violations greater than ', THRESH, ': ', NVIOL
      END IF
      DBPREC=NVIOL
      IF (ITEST.EQ.0) THEN
      CALL DECLAR( 'VIOLATIONS', 'DP', ' ', DBCOMP, DBPREC )
      DBPREC=NRMS
      CALL DECLAR( 'NUMBER', 'DP', ' ', DBCOMP, DBPREC )
      ELSE
      CALL DECLAR( 'TEST_VIOLATIONS', 'DP', ' ', DBCOMP, DBPREC )
      ENDIF
      IF (NRMS.GT.0) THEN
      DBPREC=RMS/NRMS
      DBPREC=SQRT(DBPREC)
      ELSE
      DBPREC = ZERO
      END IF
      IF (ITEST.EQ.0) THEN
      CALL DECLAR( 'RMS', 'DP', ' ', DBCOMP, DBPREC )
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
      ELSE
      CALL DECLAR( 'TEST_RMS', 'DP', ' ', DBCOMP, DBPREC )
      END IF
      WRITE(6,'(A,F8.3)') ' RMS deviation=',DBPREC
C
      RETURN
      END
C====================================================================
      SUBROUTINE CDIASS(MCMPHI,CNPHI,CIP,CJP,CKP,CLP,
     &                  CICP,CCPD,CCPE,CCPB,CCPC,CCPO,DIHCV,
     &                  IFLAGS,NATOM,X,Y,Z)
C
C Routine assigns a new dihedral angle restraint.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER MCMPHI, CNPHI
      INTEGER CIP(*), CJP(*), CKP(*), CLP(*), CICP(*), CCPD(*)
      INTEGER CCPE(*)
      DOUBLE PRECISION CCPB(*), CCPC(*), CCPO(*)
      INTEGER DIHCV(*),IFLAGS(*), NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*)
C local
      LOGICAL ERR
      INTEGER I, ISELCT, CDIHE(4)
C parameter
      DOUBLE PRECISION RAD
      PARAMETER (RAD=PI/180.0D0)
C begin
      IF (CNPHI.GE.MCMPHI) THEN
      WRITE(6,'(A)')
     & ' %CSTRAN-ERR: allocation for assignments exceeded',
     & '              increase NASSignment and run again.'
      CALL WRNDIE(-1,'CSTRAN',
     & 'exceeded allocation for assigments')
      ELSE
      CNPHI=CNPHI+1
      END IF
      ERR=.FALSE.
      DO I=1,4
      CALL SELCTA(IFLAGS,ISELCT,X,Y,Z,.TRUE.)
      IF (ISELCT.NE.1) THEN
      WRITE(6,'(A)')
     & ' %CSTRAN-ERR: selection has to contain exactly one atom.'
      ERR=.TRUE.
      END IF
      CALL MAKIND(IFLAGS,NATOM,ISELCT)
      CDIHE(I)=IFLAGS(1)
      END DO
      CIP(CNPHI)=CDIHE(1)
      CJP(CNPHI)=CDIHE(2)
      CKP(CNPHI)=CDIHE(3)
      CLP(CNPHI)=CDIHE(4)
      CALL NEXTF('force-constant=',CCPC(CNPHI))
      CALL NEXTF('minimum-angle=',CCPB(CNPHI))
      CALL NEXTF('range-angle=',CCPO(CNPHI))
      CALL NEXTI('exponent=',CCPE(CNPHI))
      CCPB(CNPHI)=CCPB(CNPHI)*RAD
      CCPO(CNPHI)=ABS(CCPO(CNPHI)*RAD)
      CCPD(CNPHI)=0
C include the dihedral in the working set
      DIHCV(CNPHI)=-1
      DO I=1,NATOM
      CICP(I)=0
      END DO
      IF (ERR) CNPHI=CNPHI-1
      RETURN
      END
C===================================================================
      SUBROUTINE ECNSTR(EC,REFX,REFY,REFZ,KCNSTR,NATOM,KCEXPN,
     &                  QPLANAR,PNORMAL)
C
C Routine computes harmonic atomic restraints
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'consta.inc'
      DOUBLE PRECISION EC, PNORMAL(*)
      DOUBLE PRECISION REFX(*), REFY(*), REFZ(*), KCNSTR(*)
      INTEGER NATOM
      INTEGER KCEXPN
      LOGICAL QPLANAR
C local
      INTEGER I
      DOUBLE PRECISION  CEXPN, E1, C, AX, AY, AZ, S, DF
      DOUBLE PRECISION  R
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      CEXPN=KCEXPN
      EC=ZERO
      DO I=1,NATOM
      C=KCNSTR(I)
      IF (C.GT.ZERO .OR. (C.LT.ZERO .AND. .NOT.QPLANAR)) THEN
      AX=X(I)-REFX(I)
      AY=Y(I)-REFY(I)
      AZ=Z(I)-REFZ(I)
      S=AX*AX+AY*AY+AZ*AZ
      IF (S.GT.RSMALL) THEN
      R=SQRT(S)
      E1=C*R**KCEXPN
      DF=CEXPN*E1/S
      EC=EC+E1
      DX(I)=DX(I)+AX*DF
      DY(I)=DY(I)+AY*DF
      DZ(I)=DZ(I)+AZ*DF
      END IF
      ELSE IF (C.LT.ZERO .AND. QPLANAR) THEN
      AX=PNORMAL(1)*(X(I)-REFX(I))
      AY=PNORMAL(2)*(Y(I)-REFY(I))
      AZ=PNORMAL(3)*(Z(I)-REFZ(I))
      R=AX+AY+AZ
      S=R*R
      IF (S.GT.RSMALL) THEN
      E1=ABS(C)*R**KCEXPN
      DF=CEXPN*E1/S
      EC=EC+E1
      DX(I)=DX(I)+R*PNORMAL(1)*DF
      DY(I)=DY(I)+R*PNORMAL(2)*DF
      DZ(I)=DZ(I)+R*PNORMAL(3)*DF
      END IF
      END IF
      END DO
      RETURN
      END
C===================================================================
      SUBROUTINE EZHARM(EC,KZCNSTR,NATOM,KZCEXP,
     &                  ZHARMAX,ZHARMIN)
C
C Routine computes harmonic z atomic restraints
C
C Author: Alexandre M.J.J. Bonvin
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'consta.inc'
      DOUBLE PRECISION EC 
      DOUBLE PRECISION KZCNSTR(*), ZHARMAX, ZHARMIN, ZHARMID
      INTEGER NATOM,NZ
      INTEGER KZCEXP
C local
      INTEGER I
      DOUBLE PRECISION  CEXPN, E1, C, AX, AY, AZ, ZZ, S, DF
      DOUBLE PRECISION  R,AVZ,AVZMIN,AVZMAX,ZLIM
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      CEXPN=KZCEXP
      EC=ZERO
      NZ=0
      AVZ=ZERO
      ZHARMID=(ZHARMAX-ZHARMIN)/2.0
      
      DO I=1,NATOM
      C=KZCNSTR(I)
      IF (C.GT.ZERO) THEN
      ZZ=Z(I)
      IF (ZZ.GT.ZHARMAX) THEN
      AZ=ZZ-ZHARMAX
      ELSE IF (ZZ.LT.ZHARMIN) THEN
      AZ=ZZ-ZHARMIN
      ELSE
      AZ=0
      ENDIF
      S=AZ*AZ
      IF (S.GT.RSMALL) THEN
      R=SQRT(S)
      E1=C*R**KZCEXP
      DF=CEXPN*E1/S
      EC=EC+E1
      DX(I)=ZERO
      DY(I)=ZERO
      DZ(I)=DZ(I)+AZ*DF
      END IF
      ELSE IF (C.LT.ZERO) THEN
      AVZ=AVZ+ZZ
      NZ=NZ+1
      END IF
      END DO
      IF (NZ.GT.0) THEN
      AVZ=AVZ/NZ
      IF (AVZ.GT.ZHARMID) THEN
      ZLIM=ZHARMAX
      ELSE IF (AVZ.LT.ZHARMID) THEN
      ZLIM=ZHARMIN
      ENDIF
C
      DO I=1,NATOM
      C=KZCNSTR(I)
      IF (C.LT.ZERO) THEN
      AZ=0
      ZZ=Z(I)
      IF (ZZ.LT.ZHARMAX) THEN
      IF (ZZ.GT.ZHARMIN) THEN
      AZ=ZZ-ZLIM
      END IF
      END IF
      S=AZ*AZ
      IF (S.GT.RSMALL) THEN
      R=SQRT(S)
      E1=-C*R**KZCEXP
      DF=CEXPN*E1/S
      EC=EC+E1
      DX(I)=ZERO
      DY(I)=ZERO
      DZ(I)=DZ(I)+AZ*DF
      END IF
      END IF      
      END DO
      END IF

      RETURN
      END
C======================================================================
      SUBROUTINE CSTRES
C
C initialize some restraint variables
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
C harmonic constraint
      KCEXPN=2
      KZCEXP=2
      QPLANAR=.FALSE.
      RETURN
      END
C
C====================================================================
      SUBROUTINE DIHCVS(DIHNUM,DIHCV,PART)
C
C Routine partitions DIHEDRAL data into PART sets.
C DIHCV will contain integer numbers between 1 and PART.
C
C Author: Axel T. Brunger
C Modifed for dihedrals: Alexandre Bonvin 12/6/95
C
C
C     IMPLICIT NONE
C I/O
      INTEGER DIHNUM, DIHCV(*)
      INTEGER PART
C local
      INTEGER I, P, NP, NTRYTOT, NRETRY
      DOUBLE PRECISION RNUM
      NTRYTOT = 0
C begin
  100 CONTINUE
      NRETRY = 0
      IF (PART.GT.0) THEN
      DO I=1,DIHNUM
      CALL GGUBFS(RNUM)
      DIHCV(I)=MAX(1,MIN(PART,INT(RNUM*PART)+1))
      END DO
C
      IF (PART .EQ. DIHNUM) THEN
      DO I=1,DIHNUM
      DIHCV(I)=I
      END DO
      END IF
C
      DO P=1,PART
      NP=0
      DO I=1,DIHNUM
      IF (DIHCV(I).EQ.P) THEN
      NP=NP+1
      END IF
      END DO
      IF (NP .EQ. 0) THEN
        NRETRY = 1
      ENDIF
      WRITE(6,'(A,I3,A,I5,A)') ' For set ',P,
     & ' there are ',NP,' dihedral angle restraints.'
      END DO
      ELSE
      WRITE(6,'(A)')
     & ' Data are not partitioned or partitioning removed.'
      DO I=1,DIHNUM
      DIHCV(I)=-1
      END DO
      END IF
      NTRYTOT = NTRYTOT + NRETRY
      IF (NRETRY.GT.0.AND.NTRYTOT.LE.10) THEN
      WRITE(6,'(A)')
     & ' Test set with 0 constraints! New trial...'
      GOTO 100
      ELSE IF (NTRYTOT .GT. 10) THEN
      CALL WRNDIE(-1,'DIHCVS',
     & 'Unable to partition the dihedral data within ten trials')
      ENDIF
C
      RETURN
      END
C
C====================================================================
      SUBROUTINE STWCV(CNPHI,IPART,CIP,CJP,CKP,
     &           CLP,CICP,CCPD,CCPE,CCPB,CCPC,CCPO,DIHCV,
     &           CWNPHI,WCIP,WCJP,WCKP,
     &           WCLP,WCCPD,WCCPE,WCCPB,WCCPC,WCCPO)
C
C store dihedral angle restraints correponding to the current working
C set for cross-validation
C Author: Alexandre Bonvin
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INTEGER CNPHI, CWNPHI, IPART
      INTEGER CIP(*), CJP(*), CKP(*), CLP(*), CICP(*), CCPD(*)
      INTEGER CCPE(*), DIHCV(*)
      INTEGER WCIP(*), WCJP(*), WCKP(*), WCLP(*), WCCPD(*)
      INTEGER WCCPE(*)
      DOUBLE PRECISION CCPB(*), CCPC(*), CCPO(*)
      DOUBLE PRECISION WCCPB(*), WCCPC(*), WCCPO(*)
C local
      INTEGER I
C
      CWNPHI=0
      DO I=1,CNPHI
      IF (DIHCV(I).NE.IPART) THEN
      CWNPHI=CWNPHI+1
      WCIP(CWNPHI)=CIP(I)
      WCJP(CWNPHI)=CJP(I)
      WCKP(CWNPHI)=CKP(I)
      WCLP(CWNPHI)=CLP(I)
      WCCPD(CWNPHI)=CCPD(I)
      WCCPE(CWNPHI)=CCPE(I)
      WCCPB(CWNPHI)=CCPB(I)
      WCCPC(CWNPHI)=CCPC(I)
      WCCPO(CWNPHI)=CCPO(I)
      END IF
      END DO
      WRITE(6,'(A,I3,A,I5,A)') ' Working set ',IPART,
     & ' containing ',CWNPHI,' dihedral angle restraints.'
C
      RETURN
      END

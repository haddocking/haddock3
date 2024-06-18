      SUBROUTINE ENERGY
C
C Main target function routine
C
C Authors: Thomas Simonson and Axel T. Brunger
C
C Uses multiple Pairs of Interacting Groups (PIGs) and atom-based parameters
C
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'param.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'funct.inc'
C     ------------------------------------------------------------------
C     FANTALIN stuff (2003)
C
      double precision fenergtot
      common /fenerg/fenergtot
C     ------------------------------------------------------------------
C local
      INTEGER I, NDIM, N
      DOUBLE PRECISION TOTAL
      DOUBLE COMPLEX DCVAL
CCC modifcation ATB 4/27/08
      DOUBLE PRECISION EDUMMY
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
      IF (TIMER.GT.0) CALL DSPCPU(' ENERGY: at entry ')
C
C initialization
      DO I=1,NENERT
      RENR(I)=ZERO
      END DO
      DO I=1,NATOM
      DX(I)=ZERO
      DY(I)=ZERO
      DZ(I)=ZERO
      END DO
C
C
C loop over all active PIGs and call routines that are
C sensitive to PIGs.
      DO N=1,NPIG
      CALL ENERG2(HEAP(IINTER(N)),N)
      DO I=1,NENERT
      RENR(I)=RENR(I)+REN2(I)
      ENDDO
      ENDDO
C
C call energy routines that are insensitive to double-selections
      IF (QENER(SSHARM)) THEN
      CALL ECNSTR(RENR(SSHARM),REFX,REFY,REFZ,KCNSTR,NATOM,KCEXPN,
     &            QPLANAR,PNORMAL)
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after harm. CONStr. term ')
      END IF
C
C Z-harmonic restraints - A. Bonvin, May 2013
C
      IF (QENER(SSZHAR)) THEN
      CALL EZHARM(RENR(SSZHAR),KZCNSTR,NATOM,KZCEXP,ZHARMAX,ZHARMIN)
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after z-harm. CONStr. term')
      END IF
C
      IF (QENER(SSNOE)) THEN
      CALL ENOE(RENR(SSNOE),'NOAN',0)
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after NOE constr. term')
      END IF
C
      IF (QENER(SSDG)) THEN
      CALL EMMDG(RENR(SSDG))
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after DG term')
      END IF
C
      IF (QENER(SSCDIH)) THEN
C
C cross-validation on dihedral angles: use only working set
C
      IF (DIHICV.GT.0) THEN
      CALL ETOR(RENR(SSCDIH),EDUMMY,HEAP(IWCIP),HEAP(IWCJP),HEAP(IWCKP),
     &          HEAP(IWCLP),HEAP(ICICP),CWNPHI,HEAP(IWCCPC),
     &          HEAP(IWCCPD),HEAP(IWCCPB),'GENERAL',
     &          HEAP(IWCCPO),HEAP(IWCCPE),'NOAN',CNSCA,ONE)
      ELSE
      CALL ETOR(RENR(SSCDIH),EDUMMY,HEAP(ICIP),HEAP(ICJP),HEAP(ICKP),
     &          HEAP(ICLP),HEAP(ICICP),CNPHI,HEAP(ICCPC),
     &          HEAP(ICCPD),HEAP(ICCPB),'GENERAL',
     &          HEAP(ICCPO),HEAP(ICCPE),'NOAN',CNSCA,ONE)
      ENDIF
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after Cons-DIHEdral term ')
      END IF
C
      IF (QENER(SSXREF)) THEN
      CALL XRENER(RENR(SSXREF))
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after XREFIne term')
      END IF
C
      IF (QENER(SSNCS)) THEN
      CALL ENENCS(RENR(SSNCS))
      IF (TIMER.GT.1) THEN
      CALL DSPCPU(' ENERGY: after non-crystallog. restraint term')
      END IF
      END IF
C
      IF (QENER(SSPLN)) THEN
      CALL ENEPLN(RENR(SSPLN))
      IF (TIMER.GT.1) THEN
      CALL DSPCPU(' ENERGY: after planarity restraint term')
      END IF
      END IF
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
C
C three-bond J coupling energy
C
      IF (QENER(SSJCOUP)) THEN
      CALL ECOUP(RENR(SSJCOUP), 'ENERGY ')
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after J coupling term')
      END IF
C
C carbon chemical shift energy
C
      IF (QENER(SSCARB)) THEN
      CALL ECSHIFT(RENR(SSCARB), 'ENERGY ')
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after carbon shift term')
      END IF
C
C proton chemical shift energy
C
      IF (QENER(SSPROT)) THEN
      CALL EPROTONSHIFT(RENR(SSPROT), 'ENERGY ')
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after proton shift term')
      END IF
C
C one bond couplings energy
C
      IF (QENER(SSONEJ)) THEN
      CALL EONEBOND(RENR(SSONEJ), 'ENERGY ')
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after one bond coup term')
      END IF
C
C Intervector projection angle energy
      IF(QENER(SSVEAN)) THEN
      CALL EVEAN(RENR(SSVEAN),'ENERGY ')
      IF (TIMER.GT.1)
     $     CALL DSPCPU(' ENERGY: after VEAN-Anisotropy term')
      ENDIF
C
C Susceptibility anisotropy energy
      IF(QENER(SSSANI)) THEN
      CALL ESANI(RENR(SSSANI),'ENERGY ')
      IF (TIMER.GT.1)
     $     CALL DSPCPU(' ENERGY: after Susc-Anisotropy term')
      ENDIF
C
C Diffusion anisotropy energy
      IF(QENER(SSDANI)) THEN
      CALL EDANI(RENR(SSDANI),'ENERGY ')
      IF (TIMER.GT.1)
     $     CALL DSPCPU(' ENERGY: after Diff-Anisotropy term')
      ENDIF
C
C pseudocontact shifts energy
      IF (QENER(SSXPCS)) THEN
      CALL EXPCS(RENR(SSXPCS),'ENERGY ')
      IF (TIMER.GT.1) 
     $     CALL DSPCPU(' ENERGY: after PCS term')
      ENDIF
C
C residual dipolar couplings energy
      IF (QENER(SSXRDC)) THEN
      CALL EXRDC(RENR(SSXRDC),'ENERGY ')
      IF (TIMER.GT.1) 
     $     CALL DSPCPU(' ENERGY: after RDC term')
      ENDIF
C
C cross-correlation rates energy
      IF (QENER(SSXCCR)) THEN
      CALL EXCCR(RENR(SSXCCR),'ENERGY ')
      IF (TIMER.GT.1) 
     $     CALL DSPCPU(' ENERGY: after CCR term')
      ENDIF
C
C RDCANGLES energy
      IF (QENER(SSXANGLE)) 
     $     THEN
      CALL EXANGLES(RENR(SSXANGLE),'ENERGY ')
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after ANG term')
      ENDIF
c     ------------------------------------------------------------------
C
C database energy
C
      IF (QENER(SSRAMA)) THEN
      CALL ERAMA(RENR(SSRAMA), 'ENERGY ')
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after torsion DB term')
      END IF
C
C bond angle database energy
C
      IF (QENER(SSANGDB)) THEN
      CALL EANGLEDB(RENR(SSANGDB), 'ENERGY ')
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after angle DB term')
      END IF
C
C Rgyration energy
C
      IF (QENER(SSCOLL)) THEN
      CALL ECOLLAPSE(RENR(SSCOLL))
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after Rg term')
      END IF

C=====================================================================
C #endif
C=====================================================================
C
      IF (TIMER.EQ.1) CALL DSPCPU(' ENERGY: after ALL terms ')
C
C declare active energy terms
      DO I = 1,NENERT
      IF ( QENER(I) ) THEN
      CALL DECLAR( ANER(I), 'DP', ' ', DCVAL, RENR(I) )
      END IF
      END DO
c     ------------------------------------------------------------------
c     FANTALIN stuff (2003)
c
      fenergtot=total
c     ------------------------------------------------------------------
C
C compute total energy, and put in the accumulation slot SSENER
      TOTAL=ZERO
      DO I=1,NENERT
      IF (QENER(I)) TOTAL=TOTAL+RENR(I)
      END DO
      RENR(SSENER)=TOTAL
      CALL DECLAR( ANER(SSENER), 'DP', ' ', DCVAL, RENR(SSENER) )
C
C compute gradient and put in the accumulation slot SSSD
      TOTAL=ZERO
      NDIM=0
      DO I=1,NATOM
      TOTAL=TOTAL+DX(I)*DX(I)+DY(I)*DY(I)+DZ(I)*DZ(I)
      NDIM=NDIM+3
      END DO
      IF (NDIM.GT.0) TOTAL=SQRT(TOTAL/NDIM)
      RENR(SSSD)=TOTAL
      CALL DECLAR( ANER(SSSD), 'DP', ' ', DCVAL, RENR(SSSD) )
C
      RETURN
      END
C
      SUBROUTINE EFLAGSYMBOLS
C
C Routine exports energy flags to symbols
C
C Author: RWGK
C ============
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'ener.inc'
C
C local
      INTEGER           I, L
      CHARACTER         VARN*10
      DOUBLE PRECISION  DPVAL
      DOUBLE COMPLEX    DCVAL
C
C begin
      DPVAL = DFLOAT(0)
      DCVAL = DCMPLX(DPVAL, DPVAL)
C
      DO I = 1, NENERT
        IF (ANER(I) .NE. ' ') THEN
          VARN = 'EFLAG.' // ANER(I)
          L = LEN(VARN)
          CALL TRIMM(VARN, L)
          IF (QENER(I)) THEN
            CALL DECLAR(VARN(1:L), 'LO', 'TRUE',  DCVAL, DPVAL)
          ELSE
            CALL DECLAR(VARN(1:L), 'LO', 'FALSE', DCVAL, DPVAL)
          END IF
        END IF
      END DO
C
      RETURN
      END
C
      SUBROUTINE ENEINI
C
C Routine initializes energy terms
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'symbol.inc'
C local
      INTEGER I
C begin
C energy term flags and descriptors
      DO I=1,NENERT
      RENR(I)=0.0D0
      ANER(I)='    '
      QENER(I)=.FALSE.
      END DO
C
C accumulation slots, only to inquire, not to manipulate
      ANER(SSTOTE)='TOTE'
      ANER(SSTOTK)='TOTK'
      ANER(SSENER)='ENER'
      ANER(SSTEMP)='TEMP'
      ANER(SSSD)='GRAD'
C
C regular slots
      ANER(SSBOND)='BOND'
      QENER(SSBOND)=.TRUE.
      ANER(SSANGL)='ANGL'
      QENER(SSANGL)=.TRUE.
      ANER(SSDIHE)='DIHE'
      QENER(SSDIHE)=.TRUE.
      ANER(SSIMPR)='IMPR'
      QENER(SSIMPR)=.TRUE.
      ANER(SSVDW)='VDW '
      QENER(SSVDW)=.TRUE.
      ANER(SSELEC)='ELEC'
      QENER(SSELEC)=.TRUE.
      ANER(SSHARM)='HARM'
      ANER(SSZHAR)='ZHAR'
      ANER(SSCDIH)='CDIH'
      ANER(SSNOE)='NOE '
      ANER(SSDG)='DG  '
      ANER(SSXREF)='XREF'
      ANER(SSPVDW)='PVDW'
      ANER(SSPELE)='PELE'
      ANER(SSNCS)='NCS '
      ANER(SSPLN)='PLAN'
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ANER(SSJCOUP)='COUP'
      ANER(SSCARB)='CARB'
      ANER(SSPROT)='PROT'
      ANER(SSONEJ)='ONEB'
      ANER(SSRAMA)='RAMA'
      ANER(SSANGDB)='ANDB'
      ANER(SSVEAN)='VEAN'
      ANER(SSSANI)='SANI'
      ANER(SSDANI)='DANI'
      ANER(SSCOLL)='COLL'
      ANER(SSXPCS)='XPCS'
      ANER(SSXRDC)='XRDC'
      ANER(SSXCCR)='XCCR'
      ANER(SSXANGLE)='XANG'
C=====================================================================
C #endif
C=====================================================================
C harmonic restraints
      PNORMAL(1)=0.0D0
      PNORMAL(2)=0.0D0
      PNORMAL(3)=0.0D0
      QPLANAR=.FALSE.
C
      EOLD=0.0D0
C
      CALL EFLAGSYMBOLS
C
      RETURN
      END
C
      SUBROUTINE EFLAGS(MODE)
C
C Routine sets flags which energy routines are to be used.
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      CHARACTER*(*) MODE
C local
      INTEGER J, NINCL
      LOGICAL TERM, FOUND, MATCH
      CHARACTER*4 BNER(NENERT)
C begin
      TERM=.FALSE.
      CALL PUSEND('FLAGS>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('FLAGS>')
      CALL MISCOM('FLAGS>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-flags')
C
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(A)') ' EFLAGS: the following energy flags are set'
      NINCL=0
      DO J=1,NENERT
      IF (ANER(J).NE.'    '.AND.QENER(J)) THEN
      NINCL=NINCL+1
      BNER(NINCL)=ANER(J)
      END IF
      END DO
      WRITE(6,'(1X,A,14(A,A))') 'EFLAGS:',(' ',BNER(J),J=1,NINCL)
      ELSE IF (WD(1:4).EQ.'INCL') THEN
      TERM=.TRUE.
      ELSE IF (WD(1:4).EQ.'EXCL') THEN
      TERM=.FALSE.
      ELSE
      FOUND=.FALSE.
      DO J=1,NENERT
      CALL EQSTWC(ANER(J),4,WD(1:4),4,1,1,MATCH)
      IF (MATCH) THEN
      QENER(J)=TERM
      FOUND=.TRUE.
      END IF
      END DO
      IF (.NOT.FOUND) CALL CHKEND('FLAGS>',DONE)
      END IF
C
C don't include the accumulation variables
C TOTE, TOTK, TEMP, ENER
      QENER(SSTOTE)=.FALSE.
      QENER(SSTOTK)=.FALSE.
      QENER(SSTEMP)=.FALSE.
      QENER(SSENER)=.FALSE.
      QENER(SSSD)=.FALSE.
C
      END IF
      END DO
      DONE=.FALSE.
C
      CALL EFLAGSYMBOLS
C
      RETURN
      END
C
      SUBROUTINE PGETE(ICYCLE)
C
C Command interpreter for stand-alone and accumulation
C requests for energy
C
C Author: Axel T. Brunger
C =======================
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INTEGER ICYCLE
C local
      INTEGER I
      LOGICAL QPRINT
      CHARACTER*4 SACCU
C
C defaults
      QPRINT=.TRUE.
      SACCU='ENER'
C
C parsing
      CALL PUSEND('ENERGY>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('ENERGY>')
      CALL MISCOM('ENERGY>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-energy')
C
C turn off automatic call to ENERGY in this case
      SACCU='    '
C
      ELSE IF (WD(1:4).EQ.'ACCU') THEN
      CALL NEXTA4('ACCUmulate=',SACCU)
      ELSE
      CALL CHKEND('ENERGY>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (SACCU.EQ.'ENER') THEN
      CALL ENERGY
      ICYCLE=ICYCLE+1
      IF (QPRINT) THEN
      WRITE(6,'(A,I6,A)')
     & ' --------------- cycle=',ICYCLE,
     & ' --------------------------------------------------'
      CALL PRINTE(RENR)
      WRITE(6,'(2A)') ' --------------------------------------------',
     & '-----------------------------------'
      END IF
C
C accumulation stuff
      IF (ICYCLE.EQ.1) THEN
      DO I=1,NENERT
      RENRA(I)=0.0D0
      RENR2A(I)=0.0D0
      END DO
      END IF
      DO I=1,NENERT
      RENRA(I)=RENRA(I)+RENR(I)
      RENR2A(I)=RENR2A(I)+RENR(I)*RENR(I)
      END DO
C
      ELSE IF (SACCU.EQ.'RESE') THEN
      ICYCLE=0
      DO I=1,NENERT
      RENRA(I)=0.0D0
      RENR2A(I)=0.0D0
      END DO
C
      ELSE IF (SACCU.EQ.'PRIN') THEN
      DO I=1,NENERT
      RENRP(I)=RENR(I)
      RENR(I)=RENRA(I)/ICYCLE
      END DO
      WRITE(6,'(A,I8,A)') ' --------------- averages for last ',
     & ICYCLE,'  steps -----------------------------'
      CALL PRINTE(RENR)
      WRITE(6,'(2A)') ' --------------------------------------------',
     &                 '-----------------------------------'
      DO I=1,NENERT
      RENR(I)=RENR2A(I)/ICYCLE-RENR(I)*RENR(I)
      IF (RENR(I).LE.RSMALL) THEN
      RENR(I)=0.0D0
      ELSE
      RENR(I)=SQRT(RENR(I))
      END IF
      END DO
      WRITE(6,'(A,I8,A)') ' --------------- rms fluctations for last ',
     & ICYCLE,' steps -----------------------'
      CALL PRINTE(RENR)
      WRITE(6,'(2A)') ' --------------------------------------------',
     &                 '-----------------------------------'
      DO I=1,NENERT
      RENR(I)=RENRP(I)
      END DO
      END IF
C
      RETURN
      END
C
      SUBROUTINE PRINTE(RENRL)
C
C See routine PRINT2 below
C
C Author: Axel T. Brunger
C =======================
C
C I/O
      DOUBLE PRECISION RENRL(*)
C begin
      CALL PRINT2(.TRUE.,RENRL)
      RETURN
      END
C
      SUBROUTINE PRINT2(QSPCL,RENRL)
C
C Routine prints energies, if QSPCL is true, then this includes the
C Etotal, the grad(E), and the perturbation term.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      LOGICAL QSPCL
      DOUBLE PRECISION RENRL(*)
      INCLUDE 'cns.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'deriv.inc'
C local
      CHARACTER*80 LINE
      CHARACTER*20 ST
      INTEGER NTERMS, J, START, STOP, STLEN
C parameters
      DOUBLE PRECISION MAXNUM
      DOUBLE PRECISION MINNUM
      PARAMETER (MAXNUM=999999.D0)
      PARAMETER (MINNUM=-99999.D0)
C begin
C
C write total energy
      LINE=' '
      LINE(2:2)='|'
      LINE(80:80)='|'
      IF (QSPCL) THEN
      IF (RENRL(SSENER).GT.MAXNUM.OR.RENRL(SSENER).LT.MINNUM) THEN
      WRITE(ST,'(E10.2)') RENRL(SSENER)
      ELSE
      WRITE(ST,'(F10.3)') RENRL(SSENER)
      END IF
      STLEN=10
      CALL TRIML(ST,STLEN)
      WRITE(LINE(3:21),'(4A)') ' Etotal =',ST(1:STLEN)
C
C write the RMS gradient
      IF (RENRL(SSSD).GT.MAXNUM.OR.RENRL(SSSD).LT.MINNUM) THEN
      WRITE(ST,'(E10.2)') RENRL(SSSD)
      ELSE
      WRITE(ST,'(F10.3)') RENRL(SSSD)
      END IF
      STLEN=10
      CALL TRIML(ST,STLEN)
      WRITE(LINE(22:40),'(4A)') ' grad(E)=',ST(1:STLEN)
C
      NTERMS=2
      ELSE
      NTERMS=0
      END IF
C
C now write the partial energy terms
      DO J=1,NENERT
      IF (QENER(J).AND.ANER(J).NE.'    ') THEN
      IF (NTERMS.GE.4) THEN
      WRITE(6,'(A80)') LINE
      LINE=' '
      LINE(2:2)='|'
      LINE(80:80)='|'
      NTERMS=0
      END IF
      NTERMS=NTERMS+1
      IF (RENRL(J).GT.MAXNUM.OR.RENRL(J).LT.MINNUM) THEN
      WRITE(ST,'(E10.2)') RENRL(J)
      ELSE
      WRITE(ST,'(F10.3)') RENRL(J)
      END IF
      STLEN=10
      CALL TRIML(ST,STLEN)
      START=(NTERMS-1)*19 +3
      STOP=NTERMS*19+2
      WRITE(LINE(START:STOP),'(4A)') ' E(',ANER(J),')=',ST(1:STLEN)
      END IF
      END DO
      IF (NTERMS.GE.1) WRITE(6,'(A80)') LINE
C
      RETURN
      END
C=====================================================================
      SUBROUTINE ENERG2(INTERE,N)
C
C Front-end routine for various partial energy routines
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'param.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INTEGER N, INTERE(*)
C local
      INTEGER I, NNPIG
      LOGICAL ERR
C parameters
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C initialization
      DO I=1,NENERT
      REN2(I)=ZERO
      RENV(I)=ZERO
      ENDDO
C
C check that all atomic coordinates are known
      CALL ATMCHK(HEAP(IINTER(N)),ERR)
      IF (.NOT.ERR) THEN
C
      IF (NBOND.GT.0.AND.QENER(SSBOND).AND.
     &    (ABS(PIGWGHT(N,SSBOND)).GE.RSMALL.OR.
     &     ABS(PIGAVWT(N,SSBOND)).GE.RSMALL)) THEN
      CALL EBOND(REN2(SSBOND),RENV(SSBOND),IB,JB,INTERE,NBOND,
     &  ACBC,ACBB,QACBC,QACBB,'NOAN',
     &  PIGWGHT(N,SSBOND),PIGAVWT(N,SSBOND),IAC,SEGID,RESID,TYPE)
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after BOND term ')
      END IF
C
      IF (NTHETA.GT.0.AND.QENER(SSANGL).AND.
     &    (ABS(PIGWGHT(N,SSANGL)).GE.RSMALL.OR.
     &     ABS(PIGAVWT(N,SSANGL)).GE.RSMALL)) THEN
      CALL EANGLE(REN2(SSANGL),RENV(SSANGL),IT,JT,KT,INTERE,
     &    NTHETA,ACTC,ACTB,ACTUC,ACTUB,QACTC,QACTB,QACTUC,QACTUB,
     &    'NOAN',PIGWGHT(N,SSANGL),PIGAVWT(N,SSANGL),
     &    IAC,SEGID,RESID,TYPE)
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after ANGLE term ')
      END IF
C
      IF (NPHI.GT.0.AND.QENER(SSDIHE).AND.
     &    (ABS(PIGWGHT(N,SSDIHE)).GE.RSMALL.OR.
     &     ABS(PIGAVWT(N,SSDIHE)).GE.RSMALL)) THEN
      CALL EDIHE(REN2(SSDIHE),RENV(SSDIHE),IP,JP,KP,LP,INTERE,
     &      NPHI,ACPC,ACPD,ACPB,
     &      QACPC,QACPD,QACPB,'NORMAL',0.0D0,0,'NOAN',
     &      PIGWGHT(N,SSDIHE),PIGAVWT(N,SSDIHE),IAC,SEGID,RESID,TYPE)
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after DIHEdral term ')
      END IF
C
      IF (NIMPHI.GT.0.AND.QENER(SSIMPR).AND.
     &    (ABS(PIGWGHT(N,SSIMPR)).GE.RSMALL.OR.
     &     ABS(PIGAVWT(N,SSIMPR)).GE.RSMALL)) THEN
      CALL EIMPR(REN2(SSIMPR),RENV(SSIMPR),IM,JM,KM,LM,INTERE,
     &    NIMPHI,ACIC,ACID,ACIB,
     &    QACIC,QACID,QACIB,'NORMAL',0.0D0,0,'NOAN',
     &    PIGWGHT(N,SSIMPR),PIGAVWT(N,SSIMPR),IAC,SEGID,RESID,TYPE)
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after IMPRoper term ')
      END IF
C
      IF (QENER(SSELEC).OR.QENER(SSVDW)  .OR.
     &    QENER(SSPELE).OR.QENER(SSPVDW)) THEN
C
C determine the last PIG that is turned on.
      NNPIG=0
      DO I=1,NPIG
      IF (ABS(PIGWGHT(I,SSVDW)) .GE.RSMALL.OR.
     &    ABS(PIGWGHT(I,SSELEC)).GE.RSMALL.OR.
     &    ABS(PIGAVWT(I,SSVDW)) .GE.RSMALL.OR.
     &    ABS(PIGAVWT(I,SSELEC)).GE.RSMALL) THEN
      NNPIG=I
      END IF
      END DO
C
      IF (ABS(PIGWGHT(N,SSVDW)) .GE.RSMALL.OR.
     &    ABS(PIGWGHT(N,SSELEC)).GE.RSMALL.OR.
     &    ABS(PIGAVWT(N,SSVDW)) .GE.RSMALL.OR.
     &    ABS(PIGAVWT(N,SSELEC)).GE.RSMALL) THEN
      CALL ENBOND(REN2(SSVDW),REN2(SSELEC),REN2(SSPVDW),REN2(SSPELE),
     &            RENV(SSVDW),RENV(SSELEC),RENV(SSPVDW),RENV(SSPELE),N,
     &            NNPIG,'NOAN')
      IF (TIMER.GT.1) CALL DSPCPU(' ENERGY: after NONBonded term ')
      END IF
      END IF
C
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE ATMCHK(INTERE,ERR)
C
C subroutine checks whether all atoms are known that are
C referenced in the interaction energy array.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INTEGER INTERE(*)
      LOGICAL ERR
C local
      INTEGER I
C begin
      ERR=.FALSE.
      DO I=1,NATOM
      IF (.NOT.INITIA(I,X,Y,Z).AND.INTERE(I).LT.9999) THEN
      ERR=.TRUE.
      WRITE(6,'(9A)')
     &  ' %ATMCHK-ERR: unknown coordinates for atom "',
     &    SEGID(I),'-',RESID(I),'-',RES(I),'-',TYPE(I),'"'
      END IF
      END DO
      IF (ERR) THEN
      CALL WRNDIE(-5,'ATMCHK','Unknown coordinates')
      END IF
      RETURN
      END

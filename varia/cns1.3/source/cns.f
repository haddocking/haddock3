      PROGRAM CNS
C
C  Crystallography & NMR System (CNS)
C
C  A program for X-ray crystallographic and solution NMR structure
C  determination.
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=+
C  +===============================================================+
C  + Authors: A.T.Brunger, P.D.Adams, G.M.Clore, W.L.DeLano,       +
C  +          P.Gros, R.W.Grosse-Kunstleve, J.-S.Jiang, J.M.Krahn  +
C  +          J.Kuszewski, M.Nilges, N.S.Pannu, R.J.Read,          +
C  +          L.M.Rice, G.F.Schroeder, T.Simonson, G.L.Warren      +
C  +===============================================================+
C  + Copyright 1997-2010 Yale University                           +
C  +===============================================================+
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IMPLICIT NONE
C
C common (global) variables and constants
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'param.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'seed.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'update.inc'
      INCLUDE 'version.inc'
C =====
C local
C =====
      INTEGER ICYCLE, I, PROMPTLEN, IMRT
      LOGICAL QOPEN, QFORM, QWRITE
      DOUBLE PRECISION DUMMY
      DOUBLE COMPLEX CDUMMY
      CHARACTER*12 USERNM
      CHARACTER*14 SYSNM
      CHARACTER*32 HOSTNM
      CHARACTER*11  DATE
      CHARACTER*8  TIME
      CHARACTER*20 PROMPT
C
C =====
C begin
C =====
C
C initialize the file handling facility
      CALL VOPEN(0,'$$INITIAL$$',' ',' ',ERROR)
C
C CHARACTER*4 array CSTACK initialization
      CALL CINIST(.FALSE.)
C
C HEAP initialization and sanity check
      CALL INITHP
      CALL HEAPVFY(ERROR)
      IF (ERROR) THEN
        CALL WRNDIE(-5, 'HEAPVFY', 'FATAL INTERNAL ERROR.')
        CALL DIE
      END IF
      CALL INITHP
C
C initialize command parsing
      CALL INIPRS
      CALL DECLAR( 'CNS_VERSION', 'ST', CNSVERSION, DUMMY, CDUMMY )
      CALL DECLAR( 'CNS_PATCH', 'ST', CNSPATCH, DUMMY, CDUMMY )
C
C initialize machine dependent variables
      CALL SETFPEPS
C
C initialize the interacting group setup
      CALL PIGRSET('INIT')
C
C initialize the nonbonded list
      CALL NBDSET('INIT')
C
C initialize the fast nonbonded setup
      CALL FNBSET('INIT')
C
C initialize the trajectory writing facility
      CALL TRTINI
C
C initialize the parameter learning facility
      CALL PRTINI
C
C torsion angle dynamics
      CALL TORINI
C
C initialize the distance geometry bounds matrix stores
      CALL DGINIT
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
C initialize j coupling refinement
      CALL COUPINIT
C
C initialize carbon chemical shift stuff
      CALL CSHIFTINIT
C
C initialize the proton chemical shift stuff
      CALL HSHIFTINIT
C
C initialize the 1-bond J coupling stuff
      CALL ONEINIT
C
C initialize the Ramachandran stuff
      CALL RAMAINIT
C
C initialize the angle database stuff
      CALL ANGLEDBINIT
C
C initial susceptibility anisotropy refinement
      CALL VEANINIT
C
C initial susceptibility anisotropy refinement
      CALL SANIINIT
C
C initial diffusion anisotropy refinement
      CALL DANIINIT
C
C PARA restraints  -----------------------------------------------------
C
C initialize Pseudocontact shifts refinement
      CALL XPCSINIT
C
C initialize Residual dipolar couplings refinement
      CALL XRDCINIT
C
C initialize Cross-correlation rates refinement
      CALL XCCRINIT
C
C initialize RDC-ANGLES refinement
      CALL XANGLESINIT
C     ------------------------------------------------------------------
C
C rotation search clustering
      CALL CLUSINI
C
C mapyard module
      CALL INIMYCOMMON(.TRUE.)
C
C psearch module
      CALL INIPSTCOMMON(.TRUE.)
c
c initialize Rg term
      CALL COLLAPSEINIT
C=====================================================================
C #endif
C=====================================================================
C set bond, angle, dihedral, improper type-based parameter update flags
      UPBOND=.TRUE.
      UPANGL=.TRUE.
      UPDIHE=.TRUE.
      UPIMPR=.TRUE.
C
C initialize scalar stores
      DO I=1,MAXSTO
      PTRSTO(I)=0
      LENSTO(I)=0
      END DO
C
C constraints
      CALL CSTRES
C
C dihedral angle restraints
      CALL CDIINI
C
C NOE tables
      CALL NOEINI
C
C crystallographic restraints
      CALL XRERES
C
C energy term flags and descriptors
      CALL ENEINI
C
C non-crystallographic symmetry operators and counters
      CALL NCSINR
      CALL NCSINS
C
C planarity restraints
      CALL PLNINR
C
C set counter for some ENERGY calls
      ICYCLE=0
C
C RTF
      CALL RTFINI
C
C MTF
      CALL MTFRES
C
C parameters
      CALL PARINI
C
C chebychev poynomials
      CALL CHEINI
C
C message level and program abort level
      WRNLEV=5
      BOMLEV=0
C
C random number generator seed
      SEED=314159.0D0
C
C display and print files are set to system output (default)
      DUNIT=6
      PUNIT=6
C
C journal file unit is set to -1 (switched off)
      PRUNIT=-1
C
C initialize timer
      CALL VINTIM
      TIMER=0
C
C initialize prompt
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      PROMPT='CNSsolve>'
C=====================================================================
C #else
C=====================================================================
C      PROMPT='CNSrefine>'
C=====================================================================
C #endif
C=====================================================================
      PROMPTLEN=20
      CALL TRIMM(PROMPT,PROMPTLEN)
C
C write header
      CALL CNSHEADER(DATE,TIME,USERNM,SYSNM,HOSTNM)
C
C Verify FFT: Make sure FFT library actually works
      CALL VFYFFT
C
C Verify numerical code: Make sure a representative routine 
C                        functions correctly
      CALL VFYNUM
C
      WRITE(6,'(1X)')
C
C command parsing
      CALL PUSEND(PROMPT(1:PROMPTLEN))
      DONE=.FALSE.
      DO WHILE (.NOT.DONE)
      CALL NEXTWD(PROMPT(1:PROMPTLEN))
      CALL MISCOM(PROMPT(1:PROMPTLEN),USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns')
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ANGL' ) THEN
      CALL READANGLEDB
      ELSE IF (WD(1:4).EQ.'CARB' ) THEN
      CALL READCSHIFT
      ELSE IF ((WD(1:4).EQ.'CONF' ).OR.(WD(1:4).EQ.'RAMA')) THEN
      CALL READRAMA
      ELSE IF (WD(1:4).EQ.'COUP' ) THEN
      CALL READCOUP
      ELSE IF (WD(1:4).EQ.'DANI' ) THEN
      CALL READDANI
      ELSE IF (WD(1:4).EQ.'ONEB' ) THEN
      CALL READONEBOND
      ELSE IF (WD(1:4).EQ.'SANI' ) THEN
      CALL READSANI
      ELSE IF (WD(1:4).EQ.'VEAN' ) THEN
      CALL READVEAN
      ELSE IF (WD(1:4).EQ.'PROT' ) THEN
      CALL READHSHIFT
      ELSE IF (WD(1:4).EQ.'MIND' ) THEN
      CALL MINDWR
      ELSE IF (WD(1:4).EQ.'ARIA' ) THEN
      CALL ARISET
c
c     Commands for PARArestraints (2003)
c
      ELSE IF (WD(1:4).EQ.'XPCS' ) THEN
      CALL READXPCS
      ELSE IF (WD(1:4).EQ.'XRDC' ) THEN
      CALL READXRDC
      ELSE IF (WD(1:4).EQ.'XCCR' ) THEN
      CALL READXCCR
      ELSE IF (WD(1:4).EQ.'XANG' ) THEN
      CALL READXANGLES
      ELSE IF (WD(1:4).EQ.'XT1D' ) THEN
      CALL XT1DIST
C=====================================================================
C #endif
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ARIA' ) THEN
      CALL ARISET
      ELSE IF (WD(1:4).EQ.'COLL' ) THEN
      CALL READCOLLAPSE
      ELSE IF (WD(1:4).EQ.'CONN' ) THEN
      CALL CONNEC
      ELSE IF (WD(1:4).EQ.'COOR' ) THEN
      CALL CORMAN
      ELSE IF (WD(1:4).EQ.'DELE' ) THEN
      CALL DELTIC
      ELSE IF (WD(1:4).EQ.'DIST' ) THEN
      CALL DISTAN
      ELSE IF (WD(1:4).EQ.'DO' ) THEN
      CALL SCDO
      ELSE IF (WD(1:4).EQ.'DUPL' ) THEN
      CALL DUPLIC
      ELSE IF (WD(1:4).EQ.'DYNA' ) THEN
      CALL DYNMCS
      ELSE IF (WD(1:4).EQ.'ENER' ) THEN
      CALL PGETE(ICYCLE)
      ELSE IF (WD(1:4).EQ.'FAST' ) THEN
      CALL FNBSET('PARSE')
      ELSE IF (WD(1:3).EQ.'FIX' ) THEN
      CALL FIX
      ELSE IF (WD(1:4).EQ.'FLAG' ) THEN
      CALL EFLAGS('PARSE')
      ELSE IF (WD(1:4).EQ.'IDEN' ) THEN
      CALL SCIDEN
      ELSE IF (WD(1:4).EQ.'IGRO' ) THEN
      CALL PIG
      ELSE IF (WD(1:4).EQ.'MINI' ) THEN
      CALL MINIMI
      ELSE IF (WD(1:4).EQ.'MMDG') THEN
      CALL MMDG
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      CALL MODES
      ELSE IF (WD(1:3).EQ.'NCS') THEN
      CALL NCSPAR
      ELSE IF (WD(1:4).EQ.'NOE ' ) THEN
      CALL NOESET
      ELSE IF (WD(1:4).EQ.'PARA' ) THEN
      CALL PARRDR
      ELSE IF (WD(1:4).EQ.'PATC' ) THEN
      CALL PATSET
      ELSE IF (WD(1:4).EQ.'PICK' ) THEN
      CALL PICK(.FALSE.)
      ELSE IF (WD(1:4).EQ.'PRIN' ) THEN
      CALL PRINT
      ELSE IF (WD(1:4).EQ.'READ' ) THEN
      CALL READD
      ELSE IF (WD(1:4).EQ.'REST' ) THEN
      CALL RSTRAN
      ELSE IF (WD(1:4).EQ.'ROTM' ) THEN
      CALL ROTMAN
      ELSE IF (WD(1:4).EQ.'SEGM' ) THEN
      CALL SEGMNT
      ELSE IF (WD(1:4).EQ.'SHOW' ) THEN
      CALL SCSHOW
      ELSE IF (WD(1:4).EQ.'STRU' ) THEN
      CALL MTFPRS
      ELSE IF (WD(1:4).EQ.'SURF' ) THEN
      CALL SURFAC
      ELSE IF (WD(1:4).EQ.'TEST' ) THEN
      CALL TESTCNS
      ELSE IF (WD(1:4).EQ.'TOPO' ) THEN
      CALL RTFRDR
      ELSE IF (WD(1:4).EQ.'WRIT' ) THEN
      CALL WRITT
      ELSE IF (WD(1:4).EQ.'XRAY' ) THEN
      CALL XRAY
      ELSE IF (WD(1:4).EQ.'?   ' ) THEN
      CALL SHOWAL
      ELSE IF (WD(1:4).EQ.'STOP'.OR.EOF ) THEN
      DONE=.TRUE.
      ENDIND=ENDIND-1
      ELSE
      CALL CHKEND(PROMPT(1:PROMPTLEN),DONE)
      IF (DONE) THEN
      WRITE(6,'(A)')
     &' %CNS-ERR: "END" not allowed. "STOP" terminates program'
      DONE=.FALSE.
      ENDIND=ENDIND+1
      END IF
      END IF
      END IF
      END DO
C
C terminate program
C
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
C free proton chemical shift restraints list
      CALL HSHIFTSHP
C free j-coupling restraints list
      CALL COUPHP
C free carbon chemical shift restraints list
      CALL CARBHP
C free onebond j-coupling restraints list
      CALL ONEHP
C free conf. torsion angle database restraints list
      CALL RAMAHP
C free diffusion anisotropy restraints list
      CALL DANIHP
C free susceptibility anisotropy restraints list
      CALL SANIHP
C free susceptibility anisotropy restraints list
      CALL VEANHP
C free conf. angle database restraints list
      CALL ANGLEDBFREE
C rotation search clustering
      CALL CLUSFRE
C=====================================================================
C #endif
C=====================================================================
C free dihedral angle restraints list
      CALL CDIHP(0)
C free NOE restraints list
      CALL NOEHP(0)
C free XRAY atom lists, reflections, and maps
      CALL XRAFRE
      CALL XRAATM(0)
      CALL XRMAPR(-1)
C
C free NCS group lists
      CALL NCSFIN
C
C free PLN group lists
      CALL PLNFIN
C
C free buffer stores
      DO I=1,MAXSTO
      IF (PTRSTO(I).NE.0) CALL FREHP(PTRSTO(I),IREAL8(LENSTO(I)))
      END DO
C
C free nonbonded lists
      DO I=1,NPIGMAX
      CALL RENBND(0,0,0,0,0,0,0,0,0,I)
      ENDDO
C
C free PIG arrays
      CALL PIGRSET('FREE')
C
C free up space for the trajectory writing facility
      CALL TRTFRE
C
C free up space for the parameter learning facility
      CALL PRTFRE
C
C free up torsion angle stuff
      CALL TORHP(0,.TRUE.)
C
C free up space for the distance geometry bounds matrices
      CALL DGRESZ(0)
C
C write info about STACK and HEAP
      CALL PRISTK
      CALL PRINHP
      CALL PRIEND
C
C close all files
      DO I=1,99
      IMRT=I
      CALL VINQRE('UNIT',WDT,WDTMAX,WDTLEN,QOPEN,QFORM,QWRITE,IMRT)
      IF (QOPEN) CALL VCLOSE(IMRT,'KEEP',ERROR)
      END DO
C
C write info about maximum dynamic memory usage and total CPU time
      CALL PRFINAL
C
      END
C
C======================================================================
C
      SUBROUTINE CNSHEADER(DATE,TIME,USERNM,SYSNM,HOSTNM)
C
      IMPLICIT NONE
C
      INCLUDE 'version.inc'
C I/O
      CHARACTER*11  DATE
      CHARACTER*8  TIME
      CHARACTER*12 USERNM
      CHARACTER*14 SYSNM
      CHARACTER*32 HOSTNM
C local
      INTEGER      HNLEN, TMP, PTRSZ, STLEN
      CHARACTER*(1) CNSPTMP
      CHARACTER*(4) ST
!$    integer omp_get_max_threads
C
C write header
      WRITE(6,'(10X,A)')
     &'============================================================'
      WRITE(6,'(10X,A,58X,A)') '|','|'
      WRITE(6,'(10X,A,12X,A,12X,A)')
     &'|','Crystallography & NMR System (CNS)','|'
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      WRITE(6,'(10X,A,25X,A,25X,A)')
     &'|','CNSsolve','|'
C=====================================================================
C #else
C=====================================================================
C      WRITE(6,'(10X,A,24X,A,25X,A)')
C     &'|','CNSrefine','|'
C=====================================================================
C #endif
C=====================================================================
      WRITE(6,'(10X,A,58X,A)') '|','|'
      WRITE(6,'(10X,A)')
     &'============================================================'
      CNSPTMP=CNSPATCH
      IF (CNSPTMP.EQ.'0') THEN
      WRITE(6,'(10X,2A)')
     &' Version: ',CNSVERSION
      ELSE
      WRITE(6,'(10X,4A)')
     &' Version: ',CNSVERSION,' at patch level ',CNSPATCH
      END IF
      WRITE(6,'(10X,A)')
     &' Status: Special UU release with Rg, paramagnetic'
      WRITE(6,'(10X,A)')
     &'         and Z-restraints (A. Bonvin, UU 2013)'
      WRITE(6,'(10X,A)')
     &'============================================================'
      WRITE(6,'(10X,A,/,10X,A,/,10X,A,/,10X,A)')
     &' Written by: A.T.Brunger, P.D.Adams, G.M.Clore, W.L.DeLano,',
     &'             P.Gros, R.W.Grosse-Kunstleve,J.-S.Jiang,J.M.Krahn,',
     &'             J.Kuszewski, M.Nilges, N.S.Pannu, R.J.Read,',
     &'             L.M.Rice, G.F.Schroeder, T.Simonson, G.L.Warren.'
      WRITE(6,'(10X,A)')
     &' Copyright (c) 1997-2010 Yale University'
      WRITE(6,'(10X,A)')
     &'============================================================'
C
      CALL VHOSTNAME(HOSTNM, HNLEN)
      CALL GETSYS(SYSNM,14,TMP)
C
      CALL CNSQPTRSZ(PTRSZ)
      PTRSZ=PTRSZ*8
      STLEN=4
      CALL ENCODI(PTRSZ,ST,4,STLEN)
C
      WRITE(6,'(10X,7A)')
     & ' Running on machine: ',HOSTNM(1:HNLEN),
     & ' (',SYSNM(1:TMP),',',ST(1:STLEN),'-bit)'
!$    WRITE(6,'(31X,A,I3,A)') 'with',OMP_GET_MAX_THREADS(),' threads'
      CALL GETNAM(USERNM,12,TMP)
      WRITE(6,'(10X,2A)')
     & ' Program started by: ',USERNM(1:TMP)
C
      CALL VTIME(TIME,8,TMP)
      CALL VDATE(DATE,11,TMP)
      WRITE(6,'(10X,4A)')
     & ' Program started at: ',TIME,' on ',DATE(1:TMP)
C
      WRITE(6,'(10X,A)')
     &'============================================================'
C
      WRITE(6,'(1X)')
C
      RETURN
      END
C
C======================================================================
C

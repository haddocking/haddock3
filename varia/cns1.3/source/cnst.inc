C data structure CNST
C ===================
C
C harmonic restraints, ic restraints, friction coefficients,
C dihedral restraints
C Pairs of Interacting Groups (PIG) information
C
C MAXA2 must match MAXA
C
      INTEGER   MAXA2
      PARAMETER (MAXA2=MAX_ATOM)
C
C
C  HARMONIC RESTRAINTS
C     KCNSTR     Force constant of atomic restraints
C     KZCNSTR    Force constant of atomic z-restraints
C     KCEXPN     Exponent of the restraint (default 2)
C     REFX       X component of reference coordinates
C     REFY       Y component of reference coordinates
C     REFZ       Z component of reference coordinates
      INTEGER KCEXPN
      DOUBLE PRECISION KCNSTR(MAXA2), REFX(MAXA2), REFY(MAXA2)
      DOUBLE PRECISION REFZ(MAXA2), KZCNSTR(MAXA2)
C
C     PNORMAL    3-d vector used for the definition of planar restrains
C     QPLANAR    logical variable that is true, if there are planar
C                constrains defined by NORMAL
      DOUBLE PRECISION PNORMAL(3)
      LOGICAL QPLANAR
C
C  RMSD array
      DOUBLE PRECISION RMSD(MAXA2)
C
C old garbage
      DOUBLE PRECISION CCBIC,CCTIC,CCPIC
C
C  Z-HARMONIC RESTRAINTS
C     ZHARMAX	Upper Z limit
C     ZHARMIN   Lower Z limit
C     KZCEPX    Exponent of the restraint (default 2)
      INTEGER KZCEXP
      DOUBLE PRECISION ZHARMAX, ZHARMIN
C
C
C  DIHEDRAL ANGLE RESTRAINTS
C     CIP   atom i
C     CJP   atom j
C     CKP   atom k
C     CLP   atom l
C     CICP  dummy lookup for torsion energy routine (has to be ID)
C     CNPHI dumber of dihedral angle restraints
C     MCMPHI maximum number of dihedral angle restraints
C     CCPB  equilibrium value for the dihedral energy
C     CCPC  Force constant for the dihedral energy
C     CCPO  offset for the dihedral energy
C     CCPE  exponent for the dihedral energy
C     CCPD  dummy (has to be zero)
C     CNSCA  overall scale factor
      INTEGER MCMPHI, CNPHI
      DOUBLE PRECISION CNSCA
C heap pointers
      INTEGER ICIP, ICJP, ICKP, ICLP
      INTEGER ICICP, ICCPD, ICCPE, ICCPB, ICCPC, ICCPO
      INTEGER LCICP
C
C  DIHEDRAL ANGLE RESTRAINTS working set for cross-validation
C  A. Bonvin 12/6/95
C
      INTEGER CWNPHI
C heap pointers
      INTEGER IWCIP, IWCJP, IWCKP, IWCLP
      INTEGER IWCICP, IWCCPD, IWCCPE, IWCCPB, IWCCPC, IWCCPO
C
C cross-validation array
      INTEGER HPDCV
C
C cross-validation test number
      INTEGER DIHICV
C
C
C  FRICTION COEFFICIENTS
C
C     FBETA      Langevin dynamics friction coefficient
      DOUBLE PRECISION FBETA(MAXA2)
C
C
C  INTERACTION ARRAY
C
C     IINTER (points to the INTERE arrays on the HEAP)
C     NPIGMAX maximum number of Pairs of Interacting Groups
C    (must match ener.inc, mtf.inc and nbonds.inc)
      INTEGER NPIGMAX
      PARAMETER (NPIGMAX=MAX_PIGS)
      INTEGER IINTER(NPIGMAX)
C NPIG number of Pairs of Interacting Groups (PIGs)
      INTEGER NPIG
C
C QPIGRST: flag indicating that interacting pairs have been set
      LOGICAL QPIGRST
C dimension of INTERE on heap :
      INTEGER XINTER
C
C COMMON BLOCKS
C +++++++++++++
C
C     double precision block
C
      COMMON /CNSTD/
     &     KCNSTR,  KZCNSTR,  REFX,    REFY,    REFZ,     CCBIC,
     &     CCTIC,   CCPIC,   FBETA,   CNSCA,   RMSD, PNORMAL,
     &     ZHARMAX, ZHARMIN
C
C     integer blocks
C
      COMMON /CNSTI/
     &     KCEXPN,  IINTER,  ICIP,     ICJP,        ICKP,
     &     ICLP,     ICICP,    CNPHI,   ICCPD,     ICCPE,
     &     NPIG,  XINTER,    ICCPB,    ICCPC,    ICCPO, MCMPHI,
     &     LCICP,    HPDCV,   DIHICV,
     &     IWCIP,    IWCJP,    IWCKP,    IWCLP,
     &     IWCICP,  IWCCPD,   IWCCPE,   IWCCPB,   IWCCPC,
     &     CWNPHI,  IWCCPO,
     &     KZCEXP
C
C     logical common block
      COMMON /CNSTL/ QPLANAR, QPIGRST
C
      SAVE /CNSTD/
      SAVE /CNSTI/
      SAVE /CNSTL/

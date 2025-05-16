C param.inc
C =========
C
C Purpose:
C
C To store the parameters used in evaluating the energy of
C the system.
C
C
C bonds
C  MAXCB Maximum number of bond parameters
C  NCB  Number of bond parameters
C  CBB  Equilibrium bond distance
C  CBC  Force constant of a bond
C  KCB1 chemical type of first bond atom
C  KCB2 chemical type of second bond atom
      INTEGER MAXCB
      PARAMETER (MAXCB=2000)
      INTEGER NCB
      DOUBLE PRECISION CBB(MAXCB), CBC(MAXCB)
      CHARACTER*4 KCB1(MAXCB), KCB2(MAXCB)
C
C angles
C  MAXCT Maximum number of bond angles parameters
C  NCT  Number of bond angles parameters
C  CTB  Equilibrium bond angle
C  CTC  Force constant for a bond angle
C  CTUB  Equilibrium Urey-Bradley distance
C  CTUC  Force constant for Urey-Bradley term
C  KCT1 chemical type of first angle atom
C  KCT2 chemical type of middle angle atom
C  KCT3 chemical type of third angle atom.
      INTEGER MAXCT
      PARAMETER (MAXCT=20000)
      INTEGER NCT
      DOUBLE PRECISION CTB(MAXCT), CTC(MAXCT), CTUB(MAXCT), CTUC(MAXCT)
      CHARACTER*4 KCT1(MAXCT), KCT2(MAXCT), KCT3(MAXCT)
C
C dihedrals
C  MAXCP Maximum number of dihedral parameters
C  NNCCP  Number of dihedral parameters
C  CPB  Phase shift (delta) for the dihedral energy
C  CPC  Force constant for the dihedral energy
C  CPD  Periodicity of the dihedral energy
C  KCP1 chemical type of first dihedral atom
C  KCP2 chemical type of second dihedral atom
C  KCP3 chemical type of third dihedral atom
C  KCP4 chemical type of fourth dihedral atom
      INTEGER MAXCP
      PARAMETER (MAXCP=100000)
      INTEGER NNCCP
      DOUBLE PRECISION CPB(MAXCP), CPC(MAXCP)
      INTEGER CPD(MAXCP)
      CHARACTER*4 KCP1(MAXCP), KCP2(MAXCP), KCP3(MAXCP), KCP4(MAXCP)
C
C impropers
C  MAXCI Maximum number of improper torsion parameters
C  NCI Number of improper torsion parameters
C  CIB Equilibrium improper torsion angle
C  CIC Force constant for the improper torsion
C  CID Periodicity of the improper energy
C  KCI1 chemical type of first improper atom
C  KCI2 chemical type of second improper atom
C  KCI3 chemical type of third improper atom
C  KCI4 chemical type of fourth improper atom
      INTEGER MAXCI
      PARAMETER (MAXCI=2000)
      INTEGER NCI
      DOUBLE PRECISION CIB(MAXCI), CIC(MAXCI)
      INTEGER CID(MAXCI)
      CHARACTER*4 KCI1(MAXCI), KCI2(MAXCI), KCI3(MAXCI), KCI4(MAXCI)
C
C lennard-jones parameters
C  MAXCN  Maximum number of non-bonding atom types.
C  NCN    Number of non-bonded atom types.
C  CNAC   chemical type of nonbonding atom
C  CNBA   "A" Lennard-Jones parameter (matrix)
C  CNBB   "B" Lennard-Jones parameter (matrix)
C  CNBA14 "A" 1-4 interaction Lennard-Jones parameter (matrix)
C  CNBB14 "B" 1-4 interaction Lennard-Jones parameter (matrix)
C  CNBVR  vdW radii for hard sphere repel potential
C  CBVR14 same for 1-4 terms
C  NLJAT  Number of atoms with atom-based Lennard-Jones parameters
      INTEGER MAXCN
      PARAMETER (MAXCN=1500)
      INTEGER NCN, NLJAT
      CHARACTER*4 CNAC(MAXCN)
      DOUBLE PRECISION CNBA(MAXCN,MAXCN), CNBB(MAXCN,MAXCN)
      DOUBLE PRECISION CNBA14(MAXCN,MAXCN), CNBB14(MAXCN,MAXCN)
      DOUBLE PRECISION CNBVR(MAXCN,MAXCN), CBVR14(MAXCN,MAXCN)
C
C
C COMMON BLOCKS
C =============
C
C     double precision block
C
      COMMON/PARA/ CBB,   CBC,    CTB,    CTC,  CTUB, CTUC, CPB, CPC,
     3             CIB,   CIC,    CNBA,  CNBB,   CNBA14, CNBB14,
     6             CNBVR, CBVR14
C
C     integer block
C
      COMMON/IPARA/ NCB, NCT, NNCCP,    CPD,
     &             NCI, CID, NCN, NLJAT
C
C     character string block
C
      COMMON/CPARA/ KCB1, KCB2,
     1             KCT1, KCT2, KCT3,
     2             KCP1, KCP2, KCP3, KCP4,
     3             KCI1, KCI2, KCI3, KCI4,
     4             CNAC
C
      SAVE /PARA/
      SAVE /IPARA/
      SAVE /CPARA/

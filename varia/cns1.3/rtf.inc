C rtf.inc
C residue topologies
C
C mass defaults
C  MAXAT2 maximum number of mass default specifications
C  IATC   number of mass default specifications
C  ATCT   chemical atom type for mass default
C  MASDEF value for mass default specification
      INTEGER MAXAT2
      PARAMETER (MAXAT2=1000)
      INTEGER IATC
      DOUBLE PRECISION MASDEF(MAXAT2)
      CHARACTER*4 ATCT(MAXAT2)
C
C RESIdues
C     MXRTRS  maximum number of residues in file
C     NRTRS   number of residues in file
C     NICM    dimension of linked list NIC
C     NIC     linked list for bonds, angles,...
C     RTRTYP  0 if normal residue, 1 if patch residue
C
      INTEGER MXRTRS, NICM
      PARAMETER (MXRTRS=50000,NICM=50)
      INTEGER NIC(NICM,MXRTRS)
      INTEGER NRTRS,RTRTYP(MXRTRS)
      CHARACTER*4 AA(MXRTRS)
C
C ATOMs
C     MXRTA    maximum number of atoms in topology file
C     FTP(*)   atom name
C     CHG(*)   charge
C     ARMASS(*) mass
C     MAC(*)   type
C     GRPR(*)  group index
C     DELAT(*) atom option (e.g. DELE, MODI, ADD)
C
      INTEGER MXRTA
      PARAMETER (MXRTA=5000)
      CHARACTER*5 FTP(MXRTA)
      CHARACTER*4 DELAT(MXRTA), MAC(MXRTA)
      DOUBLE PRECISION CHG(MXRTA), ARMASS(MXRTA)
      INTEGER GRPR(MXRTA)
C
C explicit nonbonding exclusions
C     MXRTX maximum number of exclusions
C     MXN(*) index into linked list MNB
C     MNB(*) linked list
C
      INTEGER MXRTX
      PARAMETER (MXRTX=2000)
      INTEGER MXN(MXRTA)
      CHARACTER*5 MNB(MXRTX)
C
C BONDs
C     MXRTB maximum number of bonds
C     MIB(*) first bond atom
C     MJB(*) second bond atom
C     DELBD(*) bond option (e.g. DELE, ADD )
C
      INTEGER MXRTB
      PARAMETER (MXRTB=4000)
      CHARACTER*5 MIB(MXRTB), MJB(MXRTB)
      CHARACTER*4 DELBD(MXRTB)
C
C ANGLes
C     MXRTT maximum number of angles
C     MIT(*) first angle atom
C     MJT(*) second angle atom
C     MKT(*) third angle atom
C     DELAN(*) angle option
C
      INTEGER MXRTT
      PARAMETER (MXRTT=8000)
      CHARACTER*5 MIT(MXRTT), MJT(MXRTT), MKT(MXRTT)
      CHARACTER*4 DELAN(MXRTT)
C
C DIHEdrals
C     MXRTP maximum number of dihedrals in topology file
C     MIP(*) first dihedral atom
C     MJP(*) second dihedral atom
C     MKP(*) third dihedral atom
C     MLP(*) fourth dihedral atom
C     DELPT(*) dihedral option
C
      INTEGER MXRTP
      PARAMETER (MXRTP=20000)
      CHARACTER*5 MIP(MXRTP),MJP(MXRTP),MKP(MXRTP),MLP(MXRTP)
      CHARACTER*4 DELPT(MXRTP)
C
C IMPRopers
C     MXRTI maximum number of impropers in topology file
C     MIM(*) first improper atom
C     MJM(*) second improper atom
C     MKM(*) third improper atom
C     MLM(*) fourth improper atom
C     DELMT(*) improper option
C
      INTEGER MXRTI
      PARAMETER (MXRTI=2000)
      CHARACTER*5 MIM(MXRTI),MJM(MXRTI),MKM(MXRTI),MLM(MXRTI)
      CHARACTER*4 DELMT(MXRTI)
C
C COMMON blocks
C =============
C
C     double precision block
C
      COMMON /RTF/ MASDEF, ARMASS, CHG
C
C     integer block
C
      COMMON /IRTF/ IATC, NRTRS, NIC, RTRTYP, MXN, GRPR
C
C     character string block
C
      COMMON /CRTF/ ATCT,
     1              DELAT, FTP, AA, MAC,
     2              MNB,
     3              DELBD, MIB, MJB,
     4              DELAN, MIT, MJT, MKT,
     5              DELPT, MIP, MJP, MKP, MLP,
     6              DELMT, MIM, MJM, MKM, MLM
C
      SAVE /RTF/
      SAVE /IRTF/
      SAVE /CRTF/

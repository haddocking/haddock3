C ener.inc
C
C total number of energy terms (note the free slots)
      INTEGER NENERT
C changes for pararestraints (AB, Sept. 2011)
C changes for z-harmonic restraints (AB, May 2013)
C      PARAMETER (NENERT=40)
      PARAMETER (NENERT=45)
C
C PIGWGHT weight of the Pairs of Interacting Groups in the Hamiltonian.
C PIGAVWT weight of the Pairs of Interacting Groups in the perturbing potential
      INTEGER NPIGMA3
      PARAMETER (NPIGMA3=MAX_PIGS)
      DOUBLE PRECISION PIGWGHT(NPIGMA3,NENERT), PIGAVWT(NPIGMA3,NENERT)
C TBATH bath temperature for dynamics with T-coupling
      DOUBLE PRECISION TBATH
C
C index parameters
      INTEGER SSTOTE, SSTOTK, SSENER, SSTEMP, SSBOND, SSANGL
      INTEGER SSDIHE, SSIMPR, SSVDW,  SSELEC, SSNCS,  SSNOE
      INTEGER SSHARM, SSCDIH, SSSD,   SSPLN,  SSPVDW, SSPELE
      INTEGER SSXREF, SSDG
C=====================================================================
Cif defined(CNS_SOLVE_COMPILE)
C=====================================================================
      INTEGER SSJCOUP, SSCARB, SSPROT, SSONEJ, SSRAMA
      INTEGER SSANGDB, SSSANI, SSDANI, SSVEAN, SSCOLL
      INTEGER SSZHAR
C     ------------------------------------------------------------------
C     These are the energy terms related to PARArestraints (AB, Sept. 2006)
C
      INTEGER SSXPCS, SSXRDC, SSXCCR, SSXANGLE
C=====================================================================
Cendif
C=====================================================================
C
      PARAMETER (SSTOTE=1)
      PARAMETER (SSTOTK=2)
      PARAMETER (SSENER=3)
      PARAMETER (SSTEMP=4)
      PARAMETER (SSSD=5)
      PARAMETER (SSBOND=6)
      PARAMETER (SSANGL=7)
      PARAMETER (SSDIHE=8)
      PARAMETER (SSIMPR=9)
      PARAMETER (SSVDW=10)
      PARAMETER (SSELEC=11)
      PARAMETER (SSPVDW=12)
      PARAMETER (SSPELE=13)
      PARAMETER (SSHARM=14)
      PARAMETER (SSCDIH=15)
      PARAMETER (SSNCS=16)
      PARAMETER (SSNOE=17)
      PARAMETER (SSXREF=18)
      PARAMETER (SSPLN=19)
      PARAMETER (SSDG=20)
C=====================================================================
Cif defined(CNS_SOLVE_COMPILE)
C=====================================================================
      PARAMETER (SSANGDB=30)
      PARAMETER (SSJCOUP=31)
      PARAMETER (SSCARB=32)
      PARAMETER (SSPROT=33)
      PARAMETER (SSONEJ=34)
      PARAMETER (SSRAMA=35)
      PARAMETER (SSSANI=36)
      PARAMETER (SSDANI=37)
      PARAMETER (SSVEAN=38)
      PARAMETER (SSCOLL=39)

C     ------------------------------------------------------------------
C     This is the energy term related to ZHARmomic restraints (2013)
C
      PARAMETER (SSZHAR=40)
C     ------------------------------------------------------------------
C     These are the energy terms related to PARArestraints (2003)
C
      PARAMETER (SSXPCS=41)
      PARAMETER (SSXRDC=42)
      PARAMETER (SSXCCR=43)
      PARAMETER (SSXANGLE=44)
C=====================================================================
Cendif
C=====================================================================
C
C The energy term names ANER are to be used for identification in
C command syntax usage. They and the logical flags QENER are
C initialised in ENEINI.
      CHARACTER*4 ANER(NENERT)
      DOUBLE PRECISION RENR(NENERT)
      DOUBLE PRECISION REN2(NENERT), RENV(NENERT)
      LOGICAL QENER(NENERT)
C
C the value of the previous call to PRINTE, PRINTD is held in EOLD
      DOUBLE PRECISION EOLD
C
C EXTRAE indicates whether PRINTE, PRINTD have to print an extended list
      LOGICAL EXTRAE
C
C accumulation variables
      DOUBLE PRECISION RENRP(NENERT), RENR2P(NENERT), RENRA(NENERT)
      DOUBLE PRECISION RENR2A(NENERT)
      INTEGER ISTPSA
C
C COMMON-BLOCKS
C
C     double precision block
C
      COMMON /DENER/ RENR, EOLD, RENRP, RENR2P, RENRA, RENR2A, REN2,
     &               RENV, PIGWGHT, PIGAVWT, TBATH
C
C     logical block
C
      COMMON /ENACCM/ QENER, EXTRAE
C
C     integer block
C
      COMMON /IENER/ ISTPSA
C
C     character string block
C
      COMMON /CAENR/ ANER
C
      SAVE /IENER/
      SAVE /DENER/
      SAVE /ENACCM/
      SAVE /CAENR/
C

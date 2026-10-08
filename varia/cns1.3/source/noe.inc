C =======
C noe.inc
C =======
C
C  data structure (COMMON block) for NOE restraints
C
C maximum number of protons in each selected group (dynamic
C assignment)
      INTEGER NOESMX
C
C maximum number of NOE restraints of current heap allocation
      INTEGER NOEMAX
C
C actual number of NOE restraints:
      INTEGER NOENUM
C
C number of OR restraints
      INTEGER NOEORN
C
C indices for linked lists for "OR" restraints
      INTEGER HPNORR
C
C linked list for set-I atoms (heap pointers)
      INTEGER HPNIPR, HPNILS
C
C linked list for set-J atoms (heap pointers)
      INTEGER HPNJPR, HPNJLS
C
C list for distances and +/- error estimates (heap pointers)
      INTEGER HPNDIS, HPNLOW, HPNHIG
C
C list for initial DEN distances
      INTEGER HPDENINIT
C
C individual weights on restraints
      INTEGER HPNWGH
C
C list for time-averaged distances
      INTEGER HPNRAV, HPNRRV
C list of distance evaluations since last reset
      INTEGER HPNNSP
C
C list of chemical shifts for the two selections
      INTEGER HPNPP1, HPNPP2
C
C volume estimate (as opposed to distance)
      INTEGER HPNVOL
C
C classification names and index list NOECND (heap pointer)
      INTEGER NOECMX
      PARAMETER (NOECMX=20)
      INTEGER HPNCND, NOECCN
      CHARACTER*4 NOECNM(NOECMX)
C
C list for potential and averaging options (class list)
      INTEGER NOEAVE(NOECMX), NOEPOT(NOECMX)
      INTEGER NOEM3, NOEM6, NOECEN, NOESUM, NOEBIH, NOESQW, NOESOF
      INTEGER NOESYM
      PARAMETER (NOEM6 =1, NOECEN=2, NOESUM=3,NOEM3=4)
      PARAMETER (NOEBIH=1, NOESQW=2, NOESOF=3, NOESYM=4)
      INTEGER NOEHDI,NOE3DP,NOELOG
      PARAMETER (NOEHDI=6, NOE3DP=7, NOELOG=8)
C
C list for square-well constant and exponent (class list)
      DOUBLE PRECISION NOECNS(NOECMX),NOEOFF(NOECMX),NOEMOF(NOECMX)
      INTEGER NOEEXP(NOECMX)
C
C list for asymptotic behaviour ("soft" square well)
      DOUBLE PRECISION NOESWI(NOECMX),NOEASY(NOECMX)
      INTEGER NOESOE(NOECMX)
C list for asymptotic behaviour on minus side ("soft" square well)
      DOUBLE PRECISION NOEMSW(NOECMX),NOEMAS(NOECMX)
      INTEGER NOEMSO(NOECMX)
C
C list for time-averaging of interproton distances
      INTEGER MEMTAU(NOECMX), DAVEXP(NOECMX), RAVEXP(NOECMX)
C options for time-averaging
      CHARACTER*4 TMODE(NOECMX), TFMODE(NOECMX)
C options for running average
      CHARACTER*4 RMODE(NOECMX)
C option for anaylse of violation with time-average
      CHARACTER*4 RANA
C
C maximum number of allowed interaction for ensemble-averaging
C with dimer molecules
      INTEGER MAXMOD
      PARAMETER (MAXMOD=120)
C
C option for ensemble-averaging
      CHARACTER*4 EMODE,DIENS,ENSMOD(MAXMOD,2)
      INTEGER IDIMER
      LOGICAL OCCU
C
C cross-validation array
      INTEGER HPNCV
C
C cross-validation test number
      INTEGER NOEICV
C
C high-dimensional NOE stuff.
C number of correlated noes for the high dimensional potential
C and the barrier hight
      INTEGER NOECOR(NOECMX)
      DOUBLE PRECISION NOEBAR(NOECMX)
C NOE scale (class list)
      DOUBLE PRECISION NOESCA(NOECMX)
C
C number of monomers in multimeric molecule (class list)
      INTEGER NMONO(NOECMX)
C
C exponent used for SUM average (class list)
      DOUBLE PRECISION NAVEXP(NOECMX)
C
C exponent used for OR restraints
      DOUBLE PRECISION OREXP(NOECMX)
C
C NOE ceiling, and temperature (global parameters)
      DOUBLE PRECISION NOECEI, NOETEM
C
C violation statistics
      DOUBLE PRECISION VIOTHR, VIOEXC
      INTEGER HPNVIO, VIOSTR
C
C ambiguous assignment statistics
      INTEGER HPNHGL, PNATOM, HPNMAT, HPNMA2, MAXSEL, MATDIM,
     &        AMBMAX, AMBMIN, AMBSTR, RNGMIN, RNGMAX
      DOUBLE PRECISION AMBLEV, AMBCUT
      CHARACTER*4 AMBMOD,AVEMOD,AMBFRM
C
C peak identifier
      INTEGER NPEAKI, HPNPID
C
C calibration matrix
      INTEGER CALCNT, MAXCAL, HPNCAL
      PARAMETER (MAXCAL=20)
      DOUBLE PRECISION VREFER(MAXCAL,MAXCAL), DREFER(MAXCAL,MAXCAL)
      DOUBLE PRECISION VCOUNT(MAXCAL,MAXCAL), DCOUNT(MAXCAL,MAXCAL)
      DOUBLE PRECISION CALFAC(MAXCAL,MAXCAL)
      DOUBLE PRECISION CALEXP, CALDIS
C
C pseudo atom corrections
      INTEGER HPNPSE
C diagonal intensities
      INTEGER HPNDIA
C
C
C COMMON BLOCK
C
C     double precision block
C
      COMMON /NOE/  NOECNS, NOEOFF, NOEMOF,
     &              NOESCA, NOECEI, NOETEM,
     &              NOESWI, NOEASY, NOEMSW, NOEMAS, NOEBAR,
     &              VIOTHR, VIOEXC, AMBLEV, AMBCUT,
     &              VREFER, DREFER, CALEXP, CALDIS, NAVEXP, OREXP,
     &              VCOUNT, DCOUNT, CALFAC
C
C     integer block
C
      COMMON /INOE/ NOESMX, NOEMAX, NOENUM, NOEORN,
     &              HPNORR, HPNIPR, HPNILS, HPNJPR, HPNJLS,
     &              HPNDIS, HPNLOW, HPNHIG, HPNWGH, HPNRAV, HPNRRV,
     &              HPNNSP, HPNCND, NOECCN, NOEAVE, NOEPOT,
     &              NOEEXP, NOESOE, NOEMSO, HPNCV,  MEMTAU, DAVEXP,
     &              NMONO,  RAVEXP, NOECOR, NOEICV, HPNVIO, VIOSTR,
     &              NPEAKI, HPNPID, HPNHGL, PNATOM, HPNMAT, HPNMA2,
     &              MAXSEL, MATDIM, HPNPP1, HPNPP2, HPNVOL, 
     &              HPNCAL, CALCNT, HPNPSE, AMBMIN, AMBMAX, AMBSTR,
     &              RNGMIN, RNGMAX, HPNDIA, IDIMER, HPDENINIT
C
C     character string block
C
      COMMON /CNOE/ NOECNM, TMODE, TFMODE, RMODE, AMBMOD, AVEMOD,
     &              AMBFRM, EMODE, DIENS, ENSMOD, RANA
C
C
C     logical string block
C
      COMMON /LNOE/ OCCU
C
      SAVE /NOE/
      SAVE /INOE/
      SAVE /CNOE/
      SAVE /LNOE/

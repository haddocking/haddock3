C     ------------------------------------------------------------------
C     RDC-angles stuff
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INTEGER MAXXANGLESCLASSES
      PARAMETER (MAXXANGLESCLASSES = 40)  
C     ------------------------------------------------------------------
C     Arrays that hold RDC-angles info
C     XANGLESASSNDX tells ending index of the RDC-angles arrays for each
C                   class
C     XANGLESFORCES holds K1 for each class
C     ------------------------------------------------------------------
      INTEGER XANGLESASSNDX(MAXXANGLESCLASSES)
      REAL XANGLESFORCES(MAXXANGLESCLASSES)
      CHARACTER*8 XANGLESCLASSNAMES(MAXXANGLESCLASSES)
      LOGICAL PRINTCLASS(MAXXANGLESCLASSES)
C     ------------------------------------------------------------------
C     MAXXANGLES is the number of slots set aside for RDC-angles
C                assignments
C     ANGLESNUM is the total number of RDC-angles entered
C     ------------------------------------------------------------------
      INTEGER MAXXANGLES,ANGLESNUM,NCLASSES,CURCLASS,XANGLESMX
C     ------------------------------------------------------------------
C     Pointers to arrays to hold atom numbers, observed THETA and PHI,
C     and errors
C     ------------------------------------------------------------------
      INTEGER XANGLESJPTR,XANGLESKPTR,XANGLESJLST,XANGLESKLST,
     &        XANGLESOBS1PTR,XANGLESOBS2PTR,XANGLESERRPTR,CALCXANGLESPTR
C     ------------------------------------------------------------------
C     Input modes
C     ------------------------------------------------------------------
      INTEGER MODE,NEW,UPDATE
      PARAMETER (NEW=1)
      PARAMETER (UPDATE=2)
C     ------------------------------------------------------------------
C     Parameters as set up in ETOR - Not used indeed
C     ------------------------------------------------------------------
      DOUBLE PRECISION MCONST
      PARAMETER(MCONST=0.0001D0)
      DOUBLE PRECISION EPS
      PARAMETER(EPS=0.1D0)
C     ------------------------------------------------------------------
C     Common blocks
C     ------------------------------------------------------------------
      COMMON /CXANGLES/XANGLESCLASSNAMES
      COMMON /IXANGLES1/XANGLESASSNDX,MAXXANGLES,ANGLESNUM,
     &                  CURCLASS,NCLASSES,XANGLESMX
      COMMON /IXANGLES2/XANGLESJPTR,XANGLESKPTR,XANGLESJLST,XANGLESKLST,
     &                  XANGLESOBS1PTR,XANGLESOBS2PTR,
     &                  XANGLESERRPTR,CALCXANGLESPTR,MODE
      COMMON /RXANGLES/XANGLESFORCES
      COMMON /LXANGLES/PRINTCLASS

      SAVE /CXANGLES/
      SAVE /IXANGLES1/
      SAVE /IXANGLES2/
      SAVE /RXANGLES/
      SAVE /LXANGLES/

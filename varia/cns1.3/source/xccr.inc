C     ------------------------------------------------------------------
C     Cross correlation rate stuff
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INTEGER MAXXCCRCLASSES
      PARAMETER (MAXXCCRCLASSES = 40)  
C     ------------------------------------------------------------------
C     Arrays that hold cross correlation rate info
C     XCCRASSNDX tells ending index of the CCR arrays for each class
C     XCCRFORCES holds K1 for each class
C     ------------------------------------------------------------------
      INTEGER XCCRASSNDX (MAXXCCRCLASSES)
      REAL XCCRFORCES(MAXXCCRCLASSES)
      REAL XCCRCOEF1(MAXXCCRCLASSES)
      CHARACTER*8 XCCRCLASSNAMES(MAXXCCRCLASSES)
      LOGICAL PRINTCLASS(MAXXCCRCLASSES),CCRFLG
C     ------------------------------------------------------------------
C     MAXXCCR is the number of slots set aside for CCR assignments
C     CCRNUM is the total number of CCR entered
C     ------------------------------------------------------------------
      INTEGER MAXXCCR,CCRNUM,NCLASSES,CURCLASS,XCCRMX
C     ------------------------------------------------------------------
C     Pointers to arrays to hold atom numbers, observed CCR and errors
C     ------------------------------------------------------------------
      INTEGER XCCRIPTR,XCCRJPTR,XCCRKPTR,XCCRILST,XCCRJLST,XCCRKLST,
     &        XCCROBSPTR,XCCRERRPTR,CALCXCCRPTR
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
      COMMON /CXCCR/XCCRCLASSNAMES
      COMMON /IXCCR1/XCCRASSNDX,MAXXCCR,CCRNUM,CURCLASS,NCLASSES,XCCRMX
      COMMON /IXCCR2/XCCRIPTR,XCCRJPTR,XCCRKPTR,
     &               XCCRILST,XCCRJLST,XCCRKLST,
     &               XCCROBSPTR,XCCRERRPTR,CALCXCCRPTR,MODE
      COMMON /RXCCR/XCCRFORCES,XCCRCOEF1
      COMMON /LXCCR/PRINTCLASS
      COMMON /WXCCR/CCRFLG

      SAVE /CXCCR/
      SAVE /IXCCR1/
      SAVE /IXCCR2/
      SAVE /RXCCR/
      SAVE /LXCCR/
      SAVE /WXCCR/

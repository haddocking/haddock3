C     ------------------------------------------------------------------
C     Residual dipolar coupling stuff
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INTEGER MAXXRDCCLASSES
      PARAMETER (MAXXRDCCLASSES = 40)  
C     ------------------------------------------------------------------
C     Arrays that hold residual dipolar coupling info
C     XRDCASSNDX tells ending index of the RDC arrays for each class
C     XRDCFORCES holds K1 for each class
C     ------------------------------------------------------------------
      INTEGER XRDCASSNDX(MAXXRDCCLASSES)
      REAL XRDCFORCES(MAXXRDCCLASSES)
      REAL XRDCCOEF1(MAXXRDCCLASSES),XRDCCOEF2(MAXXRDCCLASSES)
      CHARACTER*8 XRDCCLASSNAMES(MAXXRDCCLASSES)
      LOGICAL PRINTCLASS(MAXXRDCCLASSES),DIPFLG
C     ------------------------------------------------------------------
C     MAXXRDC is the number of slots set aside for RDC assignments
C     RDCNUM is the total number of RDC entered
C     ------------------------------------------------------------------
      INTEGER MAXXRDC,RDCNUM,NCLASSES,CURCLASS,XRDCMX
C     ------------------------------------------------------------------
C     Pointers to arrays to hold atom numbers, observed RDC and errors
C     ------------------------------------------------------------------
      INTEGER XRDCIPTR,XRDCJPTR,XRDCKPTR,XRDCLPTR,XRDCMPTR,XRDCNPTR,
     &        XRDCILST,XRDCJLST,XRDCKLST,XRDCLLST,XRDCMLST,XRDCNLST,
     &        XRDCOBSPTR,XRDCERRPTR,CALCXRDCPTR
C     ------------------------------------------------------------------
C     Input modes
C     ------------------------------------------------------------------
      INTEGER MODE,NEW,UPDATE
      PARAMETER(NEW = 1)
      PARAMETER(UPDATE = 2)
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
      COMMON /CXRDC/XRDCCLASSNAMES
      COMMON /IXRDC1/XRDCASSNDX,MAXXRDC,RDCNUM,CURCLASS,NCLASSES,XRDCMX
      COMMON /IXRDC2/XRDCIPTR,XRDCJPTR,XRDCKPTR,
     &               XRDCLPTR,XRDCMPTR,XRDCNPTR,
     &               XRDCILST,XRDCJLST,XRDCKLST,
     &               XRDCLLST,XRDCMLST,XRDCNLST,
     &               XRDCOBSPTR,XRDCERRPTR,CALCXRDCPTR,MODE
      COMMON /RXRDC/XRDCFORCES,XRDCCOEF1,XRDCCOEF2
      COMMON /LXRDC/PRINTCLASS
      COMMON /WXRDC/DIPFLG

      SAVE /CXRDC/
      SAVE /IXRDC1/
      SAVE /IXRDC2/
      SAVE /RXRDC/
      SAVE /LXRDC/
      SAVE /WXRDC/

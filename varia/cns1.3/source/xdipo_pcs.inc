C     ------------------------------------------------------------------
C     Pseudocontact shift stuff
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      INTEGER MAXXPCSCLASSES
      PARAMETER(MAXXPCSCLASSES=40)  
C     ------------------------------------------------------------------
C     Arrays that hold pseudocontact shift info
C     XPCSASSNDX tells ending index of the PCS arrays for each class
C     XPCSFORCES holds K1 for each class
C     ------------------------------------------------------------------
      INTEGER XPCSASSNDX(MAXXPCSCLASSES)
      REAL XPCSFORCES(MAXXPCSCLASSES)
      REAL XPCSCOEF1(MAXXPCSCLASSES),XPCSCOEF2(MAXXPCSCLASSES)
      CHARACTER*8 XPCSCLASSNAMES(MAXXPCSCLASSES)
      LOGICAL PRINTCLASS(MAXXPCSCLASSES),PCSFLG
C     ------------------------------------------------------------------
C     MAXXPCS is the number of slots set aside for PCS assignments
C     PCSNUM is the total number of PCS entered
C     ------------------------------------------------------------------
      INTEGER MAXXPCS,PCSNUM,NCLASSES,CURCLASS,XPCSMX
C     ------------------------------------------------------------------
C     Pointers to arrays to hold atom numbers, observed PCS and errors
C     ------------------------------------------------------------------
      INTEGER XPCSIPTR,XPCSJPTR,XPCSKPTR,XPCSLPTR,XPCSMPTR,XPCSNPTR,
     &        XPCSILST,XPCSJLST,XPCSKLST,XPCSLLST,XPCSMLST,XPCSNLST,
     &        XPCSOBSPTR,XPCSERRPTR,CALCXPCSPTR
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
      COMMON /CXPCS/XPCSCLASSNAMES
      COMMON /IXPCS1/XPCSASSNDX,MAXXPCS,PCSNUM,CURCLASS,NCLASSES,XPCSMX
      COMMON /IXPCS2/XPCSIPTR,XPCSJPTR,XPCSKPTR,
     &               XPCSLPTR,XPCSMPTR,XPCSNPTR,
     &               XPCSILST,XPCSJLST,XPCSKLST,
     &               XPCSLLST,XPCSMLST,XPCSNLST,
     &               XPCSOBSPTR,XPCSERRPTR,CALCXPCSPTR,MODE
      COMMON /RXPCS/XPCSFORCES,XPCSCOEF1,XPCSCOEF2
      COMMON /LXPCS/PRINTCLASS
      COMMON /WXPCS/PCSFLG

      SAVE /CXPCS/
      SAVE /IXPCS1/
      SAVE /IXPCS2/
      SAVE /RXPCS/
      SAVE /LXPCS/
      SAVE /WXPCS/

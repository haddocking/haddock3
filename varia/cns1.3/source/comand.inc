C comand.inc
C
C parsing data structure
C
C ASCII character mapping arrays
      INTEGER MASCII
      PARAMETER (MASCII=255)
      INTEGER ASCIIM(MASCII), ASCIIS(MASCII)
C
C words
      INTEGER WDLEN, WDMAX, WDTLEN, WDTMAX, WDDLEN, WDDMAX
      PARAMETER (WDMAX=WORD_SIZE,WDTMAX=WORD_SIZE,WDDMAX=WORD_SIZE)
C changed LEN(WD) to STRING_SIZE+1 to allow overflows past WDMAX, Joe Krahn, 2/24/10
      CHARACTER*(STRING_SIZE+1) WD
      CHARACTER*(WDTMAX) WDT
      CHARACTER*(WDDMAX) WDD
C
C current line (used only in NEXTWD)
      INTEGER COMLEN, COMMAX
      PARAMETER (COMMAX=STRING_SIZE)
      CHARACTER*(COMMAX) COMLYN
C
C display line
      INTEGER DISLEN, DISMAX, DUNIT, PUNIT
      PARAMETER (DISMAX=STRING_SIZE)
      CHARACTER*(DISMAX) DISLYN
C
C protocol unit
      INTEGER PRUNIT
C
C input and output filename
      CHARACTER*(STRING_SIZE) IFILE, OFILE
C
C display filename
      CHARACTER*(STRING_SIZE) DFILE
C
C current cursor and position
      CHARACTER CHR
      INTEGER CURSOR, CUROLD
C
C stack for stream FORTRAN units and "hold" lines
      INTEGER   MXSTRM, NSTRM
      PARAMETER (MXSTRM=10)
      INTEGER   ISTRM(MXSTRM), HLDLEN(MXSTRM)
      LOGICAL   QADD(MXSTRM)
      CHARACTER*(COMMAX) HLDLYN(MXSTRM)
C
C
C time info
      DOUBLE PRECISION STRCPU
      CHARACTER*8 STRTIM
      CHARACTER*11 STRDAT
C
C flags for parsing:
C QTERM: terminal interactive characteristics
C QECHO: make an echo line after reading a line
C QSUBS: make substitutions of words beginning with a "$"
C QEXPRS: accept +,-,*,/,#,~ as single character words
C QFILNM: filename parsing mode, only "=" is single character word
C QCAPIT: capitalize words (only effective when QFILNM is FALSE)
C QQUOT: indicates that last word parsed is a quoted string
C EOF: end of file encountered during parsing
C ERROR: read error encountered during parsing
C QSAVEW: save the current word, i.e. the next word is the current word
C DONE: parsing flag
C USED: parsing flag
      LOGICAL QTERM, QECHO, QSUBS, EOF, ERROR, QSAVEW, DONE, USED
      LOGICAL QEXPRS, QQUOT, QFILNM, QCAPIT
C
C QNEWLN: flag that gets set to true every time a new line is read
C (for distinguishing "end if" from " ...end \LF if.." )
      LOGICAL QNEWLN
C
C
C END stack for control statements
      INTEGER ENDMAX
      PARAMETER (ENDMAX=50)
      INTEGER ENDIND, ENDCUR(ENDMAX), ENDRET(ENDMAX), ENDCNT(ENDMAX)
      INTEGER ENDPTR(ENDMAX), ENDNAT(ENDMAX), ENDSEL(ENDMAX)
      CHARACTER*16 ENDKEY(ENDMAX)
      CHARACTER*4 ENDACT(ENDMAX), ENDLBL(ENDMAX)
C
C rotating command buffer for loops
      INTEGER BUFMAX
      PARAMETER (BUFMAX=NUM_BUFFER)
      INTEGER BUFIND, BUFFIL, BUFLEN(BUFMAX), BUFSTK
      CHARACTER*(COMMAX) BUFLYN(BUFMAX)
C
C general symbol handling variables
      INTEGER DBUFMX
      PARAMETER (DBUFMX=NUM_BUFFER)
      INTEGER DEFSCPMAX
      PARAMETER (DEFSCPMAX=9)
C array for inline scope count - used by DEFNEWSCOPE and DEFKILSCOPE
C to keep track of inlining directives
      INTEGER DEFINLCNT(DEFSCPMAX)
C when set to true, signifies that the next scope encountered is to be inlined
      LOGICAL QINLINE
C
C scopes for DEFINE and MACRO commands
      INTEGER DEFSETSCOPE
      INTEGER DEFEVLSCOPE
C current scope
      INTEGER DEFCURSCOPE
C number of variable types
      INTEGER DEFTYPMAX
      PARAMETER (DEFTYPMAX=4)
C index of the first symbol(1), parameter(2), procedure(3)
C or buffer(4) in each scope (NOTE: buffers currently only use scope 1)
      INTEGER DEFVARLIST(DEFSCPMAX,DEFTYPMAX)
C first entry in the buffered line list for each stream
      INTEGER DEFBUFLIST(MXSTRM)
C variable text and text length
      CHARACTER*(COMMAX) DEFPARTXT(DBUFMX)
      INTEGER DEFPARLEN(DBUFMX)
      INTEGER VARMAX
      PARAMETER (VARMAX=WORD_SIZE)
C variable name and name length
      CHARACTER*(VARMAX) DEFNAMTXT(DBUFMX)
      INTEGER DEFNAMLEN(DBUFMX)
C link arrays
C ...for multi-line parameters
      INTEGER DEFMULTREC(DBUFMX)
C ...for siblings
      INTEGER DEFNEXTREC(DBUFMX)
C ...for descendants
      INTEGER DEFBRANREC(DBUFMX)
C free link
      INTEGER DEFEMPTYREC
C flag indicating that define subsitutions should be performed
      LOGICAL QDEF
      LOGICAL QSCOPING
C
      INTEGER NUMVARFLG
      PARAMETER(NUMVARFLG=5)
C
C typed symbols
      CHARACTER*2 CMDTYP(DBUFMX)
      DOUBLE COMPLEX DCVALU(DBUFMX)
      DOUBLE PRECISION DPVALU(DBUFMX)
C
C error count
      INTEGER NERRPA
C
C buffered output
      INTEGER MXBFLN
      PARAMETER (MXBFLN=16384)
      INTEGER BODEST(DBUFMX)
      INTEGER BOUNIT(DBUFMX)
C
C
C-COMMON-DECLARATIONS
C
C     double precision block
C
      COMMON /DCOMND/ STRCPU, DCVALU, DPVALU
C
C     logical block
C
      COMMON /LCOMND/ QADD,    QTERM,   QECHO,   QSUBS,   EOF,
     &                ERROR,   QSAVEW,  DONE,    USED,    QEXPRS,
     &                QQUOT, QFILNM, QCAPIT, QDEF, QSCOPING, QNEWLN,
     &                QINLINE
C
C     integer block
C
      COMMON /COMAND/ WDLEN,   WDTLEN,  WDDLEN,  COMLEN,  DISLEN,
     &                CURSOR,  NSTRM,   ISTRM,   DUNIT,   PUNIT,
     &                PRUNIT,  HLDLEN,  ENDIND,  ENDCUR,  ENDRET,
     &                ENDCNT,  BUFIND,  BUFFIL,  BUFLEN,  BUFSTK,
     &                CUROLD,  NERRPA,  ENDPTR,
     &                ENDNAT,  ENDSEL,  ASCIIM,  ASCIIS,
     &                DEFSETSCOPE, DEFCURSCOPE, DEFVARLIST,
     &                DEFBUFLIST,  DEFPARLEN,   DEFNAMLEN,
     &                DEFNEXTREC,  DEFBRANREC,  DEFEMPTYREC,
     &                DEFMULTREC,  DEFEVLSCOPE, DEFINLCNT,
     &                BODEST, BOUNIT
C
C
C     character string block
C
      COMMON /CCOMND/ CHR, WD, WDT, WDD, COMLYN, DISLYN,
     &                STRTIM, STRDAT, IFILE, OFILE, DFILE,
     &                HLDLYN, ENDKEY, ENDACT, BUFLYN, ENDLBL,
     &                DEFPARTXT, DEFNAMTXT, CMDTYP
C
      SAVE /COMAND/
      SAVE /CCOMND/
      SAVE /DCOMND/
      SAVE /LCOMND/

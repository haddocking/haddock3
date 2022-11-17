      SUBROUTINE EVALTB( MODE, INDEX, STINDX, FLAGS, XRERR )
      IMPLICIT NONE
C
C     evaluate the parse tables to obtain a value for the right
C     hand side of the equation
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      LOGICAL XRERR
      INTEGER INDEX
      INTEGER FLAGS(*)
      CHARACTER*4 MODE
      INTEGER K, L, STINDX
C
Ctabtrace       WRITE(19,'(1X,A)')'INITIAL TABLE SETUP:'
Ctabtrace       WRITE(19,'(1X,A,A)')'MODE = ',MODE
Ctabtrace       CALL TBOUT
C
      K = 1
      L = K + 1
1     IF ( TEO(1) .NE. '<' .AND. ( VOP(L) .LE. VOP(K) ) ) THEN
C
         IF     ( TEO(K) .EQ. '+' ) THEN
            CALL XRADD( K, XRERR )
         ELSEIF ( TEO(K) .EQ. '-' ) THEN
            CALL XRSUB( K, XRERR )
         ELSEIF ( TEO(K) .EQ. '*' ) THEN
            CALL XRMULT( K, XRERR )
         ELSEIF ( TEO(K) .EQ. ',' ) THEN
            CALL XRCOMA( K, MODE, XRERR )
         ELSEIF ( TEO(K) .EQ. '/' ) THEN
            CALL XRDIV( K, XRERR )
         ELSEIF ( TEO(K) .EQ. '^' ) THEN
            CALL XRPOWR( K, XRERR )
         ELSEIF ( TEO(K) .EQ. '=' ) THEN
            CALL XREQU( K, MODE,  XRERR )
         ELSEIF ( TEO(K) .EQ. 'A' ) THEN
            CALL XRABS( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'M' ) THEN
            CALL XRMIN( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. '~' ) THEN
            CALL XRUNM( K, XRERR )
         ELSEIF ( TEO(K) .EQ. 'H' ) THEN
            CALL XRHEAV( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'Z' ) THEN
            CALL XRMAX( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'L' ) THEN
            CALL XRLOG( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'E' ) THEN
            CALL XREXP( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'G' ) THEN
            CALL XRGAUS( K, L, INDEX, XRERR )
         ELSEIF ( TEO(K) .EQ. 'W' ) THEN
            CALL XRMAXW( K, L, INDEX, XRERR )
CCC deleted NORM function ATB 25/12/06
CCC         ELSEIF ( TEO(K) .EQ. 'N' ) THEN
CCC            CALL XRNORM( K, L, INDEX, STINDX, FLAGS, XRERR )
         ELSEIF ( TEO(K) .EQ. 'B' ) THEN
            CALL XRCLOG( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'D' ) THEN
            CALL XRSIGN( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'F' ) THEN
            CALL XRASIN( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'J' ) THEN
            CALL XRACOS( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'S' ) THEN
            CALL XRSINE( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'C' ) THEN
            CALL XRCOS( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'T' ) THEN
            CALL XRTAN( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'Q' ) THEN
            CALL XRSQRT( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'I' ) THEN
            CALL XRINT( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'R' ) THEN
            CALL XRMOD( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'V' ) THEN
            CALL XRENCD( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'X' ) THEN
            CALL XRDECD( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'Y' ) THEN
            CALL XRRAND( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'b' ) THEN
            CALL XRNINT( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'c' ) THEN
            CALL XRCAPIT( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'a' ) THEN
            CALL XRIMOD( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'l' ) THEN
            CALL XRLEN( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 's' ) THEN
            CALL XRSUBSTR( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'i' ) THEN
            CALL XRINDEX( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'r' ) THEN
            CALL XRRINDEX( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'f' ) THEN
            CALL XRFORMAT( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'g' ) THEN
            CALL XRADJUSTL( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 't' ) THEN
            CALL XRTRIM( K, L, XRERR )
         ELSEIF ( TEO(K) .EQ. 'm' ) THEN
            CALL XRMATCH( K, L, XRERR )
         ENDIF
         IF ( XRERR ) THEN
            RETURN
         ENDIF
C
C        if last operation wasn't a noop remove it as processed
C
         IF(TEO(K) .NE. '_')THEN
            CALL TBPACK( K )
         ENDIF
Ctabtrace         CALL TBOUT
C
C        backup a position in the table for next comparison
C        if it is a noop go back another one
C
36       K = K - 1
         IF ( TEO(K) .EQ. '_' ) GOTO 36
C
C        want to compare to the next operation in the table
C        unless it is a noop
C
         L = K
37       L = L + 1
         IF ( TEO(L) .EQ. '_' ) GOTO 37
C
C        go back and make the comparison
C
         IF ( K .GT. 1 ) GOTO 1
      ELSE
C
C        continue working through the table of operators
C
38       K = K + 1
         IF ( TEO(K) .EQ. '_' ) GOTO 38
         L = K
39       L = L + 1
         IF ( TEO(L) .EQ. '_' ) GOTO 39
         GOTO 1
      ENDIF
      RETURN
      END
C
      SUBROUTINE HIER( MODE, XRERR )
      IMPLICIT NONE
C
C     this routine parses the string containing the equation and
C     stores all information in the parsing tables.
C     parameters:
C        mode -> type of calculation (real, comp, scal, phas, ampl)
C        xrerr -> returned true if this routine errs
C     parse tables:
C         teo -- is an array of the actual operators
C         vop -- holds the value of the operator in teo, this is
C                used to determine precedence when calculating
C         vareq -- holds character representations of all parts of
C                  the equation that are not operators
C         vartype -- tells whether the item in vareq is a double
C                    precision number, a double complex number, or
C                    a function
C         rvareq -- holds the value of vareq if it is a double
C                      precision number.
C         cvareq -- holds the value of vareq if it is a double
C                      complex number.
C         svareq -- holds the value of vareq if it is a string.
C
C
      CHARACTER*(*) MODE
      LOGICAL XRERR
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'timer.inc'
C     max variable length VNAMLN = WDMAX chars
      CHARACTER*(WDMAX) TMPWD
      INTEGER TWDLEN
      INTEGER K, I, INCPAR, VCPTR
      LOGICAL FOUND, INQUOT, OK
      DOUBLE PRECISION DECODF
      INTEGER PARBAL
      CHARACTER*1 LR, LSTOKN
      DOUBLE PRECISION DPVAL
      DOUBLE COMPLEX DCVAL
      CHARACTER*2 WDTYP
C
C     initialize everything
C
      CALL VINIT
      XRERR  = .FALSE.
      INQUOT = .FALSE.
      K = 0
      INCPAR = 0
      VCPTR = 0
      PARBAL = 0
      LR = 'L'
      IF ( MODE .EQ. 'EVAL' ) THEN
         XSTMAX = SVARMX
      ELSE
         XSTMAX = 4
      ENDIF
C
C     get next word; must be the left paren.
C
      CALL NEXTSL('Expression=')
C
      IF ( WD(1:4) .NE. '(   ' ) THEN
         CALL DSPERR('Expression>','Must begin with (')
         XRERR = .TRUE.
         RETURN
      ENDIF
C
C     if made it into the expression (found left paren) then continue
C     parsing until parens balance
C
101   CONTINUE
C
C     k points to next index in tables to place variable info
C
      TMPWD = WD
      TWDLEN = WDLEN
      IF ( TMPWD(1:4) .EQ. '(   '.AND..NOT.QQUOT ) THEN
         PARBAL = PARBAL + 1
         IF ( PARBAL .EQ. 1 ) THEN
            TMPWD = '>'
            TWDLEN = 1
            LSTOKN = 'O'
         ENDIF
      ELSEIF ( TMPWD(1:4) .EQ. ')   '.AND..NOT.QQUOT ) THEN
         PARBAL = PARBAL - 1
         IF ( PARBAL .EQ. 0 ) THEN
            TMPWD = '<'
            TWDLEN = 1
            LSTOKN = 'O'
         ENDIF
      ENDIF
      IF (QQUOT) THEN
C
C     it is a quoted string.
C
                 VARTYP(K) = 'ST'
                 SVAREQ(K) = ' '
                 SVAREQ(K)(1:WDLEN) = WD(1:WDLEN)
                 SVARLN(K) = WDLEN
C
C     look for operators and determine precedence
C
      ELSE IF ( TMPWD(1:4) .EQ. '>' ) THEN
         IF ( K .NE. 0 ) THEN
            CALL DSPERR('Expression>','Illegal character')
            XRERR = .TRUE.
            RETURN
         ENDIF
         K = K + 1
         VOP(K) = INCPAR
         TEO(K) = '>'
C
C        the show command 'show ( x )' is implemented by insert-
C        ing 'show=' into the command netting 'show (show=x)'
C        then when the assignment to show is made the value is printed
C
         IF ( MODE .EQ. 'SHOW' ) THEN
            VARTYP(1) = 'VT'
            VARNUM(1) = 29
            VAREQ(1) = 'SHOW                '
            TEO(2) = '='
            VOP(2) = 1
            K = 2
         ENDIF
      ELSEIF ( TMPWD(1:1) .EQ. '<' ) THEN
C
C     at end of equation check for balanced parens
C
         K = K + 1
         VOP(K) = INCPAR
         TEO(K) = '<'
      ELSEIF ( TMPWD(1:1) .EQ. '+' .OR. TMPWD(1:1) .EQ. '-' )THEN
         K = K + 1
         VOP(K) = 2 + INCPAR
         TEO(K) = TMPWD(1:1)
         LSTOKN = 'O'
      ELSEIF ( TMPWD(1:1) .EQ. '=' ) THEN
         IF ( K .GT. 1 ) THEN
            CALL DSPERR('EXPRESSION PARSER',
     &                  'TOO MANY TERMS TO LEFT OF EQUAL SIGN')
            XRERR = .TRUE.
            RETURN
         ENDIF
         LSTOKN = 'O'
         K = K + 1
         VOP(K) = 1 + INCPAR
         TEO(K) = '='
         LR = 'R'
      ELSEIF ( TMPWD(1:1) .EQ. ',' ) THEN
         K = K + 1
         VOP(K) = INCPAR
         TEO(K) = ','
         LSTOKN = 'O'
      ELSEIF ( TMPWD(1:1) .EQ. '^' ) THEN
         K = K + 1
         VOP(K) = 5 + INCPAR
         TEO(K) = '^'
         LSTOKN = 'O'
      ELSEIF ( TMPWD(1:1) .EQ. '~' ) THEN
         K = K + 1
         VOP(K) = 6 + INCPAR
         TEO(K) = '~'
         IF ( LSTOKN .EQ. 'S' ) THEN
            CALL DSPERR('EXPRESSION-EVALUATOR',
     &            'Illegal use of unary minus')
            XRERR = .TRUE.
            RETURN
         ELSE
            LSTOKN = 'O'
         ENDIF
      ELSEIF ( TMPWD(1:1) .EQ. '*' ) THEN
         K = K + 1
         CALL NEXTSL('Expression=')
         IF ( WD(1:1) .NE. '*' ) THEN
            CALL SAVEWD
            VOP(K) = 3 + INCPAR
            TEO(K) = '*'
         ELSE
            VOP(K) = 5 + INCPAR
            TEO(K) = '^'
         ENDIF
         LSTOKN = 'O'
      ELSEIF ( TMPWD(1:1) .EQ. '/' ) THEN
         K = K + 1
         VOP(K) = 3 + INCPAR
         TEO(K) = '/'
         LSTOKN = 'O'
      ELSEIF ( TMPWD(1:1) .EQ. '(' ) THEN
C
C        boost precedence level when within parens
C
         INCPAR = INCPAR + 10
      ELSEIF ( TMPWD(1:1) .EQ. ')' ) THEN
         INCPAR = INCPAR - 10
      ELSE
C
C        if not an operator assume it is a
C        variable, function, or a literal constant
C
         VAREQ(K) = TMPWD
         FOUND = .FALSE.
C
C to substitute it as a symbol
C ,ie $PI
         IF ( WD(1:1) .EQ. '$' ) THEN
               IF ( MODE .EQ. 'EVAL' .AND. LR .EQ. 'L' ) THEN
C
C in evaluate mode, the first symbol is to be assigned --
C don't replace it and store the variable name for the
C assignment.
C
                  VARTYP(K) = 'VT'
                  VAREQ(K) = WD
                  IF (WDLEN.GT.VNAMLN) THEN
      CALL DSPERR('HIER','name of symbol too long.')
                  END IF
                  FOUND = .TRUE.
               ELSE
                  TMPWD = WD
                  TWDLEN = WDLEN
                  CALL WDSUB( TMPWD,WDMAX,TWDLEN,FOUND,
     &                                    WDTYP, DPVAL, DCVAL )
      IF (.NOT.FOUND) THEN
      CALL DSPERR('WDSUB','symbol not found')
      END IF
                  IF ( FOUND ) THEN
                     VAREQ(K) = TMPWD
                     VARTYP(K) = WDTYP
                     IF ( WDTYP .EQ. 'DP' ) THEN
                        RVAREQ(K) = DPVAL
                     ELSEIF ( WDTYP .EQ. 'DC' ) THEN
                        CVAREQ(K) = DCVAL
                     ELSEIF ( WDTYP .EQ. 'ST') THEN
                        SVAREQ(K) = ' '
                        SVAREQ(K)(1:TWDLEN) = TMPWD(1:TWDLEN)
                        SVARLN(K) = TWDLEN
                     ELSEIF ( WDTYP .EQ. 'LO') THEN
                        SVAREQ(K) = ' '
                        SVAREQ(K)(1:TWDLEN) = TMPWD(1:TWDLEN)
                        SVARLN(K) = TWDLEN
                     ENDIF
                  ENDIF
               ENDIF
C
C check through VECTOR variables
         ELSE IF ( MODE .EQ. 'REAL' .OR. MODE .EQ. 'SHOW' ) THEN
            CALL SCHVAR( K, FOUND )
C
C check through XRAY-variables
         ELSE IF (MODE.NE.'EVAL') THEN
            CALL XCHVAR( K, FOUND )
         ENDIF
C
C        if not a valid variable name, check the functions
C
         IF ( .NOT. FOUND ) THEN
            CALL XRCHFX( MODE, K, FOUND )
            IF ( FOUND ) THEN
               VARTYP(K) = 'FX'
               K = K + 1
               VOP(K)  = 7 + INCPAR
               IF ( LSTOKN .EQ. 'S' ) THEN
                  CALL DSPERR('EXPRESSION-EVALUATOR',
     &                        'Function follows symbol')
                  XRERR = .TRUE.
                  RETURN
               ELSE
                  LSTOKN = 'O'
               ENDIF
            ELSE
C
C              if not a function we assume that word is a literal constant
C              Everything else is assumed to be an unquoted string.
C
C              Check first for quoted strings
C
               CALL CHKNUM(VAREQ(K),TWDLEN,OK)
C
               IF (OK) THEN
C
C              we think it is a number !! try to decode it
C
                 VARTYP(K) = 'DP'
                 RVAREQ(K) = DECODF(VAREQ(K),TWDLEN,OK)
C
C                here we provide an error message
C
                 IF (.NOT.OK) CALL DSPERR('EXPRESSION-EVALUATOR',
     &            'Error converting number')
      ELSEIF (WD(1:WDLEN).EQ.'TRUE'.OR.WD(1:WDLEN).EQ.'FALSE') THEN
C
C logical
                 VARTYP(K) = 'LO'
                 SVAREQ(K) = ' '
                 SVAREQ(K)(1:WDLEN) = WD(1:WDLEN)
                 SVARLN(K) = WDLEN
               ELSE
C
C                it must be a string.
C
                 VARTYP(K) = 'ST'
                 SVAREQ(K) = ' '
                 SVAREQ(K)(1:WDLEN) = WD(1:WDLEN)
                 SVARLN(K) = WDLEN
                 IF (WRNLEV .GE. 2 ) THEN
       WRITE(6,'(3A)') ' Assuming literal string "',WD(1:WDLEN),'"'
                 ENDIF
               ENDIF
C
               IF ( LSTOKN .EQ. 'S' ) THEN
                  CALL DSPERR('EXPRESSION-EVALUATOR',
     &                     'Illegal consecutive symbols')
                  XRERR = .TRUE.
                  RETURN
               ELSE
                  LSTOKN = 'S'
               ENDIF
            ENDIF
         ELSE
            IF ( LSTOKN .EQ. 'S' ) THEN
               CALL DSPERR('EXPRESSION-EVALUATOR',
     &            'Illegal consecutive symbols')
               XRERR = .TRUE.
               RETURN
            ELSE
               LSTOKN = 'S'
            ENDIF
         ENDIF
      ENDIF
C
C
      IF ( TMPWD(1:1) .NE. '<' ) THEN
C        if in evaluation mode get the variable name not it's value
         CALL NEXTSL('Expression=')
         GOTO 101
      ENDIF
C
C     total number of operations in equation
C
      NUMOPS = K
C
C     look for any unary minus signs, replace with ~
C
      DO I = 2,NUMOPS
         IF ( TEO(I) .EQ. '-' .AND. VARTYP(I-1) .EQ. '  '
     &                        .AND. VOP(I-1) .LT. VOP(I) + 5 )THEN
C
            TEO(I) = '~'
            VOP(I) = VOP(I) + 3
         ENDIF
      ENDDO
C
C     check to make sure there no undefined vareq or vartypes
C     exceptions for argument less fxs, allow ().
C
      DO I = 1,NUMOPS -1
         IF ( VARTYP(I) .EQ. '  ' .AND. TEO(I+1) .NE. '~' .AND.
     &        TEO(I) .NE. 'K' .AND. TEO(I) .NE. 'O' .AND.
     &        TEO(I) .NE. 'P' .AND. TEO(I) .NE. 'U' .AND.
     &        TEO(I) .NE. 'Y' ) THEN
            XRERR = .TRUE.
            CALL DSPERR('HIER>','Illegal consecutive operators')
            RETURN
         ENDIF
      ENDDO
      RETURN
      END
      SUBROUTINE SCASSN( INDEX, MODE, STINDX, XRERR )
      IMPLICIT NONE
C
C     the table has been evaluated leaving the final result in either
C     rvareq(1), svareq(1), or cvareq(1).  make the final assignment
C     to array which correlates to varnum(1).
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'update.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'vector.inc'
C
C     argument declarations
C
      INTEGER INDEX
      CHARACTER*4 MODE
      INTEGER STINDX
      LOGICAL XRERR
      INTEGER IST
C
C
      IF ( VARNUM(1) .EQ. 29 ) THEN
C
C        ********
C        * SHOW *
C        ********
C
C        special construct to allow (show=x) to just print out x
C
         IF ( VARTYP(1) .EQ. 'DP' ) THEN
            WRITE(6,'(A,G14.4)') ' SHOW = ',RVAREQ(1)
         ELSEIF ( VARTYP(1) .EQ. 'ST' ) THEN
            WRITE(6,'(A,A)')' SHOW = ',SVAREQ(1)(1:SVARLN(1))
         ELSEIF ( VARTYP(1) .EQ. 'DC' ) THEN
            WRITE(6,*)' SHOW = ',CVAREQ(1)
         ENDIF
      ELSEIF ( VARNUM(1) .EQ. 1 ) THEN
C
C        *******
C        *  B  *
C        *******
C
         WMAIN(INDEX)   = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 2 ) THEN
C
C        *******
C        *  Q  *
C        *******
C
         QMAIN(INDEX)   = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 3 ) THEN
C
C        *********
C        * BCOMP *
C        *********
C
         WCOMP(INDEX)   = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 4 ) THEN
C
C        *********
C        * QCOMP *
C        *********
C
         QCOMP(INDEX)   = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 5 ) THEN
C
C        **************************
C        * HARMONIC or CONSTRAINT *
C        **************************
C
         KCNSTR(INDEX)  = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 6 ) THEN
C
C        ********
C        * MASS *
C        ********
C
         AMASS(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 7 ) THEN
C
C        **********
C        * CHARGE *
C        **********
C
         CG(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 8 ) THEN
C
C        *********
C        * FBETA *
C        *********
C
         FBETA(INDEX)   = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 9 ) THEN
C
C        *******
C        *  X  *
C        *******
C
         X(INDEX)       = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 10 ) THEN
C
C        *******
C        *  Y  *
C        *******
C
         Y(INDEX)       = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 11 ) THEN
C
C        *******
C        *  Z  *
C        *******
C
         Z(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 12 ) THEN
C
C        *********
C        * XCOMP *
C        *********
C
         XCOMP(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 13 ) THEN
C
C        *********
C        * YCOMP *
C        *********
C
         YCOMP(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 14 ) THEN
C
C        *********
C        * ZCOMP *
C        *********
C
         ZCOMP(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 15 ) THEN
C
C        ********
C        *  DX  *
C        ********
C
         DX(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 16 ) THEN
C
C        ********
C        *  DY  *
C        ********
C
         DY(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 17 ) THEN
C
C        ********
C        *  DZ  *
C        ********
C
         DZ(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 18 ) THEN
C
C        ********
C        * REFX *
C        ********
C
         REFX(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 19 ) THEN
C
C        ********
C        * REFY *
C        ********
C
         REFY(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 20 ) THEN
C
C        ********
C        * REFZ *
C        ********
C
         REFZ(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .EQ. 21 ) THEN
C
C        ********
C        *  VX  *
C        ********
C
         XV(INDEX) = RVAREQ(1) * TIMFAC
      ELSEIF ( VARNUM(1) .EQ. 22 ) THEN
C
C        ********
C        *  VY  *
C        ********
C
         YV(INDEX) = RVAREQ(1) * TIMFAC
      ELSEIF ( VARNUM(1) .EQ. 23 ) THEN
C
C        ********
C        *  VZ  *
C        ********
C
         ZV(INDEX) = RVAREQ(1) * TIMFAC
      ELSEIF ( VARNUM(1) .EQ. 24 ) THEN
C
C        *********
C        * SEGID *
C        *********
C
         SEGID(INDEX) = SVAREQ(1)(1:4)
      ELSEIF ( VARNUM(1) .EQ. 25 ) THEN
C
C        *********
C        * RESID *
C        *********
C
         RESID(INDEX) = SVAREQ(1)(1:4)
      ELSEIF ( VARNUM(1) .EQ. 26 ) THEN
C
C        ***********
C        * RESNAME *
C        ***********
C
         RES(INDEX) = SVAREQ(1)(1:4)
      ELSEIF ( VARNUM(1) .EQ. 27 ) THEN
C
C        *****************
C        * IUPAC or NAME *
C        *****************
C
         TYPE(INDEX) = SVAREQ(1)(1:4)
      ELSEIF ( VARNUM(1) .EQ. 28 ) THEN
C
C        ************
C        * CHEMICAL *
C        ************
C
         IAC(INDEX) = SVAREQ(1)(1:4)
C set type-based parameter update flags each time when invoking
C this vector statement
      UPBOND=.TRUE.
      UPANGL=.TRUE.
      UPDIHE=.TRUE.
      UPIMPR=.TRUE.
      UPNBLK=.TRUE.
C
      ELSEIF ( VARNUM(1) .EQ. 50 ) THEN
C
C        ********
C        * RMSD *
C        ********
C
         RMSD(INDEX) = RVAREQ(1)
      ELSEIF ( VARNUM(1) .GE. 30 .AND. VARNUM(1) .LE. 38 ) THEN
C
C        *********************
C        * STORE# or RECALL# *
C        *********************
C
         IST = VARNUM(1) - 29
         CALL PUTSTO(RVAREQ(1),INDEX,IST,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 51 ) THEN
C
C        *********************
C        SCATTER_A1
C        *********************
C
         CALL PUTSCA('A1','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 52 ) THEN
C
C        *********************
C        SCATTER_A2
C        *********************
C
         CALL PUTSCA('A2','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 53 ) THEN
C
C        *********************
C        SCATTER_A3
C        *********************
C
         CALL PUTSCA('A3','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 54 ) THEN
C
C        *********************
C        SCATTER_A4
C        *********************
C
         CALL PUTSCA('A4','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 55 ) THEN
C
C        *********************
C        SCATTER_B1
C        *********************
C
         CALL PUTSCA('B1','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 56 ) THEN
C
C        *********************
C        SCATTER_B2
C        *********************
C
         CALL PUTSCA('B2','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 57 ) THEN
C
C        *********************
C        SCATTER_B3
C        *********************
C
         CALL PUTSCA('B3','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 58 ) THEN
C
C        *********************
C        SCATTER_B4
C        *********************
C
         CALL PUTSCA('B4','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 59 ) THEN
C
C        *********************
C        SCATTER_FP
C        *********************
C
         CALL PUTSCA('FP','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 60 ) THEN
C
C        *********************
C        SCATTER_FDP
C        *********************
C
         CALL PUTSCA('FDP','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 61 ) THEN
C
C        *********************
C        SCATTER_C
C        *********************
C
         CALL PUTSCA('C','PUT',RVAREQ(1),INDEX,XRERR)
      ELSEIF ( VARNUM(1) .EQ. 68 ) THEN
C
C        *********
C        * ALTID *
C        *********
C
         ALTID(INDEX) = SVAREQ(1)(1:1)
      ELSEIF ( VARNUM(1) .EQ. 69 ) THEN
C
C        **************************
C        ZHARMONIC 
C        **************************
C
         KZCNSTR(INDEX)  = RVAREQ(1)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PUTSCA(OP,MODE,R,INDEX,XRERR)
C
C Routine stores value R scatter array at index INDEX.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
C I/O
      CHARACTER*(*) OP, MODE
      DOUBLE PRECISION R
      INTEGER INDEX
      LOGICAL XRERR
C begin
C
      IF (HPATOF.EQ.0.OR.XRNATF.EQ.0) THEN
      CALL DSPERR ( 'PUTSCA',
     & 'SCATTer database undefined.' )
      XRERR = .TRUE.
      ELSE
      CALL PUTSC2(OP,MODE,R,INDEX,XRERR,XRNATF,HEAP(HPATOF),
     &            HEAP(HPINDF),XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP)
      END IF
      RETURN
      END
C===========================================================
      SUBROUTINE PUTSC2(OP,MODE,R,INDEX,XRERR,XRNATF,XRATOF,
     &                  XRINDF,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP)
C
C Routine stores/recalls value R scatter array at index INDEX.
C
C Warning: only minimal checking is done to make sure that the
C selection is compatible with the SCATTer definition.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'vector.inc'
C I/O
      CHARACTER*(*) OP, MODE
      DOUBLE PRECISION R
      INTEGER INDEX
      LOGICAL XRERR
      INTEGER XRNATF, XRATOF(*), XRINDF(*)
      INTEGER XRSM
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(XRSM), XRFDP(XRSM)
C local
      INTEGER IAT, IINDEX
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C we need to lookup INDEX
      IINDEX=0
      DO IAT=1,XRNATF
      IF (XRATOF(IAT).EQ.INDEX) IINDEX=IAT
      END DO
C
      IF ( IINDEX.EQ.0 ) THEN
      IF (MODE.EQ.'PUT') THEN
C
C in PUT mode print error message
      CALL DSPERR ( 'PUTSCA',
     & 'Selection is incompatible with SCATter definition.' )
      XRERR = .TRUE.
      ELSE
C
C in GET mode assign R to zero
      R=ZERO
      END IF
      ELSE
C
      IF (MODE.EQ.'PUT') THEN
C
      IF (OP.EQ.'A1') THEN
      XRSA(XRINDF(IINDEX),1)=R
      ELSEIF (OP.EQ.'A2') THEN
      XRSA(XRINDF(IINDEX),2)=R
      ELSEIF (OP.EQ.'A3') THEN
      XRSA(XRINDF(IINDEX),3)=R
      ELSEIF (OP.EQ.'A4') THEN
      XRSA(XRINDF(IINDEX),4)=R
      ELSEIF (OP.EQ.'B1') THEN
      XRSB(XRINDF(IINDEX),1)=R
      ELSEIF (OP.EQ.'B2') THEN
      XRSB(XRINDF(IINDEX),2)=R
      ELSEIF (OP.EQ.'B3') THEN
      XRSB(XRINDF(IINDEX),3)=R
      ELSEIF (OP.EQ.'B4') THEN
      XRSB(XRINDF(IINDEX),4)=R
      ELSEIF (OP.EQ.'C') THEN
      XRSC(XRINDF(IINDEX))=R
      ELSEIF (OP.EQ.'FP') THEN
      XRFP(XRINDF(IINDEX))=R
      ELSEIF (OP.EQ.'FDP') THEN
      XRFDP(XRINDF(IINDEX))=R
      END IF
C
      ELSEIF (MODE.EQ.'GET') THEN
C
      IF (OP.EQ.'A1') THEN
      R=XRSA(XRINDF(IINDEX),1)
      ELSEIF (OP.EQ.'A2') THEN
      R=XRSA(XRINDF(IINDEX),2)
      ELSEIF (OP.EQ.'A3') THEN
      R=XRSA(XRINDF(IINDEX),3)
      ELSEIF (OP.EQ.'A4') THEN
      R=XRSA(XRINDF(IINDEX),4)
      ELSEIF (OP.EQ.'B1') THEN
      R=XRSB(XRINDF(IINDEX),1)
      ELSEIF (OP.EQ.'B2') THEN
      R=XRSB(XRINDF(IINDEX),2)
      ELSEIF (OP.EQ.'B3') THEN
      R=XRSB(XRINDF(IINDEX),3)
      ELSEIF (OP.EQ.'B4') THEN
      R=XRSB(XRINDF(IINDEX),4)
      ELSEIF (OP.EQ.'C') THEN
      R=XRSC(XRINDF(IINDEX))
      ELSEIF (OP.EQ.'FP') THEN
      R=XRFP(XRINDF(IINDEX))
      ELSEIF (OP.EQ.'FDP') THEN
      R=XRFDP(XRINDF(IINDEX))
      END IF
C
      END IF
C
      END IF
C
      RETURN
      END
C
      SUBROUTINE PUTSTO(R,INDEX,IST,XRERR)
C
C Routine stores value R into store IST at index INDEX.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'vector.inc'
C I/O
      DOUBLE PRECISION R
      INTEGER INDEX, IST
      LOGICAL XRERR
C begin
      IF ( IST .LT. 1 .OR. IST .GT. 9 ) THEN
      CALL DSPERR ( 'PUTSTO','STORE non-existent' )
      XRERR = .TRUE.
      ELSEIF ( PTRSTO(IST) .EQ. 0 ) THEN
      CALL DSPERR ( 'PUTSTO','STORE pointer Corruption' )
      XRERR = .TRUE.
      ELSE
      CALL PUTST2(R,HEAP(PTRSTO(IST)),INDEX)
      ENDIF
      RETURN
      END
C
      SUBROUTINE PUTST2(R,STORE,INDEX)
C
C See routine PUTSTO above
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION R, STORE(*)
      INTEGER INDEX
C begin
      STORE(INDEX)=R
      RETURN
      END
C
      SUBROUTINE GETSTO(R,INDEX,IST,XRERR)
C
C Routine puts value of store IST at index INDEX into R.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'vector.inc'
C I/O
      DOUBLE PRECISION R
      INTEGER INDEX, IST
      LOGICAL XRERR
C begin
      IF ( IST .LT. 1 .OR. IST .GT. 9 ) THEN
      CALL DSPERR ( 'PUTSTO','STORE non-existent' )
      XRERR = .TRUE.
      ELSEIF ( PTRSTO(IST) .EQ. 0 ) THEN
      CALL DSPERR ( 'PUTSTO','STORE pointer Corruption' )
      XRERR = .TRUE.
      ELSE
      CALL GETST2(R,HEAP(PTRSTO(IST)),INDEX)
      ENDIF
      RETURN
      END
C
      SUBROUTINE GETST2(R,STORE,INDEX)
C
C See routine PUTSTO above
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION R, STORE(*)
      INTEGER INDEX
C begin
      R=STORE(INDEX)
      RETURN
      END
C
      SUBROUTINE SCDO
      IMPLICIT NONE
C
C     set aside some heap and call parsing routine for vector
C     do command
C
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
C
C     local variables
C
      INTEGER FLAGS
C
C     get heap space, and call routine to parse vector expression
C
      CALL NEXTWD('DO>')
      IF ( WD(1:4) .EQ. 'HELP' ) THEN
         CALL CNSHELP('cns-do')
         CALL CNSHELP('atom-expression')
         CALL CNSHELP('atom-object')
         CALL CNSHELP('atom-selection')
         RETURN
      ELSE
         CALL SAVEWD
      ENDIF
C
      FLAGS=ALLHP(INTEG4(NATOM))
      CALL SCPARS(HEAP(FLAGS),X,Y,Z)
      CALL FREHP(FLAGS,INTEG4(NATOM))
      RETURN
      END
      SUBROUTINE SCHVAR( K, FOUND )
      IMPLICIT NONE
C
C     check to see if variable is valid, if so get type
C
      LOGICAL FOUND
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INTEGER I, J, K, LAST
C
      FOUND = .FALSE.
C
C     determine the last character in the variable name.
C     scan through the valid names; if the given name is
C     long enough to uniquely identify the valid name then
C     compare the two
C
C this 20 is okay, since the constants are only 20 chars - WLD 980105
      DO J = 1,20
         IF ( VAREQ(K)(J:J) .NE. ' ' ) LAST = J
      ENDDO
C
      DO I = 1,SCNUMV
         IF ( LAST .GE. VUNIQ(I) ) THEN
            IF ( VAREQ(K)(1:LAST) .EQ. VNAMES(I)(1:LAST) ) THEN
               VAREQ(K) = VNAMES(I)
               FOUND = .TRUE.
               VARTYP(K) = VTYP(I)
               VARNUM(K) = VNUM(I)
               RETURN
            ENDIF
         ENDIF
      ENDDO
C
      RETURN
      END
      SUBROUTINE SCIDE2(FLAGS,X,Y,Z)
      IMPLICIT NONE
C
C
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'vector.inc'
C
C     parameter declarations
C
      INTEGER FLAGS(*)
      DOUBLE PRECISION X(*), Y(*), Z(*)
C
C     local variable declarations
C
      INTEGER COUNT
      INTEGER NSELCT
      INTEGER INDEX
      CHARACTER*4 MODE
      LOGICAL XRERR
C
      XRERR = .FALSE.
      MODE = 'REAL'
C
C
C     parse the equation into the tables
C
      CALL NEXTWD('IDENt>')
      IF ( WD(1:4) .EQ. 'HELP' ) THEN
         CALL CNSHELP('cns-identity')
         CALL CNSHELP('atom-object')
         CALL CNSHELP('atom-selection')
         RETURN
      ELSE
         CALL SAVEWD
      ENDIF
C
      CALL HIER( MODE, XRERR )
      IF ( XRERR ) THEN
         CALL DSPERR('IDENTIFY>',
     &               'Expression error, IDEN not performed.')
         RETURN
      ENDIF
C
      CALL SELCTA(FLAGS,NSELCT,X,Y,Z,.TRUE.)
C
C     no error during parsing, so perform the calculation in four steps
C     for each reflection point:
C          a.  the storing of tables is not needed here because
C              the tables never get evaluated (and altered), but
C              this routine also allocates heap space for stores
C          b.  store actual array values in parse tables
C          d.  store the value derived from step b in the array
C
      COUNT=0
      DO INDEX = 1,NATOM
C
C        step a) store or obtain parsing tables, set aside heap
C
         CALL TABMGT( INDEX, XRERR )
         IF ( XRERR ) GOTO 999
C
C        step b) fill values in for variables
C
         IF ( FLAGS(INDEX) .EQ. 1 ) THEN
             COUNT=COUNT+1
CCC             RVAREQ(1) = COUNT  !ATB 4/27/94
             RVAREQ(1) = INDEX
         ELSE
             RVAREQ(1) = 0
         ENDIF
C
C        step d) assign the table result to the array on left hand side
C
         CALL SCASSN( INDEX, MODE, INDEX, XRERR )
         IF ( XRERR ) GOTO 999
      ENDDO
C
      RETURN
C
C     print out error message and return
C
999   CALL DSPERR('IDENTIFY>','Assignment aborted.')
      WRITE(6,*)' Starting with array index: ',
     &          INDEX
      RETURN
      END
C
      SUBROUTINE SCIDEN
      IMPLICIT NONE
C
C     set aside some heap, and call routine to process the vector
C     iden expression
C
C input/output
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
C
C     local variables
C
      INTEGER FLAGS
C
C     begin
C
      FLAGS=ALLHP(INTEG4(NATOM))
      CALL SCIDE2(HEAP(FLAGS),X,Y,Z)
      CALL FREHP(FLAGS,INTEG4(NATOM))
      RETURN
      END
      SUBROUTINE SCPARS(FLAGS,X,Y,Z)
      IMPLICIT NONE
C
C     routine to do a non-recursive parse
C     this routine expects to receive an expression in parens.
C     ex.   (variable1 = variable2 * variable 3)
C
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'vector.inc'
C
C     parameter declarations
C
      INTEGER FLAGS(*)
      DOUBLE PRECISION X(*), Y(*), Z(*)
C
C     local variable declarations
C
      INTEGER NSELCT, INDEX
      INTEGER REFLCT
      CHARACTER*4 MODE
      LOGICAL XRERR
C
      XRERR = .FALSE.
C
      MODE = 'REAL'
C
C     parse the equation into the tables
C
      CALL HIER( MODE, XRERR )
      IF ( XRERR ) THEN
         CALL DSPERR('DO>',
     &               'Expression error, Calculation not performed.')
         RETURN
      ENDIF
C
      CALL SELCTA(FLAGS,NSELCT,X,Y,Z,.TRUE.)
C
C     no error during parsing, so perform the calculation in four steps
C     for each reflection point:
C          a.  on first pass save the table values and if we are
C              going to assign to a storeN set aside the heap space
C              to perform the operation, on subsequent passes restore
C              the tables from the saved values.
C          b.  store actual array values in parse tables
C          c.  evaluate the parse table to obtain a value for
C              the right hand side of the equation
C          d.  store the value derived from step b in the array
C
      INDEX = 0
      DO REFLCT = 1,NATOM
        IF ( FLAGS(REFLCT) .EQ. 1 ) THEN
           INDEX = INDEX + 1
C
Ctabtrace       WRITE(19,'(1X,A,I3)')'REFLECTION POINT: ',REFLCT
C
C         step a) store or obtain parsing tables, set aside heap
C
          CALL TABMGT( INDEX, XRERR )
          IF ( XRERR ) GOTO 999
C
C         step b) fill values in for variables
C
          CALL SCTBVL( MODE, REFLCT, INDEX, XRERR )
          IF ( XRERR ) GOTO 999
C
C         step c) evaluate the table
C
          CALL EVALTB( MODE, REFLCT, INDEX, FLAGS, XRERR )
          IF ( XRERR ) GOTO 999
C
C         step d) assign the table result to the array on left hand side
C
          CALL SCASSN( REFLCT, MODE, INDEX, XRERR )
          IF ( XRERR ) GOTO 999
        ENDIF
      ENDDO
      RETURN
C
C     print out error message and return
C
999   CALL DSPERR('SCPARS>','Assignment aborted.')
      WRITE(6,*)' Starting with array index: ', REFLCT
      RETURN
      END
      SUBROUTINE SCSHO2(FLAGS,X,Y,Z)
      IMPLICIT NONE
C
C     perform the show command
C
C
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'vector.inc'
C
C     parameter declarations
C
      INTEGER FLAGS(*)
      DOUBLE PRECISION X(*), Y(*), Z(*)
C
C     local variable declarations
C
      INTEGER NSELCT
      INTEGER I
      CHARACTER*4 SCMODE, MODE
      LOGICAL XRERR
      DOUBLE PRECISION A, A2
      INTEGER NCONF, STINDX
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
C
      XRERR = .FALSE.
      MODE = 'SHOW'
C
      CALL DECLAR( 'RESULT', 'ST', ' ', DBCOMP, DBPREC )
C
C     obtain optional show qualifier, default to element
C
      CALL NEXTWD('Show Qualifier=')
      IF ( WD(1:4) .EQ. '(   ' ) THEN
         SCMODE = 'ELEM'
         CALL SAVEWD
      ELSEIF ( WD(1:3) .EQ. 'SUM' .OR. WD(1:3) .EQ. 'MAX' .OR.
     &         WD(1:4) .EQ. 'ELEM' .OR. WD(1:3) .EQ. 'MIN' .OR.
     &         WD(1:3) .EQ. 'AVE' .OR. WD(1:4) .EQ. 'NORM' .OR.
     &         WD(1:3) .EQ. 'RMS' ) THEN
         SCMODE = WD(1:4)
      ELSEIF ( WD(1:4) .EQ. 'HELP' ) THEN
         CALL CNSHELP('cns-show')
         CALL CNSHELP('atom-object')
         CALL CNSHELP('atom-expression')
         CALL CNSHELP('atom-selection')
         RETURN
      ELSE
         CALL DSPERR('SHOW>','Illegal Qualifier.')
         RETURN
      ENDIF
C
C
C     parse the equation into the tables
C
      CALL HIER( MODE, XRERR)
      IF ( XRERR ) THEN
         CALL DSPERR('SHOW>',
     &               'Expression error, Show not performed.')
         RETURN
      ENDIF
C
      CALL SELCTA(FLAGS,NSELCT,X,Y,Z,.TRUE.)
C
C     no error during parsing, so perform the calculation in four steps
C     for each reflection point:
C          a.  on first pass save the table values, on subsequent
C              passes restore the tables from the saved values.
C          b.  store actual array values in parse tables
C          c.  evaluate the parse table to obtain a value for
C              the right hand side of the equation
C          d.  store the value derived from step b in the array
C
      STINDX = 0
      NCONF = 0
      DO I = 1,NATOM
        IF ( FLAGS(I) .EQ. 1 ) THEN
C
           STINDX = STINDX + 1
Ctabtrace       write(19,'(1x,a,i3)')'reflection point: ',I
C
C         step a) store or obtain parsing tables
C
          CALL TABMGT(  STINDX, XRERR )
          IF ( XRERR ) GOTO 999
C
C         step b) fill values in for variables
C
          CALL SCTBVL( MODE, I, STINDX, XRERR )
          IF ( XRERR ) GOTO 999
C
C         step c) evaluate the table
C
          CALL EVALTB( MODE, I, STINDX, FLAGS, XRERR )
          IF ( XRERR ) GOTO 999
C
C         step d) now have result of expression in _vareq(1)
C
          IF ( STINDX .EQ. 1 ) THEN
            IF ( VARTYP(1) .EQ. 'ST'  .AND. SCMODE.NE.'ELEM') THEN
              WRITE(6,'(3A)') ' SHOW: ',SCMODE,
     &            ' not valid for strings. Mode set to ELEM.'
              SCMODE='ELEM'
            END IF
          END IF
          IF ( SCMODE .EQ. 'ELEM' ) THEN
             IF ( STINDX .EQ. 1 ) THEN
                NCONF=0
             END IF
            IF (WRNLEV.GE.5) THEN
             IF ( VARTYP(1) .EQ. 'ST' ) THEN
              WRITE(6,1002) SEGID(I), RES(I),RESID(I),TYPE(I),
     &               SVAREQ(1)(1:4)
1002          FORMAT(' ( ',4(A,1X),')  ',A)
             ELSE
              WRITE(6,1000) SEGID(I), RES(I),RESID(I),TYPE(I),RVAREQ(1)
1000          FORMAT(' ( ',4(A,1X),')  ',G12.5)
             ENDIF
            END IF
            CALL DECLAR('RESULT',VARTYP(1),
     &           SVAREQ(1)(1:SVARLN(1)),CVAREQ(1),RVAREQ(1))
            NCONF=NCONF+1
          ELSEIF ( SCMODE(1:3) .EQ. 'SUM' ) THEN
             IF ( STINDX .EQ. 1 ) THEN
                DBPREC = 0.0D0
                NCONF=0
             END IF
             DBPREC = DBPREC + RVAREQ(1)
             NCONF  = NCONF  + 1
          ELSEIF ( SCMODE(1:3) .EQ. 'MAX' ) THEN
             IF ( STINDX .EQ. 1 ) THEN
                DBPREC = -R4BIG
                NCONF=0
             END IF
             IF ( DBPREC .LT. RVAREQ(1) ) DBPREC = RVAREQ(1)
             NCONF  = NCONF  + 1
          ELSEIF ( SCMODE(1:3) .EQ. 'MIN' ) THEN
             IF ( STINDX .EQ. 1 ) THEN
                DBPREC = R4BIG
                NCONF=0
             END IF
             IF ( DBPREC .GT. RVAREQ(1) ) DBPREC = RVAREQ(1)
             NCONF  = NCONF  + 1
          ELSEIF ( SCMODE(1:3) .EQ. 'AVE' ) THEN
             IF ( STINDX .EQ. 1 ) THEN
                DBPREC = 0.0D0
                NCONF  = 0
             ENDIF
             DBPREC = DBPREC + RVAREQ(1)
             NCONF  = NCONF  + 1
          ELSEIF ( SCMODE .EQ. 'NORM' ) THEN
             IF ( STINDX .EQ. 1 ) THEN
                DBPREC = 0.0D0
                NCONF  = 0
             ENDIF
             DBPREC = DBPREC + RVAREQ(1)**2
             NCONF = NCONF + 1
          ELSEIF ( SCMODE(1:3) .EQ. 'RMS' ) THEN
             IF ( STINDX .EQ. 1 ) THEN
                A      = 0.0D0
                A2 = 0.0D0
                NCONF  = 0
             ENDIF
             A = A + RVAREQ(1)
             A2 = A2 + RVAREQ(1)**2
             NCONF = NCONF + 1
          ENDIF
        ENDIF
      ENDDO
C
C     now output the results of all modes performed over a selection
C
      IF (NCONF.EQ.0) THEN
         IF (WRNLEV.GE.5) THEN
            WRITE(6,'(A)') ' SHOW: zero atoms selected'
         END IF
         XRERR = .TRUE.
         DBPREC = 0.0D0
         CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
      ELSE
         IF ( SCMODE(1:3) .EQ. 'SUM' ) THEN
           IF (WRNLEV.GE.5) THEN
           WRITE(6,9010)' SHOW: sum over selected elements = ',DBPREC
           END IF
9010       FORMAT ( A,F14.6 )
           CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
         ELSEIF ( SCMODE(1:3) .EQ. 'MAX' ) THEN
            IF (WRNLEV.GE.5) THEN
            WRITE(6,9010)' SHOW: maximum of selected elements = ',DBPREC
            END IF
            CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
         ELSEIF ( SCMODE(1:3) .EQ. 'MIN' ) THEN
            IF (WRNLEV.GE.5) THEN
            WRITE(6,9010)' SHOW: minimum of selected elements = ',DBPREC
            END IF
            CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
         ELSEIF ( SCMODE(1:3) .EQ. 'AVE' ) THEN
            DBPREC = DBPREC / NCONF
            IF (WRNLEV.GE.5) THEN
            WRITE(6,9010)' SHOW: average of selected elements = ',
     &                    DBPREC
            END IF
            CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
         ELSEIF ( SCMODE .EQ. 'NORM' ) THEN
            DBPREC = SQRT(DBPREC / NCONF)
            IF (WRNLEV.GE.5) THEN
            WRITE(6,9010)' SHOW: norm of selected elements = ',DBPREC
            END IF
            CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
         ELSEIF ( SCMODE(1:3) .EQ. 'RMS' ) THEN
             A = A / NCONF
             A2 = A2 / NCONF
             DBPREC=SQRT(A2-A**2)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,F14.6)')
     &' SHOW: rms deviation from aver of selected elements = ',DBPREC
      END IF
            CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
         ENDIF
      ENDIF
      RETURN
C
C     print out error message and return
C
999   CALL DSPERR('SHOW>','Assignment aborted')
      WRITE(6,*)' Starting with array index: ', I
      RETURN
      END
C
      SUBROUTINE SCSHOW
      IMPLICIT NONE
C
C     set aside some heap, and call routine to process the
C     show expression
C
C input/output
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
C
C     local variables
C
      INTEGER FLAGS
C
C     begin
C
      FLAGS=ALLHP(INTEG4(NATOM))
      CALL SCSHO2(HEAP(FLAGS),X,Y,Z)
      CALL FREHP(FLAGS,INTEG4(NATOM))
      RETURN
      END
      SUBROUTINE SCTBVL( MODE, INDEX, STINDX, XRERR )
      IMPLICIT NONE
C
C     for each variable in the parse table there is a double
C     precision, a double complex, and a string  variable (the three
C     valid data types).  the variable type tells the computation
C     routines which one is being used. in this routine get the value
C     for this index (array element 'index') and store it in the proper
C     variable, correct the variable type to match the mode of
C     calculation if needed.
C
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'vector.inc'
C
      INTEGER INDEX, STINDX
      CHARACTER*4 MODE
      LOGICAL XRERR
      INTEGER I, IST
C
C
      XRERR = .FALSE.
C     start with 2nd location, first is on left hand side of equation
      DO I = 2,NUMOPS
         IF ( VARNUM(I) .EQ. 1 ) THEN
C
C           *******
C           *  B  *
C           *******
C
            RVAREQ(I) = WMAIN(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 2 ) THEN
C
C           *******
C           *  Q  *
C           *******
C
            RVAREQ(I) = QMAIN(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 3 ) THEN
C
C           *********
C           * BCOMP *
C           *********
C
            RVAREQ(I) = WCOMP(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 4 ) THEN
C
C           *********
C           * QCOMP *
C           *********
C
            RVAREQ(I) = QCOMP(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 5 ) THEN
C
C           **************************
C           * HARMONIC or CONSTRAINT *
C           **************************
C
            RVAREQ(I) = KCNSTR(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 6 ) THEN
C
C           ********
C           * MASS *
C           ********
C
            RVAREQ(I) = AMASS(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 7 ) THEN
C
C           **********
C           * CHARGE *
C           **********
C
            RVAREQ(I) = CG(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 8 ) THEN
C
C           *********
C           * FBETA *
C           *********
C
            RVAREQ(I) = FBETA(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 9 ) THEN
C
C           *******
C           *  X  *
C           *******
C
            RVAREQ(I) = X(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 10 ) THEN
C
C           *******
C           *  Y  *
C           *******
C
            RVAREQ(I) = Y(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 11 ) THEN
C
C           *******
C           *  Z  *
C           *******
C
            RVAREQ(I) = Z(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 12 ) THEN
C
C           *********
C           * XCOMP *
C           *********
C
            RVAREQ(I) = XCOMP(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 13 ) THEN
C
C           *********
C           * YCOMP *
C           *********
C
            RVAREQ(I) = YCOMP(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 14 ) THEN
C
C           *********
C           * ZCOMP *
C           *********
C
            RVAREQ(I) = ZCOMP(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 15 ) THEN
C
C           ********
C           *  DX  *
C           ********
C
            RVAREQ(I) = DX(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 16 ) THEN
C
C           ********
C           *  DY  *
C           ********
C
            RVAREQ(I) = DY(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 17 ) THEN
C
C           ********
C           *  DZ  *
C           ********
C
            RVAREQ(I) = DZ(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 18 ) THEN
C
C           ********
C           * REFX *
C           ********
C
            RVAREQ(I) = REFX(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 19 ) THEN
C
C           ********
C           * REFY *
C           ********
C
            RVAREQ(I) = REFY(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 20 ) THEN
C
C           ********
C           * REFZ *
C           ********
C
            RVAREQ(I) = REFZ(INDEX)
         ELSEIF ( VARNUM(I) .EQ. 21 ) THEN
C
C           ********
C           *  VX  *
C           ********
C
            RVAREQ(I) = XV(INDEX) / TIMFAC
         ELSEIF ( VARNUM(I) .EQ. 22 ) THEN
C
C           ********
C           *  VY  *
C           ********
C
            RVAREQ(I) = YV(INDEX) / TIMFAC
         ELSEIF ( VARNUM(I) .EQ. 23 ) THEN
C
C           ********
C           *  VZ  *
C           ********
C
            RVAREQ(I) = ZV(INDEX) / TIMFAC
         ELSEIF ( VARNUM(I) .EQ. 24 ) THEN
C
C           *********
C           * SEGID *
C           *********
C
            SVAREQ(I)(1:4) = SEGID(INDEX)
            SVARLN(I) = 4
            CALL TRIMM(SVAREQ(I),SVARLN(I))
         ELSEIF ( VARNUM(I) .EQ. 25 ) THEN
C
C           *********
C           * RESID *
C           *********
C
            SVAREQ(I)(1:4) = RESID(INDEX)
            SVARLN(I) = 4
            CALL TRIMM(SVAREQ(I),SVARLN(I))
         ELSEIF ( VARNUM(I) .EQ. 26 ) THEN
C
C           ***********
C           * RESNAME *
C           ***********
C
            SVAREQ(I)(1:4) = RES(INDEX)
            SVARLN(I) = 4
            CALL TRIMM(SVAREQ(I),SVARLN(I))
         ELSEIF ( VARNUM(I) .EQ. 27 ) THEN
C
C           *****************
C           * IUPAC or NAME *
C           *****************
C
            SVAREQ(I)(1:4) = TYPE(INDEX)
            SVARLN(I) = 4
            CALL TRIMM(SVAREQ(I),SVARLN(I))
         ELSEIF ( VARNUM(I) .EQ. 28 ) THEN
C
C           ************
C           * CHEMICAL *
C           ************
C
            SVAREQ(I)(1:4) = IAC(INDEX)
            SVARLN(I) = 4
            CALL TRIMM(SVAREQ(I),SVARLN(I))
         ELSEIF ( VARNUM(I) .EQ. 50 ) THEN
C
C           ********
C           * RMSD *
C           ********
C
            RVAREQ(I) = RMSD(INDEX)
          ELSEIF ( VARNUM(I) .EQ. 29 ) THEN
C
C           ********
C           * SHOW *
C           ********
C
            CALL DSPERR('SCTBVL>','Illegal use of SHOW')
            XRERR = .TRUE.
            RETURN
         ELSEIF ( VARNUM(I).GE. 30 .AND. VARNUM(I) .LE. 38 ) THEN
C
C           *********************
C           * STORE# or RECALL# *
C           *********************
C
            IST = VARNUM(I) - 29
            CALL GETSTO(RVAREQ(I),INDEX,IST,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 51 ) THEN
C
C           *********************
C           SCATTER_A1
C           *********************
C
            CALL PUTSCA('A1','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 52 ) THEN
C
C           *********************
C           SCATTER_A2
C           *********************
C
            CALL PUTSCA('A2','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 53 ) THEN
C
C           *********************
C           SCATTER_A3
C           *********************
C
            CALL PUTSCA('A3','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 54 ) THEN
C
C           *********************
C           SCATTER_A4
C           *********************
C
            CALL PUTSCA('A4','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 55 ) THEN
C
C           *********************
C           SCATTER_B1
C           *********************
C
            CALL PUTSCA('B1','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 56 ) THEN
C
C           *********************
C           SCATTER_B2
C           *********************
C
            CALL PUTSCA('B2','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 57 ) THEN
C
C           *********************
C           SCATTER_B3
C           *********************
C
            CALL PUTSCA('B3','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 58 ) THEN
C
C           *********************
C           SCATTER_B4
C           *********************
C
            CALL PUTSCA('B4','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 59 ) THEN
C
C           *********************
C           SCATTER_FP
C           *********************
C
            CALL PUTSCA('FP','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 60 ) THEN
C
C           *********************
C           SCATTER_FDP
C           *********************
C
            CALL PUTSCA('FDP','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 61 ) THEN
C
C           *********************
C           SCATTER_C
C           *********************
C
            CALL PUTSCA('C','GET',RVAREQ(I),INDEX,XRERR)
         ELSEIF ( VARNUM(I) .EQ. 68 ) THEN
C
C           *********
C           * ALTID *
C           *********
C
            SVAREQ(I)(1:4) = ALTID(INDEX)
            SVARLN(I) = 1
            CALL TRIMM(SVAREQ(I),SVARLN(I))
        ENDIF
      ENDDO
C
      RETURN
      END
C=================================================================
      SUBROUTINE TABMGT( STINDX, XRERR )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'symbol.inc'
C
      LOGICAL XRERR
      INTEGER STINDX
      INTEGER I, IST
      INTEGER KVOP(MAXOPR)
      INTEGER KVRNUM(MAXOPR)
      CHARACTER*1 KTEO(MAXOPR)
      CHARACTER*2 KVRTYP(MAXOPR)
      CHARACTER*(VNAMLN) KVAREQ(MAXOPR)
      CHARACTER*4  KSVREQ(MAXOPR)
      INTEGER KSVRLN(MAXOPR)
      DOUBLE PRECISION KRVREQ(MAXOPR)
      DOUBLE COMPLEX   KCVREQ(MAXOPR)
      INTEGER KNMOPS
C
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
      SAVE KVOP, KTEO, KVRTYP, KVAREQ, KRVREQ, KCVREQ,
     &     KNMOPS, KSVREQ, KVRNUM, KSVRLN
C
C
C     on first pass if the vareq to the left of the equals sign
C     is a storeN or recallN, then set aside the heap space if needed
C
      IF ( STINDX .EQ. 1 ) THEN
         IF ( VARNUM(1) .GE. 30 .AND. VARNUM(1) .LE. 38 ) THEN
C
C           check to see if we are modifying the contents of the store
C
            IST = VARNUM(1) - 29
C
C see whethter we have to allocate space of the STORE
              IF ( PTRSTO(IST).EQ.0. OR.
     &             (PTRSTO(IST).NE.0.AND.LENSTO(IST).NE.NATOM) ) THEN
                  IF (PTRSTO(IST).NE.0) THEN
                     CALL FREHP( PTRSTO(IST), IREAL8(LENSTO(IST)))
                  END IF
                  PTRSTO(IST) = ALLHP(IREAL8(NATOM))
                  LENSTO(IST) = NATOM
                  CALL FILLR8(HEAP(PTRSTO(IST)),NATOM,ZERO)
            ENDIF
         ENDIF
         KNMOPS = NUMOPS
      ELSE
         NUMOPS = KNMOPS
      ENDIF
C
      DO I = 1,MAXOPR
C
C        on first pass we haven't modified the table yet, so
C        save it to temp variables so we can use it for the
C        other reflection points
C
C *todo* should change this to be a real scratch area instead
C        of saving all of these variables
C
         IF ( STINDX .EQ. 1 ) THEN
            KVOP(I)    = VOP(I)
            KTEO(I)    = TEO(I)
            KVRTYP(I) = VARTYP(I)
            KVAREQ(I)  = VAREQ(I)
            KRVREQ(I) = RVAREQ(I)
            KCVREQ(I) = CVAREQ(I)
            KSVREQ(I) = SVAREQ(I)(1:4)
            KSVRLN(I) = SVARLN(I)
            KVRNUM(I) = VARNUM(I)
         ELSE
C
C           evaluating the table destroys it, so for every reflection
C           point after the first, restore the table from saved values
C
            VOP(I)    = KVOP(I)
            TEO(I)    = KTEO(I)
            VARTYP(I) = KVRTYP(I)
            VAREQ(I)  = KVAREQ(I)
            RVAREQ(I) = KRVREQ(I)
            CVAREQ(I) = KCVREQ(I)
            SVAREQ(I)(1:4) = KSVREQ(I)
            SVARLN(I) = KSVRLN(I)
            VARNUM(I) = KVRNUM(I)
         ENDIF
      ENDDO
      RETURN
      END
      SUBROUTINE TBPACK( K )
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      INTEGER K, IJ
C
C     remove an entry from the parse tables
C
      NUMOPS = NUMOPS - 1
      IF ( K .NE. 1 ) THEN
         DO IJ = K,NUMOPS
            VAREQ(IJ) = VAREQ(IJ+1)
            VARTYP(IJ) = VARTYP(IJ+1)
            VOP(IJ) = VOP(IJ+1)
            TEO(IJ) = TEO(IJ+1)
            RVAREQ(IJ) = RVAREQ(IJ+1)
            CVAREQ(IJ) = CVAREQ(IJ+1)
            SVAREQ(IJ) = SVAREQ(IJ+1)
            SVARLN(IJ) = SVARLN(IJ+1)
         ENDDO
      ENDIF
      RETURN
      END
C
      SUBROUTINE VINIT
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
      INTEGER I
C
C     set values of variable names
C
      VNAMES(1) = 'B                   '
      VTYP(1)   = 'DP'
      VUNIQ(1)  = 1
      VNUM(1)   = 1
C
      VNAMES(2) = 'Q                   '
      VTYP(2)   = 'DP'
      VUNIQ(2)  = 1
      VNUM(2)   = 2
C
      VNAMES(3) = 'BCOMP               '
      VTYP(3)   = 'DP'
      VUNIQ(3)  = 4
      VNUM(3)   = 3
C
      VNAMES(4) = 'QCOMP               '
      VTYP(4)   = 'DP'
      VUNIQ(4)  = 4
      VNUM(4)   = 4
C
      VNAMES(5) = 'CONSTRAINT          '
      VTYP(5)   = 'DP'
      VUNIQ(5)  = 4
      VNUM(5)   = 5
C
      VNAMES(6) = 'HARMONIC            '
      VTYP(6)   = 'DP'
      VUNIQ(6)  = 4
      VNUM(6)   = 5
C
      VNAMES(7) = 'MASS                '
      VTYP(7)   = 'DP'
      VUNIQ(7)  = 4
      VNUM(7)   = 6
C
C
      VNAMES(8) = 'CHARGE              '
      VTYP(8)   = 'DP'
      VUNIQ(8)  = 4
      VNUM(8)   = 7
C
C
      VNAMES(9) = 'FBETA               '
      VTYP(9)   = 'DP'
      VUNIQ(9)  = 4
      VNUM(9)   = 8
C
C
      VNAMES(10)= 'X                   '
      VTYP(10)  = 'DP'
      VUNIQ(10) = 1
      VNUM(10)  = 9
C
C
      VNAMES(11)= 'Y                   '
      VTYP(11)  = 'DP'
      VUNIQ(11) = 1
      VNUM(11)  = 10
C
C
      VNAMES(12)= 'Z                   '
      VTYP(12)  = 'DP'
      VUNIQ(12) = 1
      VNUM(12)  = 11
C
C
      VNAMES(13)= 'XCOMP               '
      VTYP(13)  = 'DP'
      VUNIQ(13) = 4
      VNUM(13)  = 12
C
C
      VNAMES(14)= 'YCOMP               '
      VTYP(14)  = 'DP'
      VUNIQ(14) = 4
      VNUM(14)  = 13
C
C
      VNAMES(15)= 'ZCOMP               '
      VTYP(15)  = 'DP'
      VUNIQ(15) = 4
      VNUM(15)  = 14
C
C
      VNAMES(16)= 'DX                  '
      VTYP(16)  = 'DP'
      VUNIQ(16) = 2
      VNUM(16)  = 15
C
C
      VNAMES(17)= 'DY                  '
      VTYP(17)  = 'DP'
      VUNIQ(17) = 2
      VNUM(17)  = 16
C
C
      VNAMES(18)= 'DZ                  '
      VTYP(18)  = 'DP'
      VUNIQ(18) = 2
      VNUM(18)  = 17
C
C
      VNAMES(19)= 'REFX                '
      VTYP(19)  = 'DP'
      VUNIQ(19) = 4
      VNUM(19)  = 18
C
C
      VNAMES(20)= 'REFY                '
      VTYP(20)  = 'DP'
      VUNIQ(20) = 4
      VNUM(20)  = 19
C
C
      VNAMES(21)= 'REFZ                '
      VTYP(21)  = 'DP'
      VUNIQ(21) = 4
      VNUM(21)  = 20
C
C
      VNAMES(22)= 'VX                  '
      VTYP(22)  = 'DP'
      VUNIQ(22) = 2
      VNUM(22)  = 21
C
C
      VNAMES(23)= 'VY                  '
      VTYP(23)  = 'DP'
      VUNIQ(23) = 2
      VNUM(23)  = 22
C
C
      VNAMES(24)= 'VZ                  '
      VTYP(24)  = 'DP'
      VUNIQ(24) = 2
      VNUM(24)  = 23
C
C
      VNAMES(25)= 'SEGID               '
      VTYP(25)  = 'ST'
      VUNIQ(25) = 4
      VNUM(25)  = 24
C
C
      VNAMES(26)= 'RESID               '
      VTYP(26)  = 'ST'
      VUNIQ(26) = 4
      VNUM(26)  = 25
C
C
      VNAMES(27)= 'RESNAME             '
      VTYP(27)  = 'ST'
      VUNIQ(27) = 4
      VNUM(27)  = 26
C
      VNAMES(28)= 'IUPAC               '
      VTYP(28)  = 'ST'
      VUNIQ(28) = 4
      VNUM(28)  = 27
C
C
      VNAMES(29)= 'NAME                '
      VTYP(29)  = 'ST'
      VUNIQ(29) = 4
      VNUM(29)  = 27
C
C
      VNAMES(30)= 'CHEMICAL            '
      VTYP(30)  = 'ST'
      VUNIQ(30) = 4
      VNUM(30)  = 28
C
C
      VNAMES(31)= 'SHOW                '
      VTYP(31)  = 'VT'
      VUNIQ(31) = 4
      VNUM(31)  = 29
C
      VNAMES(32)= 'STORE1              '
      VTYP(32)  = 'DP'
      VUNIQ(32) = 6
      VNUM(32)  = 30
C
C
      VNAMES(33)= 'STORE2              '
      VTYP(33)  = 'DP'
      VUNIQ(33) = 6
      VNUM(33)  = 31
C
C
      VNAMES(34)= 'STORE3              '
      VTYP(34)  = 'DP'
      VUNIQ(34) = 6
      VNUM(34)  = 32
C
C
      VNAMES(35)= 'STORE4              '
      VTYP(35)  = 'DP'
      VUNIQ(35) = 6
      VNUM(35)  = 33
C
      VNAMES(36)= 'STORE5              '
      VTYP(36)  = 'DP'
      VUNIQ(36) = 6
      VNUM(36)  = 34
C
C
      VNAMES(37)= 'STORE6              '
      VTYP(37)  = 'DP'
      VUNIQ(37) = 6
      VNUM(37)  = 35
C
C
      VNAMES(38)= 'STORE7              '
      VTYP(38)  = 'DP'
      VUNIQ(38) = 6
      VNUM(38)  = 36
C
C
      VNAMES(39)= 'STORE8              '
      VTYP(39)  = 'DP'
      VUNIQ(39) = 6
      VNUM(39)  = 37
C
C
      VNAMES(40)= 'STORE9              '
      VTYP(40)  = 'DP'
      VUNIQ(40) = 6
      VNUM(40)  = 38
C
C
      VNAMES(41)= 'RECALL1             '
      VTYP(41)  = 'DP'
      VUNIQ(41) = 6
      VNUM(41)  = 30
C
C
      VNAMES(42)= 'RECALL2             '
      VTYP(42)  = 'DP'
      VUNIQ(42) = 6
      VNUM(42)  = 31
C
C
      VNAMES(43)= 'RECALL3             '
      VTYP(43)  = 'DP'
      VUNIQ(43) = 6
      VNUM(43)  = 32
C
C
      VNAMES(44)= 'RECALL4             '
      VTYP(44)  = 'DP'
      VUNIQ(44) = 6
      VNUM(44)  = 33
C
      VNAMES(45)= 'RECALL5             '
      VTYP(45)  = 'DP'
      VUNIQ(45) = 6
      VNUM(45)  = 34
C
C
      VNAMES(46)= 'RECALL6             '
      VTYP(46)  = 'DP'
      VUNIQ(46) = 6
      VNUM(46)  = 35
C
C
      VNAMES(47)= 'RECALL7             '
      VTYP(47)  = 'DP'
      VUNIQ(47) = 6
      VNUM(47)  = 36
C
C
      VNAMES(48)= 'RECALL8             '
      VTYP(48)  = 'DP'
      VUNIQ(48) = 6
      VNUM(48)  = 37
C
C
      VNAMES(49)= 'RECALL9             '
      VTYP(49)  = 'DP'
      VUNIQ(49) = 6
      VNUM(49)  = 38
C
      VNAMES(50)= 'RMSD                '
      VTYP(50)  = 'DP'
      VUNIQ(50) = 4
      VNUM(50)  = 50
C
      VNAMES(51)= 'SCATTER_A1          '
      VTYP(51)  = 'DP'
      VUNIQ(51) = 10
      VNUM(51)  = 51
C
      VNAMES(52)= 'SCATTER_A2          '
      VTYP(52)  = 'DP'
      VUNIQ(52) = 10
      VNUM(52)  = 52
C
      VNAMES(53)= 'SCATTER_A3          '
      VTYP(53)  = 'DP'
      VUNIQ(53) = 10
      VNUM(53)  = 53
C
      VNAMES(54)= 'SCATTER_A4          '
      VTYP(54)  = 'DP'
      VUNIQ(54) = 10
      VNUM(54)  = 54
C
      VNAMES(55)= 'SCATTER_B1          '
      VTYP(55)  = 'DP'
      VUNIQ(55) = 10
      VNUM(55)  = 55
C
      VNAMES(56)= 'SCATTER_B2          '
      VTYP(56)  = 'DP'
      VUNIQ(56) = 10
      VNUM(56)  = 56
C
      VNAMES(57)= 'SCATTER_B3          '
      VTYP(57)  = 'DP'
      VUNIQ(57) = 10
      VNUM(57)  = 57
C
      VNAMES(58)= 'SCATTER_B4          '
      VTYP(58)  = 'DP'
      VUNIQ(58) = 10
      VNUM(58)  = 58
C
      VNAMES(59)= 'SCATTER_FP          '
      VTYP(59)  = 'DP'
      VUNIQ(59) = 10
      VNUM(59)  = 59
C
      VNAMES(60)= 'SCATTER_FDP         '
      VTYP(60)  = 'DP'
      VUNIQ(60) = 11
      VNUM(60)  = 60
C
      VNAMES(61)= 'SCATTER_C           '
      VTYP(61)  = 'DP'
      VUNIQ(61) = 9
      VNUM(61)  = 61
C
      VNAMES(62)= 'NOE_DISTANCE        '
      VTYP(62)  = 'DP'
      VUNIQ(62) = 5
      VNUM(62)  = 62
C
      VNAMES(63)= 'NOE_LOWER-ERROR     '
      VTYP(63)  = 'DP'
      VUNIQ(63) = 5
      VNUM(63)  = 63
C
      VNAMES(64)= 'NOE_HIGHER-ERROR    '
      VTYP(64)  = 'DP'
      VUNIQ(64) = 5
      VNUM(64)  = 64
C
      VNAMES(65)= 'NOE_VOLUME          '
      VTYP(65)  = 'DP'
      VUNIQ(65) = 5
      VNUM(65)  = 65
C
      VNAMES(66)= 'NOE_WEIGHT          '
      VTYP(66)  = 'DP'
      VUNIQ(66) = 5
      VNUM(66)  = 66
C
      VNAMES(67)= 'NOE_TEST            '
      VTYP(67)  = 'DP'
      VUNIQ(67) = 5
      VNUM(67)  = 67
C
      VNAMES(68)= 'ALTID               '
      VTYP(68)  = 'ST'
      VUNIQ(68) = 4
      VNUM(68)  = 68
C
      VNAMES(69)= 'ZHARMONIC           '
      VTYP(69)  = 'DP'
      VUNIQ(69) = 4
      VNUM(69)  = 69
C
      SCNUMV = 69
C
C obsolete variable XRNUMV
      XRNUMV = 0
C
      DO I = 1,MAXOPR
         VOP(I) = 0
         TEO(I) = ' '
         VAREQ(I) = '                    '
         RVAREQ(I) = 0.0D0
         VARTYP(I) = '  '
         SVAREQ(I) = ' '
         SVARLN(I) = 0
         CVAREQ(I) = DCMPLX(0.0D0,0.0D0)
         VARNUM(I) = 0
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE TBOUT
      IMPLICIT NONE
C
C     this routine is used as a debugging tool only
C
      INCLUDE 'cns.inc'
      INCLUDE 'vector.inc'
C
      INTEGER I
C
C
      WRITE(19,100)
100   FORMAT(1X,'    VOP TEO VAREQ    TYPE        REAL          COMP')
      DO I = 1,NUMOPS
      WRITE(19,'(1X,I2,2X,I3,2X,A1,2X,A8,2X,A2,2X,F14.5,2(F14.4))')
     &I,VOP(I),TEO(I),VAREQ(I),VARTYP(I),RVAREQ(I),CVAREQ(I)
      ENDDO
      RETURN
      END
      SUBROUTINE EVAL
      IMPLICIT NONE
C
C
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'timer.inc'
C
C     local variable declarations
C
      INTEGER IDUMMY(1), IL1
      CHARACTER*4 MODE
      LOGICAL XRERR
C
      XRERR = .FALSE.
      MODE = 'EVAL'
C
C     look for help request
C
      CALL NEXTWD('EVALUATE>')
      IF ( WD(1:4) .EQ. 'HELP' ) THEN
         CALL CNSHELP('cns-evaluate')
         RETURN
      ELSE
         CALL SAVEWD
      ENDIF
C
C     parse the equation into the tables
C
      CALL HIER( MODE, XRERR )
      IF ( XRERR ) THEN
         CALL DSPERR('EVALUATE>',
     &               'Expression error, EVALUATE not performed.')
         RETURN
      ENDIF
C
C     no error during parsing, so perform the evaluation
C
C     step a) evaluate the table
C
      CALL EVALTB( MODE, 1, 1, IDUMMY, XRERR )
      IF ( XRERR ) GOTO 999
C
C     step b) assign the value to the variable
C
      CALL LDECLAR( VAREQ(1)(2:VNAMLN), VARTYP(1),
     &            SVAREQ(1)(1:SVARLN(1)), CVAREQ(1), RVAREQ(1) )
C
      IF (WRNLEV.GE.5) THEN
C
C     echo the assignment
C
      IL1=VNAMLN
      CALL TRIMM(VAREQ(1),IL1)
      IF (VARTYP(1).EQ.'ST') THEN
      WRITE(6,'(5A)')
     &     ' EVALUATE: symbol ',VAREQ(1)(1:IL1),' set to "',
     &      SVAREQ(1)(1:SVARLN(1)),'" (string)'
      ELSE IF (VARTYP(1).EQ.'DP') THEN
      WRITE(6,'(3A,G14.6,A)')
     &     ' EVALUATE: symbol ',VAREQ(1)(1:IL1),' set to ',
     &      RVAREQ(1),' (real)'
      ELSE IF (VARTYP(1).EQ.'DC') THEN
      WRITE(6,'(3A,G14.6,A,G14.6,A)')
     &     ' EVALUATE: symbol ',VAREQ(1)(1:IL1),' set to (',
     &     DBLE(CVAREQ(1)),',',DIMAG(CVAREQ(1)),') (complex)'
      ELSEIF (VARTYP(1).EQ.'LO') THEN
      WRITE(6,'(5A)')
     &     ' EVALUATE: symbol ',VAREQ(1)(1:IL1),' set to ',
     &      SVAREQ(1)(1:SVARLN(1)),' (logical)'
      END IF
      END IF
C
      RETURN
C
C     print out error message and return
C
999   CALL DSPERR('EVALUATE>','EVALUATE assignment aborted.')
      RETURN
      END
C
C =================================================================
C

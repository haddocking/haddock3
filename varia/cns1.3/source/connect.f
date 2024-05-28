      SUBROUTINE CONNEC
C
C Determines connectivity of a specified set of atoms by
C making a recursive search on the chemical bond graph given by
C IB(i),JB(i),i=1,NBOND
C
C The total number of sets is stored in symbol $mset.  Each
C set is assigned a unique id (1,...,mset).
C The individual set ids are optionally stored in the specified
C atom object.
C The size of each set is stored in the specified atom object.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
C local
      INTEGER MAXLEV
C parameter
      INTEGER IATBMX
      PARAMETER (IATBMX=12)
C pointer
      INTEGER NATBON, IATBON, SUBSET
      INTEGER LIST, VERTEX, IND, RETADR
      INTEGER ASETIDTO, ANSETTO
C begin
      MAXLEV=MAX(2,NBOND)
      NATBON=ALLHP(INTEG4(NATOM))
      IATBON=ALLHP(INTEG4(NATOM*IATBMX))
      SUBSET=ALLHP(INTEG4(NATOM))
      LIST=ALLHP(INTEG4(NATOM))
      VERTEX=ALLHP(INTEG4(MAXLEV))
      IND=ALLHP(INTEG4(MAXLEV))
      RETADR=ALLHP(INTEG4(MAXLEV))
      ASETIDTO=ALLHP(IREAL8(NATOM))
      ANSETTO=ALLHP(IREAL8(NATOM))
      CALL CONNE2(HEAP(NATBON),IATBMX,HEAP(IATBON),
     &             HEAP(SUBSET),HEAP(LIST),MAXLEV,
     &             HEAP(VERTEX),HEAP(IND),HEAP(RETADR),
     &             HEAP(ASETIDTO),HEAP(ANSETTO))
      CALL FREHP(ANSETTO,IREAL8(NATOM))
      CALL FREHP(ASETIDTO,IREAL8(NATOM))
      CALL FREHP(RETADR,INTEG4(MAXLEV))
      CALL FREHP(IND,INTEG4(MAXLEV))
      CALL FREHP(VERTEX,INTEG4(MAXLEV))
      CALL FREHP(LIST,INTEG4(NATOM))
      CALL FREHP(SUBSET,INTEG4(NATOM))
      CALL FREHP(IATBON,INTEG4(NATOM*IATBMX))
      CALL FREHP(NATBON,INTEG4(NATOM))
C
      RETURN
      END
C=====================================================================
      SUBROUTINE CONNE2(NATBON,IATBMX,IATBON,
     &                   SUBSET,LIST,MAXLEV,VERTEX,IND,RETADR,
     &                   ASETIDTO,ANSETTO)
C
C see CONNEC above
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INTEGER NATBON(*),IATBMX,IATBON(IATBMX,*)
      INTEGER SUBSET(*)
      INTEGER LIST(*), MAXLEV, VERTEX(*), IND(*), RETADR(*)
      DOUBLE PRECISION ASETIDTO(*), ANSETTO(*)
C local
      INTEGER  I, IBT, JBT, MSET, LEVEL, PASSED, NEXT, M
      INTEGER NSUBSE, SIZE
      LOGICAL QPRINT
      CHARACTER*4 SI, RI, RE, AT
      CHARACTER*(WDMAX) SETIDTO, NSETTO
      INTEGER LSETIDTO, LNSETTO
      DOUBLE PRECISION DBPREC
      DOUBLE COMPLEX DBCOMP
C parameter
      INTEGER MARK
      PARAMETER (MARK=-1)
C begin
C
C defaults
      CALL FILL4(SUBSET,NATOM,1)
      NSUBSE=NATOM
      QPRINT=.FALSE.
      SETIDTO=' '
      NSETTO=' '
C
C command parsing
      CALL PUSEND('CONNECT>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('CONNECT>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-connectivity')
C
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(SUBSET,NSUBSE,X,Y,Z,.TRUE.)
      ELSE IF (WD(1:4).EQ.'PRIN') THEN
      CALL NEXTLO('PRINT=',QPRINT)
      ELSE IF (WD(1:WDLEN).EQ.'SETIDTO') THEN
      CALL NEXTST('SETIDTO=',WD)
      CALL COPYST(SETIDTO,WDMAX,LSETIDTO,WD,WDLEN)
      ELSE IF (WD(1:WDLEN).EQ.'NSETTO') THEN
      CALL NEXTST('NSETTO=',WD)
      CALL COPYST(NSETTO,WDMAX,LNSETTO,WD,WDLEN)
      ELSE
      CALL CHKEND('CONNECT>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
C first make list of all bonds of all atoms (vertices)
      DO I=1,NATOM
      NATBON(I)=0
      END DO
      DO I=1,NBOND
      IBT=IB(I)
      JBT=JB(I)
      NATBON(IBT)=NATBON(IBT)+1
      NATBON(JBT)=NATBON(JBT)+1
      IF (NATBON(IBT).GT.IATBMX.OR.NATBON(JBT).GT.IATBMX) THEN
      CALL WRNDIE(-1,'<CONNEC>','IATBMX exceeded')
      END IF
      IATBON(NATBON(JBT),JBT)=IBT
      IATBON(NATBON(IBT),IBT)=JBT
      END DO
C
C do search in bond graph
      DO I=1,NATOM
      LIST(I)=MARK
      END DO
      LEVEL=1
      MSET=0
      I=1
6363  IF ( .NOT. (I.LE.NATOM) ) GOTO 6364
      IF (.NOT.(LIST(I).EQ.MARK.AND.SUBSET(I).EQ.1))GOTO 6365
      MSET=MSET+1
C
C *invoke PROCEDURE SEARCH (PASSED:=VERTEX)
      PASSED=I
      RETADR(LEVEL)=1
      GOTO 2000
2001  CONTINUE
C *return address
C
C
6365  CONTINUE
      I=I+1
      GOTO 6363
6364  CONTINUE
      GOTO 2999
C
C---- PROCEDURE SEARCH ( CALL BY VALUE: PASSED->VERTEX ) --------------
C     LOCAL: IND
C *head
2000  CONTINUE
      LEVEL=LEVEL+1
      IF (LEVEL.GT.MAXLEV) CALL WRNDIE(-2,'<CONNE2>','MAXLEV exceeded')
      VERTEX(LEVEL)=PASSED
C
C *begin
      LIST(VERTEX(LEVEL))=MSET
C
C * *loop head: FOR { IND(LEVEL)=1,NATBON(VERTEX(LEVEL) }
      IND(LEVEL)=0
7711  CONTINUE
C * *loop iteration
      IND(LEVEL)=IND(LEVEL)+1
      IF (IND(LEVEL).GT.NATBON(VERTEX(LEVEL))) GOTO 7722
C * *loop begin
C
      NEXT=IATBON(IND(LEVEL),VERTEX(LEVEL))
      IF (LIST(NEXT).NE.-1.AND.LIST(NEXT).NE.MSET) THEN
      CALL WRNDIE(0,'<CONNE2>','check algorithm')
      END IF
CCCCC      IF (LIST(NEXT).EQ.MARK.AND.SUBSET(NEXT).EQ.1) THEN
      IF (.NOT.(LIST(NEXT).EQ.MARK.AND.SUBSET(NEXT).EQ.1))GOTO 3567
C
C * * *invoke procedure SEARCH ( PASSED:=NEXT )
      PASSED=NEXT
      RETADR(LEVEL)=2
      GOTO 2000
2002  CONTINUE
C * * *return label
C
CCCCC     END IF
3567  CONTINUE
C
C * *end loop
      GOTO 7711
7722  CONTINUE
C * *exit loop
C
C *return to address
      LEVEL=LEVEL-1
      IF (LEVEL.LE.0) CALL WRNDIE(-1,'<CONNE2>','Level underflow')
      IF (RETADR(LEVEL).EQ.2) GOTO 2002
      IF (RETADR(LEVEL).EQ.1) GOTO 2001
      CALL WRNDIE(-2,'<CONNE2>','Unknown return address, check code')
C--END PROCEDURE SEARCH-----------------------------------------------
C
2999  CONTINUE
C
C summary statement
      WRITE(6,'(A,I6,A)')
     &' CONNECt: selected atoms form ',MSET,
     &' covalently disconnected set(s)'
C
C symbol
      DBPREC=MSET
      CALL DECLAR( 'MSET', 'DP', ' ', DBCOMP, DBPREC )
C
C extensive info
      IF (QPRINT) THEN
      DO M=1,MSET
      WRITE(6,'(A,I6,A)')
      WRITE(6,'(A,I6,A)') ' CONNECt: list of all atoms in set ',M
      DO I=1,NATOM
      IF (LIST(I).EQ.M) THEN
      CALL ATOMID(I,SI,RI,RE,AT)
      WRITE(6,'(11A)')
     & '        ','    ','  (',SI,' ',RI,' ',RE,' ',AT,')'
      END IF
      END DO
      END DO
      END IF
C
C fill atom object with set ids and set sizes if required
      IF (SETIDTO.NE.' '.OR.NSETTO.NE.' ') THEN
C
C initialize temporary arrays
      DO I=1,NATOM
      ASETIDTO(I)=0
      ANSETTO(I)=0
      END DO
C
C loop through all sets
      DO M=1,MSET
C
C determine size of set
      SIZE=0
      DO I=1,NATOM
      IF (LIST(I).EQ.M) THEN
      SIZE=SIZE+1
      END IF
      END DO
C
C set id is M and set size is SIZE.  Now fill the temporary arrays
      DO I=1,NATOM
      IF (LIST(I).EQ.M) THEN
      ASETIDTO(I)=M
      ANSETTO(I)=SIZE
      END IF
      END DO
C
      END DO
C
C now copy the temporary arrays into specified atom object
      IF (SETIDTO.NE.' ') THEN
      CALL CPATPR(SETIDTO,ASETIDTO)
      END IF
      IF (NSETTO.NE.' ') THEN
      CALL CPATPR(NSETTO,ANSETTO)
      END IF
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE CPATPR(WD,ARRAY)
C
C routine copies array into specified atom object
C
C all atom objects are supported except the scatter_fp
C and scatter_fdp arrays
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      CHARACTER*(*) WD
      DOUBLE PRECISION ARRAY(*)
C local
      INTEGER I, II
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      IF (WD(1:4).EQ.'X   ') THEN
      CALL COPYR8(ARRAY,X,NATOM)
      ELSE IF (WD(1:4).EQ.'Y   ') THEN
      CALL COPYR8(ARRAY,Y,NATOM)
      ELSE IF (WD(1:4).EQ.'Z   ') THEN
      CALL COPYR8(ARRAY,Z,NATOM)
      ELSE IF (WD(1:4).EQ.'B   ') THEN
      CALL COPYR8(ARRAY,WMAIN,NATOM)
      ELSE IF (WD(1:4).EQ.'Q   ') THEN
      CALL COPYR8(ARRAY,QMAIN,NATOM)
      ELSE IF (WD(1:4).EQ.'XCOM') THEN
      CALL COPYR8(ARRAY,XCOMP,NATOM)
      ELSE IF (WD(1:4).EQ.'YCOM') THEN
      CALL COPYR8(ARRAY,YCOMP,NATOM)
      ELSE IF (WD(1:4).EQ.'ZCOM') THEN
      CALL COPYR8(ARRAY,ZCOMP,NATOM)
      ELSE IF (WD(1:4).EQ.'BCOM') THEN
      CALL COPYR8(ARRAY,WCOMP,NATOM)
      ELSE IF (WD(1:4).EQ.'QCOM') THEN
      CALL COPYR8(ARRAY,QCOMP,NATOM)
      ELSE IF (WD(1:4).EQ.'REFX') THEN
      CALL COPYR8(ARRAY,REFX,NATOM)
      ELSE IF (WD(1:4).EQ.'REFY') THEN
      CALL COPYR8(ARRAY,REFY,NATOM)
      ELSE IF (WD(1:4).EQ.'REFZ') THEN
      CALL COPYR8(ARRAY,REFZ,NATOM)
      ELSE IF (WD(1:4).EQ.'RMSD') THEN
      CALL COPYR8(ARRAY,RMSD,NATOM)
      ELSE IF (WD(1:4).EQ.'HARM') THEN
      CALL COPYR8(ARRAY,KCNSTR,NATOM)
      ELSE IF (WD(1:4).EQ.'ZHAR') THEN
      CALL COPYR8(ARRAY,KZCNSTR,NATOM)
      ELSE IF (WD(1:4).EQ.'MASS') THEN
      CALL COPYR8(ARRAY,AMASS,NATOM)
      ELSE IF (WD(1:4).EQ.'CHAR') THEN
      CALL COPYR8(ARRAY,CG,NATOM)
      ELSE IF (WD(1:4).EQ.'DX  ') THEN
      CALL COPYR8(ARRAY,DX,NATOM)
      ELSE IF (WD(1:4).EQ.'DY  ') THEN
      CALL COPYR8(ARRAY,DY,NATOM)
      ELSE IF (WD(1:4).EQ.'DZ  ') THEN
      CALL COPYR8(ARRAY,DZ,NATOM)
C
      ELSE IF (WD(1:4).EQ.'VX  ') THEN
      DO I=1,NATOM
      XV(I)=XV(I)/TIMFAC
      END DO
      CALL COPYR8(ARRAY,XV,NATOM)
      DO I=1,NATOM
      XV(I)=XV(I)*TIMFAC
      END DO
      ELSE IF (WD(1:4).EQ.'VY  ') THEN
      DO I=1,NATOM
      YV(I)=YV(I)/TIMFAC
      END DO
      CALL COPYR8(ARRAY,YV,NATOM)
      DO I=1,NATOM
      YV(I)=YV(I)*TIMFAC
      END DO
      ELSE IF (WD(1:4).EQ.'VZ  ') THEN
      DO I=1,NATOM
      ZV(I)=ZV(I)/TIMFAC
      END DO
      CALL COPYR8(ARRAY,ZV,NATOM)
      DO I=1,NATOM
      ZV(I)=ZV(I)*TIMFAC
      END DO
      ELSE IF (WD(1:4).EQ.'FBET') THEN
      CALL COPYR8(ARRAY,FBETA,NATOM)
C
      ELSE IF (WD(1:5).EQ.'STORE') THEN
      IF (WD(6:6).EQ.'1') THEN
      II=1
      ELSEIF (WD(6:6).EQ.'2') THEN
      II=2
      ELSEIF (WD(6:6).EQ.'3') THEN
      II=3
      ELSEIF (WD(6:6).EQ.'4') THEN
      II=4
      ELSEIF (WD(6:6).EQ.'5') THEN
      II=5
      ELSEIF (WD(6:6).EQ.'6') THEN
      II=6
      ELSEIF (WD(6:6).EQ.'7') THEN
      II=7
      ELSEIF (WD(6:6).EQ.'8') THEN
      II=8
      ELSEIF (WD(6:6).EQ.'9') THEN
      II=9
      ELSE
      CALL DSPERR('CPATPR',' invalid atom object for this command')
      END IF
C
      IF ( PTRSTO(II).EQ.0) THEN
      PTRSTO(II) = ALLHP(IREAL8(NATOM))
      LENSTO(II) = NATOM
      CALL FILLR8(HEAP(PTRSTO(II)),NATOM,ZERO)
      ENDIF
      CALL COPYR8(ARRAY,HEAP(PTRSTO(II)),NATOM)
C
      ELSE
      CALL DSPERR('CPATPR',' invalid atom object for this command')
      END IF
      RETURN
      END

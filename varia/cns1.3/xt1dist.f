C     ------------------------------------------------------------------
      SUBROUTINE XT1DIST
C     ------------------------------------------------------------------
C     Reads in relaxation rate enhancement information
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE  "cns.inc"
      INCLUDE  "comand.inc"
      INCLUDE  "funct.inc"
      INCLUDE  "mtf.inc"
      INCLUDE  "heap.inc"
      INCLUDE  "numbers.inc"

      INTEGER SPTR
      CHARACTER KT1INPUT*132,KT1OUTPUT*132,KPDB*132
      DOUBLE PRECISION KT1

      SPTR=ALLHP(INTEG4(NATOM))
      CALL PUSEND('XT1D>')
862   CONTINUE
      CALL NEXTWD('XT1D>')
      CALL MISCOM('XT1D>',USED)
      IF (.NOT.USED) THEN
C     ------------------------------------------------------------------
C     Documentation
C     ------------------------------------------------------------------
          IF (WD(1:4).EQ.'HELP') THEN
              WRITE(DUNIT,'(10X,A)')
     &              ' XT1D {<T1-STATEMENT>} END ',
     &              ' <T1-STATEMENT>:== ',
     &              ' KINPut <NAME> <REAL> <NAME>',
     &              ' * InputFile Kvalue OutputFile *',
     &              ' STRUcture <NAME> <NAME> <NAME>',
     &              ' * PDBfile InputFile OutputFile *'
C     ------------------------------------------------------------------
C     Restraint conversion using a given value of K
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'KINP') THEN
              CALL NEXTFI('Input file for T1 =',KT1INPUT)
              PRINT*,'Using as input file for T1:',KT1INPUT
              CALL NEXTF('Value of K =',KT1)
              PRINT*,'Using as value of K:',KT1
              CALL NEXTFI('Output NOE-like file for T1 =',KT1OUTPUT)
              PRINT*,'Using as output NOE-like file for T1:',KT1OUTPUT
              CALL KT1TYPECALC(KT1INPUT,KT1,KT1OUTPUT)
C     ------------------------------------------------------------------
C     Restraint conversion using a structure
C     ------------------------------------------------------------------
          ELSE IF (WD(1:4).EQ.'STRU') THEN
              CALL NEXTFI('Input PDB file for T1 =',KPDB)
              PRINT*,'Using as input PDB file for T1:',KPDB
              CALL NEXTFI('Input file for T1 =',KT1INPUT)
              PRINT*,'Using as input file for T1:',KT1INPUT
              CALL NEXTFI('Output NOE-like file for T1 =',KT1OUTPUT)
              PRINT*,'Using as output NOE-like file for T1:',KT1OUTPUT
              CALL KT1PDB(KT1INPUT,KPDB,KT1OUTPUT)
C     ------------------------------------------------------------------
C     Check for END statement
C     ------------------------------------------------------------------
          ELSE
              CALL CHKEND('XT1D>',DONE)
          END IF
      END IF
      IF (.NOT.(DONE)) GOTO 862
      DONE=.FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))

      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE KT1TYPECALC (KT1INPUT,KT1,KT1OUTPUT)
C     ------------------------------------------------------------------
C     Converts T1 data to distance restraints using a given value of K
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      CHARACTER KT1INPUT*132,KT1OUTPUT*132
      CHARACTER BEGIN*3
      DOUBLE PRECISION KT1,R,RNOE
      INTEGER N1A,N2A
      CHARACTER*4 R1A,R2A
      INTEGER NLINES

      OPEN(21,FILE=KT1INPUT,STATUS='OLD',ERR=666)
      OPEN(31,FILE=KT1OUTPUT)

      NLINES=0
      DO I=1,10000 
10    READ(21,91,END=20,ERR=21) BEGIN
      NLINES=NLINES+1
      IF (BEGIN.NE.'   ') THEN
          IF (BEGIN.EQ.'KT1') THEN
              BACKSPACE 21
              READ(21,92) BEGIN,N1A,R1A,N2A,R2A,R
              RNOE=((KT1/R)**(1/6.))+1
              WRITE(31,94) ' assign (resid ',N1A,' and name ',R1A,
     &                     ') (resid ',N2A,' and name ',R2A,') ',
     &                     RNOE,RNOE,' 0'
          END IF
      END IF
20    CONTINUE        
      END DO
21    CONTINUE
      CLOSE (21)
      CLOSE (31)
      
91    FORMAT (A3)
92    FORMAT (A3,1X,I4,1X,A4,1X,I4,1X,A4,1X,F6.2)
94    FORMAT (A15,I4,A10,A4,A9,I4,A10,A4,A2,F5.1,F5.1,A2)      

      RETURN
666   PRINT*,'The file:',KT1INPUT,'does not exist...'
      RETURN
      END
C     ------------------------------------------------------------------
      SUBROUTINE KT1PDB (KT1INPUT,KPDB,KT1OUTPUT)
C     ------------------------------------------------------------------
C     Converts T1 data to distance restraints using a given structure
C
C     By Gabriele Cavallaro, Andrea Giachetti and Giacomo Parigi (2003)
C     ------------------------------------------------------------------
      CHARACTER KT1INPUT*132,KT1OUTPUT*132,KPDB*132
      CHARACTER BEGIN*4,BEGIN1*6
      DOUBLE PRECISION KT1,R(2000),RNOE,RMAX
      DOUBLE PRECISION X1,Y1,Z1,RA,X2,Y2,Z2
      INTEGER N1A(2000),N2A(2000)
      CHARACTER*4 R1A(2000),R2A(2000)
      INTEGER NLINES,NUMBERATOM,NORES,BB
      CHARACTER NAMEATOM*4

      OPEN(21,FILE=KT1INPUT,STATUS='OLD',ERR=666)

      RMAX=0
      NLINES=0
      DO I=1,10000
10    READ(21,91,END=20,ERR=21) BEGIN
      IF (BEGIN.NE.'    ') THEN
          IF (BEGIN.EQ.'KT1 ') THEN
              NLINES=NLINES+1
              BACKSPACE 21
              READ(21,92) BEGIN,N1A(NLINES),R1A(NLINES),N2A(NLINES),
     &                    R2A(NLINES),R(NLINES)
          END IF
      END IF
20    CONTINUE
      END DO
21    CONTINUE
      CLOSE(21)
    
      K=0
      DO J=1,NLINES
         OPEN(11,FILE=KPDB,STATUS='OLD',ERR=667)
         BB=0
         DO I=1,10000
            READ(11,91,END=102) BEGIN1
            IF (BEGIN1.NE.'      ') THEN
                IF (BEGIN1.EQ.'ATOM  ') THEN 
                    BACKSPACE 11  
                    READ(11,'(5X,I6,1X,A4,6X,I4,4X,3F8.3)')
     &                   NUMBERATOM,NAMEATOM,NORES,XA,YA,ZA
                    IF (NORES.EQ.N1A(J)) THEN
                        IF (NAMEATOM.EQ.R1A(J)) THEN 
                            X1=XA
                            Y1=YA
                            Z1=ZA
                            BB=BB+1
                        END IF
                    END IF
                    IF (NORES.EQ.N2A(J)) THEN
                        IF (NAMEATOM.EQ.R2A(J)) THEN
                            X2=XA
                            Y2=YA
                            Z2=ZA
                            BB=BB+1
                        END IF
                    END IF
                END IF
            END IF
         END DO
102      CONTINUE
         IF (BB.NE.2) THEN
             PRINT*,'%T1DIST-ERR: Mismatch in PDB file at line:',J
             STOP
         END IF
         RA=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
         K=(RA**6)*R(J)
         IF (K.GT.KMAX) THEN
             KMAX=K
         END IF
         CLOSE (11)
      END DO
      PRINT*,'Value of KMAX:',KMAX

      OPEN(31,FILE=KT1OUTPUT)
      DO J=1,NLINES
         RNOE=(KMAX/R(J))**(1/6.)
         WRITE(31,94) ' assign (resid ',N1A(J),' and name ',R1A(J),
     &                ') (resid ',N2A(J),' and name ',R2A(J),') ',
     &                RNOE,RNOE,' 0'
      END DO
      CLOSE(31)
      
91    FORMAT (A4)
92    FORMAT (A3,1X,I4,1X,A4,1X,I4,1X,A4,1X,F6.2)
94    FORMAT (A15,I4,A10,A4,A9,I4,A10,A4,A2,F8.3,F8.3,A2)

      RETURN
666   PRINT*,'The file:',KT1INPUT,' does not exist...'
      RETURN
667   PRINT*,'The file:',KPDB,' does not exist...'
      RETURN

      END

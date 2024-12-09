C
C File: cns.inc
C =============
C
C This include file contains parameters which define array dimensions in
C CNS. These parameters are used by other include files to define sizes
C of common-block arrays
C
C Authors: Axel T. Brunger, Paul Adams, Ralf Grosse-Kunstleve
C
C=================================================================
C
C maximum number of atoms
C NOTE: 5-char Hybrid-36 limit is 87440031
C
      INTEGER    MAX_ATOM
      PARAMETER (MAX_ATOM=1000000)
C
C=================================================================
C
C word size
C
      INTEGER    WORD_SIZE
      PARAMETER (WORD_SIZE=264)
C
C=================================================================
C
C string size
C
      INTEGER    STRING_SIZE
      PARAMETER (STRING_SIZE=264)
C
C=================================================================
C
C number of buffer lines
C
      INTEGER    NUM_BUFFER
      PARAMETER (NUM_BUFFER=65536)
C
C=================================================================
C
C maximum number of pairs of interacting groups
C
      INTEGER    MAX_PIGS
      PARAMETER (MAX_PIGS=5000)
C
C=================================================================
C
C length of xray-object names
C
      INTEGER    XRNAME_SIZE
      PARAMETER (XRNAME_SIZE=20)
C
C=================================================================
C

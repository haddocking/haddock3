C ==============
C
C collapse.fcm
C
C holds a few things for the collapse energy term
C 
C JJK 6/14/96
C ==============
      DOUBLE PRECISION KCOLLAPSE, RTARGET
      INTEGER CLPSSELPTR, NSELECTED
C
C common blocks
C
      COMMON /ICLP/ CLPSSELPTR, NSELECTED
      COMMON /RCLP/ KCOLLAPSE, RTARGET
      SAVE /ICLP/
      SAVE /RCLP/


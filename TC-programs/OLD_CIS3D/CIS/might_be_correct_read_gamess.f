C*********************************************************************
      SUBROUTINE dafrd(V, LEN, RECN, FNAME)
C
C     read LEN double precision data to V from FNAME with 
C     stringlength SLEN

      implicit none


C     INPUT/OUTPUT DATA
      integer*8   SLEN
      integer*8   LEN 
      integer*8   RECN
      double precision V(LEN)

      character*7  FNAME
C      dimension FNAME(SLEN)


C     INTERNAL DATA
      integer*8 RL                   
      integer*8 RC
      
      integer*8 IODA(950)

      double precision  DUMREC(4090)
      integer*8 DCC
      integer*8 CRNR
      integer*8 I,J

      RC = 4090 
      RL = RC*8

C      write (*,*) FNAME

C     open DAF file 
      open(30,FILE=FNAME,ACCESS='DIRECT',FORM='UNFORMATTED',
     &     STATUS='UNKNOWN',RECL=RL)
      
      
C     READ "table of contents"
      READ(30,REC=1) IODA

C     LEN=6
C      DO J=1,160
C         IF(IODA(J) .NE. -1) THEN 
C            write (*,*) J, IODA(J)
C     START READING
            
C            DCC=RC+1
C            CRNR = IODA(J+1)
C            
C            DO I=1,LEN
C               IF (DCC .GT. RC) THEN
C                  READ(30,REC=CRNR) DUMREC
C                  DCC=1
C                  CRNR=CRNR+1
C               ENDIF
C               V(I) = DUMREC(DCC)
C               DCC = DCC + 1
C               write (*,*) V(I)
C            ENDDOC

C         ENDIF
C      ENDDO

C     START READING
      
      DCC=RC+1
      CRNR = IODA(RECN+1)
      
      DO I=1,LEN
         IF (DCC .GT. RC) THEN
            READ(30,REC=CRNR) DUMREC
            DCC=1
            CRNR=CRNR+1
         ENDIF
         V(I) = DUMREC(DCC)
         DCC = DCC + 1
C        write (*,*) V(I)
      ENDDO
      
      CLOSE(30)

      END
      

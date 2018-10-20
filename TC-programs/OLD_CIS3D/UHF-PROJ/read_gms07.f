C*********************************************************************
      SUBROUTINE dafrd(V, LEN, RECN, FNAME)
C
C     read LEN double precision data to V from FNAME with 
C     stringlength SLEN

      implicit none


C     INPUT/OUTPUT DATA
      integer   SLEN
      integer   LEN 
      integer   RECN
      double precision V(LEN)

      character*512  FNAME
C      dimension FNAME(SLEN)

C     INTERNAL DATA
      integer RL                   
      integer RC
      parameter ( RC = 4090, RL = RC*8)
      
      integer IODA(950)

      double precision  DUMREC(4090)
      integer DCC
      integer CRNR
      integer I

C     open DAF file 
      open(30,FILE=FNAME,ACCESS='DIRECT',FORM='UNFORMATTED',
     &     RECL=RL)
      
      
C     READ "table of contents"
      READ(30,REC=1) IODA

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
      ENDDO
      
      CLOSE(30)

      END
      

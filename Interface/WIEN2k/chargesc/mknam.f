      SUBROUTINE MKNAM(FNAME,OLDNAM,ILOOP)
      
! create a filename with running index ILOOP from 
! a given 'parental' file name OLDNAM

      CHARACTER*180 FNAME,OLDNAM
      CHARACTER*4  ALOOP

      WRITE(ALOOP,'(I4)')ILOOP

      FNAME=OLDNAM

      DO I=LEN(ALOOP),1,-1
         IF(ALOOP(I:I).NE.' ') IFROM=I
      ENDDO

      ITO=LEN(ALOOP)

      DO 7777 I=LEN(FNAME),1,-1
         IF(FNAME(I:I).NE.' ') THEN
            FNAME(I+1:LEN(FNAME))='_' // ALOOP(IFROM:ITO)
            GOTO 7778
         ENDIF
 7777 CONTINUE
 7778 CONTINUE

      RETURN
      END



SUBROUTINE READW
  USE param
  USE struct
  USE com
	implicit real*8 (a-h,o-z)
	 
         READ(26)EF
         k=0
 10      READ(26,end=4000,err=4000)nek
         k=k+1
         READ(26)(weigh(k,nn),nn=1,nek)
         GOTO 10
 4000    CONTINUE
 78   FORMAT(7X,I4,3x,I4)
 79   FORMAT(2(F16.12))
 80   FORMAT(/,14x,F12.9,I8)
 81   FORMAT(A30)

         end


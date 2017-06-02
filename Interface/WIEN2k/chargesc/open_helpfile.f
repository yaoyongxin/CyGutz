SUBROUTINE open_helpfile(fnamehelp,jatom,fnamehelp2)
  IMPLICIT NONE
  INTEGER     , INTENT(IN)    :: jatom 
  CHARACTER*180, INTENT(IN)   :: fnamehelp
  CHARACTER*180, INTENT(OUT)  :: fnamehelp2
  INTEGER                     :: length,i,j,k
  i=jatom+30
  k=1000
  length=LEN_TRIM(fnamehelp)
  fnamehelp2=fnamehelp
  DO j=1,3
     k=k/10
     fnamehelp2(length+j:length+j)=CHAR(ICHAR("0")+i/k)
     i=i-i/k*k
  ENDDO
END SUBROUTINE open_helpfile

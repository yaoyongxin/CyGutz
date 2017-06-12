SUBROUTINE HARMON(N,X,Y,Z,LMAX2,F,DF,RI)                           
  use struk, ONLY: BR1
  IMPLICIT REAL*8 (A-H,O-Z)
  REAL*8, intent(out):: F(LMAX2+1,N),DF(LMAX2+1,N)
  INTEGER, intent(in):: N, LMAX2
  REAL*8, intent(in) :: X(N),Y(N),Z(N), RI
  ! locals
  REAL*8 :: A(3)
  LMX=LMAX2+1                                                        
  DO I=1,N    ! 1
     A(1)=X(I)*BR1(1,1)+Y(I)*BR1(1,2)+Z(I)*BR1(1,3)                    
     A(2)=X(I)*BR1(2,1)+Y(I)*BR1(2,2)+Z(I)*BR1(2,3)                    
     A(3)=X(I)*BR1(3,1)+Y(I)*BR1(3,2)+Z(I)*BR1(3,3)                    
     XM=SQRT(A(1)**2+A(2)**2+A(3)**2)                                  
     XA=RI*XM                                                          
     CALL SPHBES(LMAX2,XA,F(1,I))                                       
     CALL DVBES1(F(1,I),DF(1,I),XA,RI,LMX)                             
     DO J=1,LMX
        DF(J,I)=XM*DF(J,I)                                                
     ENDDO
  ENDDO
  RETURN                                                            
END SUBROUTINE HARMON


      SUBROUTINE DVBES1(FJ,DJ,SM,RI,NT)                                 
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
!-----X CALCULATE THE DERIVATIVES OF THE BESSEL FUNCTIONS.   X----X----X
!-----X   DJ=DFJ/DX WHERE X=SM*RI                                 X----X
!-----X                    D.D.KOELLING                      X----X----X
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DJ(*),FJ(*)                                             
      DATA ZUP/1.0D-5/
      X=SM                                                              
      IF(X.GT.ZUP) GOTO 20                                              
      DJ(1)=0.0D0     !ZERO                                                        
      DJ(2)=1.D0/3.D0 !THIRD                                                       
      DO 10 L=3,NT                                                      
   10 DJ(L)=0.0D0     !ZERO                                                        
      RETURN                                                            
   20 Q2=-1.D0/X      ! ONE/X                                                         
      Q3=Q2                                                             
      DJ(1)=-FJ(2)                                                      
      LM=1                                                              
      DO 30 L=2,NT                                                      
      Q3=Q3+Q2                                                          
      DJ(L)=FJ(LM)+Q3*FJ(L)                                             
      LM=LM+1                                                           
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               

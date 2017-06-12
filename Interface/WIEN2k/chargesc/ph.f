      REAL*8 FUNCTION PH(N)                                             
      IMPLICIT REAL*8 (A-H,O-Z)
      PH=1.0D0                                                          
      K=MOD(IABS(N),4)                                                  
      IF(K.EQ.0) GOTO 5                                                 
      PH=-1.0D0                                                         
    5 RETURN                                                            
      END                                                               

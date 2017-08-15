      REAL*8 FUNCTION T3J0(J1,J2,J3)                                    
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FACT/FCT(100)                                              
      J=J1+J2+J3                                                        
      F1=SQRT(FCT(J-2*J1+1)*FCT(J-2*J2+1)*FCT(J-2*J3+1)/FCT(J+3))       
      F2=FCT(J/2+1)/FCT(J/2-J1+1)/FCT(J/2-J2+1)/FCT(J/2-J3+1)           
      T3J0=PH(J/2)*F1*F2                                                
      RETURN                                                            
      END                                                               

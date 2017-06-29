      SUBROUTINE RINT13(A,B,X,Y,S,JATOM)                            
!                                                                       
!     PERFORM RADIAL INTEGRALS REQUIRED BY BHDK13                       
!                            D.D.KOELLING                               
        USE param
        USE struct
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
!      LOGICAL REL                                                       
      DIMENSION A(NRAD),B(NRAD),X(NRAD),Y(NRAD)                         
      D=EXP(DX(JATOM))                                                  
      CIN=1.d0/137.0359895d0**2                           
      IF(.NOT.REL) CIN=1E-22                                            
!                                                                       
      J=3-MOD(JRJ(JATOM),2)                                             
      J1=J-1                                                            
      R=R0(JATOM)*(D**(J-1))                                          
      R1=R/D                                                            
      Z4=0                                                              
      Z2=0                                                              
   10 Z4=Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))                                 
      R=R*D                                                             
      J=J+1                                                             
      IF(J.GE.JRJ(JATOM)) GOTO 20                                       
      Z2=Z2+R*(A(J)*X(J)+CIN*B(J)*Y(J))                                 
      R=R*D                                                             
      J=J+1                                                             
      GOTO 10                                                           
   20 P1=R0(JATOM)*(A(1)*X(1)+CIN*B(1)*Y(1))                          
      P2=R1*(A(J1)*X(J1)+CIN*B(J1)*Y(J1))                               
      S=2*Z2+4*Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))+P2                        
      S=(DX(JATOM)*S+P1)/3.0D0                                          
      IF(J1.GT.1) S=S+0.5D0*DX(JATOM)*(P1+P2)                           
      RETURN                                                            
      END                                                               

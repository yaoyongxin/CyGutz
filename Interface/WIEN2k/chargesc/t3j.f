      REAL*8 FUNCTION T3J (J1,J2,J3,M1,M2,M3)                           
!                                                                       
!     ANGULAR MOMENTA ENTER WITH 2X ACTUAL VALUES                       
!     ARRAY FCT MUST BE DECLARED AND FILLED GLOBALLY WITH               
!     J-FACTORIAL IN THE 2J+1 POSITION                                  
!     LARGEST  FACTORIAL NEEDED IS JA+JB+JC+1                           
!     FUNCTION NOTRI(J,K,L) RETURNS 1 IF K,L,J FORM A TRIANGLE,-1 IF NOT
!     FUNCTION PH IS NEEDED TO COMPUTE PHASES                           
!     NOTRI AND M-TEST IN MAIN PROGRAM                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FACT/FCT(100)                                              
!                                                                       
      T3J=0.                                                            
!     IF((M1+M2+M3).NE.0) GOTO 101                                      
!     IF(NOTRI(J1,J2,J3).LT.0) GOTO 102                                 
      K0=J1+J2-J3                                                       
      K1=J1-M1                                                          
      K2=J2+M2                                                          
      L1=J2-J3-M1                                                       
      L2=J1-J3+M2                                                       
      KMAX=MIN0(K1,K2,K0)+1                                             
      KMIN=MAX0(0,L1,L2)+1                                              
      C=SQRT(FCT(K0+1)*FCT(J1-J2+J3+1)*FCT(J2-J1+J3+1)/FCT(J1+J2+J3+3)   &
             *FCT(J1+M1+1)*FCT(K1+1)*FCT(J2-M2+1)*FCT(K2+1)*FCT(J3+M3+1) &
             *FCT(J3-M3+1))                                             
      DO 10 I=KMIN,KMAX,2                                               
      K=I-1                                                             
      T3J=T3J+C*PH(K)/FCT(I)/FCT(K0-K+1)/FCT(K1-K+1)/FCT(K2-K+1)         &
               /FCT(I-L1)/FCT(I-L2)                                     
   10 CONTINUE                                                          
      T3J=PH(J1-J2-M3)*T3J                                              
  100 RETURN                                                            
! 101 WRITE(6,1)                                                        
!   1 FORMAT(1H ,"M1+M2+M3.NE0")                                        
!     STOP                                                              
! 102 WRITE(6,2) J1,J2,J3                                               
!   2 FORMAT(1H ," TRIANGLE RULE NOT SATISFIED FOR (",3I4," )/2")       
!     STOP                                                              
      END                                                               

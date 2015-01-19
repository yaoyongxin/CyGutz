SUBROUTINE RINT13(REL,A,B,X,Y,S,JATOM)                            
  ! Input: A, B, X, Y, jatom
  !
  ! Output: S
  !                                                                       
  !     PERFORM RADIAL INTEGRALS REQUIRED BY BHDK13                       
  !                            D.D.KOELLING                               
  use param, ONLY: NRAD
  use struk, ONLY: DX, JRI, R0
  IMPLICIT NONE
  !
  LOGICAL, intent(in) :: REL   ! relativistic or not
  REAL*8, intent(in)  :: A(NRAD),B(NRAD),X(NRAD),Y(NRAD)
  REAL*8, intent(out) :: S
  INTEGER, intent(in) :: jatom
  !locals
  REAL*8 :: D, CIN, R, R1, P1, P2, Z, Z2, Z4
  INTEGER :: J, J1
  !
  D=EXP(DX(JATOM))                                                  
  CIN=1.d0/137.0359895d0**2                           
  IF(.NOT.REL) CIN=1E-22                                            
  !                                                                       
  J=3-MOD(JRI(JATOM),2)                                             
  J1=J-1                                                            
  R=R0(JATOM)*(D**(J-1))                                          
  R1=R/D                                                            
  Z4=0                                                              
  Z2=0                                                              

  DO
     Z4=Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J)) ! Z4 ~ R*A*X (r)
     R=R*D                                                             
     J=J+1                                                             
     IF(J.GE.JRI(JATOM)) EXIT
     Z2=Z2+R*(A(J)*X(J)+CIN*B(J)*Y(J)) ! Z2 ~ R*A*X (r+dr)                               
     R=R*D                                                             
     J=J+1                                                             
  ENDDO
  
  P1=R0(JATOM)*(A(1)*X(1)+CIN*B(1)*Y(1))
  P2=R1*(A(J1)*X(J1)+CIN*B(J1)*Y(J1))
  S=2*Z2+4*Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))+P2
  S=(DX(JATOM)*S+P1)/3.0D0
  IF(J1.GT.1) S=S+0.5D0*DX(JATOM)*(P1+P2)
  RETURN                                                            
  
END SUBROUTINE RINT13

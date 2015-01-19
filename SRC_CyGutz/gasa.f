!*********************************************************************
!  this modul contains code to calculate Clebsch-Gordan coefficients
!  and code to evaluate the integrals of three spherical (complex) harmonics
!*********************************************************************
!
      MODULE gasa
      USE gprec
!
      CONTAINS
!****************************************************************************
! function to evaluate (-1)^I
!****************************************************************************
      FUNCTION FS(I)
      USE gprec
      IMPLICIT NONE
      INTEGER I,FS
!
      FS=1-2*MOD(I+40,2)
      RETURN
!
      END FUNCTION FS
!
!*********************************************************************
      FUNCTION CLEBGO(FAC,J1,J2,J3,M1,M2,M3)
!
      IMPLICIT REAL(gq) (A-H,O-Z)
      REAL(gq) FAC(40)
      REAL(gq) CLEBGO
!
      IF(M3/=M1+M2) GO TO 2
      K1=J1+J2-J3+1
      K2=J3+J1-J2+1
      K3=J3+J2-J1+1
      K4=J1+J2+J3+2
      T= (2*J3+1)*FAC(K1)*FAC(K2)*FAC(K3)/FAC(K4)
      K1=J1+M1+1
      K2=J1-M1+1
      K3=J2+M2+1
      K4=J2-M2+1
      K5=J3+M3+1
      K6=J3-M3+1
      T=SQRT(T*FAC(K1)*FAC(K2)*FAC(K3)*FAC(K4)*FAC(K5)*FAC(K6))
      N1=MAX0(J2-J3-M1,J1-J3+M2,0)+1
      N2=MIN0(J1+J2-J3,J1-M1,J2+M2)+1
      IF(N1>N2) GO TO 2
      T1=0.0_gq
      DO M=N1,N2
         N=M-1
         K1=J1+J2-J3-N+1
         K2=J1-M1-N+1
         K3=J2+M2-N+1
         K4=J3-J2+M1+N+1
         K5=J3-J1-M2+N+1
         T1=T1+ (1+4*(N/2)-2*N)/(FAC(M)*FAC(K1)*FAC(K2)*FAC(K3) &
              &  *FAC(K4)*FAC(K5))
      ENDDO
      CLEBGO=T*T1
      RETURN
! coefficient is (0._gq,0._gq), drop back
 2    CONTINUE
      CLEBGO=0.0_gq
      RETURN
!
      END FUNCTION
!
!*********************************************************************
      FUNCTION CLEBG0(FAC,L1,L2,L3)
      USE gprec
      IMPLICIT REAL(gq) (A-H,O-Z)
      INTEGER X,P
      REAL(gq) FAC(40)
      REAL(gq) CLEBG0
!
      LT=L1+L2+L3
      P=LT/2
      IF(2*P/=LT) GO TO 1
      CLEBG0= SQRT( REAL(2*L3+1,KIND=gq)/(LT+1))
      CLEBG0=CLEBG0*FAC(P+1)/SQRT(FAC(2*P+1))
      X=P-L1
      CLEBG0=CLEBG0*SQRT(FAC(2*X+1))/FAC(X+1)
      X=P-L2
      CLEBG0=CLEBG0*SQRT(FAC(2*X+1))/FAC(X+1)
      X=P-L3
      CLEBG0=CLEBG0*SQRT(FAC(2*X+1))/FAC(X+1)
      IF(X>2*(X/2)) CLEBG0=-CLEBG0
      RETURN
! coefficient is (0._gq,0._gq), drop back
 1    CONTINUE
      CLEBG0=0.0_gq
      RETURN
      END FUNCTION
!
!************************* YLM3ST_COMPL ******************************
! calculate the integral of the product of three complex
! spherical harmonics
! i.e    Y_lm Y_l'm' Y_LM
! LMAX     max value for l and lp (maximum L is given by triagular rule
!             | l- lp | < L < | l + lp |
! YLM3     results (on exit)
!*********************************************************************
      SUBROUTINE YLM3ST_COMPL(L1,L2,GAUNT)
      USE gprec
      IMPLICIT REAL(gq) (A-H,O-Z)
      REAL(gq) GAUNT(-L1:L1,-L2:L2,ABS(L1-L2):L1+L2,-(L1+L2):L1+L2)
      REAL(gq) FAC(40)
      REAL(gq), PARAMETER :: SRPI =1.772453850905516027_gq
!
! set up table for factorials
      IMAX=39
      FAC(1)=1._gq
      DO I=1,IMAX
         FAC(I+1)= I*FAC(I)
      ENDDO
!
      GAUNT=0
      K2=(2*L1+1)*(2*L2+1)
      DO M1=-L1,L1; DO M2=-L2,L2
         M3=-(M1+M2)
!---------------------------------------------------------------------
! loop over L given by triangular rule
!---------------------------------------------------------------------
         Q1= SQRT( REAL(K2,KIND=gq)/4 )/SRPI*FS(M3)
         DO L3=L1-L2,L1+L2,2
            IF(ABS(M3)>L3) CYCLE
            T =CLEBGO(FAC(1),L1,L2,L3, M1, M2, -M3)
            T0=CLEBG0(FAC(1),L1,L2,L3)
            GAUNT(M1,M2,L3,M3)=Q1*T*T0/(SQRT( REAL(2*L3+1, KIND=gq)))
         ENDDO
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE
!
      END MODULE gasa

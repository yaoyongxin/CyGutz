!*************************************************!
! TEST DPFA_PA
!*************************************************!
      SUBROUTINE TEST_DPFA_PA()
      USE GUTIL; USE GPREC
      IMPLICIT NONE
      INTEGER,PARAMETER::N=2
      REAL(gq) A(N,N),PA(N,N),H(N,N)
!
      DATA A/0.43,0.1,0.1,0.2/
      H=0
      H(1,2)=1.D0/SQRT(2.D0); H(2,1)=H(1,2)
      CALL DPFA_PA(A,PA,H,N,DSIMIX,DPSIMIX)
      WRITE(*,*)'A='
      WRITE(*,*)A
      WRITE(*,*)'PA='
      WRITE(*,*)PA
      RETURN
!
      END SUBROUTINE TEST_DPFA_PA
!*************************************************!
! TEST ZPFA_PA
!*************************************************!
      SUBROUTINE TEST_ZPFA_PA()
      USE GUTIL; USE GPREC
      IMPLICIT NONE
      INTEGER,PARAMETER::N=2
      COMPLEX(gq) A(N,N),PA(N,N),H(N,N)
!
      DATA A/(0.43,0.0),(0.1,-0.1),(0.1,0.1),(0.2,0.0)/
      H=0
      H(1,2)=(0.D0,1.D0)/SQRT(2.D0); H(2,1)=CONJG(H(1,2))
      CALL ZPFA_PA(A,PA,H,N,DSIMIX,DPSIMIX)
      WRITE(*,*)'A='
      WRITE(*,*)A
      WRITE(*,*)'PA='
      WRITE(*,*)PA
      RETURN
!
      END SUBROUTINE TEST_ZPFA_PA

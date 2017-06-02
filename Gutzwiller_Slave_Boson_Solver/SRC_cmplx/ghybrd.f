!****************************************************************************
      SUBROUTINE GHYBRD(F,N,X,FVEC,RTOL,EPSFCN,IO)
      USE gprec
      IMPLICIT NONE
      INTEGER N,IO
      REAL(gq) X(N),FVEC(N),RTOL,EPSFCN
      EXTERNAL :: F
! LOCAL
      INTEGER INFO
!
      CALL HYBRD1(F,N,X,FVEC,RTOL,EPSFCN,INFO)
      SELECT CASE(INFO)
      CASE(0)
        WRITE(IO,'(" GHYBRD: ERROR! Improper input parameters.")')
        WRITE(0,'(" GHYBRD: ERROR! Improper input parameters.")')
      CASE(1)
        WRITE(IO,'(" GHYBRD: Success.")')
      CASE(2)
        WRITE(IO,'(" GHYBRD: WARNING! Number of calls to FCN has reached or exceeded 200*(N+1).")')
        WRITE(0,'(" GHYBRD: WARNING! Number of calls to FCN has reached or exceeded 200*(N+1).")')
      CASE(3)
        WRITE(IO,'(" GHYBRD: WARNING! TOL is too small. No further improvement in the approximate solution X is possible.")')
        WRITE(0,'(" GHYBRD: WARNING! TOL is too small. No further improvement in the approximate solution X is possible.")')
      CASE(4)
        WRITE(IO,'(" GHYBRD: WARNING! Iteration is not making good progress.")')
        WRITE(0,'(" GHYBRD: WARNING! Iteration is not making good progress.")')
      END SELECT
      RETURN
!
      END SUBROUTINE GHYBRD

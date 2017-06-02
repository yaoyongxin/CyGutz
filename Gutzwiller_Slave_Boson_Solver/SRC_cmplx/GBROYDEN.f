!
      MODULE GBROYDEN
      USE gprec
      IMPLICIT NONE
!
      TYPE BROYDEN_DATA
        INTEGER N,IT_RESET
        LOGICAL :: LFIRST=.TRUE.
        REAL(gq) :: ALPHA_MIX
        REAL(gq),ALLOCATABLE :: X0(:), F0(:), AMIX(:,:)
      END TYPE BROYDEN_DATA
!
      TYPE(BROYDEN_DATA) :: BR_DATA
!
      CONTAINS
!***************************************************************
      SUBROUTINE INIT_BROYDEN_DATA(N,ALPHA_MIX,IT_RESET)
      INTEGER,INTENT(IN) :: N,IT_RESET
      REAL(gq),INTENT(IN) :: ALPHA_MIX
! LOCAL
      INTEGER I
!
      IF(.NOT.BR_DATA%LFIRST)RETURN
      BR_DATA%N=N
      ALLOCATE(BR_DATA%X0(N), BR_DATA%F0(N), BR_DATA%AMIX(N,N))
      BR_DATA%ALPHA_MIX = ALPHA_MIX
      BR_DATA%IT_RESET = IT_RESET
      RETURN
!      
      END SUBROUTINE INIT_BROYDEN_DATA
!
!***************************************************************
      SUBROUTINE INV_BROYDEN_PREDICT(X,F)
      REAL(gq) X(BR_DATA%N),F(BR_DATA%N)
! LOCAL
      INTEGER I,J
      REAL(gq) RES,VNORM,VFAC
      REAL(gq) DX(BR_DATA%N),DF(BR_DATA%N),AV(BR_DATA%N),VA(BR_DATA%N)
      IF(BR_DATA%LFIRST)THEN 
        BR_DATA%LFIRST=.FALSE.
        BR_DATA%AMIX=0
        DO I=1,BR_DATA%N; BR_DATA%AMIX(I,I)=BR_DATA%ALPHA_MIX; ENDDO
      ELSE
        DX=X-BR_DATA%X0; DF=F-BR_DATA%F0
        AV=MATMUL(BR_DATA%AMIX,DF)
        VA=DF; VNORM=DOT_PRODUCT(DF,DF)
        VFAC = 1/MAX(VNORM,1.E-20_gq)
        AV=AV-DX
        DO I=1,BR_DATA%N; DO J=1,BR_DATA%N
          BR_DATA%AMIX(I,J) = BR_DATA%AMIX(I,J) - VFAC*AV(I)*VA(J)
        ENDDO; ENDDO
      ENDIF
! COMPUTE CORRECTION TO VIN FROM INVERSE JACOBIAN.
      AV=MATMUL(BR_DATA%AMIX,F)
! NO TOO BIG STEP
      RES=MAXVAL(ABS(AV))
      IF(RES.GT.1._gq)THEN
        RES=1._gq/RES
        AV=AV*RES
        BR_DATA%AMIX=BR_DATA%AMIX*RES
      ENDIF
! UPDATE X
      BR_DATA%X0=X; BR_DATA%F0=F
      X=X-AV
      RETURN
!
      END SUBROUTINE INV_BROYDEN_PREDICT
!
!*************************************************************
      SUBROUTINE BROYDEN_MIX(FCN,N,X,F,ALPHA,IRESET,MODE)
      INTEGER,INTENT(IN) :: N,MODE,IRESET
      REAL(gq),INTENT(INOUT) :: X(N),F(N)
      REAL(gq),INTENT(IN) ::  ALPHA
      EXTERNAL :: FCN
! LOCAL
      INTEGER IT,IFLAG,I_INNER
      REAL(gq) ALPHA_MIX
      INTEGER,PARAMETER :: I_INNER_MAX = 7
      REAL(gq),PARAMETER :: RTOL = 1.E-6_gq
      REAL(gq) F_MAX, F_MAX_BEST, X_BEST(N)
!
      F_MAX_BEST = 1.E4_gq; I_INNER=0
      ALPHA_MIX = ALPHA
      IF(MODE==-1)THEN
        CALL INIT_BROYDEN_DATA(N,ALPHA,IRESET)
      ENDIF
      DO
        IT=IT+1
        CALL FCN(N,X,F,IFLAG)
        F_MAX = MAXVAL(ABS(F))
        IF(IFLAG==-1)THEN
          IF(F_MAX<F_MAX_BEST)THEN
            RETURN
          ELSE
            X = X_BEST
            CALL FCN(N,X,F,IFLAG)
            RETURN
          ENDIF
        ENDIF
        IF(F_MAX<F_MAX_BEST)THEN
          I_INNER = 0
          F_MAX_BEST = F_MAX
          X_BEST = X
        ELSEIF(MODE==2)THEN ! SIMPLE MIX
          I_INNER = I_INNER + 1
          IF(I_INNER>I_INNER_MAX)THEN
            ALPHA_MIX = ALPHA_MIX/2
            X = X_BEST
            I_INNER = 0
          ENDIF
        ENDIF
        IF(MODE==-1)THEN ! INVERSE BROYDEN
          ! Reset broyden
          IF(MOD(IT,BR_DATA%IT_RESET)==0)THEN
            BR_DATA%LFIRST=.TRUE.
          ENDIF
          CALL INV_BROYDEN_PREDICT(X,F)
        ELSEIF(MODE==2)THEN
          X = X+ALPHA*F
        ELSE
          STOP " ERROR IN BROYDEN MIX: UNSUPPORTED MODE!"
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE BROYDEN_MIX
!
      END MODULE GBROYDEN

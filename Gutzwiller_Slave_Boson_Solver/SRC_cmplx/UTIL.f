      MODULE GUTIL
      USE gprec; USE GCONSTANT
      IMPLICIT NONE
      private
      public :: ORBITAL_SPIN_TRANS, A_TRANS, GET_YLM_CR, GET_YLM_CR_ALL,  &
              & OUT_TIME_USE, INV_APLX, WRT_A, HERMEV, INV, CHK_UNITARY, &
              & GET_HM_EXPAND, UHAU, ANMXBMM, ANNXB, TRACE_A, ATOFA, PFA_PA, &
              & H_REGULARIZE, CALC_LOEWNER, DSIMIX, DPSIMIX, &
              & FILE_NAME, int_to_str, FERMI_FUN, GAUSS_FUN, SET_RANGE, &
              & SET_LINEAR_SIMP, MAX_INTERVAL, LKEY_EXIST
!
      INTERFACE ORBITAL_SPIN_TRANS
        MODULE PROCEDURE DSPIN_BLK1_TRANS, ZSPIN_BLK2U_TRANS, DSPIN_BLK2_TRANS, &
                       & ZSPIN_BLK2_TRANS
      END INTERFACE ORBITAL_SPIN_TRANS
!
      INTERFACE A_TRANS
        MODULE PROCEDURE A4_TRANS_RC, A4_TRANS_C, ZA2_TRANS, DA2_TRANS
      END INTERFACE A_TRANS
!
      INTERFACE INV_APLX
        MODULE PROCEDURE ZINV_APLX
      END INTERFACE INV_APLX
!
      INTERFACE WRT_A
        MODULE PROCEDURE DWRT_ANN, ZWRT_ANN
      END INTERFACE WRT_A
!
      INTERFACE HERMEV
        MODULE PROCEDURE ZHEEV_, DSYEV_, ZHEEVX_, ZHEEVX_I, DSYEVX_, DSYEVX_I, &
                       & DDIAG_JZJP, ZDIAG_JZJP
      END INTERFACE HERMEV
!
      INTERFACE HERMEVD
        MODULE PROCEDURE ZHEEVD_, DSYEVD_
      END INTERFACE HERMEVD
!
      INTERFACE INV
        MODULE PROCEDURE DINV_, ZINV_
      END INTERFACE INV
!
      INTERFACE GET_HM_EXPAND
        MODULE PROCEDURE GET_DHM_EXPAND, GET_ZHM_EXPAND
      END INTERFACE GET_HM_EXPAND
!
      INTERFACE UHAU
        MODULE PROCEDURE DUHAU, ZUHAU
      END INTERFACE UHAU
!
      INTERFACE ANMXBMM
        MODULE PROCEDURE DANMXBMM, ZANMXBMM
      END INTERFACE ANMXBMM
!
      INTERFACE ANNXB
        MODULE PROCEDURE DANNXBNM, ZANNXBNM, ZANNXBNN, DANNXBNN
      END INTERFACE ANNXB
!
      INTERFACE TRACE_A
        MODULE PROCEDURE TRACE_DANN, TRACE_ZANN, TRACE_ZANND
      END INTERFACE TRACE_A
!
      INTERFACE ATOFA
        MODULE PROCEDURE ZATOFA, DATOFA, ZATOFA1, DATOFA1
      END INTERFACE ATOFA
!
      INTERFACE PFA_PA
        MODULE PROCEDURE DPFA_PA, ZPFA_PA
      END INTERFACE PFA_PA
!
      INTERFACE H_REGULARIZE
        MODULE PROCEDURE DH_REGULARIZE, ZH_REGULARIZE
      END INTERFACE H_REGULARIZE
!
      CONTAINS
!****************************************************************************
! ORBITAL FAST <-> SPIN FAST
!****************************************************************************
      SUBROUTINE DSPIN_BLK1_TRANS(A,N,LSFAST,ISO)
      INTEGER N,ISO
      REAL(gq) A(N)
      LOGICAL LSFAST
! LOCAL
      INTEGER I,IADD
      REAL(gq) U(N,N)
!
      IF(ISO.EQ.2)RETURN ! SOC, NOT NECESSARY
      U=0; IADD=0
      DO I=1,N,2; IADD=IADD+1; U(IADD,I)=1._gq; ENDDO
      DO I=2,N,2; IADD=IADD+1; U(IADD,I)=1._gq; ENDDO
      IF(LSFAST)THEN; A=MATMUL(TRANSPOSE(U),A)
      ELSE          ; A=MATMUL(U,A)
      ENDIF
      RETURN
!
      END SUBROUTINE DSPIN_BLK1_TRANS
!
!****************************************************************************
! ORBITAL FAST <-> SPIN FAST
!****************************************************************************
      SUBROUTINE ZSPIN_BLK2U_TRANS(A,N,M,LSFAST,ISO,LURIGHT)
      INTEGER N,M,ISO
      COMPLEX(gq) A(N,M)
      LOGICAL LSFAST
      LOGICAL,OPTIONAL::LURIGHT
! LOCAL
      INTEGER I,IADD,NM
      LOGICAL LUR
      COMPLEX(gq),ALLOCATABLE::U(:,:)
!
      IF(ISO.EQ.2) RETURN ! SOC, NOT NECESSARY
      LUR=.FALSE.
      IF(PRESENT(LURIGHT)) LUR=LURIGHT
      IF(LUR)THEN; NM=M; ELSE; NM=N; ENDIF
      ALLOCATE(U(NM,NM))
      U=0; IADD=0
      DO I=1,NM,2; IADD=IADD+1; U(IADD,I)=1._gq; ENDDO
      DO I=2,NM,2; IADD=IADD+1; U(IADD,I)=1._gq; ENDDO
      IF(LSFAST)THEN
        IF(.NOT.LUR)THEN
          CALL ZANNXBNM('C',U,A,N,M)
        ELSE
          CALL ZANMXBMM('N',A,U,N,M)
        ENDIF
      ELSE
        IF(.NOT.LUR)THEN
          CALL ZANNXBNM('N',U,A,N,M)
        ELSE
          CALL ZANMXBMM('C',A,U,N,M)
        ENDIF
      ENDIF
      DEALLOCATE(U)
      RETURN
!
      END SUBROUTINE ZSPIN_BLK2U_TRANS
!
!****************************************************************************
! ORBITAL FAST <-> SPIN FAST
!****************************************************************************
      SUBROUTINE DSPIN_BLK2_TRANS(A,N,LSFAST,ISO)
      INTEGER N,ISO
      REAL(gq) A(N,N)
      LOGICAL LSFAST
! LOCAL
      INTEGER I,IADD
      REAL(gq) U(N,N)
!
      IF(ISO.EQ.2)RETURN ! SOC, NOT NECESSARY
      U=0; IADD=0
      DO I=1,N,2; IADD=IADD+1; U(IADD,I)=1._gq; ENDDO
      DO I=2,N,2; IADD=IADD+1; U(IADD,I)=1._gq; ENDDO
      IF(LSFAST)THEN
        CALL DUHAU(A,U,N,N,TRUL='C',TRUR='N')
      ELSE
        CALL DUHAU(A,U,N,N,TRUL='N',TRUR='C')
      ENDIF
      RETURN
!
      END SUBROUTINE DSPIN_BLK2_TRANS
!
!****************************************************************************
! ORBITAL FAST <-> SPIN FAST
!****************************************************************************
      SUBROUTINE ZSPIN_BLK2_TRANS(A,N,LSFAST,ISO)
      INTEGER N,ISO
      COMPLEX(gq) A(N,N)
      LOGICAL LSFAST
! LOCAL
      INTEGER I,IADD
      COMPLEX(gq) U(N,N)
!
      IF(ISO.EQ.2)RETURN ! SOC, NOT NECESSARY
      U=0; IADD=0
      DO I=1,N,2; IADD=IADD+1; U(IADD,I)=1._gq; ENDDO
      DO I=2,N,2; IADD=IADD+1; U(IADD,I)=1._gq; ENDDO
      IF(LSFAST)THEN
        CALL ZUHAU(A,U,N,N,TRUL='C',TRUR='N')
      ELSE
        CALL ZUHAU(A,U,N,N,TRUL='N',TRUR='C')
      ENDIF
      RETURN
!
      END SUBROUTINE ZSPIN_BLK2_TRANS
!
!****************************************************************************
      SUBROUTINE A4_TRANS_RC(A,N,U)
      INTEGER N
      REAL(gq) A(N,N,N,N)
      COMPLEX(gq) U(N,N)
! LOCAL
      REAL(gq) RES
      COMPLEX(gq),ALLOCATABLE::B(:,:,:,:)
!
      ALLOCATE(B(N,N,N,N)); B=DCMPLX(A)
      CALL A4_TRANS_C(B,N,U)
      RES=MAXVAL(ABS(AIMAG(B)))
      IF(RES.GT.1.E-6_gq)THEN
        WRITE(0,'(" ERROR IN A4_TRANS_RC: MAX AIMAG =",F12.6)')RES; STOP
      ENDIF
      A=REAL(B,gq)
      RETURN
!
      END SUBROUTINE A4_TRANS_RC
!
!****************************************************************************
      SUBROUTINE A4_TRANS_C(A,N,U)
      INTEGER N
      COMPLEX(gq) A(N,N,N,N),U(N,N)
! LOCAL
      INTEGER I1,I2
!
      DO I1=1,N; DO I2=1,N
        CALL ZA2_TRANS(A(:,I1,:,I2),N,U,1)
      ENDDO; ENDDO
      DO I1=1,N; DO I2=1,N
        CALL ZA2_TRANS(A(I1,:,I2,:),N,U,1)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE A4_TRANS_C
!
!****************************************************************************
      SUBROUTINE ZA2_TRANS(A,N,U,MODE)
      INTEGER N,MODE
      COMPLEX(gq) A(N,N),U(N,N)
!
      IF(MODE>=0)THEN
        A=MATMUL(TRANSPOSE(CONJG(U)),MATMUL(A,U))
      ELSE
        A=MATMUL(MATMUL(U,A),TRANSPOSE(CONJG(U)))
      ENDIF
      RETURN
!
      END SUBROUTINE ZA2_TRANS
!
!****************************************************************************
      SUBROUTINE DA2_TRANS(A,N,U,MODE)
      INTEGER N,MODE
      REAL(gq) A(N,N),U(N,N)
!
      IF(MODE>=0)THEN
        A=MATMUL(TRANSPOSE(U),MATMUL(A,U))
      ELSE
        A=MATMUL(MATMUL(U,A),TRANSPOSE(U))
      ENDIF
      RETURN
!
      END SUBROUTINE DA2_TRANS
!
!****************************************************************************
! Complex to real spherical Harmonics transformation <C_Ylm | R_Ylm>
!****************************************************************************
      SUBROUTINE GET_YLM_CR(UCR,L)
      INTEGER L
      COMPLEX(gq) UCR(2*L+1,2*L+1)
! LOCAL
      INTEGER M,M_
!
      UCR=0
      DO M=1,2*L+1
        M_=M-L-1
        IF(M_>0)THEN
          UCR( M_+L+1,M)=(-1)**M_/SQRT(2._gq)
          UCR(-M_+L+1,M)=1/SQRT(2._gq)
        ELSEIF(M_==0)THEN
          UCR(L+1,L+1)=1
        ELSE
          UCR( M_+L+1,M)= CMPLX(0,1/SQRT(2._gq),gq)
          UCR(-M_+L+1,M)=-CMPLX(0,(-1)**M_/SQRT(2._gq),gq)
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE GET_YLM_CR
!
!****************************************************************************
      SUBROUTINE GET_YLM_CR_ALL(YLM_CR,LMAX)
      INTEGER LMAX
      COMPLEX(gq) YLM_CR(2*LMAX+1,2*LMAX+1,0:LMAX)
! LOCAL
      INTEGER L,MDIM
!
      YLM_CR=0
      DO L=0,LMAX
        MDIM=2*L+1
        CALL GET_YLM_CR(YLM_CR(1:MDIM,1:MDIM,L),L)
      ENDDO
      RETURN
!
      END SUBROUTINE GET_YLM_CR_ALL
!
!****************************************************************************
      SUBROUTINE OUT_TIME_USE(TASK,CPU_TIME,WAL_TIME,IO)
      CHARACTER(*) TASK
      REAL CPU_TIME,WAL_TIME
      INTEGER NH,NM,IO
      REAL,PARAMETER::MINUTE=60._4,HOUR=MINUTE*MINUTE
!
      NH=INT(CPU_TIME/HOUR  ); CPU_TIME=MOD(CPU_TIME,HOUR  )
      NM=INT(CPU_TIME/MINUTE); CPU_TIME=MOD(CPU_TIME,MINUTE)
      WRITE(IO,'(A25," takes",I5," hs ",I2," ms and ",F7.2," ss (CPU).")')TASK,NH,NM,CPU_TIME
      NH=INT(WAL_TIME/HOUR  ); WAL_TIME=MOD(WAL_TIME,HOUR  )
      NM=INT(WAL_TIME/MINUTE); WAL_TIME=MOD(WAL_TIME,MINUTE)
      WRITE(IO,'(25X," takes",I5," hs ",I2," ms and ",F7.2," ss (WALL).")')NH,NM,WAL_TIME
      RETURN
!
      END SUBROUTINE OUT_TIME_USE
!
!*************************************************************************************
! Calculate (A+X*I)^-1 given A^-1; See
! http://www.jstor.org/stable/2690437?seq=3
! 'On the Inverse of the Sum of Matrices'
! Much Slower than the standard inversion way!!!
!*************************************************************************************
      SUBROUTINE ZINV_APLX(A,B,N,X)
      INTEGER N
      COMPLEX(gq) A(N,N),B(N,N)
      REAL(gq) X
! LOCAL
      INTEGER K
      COMPLEX(gq) GK
      COMPLEX(gq) BX(N),V(N)
!
      B=A
      DO K=1,N
        GK=-1._gq/(1._gq+B(K,K)*X)
        BX=B(:,K)*X; V=B(K,:)
        CALL ZGERU(N,N,GK,BX,1,V,1,B,N)
      ENDDO
      RETURN
!
      END SUBROUTINE ZINV_APLX
!
!*************************************************************************************
      SUBROUTINE DWRT_ANN(A,N,IO)
      INTEGER N,IO
      REAL(gq) A(N,N)
      CHARACTER(20) FMT
!
      WRITE(IO,'(" REAL PART")')
      WRITE(FMT,'("(",I0,"F9.4)")') N
      WRITE(IO,FMT)REAL(A,gq)
      RETURN
!
      END SUBROUTINE DWRT_ANN
!
!*************************************************************************************
      SUBROUTINE ZWRT_ANN(A,N,IO)
      INTEGER N,IO
      COMPLEX(gq) A(N,N)
      CHARACTER(20) FMT
!
      WRITE(IO,'(" REAL PART")')
      WRITE(FMT,'("(",I0,"F9.4)")') N
      WRITE(IO,FMT)REAL(A,gq)
      WRITE(IO,'(" IMAG PART")')
      WRITE(IO,FMT)AIMAG(A)
      RETURN
!
      END SUBROUTINE ZWRT_ANN
!
!******************************************************************
      SUBROUTINE ZHEEV_(JOBZ,UPLO,A,W,N)
      INTEGER N
      CHARACTER JOBZ,UPLO
      COMPLEX(gq) A(N,N)
      REAL(gq)    W(N)
! Local
      INTEGER INFO,LWORK
      REAL(8),ALLOCATABLE :: RWORK(:)
      COMPLEX(8),ALLOCATABLE :: WORK(:)
!
      LWORK=32*N
      ALLOCATE(RWORK(max(1,3*N-2)),WORK(LWORK))
      CALL ZHEEV(JOBZ,UPLO,N,A,N,W,WORK,LWORK,RWORK,INFO)
      IF(INFO.NE.0)THEN
        WRITE(0,'(" ERROR IN ZHEEV_: INFO=",I5," N=",I5)')INFO,N; STOP
      ENDIF
      DEALLOCATE(RWORK,WORK)
      RETURN
!
      END SUBROUTINE ZHEEV_
!
!******************************************************************
      SUBROUTINE ZHEEVD_(JOBZ,UPLO,A,W,N)
      INTEGER N
      CHARACTER JOBZ,UPLO
      COMPLEX(gq) A(N,N)
      REAL(gq)    W(N)
! Local
      INTEGER INFO,LWORK,LRWORK,LIWORK
      INTEGER,ALLOCATABLE :: IWORK(:)
      REAL(8),ALLOCATABLE :: RWORK(:)
      COMPLEX(8),ALLOCATABLE :: WORK(:)
!
      LWORK=8*N+N**2; LRWORK=1+25*N+2*N**2; LIWORK=3+25*N
      ALLOCATE(RWORK(LRWORK),WORK(LWORK),IWORK(LIWORK))
      CALL ZHEEVD(JOBZ,UPLO,N,A,N,W,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
      IF(INFO.NE.0)THEN
        WRITE(0,'(" ERROR IN ZHEEVD_: INFO=",I5," N=",I5)')INFO,N; STOP
      ENDIF
      DEALLOCATE(IWORK,RWORK,WORK)
      RETURN
!
      END SUBROUTINE ZHEEVD_
!
!******************************************************************
      SUBROUTINE ZHEEVX_I(JOBZ,UPLO,A,W,N,IL,IU)
      INTEGER,INTENT(IN) :: N, IL, IU
      CHARACTER,INTENT(IN) :: JOBZ, UPLO
      COMPLEX(gq),INTENT(INOUT) :: A(N,N)
      REAL(gq),INTENT(OUT) :: W(N)
! LOCAL
      INTEGER M0
      REAL(gq) VL, VU, ABSTOL
      COMPLEX(gq),ALLOCATABLE :: Z(:,:)
!
      M0 = IU - IL + 1
      ABSTOL = 0._gq
      ALLOCATE(Z(N,M0))
      CALL ZHEEVX_(JOBZ,'I',UPLO,N,A,VL,VU,IL,IU,ABSTOL,M0,W,Z)
      A(:, 1:M0) = Z
      DEALLOCATE(Z)
      RETURN
!
      END SUBROUTINE ZHEEVX_I
!
!******************************************************************
      SUBROUTINE ZHEEVX_(JOBZ,RANGE,UPLO,N,A,VL,VU,IL,IU,ABSTOL,M0,W,Z)
      INTEGER N,M0,IL,IU
      CHARACTER JOBZ,RANGE,UPLO
      REAL(gq) VL,VU,ABSTOL,W(N)
      COMPLEX(gq) A(N,N),Z(N,M0)
! LOCAL
      INTEGER LWORK,INFO,M
      INTEGER,ALLOCATABLE::IWORK(:),IFAIL(:)
      REAL(gq),ALLOCATABLE::RWORK(:)
      COMPLEX(gq),ALLOCATABLE::WORK(:)
!
      LWORK=77*N
      ALLOCATE(WORK(LWORK),RWORK(7*N),IWORK(5*N),IFAIL(N)); WORK=0; RWORK=0; IWORK=0; IFAIL=0
      CALL ZHEEVX(JOBZ,RANGE,UPLO,N,A,N,VL,VU,IL,IU,ABSTOL,M,W,Z,N,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
      IF(INFO/=0)THEN
        WRITE(0,'(" ERROR IN ZHEEVX_: INFO=",I3," IFAIL=")')INFO
        WRITE(0,'(10I4)')IFAIL(1:M); STOP
      ENDIF
      DEALLOCATE(WORK,RWORK,IWORK,IFAIL)
      RETURN
!
      END SUBROUTINE ZHEEVX_
!
!******************************************************************
      SUBROUTINE DSYEVX_I(JOBZ,UPLO,A,W,N,IL,IU)
      INTEGER,INTENT(IN) :: N, IL, IU
      CHARACTER,INTENT(IN) :: JOBZ, UPLO
      REAL(gq),INTENT(INOUT) :: A(N,N)
      REAL(gq),INTENT(OUT) :: W(N)
! LOCAL
      INTEGER M0
      REAL(gq) VL, VU, ABSTOL
      REAL(gq),ALLOCATABLE :: Z(:,:)
!
      M0 = IU - IL + 1
      ABSTOL = 0._gq
      ALLOCATE(Z(N,M0))
      CALL DSYEVX_(JOBZ,'I',UPLO,N,A,VL,VU,IL,IU,ABSTOL,M0,W,Z)
      A(:, 1:M0) = Z
      DEALLOCATE(Z)
      RETURN
!
      END SUBROUTINE DSYEVX_I
!
!******************************************************************
      SUBROUTINE DSYEVX_(JOBZ,RANGE,UPLO,N,A,VL,VU,IL,IU,ABSTOL,M0,W,Z)
      INTEGER N,M0,IL,IU
      CHARACTER JOBZ,RANGE,UPLO
      REAL(gq) VL,VU,ABSTOL,W(N)
      REAL(gq) A(N,N),Z(N,M0)
! LOCAL
      INTEGER LWORK,INFO,M
      INTEGER,ALLOCATABLE::IWORK(:),IFAIL(:)
      REAL(gq),ALLOCATABLE::WORK(:)
!
      LWORK=77*N
      ALLOCATE(WORK(LWORK),IWORK(5*N),IFAIL(N)); WORK=0; IWORK=0; IFAIL=0
      CALL DSYEVX(JOBZ,RANGE,UPLO,N,A,N,VL,VU,IL,IU,ABSTOL,M,W,Z,N,WORK,LWORK,IWORK,IFAIL,INFO)
      IF(INFO/=0)THEN
        WRITE(0,'(" ERROR IN DSTYEVX_: INFO=",I3," IFAIL=")')INFO
        WRITE(0,'(10I4)')IFAIL(1:M); STOP
      ENDIF
      DEALLOCATE(WORK,IWORK,IFAIL)
      RETURN
!
      END SUBROUTINE DSYEVX_
!
!*************************************************************************
! BUG IN ZHEVX? For f orbital n_f=8, it gives unorthogonal eigen vectors!
!*************************************************************************
      SUBROUTINE DDIAG_JZJP(UPLO,RJ2,JZ,JP,VAL1,NDIM)
      CHARACTER UPLO
      INTEGER NDIM
      REAL(gq) JZ(NDIM,NDIM),JP(NDIM,NDIM)
      REAL(gq) RJ2,VAL1(NDIM)
! LOCAL
      INTEGER NV,MJ,IV,IV0,IV1,MULTJ
      REAL(gq) RJ,RMJ
      REAL(gq) RES
!
      RJ=SQRT(RJ2+.25_gq)-0.5_gq; MULTJ=NINT(2*RJ+1)
      NV=NDIM/MULTJ
      CALL DSYEV_('V',UPLO,JZ,VAL1,NDIM)
      RMJ=-RJ
      DO MJ=1,MULTJ-1
      RES=1/SQRT(RJ*(RJ+1)-RMJ*(RMJ+1))
      DO IV=1,NV
        IV0=IV+NV*(MJ-1); IV1=IV+NV*MJ
        CALL DGEMM('N','N',NDIM,1,NDIM,RES,JP,NDIM,JZ(:,IV0),NDIM,D0,JZ(:,IV1),NDIM) ! Jp |N,J,mJ>
        VAL1(IV1)=VAL1(IV0)+1
      ENDDO
      RMJ=RMJ+1
      ENDDO
      RETURN
!
      END SUBROUTINE DDIAG_JZJP
!
!
!*************************************************************************
! BUG IN ZHEVX? For f orbital n_f=8, it gives unorthogonal eigen vectors!
!*************************************************************************
      SUBROUTINE ZDIAG_JZJP(UPLO,RJ2,JZ,JP,VAL1,NDIM)
      CHARACTER UPLO
      INTEGER NDIM
      COMPLEX(gq) JZ(NDIM,NDIM),JP(NDIM,NDIM)
      REAL(gq) RJ2,VAL1(NDIM)
! LOCAL
      INTEGER NV,MJ,IV,IV0,IV1,MULTJ
      REAL(gq) RJ,RMJ
      COMPLEX(gq) ZES
!
      RJ=SQRT(RJ2+.25_gq)-0.5_gq; MULTJ=NINT(2*RJ+1)
      NV=NDIM/MULTJ
      CALL ZHEEV_('V',UPLO,JZ,VAL1,NDIM)
      RMJ=-RJ
      DO MJ=1,MULTJ-1
      ZES=1/SQRT(RJ*(RJ+1)-RMJ*(RMJ+1))
      DO IV=1,NV
        IV0=IV+NV*(MJ-1); IV1=IV+NV*MJ
        CALL ZGEMM('N','N',NDIM,1,NDIM,ZES,JP,NDIM,JZ(:,IV0),NDIM,Z0,JZ(:,IV1),NDIM) ! Jp |N,J,mJ>
        VAL1(IV1)=VAL1(IV0)+1
      ENDDO
      RMJ=RMJ+1
      ENDDO
      RETURN
!
      END SUBROUTINE ZDIAG_JZJP
!
!******************************************************************
! Not as robust as DSYEV, Someitmes may fail
!******************************************************************
      SUBROUTINE DSYEVD_(JOBZ,UPLO,A,W,N)
      INTEGER N
      CHARACTER JOBZ,UPLO
      REAL(gq) A(N,N),W(N)
! LOCAL
      INTEGER LWORK,LIWORK,INFO
      INTEGER,ALLOCATABLE::IWORK(:)
      REAL(gq),ALLOCATABLE::WORK(:)
!
      LWORK=1; LIWORK=1
      ALLOCATE(WORK(LWORK),IWORK(LIWORK))
      CALL DSYEVD(JOBZ,UPLO,N,A,N,W,WORK,-1,IWORK,-1,INFO )
      LWORK=NINT(WORK(1)); LIWORK=IWORK(1)
      DEALLOCATE(WORK,IWORK); ALLOCATE(WORK(LWORK),IWORK(LIWORK))
      CALL DSYEVD(JOBZ,UPLO,N,A,N,W,WORK,LWORK,IWORK,LIWORK,INFO )
      DEALLOCATE(WORK,IWORK)
      IF(INFO.NE.0)THEN
        WRITE(0,'(" ERROR IN DSYEVD_: INFO=",I5," N=",I5," TRY DSYEV!")')INFO,N
        CALL DSYEV_(JOBZ,UPLO,A,W,N)
      ENDIF
      RETURN
!
      END SUBROUTINE DSYEVD_
!
!******************************************************************
      SUBROUTINE DSYEV_(JOBZ,UPLO,A,W,N)
      INTEGER N
      CHARACTER JOBZ,UPLO
      REAL(gq) A(N,N),W(N)
! LOCAL
      INTEGER LWORK,INFO
      REAL(gq),ALLOCATABLE::WORK(:)
!
      LWORK=1; ALLOCATE(WORK(LWORK))
      CALL DSYEV(JOBZ,UPLO,N,A,N,W,WORK,-1,INFO )
      LWORK=NINT(WORK(1))
      DEALLOCATE(WORK); ALLOCATE(WORK(LWORK))
      CALL DSYEV(JOBZ,UPLO,N,A,N,W,WORK,LWORK,INFO )
      DEALLOCATE(WORK)
      IF(INFO.NE.0)THEN
        WRITE(0,'(" ERROR IN DSYEV_: INFO=",I5," N=",I5)')INFO,N; STOP
      ENDIF
      RETURN
!
      END SUBROUTINE DSYEV_
!
!******************************************************************
      SUBROUTINE ZGEEV_(JOBVR,A,W,N)
      CHARACTER JOBVR
      INTEGER N
      COMPLEX(gq) A(N,N),W(N)
! LOCAL
      INTEGER LWORK,INFO
      COMPLEX(gq) VL(1,1)
      REAL(gq) RWORK(2*N)
      COMPLEX(gq),ALLOCATABLE::VR(:,:),WORK(:)
!
      LWORK=16*N
      ALLOCATE(WORK(LWORK),VR(N,N))
      CALL ZGEEV('N',JOBVR,N,A,N,W,VL,1,VR,N,WORK,LWORK,RWORK,INFO)
      IF(INFO.NE.0)THEN
        WRITE(0,'(" ERROR IN ZGEEV_: INFO=",I5," N=",I5)')INFO,N; STOP
      ENDIF
      A=VR
      RETURN
!
      END SUBROUTINE ZGEEV_
!
!=======================================================================
! Inverse of a square matrix
!=======================================================================
      SUBROUTINE DINV_(A,NDIM)
      REAL(gq),intent(inout) :: A(ndim,ndim)
      INTEGER, intent(in)       :: ndim
      ! locals
      INTEGER    :: info, lwork, lda
      INTEGER    :: ipiv(ndim)
      REAL(gq):: work(ndim*64)
      lwork = ndim*64
      lda = ndim
!
      CALL DGETRF( ndim, ndim, A, lda, ipiv, info )
      if (info.ne.0) then
        WRITE(0,*)'dgetrf info=', info; STOP
      endif
      CALL DGETRI( ndim, A, lda, ipiv, work, lwork, info )
      if (info.ne.0) then
        WRITE(0,*)'zgetri info=', info; STOP
      endif
      RETURN
!
      END SUBROUTINE DINV_
!
!=======================================================================
      SUBROUTINE ZINV_(A,NDIM)
      COMPLEX(gq),intent(inout) :: A(ndim,ndim)
      INTEGER, intent(in)       :: ndim
      ! locals
      INTEGER    :: info, lwork, lda
      INTEGER    :: ipiv(ndim)
      COMPLEX(gq):: work(ndim*64)
      lwork = ndim*64
      lda = ndim
!
      CALL ZGETRF( ndim, ndim, A, lda, ipiv, info )
      if (info.ne.0) then
        WRITE(0,*)'zgetrf info=', info; STOP
      endif
      CALL ZGETRI( ndim, A, lda, ipiv, work, lwork, info )
      if (info.ne.0) then
        WRITE(0,*)'zgetri info=', info; STOP
      endif
      RETURN
!
      END SUBROUTINE ZINV_
!
!************************************************************
      SUBROUTINE CHK_UNITARY(A,N,MAXERR)
      INTEGER N
      COMPLEX(gq) A(N,N)
      REAL(gq) MAXERR
! LOCAL
      INTEGER I
      COMPLEX(gq) B(N,N)
!
      B=MATMUL(A,TRANSPOSE(CONJG(A)))
      DO I=1,N; B(I,I)=B(I,I)-1._gq; ENDDO
      MAXERR=MAXVAL(ABS(B))
      RETURN
!
      END SUBROUTINE CHK_UNITARY
!
!************************************************************
! Get Hermitian Matrix expansion coefficients (real)
! MODE>0 : forwatrd
!     <0 : backward
!************************************************************
      SUBROUTINE GET_DHM_EXPAND(A,B,N,NB,C,MODE,LTRANS)
      INTEGER N,NB,MODE
      REAL(gq) A(N,N),B(N,N,NB),C(NB)
      LOGICAL LTRANS
! LOCAL
      INTEGER I
      REAL(gq) B_(N,N)
      REAL(gq),EXTERNAL::DDOT
!
      IF(MODE<0)A=0
      DO I=1,NB
        IF(LTRANS)THEN; B_=TRANSPOSE(B(:,:,I)); ELSE; B_=B(:,:,I); ENDIF
        IF(MODE>0)THEN
          C(I)=DDOT(N*N,B_(1,1),1,A(1,1),1)
        ELSE
          A=A+B_*C(I)
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE GET_DHM_EXPAND
!
!************************************************************
      SUBROUTINE GET_ZHM_EXPAND(A,B,N,NB,C,MODE,LTRANS)
      INTEGER N,NB,MODE
      COMPLEX(gq) A(N,N),B(N,N,NB),C(NB)
      LOGICAL LTRANS
! LOCAL
      INTEGER I
      COMPLEX(gq) B_(N,N)
      COMPLEX(gq),EXTERNAL::ZDOTC
!
      IF(MODE<0)A=0
      DO I=1,NB
        IF(LTRANS)THEN; B_=TRANSPOSE(B(:,:,I)); ELSE; B_=B(:,:,I); ENDIF
        IF(MODE>0)THEN
          C(I)=ZDOTC(N*N,B_(1,1),1,A(1,1),1)
        ELSE
          A=A+B_*C(I)
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE GET_ZHM_EXPAND
!
!************************************************************
      SUBROUTINE ZUHAU(A,U,N,M,UHAU,TRUL,TRUR)
      INTEGER N,M
      COMPLEX(gq) A(N,N),U(N,M)
      COMPLEX(gq),OPTIONAL::UHAU(M,M)
      CHARACTER*1,OPTIONAL::TRUL,TRUR
! LOCAL
      CHARACTER*1 TUL,TUR
      COMPLEX(gq),ALLOCATABLE::A_(:,:),B_(:,:)
!
      IF(PRESENT(TRUL))THEN; TUL=TRUL; ELSE; TUL='C'; ENDIF
      IF(PRESENT(TRUR))THEN; TUR=TRUR; ELSE; TUR='N'; ENDIF
      ALLOCATE(A_(N,M)); A_=0
      CALL ZGEMM('N',TUR,N,M,N,Z1,A,N,U ,N,Z0,A_,N)
      ALLOCATE(B_(M,M)); B_=0
      CALL ZGEMM(TUL,'N',M,M,N,Z1,U,N,A_,N,Z0,B_,M)
      IF(PRESENT(UHAU))THEN
        UHAU=B_
      ELSE
        IF(N/=M)STOP ' ERROR IN ZUHAU: N/=M WHILE NOT PRESET UHAU!'
        A=B_
      ENDIF
      DEALLOCATE(A_,B_)
      RETURN
!
      END SUBROUTINE ZUHAU
!
!************************************************************
      SUBROUTINE DUHAU(A,U,N,M,UHAU,TRUL,TRUR)
      INTEGER N,M
      REAL(gq) A(N,N),U(N,M)
      REAL(gq),OPTIONAL::UHAU(M,M)
      CHARACTER*1,OPTIONAL::TRUL,TRUR
! LOCAL
      CHARACTER*1 TUL,TUR
      REAL(gq),ALLOCATABLE::A_(:,:),B_(:,:)
!
      IF(PRESENT(TRUL))THEN; TUL=TRUL; ELSE; TUL='C'; ENDIF
      IF(PRESENT(TRUR))THEN; TUR=TRUR; ELSE; TUR='N'; ENDIF
      ALLOCATE(A_(N,M)); A_=0
      CALL DGEMM('N',TUR,N,M,N,D1,A,N,U ,N,D0,A_,N)
      ALLOCATE(B_(M,M)); B_=0
      CALL DGEMM(TUL,'N',M,M,N,D1,U,N,A_,N,D0,B_,M)
      IF(PRESENT(UHAU))THEN
        UHAU=B_
      ELSE
        IF(N/=M)STOP ' ERROR IN DUHAU: N/=M WHILE NOT PRESET UHAU!'
        A=B_
      ENDIF
      DEALLOCATE(A_,B_)
      RETURN
!
      END SUBROUTINE DUHAU
!
!************************************************************
      SUBROUTINE ZANMXBMM(TRANSB,A,B,N,M)
      CHARACTER*1 TRANSB
      INTEGER N,M
      COMPLEX(gq) A(N,M),B(M,M)
! LOCAL
      COMPLEX(gq),ALLOCATABLE::A_(:,:)
!
      ALLOCATE(A_(N,M)); A_=0
      CALL ZGEMM('N',TRANSB,N,M,M,Z1,A,N,B,M,Z0,A_,N)
      A=A_
      DEALLOCATE(A_)
      RETURN
!
      END SUBROUTINE ZANMXBMM
!
!************************************************************
      SUBROUTINE DANMXBMM(TRANSB,A,B,N,M)
      CHARACTER*1 TRANSB
      INTEGER N,M
      REAL(gq) A(N,M),B(M,M)
! LOCAL
      REAL(gq),ALLOCATABLE::A_(:,:)
!
      ALLOCATE(A_(N,M)); A_=0
      CALL DGEMM('N',TRANSB,N,M,M,D1,A,N,B,M,D0,A_,N)
      A=A_
      DEALLOCATE(A_)
      RETURN
!
      END SUBROUTINE DANMXBMM
!
!************************************************************
      SUBROUTINE DANNXBNM(TRANSA,A,B,N,M)
      CHARACTER*1 TRANSA
      INTEGER N,M
      REAL(gq) A(N,N),B(N,M)
! LOCAL
      REAL(gq),ALLOCATABLE::B_(:,:)
!
      ALLOCATE(B_(N,M)); B_=0
      CALL DGEMM(TRANSA,'N',N,M,N,Z1,A,N,B,N,Z0,B_,N)
      B=B_
      DEALLOCATE(B_)
      RETURN
!
      END SUBROUTINE DANNXBNM
!
!************************************************************
      SUBROUTINE ZANNXBNM(TRANSA,A,B,N,M)
      CHARACTER*1 TRANSA
      INTEGER N,M
      COMPLEX(gq) A(N,N),B(N,M)
! LOCAL
      COMPLEX(gq),ALLOCATABLE::B_(:,:)
!
      ALLOCATE(B_(N,M)); B_=0
      CALL ZGEMM(TRANSA,'N',N,M,N,Z1,A,N,B,N,Z0,B_,N)
      B=B_
      DEALLOCATE(B_)
      RETURN
!
      END SUBROUTINE ZANNXBNM
!
!************************************************************
      SUBROUTINE ZANNXBNN(TRANSA,TRANSB,A,B,N,MODE)
      CHARACTER*1 TRANSA,TRANSB
      INTEGER N,MODE
      COMPLEX(gq) A(N,N),B(N,N)
! LOCAL
      COMPLEX(gq),ALLOCATABLE::A_(:,:)
!
      ALLOCATE(A_(N,N))
      CALL ZGEMM(TRANSA,TRANSB,N,N,N,Z1,A,N,B,N,Z0,A_,N)
      IF(MODE==1)THEN
        A=A_
      ELSE
        B=A_
      ENDIF
      DEALLOCATE(A_)
      RETURN
!
      END SUBROUTINE ZANNXBNN
!
!************************************************************
      SUBROUTINE DANNXBNN(TRANSA,TRANSB,A,B,N,MODE)
      CHARACTER*1 TRANSA,TRANSB
      INTEGER N,MODE
      REAL(gq) A(N,N),B(N,N)
! LOCAL
      REAL(gq),ALLOCATABLE::A_(:,:)
!
      ALLOCATE(A_(N,N))
      CALL DGEMM(TRANSA,TRANSB,N,N,N,Z1,A,N,B,N,Z0,A_,N)
      IF(MODE==1)THEN
        A=A_
      ELSE
        B=A_
      ENDIF
      DEALLOCATE(A_)
      RETURN
!
      END SUBROUTINE DANNXBNN
!
!************************************************************
      SUBROUTINE TRACE_DANN(A,N,TR)
      INTEGER N
      REAL(gq) A(N,N),TR
! LOCAL
      INTEGER I
!
      TR=0
      DO I=1,N; TR=TR+A(I,I); ENDDO
      RETURN
!
      END SUBROUTINE TRACE_DANN
!
!************************************************************
      SUBROUTINE TRACE_ZANN(A,N,TR)
      INTEGER N
      COMPLEX(gq) A(N,N),TR
! LOCAL
      INTEGER I
!
      TR=0
      DO I=1,N; TR=TR+A(I,I); ENDDO
      RETURN
!
      END SUBROUTINE TRACE_ZANN
!
!************************************************************
      SUBROUTINE TRACE_ZANND(A,N,TR)
      INTEGER N
      COMPLEX(gq) A(N,N)
      REAL(gq) TR
! LOCAL
      INTEGER I
!
      TR=0
      DO I=1,N; TR=TR+REAL(A(I,I),gq); ENDDO
      RETURN
!
      END SUBROUTINE TRACE_ZANND
!
!************************************************************
      SUBROUTINE ZATOFA(A,FA,N,MODE,COEF,LSYM)
      INTEGER N,MODE
      COMPLEX(gq) A(N,N),FA(N,N)
      REAL(gq) COEF
      LOGICAL :: LSYM
! LOCAL
      REAL(gq),ALLOCATABLE::WR(:)
      COMPLEX(gq),ALLOCATABLE::V(:,:),W(:),VINV(:,:)
!
      ALLOCATE(W(N),V(N,N)); W=0
      V=A
      IF(LSYM)THEN
        ALLOCATE(WR(N))
        CALL ZHEEV_('V','L',V,WR,N)
        W=WR; DEALLOCATE(WR)
        CALL ZATOFA1(FA,W,V,N,N,MODE,COEF)
      ELSE
        CALL ZGEEV_('V',V,W,N)
        ALLOCATE(VINV(N,N)); VINV=V
        CALL ZINV_(VINV,N)
        CALL ZATOFA1(FA,W,V,N,N,MODE,COEF,VINV=VINV)
      ENDIF
      DEALLOCATE(W,V)
      RETURN
!
      END SUBROUTINE ZATOFA
!
!************************************************************
      SUBROUTINE DATOFA(A,FA,N,MODE,COEF,LSYM)
      INTEGER N,MODE
      REAL(gq) A(N,N),FA(N,N)
      REAL(gq) COEF
      LOGICAL :: LSYM
! LOCAL
      REAL(gq),ALLOCATABLE::W(:)
      REAL(gq),ALLOCATABLE::V(:,:),VINV(:,:)
!
      ALLOCATE(W(N),V(N,N)); W=0
      V=A
      CALL DSYEV_('V','L',V,W,N)
      IF(LSYM)THEN
        CALL DATOFA1(FA,W,V,N,N,MODE,COEF)
      ELSE
        ALLOCATE(VINV(N,N)); VINV=V
        CALL DINV_(VINV,N)
        CALL DATOFA1(FA,W,V,N,N,MODE,COEF,VINV=VINV)
        DEALLOCATE(VINV)
      ENDIF
      DEALLOCATE(W,V)
      RETURN
!
      END SUBROUTINE DATOFA
!
!************************************************************
      SUBROUTINE DATOFA1(FA,W,V,N,M,MODE,COEF,VINV)
      INTEGER N,M,MODE
      REAL(gq) FA(M,M),W(N),V(M,N)
      REAL(gq) COEF
      REAL(gq),OPTIONAL::VINV(N,M)
! LOCAL
      INTEGER I
      REAL(gq) W_(N),VW(M,N)
!
      W_=W*COEF
      SELECT CASE(MODE)
      CASE(1) ! Exponential
        W_=EXP(W_)
      CASE(2) ! LOG
        DO I=1,N
          IF(W_(I)<RLBOUND)THEN
            W_(I)=RLBOUND
          ENDIF
          W_(I)=LOG(W_(I))
        ENDDO
      CASE(-1) ! Inverse
        DO I=1,N
          IF(ABS(W_(I))<RLBOUND)THEN
            W_(I)=RLBOUND
          ELSEIF(ABS(W_(I))>RUBOUND)THEN
            W_(I)=W_(I)/ABS(W_(I))*RUBOUND
          ENDIF
          W_(I)=1/W_(I)
        ENDDO
      CASE(-12) ! ^(-1/2)
        DO I=1,N
          IF(W_(I)<RLBOUND)THEN
            W_(I)=RLBOUND
          ELSEIF(W_(I)>RUBOUND)THEN
            W_(I)=RUBOUND
          ENDIF
          W_(I)=1/SQRT(W_(I))
        ENDDO
      CASE DEFAULT
        STOP ' ERROR IN DATOFA: UNDEFINED MODE!'
      END SELECT
      DO I=1,N; VW(:,I)=V(:,I)*W_(I); ENDDO
      IF(PRESENT(VINV))THEN
        CALL DGEMM('N','N',M,M,N,D1,VW,M,VINV,M,D0,FA,M)
      ELSE
        CALL DGEMM('N','C',M,M,N,D1,VW,M,V   ,M,D0,FA,M)
      ENDIF
      RETURN
!
      END SUBROUTINE DATOFA1
!
!************************************************************
      SUBROUTINE ZATOFA1(FA,W,V,N,M,MODE,COEF,VINV)
      INTEGER N,M,MODE
      COMPLEX(gq) FA(M,M),V(M,N),W(N)
      REAL(gq) COEF
      COMPLEX(gq),OPTIONAL::VINV(N,M)
! LOCAL
      INTEGER I
      COMPLEX(gq) W_(N)
      COMPLEX(gq) VW(M,N)
!
      W_=W*COEF
      SELECT CASE(MODE)
      CASE(1) ! Exponential
        W_=EXP(W_)
      CASE(2) ! LOG
        DO I=1,N
          IF(ABS(W_(I))<RLBOUND)THEN
            W_(I)=RLBOUND
          ENDIF
          W_(I)=LOG(W_(I))
        ENDDO
      CASE(-1) ! Inverse
        DO I=1,N
          IF(ABS(W_(I))<RLBOUND)THEN
            W_(I)=RLBOUND
          ELSEIF(ABS(W_(I))>RUBOUND)THEN
            W_(I)=W_(I)/ABS(W_(I))*RUBOUND
          ENDIF
          W_(I)=1/W_(I)
        ENDDO
      CASE(-12) ! ^(-1/2)
        DO I=1,N
          IF(ABS(W_(I))<RLBOUND)THEN
            W_(I)=RLBOUND
          ELSEIF(ABS(W_(I))>RUBOUND)THEN
            W_(I)=RUBOUND
          ENDIF
          W_(I)=1/SQRT(W_(I))
        ENDDO
      CASE DEFAULT
        STOP ' ERROR IN DATOFA: UNDEFINED MODE!'
      END SELECT
      DO I=1,N; VW(:,I)=V(:,I)*W_(I); ENDDO
      IF(PRESENT(VINV))THEN
        CALL ZGEMM('N','N',M,M,N,Z1,VW,M,VINV,M,Z0,FA,M)
      ELSE
        CALL ZGEMM('N','C',M,M,N,Z1,VW,M,V   ,M,Z0,FA,M)
      ENDIF
      RETURN
!
      END SUBROUTINE ZATOFA1
!
!*****************************************************************
! A is symmetric, ( \par F(A) / \par d_n --(Hermitian component) )
!*****************************************************************
      SUBROUTINE DPFA_PA(A,PA,H,N,F,FP)
      INTEGER N
      REAL(gq) A(N,N),PA(N,N),H(N,N)
      REAL(gq),EXTERNAL::F,FP
! LOCAL
      REAL(gq) U(N,N),LA(N,N)
      REAL(gq) W(N)
!
      U=A
      CALL DSYEV_('V','L',U,W,N)
      CALL DUHAU(H,U,N,N) ! In eigenstate reprsentation of A
      CALL CALC_LOEWNER(W,LA,N,F,FP)
      PA=LA*H
      U=TRANSPOSE(U)
      CALL DUHAU(PA,U,N,N) ! Back to original representation
      RETURN
!
      END SUBROUTINE DPFA_PA
!
!************************************************************
! A is Hermitian, ( \par F(A) / \par a_IJ )
!************************************************************
      SUBROUTINE ZPFA_PA(A,PA,H,N,F,FP)
      INTEGER N
      COMPLEX(gq) A(N,N),PA(N,N),H(N,N)
      REAL(gq),EXTERNAL::F,FP
! LOCAL
      COMPLEX(gq) U(N,N)
      REAL(gq) W(N),LA(N,N)
!
      U=A
      CALL ZHEEV_('V','L',U,W,N)
      CALL ZUHAU(H,U,N,N) ! In eigenstate reprsentation of A
      CALL CALC_LOEWNER(W,LA,N,F,FP)
      PA=LA*H
      U=TRANSPOSE(CONJG(U))
      CALL ZUHAU(PA,U,N,N) ! Back to original representation 
      RETURN
!
      END SUBROUTINE ZPFA_PA
!
!************************************************************
      SUBROUTINE DH_REGULARIZE(H,N,TOL)
      INTEGER N
      REAL(gq) H(N,N),TOL
! LOCAL
      INTEGER I
      REAL(gq) W(N),VW(N,N)
!      
      CALL DSYEV_('V','L',H,W,N)
      DO I=1,N
        W(I)=MAX(W(I),TOL)
        VW(:,I)=H(:,I)*W(I)
      ENDDO
      CALL DANNXBNN('N','T',VW,H,N,2)
      RETURN
!
      END SUBROUTINE DH_REGULARIZE
!
!************************************************************
      SUBROUTINE ZH_REGULARIZE(H,N,TOL)
      INTEGER N
      COMPLEX(gq) H(N,N)
      REAL(gq) TOL
! LOCAL
      INTEGER I
      REAL(gq) W(N)
      COMPLEX(gq) VW(N,N)
!
      CALL ZHEEV_('V','L',H,W,N)
      DO I=1,N
        W(I)=MAX(W(I),TOL)
        VW(:,I)=H(:,I)*W(I)
      ENDDO
      CALL ZANNXBNN('N','C',VW,H,N,2)
      RETURN
!
      END SUBROUTINE ZH_REGULARIZE
!
!************************************************************
! LOEWNER MATRIX OF HERMITIAN A GIVEN ITS EIGEN VALUES W.
!************************************************************
      SUBROUTINE CALC_LOEWNER(W,LA,N,F,FP)
      INTEGER N
      REAL(gq) W(N),LA(N,N)
      REAL(gq),EXTERNAL::F,FP
! LOCAL
      INTEGER I,J
!
      DO I=1,N; DO J=1,N
      IF(ABS(W(I)-W(J))<1.E-10_gq)THEN
        LA(I,J)=FP(W(I))
      ELSE
        LA(I,J)=(F(W(I))-F(W(J)))/(W(I)-W(J))
      ENDIF
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE CALC_LOEWNER
!
!************************************************************
      FUNCTION DSIMIX(X)
      REAL(gq),INTENT(IN) :: X
      REAL(gq) :: DSIMIX
! LOCAL
      REAL(gq) X_
!
      IF(ABS(X)<1.E-8_gq)THEN
        X_=1.E-8_gq
      ELSEIF(ABS(X-1)<1.E-8_gq)THEN
        X_=1-1.E-8_gq
      ELSEIF(X<0.OR.X>1)THEN
        WRITE(0,*) " FETAL ERROR IN DSIMIX: ILLEGAL X=",X; STOP
      ELSE
        X_=X
      ENDIF
      DSIMIX=SQRT(X_*(1-X_))
      RETURN
!
      END FUNCTION DSIMIX
!
!************************************************************
      FUNCTION DPSIMIX(X)
      REAL(gq),INTENT(IN) :: X
      REAL(gq) :: DPSIMIX
! LOCAL
      REAL(gq) X_
!
      IF(ABS(X)<1.E-8_gq)THEN
        X_=1.E-8_gq
      ELSEIF(ABS(X-1)<1.E-8_gq)THEN
        X_=1-1.E-8_gq
      ELSEIF(X<0.OR.X>1)THEN
        WRITE(0,*) " FETAL ERROR IN DSIMIX: ILLEGAL X=",X; STOP
      ELSE
        X_=X
      ENDIF
      DPSIMIX=(.5_gq-X)/SQRT(X_*(1-X_))
      RETURN
!
      END FUNCTION DPSIMIX
!
!************************************************************
      SUBROUTINE ZABS_ORDER(V,E,N)
      INTEGER N
      REAL(gq) E(N)
      COMPLEX(gq) V(N,N)
! LOCAL
      INTEGER I,ID(N)
      REAL(gq) EBUF(N)
      COMPLEX(gq),ALLOCATABLE::ZBUF(:,:)
!
      EBUF=ABS(E); ID=0
      CALL DSORT(N,EBUF,ID,.TRUE.)
      ALLOCATE(ZBUF(N,N))
      EBUF=E; ZBUF=V
      DO I=1,N
        E(I)=EBUF(ID(I))
        V(:,I)=ZBUF(:,ID(I))
      ENDDO
      DEALLOCATE(ZBUF)
      RETURN
!
      END SUBROUTINE ZABS_ORDER
!
!=======================================================================
! sorts RA in descending/ascending order, and rearanges an index array IB
!=======================================================================
      SUBROUTINE DSORT(N,RA,IB,LDSC)
      INTEGER N
      REAL(gq) RA(N)
      INTEGER IB(N)
      LOGICAL LDSC
! LOCAL
      REAL(gq) RRA,FAK
      INTEGER IIB,L,IR,J,I
!
      DO I=1,N; IB(I)=I; ENDDO
      IF (N<=1) RETURN
      IF(LDSC)THEN; FAK=1._gq; ELSE; FAK=-1._gq; ENDIF
      L=N/2+1; IR=N
10    CONTINUE
      IF(L.GT.1)THEN
         L=L-1
         RRA=RA(L); IIB=IB(L)
      ELSE
         RRA=RA(IR); IIB=IB(IR)
         RA(IR)=RA(1); IB(IR)=IB(1)
         IR=IR-1
         IF(IR.EQ.1)THEN
            RA(1)=RRA; IB(1)=IIB
            RETURN
         ENDIF
      ENDIF
      I=L; J=L+L
20    IF(J.LE.IR)THEN
         IF(J.LT.IR)THEN
            IF((RA(J)-RA(J+1))*FAK.GT.0)J=J+1
         ENDIF
         IF((RRA-RA(J))*FAK.GT.0)THEN
            RA(I)=RA(J); IB(I)=IB(J)
            I=J; J=J+J
         ELSE
            J=IR+1
         ENDIF
         GO TO 20
      ENDIF
      RA(I)=RRA; IB(I)=IIB
      GO TO 10
      RETURN
!
      END SUBROUTINE DSORT
!
!************************************************************
      SUBROUTINE LOCATE_MI(M,N,M1,LM1)
      INTEGER N,M(N),M1,LM1
! LOCAL
      INTEGER I
!
      LM1=0
      DO I=1,N
        IF(M1.EQ.M(I))THEN
          LM1=I; GOTO 10
        ENDIF
      ENDDO
10    CONTINUE
      RETURN
!
      END SUBROUTINE LOCATE_MI
!
!************************************************************
      FUNCTION L22L(L2)
      REAL(gq) L2,L22L
!
      L22L=SQRT(L2+0.25_gq)-0.5_gq
      RETURN
!
      END FUNCTION L22L
!
!************************************************************
      FUNCTION FILE_NAME(NI,IA1,IA2,MODE)
      INTEGER NI,IA1,IA2,MODE
      CHARACTER*77 FILE_NAME
! LOCAL
      CHARACTER*7 STR1,STR2,STR3
!
      WRITE(STR1,'(I7)')NI; WRITE(STR2,'(I7)')IA1; WRITE(STR3,'(I7)')IA2
      SELECT CASE(MODE)
      CASE(-11)
        FILE_NAME='BNDU_'//TRIM(ADJUSTL(STR1))//'.TXT'
      CASE(-4)
        FILE_NAME='KSWU_'//TRIM(ADJUSTL(STR1))//'.DAT'
      CASE(-3)
        FILE_NAME='GMPI_'//TRIM(ADJUSTL(STR1))//'.INP'
      CASE(-2)
        FILE_NAME='KSWT_'//TRIM(ADJUSTL(STR1))//'.DAT'
      CASE(-1)
        FILE_NAME='BNDU_'//TRIM(ADJUSTL(STR1))//'.INP'
      CASE(1)
        FILE_NAME='GLU2_'//TRIM(ADJUSTL(STR1))//'.TMP'
      CASE(2)
        FILE_NAME='GLEVEC_'//TRIM(ADJUSTL(STR1))//'.VSP'
      CASE(102)
        FILE_NAME='GLEVEC_'//TRIM(ADJUSTL(STR1))//'.SYM'
      CASE(202)
        FILE_NAME='GLEVEC_'//TRIM(ADJUSTL(STR1))//'.SYMTXT'
      CASE(4)
        FILE_NAME='GLSKIJ_'//TRIM(ADJUSTL(STR1))//'.VSP'
      CASE(104)
        FILE_NAME='GLSKIJ_'//TRIM(ADJUSTL(STR1))//'.SYM'
      CASE(204)
        FILE_NAME='GLSKIJ_'//TRIM(ADJUSTL(STR1))//'.SYMTXT'
      CASE(5)
        FILE_NAME='GLNVAR_KKH_'//TRIM(ADJUSTL(STR1))//'_'// &
            &TRIM(ADJUSTL(STR2))//'_'//TRIM(ADJUSTL(STR3))//'.VSP'
      CASE(6)
        FILE_NAME='GLMKKP_'//TRIM(ADJUSTL(STR1))//'_'// &
            &TRIM(ADJUSTL(STR2))//'_'//TRIM(ADJUSTL(STR3))//'.VSP'
      CASE(7)
        FILE_NAME='GLUKK_'//TRIM(ADJUSTL(STR1))//'.TMP'
      CASE(8)
        FILE_NAME='GLH_'//TRIM(ADJUSTL(STR1))//'.h5'
      CASE(11)
        FILE_NAME='GLS2_'//TRIM(ADJUSTL(STR1))//'.VSP'
      CASE(12)
        FILE_NAME='GLL2_'//TRIM(ADJUSTL(STR1))//'.VSP'
      CASE(13)
        FILE_NAME='GLJ2_'//TRIM(ADJUSTL(STR1))//'.VSP'
      CASE(15)
        FILE_NAME='GLPJC_'//TRIM(ADJUSTL(STR1))//'.OUT'
      CASE(1501)
        FILE_NAME='GLPJC01_'//TRIM(ADJUSTL(STR1))//'.OUT'
      CASE(16)
        FILE_NAME='GLEVEC_'//TRIM(ADJUSTL(STR1))//'.VSR'
      CASE(17)
        FILE_NAME='GLSECN_'//TRIM(ADJUSTL(STR1))//'.VSR'
      CASE(18)
        FILE_NAME='GLSECJ_'//TRIM(ADJUSTL(STR1))//'.VSR'
      CASE(19)
        FILE_NAME='GF_'//TRIM(ADJUSTL(STR1))//'_'//TRIM(ADJUSTL(STR2))//'_'//TRIM(ADJUSTL(STR3))//'.DAT'
      END SELECT
      RETURN
!
      END FUNCTION FILE_NAME
!
!***********************************************************************
      FUNCTION int_to_str(i)
      integer i
      character(len = 77) int_to_str
!
      write(int_to_str,*)i
      int_to_str = adjustl(int_to_str)
      return 
!
      end FUNCTION int_to_str
!
!***********************************************************************
      FUNCTION FERMI_FUN(X)
      REAL(gq) X,FERMI_FUN
!
      IF(X<-200._gq)THEN
        FERMI_FUN=1._gq
      ELSEIF(X<200._gq)THEN
        FERMI_FUN=1._gq/(EXP(X)+1._gq)
      ELSE
        FERMI_FUN=0._gq
      ENDIF
      RETURN
!
      END FUNCTION FERMI_FUN
!
!***********************************************************************
      FUNCTION GAUSS_FUN(X)
      REAL(gq) X,GAUSS_FUN
!
      IF(X<-7._gq)THEN
        GAUSS_FUN=2._gq
      ELSEIF(X<0._gq)THEN
        GAUSS_FUN=2-ERFC(-X)
      ELSEIF(X<7._gq)THEN
        GAUSS_FUN=ERFC(X)
      ELSE
        GAUSS_FUN=0._gq
      ENDIF
      GAUSS_FUN=GAUSS_FUN/2
      RETURN
!
      END FUNCTION GAUSS_FUN
!
!***********************************************************************
      FUNCTION MAX_INTERVAL(IARRAY,N)
      INTEGER N,MAX_INTERVAL,IARRAY(N)
! LOCAL
      INTEGER I
!
      MAX_INTERVAL=0
      DO I=2,N; MAX_INTERVAL=MAX(IARRAY(I)-IARRAY(I-1),MAX_INTERVAL); ENDDO
      RETURN
!
      END FUNCTION MAX_INTERVAL
!
!***********************************************************************
      SUBROUTINE SET_RANGE(X0,X1,X,N)
      INTEGER N
      REAL(gq) X0,X1,X(N)
! LOCAL
      INTEGER I
      REAL(gq) DX
!
      DX=(X1-X0)/(N-1)
      DO I=1,N
        X(I)=X0+(I-1)*DX
      ENDDO
      RETURN
!
      END SUBROUTINE SET_RANGE
!
!***********************************************************************
      SUBROUTINE SET_LINEAR_SIMP(WT,N,DX)
      INTEGER N
      REAL(8) WT(N),DX
! LOCAL
      INTEGER I
      WT=0
      IF(MOD(N,2)==0)THEN
        STOP ' ERROR IN SET_LINEAR_SIMP: ASSUME ODD POINTS!'
      ENDIF
      DO I=N,3,-2
        WT(I)  =DX/3._8+WT(I)
        WT(I-1)=DX/3._8*4
        WT(I-2)=DX/3._8
      ENDDO
      RETURN
      END SUBROUTINE SET_LINEAR_SIMP
!
!************************************************************
      FUNCTION LKEY_EXIST(LINE, KEY)
      CHARACTER(*) LINE, KEY
      LOGICAL LKEY_EXIST
!
      LKEY_EXIST=.FALSE.
      IF(INDEX(LINE, KEY)>0.AND.(INDEX(LINE,'#')<=0.OR. &
        &INDEX(LINE,'#')>INDEX(LINE, KEY)))THEN
        LKEY_EXIST=.TRUE.
      ENDIF
      RETURN
!
      END FUNCTION LKEY_EXIST
!
!************************************************************
      SUBROUTINE GET_IMIN1(I,N,IMIN)
      INTEGER N,I(N),IMIN
! LOCAL
      INTEGER J
!
      IMIN=100
      DO J=1,N; IF(I(J)<=0)CYCLE; IMIN=MIN(IMIN,I(J)); ENDDO
      RETURN
!
      END SUBROUTINE GET_IMIN1
!
!
      END MODULE GUTIL
!
!*************************************************************
! Down-hill simplex formalism
!*************************************************************
      SUBROUTINE AMOEBA(P,Y,NDIM,FTOL,FUNK,ITMAX,ITER)
      IMPLICIT REAL(8)(A-H,O-Z)
      INTEGER NDIM,ITMAX,ITER
      REAL(8) P(NDIM+1,NDIM),Y(NDIM+1),FTOL
      EXTERNAL :: FUNK
! LOCAL
      PARAMETER (ALPHA=1.0_8,BETA=0.5_8,GAMMA=2.0_8)
      REAL(8),ALLOCATABLE::PR(:),PRR(:),PBAR(:)
!
      ALLOCATE(PR(NDIM),PRR(NDIM),PBAR(NDIM)); PR=0; PRR=0; PBAR=0
      MPTS=NDIM+1
      ITER=0
1     ILO=1
      IF(Y(1).GT.Y(2))THEN
        IHI=1
        INHI=2
      ELSE
        IHI=2
        INHI=1
      ENDIF
      DO 11 I=1,MPTS
        IF(Y(I).LT.Y(ILO)) ILO=I
        IF(Y(I).GT.Y(IHI))THEN
          INHI=IHI
          IHI=I
        ELSE IF(Y(I).GT.Y(INHI))THEN
          IF(I.NE.IHI) INHI=I
        ENDIF
11    CONTINUE
!
!      RTOL=2.*ABS(Y(IHI)-Y(ILO))/MAX(1.E-12_8,(ABS(Y(IHI))+ABS(Y(ILO))))
      RTOL=ABS(Y(IHI)-Y(ILO))
!
!      WRITE(*,'(" ITER=",I6," FTOL=",E10.2," Y(IHI/ILO)=",2F20.12)')ITER,RTOL,Y(IHI),Y(ILO)
!
      IF(RTOL.LT.FTOL)THEN
!        WRITE(*,'(" ITER=",I6," FTOL=",E10.2," Y(IHI/ILO)=",2F20.12)')ITER,RTOL,Y(IHI),Y(ILO)
        RETURN
      ENDIF
      IF(ITER.EQ.ITMAX) WRITE(0,'("WARNING: Amoeba exceeding maximum iterations.")')
      ITER=ITER+1
      DO 12 J=1,NDIM
        PBAR(J)=0.
12    CONTINUE
      DO 14 I=1,MPTS
        IF(I.NE.IHI)THEN
          DO 13 J=1,NDIM
            PBAR(J)=PBAR(J)+P(I,J)
13        CONTINUE
        ENDIF
14    CONTINUE
      DO 15 J=1,NDIM
        PBAR(J)=PBAR(J)/NDIM
        PR(J)=(1.+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
15    CONTINUE
      CALL FUNK(PR,NDIM,YPR)
!      YPR=FUNK(PR)
      IF(YPR.LE.Y(ILO))THEN
        DO 16 J=1,NDIM
          PRR(J)=GAMMA*PR(J)+(1.-GAMMA)*PBAR(J)
16      CONTINUE
        CALL FUNK(PRR,NDIM,YPRR)
!        YPRR=FUNK(PRR)
        IF(YPRR.LT.Y(ILO))THEN
          DO 17 J=1,NDIM
            P(IHI,J)=PRR(J)
17        CONTINUE
          Y(IHI)=YPRR
        ELSE
          DO 18 J=1,NDIM
            P(IHI,J)=PR(J)
18        CONTINUE
          Y(IHI)=YPR
        ENDIF
      ELSE IF(YPR.GE.Y(INHI))THEN
        IF(YPR.LT.Y(IHI))THEN
          DO 19 J=1,NDIM
            P(IHI,J)=PR(J)
19        CONTINUE
          Y(IHI)=YPR
        ENDIF
        DO 21 J=1,NDIM
          PRR(J)=BETA*P(IHI,J)+(1.-BETA)*PBAR(J)
21      CONTINUE
        CALL FUNK(PRR,NDIM,YPRR)
!        YPRR=FUNK(PRR)
        IF(YPRR.LT.Y(IHI))THEN
          DO 22 J=1,NDIM
            P(IHI,J)=PRR(J)
22        CONTINUE
          Y(IHI)=YPRR
        ELSE
          DO 24 I=1,MPTS
            IF(I.NE.ILO)THEN
              DO 23 J=1,NDIM
                PR(J)=0.5*(P(I,J)+P(ILO,J))
                P(I,J)=PR(J)
23            CONTINUE
              CALL FUNK(PR,NDIM,Y(I))
!              Y(I)=FUNK(PR)
            ENDIF
24        CONTINUE
        ENDIF
      ELSE
        DO 25 J=1,NDIM
          P(IHI,J)=PR(J)
25      CONTINUE
        Y(IHI)=YPR
      ENDIF
      GO TO 1
!
      END SUBROUTINE AMOEBA

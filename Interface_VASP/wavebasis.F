#include "symbol.inc"
MODULE WAVEBASIS
      USE prec; USE wave
      IMPLICIT NONE

      TYPE WAVE_BASIS
        INTEGER LWREAD
        TYPE(wavespin) :: W
        COMPLEX(q),POINTER :: U(:,:,:,:) ! <W | eigen vector>
        COMPLEX(q),POINTER :: H_FOCK(:,:,:,:),H_CMR(:,:,:,:)
        INTEGER MODE  ! MODE=1, DUMP W_BASIS; MODE=2: USE W_BASIS
      END TYPE WAVE_BASIS

      TYPE(WAVE_BASIS) :: W_BASIS

CONTAINS

!*****************************************************
      SUBROUTINE WRT_WAVE_BARE()
      INTEGER ISP,NK,NB
      INTEGER,PARAMETER :: IU=90

      OPEN(IU,FILE='WAVEBARE',STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      DO ISP=1,W_BASIS%W%WDES%ISPIN; DO NK=1,W_BASIS%W%WDES%NKPTS
      WRITE(IU)W_BASIS%W%FERTOT(:,NK,ISP)
      DO NB=1,W_BASIS%W%WDES%NBANDS
        WRITE(IU)W_BASIS%W%CPTWFP(:,NB,NK,ISP)
      ENDDO; ENDDO; ENDDO
      CLOSE(IU)
      RETURN

      END SUBROUTINE WRT_WAVE_BARE
!
!*****************************************************
      SUBROUTINE READ_WAVE_BARE()
      INTEGER ISP,NK,NB
      INTEGER,PARAMETER :: IU=90

      W_BASIS%LWREAD=0
      OPEN(IU,FILE='WAVEBARE',STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL",ERR=100)
      DO ISP=1,W_BASIS%W%WDES%ISPIN; DO NK=1,W_BASIS%W%WDES%NKPTS
      READ(IU)W_BASIS%W%FERTOT(:,NK,ISP)
      DO NB=1,W_BASIS%W%WDES%NBANDS
        READ(IU)W_BASIS%W%CPTWFP(:,NB,NK,ISP)
      ENDDO; ENDDO; ENDDO
      CLOSE(IU)
      W_BASIS%LWREAD=1
100   CONTINUE
      RETURN

      END SUBROUTINE READ_WAVE_BARE
!
!*****************************************************
      SUBROUTINE INI_WAVE_BASIS(WDES,MODE)
      INTEGER MODE
      TYPE (wavedes),TARGET :: WDES
! LOCAL
      INTEGER NRPLWV,NBANDS,NKPTS,ISPIN,NPROD

      IF(W_BASIS%MODE/=2)RETURN
      NRPLWV=WDES%NRPLWV
      NPROD =WDES%NPROD
      NKPTS =WDES%NKPTS
      NBANDS=WDES%NBANDS
      ISPIN =WDES%ISPIN
      IF(MODE==2.OR.MODE==12)THEN
        W_BASIS%W%WDES=>WDES
        ALLOCATE(W_BASIS%W%CPTWFP(NRPLWV,NBANDS,NKPTS,ISPIN), &
                &W_BASIS%W%CPROJ (NPROD, NBANDS,NKPTS,ISPIN), &
                &W_BASIS%W%FERTOT(NBANDS,NKPTS,ISPIN), &
                &W_BASIS%W%CELEN (NBANDS,NKPTS,ISPIN), &
                &W_BASIS%H_FOCK(NBANDS,NBANDS,NKPTS,ISPIN), &
                &W_BASIS%H_CMR (NBANDS,NBANDS,NKPTS,ISPIN) )
        W_BASIS%W%FERWE=>W_BASIS%W%FERTOT
        CALL READ_WAVE_BARE()
      ENDIF
      ALLOCATE(W_BASIS%U(NBANDS,NBANDS,NKPTS,ISPIN)); W_BASIS%U=0
      RETURN

      END SUBROUTINE INI_WAVE_BASIS

!*****************************************************
      SUBROUTINE SOLVE_HALL_WB(WDES,W,INFO,GRID,LATT_CUR,NONLR_S,NONL_S,P, &
                         &CDIJ,CQIJ,SV,EXHF,LMDIM)
      USE TB_MODULE; USE wave; USE mgrid ; USE hamil; USE nonl_high; USE lattice; USE pseudo; USE base
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W
      TYPE (info_struct) INFO
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (potcar)      P(WDES%NTYP)
      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)
      INTEGER LMDIM
      REAL(q) :: EXHF
! LOCAL
      INTEGER NK,ISP

      EXHF=0
      DO ISP=1,WDES%ISPIN; DO NK=1,WDES%NKPTS
        CALL SOLVE_HK_WB(WDES,W,INFO,GRID,LATT_CUR,NONLR_S,NONL_S,P, &
                        &CDIJ,CQIJ,SV,EXHF,LMDIM,NK,ISP)
      ENDDO; ENDDO
      DO ISP=1,WDES%ISPIN; DO NK=1,WDES%NKPTS
        CALL ZANMXBMM('N',W%CPTWFP(:,:,NK,ISP),W_BASIS%U(:,:,NK,ISP),WDES%NRPLWV,WDES%NBANDS)
        CALL ZANMXBMM('N',W%CPROJ (:,:,NK,ISP),W_BASIS%U(:,:,NK,ISP),WDES%NPROD ,WDES%NBANDS)
      ENDDO; ENDDO
      IF(TB%KEY==2.OR.TB%KEY==12)THEN
      DO ISP=1,WDES%ISPIN; DO NK=1,WDES%NKPTS
        CALL ZUHAU(W_BASIS%H_CMR (:,:,NK,ISP),W_BASIS%U(:,:,NK,ISP),WDES%NBANDS,WDES%NBANDS)
        CALL ZUHAU(W_BASIS%H_FOCK(:,:,NK,ISP),W_BASIS%U(:,:,NK,ISP),WDES%NBANDS,WDES%NBANDS)
      ENDDO; ENDDO
      ENDIF
      RETURN

      END SUBROUTINE SOLVE_HALL_WB

!*****************************************************
      SUBROUTINE SOLVE_HK_WB(WDES,W,INFO,GRID,LATT_CUR,NONLR_S,NONL_S,P, &
                            &CDIJ,CQIJ,SV,EXHF,LMDIM,NK,ISP)
      USE TB_MODULE; USE wave; USE mgrid ; USE hamil; USE nonl_high; USE lattice; USE pseudo; USE base
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W
      TYPE (info_struct) INFO
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (potcar)      P(WDES%NTYP)
      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)
      INTEGER LMDIM,NK,ISP
      REAL(q) :: EXHF
! LOCAL
      INTEGER NBANDS
      COMPLEX(q),POINTER :: HK(:,:),H_FOCK(:,:),H_CMR(:,:)
      REAL(q) :: EVAL(WDES%NBANDS)
    
      NBANDS=WDES%NBANDS
      HK    =>W_BASIS%U     (:,:,NK,ISP)
      IF(TB%KEY==2.OR.TB%KEY==12)THEN
      H_FOCK=>W_BASIS%H_FOCK(:,:,NK,ISP)
      H_CMR =>W_BASIS%H_CMR (:,:,NK,ISP)
      ELSE
      H_FOCK=>W_BASIS%U     (:,:,NK,ISP)
      H_CMR =>W_BASIS%U     (:,:,NK,ISP)
      ENDIF
      CALL SET_HK_WB(WDES,W,INFO,GRID,LATT_CUR,NONLR_S,NONL_S,P, &
                    &CDIJ,CQIJ,SV,EXHF,LMDIM,NK,ISP,HK,H_FOCK,H_CMR)
      CALL ZHEEV_('V','L',HK,EVAL,NBANDS)
      W%CELEN(:,NK,ISP)=DCMPLX(EVAL)
      NULLIFY(H_CMR,H_FOCK,HK)
      RETURN
    
      END SUBROUTINE SOLVE_HK_WB

!*****************************************************
      SUBROUTINE SET_HK_WB(WDES,W,INFO,GRID,LATT_CUR,NONLR_S,NONL_S,P, &
                         &CDIJ,CQIJ,SV,EXHF,LMDIM,NK,ISP,HK,H_FOCK,H_CMR)
      USE TB_MODULE; USE wave; USE mgrid ; USE hamil; USE nonl_high; USE lattice; USE pseudo; USE base; USE fock
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W
      TYPE (info_struct) INFO    
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR 
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (potcar)      P(WDES%NTYP)
      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)
      COMPLEX(q)::HK(WDES%NBANDS,WDES%NBANDS),H_FOCK(WDES%NBANDS,WDES%NBANDS),H_CMR(WDES%NBANDS,WDES%NBANDS)
      INTEGER LMDIM,NK,ISP
      REAL(q) :: EXHF
! LOCAL
      TYPE (wavefun1) :: W1(1)
      TYPE (wavedes1) :: WDES1
      INTEGER NBANDS,NB,NPL,MB,MPLWV,NACT
      REAL(q) EVALUE(1)
      COMPLEX(q) CDCHF
      TYPE (wavefuna) :: WA
    
      IF((TB%KEY==2.OR.TB%KEY==12).AND.(.NOT.LHFCALC))THEN
        WRITE(0,'(" ERROR IN SETTING TB%KEY==2.OR.TB%KEY==12 WHILE LHFCALC=.F.")')
        STOP
      ENDIF
      IF(LHFCALC)THEN; NACT=2; ELSE; NACT=1; ENDIF
      CALL SETWDES(W%WDES,WDES1,NK)
      WDES1%NBANDS=1; NPL=WDES1%NGVECTOR
      IF (INFO%LREAL) THEN
         CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,W%WDES)
      ELSE
         CALL PHASE(W%WDES,NONL_S,NK)
      ENDIF
      NBANDS=W%WDES%NBANDS
      HK=0; H_FOCK=0; H_CMR=0; EVALUE=0
      CALL NEWWAVA(WA,WDES1,NACT)
      MPLWV=WDES1%GRID%MPLWV*WDES1%NRSPINORS
      ALLOCATE(W1(1)%CR(MPLWV))
      DO NB=1,NBANDS
        CALL SETWAV(W,W1(1),WDES1,NB,ISP)
        CALL FFTWAV_W1(W1(1))
        WA%CPTWFP=0; WA%CPROJ=0
        CALL HAMILTMU(WDES1,W1,NONLR_S,NONL_S,EVALUE, &
            &CDIJ,CQIJ,SV,ISP,ELEMENTS(WA,1,1))
        IF(LHFCALC)THEN
          CDCHF=0
          CALL FOCK_ACC(GRID,LMDIM,LATT_CUR,W, &
                       &NONLR_S,NONL_S,NK,ISP,W1(1)%NB,1, &
                       &WA%CW(:,2:2),P,CQIJ,CDCHF,W1)
          EXHF=EXHF+CDCHF
          CALL W1_DAXPY(ELEMENT(WA,2),1.0_q,ELEMENT(WA,1))
        ENDIF
        IF(TB%KEY==2.OR.TB%KEY==12)THEN
          CDCHF=0
          CALL FOCK_ACC(GRID,LMDIM,LATT_CUR,W_BASIS%W, &
                       &NONLR_S,NONL_S,NK,ISP,W1(1)%NB,1, &
                       &WA%CW(:,2:2),P,CQIJ,CDCHF,W1)
!          EXHF=EXHF-CDCHF; WRONG
        ENDIF
        DO MB=1,NB
          HK    (MB,NB)=W1_DOT(ELEMENT(W,WDES1,MB,ISP),ELEMENT(WA,1))
          IF(TB%KEY==2.OR.TB%KEY==12)THEN
          H_FOCK(MB,NB)=-W1_DOT(ELEMENT(W,WDES1,MB,ISP),ELEMENT(WA,2))
          H_CMR (MB,NB)= HK(MB,NB)
          HK    (MB,NB)= HK(MB,NB)+H_FOCK(MB,NB)
          ENDIF
          IF(MB==NB)CYCLE
          HK    (NB,MB)=CONJG(HK    (MB,NB))
          IF(TB%KEY==2.OR.TB%KEY==12)THEN
          H_FOCK(NB,MB)=CONJG(H_FOCK(MB,NB))
          H_CMR (NB,MB)=CONJG(H_CMR (MB,NB))
          ENDIF
        ENDDO
        IF(TB%KEY==2.OR.TB%KEY==12)THEN
! Now calculate CDCHF for bare Fock
        IF(ABS(W_BASIS%W%FERWE(NB,NK,ISP)/W_BASIS%W%FERWE(1,NK,ISP))<1.D-6)THEN
          CYCLE
        ENDIF
        CALL SETWAV(W_BASIS%W,W1(1),WDES1,NB,ISP)
        CALL FFTWAV_W1(W1(1))
        CDCHF=0
        CALL FOCK_ACC(GRID,LMDIM,LATT_CUR,W_BASIS%W, &
                     &NONLR_S,NONL_S,NK,ISP,W1(1)%NB,1, &
                     &WA%CW(:,2:2),P,CQIJ,CDCHF,W1)
        EXHF=EXHF-CDCHF
        ENDIF
      ENDDO
      CALL DELWAVA(WA)
      DEALLOCATE(W1(1)%CR)
      RETURN
    
      END SUBROUTINE SET_HK_WB

!************************************************************
      SUBROUTINE ZANMXBMM(TRANSB,A,B,N,M,AB)
      USE prec
      IMPLICIT NONE

      CHARACTER*1 TRANSB
      INTEGER N,M
      COMPLEX(q) A(N,M),B(M,M)
      COMPLEX(q),OPTIONAL :: AB(N,M)
! LOCAL
      COMPLEX(q),ALLOCATABLE::A_(:,:)
      COMPLEX(q),PARAMETER :: Z1=(1._q,0._q),Z0=(0._q,0._q)

      ALLOCATE(A_(N,M)); A_=0
      CALL ZGEMM('N',TRANSB,N,M,M,Z1,A,N,B,M,Z0,A_,N)
      IF(PRESENT(AB))THEN
        AB=A_
      ELSE
        A=A_
      ENDIF
      DEALLOCATE(A_)
      RETURN

      END SUBROUTINE ZANMXBMM
!
!************************************************************
      SUBROUTINE ZUHAU(A,U,N,M,UHAU,TRUL,TRUR)
      INTEGER N,M
      COMPLEX(q) A(N,N),U(N,M)
      COMPLEX(q),OPTIONAL::UHAU(M,M)
      CHARACTER*1,OPTIONAL::TRUL,TRUR
! LOCAL
      CHARACTER*1 TUL,TUR
      COMPLEX(q),ALLOCATABLE::A_(:,:),B_(:,:)
      COMPLEX(q),PARAMETER :: Z1=(1._q,0._q),Z0=(0._q,0._q)
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
        IF(N/=M)THEN
          STOP ' ERROR IN ZUHAU: N/=M WHILE NOT PRESET UHAU!'
        ENDIF
        A=B_
      ENDIF
      DEALLOCATE(A_,B_)
      RETURN
!
      END SUBROUTINE ZUHAU

END MODULE WAVEBASIS

!******************************************************************
      SUBROUTINE ZHEEV_(JOBZ,UPLO,A,W,N)
      USE prec
      IMPLICIT NONE

      INTEGER N
      CHARACTER JOBZ,UPLO
      COMPLEX(q) A(N,N)
      REAL(q)    W(N)
! Local
      INTEGER INFO,LWORK
      REAL(8),ALLOCATABLE :: RWORK(:)
      COMPLEX(8),ALLOCATABLE :: WORK(:)

      LWORK=32*N
      ALLOCATE(RWORK(max(1,3*N-2)),WORK(LWORK))
      CALL ZHEEV(JOBZ,UPLO,N,A,N,W,WORK,LWORK,RWORK,INFO)
      IF(INFO.NE.0)THEN
        WRITE(0,'(" ERROR IN ZHEEV_: INFO=",I5," N=",I5)')INFO,N; STOP
      ENDIF
      DEALLOCATE(RWORK,WORK)
      RETURN

      END SUBROUTINE ZHEEV_

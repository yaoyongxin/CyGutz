!**************************************************************************************************************************
! Copyright c 2013, The Ames Laboratory, Iowa State University, and Rutgers University*.  All rights reserved.
! 
! This software was authored by Yongxin Yao, Nicola Lanata*, Gabriel Kotliar*, Cai-Zhuang Wang, and Kai-Ming Ho, 
! at The Ames Laboratory and Rutgers University and was supported by the U.S. Department of Energy (DOE), Office of Science, 
! Basic Energy Sciences, Materials Science and Engineering Division.  
! The Ames Laboratory is operated by Iowa State University for DOE under U.S. Government contract DE-AC02-07CH11358.  
! The U.S. Government has the rights to use, reproduce, and distribute this software.  
! NEITHER THE GOVERNMENT, THE AMES LABORATORY, IOWA STATE UNIVERSITY, NOR RUTGERS UNIVERSITY MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
! If software is modified to produce derivative works, such modified software should be clearly marked, 
! so as to not confuse it with the version available from The Ames Laboratory and Rutgers University.
! 
! Additionally, redistribution and use in source and binary forms, with or without modification, 
! are permitted provided that the following conditions are met:
! 
!     Redistribution of source code must retain the above copyright notice, this list of conditions, 
! and the following disclaimer.
!     Redistribution in binary form must reproduce the above copyright notice, this list of conditions, 
! and the following disclaimer in the documentation and/or other materials provided with distribution.
!     Neither the name of The Ames Laboratory, Iowa State University, Rutgers University, the U.S. Government, 
! nor the names of its contributors may be used to endorse or promote products derived from this software 
! without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE AMES LABORATORY, IOWA STATE UNIVERSITY, RUTGERS UNIVERSITY, AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
! AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE GOVERNMENT, THE AMES LABORATORY, 
! IOWA STATE UNIVERSITY, RUTGERS UNIVERSITY, OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************************************************************************************************************
      MODULE CORRORB
      USE gprec; USE GREENFUN; USE SPARSE; USE GCONSTANT; USE GUTIL
      IMPLICIT NONE
!
! DEFINE TYPE
      TYPE SYM_MATRIX
        INTEGER,POINTER :: ID(:,:),IDL(:,:)
        INTEGER N,N1,N2,NR,IDSUM,NB
        INTEGER,POINTER :: IDP(:,:,:),DIMP(:)
        COMPLEX(gq),POINTER :: B(:,:,:) ! Matrix Basis
      END TYPE SYM_MATRIX
!
      TYPE CORR_ORB
        INTEGER DIM,DIM2,DIMSO ! Number of correlated orbital / spin-orbital / DIM*ISPIN
        COMPLEX(gq),POINTER :: NKS(:,:),NC_VAR(:,:),NC_PHY(:,:) ! Local density matrix for Kohn-Sham/ {c}/ parameter
        COMPLEX(gq),POINTER :: ISIMIX(:,:)
        REAL(gq)         :: NET,NELF1,NELF2 ! Total number of correlated electrons/ FIXED NEL FOR V_DC
        COMPLEX(gq),POINTER :: R(:,:),R0(:,:),D(:,:),D0(:,:),LA1(:,:),ETA(:,:),LA2(:,:),Z(:,:)
        COMPLEX(gq),POINTER :: RB(:,:,:),DB(:,:,:),LB1(:,:,:),ETB(:,:,:) ! BND%R,BND%LA1
        REAL(gq),POINTER :: N0R(:)
        COMPLEX(gq),POINTER :: LA2R(:) ! Symmetrically inequivalent LA2
        COMPLEX(gq),POINTER :: DR(:)
!
        TYPE(SYM_MATRIX) :: SYMG,SYMH,SYMHO ! General, Hermitian, Hermitian off-diagonal
        TYPE(ISYM_BK_MATRIX) SYMG_BK
        REAL(gq)         :: ELC,ELC_OLD ! Level center of the Correlated orbitals
        COMPLEX(gq),POINTER :: EL0(:,:) ! Local orbital levels in Kohn-Sham/Gutz
        COMPLEX(gq),POINTER :: UK (:,:,:,:) ! <Psi_k|Phi_loc>
        COMPLEX(gq),POINTER :: VK (:,:,:,:,:) ! CO part of bands VK
        COMPLEX(gq),POINTER :: HK0(:,:,:,:) ! <Phi_loc|H_LDA|Phi_loc>
        COMPLEX(gq),POINTER :: C2N(:,:) ! Unitary trans from complex Spherical Harmonics to local natural orbitals
        COMPLEX(gq),POINTER :: R2N(:,:),B2N(:,:),N2N(:,:)
        TYPE(GREEN_FUN)  :: GF
        REAL(gq)         :: VDC,EDC,EDCLA1,EDCUJ,EDCUJV ! Double counting potential
        COMPLEX(gq),POINTER :: VDC2(:,:),NPHY_FIX(:,:) ! CMR Onsite HF
        REAL(gq)         :: U,J,UB,JB,F0,F2,F4,F6 ! COULOMB PARAMETERS/Average value
        REAL(gq)         :: S2,L2,J2,SZ  ! <J^2>
        REAL(gq),POINTER :: V2AO(:,:,:,:) ! V2AO(M1,M3;M2,M4): Coulomb matrix within AO/SPIN-AO basis
        COMPLEX(gq),POINTER :: V2H(:,:,:,:) ! With unitary transformation
      END TYPE CORR_ORB
!
! SUBROUTINE
      CONTAINS
!
!******************************************************************************
      SUBROUTINE GROUP_SYM_ID(SYM)
      TYPE(SYM_MATRIX):: SYM
! LOCAL
      INTEGER I,J,N,N1,N2,IDL
!
      SYM%IDSUM=SUM(SYM%ID)
      N=SYM%N
      CALL GET_IMIN1(SYM%ID,N*N,N1)
      N2=MAXVAL(SYM%ID)
      IF(N1>N2)STOP ' ERROR IN SYM%ID!'
      SYM%N1=N1; SYM%N2=N2; SYM%NR=N2-N1+1
      ALLOCATE(SYM%IDP(2,N,SYM%NR)); SYM%IDP=0
      ALLOCATE(SYM%DIMP(SYM%NR)); SYM%DIMP=0
      ALLOCATE(SYM%IDL(N,N))
      DO I=1,N; DO J=1,N
        IF(SYM%ID(I,J)==0)THEN
          SYM%IDL(I,J)=0
        ELSE
          SYM%IDL(I,J)=SYM%ID(I,J)-N1+1
        ENDIF
      ENDDO; ENDDO
      DO I=1,N; DO J=1,N
        IDL=SYM%IDL(I,J)
        IF(IDL==0)CYCLE
        SYM%DIMP(IDL)=SYM%DIMP(IDL)+1
        SYM%IDP(1,SYM%DIMP(IDL),IDL)=I
        SYM%IDP(2,SYM%DIMP(IDL),IDL)=J
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE GROUP_SYM_ID
!
!******************************************************************************
! 1: SYMH -> Hermitin matrix basis
!******************************************************************************
      SUBROUTINE SYMH_TO_MATRIX_BASIS(SYM)
      TYPE(SYM_MATRIX):: SYM
! LOCAL
      INTEGER I,J,I1,J1,ISUM,IADD
      REAL(gq) RES
!
      SYM%NB=0
      DO I=1,SYM%NR
      IF(SYM%IDP(1,1,I)==SYM%IDP(2,1,I))THEN
        SYM%NB=SYM%NB+1
      ELSE
        SYM%NB=SYM%NB+1 *2
      ENDIF
      ENDDO ! I
      ALLOCATE(SYM%B(SYM%N,SYM%N,SYM%NB)); SYM%B=0
      ISUM=1
      DO I=1,SYM%NR
        RES=REAL(SYM%DIMP(I),gq)
        IADD=1
        IF(SYM%IDP(1,1,I)/=SYM%IDP(2,1,I))THEN
          RES=RES*2
          IADD=1 *2
        ENDIF
        RES=D1/SQRT(RES)
        DO J=1,SYM%DIMP(I)
        I1=SYM%IDP(1,J,I); J1=SYM%IDP(2,J,I)
        IF(I1==J1)THEN
          SYM%B(I1,J1,ISUM)=RES
        ELSE
          SYM%B(I1,J1,ISUM)=RES
          SYM%B(J1,I1,ISUM)=SYM%B(I1,J1,ISUM)
          SYM%B(I1,J1,ISUM+1)=ZI*RES
          SYM%B(J1,I1,ISUM+1)=CONJG(SYM%B(I1,J1,ISUM+1))
        ENDIF
        ENDDO ! J
        ISUM=ISUM+IADD
      ENDDO ! I
      RETURN
!
      END SUBROUTINE SYMH_TO_MATRIX_BASIS
!
!******************************************************************************
      SUBROUTINE OUT_CO_SYM(SYM,IO,NI,STR)
      TYPE(SYM_MATRIX):: SYM
      INTEGER IO,NI
      CHARACTER(*) STR
! LOCAL
      INTEGER I
!
      WRITE(IO,'(" NI=",I3,"  ",A5,"%NR=",I2)')NI,STR,SYM%NR
      DO I=1,SYM%NR
        WRITE(IO,'(20(I2,",",I2," | "))')SYM%IDP(:,1:SYM%DIMP(I),I)
      ENDDO
      RETURN
      END SUBROUTINE OUT_CO_SYM
!
!******************************************************************************
! MODE = 1: keep diagonal elements and lower part
!        0: keep lower part
!******************************************************************************
      SUBROUTINE REDUCE_SYM(SYM1,SYM2,MODE)
      TYPE(SYM_MATRIX):: SYM1,SYM2
      INTEGER MODE
! LOCAL
      INTEGER N,I,J,ID,KD,K,L,ISUM
      INTEGER ID_LIST(SYM1%N*SYM1%N)
!
      N=SYM1%N; SYM2%N=N; ID_LIST=0; ISUM=0
      ALLOCATE(SYM2%ID(N,N)); SYM2%ID=SYM1%ID
      DO I=1,N; DO J=I+MODE,N
        ID=SYM2%ID(I,J)
        IF(ID==0)CYCLE
        DO K=1,ISUM
          IF(ID==ID_LIST(K))GOTO 100
        ENDDO
        ISUM=ISUM+1
        ID_LIST(ISUM)=ID
        DO K=1,N; DO L=1,K+MODE-1
          KD=SYM1%ID(K,L)
          IF(KD<ID)CYCLE
          SYM2%ID(K,L)=SYM2%ID(K,L)-1
        ENDDO; ENDDO
100     CONTINUE
        SYM2%ID(I,J)=0
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE REDUCE_SYM
!
!******************************************************************************
      SUBROUTINE INI_SYM_LOC_ARRAY(CO)
      TYPE (CORR_ORB) CO
!
      ALLOCATE(CO%LA2R(CO%SYMH%NR),CO%DR(CO%SYMG%NR),CO%N0R(CO%SYMG%NR))
      RETURN
!
      END SUBROUTINE INI_SYM_LOC_ARRAY
!
!****************************************************************************
! MODE = 0: -->SYMG
!        1: -->SYMH
!****************************************************************************
      SUBROUTINE SET_SYM_LOC_ARRAY(A,AR,CO,LBACK,MODE)
      TYPE (CORR_ORB) CO
      COMPLEX(gq) A(CO%DIM2,CO%DIM2),AR(*)
      LOGICAL LBACK
      INTEGER MODE
! LOCAL
      INTEGER N2,NR,I,J,IDL
      INTEGER,POINTER::ID(:,:)
!
      N2=CO%DIM2
      IF(LBACK)THEN
        SELECT CASE(MODE)
        CASE(0); ID=>CO%SYMG%IDL
        CASE(1); ID=>CO%SYMH%IDL
        END SELECT
        A=0
        DO I=1,N2; DO J=1,N2
          IDL=ID(I,J)
          IF(IDL==0)CYCLE
          A(I,J)=AR(IDL)
          IF(MODE==0.OR.I==J)CYCLE
          A(J,I)=CONJG(AR(IDL))
        ENDDO; ENDDO
      ELSE
        SELECT CASE(MODE)
        CASE(0); ID=>CO%SYMG%IDP(:,1,:); NR=CO%SYMG%NR
        CASE(1); ID=>CO%SYMH%IDP(:,1,:); NR=CO%SYMH%NR
        END SELECT
        DO I=1,NR
          AR(I)=A(ID(1,I),ID(2,I))
        ENDDO
      ENDIF
      NULLIFY(ID)
      RETURN
!
      END SUBROUTINE SET_SYM_LOC_ARRAY
!
!******************************************************************************
      SUBROUTINE SET_SYM_LOC_N0R(CO)
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER I,IA,J
!
      CO%N0R=0
      DO I=1,CO%SYMG%NR
        IA=CO%SYMG%IDP(1,1,I)
        CO%N0R(I)=CO%NKS(IA,IA)
        DO J=2,CO%SYMG%DIMP(I)
          IA=CO%SYMG%IDP(1,J,I)
          IF(ABS(CO%N0R(I)-CO%NKS(IA,IA))>1.E-8_gq)THEN
            STOP ' ERROR IN SET_SYM_LOC_N0R!'
          ENDIF
        ENDDO
      ENDDO
      RETURN
!
      END SUBROUTINE SET_SYM_LOC_N0R
!
!******************************************************************************
      SUBROUTINE SET_U_MAT(CO,LHUB,NT,IU)
      TYPE (CORR_ORB) CO
      INTEGER LHUB,NT,IU
!
      SELECT CASE(LHUB)
      CASE(0)
        CALL READ_V2AO(CO,NT,IU)
      CASE(1)
        CALL SET_V2AO_AI(CO)
      CASE(2)
        CALL SET_V2AO_KM(CO) ! REAL HARMONICS, Kanamori
      CASE(3)
        CALL SET_V2AO_AII(CO)
      END SELECT
      CALL CALC_MEAN_UJ(CO)
      CALL SET_U2_MAT(CO)
      RETURN
!
      END SUBROUTINE SET_U_MAT
!
!******************************************************************************
      SUBROUTINE SET_U2_MAT(CO)
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER NA,NA2
!
      NA=CO%DIM
      NA2=NA*2
      ALLOCATE(CO%V2H(NA2,NA2,NA2,NA2)); CO%V2H=0 !V2H(a_isp1,b_isp2,c_isp1,d_isp2)
      CO%V2H(1:NA ,1:NA ,1:NA ,1:NA )=CO%V2AO
      CO%V2H(1+NA:,1+NA:,1+NA:,1+NA:)=CO%V2AO
      CO%V2H(1:NA ,1+NA:,1:NA ,1+NA:)=CO%V2AO
      CO%V2H(1+NA:,1:NA ,1+NA:,1:NA )=CO%V2AO
      RETURN
!
      END SUBROUTINE SET_U2_MAT
!
!*************************************************************************************
      SUBROUTINE CALC_CO_EL0(CO)
      USE BANDSTRU
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER ISYMI,ISYMF
      INTEGER IVEC,IKS,IKP,IKPL,NKP,ISYM,NSYM,I1,NASO,NA2
      REAL(gq) WTK
!
      ISYMI=SYM%IDI; ISYMF=SYM%IDF
      NSYM=ISYMF-ISYMI+1
      NASO=CO%DIMSO; NA2=CO%DIM2
      CO%EL0=0; CO%EL0=0
!
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      WTK=KPT%WT(IKP)/NSYM
      DO ISYM=1,NSYM
      CO%EL0(1:NASO,1:NASO)=CO%EL0(1:NASO,1:NASO)+CO%HK0(:,:,ISYM,IKPL)*WTK
      ENDDO ! ISYM
      ENDDO; ENDDO
!
      IF(BND%ISO.EQ.1)CO%EL0(1+NASO:,1+NASO:)=CO%EL0(1:NASO,1:NASO)
      CALL ZSPIN_BLK2_TRANS(CO%EL0,NA2,.TRUE.,BND%ISO) ! Orbital-fast to Spin-fast
      RETURN
!
      END SUBROUTINE CALC_CO_EL0
!
!*************************************************************************************
      SUBROUTINE ROT_CO_RLA1(CO)
      TYPE (CORR_ORB) CO
! LOCAL
      COMPLEX(gq) U(CO%DIM2,CO%DIM2)
      REAL(gq) W(CO%DIM2)
!
      U=-CO%NKS
      CALL ZSYM_BK_DIAG(U,W,CO%DIM2,CO%SYMG_BK)
!
      CALL ZANNXBNM('C',U,CO%R,CO%DIM2,CO%DIM2) ! U^+ R
      CALL ZUHAU(CO%LA1,U,CO%DIM2,CO%DIM2) ! U^+ \lambda U
      U=CONJG(U)
      CALL ZUHAU(CO%NKS,U,CO%DIM2,CO%DIM2)
      RETURN
!
      END SUBROUTINE ROT_CO_RLA1
!
!*************************************************************************************
      SUBROUTINE CALC_CO_EL0_UV(CO)
      USE BANDSTRU
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER ISYMI,ISYMF
      INTEGER IVEC,IKS,IKP,IKPL,NKP,ISYM,NSYM,NEMIN,NEMAX,NBANDS,I1,NAS,NA2
      REAL(gq) WTK
      CHARACTER*1 MATDESCRA(6)
      COMPLEX(gq) CM1(CO%DIMSO,CO%DIMSO)
      REAL(gq),POINTER::EK0(:)
      COMPLEX(gq),POINTER::UK(:,:)
      COMPLEX(gq),ALLOCATABLE::HUK(:,:)
!
      ISYMI=SYM%IDI; ISYMF=SYM%IDF
      NSYM=ISYMF-ISYMI+1
      NAS=CO%DIMSO; NA2=CO%DIM2
      ALLOCATE(HUK(BND%NMAXIN,NAS))
!
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      WTK=KPT%WT(IKP)/NSYM
      NEMIN=BND%NE(2,IKP); NEMAX=BND%NE(3,IKP); NBANDS=NEMAX-NEMIN+1
      EK0=>BND%EK0(NEMIN:NEMAX,IKP)
      DO ISYM=ISYMI,ISYMF
      UK=>CO%UK(1:NBANDS,:,ISYM,IKPL)
      DO I1=1,NBANDS; HUK(I1,1:NAS)=UK(I1,1:NAS)*EK0(I1); ENDDO
      CALL ZGEMM('C','N',NAS,NAS,NBANDS,Z1,UK,NBANDS,HUK(1:NBANDS,1:NAS),NBANDS,Z0,CM1,NAS)
      CO%EL0(1:NAS,1:NAS)=CO%EL0(1:NAS,1:NAS)+REAL(CM1,KIND=gq)*WTK
      ENDDO ! ISYM
      ENDDO; ENDDO
!
      DEALLOCATE(HUK)
      IF(BND%ISPO.EQ.1)CO%EL0(1+NAS:,1+NAS:)=CO%EL0(1:NAS,1:NAS)
      CALL ZSPIN_BLK2_TRANS(CO%EL0,NA2,.TRUE.,BND%ISO) ! Orbital-fast to Spin-fast
      RETURN
!
      END SUBROUTINE CALC_CO_EL0_UV
!
!****************************************************************************
      SUBROUTINE CALC_CO_NKS(CO,MODE)
      USE BANDSTRU
      TYPE (CORR_ORB) CO
      INTEGER MODE
! LOCAL
      INTEGER I1,NASO,NA2
!
      NASO=CO%DIMSO; NA2=CO%DIM2
      IF(MODE==0)THEN
        CALL CALC_CO_NKS_BND(CO)
      ELSEIF(CO%GF%WTYP==0)THEN
        IF(MODE==1)THEN
          CALL GREEN_OCC(1,CO%GF,OCC=CO%NKS(1:NASO*BND%NSPIN,1:NASO*BND%NSPIN))
        ELSE
          CALL GREEN_OCC2(CO%GF,CO%NKS(1:NASO*BND%NSPIN,1:NASO*BND%NSPIN),CO%RB)
        ENDIF
      ELSEIF(CO%GF%WTYP==1)THEN
        CALL GREEN_OCC_R1(1,1,CO%GF,CO%RB,OCC=CO%NKS(1:NASO*BND%NSPIN,1:NASO*BND%NSPIN))
      ENDIF
!
      IF(BND%ISPO.EQ.1)CO%NKS(1+NASO:,1+NASO:)=CO%NKS(1:NASO,1:NASO)
      CALL ZSPIN_BLK2_TRANS(CO%NKS,NA2,.TRUE.,BND%ISO)
      DO I1=1,NA2
        IF(REAL(CO%NKS(I1,I1),gq)<-1.E-2_gq)THEN
          WRITE(0,'(" ERROR in CALC_CO_NKS: Large negative NKS=",2F10.3)')CO%NKS(I1,I1)
          STOP
        ENDIF
        CO%NKS(I1,I1)=MAX(ABS(CO%NKS(I1,I1)),SMALL)
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_CO_NKS
!
!****************************************************************************
      SUBROUTINE CALC_CO_NKS_BND(CO)
      USE BANDSTRU
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER ISYMI,ISYMF,IB,ISYM,ISP,N1,N2
      INTEGER IVEC,IKS,IKP,IKPL,NKP,NSYM,NEMIN,NEMAX,NBANDS,NASO,NA2,MAX_IN_BANDS
      REAL(gq) WTK
      REAL(gq),POINTER        :: FERWE(:)
      COMPLEX(gq),POINTER     :: VK(:,:)
      COMPLEX(gq),ALLOCATABLE :: VF(:,:),NABK(:,:)
!
      ISYMI=SYM%IDI; ISYMF=SYM%IDF
      NSYM=ISYMF-ISYMI+1
      NASO=CO%DIMSO; NA2=CO%DIM2
      MAX_IN_BANDS=BND%NMAXIN
      ALLOCATE(VF(NASO,MAX_IN_BANDS),NABK(NASO,NASO))
!
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      NEMIN=BND%NE(2,IKP); NEMAX=BND%NE(3,IKP)
      NBANDS=NEMAX-NEMIN+1
      WTK=1._gq/NSYM/BND%RSPO
      DO ISP=1,BND%NSPIN
      FERWE=>BND%FERWE(NEMIN:NEMAX,IKP,ISP)
      N1=(ISP-1)*NASO+1; N2=ISP*NASO
      DO ISYM=1,NSYM
      VK=>CO%VK(:,1:NBANDS,ISYM,IKPL,ISP)
      DO IB=1,NBANDS; VF(:,IB)=VK(:,IB)*FERWE(IB)*WTK; ENDDO
      CALL ZGEMM('N','C',NASO,NASO,NBANDS,Z1,VF,NASO,VK,NASO,Z0,NABK,NASO)
      CO%NKS(N1:N2,N1:N2)=CO%NKS(N1:N2,N1:N2)+NABK
      ENDDO; ENDDO ! ISYM,ISP
      ENDDO; ENDDO
!
      NULLIFY(VK,FERWE); DEALLOCATE(VF,NABK)
      RETURN
!
      END SUBROUTINE CALC_CO_NKS_BND
!
!******************************************************************************
      SUBROUTINE WRITE_V2AO(CO,NT,IU)
      TYPE (CORR_ORB) CO
      INTEGER NT,IU
! LOCAL
      INTEGER I1,I2,I3,I4,NA
!
      NA=CO%DIM
      WRITE(IU,'("NT=",I3)')NT
      DO I1=1,NA; DO I2=1,NA; DO I3=1,NA; DO I4=1,I2
        IF(ABS(CO%V2AO(I1,I2,I3,I4))<1.E-10_gq)CYCLE
        WRITE(IU,'(4I3,2X,F16.10)')I1,I2,I3,I4,CO%V2AO(I1,I2,I3,I4)
      ENDDO; ENDDO; ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE WRITE_V2AO
!
!******************************************************************************
      SUBROUTINE READ_V2AO(CO,NT,IU)
      TYPE (CORR_ORB) CO
      INTEGER NT,IU
! LOCAL
      INTEGER NA,IT,I1,I2,I3,I4
      REAL(gq) RES
      CHARACTER*3 STR1,STR2
!
      NA=CO%DIM
      ALLOCATE(CO%V2AO(NA,NA,NA,NA)); CO%V2AO=0
      OPEN(IU,FILE='V2AO.INP',STATUS='OLD')
      DO 
        READ(IU,'(2A3)')STR1,STR2
        IF(STR1.NE.'NT=')CYCLE
        READ(STR2,*)IT
        IF(IT.NE.NT)CYCLE
        DO 
          READ(IU,*,END=100,ERR=100)I1,I2,I3,I4,RES
          CO%V2AO(I1,I2,I3,I4)=RES
          CO%V2AO(I1,I4,I3,I2)=RES
        ENDDO
      ENDDO
      WRITE(0,'(" NT=",I3)')NT
      STOP ' ERROR: CANNOT LOCATE V2AO!'
100   CONTINUE
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE READ_V2AO
!
!******************************************************************************
      SUBROUTINE SET_V2AO_KM(CO)
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER M1,M2,NA
      REAL(gq) U,UP,J
!
      NA=CO%DIM
      ALLOCATE(CO%V2AO(NA,NA,NA,NA)); CO%V2AO=0
!
      U=CO%U; J=CO%J; UP=U-2*J
      DO M1=1,NA; DO M2=1,NA
      IF(M1.EQ.M2)THEN
        CO%V2AO(M1,M2,M1,M2)=U
      ELSE
        CO%V2AO(M1,M2,M1,M2)=UP      ! different orb, arbitrary spin
        CO%V2AO(M1,M2,M2,M1)=J       ! different orb, same spin
        CO%V2AO(M1,M1,M2,M2)=J       ! different orb, same spin
      ENDIF
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE SET_V2AO_KM
!
!******************************************************************************
      SUBROUTINE SET_V2AO_AI(CO)
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER L,NA
      REAL(gq) F(0:6),RAT42,RAT62
!
      NA=CO%DIM
      ALLOCATE(CO%V2AO(NA,NA,NA,NA)); CO%V2AO=0
      L=(NA-1)/2
      F=0
! Set F using U and J
      SELECT CASE(L)
      CASE(1)
         F(2)=5._gq*CO%J
      CASE(2)
         RAT42=0.625_gq
         F(2)=14._gq/(1+RAT42)*CO%J
         F(4)=RAT42*F(2)
      CASE(3)
         RAT42=0.668_gq
         RAT62=0.494_gq
         F(2)=6435._gq/(286._gq+195._gq*RAT42+250._gq*RAT62)*CO%J
         F(4)=RAT42*F(2)
         F(6)=RAT62*F(2)
      END SELECT
      F(0)=CO%U
      CALL COULOMB_INT(F,CO)
      RETURN
!
      END SUBROUTINE SET_V2AO_AI
!
!******************************************************************************
      SUBROUTINE SET_V2AO_AII(CO)
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER NA
      REAL(gq) F(0:6)
!
      NA=CO%DIM
      ALLOCATE(CO%V2AO(NA,NA,NA,NA)); CO%V2AO=0
      F=0; F(0)=CO%F0; F(2)=CO%F2; F(4)=CO%F4; F(6)=CO%F6
      CALL COULOMB_INT(F,CO)
      RETURN
!
      END SUBROUTINE SET_V2AO_AII
!
!***********************************************************************
! SET COULOMB INTERACTION BETWEEN LOCAL ORBITALS USING F
!   <lm1|Y_LM|lm2><lm3|Y_LM*|lm4> 
! = <lm1|Y_LM|lm2><lm4|Y_LM|lm3>
! = Gaunt(-m1,m2,L,M)*Gaunt(-m4,m3,L,M)*(-1)^(m1+m4)
!***********************************************************************
      SUBROUTINE COULOMB_INT(F,CO)
      USE gasa
      TYPE (CORR_ORB) CO
      REAL(gq) F(0:6)
! LOCAL
      INTEGER L0,M1,M2,M3,M4,LMIN,LMAX
      INTEGER LL,MM
      REAL(gq) RES
      REAL(gq),ALLOCATABLE::GAUNT(:,:,:,:),U(:)
!
      L0=(CO%DIM-1)/2; LMIN=0; LMAX=2*L0
      ALLOCATE(GAUNT(-L0:L0,-L0:L0,LMIN:LMAX,-LMAX:LMAX),U(LMIN:LMAX))
      CALL YLM3ST_COMPL(L0,L0,GAUNT)
!
      DO M1=-L0,L0; DO M3=-L0,L0; DO M2=-L0,L0; DO M4=-L0,L0
      U=0
      DO LL=LMIN,LMAX,2; DO MM=-LL,LL
        RES=GAUNT(-M1,M2,LL,MM)*FS(M1)*GAUNT(-M4,M3,LL,MM)*FS(M4)
        IF(ABS(RES).LT.1.E-16_gq) CYCLE
        U(LL)=U(LL)+4*PI/(2*LL+1)*RES*F(LL)
      ENDDO; ENDDO
! Order <1,2>; <3,4>
      CO%V2AO(M1+L0+1,M3+L0+1,M2+L0+1,M4+L0+1)=SUM(U)
      ENDDO; ENDDO; ENDDO; ENDDO
      DEALLOCATE(GAUNT,U)
      RETURN
!
      END SUBROUTINE COULOMB_INT
!
!******************************************************************************
      SUBROUTINE CALC_MEAN_UJ(CO)
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER NA,M1,M2,ISU,ISJ
      REAL(gq) UB,JB
!
      NA=CO%DIM
      UB=0; JB=0; ISU=0; ISJ=0
      DO M1=1,NA; DO M2=1,NA
        UB=UB+CO%V2AO(M1,M2,M1,M2)
        ISU=ISU+1
        IF(M1.NE.M2)THEN
          JB=JB+CO%V2AO(M1,M2,M1,M2)-CO%V2AO(M1,M2,M2,M1); ISJ=ISJ+1
        ENDIF
      ENDDO; ENDDO
      UB=UB/ISU
      IF(ISJ.EQ.0)THEN
        JB=0
      ELSE
        JB=UB-JB/ISJ
      ENDIF
      CO%UB=UB; CO%JB=JB
      RETURN
!
      END SUBROUTINE CALC_MEAN_UJ
!
!***********************************************************************************************
! Need to replace N0 by NC_PHY for non-diagonal version
!***********************************************************************************************
      SUBROUTINE CALC_CO_VDC(CO,LDC)
      INTEGER LDC
      TYPE(CORR_ORB) CO
! LOCAL
      INTEGER I,J
      REAL(gq) UB,JB,NET
!
      CO%VDC=0
      IF(LDC.GT.0)THEN
        UB=CO%UB; JB=CO%JB
        SELECT CASE(LDC)
        CASE(1)
          NET=CO%NET
        CASE(2,3,4,12)
          NET=CO%NELF1
        END SELECT
        CO%VDC=UB*(NET-.5_gq)-JB/2*(NET-1)
      ELSEIF(LDC.EQ.-1.OR.LDC.EQ.-12)THEN
        CO%VDC2=0
        DO I=1,CO%DIM2; DO J=1,CO%DIM2
          CO%VDC2(I,J)=CO%VDC2(I,J)+SUM(CO%NPHY_FIX*(CO%V2H(I,:,J,:)-CO%V2H(I,:,:,J)))
        ENDDO; ENDDO
      ENDIF
      RETURN
!
      END SUBROUTINE CALC_CO_VDC
!
!****************************************************************************
      SUBROUTINE CALC_CO_LA(CO,LDC,MODE)
      INTEGER LDC,MODE
      TYPE(CORR_ORB) CO
! LOCAL
      INTEGER I,J,K
      COMPLEX(gq) RD(CO%DIM2,CO%DIM2)
      COMPLEX(gq) PN(CO%DIM2,CO%DIM2),H(CO%DIM2,CO%DIM2)
      COMPLEX(gq),POINTER::LAY(:,:),LAX(:,:)
      REAL(gq) RES
!
      IF(MODE==2)THEN
        LAY=>CO%LA2; LAX=>CO%LA1
      ELSE
        LAY=>CO%LA1; LAX=>CO%LA2
      ENDIF
      RD=MATMUL(CO%R,TRANSPOSE(CO%D))
      LAY=0
      DO I=1,CO%SYMH%NB
        H=CO%SYMH%B(:,:,I)
        CALL ZPFA_PA(CO%NKS,PN,H,CO%DIM2,DSIMIX,DPSIMIX) ! p f / p d_n
        RES=-SUM(PN*RD)*2
        LAY=LAY+CONJG(CO%SYMH%B(:,:,I))*RES
      ENDDO
      LAY=LAY-LAX
      DO I=1,CO%DIM2
        LAY(I,I)=LAY(I,I)-CO%VDC
      ENDDO
      IF(LDC<0)THEN
        LAY=LAY-CO%VDC2
      ENDIF
      NULLIFY(LAX,LAY)
      IF(MODE==2)THEN
        CALL SET_SYM_LOC_ARRAY(CO%LA2,CO%LA2R,CO,.FALSE.,1)
      ENDIF
      RETURN
!
      END SUBROUTINE CALC_CO_LA
!
!****************************************************************************
! ISIMIX = (N0(1-N0))^(-1/2)
!****************************************************************************
      SUBROUTINE CALC_CO_ISIMIX(CO,NKS)
      TYPE(CORR_ORB) CO
      COMPLEX(gq) NKS(CO%DIM2,CO%DIM2)
! LOCAL
      COMPLEX(gq) XN(CO%DIM2,CO%DIM2)
! 
      XN=NKS; XN=XN-MATMUL(XN,XN)
      CALL ZATOFA(XN,CO%ISIMIX,CO%DIM2,-12,D1,.TRUE.)
      RETURN
!
      END SUBROUTINE CALC_CO_ISIMIX
!
!****************************************************************************
! D = (N0(1-N0))^(-1/2) D0
!****************************************************************************
      SUBROUTINE CO_D0_TO_D(CO)
      TYPE(CORR_ORB) CO
!
      CO%D=MATMUL(CO%ISIMIX,CO%D0)
      RETURN
!
      END SUBROUTINE CO_D0_TO_D
!
!****************************************************************************
      SUBROUTINE CALC_CO_EDCLA1(CO)
      TYPE(CORR_ORB) CO
! LOCAL
      INTEGER I
!
      CO%EDCLA1=-SUM((CO%LA1+CO%ETA)*CO%NKS)
      RETURN
!
      END SUBROUTINE CALC_CO_EDCLA1
!
!****************************************************************************
      SUBROUTINE CALC_CO_EDCUJ(CO,LDC)
      INTEGER LDC
      TYPE(CORR_ORB) CO
! LOCAL
      INTEGER I,J
      REAL(gq) UB,JB,NET
!
      CO%EDCUJ=0
      IF(LDC.GT.0)THEN
        UB=CO%UB; JB=CO%JB
        NET=CO%NET
        SELECT CASE(LDC)
        CASE(1) ! Standard
          CO%EDCUJ=-(UB*NET*(NET-1)/2-JB*NET/2*(NET/2-1))
        CASE(2,12) ! Fix Nf, Nf(V_DC)
          CO%EDCUJV=-CO%VDC*NET ! Haule, change E_DC
          NET=CO%NELF1
          CO%EDCUJ =-(UB*NET*(NET-1)/2-JB*NET/2*(NET/2-1)) &
                  &-CO%VDC*(CO%NET-NET)
        CASE(3) ! Fix Nf, independent of V_DC
          NET=CO%NELF2
          CO%EDCUJ=-(UB*NET*(NET-1)/2-JB*NET/2*(NET/2-1)) &
                  &-CO%VDC*(CO%NET-NET)
        END SELECT
      ELSEIF(LDC==-1.OR.LDC==-12)THEN
        CO%EDCUJ=0
        DO I=1,CO%DIM2; DO J=I+1,CO%DIM2
          CO%EDCUJ=CO%EDCUJ-CO%NPHY_FIX(I,J)*SUM(CO%NPHY_FIX*(CO%V2H(I,:,J,:)-CO%V2H(I,:,:,J)))
        ENDDO; ENDDO
        IF(LDC==-12)THEN
          CO%EDCUJ=CO%EDCUJ-SUM((CO%NKS-CO%NPHY_FIX)*CO%VDC2)
        ENDIF
      ENDIF
      RETURN
!
      END SUBROUTINE CALC_CO_EDCUJ
!
      END MODULE CORRORB

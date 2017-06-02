!******************************************************************************
! Copyright c 2013, The Ames Laboratory, Iowa State University, and Rutgers
! University*.  All rights reserved.
! 
! This software was authored by Yongxin Yao, Nicola Lanata*, Gabriel Kotliar*,
! Cai-Zhuang Wang, and Kai-Ming Ho, at The Ames Laboratory and 
! Rutgers University and was supported by the U.S. 
! Department of Energy (DOE), Office of Science, 
! Basic Energy Sciences, Materials Science and Engineering Division.  
! The Ames Laboratory is operated by Iowa State University for DOE 
! under U.S. Government contract DE-AC02-07CH11358.  
! The U.S. Government has the rights to use, reproduce, and 
! distribute this software.  
! NEITHER THE GOVERNMENT, THE AMES LABORATORY, IOWA STATE UNIVERSITY, 
! NOR RUTGERS UNIVERSITY MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
! If software is modified to produce derivative works, 
! such modified software should be clearly marked, 
! so as to not confuse it with the version available from 
! The Ames Laboratory and Rutgers University.
! 
! Additionally, redistribution and use in source and binary forms, 
! with or without modification, 
! are permitted provided that the following conditions are met:
! 
!     Redistribution of source code must retain the above copyright notice,
!     this list of conditions, and the following disclaimer.
!
!     Redistribution in binary form must reproduce the above copyright notice,
!     this list of conditions, and the following disclaimer 
!     in the documentation and/or other materials provided with distribution.
!
!     Neither the name of The Ames Laboratory, Iowa State University, 
!     Rutgers University, the U.S. Government, nor the names of 
!     its contributors may be used to endorse or promote products derived 
!     from this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE AMES LABORATORY, IOWA STATE UNIVERSITY, 
! RUTGERS UNIVERSITY, AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
! THE IMPLIED WARRANTIES OF MERCHANTABILITY 
! AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  
! IN NO EVENT SHALL THE GOVERNMENT, THE AMES LABORATORY, 
! IOWA STATE UNIVERSITY, RUTGERS UNIVERSITY, OR CONTRIBUTORS BE LIABLE 
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
! OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!******************************************************************************
!
!
      MODULE CORRORB
      USE gprec; USE GREENFUN; USE GREENFUN2; USE SPARSE
      USE GCONSTANT; USE GUTIL
      USE GINFO
      USE GHDF5_BASE, ONLY: gh5_write, embedH_file_id, gh5_open, gh5_close
      IMPLICIT NONE
!
      TYPE CORR_ORB
        ! Number of correlated orbital / spin-orbital / DIM*ISPIN
        INTEGER DIM,DIM2,DIMSO 
        ! Local density matrix for Kohn-Sham/ {c}/ parameter
        COMPLEX(gq),POINTER :: NKS(:,:),NC_VAR(:,:),NC_PHY(:,:) 
        COMPLEX(gq),POINTER :: ISIMIX(:,:)
        ! Total number of correlated electrons/ FIXED NEL FOR V_DC
        REAL(gq)         :: NET,NELF1,NELF2 
        COMPLEX(gq),POINTER :: R(:,:),R0(:,:),D(:,:),D0(:,:), &
                                 &LA1(:,:),LA2(:,:),Z(:,:)
        COMPLEX(gq),POINTER :: S_VEC(:,:,:),L_VEC(:,:,:)
! CMR3 BEGIN
        COMPLEX(gq),POINTER :: R1(:,:),PFPR(:,:)
! CMR3 END
        COMPLEX(gq),POINTER :: RB(:,:,:),DB(:,:,:),LB1(:,:,:),ETB(:,:,:) ! BND%R,BND%LA1
        REAL(gq)      ,POINTER :: NKS_COEF(:),LA1_COEF(:),LA2_COEF(:),NCV_COEF(:) ! Expansion coefficients wrt the Hermitian matrix basis 
        COMPLEX(gq),POINTER :: R_COEF(:)
        INTEGER,POINTER :: M_INDEX(:,:,:),M_INDEX_L(:,:,:),M_INDEX_R(:,:,:)
        INTEGER :: DIM_HS,DIM_HS_L,DIM_HS_R,DIM_HS_E
        COMPLEX(gq),POINTER :: HS(:,:,:),HS_L(:,:,:),HS_R(:,:,:), &
                                 &HS_E(:,:,:)
        REAL(gq) :: ELC,ELC_OLD ! Level center of the Correlated orbitals
        COMPLEX(gq),POINTER :: EL0(:,:) ! Local orbital levels in Kohn-Sham/Gutz
        COMPLEX(gq),POINTER :: UK (:,:,:,:) ! <Psi_k|Phi_loc>
        COMPLEX(gq),POINTER :: VK (:,:,:,:,:) ! CO part of bands VK
        COMPLEX(gq),POINTER :: HK0(:,:,:,:) ! <Phi_loc|H_LDA|Phi_loc>
        ! Unitary trans from complex Spherical Harmonics to 
        ! local natural orbitals
        COMPLEX(gq),POINTER :: C2N(:,:) 
        COMPLEX(gq),POINTER :: R2N(:,:),B2N(:,:),N2N(:,:)
        TYPE(GREEN_FUN)  :: GF
        TYPE(GREEN_FUN2) :: GF2
        REAL(gq)         :: VDC,EDC,EDCLA1,EDCUJ,EDCUJV ! Double counting potential
        COMPLEX(gq),POINTER :: VDC2(:,:),NPHY_FIX(:,:) ! CMR Onsite HF
        REAL(gq)         :: U,J,UB,JB,F0,F2,F4,F6 ! COULOMB PARAMETERS/Average value
        REAL(gq)         :: S2,L2,J2 ! two-body
        REAL(gq)         :: S_VAL(3,2),L_VAL(3,2)  ! one-body
        REAL(gq),POINTER :: V2AO(:,:,:,:) ! V2AO(M1,M3;M2,M4): Coulomb matrix within AO/SPIN-AO basis
        COMPLEX(gq),POINTER :: V2H(:,:,:,:) ! With unitary transformation
      END TYPE CORR_ORB
!
! SUBROUTINE
      CONTAINS
!
!******************************************************************************
      SUBROUTINE SET_U_MAT(CO,LHUB,NT,IU)
      TYPE (CORR_ORB) CO
      INTEGER LHUB,NT,IU
!
      SELECT CASE(LHUB)
      CASE(-1)
        CALL READ_V2H(CO,NT,IU)
        RETURN
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
!******************************************************************************
      SUBROUTINE CALC_CO_EL0(CO)
      USE BANDSTRU
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER ISYMI,ISYMF
      INTEGER IVEC,IKP,IKPL,IKS,NKP,ISYM,NSYM,NASO,NA2
      REAL(gq) WTK
!
      ISYMI=SYM%IDI; ISYMF=SYM%IDF
      NSYM=ISYMF-ISYMI+1
      NASO=CO%DIMSO; NA2=CO%DIM2
      CO%EL0=0; CO%EL0=0; IKPL=0
!
      DO IVEC=1,GP%NVEC;   IF(GP%LKPVEC)THEN;     NKP=GP%KVEC(IVEC,2);   ELSE;     NKP=KPT%DIM;   ENDIF;   DO IKS=1,NKP;   IF(GP%LKPVEC)THEN;     IKP=GP%KVEC(IVEC,3)+IKS;     IF(GP%LOMP)THEN;       IKPL=IKP;     ELSE;       IKPL=IKPL+1;    ENDIF;   ELSE;     IKP=IKS;     IF(GP%LOMP)THEN;       IKPL=IKP;     ELSE;       IKPL=IKP-GP%MYRANK*KPT%DIML;     ENDIF;   ENDIF;   IF(IKPL.LE.0) CYCLE;   IF(IKPL.GT.KPT%DIML) EXIT
      WTK=KPT%WT(IKP)/NSYM
      DO ISYM=1,NSYM
      CO%EL0(1:NASO,1:NASO)=CO%EL0(1:NASO,1:NASO)+CO%HK0(:,:,ISYM,IKPL)*WTK
      ENDDO ! ISYM
      ENDDO; ENDDO
!
      IF(BND%ISO.EQ.1)CO%EL0(1+NASO:,1+NASO:)=CO%EL0(1:NASO,1:NASO)
      CALL ORBITAL_SPIN_TRANS(CO%EL0,NA2,.TRUE.,BND%ISO) ! Orbital-fast to Spin-fast
      RETURN
!
      END SUBROUTINE CALC_CO_EL0
!
!*************************************************************************************
      SUBROUTINE CALC_CO_EL0_UV(CO)
      USE BANDSTRU
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER ISYMI,ISYMF
      INTEGER IVEC,IKP,IKPL,IKS,NKP,ISYM,NSYM,NEMIN,NEMAX,NBANDS,I1,NAS,NA2
      REAL(gq) WTK
      COMPLEX(gq) CM1(CO%DIMSO,CO%DIMSO)
      REAL(gq),POINTER::EK0(:)
      COMPLEX(gq),POINTER::UK(:,:)
      COMPLEX(gq),ALLOCATABLE::HUK(:,:)
!
      ISYMI=SYM%IDI; ISYMF=SYM%IDF
      NSYM=ISYMF-ISYMI+1
      NAS=CO%DIMSO; NA2=CO%DIM2
      ALLOCATE(HUK(BND%NMAXIN,NAS))
      IKPL=0
!
      DO IVEC=1,GP%NVEC;   IF(GP%LKPVEC)THEN;     NKP=GP%KVEC(IVEC,2);   ELSE;     NKP=KPT%DIM;   ENDIF;   DO IKS=1,NKP;   IF(GP%LKPVEC)THEN;     IKP=GP%KVEC(IVEC,3)+IKS;     IF(GP%LOMP)THEN;       IKPL=IKP;     ELSE;       IKPL=IKPL+1;    ENDIF;   ELSE;     IKP=IKS;     IF(GP%LOMP)THEN;       IKPL=IKP;     ELSE;       IKPL=IKP-GP%MYRANK*KPT%DIML;     ENDIF;   ENDIF;   IF(IKPL.LE.0) CYCLE;   IF(IKPL.GT.KPT%DIML) EXIT
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
      CALL ORBITAL_SPIN_TRANS(CO%EL0,NA2,.TRUE.,BND%ISO) ! Orbital-fast to Spin-fast
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
        CALL GREEN_OCC_R1(1,1,CO%GF,OCC=CO%NKS(1:NASO*BND%NSPIN,1:NASO*BND%NSPIN))
      ENDIF
!
      IF(BND%ISPO.EQ.1)CO%NKS(1+NASO:,1+NASO:)=CO%NKS(1:NASO,1:NASO)
      CALL ORBITAL_SPIN_TRANS(CO%NKS,NA2,.TRUE.,BND%ISO)
      RETURN
!
      END SUBROUTINE CALC_CO_NKS
!
!****************************************************************************
      SUBROUTINE EVAL_CO_SL_VEC(CO,MODE)
      TYPE (CORR_ORB) CO
      INTEGER MODE
! LOCAL
      INTEGER I
      COMPLEX(gq) N_(CO%DIM2,CO%DIM2)
!
      IF(MODE==1)THEN
        N_=TRANSPOSE(CO%NKS)
      ELSE
        N_=TRANSPOSE(CO%NC_PHY)
      ENDIF
      DO I=1,3
        CO%S_VAL(I,MODE)=SUM(CO%S_VEC(:,:,I)*N_)
        CO%L_VAL(I,MODE)=SUM(CO%L_VEC(:,:,I)*N_)
      ENDDO
      RETURN
!
      END SUBROUTINE EVAL_CO_SL_VEC
!
!****************************************************************************
      SUBROUTINE CALC_CO_NKS_BND(CO)
      USE BANDSTRU
      TYPE (CORR_ORB) CO
! LOCAL
      INTEGER ISYMI,ISYMF,IB,ISYM,ISP,N1,N2
      INTEGER IVEC,IKP,IKPL,IKS,NKP,NSYM,NEMIN,NEMAX,NBANDS,NASO,NA2,MAX_IN_BANDS
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
      IKPL=0
      DO IVEC=1,GP%NVEC;   IF(GP%LKPVEC)THEN;     NKP=GP%KVEC(IVEC,2);   ELSE;     NKP=KPT%DIM;   ENDIF;   DO IKS=1,NKP;   IF(GP%LKPVEC)THEN;     IKP=GP%KVEC(IVEC,3)+IKS;     IF(GP%LOMP)THEN;       IKPL=IKP;     ELSE;       IKPL=IKPL+1;    ENDIF;   ELSE;     IKP=IKS;     IF(GP%LOMP)THEN;       IKPL=IKP;     ELSE;       IKPL=IKP-GP%MYRANK*KPT%DIML;     ENDIF;   ENDIF;   IF(IKPL.LE.0) CYCLE;   IF(IKPL.GT.KPT%DIML) EXIT
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
      CO%NKS(N1:N2,N1:N2)=CO%NKS(N1:N2,N1:N2)+TRANSPOSE(NABK)
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
      SUBROUTINE WRITE_V2H(CO,NT,IU)
      TYPE (CORR_ORB) CO
      INTEGER NT,IU
! LOCAL
      INTEGER I1,I2,I3,I4,NA2
!
      NA2=CO%DIM2
      WRITE(IU,'("NT=",I3)')NT
      DO I1=1,NA2; DO I2=1,NA2; DO I3=1,NA2; DO I4=1,NA2
        IF(ABS(CO%V2H(I1,I2,I3,I4))<1.E-10_gq)CYCLE
        WRITE(IU,'(4I3,2X,2F18.10)')I1,I2,I3,I4,REAL(CO%V2H(I1,I2,I3,I4)),AIMAG(CO%V2H(I1,I2,I3,I4))
      ENDDO; ENDDO; ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE WRITE_V2H
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
      SUBROUTINE READ_V2H(CO,NT,IU)
      TYPE (CORR_ORB) CO
      INTEGER NT,IU
! LOCAL
      INTEGER NA2,IT,I1,I2,I3,I4
      REAL(gq) RES,IES
      CHARACTER*3 STR1,STR2
!
      NA2=CO%DIM2
      ALLOCATE(CO%V2H(NA2,NA2,NA2,NA2)); CO%V2H=0
      OPEN(IU,FILE='V2H.INP',STATUS='OLD')
      DO 
        READ(IU,'(2A3)')STR1,STR2
        IF(STR1.NE.'NT=')CYCLE
        READ(STR2,*)IT
        IF(IT.NE.NT)CYCLE
        DO 
          READ(IU,*,END=100,ERR=100)I1,I2,I3,I4,RES,IES
          CO%V2H(I1,I2,I3,I4) = DCMPLX(RES, IES)
        ENDDO
      ENDDO
      STOP ' ERROR: CANNOT LOCATE V2H!'
100   CONTINUE
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE READ_V2H
!
!******************************************************************************
! Karamori convention, assuming real Harmonics.
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
      INTEGER I
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
      RD=MATMUL(CO%R1,TRANSPOSE(CO%D))
      LAY=0
      DO I=1,CO%DIM_HS_L
        H=TRANSPOSE(CO%HS_L(:,:,I))
        CALL PFA_PA(CO%NKS,PN,H,CO%DIM2,DSIMIX,DPSIMIX) ! p f / p d_n
        RES=-REAL(SUM(PN*RD),gq)*2
        LAY=LAY+CO%HS_L(:,:,I)*RES
      ENDDO
      LAY=LAY-LAX
      DO I=1,CO%DIM2
        LAY(I,I)=LAY(I,I)-CO%VDC
      ENDDO
      IF(LDC<0)THEN
        LAY=LAY-CO%VDC2
      ENDIF
      NULLIFY(LAX,LAY)
      RETURN
!
      END SUBROUTINE CALC_CO_LA
!
!****************************************************************************
      SUBROUTINE HM_EXPAND_CO(CO,FLAG,MODE)
      TYPE(CORR_ORB) CO
      CHARACTER*3,INTENT(IN)::FLAG
      INTEGER MODE
! LOCAL
      INTEGER NA2
      COMPLEX(gq),ALLOCATABLE :: COEF(:)
!
      NA2=CO%DIM2
      SELECT CASE(FLAG)
      CASE("RAA")
        ALLOCATE(COEF(CO%DIM_HS_R))
        IF(MODE<=0)COEF=CO%R_COEF
        CALL GET_HM_EXPAND(CO%R,CO%HS_R,NA2,CO%DIM_HS_R,COEF,MODE,.FALSE.)
        IF(MODE>0)THEN
          IF(GL%RMODE==2)THEN
            CO%R_COEF=COEF
          ELSE
            CO%R_COEF=REAL(COEF,gq)
          ENDIF
        ENDIF
      CASE("NKS")
        ALLOCATE(COEF(CO%DIM_HS_L))
        IF(MODE<=0)COEF=CO%NKS_COEF
        CALL GET_HM_EXPAND(CO%NKS,CO%HS_L,NA2,CO%DIM_HS_L,COEF,MODE,.TRUE.)
        IF(MODE>0)CO%NKS_COEF=REAL(COEF,gq)
      CASE("LA1")
        ALLOCATE(COEF(CO%DIM_HS_L))
        IF(MODE<=0)COEF=CO%LA1_COEF
        CALL GET_HM_EXPAND(CO%LA1,CO%HS_L,NA2,CO%DIM_HS_L,COEF,MODE,.FALSE.)
        IF(MODE>0)CO%LA1_COEF=REAL(COEF,gq)
      CASE("LA2")
        ALLOCATE(COEF(CO%DIM_HS_L))
        IF(MODE<=0)COEF=CO%LA2_COEF
        CALL GET_HM_EXPAND(CO%LA2,CO%HS_L,NA2,CO%DIM_HS_L,COEF,MODE,.FALSE.)
        IF(MODE>0)CO%LA2_COEF=REAL(COEF,gq)
      CASE("NCV")
        ALLOCATE(COEF(CO%DIM_HS_L))
        IF(MODE<=0)COEF=CO%NCV_COEF
        CALL GET_HM_EXPAND(CO%NC_VAR,CO%HS_L,NA2,CO%DIM_HS_L,COEF,MODE,.TRUE.)
        IF(MODE>0)CO%NCV_COEF=REAL(COEF,gq)
      CASE DEFAULT
        WRITE(0,'("ERROR IN HM_EXPAND_ALL: ILLEGAL FLAG=",A4)')FLAG; STOP
      END SELECT
      DEALLOCATE(COEF)
      RETURN
!
      END SUBROUTINE HM_EXPAND_CO
!
!****************************************************************************
      SUBROUTINE HM_EXPAND_COSYM(CO,FLAG)
      TYPE(CORR_ORB) CO
      CHARACTER*3,INTENT(IN)::FLAG
! LOCAL
      INTEGER NA2
      LOGICAL LTRANS
      COMPLEX(gq),ALLOCATABLE::COEF(:)
!
      NA2=CO%DIM2; LTRANS=.FALSE.
      SELECT CASE(FLAG)
      CASE("DA0")
        ALLOCATE(COEF(CO%DIM_HS_R))
        ! LTRANS=.TRUE. not necessary for the purpose of symmetrization.
        CALL GET_HM_EXPAND(CO%D0,CO%HS_R,NA2,CO%DIM_HS_R,COEF,1,LTRANS)
        CALL GET_HM_EXPAND(CO%D0,CO%HS_R,NA2,CO%DIM_HS_R,COEF,-1,LTRANS)
      CASE("EL0")
        ALLOCATE(COEF(CO%DIM_HS_E))
        CALL GET_HM_EXPAND(CO%EL0,CO%HS_E,NA2,CO%DIM_HS_E,COEF,1,LTRANS)
        CALL GET_HM_EXPAND(CO%EL0,CO%HS_E,NA2,CO%DIM_HS_E,COEF,-1,LTRANS)
      CASE("NCP")
        ALLOCATE(COEF(CO%DIM_HS))
        CALL GET_HM_EXPAND(CO%NC_PHY,CO%HS,NA2,CO%DIM_HS,COEF,1,LTRANS)
        CALL GET_HM_EXPAND(CO%NC_PHY,CO%HS,NA2,CO%DIM_HS,COEF,-1,LTRANS)
      CASE DEFAULT
        WRITE(0,'("ERROR IN HM_EXPAND_ALL: ILLEGAL FLAG=",A4)')FLAG; STOP
      END SELECT
      DEALLOCATE(COEF)
      RETURN
!
      END SUBROUTINE HM_EXPAND_COSYM
!
!****************************************************************************
! ISIMIX = (NKS(1-NKS))^(-1/2)
!****************************************************************************
      SUBROUTINE CALC_CO_ISIMIX(CO,NKS)
      TYPE(CORR_ORB) CO
      COMPLEX(gq) NKS(CO%DIM2,CO%DIM2)
! LOCAL
      COMPLEX(gq) XN(CO%DIM2,CO%DIM2)
! 
      XN=NKS; XN=XN-MATMUL(XN,XN)
      CALL ATOFA(XN,CO%ISIMIX,CO%DIM2,-12,D1,.TRUE.)
      RETURN
!
      END SUBROUTINE CALC_CO_ISIMIX
!
!****************************************************************************
! D = (NKS(1-NKS))^(-1/2) D0
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
!
      CO%EDCLA1=-REAL(SUM(CO%LA1*CO%NKS),gq) ! Missing transposition?
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
          CO%EDCUJ=CO%EDCUJ-REAL(CO%NPHY_FIX(I,J)*SUM(CO%NPHY_FIX*(CO%V2H(I,:,J,:)-CO%V2H(I,:,:,J))),gq)
        ENDDO; ENDDO
        IF(LDC==-12)THEN
          CO%EDCUJ=CO%EDCUJ-REAL(SUM((CO%NKS-CO%NPHY_FIX)*TRANSPOSE(CO%VDC2)),gq)
        ENDIF
      ENDIF
      RETURN
!
      END SUBROUTINE CALC_CO_EDCUJ
!
!****************************************************************************
      SUBROUTINE DUMP_EMBED_HAMILTONIAN(NI, CO)
      INTEGER NI
      TYPE(CORR_ORB) CO
! LOCAL
      INTEGER I,J,K,L,ISUM,NDIM
      COMPLEX(gq),ALLOCATABLE :: X_LIST(:)
      INTEGER,ALLOCATABLE :: I_LIST(:),J_LIST(:),K_LIST(:),L_LIST(:)
      REAL(gq),PARAMETER::RTOL=1.E-10_gq
! 
      CALL gh5_open(FILE_NAME(NI,0,0,8), embedH_file_id)
      CALL gh5_write(CO%DIM2*2, "/dim_cf", embedH_file_id)
      ! Two-body term V2H
      ISUM=0
      NDIM=CO%DIM2**4
      ALLOCATE(I_LIST(NDIM),J_LIST(NDIM),K_LIST(NDIM),L_LIST(NDIM), &
          &X_LIST(NDIM))
      I_LIST(1)=0; J_LIST(1)=0; K_LIST(1)=0; L_LIST(1)=0
      DO I=1,CO%DIM2; DO J=1,CO%DIM2; DO K=1,CO%DIM2; DO L=1,CO%DIM2
        IF(ABS(CO%V2H(I,J,K,L))<RTOL.OR.I==J.OR.K==L)CYCLE
        ISUM=ISUM+1
        I_LIST(ISUM)=I; J_LIST(ISUM)=J; K_LIST(ISUM)=K; L_LIST(ISUM)=L
      ENDDO; ENDDO; ENDDO; ENDDO
      DO I=1,ISUM
        X_LIST(I)=CO%V2H(I_LIST(I),J_LIST(I),K_LIST(I),L_LIST(I))
      ENDDO
      ! convert from one-based to zero-based
      I_LIST=I_LIST-1; J_LIST=J_LIST-1; K_LIST=K_LIST-1; L_LIST=L_LIST-1
      CALL gh5_write(X_LIST, MAX(ISUM,1), "/U_list", embedH_file_id)
      CALL gh5_write(I_LIST, MAX(ISUM,1), "/Ui_list", embedH_file_id)
      CALL gh5_write(J_LIST, MAX(ISUM,1), "/Uj_list", embedH_file_id)
      CALL gh5_write(K_LIST, MAX(ISUM,1), "/Uk_list", embedH_file_id)
      CALL gh5_write(L_LIST, MAX(ISUM,1), "/Ul_list", embedH_file_id)
      ! Or the 4-d array
      CALL gh5_write(CO%V2H,CO%DIM2,CO%DIM2,CO%DIM2,CO%DIM2,"/V2H", &
          &embedH_file_id)
      ! one-body term EL
      ISUM=0; X_LIST(1)=0
      DO I=1,CO%DIM2; DO J=1,CO%DIM2
        IF(ABS(CO%EL0(I,J))<RTOL)CYCLE
        ISUM=ISUM+1
        I_LIST(ISUM)=I-1; J_LIST(ISUM)=J-1
        X_LIST(ISUM)=CO%EL0(I,J)
        ENDDO; ENDDO
      CALL gh5_write(X_LIST, MAX(ISUM,1), "/e_list", embedH_file_id)
      CALL gh5_write(I_LIST, MAX(ISUM,1), "/ei_list", embedH_file_id)
      CALL gh5_write(J_LIST, MAX(ISUM,1), "/ej_list", embedH_file_id)
      ! Or the 2-d array
      CALL gh5_write(CO%EL0,CO%DIM2,CO%DIM2,"/H1E",embedH_file_id)
      ! D
      ISUM=0; X_LIST(1)=0
      DO I=1,CO%DIM2; DO J=1,CO%DIM2
        IF(ABS(CO%D(I,J))<RTOL)CYCLE
        ISUM=ISUM+1
        I_LIST(ISUM)=I-1; J_LIST(ISUM)=J-1
        X_LIST(ISUM)=CO%D(I,J)
      ENDDO; ENDDO
      CALL gh5_write(X_LIST, MAX(ISUM,1), "/D_list", embedH_file_id)
      CALL gh5_write(I_LIST, MAX(ISUM,1), "/Di_list", embedH_file_id)
      CALL gh5_write(J_LIST, MAX(ISUM,1), "/Dj_list", embedH_file_id)
      ! Or the 2-d array
      CALL gh5_write(CO%D,CO%DIM2,CO%DIM2,"/D",embedH_file_id)
      ! lambda_c
      ISUM=0; X_LIST(1)=0
      DO I=1,CO%DIM2; DO J=1,CO%DIM2
        IF(ABS(CO%LA2(I,J))<RTOL)CYCLE
        ISUM=ISUM+1
        I_LIST(ISUM)=I-1; J_LIST(ISUM)=J-1
        X_LIST(ISUM)=CO%LA2(I,J)
      ENDDO; ENDDO
      CALL gh5_write(X_LIST, MAX(ISUM,1), "/lambda_list", embedH_file_id)
      CALL gh5_write(I_LIST, MAX(ISUM,1), "/lambdai_list", embedH_file_id)
      CALL gh5_write(J_LIST, MAX(ISUM,1), "/lambdaj_list", embedH_file_id)
      ! Or the 2-d array
      CALL gh5_write(CO%LA2,CO%DIM2,CO%DIM2,"/Lambda_c",embedH_file_id)
      CALL gh5_close(embedH_file_id)
      DEALLOCATE(I_LIST,J_LIST,K_LIST,L_LIST,X_LIST)
      RETURN
!
      END SUBROUTINE DUMP_EMBED_HAMILTONIAN
!
!
      END MODULE CORRORB

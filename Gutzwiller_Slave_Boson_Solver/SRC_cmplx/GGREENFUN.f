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
      MODULE GREENFUN
      USE GPREC; USE GMPI
      USE GUTIL
      USE GCONSTANT
      IMPLICIT NONE
!
      TYPE GREEN_FUN
        REAL(gq) T,TEV,EWIN
        REAL(gq) ETA ! Zero+
        INTEGER NW,NBMAX,NSPIN,WTYP,NLR
        REAL(gq) WMAX,WMIN,MU(2)
        COMPLEX(gq),POINTER :: W(:) ! Frequency 
        COMPLEX(gq),POINTER :: G(:,:,:,:) ! Green function
        COMPLEX(gq),POINTER :: M(:,:,:,:) ! Self-energy
        COMPLEX(gq),POINTER :: D(:,:,:,:,:) ! Hybridization function
        COMPLEX(gq),POINTER :: D1(:) ! 1-band model
        REAL(gq) :: DS,TSP ! Half-band width, impurity-n.n. lead site hopping
        COMPLEX(gq),POINTER :: GM(:,:,:,:,:) ! \Gamma
        REAL(gq),POINTER :: WT(:)
      END TYPE GREEN_FUN
!
      TYPE(GREEN_FUN)::GF
!      
!
      CONTAINS
!
!******************************************************************************
      SUBROUTINE GF_INI(LGREEN,LUNIT,LMODEL,IU,IO)
      USE GCONSTANT
      INTEGER LGREEN,LUNIT,LMODEL,IU,IO
! LOCAL
      REAL(gq) XMIN,XMAX
      CHARACTER STR*2
!
      IF(LGREEN==0)RETURN
      IF(LUNIT==0)THEN ! eV UNIT
        XMIN=-30._gq; STR='eV'
        IF(GF%TEV>=0)THEN
          GF%T=GF%TEV
        ELSE
          GF%T=GF%T*KTOEV
        ENDIF
      ELSE
        GF%T=GF%T*KTORY; XMIN=-3._gq; STR='Ry'
      ENDIF
      XMAX=-XMIN
      IF(GF%ETA<=1.E-16_gq)THEN
        GF%ETA=GF%T
      ENDIF
      IF(GF%WTYP==0)THEN ! Complex frequency
        CALL READ_FREQ(IU)
      ELSEIF(GF%WTYP==1)THEN ! Real frequency
        CALL GEN_REAL_FREQ()
      ENDIF
      IF(GF%WTYP==1)THEN
        IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" GF: WMIN=",F12.6," WMAX=",F12.6," DS=",F12.6," TSP=",F12.6)') &
              &GF%WMIN,GF%WMAX,GF%DS,GF%TSP
      ENDIF
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'("   T=",E8.1,A3,", "," NOMEGA=",I7)')GF%T,STR,GF%NW
      GF%NLR=1  ! Default right-side lead
      IF(LMODEL==102)GF%NLR=2
      IF(GF%WTYP==0)THEN
        GF%W=GF%W*GF%T*2
        GF%EWIN=GF%NW*2*GF%T
      ENDIF
      RETURN
!
      END SUBROUTINE GF_INI
!
!******************************************************************************
      SUBROUTINE GF_HYBRD_INI()
!      
      ALLOCATE(GF%D1(GF%NW))
      CALL GET_BARE_HYBRIZ1(GF%D1,GF%DS)
      ALLOCATE(GF%D(GF%NBMAX,GF%NBMAX,GF%NW,GF%NSPIN,GF%NLR))
      IF(GF%WTYP==1)THEN
        ALLOCATE(GF%GM(GF%NBMAX,GF%NBMAX,GF%NW,GF%NSPIN,GF%NLR))
      ENDIF
      RETURN
!
      END SUBROUTINE GF_HYBRD_INI
!
!******************************************************************************
      SUBROUTINE READ_FREQ(IU)
      INTEGER IU
! LOCAL
      INTEGER IW
      REAL(gq) X,Y
!
      OPEN(IU,FILE='GFW.INP',STATUS='OLD')
      READ(IU,*)GF%NW
      ALLOCATE(GF%W(GF%NW)); GF%W=0
      DO IW=1,GF%NW,2
        READ(IU,*)X,Y
        GF%W(IW)=CMPLX(X,Y,KIND=gq)
        GF%W(IW+1)=-GF%W(IW)
      ENDDO
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE READ_FREQ
!
!******************************************************************************
      SUBROUTINE GEN_REAL_FREQ()
      REAL(gq) DW,OMG(GF%NW)
!
      CALL SET_RANGE(GF%WMIN,GF%WMAX,OMG,GF%NW)
      ALLOCATE(GF%W(GF%NW))
      GF%W=OMG+GF%ETA*ZI
      DW=REAL(GF%W(2)-GF%W(1),gq)
      ALLOCATE(GF%WT(GF%NW))
      CALL SET_LINEAR_SIMP(GF%WT,GF%NW,DW)
      RETURN
!
      END SUBROUTINE GEN_REAL_FREQ
!
!******************************************************************************
! Using local quasi-particle Green function
!******************************************************************************
      SUBROUTINE GREEN_OCC(MODE,GFL,NELE,OCC)
      INTEGER MODE
      TYPE(GREEN_FUN) :: GFL
      REAL(gq),OPTIONAL :: NELE
      COMPLEX(gq),OPTIONAL :: OCC(GFL%NBMAX*GFL%NSPIN,GFL%NBMAX*GFL%NSPIN)
! LOCAL
      INTEGER IW,IB,NDIM,NSPIN,ISP,N1,N2
      COMPLEX(gq),ALLOCATABLE :: MOCC(:,:,:)
!
      ! 'GREEN_OCC'
      NDIM=GFL%NBMAX; NSPIN=GFL%NSPIN
      ALLOCATE(MOCC(NDIM,NDIM,NSPIN)); MOCC=0
      DO IW=1,GFL%NW; MOCC=MOCC+GFL%G(:,:,IW,:); ENDDO
      MOCC=MOCC*GFL%T
      DO IB=1,NDIM; MOCC(IB,IB,:)=MOCC(IB,IB,:)+0.5_gq; ENDDO
      IF(MODE==0)THEN
        NELE=0
        DO IB=1,NDIM; NELE=NELE+REAL(SUM(MOCC(IB,IB,:)),gq); ENDDO
      ELSEIF(MODE==1)THEN
        DO ISP=1,GFL%NSPIN
          N1=(ISP-1)*NDIM+1; N2=ISP*NDIM
          OCC(N1:N2,N1:N2)=TRANSPOSE(MOCC(:,:,ISP))
        ENDDO
      ENDIF
      DEALLOCATE(MOCC)
      ! 'GREEN_OCC'
      RETURN
!
      END SUBROUTINE GREEN_OCC
!******************************************************************************
! Real frequency
! MODE = 0: FERMI_WEIGHT == 1, TEST NORMALIZATION
!        1: USUAL FERMI_WEIGHT
!******************************************************************************
      SUBROUTINE GREEN_OCC_R1(MODE,WMODE,GFL,NELE,OCC,NELET)
      TYPE(GREEN_FUN) :: GFL
      INTEGER MODE,WMODE
      REAL(gq),OPTIONAL :: NELE(GFL%NBMAX,GFL%NSPIN),NELET
      COMPLEX(gq),OPTIONAL :: OCC(GFL%NBMAX*GFL%NSPIN,GFL%NBMAX*GFL%NSPIN)
! LOCAL
      INTEGER ISP,IW,N1,N2,ILR
      REAL(gq) OMG,FW
      COMPLEX(gq),ALLOCATABLE :: MOCC(:,:),MOCC1(:,:)
!
      IF(PRESENT(NELE ))NELE =0
      IF(PRESENT(NELET))NELET=0
      IF(PRESENT(OCC))OCC=0
      FW=1._gq
      ALLOCATE(MOCC(GFL%NBMAX,GFL%NBMAX),MOCC1(GFL%NBMAX,GFL%NBMAX))
      DO ILR=1,GFL%NLR
      DO ISP=1,GFL%NSPIN
      MOCC=0
      DO IW=1,GFL%NW
        OMG=REAL(GFL%W(IW),gq)
        IF(MODE/=0)THEN
          IF(GFL%T>0)THEN
            OMG=(OMG-GFL%MU(ILR))/GFL%T
            FW=FERMI_FUN(OMG)
          ELSE
            IF(OMG<GFL%MU(ILR))THEN
              FW=1._gq
            ELSEIF(ABS(OMG-GFL%MU(ILR))<1.E-16_gq)THEN
              FW=.5_gq
            ELSE
              FW=0._gq
            ENDIF
          ENDIF
        ENDIF
        IF(FW<1.E-16_gq)CYCLE
        IF(WMODE==1)THEN
          MOCC1=MATMUL(GFL%G(:,:,IW,ISP),MATMUL(GFL%GM(:,:,IW,ISP,ILR),TRANSPOSE(CONJG(GFL%G(:,:,IW,ISP)))))
        ELSE
          MOCC1=MATMUL(GFL%G(:,:,IW,ISP),MATMUL(&
                     &(TRANSPOSE(CONJG(GFL%D(:,:,IW,ISP,ILR)))-GFL%D(:,:,IW,ISP,ILR)), &
                     & TRANSPOSE(CONJG(GFL%G(:,:,IW,ISP)))))/ZITPI
        ENDIF
        MOCC=MOCC+MOCC1*FW*GF%WT(IW)
      ENDDO ! IW
      IF(PRESENT(OCC))THEN
        N1=(ISP-1)*GFL%NBMAX+1; N2=ISP*GFL%NBMAX
        OCC(N1:N2,N1:N2)=OCC(N1:N2,N1:N2)+TRANSPOSE(MOCC)
      ENDIF
      IF(PRESENT(NELE))THEN
        DO N1=1,GFL%NBMAX
          NELE(N1,ISP)=NELE(N1,ISP)+REAL(MOCC(N1,N1),gq)
        ENDDO
      ENDIF
      IF(PRESENT(NELET))THEN
        DO N1=1,GFL%NBMAX
          NELET=NELET+REAL(MOCC(N1,N1),gq)
        ENDDO
      ENDIF
      ENDDO ! ISP
      ENDDO ! ILR
      RETURN
!
      END SUBROUTINE GREEN_OCC_R1
!
!******************************************************************************
! Using coherent part of the local Green function, to be revised.
!******************************************************************************
      SUBROUTINE GREEN_OCC2(GFL,OCC,R)
      TYPE(GREEN_FUN) :: GFL
      COMPLEX(gq) :: OCC(GFL%NBMAX*GFL%NSPIN,GFL%NBMAX*GFL%NSPIN),R(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN)
! LOCAL
      INTEGER IW,IB,NDIM,NSPIN,ISP,N1,N2
      COMPLEX(gq),ALLOCATABLE :: MOCC(:,:,:),RINV(:,:)
!
      ! 'GREEN_OCC2'
      NDIM=GFL%NBMAX; NSPIN=GFL%NSPIN
      ALLOCATE(MOCC(NDIM,NDIM,NSPIN)); MOCC=0
!
      DO IW=1,GFL%NW; MOCC=MOCC+GFL%G(:,:,IW,:); ENDDO
      MOCC=MOCC*GFL%T
!
      ALLOCATE(RINV(NDIM,NDIM))
      DO ISP=1,GFL%NSPIN
        RINV=R(:,:,ISP)
        CALL INV(RINV,NDIM)
        CALL ANNXB('C','N',RINV,MOCC(:,:,ISP),NDIM,2)
        CALL ANNXB('N','N',MOCC(:,:,ISP),RINV,NDIM,2) ! 1/R^H T sum{G} 1/R
         DO IB=1,NDIM; MOCC(IB,IB,ISP)=MOCC(IB,IB,ISP)+0.5_gq; ENDDO
      ENDDO ! ISP
      DO ISP=1,GFL%NSPIN
        N1=(ISP-1)*NDIM+1; N2=ISP*NDIM
        OCC(N1:N2,N1:N2)=TRANSPOSE(MOCC(:,:,ISP))
      ENDDO
      DEALLOCATE(MOCC)
      ! 'GREEN_OCC2'
      RETURN
!
      END SUBROUTINE GREEN_OCC2
!
!******************************************************************************
      SUBROUTINE GREEN_DA(GFL,MU,R,D)
      TYPE(GREEN_FUN) :: GFL
      COMPLEX(gq) :: MU(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN), &
                       & R(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN), &
                       & D(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN)
! LOCAL
      INTEGER IW,IB,NDIM,ISP
      COMPLEX(gq) TMP(GFL%NBMAX,GFL%NBMAX)
      COMPLEX(gq) RD(GFL%NBMAX,GFL%NBMAX),RINV(GFL%NBMAX,GFL%NBMAX)
!
      ! 'GREEN_DA'
      NDIM=GFL%NBMAX
      DO ISP=1,GFL%NSPIN
        RD=0
        DO IW=1,GFL%NW
          TMP=-MU(:,:,ISP)
          DO IB=1,NDIM; TMP(IB,IB)=TMP(IB,IB)+GFL%W(IW); ENDDO ! W-lambda-eta+Ef
          CALL ANNXB('N','N',TMP,GFL%G(:,:,IW,ISP),NDIM,1)
          RD=RD+TMP
        ENDDO
        DO IB=1,NDIM; RD(IB,IB)=RD(IB,IB)-GFL%NW; ENDDO
        RD=RD*GFL%T
        RINV=R(:,:,ISP)
        CALL INV(RINV,NDIM)
        CALL ANNXB('N','N',RINV,RD,NDIM,2)
        D(:,:,ISP)=TRANSPOSE(RD)
      ENDDO
      ! 'GREEN_DA'
      RETURN
!
      END SUBROUTINE GREEN_DA
!
!******************************************************************************
      SUBROUTINE GREEN_DA_R1(GFL,R,D,WMODE)
      TYPE(GREEN_FUN) :: GFL
      COMPLEX(gq) :: R(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN), &
                       &D(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN)
      INTEGER WMODE ! Various ways for the same calculation.
! LOCAL
      INTEGER ISP,IW,I,ILR
      REAL(gq) OMG,FW
      COMPLEX(gq) RD(GFL%NBMAX,GFL%NBMAX),RINV(GFL%NBMAX,GFL%NBMAX)
      COMPLEX(gq) TMP(GFL%NBMAX,GFL%NBMAX)
!
      D=0
      DO ILR=1,GFL%NLR
      DO ISP=1,GFL%NSPIN
      RD=0
      DO IW=1,GFL%NW
        OMG=REAL(GFL%W(IW),gq)
        IF(GFL%T>0)THEN
          OMG=(OMG-GFL%MU(ILR))/GFL%T
          FW=FERMI_FUN(OMG)
        ELSE
          IF(OMG<GFL%MU(ILR))THEN
            FW=1._gq
          ELSEIF(ABS(OMG-GFL%MU(ILR))<1.E-16_gq)THEN
            FW=.5_gq
          ELSE
            FW=0._gq
          ENDIF
        ENDIF
        IF(FW<1.E-16_gq)CYCLE
        IF(WMODE==1.OR.WMODE==2)THEN
          TMP=0; DO I=1,GFL%NLR; TMP=TMP+GFL%D(:,:,IW,ISP,I); ENDDO ! Total \Delta
          TMP=MATMUL(TMP,GFL%G(:,:,IW,ISP))
          DO I=1,GFL%NBMAX; TMP(I,I)=TMP(I,I)+1._gq; ENDDO
          IF(WMODE==1)THEN
            TMP=MATMUL(TMP,MATMUL(GFL%GM(:,:,IW,ISP,ILR),TRANSPOSE(CONJG(GFL%G(:,:,IW,ISP)))))
            RD=RD+TMP*FW*GF%WT(IW)
          ELSE
            TMP=MATMUL(TMP,MATMUL( &
               &(TRANSPOSE(CONJG(GFL%D(:,:,IW,ISP,ILR)))-GFL%D(:,:,IW,ISP,ILR)), &
               &TRANSPOSE(CONJG(GFL%G(:,:,IW,ISP)))))
            RD=RD+TMP*FW*GF%WT(IW)/ZITPI
          ENDIF
        ELSE
          TMP=-MATMUL(GFL%D(:,:,IW,ISP,ILR),GFL%G(:,:,IW,ISP)) &
           &+MATMUL(TRANSPOSE(CONJG(GFL%D(:,:,IW,ISP,ILR))),TRANSPOSE(CONJG(GFL%G(:,:,IW,ISP))))
           RD=RD+TMP*FW*GF%WT(IW)/ZITPI
        ENDIF
      ENDDO ! IW
      RINV=R(:,:,ISP)
      CALL INV(RINV,GFL%NBMAX)
      CALL ANNXB('N','N',RINV,RD,GFL%NBMAX,2)
      D(:,:,ISP)=D(:,:,ISP)+TRANSPOSE(RD)
      ENDDO ! ISP
      ENDDO ! ILR
      RETURN
!
      END SUBROUTINE GREEN_DA_R1
!
!******************************************************************************
      SUBROUTINE GREEN_DA2(GFL,R,D,EF)
      TYPE(GREEN_FUN) :: GFL
      COMPLEX(gq) :: R(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN), &
                       &D(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN)
      REAL(gq) :: EF
! LOCAL
      INTEGER IW,IB,NDIM,ISP
      COMPLEX(gq) RD(GFL%NBMAX,GFL%NBMAX),RINV(GFL%NBMAX,GFL%NBMAX)
      COMPLEX(gq) TMP(GFL%NBMAX,GFL%NBMAX)
!
      ! 'GREEN_DA2'
      NDIM=GFL%NBMAX
      DO ISP=1,GFL%NSPIN
        RD=0
        DO IW=1,GFL%NW
          TMP=-GFL%M(:,:,IW,ISP)
          DO IB=1,NDIM; TMP(IB,IB)=TMP(IB,IB)+GFL%W(IW)+EF; ENDDO ! z + Ef - Sigma
          CALL ANNXB('N','N',TMP,GFL%G(:,:,IW,ISP),NDIM,1)
          RD=RD+TMP
        ENDDO
        DO IB=1,NDIM; RD(IB,IB)=RD(IB,IB)-GF%NW; ENDDO ! -PI_i
        RD=RD*GF%T
        RINV=R(:,:,ISP)
        CALL INV(RINV,NDIM)
        CALL ANNXB('N','N',RINV,RD,NDIM,2)
        D(:,:,ISP)=TRANSPOSE(RD)
      ENDDO
      ! 'GREEN_DA2'
      RETURN
!
      END SUBROUTINE GREEN_DA2
!
!=============================================================================
      SUBROUTINE PLOT_GF(IU)
      INTEGER IU
! LOCAL
      INTEGER I1,I2,IW
!
      IF(GP%MYRANK.NE.GP%MASTER) RETURN
      DO I1=1,GF%NBMAX; DO I2=1,GF%NBMAX
      OPEN(IU,FILE=FILE_NAME(0,I1,I2,19),STATUS='REPLACE')
      DO IW=1,GF%NW
        IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IU,'(F16.9,2F12.6)')REAL(GF%W(IW),gq),REAL(GF%G(I1,I2,IW,1),gq),AIMAG(GF%G(I1,I2,IW,1))
      ENDDO
      CLOSE(IU)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE PLOT_GF
!
!=============================================================================
      SUBROUTINE CALC_GF_1K(HK,GK,N1,N2,MU,OMG,ETOP,N3,SIGMA)
      INTEGER N1,N2,N3
      REAL(gq) MU,ETOP
      COMPLEX(gq) OMG
      COMPLEX(gq) HK(N1,N1),GK(N2,N2)
      COMPLEX(gq),OPTIONAL :: SIGMA(N3,N3)
! LOCAL
      INTEGER I
!
      GK=0; GK(1:N1,1:N1)=-HK
      DO I=N1+1,N2; GK(I,I)=-ETOP; ENDDO
      DO I=1,N2; GK(I,I)=GK(I,I)+OMG+MU; ENDDO
      IF(PRESENT(SIGMA))THEN
        GK(1:N3,1:N3)=GK(1:N3,1:N3)-SIGMA
      ENDIF
      CALL INV(GK,N2)
      RETURN
!
      END SUBROUTINE CALC_GF_1K
!
!=============================================================================
      SUBROUTINE CALC_GF_1KS(GK0,GK,N,MU)
      USE GUTIL
      INTEGER N
      REAL(gq) MU
      COMPLEX(gq) GK0(N,N),GK(N,N)
!
      CALL INV_APLX(GK0,GK,N,MU)
      RETURN
!
      END SUBROUTINE CALC_GF_1KS
!
!=============================================================================
      SUBROUTINE CALC_SELF_ENERGY(GFL,R,LA,EF)
      TYPE(GREEN_FUN) :: GFL
      COMPLEX(gq) LA(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN),R(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN)
      REAL(gq) :: EF
! LOCAL
      INTEGER ISP,IW,IB
      COMPLEX(gq) RHR(GFL%NBMAX,GFL%NBMAX),RINV(GFL%NBMAX,GFL%NBMAX),RILARIH(GFL%NBMAX,GFL%NBMAX)
!
      ! 'CALC_SELF_ENERGY'
      DO ISP=1,GFL%NSPIN
      RINV=R(:,:,ISP)
      CALL INV(RINV,GFL%NBMAX) ! 1/R
      RILARIH=LA(:,:,ISP)
      CALL ANNXB('N','N',RINV,RILARIH,GFL%NBMAX,2)
      CALL ANNXB('N','C',RILARIH,RINV,GFL%NBMAX,1) ! 1/R (La + eta) 1/R^H
      RHR=R(:,:,ISP)
      CALL ANNXB('C','N',R(:,:,ISP),RHR,GFL%NBMAX,2) ! R^H R
      RHR=-RHR
      DO IB=1,GFL%NBMAX; RHR(IB,IB)=RHR(IB,IB)+1._gq; ENDDO ! 1-R^H R
      CALL ANNXB('N','N',RINV,RHR,GFL%NBMAX,2)
      CALL ANNXB('N','C',RHR,RINV,GFL%NBMAX,1) ! 1/R (1-R^h R) 1/R^H
      DO IW=1,GFL%NW
        GFL%M(:,:,IW,ISP)=RILARIH-(GFL%W(IW)+EF)*RHR
      ENDDO
      ENDDO ! ISP
      ! 'CALC_SELF_ENERGY'
      RETURN
!
      END SUBROUTINE CALC_SELF_ENERGY
!
!******************************************************************************
! R \Gamma R^+; renormalized semi-circular density of states
!******************************************************************************
      SUBROUTINE GET_RENORM_GAMMA1(GM,R,N,JN)
      INTEGER N,JN
      COMPLEX(gq) R(N,N,GF%NSPIN),GM(N,N,GF%NW,GF%NSPIN)
! LOCAL
      INTEGER IW,I,ISP
      REAL(gq) OMG,FAKT
!
      GM=0; FAKT=(GF%TSP/GF%DS*2)**2
      DO ISP=1,GF%NSPIN; DO IW=1,GF%NW
      OMG=REAL(GF%W(IW),gq)
      IF(ABS(OMG)>=GF%DS)CYCLE
      DO I=1,N
        GM(I,:,IW,ISP)=R(I,JN,ISP)*SQRT((GF%DS)**2-OMG**2)*CONJG(R(:,JN,ISP))/2/PI*FAKT
        GM(I,I,IW,ISP)=GM(I,I,IW,ISP)+GF%ETA/PI
      ENDDO; ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE GET_RENORM_GAMMA1
!
!******************************************************************************
      SUBROUTINE GET_RENORM_HYBRIZ1(DELTA,R,RDR,N,JN)
      INTEGER N,JN
      COMPLEX(gq) DELTA(GF%NW),RDR(N,N,GF%NW,GF%NSPIN) ! DELTA(1,1,NW)
      COMPLEX(gq) R(N,N,GF%NSPIN)
! LOCAL
      INTEGER IW,I,ISP
! 
      DO ISP=1,GF%NSPIN; DO IW=1,GF%NW; DO I=1,N
        RDR(I,:,IW,ISP)=R(I,JN,ISP)*DELTA(IW)*CONJG(R(:,JN,ISP)) ! Right lead
      ENDDO; ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE GET_RENORM_HYBRIZ1
!
!******************************************************************************
!  Linear one-band nearest neighbour hopping; Ds=2ts
!******************************************************************************
      SUBROUTINE GET_BARE_HYBRIZ1(DELTA,DS)
      COMPLEX(gq) DELTA(GF%NW)
      REAL(gq) DS ! Half width of the semicircular DOS = 2ts
! LOCAL
      INTEGER I
      REAL(gq) SGN
      COMPLEX(gq) ZTMP,FAKT
!
      FAKT=(GF%TSP/GF%DS*2)**2  ! (t'/t)^2
      DO I=1,GF%NW
        ZTMP=GF%W(I)**2-DS**2
        SGN=REAL(GF%W(I),gq)*AIMAG(GF%W(I))
        IF(SGN>=0)THEN; SGN=1._gq; ELSE; SGN=-1._gq; ENDIF
        ZTMP=SQRT(ZTMP)*SGN
        DELTA(I)=(GF%W(I)-ZTMP)/2*FAKT
      ENDDO
      RETURN
!
      END SUBROUTINE GET_BARE_HYBRIZ1
!
!******************************************************************************
      SUBROUTINE CHK_W_ETA(IO)
      INTEGER IO
! LOCAL
      INTEGER IW
      REAL(gq) RINT
!
      RINT=0
      DO IW=1,GF%NW
        RINT=RINT-GF%WT(IW)*AIMAG(1._gq/(GF%W(IW)+GF%ETA*ZI))
      ENDDO
      RINT=RINT/PI
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" GF%ETA=",E10.2," NW=",I10," NORMALIZATION(1) CHECK=",F12.8)')GF%ETA,GF%NW,RINT
      RETURN
      END SUBROUTINE CHK_W_ETA
!
!******************************************************************************
      SUBROUTINE TEST_GFERMI(XMIN,XMAX,IU)
      REAL(gq) XMIN,XMAX
      INTEGER IU
! LOCAL
      INTEGER  I
      INTEGER,PARAMETER :: NSTEP=1000
      REAL(gq) X,DELTA
      COMPLEX(gq) Y
!
      ! 'TEST_GFERMI'
      DELTA=(XMAX-XMIN)/NSTEP
      OPEN(IU,FILE='TESTGF.OUT',STATUS='REPLACE')
      X=XMIN
      DO I=1,NSTEP
        X=X+DELTA
        Y=GFERMI(X)
        IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IU,*)X,REAL(Y,KIND=gq),AIMAG(Y)
      ENDDO
      ! 'TEST_GFERMI'
      RETURN
!
      END SUBROUTINE TEST_GFERMI
!
      FUNCTION GFERMI(X)
      INTEGER IW
      REAL(gq) X
      COMPLEX(gq) GFERMI
!
      GFERMI=0
      DO IW=1,GF%NW; GFERMI=GFERMI+1/(GF%W(IW)-X); ENDDO
      GFERMI=GFERMI*GF%T+0.5_gq
      RETURN
!
      END FUNCTION GFERMI
!
      END MODULE GREENFUN

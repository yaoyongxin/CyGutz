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
      MODULE GREENFUN2
      USE GPREC; USE GMPI; USE GUTIL; USE GCONSTANT
      IMPLICIT NONE
!
      TYPE GREEN_FUN2
        REAL(gq) T,TEV,EWIN
        REAL(gq) ETA ! Zero+
        INTEGER NW,NBMAX,NSPIN,NLR,NV       ! NV: number of function values to be obtained by integration.
        REAL(gq) WMAX,WMIN,MU(2)
        COMPLEX(gq),POINTER :: G(:,:,:)     ! Green function
        COMPLEX(gq),POINTER :: D(:,:,:,:)   ! Hybridization function, R * D1 * R^{\dagger}
        COMPLEX(gq),POINTER :: D1R(:,:,:,:) ! D1 * R^{\dagger}
        COMPLEX(gq) :: D1 ! 1-band model
        REAL(gq) :: DS,TSP ! Half-band width, impurity-n.n. lead site hopping
        COMPLEX(gq),POINTER :: GM(:,:,:,:)   ! R \Gamma_1 R^{\dagger}  
        COMPLEX(gq),POINTER :: GM1R(:,:,:,:) !   \Gamma_1 R^{\dagger}
      END TYPE GREEN_FUN2
!
      TYPE(GREEN_FUN2)::GF2
!      
!
      CONTAINS
!
!******************************************************************************
      SUBROUTINE GF2_INI(LUNIT, LMODEL, IO)
      INTEGER,INTENT(IN) :: LUNIT, LMODEL, IO
! LOCAL
      CHARACTER STR*2
!
      IF(LUNIT==0)THEN ! eV UNIT
        STR='eV'
        IF(GF2%TEV>=0)THEN
          GF2%T=GF2%TEV
        ELSE
          GF2%T=GF2%T*KTOEV
        ENDIF
      ELSE
        GF2%T=GF2%T*KTORY; STR='Ry'
      ENDIF
      IF(GF2%ETA<=1.E-16_gq)THEN
        GF2%ETA=GF2%T
      ENDIF
      GF2%NLR=1  ! Default right-side lead
      IF(LMODEL==102)GF2%NLR=2
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'("   T=",E8.1,A3,", "," Adaptive frequency sampling.")')GF2%T,STR
      RETURN
!
      END SUBROUTINE GF2_INI
!
!******************************************************************************
      SUBROUTINE GF2_HYBRD_INI()
!      
      ALLOCATE(GF2%D(GF2%NBMAX,GF2%NBMAX,GF2%NSPIN,GF2%NLR))
      ALLOCATE(GF2%D1R(GF2%NBMAX,GF2%NBMAX,GF2%NSPIN,GF2%NLR))
      ALLOCATE(GF2%GM(GF2%NBMAX,GF2%NBMAX,GF2%NSPIN,GF2%NLR))
      ALLOCATE(GF2%GM1R(GF2%NBMAX,GF2%NBMAX,GF2%NSPIN,GF2%NLR))
      RETURN
!
      END SUBROUTINE GF2_HYBRD_INI
!
!******************************************************************************
      SUBROUTINE GREEN2_OCC_R1(MODE,WMODE,GFL,OMG,OCC,NELE)
      TYPE(GREEN_FUN2) :: GFL
      INTEGER,INTENT(IN) :: MODE,WMODE
      COMPLEX(gq),INTENT(IN) :: OMG
      COMPLEX(gq),OPTIONAL,INTENT(OUT) :: OCC(GFL%NBMAX*GFL%NSPIN,GFL%NBMAX*GFL%NSPIN)
      REAL(gq),OPTIONAL,INTENT(OUT) :: NELE(GFL%NBMAX,GFL%NSPIN)
! LOCAL
      INTEGER ISP,N1,N2,ILR
      REAL(gq) FW, ROMG
      COMPLEX(gq),ALLOCATABLE :: MOCC(:,:)
!
      ALLOCATE(MOCC(GFL%NBMAX,GFL%NBMAX))
      DO ILR=1,GFL%NLR
      IF(MODE/=0)THEN
        ROMG = REAL(OMG)
        IF(GFL%T>0)THEN
          ROMG=(ROMG-GFL%MU(ILR))/GFL%T
          FW=FERMI_FUN(ROMG)
        ELSE
          IF(ROMG<GFL%MU(ILR))THEN
            FW=1._gq
          ELSEIF(ABS(ROMG-GFL%MU(ILR))<1.E-16_gq)THEN
            FW=.5_gq
          ELSE
            FW=0._gq
          ENDIF
        ENDIF
      ELSE
        FW = 1._gq
      ENDIF
      IF(FW<1.E-16_gq)CYCLE
      DO ISP=1,GFL%NSPIN
      MOCC=0
      IF(WMODE==1)THEN
        MOCC=MATMUL(GFL%G(:,:,ISP),MATMUL(GFL%GM(:,:,ISP,ILR),TRANSPOSE(CONJG(GFL%G(:,:,ISP)))))
      ELSE
        MOCC=MATMUL(GFL%G(:,:,ISP),MATMUL(&
                   &(TRANSPOSE(CONJG(GFL%D(:,:,ISP,ILR)))-GFL%D(:,:,ISP,ILR)), &
                   & TRANSPOSE(CONJG(GFL%G(:,:,ISP)))))/ZITPI
      ENDIF
      MOCC=MOCC*FW
      IF(PRESENT(OCC))THEN
        N1=(ISP-1)*GFL%NBMAX+1; N2=ISP*GFL%NBMAX
        OCC(N1:N2,N1:N2)=TRANSPOSE(MOCC)
      ENDIF
      IF(PRESENT(NELE))THEN
        DO N1=1,GFL%NBMAX
          NELE(N1,ISP)=REAL(MOCC(N1,N1),gq)
        ENDDO
      ENDIF
      ENDDO ! ISP
      ENDDO ! ILR
      RETURN
!
      END SUBROUTINE GREEN2_OCC_R1
!
!******************************************************************************
      SUBROUTINE GREEN2_DA_R1(GFL,OMG,R,D,WMODE)
      TYPE(GREEN_FUN2) :: GFL
      COMPLEX(gq),INTENT(IN) :: R(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN)
      COMPLEX(gq),INTENT(OUT) :: D(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN)
      INTEGER,INTENT(IN) :: WMODE ! Various ways for the same calculation.
      COMPLEX(gq),INTENT(IN) :: OMG
! LOCAL
      INTEGER ISP,I,ILR
      REAL(gq) FW, ROMG
      COMPLEX(gq) RD(GFL%NBMAX,GFL%NBMAX),RINV(GFL%NBMAX,GFL%NBMAX)
      COMPLEX(gq) TMP(GFL%NBMAX,GFL%NBMAX)
!
      D=0
      DO ILR=1,GFL%NLR
      ROMG = REAL(OMG)
      IF(GFL%T>0)THEN
        ROMG=(ROMG-GFL%MU(ILR))/GFL%T
        FW=FERMI_FUN(ROMG)
      ELSE
        IF(ROMG<GFL%MU(ILR))THEN
          FW=1._gq
        ELSEIF(ABS(ROMG-GFL%MU(ILR))<1.E-16_gq)THEN
          FW=.5_gq
        ELSE
          FW=0._gq
        ENDIF
      ENDIF
      IF(FW<1.E-16_gq)CYCLE
      DO ISP=1,GFL%NSPIN
        IF(WMODE==1.OR.WMODE==2)THEN
          TMP=0; DO I=1,GFL%NLR; TMP=TMP+GFL%D(:,:,ISP,I); ENDDO ! Total \Delta
          TMP=MATMUL(TMP,GFL%G(:,:,ISP))
          ! 1+\Delta*G
          DO I=1,GFL%NBMAX; TMP(I,I)=TMP(I,I)+1._gq; ENDDO
          IF(WMODE==1)THEN
            ! (1 + \Delta*G) * \Gamma * G^{\dagger}
            TMP=MATMUL(TMP,MATMUL(GFL%GM(:,:,ISP,ILR),TRANSPOSE(CONJG(GFL%G(:,:,ISP)))))
            RD=TMP*FW
          ELSE
            ! (1 + \Delta*G) * (\Delta^{\dagger} - \Delta) * G^{\dagger}
            TMP=MATMUL(TMP,MATMUL( &
               &(TRANSPOSE(CONJG(GFL%D(:,:,ISP,ILR)))-GFL%D(:,:,ISP,ILR)), &
               &TRANSPOSE(CONJG(GFL%G(:,:,ISP)))))
            RD=TMP*FW/ZITPI
          ENDIF
        ELSE
          TMP=-MATMUL(GFL%D(:,:,ISP,ILR),GFL%G(:,:,ISP)) &
           &+MATMUL(TRANSPOSE(CONJG(GFL%D(:,:,ISP,ILR))),TRANSPOSE(CONJG(GFL%G(:,:,ISP))))
           RD=TMP*FW/ZITPI
        ENDIF
        RINV=R(:,:,ISP)
        CALL INV(RINV,GFL%NBMAX)
        CALL ANNXB('N','N',RINV,RD,GFL%NBMAX,2)
        D(:,:,ISP)=TRANSPOSE(RD)
      ENDDO ! ISP
      ENDDO ! ILR
      RETURN
!
      END SUBROUTINE GREEN2_DA_R1
!
!******************************************************************************
      SUBROUTINE GREEN2_DA_R2(GFL,OMG,R,D)
      TYPE(GREEN_FUN2) :: GFL
      COMPLEX(gq),INTENT(IN) :: R(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN)
      COMPLEX(gq),INTENT(OUT) :: D(GFL%NBMAX,GFL%NBMAX,GFL%NSPIN)
      COMPLEX(gq),INTENT(IN) :: OMG
! LOCAL
      INTEGER ISP,I,J,ILR,IDX(GFL%NBMAX),ISUM
      REAL(gq) FW, ROMG
      COMPLEX(gq) RD(GFL%NBMAX,GFL%NBMAX), &
                    &RINV(GFL%NBMAX,GFL%NBMAX),RINV_(GFL%NBMAX,GFL%NBMAX)
      COMPLEX(gq) TMP(GFL%NBMAX,GFL%NBMAX)
!
      D=0
      DO ILR=1,GFL%NLR
      ROMG = REAL(OMG)
      IF(GFL%T>0)THEN
        ROMG=(ROMG-GFL%MU(ILR))/GFL%T
        FW=FERMI_FUN(ROMG)
      ELSE
        IF(ROMG<GFL%MU(ILR))THEN
          FW=1._gq
        ELSEIF(ABS(ROMG-GFL%MU(ILR))<1.E-16_gq)THEN
          FW=.5_gq
        ELSE
          FW=0._gq
        ENDIF
      ENDIF
      IF(FW<1.E-16_gq)CYCLE
      DO ISP=1,GFL%NSPIN
        ! Total \Delta_1 * R^{\dagger}
        TMP=0; DO I=1,GFL%NLR; TMP=TMP+GFL%D1R(:,:,ISP,I); ENDDO
        TMP=MATMUL(TMP,GFL%G(:,:,ISP)) ! \Delta_1 * R^{\dagger} * G
        ! \Delta_1 * R^{\dagger} * G * \Gamma * G^{\dagger}
        TMP=MATMUL(TMP,MATMUL(GFL%GM(:,:,ISP,ILR),TRANSPOSE(CONJG(GFL%G(:,:,ISP)))))
        TMP=TMP+MATMUL(GFL%GM1R(:,:,ISP,ILR),TRANSPOSE(CONJG(GFL%G(:,:,ISP))))
        RD=TMP*FW
        D(:,:,ISP)=TRANSPOSE(RD)
      ENDDO ! ISP
      ENDDO ! ILR
      RETURN
!
      END SUBROUTINE GREEN2_DA_R2
!
!******************************************************************************
! R \Gamma R^+; renormalized semi-circular density of states
!******************************************************************************
      SUBROUTINE GET_RENORM_GAMMA21(GM,OMG,R,N,JN)
      INTEGER,INTENT(IN) :: N,JN
      COMPLEX(gq),INTENT(IN) :: OMG
      COMPLEX(gq),INTENT(IN) :: R(N,N,GF2%NSPIN)
      COMPLEX(gq),INTENT(OUT) :: GM(N,N,GF2%NSPIN)
! LOCAL
      INTEGER I,ISP
      REAL(gq) FAKT,ROMG
!
      GM=0; FAKT=(GF2%TSP/GF2%DS*2)**2
      ROMG = REAL(OMG)
      IF(ABS(ROMG)>=GF2%DS)RETURN
      DO ISP=1,GF2%NSPIN; DO I=1,N
        GM(I,:,ISP)=R(I,JN,ISP)*SQRT((GF2%DS)**2-ROMG**2)*CONJG(R(:,JN,ISP))/2/PI*FAKT
        GM(I,I,ISP)=GM(I,I,ISP)+GF2%ETA/PI
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE GET_RENORM_GAMMA21
!
!******************************************************************************
! \Gamma R^+; 1-side-renormalized semi-circular density of states
!******************************************************************************
      SUBROUTINE GET_RENORM_GAMMA11(GM,OMG,R,N,JN)
      INTEGER,INTENT(IN) :: N,JN
      COMPLEX(gq),INTENT(IN) :: OMG
      COMPLEX(gq),INTENT(IN) :: R(N,N,GF2%NSPIN)
      COMPLEX(gq),INTENT(OUT) :: GM(N,N,GF2%NSPIN)
! LOCAL
      INTEGER I,ISP
      REAL(gq) FAKT,ROMG
!
      GM=0; FAKT=(GF2%TSP/GF2%DS*2)**2
      ROMG = REAL(OMG)
      IF(ABS(ROMG)>=GF2%DS)RETURN
      DO ISP=1,GF2%NSPIN
        GM(JN,:,ISP)=SQRT((GF2%DS)**2-ROMG**2)*CONJG(R(:,JN,ISP))/2/PI*FAKT
        DO I=1,N
          GM(I,I,ISP)=GM(I,I,ISP)+GF2%ETA/PI
        ENDDO
      ENDDO
      RETURN
!
      END SUBROUTINE GET_RENORM_GAMMA11
!
!******************************************************************************
      SUBROUTINE GET_RENORM_HYBRIZ21(DELTA,R,RDR,N,JN)
      INTEGER N,JN
      COMPLEX(gq) DELTA,RDR(N,N,GF2%NSPIN) ! DELTA(1,1,NW)
      COMPLEX(gq) R(N,N,GF2%NSPIN)
! LOCAL
      INTEGER I,ISP
! 
      DO ISP=1,GF2%NSPIN; DO I=1,N
        RDR(I,:,ISP)=R(I,JN,ISP)*DELTA*CONJG(R(:,JN,ISP)) ! Right lead
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE GET_RENORM_HYBRIZ21
!
!******************************************************************************
      SUBROUTINE GET_RENORM_HYBRIZ11(DELTA,R,DR,N,JN)
      INTEGER N,JN
      COMPLEX(gq) DELTA,DR(N,N,GF2%NSPIN) ! DELTA(1,1,NW)
      COMPLEX(gq) R(N,N,GF2%NSPIN)
! LOCAL
      INTEGER ISP
! 
      DO ISP=1,GF2%NSPIN
        DR(JN,:,ISP)=DELTA*CONJG(R(:,JN,ISP)) ! Right lead
      ENDDO
      RETURN
!
      END SUBROUTINE GET_RENORM_HYBRIZ11
!
!******************************************************************************
!  Linear one-band nearest neighbour hopping; Ds=2ts
!******************************************************************************
      SUBROUTINE GET_BARE_HYBRIZ21(DELTA,OMG,DS)
      COMPLEX(gq) DELTA
      REAL(gq) DS ! Half width of the semicircular DOS = 2ts
      COMPLEX(gq) OMG
! LOCAL
      REAL(gq) SGN
      COMPLEX(gq) ZTMP,FAKT
!
      FAKT=(GF2%TSP/GF2%DS*2)**2  ! (t'/t)^2
      ZTMP=OMG**2-DS**2
      SGN=REAL(OMG,gq)*AIMAG(OMG)
      IF(SGN>=0)THEN; SGN=1._gq; ELSE; SGN=-1._gq; ENDIF
      ZTMP=SQRT(ZTMP)*SGN
      DELTA=(OMG-ZTMP)/2*FAKT
      RETURN
!
      END SUBROUTINE GET_BARE_HYBRIZ21
!
!
      END MODULE GREENFUN2

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
      MODULE WAREHOUSE
      USE gprec; USE GMPI; USE GUTIL
      USE GCONSTANT, ONLY: ZI
!
      IMPLICIT NONE
! DEFINE TYPE
      TYPE FROZEN
        INTEGER :: NSORB
        REAL(gq) :: NELECT
        INTEGER,POINTER :: IDX_ORB(:)
      END TYPE FROZEN
!
      TYPE MATRIX_BASIS
        INTEGER :: DIMHST=0,DIMHSMAX=0,DIMHST_RED
        ! Dimension of the Matrix basis for each atom
        INTEGER,POINTER :: DIM_HS(:) 
        COMPLEX(gq),POINTER :: HS(:,:,:,:)=>NULL() ! Matrix basis
        INTEGER,POINTER :: M_STRUCT(:,:,:)=>NULL() ! self-energy structure
        INTEGER,POINTER :: IDX_RED2FULL(:)=>NULL()
        INTEGER,POINTER :: M_INDEX(:,:,:,:)=>NULL()
      END TYPE MATRIX_BASIS
!
      TYPE WARE_HOUSE
        INTEGER NIONS,NTYP
        INTEGER NA2MAX,NASOMAX,NASOTOT
        LOGICAL :: LCALC_LS=.FALSE.
        COMPLEX(gq),POINTER :: C2N(:,:,:)  ! Complex Harmonics to "Natural" basis
        COMPLEX(gq),POINTER :: R2N(:,:,:),B2N(:,:,:)  ! Original basis (can be real Harmonics) to "Natural" basis
        COMPLEX(gq),POINTER :: N2N(:,:,:) ! Additional rotations
        COMPLEX(gq),POINTER :: R(:,:,:),R0(:,:,:),D(:,:,:),D0(:,:,:),LA1(:,:,:),LA2(:,:,:),Z(:,:,:)
        COMPLEX(gq),POINTER :: COPY(:,:,:)
        COMPLEX(gq),POINTER :: S_VEC(:,:,:,:)=>NULL(),L_VEC(:,:,:,:)=>NULL()
! CMR3 BEGIN
        COMPLEX(gq),POINTER :: R1(:,:,:),PFPR(:,:,:) 
! CMR3 END R1 is R, R is F(R)
        COMPLEX(gq),POINTER :: VDC2(:,:,:),NPHY_FIX(:,:,:)
        REAL(gq),POINTER :: TFC(:,:),TFF(:,:,:) ! f-c and f-f compoments
        REAL(gq) :: TCC,EF
        COMPLEX(gq),POINTER :: NKS(:,:,:),NC_VAR(:,:,:),NC_PHY(:,:,:),EL0(:,:,:)
        REAL(gq),POINTER :: NKS_COEF(:),LA1_COEF(:),LA2_COEF(:),NCV_COEF(:)
        COMPLEX(gq),POINTER :: R_COEF(:)
        COMPLEX(gq),POINTER :: ISIMIX(:,:,:) ! 1/SQRT(NKS(1-NKS))
        TYPE (FROZEN),POINTER :: FZ(:)
        ! Hermitian matrix basis, (Hermitian) maxtrix basis with 
        ! selective Mott localization.
        TYPE (MATRIX_BASIS) :: HM, HM_L, HM_R, HM_E
      END TYPE WARE_HOUSE
! VARIABLE
      TYPE (WARE_HOUSE),SAVE :: WH
!
! SUBROUTINE
      CONTAINS
!
!******************************************************************************
      SUBROUTINE ALLOC_WAREHOUSE(NA2MAX,NIONS)
      INTEGER NA2MAX,NIONS
!
      WH%NIONS=NIONS; WH%NA2MAX=NA2MAX
      ALLOCATE(WH%R     (NA2MAX,NA2MAX,NIONS), &
              &WH%R1    (NA2MAX,NA2MAX,NIONS), &
              &WH%PFPR  (NA2MAX,NA2MAX,NIONS), &
              &WH%R0    (NA2MAX,NA2MAX,NIONS), & 
              &WH%C2N   (NA2MAX,NA2MAX,NIONS), &
              &WH%R2N   (NA2MAX,NA2MAX,NIONS), &
              &WH%N2N   (NA2MAX,NA2MAX,NIONS), &
              &WH%Z     (NA2MAX,NA2MAX,NIONS), &
              &WH%D     (NA2MAX,NA2MAX,NIONS), &
              &WH%D0    (NA2MAX,NA2MAX,NIONS), &
              &WH%LA1   (NA2MAX,NA2MAX,NIONS), &
              &WH%LA2   (NA2MAX,NA2MAX,NIONS), &
              &WH%NKS   (NA2MAX,NA2MAX,NIONS), &
              &WH%ISIMIX(NA2MAX,NA2MAX,NIONS), &
              &WH%NC_VAR(NA2MAX,NA2MAX,NIONS), &
              &WH%NC_PHY(NA2MAX,NA2MAX,NIONS), &
              &WH%COPY  (NA2MAX,NA2MAX,NIONS), &
              &WH%VDC2  (NA2MAX,NA2MAX,NIONS), &
              &WH%NPHY_FIX(NA2MAX,NA2MAX,NIONS), &
              &WH%EL0   (NA2MAX,NA2MAX,NIONS))
      WH%LA1=0; WH%LA2=0; WH%EL0=0; WH%VDC2=0
      WH%NPHY_FIX(1,1,1)=2._gq
      RETURN
!
      END SUBROUTINE ALLOC_WAREHOUSE
!
!****************************************************************************
      SUBROUTINE WRT_WH_MATRIX_LIST(IU, FNAME, A)
      INTEGER,INTENT(IN)::IU
      CHARACTER(*) FNAME
      COMPLEX(gq) A(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
! LOCAL
      INTEGER NI
      CHARACTER(40) FMT
      REAL(gq) IBUF(WH%NA2MAX, WH%NA2MAX)
!
      IF(GP%MYRANK.NE.GP%MASTER) RETURN
      IBUF = 0
      WRITE(FMT,'(A,I3,A)')"(",WH%NA2MAX,"F22.16)"
      OPEN(IU,FILE=FNAME,STATUS='REPLACE')
      DO NI = 1,WH%NIONS
        WRITE(IU,FMT) REAL(A(:,:,NI))
        WRITE(IU,FMT)AIMAG(A(:,:,NI))
      ENDDO
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE WRT_WH_MATRIX_LIST
!
!****************************************************************************
      SUBROUTINE READ_WH_MATRIX_LIST(IU, FNAME, A)
      INTEGER IU
      CHARACTER(*) FNAME
      COMPLEX(gq) A(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
! LOCAL
      INTEGER NI
      REAL(gq) RBUF(WH%NA2MAX,WH%NA2MAX),IBUF(WH%NA2MAX,WH%NA2MAX)
!
      OPEN(IU,FILE=FNAME,STATUS='OLD')
      DO NI = 1,WH%NIONS
        READ(IU,*) RBUF; READ(IU,*) IBUF
        A(:,:,NI) = RBUF + ZI*IBUF
      ENDDO
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE READ_WH_MATRIX_LIST
!
!****************************************************************************
      SUBROUTINE READ_WH_MATRIX_STRUCT(IU, FNAME, MB)
      INTEGER IU
      CHARACTER(*) FNAME
      TYPE (MATRIX_BASIS) :: MB
! LOCAL
      INTEGER NI,NA2,NIP
!
      OPEN(IU,FILE=FNAME,STATUS='OLD',ERR=100)
      ALLOCATE(MB%M_STRUCT(WH%NA2MAX,WH%NA2MAX,WH%NIONS))
      DO NI = 1,WH%NIONS
        READ(IU,*) NIP, NA2
        IF(NIP/=NI)THEN
          STOP " ERROR IN READ_WH_SIGMA_STRUCT: CORRUPTED FILE WH_SIGMA_STRUCT.INP!"
        ENDIF
        READ(IU,*) MB%M_STRUCT(1:NA2,1:NA2,NI)
      ENDDO
      CLOSE(IU)
100   RETURN
!
      END SUBROUTINE READ_WH_MATRIX_STRUCT
!
!****************************************************************************
      SUBROUTINE WRT_WH_RLNEF(IU, FNAME)
      INTEGER IU
      character(*) FNAME
!
      IF(GP%MYRANK.NE.GP%MASTER) RETURN
      OPEN(IU,FILE=FNAME,STATUS='REPLACE')
      WRITE(IU,*) WH%NA2MAX,WH%NIONS
      CALL ZWRT_WH_RLNEF1(IU,WH%R)
      CALL ZWRT_WH_RLNEF1(IU,WH%LA1)
      CALL ZWRT_WH_RLNEF1(IU,WH%NKS)
      CALL ZWRT_WH_RLNEF1(IU,WH%EL0)
      WRITE(IU,*) WH%EF
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE WRT_WH_RLNEF
!
!****************************************************************************
      SUBROUTINE DWRT_WH_RLNEF1(IU,X)
      INTEGER IU
      REAL(gq) X(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
! LOCAL
      INTEGER NI
      REAL(gq) IBUF(WH%NA2MAX,WH%NA2MAX)
      CHARACTER(40) FMT
!
      IBUF=0
      WRITE(FMT,'(A,I3,A)')"(",WH%NA2MAX,"F22.16)"
      DO NI=1,WH%NIONS
        WRITE(IU,FMT) REAL(X(:,:,NI))
        WRITE(IU,FMT)IBUF
      ENDDO
      RETURN
!
      END SUBROUTINE DWRT_WH_RLNEF1
!
!****************************************************************************
      SUBROUTINE ZWRT_WH_RLNEF1(IU,X)
      INTEGER IU
      COMPLEX(gq) X(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
! LOCAL
      INTEGER NI
      CHARACTER(40) FMT
!
      WRITE(FMT,'(A,I3,A)')"(",WH%NA2MAX,"F22.16)"
      DO NI=1,WH%NIONS
        WRITE(IU,FMT) REAL(X(:,:,NI))
        WRITE(IU,FMT)AIMAG(X(:,:,NI))
      ENDDO
      RETURN
!
      END SUBROUTINE ZWRT_WH_RLNEF1
!
!****************************************************************************
      SUBROUTINE WRT_WH_RLNEF_ORIG(IU)
      INTEGER IU
! LOCAL
      INTEGER NI
      COMPLEX(gq) ZBUF(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
!
      IF(GP%MYRANK.NE.GP%MASTER) RETURN
      OPEN(IU,FILE='WH_RLNEF.OUT_ORIG',STATUS='REPLACE')
      WRITE(IU,*) WH%NA2MAX,WH%NIONS
      DO NI=1,WH%NIONS
        ZBUF(:,:,NI)=MATMUL(WH%N2N(:,:,NI),MATMUL(WH%R(:,:,NI), &
                    &TRANSPOSE(CONJG(WH%N2N(:,:,NI)))))
      ENDDO
      CALL ZWRT_WH_RLNEF1(IU,ZBUF)
      DO NI=1,WH%NIONS
        ZBUF(:,:,NI)=MATMUL(WH%N2N(:,:,NI),MATMUL(WH%LA1(:,:,NI), &
                    &TRANSPOSE(CONJG(WH%N2N(:,:,NI)))))
      ENDDO
      CALL ZWRT_WH_RLNEF1(IU,ZBUF)
      DO NI=1,WH%NIONS
        ZBUF(:,:,NI)=MATMUL(TRANSPOSE(WH%N2N(:,:,NI)),MATMUL(WH%NKS(:,:,NI), &
                    &CONJG(WH%N2N(:,:,NI))))
      ENDDO
      CALL ZWRT_WH_RLNEF1(IU,ZBUF)
      DO NI=1,WH%NIONS
        ZBUF(:,:,NI)=MATMUL(WH%N2N(:,:,NI),MATMUL(WH%EL0(:,:,NI), &
                    &TRANSPOSE(CONJG(WH%N2N(:,:,NI)))))
      ENDDO
      CALL ZWRT_WH_RLNEF1(IU,ZBUF)
      WRITE(IU,*) WH%EF
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE WRT_WH_RLNEF_ORIG
!
!****************************************************************************
      SUBROUTINE READ_WH_RLNEF(IU)
      INTEGER IU
! LOCAL
      INTEGER NA2MAX,NIONS
!
      OPEN(IU,FILE='WH_RLNEF.INP',STATUS='OLD')
      READ(IU,*)NA2MAX,NIONS
      IF((NA2MAX/=WH%NA2MAX).OR.(NIONS/=WH%NIONS))THEN
        WRITE(0,'(" ERROR IN READ_WH_RLNEF: NA2MAX=",2I3," NIONS=",2I3)')NA2MAX,WH%NA2MAX,NIONS,WH%NIONS; STOP
      ENDIF
      CALL READ_WH_RLNEF1(IU,WH%R)
      CALL READ_WH_RLNEF1(IU,WH%LA1)
      CALL READ_WH_RLNEF1(IU,WH%NKS)
      WH%COPY(1,1,1)=1.E5_gq
      CALL READ_WH_RLNEF1(IU,WH%COPY)
      READ(IU,*,ERR=100,END=100)WH%EF
100   CLOSE(IU)
      RETURN
!
      END SUBROUTINE READ_WH_RLNEF
!
!****************************************************************************
      SUBROUTINE READ_WH_RLNEF1(IU,X)
      INTEGER IU
      COMPLEX(gq) X(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
! LOCAL
      INTEGER NI
      REAL(gq) RBUF(WH%NA2MAX,WH%NA2MAX),IBUF(WH%NA2MAX,WH%NA2MAX)
!
      DO NI=1,WH%NIONS
      READ(IU,*,ERR=100,END=100)RBUF; READ(IU,*,ERR=100,END=100)IBUF
        X(:,:,NI)=RBUF+ZI*IBUF
      ENDDO
100   RETURN
!
      END SUBROUTINE READ_WH_RLNEF1
!
!****************************************************************************
! READ [S_x, S_y, S_z] and [L_x, L_y, L_z]
!****************************************************************************
      SUBROUTINE READ_WH_SL_VEC(IU)
      INTEGER IU
! LOCAL
      INTEGER I
      LOGICAL LEXIST
!
      OPEN(IU,FILE='WH_SL_VEC.INP',STATUS='OLD',ERR=100)
      WH%LCALC_LS=.TRUE.
      ALLOCATE(WH%S_VEC(WH%NA2MAX,WH%NA2MAX,3,WH%NIONS))   
      ALLOCATE(WH%L_VEC(WH%NA2MAX,WH%NA2MAX,3,WH%NIONS))   
      DO I=1,3
        CALL READ_WH_RLNEF1(IU,WH%S_VEC(:,:,I,:))
      ENDDO
      DO I=1,3
        CALL READ_WH_RLNEF1(IU,WH%L_VEC(:,:,I,:))
      ENDDO
      CLOSE(IU)
      RETURN
100   CONTINUE
      WH%LCALC_LS=.FALSE.
      RETURN
!
      END SUBROUTINE READ_WH_SL_VEC
!
!****************************************************************************
      SUBROUTINE READ_WH_HS(IU, FNAME, MB)
      INTEGER IU
      CHARACTER(*) FNAME
      TYPE (MATRIX_BASIS) :: MB
! LOCAL
      INTEGER NA2,DIM_HS,NIONS,NI,IH,IADD,DIMHSMAX
      REAL(gq) RBUF(WH%NA2MAX,WH%NA2MAX),IBUF(WH%NA2MAX,WH%NA2MAX)
!
      OPEN(IU,FILE=FNAME,STATUS='OLD')
      ALLOCATE(MB%DIM_HS(WH%NIONS))
      IADD=0
      DO WHILE(IADD<2)
        READ(IU,*)NIONS
        IF(NIONS /= WH%NIONS)THEN
          WRITE(0,'(" INCONSISTENCY ERROR: NIONS_FILE=",I5," WHILE WH%NIONS=",I5)')NIONS,WH%NIONS; STOP
        ENDIF
        DO NI=1,NIONS
          READ(IU,*)NA2,DIM_HS
          MB%DIM_HS(NI)=DIM_HS
          DO IH=1,DIM_HS
            READ(IU,*)RBUF(1:NA2,1:NA2)
            READ(IU,*)IBUF(1:NA2,1:NA2)
            IF(IADD>0) MB%HS(1:NA2,1:NA2,IH,NI)=RBUF(1:NA2,1:NA2)+ZI*IBUF(1:NA2,1:NA2)
          ENDDO
        ENDDO
        IF(IADD == 0)THEN
          DIMHSMAX = MAXVAL(MB%DIM_HS)
          ALLOCATE(MB%HS(WH%NA2MAX,WH%NA2MAX,DIMHSMAX,NIONS)); MB%HS=0
          REWIND(IU)
        ENDIF
        IADD = IADD+1
      ENDDO
      CLOSE(IU)
      MB%DIMHSMAX=MAXVAL(MB%DIM_HS)
      RETURN
!
      END SUBROUTINE READ_WH_HS
!
!****************************************************************************
      SUBROUTINE SET_WH_HS(IU, FNAME, MB, MA)
      INTEGER,INTENT(IN) :: IU
      CHARACTER(*),INTENT(IN) :: FNAME
      TYPE(MATRIX_BASIS) :: MB
      TYPE(MATRIX_BASIS),OPTIONAL :: MA
!
      LOGICAL LEXIST
!
      INQUIRE(FILE=FNAME,EXIST=LEXIST)
      IF(LEXIST)THEN
        CALL READ_WH_HS(IU, FNAME, MB)
      ELSEIF(PRESENT(MA))THEN
        MB%DIM_HS=>MA%DIM_HS
        MB%HS=>MA%HS
      ENDIF
      RETURN
!
      END SUBROUTINE SET_WH_HS
!
!****************************************************************************
! Same format as GUTZ4.INP
!****************************************************************************
      SUBROUTINE SET_WH_N2N(IU,IO)
      INTEGER IU,IO
! LOCAL
      INTEGER I,NA2MAX,NIONS
      REAL(gq) RN2N(WH%NA2MAX,WH%NA2MAX),IN2N(WH%NA2MAX,WH%NA2MAX)
!
      OPEN(IU,FILE='WH_N2N.INP',STATUS='OLD',ERR=100)
      READ(IU,*)NA2MAX,NIONS
      IF(NA2MAX/=WH%NA2MAX.OR.NIONS/=WH%NIONS)THEN
        STOP ' ERROR IN SET_WH_N2N: DIMENSION NOT MATCH!'
      ENDIF
      DO I=1,NIONS
        READ(IU,*)RN2N; READ(IU,*)IN2N
        WH%N2N(:,:,I)=RN2N+ZI*IN2N
      ENDDO
      CLOSE(IU)
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" WH_N2N.INP READ IN.")')
      RETURN
100   CONTINUE
      WH%N2N=0
      DO I=1,WH%NA2MAX; WH%N2N(I,I,:)=1._gq; ENDDO
      RETURN
!
      END SUBROUTINE SET_WH_N2N
!
!****************************************************************************
! FROZEN SPIN-ORBITAL INFORMATION
!****************************************************************************
      SUBROUTINE READ_FROZEN(IU)
      INTEGER IU
! LOCAL
      INTEGER NI,NI_
      LOGICAL LEXIST
!
      INQUIRE(FILE='FROZEN.INP',EXIST=LEXIST)
      ALLOCATE(WH%FZ(WH%NIONS))
      IF(.NOT.LEXIST)THEN
        WH%FZ(:)%NSORB=0; WH%FZ(:)%NELECT=0
        RETURN
      ENDIF
      OPEN(IU,FILE='FROZEN.INP',STATUS='OLD')
      DO NI=1,WH%NIONS
        READ(IU,*)NI_
        IF(NI/=NI_)THEN
          STOP " ERROR IN FROZEN.INP: NI/=NI_!"
        ENDIF
        READ(IU,*)WH%FZ(NI)%NSORB,WH%FZ(NI)%NELECT
        IF(WH%FZ(NI)%NSORB>0)THEN
          ALLOCATE(WH%FZ(NI)%IDX_ORB(WH%FZ(NI)%NSORB))
          READ(IU,*)WH%FZ(NI)%IDX_ORB
        ENDIF
      ENDDO
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE READ_FROZEN
!
!****************************************************************************
      SUBROUTINE OUT_FROZEN(IO)
      INTEGER IO
! LOCAL
      INTEGER NI
!
      IF(GP%MYRANK.NE.GP%MASTER) RETURN
      WRITE(IO,'(" FROZEN ORBITAL INFO:")')
      DO NI=1,WH%NIONS
        WRITE(IO,'(" NI=",I3," NSORB=",I2," NELECT=",F6.1)')NI,WH%FZ(NI)%NSORB,WH%FZ(NI)%NELECT
        IF(WH%FZ(NI)%NSORB>0)THEN
          WRITE(IO,'("    IDX_ORB=",14I3)')WH%FZ(NI)%IDX_ORB
        ENDIF
      ENDDO
      RETURN
!      
      END SUBROUTINE OUT_FROZEN
!      
!****************************************************************************
      SUBROUTINE MODIFY_R_LA1_FROZEN()
      INTEGER NI,I,I_
!
      DO NI=1,WH%NIONS; DO I=1,WH%FZ(NI)%NSORB
        I_=WH%FZ(NI)%IDX_ORB(I)
        WH%R(I_,:,NI)=0
        WH%LA1(I_,I_,NI)=30._gq
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE MODIFY_R_LA1_FROZEN
!
!****************************************************************************
! CHECK WH_HS
!****************************************************************************
      SUBROUTINE CHK_WH_HS(IO,MB)
      INTEGER,INTENT(IN)::IO
      TYPE (MATRIX_BASIS) :: MB
! LOCAL
      INTEGER NI,I1,DIMHS
      COMPLEX(gq) COEF(MB%DIMHSMAX)
!
      IF(GP%MYRANK.NE.GP%MASTER) RETURN
      WRITE(IO,'(" CHK_WH_HS WITH DIMHST = ",I5)')MB%DIMHST
      DO NI = 1, WH%NIONS
      DIMHS = MB%DIM_HS(NI)
      DO I1 = 1, DIMHS
        CALL GET_HM_EXPAND(MB%HS(:,:,I1,NI),MB%HS(:,:,:,NI),WH%NA2MAX,DIMHS,COEF,1,.FALSE.)
        WRITE(IO,'(" NI=",I3," I_HS=",I4," COEF:")')NI,I1
        WRITE(IO,'(5(3X,2F7.3))')COEF(1:DIMHS)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE CHK_WH_HS
!
!
      END MODULE WAREHOUSE

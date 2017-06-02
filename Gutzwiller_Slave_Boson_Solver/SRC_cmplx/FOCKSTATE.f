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
      MODULE FOCKSTATE
      USE gprec; USE SPARSE; USE GMPI
      IMPLICIT NONE
!
! DEFINE TYPES
      TYPE SECTOR1
        INTEGER DIM
        INTEGER,POINTER :: ID(:)       ! SECTOR START/END INDEX DIM=(N+1)
        REAL(gq),POINTER:: VAL(:)
        INTEGER,POINTER :: ID_FROZEN(:)=>NULL()
        INTEGER,POINTER :: I_BS(:)=>NULL()
      END TYPE SECTOR1
!
      TYPE FOCK_STATE
        INTEGER         :: DIM,IDIM    ! DIMENSION OF FOCK STATES (BS) AND INVERSION-BASIS (IBS)
        INTEGER         :: NA,NA2      ! TOTAL (SPIN) ORBITAL DIM
        LOGICAL         :: LFZ1SET     ! ONLY PICK ONE SET OF OCCUPIED MOTT LOCALIZED ORBITALS
        INTEGER         :: OCC(2)      ! OCCUPATION RANGE FOR FOCK STATES BASIS, OCC2(1)=MAX(OCC(1)-2,0)
        INTEGER,POINTER :: ID_FSC_N(:) ! N-Block index for C
        TYPE(SECTOR1)   :: SEC_N       ! SECTOR ACCORDING TO OCC-N
        INTEGER,POINTER :: BS(:),IBS(:)! FOCK STATE LABEL,(Full basis, INVERSE TABLE)
        INTEGER,POINTER :: BS_JSZ(:)    ! Sz of each Fock states (assuming up,dn,up,dn order)
                                       ! Or Jz assuming J-1/2 and J+1/2 block. 
        TYPE(DCOO_MATRIX),POINTER :: C(:)  ! c sparse matrix representation
      END TYPE FOCK_STATE
!
      CONTAINS
! SUBROUTINES
!****************************************************************************
      SUBROUTINE INI_FOCKSTATE(FS,N_FZ,LGPRJ,IDX_FZ,NE_FZ)
      TYPE (FOCK_STATE) :: FS
      INTEGER,INTENT(IN) :: N_FZ,LGPRJ
      INTEGER,OPTIONAL :: NE_FZ,IDX_FZ(N_FZ)
! LOCAL
      INTEGER I,I1,I2,NA2,NEU,NEL
!
      NEL=FS%OCC(1); NEU=FS%OCC(2)
      IF(NEL.LT.0.OR.NEL.GT.NEU.OR.NEU.GT.FS%NA2) STOP ' ERROR: ILLEGAL FS%OCC!'
      NA2=FS%NA2
      CALL NCHM(NA2,NEL-1,I1,.TRUE.)
      CALL NCHM(NA2,NEU  ,I2,.TRUE.)
      FS%DIM=I2-I1; FS%IDIM=2**NA2
      ALLOCATE(FS%BS(FS%DIM),FS%IBS(FS%IDIM)); FS%BS=0; FS%IBS=0
      CALL SET_FS_BS(FS)
      ALLOCATE(FS%SEC_N%ID_FROZEN(FS%SEC_N%DIM+1)); FS%SEC_N%ID_FROZEN(1)=1
      IF(N_FZ>0)THEN
        CALL MODIFY_FS_BS_FROZEN(FS,N_FZ,IDX_FZ,NE_FZ)
      ELSE
        FS%SEC_N%ID_FROZEN=FS%SEC_N%ID
        ALLOCATE(FS%SEC_N%I_BS(FS%DIM))
        DO I=1,FS%DIM; FS%SEC_N%I_BS(I)=I; ENDDO
      ENDIF
      DO I=1,FS%DIM; FS%IBS(FS%BS(I)+1)=I; ENDDO
      CALL SET_FS_CMAT(FS)
      CALL SET_FS_BS_JSZ(FS,LGPRJ)
      RETURN
!
      END SUBROUTINE INI_FOCKSTATE
!
!=============================================================================
! SET UP FOCK STATE BASIS LABEL
!=============================================================================
      SUBROUTINE SET_FS_BS(FS)
      TYPE (FOCK_STATE) FS
! LOCAL
      INTEGER I,N1,NA2,NEL,NEU,NBASE
      INTEGER,ALLOCATABLE :: IDX(:)
!
      NA2=FS%NA2
      NEL=FS%OCC(1); NEU=FS%OCC(2)
      FS%SEC_N%DIM=NEU-NEL+1
      ALLOCATE(FS%SEC_N%ID(FS%SEC_N%DIM+1),FS%SEC_N%VAL(FS%SEC_N%DIM))
      FS%SEC_N%ID=0; FS%SEC_N%VAL=0
      ALLOCATE(IDX(NEL:NEU)); IDX=0
      CALL NCHM(NA2,NEL-1,NBASE,.TRUE.)
      DO I=NEL,NEU
        CALL NCHM(NA2,I,IDX(I),.TRUE.)
      ENDDO
      IDX=IDX-NBASE
!
      FS%SEC_N%ID(1)=1
      DO I=2,FS%SEC_N%DIM+1
        FS%SEC_N%ID(I)=IDX(NEL+I-2)+1
      ENDDO
      IDX=FS%SEC_N%ID(1:FS%SEC_N%DIM)-1
      DO I=NEL,NEU
        FS%SEC_N%VAL(I-NEL+1)=REAL(I,gq)
      ENDDO
!
      DO I=0,2**NA2-1
        CALL SUM_OCC(I,0,NA2-1,N1)
        IF(N1.LT.NEL.OR.N1.GT.NEU)GOTO 100
        IDX(N1)=IDX(N1)+1
        FS%BS(IDX(N1))=I
100     CONTINUE
      ENDDO ! I
      IF(IDX(NEU).NE.FS%DIM)THEN
        WRITE(0,'(" ERROR IN SET_RFK_BS: IDX(NEU) vs FS%DIM=",2I9)')IDX(NEU),FS%DIM
        STOP
      ENDIF
      DEALLOCATE(IDX)
      RETURN
!
      END SUBROUTINE SET_FS_BS
!
!****************************************************************************
      SUBROUTINE MODIFY_FS_BS_FROZEN(FS,N_FZ,IDX_FZ,NE_FZ)
      TYPE (FOCK_STATE) FS
      INTEGER N_FZ,NE_FZ,IDX_FZ(N_FZ)
! LOCAL
      INTEGER I,N,ITMP,NOCC,NOCC_,M
      INTEGER FS_BS(FS%DIM,0:FS%SEC_N%DIM),N_BS(0:FS%SEC_N%DIM)
      INTEGER I_BS(FS%DIM)
!
      N_BS=0
      DO N=1,FS%SEC_N%DIM
      FS%SEC_N%ID_FROZEN(N+1)=FS%SEC_N%ID_FROZEN(N)
      M=FS%SEC_N%ID(N+1)-FS%SEC_N%ID(N)
      IF(FS%SEC_N%VAL(N)+0.1<NE_FZ)THEN
        FS_BS(N_BS(0)+1:N_BS(0)+M,0)=FS%BS(FS%SEC_N%ID(N):FS%SEC_N%ID(N+1)-1)
        N_BS(0)=N_BS(0)+M
      ELSE
        N_BS(N)=M
        FS_BS(1:M,N)=FS%BS(FS%SEC_N%ID(N):FS%SEC_N%ID(N+1)-1)
        DO I=FS%SEC_N%ID(N),FS%SEC_N%ID(N+1)-1
          ITMP=FS%BS(I)
          CALL SUM_OCC_LIST(ITMP,N_FZ,IDX_FZ,NOCC)
          IF(FS%LFZ1SET)THEN
            CALL SUM_OCC_LIST(ITMP,NE_FZ,IDX_FZ(1:NE_FZ),NOCC_)
            IF(NOCC_==NOCC)THEN
              NOCC=NOCC_
            ELSE
              NOCC=NE_FZ+1
            ENDIF
          ENDIF
          IF(NOCC==NE_FZ)THEN
            I_BS(FS%SEC_N%ID_FROZEN(N+1))=I
            FS%SEC_N%ID_FROZEN(N+1)=FS%SEC_N%ID_FROZEN(N+1)+1
          ENDIF
        ENDDO
      ENDIF
      ENDDO
      ALLOCATE(FS%SEC_N%I_BS(FS%SEC_N%ID_FROZEN(FS%SEC_N%DIM+1)-1))
      FS%SEC_N%I_BS=I_BS(1:FS%SEC_N%ID_FROZEN(FS%SEC_N%DIM+1)-1)
! update 
      FS%BS(1:N_BS(0))=FS_BS(1:N_BS(0),0)
      FS%SEC_N%ID(1)=N_BS(0)+1
      DO N=1,FS%SEC_N%DIM
        FS%SEC_N%ID(N+1)=FS%SEC_N%ID(N)+N_BS(N)
        FS%BS(FS%SEC_N%ID(N):FS%SEC_N%ID(N+1)-1)=FS_BS(1:N_BS(N),N)
      ENDDO
      RETURN
!
      END SUBROUTINE MODIFY_FS_BS_FROZEN
!
!****************************************************************************
      SUBROUTINE SET_FS_CMAT(FS)
      TYPE (FOCK_STATE) FS
! LOCAL
      INTEGER NA2,FSDIM,I,JF,ITMP,NSUM,JFP,ISIG,NNZ
!
      NA2=FS%NA2; FSDIM=FS%DIM
      ALLOCATE(FS%C(NA2))
      CALL SET_FSC_NNZ(FS,NNZ)
! C, COO format
      DO I=1,NA2
      CALL ALLOC_DCOO(FS%C(I),NNZ,FSDIM,FSDIM)
      NSUM=0
      DO JF=1,FSDIM
        ITMP=FS%BS(JF)
        ISIG=1
        CALL A_STAT(ITMP,I-1,.TRUE.,ISIG) ! <JF|C|JFP>
        JFP=FS%IBS(ITMP+1)
        IF(JFP.GT.0.AND.ISIG.NE.0)THEN
          NSUM=NSUM+1
          FS%C(I)%I(NSUM)=JF; FS%C(I)%J(NSUM)=JFP
          FS%C(I)%A(NSUM)=ISIG
        ENDIF
      ENDDO !JF
      IF(NSUM/=NNZ)STOP ' ERROR: CHECK C MATRIX!'
      ENDDO ! I
      RETURN
!
      END SUBROUTINE SET_FS_CMAT
!
!****************************************************************************
      SUBROUTINE SET_FS_BS_JSZ(FS,LGPRJ)
      TYPE (FOCK_STATE) FS
      INTEGER,INTENT(IN)::LGPRJ
! LOCAL
      INTEGER I,J,NA2,ITMP,LL,NBASE
!
      ALLOCATE(FS%BS_JSZ(FS%DIM)); FS%BS_JSZ=0
      NA2=FS%NA2
      DO I=1,FS%DIM
        ITMP=FS%BS(I)
        FS%BS_JSZ(I)=0
        IF(LGPRJ==15)THEN
          LL=(NA2-2)/4
          NBASE=0
          DO J=NBASE+1,NBASE+2*LL
            IF(BTEST(ITMP,J-1))FS%BS_JSZ(I)=FS%BS_JSZ(I)+2*(J-NBASE-1)-5
          ENDDO
          NBASE=NBASE+2*LL
          DO J=NBASE+1,NBASE+2*LL+2
            IF(BTEST(ITMP,J-1))FS%BS_JSZ(I)=FS%BS_JSZ(I)+2*(J-NBASE-1)-7
          ENDDO
        ELSE
          DO J=1,NA2,2
            IF(BTEST(ITMP,J-1))FS%BS_JSZ(I)=FS%BS_JSZ(I)+1
            IF(BTEST(ITMP,J  ))FS%BS_JSZ(I)=FS%BS_JSZ(I)-1
          ENDDO
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE SET_FS_BS_JSZ
!
!****************************************************************************
! <F_I | c_{IR}^+ c_IL  | F_J>
!****************************************************************************
      SUBROUTINE CALC_NIJ_FOCK(A,FS,N1,N2,IR,IL)
      TYPE(FOCK_STATE),INTENT(IN)::FS
      INTEGER N1,N2,IR,IL
      TYPE (ZCSR_MATRIX),INTENT(OUT)::A
! LOCAL
      INTEGER I,BS,IBS,SGN,NNZ,NFK
!
      NFK=N2-N1+1
      NNZ=NFK
      CALL ALLOC_ZCSR(A,NNZ,NFK,NFK)
      A%I(1)=1
      DO I=1,NFK
        A%I(I+1)=A%I(I)
        BS=FS%BS(I+N1-1); SGN=1
        CALL A_STAT(BS,IR-1,.FALSE.,SGN) ! < F_I | c_{IR}^+
        IF(SGN==0)CYCLE
        CALL A_STAT(BS,IL-1,.TRUE.,SGN) ! < F_I | c_{IR}^+ c_{IL}
        IF(SGN==0)CYCLE
        IBS=FS%IBS(BS+1)-N1+1
        A%J(A%I(I+1))=IBS
        A%A(A%I(I+1))=SGN
        A%I(I+1)=A%I(I+1)+1
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_NIJ_FOCK
!
!****************************************************************************
      SUBROUTINE SET_FSC_NNZ(FS,NNZ)
      TYPE (FOCK_STATE) FS
      INTEGER NNZ
! LOCAL
      INTEGER I,NCM,NEL,NEU,ISUM
!
      NEL=FS%OCC(1); NEU=FS%OCC(2)
      ALLOCATE(FS%ID_FSC_N(NEU-NEL+2))
      NNZ=0; ISUM=2; FS%ID_FSC_N(1:2)=1 ! According to right |F>
      DO I=NEL+1,NEU
        CALL NCHM(FS%NA2-1,I-1,NCM,.FALSE.)
        NNZ=NNZ+NCM
        ISUM=ISUM+1
        FS%ID_FSC_N(ISUM)=NNZ+1
      ENDDO
      RETURN
!
      END SUBROUTINE SET_FSC_NNZ
!
!=============================================================================
! CALC. AI+ / AI |  > ON EXIT: INTEGER
! ISIG: SIGN FOR THIS ACTION (ANTISYMMETRY)
!=============================================================================
      SUBROUTINE A_STAT(N,M,LACT,ISIG)
      INTEGER N,M,FS,I,OCC,ISIG ! M NEEDS TO BE SHIFTTED BY ONE.
      LOGICAL LACT,LBTEST
!
      FS(I)=1-2*MOD(I+20,2) ! CALC. (-1)^I, IN THE DECLARATION AREA!
      LBTEST=BTEST(N,M)
      IF(LACT.EQV.LBTEST) THEN
        ISIG=0
        RETURN
      ENDIF
!
      CALL SUM_OCC(N,0,M-1,OCC)
      ISIG=ISIG*FS(OCC)
      IF(LACT) THEN
        N=IBSET(N,M)
      ELSE
        N=IBCLR(N,M)
      ENDIF
      RETURN
!
      END SUBROUTINE A_STAT
!
!=============================================================================
! GET THE NUMBER OF OCCUPIED STATES IN RANGE OF [M1,M2]
!=============================================================================
      SUBROUTINE SUM_OCC(N,M1,M2,OCC)
      INTEGER N,M1,M2,OCC
      INTEGER I
!
      OCC=0
      DO I=M1,M2
        IF(BTEST(N,I)) OCC=OCC+1
      ENDDO
      RETURN
!
      END SUBROUTINE SUM_OCC
!
!=============================================================================
! ORBITAL INDEX LIST, one-based
!=============================================================================
      SUBROUTINE SUM_OCC_LIST(N,M,LIST,OCC)
      INTEGER N,M,OCC,LIST(M)
      INTEGER I
!
      OCC=0
      DO I=1,M
        IF(BTEST(N,LIST(I)-1)) OCC=OCC+1
      ENDDO
      RETURN
!
      END SUBROUTINE SUM_OCC_LIST
!
!=============================================================================
! CALCULATE N-CHOOSE-I, I=0,...,M  ( N )
!                                  ( M )
!=============================================================================
      SUBROUTINE NCHM(N,M,NCM1,TSUM)
      INTEGER N,M,NCM1
      LOGICAL TSUM
! LOCAL
      REAL(gq),ALLOCATABLE :: FAC(:) ! STORE LARGER VALUE(TEMP)
      INTEGER,SAVE:: IMAX=-30
      REAL(gq),SAVE:: NCM(30,0:30)
      INTEGER I,J,M0
!
      NCM1=0
      IF(N.GT.ABS(IMAX)) STOP 'INCREASE IMAX!'
      IF(M.GT.N.OR.M.LT.0.OR.N.LE.0)RETURN
! INITIALIZE ARRAY NCM
      IF(IMAX.LT.0)THEN
        IMAX=-IMAX
        ALLOCATE(FAC(0:IMAX+1))
        NCM=0; FAC(0)=1._gq
        DO I=0,IMAX
          FAC(I+1)=(I+1)*FAC(I)
        ENDDO
        DO I=1,IMAX; DO J=0,I
          NCM(I,J)=FAC(I)/FAC(J)/FAC(I-J)
        ENDDO; ENDDO
        DEALLOCATE(FAC)
      ENDIF ! NCM(1,1).EQ.0
!
      M0=M
      IF(TSUM) M0=0
      DO I=M0,M
        NCM1=NCM1+NINT(NCM(N,I))
      ENDDO
!
      END SUBROUTINE NCHM
!
!=============================================================================
      SUBROUTINE WRT_SEC1(SEC,IU)
      TYPE(SECTOR1) SEC
      INTEGER IU
!
      WRITE(IU)SEC%DIM; WRITE(IU)SEC%ID(1:SEC%DIM+1); WRITE(IU)SEC%VAL(1:SEC%DIM)
      RETURN
!
      END SUBROUTINE WRT_SEC1
!
!=============================================================================
      SUBROUTINE READ_SEC1(SEC,IU)
      TYPE(SECTOR1) SEC
      INTEGER IU
!
      READ(IU)SEC%DIM
      ALLOCATE(SEC%ID(SEC%DIM+1),SEC%VAL(SEC%DIM))
      READ(IU)SEC%ID; READ(IU)SEC%VAL
      RETURN
!
      END SUBROUTINE READ_SEC1
!
!=============================================================================
      SUBROUTINE READ_SEC1_TXT(SEC,IU)
      TYPE(SECTOR1) SEC
      INTEGER IU
!
      READ(IU,*)SEC%DIM
      ALLOCATE(SEC%ID(SEC%DIM+1),SEC%VAL(SEC%DIM))
      READ(IU,*)SEC%ID; READ(IU,*)SEC%VAL
      RETURN
!
      END SUBROUTINE READ_SEC1_TXT
!
!****************************************************************************
      SUBROUTINE LOAD_SEC1(SEC,NI,IFILE)
      TYPE(SECTOR1) SEC
      INTEGER NI,IFILE
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,0,IFILE))),STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      CALL READ_SEC1(SEC,GIU)
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE LOAD_SEC1
!
!****************************************************************************
      SUBROUTINE DUMP_SEC1(SEC,NI,IFILE)
      TYPE(SECTOR1) SEC
      INTEGER NI,IFILE
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,0,IFILE))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      CALL WRT_SEC1(SEC,GIU)
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE DUMP_SEC1
!
!=============================================================================
      SUBROUTINE SEC1_SIMPLIFY(SEC)
      TYPE(SECTOR1) SEC
! LOCAL
      INTEGER I,IADD
!
      IADD=0
      DO I=1,SEC%DIM
        IF(SEC%ID(I+1)<=SEC%ID(I))CYCLE
        IADD=IADD+1
        SEC%ID (IADD)=SEC%ID (I)
        SEC%VAL(IADD)=SEC%VAL(I)
      ENDDO
      SEC%ID(IADD+1)=SEC%ID(SEC%DIM+1)
      SEC%DIM=IADD
      RETURN
!
      END SUBROUTINE SEC1_SIMPLIFY
!
      END MODULE FOCKSTATE
!

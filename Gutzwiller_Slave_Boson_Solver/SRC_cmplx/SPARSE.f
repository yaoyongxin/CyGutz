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
      MODULE SPARSE
      USE gprec; USE gconstant; USE GMPI; USE GUTIL
      IMPLICIT NONE
! DEFINE TYPES
      TYPE DCSR_MATRIX
        INTEGER :: NROW=0,NCOL=0
        REAL(gq),POINTER :: A(:)
        INTEGER,POINTER :: I(:),J(:)
      END TYPE DCSR_MATRIX
!
      TYPE ZCSR_MATRIX 
        INTEGER :: NROW=0,NCOL=0
        COMPLEX(gq),POINTER :: A(:)
        INTEGER,POINTER :: I(:),J(:)
      END TYPE ZCSR_MATRIX
!
      TYPE ZCOO_MATRIX
        INTEGER :: NROW=0,NCOL=0,NNZ=0
        COMPLEX(gq),POINTER :: A(:)
        INTEGER,POINTER :: I(:),J(:)
      END TYPE ZCOO_MATRIX
!
      TYPE DCOO_MATRIX
        INTEGER :: NROW=0,NCOL=0,NNZ=0
        REAL(gq),POINTER :: A(:)
        INTEGER,POINTER :: I(:),J(:)
      END TYPE DCOO_MATRIX
!
      TYPE DCSR_VECTOR 
        INTEGER :: DIM=0
        REAL(gq),POINTER :: A(:)
        INTEGER,POINTER :: I(:)
      END TYPE DCSR_VECTOR
!
      TYPE DBD_SUBMAT
        INTEGER :: NROW=0,NCOL=0
        LOGICAL LZERO,LSPARSE
        REAL(gq),POINTER :: A(:,:)=>NULL()
        TYPE(DCSR_MATRIX) :: ACSR
      END TYPE DBD_SUBMAT
!
      TYPE DBD_MATRIX
        INTEGER :: NBK=0,DIM=0,JSHIFT ! JSHIFT=-1,0, or 1
        TYPE(DBD_SUBMAT),POINTER :: BK(:)=>NULL()
      END TYPE DBD_MATRIX
!
      TYPE ZBD_SUBMAT
        INTEGER :: NROW=0,NCOL=0
        LOGICAL LZERO,LSPARSE
        COMPLEX(gq),POINTER :: A(:,:)=>NULL()
        TYPE(ZCSR_MATRIX) :: ACSR
      END TYPE ZBD_SUBMAT
!
      TYPE ZBD_MATRIX
        INTEGER :: NBK=0,DIM=0,JSHIFT ! JSHIFT=-1,0, or 1
        TYPE(ZBD_SUBMAT),POINTER :: BK(:)=>NULL()
      END TYPE ZBD_MATRIX
!
      TYPE IBD_SUBMAT
        INTEGER :: NROW=0,NCOL=0
        INTEGER,POINTER :: A(:,:)=>NULL()
      END TYPE IBD_SUBMAT
!
      TYPE IBD_MATRIX
        INTEGER :: NBK=0,DIM=0
        TYPE(IBD_SUBMAT),POINTER :: BK(:)=>NULL()
      END TYPE IBD_MATRIX
!
      TYPE ISYM_BK_SUBMAT
        INTEGER :: N=0,IMAP=0
        INTEGER,POINTER :: IJ(:,:,:)
      END TYPE ISYM_BK_SUBMAT
!
      TYPE ISYM_BK_MATRIX
        INTEGER :: NBK=0,DIM=0
        TYPE(ISYM_BK_SUBMAT),POINTER :: BK(:)
      END TYPE ISYM_BK_MATRIX
!
      CONTAINS
!
! ---------------------------DSM-----------------------------------
      SUBROUTINE ALLOC_DCSR(A,NNZ,NROW,NCOL)
      INTEGER NNZ,NROW,NCOL
      TYPE(DCSR_MATRIX)::A
!
      !(*,'(" NROW=",I8," NCOL=",I8," NNZ=",I10)')NROW,NCOL,NNZ
      ALLOCATE(A%A(NNZ),A%J(NNZ),A%I(NROW+1))
      A%NROW=NROW; A%NCOL=NCOL
      A%A=0; A%J=0; A%I=0
      GMEM_SIZE=GMEM_SIZE+(REAL(SIZE(A%A),gq)*8+REAL(SIZE(A%J),gq)*4)*1.E-9_gq
      ! 'ALLOC_DCSR'
      RETURN
!
      END SUBROUTINE ALLOC_DCSR
!
!******************************************************************
      SUBROUTINE DCSR_WRT(A,IU)
      TYPE(DCSR_MATRIX)::A
      INTEGER IU
! LOCAL
      INTEGER NNZ1,NNZ2,NNZ3
!
      NNZ1=UBOUND(A%A,1); NNZ2=UBOUND(A%J,1); NNZ3=A%I(A%NROW+1)-1
      IF(NNZ1.NE.NNZ2.OR.NNZ1.NE.NNZ3)THEN
        WRITE(0,'(" ERROR: INCONSISTENT NNZ1/2/3=",3I7)')NNZ1,NNZ2,NNZ3
        STOP
      ENDIF
      WRITE(IU)A%NROW,A%NCOL
      WRITE(IU)A%I; WRITE(IU)A%A; WRITE(IU)A%J
      RETURN
!
      END SUBROUTINE DCSR_WRT
!
!******************************************************************
      SUBROUTINE ZCSR_WRT(A,IU)
      TYPE(ZCSR_MATRIX)::A
      INTEGER IU
! LOCAL
      INTEGER NNZ1,NNZ2,NNZ3
!
      NNZ1=UBOUND(A%A,1); NNZ2=UBOUND(A%J,1); NNZ3=A%I(A%NROW+1)-1
      IF(NNZ1.NE.NNZ2.OR.NNZ1.NE.NNZ3)THEN
        WRITE(0,'(" ERROR: INCONSISTENT NNZ1/2/3=",3I7)')NNZ1,NNZ2,NNZ3
        STOP
      ENDIF
      WRITE(IU)A%NROW,A%NCOL
      WRITE(IU)A%I; WRITE(IU)A%A; WRITE(IU)A%J
      RETURN
!
      END SUBROUTINE ZCSR_WRT
!
!*****************************************************************************
      SUBROUTINE DCSR_DUMP(A,NI,IA1,IA2,MODE)
      TYPE(DCSR_MATRIX)::A
      INTEGER NI,IA1,IA2,MODE
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,IA1,IA2,MODE))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      CALL DCSR_WRT(A,GIU)
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE DCSR_DUMP
!
!*****************************************************************************
      SUBROUTINE ZCSR_DUMP(A,NI,IA1,IA2,MODE)
      TYPE(ZCSR_MATRIX)::A
      INTEGER NI,IA1,IA2,MODE
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,IA1,IA2,MODE))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      CALL ZCSR_WRT(A,GIU)
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE ZCSR_DUMP
!
!*****************************************************************************
      SUBROUTINE DCSR2_DUMP(A,N,NI,MODE)
      INTEGER N,NI,MODE
      TYPE(DCSR_MATRIX)::A(N,N)
! LOCAL
      INTEGER IA1,IA2
!
      DO IA1=1,N; DO IA2=1,N
        IF(A(IA1,IA2)%NROW<=0)CYCLE
        CALL DCSR_DUMP(A(IA1,IA2),NI,IA1,IA2,MODE)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE DCSR2_DUMP
!
!*****************************************************************************
      SUBROUTINE ZCSR2_DUMP(A,N,NI,MODE)
      INTEGER N,NI,MODE
      TYPE(ZCSR_MATRIX)::A(N,N)
! LOCAL
      INTEGER IA1,IA2
!
      DO IA1=1,N; DO IA2=1,N
        IF(A(IA1,IA2)%NROW<=0)CYCLE
        CALL ZCSR_DUMP(A(IA1,IA2),NI,IA1,IA2,MODE)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE ZCSR2_DUMP
!
!******************************************************************
      SUBROUTINE DCSR_READ(A,IU)
      TYPE(DCSR_MATRIX)::A
      INTEGER IU,NNZ
!
      READ(IU)A%NROW,A%NCOL
      ALLOCATE(A%I(A%NROW+1))
      READ(IU)A%I
      NNZ=A%I(A%NROW+1)-1; ALLOCATE(A%J(NNZ),A%A(NNZ))
      READ(IU)A%A; READ(IU)A%J
      GMEM_SIZE=GMEM_SIZE+(REAL(SIZE(A%A),gq)*8+REAL(SIZE(A%J),gq)*4)*1.E-9_gq
      ! 'DCSR_READ'
      RETURN
!
      END SUBROUTINE DCSR_READ
!
!******************************************************************
      SUBROUTINE ZCSR_READ(A,IU)
      TYPE(ZCSR_MATRIX)::A
      INTEGER IU,NNZ
!
      READ(IU)A%NROW,A%NCOL
      ALLOCATE(A%I(A%NROW+1))
      READ(IU)A%I
      NNZ=A%I(A%NROW+1)-1; ALLOCATE(A%J(NNZ),A%A(NNZ))
      READ(IU)A%A; READ(IU)A%J
      GMEM_SIZE=GMEM_SIZE+(REAL(SIZE(A%A),gq)*16+REAL(SIZE(A%J),gq)*4)*1.E-9_gq
      ! 'ZCSR_READ'
      RETURN
!
      END SUBROUTINE ZCSR_READ
!
!*****************************************************************************
      SUBROUTINE DCSR_LOAD(A,NI,IA1,IA2,MODE)
      TYPE(DCSR_MATRIX)::A
      INTEGER NI,IA1,IA2,MODE
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,IA1,IA2,MODE))),STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      CALL DCSR_READ(A,GIU)
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE DCSR_LOAD
!
!*****************************************************************************
      SUBROUTINE DCSR2_LOAD(A,N,NI,MODE)
      INTEGER N,NI,MODE
      TYPE(DCSR_MATRIX)::A(N,N)
! LOCAL
      INTEGER IA1,IA2
!
      DO IA1=1,N; DO IA2=1,N
        IF(A(IA1,IA2)%NROW<=0)CYCLE
        CALL DCSR_LOAD(A(IA1,IA2),NI,IA1,IA2,MODE)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE DCSR2_LOAD
!
!*****************************************************************************
      SUBROUTINE ZCSR_LOAD(A,NI,IA1,IA2,MODE)
      TYPE(ZCSR_MATRIX)::A
      INTEGER NI,IA1,IA2,MODE
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,IA1,IA2,MODE))),STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      CALL ZCSR_READ(A,GIU)
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE ZCSR_LOAD
!
!*****************************************************************************
      SUBROUTINE ZCSR2_LOAD(A,N,NI,MODE)
      INTEGER N,NI,MODE
      TYPE(ZCSR_MATRIX)::A(N,N)
! LOCAL
      INTEGER IA1,IA2
!
      DO IA1=1,N; DO IA2=1,N
        IF(A(IA1,IA2)%NROW<=0)CYCLE
        CALL ZCSR_LOAD(A(IA1,IA2),NI,IA1,IA2,MODE)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE ZCSR2_LOAD
!
!******************************************************************
      SUBROUTINE WRT_DCOO(A,IU)
      TYPE(DCOO_MATRIX)::A
      INTEGER IU
!
      WRITE(IU)A%NROW,A%NCOL,A%NNZ
      WRITE(IU)A%I; WRITE(IU)A%J; WRITE(IU)A%A
      RETURN
!
      END SUBROUTINE WRT_DCOO
!
!******************************************************************
      SUBROUTINE READ_DCOO(A,IU)
      TYPE(DCOO_MATRIX)::A
      INTEGER IU
!
      READ(IU)A%NROW,A%NCOL,A%NNZ
      ALLOCATE(A%I(A%NNZ),A%J(A%NNZ),A%A(A%NNZ))
      READ(IU)A%I; READ(IU)A%J; READ(IU)A%A
      GMEM_SIZE=GMEM_SIZE+REAL(A%NNZ,gq)*24*1.E-9_gq
      ! 'READ_DCOO'
      RETURN
!
      END SUBROUTINE READ_DCOO
!
!******************************************************************
      SUBROUTINE READ_DCOO_TXT(A,IU)
      TYPE(DCOO_MATRIX)::A
      INTEGER IU
!
      READ(IU,*)A%NROW,A%NCOL,A%NNZ
      ALLOCATE(A%I(A%NNZ),A%J(A%NNZ),A%A(A%NNZ))
      READ(IU,*)A%I; READ(IU,*)A%J; READ(IU,*)A%A
      GMEM_SIZE=GMEM_SIZE+REAL(A%NNZ,gq)*24*1.E-9_gq
      ! 'READ_DCOO_TXT'
      RETURN
!
      END SUBROUTINE READ_DCOO_TXT
!
!******************************************************************
      SUBROUTINE ALLOC_DCSR_DIAG(A,N)
      INTEGER N
      TYPE(DCSR_MATRIX)::A
! LOCAL
      INTEGER I
!
      ALLOCATE(A%A(N),A%J(N),A%I(N+1))
      A%NROW=N; A%NCOL=N; A%A=0
      DO I=1,N; A%J(I)=I; A%I(I)=I; ENDDO
      A%I(N+1)=N+1
      RETURN
!
      END SUBROUTINE ALLOC_DCSR_DIAG
!
!******************************************************************
      SUBROUTINE DEALLOC_DCSR(A)
      TYPE(DCSR_MATRIX)::A
!
      GMEM_SIZE=GMEM_SIZE-(REAL(SIZE(A%A),gq)*8+REAL(SIZE(A%J),gq)*4)*1.E-9_gq
      DEALLOCATE(A%A,A%J,A%I); A%NROW=0; A%NCOL=0 ! DEALLOCATE IF MEMORY IS NOT ASSOCIATED WITH OTHERS
      ! 'DEALLOC_DCSR'
      RETURN
!
      END SUBROUTINE DEALLOC_DCSR
!
!******************************************************************
      SUBROUTINE ALLOC_DSV(V,N)
      INTEGER N
      TYPE(DCSR_VECTOR)::V
!
      ALLOCATE(V%A(N),V%I(N))
      V%DIM=N; V%A=0; V%I=0
      RETURN
!
      END SUBROUTINE ALLOC_DSV
!
!******************************************************************
      SUBROUTINE DEALLOC_DSV(V)
      TYPE(DCSR_VECTOR)::V
!
      DEALLOCATE(V%A,V%I); V%DIM=0
      RETURN
!
      END SUBROUTINE DEALLOC_DSV
!
!******************************************************************
      SUBROUTINE DCSR_SIZE(A,SIZEA)
      TYPE(DCSR_MATRIX)::A
      REAL(gq) SIZEA
! LOCAL
      INTEGER NNZ
!
      NNZ=A%I(A%NROW+1)-1
      SIZEA=REAL(NNZ,gq)*8._gq*2*1.E-9_gq
      RETURN
!
      END SUBROUTINE DCSR_SIZE
!
!******************************************************************
      FUNCTION DCSR_MAX_OFF_DIAG(A)
      TYPE(DCSR_MATRIX)::A
      REAL(gq) DCSR_MAX_OFF_DIAG
! LOCAL
      INTEGER IROW,I,ICOL
!
      DCSR_MAX_OFF_DIAG=0
      DO IROW=1,A%NROW; DO I=A%I(IROW),A%I(IROW+1)-1
        ICOL=A%J(I)
        IF(IROW.EQ.ICOL)CYCLE
        IF(ABS(A%A(I)).LE.ABS(DCSR_MAX_OFF_DIAG))CYCLE
       DCSR_MAX_OFF_DIAG=A%A(I)
      ENDDO; ENDDO
      RETURN
!
      END FUNCTION DCSR_MAX_OFF_DIAG
!
!*****************************************************************************
      SUBROUTINE SM_DCSRDCSC(A,AT)
      TYPE(DCSR_MATRIX)::A
      TYPE(DCSR_MATRIX),OPTIONAL::AT
! LOCAL
      INTEGER NROW,NCOL,NNZ
      TYPE(DCSR_MATRIX)::B
!
      NROW=A%NROW; NCOL=A%NCOL
      NNZ=A%I(NROW+1)-1
      CALL ALLOC_DCSR(B,NNZ,NCOL,NROW) ! Transpose
      CALL DCSRCSC2(NROW,NCOL,1,1,A%A,A%J,A%I,B%A,B%J,B%I)
      IF(PRESENT(AT))THEN; CALL DCSR_COPY(B,AT)
      ELSE; CALL DEALLOC_DCSR(A); CALL DCSR_COPY(B,A)
      ENDIF
      CALL DEALLOC_DCSR(B)
      RETURN
!
      END SUBROUTINE SM_DCSRDCSC
!
!*****************************************************************************
! Take out sub-matrix (N1:N2,M1:M2) of a BLOCK DIAGONAL CSR sparse matrix
!*****************************************************************************
      SUBROUTINE DCSR_DIA_SUBMAT(A,ABLK,N1,N2,M1,M2)
      INTEGER N1,M1,N2,M2
      TYPE(DCSR_MATRIX)::A,ABLK
! LOCAL
      INTEGER NROW,NCOL,NNZ
!
      NROW=N2-N1+1; NCOL=M2-M1+1
      NNZ=A%I(N2+1)-A%I(N1)
      CALL ALLOC_DCSR(ABLK,NNZ,NROW,NCOL)
      ABLK%A(1:NNZ)=A%A(A%I(N1):A%I(N2+1)-1)
      ABLK%J(1:NNZ)=A%J(A%I(N1):A%I(N2+1)-1)-N1+1
      ABLK%I(1:NROW+1)=A%I(N1:N2+1)-A%I(N1)+1
      IF(MINVAL(ABLK%J).LE.0.OR.MAXVAL(ABLK%J).GT.NCOL)THEN
        STOP 'FETAL ERROR IN DCSR_DIA_SUBMAT: CHECK A IN DCSR_SUBMAT--NOT BLOCK DIAGONAL!'
      ENDIF
      RETURN
!
      END SUBROUTINE DCSR_DIA_SUBMAT
!
!*****************************************************************************
      SUBROUTINE DCSR_SUBMAT(A,N1,N2,M1,M2,ABLK)
      INTEGER N1,M1,N2,M2
      TYPE(DCSR_MATRIX)::A
      TYPE(DCSR_MATRIX),OPTIONAL::ABLK
! LOCAL
      INTEGER NROW,NCOL,NNZ
      TYPE(DCSR_MATRIX)::BUF
!
      NROW=N2-N1+1; NCOL=M2-M1+1
      NNZ=A%I(N2+1)-A%I(N1)
      CALL ALLOC_DCSR(BUF,NNZ,NROW,NCOL)
      CALL DSUBMAT(A%NROW,1,N1,N2,M1,M2,A%A,A%J,A%I,NROW,NCOL,BUF%A,BUF%J,BUF%I)
      IF(PRESENT(ABLK))THEN
        CALL DCSR_COPY(BUF,ABLK)
      ELSE
        CALL DEALLOC_DCSR(A)
        CALL DCSR_COPY(BUF,A)
      ENDIF
      CALL DEALLOC_DCSR(BUF)
      RETURN
!
      END SUBROUTINE DCSR_SUBMAT
!
!*****************************************************************************
      SUBROUTINE DCSR_COPY(A,B)
      TYPE(DCSR_MATRIX)::A,B
! LOCAL
      INTEGER NROW,NCOL,NNZ
!
      NROW=A%NROW; NCOL=A%NCOL; NNZ=A%I(NROW+1)-1
      CALL ALLOC_DCSR(B,NNZ,NROW,NCOL)
      CALL DCSR_COPY1(A,B)
      RETURN
!
      END SUBROUTINE DCSR_COPY
!
!*****************************************************************************
      SUBROUTINE DCSR_COPY1(A,B)
      TYPE(DCSR_MATRIX)::A,B
! LOCAL
      INTEGER NROW,NNZ
!
      NROW=A%NROW; NNZ=A%I(NROW+1)-1
      B%NROW=A%NROW; B%NCOL=A%NCOL
      B%A=A%A(1:NNZ); B%J=A%J(1:NNZ); B%I=A%I
      RETURN
!
      END SUBROUTINE DCSR_COPY1
!
!*****************************************************************************
      SUBROUTINE DCSR_COMPACT(A)
      TYPE(DCSR_MATRIX)::A
! LOCAL
      TYPE(DCSR_MATRIX)::BUF
!
      CALL DCSR_COPY(A,BUF)
      CALL DEALLOC_DCSR(A)
      CALL DCSR_COPY(BUF,A)
      CALL DEALLOC_DCSR(BUF)
      RETURN
!
      END SUBROUTINE DCSR_COMPACT
!
!*****************************************************************************
      SUBROUTINE GET_ISYM_BK_MATRIX(IJ,IJBK,N)
      INTEGER N,IJ(N,N)
      TYPE(ISYM_BK_MATRIX) IJBK
! LOCAL
      INTEGER IMAX,I,J,K,JD,KD,NN,ISUM,IJ_BK(N,N,N),IJ_BK1(2,N,N)
      LOGICAL LDONE(N)
!
      IJ_BK=0; ISUM=0; LDONE=.FALSE.
      DO I=1,N
      IF(LDONE(I))CYCLE
      ISUM=ISUM+1
      DO J=I,N
      IF(IJ(J,I)==0)CYCLE
      IF(LDONE(J))THEN
        STOP ' ERROR-I IN GET_ISYM_BK_MATRIX!'
      ELSE
        LDONE(J)=.TRUE.
      ENDIF
      DO K=I,N
      IF(IJ(J,K)==0)CYCLE
      IJ_BK(J,K,ISUM)=IJ(J,K)
      ENDDO; ENDDO; ENDDO
!
      IJBK%NBK=ISUM; IJBK%DIM=N
      ALLOCATE(IJBK%BK(ISUM))
      DO I=1,IJBK%NBK
      IJ_BK1=0; JD=0
      DO J=1,N
      IF(SUM(IJ_BK(J,:,I))==0)THEN
        JD=JD+1; CYCLE
      ENDIF
      KD=0
      DO K=1,N
      IF(SUM(IJ_BK(:,K,I))==0)THEN
        KD=KD+1; CYCLE
      ENDIF
      IJ_BK1(1,J-JD,K-KD)=J; IJ_BK1(2,J-JD,K-KD)=K
      ENDDO ! K
      ENDDO ! J
      NN=N-JD
      IJBK%BK(I)%N=NN; ALLOCATE(IJBK%BK(I)%IJ(2,NN,NN))
      IJBK%BK(I)%IJ=IJ_BK1(:,1:NN,1:NN)
      IJBK%BK(I)%IMAP=I
      IMAX=MAXVAL(IJ_BK(:,:,I))
      DO J=1,I-1
        IF(MAXVAL(IJ_BK(:,:,J))==IMAX)THEN
          IJBK%BK(I)%IMAP=J
        ENDIF
      ENDDO ! J
      ENDDO ! I
!
      IF(SUM(IJBK%BK(:)%N)/=IJBK%DIM)THEN
        STOP ' ERROR-II IN GET_ISYM_BK_MATRIX!'
      ENDIF
      RETURN
!
      END SUBROUTINE GET_ISYM_BK_MATRIX
!
!*****************************************************************************
      SUBROUTINE DSYM_BK_DIAG(A,W,N,IJBK)
      INTEGER N
      TYPE(ISYM_BK_MATRIX) IJBK
      REAL(gq) A(N,N),W(N)
! LOCAL
      INTEGER IBK,I,J,NN,IBASE,IMAP
      REAL(gq),ALLOCATABLE :: ABK(:,:,:),W_(:,:)
!
      NN=MAXVAL(IJBK%BK(:)%N)
      ALLOCATE(ABK(NN,NN,IJBK%NBK),W_(NN,IJBK%NBK)); ABK=0; W_=0
      ABK=0; IBASE=0
      DO IBK=1,IJBK%NBK
      NN=IJBK%BK(IBK)%N
      IMAP=IJBK%BK(IBK)%IMAP
      IF(IMAP==IBK)THEN
        DO I=1,NN; DO J=1,NN
          ABK(I,J,IBK)=A(IJBK%BK(IBK)%IJ(1,I,J),IJBK%BK(IBK)%IJ(2,I,J))
        ENDDO; ENDDO ! I,J
        CALL HERMEV('V','L',ABK(1:NN,1:NN,IBK),W_(1:NN,IBK),NN)
      ENDIF
      W(IBASE+1:IBASE+NN)=W_(1:NN,IMAP)
      IBASE=IBASE+NN
      DO I=1,NN; DO J=1,NN
        A(IJBK%BK(IBK)%IJ(1,I,J),IJBK%BK(IBK)%IJ(2,I,J))=ABK(I,J,IMAP)
      ENDDO; ENDDO ! I,J
      ENDDO ! IBK
      DEALLOCATE(ABK,W_)
      RETURN
!
      END SUBROUTINE DSYM_BK_DIAG
!
!*****************************************************************************
      SUBROUTINE ZSYM_BK_DIAG(A,W,N,IJBK)
      INTEGER N
      TYPE(ISYM_BK_MATRIX) IJBK
      COMPLEX(gq) A(N,N)
      REAL(gq) W(N)
! LOCAL
      INTEGER IBK,I,J,NN,IBASE,IMAP
      COMPLEX(gq),ALLOCATABLE :: ABK(:,:,:)
      REAL(gq),ALLOCATABLE :: W_(:,:)
!
      NN=MAXVAL(IJBK%BK(:)%N)
      ALLOCATE(ABK(NN,NN,IJBK%NBK),W_(NN,IJBK%NBK)); ABK=0; W_=0
      ABK=0; IBASE=0
      DO IBK=1,IJBK%NBK
      NN=IJBK%BK(IBK)%N
      IMAP=IJBK%BK(IBK)%IMAP
      IF(IMAP==IBK)THEN
        DO I=1,NN; DO J=1,NN
          ABK(I,J,IBK)=A(IJBK%BK(IBK)%IJ(1,I,J),IJBK%BK(IBK)%IJ(2,I,J))
        ENDDO; ENDDO ! I,J
        CALL HERMEV('V','L',ABK(1:NN,1:NN,IBK),W_(1:NN,IBK),NN)
      ENDIF
      W(IBASE+1:IBASE+NN)=W_(1:NN,IMAP)
      IBASE=IBASE+NN
      DO I=1,NN; DO J=1,NN
        A(IJBK%BK(IBK)%IJ(1,I,J),IJBK%BK(IBK)%IJ(2,I,J))=ABK(I,J,IMAP)
      ENDDO; ENDDO ! I,J
      ENDDO ! IBK
      DEALLOCATE(ABK,W_)
      RETURN
!
      END SUBROUTINE ZSYM_BK_DIAG
!
!*****************************************************************************
      SUBROUTINE ZCSR_COMPACT(A)
      TYPE(ZCSR_MATRIX)::A
! LOCAL
      TYPE(ZCSR_MATRIX)::BUF
!
      CALL ZCSR_COPY(A,BUF)
      CALL DEALLOC_ZCSR(A)
      CALL ZCSR_COPY(BUF,A)
      CALL DEALLOC_ZCSR(BUF)
      RETURN
!
      END SUBROUTINE ZCSR_COMPACT
!
!*********************************************************************
! Direct sum of two csr matrices to form a block diagonal matrix
!*********************************************************************
      SUBROUTINE DCSR_DSUM_AB(A,B,AB,LFIRST)
      TYPE(DCSR_MATRIX)::A,B
      TYPE(DCSR_MATRIX),OPTIONAL::AB
      LOGICAL,OPTIONAL::LFIRST
! LOCAL
      INTEGER NROW1,NROW2,NCOL1,NCOL2,NNZ1,NNZ2
      INTEGER NROW,NCOL,NNZ
      TYPE(DCSR_MATRIX)::C
!
      IF(PRESENT(LFIRST))THEN
        IF(LFIRST)THEN
          IF(PRESENT(AB))THEN
            CALL DCSR_COPY(B,AB)
          ELSE
            CALL DCSR_COPY(B,A)
          ENDIF
          RETURN
        ENDIF
      ENDIF
      NROW1=A%NROW; NCOL1=A%NCOL; NNZ1=A%I(NROW1+1)-1
      NROW2=B%NROW; NCOL2=B%NCOL; NNZ2=B%I(NROW2+1)-1
      NROW=NROW1+NROW2; NCOL=NCOL1+NCOL2; NNZ=NNZ1+NNZ2
      CALL ALLOC_DCSR(C,NNZ,NROW,NCOL)
      C%A(1:NNZ1) =A%A; C%J(1:NNZ1) =A%J      ; C%I(1:NROW1) =A%I(1:NROW1)
      C%A(1+NNZ1:)=B%A; C%J(1+NNZ1:)=B%J+NCOL1; C%I(1+NROW1:)=B%I+NNZ1
      IF(PRESENT(AB))THEN; CALL DCSR_COPY(C,AB)
      ELSE; CALL DEALLOC_DCSR(A); CALL DCSR_COPY(C,A)
      ENDIF
      CALL DEALLOC_DCSR(C)
      RETURN
!
      END SUBROUTINE DCSR_DSUM_AB
!
!******************************************************************
      SUBROUTINE DCSR_DSUM_DNS(ACSR,ADNS,N,M,UPLO,LFIRST)
      TYPE(DCSR_MATRIX)::ACSR
      INTEGER N,M
      REAL(gq)ADNS(N,M)
      CHARACTER*1 UPLO
      LOGICAL LFIRST
! LOCAL
      TYPE(DCSR_MATRIX)::BUF
!
      IF(LFIRST)THEN
        CALL SM_DDNSDCSR(ADNS,N,M,ACSR,UPLO)
      ELSE
        CALL SM_DDNSDCSR(ADNS,N,M,BUF,UPLO)
        CALL DCSR_DSUM_AB(ACSR,BUF)
        CALL DEALLOC_DCSR(BUF)
      ENDIF
      RETURN
!
      END SUBROUTINE DCSR_DSUM_DNS
!
!*********************************************************************
      SUBROUTINE DCSR_AMUB_SK(TRANS,A,B,C,LB)
      CHARACTER*1 TRANS
      TYPE(DCSR_MATRIX)::A,B
      TYPE(DCSR_MATRIX),OPTIONAL::C
      INTEGER,OPTIONAL::LB
! LOCAL
      INTEGER NNZ,IERR,NROW1,NCOL1,NROW2,NCOL2
      INTEGER,ALLOCATABLE::IW(:)
      TYPE(DCSR_MATRIX)::A_,D
!
      IF(TRANS.EQ.'N'.OR.TRANS.EQ.'n')THEN
        CALL DCSR_COPY(A,A_)
      ELSE
        CALL SM_DCSRDCSC(A,AT=A_)
      ENDIF
      NROW1=A_%NROW; NCOL1=A_%NCOL
      NROW2=B %NROW; NCOL2=B %NCOL
      ALLOCATE(IW(NCOL1)); IW=0
      IF(NCOL1.NE.NROW2)THEN
        WRITE(0,'(" NCOL1 vs NROW2:",2I6)')NCOL1,NROW2
        WRITE(0,'(" NROW1 vs NCOL2:",2I6)')NROW1,NCOL2
        STOP' ERROR: NCOL1.NE.NROW2 IN DCSR_AMUB_SK!'
      ENDIF
      NNZ=MIN(MAXNNZ,NROW1*NCOL2)
      CALL ALLOC_DCSR(D,NNZ,NROW1,NCOL2)
      IERR=0
      CALL DAMUB(NROW1,NCOL1,1,A_%A,A_%J,A_%I,B%A,B%J,B%I,D%A,D%J,D%I,NNZ,IW,IERR)
      DEALLOCATE(IW)
      CALL DEALLOC_DCSR(A_)
      IF(IERR.NE.0)STOP' ERROR IN DCSR_AMUB_SK!'
      IF(PRESENT(C))THEN
        CALL DCSR_COPY(D,C)
      ELSEIF(PRESENT(LB))THEN
        CALL DEALLOC_DCSR(B)
        CALL DCSR_COPY(D,B)
      ELSE
        CALL DEALLOC_DCSR(A)
        CALL DCSR_COPY(D,A)
      ENDIF
      CALL DEALLOC_DCSR(D)
      RETURN
!
      END SUBROUTINE DCSR_AMUB_SK
!
!*********************************************************************
      SUBROUTINE DCSRMUDEN_SK(ACSR,B,AB,NC)
      TYPE(DCSR_MATRIX)::ACSR
      INTEGER NC
      REAL(gq) B(ACSR%NCOL,NC),AB(ACSR%NROW,NC)
! LOCAL
      INTEGER I,J,JL
!
      AB=0
      DO I=1,ACSR%NROW; DO J=ACSR%I(I),ACSR%I(I+1)-1
        JL=ACSR%J(J)
        CALL DAXPY(NC,ACSR%A(J),B(JL,:),1,AB(I,:),1)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE DCSRMUDEN_SK
!
!*********************************************************************
      SUBROUTINE ZCSRMUDEN_SK(ACSR,B,AB,NC)
      TYPE(ZCSR_MATRIX)::ACSR
      INTEGER NC
      COMPLEX(gq) B(ACSR%NCOL,NC),AB(ACSR%NROW,NC)
! LOCAL
      INTEGER I,J,JL
!
      AB=0
      DO I=1,ACSR%NROW; DO J=ACSR%I(I),ACSR%I(I+1)-1
        JL=ACSR%J(J)
        CALL ZAXPY(NC,ACSR%A(J),B(JL,:),1,AB(I,:),1)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE ZCSRMUDEN_SK
!
!*********************************************************************
      SUBROUTINE DCSR_APLSB_SK(TRANS,A,S,B,C)
      CHARACTER*1 TRANS
      REAL(gq) S
      TYPE(DCSR_MATRIX)::A,B
      TYPE(DCSR_MATRIX),OPTIONAL::C
! LOCAL
      INTEGER NNZ,IERR,NNZ1,NNZ2,NROW1,NCOL1,NROW2,NCOL2
      INTEGER,ALLOCATABLE::IW(:)
      TYPE(DCSR_MATRIX)::B_,D
!
      IF(TRANS.EQ.'N'.OR.TRANS.EQ.'n')THEN
        CALL DCSR_COPY(B,B_)
      ELSE
        CALL SM_DCSRDCSC(B,AT=B_)
      ENDIF
      NROW1=A %NROW; NCOL1=A %NCOL
      NROW2=B_%NROW; NCOL2=B_%NCOL
      IF(NROW1.NE.NROW2)STOP' ERROR: NROW1.NE.NROW2 IN DCSR_APLSB_SK!'
      IF(NCOL1.NE.NCOL2)STOP' ERROR: NCOL1.NE.NCOL2 IN DCSR_APLSB_SK!'
      ALLOCATE(IW(NCOL1)); IW=0
      NNZ1=A%I(NROW1+1)-1; NNZ2=B_%I(NROW2+1)-1
      NNZ=MIN(MAXNNZ,NNZ1+NNZ2)
      CALL ALLOC_DCSR(D,NNZ,NROW1,NCOL1)
      IERR=0
      CALL DAPLSB(NROW1,NCOL1,A%A,A%J,A%I,S,B_%A,B_%J,B_%I,D%A,D%J,D%I,NNZ,IW,IERR)
      IF(IERR.NE.0)STOP' ERROR IN DCSR_APLSB_SK!'
      CALL DEALLOC_DCSR(B_)
      DEALLOCATE(IW)
      IF(PRESENT(C))THEN
        CALL DCSR_COPY(D,C)
      ELSE
        CALL DEALLOC_DCSR(A)
        CALL DCSR_COPY(D,A)
      ENDIF
      CALL DEALLOC_DCSR(D)
      RETURN
!
      END SUBROUTINE DCSR_APLSB_SK
!
!*****************************************************************************
      SUBROUTINE DBM_APLSB(A,S,B)
      TYPE(DBD_MATRIX),INTENT(INOUT)::A
      REAL(gq),INTENT(IN)::S
      TYPE(DBD_MATRIX),INTENT(IN)::B
! LOCAL
      INTEGER I
!
      DO I=1,A%NBK
        IF(A%BK(I)%LZERO)CYCLE
        IF(A%BK(I)%LSPARSE)THEN
          CALL DCSR_APLSB_SK('N',A%BK(I)%ACSR,S,B%BK(I)%ACSR)
        ELSE
          A%BK(I)%A=A%BK(I)%A+S*B%BK(I)%A
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE DBM_APLSB
!
!*********************************************************************
      SUBROUTINE DCSR_GESAMUX_SK(TRANS,S,A,X,AX)
      CHARACTER*1 TRANS
      TYPE(DCSR_MATRIX)::A
      REAL(gq) S,X(*),AX(*)
! LOCAL
      INTEGER N2
      REAL(gq),ALLOCATABLE::AX_(:)
!
      IF(TRANS=='T'.OR.TRANS=='t'.OR.TRANS=='C'.OR.TRANS=='c')THEN
        N2=A%NCOL
      ELSE
        N2=A%NROW
      ENDIF
      ALLOCATE(AX_(N2)); AX_=0
      CALL DCSR_GEAMUX_SK(TRANS,A,X,AX_)
      CALL DAXPY(N2,S,AX_,1,AX,1)
      DEALLOCATE(AX_)
      RETURN
!
      END SUBROUTINE DCSR_GESAMUX_SK
!
!*********************************************************************
      SUBROUTINE ZCSR_GESAMUX_SK(TRANS,S,A,X,AX)
      CHARACTER*1 TRANS
      TYPE(ZCSR_MATRIX)::A
      COMPLEX(gq) S,X(*),AX(*)
! LOCAL
      INTEGER N2
      COMPLEX(gq),ALLOCATABLE::AX_(:)
!
      IF(TRANS=='T'.OR.TRANS=='t'.OR.TRANS=='C'.OR.TRANS=='c')THEN
        N2=A%NCOL
      ELSE
        N2=A%NROW
      ENDIF
      ALLOCATE(AX_(N2)); AX_=0
      CALL ZCSR_GEAMUX_SK(TRANS,A,X,AX_)
      CALL ZAXPY(N2,S,AX_,1,AX,1)
      DEALLOCATE(AX_)
      RETURN
!
      END SUBROUTINE ZCSR_GESAMUX_SK
!
!*********************************************************************
      SUBROUTINE DCSR_GEAMUX_SK(TRANS,A,X,AX)
      CHARACTER*1 TRANS
      TYPE(DCSR_MATRIX)::A
      REAL(gq) X(*),AX(*)
!
      IF(TRANS=='T'.OR.TRANS=='t'.OR.TRANS=='C'.OR.TRANS=='c')THEN
        CALL DATMUXR(A%NCOL,A%NROW,X,AX,A%A,A%J,A%I)
      ELSE
        CALL DAMUX(A%NROW,X,AX,A%A,A%J,A%I)
      ENDIF
      RETURN
!
      END SUBROUTINE DCSR_GEAMUX_SK
!
!*********************************************************************
      SUBROUTINE ZCSR_GEAMUX_SK(TRANS,A,X,AX)
      CHARACTER*1 TRANS
      TYPE(ZCSR_MATRIX)::A
      COMPLEX(gq) X(*),AX(*)
!
      IF(TRANS=='T'.OR.TRANS=='t')THEN
        CALL ZATMUXR(A%NCOL,A%NROW,X,AX,A%A,A%J,A%I)
      ELSEIF(TRANS=='C'.OR.TRANS=='c')THEN
        CALL ZAHMUXR(A%NCOL,A%NROW,X,AX,A%A,A%J,A%I)
      ELSE
        CALL ZAMUX(A%NROW,X,AX,A%A,A%J,A%I)
      ENDIF
      RETURN
!
      END SUBROUTINE ZCSR_GEAMUX_SK
!
!*********************************************************************
      SUBROUTINE DCSR_SYAMUX_SK(UPLO,A,X,AX)
      CHARACTER*1 UPLO
      TYPE(DCSR_MATRIX)::A
      REAL(gq) X(A%NCOL),AX(A%NROW)
! LOCAL
      INTEGER I,K,J
!
      AX=0
      DO I=1,A%NROW; DO K=A%I(I),A%I(I+1)-1
          J=A%J(K)
          AX(I)=AX(I)+A%A(K)*X(J)
          IF(J.EQ.I)CYCLE
          IF(J.GT.I.AND.(UPLO=='U'.OR.UPLO=='u'))THEN
            AX(J)=AX(J)+A%A(K)*X(I)
          ELSEIF(J.LT.I.AND.(UPLO=='L'.OR.UPLO=='l'))THEN
            AX(J)=AX(J)+A%A(K)*X(I)
          ENDIF
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE DCSR_SYAMUX_SK
!
!*********************************************************************
      SUBROUTINE ZCSR_SYAMUX_SK(UPLO,A,X,AX)
      CHARACTER*1 UPLO
      TYPE(ZCSR_MATRIX)::A
      COMPLEX(gq) X(A%NCOL),AX(A%NROW)
! LOCAL
      INTEGER I,K,J
!
      AX=0
      DO I=1,A%NROW; DO K=A%I(I),A%I(I+1)-1
          J=A%J(K)
          AX(I)=AX(I)+A%A(K)*X(J)
          IF(J.EQ.I)CYCLE
          IF(J.GT.I.AND.(UPLO=='U'.OR.UPLO=='u'))THEN
            AX(J)=AX(J)+CONJG(A%A(K))*X(I)
          ELSEIF(J.LT.I.AND.(UPLO=='L'.OR.UPLO=='l'))THEN
            AX(J)=AX(J)+CONJG(A%A(K))*X(I)
          ENDIF
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE ZCSR_SYAMUX_SK
!
!******************************************************************
      SUBROUTINE SM_DDNSDCSR(ADNS,NROW,NCOL,ACSR,UPLO)
      INTEGER NROW,NCOL
      REAL(gq) ADNS(NROW,NCOL)
      TYPE(DCSR_MATRIX)::ACSR
      CHARACTER*1 UPLO
! LOCAL
      INTEGER NNZ,I1,IERR
      TYPE(DCSR_MATRIX)::BUF
!
      IF(UPLO.EQ.'U'.OR.UPLO.EQ.'u')THEN
        DO I1=1,NCOL; ADNS(I1+1:NROW,I1)=0; ENDDO
      ELSEIF(UPLO.EQ.'L'.OR.UPLO.EQ.'l')THEN
        DO I1=1,NCOL; ADNS(1:I1-1,I1)=0; ENDDO
      ENDIF
      NNZ=MIN(MAXNNZ,NROW*NCOL)
      CALL ALLOC_DCSR(BUF,NNZ,NROW,NCOL)
      IERR=0
      CALL DDNSCSR(NROW,NCOL,NNZ,ADNS,NROW,BUF%A,BUF%J,BUF%I,IERR)
      CALL DCSR_COPY(BUF,ACSR)
      CALL DEALLOC_DCSR(BUF)
      RETURN
!
      END SUBROUTINE SM_DDNSDCSR
!
!******************************************************************
      SUBROUTINE SM_ZCSRZDNS(ACSR,ADNS)
      TYPE(ZCSR_MATRIX)::ACSR
      COMPLEX(gq) ADNS(ACSR%NROW,ACSR%NCOL)
! LOCAL
      INTEGER IERR
!
      CALL ZCSRDNS(ACSR%NROW,ACSR%NCOL,ACSR%A,ACSR%J,ACSR%I,ADNS,ACSR%NROW,IERR)
      IF(IERR.NE.0)STOP' ERROR IN SM_ZCSRZDNS!'
      RETURN
!
      END SUBROUTINE SM_ZCSRZDNS
!
!******************************************************************
      SUBROUTINE SM_DCSRDDNS(ACSR,ADNS)
      TYPE(DCSR_MATRIX)::ACSR
      REAL(gq) ADNS(ACSR%NROW,ACSR%NCOL)
! LOCAL
      INTEGER IERR
!
      CALL DCSRDNS(ACSR%NROW,ACSR%NCOL,ACSR%A,ACSR%J,ACSR%I,ADNS,ACSR%NROW,IERR)
      IF(IERR.NE.0)STOP' ERROR IN SM_DCSRDDNS!'
      RETURN
!
      END SUBROUTINE SM_DCSRDDNS
!
!******************************************************************
      SUBROUTINE DCSR_GETDIA(JOB,A,DIAG,IOFF,LADD)
      TYPE(DCSR_MATRIX)::A
      INTEGER JOB,IOFF
      REAL(gq) DIAG(A%NROW)
      LOGICAL LADD
! LOCAL
      INTEGER NNZ,IDIAG(A%NROW)
      REAL(gq) DIAG_(A%NROW)
!
      DIAG_=0      
      CALL DGETDIA(A%NROW,A%NCOL,JOB,A%A,A%J,A%I,NNZ,DIAG_,IDIAG,IOFF)
      IF(LADD)THEN
        DIAG=DIAG+DIAG_
      ELSE
        DIAG=DIAG_
      ENDIF
      RETURN
!
      END SUBROUTINE DCSR_GETDIA
!
!******************************************************************
      SUBROUTINE  TRACE_DCSR(A,TR)
      TYPE(DCSR_MATRIX)::A
      REAL(gq) TR
! LOCAL
      INTEGER I,J
!
      TR=0
      DO I=1,A%NROW; DO J=A%I(I),A%I(I+1)-1
        IF(A%J(J)/=I)CYCLE
        TR=TR+A%A(J); EXIT
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE  TRACE_DCSR
!
!------------------------------ZSM----------------------------------
      SUBROUTINE ALLOC_ZCSR(A,NNZ,NROW,NCOL)
      INTEGER NNZ,NROW,NCOL
      TYPE(ZCSR_MATRIX)::A
!
      !(*,'(" NROW=",I8," NCOL=",I8," NNZ=",I10)')NROW,NCOL,NNZ
      ALLOCATE(A%A(NNZ),A%J(NNZ),A%I(NROW+1))
      A%NROW=NROW; A%NCOL=NCOL
      A%A=0; A%J=0; A%I=0
      GMEM_SIZE=GMEM_SIZE+(REAL(SIZE(A%A),gq)*16+REAL(SIZE(A%J),gq)*4)*1.E-9_gq
      ! 'ALLOC_ZCSR'
      RETURN
!
      END SUBROUTINE ALLOC_ZCSR
!
!******************************************************************
      SUBROUTINE DEALLOC_ZCSR(A)
      TYPE(ZCSR_MATRIX)::A
!
      GMEM_SIZE=GMEM_SIZE-(REAL(SIZE(A%A),gq)*16+REAL(SIZE(A%J),gq)*4)*1.E-9_gq
      DEALLOCATE(A%A,A%J,A%I); A%NROW=0; A%NCOL=0 ! DEALLOCATE IF MEMORY IS NOT ASSOCIATED WITH OTHERS
      ! 'DEALLOC_ZCSR'
      RETURN
!
      END SUBROUTINE DEALLOC_ZCSR
!
!******************************************************************
      SUBROUTINE ALLOC_ZCSR_DIAG(A,N)
      INTEGER N
      TYPE(ZCSR_MATRIX)::A
! LOCAL
      INTEGER I
!
      ALLOCATE(A%A(N),A%J(N),A%I(N+1))
      A%NROW=N; A%NCOL=N; A%A=0
      DO I=1,N; A%J(I)=I; A%I(I)=I; ENDDO
      A%I(N+1)=N+1
      RETURN
!
      END SUBROUTINE ALLOC_ZCSR_DIAG
!
!******************************************************************
      SUBROUTINE ALLOC_DCSR_KKP(A,ID,N,MODE,MODE2)
      TYPE(DCSR_MATRIX)::A
      INTEGER N,ID(N+1),MODE,MODE2
! LOCAL
      INTEGER NNZ,J,NNZL,NPHIK
!
      NNZ=0; NPHIK=ID(N+1)-1
      DO J=1,N
        NNZL=ID(J+1)-ID(J)
        IF(MODE>0)THEN
          NNZ=NNZ+NNZL*MODE
        ELSE
          IF(ABS(MODE2)==1)THEN
            NNZ=NNZ+(NNZL+1)*NNZL/2
          ELSE
            NNZ=NNZ+NNZL*NNZL
          ENDIF
        ENDIF
      ENDDO
      CALL ALLOC_DCSR(A,NNZ,NPHIK,NPHIK)
      RETURN
!
      END SUBROUTINE ALLOC_DCSR_KKP
!
!******************************************************************
      SUBROUTINE ALLOC_ZCSR_KKP(A,ID,N,MODE,MODE2)
      TYPE(ZCSR_MATRIX)::A
      INTEGER N,ID(N+1),MODE,MODE2
! LOCAL
      INTEGER NNZ,J,NNZL,NPHIK
!
      NNZ=0; NPHIK=ID(N+1)-1
      DO J=1,N
        NNZL=ID(J+1)-ID(J)
        IF(MODE>0)THEN
          NNZ=NNZ+NNZL*MODE
        ELSE
          IF(ABS(MODE2)==1)THEN
            NNZ=NNZ+(NNZL+1)*NNZL/2
          ELSE
            NNZ=NNZ+NNZL*NNZL
          ENDIF
        ENDIF
      ENDDO
      CALL ALLOC_ZCSR(A,NNZ,NPHIK,NPHIK)
      RETURN
!
      END SUBROUTINE ALLOC_ZCSR_KKP
!
!******************************************************************
! non-Hermitian, store full matrix, form occ: <N+1, N>
!******************************************************************
      SUBROUTINE ALLOC_ZCSR_KKP1(A,ID,N)
      TYPE(ZCSR_MATRIX)::A
      INTEGER N,ID(N+1)
! LOCAL
      INTEGER NNZ,J,NNZL,NNZR,NPHIK
!
      NNZ=0; NPHIK=ID(N+1)-1
      DO J=2,N
        NNZL=ID(J+1)-ID(J); NNZR=ID(J)-ID(J-1)
        NNZ=NNZ+NNZL*NNZR
      ENDDO
      CALL ALLOC_ZCSR(A,NNZ,NPHIK,NPHIK)
      RETURN
!
      END SUBROUTINE ALLOC_ZCSR_KKP1
!
!******************************************************************
      SUBROUTINE ALLOC_ZCOO(A,NNZ,NROW,NCOL)
      TYPE(ZCOO_MATRIX)::A
      INTEGER NROW,NCOL,NNZ
!
      !(*,'(" NROW=",I8," NCOL=",I8," NNZ=",I10)')NROW,NCOL,NNZ
      A%NNZ=NNZ; A%NROW=NROW; A%NCOL=NCOL
      ALLOCATE(A%A(NNZ),A%I(NNZ),A%J(NNZ))
      A%A=0; A%I=0; A%J=0
      GMEM_SIZE=GMEM_SIZE+REAL(SIZE(A%A),gq)*48*1.E-9_gq
      ! 'ALLOC_ZCOO'
      RETURN
!
      END SUBROUTINE ALLOC_ZCOO
!
!******************************************************************
      SUBROUTINE ALLOC_DCOO(A,NNZ,NROW,NCOL)
      TYPE(DCOO_MATRIX)::A
      INTEGER NROW,NCOL,NNZ
!
      !(*,'(" NROW=",I8," NCOL=",I8," NNZ=",I10)')NROW,NCOL,NNZ
      A%NNZ=NNZ; A%NROW=NROW; A%NCOL=NCOL
      ALLOCATE(A%A(NNZ),A%I(NNZ),A%J(NNZ))
      A%A=0; A%I=0; A%J=0
      GMEM_SIZE=GMEM_SIZE+REAL(SIZE(A%A),gq)*32*1.E-9_gq
      ! 'ALLOC_DCOO'
      RETURN
!
      END SUBROUTINE ALLOC_DCOO
!
!******************************************************************
      SUBROUTINE DBM_TRANSPOSE(A)
      TYPE(DBD_MATRIX) :: A
! LOCAL
      INTEGER I
!
      IF(A%JSHIFT/=0)STOP' ERROR IN DBM_TRANSPOSE: ONLY FOR BLOCK DIAGONAL MATRIX!'
      DO I=1,A%NBK
        A%BK(I)%A=TRANSPOSE(A%BK(I)%A)
      ENDDO
      RETURN
!
      END SUBROUTINE DBM_TRANSPOSE
!
!******************************************************************
      SUBROUTINE DEALLOC_ZCOO(A)
      TYPE(ZCOO_MATRIX)::A
!
      A%NNZ=0; A%NROW=0; A%NCOL=0
      GMEM_SIZE=GMEM_SIZE-REAL(SIZE(A%A),gq)*48*1.E-9_gq
      DEALLOCATE(A%A,A%I,A%J)
      ! 'DEALLOC_ZCOO'
      RETURN
!
      END SUBROUTINE DEALLOC_ZCOO
!
!******************************************************************
      SUBROUTINE DEALLOC_DCOO(A)
      TYPE(DCOO_MATRIX)::A
!
      A%NNZ=0; A%NROW=0; A%NCOL=0
      GMEM_SIZE=GMEM_SIZE-REAL(SIZE(A%A),gq)*32*1.E-9_gq
      DEALLOCATE(A%A,A%I,A%J)
      ! 'DEALLOC_DCOO'
      RETURN
!
      END SUBROUTINE DEALLOC_DCOO
!
!*****************************************************************************
      SUBROUTINE ZCSR_COPY(A,B)
      TYPE(ZCSR_MATRIX)::A,B
! LOCAL
      INTEGER NROW,NCOL,NNZ
!
      NROW=A%NROW; NCOL=A%NCOL; NNZ=A%I(NROW+1)-1
      CALL ALLOC_ZCSR(B,NNZ,NROW,NCOL)
      CALL ZCSR_COPY1(A,B)
      RETURN
!
      END SUBROUTINE ZCSR_COPY
!
!*****************************************************************************
      SUBROUTINE ZCSR_COPY1(A,B)
      TYPE(ZCSR_MATRIX)::A,B
! LOCAL
      INTEGER NROW,NNZ
!
      NROW=A%NROW; NNZ=A%I(NROW+1)-1
      B%NROW=A%NROW; B%NCOL=A%NCOL
      B%A=A%A(1:NNZ); B%J=A%J(1:NNZ); B%I=A%I
      RETURN
!
      END SUBROUTINE ZCSR_COPY1
!
!*****************************************************************************
      SUBROUTINE ZCOO_COPY(A,B)
      TYPE(ZCOO_MATRIX)::A,B
! LOCAL
      INTEGER NROW,NCOL,NNZ
!
      NROW=A%NROW; NCOL=A%NCOL; NNZ=A%NNZ
      CALL ALLOC_ZCOO(B,NNZ,NROW,NCOL)
      CALL ZCOO_COPY1(A,B)
      RETURN
!
      END SUBROUTINE ZCOO_COPY
!
!*****************************************************************************
      SUBROUTINE ZCOO_COPY1(A,B)
      TYPE(ZCOO_MATRIX)::A,B
! LOCAL
      INTEGER NNZ
!
      B%NROW=A%NROW; B%NCOL=A%NCOL; B%NNZ=A%NNZ; NNZ=A%NNZ
      B%A=A%A(1:NNZ); B%I=A%I(1:NNZ); B%J=A%J(1:NNZ)
      RETURN
!
      END SUBROUTINE ZCOO_COPY1
!
!*********************************************************************
      SUBROUTINE ZCSR_APLSB_SK(TRANS,A,S,B,C)
      CHARACTER*1 TRANS
      COMPLEX(gq) S
      TYPE(ZCSR_MATRIX)::A,B
      TYPE(ZCSR_MATRIX),OPTIONAL::C
! LOCAL
      INTEGER NNZ,IERR,NNZ1,NNZ2,NROW1,NCOL1,NROW2,NCOL2
      INTEGER,ALLOCATABLE::IW(:)
      TYPE(ZCSR_MATRIX)::B_,D
!
      IF(TRANS.EQ.'N'.OR.TRANS.EQ.'n')THEN
        CALL ZCSR_COPY(B,B_)
      ELSE
        CALL SM_ZCSRZCSC(B,AT=B_)
        IF(TRANS.EQ.'C'.OR.TRANS.EQ.'c') B_%A=CONJG(B_%A)
      ENDIF
      NROW1=A %NROW; NCOL1=A %NCOL
      NROW2=B_%NROW; NCOL2=B_%NCOL
      IF(NROW1.NE.NROW2)STOP' ERROR: NROW1.NE.NROW2 IN ZCSR_APLSB_SK!'
      IF(NCOL1.NE.NCOL2)STOP' ERROR: NCOL1.NE.NCOL2 IN ZCSR_APLSB_SK!'
      ALLOCATE(IW(NCOL1)); IW=0
      NNZ1=A%I(NROW1+1)-1; NNZ2=B_%I(NROW2+1)-1
      NNZ=MIN(MAXNNZ,NNZ1+NNZ2)
      CALL ALLOC_ZCSR(D,NNZ,NROW1,NCOL1)
      IERR=0
      CALL ZAPLSB(NROW1,NCOL1,A%A,A%J,A%I,S,B_%A,B_%J,B_%I,D%A,D%J,D%I,NNZ,IW,IERR)
      IF(IERR.NE.0)STOP' ERROR IN DCSR_APLSB_SK!'
      CALL DEALLOC_ZCSR(B_)
      DEALLOCATE(IW)
      IF(PRESENT(C))THEN
        CALL ZCSR_COPY(D,C)
      ELSE
        CALL DEALLOC_ZCSR(A)
        CALL ZCSR_COPY(D,A)
      ENDIF
      CALL DEALLOC_ZCSR(D)
      RETURN
!
      END SUBROUTINE ZCSR_APLSB_SK
!
!*****************************************************************************
      SUBROUTINE ZBM_APLSB(A,S,B)
      TYPE(ZBD_MATRIX),INTENT(INOUT)::A
      COMPLEX(gq),INTENT(IN)::S
      TYPE(ZBD_MATRIX),INTENT(IN)::B
! LOCAL
      INTEGER I
!
      DO I=1,A%NBK
        IF(A%BK(I)%LZERO)CYCLE
        IF(A%BK(I)%LSPARSE)THEN
          CALL ZCSR_APLSB_SK('N',A%BK(I)%ACSR,S,B%BK(I)%ACSR)
        ELSE
          A%BK(I)%A=A%BK(I)%A+S*B%BK(I)%A
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE ZBM_APLSB
!
!*****************************************************************************
      SUBROUTINE SM_ZCSRZCSC(A,AT)
      TYPE(ZCSR_MATRIX)::A
      TYPE(ZCSR_MATRIX),OPTIONAL::AT
! LOCAL
      INTEGER NROW,NCOL,NNZ
      TYPE(ZCSR_MATRIX)::B
!
      NROW=A%NROW; NCOL=A%NCOL
      NNZ=A%I(NROW+1)-1
      CALL ALLOC_ZCSR(B,NNZ,NCOL,NROW) ! Transpose
      CALL ZCSRCSC2(NROW,NCOL,1,1,A%A,A%J,A%I,B%A,B%J,B%I)
      IF(PRESENT(AT))THEN; CALL ZCSR_COPY(B,AT)
      ELSE; CALL DEALLOC_ZCSR(A); CALL ZCSR_COPY(B,A)
      ENDIF
      CALL DEALLOC_ZCSR(B)
      RETURN
!
      END SUBROUTINE SM_ZCSRZCSC
!
!*********************************************************************
      SUBROUTINE ZCSR_AMUB_SK(TRANS,A,B,C,LB)
      CHARACTER*1 TRANS
      TYPE(ZCSR_MATRIX)::A,B
      TYPE(ZCSR_MATRIX),OPTIONAL::C
      INTEGER,OPTIONAL::LB
! LOCAL
      INTEGER NNZ,IERR,NROW1,NCOL1,NROW2,NCOL2
      INTEGER,ALLOCATABLE::IW(:)
      TYPE(ZCSR_MATRIX)::A_,D
!
      IF(TRANS.EQ.'N'.OR.TRANS.EQ.'n')THEN
        CALL ZCSR_COPY(A,A_)
      ELSE
        CALL SM_ZCSRZCSC(A,AT=A_)
        IF(TRANS.EQ.'C'.OR.TRANS.EQ.'c')THEN
          A_%A=CONJG(A_%A)
        ENDIF
      ENDIF
      NROW1=A_%NROW; NCOL1=A_%NCOL
      NROW2=B %NROW; NCOL2=B %NCOL
      ALLOCATE(IW(NCOL1)); IW=0
      IF(NCOL1.NE.NROW2)THEN
        WRITE(0,'(" NCOL1 vs NROW2:",2I6)')NCOL1,NROW2
        WRITE(0,'(" NROW1 vs NCOL2:",2I6)')NROW1,NCOL2
        STOP' ERROR: NCOL1.NE.NROW2 IN ZCSR_AMUB_SK!'
      ENDIF
      NNZ=MIN(MAXNNZ,NROW1*NCOL2)
      CALL ALLOC_ZCSR(D,NNZ,NROW1,NCOL2)
      IERR=0
      CALL ZAMUB(NROW1,NCOL1,1,A_%A,A_%J,A_%I,B%A,B%J,B%I,D%A,D%J,D%I,NNZ,IW,IERR)
      DEALLOCATE(IW)
      CALL DEALLOC_ZCSR(A_)
      IF(IERR.NE.0)STOP' ERROR IN ZCSR_AMUB_SK!'
      IF(PRESENT(C))THEN
        CALL ZCSR_COPY(D,C)
      ELSEIF(PRESENT(LB))THEN
        CALL DEALLOC_ZCSR(B)
        CALL ZCSR_COPY(D,B)
      ELSE
        CALL DEALLOC_ZCSR(A)
        CALL ZCSR_COPY(D,A)
      ENDIF
      CALL DEALLOC_ZCSR(D)
      RETURN
!
      END SUBROUTINE ZCSR_AMUB_SK
!
!******************************************************************
      SUBROUTINE ZCSR_GETDIA(JOB,A,DIAG,IOFF,LADD)
      TYPE(DCSR_MATRIX)::A
      INTEGER JOB,IOFF
      COMPLEX(gq) DIAG(A%NROW)
      LOGICAL LADD
! LOCAL
      INTEGER NNZ,IDIAG(A%NROW)
      COMPLEX(gq) DIAG_(A%NROW)
!
      DIAG_=0
      CALL ZGETDIA(A%NROW,A%NCOL,JOB,A%A,A%J,A%I,NNZ,DIAG_,IDIAG,IOFF)
      IF(LADD)THEN
        DIAG=DIAG+DIAG_
      ELSE
        DIAG=DIAG_
      ENDIF
      RETURN
!
      END SUBROUTINE ZCSR_GETDIA
!
!******************************************************************
      SUBROUTINE SM_ZDNSZCSR(ADNS,NROW,NCOL,ACSR,UPLO)
      INTEGER NROW,NCOL
      COMPLEX(gq) ADNS(NROW,NCOL)
      TYPE(ZCSR_MATRIX)::ACSR
      CHARACTER*1 UPLO
! LOCAL
      INTEGER NNZ,I1,IERR
      TYPE(ZCSR_MATRIX)::BUF
!
      IF(UPLO.EQ.'U'.OR.UPLO.EQ.'u')THEN
        DO I1=1,NCOL; ADNS(I1+1:NROW,I1)=0; ENDDO
      ELSEIF(UPLO.EQ.'L'.OR.UPLO.EQ.'l')THEN
        DO I1=1,NCOL; ADNS(1:I1-1,I1)=0; ENDDO
      ENDIF
      NNZ=MIN(MAXNNZ,NROW*NCOL)
      CALL ALLOC_ZCSR(BUF,NNZ,NROW,NCOL)
      IERR=0
      CALL ZDNSCSR(NROW,NCOL,NNZ,ADNS,NROW,BUF%A,BUF%J,BUF%I,IERR)
      NNZ=BUF%I(NROW+1)-1
      CALL ZCSR_COPY(BUF,ACSR)
      CALL DEALLOC_ZCSR(BUF)
      RETURN
!
      END SUBROUTINE SM_ZDNSZCSR
!
!******************************************************************
      SUBROUTINE SM_ZCOOZCSR(ACOO,BCSR)
      TYPE(ZCOO_MATRIX)::ACOO
      TYPE(ZCSR_MATRIX)::BCSR
! LOCAL
      INTEGER NROW,NCOL,NNZ
      INTEGER,ALLOCATABLE::IR(:)
!
      NROW=ACOO%NROW; NCOL=ACOO%NCOL; NNZ=ACOO%NNZ
      ALLOCATE(IR(NNZ)); IR=ACOO%I
      CALL ALLOC_ZCSR(BCSR,NNZ,NROW,NCOL)
      CALL ZCOOCSR(NROW,NNZ,ACOO%A,IR,ACOO%J,BCSR%A,BCSR%J,BCSR%I)
      DEALLOCATE(IR)
      RETURN
!
      END SUBROUTINE SM_ZCOOZCSR
!
!******************************************************************
      SUBROUTINE SM_DCOOZCSR(ACOO,BCSR)
      TYPE(DCOO_MATRIX)::ACOO
      TYPE(ZCSR_MATRIX)::BCSR
! LOCAL
      INTEGER NROW,NCOL,NNZ
      INTEGER,ALLOCATABLE::IR(:)
      COMPLEX(gq),ALLOCATABLE::A(:)
!
      NROW=ACOO%NROW; NCOL=ACOO%NCOL; NNZ=ACOO%NNZ
      ALLOCATE(IR(NNZ),A(NNZ)); IR=ACOO%I; A=DCMPLX(ACOO%A)
      CALL ALLOC_ZCSR(BCSR,NNZ,NROW,NCOL)
      CALL ZCOOCSR(NROW,NNZ,A,IR,ACOO%J,BCSR%A,BCSR%J,BCSR%I)
      DEALLOCATE(IR,A)
      RETURN
!
      END SUBROUTINE SM_DCOOZCSR
!
!*********************************************************************
! Direct sum of two csr matrices to form a "block diagonal" matrix
!*********************************************************************
      SUBROUTINE ZCSR_DSUM_AB(A,B,AB,LFIRST)
      TYPE(ZCSR_MATRIX)::A,B
      TYPE(ZCSR_MATRIX),OPTIONAL::AB
      LOGICAL,OPTIONAL::LFIRST
! LOCAL
      INTEGER NROW1,NROW2,NCOL1,NCOL2,NNZ1,NNZ2
      INTEGER NROW,NCOL,NNZ
      TYPE(ZCSR_MATRIX)::C
!
      IF(PRESENT(LFIRST))THEN
        IF(LFIRST)THEN
          IF(PRESENT(AB))THEN
            CALL ZCSR_COPY(B,AB)
          ELSE
            CALL ZCSR_COPY(B,A)
          ENDIF
          RETURN
        ENDIF
      ENDIF
      NROW1=A%NROW; NCOL1=A%NCOL; NNZ1=A%I(NROW1+1)-1
      NROW2=B%NROW; NCOL2=B%NCOL; NNZ2=B%I(NROW2+1)-1
      NROW=NROW1+NROW2; NCOL=NCOL1+NCOL2; NNZ=NNZ1+NNZ2
      CALL ALLOC_ZCSR(C,NNZ,NROW,NCOL)
      C%A(1:NNZ1) =A%A; C%J(1:NNZ1) =A%J      ; C%I(1:NROW1) =A%I(1:NROW1)
      C%A(1+NNZ1:)=B%A; C%J(1+NNZ1:)=B%J+NCOL1; C%I(1+NROW1:)=B%I+NNZ1
      IF(PRESENT(AB))THEN; CALL ZCSR_COPY(C,AB)
      ELSE; CALL DEALLOC_ZCSR(A); CALL ZCSR_COPY(C,A)
      ENDIF
      CALL DEALLOC_ZCSR(C)
      RETURN
!
      END SUBROUTINE ZCSR_DSUM_AB
!
!
!******************************************************************
      SUBROUTINE ZCSR_DSUM_DNS(ACSR,ADNS,N,M,UPLO,LFIRST)
      TYPE(ZCSR_MATRIX)::ACSR
      INTEGER N,M
      COMPLEX(gq)ADNS(N,M)
      CHARACTER*1 UPLO
      LOGICAL LFIRST
! LOCAL
      TYPE(ZCSR_MATRIX)::BUF
!
      IF(LFIRST)THEN
        CALL SM_ZDNSZCSR(ADNS,N,M,ACSR,UPLO)
      ELSE
        CALL SM_ZDNSZCSR(ADNS,N,M,BUF,UPLO)
        CALL ZCSR_DSUM_AB(ACSR,BUF)
        CALL DEALLOC_ZCSR(BUF)
      ENDIF
      RETURN
!
      END SUBROUTINE ZCSR_DSUM_DNS
!
!*****************************************************************************
! Take out sub-matrix (N1:N2,M1:M2) of a BLOCK DIAGONAL CSR sparse matrix
!*****************************************************************************
      SUBROUTINE ZCSR_DIA_SUBMAT(A,ABLK,N1,N2,M1,M2)
      INTEGER N1,M1,N2,M2
      TYPE(ZCSR_MATRIX)::A,ABLK
! LOCAL
      INTEGER NROW,NCOL,NNZ
!
      NROW=N2-N1+1; NCOL=M2-M1+1
      NNZ=A%I(N2+1)-A%I(N1)
      CALL ALLOC_ZCSR(ABLK,NNZ,NROW,NCOL)
      ABLK%A(1:NNZ)=A%A(A%I(N1):A%I(N2+1)-1)
      ABLK%J(1:NNZ)=A%J(A%I(N1):A%I(N2+1)-1)-M1+1
      ABLK%I(1:NROW+1)=A%I(N1:N2+1)-A%I(N1)+1
      IF(MINVAL(ABLK%J).LE.0.OR.MAXVAL(ABLK%J).GT.NCOL)THEN
        WRITE(0,'(" NCOL=",I10)')NCOL
        WRITE(0,'(" ABLK%J=")'); WRITE(0,'(10I10)')ABLK%J
        WRITE(0,'(" ABLK%A=")'); WRITE(0,'(5F20.12)')ABLK%A
        STOP 'FETAL ERROR IN ZCSR_DIA_SUBMAT: A NOT BLOCK DIAGONAL!'
      ENDIF
      RETURN
!
      END SUBROUTINE ZCSR_DIA_SUBMAT
!
!*****************************************************************************
      SUBROUTINE ZCSR_DIADNS_SUBMAT(ACSR,ABLK,N1,N2,M1,M2)
      INTEGER N1,M1,N2,M2
      TYPE(ZCSR_MATRIX)::ACSR
      COMPLEX(gq)ABLK(N2-N1+1,M2-M1+1)
! LOCAL
      TYPE(ZCSR_MATRIX)::BUF
!
      CALL ZCSR_DIA_SUBMAT(ACSR,BUF,N1,N2,M1,M2)
      CALL SM_ZCSRZDNS(BUF,ABLK)
      RETURN
!
      END SUBROUTINE ZCSR_DIADNS_SUBMAT
!
!*****************************************************************************
      SUBROUTINE ZCSR_SUBMAT(A,N1,N2,M1,M2,ABLK)
      INTEGER N1,M1,N2,M2
      TYPE(ZCSR_MATRIX)::A
      TYPE(ZCSR_MATRIX),OPTIONAL::ABLK
! LOCAL
      INTEGER NROW,NCOL,NNZ
      TYPE(ZCSR_MATRIX)::BUF
!
      NROW=N2-N1+1; NCOL=M2-M1+1
      NNZ=A%I(N2+1)-A%I(N1)
      CALL ALLOC_ZCSR(BUF,NNZ,NROW,NCOL)
      CALL ZSUBMAT(A%NROW,1,N1,N2,M1,M2,A%A,A%J,A%I,NROW,NCOL,BUF%A,BUF%J,BUF%I)
      IF(PRESENT(ABLK))THEN
        CALL ZCSR_COPY(BUF,ABLK)
      ELSE
        CALL DEALLOC_ZCSR(A)
        CALL ZCSR_COPY(BUF,A)
      ENDIF
      CALL DEALLOC_ZCSR(BUF)
      RETURN
!
      END SUBROUTINE ZCSR_SUBMAT
!
!*****************************************************************************
      SUBROUTINE ZCSR_UHAU(A,U)
      TYPE(ZCSR_MATRIX)::A,U
! LOCAL
      TYPE(ZCSR_MATRIX)::BUF
!
      ! "ZCSR_UHAU"
      CALL ZCSR_AMUB_SK('N',A,U,C=BUF)
      CALL DEALLOC_ZCSR(A)
      CALL ZCSR_AMUB_SK('C',U,BUF,C=A)
      CALL DEALLOC_ZCSR(BUF)
      ! "ZCSR_UHAU"
      RETURN
!
      END SUBROUTINE ZCSR_UHAU
!
!******************************************************************
      SUBROUTINE ZCSR_SIMPLIFY(A,SMALL,UPLO)
      TYPE(ZCSR_MATRIX)::A
      REAL(gq) SMALL
      CHARACTER*1 UPLO
! LOCAL
      INTEGER NROW,NNZ,IERR
      TYPE(ZCSR_MATRIX)::BUF
!
      NROW=A%NROW; NNZ=A%I(NROW+1)-1
      IERR=0
      CALL ZFILTER(NROW,1,SMALL,A%A,A%J,A%I,A%A,A%J,A%I,NNZ,IERR) ! in place
      IF(IERR.NE.0)STOP' ERROR IN DCSR_SIMPLIFY!'
      IF(UPLO.EQ.'U'.OR.UPLO.EQ.'u')THEN
        CALL ZGETU(NROW,A%A,A%J,A%I,A%A,A%J,A%I) ! in place
      ELSEIF(UPLO.EQ.'L'.OR.UPLO.EQ.'l')THEN
        CALL ZGETL(NROW,A%A,A%J,A%I,A%A,A%J,A%I) ! in place
      ENDIF
      CALL ZCSR_COPY(A,BUF)
      CALL DEALLOC_ZCSR(A)
      CALL ZCSR_COPY(BUF,A)
      CALL DEALLOC_ZCSR(BUF)
      RETURN
!
      END SUBROUTINE ZCSR_SIMPLIFY
!
!******************************************************************
      SUBROUTINE ZCSR_SHRINK(A)
      TYPE(ZCSR_MATRIX)::A
! LOCAL
      TYPE(ZCSR_MATRIX)::BUF
!
      CALL ZCSR_COPY(A,BUF)
      CALL DEALLOC_ZCSR(A)
      CALL ZCSR_COPY(BUF,A)
      CALL DEALLOC_ZCSR(BUF)
      RETURN
!
      END SUBROUTINE ZCSR_SHRINK
!
!******************************************************************
      SUBROUTINE DCSR_SET_ROW(A, I_ROW, A_ROW)
      TYPE(DCSR_MATRIX),INTENT(INOUT) :: A
      INTEGER, INTENT(IN) :: I_ROW
      REAL(gq), INTENT(IN) :: A_ROW(A%NCOL)
! LOCAL
      INTEGER J
!
      IF(I_ROW == 1) A%I(1) = 1
      A%I(I_ROW + 1) = A%I(I_ROW)
      DO J = 1, A%NCOL
        IF(ABS(A_ROW(J)) < SMALL) CYCLE
        A%A(A%I(I_ROW + 1)) = A_ROW(J)
        A%J(A%I(I_ROW + 1)) = J
        A%I(I_ROW + 1) = A%I(I_ROW + 1) + 1
      ENDDO
      RETURN
!
      END SUBROUTINE DCSR_SET_ROW
!
!******************************************************************
      SUBROUTINE ZCSR_SET_ROW(A, I_ROW, A_ROW)
      TYPE(ZCSR_MATRIX),INTENT(INOUT) :: A
      INTEGER, INTENT(IN) :: I_ROW
      COMPLEX(gq), INTENT(IN) :: A_ROW(A%NCOL)
! LOCAL
      INTEGER J
!
      IF(I_ROW == 1) A%I(1) = 1
      A%I(I_ROW + 1) = A%I(I_ROW)
      DO J = 1, A%NCOL
        IF(ABS(A_ROW(J)) < SMALL) CYCLE
        A%A(A%I(I_ROW + 1)) = A_ROW(J)
        A%J(A%I(I_ROW + 1)) = J
        A%I(I_ROW + 1) = A%I(I_ROW + 1) + 1
      ENDDO
      RETURN
!
      END SUBROUTINE ZCSR_SET_ROW
!
!******************************************************************
      SUBROUTINE DCSR_SIMPLIFY(A,SMALL,UPLO)
      TYPE(DCSR_MATRIX)::A
      REAL(gq) SMALL
      CHARACTER*1 UPLO
! LOCAL
      INTEGER NROW,NNZ,IERR
      TYPE(DCSR_MATRIX)::BUF
!
      NROW=A%NROW; NNZ=A%I(NROW+1)-1
      IERR=0
      CALL DFILTER(NROW,1,SMALL,A%A,A%J,A%I,A%A,A%J,A%I,NNZ,IERR) ! in place
      IF(IERR.NE.0)STOP' ERROR IN DCSR_SIMPLIFY!'
      IF(UPLO.EQ.'U'.OR.UPLO.EQ.'u')THEN
        CALL DGETU(NROW,A%A,A%J,A%I,A%A,A%J,A%I) ! in place
      ELSEIF(UPLO.EQ.'L'.OR.UPLO.EQ.'l')THEN
        CALL DGETL(NROW,A%A,A%J,A%I,A%A,A%J,A%I) ! in place
      ENDIF
      CALL DCSR_COPY(A,BUF)
      CALL DEALLOC_DCSR(A)
      CALL DCSR_COPY(BUF,A)
      CALL DEALLOC_DCSR(BUF)
      RETURN
!
      END SUBROUTINE DCSR_SIMPLIFY
!
!******************************************************************
      SUBROUTINE DCSR_SHRINK(A)
      TYPE(DCSR_MATRIX)::A
! LOCAL
      TYPE(DCSR_MATRIX)::BUF
!
      CALL DCSR_COPY(A,BUF)
      CALL DEALLOC_DCSR(A)
      CALL DCSR_COPY(BUF,A)
      CALL DEALLOC_DCSR(BUF)
      RETURN
!
      END SUBROUTINE DCSR_SHRINK
!
!****************************************************************************
! MODE > 0: Trace( S(k)^\dagger U SP(k') )
!     else: Trace( SP(k')       U S(k)^\dagger )
! Index order: A_J1,I1 * U_I1,I2 * B_I2,J2
!****************************************************************************
      SUBROUTINE TR_DCOODDNSDCOO(TR,S,A,SP,N,IBASE,MODE)
      TYPE(DCOO_MATRIX)::S,SP
      INTEGER N,IBASE,MODE
      REAL(gq) TR,A(N,N)
! LOCAL
      INTEGER INZ,INZP,I1,J1,I2,J2
!
      TR=0
      DO INZ=1,S%NNZ
      IF(MODE>0)THEN
        J1=S%J(INZ); I1=S%I(INZ)
      ELSE
        J2=S%I(INZ); I2=S%J(INZ)
      ENDIF
      DO INZP=1,SP%NNZ
      IF(MODE>0)THEN
        J2=SP%J(INZP); I2=SP%I(INZP)
      ELSE
        J1=SP%I(INZP); I1=SP%J(INZP)
      ENDIF
      IF(J1/=J2)CYCLE
      TR=TR+S%A(INZ)*A(I1-IBASE,I2-IBASE)*SP%A(INZP)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE TR_DCOODDNSDCOO
!
!****************************************************************************
      SUBROUTINE TR_DCOOZDNSDCOO(TR,S,A,SP,N,IBASE,MODE)
      TYPE(DCOO_MATRIX)::S,SP
      INTEGER N,IBASE,MODE
      COMPLEX(gq) TR,A(N,N)
! LOCAL
      INTEGER INZ,INZP,I1,J1,I2,J2
!
      TR=0
      DO INZ=1,S%NNZ
      IF(MODE>0)THEN
        J1=S%J(INZ); I1=S%I(INZ)
      ELSE
        J2=S%I(INZ); I2=S%J(INZ)
      ENDIF
      DO INZP=1,SP%NNZ
      IF(MODE>0)THEN
        J2=SP%J(INZP); I2=SP%I(INZP)
      ELSE
        J1=SP%I(INZP); I1=SP%J(INZP)
      ENDIF
      IF(J1/=J2)CYCLE
      TR=TR+S%A(INZ)*A(I1-IBASE,I2-IBASE)*SP%A(INZP)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE TR_DCOOZDNSDCOO
!
!****************************************************************************
! Trace( S(k)^\dagger f^\dagger SP(k') f )
! Index order: S*_J1,I1  F2^+_I1,I2  S_I2,J2  F1_J2,J1
!****************************************************************************
      SUBROUTINE TR_DCOODDNSDCOODDNS(TR,S,F2,SP,F1,N1,N2,IBASE,JBASE)
      TYPE(DCOO_MATRIX)::S,SP
      INTEGER N1,N2,IBASE,JBASE
      REAL(gq) TR,F1(N1,N2),F2(N1,N2)
! LOCAL
      INTEGER INZ,INZP,I1,J1,I2,J2
!
      TR=0
      DO INZ=1,S%NNZ
      J1=S%J(INZ); I1=S%I(INZ)
      DO INZP=1,SP%NNZ
      J2=SP%J(INZP); I2=SP%I(INZP)
      TR=TR+S%A(INZ)*F2(I2-IBASE,I1-JBASE)*SP%A(INZP)*F1(J2-IBASE,J1-JBASE)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE TR_DCOODDNSDCOODDNS
!
!****************************************************************************
      SUBROUTINE TR_DCOOZDNSDCOOZDNS(TR,S,F2,SP,F1,N1,N2,IBASE,JBASE)
      TYPE(DCOO_MATRIX)::S,SP
      INTEGER N1,N2,IBASE,JBASE
      COMPLEX(gq) TR,F1(N1,N2),F2(N1,N2)
! LOCAL
      INTEGER INZ,INZP,I1,J1,I2,J2
!
      TR=0
      DO INZ=1,S%NNZ
      J1=S%J(INZ); I1=S%I(INZ)
      DO INZP=1,SP%NNZ
      J2=SP%J(INZP); I2=SP%I(INZP)
      TR=TR+S%A(INZ)*CONJG(F2(I2-IBASE,I1-JBASE))*SP%A(INZP)*F1(J2-IBASE,J1-JBASE)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE TR_DCOOZDNSDCOOZDNS
!
!*****************************************************************************
      SUBROUTINE ZCSR_BLK_MAXODE(A,NBK,ID,MAXODE)
      INTEGER NBK,ID(NBK+1)
      REAL(gq) MAXODE,RES
      TYPE(ZCSR_MATRIX)::A
! LOCAL
      INTEGER IBK,I,J,IJ,N1,N2
!
      DO IBK=1,NBK; N1=ID(IBK); N2=ID(IBK+1)-1
      DO I=N1,N2; DO J=A%I(I),A%I(I+1)-1
        IJ=A%J(J)
        IF(IJ>=N1.AND.IJ<=N2)CYCLE
        RES=ABS(A%A(J))
        IF(MAXODE>RES)CYCLE
        MAXODE=RES
      ENDDO; ENDDO; ENDDO ! IBK
      RETURN
!
      END SUBROUTINE ZCSR_BLK_MAXODE
!
!*********************************************************************
      SUBROUTINE DBMTODCSR(AB,AS)
      TYPE(DBD_MATRIX)::AB
      TYPE(DCSR_MATRIX)::AS
! LOCAL
      INTEGER I,IMIN,IMAX,NDROW,NDCOL
!
      CALL DBM_CYC_BOUND(AB,IMIN,IMAX,NDROW,NDCOL)
      DO I=IMIN,IMAX
        IF(AB%BK(I)%LSPARSE)THEN
          CALL DCSR_DSUM_AB(AS,AB%BK(I)%ACSR,LFIRST=(I==IMIN))
        ELSE
          CALL DCSR_DSUM_DNS(AS,AB%BK(I)%A,AB%BK(I)%NROW,AB%BK(I)%NCOL,'N',I==IMIN)
        ENDIF
      ENDDO
      CALL DCSR_ROWCOL_SHIFT(AS,NDROW,NDCOL,AB%JSHIFT)
      RETURN
!
      END SUBROUTINE DBMTODCSR
!
!*********************************************************************
      SUBROUTINE DBMTOBKDCSR(AB,RCUT_)
      TYPE(DBD_MATRIX)::AB
      REAL(gq),OPTIONAL::RCUT_
! LOCAL
      INTEGER I,IMIN,IMAX
      REAL(gq) RCUT
!
      RCUT=1.E-16_gq
      IF(PRESENT(RCUT_))RCUT=RCUT_
      SELECT CASE(AB%JSHIFT)
      CASE(-1); IMIN=2; IMAX=AB%NBK
      CASE( 0); IMIN=1; IMAX=AB%NBK
      CASE( 1); IMIN=1; IMAX=AB%NBK-1
      END SELECT
      DO I=IMIN,IMAX
        AB%BK(I)%LSPARSE=.TRUE.
        CALL SM_DDNSDCSR(AB%BK(I)%A,AB%BK(I)%NROW,AB%BK(I)%NCOL,AB%BK(I)%ACSR,'A')
        DEALLOCATE(AB%BK(I)%A)
        CALL DCSR_SIMPLIFY(AB%BK(I)%ACSR,RCUT,'N')
        GMEM_SIZE=GMEM_SIZE-REAL(SIZE(AB%BK(I)%A),gq)*8*1.E-9_gq
      ENDDO
      ! 'DBMTOBKDCSR'
      RETURN
!
      END SUBROUTINE DBMTOBKDCSR
!
!*********************************************************************
      SUBROUTINE ZBMTOBKZCSR(AB,RCUT_)
      TYPE(ZBD_MATRIX)::AB
      REAL(gq),OPTIONAL::RCUT_
! LOCAL
      INTEGER I,IMIN,IMAX
      REAL(gq) RCUT
!
      RCUT=1.E-16_gq
      IF(PRESENT(RCUT_))RCUT=RCUT_
      SELECT CASE(AB%JSHIFT)
      CASE(-1); IMIN=2; IMAX=AB%NBK
      CASE( 0); IMIN=1; IMAX=AB%NBK
      CASE( 1); IMIN=1; IMAX=AB%NBK-1
      END SELECT
      DO I=IMIN,IMAX
        AB%BK(I)%LSPARSE=.TRUE.
        CALL SM_ZDNSZCSR(AB%BK(I)%A,AB%BK(I)%NROW,AB%BK(I)%NCOL,AB%BK(I)%ACSR,'A')
        DEALLOCATE(AB%BK(I)%A)
        CALL ZCSR_SIMPLIFY(AB%BK(I)%ACSR,RCUT,'N')
        GMEM_SIZE=GMEM_SIZE-REAL(SIZE(AB%BK(I)%A),gq)*16*1.E-9_gq
      ENDDO
      ! 'ZBMTOBKZCSR'
      RETURN
!
      END SUBROUTINE ZBMTOBKZCSR
!
!*********************************************************************
      SUBROUTINE ZBMTODCSR(AB,AS)
      TYPE(ZBD_MATRIX)::AB
      TYPE(DCSR_MATRIX)::AS
! LOCAL
      INTEGER I,IMIN,IMAX,NDROW,NDCOL,NROW,NCOL
      REAL(gq),ALLOCATABLE::B(:,:)
!
      CALL ZBM_CYC_BOUND(AB,IMIN,IMAX,NDROW,NDCOL)
      DO I=IMIN,IMAX
        NROW=AB%BK(I)%NROW; NCOL=AB%BK(I)%NCOL
        ALLOCATE(B(NROW,NCOL)); B=REAL(AB%BK(I)%A,gq)
        CALL DCSR_DSUM_DNS(AS,B,NROW,NCOL,'N',I==IMIN)
        DEALLOCATE(B)
      ENDDO
      CALL DCSR_ROWCOL_SHIFT(AS,NDROW,NDCOL,AB%JSHIFT)
      RETURN
!
      END SUBROUTINE ZBMTODCSR
!
!*********************************************************************
      SUBROUTINE ZBMTOZCSR(AB,AS)
      TYPE(ZBD_MATRIX)::AB
      TYPE(ZCSR_MATRIX)::AS
! LOCAL
      INTEGER I,IMIN,IMAX,NDROW,NDCOL
!
      CALL ZBM_CYC_BOUND(AB,IMIN,IMAX,NDROW,NDCOL)
      DO I=IMIN,IMAX
        IF(AB%BK(I)%LSPARSE)THEN
          CALL ZCSR_DSUM_AB(AS,AB%BK(I)%ACSR,LFIRST=(I==IMIN))
        ELSE
          CALL ZCSR_DSUM_DNS(AS,AB%BK(I)%A,AB%BK(I)%NROW,AB%BK(I)%NCOL,'N',I==IMIN)
        ENDIF
      ENDDO
      CALL ZCSR_ROWCOL_SHIFT(AS,NDROW,NDCOL,AB%JSHIFT)
      RETURN
!
      END SUBROUTINE ZBMTOZCSR
!
!*********************************************************************
      SUBROUTINE ZBM_CYC_BOUND(AB,IMIN,IMAX,NDROW,NDCOL)
      TYPE(ZBD_MATRIX)::AB
      INTEGER IMIN,IMAX,NDROW,NDCOL
!
      SELECT CASE(AB%JSHIFT)
      CASE(-1)
        IMIN=2; IMAX=AB%NBK; NDROW=AB%BK(1)%NROW; NDCOL=AB%BK(AB%NBK)%NROW
      CASE( 0)
        IMIN=1; IMAX=AB%NBK; NDROW=0; NDCOL=0
      CASE( 1)
        IMIN=1; IMAX=AB%NBK-1; NDROW=AB%BK(AB%NBK)%NCOL; NDCOL=AB%BK(1)%NCOL
      END SELECT
      RETURN
!
      END SUBROUTINE ZBM_CYC_BOUND
!
!*********************************************************************
      SUBROUTINE DBM_CYC_BOUND(AB,IMIN,IMAX,NDROW,NDCOL)
      TYPE(DBD_MATRIX)::AB
      INTEGER IMIN,IMAX,NDROW,NDCOL
!
      SELECT CASE(AB%JSHIFT)
      CASE(-1)
        IMIN=2; IMAX=AB%NBK; NDROW=AB%BK(1)%NROW; NDCOL=AB%BK(AB%NBK)%NROW
      CASE( 0)
        IMIN=1; IMAX=AB%NBK; NDROW=0; NDCOL=0
      CASE( 1)
        IMIN=1; IMAX=AB%NBK-1; NDROW=AB%BK(AB%NBK)%NCOL; NDCOL=AB%BK(1)%NCOL
      END SELECT
      RETURN
!
      END SUBROUTINE DBM_CYC_BOUND
!
!*********************************************************************
      SUBROUTINE DBM_SIMPLIFY(A,NBK,ID)
      TYPE(DBD_MATRIX)::A
      INTEGER NBK,ID(NBK+1)
! LOCAL
      INTEGER IBK_A,IBK_B,IBASE_A,IBASE_B,N1,N2,I
      REAL(gq) MAXERR
      TYPE(DBD_MATRIX)::B
!
      IF(A%JSHIFT/=0)STOP ' ERROR IN DBM_SIMPLIFY: ONLY FOR DIAGONAl BLOCK MATRIX!'
      IF(A%NBK==NBK)RETURN
      CALL ALLOC_DBM(B,NBK,ID,0,0)
      IBASE_A=0; IBASE_B=0; IBK_A=1
      DO IBK_B=1,NBK
        DO I=IBK_A,A%NBK
          IF(A%BK(I)%LZERO)CYCLE
          IF(IBASE_B>=IBASE_A+A%BK(I)%NROW)THEN
            IBASE_A=IBASE_A+A%BK(I)%NROW; CYCLE
          ENDIF
          EXIT
        ENDDO
        IBK_A=I
        N1=IBASE_B-IBASE_A+1; N2=N1+B%BK(IBK_B)%NROW-1
        IF(N2>IBASE_A+A%BK(IBK_A)%NROW)STOP ' ERROR-I IN DBM_SIMPLIFY!'
        MAXERR=MAXVAL(ABS(A%BK(IBK_A)%A(:N1-1,N1:N2)))
        MAXERR=MAX(MAXERR,MAXVAL(ABS(A%BK(IBK_A)%A(N2+1:,N1:N2))))
        MAXERR=MAX(MAXERR,MAXVAL(ABS(A%BK(IBK_A)%A(N1:N2,:N1-1))))
        MAXERR=MAX(MAXERR,MAXVAL(ABS(A%BK(IBK_A)%A(N1:N2,N2+1:))))
        IF(MAXERR>1.E-10_gq)STOP ' ERROR-II IN DBM_SIMPLIFY!'
        B%BK(IBK_B)%A=A%BK(IBK_A)%A(N1:N2,N1:N2)
        IBASE_B=IBASE_B+B%BK(IBK_B)%NROW
      ENDDO
      CALL DEALLOC_DBM(A)
      CALL DBM_COPY(B,A,1)
      CALL DEALLOC_DBM(B)
      RETURN
!
      END SUBROUTINE DBM_SIMPLIFY
!
!*********************************************************************
      SUBROUTINE ZBM_SIMPLIFY(A,NBK,ID)
      TYPE(ZBD_MATRIX)::A
      INTEGER NBK,ID(NBK+1)
! LOCAL
      INTEGER IBK_A,IBK_B,IBASE_A,IBASE_B,N1,N2,I
      REAL(gq) MAXERR
      TYPE(ZBD_MATRIX)::B
!
      IF(A%JSHIFT/=0)STOP ' ERROR IN ZBM_SIMPLIFY: ONLY FOR DIAGONAl BLOCK MATRIX!'
      IF(A%NBK==NBK)RETURN
      CALL ALLOC_ZBM(B,NBK,ID,0,0)
      IBASE_A=0; IBASE_B=0; IBK_A=1
      DO IBK_B=1,NBK
        DO I=IBK_A,A%NBK
          IF(A%BK(I)%LZERO)CYCLE
          IF(IBASE_B>=IBASE_A+A%BK(I)%NROW)THEN
            IBASE_A=IBASE_A+A%BK(I)%NROW; CYCLE
          ENDIF
          EXIT
        ENDDO
        IBK_A=I
        N1=IBASE_B-IBASE_A+1; N2=N1+B%BK(IBK_B)%NROW-1
        IF(N2>IBASE_A+A%BK(IBK_A)%NROW)STOP ' ERROR-I IN ZBM_SIMPLIFY!'
        MAXERR=MAXVAL(ABS(A%BK(IBK_A)%A(:N1-1,N1:N2)))
        MAXERR=MAX(MAXERR,MAXVAL(ABS(A%BK(IBK_A)%A(N2+1:,N1:N2))))
        MAXERR=MAX(MAXERR,MAXVAL(ABS(A%BK(IBK_A)%A(N1:N2,:N1-1))))
        MAXERR=MAX(MAXERR,MAXVAL(ABS(A%BK(IBK_A)%A(N1:N2,N2+1:))))
        IF(MAXERR>1.E-10_gq)STOP ' ERROR-II IN ZBM_SIMPLIFY!'
        B%BK(IBK_B)%A=A%BK(IBK_A)%A(N1:N2,N1:N2)
        IBASE_B=IBASE_B+B%BK(IBK_B)%NROW
      ENDDO
      CALL DEALLOC_ZBM(A)
      CALL ZBM_COPY(B,A)
      CALL DEALLOC_ZBM(B)
      RETURN
!
      END SUBROUTINE ZBM_SIMPLIFY
!
!*****************************************************************************
      SUBROUTINE DBM_COPY(A,B,MODE)
      TYPE(DBD_MATRIX)::A,B
      INTEGER MODE
! LOCAL
      INTEGER I
!
      B%NBK=A%NBK; B%DIM=A%DIM; B%JSHIFT=A%JSHIFT
      ALLOCATE(B%BK(B%NBK))
      DO I=1,B%NBK
        B%BK(I)%NROW=A%BK(I)%NROW; B%BK(I)%NCOL=A%BK(I)%NCOL
        B%BK(I)%LZERO  =A%BK(I)%LZERO
        B%BK(I)%LSPARSE=A%BK(I)%LSPARSE
        IF(B%BK(I)%LZERO)CYCLE
        ALLOCATE(B%BK(I)%A(B%BK(I)%NROW,B%BK(I)%NCOL))
        IF(MODE>0)THEN
          B%BK(I)%A=A%BK(I)%A
        ELSE
           B%BK(I)%A=0
        ENDIF
        GMEM_SIZE=GMEM_SIZE+REAL(SIZE(A%BK(I)%A),gq)*8*1.E-9_gq
      ENDDO
      RETURN
!
      END SUBROUTINE DBM_COPY
!
!*****************************************************************************
      SUBROUTINE DBM_TRACE(A,TR)
      TYPE(DBD_MATRIX)::A
      REAL(gq) TR
! LOCAL
      INTEGER IBK
      REAL(gq) TR_
!
      TR=0
      IF(A%JSHIFT/=0)RETURN
      DO IBK=1,A%NBK
        IF(A%BK(IBK)%LSPARSE)THEN
          CALL TRACE_DCSR(A%BK(IBK)%ACSR,TR_)
        ELSE
          CALL TRACE_A(A%BK(IBK)%A,A%BK(IBK)%NROW,TR_)
        ENDIF
        TR=TR+TR_
      ENDDO
      RETURN
!
      END SUBROUTINE DBM_TRACE
!
!*********************************************************************
      SUBROUTINE DBM_EIGENSYS(A,EV,MODE)
      TYPE(DBD_MATRIX)::A
      REAL(gq) EV(A%DIM)
      INTEGER MODE
! LOCAL
      INTEGER I,N1,NBASE
!
      EV=0; NBASE=1
      DO I=1,A%NBK
        IF(A%BK(I)%LZERO)THEN
          CYCLE
        ENDIF
        N1=A%BK(I)%NROW
        IF(N1/=A%BK(I)%NCOL)THEN
          WRITE(0,'(" A%BK(I)%NROW/NCOL=",2I8)')N1,A%BK(I)%NCOL
          STOP ' ERROR: A%BK(I)%A NOT SQUARE MATRIX!'
        ENDIF
        IF(MODE==-1)THEN
          A%BK(I)%A=-A%BK(I)%A
        ENDIF
        CALL HERMEV('V','L',A%BK(I)%A,EV(NBASE:NBASE+N1-1),N1)
        IF(MODE==-1)THEN
          EV(NBASE:NBASE+N1-1)=-EV(NBASE:NBASE+N1-1)
        ENDIF
        NBASE=NBASE+N1
      ENDDO
      RETURN
!
      END SUBROUTINE DBM_EIGENSYS
!
!*********************************************************************
      SUBROUTINE ZBM_EIGENSYS(A,EV,MODE)
      TYPE(ZBD_MATRIX)::A
      REAL(gq) EV(A%DIM)
      INTEGER MODE
! LOCAL
      INTEGER I,N1,NBASE
!
      EV=0; NBASE=1
      DO I=1,A%NBK
        IF(A%BK(I)%LZERO)THEN
          CYCLE
        ENDIF
        N1=A%BK(I)%NROW
        IF(N1/=A%BK(I)%NCOL)THEN
          WRITE(0,'(" A%BK(I)%NROW/NCOL=",2I8)')N1,A%BK(I)%NCOL
          STOP ' ERROR: A%BK(I)%A NOT SQUARE MATRIX!'
        ENDIF
        IF(MODE==-1)THEN
          A%BK(I)%A=-A%BK(I)%A
        ENDIF
        CALL HERMEV('V','L',A%BK(I)%A,EV(NBASE:NBASE+N1-1),N1)
        IF(MODE==-1)THEN
          EV(NBASE:NBASE+N1-1)=-EV(NBASE:NBASE+N1-1)
        ENDIF
        NBASE=NBASE+N1
      ENDDO
      RETURN
!
      END SUBROUTINE ZBM_EIGENSYS
!
!*********************************************************************
      SUBROUTINE DBM_AXB(A,B)
      TYPE(DBD_MATRIX)::A,B
! LOCAL
      INTEGER IBK_A,IBK_B,NROWB_SUM,NCOLB_SUM,NROWA,NROWB,NCOLB
      REAL(gq),ALLOCATABLE::AB(:,:)
!
      NROWB_SUM=0; NCOLB_SUM=0; IBK_A=1
      DO IBK_B=1,B%NBK
        NROWB=B%BK(IBK_B)%NROW; NCOLB=B%BK(IBK_B)%NCOL
        NROWA=A%BK(IBK_A)%NROW
        IF(NCOLB>0)THEN
          ALLOCATE(AB(NROWA,NCOLB))
          CALL DGEMM('N','N',NROWA,NCOLB,NROWB,D1,A%BK(IBK_A)%A(:,NROWB_SUM+1:NROWB_SUM+NROWB), &
                    &NROWA,B%BK(IBK_B)%A(:,1:NCOLB),NROWB,D0,AB,NROWA)
          A%BK(IBK_A)%A(:,NCOLB_SUM+1:NCOLB_SUM+NCOLB)=AB
          DEALLOCATE(AB)
        ENDIF
        NROWB_SUM=NROWB_SUM+NROWB; NCOLB_SUM=NCOLB_SUM+NCOLB
        IF(NROWB_SUM==NROWA)THEN
          A%BK(IBK_A)%NCOL=NCOLB_SUM
          IBK_A=IBK_A+1; NROWB_SUM=0; NCOLB_SUM=0
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE DBM_AXB
!
!*********************************************************************
      SUBROUTINE ZBM_AXB(A,B)
      TYPE(ZBD_MATRIX)::A,B
! LOCAL
      INTEGER IBK_A,IBK_B,NROWB_SUM,NCOLB_SUM,NROWA,NROWB,NCOLB
      COMPLEX(gq),ALLOCATABLE::AB(:,:)
!
      NROWB_SUM=0; NCOLB_SUM=0; IBK_A=1
      DO IBK_B=1,B%NBK
        NROWB=B%BK(IBK_B)%NROW; NCOLB=B%BK(IBK_B)%NCOL
        NROWA=A%BK(IBK_A)%NROW
        IF(NCOLB>0)THEN
          ALLOCATE(AB(NROWA,NCOLB))
          CALL ZGEMM('N','N',NROWA,NCOLB,NROWB,Z1,A%BK(IBK_A)%A(:,NROWB_SUM+1:NROWB_SUM+NROWB), &
                    &NROWA,B%BK(IBK_B)%A(:,1:NCOLB),NROWB,Z0,AB,NROWA)
          A%BK(IBK_A)%A(:,NCOLB_SUM+1:NCOLB_SUM+NCOLB)=AB
          DEALLOCATE(AB)
        ENDIF
        NROWB_SUM=NROWB_SUM+NROWB; NCOLB_SUM=NCOLB_SUM+NCOLB
        IF(NROWB_SUM==NROWA)THEN
          A%BK(IBK_A)%NCOL=NCOLB_SUM
          IBK_A=IBK_A+1; NROWB_SUM=0; NCOLB_SUM=0
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE ZBM_AXB
!
!*********************************************************************
      SUBROUTINE ZCSR_ROWCOL_SHIFT(AS,NDROW,NDCOL,JSHIFT)
      TYPE(ZCSR_MATRIX)::AS
      INTEGER NDROW,NDCOL,JSHIFT
! LOCAL
      INTEGER NROW
      INTEGER,ALLOCATABLE::IA(:)
!
      IF(JSHIFT==0)RETURN
      ALLOCATE(IA(AS%NROW+1)); IA=AS%I
      NROW=AS%NROW+NDROW
      DEALLOCATE(AS%I); ALLOCATE(AS%I(NROW+1))
      SELECT CASE(JSHIFT)
      CASE(-1)
        AS%I(1:NDROW)=1; AS%I(1+NDROW:)=IA
      CASE( 1)
        AS%I(1:NROW+1)=IA; AS%I(2+NROW:)=AS%I(NROW+1)
        AS%J=AS%J+NDCOL
      END SELECT
      DEALLOCATE(IA)
      AS%NROW=NROW; AS%NCOL=AS%NCOL+NDCOL
      RETURN
!
      END SUBROUTINE ZCSR_ROWCOL_SHIFT
!
!*********************************************************************
      SUBROUTINE DCSR_ROWCOL_SHIFT(AS,NDROW,NDCOL,JSHIFT)
      TYPE(DCSR_MATRIX)::AS
      INTEGER NDROW,NDCOL,JSHIFT
! LOCAL
      INTEGER NROW
      INTEGER,ALLOCATABLE::IA(:)
!
      IF(JSHIFT==0)RETURN
      ALLOCATE(IA(AS%NROW+1)); IA=AS%I
      NROW=AS%NROW+NDROW
      DEALLOCATE(AS%I); ALLOCATE(AS%I(NROW+1))
      SELECT CASE(JSHIFT)
      CASE(-1)
        AS%I(1:NDROW)=1; AS%I(1+NDROW:)=IA
      CASE( 1)
        AS%I(1:NROW+1)=IA; AS%I(2+NROW:)=AS%I(NROW+1)
        AS%J=AS%J+NDCOL
      END SELECT
      DEALLOCATE(IA)
      AS%NROW=NROW; AS%NCOL=AS%NCOL+NDCOL
      RETURN
!
      END SUBROUTINE DCSR_ROWCOL_SHIFT
!
!*********************************************************************
      SUBROUTINE ZCSRTOZBM(AS,AB,NBK,ID,JSHIFT)
      TYPE(ZBD_MATRIX)::AB
      TYPE(ZCSR_MATRIX)::AS
      INTEGER NBK,ID(NBK+1),JSHIFT
! LOCAL
      INTEGER I,N1,N2,M1,M2,IBASE
!
      CALL ALLOC_ZBM(AB,NBK,ID,JSHIFT,0)
      IBASE=ID(1)-1
      DO I=1,NBK
        IF(AB%BK(I)%LZERO)CYCLE
        N1=ID(I)-IBASE;        N2=ID(I+1)-1-IBASE
        M1=ID(I+JSHIFT)-IBASE; M2=ID(I+JSHIFT+1)-1-IBASE
        CALL ZCSR_DIADNS_SUBMAT(AS,AB%BK(I)%A,N1,N2,M1,M2)
      ENDDO
      RETURN
!
      END SUBROUTINE ZCSRTOZBM
!
!*********************************************************************
      SUBROUTINE ZCSRTODBM(AS,AB,NBK,ID,JSHIFT)
      TYPE(DBD_MATRIX)::AB
      TYPE(ZCSR_MATRIX)::AS
      INTEGER NBK,ID(NBK+1),JSHIFT
! LOCAL
      INTEGER I,N1,N2,M1,M2,IBASE
      REAL(gq) MAXIMAG
      COMPLEX(gq),ALLOCATABLE::ABLK(:,:)
!
      CALL ALLOC_DBM(AB,NBK,ID,JSHIFT,0)
      IBASE=ID(1)-1
      DO I=1,NBK
        IF(AB%BK(I)%LZERO)CYCLE
        N1=ID(I)-IBASE;        N2=ID(I+1)-1-IBASE
        M1=ID(I+JSHIFT)-IBASE; M2=ID(I+JSHIFT+1)-1-IBASE
        ALLOCATE(ABLK(N2-N1+1,M2-M1+1)); ABLK=0
        CALL ZCSR_DIADNS_SUBMAT(AS,ABLK,N1,N2,M1,M2)
        MAXIMAG=MAXVAL(ABS(AIMAG(ABLK)))
        IF(MAXIMAG>1.E-10_gq)THEN
          STOP ' ERROR: BD_MATRIX HAS IMAGINARY PART, USE COMPLEX VERSION!'
        ENDIF
        AB%BK(I)%A=REAL(ABLK,gq)
        DEALLOCATE(ABLK)
      ENDDO
      RETURN
!
      END SUBROUTINE ZCSRTODBM
!
!*****************************************************************************
! MODE=0: LSPARSE=.FALSE.
!     >0: LSPARSE=.TRUE., NNZ given
!*****************************************************************************
      SUBROUTINE ALLOC_ZBM(A,NBK,ID,JSHIFT,MODE)
      TYPE(ZBD_MATRIX)::A
      INTEGER NBK,ID(NBK+1),JSHIFT,MODE
! LOCAL
      INTEGER I,N1,N2,M1,M2,IMIN,IMAX,NROW,NCOL,NNZ
!
      !(*,'(" NBK=",I4)')NBK
      A%NBK=NBK; A%JSHIFT=JSHIFT; A%DIM=ID(NBK+1)-ID(1)
      SELECT CASE(JSHIFT)
      CASE(-1); IMIN=2; IMAX=NBK
      CASE( 0); IMIN=1; IMAX=NBK
      CASE( 1); IMIN=1; IMAX=NBK-1
      END SELECT
      ALLOCATE(A%BK(NBK))
      DO I=1,NBK
        N1=ID(I); N2=ID(I+1)-1; NROW=N2-N1+1
        A%BK(I)%NROW=NROW; A%BK(I)%LZERO=.FALSE.
        IF(MODE>0)THEN
          A%BK(I)%LSPARSE=.TRUE.
        ELSE
          A%BK(I)%LSPARSE=.FALSE.
        ENDIF
        IF(I<IMIN.OR.I>IMAX)THEN
          A%BK(I)%NCOL=NROW; A%BK(I)%LZERO=.TRUE.; CYCLE
        ENDIF
        M1=ID(I+JSHIFT); M2=ID(I+JSHIFT+1)-1; NCOL=M2-M1+1
        A%BK(I)%NCOL=NCOL
        IF(MODE>0)THEN
          NNZ=MIN(NROW,NCOL)*MODE
          CALL ALLOC_ZCSR(A%BK(I)%ACSR,NNZ,NROW,NCOL)
        ELSE
          ALLOCATE(A%BK(I)%A(NROW,NCOL)); A%BK(I)%A=0
          GMEM_SIZE=GMEM_SIZE+REAL(SIZE(A%BK(I)%A),gq)*16*1.E-9_gq
        ENDIF
      ENDDO
      ! 'ALLOC_ZBM'
      RETURN
!
      END SUBROUTINE ALLOC_ZBM
!
!*****************************************************************************
! MODE=0: LSPARSE=.FALSE.
!     >0: LSPARSE=.TRUE., NNZ given
!*****************************************************************************
      SUBROUTINE ALLOC_DBM(A,NBK,ID,JSHIFT,MODE)
      TYPE(DBD_MATRIX)::A
      INTEGER NBK,ID(NBK+1),JSHIFT,MODE
! LOCAL
      INTEGER I,N1,N2,M1,M2,IMIN,IMAX,NROW,NCOL,NNZ
!
      !(*,'(" NBK=",I4)')NBK
      A%NBK=NBK; A%JSHIFT=JSHIFT; A%DIM=ID(NBK+1)-ID(1)
      SELECT CASE(JSHIFT)
      CASE(-1); IMIN=2; IMAX=NBK
      CASE( 0); IMIN=1; IMAX=NBK
      CASE( 1); IMIN=1; IMAX=NBK-1
      END SELECT
      ALLOCATE(A%BK(NBK))
      DO I=1,NBK
        N1=ID(I); N2=ID(I+1)-1; NROW=N2-N1+1
        A%BK(I)%NROW=NROW; A%BK(I)%LZERO=.FALSE.
        IF(MODE>0)THEN
          A%BK(I)%LSPARSE=.TRUE.
        ELSE
          A%BK(I)%LSPARSE=.FALSE.
        ENDIF
        IF(I<IMIN.OR.I>IMAX)THEN
          A%BK(I)%NCOL=NROW; A%BK(I)%LZERO=.TRUE.; CYCLE
        ENDIF
        M1=ID(I+JSHIFT); M2=ID(I+JSHIFT+1)-1; NCOL=M2-M1+1
        A%BK(I)%NCOL=NCOL
        IF(NROW==0.OR.NCOL==0)THEN
          A%BK(I)%LZERO=.TRUE.; CYCLE
        ENDIF
        IF(MODE>0)THEN
          NNZ=MIN(NROW,NCOL)*MODE
          CALL ALLOC_DCSR(A%BK(I)%ACSR,NNZ,NROW,NCOL)
        ELSE
          ALLOCATE(A%BK(I)%A(NROW,NCOL)); A%BK(I)%A=0
          GMEM_SIZE=GMEM_SIZE+REAL(SIZE(A%BK(I)%A),gq)*8*1.E-9_gq
        ENDIF
      ENDDO
      ! 'ALLOC_DBM'
      RETURN
!
      END SUBROUTINE ALLOC_DBM
!
!*****************************************************************************
      SUBROUTINE DEALLOC_ZBM(A)
      TYPE(ZBD_MATRIX)::A
! LOCAL
      INTEGER I
!
      DO I=1,A%NBK
        IF(A%BK(I)%LZERO)CYCLE
        IF(A%BK(I)%LSPARSE)THEN
          CALL DEALLOC_ZCSR(A%BK(I)%ACSR)
        ELSE
          GMEM_SIZE=GMEM_SIZE-REAL(SIZE(A%BK(I)%A),gq)*16*1.E-9_gq
          DEALLOCATE(A%BK(I)%A)
        ENDIF
      ENDDO
      A%NBK=0; A%DIM=0
      DEALLOCATE(A%BK)
      ! 'DEALLOC_ZBM'
      RETURN
!
      END SUBROUTINE DEALLOC_ZBM
!
!*****************************************************************************
      SUBROUTINE DEALLOC_DBM(A)
      TYPE(DBD_MATRIX)::A
! LOCAL
      INTEGER I
!
      DO I=1,A%NBK
        IF(A%BK(I)%LZERO)CYCLE
        IF(A%BK(I)%LSPARSE)THEN
          CALL DEALLOC_DCSR(A%BK(I)%ACSR)
        ELSE
          GMEM_SIZE=GMEM_SIZE-REAL(SIZE(A%BK(I)%A),gq)*8*1.E-9_gq
          DEALLOCATE(A%BK(I)%A)
        ENDIF
      ENDDO
      A%NBK=0; A%DIM=0
      DEALLOCATE(A%BK)
      ! 'DEALLOC_DBM'
      RETURN
!
      END SUBROUTINE DEALLOC_DBM
!
!*****************************************************************************
      SUBROUTINE SET_IDEN_DBM(A)
      TYPE(DBD_MATRIX)::A
! LOCAL
      INTEGER I,J
!
      IF(A%JSHIFT.NE.0)THEN
        STOP 'FETAL ERROR: OFF-DIAGONAL MATRIX CANNOT BE SET IDENTITY!'
      ENDIF
      DO I=1,A%NBK
        IF(A%BK(I)%LZERO)CYCLE
        DO J=1,A%BK(I)%NROW
          A%BK(I)%A(J,J)=1._gq
        ENDDO
      ENDDO
      RETURN
!
      END SUBROUTINE SET_IDEN_DBM
!
!*****************************************************************************
      SUBROUTINE SET_IDEN_ZBM(A)
      TYPE(ZBD_MATRIX)::A
! LOCAL
      INTEGER I,J
!
      IF(A%JSHIFT.NE.0)THEN
        STOP 'FETAL ERROR: OFF-DIAGONAL MATRIX CANNOT BE SET IDENTITY!'
      ENDIF
      DO I=1,A%NBK
        IF(A%BK(I)%LZERO)CYCLE
        DO J=1,A%BK(I)%NROW
          A%BK(I)%A(J,J)=1._gq
        ENDDO
      ENDDO
      RETURN
!
      END SUBROUTINE SET_IDEN_ZBM
!
!*****************************************************************************
      SUBROUTINE ZBM_COPY(A,B)
      TYPE(ZBD_MATRIX)::A,B
! LOCAL
      INTEGER I,NBK,NROW,NCOL
!
      NBK=A%NBK
      B%NBK=NBK; B%JSHIFT=A%JSHIFT; B%DIM=A%DIM
      ALLOCATE(B%BK(NBK))
      DO I=1,NBK
        NROW=A%BK(I)%NROW; NCOL=A%BK(I)%NCOL
        B%BK(I)%NROW=NROW; B%BK(I)%NCOL=NCOL
        B%BK(I)%LZERO  =A%BK(I)%LZERO
        B%BK(I)%LSPARSE=A%BK(I)%LSPARSE
        IF(A%BK(I)%LZERO)CYCLE
        ALLOCATE(B%BK(I)%A(NROW,NCOL))
        B%BK(I)%A=A%BK(I)%A
        GMEM_SIZE=GMEM_SIZE+REAL(SIZE(B%BK(I)%A),gq)*16*1.E-9_gq
      ENDDO
      ! 'ZBM_COPY'
      RETURN
!
      END SUBROUTINE ZBM_COPY
!
!*****************************************************************************
      SUBROUTINE ZBM_COPY1(A,B)
      TYPE(ZBD_MATRIX)::A,B
! LOCAL
      INTEGER I
!
      DO I=1,A%NBK
        IF(A%BK(I)%LZERO)CYCLE
        B%BK(I)%A=A%BK(I)%A
      ENDDO
      RETURN
!
      END SUBROUTINE ZBM_COPY1
!
!*****************************************************************************
      SUBROUTINE ZBM_LINK(A,B,IBK)
      TYPE(ZBD_MATRIX)::A,B
      INTEGER IBK
! LOCAL
      INTEGER IBASE
!
      IBASE=IBK-1
      B%NBK=A%NBK-IBASE; B%JSHIFT=A%JSHIFT; B%DIM=A%DIM-SUM(A%BK(1:IBASE)%NROW)
      B%BK=>A%BK(IBASE+1:)
      RETURN
!
      END SUBROUTINE ZBM_LINK
!
!*****************************************************************************
      SUBROUTINE ZBM_DUMP(A,NI,IA1,IA2,MODE)
      TYPE(ZBD_MATRIX),INTENT(IN)::A
      INTEGER NI,IA1,IA2,MODE
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,IA1,IA2,MODE))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      CALL ZBM_WRT(A,GIU)
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE ZBM_DUMP
!
!*****************************************************************************
      SUBROUTINE ZBM2_DUMP(A,N,NI,MODE)
      INTEGER N,NI,MODE
      TYPE(ZBD_MATRIX),INTENT(IN)::A(N,N)
! LOCAL
      INTEGER IA1,IA2
!
      DO IA1=1,N; DO IA2=1,N
        IF(A(IA1,IA2)%DIM<=0)CYCLE  ! No use.
        CALL ZBM_DUMP(A(IA1,IA2),NI,IA1,IA2,MODE)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE ZBM2_DUMP
!
!*****************************************************************************
      SUBROUTINE ZBM_LOAD(A,NI,IA1,IA2,MODE)
      TYPE(ZBD_MATRIX)::A
      INTEGER NI,IA1,IA2,MODE,INFO
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,IA1,IA2,MODE))),STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL",ERR=100)
      INFO=0
      CALL ZBM_READ(A,GIU,INFO)
      IF(INFO/=0)THEN
        REWIND(GIU); CALL ZBM_READ(A,GIU,INFO)
      ENDIF
      CLOSE(GIU)
      RETURN
100   CONTINUE
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,IA1,IA2,MODE+100))),STATUS='OLD')
      CALL ZBM_READ_TXT(A,GIU)
      CLOSE(GIU)
      CALL ZBM_DUMP(A,NI,IA1,IA2,MODE)
      RETURN
!
      END SUBROUTINE ZBM_LOAD
!
!*****************************************************************************
      SUBROUTINE ZBM2_LOAD(A,N,NI,MODE)
      INTEGER N,NI,MODE
      TYPE(ZBD_MATRIX)::A(N,N)
! LOCAL
      INTEGER IA1,IA2
!
      DO IA1=1,N; DO IA2=1,N
        IF(A(IA1,IA2)%DIM<=0)CYCLE
        CALL ZBM_LOAD(A(IA1,IA2),NI,IA1,IA2,MODE); 
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE ZBM2_LOAD
!
!*****************************************************************************
      SUBROUTINE ZBM_UHAU(A,U)
      TYPE(ZBD_MATRIX)::A,U
! LOCAL
      INTEGER I,NROW,NCOL
      COMPLEX(gq),ALLOCATABLE::AU(:,:)
      !
!
      ! "ZBM_UHAU"
      IF(A%JSHIFT/=0)THEN
        STOP ' ERROR IN ZBM_UHAU: A%JSHIFT/=0!'
      ENDIF
      IF(U%JSHIFT/=0)THEN
        STOP ' ERROR IN ZBM_UHAU: U%JSHIFT/=0!'
      ENDIF
      IF(A%NBK/=U%NBK)THEN
        STOP ' ERROR IN ZBM_UHAU: A%NBK/=U%NBK!'
      ENDIF
      DO I=1,U%NBK
        NROW=U%BK(I)%NROW; NCOL=U%BK(I)%NCOL
        IF(NROW/=A%BK(I)%NROW)THEN
          STOP ' ERROR IN ZBM_UHAU: A%BK(I)%NROW/=U%BK(I)%NROW!'
        ENDIF
        IF(NROW/=A%BK(I)%NCOL)THEN
          STOP ' ERROR IN ZBM_UHAU: A%BK(I)%NCOL/=U%BK(I)%NROW!'
        ENDIF
        IF(A%BK(I)%LZERO)CYCLE
        IF(NCOL==0.OR.NROW==0)THEN
          A%BK(I)%LZERO=.TRUE.; A%BK(I)%NROW=0; A%BK(I)%NCOL=0
          DEALLOCATE(A%BK(I)%A); CYCLE
        ENDIF
        !(0,'(" IBK=",I3," NROW/NCOL=",2I6)')I,NROW,NCOL
        ALLOCATE(AU(NROW,NCOL)); AU=0
        !
        CALL ZGEMM('N','N',NROW,NCOL,NROW,Z1,A%BK(I)%A,NROW,U%BK(I)%A,NROW,Z0,AU,NROW)
        IF(NROW/=NCOL)THEN
          DEALLOCATE(A%BK(I)%A); ALLOCATE(A%BK(I)%A(NCOL,NCOL))
          A%BK(I)%NROW=NCOL; A%BK(I)%NCOL=NCOL
        ENDIF
        CALL ZGEMM('C','N',NCOL,NCOL,NROW,Z1,U%BK(I)%A,NROW,AU,NROW,Z0,A%BK(I)%A,NCOL)
        !
        DEALLOCATE(AU)
      ENDDO
      ! "ZBM_UHAU"
      RETURN
!
      END SUBROUTINE ZBM_UHAU
!
!*****************************************************************************
      SUBROUTINE ZBM_GEMM(TRANSA,TRANSB,A,B,LB)
      CHARACTER*1 TRANSA,TRANSB
      TYPE(ZBD_MATRIX)::A,B
      INTEGER,OPTIONAL::LB
! LOCAL
      INTEGER I,NROW,NCOL
      COMPLEX(gq),ALLOCATABLE::AB(:,:)
!
      IF(A%JSHIFT/=0.OR.B%JSHIFT/=0)STOP 'ERROR IN ZBM_GEMM: U%JSHIFT/=0!'
      DO I=1,A%NBK
        NROW=A%BK(I)%NROW; NCOL=A%BK(I)%NCOL
        IF(NROW<=0.OR.NCOL<=0)CYCLE
        ALLOCATE(AB(NROW,NCOL)); AB=0
        CALL ZGEMM(TRANSA,TRANSB,NROW,NCOL,NCOL,Z1,A%BK(I)%A,NROW,B%BK(I)%A,NROW,Z0,AB,NROW) 
        IF(PRESENT(LB))THEN
          B%BK(I)%A=AB
        ELSE
          A%BK(I)%A=AB
        ENDIF
        DEALLOCATE(AB)
      ENDDO
      RETURN      
!
      END SUBROUTINE ZBM_GEMM
!
!*****************************************************************************
      SUBROUTINE DBM_TR_RHOA(RHO,A,VAL)
      TYPE(DBD_MATRIX)::RHO,A
      REAL(gq) VAL
! LOCAL
      INTEGER IBK_R,IBK_A,IDIF,ISUM_A,IBASE_R,I,NROW_R
!
      VAL=0
      IBK_A=1; ISUM_A=A%BK(1)%NROW+1; IBASE_R=1
      DO IBK_R=1,RHO%NBK
      IF(IBASE_R>=ISUM_A)THEN
        IBK_A=IBK_A+1
        ISUM_A=ISUM_A+A%BK(IBK_A)%NROW
      ENDIF
      IDIF=IBASE_R-(ISUM_A-A%BK(IBK_A)%NROW)
      NROW_R=RHO%BK(IBK_R)%NROW
      DO I=1,NROW_R
        VAL=VAL+SUM(RHO%BK(IBK_R)%A(I,:)*A%BK(IBK_A)%A(IDIF+1:IDIF+NROW_R,IDIF+I))
      ENDDO
      IBASE_R=IBASE_R+NROW_R
      ENDDO
      RETURN
!
      END SUBROUTINE DBM_TR_RHOA
!
!*****************************************************************************
      SUBROUTINE ZBM_TR_RHOA(RHO,A,VAL)
      TYPE(ZBD_MATRIX)::RHO,A
      COMPLEX(gq) VAL
! LOCAL
      INTEGER IBK_R,IBK_A,IDIF,ISUM_A,IBASE_R,I,NROW_R
!
      VAL=0
      IBK_A=1; ISUM_A=A%BK(1)%NROW+1; IBASE_R=1
      DO IBK_R=1,RHO%NBK
      IF(IBASE_R>=ISUM_A)THEN
        IBK_A=IBK_A+1
        ISUM_A=ISUM_A+A%BK(IBK_A)%NROW
      ENDIF
      IDIF=IBASE_R-(ISUM_A-A%BK(IBK_A)%NROW)
      NROW_R=RHO%BK(IBK_R)%NROW
      DO I=1,NROW_R
        VAL=VAL+SUM(RHO%BK(IBK_R)%A(I,:)*A%BK(IBK_A)%A(IDIF+1:IDIF+NROW_R,IDIF+I))
      ENDDO
      IBASE_R=IBASE_R+NROW_R
      ENDDO
      RETURN
!
      END SUBROUTINE ZBM_TR_RHOA
!
!*****************************************************************************
      SUBROUTINE ZBM_AMUX(A,X,LHM)
      TYPE(ZBD_MATRIX)::A
      COMPLEX(gq),TARGET::X(A%DIM)
      LOGICAL LHM ! Hermitian matrix
! LOCAL
      INTEGER IBK,IBASE,NROW,NCOL,I1,I2,J1,J2
      COMPLEX(gq),ALLOCATABLE::AX(:)
!
      ALLOCATE(AX(A%DIM)); AX=0
      IBASE=0
      DO IBK=1,A%NBK
      NROW=A%BK(IBK)%NROW; NCOL=A%BK(IBK)%NCOL
      IF(A%BK(IBK)%LZERO)GOTO 100
      I1=IBASE+1; I2=I1+NROW-1
      IF(I1>I2)GOTO 100
      J1=IBASE+A%JSHIFT*NCOL+1; J2=J1+NCOL-1
      IF(J1>J2.OR.J1<1.OR.J1>A%DIM)GOTO 100
      IF(A%BK(IBK)%LSPARSE)THEN
        CALL ZCSR_GEAMUX_SK('N',A%BK(IBK)%ACSR,X(J1:J2),AX(I1:I2))
        IF(LHM.AND.A%JSHIFT/=0)THEN
        CALL ZCSR_GEAMUX_SK('C',A%BK(IBK)%ACSR,X(I1:I2),AX(J1:J2))
        ENDIF
      ELSE
        CALL ZGEMV('N',NROW,NCOL,Z1,A%BK(IBK)%A,NROW,X(J1:J2),1,Z1,AX(I1:I2),1) ! M
        IF(LHM.AND.A%JSHIFT/=0)THEN
        CALL ZGEMV('C',NROW,NCOL,Z1,A%BK(IBK)%A,NROW,X(I1:I2),1,Z1,AX(J1:J2),1) ! M^H
        ENDIF
      ENDIF
100   CONTINUE
      IBASE=IBASE+NROW
      ENDDO
      X=AX
      DEALLOCATE(AX)
      RETURN
!
      END SUBROUTINE ZBM_AMUX
!
!*****************************************************************************
      SUBROUTINE ZBM_TR_RHOUHVU(RHO,U,V,VAL)
      TYPE(ZBD_MATRIX)::RHO,U
      TYPE(ZCSR_MATRIX) :: V
      COMPLEX(gq) VAL
! LOCAL
      INTEGER IBK_R,IBK_A,IDIF,IBASE_AL,IBASE_AR,IBASE_R,I,NROW_R,NROW,NCOL,N1,N2
      LOGICAL LDO
      TYPE(ZCSR_MATRIX) :: NCSR
      COMPLEX(gq),ALLOCATABLE::UHVU(:,:),VU(:,:)
!
      VAL=0
      IBK_A=1; IBASE_AL=1; IBASE_AR=1; IBASE_R=1; LDO=.TRUE.
      DO IBK_R=1,RHO%NBK
      IF(IBASE_R>=IBASE_AR+U%BK(IBK_A)%NCOL)THEN
        DO I=IBK_A,U%NBK-1
          IBASE_AL=IBASE_AL+U%BK(I)%NROW
          IBASE_AR=IBASE_AR+U%BK(I)%NCOL
          IF(U%BK(I+1)%LZERO)THEN
            CYCLE
          ELSE
            IBK_A=I+1
            IF(IBASE_R>=IBASE_AR+U%BK(IBK_A)%NCOL)THEN
              STOP ' ERROR IN ZBM_TR_RHOUHVU: RHO, U NOT MATCH!'
            ENDIF
            EXIT
          ENDIF
        ENDDO
        LDO=.TRUE.
      ENDIF
      IF(LDO)THEN
        NROW=U%BK(IBK_A)%NROW; NCOL=U%BK(IBK_A)%NCOL
        ALLOCATE(VU(NROW,NCOL))
        N1=IBASE_AL; N2=IBASE_AL+NROW-1
        CALL ZCSR_DIA_SUBMAT(V,NCSR,N1,N2,N1,N2)
        CALL ZCSRMUDEN_SK(NCSR,U%BK(IBK_A)%A,VU,NCOL)
        CALL DEALLOC_ZCSR(NCSR)
        IF(IBK_A>1) DEALLOCATE(UHVU)
        ALLOCATE(UHVU(NCOL,NCOL)); UHVU=0
        CALL ZGEMM('C','N',NCOL,NCOL,NROW,Z1,U%BK(IBK_A)%A,NROW,VU,NROW,Z0,UHVU,NCOL)
        DEALLOCATE(VU)
        LDO=.FALSE.
      ENDIF
      IDIF=IBASE_R-IBASE_AR
      NROW_R=RHO%BK(IBK_R)%NROW
      DO I=1,NROW_R
        VAL=VAL+SUM(RHO%BK(IBK_R)%A(I,:)*UHVU(IDIF+1:IDIF+NROW_R,IDIF+I))
      ENDDO
      IBASE_R=IBASE_R+NROW_R
      ENDDO
      DEALLOCATE(UHVU)
      RETURN
!
      END SUBROUTINE ZBM_TR_RHOUHVU
!
!*****************************************************************************
! MODE>0: A->ABLK; MODE<=0: ABLK->A
!*****************************************************************************
      SUBROUTINE ZBM_DIA_SUBMAT(A,ABLK,N1,N2,M1,M2,MODE)
      INTEGER N1,M1,N2,M2
      TYPE(ZBD_MATRIX)::A
      COMPLEX(gq)ABLK(N2-N1+1,M2-M1+1)
      INTEGER MODE
! LOCAL
      INTEGER I,J1,J2,IBASE,NROW,N1_,N2_,M1_,M2_
      REAL(gq) ERR
      REAL(gq),PARAMETER::RCUT=1.E-12_gq
!
      IF(A%JSHIFT/=0)THEN
        STOP 'ERROR IN ZBM_DIA_SUBMAT: U%JSHIFT/=0!'
      ENDIF
      IBASE=0; ERR=0
      DO I=1,A%NBK
        IF(MODE==-2)THEN
          NROW=A%BK(I)%NCOL  ! ANMXBMM
        ELSE
          NROW=A%BK(I)%NROW
        ENDIF
        N1_=N1-IBASE; N2_=N2-IBASE; M1_=M1-IBASE; M2_=M2-IBASE
        IF(N1_>NROW)THEN
          IBASE=IBASE+NROW; CYCLE
        ENDIF
        IF(N2_>NROW)EXIT
        IF(MODE==1.OR.MODE==2)THEN
          ABLK=A%BK(I)%A(N1_:N2_,M1_:M2_)
          IF(MODE==1)THEN
          DO J1=N1_,N2_; DO J2=1,NROW
            IF(J2>=N1_.AND.J2<=N2_)CYCLE
            ERR=MAX(ERR,ABS(A%BK(I)%A(J1,J2)))
            ERR=MAX(ERR,ABS(A%BK(I)%A(J2,J1)))
          ENDDO; ENDDO
          IF(ERR.GT.RCUT)THEN
            WRITE(0,'(" ERROR IN ZBM_GET_DIA_SUBMAT--NOT BLOCK DIAGONAL! ERR=",E10.2)')ERR
            STOP
          ENDIF
          ENDIF
        ELSEIF(MODE==-1)THEN
          A%BK(I)%A(N1_:N2_,M1_:M2_)=ABLK
        ELSEIF(MODE==-2)THEN
          CALL ANMXBMM('N',A%BK(I)%A(:,M1_:M2_),ABLK,A%BK(I)%NROW,N2-N1+1)
        ELSE
          STOP 'ERROR INI ZBM_DIA_SUBMAT: ILLEGAL MODE!'
        ENDIF
        GOTO 100
      ENDDO
      STOP 'ERROR IN ZBM_DIA_SUBMAT: SUB_DIAMAT NOT IN A SINGLE BLOCK!'
100   CONTINUE
      RETURN
!
      END SUBROUTINE ZBM_DIA_SUBMAT
!
!*****************************************************************************
      SUBROUTINE ZBM_SUBMAT_IJBASE(A,N1,N2,M1,M2,IBK,IRBASE,JCBASE)
      INTEGER N1,M1,N2,M2,IRBASE,JCBASE,IBK
      TYPE(ZBD_MATRIX)::A
! LOCAL
      INTEGER I,NROW
!
      IRBASE=0
      DO I=1,A%NBK
        NROW=A%BK(I)%NROW
        IF(IRBASE+NROW<N1)THEN
          IRBASE=IRBASE+NROW; CYCLE
        ENDIF
        IF(IRBASE+NROW<N2)THEN
          STOP ' ERROR IN ZBM_SUBMAT_IDX: SUBMAT NOT IN BLOCK!'
        ENDIF
        SELECT CASE(A%JSHIFT)
        CASE(-1)
          JCBASE=IRBASE-A%BK(I-1)%NROW
          IF(M2>JCBASE+A%BK(I-1)%NROW)THEN
            STOP ' ERROR IN ZBM_SUBMAT_IDX: M2>JCBASE+A%BK(I-1)%NROW!'
          ENDIF
        CASE(0)
          JCBASE=IRBASE
          IF(M2>JCBASE+NROW)THEN
            STOP ' ERROR IN ZBM_SUBMAT_IDX: M2>JCBASE+NROW!'
          ENDIF
        CASE(1)
          JCBASE=IRBASE+NROW
          IF(M2>JCBASE+JCBASE+A%BK(I+1)%NROW)THEN
            STOP ' ERROR IN ZBM_SUBMAT_IDX: M2>JCBASE+JCBASE+A%BK(I+1)%NROW!'
          ENDIF
        END SELECT
        IF(M1<=JCBASE)THEN
          STOP ' ERROR IN ZBM_SUBMAT_IDX: M1<=JCBASE!'
        ENDIF
        IBK=I
        RETURN
      ENDDO
      STOP 'ERROR IN ZBM_SUBMAT_IDX: CAN NOT LOCATE BLOCK!'
      RETURN
!
      END SUBROUTINE ZBM_SUBMAT_IJBASE
!
!*********************************************************************
      SUBROUTINE IBM_WRT(A,IU)
      TYPE(IBD_MATRIX)::A
      INTEGER IU
! LOCAL
      INTEGER I
!
      WRITE(IU)A%NBK,A%DIM
      DO I=1,A%NBK
        WRITE(IU)A%BK(I)%NROW,A%BK(I)%NCOL
        WRITE(IU)A%BK(I)%A(:,1:A%BK(I)%NCOL)
      ENDDO
      RETURN
!
      END SUBROUTINE IBM_WRT
!
!*********************************************************************
      SUBROUTINE DBM_WRT(A,IU)
      TYPE(DBD_MATRIX)::A
      INTEGER IU
! LOCAL
      INTEGER I
!
      WRITE(IU)A%NBK,A%DIM,A%JSHIFT
      DO I=1,A%NBK
        WRITE(IU)A%BK(I)%NROW,A%BK(I)%NCOL,A%BK(I)%LZERO,A%BK(I)%LSPARSE
        IF(A%BK(I)%LZERO)CYCLE
        IF(A%BK(I)%LSPARSE)THEN
          CALL DCSR_WRT(A%BK(I)%ACSR,IU)
        ELSE
          WRITE(IU)A%BK(I)%A(:,1:A%BK(I)%NCOL)
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE DBM_WRT
!
!*********************************************************************
      SUBROUTINE ZBM_WRT(A,IU)
      TYPE(ZBD_MATRIX)::A
      INTEGER IU
! LOCAL
      INTEGER I
!
      WRITE(IU)A%NBK,A%DIM,A%JSHIFT
      DO I=1,A%NBK
        WRITE(IU)A%BK(I)%NROW,A%BK(I)%NCOL,A%BK(I)%LZERO,A%BK(I)%LSPARSE
        IF(A%BK(I)%LZERO)CYCLE
        IF(A%BK(I)%LSPARSE)THEN
          CALL ZCSR_WRT(A%BK(I)%ACSR,IU)
        ELSE
          WRITE(IU)A%BK(I)%A(:,1:A%BK(I)%NCOL)
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE ZBM_WRT
!
!*********************************************************************
      SUBROUTINE DBM_WRT2(A,IU)
      TYPE(DBD_MATRIX)::A
      INTEGER IU
! LOCAL
      CHARACTER(50) FMT
      INTEGER I,J
!
      WRITE(IU,'(" A%NBK=",I3," A%DIM=",I10," A%JSHIFT=",I2)')A%NBK,A%DIM,A%JSHIFT
      DO I=1,A%NBK
        WRITE(IU,'(" BK_NROW=",I8," NCOL=",I8," LZERO=",L2," LSPARSE=",L2)')A%BK(I)%NROW,A%BK(I)%NCOL,A%BK(I)%LZERO,A%BK(I)%LSPARSE
        IF(A%BK(I)%LZERO)CYCLE
        IF(.NOT.A%BK(I)%LSPARSE)THEN
          WRITE(FMT,'("(",I0,"F7.3)")')A%BK(I)%NCOL
          WRITE(IU,FMT)(A%BK(I)%A(J,1:A%BK(I)%NCOL),J=1,A%BK(I)%NROW)
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE DBM_WRT2
!
!*********************************************************************
      SUBROUTINE ZBM_WRT2(A,IU)
      TYPE(ZBD_MATRIX)::A
      INTEGER IU
! LOCAL
      CHARACTER(50) FMT
      INTEGER I,J
!
      WRITE(IU,'(" A%NBK=",I3," A%DIM=",I10," A%JSHIFT=",I2)')A%NBK,A%DIM,A%JSHIFT
      DO I=1,A%NBK
        WRITE(IU,'(" BK_NROW=",I8," NCOL=",I8," LZERO=",L2," LSPARSE=",L2)')A%BK(I)%NROW,A%BK(I)%NCOL,A%BK(I)%LZERO,A%BK(I)%LSPARSE
        IF(A%BK(I)%LZERO)CYCLE
        IF(.NOT.A%BK(I)%LSPARSE)THEN
          WRITE(FMT,'("(",I0,"F7.3)")')A%BK(I)%NCOL
          WRITE(IU,'(" REAL PART")')
          WRITE(IU,FMT)(REAL(A%BK(I)%A(J,:),gq),J=1,A%BK(I)%NROW)
          WRITE(IU,'(" IMAG PART")')
          WRITE(IU,FMT)(AIMAG(A%BK(I)%A(J,:)),J=1,A%BK(I)%NROW)
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE ZBM_WRT2
!
!*********************************************************************
      SUBROUTINE IBM_READ(A,IU)
      TYPE(IBD_MATRIX)::A
      INTEGER IU
! LOCAL
      INTEGER I
!
      READ(IU)A%NBK,A%DIM
      ALLOCATE(A%BK(A%NBK))
      DO I=1,A%NBK
        READ(IU)A%BK(I)%NROW,A%BK(I)%NCOL
        ALLOCATE(A%BK(I)%A(A%BK(I)%NROW,A%BK(I)%NCOL))
        READ(IU)A%BK(I)%A
        GMEM_SIZE=GMEM_SIZE+REAL(SIZE(A%BK(I)%A),gq)*8*1.E-9_gq
      ENDDO
      ! 'IBM_READ'
      RETURN
      RETURN
!
      END SUBROUTINE IBM_READ
!
!*********************************************************************
      SUBROUTINE DBM_READ(A,IU,INFO)
      TYPE(DBD_MATRIX)::A
      INTEGER IU,INFO
! LOCAL
      INTEGER I
!
      READ(IU)A%NBK,A%DIM,A%JSHIFT
      IF(INFO==0)THEN
        ALLOCATE(A%BK(A%NBK))
      ENDIF
      DO I=1,A%NBK
        IF(INFO==0)THEN
          READ(IU,ERR=100)A%BK(I)%NROW,A%BK(I)%NCOL,A%BK(I)%LZERO,A%BK(I)%LSPARSE
        ELSE
          READ(IU,ERR=100)A%BK(I)%NROW,A%BK(I)%NCOL,A%BK(I)%LZERO; A%BK(I)%LSPARSE=.FALSE.
        ENDIF
        IF(A%BK(I)%LZERO)CYCLE
        IF(A%BK(I)%LSPARSE)THEN
          CALL DCSR_READ(A%BK(I)%ACSR,IU)
        ELSE
          ALLOCATE(A%BK(I)%A(A%BK(I)%NROW,A%BK(I)%NCOL))
          READ(IU)A%BK(I)%A
          GMEM_SIZE=GMEM_SIZE+REAL(SIZE(A%BK(I)%A),gq)*8*1.E-9_gq
        ENDIF
      ENDDO
      ! 'DBM_READ'
      RETURN
100   INFO=1
      RETURN
!
      END SUBROUTINE DBM_READ
!
!*********************************************************************
      SUBROUTINE ZBM_READ(A,IU,INFO)
      TYPE(ZBD_MATRIX)::A
      INTEGER IU,INFO
! LOCAL
      INTEGER I
!
      READ(IU)A%NBK,A%DIM,A%JSHIFT
      IF(INFO==0)THEN
        ALLOCATE(A%BK(A%NBK))
      ENDIF
      DO I=1,A%NBK
        IF(INFO==0)THEN
          READ(IU,ERR=100)A%BK(I)%NROW,A%BK(I)%NCOL,A%BK(I)%LZERO,A%BK(I)%LSPARSE
        ELSE
          READ(IU,ERR=100)A%BK(I)%NROW,A%BK(I)%NCOL,A%BK(I)%LZERO; A%BK(I)%LSPARSE=.FALSE.
        ENDIF
        IF(A%BK(I)%LZERO)CYCLE
        IF(A%BK(I)%LSPARSE)THEN
          CALL ZCSR_READ(A%BK(I)%ACSR,IU)
        ELSE
          ALLOCATE(A%BK(I)%A(A%BK(I)%NROW,A%BK(I)%NCOL))
          READ(IU)A%BK(I)%A
          GMEM_SIZE=GMEM_SIZE+REAL(SIZE(A%BK(I)%A),gq)*8*1.E-9_gq
        ENDIF
      ENDDO
      ! 'ZBM_READ'
      RETURN
100   INFO=1
      RETURN
!
      END SUBROUTINE ZBM_READ
!
!*********************************************************************
      SUBROUTINE ZBM_READ_TXT(A,IU)
      TYPE(ZBD_MATRIX)::A
      INTEGER IU
! LOCAL
      INTEGER I
!
      READ(IU,*)A%NBK,A%DIM,A%JSHIFT
      ALLOCATE(A%BK(A%NBK))
      DO I=1,A%NBK
        READ(IU,*)A%BK(I)%NROW,A%BK(I)%NCOL,A%BK(I)%LZERO; A%BK(I)%LSPARSE=.FALSE.
        IF(A%BK(I)%LZERO)CYCLE
        ALLOCATE(A%BK(I)%A(A%BK(I)%NROW,A%BK(I)%NCOL))
        READ(IU,*)A%BK(I)%A
        GMEM_SIZE=GMEM_SIZE+REAL(SIZE(A%BK(I)%A),gq)*8*1.E-9_gq
      ENDDO
      ! 'ZBM_READ_TXT'
      RETURN
!
      END SUBROUTINE ZBM_READ_TXT
!
!*********************************************************************
      SUBROUTINE DBM_READ_TXT(A,IU)
      TYPE(DBD_MATRIX)::A
      INTEGER IU
! LOCAL
      INTEGER I
!
      READ(IU,*)A%NBK,A%DIM,A%JSHIFT
      ALLOCATE(A%BK(A%NBK))
      DO I=1,A%NBK
        READ(IU,*)A%BK(I)%NROW,A%BK(I)%NCOL,A%BK(I)%LZERO; A%BK(I)%LSPARSE=.FALSE.
        IF(A%BK(I)%LZERO)CYCLE
        ALLOCATE(A%BK(I)%A(A%BK(I)%NROW,A%BK(I)%NCOL))
        READ(IU,*)A%BK(I)%A
        GMEM_SIZE=GMEM_SIZE+REAL(SIZE(A%BK(I)%A),gq)*8*1.E-9_gq
      ENDDO
      ! 'DBM_READ_TXT'
      RETURN
!
      END SUBROUTINE DBM_READ_TXT
!
!*********************************************************************
! Y = A.H * X
!*********************************************************************
      SUBROUTINE ZAHMUXR(M,N,X,Y,A,JA,IA)
      INTEGER M,N,IA(N+1),JA(*)
      COMPLEX(gq) X(N),Y(M),A(*)
! LOCAL
      INTEGER I,K
!
      Y=0
      DO I=1,N; DO K=IA(I),IA(I+1)-1
        Y(JA(K))=Y(JA(K))+X(I)*CONJG(A(K))
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE ZAHMUXR
!
!*****************************************************************************
      SUBROUTINE ZBM2_TRACE(A,B,TR)
      TYPE(ZBD_MATRIX)::A,B
      COMPLEX(gq) TR
! LOCAL
      INTEGER IBK
      COMPLEX(gq) TR_
!
      TR=0
      IF(A%JSHIFT/=0)RETURN
      IF(A%NBK/=B%NBK)THEN
        WRITE(0,'(" WARNING: ZBM2_TRACE FAILED SINCE A%NBK/=B%NBK!", I0, " VS ", I0)')A%NBK,B%NBK
        TR=-99._gq; RETURN
      ENDIF
      DO IBK=1,A%NBK
        IF(B%BK(IBK)%LSPARSE)THEN
          WRITE(0,'(" WARNING: ZBM2_TRACE FAILED SINCE B%BK(IBK)%LSPARSE==.TRUE.!")')
          TR=-99._gq; RETURN
        ENDIF
        IF(A%BK(IBK)%LSPARSE)THEN
          CALL TRACE_ZCSRDNS(A%BK(IBK)%ACSR,B%BK(IBK)%A,TR_)
        ELSE
          TR_=SUM(A%BK(IBK)%A*TRANSPOSE(B%BK(IBK)%A))
        ENDIF
        TR=TR+TR_
      ENDDO
      RETURN
!
      END SUBROUTINE ZBM2_TRACE
!
!*****************************************************************************
      SUBROUTINE DBM2_TRACE(A,B,TR)
      TYPE(DBD_MATRIX)::A,B
      REAL(gq) TR
! LOCAL
      INTEGER IBK
      REAL(gq) TR_
!
      TR=0
      IF(A%JSHIFT/=0)RETURN
      IF(A%NBK/=B%NBK)THEN
        WRITE(0,'(" WARNING: ZBM2_TRACE FAILED SINCE A%NBK/=B%NBK!", I0, " VS ", I0)')A%NBK,B%NBK
        TR=-99._gq; RETURN
      ENDIF
      DO IBK=1,A%NBK
        IF(B%BK(IBK)%LSPARSE)THEN
          WRITE(0,'(" WARNING: ZBM2_TRACE FAILED SINCE B%BK(IBK)%LSPARSE==.TRUE.!")')
          TR=-99._gq; RETURN
        ENDIF
        IF(A%BK(IBK)%LSPARSE)THEN
          CALL TRACE_DCSRDNS(A%BK(IBK)%ACSR,B%BK(IBK)%A,TR_)
        ELSE
          TR_=SUM(A%BK(IBK)%A*TRANSPOSE(B%BK(IBK)%A))
        ENDIF
        TR=TR+TR_
      ENDDO
      RETURN
!
      END SUBROUTINE DBM2_TRACE
!
!*****************************************************************************
      SUBROUTINE TRACE_ZCSRDNS(ACSR,B,TR)
      TYPE(ZCSR_MATRIX),INTENT(IN) :: ACSR
      COMPLEX(gq),INTENT(IN) :: B(ACSR%NROW,ACSR%NCOL)
      COMPLEX(gq),INTENT(OUT) :: TR
! LOCAL
      INTEGER I,J,J_
!
      TR=0
      DO I=1,ACSR%NROW; DO J_=ACSR%I(I),ACSR%I(I+1)-1
        J=ACSR%J(J_)
        TR=TR+ACSR%A(J_)*B(J,I)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE TRACE_ZCSRDNS
!
!*****************************************************************************
      SUBROUTINE TRACE_DCSRDNS(ACSR,B,TR)
      TYPE(DCSR_MATRIX),INTENT(IN) :: ACSR
      REAL(gq),INTENT(IN) :: B(ACSR%NROW,ACSR%NCOL)
      REAL(gq),INTENT(OUT) :: TR
! LOCAL
      INTEGER I,J,J_
!
      TR=0
      DO I=1,ACSR%NROW; DO J_=ACSR%I(I),ACSR%I(I+1)-1
        J=ACSR%J(J_)
        TR=TR+ACSR%A(J_)*B(J,I)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE TRACE_DCSRDNS
!
!
      END MODULE SPARSE

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
!********************************************************************************
      SUBROUTINE INI_GMPI_(IU)
      USE GUTZ
      USE MPI
      IMPLICIT NONE
      INTEGER IU
! LOCAL
      INTEGER MYRANK,NPROCS,MASTER,NVEC,IERR,MYRANK_,NPROCS_
      INTEGER,ALLOCATABLE::KVEC(:,:)
      LOGICAL LKPVEC
!
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
      OPEN(IU,FILE=TRIM(ADJUSTL(FILE_NAME(MYRANK,0,-3))),STATUS='OLD')
      READ(IU,*)MYRANK_,NPROCS_,MASTER,LKPVEC,NVEC
      IF(MYRANK_/=MYRANK)THEN
        WRITE(0,'(" ERROR: MYRANK VS MYRANK-IN",2I3)')MYRANK,MYRANK_; STOP
      ENDIF
      IF(NPROCS_/=NPROCS)THEN
        WRITE(0,'(" ERROR: NPROCS VS NPROCS-IN",2I3)')NPROCS,NPROCS_; STOP
      ENDIF
      IF(LKPVEC)THEN
        ALLOCATE(KVEC(NVEC,3))
        READ(IU,*)KVEC
      ENDIF
!
      IF(LKPVEC)THEN
        CALL INI_GMPI(MYRANK,NPROCS,MASTER,LKPVEC,NVEC,KVEC)
      ELSE
        CALL INI_GMPI(MYRANK,NPROCS,MASTER,LKPVEC,NVEC)
      ENDIF
      RETURN
!
      END SUBROUTINE INI_GMPI_
!
!********************************************************************************
      SUBROUTINE GUTZ1_INI_(IU,IO)
      USE GUTZ
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      INTEGER NAT,LUNIT
!
      OPEN(IU,FILE='GUTZ1.INP',STATUS='OLD')
      READ(IU,*)NAT,LUNIT
      CLOSE(IU)
!
      CALL GUTZ1_INI(NAT,IU,IO,LUNIT)
      RETURN
!
      END SUBROUTINE GUTZ1_INI_
!
!********************************************************************************
      SUBROUTINE GUTZ2_SET_NISO_()
      USE GUTZ
      IMPLICIT NONE
      INTEGER ISO,ISPIN_IN,NBMAX,NKPT
!
      OPEN(GL%IU,FILE='GUTZ2.INP',STATUS='OLD')
      READ(GL%IU,*)ISO,ISPIN_IN,NBMAX,NKPT
      CLOSE(GL%IU)
!
      CALL GUTZ2_SET_NISO(ISO,ISPIN_IN,NBMAX,NKPT)
      RETURN
!
      END SUBROUTINE GUTZ2_SET_NISO_
!
!********************************************************************************
      SUBROUTINE GUTZ3_LOCATE_SYM_IE_()
      USE GUTZ
      IMPLICIT NONE
      INTEGER  IORD,I
      INTEGER,ALLOCATABLE :: IZ(:,:,:)
      REAL(gq),ALLOCATABLE :: TAU(:,:)
!
      OPEN(GL%IU,FILE='GUTZ3.INP',STATUS='OLD')
      READ(GL%IU,*)IORD
      IF(IORD>0)THEN
        ALLOCATE(IZ(3,3,IORD),TAU(3,IORD))
        READ(GL%IU,*)(IZ(:,:,I),TAU(:,I),I=1,IORD)
      ENDIF
      CLOSE(GL%IU)
!
      IF(IORD>0)THEN
        CALL GUTZ3_LOCATE_SYM_IE(IORD,IZ,TAU)
      ELSE
        CALL GUTZ3_LOCATE_SYM_IE(IORD)
      ENDIF
      RETURN
!
      END SUBROUTINE GUTZ3_LOCATE_SYM_IE_
!
!********************************************************************************
      SUBROUTINE GUTZ4_SET_C2N_UH_()
      USE GUTZ
      IMPLICIT NONE
      INTEGER NASOMAX,NIONS
      COMPLEX(gq),ALLOCATABLE::C2N(:,:,:)
!
      OPEN(GL%IU,FILE='GUTZ4.INP',STATUS='OLD')
      READ(GL%IU,*)NASOMAX,NIONS
      ALLOCATE(C2N(NASOMAX,NASOMAX,NIONS))
      READ(GL%IU,*)C2N
      CLOSE(GL%IU)
!
      CALL GUTZ4_SET_C2N_UH(C2N,NASOMAX,NIONS)
      RETURN
!
      END SUBROUTINE GUTZ4_SET_C2N_UH_
!
!********************************************************************************
      SUBROUTINE GUTZ5_INI_CORR_UVEK_()
      USE GUTZ
      IMPLICIT NONE
      INTEGER NKP,NKPL,NKPP,ISMEAR
      REAL(gq) NELET,DELTA
      INTEGER NI,NBASE,NASO,I
      INTEGER ,ALLOCATABLE::NE(:,:)
      REAL(gq),ALLOCATABLE::WT(:),KX(:),KY(:),KZ(:)
      CHARACTER*10,ALLOCATABLE::KNAME(:)
!
      KPT%FILE_NAME=''; NKPP=0
      OPEN(GL%IU,FILE='GUTZ5.INP',STATUS='OLD')
      READ(GL%IU,*)NKP
      ALLOCATE(WT(MAX(1,NKP)),NE(3,MAX(1,NKP)))
      IF(NKP<=0)THEN ! Impurity model
        NKP=1; WT=0; ISMEAR=-1; DELTA=0.01_gq
        NE=0; NKPP=2
        READ(GL%IU,*)NELET
        IF(GL%LMODEL==0)THEN
          STOP ' ERROR IN GUTZ5_INI_CORR_UVEK_: NKP<0 WHILE GL%LMODEL=0!'
        ENDIF
      ELSE
        READ(GL%IU,*)WT
        READ(GL%IU,*)ISMEAR,DELTA
        READ(GL%IU,*)NELET
        READ(GL%IU,*)NE
        READ(GL%IU,*,ERR=100,END=100)KPT%FILE_NAME
        READ(GL%IU,*)NKPP
        IF(NKPP==NKP)THEN
          ALLOCATE(KX(NKP),KY(NKP),KZ(NKP),KNAME(NKP))
          READ(GL%IU,*)(KX(I),KY(I),KZ(I),KNAME(I),I=1,NKP)
        ENDIF
100     CONTINUE
      ENDIF
      CLOSE(GL%IU)
!
      KPT%ISMEAR=ISMEAR; KPT%DELTA=DELTA
      CALL SET_KPT_ICOR(GL%IO)
      IF(GP%LKPVEC)THEN
        NKPL=SUM(GP%KVEC(:,2))
      ELSE
        NKPL=CEILING(NKP/REAL(GP%NPROCS,gq))
      ENDIF
      IF(GL%LMODEL/=0)THEN
        NKPL=0
      ENDIF
      IF(NKPP==NKP)THEN
        CALL GUTZ5_INI_CORR_UVEK(NKP,NKPL,WT,NELET,NE,KX,KY,KZ,KNAME)
      ELSE
        CALL GUTZ5_INI_CORR_UVEK(NKP,NKPL,WT,NELET,NE)
      ENDIF
      IF(GL%LMODEL==0)THEN
        CALL READ_BNDU(GL%IU,GL%LBNDU)
        IF(GL%LUNIT==0)THEN
          CALL ROTATE_BNDU('N',0) ! VASP CONVENTION, POST ROTATE
        ELSE
          CALL ROTATE_BNDU('N',1) ! Apply additional rotation.
        ENDIF
      ENDIF
      RETURN
!
      END SUBROUTINE GUTZ5_INI_CORR_UVEK_
!
!

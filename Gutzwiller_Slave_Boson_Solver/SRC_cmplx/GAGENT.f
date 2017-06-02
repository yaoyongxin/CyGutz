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
!********************************************************************************
      SUBROUTINE INI_GMPI_(IU,IO)
      USE GUTZ
      USE MPI
      IMPLICIT NONE
      INTEGER IU,IO
! LOCAL
      INTEGER MYRANK,NPROCS,MASTER,NVEC,IERR,I,MYRANK_,NPROCS_
      INTEGER,ALLOCATABLE::KVEC(:,:)
      LOGICAL LKPVEC,LOMP
      CHARACTER*77 F_NAME
!
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
      INQUIRE(FILE='GOMP.INP',EXIST=LOMP)
      IF(LOMP)THEN
        F_NAME='GOMP.INP'
      ELSE
        F_NAME=FILE_NAME(MYRANK,0,0,-3)
      ENDIF
      ! GMPI_0.INP for myrank = 0.
      OPEN(IU,FILE=TRIM(ADJUSTL(F_NAME)),STATUS='OLD')
      ! my rank of the current processor, total number of processors. 
      ! the master processor
      ! LKPVEC = true and NVEC: Wien2k k-ponits parallelization
      READ(IU,*)MYRANK_,NPROCS_,MASTER,LKPVEC,NVEC
      IF(LOMP)THEN
        CALL OMP_SET_NUM_THREADS(NPROCS_)
        NPROCS=NPROCS_
        WRITE(IO,'(" OMP_NUM_THREADS=",I2)')NPROCS_
      ELSE
        CALL OMP_SET_NUM_THREADS(1)
      ENDIF
      IF(MYRANK_/=MYRANK)THEN
        WRITE(0,'(" ERROR: MYRANK VS MYRANK-IN",2I3)')MYRANK,MYRANK_; STOP
      ENDIF
      IF(NPROCS_/=NPROCS)THEN
        WRITE(0,'(" ERROR: NPROCS VS NPROCS-IN",2I3)')NPROCS,NPROCS_; STOP
      ENDIF
      IF(LKPVEC)THEN ! Wien2k k-ponits parallelization
        ALLOCATE(KVEC(NVEC,3))
        READ(IU,*) (KVEC(I,:),I=1,NVEC)
      ENDIF
!
      IF(LKPVEC)THEN
        CALL INI_GMPI(MYRANK,NPROCS,MASTER,LOMP,LKPVEC,NVEC,KVEC)
      ELSE
        CALL INI_GMPI(MYRANK,NPROCS,MASTER,LOMP,LKPVEC,NVEC)
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
      ! NAT: number of inequivalent atoms in case.struct.
      ! LUNIT = 1: Wien2k convention, Rydberg/Bohr.
      ! LUNIT = 0: VASP convention, eV/A.
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
      ! ISO = 1: No spin-orbit; 2: spin-orbit
      ! ISPIN_IN (for the input band dispersion) = 1: no spin-splitting; 2: spin splitting
      ! NBMAX: Maximal number of bands (for one k-point) involved over all the k-points.
      ! NKPT: number of k-points
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
      ! Wien2k symmetry operations.
      ! IORD: number of symmetry operations
      READ(GL%IU,*)IORD
      IF(IORD>0)THEN
        ! Rotations + translations
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
      USE GUTZ; USE GCONSTANT
      IMPLICIT NONE
      INTEGER NASOMAX,NIONS,I,J
      COMPLEX(gq),ALLOCATABLE::C2N(:,:,:)
      REAL(gq),ALLOCATABLE::RA(:,:),IA(:,:)
!
      OPEN(GL%IU,FILE='GUTZ4.INP',STATUS='OLD')
      ! Maximal number of spin-orbitals (for each atom) over all the correlated atoms.
      READ(GL%IU,*)NASOMAX,NIONS
      ! C2N = < complex spherical harmonics | rotated basis in case.indmfl >
      ALLOCATE(C2N(NASOMAX,NASOMAX,NIONS)); C2N = 0
      ALLOCATE(RA(NASOMAX,NASOMAX),IA(NASOMAX,NASOMAX))
      DO I=1,NIONS
        READ(GL%IU,*, END = 100, ERR = 100)RA
        READ(GL%IU,*, END = 100, ERR = 100)IA
        C2N(:,:,I)=RA+ZI*IA
      ENDDO
      GOTO 101
100   CONTINUE
      DO I=1,NIONS; DO J=1,NASOMAX
        C2N(J, J, I) = 1._gq
      ENDDO; ENDDO
101   CONTINUE
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
      ! NKP: number of k-points
      READ(GL%IU,*)NKP
      ! WT: k-point weight
      !  NE(1:3): total number of bands, start index and end index for correlated bands.
      ALLOCATE(WT(MAX(1,NKP)),NE(3,MAX(1,NKP)))
      IF(NKP<=0)THEN ! Impurity model
        NKP=1; WT=0; ISMEAR=-1; DELTA=0.01_gq
        NE=0; NKPP=2
        ! NELET: number of electrons.
        READ(GL%IU,*)NELET
        IF(GL%LMODEL==0)THEN
          STOP ' ERROR IN GUTZ5_INI_CORR_UVEK_: NKP<0 WHILE GL%LMODEL=0!'
        ENDIF
      ELSE
        READ(GL%IU,*)WT
        ! ISMEAR = -5: tetrahedron BZ integral; -1: fermi smearing; 0: Gaussian smearing
        ! DELTA: broadening factor
        READ(GL%IU,*)ISMEAR,DELTA
        READ(GL%IU,*)NELET
        READ(GL%IU,*)NE
        ! case.kgen file name for wien2k.
        READ(GL%IU,*,ERR=100,END=100)KPT%FILE_NAME
        READ(GL%IU,*)NKPP
        IF(NKPP==NKP)THEN
          ALLOCATE(KX(NKP),KY(NKP),KZ(NKP),KNAME(NKP))
          DO I=1,NKP
            READ(GL%IU,'(3F16.12,A10)')KX(I),KY(I),KZ(I),KNAME(I)
          ENDDO
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
        CALL PROC_BNDU(GL%LBNDU)
        CALL ROTATE_BNDU('N')
      ENDIF
      RETURN
!
      END SUBROUTINE GUTZ5_INI_CORR_UVEK_
!

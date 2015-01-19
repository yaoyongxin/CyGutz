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
      MODULE WAREHOUSE
      USE gprec; USE GMPI; USE GUTIL
      IMPLICIT NONE
! DEFINE TYPE
      TYPE WARE_HOUSE
        INTEGER DIMX1,DIMX2,DIMX3,DIMXT,NIONS,NTYP,DIMXG,DIMXH,DIMXHO
        INTEGER NA2MAX,NASOMAX,NA2TOT,NASOTOT
        COMPLEX(gq),POINTER :: C2N(:,:,:)  ! Complex Harmonics to "Natural" basis
        COMPLEX(gq),POINTER :: R2N(:,:,:),B2N(:,:,:)  ! Original basis (can be real Harmonics) to "Natural" basis
        COMPLEX(gq),POINTER :: N2N(:,:,:) ! Additional rotations
        COMPLEX(gq),ALLOCATABLE :: X(:),XNC(:),XN0(:) ! Fixed point
        COMPLEX(gq),POINTER :: R(:,:,:),R0(:,:,:),D(:,:,:),D0(:,:,:),LA1(:,:,:),LA2(:,:,:),ETA(:,:,:),Z(:,:,:)
        COMPLEX(gq),POINTER :: VDC2(:,:,:),NPHY_FIX(:,:,:)
        REAL(gq),POINTER :: TFC(:,:),TFF(:,:,:) ! f-c and f-f compoments
        REAL(gq) :: TCC,EF
        COMPLEX(gq),POINTER :: NKS(:,:,:),NC_VAR(:,:,:),NC_PHY(:,:,:),EL0(:,:,:)
        COMPLEX(gq),POINTER :: ISIMIX(:,:,:) ! 1/SQRT(NKS(1-NKS))
        INTEGER,ALLOCATABLE :: LXREAL(:)
      END TYPE WARE_HOUSE
! VARIABLE
      TYPE (WARE_HOUSE),SAVE :: WH
!
! SUBROUTINE
      CONTAINS
!
!*************************************************************************************
      SUBROUTINE ALLOC_WAREHOUSE(NA2MAX,NIONS,NA2TOT)
      INTEGER NA2MAX,NIONS,NA2TOT
!
      WH%NIONS=NIONS; WH%NA2MAX=NA2MAX; WH%NA2TOT=NA2TOT
      ALLOCATE(WH%R     (NA2MAX,NA2MAX,NIONS), &
              &WH%R0    (NA2MAX,NA2MAX,NIONS), & 
              &WH%C2N   (NA2MAX,NA2MAX,NIONS), &
              &WH%R2N   (NA2MAX,NA2MAX,NIONS), &
              &WH%N2N   (NA2MAX,NA2MAX,NIONS), &
              &WH%Z     (NA2MAX,NA2MAX,NIONS), &
              &WH%D     (NA2MAX,NA2MAX,NIONS), &
              &WH%D0    (NA2MAX,NA2MAX,NIONS), &
              &WH%LA1   (NA2MAX,NA2MAX,NIONS), &
              &WH%LA2   (NA2MAX,NA2MAX,NIONS), &
              &WH%ETA   (NA2MAX,NA2MAX,NIONS), &
              &WH%NKS   (NA2MAX,NA2MAX,NIONS), &
              &WH%ISIMIX(NA2MAX,NA2MAX,NIONS), &
              &WH%NC_VAR(NA2MAX,NA2MAX,NIONS), &
              &WH%NC_PHY(NA2MAX,NA2MAX,NIONS), &
              &WH%VDC2  (NA2MAX,NA2MAX,NIONS), &
              &WH%NPHY_FIX(NA2MAX,NA2MAX,NIONS), &
              &WH%EL0   (NA2MAX,NA2MAX,NIONS)  )
      WH%LA1=0; WH%LA2=0; WH%ETA=0; WH%EL0=0; WH%VDC2=0
      WH%NPHY_FIX(1,1,1)=2._gq
      RETURN
!
      END SUBROUTINE ALLOC_WAREHOUSE
!
!****************************************************************************
      SUBROUTINE WRT_WH_X(IU)
      INTEGER IU,I
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      OPEN(IU,FILE='GLX.OUT',STATUS='REPLACE')
      WRITE(IU,*)WH%X(1:WH%DIMX1); WRITE(IU,*)
      WRITE(IU,*)WH%X(1+WH%DIMX1:WH%DIMX2); WRITE(IU,*)
      WRITE(IU,*)WH%X(1+WH%DIMX2:WH%DIMX3); WRITE(IU,*)
      WRITE(IU,*)WH%X(1+WH%DIMX3:WH%DIMXT); WRITE(IU,*)
      WRITE(IU,*)WH%EF
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE WRT_WH_X
!
!****************************************************************************
      SUBROUTINE WRT_WHEL0(IU)
      INTEGER IU
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      OPEN(IU,FILE='WHEL0.OUT',STATUS='REPLACE')
      WRITE(IU,'("WH%EL0:" )'); WRITE(IU,*)WH%EL0
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE WRT_WHEL0
!
!****************************************************************************
      SUBROUTINE READ_WHEL0(IU)
      INTEGER IU
!
      OPEN(IU,FILE='WHEL0.INP',STATUS='OLD')
      READ(IU,*)WH%EL0
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE READ_WHEL0
!
!****************************************************************************
      SUBROUTINE WRT_WH_RLNEF(IU)
      INTEGER IU
! LOCAL
      INTEGER NI
      COMPLEX(gq) ZBUF(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
      CHARACTER(40) FMT
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      WRITE(FMT,'(A,I3,A)')"(",WH%NA2MAX,'("(",F18.12,",",F18.12,")"))'
      OPEN(IU,FILE='WH_RLNEF.OUT',STATUS='REPLACE')
      WRITE(IU,*) WH%NA2MAX,WH%NIONS
      ZBUF=WH%R;   DO NI=1,WH%NIONS; WRITE(IU,FMT)ZBUF(:,:,NI); WRITE(IU,*); ENDDO
      WRITE(IU,*)
      ZBUF=WH%LA1; DO NI=1,WH%NIONS; WRITE(IU,FMT)ZBUF(:,:,NI); WRITE(IU,*); ENDDO
      WRITE(IU,*)
      ZBUF=WH%NKS; DO NI=1,WH%NIONS; WRITE(IU,FMT)ZBUF(:,:,NI); WRITE(IU,*); ENDDO
      WRITE(IU,*)
      ZBUF=WH%ETA; DO NI=1,WH%NIONS; WRITE(IU,FMT)ZBUF(:,:,NI); WRITE(IU,*); ENDDO
      WRITE(IU,*)
      WRITE(IU,*)WH%EF
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE WRT_WH_RLNEF
!
!****************************************************************************
      SUBROUTINE READ_WH_RLNEF(IU)
      INTEGER IU
! LOCAL
      INTEGER NA2MAX,NIONS
      COMPLEX(gq),ALLOCATABLE :: ZBUF(:,:,:)
!
      OPEN(IU,FILE='WH_RLNEF.INP',STATUS='OLD')
      READ(IU,*)NA2MAX,NIONS
      ALLOCATE(ZBUF(NA2MAX,NA2MAX,NIONS))
      READ(IU,*)ZBUF; WH%R  (1:NA2MAX,1:NA2MAX,1:NIONS)=ZBUF
      READ(IU,*)ZBUF; WH%LA1(1:NA2MAX,1:NA2MAX,1:NIONS)=ZBUF
      READ(IU,*,ERR=100,END=100)ZBUF; WH%NKS(1:NA2MAX,1:NA2MAX,1:NIONS)=ZBUF
      READ(IU,*,ERR=100,END=100)ZBUF; WH%ETA(1:NA2MAX,1:NA2MAX,1:NIONS)=ZBUF
      READ(IU,*,ERR=100,END=100)WH%EF
100   CONTINUE
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE READ_WH_RLNEF
!
!****************************************************************************
! Same format as GUTZ4.INP
!****************************************************************************
      SUBROUTINE SET_WH_N2N(IU,IO)
      INTEGER IU,IO
! LOCAL
      INTEGER I,NA2MAX,NIONS
!
      OPEN(IU,FILE='WH_N2N.INP',STATUS='OLD',ERR=100)
      READ(IU,*)NA2MAX,NIONS
      IF(NA2MAX/=WH%NA2MAX.OR.NIONS/=WH%NIONS)THEN
        STOP ' ERROR IN SET_WH_N2N: DIMENSION NOT MATCH!'
      ENDIF
      READ(IU,*)WH%N2N
      CLOSE(IU)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(IO,'(" WH_N2N.INP READ IN.")')
      RETURN
100   CONTINUE
      WH%N2N=0
      DO I=1,WH%NA2MAX; WH%N2N(I,I,:)=1._gq; ENDDO
      RETURN
!
      END SUBROUTINE SET_WH_N2N
!
!
      END MODULE WAREHOUSE

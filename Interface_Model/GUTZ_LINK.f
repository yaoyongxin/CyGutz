!****************************************************************************
      SUBROUTINE GMPI_WRT()
      IMPLICIT NONE
      INTEGER,PARAMETER :: IU=90

      OPEN(IU,FILE='GMPI_0.INP',STATUS='REPLACE')
      WRITE(IU,*)0,1,0,.FALSE.,1 
      WRITE(IU,*)' MYRANK,NPROCS,MASTER,VECTOR_PARA,NVECTOR,VECTORS'
      CLOSE(IU)
      RETURN

      END SUBROUTINE GMPI_WRT

!****************************************************************************
      SUBROUTINE GUTZ1_WRT(NAT)
      IMPLICIT NONE
      INTEGER NAT
! LOCAL
      INTEGER,PARAMETER :: IU=90

      OPEN(IU,FILE='GUTZ1.INP',STATUS='REPLACE')
      WRITE(IU,*)NAT,0
      WRITE(IU,'(" NAT,LUNIT!")')
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ1_WRT

!****************************************************************************
      SUBROUTINE GUTZ2_WRT(ISO,ISO2,NUME,NKPT)
      IMPLICIT NONE
      INTEGER ISO,ISO2,NUME,NKPT
! LOCAL
      INTEGER,PARAMETER :: IU=90

      OPEN(IU,FILE='GUTZ2.INP',STATUS='REPLACE')
      WRITE(IU,*)ISO,ISO2,NUME,NKPT,1.D0
      WRITE(IU,'(" ISO,ISO2,NUME,NKPT")')
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ2_WRT

!****************************************************************************
      SUBROUTINE GUTZ3_WRT()
      IMPLICIT NONE
      INTEGER,PARAMETER :: IORD=1
      INTEGER IZ(3,3,IORD)
      REAL(8) TAU(3,IORD)
! LOCAL
      INTEGER I
      INTEGER,PARAMETER :: IU=90

      IZ=0; TAU=0
      DO I=1,3; IZ(I,I,IORD)=1; ENDDO
      OPEN(IU,FILE='GUTZ3.INP',STATUS='REPLACE')
      WRITE(IU,*)IORD
      WRITE(IU,'(3(3I6,/),3F6.2,/)')(IZ(:,:,I),TAU(:,I),I=1,IORD)
      WRITE(IU,*)' ! (IZ(:,:,I),TAU(:,I),I=1,IORD)'
      RETURN

      END SUBROUTINE GUTZ3_WRT

!*****************************************************************
! Spherical Harmonics->Projector
!*****************************************************************
      SUBROUTINE GUTZ4_WRT(C2N,MAXDIM2,NORBITALS)
      IMPLICIT NONE
      INTEGER MAXDIM2,NORBITALS
      COMPLEX(8) C2N(MAXDIM2,MAXDIM2,NORBITALS)
! LOCAL
      INTEGER,PARAMETER :: IU=90

      OPEN(IU,FILE='GUTZ4.INP',STATUS='REPLACE')
      WRITE(IU,*)MAXDIM2,NORBITALS
      WRITE(IU,*)C2N
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ4_WRT

!****************************************************************************
      SUBROUTINE GUTZ5_WRT(NKP,WT,NELET,NE,EFMOD,DELTA)
      IMPLICIT NONE
      INTEGER NKP
      REAL(8) WT(NKP),NELET,DELTA
      INTEGER :: NE(3,NKP)
      CHARACTER*5::EFMOD
! LOCAL
      INTEGER I,ISMEAR
      INTEGER,PARAMETER :: IU=90

      OPEN(IU,FILE='GUTZ5.INP',STATUS='REPLACE')
      WRITE(IU,*)NKP
      WRITE(IU,*)WT
      IF(EFMOD=='GAUSS')THEN
        ISMEAR=0
      ELSEIF(EFMOD=='TETRA')THEN
        ISMEAR=-5
      ELSE
        WRITE(0,'(" EFMOD=",A6)')EFMOD; STOP ' ERROR: UNSUPPORTED EFMOD!'
      ENDIF
      WRITE(IU,*)ISMEAR,DELTA
      WRITE(IU,*)NELET
      WRITE(IU,*)NE
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ5_WRT


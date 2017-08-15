!*************************************************!
! COPYRIGHT BY YONGXIN YAO AND NICOLA LANATA 2012 !
!*************************************************!
      MODULE GMPI
      USE gprec
      IMPLICIT NONE
#ifdef MPI
  include 'mpif.h'
#endif
      TYPE G_MPI
        INTEGER NPROCS
        INTEGER MYRANK
        CHARACTER*10 CPUID ! MYRANK in string representation
        INTEGER MASTER
        INTEGER NIONL
! K-parallel of Wien2K
        LOGICAL LKPVEC
        INTEGER NVEC
        INTEGER,POINTER::KVEC(:,:) ! NKP, NKSTART0
      END TYPE G_MPI

      TYPE(G_MPI) GP

      CONTAINS

!**************************************************
      SUBROUTINE INI_GMPI(MYRANK,NPROCS,MASTER,LKPVEC,NVEC,KVEC)
      INTEGER MYRANK,NPROCS,MASTER,NVEC
      INTEGER,OPTIONAL::KVEC(NVEC,3)
      LOGICAL LKPVEC

      GP%NPROCS=NPROCS; GP%MYRANK=MYRANK; GP%MASTER=MASTER
      WRITE(GP%CPUID,'(I10)')GP%MYRANK
      GP%LKPVEC=LKPVEC; GP%NVEC=NVEC
      IF(PRESENT(KVEC))THEN
        ALLOCATE(GP%KVEC(NVEC,3)); GP%KVEC=KVEC
      ENDIF
      RETURN

      END SUBROUTINE INI_GMPI


#ifdef MPI
!**************************************************
      SUBROUTINE GMPI_Barrier()
      INTEGER IERR

      CALL MPI_Barrier(MPI_COMM_WORLD,IERR)

      END SUBROUTINE GMPI_Barrier

!**************************************************
      SUBROUTINE IMAX1_MASTER_MPI(N)
      INTEGER N
! LOCAL
      INTEGER I(1)

      I(1)=N; CALL IMAX_MASTER_MPI(I,1); N=I(1)
      RETURN

      END SUBROUTINE IMAX1_MASTER_MPI

!**************************************************
      SUBROUTINE IMAX_MASTER_MPI(I,N)
      INTEGER N,I(N)
! LOCAL
      INTEGER IERR,MAXI(N)

      CALL MPI_REDUCE(I,MAXI,N,MPI_INTEGER,MPI_MAX,GP%MASTER,MPI_COMM_WORLD,IERR)
      IF(GP%MYRANK.EQ.GP%MASTER)I=MAXI
      RETURN

      END SUBROUTINE IMAX_MASTER_MPI

!**************************************************
      SUBROUTINE DMAX1_MASTER_MPI(A)
      REAL(gq) A
! LOCAL
      REAL(gq) B(1)

      B(1)=A; CALL DMAX_MASTER_MPI(B,1); A=B(1)
      RETURN

      END SUBROUTINE DMAX1_MASTER_MPI

!**************************************************
      SUBROUTINE DMAX_MASTER_MPI(A,N)
      INTEGER N
      REAL(gq) A(N)
! LOCAL
      INTEGER IERR
      REAL(gq) MAXA(N)

      CALL MPI_REDUCE(A,MAXA,N,MPI_DOUBLE_PRECISION,MPI_MAX,GP%MASTER,MPI_COMM_WORLD,IERR)
      IF(GP%MYRANK.EQ.GP%MASTER)A=MAXA
      RETURN

      END SUBROUTINE DMAX_MASTER_MPI

!**************************************************
      SUBROUTINE IMIN1_MASTER_MPI(N)
      INTEGER N
! LOCAL
      INTEGER M(1)

      M(1)=N; CALL IMIN_MASTER_MPI(M,1); N=M(1)
      RETURN

      END SUBROUTINE IMIN1_MASTER_MPI

!**************************************************
      SUBROUTINE IMIN_MASTER_MPI(I,N)
      INTEGER N,I(N)
! LOCAL
      INTEGER IERR,MINI(N)

      CALL MPI_REDUCE(I,MINI,N,MPI_INTEGER,MPI_MIN,GP%MASTER,MPI_COMM_WORLD,IERR)
      IF(GP%MYRANK.EQ.GP%MASTER)I=MINI
      RETURN

      END SUBROUTINE IMIN_MASTER_MPI

!**************************************************
      SUBROUTINE ISUM1_MASTER_MPI(I)
      INTEGER I
! LOCAL
      INTEGER J(1)

      J(1)=I; CALL ISUM_MASTER_MPI(J,1); I=J(1)
      RETURN

      END SUBROUTINE ISUM1_MASTER_MPI

!**************************************************
      SUBROUTINE ISUM_MASTER_MPI(I,N)
      INTEGER N,I(N)
! LOCAL
      INTEGER IERR
      INTEGER,ALLOCATABLE::J(:)

      ALLOCATE(J(N)); J=0
      CALL MPI_REDUCE(I,J,N,MPI_INTEGER,MPI_SUM,GP%MASTER,MPI_COMM_WORLD,IERR)
      IF(GP%MYRANK.EQ.GP%MASTER)I=J
      DEALLOCATE(J)
      RETURN

      END SUBROUTINE ISUM_MASTER_MPI

!**************************************************
      SUBROUTINE DSUM1_MASTER_MPI(A)
      REAL(gq) A
! LOCAL
      REAL(gq) B(1)

      B(1)=A; CALL DSUM_MASTER_MPI(B,1); A=B(1)
      RETURN

      END SUBROUTINE DSUM1_MASTER_MPI

!**************************************************
      SUBROUTINE DSUM_MASTER_MPI(A,N)
      INTEGER N
      REAL(gq) A(N)
! LOCAL
      INTEGER IERR
      REAL(gq),ALLOCATABLE::B(:)

      ALLOCATE(B(N)); B=0
      CALL MPI_REDUCE(A,B,N,MPI_DOUBLE_PRECISION,MPI_SUM,GP%MASTER,MPI_COMM_WORLD,IERR)
      IF(GP%MYRANK.EQ.GP%MASTER)A=B
      DEALLOCATE(B)
      RETURN

      END SUBROUTINE DSUM_MASTER_MPI

!**************************************************
      SUBROUTINE CSUM1_MASTER_MPI(A)
      COMPLEX(gq) A
! LOCAL
      COMPLEX(gq) B(1)

      B(1)=A; CALL CSUM_MASTER_MPI(B,1); A=B(1)
      RETURN

      END SUBROUTINE CSUM1_MASTER_MPI

!**************************************************
      SUBROUTINE CSUM_MASTER_MPI(A,N)
      INTEGER N
      COMPLEX(gq) A(N)
! LOCAL
      INTEGER IERR
      COMPLEX(gq),ALLOCATABLE::B(:)

      ALLOCATE(B(N)); B=0
      CALL MPI_REDUCE(A,B,N,MPI_DOUBLE_COMPLEX,MPI_SUM,GP%MASTER,MPI_COMM_WORLD,IERR)
      IF(GP%MYRANK.EQ.GP%MASTER)A=B
      DEALLOCATE(B)
      RETURN

      END SUBROUTINE CSUM_MASTER_MPI

!**************************************************
      SUBROUTINE IMAX1_ALL_MPI(I)
      INTEGER I
! LOCAL
      INTEGER J(1)

      J(1)=I; CALL IMAX_ALL_MPI(J,1); I=J(1)
      RETURN

      END SUBROUTINE IMAX1_ALL_MPI

!**************************************************
      SUBROUTINE IMAX_ALL_MPI(I,N)
      INTEGER N,I(N)
! LOCAL
      INTEGER IERR,MAXI(N)

      CALL MPI_ALLREDUCE(I,MAXI,N,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
      I=MAXI
      RETURN

      END SUBROUTINE IMAX_ALL_MPI

!**************************************************
      SUBROUTINE IMIN1_ALL_MPI(I)
      INTEGER I
! LOCAL
      INTEGER J(1)

      J(1)=I; CALL IMIN_ALL_MPI(J,1); I=J(1)
      RETURN

      END SUBROUTINE IMIN1_ALL_MPI

!**************************************************
      SUBROUTINE IMIN_ALL_MPI(I,N)
      INTEGER N,I(N)
! LOCAL
      INTEGER IERR,MINI(N)

      CALL MPI_ALLREDUCE(I,MINI,N,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,IERR)
      I=MINI
      RETURN

      END SUBROUTINE IMIN_ALL_MPI

!**************************************************
      SUBROUTINE DMAX1_ALL_MPI(A)
      REAL(gq) A
! LOCAL
      REAL(gq) B(1)

      B(1)=A; CALL DMAX_ALL_MPI(B,1); A=B(1)
      RETURN

      END SUBROUTINE DMAX1_ALL_MPI

!**************************************************
      SUBROUTINE DMAX_ALL_MPI(A,N)
      INTEGER N
      REAL(gq) A(N)
! LOCAL
      INTEGER IERR
      REAL(gq) MAXA(N)

      CALL MPI_ALLREDUCE(A,MAXA,N,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
      A=MAXA
      RETURN

      END SUBROUTINE DMAX_ALL_MPI

!**************************************************
      SUBROUTINE ISUM1_ALL_MPI(I)
      INTEGER I
! LOCAL
      INTEGER J(1)

      J(1)=I; CALL ISUM_ALL_MPI(J,1); I=J(1)
      RETURN

      END SUBROUTINE ISUM1_ALL_MPI

!**************************************************
      SUBROUTINE ISUM_ALL_MPI(I,N)
      INTEGER N,I(N)
! LOCAL
      INTEGER IERR
      INTEGER,ALLOCATABLE::J(:)

      ALLOCATE(J(N)); J=0
      CALL MPI_ALLREDUCE(I,J,N,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      I=J; DEALLOCATE(J)
      RETURN

      END SUBROUTINE ISUM_ALL_MPI

!**************************************************
      SUBROUTINE DSUM1_ALL_MPI(A)
      REAL(gq) A
! LOCAL
      REAL(gq) B(1)

      B(1)=A; CALL DSUM_ALL_MPI(B,1); A=B(1)
      RETURN

      END SUBROUTINE DSUM1_ALL_MPI

!**************************************************
      SUBROUTINE DSUM_ALL_MPI(A,N)
      INTEGER N
      REAL(gq) A(N)
! LOCAL
      INTEGER IERR
      REAL(gq),ALLOCATABLE::B(:)

      ALLOCATE(B(N)); B=0
      CALL MPI_ALLREDUCE(A,B,N,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
      A=B; DEALLOCATE(B)
      RETURN

      END SUBROUTINE DSUM_ALL_MPI

!**************************************************
      SUBROUTINE CSUM1_ALL_MPI(A)
      COMPLEX(gq) A
! LOCAL
      COMPLEX(gq) B(1)

      B(1)=A; CALL CSUM_ALL_MPI(B,1); A=B(1)
      RETURN

      END SUBROUTINE CSUM1_ALL_MPI

!**************************************************
      SUBROUTINE CSUM_ALL_MPI(A,N)
      INTEGER N
      COMPLEX(gq) A(N)
! LOCAL
      INTEGER IERR
      COMPLEX(gq),ALLOCATABLE::B(:)

      ALLOCATE(B(N)); B=0
      CALL MPI_ALLREDUCE(A,B,N,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
      A=B; DEALLOCATE(B)
      RETURN

      END SUBROUTINE CSUM_ALL_MPI

#endif

      END MODULE GMPI

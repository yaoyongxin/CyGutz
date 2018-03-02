!**************************************************!
!******* § Copyright by Kristjan Haule, 2002 ******!
!**************************************************!
MODULE com_mpi
!!! This module should contain everything connected with paralelezation 
!!! of the code with MPI

#ifdef MPI
  include 'mpif.h'
#endif

  INTEGER :: myrank    ! processor ID
  INTEGER :: nprocs    ! # of all processors awailable
  INTEGER :: ierr      ! returned error code
  INTEGER :: master    ! # of master processor
  CHARACTER*3 :: cpuID ! number of cpu in string representation
  INTEGER, PARAMETER :: clen = 100      ! Length of characters in command-line
  CHARACTER*100, ALLOCATABLE :: argv(:) ! Command-line arguments
  INTEGER      :: nargs                 ! Number of command line arguments
  LOGICAL      :: vector_para           ! is this parallel Wien2K run?
  INTEGER      :: vectors(20,3)        ! which vector files should be read
  INTEGER      :: nvector
  CHARACTER*200:: VECFN(4)
  CHARACTER*200:: fvectors(20,4)
  INTEGER      :: pr_proc, pr_procr
  INTEGER, ALLOCATABLE :: pr_procs(:)
CONTAINS

#ifdef MPI
! What needs to be done for parallel job
  SUBROUTINE start_MPI
    IMPLICIT NONE
    INTEGER :: iargc ! external function gives number of command line arguments
    INTEGER :: j
    ! getting in contact with MPI
    CALL MPI_INIT( ierr )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr)
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr)
!    PRINT *,'nprocs=',nprocs,'myrank =',myrank 
    master = 0
    write(cpuID,'(I3)') myrank
    ! Get command-line arguments
    IF (myrank .EQ. master) THEN
       nargs = iargc()
       IF (nargs .GT. 4) nargs = nargs-4  ! Seems that MPI adds 4 additional arguments which we
       ALLOCATE (argv(nargs))                                 ! wouldn't like to parse
       WRITE(*,'(A,I2)') 'nargs=', nargs
       DO j=1,nargs
          CALL getarg(j, argv(j))
          WRITE(*,'(A,A)') 'argi=', TRIM(argv(j))
       ENDDO
    ENDIF
    ! Send the number of arguments to other nodes
    CALL MPI_BCAST(nargs, 1, MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    IF (myrank .NE. master) THEN
       ALLOCATE (argv(nargs))  ! Only now we can allocate correct size of array
    ENDIF
    ! Send all arguments to other nodes
    CALL MPI_BCAST(argv, nargs*clen, MPI_CHARACTER, master, MPI_COMM_WORLD,ierr)
  END SUBROUTINE start_MPI

  SUBROUTINE stop_MPI
    CALL MPI_FINALIZE(ierr)
  ENDSUBROUTINE stop_MPI

  SUBROUTINE FilenameMPI(infout)
    CHARACTER(LEN=*) :: infout
    infout    = TRIM(infout)//"."//trim(ADJUSTL(cpuID))
  ENDSUBROUTINE FilenameMPI

  SUBROUTINE FindMaxK_MPI(kmax)
    INTEGER, intent(inout) :: kmax(3)
    INTEGER :: tkmax(3)
    CALL MPI_ALLREDUCE(kmax, tkmax, 3, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    kmax(:) = tkmax(:)
  END SUBROUTINE FindMaxK_MPI

  SUBROUTINE FilenameMPI2(infout)
    CHARACTER(LEN=*) :: infout
    infout    = TRIM(infout)//trim(ADJUSTL(cpuID))
  ENDSUBROUTINE FilenameMPI2

  SUBROUTINE Gather_procs(pr_procr, pr_procs, nprocs)
    IMPLICIT NONE
    INTEGER, intent(in) :: pr_procr, nprocs
    INTEGER, intent(out):: pr_procs(nprocs)
    INTEGER :: ier
    CALL MPI_GATHER(pr_procr, 1, MPI_INTEGER, pr_procs, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER pr_procs', ierr
  END SUBROUTINE Gather_procs

  SUBROUTINE Reduce_MPI(xwt, w_RHOLM, w_vRHOLM, w_xwt1, sumfft, vsumfft, NRAD, LM_MAX, nat, iff1, iff2, iff3)
    IMPLICIT NONE
    REAL*8, intent(inout) :: xwt
    REAL*8, intent(inout) :: w_RHOLM(1:NRAD,1:LM_MAX,1:nat), w_vRHOLM(1:NRAD,1:LM_MAX,1:nat)
    REAL*8, intent(inout) :: w_xwt1(0:21,1:nat)
    COMPLEX*16, intent(inout) :: sumfft(iff1,iff2,iff3), vsumfft(iff1,iff2,iff3)
    INTEGER, intent(in)   :: NRAD, LM_MAX, nat, iff1, iff2, iff3
    !
    REAL*8 :: w_xwt
    REAL*8, allocatable :: ww_RHOLM(:,:,:), ww_vRHOLM(:,:,:)
    REAL*8, allocatable :: ww_xwt1(:,:)
    COMPLEX*16, allocatable :: w_sumfft(:,:,:), w_vsumfft(:,:,:)

    CALL MPI_REDUCE(xwt, w_xwt, 1,               MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 1', ierr
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 2', ierr
    
    if (myrank.eq.master) then
       xwt = w_xwt
    endif

    ALLOCATE( ww_RHOLM(1:NRAD,1:LM_MAX,1:nat) )
    CALL MPI_REDUCE(w_RHOLM, ww_RHOLM, NRAD*LM_MAX*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 4', ierr
    if (myrank.eq.master)w_RHOLM(:,:,:) = ww_RHOLM(:,:,:)
    DEALLOCATE( ww_RHOLM )

    ALLOCATE( ww_vRHOLM(1:NRAD,1:LM_MAX,1:nat) )
    CALL MPI_REDUCE(w_vRHOLM, ww_vRHOLM, NRAD*LM_MAX*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 5', ierr
    if (myrank.eq.master)w_vRHOLM(:,:,:) = ww_vRHOLM(:,:,:)
    DEALLOCATE( ww_vRHOLM )

    ALLOCATE( ww_xwt1(0:21,1:nat) )
    CALL MPI_REDUCE(w_xwt1, ww_xwt1, 22*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 6', ierr
    if (myrank.eq.master)w_xwt1(:,:) = ww_xwt1(:,:)
    DEALLOCATE( ww_xwt1 )

    ALLOCATE( w_sumfft(iff1,iff2,iff3) )
    CALL MPI_REDUCE(sumfft, w_sumfft,   iff1*iff2*iff3, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 11', ierr
    if (myrank.eq.master)sumfft(:,:,:) = w_sumfft(:,:,:)
    DEALLOCATE( w_sumfft )
    
    ALLOCATE( w_vsumfft(iff1,iff2,iff3) )
    CALL MPI_REDUCE(vsumfft, w_vsumfft,   iff1*iff2*iff3, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 12', ierr
    if (myrank.eq.master)vsumfft(:,:,:) = w_vsumfft(:,:,:)
    DEALLOCATE( w_vsumfft )
  END SUBROUTINE Reduce_MPI
  
#else

! What needs to be done for serial job

  SUBROUTINE start_MPI
    IMPLICIT NONE
    INTEGER :: iargc ! external function gives number of command line arguments
    INTEGER :: j
    myrank=0
    master=0
    nprocs=1
    ! Get command-line arguments
    nargs = iargc()
    ALLOCATE (argv(nargs))
    DO j=1,nargs
       CALL getarg(j, argv(j))
    ENDDO
    cpuID='0'

    print *, 'nprocs=', nprocs
    
  END SUBROUTINE start_MPI

  SUBROUTINE stop_MPI
  ENDSUBROUTINE stop_MPI

  SUBROUTINE FilenameMPI2(infout)
    CHARACTER(LEN=*) :: infout
  ENDSUBROUTINE FilenameMPI2

  SUBROUTINE FilenameMPI(infout)
    CHARACTER(LEN=*) :: infout
  ENDSUBROUTINE FilenameMPI

  SUBROUTINE Reduce_MPI(xwt, w_RHOLM, w_vRHOLM, w_xwt1, sumfft, vsumfft, &
          &NRAD, LM_MAX, nat, iff1, iff2, iff3)
    IMPLICIT NONE
    REAL*8, intent(inout) :: xwt
    REAL*8, intent(inout) :: w_RHOLM(1:NRAD,1:LM_MAX,1:nat), w_vRHOLM(1:NRAD,1:LM_MAX,1:nat)
    REAL*8, intent(inout) :: w_xwt1(0:21,1:nat)
    COMPLEX*16, intent(inout) :: sumfft(iff1,iff2,iff3), vsumfft(iff1,iff2,iff3)
    INTEGER, intent(in)   :: NRAD, LM_MAX, nat, iff1, iff2, iff3
  END SUBROUTINE Reduce_MPI

  SUBROUTINE FindMaxK_MPI(kmax)
    INTEGER, intent(inout) :: kmax(3)
  END SUBROUTINE FindMaxK_MPI

  SUBROUTINE Gather_procs(pr_procr, pr_procs, nprocs)
    IMPLICIT NONE
    INTEGER, intent(in) :: pr_procr, nprocs
    INTEGER, intent(out):: pr_procs(nprocs)
    pr_procs(1) = pr_procr
  END SUBROUTINE Gather_procs

#endif
END MODULE com_mpi

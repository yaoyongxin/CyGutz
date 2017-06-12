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

  SUBROUTINE FindMax_MPI(max_bands, nbandsk, pr_proc)
    INTEGER, intent(out) :: max_bands
    INTEGER, intent(in)  :: pr_proc
    INTEGER, intent(in)  :: nbandsk(pr_proc)
    ! locals
    INTEGER :: maxb
    maxb=1
    DO i=1,pr_proc
       maxb = max(maxb,nbandsk(i))
    ENDDO
    !print *, 'max_bands on', myrank, '=', maxb
    CALL MPI_ALLREDUCE(maxb, max_bands, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    !print *, 'total number of maxbands=', max_bands
    ! MPI_ALLREDUCE(sendbuf, recbuf, count, MPI_TYPE, MPI_OP, comm, ierr)
  END SUBROUTINE FindMax_MPI

  SUBROUTINE FindMaxK_MPI(kmax)
    INTEGER, intent(inout) :: kmax(3)
    INTEGER :: tkmax(3)
    CALL MPI_ALLREDUCE(kmax, tkmax, 3, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    kmax(:) = tkmax(:)
  END SUBROUTINE FindMaxK_MPI

  SUBROUTINE Reduce_MPI0(Olapm, maxdim, ncix, max_nbands)
    IMPLICIT NONE
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix)
    INTEGER, intent(inout) :: max_nbands
    INTEGER, intent(in) :: maxdim, ncix
    !
    COMPLEX*16 :: w_Olapm(maxdim,maxdim,ncix)
    INTEGER :: w_max_nbands
    
    CALL MPI_ALLREDUCE(Olapm, w_Olapm, maxdim*maxdim*ncix, MPI_DOUBLE_COMPLEX,  MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 1', ierr
    Olapm(:,:,:) = w_Olapm(:,:,:)

    CALL MPI_ALLREDUCE(max_nbands, w_max_nbands, 1, MPI_INTEGER,  MPI_MAX, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 2', ierr
    max_nbands = w_max_nbands
  END SUBROUTINE Reduce_MPI0

  SUBROUTINE FilenameMPI2(infout)
    CHARACTER(LEN=*) :: infout
    infout    = TRIM(infout)//trim(ADJUSTL(cpuID))
  ENDSUBROUTINE FilenameMPI2

  SUBROUTINE Reduce_MPI1(pr_proc, nomega, nume, Ebmax, Ebmin, weib, gtot, numex, recomputeEF)
    use muzero, ONLY: nkp, wgh, nemm, zEk, Ek, max_nbands
    IMPLICIT NONE
    INTEGER, intent(in)   :: pr_proc, nomega, nume, numex, recomputeEF
    REAL*8, intent(inout) :: Ebmin(numex), Ebmax(numex), weib(numex)
    COMPLEX*16, intent(inout) :: gtot(nomega)
    !
    INTEGER :: w_nkp
    INTEGER, allocatable    :: w_nemm(:,:)
    COMPLEX*16, allocatable :: w_zEk(:,:,:), w_gtot(:)
    REAL*8, allocatable     :: w_Ek(:,:), w_wgh(:), w_Ebmin(:), w_Ebmax(:), w_eib(:)
    INTEGER :: i, ikp, iom
    
    CALL MPI_REDUCE(nkp, w_nkp, 1, MPI_INTEGER, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 1', ierr
    if (myrank.eq.master) then
       nkp = w_nkp
    endif

    ALLOCATE( w_wgh(pr_proc*nprocs) )
    CALL MPI_GATHER(wgh, pr_proc, MPI_REAL8, w_wgh, pr_proc, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 4', ierr
    if (myrank.eq.master) then
       DEALLOCATE( wgh )
       ALLOCATE( wgh(nkp) )
       wgh(:nkp) = w_wgh(:nkp)
    endif
    DEALLOCATE( w_wgh )

    ALLOCATE( w_nemm(3,pr_proc*nprocs) )
    CALL MPI_GATHER(nemm, pr_proc*3, MPI_INTEGER, w_nemm, pr_proc*3, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 1', ierr
    if (myrank.eq.master) then
       DEALLOCATE( nemm )
       ALLOCATE( nemm(3,nkp) )
       nemm(:,:nkp) = w_nemm(:,:nkp)
    endif
    DEALLOCATE( w_nemm )

    if (myrank.eq.master) print *, 'Data-size pr_processor[MBy]=', 16.*max_nbands*nomega*pr_proc/(1024**2)
    if (myrank.eq.master) print *, 'Total datasize is[MBy]=', 16.*max_nbands*nomega*pr_proc*nprocs/(1024**2)
    if (myrank.eq.master) WRITE(*,'(A,I4,1x,A,I4,1x,A,I4,1x,A,I4,1x)') 'max_nbands=', max_nbands, 'nomega=', nomega, 'pr_proc=', pr_proc, 'nprocs=', nprocs
    
    ALLOCATE( w_zEk(max_nbands,nomega,pr_proc*nprocs) )
    CALL MPI_GATHER(zEk, max_nbands*nomega*pr_proc, MPI_DOUBLE_COMPLEX, w_zEk, max_nbands*nomega*pr_proc, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 2', ierr
    if (myrank.eq.master) then
       DEALLOCATE( zEk )
       ALLOCATE( zEk(max_nbands,nomega,nkp) )
       zEk(:,:,:nkp) = w_zEk(:,:,:nkp)
    endif
    DEALLOCATE( w_zEk )
    
    ALLOCATE( w_Ek(nume,pr_proc*nprocs) )
    CALL MPI_GATHER(Ek, nume*pr_proc, MPI_REAL8, w_Ek, nume*pr_proc, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 3', ierr
    if (myrank.eq.master) then
       DEALLOCATE( Ek )
       ALLOCATE( Ek(nume,nkp) )
       Ek(:,:nkp) = w_Ek(:,:nkp)
    endif
    DEALLOCATE( w_Ek )

    IF (recomputeEF.GT.1) THEN
       ALLOCATE( w_gtot(nomega) )
       CALL MPI_REDUCE(gtot,w_gtot, nomega, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
       if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 4', ierr
       if (myrank.eq.master) then
          gtot(:) = w_gtot(:)
       endif
       DEALLOCATE( w_gtot )
       
       ALLOCATE( w_Ebmax(numex) )
       CALL MPI_REDUCE(Ebmax,w_Ebmax, numex, MPI_DOUBLE_PRECISION, MPI_MAX, master, MPI_COMM_WORLD, ierr)
       if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 5', ierr
       if (myrank.eq.master) then
          Ebmax(:) = w_Ebmax(:)
       endif
       DEALLOCATE( w_Ebmax )
       
       ALLOCATE( w_Ebmin(numex) )
       CALL MPI_REDUCE(Ebmin,w_Ebmin, numex, MPI_DOUBLE_PRECISION, MPI_MIN, master, MPI_COMM_WORLD, ierr)
       if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 6', ierr
       if (myrank.eq.master) then
          Ebmin(:) = w_Ebmin(:)
       endif
       DEALLOCATE( w_Ebmin )
       
       ALLOCATE( w_eib(numex) )
       CALL MPI_REDUCE(weib,w_eib, numex, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, ierr)
       if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 7', ierr
       if (myrank.eq.master) then
          weib(:) = w_eib(:)
       endif
       DEALLOCATE( w_eib )
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  END SUBROUTINE Reduce_MPI1

  SUBROUTINE Reduce_MPI2(nume,nkpt,sumw0)
    IMPLICIT NONE
    INTEGER :: nume,nkpt
    REAL*8  :: sumw0
    !
    REAL*8  :: sumw0_final
    INTEGER :: nkpt_final,nume_final

    if (vector_para) then
        CALL MPI_ALLREDUCE(nume, nume_final, 1, MPI_INTEGER,  MPI_SUM, MPI_COMM_WORLD, ierr)
        if (ierr.ne.0) print *, 'ERROR in MPI_ALLREDUCE 0', ierr
        nume=nume_final
        CALL MPI_ALLREDUCE(nkpt, nkpt_final, 1, MPI_INTEGER,  MPI_SUM, MPI_COMM_WORLD, ierr)
        if (ierr.ne.0) print *, 'ERROR in MPI_ALLREDUCE 0', ierr
        nkpt=nkpt_final
        CALL MPI_ALLREDUCE(sumw0, sumw0_final, 1, MPI_REAL8,  MPI_SUM, MPI_COMM_WORLD, ierr)
        if (ierr.ne.0) print *, 'ERROR in MPI_ALLREDUCE 1', ierr
        sumw0=sumw0_final
    endif

  END SUBROUTINE Reduce_MPI2

  SUBROUTINE Gather_procs(pr_procr, pr_procs, nprocs)
    IMPLICIT NONE
    INTEGER, intent(in) :: pr_procr, nprocs
    INTEGER, intent(out):: pr_procs(nprocs)
    INTEGER :: ier
    CALL MPI_GATHER(pr_procr, 1, MPI_INTEGER, pr_procs, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER pr_procs', ierr
  END SUBROUTINE Gather_procs

  SUBROUTINE Bcast_MPI(EF)
    IMPLICIT NONE
    REAL*8, intent(inout) :: EF
    CALL MPI_BCAST(EF, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST 1', ierr
  END SUBROUTINE Bcast_MPI
  
  SUBROUTINE Reduce_MPI(xwt, etot, w_RHOLM, w_vRHOLM, w_xwt1, w_xwteh, w_xwtel, w_xwt1h, w_xwt1l, sumfft, vsumfft, nomega, NRAD, LM_MAX, nat, iff1, iff2, iff3)
    IMPLICIT NONE
    REAL*8, intent(inout) :: xwt, etot
    REAL*8, intent(inout) :: w_RHOLM(1:NRAD,1:LM_MAX,1:nat), w_vRHOLM(1:NRAD,1:LM_MAX,1:nat)
    REAL*8, intent(inout) :: w_xwt1(0:21,1:nat), w_xwteh(0:3,1:nat), w_xwtel(0:3,1:nat), w_xwt1h(0:3,1:nat), w_xwt1l(0:3,1:nat)
    COMPLEX*16, intent(inout) :: sumfft(iff1,iff2,iff3), vsumfft(iff1,iff2,iff3)
    INTEGER, intent(in)   :: nomega, NRAD, LM_MAX, nat, iff1, iff2, iff3
    !
    REAL*8 :: w_xwt, w_etot
    REAL*8, allocatable :: ww_RHOLM(:,:,:), ww_vRHOLM(:,:,:)
    REAL*8, allocatable :: ww_xwt1(:,:), ww_xwteh(:,:), ww_xwtel(:,:), ww_xwt1h(:,:), ww_xwt1l(:,:)
    COMPLEX*16, allocatable :: w_sumfft(:,:,:), w_vsumfft(:,:,:)

    CALL MPI_REDUCE(xwt, w_xwt, 1,               MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 1', ierr
    CALL MPI_REDUCE(etot,w_etot,1,               MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 2', ierr
    
    if (myrank.eq.master) then
       xwt = w_xwt
       etot=w_etot
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

    ALLOCATE( ww_xwteh(0:3,1:nat) )
    CALL MPI_REDUCE(w_xwteh, ww_xwteh, 4*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 7', ierr
    if (myrank.eq.master)w_xwteh(:,:) = ww_xwteh(:,:)
    DEALLOCATE( ww_xwteh )

    ALLOCATE( ww_xwtel(0:3,1:nat) )
    CALL MPI_REDUCE(w_xwtel, ww_xwtel, 4*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 8', ierr
    if (myrank.eq.master)w_xwtel(:,:) = ww_xwtel(:,:)
    DEALLOCATE( ww_xwtel )

    ALLOCATE( ww_xwt1h(0:3,1:nat) )
    CALL MPI_REDUCE(w_xwt1h, ww_xwt1h, 4*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 9', ierr
    if (myrank.eq.master)w_xwt1h(:,:) = ww_xwt1h(:,:)
    DEALLOCATE( ww_xwt1h )

    ALLOCATE( ww_xwt1l(0:3,1:nat) )
    CALL MPI_REDUCE(w_xwt1l, ww_xwt1l, 4*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 10', ierr
    if (myrank.eq.master)w_xwt1l(:,:) = ww_xwt1l(:,:)
    DEALLOCATE( ww_xwt1l )

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

  SUBROUTINE FindMax_MPI(max_bands, nbandsk, pr_proc)
    INTEGER, intent(out) :: max_bands
    INTEGER, intent(in)  :: pr_proc
    INTEGER, intent(in)  :: nbandsk(pr_proc)
    ! locals
    INTEGER :: maxb
    maxb=1
    DO i=1,pr_proc
       maxb = max(maxb,nbandsk(i))
    ENDDO
    max_bands = maxb
  END SUBROUTINE FindMax_MPI

  SUBROUTINE Reduce_MPI0(Olapm, maxdim, ncix, max_nbands)
    IMPLICIT NONE
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix)
    INTEGER, intent(in) :: maxdim, ncix, max_nbands
  END SUBROUTINE Reduce_MPI0
  
  SUBROUTINE Reduce_MPI1(pr_proc, nomega, nume, Ebmax, Ebmin, weib, gtot, numex, recomputeEF)
    use muzero, ONLY: nkp, wgh, nemm, zEk, Ek, max_nbands
    IMPLICIT NONE
    INTEGER, intent(in) :: pr_proc, nomega, nume, numex, recomputeEF
    REAL*8, intent(inout) :: Ebmin(numex), Ebmax(numex), weib(numex)
    COMPLEX*16, intent(inout) :: gtot(nomega)
  END SUBROUTINE Reduce_MPI1

  SUBROUTINE Reduce_MPI(xwt, etot, w_RHOLM, w_vRHOLM, w_xwt1, w_xwteh, w_xwtel, w_xwt1h, w_xwt1l, sumfft, vsumfft, nomega, NRAD, LM_MAX, nat, iff1, iff2, iff3)
    IMPLICIT NONE
    REAL*8, intent(inout) :: xwt, etot
    REAL*8, intent(inout) :: w_RHOLM(1:NRAD,1:LM_MAX,1:nat), w_vRHOLM(1:NRAD,1:LM_MAX,1:nat)
    REAL*8, intent(inout) :: w_xwt1(0:21,1:nat), w_xwteh(0:3,1:nat), w_xwtel(0:3,1:nat), w_xwt1h(0:3,1:nat), w_xwt1l(0:3,1:nat)
    COMPLEX*16, intent(inout) :: sumfft(iff1,iff2,iff3), vsumfft(iff1,iff2,iff3)
    INTEGER, intent(in)   :: nomega, NRAD, LM_MAX, nat, iff1, iff2, iff3
  END SUBROUTINE Reduce_MPI

  SUBROUTINE FindMaxK_MPI(kmax)
    INTEGER, intent(inout) :: kmax(3)
  END SUBROUTINE FindMaxK_MPI

  SUBROUTINE Bcast_MPI(EF)
    IMPLICIT NONE
    REAL*8, intent(inout) :: EF
  END SUBROUTINE Bcast_MPI

  SUBROUTINE Gather_procs(pr_procr, pr_procs, nprocs)
    IMPLICIT NONE
    INTEGER, intent(in) :: pr_procr, nprocs
    INTEGER, intent(out):: pr_procs(nprocs)
    pr_procs(1) = pr_procr
  END SUBROUTINE Gather_procs

#endif
END MODULE com_mpi

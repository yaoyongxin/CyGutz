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
  
  SUBROUTINE Reduce_MPI(gloc, gtot, gmloc, Olapm, Eimpm, norbitals, nomega, maxdim, ncix)                                                                     
    IMPLICIT NONE                                                                                                                                             
    COMPLEX*16, intent(inout) :: gloc(norbitals,nomega)                                                                                                       
    COMPLEX*16, intent(inout) :: gtot(nomega)                                                                                                                 
    COMPLEX*16, intent(inout) :: gmloc(maxdim, maxdim, ncix, nomega)                                                                                          
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix), Eimpm(maxdim,maxdim,ncix)                                                                         
    INTEGER, intent(in)       :: norbitals, nomega, maxdim, ncix                                                                                              
    ! locals                                                                                                                                                  
    COMPLEX*16, allocatable :: wgtot(:)                                                                                                                       
    COMPLEX*16, allocatable :: wgloc(:,:)                                                                                                                     
    COMPLEX*16, allocatable :: wOlapm(:,:,:), wEimpm(:,:,:)                                                                                                   
    COMPLEX*16, allocatable :: wgmloc(:,:,:,:)                                                                                                                
                                                                                                                                                              
    if (myrank.eq.master) then                                                                                                                                
       ALLOCATE( wgtot(nomega), wgloc(norbitals,nomega), wOlapm(maxdim,maxdim,ncix), wEimpm(maxdim,maxdim,ncix) )                                             
    endif
    CALL MPI_REDUCE(gtot,  wgtot,  nomega,           MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)                                               
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 1', ierr
    CALL MPI_REDUCE(gloc,  wgloc,  norbitals*nomega, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)                                               
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 2', ierr
    CALL MPI_REDUCE(Olapm, wOlapm, maxdim*maxdim*ncix, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)                                             
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 3', ierr
    CALL MPI_REDUCE(Eimpm, wEimpm, maxdim*maxdim*ncix, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)                                             
    if (myrank.eq.master) then                                                                                                                                
       gtot = wgtot                                                                                                                                           
       gloc = wgloc                                                                                                                                           
       Olapm = wOlapm                                                                                                                                         
       Eimpm = wEimpm                                                                                                                                         
       DEALLOCATE( wgtot, wgloc, wOlapm, wEimpm )                                                                                                             
    endif
                                                                                                                                                              
    if (myrank.eq.master) then                                                                                                                                
       ALLOCATE( wgmloc(maxdim, maxdim, ncix, nomega) )                                                                                                       
    endif
    CALL MPI_REDUCE(gmloc, wgmloc, maxdim*maxdim*ncix*nomega, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)                                      
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 4', ierr
    if (myrank.eq.master) then                                                                                                                                
       gmloc = wgmloc                                                                                                                                         
       DEALLOCATE( wgmloc )                                                                                                                                   
    endif
    !MPI_REDUCE(SNDBUF, RECVBUF, COUNT, DATATYPE, OPERATOR, ROOT, COMM, IERR)                                                                                 
    ! OPERATOR = MPI_SUM                                                                                                                                      
  END SUBROUTINE Reduce_MPI                                                                                                                                   

  SUBROUTINE AllReduce_MPI(Olapm, maxdim, ncix)                                                                     
    IMPLICIT NONE                                                                                                                                             
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix)
    INTEGER, intent(in)       :: maxdim, ncix                                                                                              
    ! locals                                                                                                                                                  
    COMPLEX*16, allocatable :: wOlapm(:,:,:)
    
    ALLOCATE( wOlapm(maxdim,maxdim,ncix) )
    
    CALL MPI_ALLREDUCE(Olapm, wOlapm, maxdim*maxdim*ncix, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 3', ierr
    Olapm = wOlapm
    DEALLOCATE( wOlapm )
    
    !MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
    !<type> SENDBUF(*), RECVBUF(*)
    !INTEGER COUNT, DATATYPE, OP, COMM, IERROR 

  END SUBROUTINE AllReduce_MPI

     
  SUBROUTINE Gather_MPI(tEk, tnbands, tnemin, Ekp, nbandsk, nemink, pr_proc, nprocs, numkpt, nume, max_bands, nomega, nsymop)
    IMPLICIT NONE
    INTEGER, intent(in)     :: pr_proc, nprocs, numkpt, nume, nomega, nsymop, max_bands
    COMPLEX*16, intent(out) :: tEk(nomega,nsymop,pr_proc*nprocs,max_bands)
    INTEGER, intent(out)    :: tnbands(pr_proc*nprocs), tnemin(pr_proc*nprocs)
    COMPLEX*16, intent(in)  :: Ekp(nomega,nsymop,pr_proc,nume)
    INTEGER, intent(in)     :: nbandsk(pr_proc), nemink(pr_proc)
    ! locals
    INTEGER :: i
    
    CALL MPI_GATHER(nbandsk, pr_proc, MPI_INTEGER, tnbands, pr_proc, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 1', ierr
    
    CALL MPI_GATHER(nemink,  pr_proc, MPI_INTEGER, tnemin,  pr_proc, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 2', ierr
    
    DO i=1,max_bands
       CALL MPI_GATHER(Ekp(:,:,:,i),  nomega*nsymop*pr_proc, MPI_DOUBLE_COMPLEX, tEk(:,:,:,i),  nomega*nsymop*pr_proc, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
       if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 3', ierr
    enddo
    !MPI_GATHER(SNDBUF, SCOUNT, DATATYPE, RECVBUF, RCOUNT, RDATATYPE, ROOT, COMM, IERR)
    ! MPI_DOUBLE_COMPLEX
    ! MPI_2DOUBLE_PRECISION
  END SUBROUTINE Gather_MPI

  SUBROUTINE FilenameMPI2(infout)
    CHARACTER(LEN=*) :: infout
    infout    = TRIM(infout)//trim(ADJUSTL(cpuID))
  ENDSUBROUTINE FilenameMPI2


!****************************************************************************
      subroutine zsum_master_mpi(a,n)
      integer n
      complex(8) a(n)
! local
      integer ierr
      complex(8) b(n)

      b=0
      call mpi_reduce(a,b,n,mpi_double_complex,mpi_sum,master, &
            &mpi_comm_world,ierr)
      a=b
      return

      end subroutine zsum_master_mpi


!****************************************************************************
      SUBROUTINE DSUM_ALL_MPI(A,N)
      INTEGER N
      REAL(8) A(N)
! LOCAL
      INTEGER IERR
      REAL(8),ALLOCATABLE::B(:)

      ALLOCATE(B(N)); B=0
      CALL MPI_ALLREDUCE(A,B,N,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
      A=B; DEALLOCATE(B)
      RETURN

      END SUBROUTINE DSUM_ALL_MPI

!**************************************************
      SUBROUTINE DSUM_MASTER_MPI(A,N)
      INTEGER N
      REAL(8) A(N)
! LOCAL
      INTEGER IERR
      REAL(8),ALLOCATABLE::B(:)

      ALLOCATE(B(N)); B=0
      CALL MPI_REDUCE(A,B,N,MPI_DOUBLE_PRECISION,MPI_SUM,MASTER,MPI_COMM_WORLD,IERR)
      IF(MYRANK.EQ.MASTER)A=B
      DEALLOCATE(B)
      RETURN

      END SUBROUTINE DSUM_MASTER_MPI

!*****************************************************************************
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
      SUBROUTINE ISUM_MASTER_MPI(I,N)
      INTEGER N,I(N)
! LOCAL
      INTEGER IERR
      INTEGER,ALLOCATABLE::J(:)

      ALLOCATE(J(N)); J=0
      CALL MPI_REDUCE(I,J,N,MPI_INTEGER,MPI_SUM,MASTER,MPI_COMM_WORLD,IERR)
      IF(MYRANK.EQ.MASTER)I=J
      DEALLOCATE(J)
      RETURN

      END SUBROUTINE ISUM_MASTER_MPI

!*****************************************************************************
      SUBROUTINE IMAX1_ALL_MPI(I)
      INTEGER I
! LOCAL
      INTEGER J(1)

      J(1)=I; CALL IMAX_ALL_MPI(J,1); I=J(1)
      RETURN

      END SUBROUTINE IMAX1_ALL_MPI

!********************************************************************************************
      SUBROUTINE IMAX_ALL_MPI(I,N)
      INTEGER N,I(N)
! LOCAL
      INTEGER IERR,MAXI(N)

      CALL MPI_ALLREDUCE(I,MAXI,N,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
      I=MAXI
      RETURN

      END SUBROUTINE IMAX_ALL_MPI
 
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
  END SUBROUTINE start_MPI

  SUBROUTINE stop_MPI
  ENDSUBROUTINE stop_MPI

  SUBROUTINE FilenameMPI(infout)
    CHARACTER(LEN=*) :: infout
  ENDSUBROUTINE FilenameMPI

  SUBROUTINE FilenameMPI2(infout)
    CHARACTER(LEN=*) :: infout
    infout    = TRIM(infout)//trim(ADJUSTL(cpuID))
  ENDSUBROUTINE FilenameMPI2

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

  SUBROUTINE Reduce_MPI(gloc, gtot, gmloc, Olapm, Eimpm, norbitals, nomega, maxdim, ncix)
    IMPLICIT NONE
    COMPLEX*16, intent(inout) :: gloc(norbitals,nomega)
    COMPLEX*16, intent(inout) :: gtot(nomega)
    COMPLEX*16, intent(inout) :: gmloc(maxdim, maxdim, ncix, nomega)
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix), Eimpm(maxdim,maxdim,ncix)
    INTEGER, intent(in)       :: norbitals, nomega, maxdim, ncix
    ! locals
  END SUBROUTINE Reduce_MPI

  SUBROUTINE AllReduce_MPI(Olapm, maxdim, ncix)                                                                     
    IMPLICIT NONE                                                                                                                                             
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix)
    INTEGER, intent(in)       :: maxdim, ncix                                                                                              
  END SUBROUTINE AllReduce_MPI
  
  SUBROUTINE Gather_MPI(tEk, tnbands, tnemin, Ekp, nbandsk, nemink, pr_proc, nprocs, numkpt, nume, max_bands, nomega, nsymop)
    IMPLICIT NONE
    INTEGER, intent(in)     :: pr_proc, nprocs, numkpt, nume, nomega, nsymop, max_bands
    COMPLEX*16, intent(out) :: tEk(nomega,nsymop,pr_proc*nprocs,max_bands)
    INTEGER, intent(out)    :: tnbands(pr_proc*nprocs), tnemin(pr_proc*nprocs)
    COMPLEX*16, intent(in)  :: Ekp(nomega,nsymop,pr_proc,nume)
    INTEGER, intent(in)     :: nbandsk(pr_proc), nemink(pr_proc)
    tnbands = nbandsk
    tnemin = nemink
    tEk(:,:,:,:) = Ekp(:,:,:,:max_bands)
  END SUBROUTINE Gather_MPI

#endif
END MODULE com_mpi

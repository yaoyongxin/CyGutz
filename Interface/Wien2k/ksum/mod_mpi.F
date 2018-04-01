!**************************************************!
!******* § copyright by kristjan haule, 2002 ******!
!**************************************************!
module com_mpi

#ifdef mpi_mode
  include 'mpif.h'
#endif

  integer :: myrank    ! processor id
  integer :: nprocs    ! # of all processors awailable
  integer :: ierr      ! returned error code
  integer :: master    ! # of master processor
  character*3 :: cpuid ! number of cpu in string representation
  integer, parameter :: clen = 100      ! length of characters in command-line
  character*100, allocatable :: argv(:) ! command-line arguments
  integer      :: nargs                 ! number of command line arguments
  logical      :: vector_para           ! is this parallel wien2k run?
  integer      :: vectors(20,3)        ! which vector files should be read
  integer      :: nvector
  character*200:: vecfn(4)
  character*200:: fvectors(20,6)
contains

#ifdef mpi_mode
! what needs to be done for parallel job
  subroutine start_mpi
    implicit none
    integer :: iargc ! external function gives number of command line arguments
    integer :: j
    ! getting in contact with mpi
    call mpi_init( ierr )
    call mpi_comm_size( mpi_comm_world, nprocs, ierr)
    call mpi_comm_rank( mpi_comm_world, myrank, ierr)
!    print *,'nprocs=',nprocs,'myrank =',myrank 
    master = 0
    write(cpuid,'(i3)') myrank
    ! get command-line arguments
    if (myrank .eq. master) then
       nargs = iargc()
       if (nargs .gt. 4) nargs = nargs-4  ! seems that mpi adds 4 additional arguments which we
       allocate (argv(nargs))                                 ! wouldn't like to parse
       write(*,'(a,i2)') 'nargs=', nargs
       do j=1,nargs
          call getarg(j, argv(j))
          write(*,'(a,a)') 'argi=', trim(argv(j))
       enddo
    endif
    ! send the number of arguments to other nodes
    call mpi_bcast(nargs, 1, mpi_integer, master, mpi_comm_world,ierr)
    if (myrank .ne. master) then
       allocate (argv(nargs))  ! only now we can allocate correct size of array
    endif
    ! send all arguments to other nodes
    call mpi_bcast(argv, nargs*clen, mpi_character, master, mpi_comm_world,ierr)
  end subroutine start_mpi

  subroutine stop_mpi
    call mpi_finalize(ierr)
  endsubroutine stop_mpi

  subroutine filenamempi(infout)
    character(len=*) :: infout
    infout    = trim(infout)//"."//trim(adjustl(cpuid))
  endsubroutine filenamempi

  subroutine findmax_mpi(max_bands, nbandsk, pr_proc)
    integer, intent(out) :: max_bands
    integer, intent(in)  :: pr_proc
    integer, intent(in)  :: nbandsk(pr_proc)
    ! locals
    integer :: maxb
    maxb=1
    do i=1,pr_proc
       maxb = max(maxb,nbandsk(i))
    enddo
    !print *, 'max_bands on', myrank, '=', maxb
    call mpi_allreduce(maxb, max_bands, 1, mpi_integer, mpi_max, mpi_comm_world, ierr)
    !print *, 'total number of maxbands=', max_bands
    ! mpi_allreduce(sendbuf, recbuf, count, mpi_type, mpi_op, comm, ierr)
  end subroutine findmax_mpi
  
  subroutine filenamempi2(infout)
    character(len=*) :: infout
    infout    = trim(infout)//trim(adjustl(cpuid))
  endsubroutine filenamempi2


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
      subroutine dsum_all_mpi(a,n)
      integer n
      real(8) a(n)
! local
      integer ierr
      real(8),allocatable::b(:)

      allocate(b(n)); b=0
      call mpi_allreduce(a,b,n,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
      a=b; deallocate(b)
      return

      end subroutine dsum_all_mpi

!**************************************************
      subroutine dsum_master_mpi(a,n)
      integer n
      real(8) a(n)
! local
      integer ierr
      real(8),allocatable::b(:)

      allocate(b(n)); b=0
      call mpi_reduce(a,b,n,mpi_double_precision,mpi_sum,master,mpi_comm_world,ierr)
      if(myrank.eq.master)a=b
      deallocate(b)
      return

      end subroutine dsum_master_mpi

!*****************************************************************************
      subroutine isum_all_mpi(i,n)
      integer n,i(n)
! local
      integer ierr
      integer,allocatable::j(:)

      allocate(j(n)); j=0
      call mpi_allreduce(i,j,n,mpi_integer,mpi_sum,mpi_comm_world,ierr)
      i=j; deallocate(j)
      return

      end subroutine isum_all_mpi

!**************************************************
      subroutine isum_master_mpi(i,n)
      integer n,i(n)
! local
      integer ierr
      integer,allocatable::j(:)

      allocate(j(n)); j=0
      call mpi_reduce(i,j,n,mpi_integer,mpi_sum,master,mpi_comm_world,ierr)
      if(myrank.eq.master)i=j
      deallocate(j)
      return

      end subroutine isum_master_mpi

!*****************************************************************************
      subroutine imax1_all_mpi(i)
      integer i
! local
      integer j(1)

      j(1)=i; call imax_all_mpi(j,1); i=j(1)
      return

      end subroutine imax1_all_mpi

!********************************************************************************************
      subroutine imax_all_mpi(i,n)
      integer n,i(n)
! local
      integer ierr,maxi(n)

      call mpi_allreduce(i,maxi,n,mpi_integer,mpi_max,mpi_comm_world,ierr)
      i=maxi
      return

      end subroutine imax_all_mpi
 
#else

! what needs to be done for serial job

  subroutine start_mpi
    implicit none
    integer :: iargc ! external function gives number of command line arguments
    integer :: j
    myrank=0
    master=0
    nprocs=1
    ! get command-line arguments
    nargs = iargc()
    allocate (argv(nargs))
    do j=1,nargs
       call getarg(j, argv(j))
    enddo
    cpuid='0'
  end subroutine start_mpi

  subroutine stop_mpi
  endsubroutine stop_mpi

  subroutine filenamempi(infout)
    character(len=*) :: infout
  endsubroutine filenamempi

  subroutine filenamempi2(infout)
    character(len=*) :: infout
    infout    = trim(infout)//trim(adjustl(cpuid))
  endsubroutine filenamempi2

  subroutine findmax_mpi(max_bands, nbandsk, pr_proc)
    integer, intent(out) :: max_bands
    integer, intent(in)  :: pr_proc
    integer, intent(in)  :: nbandsk(pr_proc)
    ! locals
    integer :: maxb
    maxb=1
    do i=1,pr_proc
       maxb = max(maxb,nbandsk(i))
    enddo
    max_bands = maxb
  end subroutine findmax_mpi

#endif
end module com_mpi

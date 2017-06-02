program vectorso2hdf5
    use gprec
    use ghdf5_base
    use hdf5
    implicit none

    character*10 fname
    integer nat
    integer,parameter::iu=7

    call gh5_init()
    call get_command_argument(1, fname)

    open(iu,file=trim(adjustl(fname))//'.struct',status='old')
    read(iu,*)
    read(iu,'(27X,I3)')nat
    close(iu)

    call filetohdf5(trim(adjustl(fname))//'.vectorso', &
            &trim(adjustl(fname))//"_vectorso.h5",nat,iu)
    call filetohdf5(trim(adjustl(fname))//'.vectorsodn', &
            &trim(adjustl(fname))//"_vectorsodn.h5",nat,iu)
    call gh5_end()


    contains

    subroutine filetohdf5(fname,fname_h5,nat,iu)
    character(*)fname,fname_h5
    integer nat,iu

    integer nv,ne,i,ikp
    real(8) sk(3),weight
    character*10 bname
    integer,parameter::lomax=3,nloat=3,lmax2=10
    real(8),allocatable::emist(:,:),elo(:,:,:),ee(:)
    complex(8),allocatable::a(:,:)
    integer,allocatable::kv(:,:),num(:)

    call gh5_open_w(fname_h5,f_id)
    open(iu,file=fname,form='unformatted',status='old',access='sequential')

    allocate(emist(0:lmax2,nat),elo(0:lomax,nloat,nat))
    do i=1,nat
        read(iu)emist(0:lmax2,i); read(iu)elo(0:lomax,:,i)
    enddo
    call gh5_write(emist,lmax2+1,nat,'/emist_list',f_id)
    call gh5_write(elo,lomax+1,nloat,nat,'/elo_list',f_id)
    deallocate(emist,elo)

    ikp=0
    do 
        read(iu,end=101)sk,bname,nv,ne,weight
        ikp=ikp+1
        call gh5_create_group("/IKP_"//trim(int_to_str(ikp)),f_id)
        call gh5_write(sk,3,"/IKP_"//trim(int_to_str(ikp))//"/sk",f_id)
        call gh5_write(bname,10,"/IKP_"//trim(int_to_str(ikp))//"/bname",f_id)
        call gh5_write(nv,"/IKP_"//trim(int_to_str(ikp))//"/nv",f_id)
        call gh5_write(ne,"/IKP_"//trim(int_to_str(ikp))//"/ne",f_id)
        call gh5_write(weight,"/IKP_"//trim(int_to_str(ikp))//"/weight",f_id)

        allocate(kv(3,nv))
        read(iu)kv
        call gh5_write(kv,3,nv,"/IKP_"//trim(int_to_str(ikp))//"/KV",f_id)
        deallocate(kv)
        
        allocate(num(ne),ee(ne),a(nv,ne))
        do i=1,ne
            read(iu)num(i),ee(i)
            read(iu)a(:,i)
        enddo
        call gh5_write(num,ne,"/IKP_"//trim(int_to_str(ikp))//"/num",f_id)
        call gh5_write(ee,ne,"/IKP_"//trim(int_to_str(ikp))//"/ee",f_id)
        call gh5_write(a,nv,ne,"/IKP_"//trim(int_to_str(ikp))//"/A",f_id)
        deallocate(num,ee,a)
    enddo

101 continue
    call gh5_write(ikp,"/kptdim",f_id)
    close(iu)
    call gh5_close(f_id)
    return

    end subroutine filetohdf5


    function int_to_str(i)
    integer i
    character(len = 77) int_to_str

    write(int_to_str,*)i
    int_to_str = adjustl(int_to_str)
    return

    end function int_to_str

end program vectorso2hdf5

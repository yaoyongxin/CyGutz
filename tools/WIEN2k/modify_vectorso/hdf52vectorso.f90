program hdf52vectorso
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

    call hdf5tofile(trim(adjustl(fname))//'.vectorso', &
            &trim(adjustl(fname))//"_vectorso.h5",nat,iu)
    call hdf5tofile(trim(adjustl(fname))//'.vectorsodn', &
            &trim(adjustl(fname))//"_vectorsodn.h5",nat,iu)
    call gh5_end()


    contains

    subroutine hdf5tofile(fname,fname_h5,nat,iu)
    character(*)fname,fname_h5
    integer nat,iu

    integer nv,ne,i,ikp,nkpts
    real(8) sk(3),weight
    character*10 bname
    integer,parameter::lomax=3,nloat=3,lmax2=10
    real(8),allocatable::emist(:,:),elo(:,:,:),ee(:)
    complex(8),allocatable::a(:,:)
    integer,allocatable::kv(:,:),num(:)

    call gh5_open_r(fname_h5,f_id)
    open(iu,file=fname,form='unformatted',status='replace',access='sequential')

    allocate(emist(0:lmax2,nat),elo(0:lomax,nloat,nat))
    call gh5_read(emist,lmax2+1,nat,'/emist_list',f_id)
    call gh5_read(elo,lomax+1,nloat,nat,'/elo_list',f_id)
    do i=1,nat
        write(iu)emist(0:lmax2,i); write(iu)elo(0:lomax,:,i)
    enddo
    deallocate(emist,elo)

    call gh5_read(nkpts,"/kptdim",f_id)
    do ikp=1,nkpts
        call gh5_read(sk,3,"/IKP_"//trim(int_to_str(ikp))//"/sk",f_id)
        call gh5_read(bname,10,"/IKP_"//trim(int_to_str(ikp))//"/bname",f_id)
        call gh5_read(nv,"/IKP_"//trim(int_to_str(ikp))//"/nv",f_id)
        call gh5_read(ne,"/IKP_"//trim(int_to_str(ikp))//"/ne",f_id)
        call gh5_read(weight,"/IKP_"//trim(int_to_str(ikp))//"/weight",f_id)
        write(iu)sk,bname,nv,ne,weight

        allocate(kv(3,nv))
        call gh5_read(kv,3,nv,"/IKP_"//trim(int_to_str(ikp))//"/KV",f_id)
        write(iu)kv
        deallocate(kv)
        
        allocate(num(ne),ee(ne),a(nv,ne))
        call gh5_read(num,ne,"/IKP_"//trim(int_to_str(ikp))//"/num",f_id)
        call gh5_read(ee,ne,"/IKP_"//trim(int_to_str(ikp))//"/eep",f_id)
        call gh5_read(a,nv,ne,"/IKP_"//trim(int_to_str(ikp))//"/AP",f_id)
        do i=1,ne
            write(iu)num(i),ee(i)
            write(iu)a(:,i)
        enddo
        deallocate(num,ee,a)
    enddo

101 continue
    close(iu)
    call gh5_close(f_id)
    return

    end subroutine hdf5tofile


    function int_to_str(i)
    integer i
    character(len = 77) int_to_str

    write(int_to_str,*)i
    int_to_str = adjustl(int_to_str)
    return

    end function int_to_str

end program hdf52vectorso

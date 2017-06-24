program testrotylm
    implicit none
    integer l,i,iu
    real(8) rot(3,3)

    iu=8
    open(iu,file='param.inp',status='old')
    read(iu,*)l
    do i=1,3
        read(iu,*)rot(i,:)
    enddo
    close(iu)
    call get_rotmat_ylm(rot, l)

end program testrotylm

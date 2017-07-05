program testlxyz
    implicit none
    integer l

    open(7,file='param.inp',status='old')
    read(7,*)l
    call SETUP_L(l)


end program testlxyz

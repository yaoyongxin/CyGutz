program testcoulombmatrix
    integer l
    real(8) u,j

    open(7,file='param.inp',status='old')
    read(7,*)l,u,j
    close(7)

    call Vcalc(l,u,j)

end program testcoulombmatrix

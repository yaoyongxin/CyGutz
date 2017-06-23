program testylm
    integer lmax
    real(8) v(3)
   
    open(7,file='param.inp',status='old')
    read(7,*)lmax
    read(7,*)v
    call ylm(v,lmax)


end program testylm

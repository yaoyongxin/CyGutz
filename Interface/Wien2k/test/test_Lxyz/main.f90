program testlxyz
    implicit none
    integer l

    open(7,file='param.inp',status='old')
    read(7,*)l
    call calc_lxyz(l)


    contains

    subroutine calc_lxyz(l)
    integer l

    integer ml1,ml2,j1,j2,j3
    complex(8) dima
    complex(8) lx(-l:l,-l:l),ly(-l:l,-l:l),lz(-l:l,-l:l), &
            &lp(-l:l,-l:l),ln(-l:l,-l:l),tmp(-l:l,-l:l)

    lx=0; ly=0; lz=0
    dima=(0.0d0,1.0d0)
    do ml1=-L,L
        Lz(ml1,ml1)=dfloat(ml1)
        do ml2=-L,L
            if((ml1+1).eq.ml2)then
                j1=L-ml2+1
                j2=L+ml2
                j3=j1*j2
                Lx(ml1,ml2)=sqrt(dfloat(j3))/2.d0
                Ly(ml1,ml2)=sqrt(dfloat(j3))/(2.d0*dima)
            endif
            if((ml1-1).eq.ml2)then
                j1=L+ml2+1
                j2=L-ml2
                j3=j1*j2
                Lx(ml1,ml2)=sqrt(dfloat(j3))/2.d0
                Ly(ml1,ml2)=-sqrt(dfloat(j3))/(2.d0*dima)
            endif
        enddo
    enddo
    
    lp=lx+ly*dima
    ln=lx-ly*dima

    ! Check commutation.
    tmp = (matmul(lx,ly)-matmul(ly,lx))/dima

    open(8,file='w_lxyz.txt',status='replace')
    call write_matrix(lx(-l:l,-l:l),2*l+1,8)
    call write_matrix(ly(-l:l,-l:l),2*l+1,8)
    call write_matrix(lz(-l:l,-l:l),2*l+1,8)
    call write_matrix(lp(-l:l,-l:l),2*l+1,8)
    call write_matrix(ln(-l:l,-l:l),2*l+1,8)
    call write_matrix(tmp(-l:l,-l:l),2*l+1,8)
    close(8)

    end subroutine calc_lxyz


    subroutine write_matrix(a,n,iu)
    integer n,iu
    complex(8) a(n,n)

    integer m1,m2

    do m1=1,n
        do m2=1,n
            write(iu,'(f12.6)',advance='no')real(a(m1,m2))
        enddo
        write(iu,*)
    enddo
    write(iu,*)
    do m1=1,n
        do m2=1,n
            write(iu,'(f12.6)',advance='no')aimag(a(m1,m2))
        enddo
        write(iu,*)
    enddo
    write(iu,*)
    return

    end subroutine write_matrix


end program testlxyz

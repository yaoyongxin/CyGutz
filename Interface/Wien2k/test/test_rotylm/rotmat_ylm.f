! given a proper rotation matrix in 3d space and angular momentum l,
! it prints the rotation matrix in complex spherical harmonics basis.
  
subroutine get_rotmat_ylm(rot,l)
    implicit none
    integer::l
    real(8)::rot(3,3)
      
    integer dimlm
    real(8) a,b,c,det
    complex(8) rotylm(2*l+1,2*l+1)
      
    det=1 ! proper rotation
    call euler(rot,a,b,c)
    dimlm=2*l+1
    call dmat(l,a,b,c,det,rotylm,dimlm)  
    call wrt_rotylm(rotylm,dimlm)
    return 
      
end subroutine get_rotmat_ylm
 

subroutine wrt_rotylm(rotylm,dimlm)
    implicit none
    integer dimlm
    complex(8) rotylm(dimlm,dimlm)

    integer::iu=8
    integer i,j

    open(iu,file='w_rotylm.txt',status='replace')
    do i=1,dimlm
        do j=1,dimlm
            write(iu,'(2f10.6,2x)',advance='no')rotylm(i,j)
        enddo
        write(iu,*)
    enddo
    close(iu)
    return

end subroutine wrt_rotylm

  
subroutine euler(rot,a,b,c)
    !*************************************************************************
    !*************************************************************************
    ! %% this subroutine calculates the euler angles a, b and c of rot.  %%
    ! %% the result are stored in a,b,c. (same as in src_lapwdm/euler.f) %%
    !*************************************************************************
    !*************************************************************************
      
    implicit none
    real(kind=8) :: a,aa,b,bb,c,cc,zero,pi,y_norm
    real(kind=8), dimension(3,3) :: rot, rot_temp
    real(kind=8), dimension(3) :: z,zz,y,yy,yyy,pom,x,xx
    integer :: i,j
    ! definition of the constants
    zero=0d0
    pi=acos(-1d0)
    ! definition of rot_temp=id
    do i=1,3
        do j=1,3
            rot_temp(i,j)=0
            if (i.eq.j) rot_temp(i,i)=1
        enddo
    enddo
    ! initialization of y=e_y, z=e_z, yyy and zz
    do j=1,3
        y(j)=rot_temp(j,2)
        yyy(j)=rot(j,2)
        z(j)=rot_temp(j,3)
        zz(j)=rot(j,3)
    enddo
    ! calculation of yy
    call vecprod(z,zz,yy)
    y_norm=dsqrt(dot_product(yy,yy))
    if (y_norm.lt.1d-10) then
        ! if yy=0, this implies that b is zero or pi
        if (abs(dot_product(y,yyy)).gt.1d0) then
            aa=dot_product(y,yyy)/abs(dot_product(y,yyy))
            a=acos(aa)
        else
            a=acos(dot_product(y,yyy))
        endif
        !
        if (dot_product(z,zz).gt.zero) then
            c=zero
            b=zero
            if (yyy(1).gt.zero) a=2*pi-a
        else
            c=a
            a=zero
            b=pi
            if (yyy(1).lt.zero) c=2*pi-c
        endif
    else
        ! if yy is not 0, then b belongs to ]0,pi[
        do j=1,3
            yy(j)=yy(j)/y_norm
        enddo
        !
        aa=dot_product(y,yy)
        bb=dot_product(z,zz)
        cc=dot_product(yy,yyy)
        if (abs(aa).gt.1d0) aa=aa/abs(aa)
        if (abs(bb).gt.1d0) bb=bb/abs(bb)
        if (abs(cc).gt.1d0) cc=cc/abs(cc)
        b=acos(bb)
        a=acos(aa)
        c=acos(cc)
        if (yy(1).gt.zero) a=2*pi-a
        call vecprod(yy,yyy,pom)
        if (dot_product(pom,zz).lt.zero) c=2*pi-c
    endif
      
end subroutine euler
  
  
subroutine vecprod(a,b,c)
    !*************************************************************************
    !*************************************************************************
    ! %% this subroutine calculates the vector product of a and b.       %%
    ! %% the result is stored in c. (same as in src_lapwdm/euler.f)      %%
    !*************************************************************************
    !*************************************************************************
      
    implicit none
    real(kind=8), dimension(3) :: a,b,c
      
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
    return
      
end subroutine vecprod
  
  
  
subroutine dmat(l,a,b,c,det,dd,length)
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%                                                                 %%
    ! %% this subroutine computes the inverse of the matrix of the       %%
    ! %% representation of size (2*l+1) associated to the rotation       %%
    ! %% described by (a,b,c) angles in euler description and with       %%
    ! %% determinant det.                                                %%
    ! %% the obtained matrix is put in the variable dd.                  %%
    ! %%                                                                 %%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    ! definiton of the variables :
    ! ----------------------------
    implicit real*8 (a-h,o-z)
    integer l,m,n,ifac,length
    complex*16 izero,imag, dd
    dimension dd(length,length)

    imag=(0d0,1d0)
    izero=(0d0,0d0)
    pi=acos(-1d0)
      
    do m=-l,l
        do n=-l,l
            call d_matrix(l,m,n,b,dm)
            if (det.lt.-0.5) then
                dd(l+m+1,n+l+1)=(-1)**l*cdexp(imag*n*a) &
                        &*cdexp(imag*m*c)*dm
            else
                dd(l+m+1,n+l+1)=cdexp(imag*n*a) &
                        &*cdexp(imag*m*c)*dm
            end if
        end do
    end do
    return
      
end subroutine dmat
  
  
subroutine d_matrix(l,m,n,b,dm)
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%                                                                 %%
    ! %% this subroutine is called by the subroutine dmat to compute the %%
    ! %% the value of the coefficient dm.                                %%
    ! %%                                                                 %%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    ! definiton of the variables :
    ! ----------------------------
    implicit real*8 (a-h,o-z)
    integer l,m,n,t
      
    sum=0d0
      
    f1=dfloat(ifac(l+m)*ifac(l-m))/&
            &dfloat(ifac(l+n)*ifac(l-n))
      
    do t=0,2*l
        if ((l-m-t).ge.0.and.(l-n-t).ge.0.and.(t+n+m).ge.0) then
            ! general factor
            f2=dfloat(ifac(l+n)*ifac(l-n))/dfloat(ifac(l-m-t) &
                    &*ifac(m+n+t)*ifac(l-n-t)*ifac(t))
            ! factor with sin(b/2)
            if ((2*l-m-n-2*t).eq.0) then
                f3=1.
            else
                f3=(sin(b/2))**(2*l-m-n-2*t)
            end if
            ! factor with cos(b/2)
            if ((2*t+n+m).eq.0) then
                f4=1.
            else
                f4=(cos(b/2))**(2*t+n+m)
            end if
            sum=sum+(-1)**(l-m-t)*f2*f3*f4
        end if
    end do
    dm=sqrt(f1)*sum
    return

end subroutine d_matrix
  
  
integer function ifac(n)
    integer n
  
    ! definiton of the variables :
    if (n.eq.0) then
        ifac=1
    else
        ifac=1
        do j=1,n
            ifac=ifac*j
        end do
    end if
    return 
end function ifac
  

	Subroutine dmat(l,a,b,c,det,DD)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMPLEX*16 izero,imag,dd
	dimension DD(7,7)                       
        INTEGER l,m,n,ifac
	imag=(0d0,1d0)
	izero=(0d0,0d0)
	pi=acos(-1d0)

	do m=-l,l
	do n=-l,l
	call d_matrix(l,m,n,b,dm)
	if (det.lt.-0.5) then
        dd(l+m+1,n+l+1)=(-1)**l*cdexp(imag*n*a)*cdexp(imag*m*c)*dm
	else
        dd(l+m+1,n+l+1)=cdexp(imag*n*a)*cdexp(imag*m*c)*dm
	end if
 3      format(2I3,2f10.6)
	end do
	end do
	do j=1,2*l+1
	end do
 5      format(7(2f6.3,1X))

	end

!_______________________________________________________



	Subroutine d_matrix(l,m,n,b,dm)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER l,m,n,t

 	sum=0d0

        f1=dfloat(ifac(l+m)*ifac(l-m))/ &
           dfloat(ifac(l+n)*ifac(l-n))


	do t=0,2*l

	if ((l-m-t).ge.0.AND.(l-n-t).ge.0.AND.(t+n+m).ge.0) then

	f2=dfloat(ifac(l+n)*ifac(l-n))/dfloat(ifac(l-m-t) &
          *ifac(m+n+t)*ifac(l-n-t)*ifac(t))

	if ((2*l-m-n-2*t).eq.0) then
        f3=1.
	else
	f3=(sin(b/2))**(2*l-m-n-2*t)
	end if
	if ((2*t+n+m).eq.0) then
	f4=1.
	else
	f4=(cos(b/2))**(2*t+n+m)
	end if

!	write(12,*)f1,f2,f3,f4
	sum=sum+(-1)**(l-m-t)*f2*f3*f4
        end if
	end do

	dm=sqrt(f1)*sum
	end

!__________________________________________________________

	Integer Function ifac(n)

	if (n.eq.0) then
	ifac=1
	else
	ifac=1
	do j=1,n
	ifac=ifac*j
	end do
	end if
	end

!___________________________________________________________
	  
      

























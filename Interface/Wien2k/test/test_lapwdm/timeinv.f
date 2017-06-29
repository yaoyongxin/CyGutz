	SUBROUTINE TIMEINV(TMAT,l)
	IMPLICIT REAL*8 (A-H,O-Z)
	COMPLEX*16 TMAT(7,7),TMP(7,7),SUM,CZERO
	DIMENSION TINV(7,7)
	DATA CZERO /(0.D0,0.D0)/, ZERO /0.D0/

	N=2*l+1
        do i=-l,l
        do j=-l,l
        tinv(i+L+1,j+L+1)=zero
        if (i.eq.(-j)) tinv(i+L+1,j+L+1)=(-1)**i
        end do
        end do

	do i=1,N
	do j=1,N
 	tmp(i,j)=tmat(i,j)
	end do
	end do


	do i=1,N
	do j=1,N
	sum=czero
	do k=1,N
	sum=sum+tinv(i,k)*tmp(k,j)
	end do
	tmat(i,j)=sum
	end do
        end do


	end


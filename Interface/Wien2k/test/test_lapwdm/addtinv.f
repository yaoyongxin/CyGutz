	SUBROUTINE ADDTINV(ll,USYM)

        IMPLICIT REAL*8 (A-H,O-Z)
        COMPLEX*16 T(7,7),tmp(7,7),usym(7,7)
        COMPLEX*16 sum,cone,czero

	DATA cone/(1.d0,0.d0)/,czero/(0.d0,0.d0)/

	n=2*ll+1
	do i=1,n
	do j=1,n
	T(i,j)=czero
	end do
	T(i,i)=cone
	end do
 5      format(7f8.4)
	call timeinv(t,ll)

	do i=1,n
	do j=1,n
	sum=czero
	do k=1,n
	do l=1,n
	sum=sum+t(k,i)*t(l,j)*conjg(usym(k,l))
	end do
	end do
	tmp(i,j)=sum
	end do
	end do

	do i=1,n
	do j=1,n
	usym(i,j)=(usym(i,j)+tmp(i,j))/2
	end do
	end do

	end
		

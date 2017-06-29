	SUBROUTINE ADDTINVSO(ll,USYM)

        IMPLICIT REAL*8 (A-H,O-Z)
        COMPLEX*16 T(7,7),tmp(7,7,3),usym(7,7,3)
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
         write(71,*)'sym'
        do i=1,3
        write(71,5)(real(usym(i,j,2)),j=1,3)
        end do
	call timeinv(t,ll)

	do ii=1,3
	do i=1,n
	do j=1,n
	sum=czero
	do k=1,n
	do l=1,n
	sum=sum+t(k,i)*t(l,j)*conjg(usym(k,l,ii))
	end do
	end do
	tmp(i,j,ii)=sum
	end do
	end do
        end do

	write(71,*)'tmp'
	do i=1,3
	write(71,5)(real(tmp(i,j,2)),j=1,3)
	end do

	do i=1,n
	do j=1,n
	usym(i,j,1)=(usym(i,j,1)+tmp(i,j,3))/2
        usym(i,j,3)=(usym(i,j,3)+tmp(i,j,1))/2
        usym(i,j,2)=(usym(i,j,2)-conjg(tmp(j,i,2)))/2
	end do
	end do

	end
		

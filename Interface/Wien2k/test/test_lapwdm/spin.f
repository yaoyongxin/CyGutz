SUBROUTINE SPIN(snn,sx,sy,sz)
  USE param
  USE sporb
        IMPLICIT REAL*8 (A-H,O-Z)
        COMPLEX*16 stm(2,2),CZERO,CONE,snn(3),dmat(2,2)
        COMPLEX*16 SUM,DNEW(2,2)
      DATA CZERO/(0.D0,0.D0)/,CONE/(0.d0,1.d0)/


	CALL EULER(rotinv,a,b,c)
 9	format (3(f14.8,5x))
!write(6,9)((rotinv(i,j),j=1,3),i=1,3)
!write(6,*)'spin abc:',a,b,c
	
	stm(1,1)=cos(b/2)*exp(cone/2*(a+c))
	stm(1,2)=sin(b/2)*exp(cone/2*(c-a))
	stm(2,1)=-sin(b/2)*exp(cone/2*(a-c))
	stm(2,2)=cos(b/2)*exp(-cone/2*(a+c))
!write(6,*)'TRANSF matrix'
!write(6,8)((stm(i,j),i=1,2),j=1,2)

	dmat(1,1)=snn(1)
	dmat(1,2)=snn(2)
	dmat(2,1)=conjg(snn(2))
	dmat(2,2)=snn(3)
!       write(6,*)'SMAT spincoord:'
!       write(6,8)((dmat(i,j),j=1,2),i=1,2)
 8	format (2(2f14.8,5x))
	do i=1,2
	do j=1,2
	sum=czero
	do k=1,2
	do l=1,2
	sum=conjg(stm(i,k))*dmat(k,l)*stm(j,l)+sum
	end do
	end do
	dnew(i,j)=sum
	end do
	end do
!       write(6,*)'SMAT globa:'
!       write(6,8)((dnew(i,j),j=1,2),i=1,2)

	sx=2.d0*dble(dnew(1,2))
	sy=2.d0*aimag(dnew(1,2))
	sz=real(dnew(1,1)-dnew(2,2))

	end

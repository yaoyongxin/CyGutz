SUBROUTINE SYMMETRIZATION(ICASE,LL,usym,umat,IPH)
  USE param
  USE struct
  USE rotat
  USE case
  USE sporb
  USE matpdf
	IMPLICIT REAL*8 (A-H,O-Z)

        COMPLEX*16 UMAT,USYM,czero,sum, &
                   IMAG,upom,ss
	DIMENSION  &
                  USYM(2*LXDOS+1,2*LXDOS+1), &
                  umat(2*lxdos+1,2*lxdos+1,ndif), &
                  upom(2*lxdos+1,2*lxdos+1)

      parameter ( imag=cmplx(0.D0,1.D0))
!      DATA imag /(0.,1.)/
	jatom=iatom(icase)
        write(6,*)'SYMMETRIZATION OF THE DENSITY MATRIX:',lL
        write(6,*)(itr(i,jatom),i=1,iord)

 333    format(5(2e12.5,1x))
	n=2*LL+1
	czero=(0d0,0d0)

    write(0,'(" phase factor:")')
    write(0,*)phase

    write(0,'(" det factor:")')
    write(0,*)det
   

!________________________________________
	if(ll.eq.0) then
	sum=czero
	do mu=1,mult(jatom)
	sum=sum+umat(1,1,mu)
	end do
	usym(1,1)=sum/mult(jatom)

!_________________________________________
        else if (ll.eq.1) then
	do i=1,n
	do j=1,n
	usym(i,j)=czero
	end do
	end do

        do 100 II=1,iord
        do i=1,n
        do j=1,n
        sum=czero
        do k=1,n
        do l=1,n
	ss=matp(k,i,ii)* &
        exp(iph*imag*phase(ii))* &
        conjg(matp(l,j,ii))

	if (det(ii).lt.-.5) then
        sum=sum+ss*conjg(umat(k,l,itr(ii,jatom)))
	else
        sum=sum+ss*umat(k,l,itr(ii,jatom)) 
	end if

        end do
        end do
        usym(i,j)=usym(i,j)+sum/iord
        end do
        end do
 100    continue

!_________________________________________
        elseif (ll.eq.2) then
	do i=1,n
	do j=1,n
	usym(i,j)=czero
	enddo
	enddo

        do 200 ii=1,iord
	do i=1,n
	do j=1,n
	sum=czero
	do k=1,n
	do l=1,n
	ss=matd(k,i,ii)* &
        exp(iph*imag*phase(ii))* &
        conjg(matd(l,j,ii))

        if (det(ii).lt.-.5) then
        sum=sum+ss*conjg(umat(k,l,itr(ii,jatom)))
        else
        sum=sum+ss*umat(k,l,itr(ii,jatom))
        end if

        end do
        end do
        usym(i,j)=usym(i,j)+sum/iord
        end do
        end do
 200    continue

!__________________________________________________
        elseif (ll.eq.3) then
        do i=1,n
        do j=1,n
        usym(i,j)=czero
        end do 
        end do 

        do 300 ii=1,iord
        do i=1,n
        do j=1,n
        sum=czero
        do k=1,n
        do l=1,n
	ss=matf(k,i,ii)* &
        exp(iph*imag*phase(ii))* &
        conjg(matf(l,j,ii))

        if (det(ii).lt.-.5) then
        sum=sum+ss*conjg(umat(k,l,itr(ii,jatom)))
        else
        sum=sum+ss*umat(k,l,itr(ii,jatom))
        end if

        end do
        end do
        usym(i,j)=usym(i,j)+sum/iord
        upom(i,j)=sum
        end do
        end do
 5      format(7f10.6)
 300    continue

	end if
	end

	



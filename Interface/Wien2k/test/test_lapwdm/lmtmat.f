SUBROUTINE lmtmat
  USE param
  USE sporb
  USE matpdf
  USE symop
  USE struct
	IMPLICIT REAL*8 (A-H,O-Z)
        COMPLEX*16 tmat,tmp

	dimension sym(3,3),tmat(7,7), &
        tmp(7,7),t1(3,3),t2(3,3),t3(3,3)

	pi=acos(-1.d0)
        if (iprint.ge.1) then
	write(6,*)
	write(6,*)'SYMMETRY OPERATIONS IN Cartesian coord. system'
	end if
	
	if (iprint.ge.2) then
	write(6,*)'TRANSFORMATION MATRICES IN P,D,F REPRESENTATION:'
        end if

!       a=5.0
!       b=3.0
!       c=3.0
!       t1(1,1)=cos(a)
!t1(1,2)=sin(a)
!t1(1,3)=0
!t1(2,1)=-sin(a)
!t1(2,2)=cos(a)
!t1(2,3)=0
!t1(3,1)=0
!t1(3,2)=0
!t1(3,3)=1
!       t2(1,1)=cos(b)
!       t2(1,2)=0
!       t2(1,3)=sin(b)
!       t2(2,1)=0
!       t2(2,2)=1
!       t2(2,3)=0
!       t2(3,1)=-sin(b)
!       t2(3,2)=0
!       t2(3,3)=cos(b)
!       t3(1,1)=cos(c)
!       t3(1,2)=sin(c)
!       t3(1,3)=0
!       t3(2,1)=-sin(c)
!       t3(2,2)=cos(c)
!       t3(2,3)=0
!       t3(3,1)=0
!       t3(3,2)=0
!       t3(3,3)=1

!do i=1,3
!do j=1,3
!sum=0.d0
!do k=1,3
!do l=1,3
!sum=sum+t1(i,k)*t2(k,l)*t3(l,j)
!end do
!end do
!       smat(i,j,1)=sum
!end do
!end do

        do 100 iio=1,iord
	do iii=1,3
	do jjj=1,3
	sym(iii,jjj)=smat(iii,jjj,iio)
	end do
	end do
	if (iprint.ge.2) then
        write(6,*)'SYM. OP. no.: ',iio
	do jj=1,3
	write(6,5)(sym(jj,ii),ii=1,3)
	end do
	end if
!____________________________________________________________
        dd =sym(1,1)*(sym(2,2)*sym(3,3)- &
                         sym(3,2)*sym(2,3))- &
             sym(1,2)*(sym(2,1)*sym(3,3)- &
                         sym(3,1)*sym(2,3))+ &
             sym(1,3)*(sym(2,1)*sym(3,2)- &
                         sym(2,2)*sym(3,1))
!_____________________________________________________________
	
	if (dd.lt.0d0) then
        do i=1,3
        do j=1,3
        sym(i,j)=-sym(i,j)
        end do
        end do
        end if
 
	
	call euler(sym,a,b,c)
        if (iprint.ge.1) then
   write(6,3)iio,a*180/pi,b*180/pi,c*180/pi
 3  format(i2,2x,'euler angles: a,b,c ',3f6.1)
        end if
	do 50 l=1,3
!        call null(tmat)
        tmat=(0.d0,0.d0)
	call dmat(l,a,b,c,dd,tmat)
	if (det(iio).lt.-.5) call timeinv(tmat,l)
        if (iprint.ge.2) then
	if (l.eq.1) write(6,*)'p-representaion:'
        if (l.eq.2) write(6,*)'d-representaion:'
        if (l.eq.3) write(6,*)'f-representaion:'
        do m=-l,l
         write(6,5)(tmat(l+m+1,l+n+1),n=-l,l)
        end do
        end if
 5      format(7(2f6.3,1X))

	if (l.eq.1) then
	do i=1,3
	do j=1,3
	  matp(i,j,iio)=tmat(i,j)
	end do
	end do
	else if (l.eq.2) then
        do i=1,5
        do j=1,5
          matd(i,j,iio)=tmat(i,j)
        end do
        end do
        else 
        do i=1,7
        do j=1,7
          matf(i,j,iio)=tmat(i,j)
        end do 
        end do
	end if
 50     continue
 100    continue
	end


	subroutine null(a)
	complex*16 a
	dimension a(7,7)
        a=(0d0,0d0)
	end


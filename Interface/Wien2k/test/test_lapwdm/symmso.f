SUBROUTINE SYMMSO(ICASE,LL,usym,umat,IPH)
  USE param
  USE struct
  USE rotat
  USE case
  USE sporb
  USE matpdf
  USE symop
	IMPLICIT REAL*8 (A-H,O-Z)

        COMPLEX*16 UMAT,USYM,czero,sum, &
                   IMAG,ss,utmp
	DIMENSION  &
                  USYM(2*LXDOS+1,2*LXDOS+1,3), &
                  umat(2*lxdos+1,2*lxdos+1,ndif,3), &
                  utmp(2*lxdos+1,2*lxdos+1,ndif,2,2), &
                  index(2,2)

      DATA imag /(0.,1.)/
	jatom=iatom(icase)
        write(6,*)'SYMMETRIZATION OF THE DENSITY MATRIX:',lL
        write(6,*)(itr(i,jatom),i=1,iord)

 333    format(5(2e12.5,1x))
	n=2*LL+1
	czero=(0d0,0d0)

        do 55 mu=1,ndif
	do 55 i=1,n
	do 55 j=1,n
        utmp(i,j,mu,1,1)=umat(i,j,mu,1)
	utmp(i,j,mu,1,2)=umat(i,j,mu,2)
	utmp(i,j,mu,2,1)=dconjg(umat(j,i,mu,2))
        utmp(i,j,mu,2,2)=umat(i,j,mu,3)
 55     continue
	index(1,1)=1
	index(1,2)=2
	index(2,1)=0
        index(2,2)=3
        
	
 5      format(7f10.4)
!________________________________________
        if(ll.eq.0) then
        sum=czero
        do mu=1,mult(jatom)
        sum=sum+umat(1,1,mu,1)
        end do
        usym(1,1,1)=sum/mult(jatom)
        usym(1,1,3)=usym(1,1,1)
	usym(1,1,2)=czero
!_________________________________________
        else if (ll.eq.1) then
	do i=1,n
	do j=1,n
        do ii=1,3
	usym(i,j,ii)=czero
	end do
	end do
        end do

        do 100 II=1,iord
        do isi=1,2
	do isj=1,2
        do i=1,n
        do j=1,n
        sum=czero
        do isk=1,2
	do isl=1,2
        do k=1,n
        do l=1,n
	ss=matp(k,i,ii)* &
        dconjg(matp(l,j,ii))* &
        spmt(isk,isi,ii)* &
        dconjg(spmt(isl,isj,ii))
        sum=sum+ss*utmp(k,l,itr(ii,jatom),isk,isl) 

        end do
        end do
        end do
        end do
        usym(i,j,index(isi,isj))=usym(i,j,index(isi,isj))+sum/iord
        end do
        end do
        end do
        end do

 100    continue

!_________________________________________
        else if (ll.eq.2) then
        do i=1,n
        do j=1,n
        do ii=1,3
        usym(i,j,ii)=czero
        end do
        end do
        end do

        do 200 II=1,iord
        do isi=1,2
        do isj=1,2
        do i=1,n
        do j=1,n
        sum=czero
        do isk=1,2
        do isl=1,2
        do k=1,n
        do l=1,n
        ss=matd(k,i,ii)* &
        dconjg(matd(l,j,ii))* &
        spmt(isk,isi,ii)* &
        dconjg(spmt(isl,isj,ii))
        sum=sum+ss*utmp(k,l,itr(ii,jatom),isk,isl)

        end do
        end do
        end do
        end do
        usym(i,j,index(isi,isj))=usym(i,j,index(isi,isj))+sum/iord
        end do
        end do
        end do
        end do

 200    continue

!_________________________________________
        else if (ll.eq.3) then
        do i=1,n
        do j=1,n
        do ii=1,3
        usym(i,j,ii)=czero
        end do
        end do
        end do

        do 300 II=1,iord
        do isi=1,2
        do isj=1,2
        do i=1,n
        do j=1,n
        sum=czero
        do isk=1,2
        do isl=1,2
        do k=1,n
        do l=1,n
        ss=matf(k,i,ii)* &
        dconjg(matf(l,j,ii))* &
        spmt(isk,isi,ii)* &
        dconjg(spmt(isl,isj,ii))
        sum=sum+ss*utmp(k,l,itr(ii,jatom),isk,isl)

        end do
        end do
        end do
        end do
        usym(i,j,index(isi,isj))=usym(i,j,index(isi,isj))+sum/iord
        end do
        end do
        end do
        end do

 300    continue

	end if
	end

	



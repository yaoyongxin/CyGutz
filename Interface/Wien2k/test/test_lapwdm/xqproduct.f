      subroutine XQproduct(L,X,Usym,iprint)
! XQproduct calculates product                               
! Z=(X+)QTL, where X is the X-operator, X+ its herm.
! conjugate, QTL is symmetrized population matrix in global co-ordinates      
      USE param
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 X(0:LMAX2,-LMAX2:LMAX2,-LMAX2:LMAX2,2,2) 
      complex*16 Usym(7,7,3),czero
      complex*16 Y(14,14),YQTL(14,14),Z(14,14),Zsym(7,7,3)
      complex*16 dz
      czero=(0.d0,0.d0)
      n=2*L+1
      DO 11 MS1=1,2  ! matrix of X-operator
      DO 11 MS2=1,2
      M1=0
      DO 11 ML1=-L,L
      M1=M1+1
      J1=(2*L+1)*(MS1-1)+M1
      M2=0
      DO 11 ML2=-L,L
      M2=M2+1
      J2=(2*L+1)*(MS2-1)+M2
      Y(j1,j2)=X(L,ml1,ml2,ms1,ms2)
   11 CONTINUE
      do 12 i=1,n  ! symmetrized population matrix
       do 12 j=1,n
        Yqtl(i,j)=Usym(i,j,1)
        Yqtl(i,j+n)=Usym(i,j,2)
        Yqtl(i+n,j)=dconjg(Usym(j,i,2))
        Yqtl(i+n,j+n)=Usym(i,j,3)
  12  continue        
! matrix product Z=[hermconj(Yqtl)*Y + Y*hermconjg(Yqtl)]
      do i=1,2*n
       do j=1,2*n
        Z(i,j)=czero
        do k=1,2*n
         Z(i,j)=Z(i,j)+(Y(k,i)*dconjg(Yqtl(k,j))+dconjg(Yqtl(i,k))*Y(j,k))/2.
        enddo
       enddo
      enddo
! decompose Z on upup, updn, dndn blocks
      do i=1,n
       do j=1,n
        Usym(i,j,1)=Z(i,j)
        Usym(i,j,2)=Z(i+n,j)
        Usym(i,j,3)=Z(i+n,j+n)
       enddo
      enddo
      if(iprint.ge.1)then
       write(6,*)' X matrix',L,n
       call printx(2*n,Y)
       write(6,*)' QTL matrix',L,n
       call printx(2*n,Yqtl)
       write(6,*)' Z= [X(Q+) + (Q+)X]/2 matrix',L,n
       call printx(2*n,Z)
      endif
      return
      end

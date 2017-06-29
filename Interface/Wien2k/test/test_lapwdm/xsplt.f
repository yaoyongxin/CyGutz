SUBROUTINE XSPLT(nemin,nemax,mu,L,nd,ISPIN,weight)    
  USE param
  USE xxa
  USE xdos
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION WEIGHT(NUME)                            
      COMMON /RINTEG/  RI_MAT(0:lmax2,nrf,nrf,2,2)
!--------------------------------------------------------------------- 
	do 100 num=nemin,nemax
        ii=0
        do 100 is1=1,ISPIN
	do 100 is2=is1,ISPIN
        ii=ii+1
        do 100 irf1=1,nrf
        do 100 irf2=1,nrf
        do 100 ly=1,2*l+1
          do 100 lpy=1,2*l+1
       ly1=ly+L*L
       lpy1=lpy+L*L
          xqtl(ly,lpy,mu,ii,nd)=xqtl(ly,lpy,mu,ii,nd)+ &
        dconjg(alm(ly1,num,mu,irf1,is1))*alm(lpy1,num,mu,irf2,is2)* &
        ri_mat(L,irf1,irf2,is1,is2)*weight(num)
 100   continue
      RETURN                                                            
      END                                                               

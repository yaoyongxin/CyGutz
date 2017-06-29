      subroutine couplx(X,L,theta,phi)      
! THIS ROUTINE CALCULATES THE MATRIX ELEMENT OF
! spin and orbital momntum operators in crystal systems
! (Lx,Ly,Lz), (Sx,Sy,Sz) and system with magnetization as
! quantization axis (Lksi,Leta,Ldzeta), (Sksi,Seta,Sdzeta)
! then calls subroutine XOPER in which various matrices of
! various X(L,S) operators are calculated by multiplying
! Salpha and Lalpha matrices. P. Novak 2006.
      USE param
      USE struct 
      IMPLICIT REAL*8(A-H,O-Z)
      common /aver/ krad,kls,cx(-20:20,20),iprx
      COMPLEX*16 czero,dima
      complex*16 Lx(-LMAX2:LMAX2,-LMAX2:LMAX2)    &
                ,Ly(-LMAX2:LMAX2,-LMAX2:LMAX2)    &
                ,Lz(-LMAX2:LMAX2,-LMAX2:LMAX2)    &
                ,Lksi(-LMAX2:LMAX2,-LMAX2:LMAX2)  &
                ,Leta(-LMAX2:LMAX2,-LMAX2:LMAX2)  &
                ,Ldzeta(-LMAX2:LMAX2,-LMAX2:LMAX2)
      complex*16 sksi(2,2),seta(2,2),sdzeta(2,2)   &
                ,sx(2,2),sy(2,2),sz(2,2)           &
                ,dsx(14,14),dsy(14,14),dsz(14,14)  &
                ,dsksi(14,14),dseta(14,14),dsdzeta(14,14)  &
                ,dLx(14,14),dLy(14,14),dLz(14,14)  &
                ,dLksi(14,14),dLeta(14,14),dLdzeta(14,14)  &
                ,Bdip(14,14),bz2(14,14)   &
                ,soc(14,14),X(0:LMAX2,-LMAX2:LMAX2,-LMAX2:LMAX2,2,2) 
      common /spinorb/ dsx,dsy,dsz,dsksi,dseta,dsdzeta  &
                ,dLx,dLy,dLz,dLksi,dLeta,dLdzeta  
      real*8 SM
      integer S2
      dimension rot(3,3)
!**********************************************************************
      czero=(0.d0,0.d0)
      dima=(0.0d0,1.0d0)
      S2=1
!  spin variables in coordinate system with diagonal s_dzeta (dzeta||M)
      do 14 ms=1,S2+1
      do 14 msp=1,S2+1
      sksi(msp,ms)=czero
      seta(msp,ms)=czero
      sdzeta(msp,ms)=czero
14    continue
        S=dfloat(S2)/2.d0
         do MS1=1,S2+1
          SM=-S-1.d0
           do MS2=1,S2+1
            SM=SM+1.d0
            asp=0.d0
            asm=0.d0  
            sdzeta(MS1,MS2)=czero
            if((MS1+1).eq.MS2)then
            asp=sqrt((S-SM+1.d0)*(S+SM))
            sksi(MS1,MS2)=asp/2.d0
            seta(MS1,MS2)=-asp/(2.d0*dima)
            endif
            if((MS1-1).eq.MS2)then
            asm=sqrt((S+SM+1.d0)*(S-SM))
            sksi(MS1,MS2)=asm/2.d0
            seta(MS1,MS2)=asm/(2.d0*dima)
            endif
            if(MS1.eq.MS2)sdzeta(MS1,MS2)=SM
            enddo
           enddo
      if(iprx.gt.1)then
      write(6,*)' ********* sksi real part*********'
      do 508 ms=1,S2+1
508   write(6,602)(dble(sksi(ms,ms1)),ms1=1,S2+1)
      write(6,*)' ********* sksi imag part*********'
      do 518 ms=1,S2+1
518   write(6,602)(dimag(sksi(ms,ms1)),ms1=1,S2+1)
      write(6,*)' ********* seta real part*********'
      do 509 ms=1,S2+1
509   write(6,602)(dble(seta(ms,ms1)),ms1=1,S2+1)
      write(6,*)' ********* seta imag part*********'
      do 519 ms=1,S2+1
519   write(6,602)(dimag(seta(ms,ms1)),ms1=1,S2+1)
      write(6,*)' ********* sdzeta *********'
      do 510 ms=1,S2+1
510   write(6,602)(dble(sdzeta(ms,ms1)),ms1=1,S2+1)
!  set the rotation matrix 
      write(6,550)theta,phi
550   format(' theta=',f6.1,' phi=',f6.1)
      endif
! axis ksi in (x,y,z) system
      rot(1,1)=cos(phi)*cos(theta)
      rot(2,1)=sin(phi)*cos(theta)
      rot(3,1)=-sin(theta)
!     write(6,*)' ksi=',(rot(j,1),j=1,3)
! axis dzeta in (x,y,z) system
      rot(1,3)=cos(phi)*sin(theta)
      rot(2,3)=sin(phi)*sin(theta)
      rot(3,3)=cos(theta)
!     write(6,*)' dzt=',(rot(j,3),j=1,3)
! axis eta in (x,y,z) system, eta=[dzeta,ksi]
      rot(1,2)= rot(2,3)*rot(3,1)-rot(3,3)*rot(2,1)  
      rot(2,2)=-rot(1,3)*rot(3,1)+rot(3,3)*rot(1,1) 
      rot(3,2)= rot(1,3)*rot(2,1)-rot(2,3)*rot(1,1)
!     write(6,*)' eta=',(rot(j,2),j=1,3)
      if(iprx.gt.1)then
      write(6,*)' ********* rotation matrix *********'
      do 512 i=1,3
512   write(6,602)(rot(i,i1),i1=1,3)
      endif
!  express sx, sy, sz through sksi, seta, sdzeta
      do 15 msp=1,S2+1
      do 15 ms=1,S2+1
      sx(msp,ms)=rot(1,1)*sksi(msp,ms)+rot(1,2)*seta(msp,ms)+  &
                 rot(1,3)*sdzeta(msp,ms)
      sy(msp,ms)=rot(2,1)*sksi(msp,ms)+rot(2,2)*seta(msp,ms)+  &
                 rot(2,3)*sdzeta(msp,ms)
      sz(msp,ms)=rot(3,1)*sksi(msp,ms)+rot(3,2)*seta(msp,ms)+  &
                 rot(3,3)*sdzeta(msp,ms)
15    continue
      if(iprx.gt.1)then
      write(6,*)' ********* sx real part*********'
      do 608 ms=1,S2+1
608   write(6,602)(dble(sx(ms,ms1)),ms1=1,S2+1)
      write(6,*)' ********* sx imag part*********'
      do 618 ms=1,S2+1
618   write(6,602)(dimag(sx(ms,ms1)),ms1=1,S2+1)
      write(6,*)' ********* sy real part*********'
      do 609 ms=1,S2+1
609   write(6,602)(dble(sy(ms,ms1)),ms1=1,S2+1)
      write(6,*)' ********* sy imag part*********'
      do 619 ms=1,S2+1
619   write(6,602)(dimag(sy(ms,ms1)),ms1=1,S2+1)
      write(6,*)' ********* sz real part*********'
      do 520 ms=1,S2+1
520   write(6,602)(dble(sz(ms,ms1)),ms1=1,S2+1)
      write(6,*)' ********* sz imag part*********'
      do 521 ms=1,S2+1
521   write(6,602)(dimag(sz(ms,ms1)),ms1=1,S2+1)
      endif
      do 13 m=-L,L
      do 13 mp=-L,L
      Lx(mp,m)=czero
      Ly(mp,m)=czero
      Lz(mp,m)=czero
      Lksi(mp,m)=czero
      Leta(mp,m)=czero
      Ldzeta(mp,m)=czero
13    continue
! components of the angular momentum, in crystal co-ordinate system
      do 12 ml1=-L,L
      Lz(ml1,ml1)=dfloat(ml1)
       do 12 ml2=-L,L
         alp=0.d0
         alm=0.d0
         if((ml1+1).eq.ml2)then
         j1=L-ml2+1
         j2=L+ml2
         j3=j1*j2
         Lx(ml1,ml2)=sqrt(float(j3))/2.d0
         Ly(ml1,ml2)=sqrt(float(j3))/(2.d0*dima)
         endif
         if((ml1-1).eq.ml2)then
         j1=L+ml2+1
         j2=L-ml2
         j3=j1*j2
         Lx(ml1,ml2)=sqrt(float(j3))/2.d0
         Ly(ml1,ml2)=-sqrt(float(j3))/(2.d0*dima)
         endif
12    continue
! components of the angular momentum, in (ksi,eta,dzeta) co-ordinate system
      do 16 ml1=-L,L
       do 16 ml2=-L,L
        Lksi(ml1,ml2)=rot(1,1)*Lx(ml1,ml2)+rot(2,1)*Ly(ml1,ml2)+ &
                        rot(3,1)*Lz(ml1,ml2)
        Leta(ml1,ml2)=rot(1,2)*Lx(ml1,ml2)+rot(2,2)*Ly(ml1,ml2)+ &
                        rot(3,2)*Lz(ml1,ml2)
        Ldzeta(ml1,ml2)=rot(1,3)*Lx(ml1,ml2)+rot(2,3)*Ly(ml1,ml2)+ &
                        rot(3,3)*Lz(ml1,ml2)
16    continue
      if(iprx.gt.2)then
      write(6,*)' ********* Lx real part*********'
      do 500 m=-L,L
500   write(6,602)(dble(Lx(m,m1)),m1=-L,L)
      write(6,*)' ********* Lx imag part*********'
      do 530 m=-L,L
530   write(6,602)(dble(Lx(m,m1)),m1=-L,L)
      write(6,*)' ********* Ly real part*********'
      do 900 m=-L,L
900   write(6,602)(dble(Ly(m,m1)),m1=-L,L)
      write(6,*)' ********* Ly imag part*********'
      do 930 m=-L,L
930   write(6,602)(dimag(Ly(m,m1)),m1=-L,L)
602   format(10f7.3)
      write(6,*)' ********* Lz *********'
      do 502 m=-L,L
502   write(6,602)(dble(Lz(m,m1)),m1=-L,L)
      endif
! Spin and orbital momentum matrices in (2S+1)*(2L+1) space. Ordering is
! mS=-1/2,mL=-L, (-1/2,-L+1)...(-1/2,L),(1/2,-L)......(1/2,L)
      DO 11 MS1=1,S2+1
      DO 11 MS2=1,S2+1
      M1=0
      DO 11 ML1=-L,L
      M1=M1+1
      J1=(2*L+1)*(MS1-1)+M1
      M2=0
      DO 11 ML2=-L,L
      M2=M2+1
      J2=(2*L+1)*(MS2-1)+M2
      dsx(j1,j2)=czero
      dsy(j1,j2)=czero
      dsz(j1,j2)=czero
      dsksi(j1,j2)=czero
      dseta(j1,j2)=czero
      dsdzeta(j1,j2)=czero
      dlx(j1,j2)=czero
      dly(j1,j2)=czero
      dlz(j1,j2)=czero
      dlksi(j1,j2)=czero
      dleta(j1,j2)=czero
      dldzeta(j1,j2)=czero
      if(ML1.eq.ML2)then
      dsx(j1,j2)=sx(MS1,MS2)
      dsy(j1,j2)=sy(MS1,MS2)
      dsz(j1,j2)=sz(MS1,MS2)
      dsksi(j1,j2)=sksi(MS1,MS2)
      dseta(j1,j2)=seta(MS1,MS2)
      dsdzeta(j1,j2)=sdzeta(MS1,MS2)
      endif
      if(MS1.eq.MS2)then
      dLx(j1,j2)=Lx(ML1,ML2)
      dLy(j1,j2)=Ly(ML1,ML2)
      dLz(j1,j2)=Lz(ML1,ML2)
      dLksi(j1,j2)=Lksi(ML1,ML2)
      dLeta(j1,j2)=Leta(ML1,ML2)
      dLdzeta(j1,j2)=Ldzeta(ML1,ML2)
      endif
   11 CONTINUE
      n=(2*L+1)*(S2+1)
      if(iprx.gt.1)then
        write(6,*)
        write(6,*)' Matrix elements Lx '
        call printx(n,dlx)
        write(6,*)
        write(6,*)' Matrix elements Ly '
        call printx(n,dly)
        write(6,*)
        write(6,*)' Matrix elements Lz '
        call printx(n,dlz)
        write(6,*)
        write(6,*)' Matrix elements Sx'
        call printx(n,dsx)
        write(6,*)
        write(6,*)' Matrix elements Sy'
        call printx(n,dsy)
        write(6,*)
        write(6,*)' Matrix elements Sz'
        call printx(n,dsz)
       endif
       call Xoper(L,S2,X)
      return
      end

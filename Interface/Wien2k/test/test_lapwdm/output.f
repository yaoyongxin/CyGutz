SUBROUTINE OUTPUT(ICASE,LL,USYM,IBL,COUP,MAGN,SO)
  USE param
  USE struct
  USE case
        IMPLICIT REAL*8 (A-H,O-Z)
        COMPLEX*16 COUP(0:LMAX2,-LMAX2:LMAX2,-LMAX2:LMAX2,2,2)
        common /aver/ krad,kls,cx(-20:20,20),iprx
        COMMON/ANGL/     THETA,PHI

        COMPLEX*16 USYM,DENSMAT,SS(3),czero
        complex*16 uhelp(2*LXDOS+1,2*LXDOS+1,3)
        dimension xsum(3),trace(3)
        DIMENSION USYM(2*LXDOS+1,2*LXDOS+1,3), &
                  DENSMAT(2*LXDOS+1,2*LXDOS+1), &
                  torb(3),index(3)
        CHARACTER*4   :: spsp
      czero=(0.d0,0.d0)
      N=2*LL+1
        index(1)=7
        index(2)=11
        index(3)=8

      if ((iprint.ge.1).or.(iprx.ge.1)) then
       write(6,*)'symmetrized matrix at repre. atom in global coord.'
       DO II=1,IBL
         DO I=1,N
            WRITE(6,567)(REAL(usym(I,J,II)),J=1,N)
         end do
         WRITE(6,556)
         DO I=1,N
            WRITE(6,567)(dimag(usym(I,J,II)),J=1,N)
         end do
         write(6,*)
       END DO
      endif
      if((krad.ne.0).and.(iprx.gt.0))then
       DO II=1,2  
        write(6,*)' X-operator real part, block',ii
         DO I=1,N
          mi=i-ll-1
            WRITE(6,567)(REAL(coup(LL,mi,J-LL-1,II,II)),J=1,N)
         end do
         write(6,*)
         write(6,*)' X-operator imag. part, block',ii
         DO I=1,N
          mi=i-ll-1
            WRITE(6,567)(dimag(coup(LL,mi,J-LL-1,II,II)),J=1,N)
         end do
         write(6,*)
       END DO
      endif
! for <|X|> calculation replace usym by [usym*coup+coup*usym]/2
      if(krad.ne.0)then
       call xqproduct(LL,COUP,Usym,iprx)
       write(6,*)' product Usym*Coupl'
       DO II=1,IBL
         DO I=1,N
            WRITE(6,567)(REAL(usym(I,J,II)),J=1,N)
         end do
         WRITE(6,556)
         DO I=1,N
            WRITE(6,567)(dimag(usym(I,J,II)),J=1,N)
         end do
         write(6,*)
       END DO
      endif
       do i=1,3
        torb(i)=0.d0
        xsum(i)=0.d0
      end do
      if(krad.ne.0)then
       write(6,*)'DENSITY MATRIX IN LOCAL COORD. SYSTEM of repre. atom'
       write(6,*)'spin coordinate system according to case.inso'
      endif
      write(6,*)
  DO 999 II=1,IBL
      IPH=(3-II)*(II-1)
      CALL LOCMAT(ICASE,LL,USYM(1,1,II),DENSMAT)
       write(6,*)' product Usym*Coupl after densmat'
         DO I=1,N
            WRITE(6,567)(REAL(usym(I,J,II)),J=1,N)
         end do
         WRITE(6,556)
         DO I=1,N
            WRITE(6,567)(dimag(usym(I,J,II)),J=1,N)
         end do
         write(6,*)
      if(krad.ne.0)then
       do i=1,n
        ml=i-(2*ll+1)
        xsum(ii)=xsum(ii)+densmat(i,i)
       enddo
      else
      write(6,*)'_______________________________________'
      SS(II)=(0.d0,0.d0)
567  format(7x,7f9.5)
      if (ii.eq.1) then
         spsp='UPUP'
      else if (ii.eq.2) then
         spsp='UPDN'
      else
         spsp='DNDN'
      end if
      write(6,*) spsp,' block:'
      WRITE(21,557) spsp,ll
      write(6,*)'REAL:'
      DO I=1,N
        WRITE(6,567)(REAL(DENSMAT(I,J)),J=1,N)
        WRITE(21,567)(REAL(DENSMAT(I,J)),J=1,N)
        SS(II)=SS(II)+DENSMAT(I,I)
      END DO
        WRITE(6,*)
      write(6,*)'IMAG:'
      write(21,556)
      DO I=1,N
        WRITE(6,567)(DIMAG(DENSMAT(I,J)),J=1,N)
        WRITE(21,567)(DIMAG(DENSMAT(I,J)),J=1,N)
      END DO
 WRITE(6,*)
 WRITE(21,'(/,'':TRA'',i3.3,'':  TRACE of '',a4,'' MATRIX='',2f10.5)')iatom(icase),spsp,SS(II)
 WRITE(6,'(/,'':TRA'',i3.3,'':  TRACE of '',a4,'' MATRIX='',2f10.5)')iatom(icase),spsp,SS(II)
      
      if (IPH.eq.0) then
      CALL ORB(LL,USYM(1,1,II),tx,ty,tz)
      torb(1)=torb(1)+tx
      torb(2)=torb(2)+ty
      torb(3)=torb(3)+tz
   WRITE(21,'('':POM'',i3.3,a2,'' ORBITAL MOMENT in global orthog. system='',3f9.5)') &
         iatom(icase),spsp,tx,ty,tz
   WRITE(6,'('':POM'',i3.3,a2,'':Partial ORBITAL MOMENT in global orthog. system='',3f9.5)') &
         iatom(icase),spsp,tx,ty,tz
      end if

      write(6,*)
      write(6,*)
      write(21,*)
      write(21,*)

! output of density matrix
      
      write(index(ii),100)iatom(icase)
100   format(i5,' atom density matrix')
      write(index(ii),102)ll,tx,ty,tz
102   format(i5,3f10.6, ' L, Lx,Ly,Lz in global orthogonal system')
      do i=1,n
      write(index(ii),101)(densmat(i,j),j=1,n)
      enddo
101   format(2(2es16.8,2x))
!......end special output
      endif

999 CONTINUE
      if(krad.ne.0)then
       xsum(1)=cx(krad,kls)*xsum(1)
       xsum(3)=cx(krad,kls)*xsum(3)
       xtot=xsum(1)+xsum(3)
       if(abs(krad).lt.11)then
        write(6,103)iatom(icase),ll,xsum(1),xsum(3),xtot
        write(21,103)iatom(icase),ll,xsum(1),xsum(3),xtot
       else
        dl=dfloat(2*ll+1)
        write(6,104)iatom(icase),ll,xsum(1)/dl,xsum(3)/dl
        write(21,104)iatom(icase),ll,xsum(1)/dl,xsum(3)/dl
       endif
103   format(':XOP',i3,i3,4f12.5)
104   format(':XOP',i3,i3,4e12.5)
      else
        CALL SPIN(SS,Sx,sy,sz)
        tM=sin(theta)*(cos(phi)*torb(1)+sin(phi)*torb(2))+cos(theta)*torb(3)
        ssM=sin(theta)*(cos(phi)*sX+sin(phi)*sy)+cos(theta)*sz
        write(21,*)
        write(6,*)
        WRITE(6,379)iatom(icase),torb(1),torb(2),torb(3),tM
        WRITE(6,378)iatom(icase),SX,SY,SZ,ssm
        WRITE(21,379)iatom(icase),torb(1),torb(2),torb(3),tM
        WRITE(21,378)iatom(icase),SX,SY,SZ,ssm
 379    format(':ORB',i3.3,':  ORBITAL MOMENT:',3f9.5,' PROJECTION ON M' &
        ,f9.5)
 378    format(':SPI',i3.3,':  SPIN MOMENT:',3f10.5,' PROJECTION ON M' &
        ,f9.5)
557   format(' Density matrix ',a4,' block, real part.  L=',i2)
556   format(' Density matrix, imag part')
      endif

      END

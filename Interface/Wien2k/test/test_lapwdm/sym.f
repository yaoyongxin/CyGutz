!..... progaram sym is used for calculation of transformation parameters
!..... of so wavefunctions. It transforms the symmetry operations
!..... to frame of reference z||M. Variable det=1 ... M not reversed
!..... det=-1 ... M reversed ( must by later coupled to time  inversion ).
!.... Angle phase(I) is used for spin part of transformation matrix 
!..... acting of spinor function.
       
SUBROUTINE SYM(MAGN)
  USE param
  USE struct
  USE sporb
  USE symop
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL MAGN
      COMMON/ANGL/     THETA,PHI

      DIMENSION rot(3,3),trans(3,3)
      COMPLEX*16 imag
      parameter (imag=(0.d0,1.d0))
!     parameter (pi=3.141592654d0) 
      pi = acos(-1.D0)
      ct=cos(theta)
      cf=cos(phi)
      st=sin(theta)
      sf=sin(phi)
      rot(1,1)=cf*ct
      rot(1,2)=-sf
      rot(1,3)=cf*st
      rot(2,1)=sf*ct
      rot(2,2)=cf
      rot(2,3)=sf*st
      rot(3,1)=-st
      rot(3,2)=0
      rot(3,3)=ct
      if (iprint.ge.1) then
      write(6,*)'SPIN COORD SYSTEM:'
      write(6,5)((rot(i,j),j=1,3),i=1,3)
      write(6,*)
      end if
      call INVERSSYMDEF(rot,rotinv)
  write(6,*)'SYMMETRY OPERATIONS IN SPIN COORD. SYSTEM'
  write(6,*)
      do 100 i=1,IORD
         call transform(trans,rotinv,smat(1,1,i),rot)
 5       format(3f12.6)
 	if (iprint.ge.2) then
	WRITE(6,*)'SYMMETRY OPERATION: ',I
	write(6,*)'in global cartesian coordinate system:'
	write(6,5)((smat(ii,j,i),j=1,3),ii=1,3)
        write(6,*)
	write(6,*)'in spin coordinate system:'
        write(6,5)((trans(ii,jj),jj=1,3),ii=1,3)
        write(6,*)
	end if

!.... det(i)=-1 operation must be coupled to time inversion
!.... det(i)=1 operation must not be coupled to time inversion

         det(i)=trans(1,1)*trans(2,2)-trans(1,2)*trans(2,1)

!...DD=-1 operation(i)=rotation(i)*(-I)

         DD=trans(1,1)*(trans(2,2)*trans(3,3)-trans(2,3)*trans(3,2))-&
            trans(1,2)*(trans(2,1)*trans(3,3)-trans(2,3)*trans(3,1))+&
            trans(1,3)*(trans(2,1)*trans(3,2)-trans(2,2)*trans(3,1))
         if (DD.lt.-0.5) then
         do ii=1,3
         do jj=1,3
         trans(ii,jj)=-trans(ii,jj)
         end do
         end do
         end if

!...Euler angles calculated for the rotation(i)
!...phase(i) determines the phase shift between 
!...the spin upand spin dn components (relevant only
!... for updn part of the density matrix

         call euler(trans,a,b,c)
         if (iprint.ge.1) then
   write(6,3)i,a*180/pi,b*180/pi,c*180/pi
         end if
         if (MAGN) then
         if (det(i).gt.0.5) then
         phase(i)=a+c
         else
         phase(i)=a-c
         end if 

	 if (iprint.ge.1) then
	 write(6,*)'# of operation, phase, det:'
         write(6,*)i,phase(i),det(I)
   3     format(i2,1x,'Euler angles: a,b,c:  ',3f7.1)
         write(6,*)
	 end if
	 if (abs(1.-abs(det(i))).gt.1d-2) then
 	 write(6,*)'symm. operation ',i,' so-det=',det(i)
	 STOP
	 end if

         else
         det(i)=1.d0
	 spmt(1,1,i)=exp(imag*(a+c)/2.)*cos(b/2.)
         spmt(1,2,i)=exp(-imag*(a-c)/2.)*sin(b/2.)
	 spmt(2,1,i)=-dconjg(spmt(1,2,i))
	 spmt(2,2,i)=dconjg(spmt(1,1,i))
         if (iprint.ge.2) then
	 do ii=1,2
	 write(6,6)(spmt(ii,jj,i),jj=1,2)
	 end do
	 end if
 6      format(7(2f6.3,1X))
	 end if
         
	 
 100  continue
      end
      
      subroutine transform(T,Pinv,A,P)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION T(3,3),P(3,3),A(3,3),Pinv(3,3)
      
      do i=1,3
         do j=1,3
            sum=0
            do k=1,3
               do l=1,3
                  sum=sum+Pinv(i,k)*A(k,l)*P(l,j)
               end do
            end do
            T(i,j)=sum
         end do
!       write(6,*)'Transf. matrix:',(T(i,k),k=1,3)
      end do
      return
      end

        real*8 function ACOSS(x)
        implicit real*8 (A-H,O-Z)

        if ((abs(x)-1.).gt.1d-4) then
        write(6,*)'x=',x
        stop 'ACOSS ERROR'
        end if
        if (x.ge.1) then
        acoss=0.0
        else if (x.le.-1.) then
        acoss=acos(-1.0)
        else
        acoss=acos(x)
        end if
        end


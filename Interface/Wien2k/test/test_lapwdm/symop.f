      SUBROUTINE symoper
!
!     transforms symmetri matrices to cartesian system
!
!     for .orhto. systems   smat == imat
! 
!     for not.ortho. systems smat == BR1 . imat . BR1^-1 
        USE param
        USE struct
        USE symop
!
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 BR1in(3,3)
      LOGICAL          ORTHO
      COMMON /ORTH/   ORTHO
      COMMON /GENER/  BR1(3,3),BR2(3,3)

!.........inverssymdef is a subroutine in sph-UP.frc and
!.........calculates the inverse of an 3*3 matrix.......
!	write(6,*)'symop called'
      CALL INVERSSYMDEF(BR1,BR1in)

!       write(6,*)'invers end',iord
! 
!.......define symmetry matrices for kartesian system
!
	ior=iord
        do i=1,IORD
          do j=1,3
          do k=1,3
          smat(j,k,i)=0. 
          end do
          end do
        end do

 
        do i=1,IORD
        write(78,*)i
        do ii=1,3
!        write(78,*)(iz(ii,jj,i),jj=1,3)
        end do
          do j=1,3
          do k=1,3
            do l=1,3
            do m=1,3
          smat(j,k,i)= smat(j,k,i) +  &
               BR1(j,l) * iz(l,m,i) * BR1in(m,k)
            end do
            end do
          end do
          end do
        end do
! 
!.......define symmetry matrices for kartesian system
!
       return
 11   FORMAT(3(3I2,F10.5/))
 12   FORMAT(3(3f8.5/))

       end

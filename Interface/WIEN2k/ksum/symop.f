SUBROUTINE symoper
  !
  !     transforms symmetry matrices to cartesian system
  !
  !     for .ortho. systems   opimat == iz
  ! 
  !     for not.ortho. systems opimat == BR1 . iz . BR1^-1 
  !
  USE param
  USE struct
  USE sym2
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
  !.......define symmetry matrices for cartesian system
  !
  ior=iord
  do i=1,IORD
     do j=1,3
        do k=1,3
           opimat(j,k,i)=0. 
        end do
     end do
  end do
 
  do i=1,IORD
     do j=1,3
        do k=1,3
           do l=1,3
              do m=1,3
                 opimat(j,k,i)= opimat(j,k,i) + BR1(j,l) * iz(l,m,i) * BR1in(m,k)
              end do
           end do
        end do
     end do
  end do
  
  return
11 FORMAT(3(3I2,F10.5/))
12 FORMAT(3(3f8.5/))
end SUBROUTINE symoper

        subroutine SymmRot(rotloc,NAT)

        implicit real*8 (a-h,o-z)
!
!---------------------------------------------------------------------  
!  Patch up symmetry issues of rotation matrices
!  L. D. Marks, July 2010
      DOUBLE PRECISION rotloc(3,3,NAT), rold(3,3,NAT)
      rold=rotloc
!
!     Test various integer combinations -- this will catch most things
      DO IRTest = 2,6
        RotTest = 1.D0/sqrt(DBLE(IRTEST))
        do jatom=1,NAT
        do j1=1,3
        do j2=1,3
                test1=abs(rotloc(j2,j1,jatom))-RotTest
                if(abs(test1).lt.1d-6)rotloc(j2,j1,jatom)=sign(RotTest,rotloc(j2,j1,jatom))
        enddo
        enddo
        enddo
     ENDDO
        RotTest = sqrt(0.75D0)
        do jatom=1,NAT
        do j1=1,3
        do j2=1,3
                test1=abs(rotloc(j2,j1,jatom))-RotTest
                if(abs(test1).lt.1d-6)rotloc(j2,j1,jatom)=sign(RotTest,rotloc(j2,j1,jatom))
        enddo
        enddo
        enddo

        RotTest = sqrt(2.D0/3.D0)
        do jatom=1,NAT
        do j1=1,3
        do j2=1,3
                test1=abs(rotloc(j2,j1,jatom))-RotTest
                if(abs(test1).lt.1d-6)rotloc(j2,j1,jatom)=sign(RotTest,rotloc(j2,j1,jatom))
        enddo
        enddo
        enddo

     do Jatom=1,NAT
       write(76,*)rotloc(1:3,1:3,jatom)
       write(75,*)rold  (1:3,1:3,jatom)
     enddo

      return
      end

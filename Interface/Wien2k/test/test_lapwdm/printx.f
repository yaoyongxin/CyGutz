      subroutine printx(n,x)
      complex*16 x(14,14)
            write(6,*)' Real part'
      do j1=1,n
       write(6,629)(dble(x(j1,j2)),j2=1,n)
      enddo
            write(6,*)' Imaginary part'
      do j1=1,n
       write(6,629)(dimag(x(j1,j2)),j2=1,n)
      enddo
        write(6,*)
629   format(14f7.3)
      return
      end

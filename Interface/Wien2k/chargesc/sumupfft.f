SUBROUTINE sumupfft(ifft1,ifft2,ifft3,w,fft,sumfft)
  !
  ! Accumulates fft fields
  ! GM 6/6-00
  IMPLICIT NONE
  INTEGER, intent(in) :: ifft1,ifft2,ifft3
  REAL*8, intent(in)  :: w
  COMPLEX*16, intent(in) :: fft(ifft1,ifft2,ifft3)
  COMPLEX*16, intent(out):: sumfft(ifft1,ifft2,ifft3)
  ! locals
  INTEGER    :: i1,i2,i3
  INTRINSIC DCONJG
  do i3=1,ifft3
     do i2=1,ifft2
        do i1=1,ifft1
           sumfft(i1,i2,i3)=sumfft(i1,i2,i3)+ w*DBLE(FFT(i1,i2,i3)*conjg(FFT(i1,i2,i3)))
        enddo
     enddo
  enddo
  return
end SUBROUTINE sumupfft



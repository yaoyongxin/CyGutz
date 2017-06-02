SUBROUTINE setfft1(nkk,ifft1,ifft2,ifft3,a,fft,keigen)
  !
  ! Sets the fft-fields from the eigenvectors a and k-vectors
  ! GM 6/6-00
  USE defs
  USE param
  IMPLICIT NONE
  !  !_REAL     REAL*8     :: A(nmat)
  !  !_COMPLEX  COMPLEX*16 :: A(nmat)
  COMPLEX*16, intent(in):: A(nmat)
  INTEGER, intent(in)   :: nkk,ifft1,ifft2,ifft3
  INTEGER, intent(in)   :: KEIGEN(3,NMAT)
  COMPLEX*16,intent(out):: FFT(IFFT1,IFFT2,IFFT3)
  ! locals
  INTEGER             :: I, I1, I2, I3
  REAL*8              :: TPI
  
  fft=zeroc
  DO i=1,nkk                                                      
     i1 = keigen(1,I)
     i2 = keigen(2,I)
     i3 = keigen(3,I)
     IF(I1.LT.0) I1 = I1+IFFT1                                         
     IF(I2.LT.0) I2 = I2+IFFT2                                         
     IF(I3.LT.0) I3 = I3+IFFT3
     FFT(I1+1,I2+1,I3+1)=a(i)
  ENDDO
END SUBROUTINE setfft1

SUBROUTINE setfft1c(nkk,ifft1,ifft2,ifft3,a,fft,keigen)
  !
  ! Sets the fft-fields from the eigenvectors a and k-vectors
  !
  USE defs
  USE param
  IMPLICIT NONE
  INTEGER, intent(in)    :: nkk,ifft1,ifft2,ifft3
  COMPLEX*16, intent(in) :: A(nkk)
  INTEGER, intent(in)    :: KEIGEN(3,NMAT)
  COMPLEX*16, intent(out):: FFT(IFFT1,IFFT2,IFFT3)
  ! locals
  INTEGER             :: I, I1, I2, I3
  REAL*8              :: TPI
  !
  fft=zeroc
  DO i=1,nkk                                                      
     i1 = keigen(1,I)
     i2 = keigen(2,I)
     i3 = keigen(3,I)
     IF(I1.LT.0) I1 = I1+IFFT1
     IF(I2.LT.0) I2 = I2+IFFT2
     IF(I3.LT.0) I3 = I3+IFFT3
     FFT(I1+1,I2+1,I3+1)=a(i)
  ENDDO
END SUBROUTINE setfft1c

SUBROUTINE GETFFT(nwave,ifft1,ifft2,ifft3,rho1,fft,kmax)
  !
  !.....sets rho1 from fft-field
  !
  USE param
  use reclat
  IMPLICIT NONE
  INTEGER     :: ifft1,ifft2,ifft3
  INTEGER     :: G(3),NST,STG(3,NSYM),IND(NSYM),kmax(3)
  INTEGER     :: j,jj,index,i1,i2,i3,nwave
  COMPLEX*16  :: FFT(IFFT1,IFFT2,IFFT3),rho1(*)          
  COMPLEX*16  :: TAUP(NSYM)
  index=0
  do j=1,nwave
     g(1)=kzz(1,j)
     g(2)=kzz(2,j)
     g(3)=kzz(3,j)
     call STERN(G,NST,STG,TAUP)
     do jj=1,NST
        index=index+1
        I1=stg(1,jj)
        I2=stg(2,jj)
        I3=stg(3,jj)
        IF(IABS(i1).GT.2*kmax(1)) EXIT ! GOTO 22
        IF(IABS(i2).GT.2*kmax(2)) EXIT ! GOTO 22
        IF(IABS(i3).GT.2*kmax(3)) EXIT ! GOTO 22
        IF(I1.LT.0) I1=I1+IFFT1                                         
        IF(I2.LT.0) I2=I2+IFFT2                                         
        IF(I3.LT.0) I3=I3+IFFT3                                         
        rho1(index)=FFT(I1+1,I2+1,I3+1)                  
     enddo
  enddo
  return
END SUBROUTINE GETFFT






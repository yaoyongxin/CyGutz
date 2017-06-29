SUBROUTINE inispl                                          
  USE param
  USE xdos
  IMPLICIT REAL*8 (A-H,O-Z)
  !
  xqtl(1:2*lxdos+1,1:2*lxdos+1,1:ndif,1:3,1:lmax2+1)=(0.d0,0.d0)
END SUBROUTINE inispl

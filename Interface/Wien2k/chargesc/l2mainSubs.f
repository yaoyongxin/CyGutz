SUBROUTINE Read_Vec_Spin(ikp, E, wgh, As, As_lo, kx, ky, kz, kxlo, kylo, kzlo, bkx, bky, bkz, bkxlo, bkylo, bkzlo, nemin, nemax, more_kpoints, n0, emin, emax, iso, nmat, nume, nnlo)
  use com, ONLY: weigh
  use dmf, ONLY: Qcomplex
  IMPLICIT NONE
  INTEGER, intent(in)  :: ikp
  REAL*8,  intent(out) :: E(nume), wgh !, weight(nume)
  !  !_REAL REAL*8,intent(out)        :: As(nmat,nume,iso)
  !  !_REAL REAL*8,intent(out)        :: As_lo(nnlo,nume,iso)
  !  !_COMPLEX COMPLEX*16,intent(out) :: As(nmat,nume,iso)
  !  !_COMPLEX COMPLEX*16,intent(out) :: As_lo(nnlo,nume,iso)
  !
  COMPLEX*16,intent(out) :: As(nmat,nume,iso)
  COMPLEX*16,intent(out) :: As_lo(nnlo,nume,iso)
  INTEGER, intent(out) :: kx(nmat), ky(nmat), kz(nmat), kxlo(nnlo), kylo(nnlo), kzlo(nnlo)
  REAL*8,  intent(out) :: bkx(nmat), bky(nmat), bkz(nmat), bkxlo(nnlo), bkylo(nnlo), bkzlo(nnlo)
  INTEGER, intent(out) :: nemin, nemax
  LOGICAL, intent(out) :: more_kpoints
  INTEGER, intent(out) :: n0
  REAL*8,  intent(in)  :: emin, emax
  INTEGER, intent(in)  :: iso, nmat, nume, nnlo
  ! locals
  REAL*8  :: As_tmp(nmat)
  INTEGER :: i, is, itape, ios, ne, num, nlov
  REAL*8  :: s,t,z, vnorm1
  CHARACTER*10 :: bname
  
  DO is=1,iso
     itape=8+is
     READ(itape,IOSTAT=ios) s,t,z,bname,n0,ne,wgh
     more_kpoints=.FALSE.
     IF (ios /= 0) EXIT
     more_kpoints=.TRUE.
     
     READ(itape) (kx(i),ky(i),kz(i),i=1,n0-nnlo),(kxlo(i),kylo(i),kzlo(i),i=1,nnlo)
     
     DO i=1,n0-nnlo
        bkx(i)=(s+kx(i))
        bky(i)=(t+ky(i))
        bkz(i)=(z+kz(i))
     ENDDO
     DO i=1,nnlo
        bkxlo(i)=(s+kxlo(i))
        bkylo(i)=(t+kylo(i))
        bkzlo(i)=(z+kzlo(i))
     ENDDO
     nemin=1
     nemax=0
     DO
        READ(itape) num,E(num)
        if (Qcomplex) then
           READ(itape) (As(i,num,is),i=1,n0)
        else
           READ(itape) (As_tmp(i),i=1,n0)
           As(:n0,num,is) = As_tmp(:n0)
        endif
        As_lo(:nnlo,num,is) = As(n0-nnlo+1:n0,num,is)
        IF(e(num).LT.emin) nemin=nemin+1                           
        IF(e(num).LE.emax) nemax=nemax+1     
        IF(num.EQ.ne) EXIT
     ENDDO
  ENDDO     
  RETURN
END SUBROUTINE Read_Vec_Spin

SUBROUTINE DMFT_WEIGHTS(zw2, Aweight, nbands, debug)
  IMPLICIT NONE
  REAL*8,     intent(out)  :: zw2(nbands)
  COMPLEX*16, intent(inout):: Aweight(nbands,nbands)
  INTEGER,    intent(in)   :: nbands
  LOGICAL,    intent(in)   :: debug
  ! locals
  INTEGER    :: imx(nbands), timx, i, j, l, n
  INTEGER    :: idxarr(nbands)
  COMPLEX*16 :: zweight(nbands,nbands), tt(nbands), cc
  INTEGER    ::  lwork, lrwork, liwork
  complex*16, allocatable :: work(:)
  real*8,     allocatable :: rwork(:)
  integer,    allocatable :: iwork(:)
  integer :: info
  
  lwork = 2*nbands+nbands*nbands
  lrwork =  5*nbands + 2*nbands*nbands + 1
  liwork = 3+5*nbands
  ALLOCATE ( work(lwork), rwork(lrwork), iwork(liwork) )
     
  ! Aweight = Aw . zw2 . Aw^+
  zweight(:,:) = Aweight(:,:)
  Aweight(:,:) = -Aweight(:,:)

  CALL ZHEEVD('V','U', nbands, Aweight, nbands, zw2, work, lwork, rwork, lrwork, iwork, liwork, info )

  if (info .ne. 0) then
     print *, 'Diagonalization of weights failed. Info-zheevd=', info
  endif
  
  ! We defines zw2 with minus before to sort eigenvalues from largest to smallest
  ! To make them positive, change sign here
  zw2 = -zw2

  CALL eig_order_abs_val(zw2, idxarr, nbands)

  CALL permute_eigensystem1(idxarr, zw2, Aweight, nbands)
     
  ! Sort degenerate eigenvalues to make it similar to identity
  do i=1,nbands
     timx=1
     do l=2,nbands
        if (abs(Aweight(l,i)).gt.abs(Aweight(timx,i))) timx=l
     enddo
     imx(i)=timx
  enddo
  do i=1,nbands
     do j=i+1,nbands
        if (abs(zw2(i)-zw2(j)).LT.1e-9 .and. imx(i).gt.imx(j)) then
           tt(:) = Aweight(:,i)
           Aweight(:,i) = Aweight(:,j)
           Aweight(:,j) = tt(:)
        endif
     enddo
  enddo
  
  if (debug) then
     ! checking exact diagonalization
     do i=1,nbands
        do j=1,nbands
           cc=0
           do n=1,nbands
              cc = cc + Aweight(i,n)*zw2(n)*conjg(Aweight(j,n))
           enddo
           if (abs(zweight(i,j)-cc).GT.1e-8) print *, 'd',i,j,abs(zweight(i,j)-cc)
        enddo
     enddo
  endif
  DEALLOCATE ( work, rwork, iwork )
  
END SUBROUTINE DMFT_WEIGHTS

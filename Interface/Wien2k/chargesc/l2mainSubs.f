SUBROUTINE Read_Vec_Spin(ikp, E, As, As_lo, kx, ky, kz, vnorm, &
        &kxlo, kylo, kzlo, bkx, bky, bkz, bkxlo, bkylo, bkzlo, &
        &more_kpoints, n0, emin, nemin, iso, nmat, nume, nnlo)
  use dmf, ONLY: Qcomplex
  IMPLICIT NONE
  INTEGER, intent(in)  :: ikp
  REAL*8,  intent(out) :: E(nume), vnorm(nume)
  COMPLEX*16,intent(out) :: As(nmat,nume,iso)
  COMPLEX*16,intent(out) :: As_lo(nnlo,nume,iso)
  INTEGER, intent(out) :: kx(nmat), ky(nmat), kz(nmat), kxlo(nnlo), kylo(nnlo), kzlo(nnlo)
  REAL*8,  intent(out) :: bkx(nmat), bky(nmat), bkz(nmat), bkxlo(nnlo), bkylo(nnlo), bkzlo(nnlo)
  LOGICAL, intent(out) :: more_kpoints
  INTEGER, intent(out) :: n0, nemin
  REAL*8,  intent(in)  :: emin
  INTEGER, intent(in)  :: iso, nmat, nume, nnlo
  ! locals
  REAL*8  :: As_tmp(nmat)
  INTEGER :: i, is, itape, ios, ne, num, nlov
  REAL*8  :: s,t,z, vnorm1, wgh
  CHARACTER*10 :: bname
  
  DO is=1,iso
     itape=8+is
     READ(itape,IOSTAT=ios) s,t,z,bname,n0,ne,wgh
     more_kpoints=.FALSE.
     IF (ios /= 0) EXIT
     more_kpoints=.TRUE.
     
     READ(itape) (kx(i),ky(i),kz(i),i=1,n0-nnlo), &
         &(kxlo(i),kylo(i),kzlo(i),i=1,nnlo)
     
     DO i=1,n0-nnlo
        bkx(i)=s+kx(i)
        bky(i)=t+ky(i)
        bkz(i)=z+kz(i)
     ENDDO
     DO i=1,nnlo
        bkxlo(i)=s+kxlo(i)
        bkylo(i)=t+kylo(i)
        bkzlo(i)=z+kzlo(i)
     ENDDO
     nemin=1
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
        IF(num.EQ.ne) EXIT
     ENDDO
  ENDDO    

  READ(12,202,iostat=ios) (vnorm(i),i=1,ne)
  IF (ios /= 0) THEN
      vnorm=1.d0
  ENDIF

  RETURN

202 FORMAT(4e20.12)
END SUBROUTINE Read_Vec_Spin


SUBROUTINE DMFT_WEIGHTS(zw2, Aweight, nbands)
    IMPLICIT NONE
    REAL*8,     intent(out)  :: zw2(nbands)
    COMPLEX*16, intent(inout):: Aweight(nbands,nbands)
    INTEGER,    intent(in)   :: nbands
    ! locals
    INTEGER    ::  lwork, lrwork
    INTEGER    :: idxarr(nbands)
    complex*16, allocatable :: work(:)
    real*8,     allocatable :: rwork(:)
    integer :: info
    
    lwork = nbands+nbands*nbands
    lrwork =  3*nbands
    ALLOCATE(work(lwork),rwork(lrwork))
       
    Aweight = -Aweight
    call zheev("v","u",nbands,Aweight,nbands,zw2,work,lwork,rwork,info)
    if (info .ne. 0) then
       write(0,"(A,I0)")'Diagonalization of weights failed. Info-zheevd=',info
       stop
    endif
    ! change sign back
    zw2 = -zw2
    DEALLOCATE(work, rwork)
    CALL eig_order_abs_val(zw2, idxarr, nbands)
    CALL permute_eigensystem1(idxarr, zw2, Aweight, nbands)
    return
  
END SUBROUTINE DMFT_WEIGHTS


SUBROUTINE eig_order_abs_val(ev, idxarr, ndim)
  IMPLICIT NONE
!!!-----------------------------------------------------------------!!!
!!! This routine sorts complex eigenvalues of a matrix according to !!!
!!! its real parts with the smallest in the first slot and reorders !!!
!!! the matrices of left (row) and right (column) eigenvectors in a !!!
!!! corresponding manner.                                           !!!
!!!-----------------------------------------------------------------!!!
  !---------- Passed variables ----------
  real*8, intent(in)   :: ev(ndim)         ! Array of eigenvalues
  INTEGER, intent(out) :: idxarr(ndim)     ! Index array which gives proper     order
  INTEGER :: ndim                            ! Dimension of matrices
  !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
  !---------- Parameters ----------
  REAL*8, PARAMETER :: maxval = 1000.d0
  !---------- Local variables ----------
  LOGICAL, ALLOCATABLE :: sorted(:)
  REAL*8,  ALLOCATABLE :: sortonr(:)
  INTEGER :: p
  INTEGER :: q
  INTEGER :: idx
  REAL*8  :: min
  !---------- Allocate dynamic memory storage ----------
  ALLOCATE(sortonr(ndim), sorted(ndim))
  !---------- Initialize arrays ----------
  idxarr = 0
  sortonr = -DBLE(abs(ev))
  sorted = .FALSE.
  !---------- Create index array for real value ----------
  sorted = .FALSE.
  DO p = 1,ndim
     min = maxval
     DO q = 1,ndim
        IF(.NOT.sorted(q).AND.min.GT.sortonr(q)) THEN
           min = sortonr(q)
           idx = q
        ENDIF
     ENDDO
     idxarr(p) = idx
     sorted(idx) = .TRUE.
  ENDDO
  DEALLOCATE(sortonr, sorted)
  RETURN
END SUBROUTINE eig_order_abs_val


SUBROUTINE permute_eigensystem1(idxarr, ev, evr, ndim)
  IMPLICIT NONE
  !---------- Passed variables ----------
  INTEGER, intent(in)       :: idxarr(ndim)     ! Index array which gives       proper order
  real*8, intent(inout)     :: ev(ndim)         ! Array of eigenvalues
  COMPLEX*16, intent(inout) :: evr(ndim,ndim)   ! Matrix of right eigenvectors  (column)
  INTEGER :: ndim                               ! Dimension of matrices
  !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
  !---------- Local variables ------------------
  INTEGER :: p
  COMPLEX*16, ALLOCATABLE :: eval(:)
  COMPLEX*16, ALLOCATABLE :: evec(:,:)
  ALLOCATE(eval(ndim), evec(ndim,ndim))
  !---------- Permute the eigenvalues ----------
  DO p = 1,ndim
     eval(p) = ev(idxarr(p))
  ENDDO
  ev = eval
  !---------- Permute the right eigenvectors ----------
  DO p = 1,ndim
     evec(:,p) = evr(:,idxarr(p))
  ENDDO
  evr = evec
  !---------- Deallocate dynamic memory storage ----------
  DEALLOCATE(eval, evec)
  RETURN
END SUBROUTINE permute_eigensystem1


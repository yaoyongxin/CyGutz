subroutine read_vec_spin(ikp, e, as, as_lo, kx, ky, kz, vnorm, &
        &kxlo, kylo, kzlo, bkx, bky, bkz, bkxlo, bkylo, bkzlo, &
        &more_kpoints, n0, emin, nemin, iso, nmat, nume, nnlo)
  use dmf, only: qcomplex
  implicit none
  integer, intent(in)  :: ikp
  real*8,  intent(out) :: e(nume), vnorm(nume)
  complex*16,intent(out) :: as(nmat,nume,iso)
  complex*16,intent(out) :: as_lo(nnlo,nume,iso)
  integer, intent(out) :: kx(nmat), ky(nmat), kz(nmat), kxlo(nnlo), kylo(nnlo), kzlo(nnlo)
  real*8,  intent(out) :: bkx(nmat), bky(nmat), bkz(nmat), bkxlo(nnlo), bkylo(nnlo), bkzlo(nnlo)
  logical, intent(out) :: more_kpoints
  integer, intent(out) :: n0, nemin
  real*8,  intent(in)  :: emin
  integer, intent(in)  :: iso, nmat, nume, nnlo
  ! locals
  real*8  :: as_tmp(nmat)
  integer :: i, is, itape, ios, ne, num, nlov
  real*8  :: s,t,z, wgh
  character*10 :: bname
  
  do is=1,iso
     itape=8+is
     read(itape,iostat=ios) s,t,z,bname,n0,ne,wgh
     more_kpoints=.false.
     if (ios /= 0) exit
     more_kpoints=.true.
     
     read(itape) (kx(i),ky(i),kz(i),i=1,n0-nnlo), &
         &(kxlo(i),kylo(i),kzlo(i),i=1,nnlo)
     
     do i=1,n0-nnlo
        bkx(i)=s+kx(i)
        bky(i)=t+ky(i)
        bkz(i)=z+kz(i)
     enddo
     do i=1,nnlo
        bkxlo(i)=s+kxlo(i)
        bkylo(i)=t+kylo(i)
        bkzlo(i)=z+kzlo(i)
     enddo
     nemin=1
     do
        read(itape) num,e(num)
        if (qcomplex) then
           read(itape) (as(i,num,is),i=1,n0)
        else
           read(itape) (as_tmp(i),i=1,n0)
           as(:n0,num,is) = as_tmp(:n0)
        endif
        as_lo(:nnlo,num,is) = as(n0-nnlo+1:n0,num,is)
        if(e(num).lt.emin) nemin=nemin+1
        if(num.eq.ne) exit
     enddo
  enddo    

  read(12,202,iostat=ios) (vnorm(i),i=1,ne)
  if (ios /= 0) then
      vnorm=1.d0
  endif

  return

202 format(4e20.12)
end subroutine read_vec_spin


subroutine dmft_weights(zw2, aweight, nbands)
    implicit none
    real*8,     intent(out)  :: zw2(nbands)
    complex*16, intent(inout):: aweight(nbands,nbands)
    integer,    intent(in)   :: nbands
    ! locals
    integer    ::  lwork, lrwork
    integer    :: idxarr(nbands)
    complex*16, allocatable :: work(:)
    real*8,     allocatable :: rwork(:)
    integer :: info
    
    lwork = nbands+nbands*nbands
    lrwork =  3*nbands
    allocate(work(lwork),rwork(lrwork))
       
    aweight = -aweight
    call zheev("v","u",nbands,aweight,nbands,zw2,work,lwork,rwork,info)
    if (info .ne. 0) then
       write(0,"(a,i0)")'diagonalization of weights failed. info-zheevd=',info
       stop
    endif
    ! change sign back
    zw2 = -zw2
    deallocate(work, rwork)
    call eig_order_abs_val(zw2, idxarr, nbands)
    call permute_eigensystem1(idxarr, zw2, aweight, nbands)
    return
  
end subroutine dmft_weights


subroutine eig_order_abs_val(ev, idxarr, ndim)
  implicit none
!!!-----------------------------------------------------------------!!!
!!! this routine sorts complex eigenvalues of a matrix according to !!!
!!! its real parts with the smallest in the first slot and reorders !!!
!!! the matrices of left (row) and right (column) eigenvectors in a !!!
!!! corresponding manner.                                           !!!
!!!-----------------------------------------------------------------!!!
  !---------- passed variables ----------
  real*8, intent(in)   :: ev(ndim)         ! array of eigenvalues
  integer, intent(out) :: idxarr(ndim)     ! index array which gives proper     order
  integer :: ndim                            ! dimension of matrices
  !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
  !---------- parameters ----------
  real*8, parameter :: maxval = 1000.d0
  !---------- local variables ----------
  logical, allocatable :: sorted(:)
  real*8,  allocatable :: sortonr(:)
  integer :: p
  integer :: q
  integer :: idx
  real*8  :: min
  !---------- allocate dynamic memory storage ----------
  allocate(sortonr(ndim), sorted(ndim))
  !---------- initialize arrays ----------
  idxarr = 0
  sortonr = -dble(abs(ev))
  sorted = .false.
  !---------- create index array for real value ----------
  sorted = .false.
  do p = 1,ndim
     min = maxval
     do q = 1,ndim
        if(.not.sorted(q).and.min.gt.sortonr(q)) then
           min = sortonr(q)
           idx = q
        endif
     enddo
     idxarr(p) = idx
     sorted(idx) = .true.
  enddo
  deallocate(sortonr, sorted)
  return
end subroutine eig_order_abs_val


subroutine permute_eigensystem1(idxarr, ev, evr, ndim)
  implicit none
  !---------- passed variables ----------
  integer, intent(in)       :: idxarr(ndim)     ! index array which gives       proper order
  real*8, intent(inout)     :: ev(ndim)         ! array of eigenvalues
  complex*16, intent(inout) :: evr(ndim,ndim)   ! matrix of right eigenvectors  (column)
  integer :: ndim                               ! dimension of matrices
  !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
  !---------- local variables ------------------
  integer :: p
  complex*16, allocatable :: eval(:)
  complex*16, allocatable :: evec(:,:)
  allocate(eval(ndim), evec(ndim,ndim))
  !---------- permute the eigenvalues ----------
  do p = 1,ndim
     eval(p) = ev(idxarr(p))
  enddo
  ev = eval
  !---------- permute the right eigenvectors ----------
  do p = 1,ndim
     evec(:,p) = evr(:,idxarr(p))
  enddo
  evr = evec
  !---------- deallocate dynamic memory storage ----------
  deallocate(eval, evec)
  return
end subroutine permute_eigensystem1


!find max reciprocal lattice vectors
subroutine get_kmax(kmax,nat,nnlo)
    use param,only:fh_vec
    use com_mpi,only:nvector,vector_para,fvectors,findmaxk_mpi,vectors
    use dmf,only: qcomplex
    implicit none
    integer,intent(in)::nnlo,nat
    integer,intent(out)::kmax(3)

    integer :: ivector,itmp,n,ne,i,j,iks,ios,num
    integer,allocatable::keigen(:,:)
    real(8) :: rtmp
    complex(8) :: ztmp
    character :: stmp*10

    kmax=0
    do ivector=1,nvector
        if (vector_para) then
        open(fh_vec,file=fvectors(ivector,1), &
                &status='old',form='unformatted')
        else
            rewind(fh_vec)
        endif
        do i=1,nat
            read(fh_vec) rtmp; read(fh_vec) rtmp
        enddo
        num=0
        do iks=1,vectors(ivector,2) 
            read(fh_vec,iostat=ios) rtmp,rtmp,rtmp,stmp,n,ne
            if(ios/=0)exit
            if(allocated(keigen))then
                if(num<n)then
                    num=n+7
                    deallocate(keigen)
                    allocate(keigen(3,num))
                endif
            else
                num=n+7
                allocate(keigen(3,num))
            endif
            read(fh_vec) ((keigen(i,j),i=1,3),j=1,n)
            do i=n-nnlo,1,-1
                do j=1,3
                    kmax(j)=max(keigen(j,i),kmax(j))
                enddo
            enddo
            ! skip-read the remaining
            num=0
            do while(num/=ne)
                read(fh_vec)num,rtmp
                if(qcomplex)then
                    read(fh_vec)(ztmp, i=1,n)
                else
                    read(fh_vec)(rtmp, i=1,n)
                endif
            enddo
        enddo
        if (vector_para) then
            close(fh_vec)
        else
            rewind(fh_vec)
        endif
    enddo
    deallocate(keigen)
    call findmaxk_mpi(kmax)
    return

end subroutine get_kmax

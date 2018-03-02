!*****************************************************************************
! set bnd_ne: correlated band interval
!*****************************************************************************
subroutine gutz_set_ef_lda(nume,nkpt,projector)
    use com, only: emin,emax,ef,iso,ispin_pol,isospin,bnd_ne,bnd_ek
    use struct, only: nat
    use com_mpi,only: nvector,vector_para,fvectors,vectors, &

#ifdef MPI
            &dsum_all_mpi,isum_all_mpi, &
#endif

            &myrank,master
    implicit none
    integer nume,nkpt,projector
      
    integer k,ikp,n,nen,i,ii,ios,num,nemin,nemax
    integer ivector,itape,is
    real(8) ss,t,z,e1
    logical linef
    character*10 kname
    character*200   fname
    real(8),parameter :: rytoev=13.60569193_8
    integer,parameter :: iu=90

    allocate(bnd_ek(nume,nkpt,isospin))
    bnd_ek=3.d0 ! Rydberg
    allocate(bnd_ne(3,nkpt,isospin)); bnd_ne=0
      
    ! Use only one spin-channel
    do ivector=1,nvector
        if (vector_para) then
            do is=1,isospin    !------ over up/dn ---------------------!
                fname = fvectors(ivector, is+2)
                itape=61-is
                open(itape,file=fname,status='old',form='formatted')
            enddo
        endif
        do is=1,isospin
            itape = 61-is
            do i=1,nat
                read(itape,*); read(itape,*)
            enddo
            ios=0; k=0
            do while (ios == 0)
                read(itape,'(3e19.12,a10,2i6)',iostat=ios) ss,t,z,kname,n,nen
                if (ios /= 0) cycle
                k=k+1
                if (vector_para) then
                    ikp = vectors(ivector,3)+k
                else
                    ikp=k
                endif
                bnd_ne(1,ikp,is)=nen ! global
                do ii=1,nen
                    read(itape,*) num,e1
                    bnd_ek(ii,ikp,is)=e1
                enddo
            enddo
        end do
        do is=1,isospin    !------ over up/dn ---------------------!
            itape=61-is
            if (vector_para) then
                close(itape)
            else
                rewind(itape)
            endif
        end do
    end do
#ifdef MPI
    call dsum_all_mpi(bnd_ek,nume*nkpt*isospin)
    call isum_all_mpi(bnd_ne,3*nkpt*isospin)
#endif
    inquire(file='EFLDA.INP',exist=linef)
    if(linef)then
        call gutz_read_fermi(ef,iu)
        if (myrank.eq.master) &
                &write(6,'(" fermi level read in from EFLDA.INP")')
    else
        call gutz_fermi(iso,isospin)
        call gutz_wrt_fermi(ef,iu)
    endif
    if (myrank.eq.master) write(6,'(" fermi-lda=",f8.3," ryd,", &
            &f8.3," ev.")')ef,ef*rytoev
      
    emin=emin+ef; emax=emax+ef
    if (myrank.eq.master)then
        write(6,'(" absolute emin=",f8.3," ryd.",f8.3, &
                &" ev.")')emin,emin*rytoev
        write(6,'(" absolute emax=",f8.3," ryd.",f8.3, &
                &" ev.")')emax,emax*rytoev
    endif
     
    do is=1,isospin
        do k=1,nkpt
            nemin=1; nemax=0
            do ii=1,bnd_ne(1,k,is)
                if(bnd_ek(ii,k,is).lt.emin) nemin=nemin+1
                if(bnd_ek(ii,k,is).lt.emax) nemax=nemax+1
            enddo ! ii
            bnd_ne(2,k,is)=nemin; bnd_ne(3,k,is)=nemax
        enddo ! k
    enddo
    if(projector>0)then  ! keep same number of bands over k
        do is=1,isospin
            nemin=minval(bnd_ne(2,:,is)); bnd_ne(2,:,is)=nemin
            nemax=maxval(bnd_ne(3,:,is))
            do k=1,nkpt
                bnd_ne(3,k,is)=min(nemax,bnd_ne(1,k,is))
            enddo
        enddo
    else
        projector=-projector
    endif
    deallocate(bnd_ek)
    return
      
end subroutine gutz_set_ef_lda


subroutine set_gfile(deffn)
    use com,only: gfile
    implicit none
    character*(*)deffn

    integer idx

    idx = index(deffn,".def")
    if(deffn(idx-2:idx-1)=="dn")then
        gfile=2
    endif
    return
    
end subroutine set_gfile


subroutine gutz_fermi(iso,isospin)
    ! nelec is the number of electrons in this compound
    ! ef is the fermi-energy to be calculated
    use param,   only: nume,nkpt
    use com_mpi, only: myrank,master
    use com,only : efmod,ef,delta,bnd_ne,bnd_ek
    use bandm
    implicit none
    integer iso,isospin
      
    integer ik,ib,isp
    real(8) e1
      
    call init_bandm(nume)
     
    do isp=1,isospin
        do ik=1,nkpt
            do ib=1,bnd_ne(1,ik,isp)
                e1=bnd_ek(ib,ik,isp)
                if(e1.gt.ebmax(ib)) ebmax(ib)=e1
                if(e1.lt.ebmin(ib)) ebmin(ib)=e1
            enddo
        enddo
    enddo

    if(efmod.eq.'GAUSS'.or.efmod.eq.'TEMPS') then 
        ! just to get the energy window
        !         ef=ef-int(ef)
        !         if(ef.gt.0.1d0) ef=0.002d0
        ef=delta
        if (myrank.eq.master) &
                &write(6,'(" gauss-smearing with ef = ",f10.5,4h ry )') ef
        call gutz_fermig(iso,isospin)
    elseif(efmod.eq.'TETRA') then
        !     tetraeder method of bloechl (for ef.lt.100  corr is used!!)
        call gutz_fermi_tetra(iso,isospin)
        if (myrank.eq.master) &
                &write(6,'(" tetrahedron with ef = ",f10.5,4h ry )') ef
    else
        write(0,'(" efmod=",a5)')efmod
        stop 'not implimented yet!'
    end if
    return
      
end subroutine gutz_fermi
  
!*****************************************************************************
subroutine gutz_fermi_tetra(iso,isospin)
    ! nelec is the number of electrons in this compound
    ! ef is the fermi-energy to be calculated
    use com_mpi, only: myrank,master
    use param,only: nkpt,nume
    use com,only : ef,elecn,bnd_ne,bnd_ek
      
    implicit real*8 (a-h,o-z)
    parameter (nw=250000)
    integer          :: iw(nw)
    real(8),allocatable :: weight_(:),e_(:,:)
      
    nemax=maxval(bnd_ne(1,:,:))
    elecn=elecn-1.d-10
    !para begin
    ! ensure to mimick the proper vector file
    ! if running on multiple processors
    if (myrank.eq.master) write(6,*)'call eord...'
    allocate(e_(nemax*isospin,nkpt)); e_=0
    call eord(e_,bnd_ek,nemax,nkpt,isospin)
    if (myrank.eq.master) write(6,*)'call dos...'
    allocate(weight_(nemax*isospin*nkpt)); weight_=0
    call dos(nemax*isospin,nkpt,e_,weight_,elecn/2.d0*iso*isospin,d,ef,iw,nw)
    if (myrank.eq.master) write(6,*)'call eweigh...'
    elecn=elecn+1.d-10
    deallocate(e_,weight_)
    return
      
end subroutine gutz_fermi_tetra
  
!*****************************************************************************
subroutine gutz_fermig(iso,isospin)
    use com_mpi, only: myrank,master
    use param,only: nkpt,nume
    use kpts,only: weight=>mweight
    use com,only : ef,elecn,bnd_ek
      
    implicit real*8 (a-h,o-z)
    real*8 ztot,del1
    integer i,j,k,isp
    real(8),allocatable :: we(:,:,:)
      
    ztot=elecn
    del1=ef
    allocate(we(nume,nkpt,isospin)); we=0
    call fermi_gs(ztot,del1,iso,isospin,bnd_ek,we,nume,weight,nkpt, &
            &3.9d0,1.d-4,ef)
    deallocate(we)
    return
    154   format(/,9x,'bandenergy correction:',f10.6)
      
end subroutine gutz_fermig
  
!*****************************************************************************
subroutine fermi_gs(ztot,delta,iso,isospin,e,we,nbmax,wtk,nkpt,&
            &ecut,estep,ef)
    implicit none
    integer,intent(in)::nbmax,nkpt,iso,isospin
    real(8),intent(in)::ztot,delta,ecut,estep,e(nbmax,nkpt,isospin),wtk(nkpt)
    real(8),intent(out)::we(nbmax,nkpt,isospin),ef
      
    real(8) efermi,nel,nel0,dnnt,ef0,ef1,de,wt,fac,e1
    integer ib,ikp,isp,nnt,nc,rso,rspo
    real(8),parameter::tol=1.e-6_8, pi=3.141592653589793238_8
    integer,parameter::ncmax=10000
      
    rso=3-iso; rspo=3-max(iso,isospin)
    nnt=ztot/rso; dnnt=ztot/rso-nnt
    if(nnt==0)then
        ef0=e(1,1,1)
    else
        ef0=e(nnt,1,1)+dnnt*(e(nnt+1,1,1)-e(nnt,1,1))
    endif
    efermi=ef0
    do nc=1,ncmax
        we=0
        do isp=1,isospin
            do ikp=1,nkpt
                do ib=1,nbmax
                    if(e(ib,ikp,isp)>ecut)cycle
                    e1=e(ib,ikp,isp)
                    de=(e1-efermi)/delta
                    if(de<-3._8)then
                        wt=2._8
                    elseif(de<0._8)then
                        wt=2-erfc(-de)
                    elseif(de<3._8)then
                        wt=erfc(de)
                    else
                        wt=0._8
                    endif
                    wt=wt/2*rspo
                    we(ib,ikp,isp)=wt*wtk(ikp)
                enddo
            enddo
        enddo
        nel=sum(we)
        ! special treatment of the first iteration step
        if(nc==1)then
            efermi=efermi+estep
            nel0=nel
            cycle
        endif
        ! all other steps
        if(abs(ztot-nel)<tol)goto 100 ! located ef
        ef1=efermi
        fac=(nel0-nel)/(ztot-nel)
        if(abs(fac)>.1_8)then
            efermi=efermi+(ef0-efermi)/fac
            ef0=ef1; nel0=nel
        else
            ef0=ef1; nel0=nel
            if((ztot-nel)<0)then
                efermi=efermi-estep
            else
                efermi=efermi+estep
            endif
        endif
    enddo ! nc
      
    write(0,*) 'error in fermi_gs: fermilevel not converged!'; stop
    100   continue
    ef=efermi
    return
      
end subroutine fermi_gs
 

!*************************************************************************
! locate the index of identity operator
!*************************************************************************
function get_sym_ie(iord,iz,tau)
    integer,intent(in) :: iord, iz(3,3,iord)
    real(8),intent(in) :: tau(3,iord)
    integer :: get_sym_ie
      
    integer isym,i,j
      
sy: do isym=1,iord
        do i=1,3
            do j=1,3
                if(i.eq.j)then
                    if(iz(i,j,isym).ne.1) cycle sy
                else
                    if(iz(i,j,isym).ne.0) cycle sy
                endif
            enddo
        enddo
        do i=1,3
            if(abs(tau(i,isym))>1.d-16) cycle sy
        enddo
        get_sym_ie=isym
        return
    enddo sy
    stop ' error: failed to locate identity operation!'
      
end function get_sym_ie


!*****************************************************************************
! spherical harmonics->projector
!*****************************************************************************
subroutine gutz4_wrt(cfx,maxdim2,norbitals,nindo)
    use ghdf5_base
    use com,only: gfile
    implicit none
    integer maxdim2,norbitals,nindo(norbitals)
    complex(8) cfx(maxdim2,maxdim2,norbitals,norbitals)
      
    integer i1,i2,nindo_
    complex(8) c2n(maxdim2,maxdim2)
    character str*7
    
    if(gfile==2)return
    do i1=1,norbitals
        do i2=1,norbitals
            if(i1.eq.i2)cycle
            if(sum(abs(cfx(:,:,i1,i2))).lt.1.e-6_8)cycle
            write(*,'(" not implimented: cfx has non-zero cross terms &
                    &over norbitals!")'); stop
        enddo
        write(str,'(I7)')i1
        call gh5_create_group("/IMPURITY_"//trim(adjustl(str)),log_file_id)
    enddo
     
    do i1=1,norbitals
        write(str,'(I7)')i1
        nindo_=nindo(i1)
        c2n=0
        c2n(1:nindo_,1:nindo_)=transpose(conjg( &
                &cfx(1:nindo_,1:nindo_,i1,i1)))
        call gh5_write(c2n(1:nindo_,1:nindo_),nindo_,nindo_,"/IMPURITY_"//&
                &trim(adjustl(str))//&
                &"/T_SPHERICALHARMONICS_TO_LOCALBASIS1", &
                &log_file_id)
    enddo
    return
      
end subroutine gutz4_wrt
  
!*****************************************************************************
subroutine gutz5_wrt(nkp,wt,nelet,ne,kx,ky,kz,kname,efmod,delta)
    use ghdf5_base
    use com,only: gfile
    implicit none
    integer nkp
    real(8) wt(nkp),nelet,delta
    integer :: ne(3,nkp)
    real(8) :: kx(nkp),ky(nkp),kz(nkp)
    character*10::kname(nkp)
    character::efmod*5,str1*1
    common /gener/ br1(3,3), br2(3,3)
    real(8) :: br1, br2
      
    integer i,ismear
    character::file_name*512
      
    if(efmod=='GAUSS')then
        ismear=0
    elseif(efmod=='TETRA')then
        ismear=-5
    elseif(efmod=='TEMPS')then
        ismear=-1
    else
        write(0,'(" efmod=",a6)')efmod; stop ' error: unsupported efmod!'
    endif
    
    if(gfile==1)then
        call gh5_write(nkp,'/kptdim',log_file_id)
        call gh5_write(wt,nkp,'/kptwt',log_file_id)
        call gh5_write(ismear,'/ismear',log_file_id)
        call gh5_write(delta,'/delta',log_file_id)
        call gh5_write(nelet,'/nelectron',log_file_id)
        inquire(unit=14,name=file_name)
        file_name=file_name(index(file_name,'/',.true.)+1:)
        call gh5_write(file_name,512,'/kptfname',log_file_id)
        call gh5_write(kname,10,nkp,'/kptname',log_file_id)
        call gh5_write(kx,nkp,'/kptx',log_file_id)
        call gh5_write(ky,nkp,'/kpty',log_file_id)
        call gh5_write(kz,nkp,'/kptz',log_file_id)
        call gh5_write(br2,3,3,'/recip_prim_vec',log_file_id)
    endif
    write(str1,'(i1)')gfile
    call gh5_write(ne,3,nkp,'/NE_LIST_SPIN'//str1,log_file_id)
    return
      
end subroutine gutz5_wrt
  
!*****************************************************************************
subroutine gmpi_wrt()
    use com_mpi
    use com, only: gfile
    use ghdf5_base
    implicit none
      
    integer i
    character*7 str1
    
    if(nprocs==1.or.gfile==2)return
    write(str1,'(i7)')myrank
    call gh5_open_w('GMPI_'//trim(adjustl(str1))//'.h5', f_id)
    call gh5_write(myrank,'/myrank',f_id)
    call gh5_write(nprocs,'/nprocs',f_id)
    call gh5_write(master,'/master',f_id)
    call gh5_write(nvector,'/nvec',f_id)
    if(vector_para)then
        call gh5_write(vectors(1:nvector,1:3),nvector,3,'/KVEC',f_id)
    else
        stop " Only vector_para is supported!"
    endif
    call gh5_close(f_id)
    return
      
end subroutine gmpi_wrt
  
!*****************************************************************************
subroutine gutz_read_fermi(ef,iu)
    use com, only: gfile
    implicit none
    integer iu,i
    real(8) ef
      
    open(iu,file='EFLDA.INP',status='old')
    do i=1,gfile
        read(iu,*)ef
    enddo
    close(iu)
    return
      
end subroutine gutz_read_fermi
  
!*****************************************************************************
subroutine gutz_wrt_fermi(ef,iu)
    use com, only: gfile
    implicit none
    integer iu
    real(8) ef
    
    if(gfile==1)then
        open(iu,file='EFLDA.OUT',status='replace')
    else
        open(iu,file='EFLDA.OUT',status='old',access="append")
    endif
    write(iu,*)ef,' ! ef_lda'
    close(iu)
    return
      
end subroutine gutz_wrt_fermi


! Assume m>=n
subroutine zlowdin_orthogonalization(a,m,n)
    implicit none
    integer,intent(in)::m,n
    complex(8),intent(inout)::a(m,n)

    integer lwork,info
    complex(8) u(m,m),vt(n,n),work(8*m)
    real(8) s(n),rwork(5*n)
    complex(8),parameter::z1=cmplx(1.d0,0.d0),z0=cmplx(0.d0,0.d0)

    if(m<n)then
        stop ' error in zlowdin_orthogonalization: m<n!'
    endif
    lwork=8*m
    call zgesvd('a','a',m,n,a,m,s,u,m,vt,n,work,lwork,rwork,info)
    if(info/=0)then
        stop ' error in lowdin_orth_zmatrix: zgesvd failed!'
    endif
    call zgemm('n','n',m,n,n,z1,u(:,1:n),m,vt,n,z0,a,m)
    return

end subroutine zlowdin_orthogonalization


! Assume m>n
subroutine zappend_complementary_space(a,m,n)
    implicit none
    integer,intent(in)::m,n
    complex(8),intent(inout)::a(m,m)

    integer lwork,info
    real(8) w(m),rwork(3*m)
    complex(8),parameter::z1=cmplx(1.d0,0.d0),z0=cmplx(0.d0,0.d0)
    complex(8) aa(m,m),work(8*m)

    if(m<n)then
        stop ' error in zappend_complement_space: m<n!'
    elseif(m==n)then
        ! no complementary subspace
        return
    endif
    call zgemm('n','c',m,m,n,z1,a(:,1:n),m,a(:,1:n),m,z0,aa,m)
    lwork=8*m
    call zheev('v','l',m,aa,m,w,work,lwork,rwork,info)
    if(info/=0)then
        stop ' error in zappend_complement_space: zheev failed!'
    endif
    a(:,n+1:)=aa(:,1:m-n)
    return

end subroutine zappend_complementary_space


subroutine zget_hamilt_from_ev(v,w,m)
    implicit none
    integer,intent(in) :: m
    complex(8) v(m,m)
    real(8) w(m)

    integer i
    complex(8) vt(m,m),vte(m,m)
    complex(8),parameter::z1=cmplx(1.d0,0.d0),z0=cmplx(0.d0,0.d0)

    ! transform to <local orbitals | psi>
    vt=transpose(conjg(v)) 
    forall(i=1:m)vte(:,i)=vt(:,i)*w(i)
    call zgemm('n','c',m,m,m,z1,vte,m,vt,m,z0,v,m)
    return

end subroutine zget_hamilt_from_ev

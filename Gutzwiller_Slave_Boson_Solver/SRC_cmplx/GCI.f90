!******************************************************************************
! Copyright c 2013, The Ames Laboratory, Iowa State University, and Rutgers
! University*.  All rights reserved.
!
! This software was authored by Yongxin Yao, Nicola Lanata*, Gabriel Kotliar*,
! Cai-Zhuang Wang, and Kai-Ming Ho, at The Ames Laboratory and
! Rutgers University and was supported by the U.S.
! Department of Energy (DOE), Office of Science,
! Basic Energy Sciences, Materials Science and Engineering Division.
! The Ames Laboratory is operated by Iowa State University for DOE
! under U.S. Government contract DE-AC02-07CH11358.
! The U.S. Government has the rights to use, reproduce, and
! distribute this software.
! NEITHER THE GOVERNMENT, THE AMES LABORATORY, IOWA STATE UNIVERSITY,
! NOR RUTGERS UNIVERSITY MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
! If software is modified to produce derivative works,
! such modified software should be clearly marked,
! so as to not confuse it with the version available from
! The Ames Laboratory and Rutgers University.
!
! Additionally, redistribution and use in source and binary forms,
! with or without modification,
! are permitted provided that the following conditions are met:
!
!     Redistribution of source code must retain the above copyright notice,
!     this list of conditions, and the following disclaimer.
!
!     Redistribution in binary form must reproduce the above copyright notice,
!     this list of conditions, and the following disclaimer
!     in the documentation and/or other materials provided with distribution.
!
!     Neither the name of The Ames Laboratory, Iowa State University,
!     Rutgers University, the U.S. Government, nor the names of
!     its contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE AMES LABORATORY, IOWA STATE UNIVERSITY,
! RUTGERS UNIVERSITY, AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY
! AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE GOVERNMENT, THE AMES LABORATORY,
! IOWA STATE UNIVERSITY, RUTGERS UNIVERSITY, OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
! OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!******************************************************************************

module gci
    ! optimized for CISD.
    use gprec
    use gconstant
    use gutil
    implicit none

    type icscn_matrix
        integer::ncol=0
        integer,pointer::p(:,:),i(:)=>null(),j(:),phase(:)
    end type icscn_matrix

    type gci_mem
        integer::norb,norb2,nstates=0,ci_order=2
        type(icscn_matrix)::h1,h2
        integer,pointer::bs(:),ibs(:),idx(:)
    end type gci_mem

    integer,allocatable::mem_binom(:,:)
    integer::gci_nstates
    complex(q),allocatable::gh_ci(:)
    integer,parameter,dimension(7)::nnz_pqrs=(/132,9504,134316,857856,3549300,&
            &11202912,29407644/)
    integer,parameter,dimension(7)::nnz_pq=(/24,564,4488,19848,62760,160284,&
            &353304/)
    type(gci_mem),allocatable::gci_mem_list(:)

    private
    public::gci_mem,gci_init,set_full_fock_states_ci,solve_hembed_ci, &
            &gci_nstates,gh_ci,set_h1_ci,calc_dm_ci,gci_mem_list, &
            &init_gci_mem_list,solve_hembed_ci_drive

    contains


    subroutine gci_init(n)
    integer,intent(in)::n

    allocate(mem_binom(n,0:n))
    call set_binoms(n,mem_binom)
    return

    end subroutine gci_init


    subroutine init_gci_mem_list(norb_list,ntyp,io)
    integer,intent(in)::io,ntyp,norb_list(ntyp)
    
    integer i

    allocate(gci_mem_list(ntyp))
    do i=1,ntyp
        call set_full_fock_states_ci(gci_mem_list(i),norb_list(i))
        call set_h1_ci(gci_mem_list(i))
        if(io>0)then
            write(io,'(" ityp_local_imp = ",i2," dim-phi = ",i0)')i,&
                    &gci_mem_list(i)%nstates
        endif
    enddo
    return

    end subroutine init_gci_mem_list


    subroutine set_full_fock_states_ci(dmem,norb,ci_order)
    type(gci_mem),intent(inout)::dmem
    integer,intent(in)::norb
    integer,optional,intent(in)::ci_order

    integer i

    if(.not.associated(dmem%bs))then
        dmem%norb=norb; dmem%norb2=norb*2
        call set_fock_state_indices(norb,mem_binom,dmem%idx,dmem%bs,dmem%ibs)
    endif

    if(present(ci_order))then
        if(dmem%ci_order/=ci_order)then
            if(associated(dmem%h1%i)) call reset_icscn_matrix(dmem%h1)
            if(associated(dmem%h2%i)) call reset_icscn_matrix(dmem%h2)
            if(allocated(gh_ci)) deallocate(gh_ci)
            dmem%ci_order=ci_order
            dmem%nstates=0
        endif
    endif
    if(dmem%nstates==0)then
        do i=0,dmem%ci_order
            dmem%nstates=dmem%nstates+(dmem%idx(i+1)-dmem%idx(i))**2
        enddo
    endif
    return

    end subroutine set_full_fock_states_ci


    subroutine solve_hembed_ci_drive(dmem,h1e,v2e,dm,etot)
    type(gci_mem),intent(inout)::dmem
    complex(q),intent(in)::h1e(dmem%norb2,dmem%norb2), &
            &v2e(dmem%norb2,dmem%norb2,dmem%norb2,dmem%norb2)
    complex(q),intent(inout)::dm(dmem%norb2,dmem%norb2)
    real(q),intent(out)::etot

    integer method
    complex(q) eri(dmem%norb2,dmem%norb2,dmem%norb2,dmem%norb2), &
            &v(dmem%nstates),ao2mo(dmem%norb2,dmem%norb2)
    real(q) w(dmem%norb2)

    ! dm_ij = <c^\dagger_i c_j>
    ao2mo=-transpose(dm)
    call hermev('v','l',ao2mo,w,dmem%norb2)
    call get_eri(v2e,h1e,eri,dmem%norb2)
    call a_trans(eri,dmem%norb2,ao2mo)
    if(dmem%norb2<=8)then
        method=0
    else
        method=1
    endif
    call solve_hembed_ci(dmem,eri,v,w(1:2),method=method)
    etot=w(1)
    call calc_dm_ci(dmem,v,dm)
    ao2mo=conjg(ao2mo)
    call uhau(dm,ao2mo,dmem%norb2,dmem%norb2,trul='n',trur='c')
    return

    end subroutine solve_hembed_ci_drive


    subroutine solve_hembed_ci(dmem,eri,v,w,method)
    type(gci_mem),intent(inout)::dmem
    complex(q),intent(in)::eri(dmem%norb2,dmem%norb2,dmem%norb2,dmem%norb2)
    complex(q),intent(inout)::v(dmem%nstates)
    integer,intent(in),optional::method
    real(q),intent(inout)::w(2)
    integer::i_method=1
    external::av_gci

    if(present(method))then
        i_method=method
    endif

    call set_hembed_ci(dmem,eri)
    gci_nstates=dmem%nstates
    if(i_method==0)then
        call lapack_diag_ci(v,dmem%nstates,w)
    else
        call primme_diag(v,dmem%nstates,w,av_gci)
    endif
    return

    end subroutine solve_hembed_ci


    !*************************************************************************
    subroutine lapack_diag_ci(v,n,w)
    integer,intent(in) :: n
    complex(q),intent(out) :: v(n)

    real(q),intent(out) :: w(2)

    integer i
    real(q) w_(n)
    complex(q) z(n,n),v1(n)

    call hermev('v','u',gh_ci,w_,z,n)
    v=z(:,1)
    w=w_(:2)
    return

    end subroutine lapack_diag_ci


    subroutine set_hembed_ci(dmem,eri)
    type(gci_mem),intent(inout)::dmem
    complex(q),intent(in)::eri(dmem%norb2,dmem%norb2,dmem%norb2,dmem%norb2)

    integer i,j1,j2,ij,p,q,r,s,itmp(0:4),i_sign(4),n_hole(4),mode,nnz
    integer ibase(0:dmem%ci_order),ibs1,ibs2,mask1,nbs,istate,jstate

    if(associated(dmem%h2%i))then
        call set_hembed_ci1(dmem,eri)
        return
    endif

    if(dmem%h2%ncol==0)then
        dmem%h2%ncol=dmem%nstates
        allocate(dmem%h2%j(dmem%nstates+1))
        if(dmem%norb/2>size(nnz_pqrs).or.dmem%ci_order/=2)then
            mode=0
        else
            mode=1; nnz=nnz_pqrs(dmem%norb/2)
        endif
    else
        mode=1
        nnz=dmem%h2%j(dmem%h2%ncol+1)-1
    endif

    if(mode==1)then
        allocate(dmem%h2%p(4,nnz),dmem%h2%i(nnz),dmem%h2%phase(nnz))
    endif
    nbs=ishft(1,dmem%norb); mask1=nbs-1
    ibase=0
    do i=0,dmem%ci_order-1
        ibase(i+1)=ibase(i)+(dmem%idx(i+1)-dmem%idx(i))**2
    enddo
    ! Only save upper triangular part.
    if(.not.allocated(gh_ci)) &
            &allocate(gh_ci(dmem%nstates*(dmem%nstates+1)/2))
    gh_ci=0
    jstate=0; nnz=1
    do i=0,dmem%ci_order
    ! \psi_0 occupied part
    do j1=dmem%idx(dmem%norb-i+1)-1,dmem%idx(dmem%norb-i),-1
    ! \psi_0 unoccupied part
    do j2=dmem%idx(i),dmem%idx(i+1)-1
        jstate=jstate+1
        dmem%h2%j(jstate)=nnz
        itmp(0)=ior(dmem%bs(j1),ishft(dmem%bs(j2),dmem%norb))
        do s=1,dmem%norb2
            itmp(1)=itmp(0); i_sign(1)=1
            ! s^-
            call act_state(itmp(1),s-1,.false.,i_sign(1))
            if(i_sign(1)==0)cycle
            if(s<=dmem%norb)then
                n_hole(1)=i+1
            else
                n_hole(1)=i
            endif
            do r=1,dmem%norb2
                itmp(2)=itmp(1); i_sign(2)=i_sign(1)
                ! r^\dagger
                call act_state(itmp(2),r-1,.true.,i_sign(2))
                if(i_sign(2)==0)cycle
                if(r<=dmem%norb)then
                    n_hole(2)=n_hole(1)-1
                else
                    n_hole(2)=n_hole(1)
                endif
                do q=1,dmem%norb2
                    itmp(3)=itmp(2); i_sign(3)=i_sign(2)
                    ! q^-
                    call act_state(itmp(3),q-1,.false.,i_sign(3))
                    if(i_sign(3)==0)cycle
                    if(q<=dmem%norb)then
                        n_hole(3)=n_hole(2)+1
                    else
                        n_hole(3)=n_hole(2)
                    endif
                    do p=1,dmem%norb2
                        itmp(4)=itmp(3); i_sign(4)=i_sign(3)
                        ! p^\dagger
                        call act_state(itmp(4),p-1,.true.,i_sign(4))
                        if(i_sign(4)==0)cycle
                        if(p<=dmem%norb)then
                            n_hole(4)=n_hole(3)-1
                        else
                            n_hole(4)=n_hole(3)
                        endif
                        ! upper triangular part
                        if(n_hole(4)>i)cycle
                        ibs1=dmem%ibs(iand(itmp(4),mask1))
                        if(ibs1<j1)cycle
                        ibs2=dmem%ibs(ishft(itmp(4),-dmem%norb))
                        if(ibs1==j1.and.ibs2>j2)cycle
                        istate=ibase(n_hole(4))+ &
                                &(dmem%idx(dmem%norb-n_hole(4)+1)-1-ibs1)* &
                                &(dmem%idx(n_hole(4)+1)-dmem%idx(n_hole(4)))+&
                                &ibs2-dmem%idx(n_hole(4))+1
                        ij=istate+(jstate-1)*jstate/2
                        gh_ci(ij)=gh_ci(ij)+eri(p,q,r,s)*i_sign(4)
                        if(mode==1)then
                            dmem%h2%p(:,nnz)=(/p,q,r,s/)
                            dmem%h2%i(nnz)=istate
                            dmem%h2%phase(nnz)=i_sign(4)
                        endif
                        nnz=nnz+1
                    enddo
                enddo
            enddo
        enddo
    enddo; enddo; enddo
    dmem%h2%j(dmem%h2%ncol+1)=nnz
    return

    end subroutine set_hembed_ci


    subroutine set_hembed_ci1(dmem,eri)
    type(gci_mem),intent(inout)::dmem
    complex(q),intent(in)::eri(dmem%norb2,dmem%norb2,dmem%norb2,dmem%norb2)
    
    integer i,istate,jstate,ij

    gh_ci=0
    do jstate=1,dmem%h2%ncol
        do i=dmem%h2%j(jstate),dmem%h2%j(jstate+1)-1
            istate=dmem%h2%i(i)
            ij=istate+(jstate-1)*jstate/2
            gh_ci(ij)=gh_ci(ij)+&
                    &eri(dmem%h2%p(1,i),dmem%h2%p(2,i),dmem%h2%p(3,i),&
                    &dmem%h2%p(4,i))*dmem%h2%phase(i)
        enddo
    enddo
    return

    end subroutine set_hembed_ci1


    subroutine set_h1_ci(dmem)
    type(gci_mem),intent(inout)::dmem

    integer i,j1,j2,p,q,itmp(0:2),i_sign(2),n_hole(2),mode,nnz
    integer ibase(0:dmem%ci_order),ibs1,ibs2,mask1,nbs,istate,jstate

    do 
    if(dmem%h1%ncol==0)then
        dmem%h1%ncol=dmem%nstates
        allocate(dmem%h1%j(dmem%nstates+1))
        if(dmem%norb/2>size(nnz_pq).or.dmem%ci_order/=2)then
            mode=0
        else
            mode=1; nnz=nnz_pq(dmem%norb/2)
        endif
    else
        mode=1
        nnz=dmem%h1%j(dmem%h1%ncol+1)-1
    endif

    if(mode==1)then
        allocate(dmem%h1%p(2,nnz),dmem%h1%i(nnz),dmem%h1%phase(nnz))
    endif
    nbs=ishft(1,dmem%norb); mask1=nbs-1
    ibase=0
    do i=0,dmem%ci_order-1
        ibase(i+1)=ibase(i)+(dmem%idx(i+1)-dmem%idx(i))**2
    enddo
    ! Only save upper triangular part.
    jstate=0; nnz=1
    do i=0,dmem%ci_order
    ! \psi_0 occupied part
    do j1=dmem%idx(dmem%norb-i+1)-1,dmem%idx(dmem%norb-i),-1
    ! \psi_0 unoccupied part
    do j2=dmem%idx(i),dmem%idx(i+1)-1
        jstate=jstate+1
        dmem%h1%j(jstate)=nnz
        itmp(0)=ior(dmem%bs(j1),ishft(dmem%bs(j2),dmem%norb))
        do q=1,dmem%norb2
            itmp(1)=itmp(0); i_sign(1)=1
            ! q^-
            call act_state(itmp(1),q-1,.false.,i_sign(1))
            if(i_sign(1)==0)cycle
            if(q<=dmem%norb)then
                n_hole(1)=i+1
            else
                n_hole(1)=i
            endif
            do p=1,dmem%norb2
                itmp(2)=itmp(1); i_sign(2)=i_sign(1)
                ! p^\dagger
                call act_state(itmp(2),p-1,.true.,i_sign(2))
                if(i_sign(2)==0)cycle
                if(p<=dmem%norb)then
                    n_hole(2)=n_hole(1)-1
                else
                    n_hole(2)=n_hole(1)
                endif
                ! upper triangular part
                if(n_hole(2)>i)cycle
                ibs1=dmem%ibs(iand(itmp(2),mask1))
                if(ibs1<j1)cycle
                ibs2=dmem%ibs(ishft(itmp(2),-dmem%norb))
                if(ibs1==j1.and.ibs2>j2)cycle
                istate=ibase(n_hole(2))+ &
                        &(dmem%idx(dmem%norb-n_hole(2)+1)-1-ibs1)* &
                        &(dmem%idx(n_hole(2)+1)-dmem%idx(n_hole(2)))+&
                        &ibs2-dmem%idx(n_hole(2))+1
                if(mode==1)then
                    dmem%h1%p(:,nnz)=(/p,q/)
                    dmem%h1%i(nnz)=istate
                    dmem%h1%phase(nnz)=i_sign(2)
                endif
                nnz=nnz+1
            enddo
        enddo
    enddo; enddo; enddo
    dmem%h1%j(dmem%h1%ncol+1)=nnz
    if(mode==1)return
    enddo

    end subroutine set_h1_ci


    subroutine calc_dm_ci(dmem,v,dm)
    type(gci_mem),intent(in)::dmem
    complex(q),intent(in)::v(dmem%nstates) 
    complex(q),intent(out)::dm(dmem%norb2,dmem%norb2)

    complex(q) vp(dmem%nstates,dmem%norb2*(dmem%norb2+1)/2), &
            &dm1(dmem%norb2*(dmem%norb2+1)/2)
    integer i,istate,jstate,p,q,pq

    vp=0
    do jstate=1,dmem%nstates
        if(abs(v(jstate))<1.d-12)cycle
        do i=dmem%h1%j(jstate),dmem%h1%j(jstate+1)-1
            istate=dmem%h1%i(i)
            p=dmem%h1%p(1,i); q=dmem%h1%p(2,i)
            if(p<=q)then
                pq=p+(q-1)*q/2
                vp(istate,pq)=vp(istate,pq)+v(jstate)*dmem%h1%phase(i)
            else
                if(abs(v(istate))>1.d-12)then
                    pq=q+(p-1)*p/2
                    vp(jstate,pq)=vp(jstate,pq)+v(istate)*dmem%h1%phase(i)
                endif
            endif
        enddo
    enddo

    call zgemm('c','n',1,dmem%norb2*(dmem%norb2+1)/2,dmem%nstates,z1, &
            &v(1),dmem%nstates,vp(1,1),dmem%nstates,z0,dm1(1),1)
    do p=1,dmem%norb2; do q=p,dmem%norb2
        pq=p+(q-1)*q/2
        dm(p,q)=dm1(pq)
        dm(q,p)=conjg(dm1(pq))
    enddo; enddo
    return

    end subroutine calc_dm_ci


    subroutine reset_icscn_matrix(h)
    type(icscn_matrix),intent(inout)::h

    h%ncol=0
    deallocate(h%p,h%j,h%i,h%phase)
    nullify(h%p,h%j,h%i,h%phase)
    return

    end subroutine reset_icscn_matrix


end module gci


subroutine av_gci(v1,v2,k,primme)
    use gprec
    use gci
    use gconstant
    implicit none
    integer,intent(in)::k
    integer(8),intent(in)::primme
    complex(q),intent(in)::v1(*)
    complex(q),intent(out)::v2(*)

    integer n,i

    n=gci_nstates
    do i=1,k
        call zhpmv('u',n,z1,gh_ci,v1(n*(i-1)+1),1,z0,v2(n*(i-1)+1),1)
    enddo
    return

end subroutine av_gci

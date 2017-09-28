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

module gdiis
    use gprec
    use gutil
    use gconstant
    implicit none
      
    type diis_mem
        integer::head=0,space=8,min_space=1,numx=0
        integer::dimx
        logical::x_err_provided=.false.
        complex(q),pointer::x_err(:,:),x_vec(:,:),x_prev(:)=>null(), &
                &h(:,:),hinv(:)
    end type diis_mem

    contains


    subroutine gdiis_solve(f,n,x,fvec)
    integer,intent(in) :: n
    complex(q),intent(inout) :: x(n),fvec(n)
    external :: f

    type(diis_mem) :: dmem
    integer iflag

    call diis_vector_init(dmem,n)
    do 
        call f(n,x,fvec,iflag)
        if(iflag==-1)return
        call gdiis_update(dmem,x,fvec)
    enddo
    return

    end subroutine gdiis_solve

      
    subroutine gdiis_update(dmem,x,xerr)
    type(diis_mem),intent(inout)::dmem
    complex(q),intent(inout)::x(dmem%dimx)
    complex(q),intent(in),optional::xerr(dmem%dimx)

    integer i
    character*75 stmp
    complex(q),pointer::hinv(:,:)
    complex(q),external::zdotc

    ! For unknown ifort reason, it fixes.
    write(stmp,*)x(1)

    if(present(xerr))then
        if(dmem%head==dmem%space)then
            dmem%head=1
        else
            dmem%head=dmem%head+1
        endif
        dmem%x_err(:,dmem%head)=xerr(1:dmem%dimx)
        dmem%x_err_provided=.true.
    endif
      
    if(dmem%x_err_provided)then
        dmem%x_vec(:,dmem%head)=x(1:dmem%dimx)
        if(dmem%numx<dmem%space)then
            dmem%numx=dmem%numx+1
        endif
    elseif(.not.associated(dmem%x_prev))then
        allocate(dmem%x_prev(dmem%dimx))
        dmem%x_prev=x(1:dmem%dimx)
        return
    else
        if(dmem%head==dmem%space)then
            dmem%head=1
        else
            dmem%head=dmem%head+1
        endif
        dmem%x_vec(:,dmem%head)=x(1:dmem%dimx)
        dmem%x_err(:,dmem%head)=x(1:dmem%dimx)-dmem%x_prev
        if(dmem%numx<dmem%space)then
            dmem%numx=dmem%numx+1
        endif
    endif
      
    ! update h_{ij} = <v_err_i | v_err_j>
    do i=1,dmem%numx
        dmem%h(i,dmem%head)=zdotc(dmem%dimx,dmem%x_err(1,i),1, &
                &dmem%x_err(1,dmem%head),1)
        dmem%h(dmem%head,i)=conjg(dmem%h(i,dmem%head))
    enddo
      
    if(dmem%numx<dmem%min_space)return
      
    ! interpolate the solution
    ! ( 0   1  ...  1      1  ) ( lambda )   (1)
    ! ( 1  h1,1...h1,n-1 h1,n ) ( c1     )   (0)
    ! (   ................    ) ( .....  ) = (.)
    ! ( 1  hn,1...hn,n-1 hn,n ) ( cn     )   (0)

    hinv(1:dmem%numx+1,1:dmem%numx+1)=>dmem%hinv(1:(1+dmem%numx)**2)
    call atofa(dmem%h(0:dmem%numx,0:dmem%numx), &
            &hinv,dmem%numx+1,-1,d1,.true.)

    ! hinv(2:dmem%numx,1) is the solution.
    call zgemm('n','n',dmem%dimx,1,dmem%numx,z1,dmem%x_vec(1,1),dmem%dimx,&
            &hinv(2,1),dmem%numx,z0,x(1),dmem%dimx)
    nullify(hinv)
    return
    
    end subroutine gdiis_update
      
      
    subroutine diis_vector_init(dmem,dimx)
    integer,intent(in)::dimx
    type(diis_mem),intent(out)::dmem
      
    dmem%dimx=dimx
    allocate(dmem%x_err(dimx,dmem%space),dmem%x_vec(dimx,dmem%space), &
            &dmem%h(0:dmem%space,0:dmem%space), &
            &dmem%hinv((1+dmem%space)**2))
    dmem%h=0; dmem%h(0,1:)=1._q; dmem%h(1:,0)=1._q
    return
      
    end subroutine diis_vector_init
      
      
end module gdiis

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

module dcstd
    use gprec
    use gconstant
    use warehouse
    use ghdf5_base
    implicit none

    type dc_std
        integer::mode=1
        real(q),allocatable::e(:),vpot(:),u_avg(:),j_avg(:),nelf(:)
    end type dc_std

    type(dc_std)::dc


    contains


    subroutine init_dc_std(io)
    integer,intent(in)::io

    allocate(dc%vpot(wh%num_imp),dc%u_avg(wh%num_imp),dc%j_avg(wh%num_imp), &
            &dc%nelf(wh%num_imp),dc%e(wh%num_imp))
    dc%u_avg=0; dc%j_avg=0; dc%nelf=0
    call set_nelf_list(io)
    return

    end subroutine init_dc_std


    subroutine set_nelf_list(io)
    integer,intent(in)::io

    logical lexist

    call gh5_open_r('GPARAM.h5',f_id)
    call gh5_read(dc%mode,'/dc_mode',f_id)
    if(dc%mode>0)then
        call gh5_read(dc%u_avg,wh%num_imp,'/dc_u_avg',f_id)
        call gh5_read(dc%j_avg,wh%num_imp,'/dc_j_avg',f_id)
    endif
    if(dc%mode>1)then
        call h5lexists_f(f_id,'/dc_nelf_list',lexist,gh5_err)
        if(lexist)then
            call gh5_read(dc%nelf,wh%num_imp,'/dc_nelf_list',f_id)
        else
            dc%nelf=wh%co(:)%net
        endif
    endif
    call gh5_close(f_id)

    if(io>0)then
        select case(dc%mode)
        case(0)
            write(io,'(" no double counting.")')
        case(1)
            write(io,'(" standard fully localized limit (FLL) dc.")')
        case(2)
            write(io,'(" fixed double counting.")')
        case(12)
            write(io,'(" standard FLL dc with dc only updated at the &
                    &electron density cycles.")')
        case default
            stop ' error in GDC_NELF.INP: dc%mode not defined!'
        end select

        if(dc%mode>1)then
            write(io,'(" input nelf:")')
            write(io,'(8x,5(f8.3,2x))')dc%nelf
        endif

        if(dc%mode>0)then
            write(io,'(" average hubbard u list:")')
            write(io,'(8x,5(f8.3,2x))')dc%u_avg
            write(io,'(" average hund j list:")')
            write(io,'(8x,5(f8.3,2x))')dc%j_avg
        endif
    endif
    return

    end subroutine set_nelf_list


    subroutine update_nelf_list(io)
    integer,intent(in)::io

    integer idx(1)
    real(q) ndiff(wh%num_imp)
    character(20) fmt


    ndiff=wh%co(:)%net-dc%nelf
    idx=maxloc(abs(ndiff))

    if(io>0)then
        write(io,'(" max nelf diff=",f14.8)')ndiff(idx(1))
    endif

    if(io>0)then
        call gh5_open_w('GDC_NELF_OUT.h5',f_id)
        call gh5_write(dc%nelf,wh%num_imp,'/dc_nelf_list_inp',f_id)
        dc%nelf=wh%co(:)%net
        call gh5_write(dc%nelf,wh%num_imp,'/dc_nelf_list_out',f_id)
        call gh5_close(f_id)
    endif
    return

    end subroutine update_nelf_list


    subroutine calc_vdc_list()
    integer i,j

    if(dc%mode==1)then
        dc%nelf=wh%co(:)%net
    endif
    dc%vpot=dc%u_avg*(dc%nelf-.5_q)-dc%j_avg/2*(dc%nelf-1)

    ! add to co%la2
    wh%la2=0
    do i=1,wh%num_imp
        do j=1,wh%na2_imp(i)
            wh%co(i)%la2(j,j)=-dc%vpot(i)
        enddo
    enddo
    return

    end subroutine calc_vdc_list


    subroutine add_vdc_to_la1_list()

    call calc_vdc_list()
    wh%la1=wh%la1+wh%la2
    return

    end subroutine add_vdc_to_la1_list


    subroutine calc_edc_list()

    if(dc%mode==1)then
        dc%nelf=wh%co(:)%net
    endif
    dc%e=dc%u_avg*dc%nelf*(dc%nelf-1)/2-dc%j_avg*dc%nelf/2*(dc%nelf/2-1)+ &
            &+dc%vpot*(wh%co(:)%net-dc%nelf)
    return
    
    end subroutine calc_edc_list


    subroutine output_energies_dc(io)
    integer,intent(in)::io

    integer i

    write(io,'(" impurity-wise interaction dc energy:")')
    write(io,'(4x,5f14.7)')(dc%e(i), i=1,wh%num_imp)
    write(io,'(" total U-dc energy = ",f0.7)')sum(dc%e)
    return

    end subroutine output_energies_dc



end module dcstd

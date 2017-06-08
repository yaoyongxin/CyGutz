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

program cygutz

    use ghdf5_base, only: gh5_init, gh5_end
#ifdef mpi_mode
    use mpi
#endif
    use gparam
    use bandstru
    use warehouse
    use dcstd
    use gmpi
    use gkernel
    use psave
    implicit none
    integer ierr
    real ta1,ta2,tb1,tb2
    integer tib1,tib2,tirate
    external::g_fcn_rl
      
    call cpu_time(ta1)
    call system_clock(tib1,tirate)
    tb1=real(tib1,4)/real(tirate,4)
#ifdef mpi_mode
    call mpi_init(ierr)
#endif
    call gh5_init()
    call init_gmpi()
    if(gp%io>0)then
        open(gp%io,file='GUTZ.LOG',status='replace')
#ifdef with_compile_date
        write(gp%io,'(" cygutz (atomic rydberg units) vmpi_5.0 built on", &
                &a12)')__DATE__
#else
        write(gp%io,'(" cygutz (atomic rydberg units) vmpi_5.0")')
#endif
    endif

    call set_gparam(gp%io)
    call init_warehouse(gp%io)
    call set_bnd_info(gp%io)
    call set_gmpi() !< kpt%diml was set here.
    call alloc_bnd()
    call read_bare_hamiltonian()
    call rotate_bare_hamiltonian()

    ! correlated block energy window.
    call calc_corr_ebwidth(gp%io)

    ! check the bare band dispersions.
    call calc_band_all(gp%io)
    call gutz_fermi(gp%io)
    call calc_nks()
    call calc_nks_pp(gp%io)

    ! single out the local one-body part.
    call rm_h1e_from_bare_hamiltonian()

    ! double counting
    call init_dc_std(gp%io)

    ! run the kernel
    call init_gkernel(gp%io)
    call g_newton_solver(gp%io,g_fcn_rl)

    call postsave()

    ! update electron density
    if(g_updaterho>0)then
        call calc_rnrl()
        call map_wh_bnd_matrix(wh%nc_phy,bnd%nc_phy,.false.)
        call calc_kswt(gp%io)
        call gh5_wrt_kswt()
    endif

    if(gp%io>0)then
        call output_energies(gp%io)
        call gh5_wrt_wh_rl('WH_RL_OUT.h5')
        call update_nelf_list(gp%io)
    endif

    call cpu_time(ta2)
    call system_clock(tib2,tirate)
    tb2=real(tib2,4)/real(tirate,4)
    if(gp%io>0)then
        call out_time_use('total',ta2-ta1,tb2-tb1,gp%io)
        close(gp%io)
    endif
    call gh5_end()
#ifdef mpi_mode
    call mpi_finalize(ierr)
#endif
      
end program cygutz

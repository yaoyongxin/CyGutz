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
!
!
!
      MODULE GINFO
      USE gprec
      IMPLICIT NONE
!
      TYPE GL_INFO
        INTEGER IO,IU
        INTEGER ITER1,ITER2,ITER_RHO,ITER_GF,ITER_LA1,ITER_LA2
        INTEGER NASOMAX,NA2MAX,NAT_TOT,NEV_F,NCV_F,NCV_P
        INTEGER LDC,LSCF,LGREEN,LGPRJ,LDIAPJ,LHUB,LUNIT,LEBANDCOMP, &
                &LEL0,LBSCODE,LEFERMI
        INTEGER LPHISYM 
        ! switch for symmetrization of \phi. 0: off; 
        !                                   -1: generate; 
        !                                    1: use
        ! LBSCODE=0: eV and real Harmonics (VASP); 
        ! 1: Ryd and complex Harmonics (Wien2K)
        INTEGER LENTANGLES,LSOLVER,LEIGV,LPJC,LBNDU,LB2N
        INTEGER LCHKLOC,LV2AO,LHARMONICS,LNEWTON
        ! LHARMONICS = 0: Complex Harmonics
        !              1: Single real cubic Harmonics
        INTEGER LMODEL,LENSEMBLE,LPLTGF,LIADBOCC,LCURRENT
        INTEGER NMAX_ITER,NMAX_INNER,LMCFLY,LXGUESS
        INTEGER LFZ1SET,NB_RESET
        REAL(gq) RTOL,RCUT_HS
        REAL(gq) MAXERROR,DCV_ERROR
        INTEGER RMODE
! CMR3 begin
        INTEGER LFUNR
! CMR3 end
        REAL(gq) DCMIX_A,RHO_CUT,RMIX_A
        INTEGER LH5WRT,LA1ADJUST
      END TYPE GL_INFO
!
      TYPE (GL_INFO) ,SAVE :: GL
      INTEGER GL_NI
      !$OMP THREADPRIVATE(GL_NI)
!
!
      END MODULE GINFO

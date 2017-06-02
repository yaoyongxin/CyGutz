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
      PROGRAM CYGUTZ
      USE GUTZ
      USE GHDF5_BASE, ONLY: gh5_init, gh5_end
      USE MPI
      IMPLICIT NONE
      INTEGER IERR
      REAL :: TA1,TA2,TB1,TB2; INTEGER TIB1,TIB2,TIRATE
!
      CALL MPI_INIT(IERR)
      CALL CPU_TIME(TA1); CALL SYSTEM_CLOCK(TIB1,TIRATE); TB1=REAL(TIB1,4)/REAL(TIRATE,4)
      OPEN(UNIT=6,FILE='GUTZ.LOG',STATUS='REPLACE')
      CALL INI_GMPI_(90,6)
      IF(GP%MYRANK.EQ.GP%MASTER) CALL gh5_init("glog.h5")
      CALL GUTZ1_INI_(90,6)
      CALL GUTZ2_SET_NISO_()
      CALL GUTZ3_LOCATE_SYM_IE_()
      CALL GUTZ4_SET_C2N_UH_()
      CALL GUTZ5_INI_CORR_UVEK_()
      IF(KPT%FILE_NAME/='')THEN
        OPEN(UNIT=IU_KGEN,FILE=TRIM(ADJUSTL(KPT%FILE_NAME)),STATUS='OLD')
      ENDIF
      CALL GUTZ6_SOLVE()
      IF(GL%LSCF==2)THEN
        GOTO 100
      ENDIF
      IF(GL%LGREEN==0)THEN
        CALL CALC_KSWT(0)
        IF(GL%LBSCODE==1)THEN
          CALL WRT_KSWT()
        ELSEIF(GL%LBSCODE==0)THEN
          CALL DIAG_KSWT()
          CALL WRT0_KSWT(-2)
          IF(GL%LBNDU==1)THEN
            CALL CALC_KSWT(1)
            CALL DIAG_KSWT()
            CALL WRT0_KSWT(-4)
          ENDIF
        ENDIF
      ENDIF
100   CONTINUE
      CALL CPU_TIME(TA2); CALL SYSTEM_CLOCK(TIB2,TIRATE); TB2=REAL(TIB2,4)/REAL(TIRATE,4)
      IF(GP%MYRANK.EQ.GP%MASTER) CALL OUT_TIME_USE('TOTAL',TA2-TA1,TB2-TB1,GL%IO)
      IF(KPT%FILE_NAME/='')THEN
        CLOSE(IU_KGEN)
      ENDIF
      CLOSE(6)
      IF(GP%MYRANK.EQ.GP%MASTER) CALL gh5_end()
      CALL MPI_FINALIZE(IERR)
!
      END PROGRAM CYGUTZ

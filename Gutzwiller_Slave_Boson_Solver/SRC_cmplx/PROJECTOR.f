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
      MODULE GPROJECTOR
      USE gprec; USE SPARSE; USE FOCKSTATE; USE CORRORB
      USE LOCALHAMIL; USE gconstant
      USE GHDF5, ONLY: gh5_write_compound
      USE GHDF5_BASE, ONLY: gh5_write, log_file_id
      USE GINFO
      IMPLICIT NONE
!
! DEFINE TYPES
      TYPE PROJ
        INTEGER :: LMCFLY
        TYPE(ZCSR_MATRIX) :: FHL  ! F(U+M+N_var+N_phy)
        INTEGER :: N_PHIK,NEV_F,NCV_F,LEIGV,NCV_P
        TYPE(SECTOR1)    :: SEC_N
        INTEGER,POINTER :: ID_PHIK_N(:),ID_PHIK_J(:)
        TYPE(DCOO_MATRIX),POINTER :: SKIJ(:)
        ! begin \phi symmetrization.
        TYPE(ZBD_MATRIX) :: EVEC_SYM 
        TYPE(SECTOR1)    :: SEC_N_SYM
        INTEGER :: N_PHIK_SYM
        INTEGER,POINTER :: ID_PHIK_N_SYM(:)=>NULL()
        TYPE(DCOO_MATRIX),POINTER :: SKIJ_SYM(:)=>NULL()
        ! end \phi symmetrization.
        TYPE(IBD_MATRIX) :: ID_SKIJ ! Inverse table of SKIJ
        TYPE(ZCSR_MATRIX),POINTER :: NVAR_KKP(:,:)
        TYPE(ZBD_MATRIX),POINTER :: MKKP(:,:)
        TYPE(ZCSR_MATRIX) :: U_KKH
        REAL(gq)         :: EVL(2) ! Two lowest eigen values of HL
        COMPLEX(gq),POINTER :: C(:)=>NULL()   ! LOCAL CONFIGURATION PROBABILITY OR {c} in the notes
        COMPLEX(gq),POINTER :: C_BEST(:)=>NULL()
        INTEGER :: ITER
        INTEGER :: LENTANGLES,LSCF
        LOGICAL :: LSAVE_MNSV,LSET_MNSV
      END TYPE PROJ
!
      CONTAINS
!
!****************************************************************************
! CALC ATOM-TYPE DEPENDENT VARIATIONAL SETUP
!****************************************************************************
      SUBROUTINE SET_PROJ_TYPE(CO,FS,HL,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER MODE
!
      CALL MAP_HL_FS(FS,HL)
      PJ%SEC_N=HL%SEC_N
      IF((HL%LGPRJ/=21).AND.(HL%LGPRJ/=11).AND.(HL%LGPRJ/=14) &
          &.AND.(HL%LGPRJ/=99).AND.(HL%LGPRJ/=15))THEN
        CALL ALLOC_ZBM(HL%EVEC,HL%SEC_N%DIM,HL%SEC_N%ID,0,0)
      ENDIF
!
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15.OR.HL%LGPRJ==99)THEN 
        HL%EVEC%DIM=HL%DIM; HL%SEC_J=HL%SEC_N
        IF(HL%LGPRJ==11.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15)THEN
          HL%EVEC%NBK=-1
        ENDIF
      ELSEIF(HL%LGPRJ==12)THEN
        CALL SET_IDEN_ZBM(HL%EVEC)
        HL%SEC_J=HL%SEC_N
      ELSE
        CALL SET_HL_SLJ2(CO,FS,HL)
        IF(HL%LGPRJ/=21)THEN
          CALL DIAG_HL_SLJ2(HL)
        ELSE
          CALL ZBM_LOAD(HL%EVEC,HL%NI,0,0,16)
          CALL LOAD_SEC1(HL%SEC_J,HL%NI,18)
        ENDIF
        SELECT CASE(HL%LGPRJ)
        CASE(1,21)
          CALL SET_VNJ_JZ(HL%JZ,HL%JP,HL%SEC_J,HL%SEC_JZ,HL)
        CASE(2)
        CASE(3)
          CALL SET_VNJ_JZ(HL%SZ,HL%SP,HL%SEC_S,HL%SEC_SZ,HL)
          HL%SEC_J=HL%SEC_S; HL%SEC_JZ=HL%SEC_SZ
        CASE(4)
          HL%SEC_J=HL%SEC_SZ
        CASE DEFAULT
          STOP ' ERROR: ILLEGAL LGPRJ!'
        END SELECT
      ENDIF
      IF(PJ%LSAVE_MNSV)THEN
        IF(HL%EVEC%NBK>0.AND.GL%LH5WRT>=1)THEN
          IF(GP%MYRANK.EQ.GP%MASTER) CALL gh5_write_compound(HL%EVEC,"/Impurity_"// &
              &TRIM(int_to_str(HL%NI))//"/FockToCurrentBasis", log_file_id)
        ENDIF
        IF(GP%MYRANK.EQ.GP%MASTER) CALL gh5_write(HL%SEC_J%ID , HL%SEC_J%DIM+1, "/Impurity_"// &
            &TRIM(int_to_str(HL%NI))//"/Sec_ID", log_file_id)
        IF(GP%MYRANK.EQ.GP%MASTER) CALL gh5_write(HL%SEC_J%VAL, HL%SEC_J%DIM, "/Impurity_"// &
            &TRIM(int_to_str(HL%NI))//"/Sec_VAL", log_file_id)
      ENDIF
      IF(PJ%LSAVE_MNSV.AND.(HL%LGPRJ/=11).AND.(HL%LGPRJ/=14) &
          &.AND.(HL%LGPRJ/=99).AND.(HL%LGPRJ/=15))THEN
        CALL ZBM_DUMP(HL%EVEC,HL%NT,0,0,2)
      ENDIF
      IF(PJ%LSAVE_MNSV.AND.(GL%LPHISYM==-1))THEN
        CALL ZBM_DUMP(HL%EVEC,HL%NT,0,0,102)
      ELSEIF(GL%LPHISYM==1)THEN
        CALL ZBM_LOAD(PJ%EVEC_SYM,HL%NT,0,0,102)
      ENDIF
      MODE=HL%LDIAPJ
      IF(HL%LGPRJ==14.OR.HL%LGPRJ==15)MODE=2
      IF(HL%LGPRJ/=99)THEN
        CALL SET_SKIJ(FS,HL,PJ,MODE)
      ENDIF
      IF(HL%LGPRJ/=99)THEN
        CALL SET_NVAR_KKP(CO,FS,HL,PJ)
        CALL SET_MKKP(CO,FS,HL,PJ)
        IF(PJ%LSAVE_MNSV)THEN
          CALL SKIJ_DUMP(CO,PJ,HL)
          CALL ZCSR2_DUMP(PJ%NVAR_KKP,CO%DIM2,HL%NT,5)
          CALL ZBM2_DUMP(PJ%MKKP,CO%DIM2,HL%NT,6)
          IF(GL%LH5WRT>=3)THEN
            IF(GP%MYRANK.EQ.GP%MASTER) CALL GH5_SKIJ_DUMP(CO,PJ,HL)
          ENDIF
        ENDIF
      ENDIF
      IF(GL%LPHISYM==-1)THEN
        PJ%SEC_N_SYM = HL%SEC_N
        PJ%ID_PHIK_N_SYM => PJ%ID_PHIK_N
        PJ%SKIJ_SYM => PJ%SKIJ
        PJ%N_PHIK_SYM = PJ%N_PHIK
        CALL SKIJ_SYM_DUMP(PJ,HL%NT)
      ELSEIF(GL%LPHISYM==1)THEN
        CALL SKIJ_SYM_LOAD(PJ,HL%NT)
      ENDIF
      RETURN
!
      END SUBROUTINE SET_PROJ_TYPE
!
!****************************************************************************
! CALC ATOM-DEPENDENT VARIATIONAL SETUP
!****************************************************************************
      SUBROUTINE SET_PROJ_ATOM(CO,FS,HL,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      LOGICAL LSKIPU
!
      IF(PJ%LSCF==1)THEN
        LSKIPU=.FALSE.
      ELSE
        INQUIRE(FILE=TRIM(ADJUSTL(FILE_NAME(HL%NI,0,0,1))),EXIST=LSKIPU)
      ENDIF
      IF(HL%LGPRJ == 99) THEN
        CALL SET_HLOC_SPARSE(CO,FS,HL,LSKIPU)
      ELSE
        CALL SET_HLOC(CO,FS,HL,LSKIPU)
      ENDIF
      IF(GL%LH5WRT >= 1)THEN
        IF(GP%MYRANK.EQ.GP%MASTER) CALL gh5_write_compound(HL%H,"/Impurity_"//TRIM(int_to_str(HL%NI))// &
            &"/H", log_file_id)
      ENDIF
      IF(HL%LGPRJ == 99) RETURN
      IF(.NOT.LSKIPU)THEN
        IF((HL%LGPRJ/=11).AND.(HL%LGPRJ/=14).AND.(HL%LGPRJ/=15))THEN
          CALL ZBM_UHAU(HL%U,HL%EVEC)
          CALL ZBM_SIMPLIFY(HL%U,HL%SEC_J%DIM,HL%SEC_J%ID)
        ENDIF
        CALL ZBM_DUMP(HL%U,HL%NI,0,0,1)
        IF(PJ%LMCFLY>0)THEN
          CALL DEALLOC_ZBM(HL%U)
        ENDIF
      ELSEIF(PJ%LMCFLY<0)THEN
        CALL ZBM_LOAD(HL%U,HL%NI,0,0,1)
      ENDIF
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15)THEN
        CALL ZBMTOBKZCSR(HL%H)
        IF(PJ%LMCFLY==2)THEN
          PJ%U_KKH%NROW=0
        ELSE
          CALL SET_UKKH_FOCK(PJ%U_KKH,HL%H,HL,PJ)
        ENDIF
      ELSE
        CALL ZBM_UHAU(HL%H,HL%EVEC)
        CALL ZBM_SIMPLIFY(HL%H,HL%SEC_J%DIM,HL%SEC_J%ID)
        CALL SET_UKKH(PJ%U_KKH,HL%H,HL,PJ)
      ENDIF
      IF(.NOT.((HL%LGPRJ==11.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15) &
          &.AND.PJ%LMCFLY==2))THEN
        CALL DEALLOC_ZBM(HL%H)
      ENDIF
      IF(ASSOCIATED(PJ%C))THEN
        DEALLOCATE(PJ%C, PJ%C_BEST)
      ENDIF
      ALLOCATE(PJ%C(PJ%N_PHIK), PJ%C_BEST(PJ%N_PHIK)) ! Allocate c_{k}
      PJ%C=0; PJ%C_BEST(1) = 1.E2_gq
      RETURN
!
      END SUBROUTINE SET_PROJ_ATOM
!
!****************************************************************************
      SUBROUTINE SKIJ_DUMP(CO,PJ,HL)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      TYPE(LOCAL_HAMIL) HL
! LOCAL
      INTEGER I
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(HL%NT,0,0,4))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      WRITE(GIU)CO%ELC
      CALL WRT_SEC1(HL%SEC_J,GIU)
      WRITE(GIU)HL%SEC_N%DIM; WRITE(GIU)PJ%ID_PHIK_N
      WRITE(GIU)HL%SEC_J%DIM; WRITE(GIU)PJ%ID_PHIK_J
      WRITE(GIU)PJ%N_PHIK
      DO I=1,PJ%N_PHIK; CALL WRT_DCOO(PJ%SKIJ(I),GIU); ENDDO
      CALL IBM_WRT(PJ%ID_SKIJ,GIU)
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE SKIJ_DUMP
!
!****************************************************************************
      SUBROUTINE GH5_SKIJ_DUMP(CO,PJ,HL)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      TYPE(LOCAL_HAMIL) HL
! LOCAL
      TYPE(ZCSR_MATRIX) BCSR
      INTEGER I
!
      CALL gh5_write(PJ%N_PHIK,"/Impurity_"//TRIM(int_to_str(HL%NI))// &
          &"/NUM_HMAT_BASIS",log_file_id)
      DO I=1,PJ%N_PHIK
        CALL SM_DCOOZCSR(PJ%SKIJ(I),BCSR)
        CALL gh5_write_compound(BCSR,"/Impurity_"// &
            &TRIM(int_to_str(HL%NI))//"/HMAT_BASIS_"// &
            &TRIM(int_to_str(I)), log_file_id)
        CALL DEALLOC_ZCSR(BCSR)
      ENDDO
      RETURN
!
      END SUBROUTINE GH5_SKIJ_DUMP     
!
!****************************************************************************
      SUBROUTINE SKIJ_SYM_DUMP(PJ,NI)
      TYPE(PROJ)        PJ
      INTEGER,INTENT(IN) :: NI
! LOCAL
      INTEGER I
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,0,104))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      CALL WRT_SEC1(PJ%SEC_N_SYM,GIU)
      WRITE(GIU)PJ%ID_PHIK_N_SYM; WRITE(GIU)PJ%N_PHIK_SYM
      DO I=1,PJ%N_PHIK_SYM; CALL WRT_DCOO(PJ%SKIJ_SYM(I),GIU); ENDDO
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE SKIJ_SYM_DUMP
!
!****************************************************************************
      SUBROUTINE SKIJ_LOAD(CO,PJ,HL)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      TYPE(LOCAL_HAMIL) HL
! LOCAL
      INTEGER N,I
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(HL%NTMAP,0,0,4))),STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      READ(GIU)CO%ELC_OLD
      CALL READ_SEC1(HL%SEC_J,GIU)
      READ(GIU)N
      IF(N/=HL%SEC_N%DIM)THEN
        WRITE(0,'(" ERROR IN SKIJ_LOAD: N vs HL%SEC_N%DIM=",2I8)')N,HL%SEC_N%DIM; STOP
      ENDIF
      ALLOCATE(PJ%ID_PHIK_N(N+1)); READ(GIU)PJ%ID_PHIK_N
      READ(GIU)N
      IF(N/=HL%SEC_J%DIM)THEN
        WRITE(0,'(" ERROR IN SKIJ_LOAD: N vs HL%SEC_J%DIM=",2I8)')N,HL%SEC_J%DIM; STOP
      ENDIF
      ALLOCATE(PJ%ID_PHIK_J(N+1)); READ(GIU)PJ%ID_PHIK_J
      READ(GIU)PJ%N_PHIK
      ALLOCATE(PJ%SKIJ(PJ%N_PHIK))
      DO I=1,PJ%N_PHIK; CALL READ_DCOO(PJ%SKIJ(I),GIU); ENDDO
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15) &
          &CALL IBM_READ(PJ%ID_SKIJ,GIU)
      CLOSE(GIU)
      !(*,'(" N_PHIK=",I10)')PJ%N_PHIK
      IF(PJ%LMCFLY==0)THEN
        IF(PJ%N_PHIK>1000.OR.HL%DIM>1000)THEN
          PJ%LMCFLY= 2
        ELSE
          PJ%LMCFLY=-1
        ENDIF
      ENDIF
      IF(PJ%LEIGV/=0 .AND. PJ%LEIGV<=1)THEN
        PJ%LMCFLY=-1
      ENDIF
      RETURN
!
      END SUBROUTINE SKIJ_LOAD
!
!****************************************************************************
      SUBROUTINE SKIJ_SYM_LOAD(PJ,NI)
      TYPE(PROJ)        PJ
      INTEGER,INTENT(IN) :: NI
! LOCAL
      INTEGER I
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,0,104))),STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL",ERR=100)
      CALL READ_SEC1(PJ%SEC_N_SYM,GIU)
      ALLOCATE(PJ%ID_PHIK_N_SYM(PJ%SEC_N_SYM%DIM+1)); READ(GIU)PJ%ID_PHIK_N_SYM
      READ(GIU)PJ%N_PHIK_SYM
      ALLOCATE(PJ%SKIJ_SYM(PJ%N_PHIK_SYM))
      DO I=1,PJ%N_PHIK_SYM; CALL READ_DCOO(PJ%SKIJ_SYM(I),GIU); ENDDO
      CLOSE(GIU)
      RETURN
100   CONTINUE
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,0,204))),STATUS='OLD')
      CALL READ_SEC1_TXT(PJ%SEC_N_SYM,GIU)
      ALLOCATE(PJ%ID_PHIK_N_SYM(PJ%SEC_N_SYM%DIM+1)); READ(GIU,*)PJ%ID_PHIK_N_SYM
      READ(GIU,*)PJ%N_PHIK_SYM
      ALLOCATE(PJ%SKIJ_SYM(PJ%N_PHIK_SYM))
      DO I=1,PJ%N_PHIK_SYM; CALL READ_DCOO_TXT(PJ%SKIJ_SYM(I),GIU); ENDDO
      CLOSE(GIU)
      CALL SKIJ_SYM_DUMP(PJ,NI)
      RETURN
!
      END SUBROUTINE SKIJ_SYM_LOAD
!
!****************************************************************************
      SUBROUTINE CLR_PJ(CO,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER I1,I2
!
      DO I1=1,CO%DIM2; DO I2=1,CO%DIM2
        IF(PJ%MKKP(I1,I2)%DIM<=0)CYCLE
        CALL DEALLOC_ZBM(PJ%MKKP(I1,I2))
      ENDDO; ENDDO
      DO I1=1,CO%DIM2; DO I2=I1,CO%DIM2
        IF(PJ%NVAR_KKP(I1,I2)%NROW<=0)CYCLE
        CALL DEALLOC_ZCSR(PJ%NVAR_KKP(I1,I2))
      ENDDO; ENDDO
      DEALLOCATE(PJ%MKKP,PJ%NVAR_KKP)
      RETURN
!
      END SUBROUTINE CLR_PJ
!
!****************************************************************************
! READ ATOM-TYPE DEPENDENT VARIATIONAL SETUP
!****************************************************************************
      SUBROUTINE READ_PROJ_TYPE(IO,CO,FS,HL,PJ)
      INTEGER IO
      TYPE(CORR_ORB)    CO
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER IA1,IA2
!
      CALL MAP_HL_FS(FS,HL)
      PJ%SEC_N=HL%SEC_N
      CALL SKIJ_LOAD(CO,PJ,HL)
      IF(GL%LPHISYM==1)THEN
        CALL SKIJ_SYM_LOAD(PJ,HL%NTMAP)
      ENDIF
      WRITE(IO,'(" CO%ELC VS ELC_OLD:",2F8.3)')CO%ELC,CO%ELC_OLD
!
      ALLOCATE(PJ%NVAR_KKP(CO%DIM2,CO%DIM2))
      ALLOCATE(PJ%MKKP(CO%DIM2,CO%DIM2))
      PJ%NVAR_KKP(:,:)%NROW=-1; PJ%MKKP(:,:)%DIM=-1
      DO IA1=1,CO%DIM2; DO IA2=1,CO%DIM2
        IF(CO%M_INDEX_R(1,IA1,IA2)==IA1.AND.CO%M_INDEX_R(2,IA1,IA2)==IA2)THEN
          PJ%MKKP(IA1,IA2)%DIM=1
        ENDIF
        IF(CO%M_INDEX_L(1,IA1,IA2)==IA1.AND.CO%M_INDEX_L(2,IA1,IA2)==IA2.AND.IA1<=IA2)THEN
          PJ%NVAR_KKP(IA1,IA2)%NROW=1 ! Lower triangle enough due to hermitian.
        ENDIF
      ENDDO; ENDDO
!
      CALL ZCSR2_LOAD(PJ%NVAR_KKP,CO%DIM2,HL%NTMAP,5)
      CALL ZBM2_LOAD(PJ%MKKP,CO%DIM2,HL%NTMAP,6)
!
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15)THEN
        HL%EVEC%NBK=-1; HL%EVEC%DIM=HL%DIM
      ELSE
        CALL ZBM_LOAD(HL%EVEC,HL%NTMAP,0,0,2)
      ENDIF
      IF(GL%LPHISYM==1)THEN
        CALL ZBM_LOAD(PJ%EVEC_SYM,HL%NTMAP,0,0,102)
      ENDIF
      IF(PJ%LSAVE_MNSV)THEN
        IF((HL%EVEC%NBK>0).AND.(GL%LH5WRT>=1))THEN
          IF(GP%MYRANK.EQ.GP%MASTER) CALL gh5_write_compound(HL%EVEC,"/Impurity_"// &
              &TRIM(int_to_str(HL%NI))//"/FockToCurrentBasis", log_file_id)
        ENDIF
        IF(GP%MYRANK.EQ.GP%MASTER) CALL gh5_write(HL%SEC_J%ID, HL%SEC_J%DIM+1, "/Impurity_"// &
            &TRIM(int_to_str(HL%NI))//"/Sec_ID", log_file_id)
        IF(GP%MYRANK.EQ.GP%MASTER) CALL gh5_write(HL%SEC_J%VAL, HL%SEC_J%DIM, "/Impurity_"// &
            &TRIM(int_to_str(HL%NI))//"/Sec_VAL", log_file_id)
      ENDIF
      RETURN
!
      END SUBROUTINE READ_PROJ_TYPE
!
!****************************************************************************
      SUBROUTINE CALC_PJ_LC(CO,HL,PJ,IO,LSCF)
      TYPE(CORR_ORB)    CO
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      INTEGER IO,LSCF
! LOCAL
      LOGICAL :: LCINP
!
      IF(LSCF==105.OR.LSCF==107)THEN
        INQUIRE(FILE=TRIM(ADJUSTL(FILE_NAME(HL%NI,0,0,15))),EXIST=LCINP)
        IF(LCINP)THEN
          CALL LOAD_PJ_C(PJ,HL%NI)
          PJ%ITER=-1; PJ%EVL=0
          GOTO 100
        ENDIF
      ENDIF
      CALL CALC_FHL(CO,PJ)
      CALL DIAG_FHL(PJ,IO)
      IF(LSCF.EQ.-11)THEN
        CALL LOAD_PJ_C01(PJ,HL%NI)
      ENDIF
      IF(PJ%LENTANGLES==10)THEN
        CALL OPT_PHIK_PJENS_DS(PJ)
      ENDIF
      IF(PJ%LMCFLY<0)THEN
        CALL DEALLOC_ZCSR(PJ%FHL)
      ENDIF
100   CONTINUE
      CALL CALC_CZCSRC2(PJ%NVAR_KKP,CO%NC_VAR,CO,PJ,-1)
      RETURN
!
      END SUBROUTINE CALC_PJ_LC
!
!****************************************************************************
      SUBROUTINE CALC_PJ_NCVAR(CO,HL,PJ,IO)
      TYPE(CORR_ORB)    CO
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      INTEGER IO
!
      CALL CALC_FHL(CO,PJ)
      CALL DIAG_FHL(PJ,IO)
      IF(PJ%LMCFLY<0)THEN
        CALL DEALLOC_ZCSR(PJ%FHL)
      ENDIF
      CALL CALC_CZCSRC2(PJ%NVAR_KKP,CO%NC_VAR,CO,PJ,-1)
      RETURN
!
      END SUBROUTINE CALC_PJ_NCVAR
!
!****************************************************************************
      SUBROUTINE CALC_R(CO,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
!
      CALL CALC_CZBMC_COSYM(PJ%MKKP,CO%R0,CO,PJ)
      CO%R=MATMUL(TRANSPOSE(CO%ISIMIX),CO%R0)
      RETURN
!
      END SUBROUTINE CALC_R
!
!******************************************************************************
      SUBROUTINE SET_SKIJ(FS,HL,PJ,MODE)
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      INTEGER MODE
! LOCAL
      INTEGER NPHIK
!
      ALLOCATE(PJ%ID_PHIK_N(HL%SEC_N%DIM+1))
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15)THEN
        CALL GET_PHIK_DIM_N(HL%SEC_N,PJ%ID_PHIK_N,NPHIK,MODE,FS%BS_JSZ(FS%DIM-HL%DIM+1:),HL%DIM)
        PJ%ID_PHIK_J=>PJ%ID_PHIK_N
      ELSE
        ALLOCATE(PJ%ID_PHIK_J(HL%SEC_J%DIM+1))
        CALL GET_PHIK_DIM_NJ(HL%SEC_N,HL%SEC_J,PJ%ID_PHIK_N,PJ%ID_PHIK_J,NPHIK,HL%LGPRJ,MODE)
      ENDIF
      ALLOCATE(PJ%SKIJ(NPHIK))
      PJ%N_PHIK=NPHIK
      !(*,'(" N_PHIK=",I10)')PJ%N_PHIK
      IF(PJ%LMCFLY==0)THEN
        IF(PJ%N_PHIK>1000.OR.HL%DIM>1000)THEN
          PJ%LMCFLY= 2
        ELSE
          PJ%LMCFLY=-1
        ENDIF
      ENDIF
      IF(PJ%LEIGV/=0 .AND. PJ%LEIGV<=1)THEN
        PJ%LMCFLY=-1
      ENDIF
      IF(HL%LGPRJ==2.OR.HL%LGPRJ==4.OR.HL%LGPRJ==11.OR.HL%LGPRJ==12 &
          &.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15)THEN
        CALL SET_SKIJ_N(PJ%SKIJ,HL%SEC_J,NPHIK,MODE,PJ%ID_SKIJ,FS%BS_JSZ(FS%DIM-HL%DIM+1:))
      ELSE
        CALL SET_SKIJ_J(PJ%SKIJ,HL%SEC_J,NPHIK,MODE)
      ENDIF
      RETURN
!
      END SUBROUTINE SET_SKIJ
!
!****************************************************************************
      SUBROUTINE GET_PHIK_DIM_N(SEC_N,ID_PHIK_N,NPHIK,MODE,BS_JSZ,NBS)
      TYPE(SECTOR1)::SEC_N
      INTEGER NPHIK,NBS,MODE,ID_PHIK_N(SEC_N%DIM+1),BS_JSZ(NBS)
! LOCAL
      INTEGER I,N,IBASE,NPH1
!
      NPHIK=0; ID_PHIK_N(1)=1; IBASE=0
      DO I=1,SEC_N%DIM
        N=SEC_N%ID(I+1)-SEC_N%ID(I)
        IF(MODE==1)THEN
          NPHIK=NPHIK+N
        ELSEIF(MODE==2)THEN
          CALL GET_NPHIK_NBJSZ(NPH1,BS_JSZ(IBASE+1:IBASE+N),N)
          NPHIK=NPHIK+NPH1
        ELSE
          NPHIK=NPHIK+N*(SEC_N%ID_FROZEN(I+1)-SEC_N%ID_FROZEN(I))
        ENDIF
        ID_PHIK_N(I+1)=NPHIK+1
        IBASE=IBASE+N
      ENDDO
      RETURN
!
      END SUBROUTINE GET_PHIK_DIM_N
!
!****************************************************************************
      SUBROUTINE GET_PHIK_DIM_NJ(SEC_N,SEC_J,ID_PHIK_N,ID_PHIK_J,NPHIK,LGPRJ,MODE)
      TYPE(SECTOR1)::SEC_N,SEC_J
      INTEGER NPHIK,LGPRJ,MODE,ID_PHIK_N(SEC_N%DIM+1),ID_PHIK_J(SEC_J%DIM+1)
! LOCAL
      INTEGER I,J,N,MULJ,INBK
      REAL(gq) RJ
!
      NPHIK=0; INBK=1; ID_PHIK_N=0; ID_PHIK_N(1)=1; ID_PHIK_J=0; ID_PHIK_J(1)=1
      DO I=1,SEC_J%DIM
        IF(LGPRJ==1.OR.LGPRJ==3.OR.LGPRJ==21)THEN
          RJ=SQRT(SEC_J%VAL(I)+0.25_gq)-0.5_gq; MULJ=NINT(2*RJ+1)
        ELSE
          MULJ=1
        ENDIF
        N=(SEC_J%ID(I+1)-SEC_J%ID(I))/MULJ ! Number of repeated irreducible representation
        IF(MODE==1)THEN
          NPHIK=NPHIK+N
        ELSE
          NPHIK=NPHIK+N*N
        ENDIF
        DO J=INBK,SEC_N%DIM
          IF(SEC_J%ID(I)<SEC_N%ID(J+1))EXIT
          IF(ID_PHIK_N(J+1)==0)THEN
            ID_PHIK_N(J+1)=ID_PHIK_N(J)
          ENDIF
        ENDDO
        INBK=J
        ID_PHIK_N(INBK+1)=NPHIK+1
        ID_PHIK_J(I+1)=NPHIK+1
      ENDDO
      RETURN
!
      END SUBROUTINE GET_PHIK_DIM_NJ
!
!****************************************************************************
      SUBROUTINE SET_SKIJ_N(SKIJ,SEC_N,NPHIK,MODE,ID_SKIJ,BS_JSZ)
      INTEGER NPHIK,MODE,BS_JSZ(*)
      TYPE(DCOO_MATRIX) :: SKIJ(NPHIK)
      TYPE(IBD_MATRIX) :: ID_SKIJ ! Inverse table of SKIJ
      TYPE(SECTOR1)::SEC_N
! LOCAL
      INTEGER I,N,N1,N2,N2_,IK,IBASE,MULN,MN,NDIM,SZ1,SZ2
      REAL(gq) RES
!
      IK=0; NDIM=SEC_N%ID(SEC_N%DIM+1)-1; MULN=1; RES=1._gq
      ID_SKIJ%NBK=SEC_N%DIM
      ALLOCATE(ID_SKIJ%BK(SEC_N%DIM))
      DO I=1,SEC_N%DIM
        IBASE=SEC_N%ID(I)-1
        N=SEC_N%ID(I+1)-SEC_N%ID(I)
        ID_SKIJ%BK(I)%NROW=N; ID_SKIJ%BK(I)%NCOL=N
        ALLOCATE(ID_SKIJ%BK(I)%A(N,N)); ID_SKIJ%BK(I)%A=0
        DO N1=1,N
        SZ1=BS_JSZ(IBASE+N1)
        DO N2=SEC_N%ID_FROZEN(I),SEC_N%ID_FROZEN(I+1)-1
          N2_=SEC_N%I_BS(N2)-IBASE
          IF(MODE==1)THEN
            IF(N2_/=N1)CYCLE
          ELSEIF(MODE==2)THEN
            SZ2=BS_JSZ(IBASE+N2_)
            IF(SZ1/=SZ2)CYCLE
          ENDIF
          IK=IK+1
          CALL ALLOC_DCOO(SKIJ(IK),MULN,NDIM,NDIM)
          SKIJ(IK)%A=RES
          DO MN=1,MULN
            SKIJ(IK)%I(MN)=IBASE+N1+(MN-1)*N
            SKIJ(IK)%J(MN)=IBASE+N2_+(MN-1)*N
          ENDDO
          ID_SKIJ%BK(I)%A(N1,N2_)=IK
        ENDDO; ENDDO
      ENDDO
      ID_SKIJ%DIM=SUM(ID_SKIJ%BK(:)%NROW)
      RETURN
!
      END SUBROUTINE SET_SKIJ_N
!
!****************************************************************************
      SUBROUTINE SET_SKIJ_J(SKIJ,SEC_J,NPHIK,MODE)
      INTEGER NPHIK,MODE
      TYPE(DCOO_MATRIX) :: SKIJ(NPHIK)
      TYPE(SECTOR1)::SEC_J
! LOCAL
      INTEGER I,N,N1,N2,IK,IBASE,MULJ,MJ,NDIM,NSTR,NEND
      REAL(gq) RES,RJ
!
      IK=0; NDIM=SEC_J%ID(SEC_J%DIM+1)-1
      DO I=1,SEC_J%DIM
        IBASE=SEC_J%ID(I)-1
        RJ=SQRT(SEC_J%VAL(I)+0.25_gq)-0.5_gq; MULJ=NINT(2*RJ+1)
        N=(SEC_J%ID(I+1)-SEC_J%ID(I))/MULJ ! Number of repeated irreducible representation
        RES=SQRT(1._gq/MULJ)
        DO N1=1,N
        IF(MODE==1)THEN  ! Diagonal 
          NSTR=N1; NEND=N1
        ELSE             ! Full
          NSTR=1; NEND=N
        ENDIF
        DO N2=NSTR,NEND
          IK=IK+1
          CALL ALLOC_DCOO(SKIJ(IK),MULJ,NDIM,NDIM)
          SKIJ(IK)%A=RES
          DO MJ=1,MULJ
            SKIJ(IK)%I(MJ)=IBASE+N1+(MJ-1)*N
            SKIJ(IK)%J(MJ)=IBASE+N2+(MJ-1)*N
          ENDDO
        ENDDO; ENDDO
      ENDDO
      RETURN
!
      END SUBROUTINE SET_SKIJ_J
!
!****************************************************************************
      SUBROUTINE GET_NPHIK_NBJSZ(NPHIK,BS_JSZ,N)
      INTEGER NPHIK,N,BS_JSZ(N)
! LOCAL
      INTEGER I,IMIN,IMAX
      INTEGER,ALLOCATABLE::ICOUNT(:)
!
      NPHIK=0
      IMIN=MINVAL(BS_JSZ); IMAX=MAXVAL(BS_JSZ)
      ALLOCATE(ICOUNT(IMIN:IMAX)); ICOUNT=0
      DO I=1,N; ICOUNT(BS_JSZ(I))=ICOUNT(BS_JSZ(I))+1; ENDDO
      DO I=IMIN,IMAX; NPHIK=NPHIK+ICOUNT(I)**2; ENDDO
      RETURN
!
      END SUBROUTINE GET_NPHIK_NBJSZ
!
!****************************************************************************
! Only include U+N parts, M will be evaluated on fly
!****************************************************************************
      SUBROUTINE CALC_FHL(CO,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      !
!
      IF(PJ%LMCFLY>0)RETURN
      !"CALC_FHL"
      !
      CALL CALC_U_FHL(PJ%FHL,PJ)
      CALL CALC_N_FHL(PJ%FHL,CO,PJ)
      CALL CALC_M_FHL(PJ%FHL,CO,PJ)
      CALL ZCSR_SIMPLIFY(PJ%FHL,1.E-16_gq,'L')
      !
      !"CALC_FHL"
      RETURN
!
      END SUBROUTINE CALC_FHL
!
!****************************************************************************
! Matrix element contribution from U
!****************************************************************************
      SUBROUTINE CALC_U_FHL(A,PJ)
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX) A
!
      CALL ZCSR_COPY(PJ%U_KKH,A)
      RETURN
!
      END SUBROUTINE CALC_U_FHL
!
!****************************************************************************
! Matrix element contribution from N
!****************************************************************************
      SUBROUTINE CALC_N_FHL(A,CO,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX) A
! LOCAL
      INTEGER I1,I2
      REAL(gq),PARAMETER::RTOL=1.E-10_gq
      COMPLEX(gq) ZES
      !
!
      !
      DO I1=1,CO%DIM2; DO I2=I1,CO%DIM2
        IF(PJ%NVAR_KKP(I1,I2)%NROW<=0)CYCLE
        ZES=CO%LA2(I1,I2)
        IF(ABS(ZES)<RTOL)CYCLE
        CALL ZCSR_APLSB_SK('N',A,ZES,PJ%NVAR_KKP(I1,I2))
        IF(I1==I2)CYCLE
        ZES=CONJG(ZES)
        CALL ZCSR_APLSB_SK('C',A,ZES,PJ%NVAR_KKP(I1,I2)) ! h.c.
      ENDDO; ENDDO
      !
      RETURN
!
      END SUBROUTINE CALC_N_FHL
!
!****************************************************************************
! Matrix element contribution from M
!****************************************************************************
      SUBROUTINE CALC_M_FHL(A,CO,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX) A
! LOCAL
      INTEGER I1,I2
      TYPE(ZCSR_MATRIX) MKKP
      REAL(gq),PARAMETER::RTOL=1.E-10_gq
      !
!
      !
      DO I1=1,CO%DIM2; DO I2=1,CO%DIM2
        IF(PJ%MKKP(I1,I2)%DIM<=0)CYCLE
        IF(ABS(CO%D(I1,I2))<RTOL)CYCLE
        CALL ZBMTOZCSR(PJ%MKKP(I1,I2),MKKP)
        CALL ZCSR_APLSB_SK('N',A,CO%D(I1,I2),MKKP)
        CALL DEALLOC_ZCSR(MKKP)
      ENDDO; ENDDO
      !
      RETURN
!
      END SUBROUTINE CALC_M_FHL
!
!****************************************************************************
      SUBROUTINE DIAG_FHL(PJ,IO)
      use gprimme
      TYPE(PROJ)       PJ
      INTEGER IO
! LOCAL
      INTEGER NDIM,NEV,NCV,NCONV,I_FIX,NEST
      COMPLEX(gq) SIGMA
      REAL(gq) ABSTOL, SIGMA_R, NORM
      REAL(gq),ALLOCATABLE::WE(:)
      COMPLEX(gq),ALLOCATABLE::ZE(:,:)
      EXTERNAL :: PHIK_AV, PHIK_AV_M
!
      !"DIAG_FHL"
      SIGMA=0; SIGMA_R=0; ABSTOL=0; I_FIX=0
      NDIM=PJ%N_PHIK; NEV=MIN(PJ%NEV_F,NDIM-2); NEV=MAX(NEV,1)
      IF(PJ%NCV_F==0)THEN
        NCV=MIN(NDIM,36*NEV)
      ELSE
        NCV=MIN(PJ%NCV_F,NDIM)
      ENDIF
      PJ%ITER=0
      IF(PJ%LMCFLY>0 .AND. PJ%LEIGV==0)THEN
        PJ%LEIGV=2
      ENDIF
      IF(PJ%LEIGV==0)THEN
        IF(NDIM.GT.100)THEN
          PJ%LEIGV=2
        ELSE
          PJ%LEIGV=1
        ENDIF
      ENDIF
      IF(PJ%LEIGV > 1)THEN
        ALLOCATE(WE(NEV),ZE(NDIM,NEV))
100     CONTINUE
        WE=0; ZE=0
        IF(PJ%LEIGV==2)THEN
          CALL DSDRV1(NDIM *2,NEV,NCV,'I','SA',ABSTOL,SIGMA_R,.TRUE.,WE,ZE,3000,NCONV,PHIK_AV)
        ELSEIF(PJ%LEIGV==5)THEN
! It canot converge in some cases.
!          IF(REAL(PJ%C_BEST(1)) > 10._gq)THEN
            NEST = 0
!          ELSE
!            NEST = 1
!            ZE(:,1) = PJ%C_BEST
!          ENDIF
          ! Need this high accuracy for the solver.
          CALL Zprimme(PJ%NCV_P, 1, 3000, 1.E-10_gq, 1, 0, 0, NDIM, NEST, WE(1), ZE(1,1), PHIK_AV_M)
          NCONV = 1
        ELSE
          CALL ZHDRV1(NDIM,NEV,NCV,'I','SR',ABSTOL,SIGMA,.TRUE.,WE,ZE,3000,NCONV,PHIK_AV)
        ENDIF
        IF(NCONV.EQ.0)THEN
          IF(I_FIX < 7)THEN
            NCV = MIN(NCV + 100, NDIM - 1)
            I_FIX = I_FIX + 1; GOTO 100 
          ENDIF
          IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" SEVERE WARNING: ZHDRV1--No converged Ritz values!")')
          IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(0 ,'(" SEVERE WARNING: ZHDRV1--No converged Ritz values!")')
        ELSEIF(NCONV.EQ.1)THEN
          IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" WARNING-ZHDRV1: Only 1 converged Ritz value!")')
        ENDIF
      ELSE
        ALLOCATE(ZE(NDIM,NDIM),WE(NDIM)); ZE=0; WE=0
        CALL SM_ZCSRZDNS(PJ%FHL,ZE)
        IF (PJ%LEIGV == 1) THEN
          CALL HERMEV('V','L',ZE,WE,NDIM)
        ELSE
          CALL HERMEV('V','L',ZE,WE,NDIM,1,2)
        ENDIF
      ENDIF
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" PJ%C_OLD.PJ%C  =",2F12.8)')DOT_PRODUCT(PJ%C,ZE(:,1))
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" PJ%C VAL CHANGE=",2F12.8)')PJ%EVL-WE(1:2)
      PJ%C=ZE(:,1)
      PJ%EVL=WE(1:2)
      DEALLOCATE(WE,ZE)
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" DIAG_FHL ITER=",I5," PJ%EVL=",2F20.8)')PJ%ITER,PJ%EVL
      IF(ABS(PJ%EVL(1)-PJ%EVL(2))<1.E-6_gq)THEN
        IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" WARNING: ALMOST DEGENERATE GROUND STATE!")')
      ENDIF
      IF(GL%LPHISYM==1)THEN
        CALL SYM_SKIJ(PJ,NORM)
        IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" SYM_SKIJ NORM = ", F16.10)')NORM
      ENDIF
      !"DIAG_FHL"
      RETURN
!
      END SUBROUTINE DIAG_FHL
!
!****************************************************************************
! For CO%NC_VAR
!****************************************************************************
      SUBROUTINE CALC_CZCSRC2(A,X,CO,PJ,MODE)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX)::A(CO%DIM2,CO%DIM2)
      COMPLEX(gq)     ::X(CO%DIM2,CO%DIM2)
      INTEGER            MODE
! LOCAL
      INTEGER I1,I2,I1P,I2P
      COMPLEX(gq) ZES
!
      X=0
      DO I1=1,CO%DIM2; DO I2=I1,CO%DIM2
        I1P=CO%M_INDEX_L(1,I1,I2); I2P=CO%M_INDEX_L(2,I1,I2)
        IF(I1P==0)CYCLE
        IF(I1P==I1.AND.I2P==I2)THEN
          CALL CALC_CZCSRCP(PJ%C,PJ%C,A(I1,I2),ZES,MODE)
          X(I1,I2)=ZES/CO%M_INDEX_L(3,I1,I2)
        ELSE
          X(I1,I2)=X(I1P,I2P)
        ENDIF
        IF(I1/=I2)THEN
          X(I2,I1)=CONJG(X(I1,I2))
        ENDIF
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE CALC_CZCSRC2
!
!****************************************************************************
! For CO%R0, R0_{ij} = <\phi | M_{ij} |\phi >
!****************************************************************************
      SUBROUTINE CALC_CZBMC_COSYM(A,X,CO,PJ)
      TYPE(CORR_ORB),INTENT(IN)::CO
      TYPE(PROJ),INTENT(IN)::PJ
      TYPE(ZBD_MATRIX),INTENT(IN)::A(CO%DIM2,CO%DIM2)
      COMPLEX(gq),INTENT(OUT)::X(CO%DIM2,CO%DIM2)
! LOCAL
      INTEGER IA1,IA2
      COMPLEX(gq) ZES
!
      X=0
      DO IA1=1,CO%DIM2; DO IA2=1,CO%DIM2
        IF(A(IA1,IA2)%DIM<=0)CYCLE
        CALL CALC_CZBMC(PJ%C,A(IA1,IA2),ZES,.FALSE.)
        X(IA1,IA2)=ZES
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE CALC_CZBMC_COSYM
!
!****************************************************************************
      SUBROUTINE CALC_EGAMM(HL,PJ)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      COMPLEX(gq) ZES
!
      IF(PJ%U_KKH%NROW>0)THEN
        CALL CALC_CZCSRCP(PJ%C,PJ%C,PJ%U_KKH,ZES,1)
        HL%EGAMM=REAL(ZES,gq)
      ELSE
        CALL CALC_EGAMM_HL(HL)
      ENDIF
      RETURN
!
      END SUBROUTINE CALC_EGAMM
!
!****************************************************************************
      SUBROUTINE CALC_PJ_RHO(HL,PJ,MODE)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      INTEGER MODE
! LOCAL
      INTEGER NBK,IBK,IBASE,K,IJ,I,J
      INTEGER,POINTER::ID(:)
      COMPLEX(gq),EXTERNAL::ZDOTC
!
      NBK=HL%SEC_J%DIM; ID=>HL%SEC_J%ID
      IF(HL%RHO%NBK<=0)THEN
        CALL ALLOC_ZBM(HL%RHO,NBK,ID,0,0)
      ENDIF
      DO IBK=1,NBK
      IBASE=HL%SEC_J%ID(IBK)-1; HL%RHO%BK(IBK)%A=0
      DO K=PJ%ID_PHIK_J(IBK),PJ%ID_PHIK_J(IBK+1)-1
      DO IJ=1,PJ%SKIJ(K)%NNZ
        I=PJ%SKIJ(K)%I(IJ)-IBASE; J=PJ%SKIJ(K)%J(IJ)-IBASE
        HL%RHO%BK(IBK)%A(I,J)=HL%RHO%BK(IBK)%A(I,J)+PJ%SKIJ(K)%A(IJ)*PJ%C(K)
      ENDDO; ENDDO; ENDDO
      IF(GL%LH5WRT>=1.AND.MODE>0)THEN
        IF(GP%MYRANK.EQ.GP%MASTER) CALL gh5_write_compound(HL%RHO,"/Impurity_"// &
            &TRIM(int_to_str(HL%NI))//"/phi", log_file_id)
      ENDIF
      CALL ZBM_GEMM('N','C',HL%RHO,HL%RHO) ! phi phi*
      NULLIFY(ID)
      RETURN
!
      END SUBROUTINE CALC_PJ_RHO
!
!****************************************************************************
      SUBROUTINE SET_NVAR_KKP(CO,FS,HL,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER I1,I2,I1P,I2P
      TYPE(ZCSR_MATRIX) :: NCSR
      TYPE(ZCSR_MATRIX) :: NKKP
      !
!
      ! "SET_NVAR_KKP"
      !
      ALLOCATE(PJ%NVAR_KKP(CO%DIM2,CO%DIM2))
      PJ%NVAR_KKP%NROW=-1
      DO I1=1,CO%DIM2; DO I2=I1,CO%DIM2
        I1P=CO%M_INDEX_L(1,I1,I2); I2P=CO%M_INDEX_L(2,I1,I2)
        IF(I1P==0)CYCLE
        CALL CALC_NIJ_FOCK(NCSR,FS,FS%DIM-HL%DIM+1,FS%DIM,I1,I2)
        IF(I1P==I1.AND.I2P==I2)THEN
          CALL SET_NKKP(PJ%NVAR_KKP(I1,I2),NCSR,HL,PJ,-1)
        ELSE
          CALL SET_NKKP(NKKP,NCSR,HL,PJ,-1)
          CALL ZCSR_APLSB_SK('N',PJ%NVAR_KKP(I1P,I2P),Z1,NKKP)
          CALL DEALLOC_ZCSR(NKKP)
        ENDIF
        CALL DEALLOC_ZCSR(NCSR)
      ENDDO; ENDDO
      !
      ! "SET_NVAR_KKP"
      RETURN
!
      END SUBROUTINE SET_NVAR_KKP
!
!****************************************************************************
! U tensor; assume block diagonal in |N> and Hermitian
!****************************************************************************
      SUBROUTINE SET_UKKH_FOCK(UKKH,U,HL,PJ)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX):: UKKH
      TYPE(ZBD_MATRIX) :: U
      !
! LOCAL
      INTEGER NNZ
!
      ! "SET_UKKH_FOCK"
      !
      CALL COUNT_NNZ_UKKH_FOCK(U,HL,PJ,NNZ)
      CALL ALLOC_ZCSR(UKKH,NNZ,PJ%N_PHIK,PJ%N_PHIK)
      CALL CALC_UKKH_FOCK(UKKH,U,HL,PJ)
      !
      ! "SET_UKKH_FOCK"
      RETURN
!
      END SUBROUTINE SET_UKKH_FOCK
!
!****************************************************************************
! U tensor; assume block diagonal in |N,J> and Hermitian
!****************************************************************************
      SUBROUTINE SET_UKKH(UKKH,U,HL,PJ)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX):: UKKH
      TYPE(ZBD_MATRIX) :: U
      !
!
      ! "SET_UKKH"
      !
      CALL ALLOC_ZCSR_KKP(UKKH,PJ%ID_PHIK_J,HL%SEC_J%DIM,0,1)
      CALL CALC_UKKH(UKKH,U,HL,PJ)
      CALL ZCSR_SIMPLIFY(UKKH,1.E-16_gq,'L')
      !
      ! "SET_UKKH"
      RETURN
!
      END SUBROUTINE SET_UKKH
!
!****************************************************************************
! Trace( S(k)^\dagger U SP(k') )
! |G2><G1| U |G3><G2|
!****************************************************************************
      SUBROUTINE CALC_UKKH_FOCK(UKKH,U,HL,PJ)
      TYPE(ZCSR_MATRIX):: UKKH
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(ZBD_MATRIX) :: U
! LOCAL
      INTEGER I,J,K,N,INZ,IR,JR,NFK,JPHIK
!
      INZ=0; UKKH%I(1)=1
      DO N=1,HL%SEC_N%DIM
      NFK=HL%SEC_N%ID(N+1)-HL%SEC_N%ID(N)
      DO K=PJ%ID_PHIK_N(N),PJ%ID_PHIK_N(N+1)-1
      IR=PJ%SKIJ(K)%I(1)-HL%SEC_N%ID(N)+1
      JR=PJ%SKIJ(K)%J(1)-HL%SEC_N%ID(N)+1
      DO I=U%BK(N)%ACSR%I(IR),U%BK(N)%ACSR%I(IR+1)-1
        J=U%BK(N)%ACSR%J(I)
        JPHIK=PJ%ID_SKIJ%BK(N)%A(J,JR)
        IF(JPHIK<=0.OR.JPHIK>K)CYCLE ! Lower triangle Hermitian part
        INZ=INZ+1
        UKKH%J(INZ)=JPHIK
        UKKH%A(INZ)=U%BK(N)%ACSR%A(I)
      ENDDO ! I
      UKKH%I(K+1)=INZ+1
      ENDDO; ENDDO ! K,N
      RETURN
!
      END SUBROUTINE CALC_UKKH_FOCK
!
!****************************************************************************
      SUBROUTINE COUNT_NNZ_UKKH_FOCK(U,HL,PJ,INZ)
      TYPE(LOCAL_HAMIL),INTENT(IN) :: HL
      TYPE(PROJ),INTENT(IN) :: PJ
      TYPE(ZBD_MATRIX),INTENT(IN) :: U
      INTEGER,INTENT(OUT) :: INZ
! LOCAL
      INTEGER I,J,K,N,IR,JR,NFK,JPHIK
!
      INZ=0
      DO N=1,HL%SEC_N%DIM
      NFK=HL%SEC_N%ID(N+1)-HL%SEC_N%ID(N)
      DO K=PJ%ID_PHIK_N(N),PJ%ID_PHIK_N(N+1)-1
      IR=PJ%SKIJ(K)%I(1)-HL%SEC_N%ID(N)+1
      JR=PJ%SKIJ(K)%J(1)-HL%SEC_N%ID(N)+1
      DO I=U%BK(N)%ACSR%I(IR),U%BK(N)%ACSR%I(IR+1)-1
        J=U%BK(N)%ACSR%J(I)
        JPHIK=PJ%ID_SKIJ%BK(N)%A(J,JR)
        IF(JPHIK<=0.OR.JPHIK>K)CYCLE ! Lower triangle Hermitian part
        INZ=INZ+1
      ENDDO ! I
      ENDDO; ENDDO ! K,N
      RETURN
!
      END SUBROUTINE COUNT_NNZ_UKKH_FOCK
!
!****************************************************************************
! Trace( S(k)^\dagger U SP(k') )
! |G2><G1| U |G3><G2|
!****************************************************************************
      SUBROUTINE APPLY_UKKH_FOCK(HL,PJ,X,HX)
      TYPE(LOCAL_HAMIL),INTENT(IN) :: HL
      TYPE(PROJ),INTENT(IN)        :: PJ
      COMPLEX(gq),INTENT(IN) :: X(PJ%N_PHIK)
      COMPLEX(gq),INTENT(OUT) :: HX(PJ%N_PHIK)
! LOCAL
      INTEGER I,J,K,N,IR,JR,NFK,JPHIK
!
      HX=0
      DO N=1,HL%SEC_N%DIM
      NFK=HL%SEC_N%ID(N+1)-HL%SEC_N%ID(N)
      DO K=PJ%ID_PHIK_N(N),PJ%ID_PHIK_N(N+1)-1
      IF(ABS(X(K))<SMALL)CYCLE
      IR=PJ%SKIJ(K)%I(1)-HL%SEC_N%ID(N)+1
      JR=PJ%SKIJ(K)%J(1)-HL%SEC_N%ID(N)+1
      DO I=HL%H%BK(N)%ACSR%I(IR),HL%H%BK(N)%ACSR%I(IR+1)-1
        J=HL%H%BK(N)%ACSR%J(I)
        JPHIK=PJ%ID_SKIJ%BK(N)%A(J,JR)
        IF(JPHIK<=0)CYCLE 
        HX(JPHIK)=HX(JPHIK)+CONJG(HL%H%BK(N)%ACSR%A(I))*X(K)
      ENDDO ! I
      ENDDO; ENDDO ! K,N
      RETURN
!
      END SUBROUTINE APPLY_UKKH_FOCK
!
!****************************************************************************
! Trace( S(k)^\dagger U SP(k') )
!****************************************************************************
      SUBROUTINE CALC_UKKH(UKKH,U,HL,PJ)
      TYPE(ZCSR_MATRIX):: UKKH
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(ZBD_MATRIX) :: U
! LOCAL
      INTEGER K,KP,J,IBASE,ISUM
!
      ISUM=1; IBASE=0; UKKH%I(1)=1
      DO J=1,HL%SEC_J%DIM
      DO K =PJ%ID_PHIK_J(J),PJ%ID_PHIK_J(J+1)-1
      DO KP=PJ%ID_PHIK_J(J),K
        CALL TR_DCOOZDNSDCOO(UKKH%A(ISUM),PJ%SKIJ(K),U%BK(J)%A,PJ%SKIJ(KP),U%BK(J)%NROW,IBASE,1)
        IF(ABS(UKKH%A(ISUM))<1.E-16_gq)CYCLE
        UKKH%J(ISUM)=KP
        ISUM=ISUM+1
      ENDDO
      UKKH%I(K+1)=ISUM
      ENDDO
      IBASE=IBASE+U%BK(J)%NROW
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_UKKH
!
!****************************************************************************
! MODE = -1: N_VAR; 1: N_PHY
!****************************************************************************
      SUBROUTINE SET_NKKP(NKKP,NCSR,HL,PJ,MODE)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX) :: NKKP
      TYPE(ZCSR_MATRIX) :: NCSR
      INTEGER MODE
! LOCAL
      !
!
      ! "SET_NKKP"
      !
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15)THEN
        CALL ALLOC_ZCSR_KKP(NKKP,PJ%ID_PHIK_J,HL%SEC_J%DIM,1,MODE)
        CALL CALC_NKKP_FOCK(NKKP,NCSR,HL,PJ,MODE)
      ELSE
        CALL ALLOC_ZCSR_KKP(NKKP,PJ%ID_PHIK_J,HL%SEC_J%DIM,0,MODE)
        CALL CALC_NKKP(NKKP,NCSR,HL,PJ,MODE)
      ENDIF
      CALL ZCSR_SHRINK(NKKP)
      !
      ! "SET_NKKP"
      RETURN
!
      END SUBROUTINE SET_NKKP
!
!****************************************************************************
! MODE > 0: Trace( S^\dagger NCSR S ) -- |G2><G1| N |G3><G2|
!     else: Trace( S^\dagger S NCSR ) -- |G2><G1| |G1><G3| N |G2>
!****************************************************************************
      SUBROUTINE CALC_NKKP_FOCK(NKKP,NCSR,HL,PJ,MODE)
      TYPE(ZCSR_MATRIX):: NKKP
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX) :: NCSR
      INTEGER MODE
! LOCAL
      INTEGER N,K,KP,ISUM,IFK1,IFK2,IFK3,J
      !
!
      ! "CALC_NKKP_FOCK"
      !
      IF(MODE<=0)CALL SM_ZCSRZCSC(NCSR) ! <F_j| N |F_i>
      NKKP%I(1)=1; ISUM=1
      DO N=1,HL%SEC_N%DIM
      DO K=PJ%ID_PHIK_N(N),PJ%ID_PHIK_N(N+1)-1
      IFK1=PJ%SKIJ(K)%I(1); IFK2=PJ%SKIJ(K)%J(1)
      IF(MODE>0)THEN
        DO J=NCSR%I(IFK1),NCSR%I(IFK1+1)-1
          IFK3=NCSR%J(J)
          KP=PJ%ID_SKIJ%BK(N)%A(IFK3-HL%SEC_N%ID(N)+1,IFK2-HL%SEC_N%ID(N)+1)
          IF(KP<=0)CYCLE
          NKKP%J(ISUM)=KP
          NKKP%A(ISUM)=NCSR%A(J)
          ISUM=ISUM+1
        ENDDO
      ELSE
        DO J=NCSR%I(IFK2),NCSR%I(IFK2+1)-1
          IFK3=NCSR%J(J)
          KP=PJ%ID_SKIJ%BK(N)%A(IFK1-HL%SEC_N%ID(N)+1,IFK3-HL%SEC_N%ID(N)+1)
          IF(KP<=0)CYCLE
          NKKP%J(ISUM)=KP
          NKKP%A(ISUM)=NCSR%A(J)
          ISUM=ISUM+1
        ENDDO
      ENDIF
      NKKP%I(K+1)=ISUM
      ENDDO; ENDDO
      !
      ! "CALC_NKKP_FOCK"
      RETURN
!
      END SUBROUTINE CALC_NKKP_FOCK
!
!****************************************************************************
! MODE > 0: Trace( S^\dagger U* VN U SP )
!     else: Trace( SP        U* VN U S^\dagger )
!****************************************************************************
      SUBROUTINE CALC_NKKP(NKKP,NCSR,HL,PJ,MODE)
      TYPE(ZCSR_MATRIX):: NKKP
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX) :: NCSR
      INTEGER MODE
! LOCAL
      INTEGER K,KP,I,J,IBK,IBASE,JBASE,ISUM,NROW,NCOL,N1,N2
      LOGICAL LDO
      TYPE(ZCSR_MATRIX) :: NCSR_
      COMPLEX(gq),ALLOCATABLE::UHVU(:,:),VU(:,:)
      !
!
      !
      IBK=1; ISUM=1; LDO=.TRUE.; IBASE=0
      NKKP%I(1)=1
      DO J=1,HL%SEC_J%DIM
      IF(HL%SEC_J%ID(J)>=HL%SEC_N%ID(IBK+1))THEN
        DO I=IBK+1,HL%EVEC%NBK
          IF(HL%SEC_N%ID(I+1)>HL%SEC_J%ID(J))THEN
            IBK=I; LDO=.TRUE.; EXIT
          ENDIF
          IBASE=IBASE+HL%EVEC%BK(I)%NROW
        ENDDO
      ENDIF 
      IF(LDO)THEN
        NROW=HL%EVEC%BK(IBK)%NROW; NCOL=HL%EVEC%BK(IBK)%NCOL
        ALLOCATE(VU(NROW,NCOL))
        N1=HL%SEC_N%ID(IBK); N2=HL%SEC_N%ID(IBK+1)-1
        CALL ZCSR_DIA_SUBMAT(NCSR,NCSR_,N1,N2,N1,N2)
        CALL ZCSRMUDEN_SK(NCSR_,HL%EVEC%BK(IBK)%A,VU,NCOL)
        CALL DEALLOC_ZCSR(NCSR_)
        IF(J>1) DEALLOCATE(UHVU)
        ALLOCATE(UHVU(NCOL,NCOL)); UHVU=0
        CALL ZGEMM('C','N',NCOL,NCOL,NROW,Z1,HL%EVEC%BK(IBK)%A,NROW,VU,NROW,Z0,UHVU,NCOL)
        DEALLOCATE(VU)
        IBASE=IBASE+NROW
        LDO=.FALSE.
      ENDIF
      JBASE=HL%SEC_N%ID(IBK)-1
      DO K =PJ%ID_PHIK_J(J),PJ%ID_PHIK_J(J+1)-1
      DO KP=PJ%ID_PHIK_J(J),PJ%ID_PHIK_J(J+1)-1
        CALL TR_DCOOZDNSDCOO(NKKP%A(ISUM),PJ%SKIJ(K),UHVU,PJ%SKIJ(KP),NCOL,JBASE,MODE)
        IF(ABS(NKKP%A(ISUM))<1.E-16_gq)CYCLE
        NKKP%J(ISUM)=KP
        ISUM=ISUM+1
      ENDDO ! KP
      NKKP%I(K+1)=ISUM
      ENDDO; ENDDO ! J
      DEALLOCATE(UHVU)
      !
      RETURN
!
      END SUBROUTINE CALC_NKKP
!
!****************************************************************************
      SUBROUTINE SET_MKKP(CO,FS,HL,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER IA1,IA2,MODE,IA1P,IA2P
      TYPE(ZBD_MATRIX) :: MKKP
      !
!
      ! "SET_MKKP"
      !
      ALLOCATE(PJ%MKKP(CO%DIM2,CO%DIM2)); PJ%MKKP(:,:)%DIM=-1
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15)THEN
        MODE=1
      ELSE
        MODE=0
      ENDIF
      DO IA1=1,CO%DIM2; DO IA2=1,CO%DIM2
      IA1P=CO%M_INDEX_R(1,IA1,IA2); IA2P=CO%M_INDEX_R(2,IA1,IA2)
      IF(IA1P==0)CYCLE
      IF(IA1P==IA1.AND.IA2P==IA2)THEN
        CALL ALLOC_ZBM(PJ%MKKP(IA1,IA2),HL%SEC_N%DIM,PJ%ID_PHIK_N,-1,MODE)
      ENDIF
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14.OR.HL%LGPRJ==15)THEN
        IF(IA1P==IA1.AND.IA2P==IA2)THEN
          CALL CALC_MKKP_FOCK(PJ%MKKP(IA1,IA2),IA1,IA2,HL,PJ,FS)
        ELSE
          CALL ALLOC_ZBM(MKKP,HL%SEC_N%DIM,PJ%ID_PHIK_N,-1,MODE)
          CALL CALC_MKKP_FOCK(MKKP,IA1,IA2,HL,PJ,FS)
          CALL ZBM_APLSB(PJ%MKKP(IA1P,IA2P),Z1,MKKP)
          CALL DEALLOC_ZBM(MKKP)
        ENDIF
      ELSE
        CALL CALC_MKKP(PJ%MKKP(IA1P,IA2P),IA1,IA2,HL,PJ,FS)
      ENDIF
      IF(HL%LGPRJ==4.OR.HL%LGPRJ==12)THEN
        CALL ZBMTOBKZCSR(PJ%MKKP(IA1,IA2))
      ENDIF
      ENDDO; ENDDO ! IA1,IA2
      !
      ! "SET_MKKP"
      RETURN
!
      END SUBROUTINE SET_MKKP
!
!****************************************************************************
! M_{IA1,IA2} =
! Trace( S(k)^\dagger f^\dagger_{IA2}  SP(k') f_{IA1} )
!       |G2><G1| f+_{IA2} |G3><G4| f_{IA1} |G2><G2|
!****************************************************************************
      SUBROUTINE CALC_MKKP_FOCK(MKKP,IA1,IA2,HL,PJ,FS)
      TYPE(ZBD_MATRIX) :: MKKP
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(FOCK_STATE)  FS
      INTEGER IA1,IA2
! LOCAL
      INTEGER N,K,KP,KO,INZ,BS1,BS2,BS3,BS4,SGN1,SGN2,NBASE,IFK3,IFK4
      REAL(gq) RES
      !
!
      ! "CALC_MKKP_FOCK"
      !
      NBASE=FS%DIM-HL%DIM
      DO N=2,HL%SEC_N%DIM
      MKKP%BK(N)%ACSR%I(1)=1; INZ=0
      DO K=PJ%ID_PHIK_N(N),PJ%ID_PHIK_N(N+1)-1
        KO=K-PJ%ID_PHIK_N(N)+1
        BS1=FS%BS(PJ%SKIJ(K)%I(1)+NBASE)
        BS2=FS%BS(PJ%SKIJ(K)%J(1)+NBASE)
        BS3=BS1; SGN1=1
        CALL A_STAT(BS3,IA2-1,.FALSE.,SGN1) ! <G1 | f+_{IA2} 
        IF(SGN1==0)GOTO 100
        IFK3=FS%IBS(BS3+1)-NBASE-HL%SEC_N%ID(N-1)+1
        IF(IFK3<=0)GOTO 100
        BS4=BS2; SGN2=1
        CALL A_STAT(BS4,IA1-1,.FALSE.,SGN2) ! f_{IA1} |G2>
        IF(SGN2==0)GOTO 100
        IFK4=FS%IBS(BS4+1)-NBASE-HL%SEC_N%ID(N-1)+1
        IF(IFK4<=0)GOTO 100
        KP=PJ%ID_SKIJ%BK(N-1)%A(IFK3,IFK4)-PJ%ID_PHIK_N(N-1)+1
        IF(KP<=0)GOTO 100
        RES=PJ%SKIJ(K)%A(1)*SGN1*SGN2*PJ%SKIJ(KP+PJ%ID_PHIK_N(N-1)-1)%A(1)
        INZ=INZ+1
        MKKP%BK(N)%ACSR%J(INZ)=KP
        MKKP%BK(N)%ACSR%A(INZ)=RES
100     MKKP%BK(N)%ACSR%I(KO+1)=INZ+1
      ENDDO ! K
      CALL ZCSR_COMPACT(MKKP%BK(N)%ACSR)
      ENDDO ! N
      !
      ! "CALC_MKKP_FOCK"
      RETURN
!
      END SUBROUTINE CALC_MKKP_FOCK
!
!****************************************************************************
! Trace( S(k)^\dagger f^\dagger SP(k') f )
!****************************************************************************
      SUBROUTINE CALC_MKKP(MKKP,IA1,IA2,HL,PJ,FS)
      TYPE(ZBD_MATRIX) :: MKKP
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(FOCK_STATE)  FS
      INTEGER IA1,IA2
! LOCAL
      INTEGER N,K,KP,NROWL,NROWR,NCOLL,NCOLR,I,J,IBASE,JBASE,ND_OCC
      COMPLEX(gq) ZES
      COMPLEX(gq),POINTER::CU(:,:),UHC1U(:,:),UHC2U(:,:)
      !
!
      ! "CALC_MKKP"
      !
      ND_OCC=HL%OCC(1)-FS%OCC(1)
      DO N=2,HL%SEC_N%DIM
      IF(MKKP%BK(N)%LZERO)CYCLE
      NROWL=FS%SEC_N%ID(N+ND_OCC  )-FS%SEC_N%ID(N-1+ND_OCC)
      NCOLL=HL%SEC_N%ID(N  )-HL%SEC_N%ID(N-1); IBASE=FS%SEC_N%ID(N-1+ND_OCC)-1
      NROWR=FS%SEC_N%ID(N+1+ND_OCC)-FS%SEC_N%ID(N+ND_OCC  )
      NCOLR=HL%SEC_N%ID(N+1)-HL%SEC_N%ID(N  ); JBASE=FS%SEC_N%ID(N+ND_OCC  )-1
      ALLOCATE(CU(NROWL,NCOLR)); CU=0
      DO K=HL%ID_FSC_N(N),HL%ID_FSC_N(N+1)-1
        I=FS%C(IA1)%I(K)-IBASE; J=FS%C(IA1)%J(K)-JBASE
        CU(I,:)=FS%C(IA1)%A(K)*HL%EVEC%BK(N)%A(J,:)
      ENDDO
      ALLOCATE(UHC1U(NCOLL,NCOLR)); UHC1U=0
      CALL ZGEMM('C','N',NCOLL,NCOLR,NROWL,Z1,HL%EVEC%BK(N-1)%A,NROWL,CU,NROWL,Z0,UHC1U,NCOLL)
      IF(IA2==IA1)THEN
        UHC2U=>UHC1U
      ELSE
        CU=0
        DO K=HL%ID_FSC_N(N),HL%ID_FSC_N(N+1)-1
          I=FS%C(IA2)%I(K)-IBASE; J=FS%C(IA2)%J(K)-JBASE
          CU(I,:)=FS%C(IA2)%A(K)*HL%EVEC%BK(N)%A(J,:)
        ENDDO
        ALLOCATE(UHC2U(NCOLL,NCOLR)); UHC2U=0
        CALL ZGEMM('C','N',NCOLL,NCOLR,NROWL,Z1,HL%EVEC%BK(N-1)%A,NROWL,CU,NROWL,Z0,UHC2U,NCOLL)
      ENDIF
      DEALLOCATE(CU)
      IBASE=HL%SEC_N%ID(N-1)-1; JBASE=HL%SEC_N%ID(N)-1
      DO K =PJ%ID_PHIK_N(N  ),PJ%ID_PHIK_N(N+1)-1; I=K -PJ%ID_PHIK_N(N  )+1
      DO KP=PJ%ID_PHIK_N(N-1),PJ%ID_PHIK_N(N  )-1; J=KP-PJ%ID_PHIK_N(N-1)+1
        CALL TR_DCOOZDNSDCOOZDNS(ZES,PJ%SKIJ(K),UHC2U,PJ%SKIJ(KP),UHC1U,NCOLL,NCOLR,IBASE,JBASE)
        MKKP%BK(N)%A(I,J)=MKKP%BK(N)%A(I,J)+ZES
      ENDDO; ENDDO
      IF(IA2==IA1)THEN
        NULLIFY(UHC2U)
      ELSE
        DEALLOCATE(UHC2U)
      ENDIF
      DEALLOCATE(UHC1U)
      ENDDO ! N
      !
      ! "CALC_MKKP"
      RETURN
!
      END SUBROUTINE CALC_MKKP
!
!****************************************************************************
! Symmetrize \phi_k or equivalently, skij
!****************************************************************************
      SUBROUTINE SYM_SKIJ(PJ,NORM)
      TYPE(PROJ) :: PJ
      REAL(gq) NORM
! LOCAL
      COMPLEX(gq),ALLOCATABLE::PHI(:,:)
      INTEGER I,J,NBASE,N1,MBASE,M1
!
      NBASE=PJ%SEC_N%ID(1)-1
      MBASE=PJ%SEC_N_SYM%ID(1)-1
      DO I=1,PJ%SEC_N%DIM
        N1=PJ%SEC_N%ID(I+1)-PJ%SEC_N%ID(I)
        IF(N1>0)THEN
          ALLOCATE(PHI(N1,N1))
          CALL SKIJ_TO_PHI(PJ,PJ%ID_PHIK_N(I),PJ%ID_PHIK_N(I+1),PHI,N1,NBASE)
          CALL UHAU(PHI,PJ%EVEC_SYM%BK(I)%A,N1,N1)
          CALL SYM_PHI_BK(PHI,N1,PJ%SKIJ_SYM,PJ%ID_PHIK_N_SYM(I),PJ%ID_PHIK_N_SYM(I+1),MBASE)
          CALL UHAU(PHI,PJ%EVEC_SYM%BK(I)%A,N1,N1,TRUL='N',TRUR='C')
          CALL SYM_PHI_BK(PHI,N1,PJ%SKIJ,PJ%ID_PHIK_N(I),PJ%ID_PHIK_N(I+1),NBASE,COEF=PJ%C)
          DEALLOCATE(PHI)
          NBASE=NBASE+N1
        ENDIF
        MBASE=MBASE+PJ%SEC_N_SYM%ID(I+1)-PJ%SEC_N_SYM%ID(I)
      ENDDO
!
      NORM=SQRT(DOT_PRODUCT(PJ%C,PJ%C))
      PJ%C=PJ%C/NORM
      RETURN
!      
      END SUBROUTINE SYM_SKIJ
!
!****************************************************************************
      SUBROUTINE SKIJ_TO_PHI(PJ,N1,N2,PHI,NP,NBASE)
      TYPE(PROJ) :: PJ
      INTEGER N1,N2,NP,NBASE
      COMPLEX(gq) PHI(NP,NP)
! LOCAL
      INTEGER I,J,I_,J_
!
      PHI=0
      DO I=N1,N2-1; DO J=1,PJ%SKIJ(I)%NNZ
        I_=PJ%SKIJ(I)%I(J)-NBASE
        J_=PJ%SKIJ(I)%J(J)-NBASE
        PHI(I_,J_)=PHI(I_,J_)+PJ%C(I)*PJ%SKIJ(I)%A(J)
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE SKIJ_TO_PHI
!
!****************************************************************************
      SUBROUTINE SYM_PHI_BK(PHI,NP,SKIJ,N1,N2,NBASE,COEF)
      INTEGER NP,N1,N2,NBASE
      COMPLEX(gq) PHI(NP,NP)
      TYPE(DCOO_MATRIX)::SKIJ(:)
      COMPLEX(gq),TARGET,OPTIONAL::COEF(:)
! LOCAL
      INTEGER MODE,I,J,I_,J_
      COMPLEX(gq),POINTER::C(:)
!
      IF(PRESENT(COEF))THEN
        MODE=1; C=>COEF
      ELSE
        MODE=0
        ALLOCATE(C(N2))
      ENDIF
      DO I=N1,N2-1
      C(I)=0
      DO J=1,SKIJ(I)%NNZ
        I_=SKIJ(I)%I(J)-NBASE; J_=SKIJ(I)%J(J)-NBASE
        C(I)=C(I)+SKIJ(I)%A(J)*PHI(I_,J_)
      ENDDO; ENDDO
      IF(MODE==1)THEN
        NULLIFY(C); RETURN
      ENDIF
      PHI=0
      DO I=N1,N2-1; DO J=1,SKIJ(I)%NNZ
        I_=SKIJ(I)%I(J)-NBASE; J_=SKIJ(I)%J(J)-NBASE
        PHI(I_,J_)=PHI(I_,J_)+SKIJ(I)%A(J)*C(I)
      ENDDO; ENDDO
      DEALLOCATE(C)
      RETURN
!
      END SUBROUTINE SYM_PHI_BK
!      
!****************************************************************************
! Contribution of N,M matrix
!****************************************************************************
      SUBROUTINE ACT_ZBMC(N,V1,V2,CO,HL,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      INTEGER N
      COMPLEX(gq) V1(N),V2(N)
! LOCAL
      INTEGER IA1,IA2,IBK,NROW,NCOL,I1,I2,J1,J2
      REAL(gq),PARAMETER::RTOL=1.E-10_gq
      COMPLEX(gq) ZES
      !
!
      !
      DO IA1=1,CO%DIM2; DO IA2=IA1,CO%DIM2
      IF(PJ%NVAR_KKP(IA1,IA2)%NROW<=0)CYCLE
      ZES=CO%LA2(IA1,IA2)
      IF(ABS(ZES)<RTOL)CYCLE
      CALL ZCSR_GESAMUX_SK('N',ZES,PJ%NVAR_KKP(IA1,IA2),V1,V2)
      IF(IA1/=IA2)THEN ! H.C.
        ZES=CONJG(ZES)
        CALL ZCSR_GESAMUX_SK('C',ZES,PJ%NVAR_KKP(IA1,IA2),V1,V2)
      ENDIF
      ENDDO; ENDDO
!
      DO IA1=1,CO%DIM2; DO IA2=1,CO%DIM2
      IF(PJ%MKKP(IA1,IA2)%DIM<=0)CYCLE
      ZES=CO%D(IA1,IA2)
      IF(ABS(ZES)<RTOL)CYCLE
      DO IBK=2,HL%SEC_N%DIM
        IF(PJ%MKKP(IA1,IA2)%BK(IBK)%LZERO)CYCLE
        I1=PJ%ID_PHIK_N(IBK); J1=PJ%ID_PHIK_N(IBK-1)
        NROW=PJ%ID_PHIK_N(IBK+1)-I1; NCOL=I1-J1
        I2=I1+NROW-1; J2=J1+NCOL-1
        IF(I1>I2.OR.J1>J2)CYCLE
        ZES=CO%D(IA1,IA2)
        IF(PJ%MKKP(IA1,IA2)%BK(IBK)%LSPARSE)THEN
          CALL ZCSR_GESAMUX_SK('N',ZES,PJ%MKKP(IA1,IA2)%BK(IBK)%ACSR,V1(J1:J2),V2(I1:I2))
          ZES=CONJG(ZES)
          CALL ZCSR_GESAMUX_SK('C',ZES,PJ%MKKP(IA1,IA2)%BK(IBK)%ACSR,V1(I1:I2),V2(J1:J2))
        ELSE
          CALL ZGEMV('N',NROW,NCOL,ZES,PJ%MKKP(IA1,IA2)%BK(IBK)%A,NROW,V1(J1:J2),1,Z1,V2(I1:I2),1) ! M
          ZES=CONJG(ZES)
          CALL ZGEMV('C',NROW,NCOL,ZES,PJ%MKKP(IA1,IA2)%BK(IBK)%A,NROW,V1(I1:I2),1,Z1,V2(J1:J2),1) ! M^H
        ENDIF
      ENDDO; ENDDO; ENDDO
      !
      RETURN
!
      END SUBROUTINE ACT_ZBMC
!
!****************************************************************************
      SUBROUTINE DUMP_PJ_C(PJ,NI)
      TYPE(PROJ)        PJ
      INTEGER NI
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,0,15))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      WRITE(GIU)PJ%C
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE DUMP_PJ_C
!
!****************************************************************************
      SUBROUTINE LOAD_PJ_C(PJ,NI)
      TYPE(PROJ)        PJ
      INTEGER NI
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,0,15))),STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      READ(GIU)PJ%C
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE LOAD_PJ_C
!
!****************************************************************************
      SUBROUTINE DUMP_PJ_C01(PJ,NI)
      TYPE(PROJ)        PJ
      INTEGER NI
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,0,1501))),STATUS='REPLACE')
      WRITE(GIU,'(F20.12)')PJ%C
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE DUMP_PJ_C01
!
!****************************************************************************
      SUBROUTINE LOAD_PJ_C01(PJ,NI)
      TYPE(PROJ)        PJ
      INTEGER NI
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,0,1501))),STATUS='OLD')
      READ(GIU,*)PJ%C
      CLOSE(GIU)
      RETURN
!
      END SUBROUTINE LOAD_PJ_C01
!
!****************************************************************************
      SUBROUTINE OPT_PHIK_PJENS_DS(PJ)
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER I,ITMAX,ITER
      REAL(gq) :: FTOL
      REAL(gq) :: Y(PJ%N_PHIK+1)
      COMPLEX(gq),ALLOCATABLE :: V(:,:)
      EXTERNAL::GUTZ_FCN_PJS
!
      ALLOCATE(V(PJ%N_PHIK+1,PJ%N_PHIK))
      DO I=1,PJ%N_PHIK+1
        V(I,:)=PJ%C
        IF(I>1)THEN
          V(I,I-1)=V(I,I-1)+1.E-3_gq
        ENDIF
        CALL GUTZ_FCN_PJS(V(I,:),PJ%N_PHIK,Y(I))
      ENDDO
      FTOL=1.E-8_gq; ITMAX=PJ%N_PHIK*1000
      CALL AMOEBA(V,Y,PJ%N_PHIK,FTOL,GUTZ_FCN_PJS,ITMAX,ITER)
      RETURN
!
      END SUBROUTINE OPT_PHIK_PJENS_DS
!
!
      END MODULE GPROJECTOR

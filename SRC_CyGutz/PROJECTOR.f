!**************************************************************************************************************************
! Copyright c 2013, The Ames Laboratory, Iowa State University, and Rutgers University*.  All rights reserved.
! 
! This software was authored by Yongxin Yao, Nicola Lanata*, Gabriel Kotliar*, Cai-Zhuang Wang, and Kai-Ming Ho, 
! at The Ames Laboratory and Rutgers University and was supported by the U.S. Department of Energy (DOE), Office of Science, 
! Basic Energy Sciences, Materials Science and Engineering Division.  
! The Ames Laboratory is operated by Iowa State University for DOE under U.S. Government contract DE-AC02-07CH11358.  
! The U.S. Government has the rights to use, reproduce, and distribute this software.  
! NEITHER THE GOVERNMENT, THE AMES LABORATORY, IOWA STATE UNIVERSITY, NOR RUTGERS UNIVERSITY MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
! If software is modified to produce derivative works, such modified software should be clearly marked, 
! so as to not confuse it with the version available from The Ames Laboratory and Rutgers University.
! 
! Additionally, redistribution and use in source and binary forms, with or without modification, 
! are permitted provided that the following conditions are met:
! 
!     Redistribution of source code must retain the above copyright notice, this list of conditions, 
! and the following disclaimer.
!     Redistribution in binary form must reproduce the above copyright notice, this list of conditions, 
! and the following disclaimer in the documentation and/or other materials provided with distribution.
!     Neither the name of The Ames Laboratory, Iowa State University, Rutgers University, the U.S. Government, 
! nor the names of its contributors may be used to endorse or promote products derived from this software 
! without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE AMES LABORATORY, IOWA STATE UNIVERSITY, RUTGERS UNIVERSITY, AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
! AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE GOVERNMENT, THE AMES LABORATORY, 
! IOWA STATE UNIVERSITY, RUTGERS UNIVERSITY, OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************************************************************************************************************
      MODULE GPROJECTOR
      USE gprec; USE SPARSE; USE FOCKSTATE; USE CORRORB; USE LOCALHAMIL; USE gconstant
      IMPLICIT NONE
!
! DEFINE TYPES
      TYPE PROJ
        LOGICAL         :: LMCFLY
        TYPE(ZCSR_MATRIX) :: FHL  ! F(U+M+N_var+N_phy)
        INTEGER :: N_PHIK,NEV_F,NCV_F,LEIGV
        INTEGER,POINTER :: ID_PHIK_N(:),ID_PHIK_J(:)
        TYPE(DCOO_MATRIX),ALLOCATABLE :: SKIJ(:)
        TYPE(IBD_MATRIX) :: ID_SKIJ ! Inverse table of SKIJ
        TYPE(ZCSR_MATRIX),ALLOCATABLE :: NVAR_KKP(:)
        TYPE(ZBD_MATRIX),ALLOCATABLE :: MKKP(:)
        TYPE(ZCSR_MATRIX) :: U_KKH
        REAL(gq)         :: EVL(2) ! Two lowest eigen values of HL
        REAL(gq)         :: EGAMM,EPOT2  ! ONSITE ENERGY/COULOMB POTENTIAL ENERGY
        COMPLEX(gq),POINTER :: C(:)   ! LOCAL CONFIGURATION PROBABILITY OR {c} in the notes
        TYPE(ZBD_MATRIX) :: RHO
        COMPLEX(gq),POINTER :: RHO0(:) ! Fock states mean-field probabilities
        REAL(gq),POINTER :: RHO_EV(:)
        REAL(gq) PJ_ENS ! Projection entropy
        REAL(gq),POINTER :: CW_N(:) ! Reduced configuration weight 
        INTEGER :: ITER
        INTEGER :: LENTANGLES
      END TYPE PROJ
!
      CONTAINS
!
!****************************************************************************
      SUBROUTINE SET_PROJ(CO,FS,HL,PJ,ITER_RHO)
      TYPE(CORR_ORB)    CO
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      INTEGER ITER_RHO
! LOCAL
      INTEGER MODE
!
      CALL MAP_HL_FS(FS,HL)
      IF((HL%LGPRJ/=21).AND.(HL%LGPRJ/=11).AND.(HL%LGPRJ/=14))THEN
        CALL ALLOC_ZBM(HL%EVEC,HL%SEC_N%DIM,HL%SEC_N%ID,0,0)
      ENDIF
!
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)THEN 
        HL%EVEC%NBK=-1; HL%EVEC%DIM=HL%DIM
        HL%SEC_J=HL%SEC_N
      ELSEIF(HL%LGPRJ==12)THEN
        CALL SET_IDEN_ZBM(HL%EVEC)
        HL%SEC_J=HL%SEC_N
      ELSE
        CALL SET_HL_SLJ2(CO,FS,HL)
        IF(HL%LGPRJ/=21)THEN
          CALL DIAG_HL_SLJ2(HL)
        ELSE
          CALL ZBM_LOAD(HL%EVEC,HL%NI,0,16)
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
      IF((HL%LGPRJ/=11).AND.(HL%LGPRJ/=14))CALL ZBM_DUMP(HL%EVEC,HL%NI,0,2)
      MODE=HL%LDIAPJ
      IF(HL%LGPRJ==14)MODE=2
      CALL SET_SKIJ(FS,HL,PJ,ITER_RHO,MODE)
      CALL SKIJ_DUMP(CO,PJ,HL,4)
      CALL SET_NVAR_KKP(CO,FS,HL,PJ)
      CALL ZCSR2_DUMP(PJ%NVAR_KKP,CO%SYMH%NR,HL%NI,5)
      CALL SET_MKKP(CO,FS,HL,PJ)
      CALL ZBM2_DUMP(PJ%MKKP,CO%SYMG%NR,HL%NI,6)
      CALL SET_HLOC(CO,FS,HL,.FALSE.)
!     CALL ZBM_WRT2(HL%H,6)
      IF((HL%LGPRJ/=11).AND.(HL%LGPRJ/=14))THEN
        CALL ZBM_UHAU(HL%U,HL%EVEC)
        CALL ZBM_SIMPLIFY(HL%U,HL%SEC_J%DIM,HL%SEC_J%ID)
      ENDIF
      CALL ZBM_DUMP(HL%U,HL%NI,0,1)
      IF(PJ%LMCFLY)THEN
        CALL DEALLOC_ZBM(HL%U)
      ENDIF
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)THEN
        CALL ZBMTOBKZCSR(HL%H)
        CALL SET_UKKH_FOCK(PJ%U_KKH,HL%H,HL,PJ)
      ELSE
        CALL ZBM_UHAU(HL%H,HL%EVEC)
        CALL ZBM_SIMPLIFY(HL%H,HL%SEC_J%DIM,HL%SEC_J%ID)
        CALL SET_UKKH(PJ%U_KKH,HL%H,HL,PJ)
      ENDIF
      CALL DEALLOC_ZBM(HL%H)
      IF(PJ%LMCFLY.AND.(HL%LGPRJ/=11).AND.(HL%LGPRJ/=14))THEN
        CALL DEALLOC_ZBM(HL%EVEC)
      ENDIF
      IF(ITER_RHO>0) DEALLOCATE(PJ%C)
      ALLOCATE(PJ%C(PJ%N_PHIK)) ! Allocate c_{k}
      PJ%C=0
      RETURN
!
      END SUBROUTINE SET_PROJ
!
!****************************************************************************
      SUBROUTINE SKIJ_DUMP(CO,PJ,HL,MODE)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      TYPE(LOCAL_HAMIL) HL
      INTEGER MODE
! LOCAL
      INTEGER I
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(HL%NI,0,MODE))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
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
      SUBROUTINE SKIJ_LOAD(CO,PJ,HL,MODE)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      TYPE(LOCAL_HAMIL) HL
      INTEGER MODE
! LOCAL
      INTEGER N,I
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(HL%NI,0,MODE))),STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
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
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)CALL IBM_READ(PJ%ID_SKIJ,GIU)
      CLOSE(GIU)
      !(*,'(" N_PHIK=",I10)')PJ%N_PHIK
      IF(.NOT.PJ%LMCFLY.AND.(PJ%N_PHIK>1000.OR.HL%DIM>1000))THEN
        PJ%LMCFLY=.TRUE.
      ENDIF
      RETURN
!
      END SUBROUTINE SKIJ_LOAD
!
!****************************************************************************
      SUBROUTINE CLR_PJ(CO,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER I
!
      DO I=1,CO%SYMG%NR
        CALL DEALLOC_ZBM(PJ%MKKP(I))
      ENDDO
      DO I=1,CO%SYMH%NR
        CALL DEALLOC_ZCSR(PJ%NVAR_KKP(I))
      ENDDO
      DEALLOCATE(PJ%MKKP,PJ%NVAR_KKP)
      RETURN
!
      END SUBROUTINE CLR_PJ
!
!****************************************************************************
      SUBROUTINE READ_PROJ(IO,CO,FS,HL,PJ,ITER_RHO,LSCF)
      INTEGER IO,ITER_RHO,LSCF
      TYPE(CORR_ORB)    CO
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      LOGICAL LSKIPU
!
      CALL MAP_HL_FS(FS,HL)
      CALL SKIJ_LOAD(CO,PJ,HL,4)
      WRITE(IO,'(" CO%ELC VS ELC_OLD:",2F8.3)')CO%ELC,CO%ELC_OLD
!
      ALLOCATE(PJ%NVAR_KKP(CO%SYMH%NR))
      CALL ZCSR2_LOAD(PJ%NVAR_KKP,CO%SYMH%NR,HL%NI,5)
      ALLOCATE(PJ%MKKP(CO%SYMG%NR))
      CALL ZBM2_LOAD(PJ%MKKP,CO%SYMG%NR,HL%NI,6)
!
      INQUIRE(FILE=TRIM(ADJUSTL(FILE_NAME(HL%NI,0,1))),EXIST=LSKIPU)
      CALL SET_HLOC(CO,FS,HL,LSKIPU)
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)THEN
        HL%EVEC%NBK=-1; HL%EVEC%DIM=HL%DIM
      ELSE
        CALL ZBM_LOAD(HL%EVEC,HL%NI,0,2)
      ENDIF
      IF(.NOT.LSKIPU)THEN
        IF(HL%LGPRJ/=11.OR.HL%LGPRJ/=14)THEN
          CALL ZBM_UHAU(HL%U,HL%EVEC)
          CALL ZBM_SIMPLIFY(HL%U,HL%SEC_J%DIM,HL%SEC_J%ID)
        ENDIF
        CALL ZBM_DUMP(HL%U,HL%NI,0,1)
        IF(PJ%LMCFLY) CALL DEALLOC_ZBM(HL%U)
      ENDIF
      IF((.NOT.PJ%LMCFLY).AND.LSKIPU)THEN
        CALL ZBM_LOAD(HL%U,HL%NI,0,1)
      ENDIF
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)THEN
        CALL ZBMTOBKZCSR(HL%H)
        CALL SET_UKKH_FOCK(PJ%U_KKH,HL%H,HL,PJ)
      ELSE
        CALL ZBM_UHAU(HL%H,HL%EVEC)
        CALL ZBM_SIMPLIFY(HL%H,HL%SEC_J%DIM,HL%SEC_J%ID)
        CALL SET_UKKH(PJ%U_KKH,HL%H,HL,PJ)
      ENDIF
      IF(LSCF==106)THEN
        CALL ZBM_DUMP(HL%H,HL%NI,0,8)
        CALL CALC_DUMP_HL_NAS(CO,HL,FS,HL%NI)
        CALL OUT_SEC(HL)
      ENDIF
      CALL DEALLOC_ZBM(HL%H)
!
      IF(PJ%LMCFLY.AND.(HL%LGPRJ/=11).AND.(HL%LGPRJ/=14))CALL DEALLOC_ZBM(HL%EVEC)
      IF(ITER_RHO>0)DEALLOCATE(PJ%C)
      ALLOCATE(PJ%C(PJ%N_PHIK)) ! Allocate c_{k}
      PJ%C=0
      RETURN
!
      END SUBROUTINE READ_PROJ
!
!****************************************************************************
      SUBROUTINE CALC_PJ_LC(CO,HL,PJ,IO,LSCF)
      TYPE(CORR_ORB)    CO
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      INTEGER IO,LSCF
! LOCAL
      INTEGER I
      LOGICAL :: LCINP
!
      IF(LSCF==105.OR.LSCF==106.OR.LSCF==107)THEN
        INQUIRE(FILE=TRIM(ADJUSTL(FILE_NAME(HL%NI,0,15))),EXIST=LCINP)
        IF(LCINP)THEN
          CALL LOAD_PJ_C(PJ,HL%NI)
          PJ%ITER=-1; PJ%EVL=0
          GOTO 100
        ENDIF
      ENDIF
      CALL CALC_FHL(CO,HL,PJ)
      CALL DIAG_FHL(PJ,IO)
      IF(LSCF.EQ.-11)THEN
        CALL LOAD_PJ_C01(PJ,HL%NI)
      ENDIF
      IF(PJ%LENTANGLES==10)THEN
        CALL OPT_PHIK_PJENS_DS(PJ)
      ENDIF
      IF(.NOT.PJ%LMCFLY)THEN
        CALL DEALLOC_ZCSR(PJ%FHL)
      ENDIF
100   CONTINUE
      CALL CALC_CZCSRC_COSYM(PJ%NVAR_KKP,CO%NC_VAR,CO,PJ,-1)
      RETURN
!
      END SUBROUTINE CALC_PJ_LC
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
      SUBROUTINE SET_SKIJ(FS,HL,PJ,ITER_RHO,MODE)
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      INTEGER ITER_RHO,MODE
! LOCAL
      INTEGER NPHIK
!
      ALLOCATE(PJ%ID_PHIK_N(HL%SEC_N%DIM+1))
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)THEN
        CALL GET_PHIK_DIM_N(HL%SEC_N,PJ%ID_PHIK_N,NPHIK,MODE,FS%BS_SZ(FS%DIM-HL%DIM+1:),HL%DIM)
        PJ%ID_PHIK_J=>PJ%ID_PHIK_N
      ELSE
        ALLOCATE(PJ%ID_PHIK_J(HL%SEC_J%DIM+1))
        CALL GET_PHIK_DIM_NJ(HL%SEC_N,HL%SEC_J,PJ%ID_PHIK_N,PJ%ID_PHIK_J,NPHIK,HL%LGPRJ,MODE)
      ENDIF
      ALLOCATE(PJ%SKIJ(NPHIK))
      PJ%N_PHIK=NPHIK
      !(*,'(" N_PHIK=",I10)')PJ%N_PHIK
      IF(.NOT.PJ%LMCFLY.AND.(PJ%N_PHIK>1000.OR.HL%DIM>1000))THEN
        PJ%LMCFLY=.TRUE.
      ENDIF
      IF(HL%LGPRJ==2.OR.HL%LGPRJ==4.OR.HL%LGPRJ==11.OR.HL%LGPRJ==12.OR.HL%LGPRJ==14)THEN
        CALL SET_SKIJ_N(PJ%SKIJ,HL%SEC_J,NPHIK,MODE,PJ%ID_SKIJ,FS%BS_SZ(FS%DIM-HL%DIM+1:))
      ELSE
        CALL SET_SKIJ_J(PJ%SKIJ,HL%SEC_J,NPHIK,MODE)
      ENDIF
      RETURN
!
      END SUBROUTINE SET_SKIJ
!
!****************************************************************************
      SUBROUTINE GET_PHIK_DIM_N(SEC_N,ID_PHIK_N,NPHIK,MODE,BS_SZ,NBS)
      TYPE(SECTOR1)::SEC_N
      INTEGER NPHIK,NBS,MODE,ID_PHIK_N(SEC_N%DIM+1),BS_SZ(NBS)
! LOCAL
      INTEGER I,N,IBASE,NPH1
!
      NPHIK=0; ID_PHIK_N(1)=1; IBASE=0
      DO I=1,SEC_N%DIM
        N=SEC_N%ID(I+1)-SEC_N%ID(I)
        IF(MODE==1)THEN
          NPHIK=NPHIK+N
        ELSEIF(MODE==2)THEN
          CALL GET_NPHIK_NBSZ(NPH1,BS_SZ(IBASE+1:IBASE+N),N)
          NPHIK=NPHIK+NPH1
        ELSE
          NPHIK=NPHIK+N*N
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
      SUBROUTINE SET_SKIJ_N(SKIJ,SEC_N,NPHIK,MODE,ID_SKIJ,BS_SZ)
      INTEGER NPHIK,MODE,BS_SZ(*)
      TYPE(DCOO_MATRIX) :: SKIJ(NPHIK)
      TYPE(IBD_MATRIX) :: ID_SKIJ ! Inverse table of SKIJ
      TYPE(SECTOR1)::SEC_N
! LOCAL
      INTEGER I,N,N1,N2,IK,IBASE,MULN,MN,NDIM,NSTR,NEND,SZ1,SZ2
      REAL(gq) RES
!
      IK=0; NDIM=SEC_N%ID(SEC_N%DIM+1)-1
      ID_SKIJ%NBK=SEC_N%DIM
      ALLOCATE(ID_SKIJ%BK(SEC_N%DIM))
      DO I=1,SEC_N%DIM
        IBASE=SEC_N%ID(I)-1
        N=SEC_N%ID(I+1)-SEC_N%ID(I)
        RES=1._gq; MULN=1
        ID_SKIJ%BK(I)%NROW=N; ID_SKIJ%BK(I)%NCOL=N
        ALLOCATE(ID_SKIJ%BK(I)%A(N,N)); ID_SKIJ%BK(I)%A=0
        DO N1=1,N
        IF(MODE==1)THEN  ! Diagonal
          NSTR=N1; NEND=N1
        ELSE             ! Full
          NSTR=1; NEND=N
          SZ1=BS_SZ(IBASE+N1)
        ENDIF
        DO N2=NSTR,NEND
          IF(MODE==2)THEN
            SZ2=BS_SZ(IBASE+N2)
            IF(SZ1/=SZ2)CYCLE
          ENDIF
          IK=IK+1
          CALL ALLOC_DCOO(SKIJ(IK),MULN,NDIM,NDIM)
          SKIJ(IK)%A=RES
          DO MN=1,MULN
            SKIJ(IK)%I(MN)=IBASE+N1+(MN-1)*N
            SKIJ(IK)%J(MN)=IBASE+N2+(MN-1)*N
          ENDDO
          ID_SKIJ%BK(I)%A(N1,N2)=IK
        ENDDO; ENDDO
      ENDDO
      ID_SKIJ%DIM=SUM( ID_SKIJ%BK(:)%NROW)
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
      SUBROUTINE GET_NPHIK_NBSZ(NPHIK,BS_SZ,N)
      INTEGER NPHIK,N,BS_SZ(N)
! LOCAL
      INTEGER I,IMIN,IMAX
      INTEGER,ALLOCATABLE::ICOUNT(:)
!
      NPHIK=0
      IMIN=MINVAL(BS_SZ); IMAX=MAXVAL(BS_SZ)
      ALLOCATE(ICOUNT(IMIN:IMAX)); ICOUNT=0
      DO I=1,N; ICOUNT(BS_SZ(I))=ICOUNT(BS_SZ(I))+1; ENDDO
      DO I=IMIN,IMAX; NPHIK=NPHIK+ICOUNT(I)**2; ENDDO
      RETURN
!
      END SUBROUTINE GET_NPHIK_NBSZ
!
!****************************************************************************
      SUBROUTINE CHK_COMMU_PHI_J(PJ,J)
      TYPE(PROJ)        PJ
      TYPE(ZBD_MATRIX) J
! LOCAL
      INTEGER I
!
      DO I=1,PJ%N_PHIK
        CALL COMMU_DCOOZBM1(PJ%SKIJ(I),J)
      ENDDO
      RETURN
!
      END SUBROUTINE CHK_COMMU_PHI_J
!
!****************************************************************************
! Only include U+N parts, M will be evaluated on fly
!****************************************************************************
      SUBROUTINE CALC_FHL(CO,HL,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      !
!
      IF(PJ%LMCFLY)RETURN
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
      INTEGER I
      COMPLEX(gq) ZES
      !
!
      !
      DO I=1,CO%SYMH%NR
        ZES=CO%LA2R(I)
        IF(ABS(ZES)<1.E-10_gq)CYCLE
        CALL ZCSR_APLSB_SK('N',A,ZES,PJ%NVAR_KKP(I))
        IF(CO%SYMH%IDP(1,1,I)/=CO%SYMH%IDP(2,1,I))THEN ! Off-diagonal
          ZES=CONJG(ZES)
          CALL ZCSR_APLSB_SK('C',A,ZES,PJ%NVAR_KKP(I))
        ENDIF
      ENDDO
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
      INTEGER I
      TYPE(ZCSR_MATRIX) MKKP
      !
!
      !
      DO I=1,CO%SYMG%NR
        CALL ZBMTOZCSR(PJ%MKKP(I),MKKP)
        CALL ZCSR_APLSB_SK('N',A,CO%DR(I),MKKP)
        CALL DEALLOC_ZCSR(MKKP)
      ENDDO
      !
      RETURN
!
      END SUBROUTINE CALC_M_FHL
!
!****************************************************************************
      SUBROUTINE DIAG_FHL(PJ,IO)
      TYPE(PROJ)       PJ
      INTEGER IO
! LOCAL
      INTEGER NDIM,NEV,NCV,NCONV
      COMPLEX(gq) SIGMA
      REAL(gq) ABSTOL
      REAL(gq),ALLOCATABLE::WE(:)
      COMPLEX(gq),ALLOCATABLE::ZE(:,:)
      EXTERNAL::PHIK_AV
!
      !"DIAG_FHL"
      SIGMA=0; ABSTOL=0
      NDIM=PJ%U_KKH%NROW; NEV=MIN(PJ%NEV_F,NDIM-1)
      IF(PJ%NCV_F==0)THEN
        NCV=MIN(NDIM,36*NEV)
      ELSE
        NCV=PJ%NCV_F
      ENDIF
      PJ%ITER=0
      IF(PJ%LMCFLY)THEN
        PJ%LEIGV=2
      ENDIF
      IF(PJ%LEIGV==0)THEN
        IF(NDIM.GT.100)THEN
          PJ%LEIGV=2
        ELSE
          PJ%LEIGV=1
        ENDIF
      ENDIF
      IF(PJ%LEIGV==2)THEN
        ALLOCATE(WE(NEV),ZE(NDIM,NEV)); WE=0; ZE=0
        CALL ZHDRV1(NDIM,NEV,NCV,'I','SR',ABSTOL,SIGMA,.TRUE.,WE,ZE,1000,NCONV,PHIK_AV)
        IF(NCONV.EQ.0)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(IO,'(" SEVERE WARNING: ZHDRV1--No converged Ritz values!")')
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(0 ,'(" SEVERE WARNING: ZHDRV1--No converged Ritz values!")')
        ELSEIF(NCONV.EQ.1)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(IO,'(" WARNING-ZHDRV1: Only 1 converged Ritz value!")')
        ENDIF
      ELSE
        ALLOCATE(ZE(NDIM,NDIM),WE(NDIM)); ZE=0; WE=0
        CALL SM_ZCSRZDNS(PJ%FHL,ZE)
        CALL ZHEEV_('V','L',ZE,WE,NDIM)
      ENDIF
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(IO,'(" PJ%C_OLD.PJ%C  =",2F12.8)')DOT_PRODUCT(PJ%C,ZE(:,1))
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(IO,'(" PJ%C VAL CHANGE=",2F12.8)')PJ%EVL-WE(1:2)
      PJ%C=ZE(:,1)
      PJ%EVL=WE(1:2)
      DEALLOCATE(WE,ZE)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(IO,'(" DIAG_FHL ITER=",I5," PJ%EVL=",2F20.8)')PJ%ITER,PJ%EVL
      IF(ABS(PJ%EVL(1)-PJ%EVL(2))<1.E-6_gq)THEN
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(IO,'(" WARNING: ALMOST DEGENERATE GROUND STATE!")')
      ENDIF
      !"DIAG_FHL"
      RETURN
!
      END SUBROUTINE DIAG_FHL
!
!****************************************************************************
! For CO%NC_VAR, etc
!****************************************************************************
      SUBROUTINE CALC_CZCSRC_COSYM(A,X,CO,PJ,MODE,N0R)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX)::A(CO%SYMH%NR)
      COMPLEX(gq)     ::X(CO%DIM2,CO%DIM2)
      INTEGER            MODE
      REAL(gq),OPTIONAL::N0R(CO%DIM2)
! LOCAL
      INTEGER DIMP,I
      COMPLEX(gq) ZES
      COMPLEX(gq) XR(CO%SYMH%NR)
!
      DO I=1,CO%SYMH%NR
        DIMP=CO%SYMH%DIMP(I)
        CALL CALC_CZCSRC(PJ%C,A(I),ZES,MODE)
        XR(I)=ZES/DIMP
        IF(PRESENT(N0R))XR(I)= XR(I)/SQRT(N0R(I)*(1-N0R(I)))
      ENDDO
      CALL SET_SYM_LOC_ARRAY(X,XR,CO,.TRUE.,1)
      RETURN
!
      END SUBROUTINE CALC_CZCSRC_COSYM
!
!****************************************************************************
! For CO%R0, R0_{ij} = <\phi | c_j+ f_i |\phi >
!****************************************************************************
      SUBROUTINE CALC_CZBMC_COSYM(A,X,CO,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
      TYPE(ZBD_MATRIX)::A(CO%SYMG%NR)
      COMPLEX(gq)    ::X(CO%DIM2,CO%DIM2)
! LOCAL
      INTEGER DIMP,I
      COMPLEX(gq) ZES
      COMPLEX(gq) XR(CO%SYMG%NR)
!
      DO I=1,CO%SYMG%NR
        DIMP=CO%SYMG%DIMP(I)
        CALL CALC_CZBMC(PJ%C,A(I),ZES,.FALSE.)
        XR(I)=ZES/DIMP
      ENDDO
      CALL SET_SYM_LOC_ARRAY(X,XR,CO,.TRUE.,0)
      RETURN
!
      END SUBROUTINE CALC_CZBMC_COSYM
!
!****************************************************************************
! CO%NC_PHY
!****************************************************************************
      SUBROUTINE CALC_NCPHY_COSYM(FS,HL,CO,PJ)
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(CORR_ORB)    CO
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER N1,I,N,J,NBASE,K,K_
      COMPLEX(gq) XR(CO%SYMH%NR)
      TYPE(ZCSR_MATRIX) :: NCSR
!
      DO I=1,CO%SYMH%NR
        N1=CO%SYMH%DIMP(I)
        CALL CALC_NIJ_SYMGP(NCSR,FS,FS%DIM-HL%DIM+1,FS%DIM,CO%SYMH%IDP(:,1:N1,I),N1)
        IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)THEN
          XR(I)=0; NBASE=0
          DO N=1,PJ%RHO%NBK
          DO J=1,PJ%RHO%BK(N)%NROW ! N_jk \Rho_kj
          DO K_=NCSR%I(J+NBASE),NCSR%I(J+NBASE+1)-1
            K=NCSR%J(K_)-NBASE
            XR(I)=XR(I)+NCSR%A(K_)*PJ%RHO%BK(N)%A(K,J)
          ENDDO; ENDDO
          NBASE=NBASE+PJ%RHO%BK(N)%NROW
          ENDDO
        ELSE
          CALL ZBM_TR_RHOUHVU(PJ%RHO,HL%EVEC,NCSR,XR(I))
        ENDIF
        CALL DEALLOC_ZCSR(NCSR)
        XR(I)=XR(I)/N1
      ENDDO
      CALL SET_SYM_LOC_ARRAY(CO%NC_PHY,XR,CO,.TRUE.,1)
      RETURN
!
      END SUBROUTINE CALC_NCPHY_COSYM
!
!****************************************************************************
      SUBROUTINE CALC_EGAMM(PJ)
      TYPE(PROJ)        PJ
! LOCAL
      COMPLEX(gq) ZES
!
      CALL CALC_CZCSRC(PJ%C,PJ%U_KKH,ZES,1)
      PJ%EGAMM=ZES
      RETURN
!
      END SUBROUTINE CALC_EGAMM
!
!****************************************************************************
      SUBROUTINE CALC_UPOT(HL,PJ)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      COMPLEX(gq) RES
!
      !"CALC_UPOT"
      CALL ZBM_TR_RHOA(PJ%RHO,HL%U,RES)
      PJ%EPOT2=RES
      !"CALC_UPOT"
      RETURN
!
      END SUBROUTINE CALC_UPOT
!
!****************************************************************************
      SUBROUTINE CALC_PJ_RHO(HL,PJ)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER NBK,IBK,IBASE,K,IJ,I,J
      INTEGER,POINTER::ID(:)
      COMPLEX(gq),EXTERNAL::ZDOTC
!
      NBK=HL%SEC_J%DIM; ID=>HL%SEC_J%ID
      IF(PJ%RHO%NBK<=0)THEN
        CALL ALLOC_ZBM(PJ%RHO,NBK,ID,0,0)
      ENDIF
      DO IBK=1,NBK
      IBASE=HL%SEC_J%ID(IBK)-1; PJ%RHO%BK(IBK)%A=0
      DO K=PJ%ID_PHIK_J(IBK),PJ%ID_PHIK_J(IBK+1)-1
      DO IJ=1,PJ%SKIJ(K)%NNZ
        I=PJ%SKIJ(K)%I(IJ)-IBASE; J=PJ%SKIJ(K)%J(IJ)-IBASE
        PJ%RHO%BK(IBK)%A(I,J)=PJ%RHO%BK(IBK)%A(I,J)+PJ%SKIJ(K)%A(IJ)*PJ%C(K)
      ENDDO; ENDDO; ENDDO
      CALL ZBM_GEMM('N','C',PJ%RHO,PJ%RHO) ! phi phi*
      NULLIFY(ID)
      RETURN
!
      END SUBROUTINE CALC_PJ_RHO
!
!****************************************************************************
      SUBROUTINE OUT_LOCAL_DBOCC1(PJ,FS,HL,IO,IA)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(FOCK_STATE)  FS
      INTEGER IO,IA
! LOCAL
      INTEGER I,IMIN,IMAX
      COMPLEX(gq) RES
      TYPE(ZBD_MATRIX) :: OP
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      IF(IA<=0)RETURN
      IF(IA<=FS%NA)THEN
        IMIN=IA; IMAX=IA
      ELSE
        IMIN=1; IMAX=FS%NA
      ENDIF
      CALL ALLOC_ZBM(OP,HL%SEF_N%DIM,HL%SEF_N%ID,0,0)
      DO I=IMIN,IMAX
        CALL CALC_LOCAL_OPR(FS,HL,OP,I,2)
        CALL ZBM_TR_RHOA(PJ%RHO,OP,RES)
        WRITE(IO,'(" I_ORBITAL=",I3," DOUBLE OCCUPANCY=",2F12.6)')I,RES
      ENDDO
      CALL DEALLOC_ZBM(OP)
      RETURN
!
      END SUBROUTINE OUT_LOCAL_DBOCC1
!
!****************************************************************************
! Tr(\Rho \sum_{\sigma}{c_{1\sigma}^{\dagger} c_{0 \sigma} - h.c.})
!****************************************************************************
      SUBROUTINE OUT_LOCAL_CURRENT1(PJ,FS,HL,IO,IA)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(FOCK_STATE)  FS
      INTEGER IO,IA
! LOCAL
      COMPLEX(gq) RES
      TYPE(ZBD_MATRIX) :: OP
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      IF(IA<=1)RETURN
      CALL ALLOC_ZBM(OP,HL%SEF_N%DIM,HL%SEF_N%ID,0,0)
      CALL CALC_LOCAL_OPR(FS,HL,OP,IA,1)
      CALL ZBM_TR_RHOA(PJ%RHO,OP,RES)
      CALL DEALLOC_ZBM(OP)
      WRITE(IO,'(" ExpValCurrentOpt=",2F12.6)')RES
      RETURN
!
      END SUBROUTINE OUT_LOCAL_CURRENT1
!
!****************************************************************************
! MODE=2: Single orbital double occupancy operator--n_{i up} n_{i dn}
!      1: Current operator: c_{1 }^{+} c_{0} -c_{0}^{+} c_{1}
!****************************************************************************
      SUBROUTINE CALC_LOCAL_OPR(FS,HL,OP,IA,MODE)
      TYPE(FOCK_STATE)   :: FS
      TYPE(LOCAL_HAMIL)  :: HL
      TYPE(ZBD_MATRIX) :: OP
      INTEGER IA,MODE
! LOCAL
      INTEGER NBK,ISEC,N1,NBASE,IFS,BS,BS1,SGN,ISP,IBS
      INTEGER,POINTER::ID(:)
!
      NBK=HL%SEF_N%DIM; ID=>HL%SEF_N%ID
      DO ISEC=1,NBK
      NBASE=ID(ISEC)-1+FS%DIM-HL%DIM
      N1=ID(ISEC+1)-ID(ISEC)
      OP%BK(ISEC)%A=0
      DO IFS=1,N1
        BS=FS%BS(NBASE+IFS)
        IF(MODE==2)THEN
          IF(.NOT.BTEST(BS,2*IA-2))CYCLE
          IF(.NOT.BTEST(BS,2*IA-1))CYCLE
          OP%BK(ISEC)%A(IFS,IFS)=1._gq
        ELSE
          DO ISP=1,2
            SGN=1; BS1=BS
            CALL A_STAT(BS1,2*IA-3+ISP,.FALSE.,SGN)
            IF(SGN==0)GOTO 101
            CALL A_STAT(BS1,2*IA-5+ISP,.TRUE.,SGN)
            IF(SGN==0)GOTO 101
            IBS=FS%IBS(BS1+1)-NBASE
            OP%BK(ISEC)%A(IBS,IFS)=OP%BK(ISEC)%A(IBS,IFS)+REAL(SGN,gq)
101         SGN=1; BS1=BS
            CALL A_STAT(BS1,2*IA-5+ISP,.FALSE.,SGN)
            IF(SGN==0)CYCLE
            CALL A_STAT(BS1,2*IA-3+ISP,.TRUE.,SGN)
            IF(SGN==0)CYCLE
            IBS=FS%IBS(BS1+1)-NBASE
            OP%BK(ISEC)%A(IBS,IFS)=OP%BK(ISEC)%A(IBS,IFS)-REAL(SGN,gq)
          ENDDO
        ENDIF
      ENDDO; ENDDO
      IF(HL%LGPRJ/=11.AND.HL%LGPRJ/=14)THEN
        IF(HL%EVEC%NBK<=0)THEN
          CALL ZBM_LOAD(HL%EVEC,HL%NI,0,2)
        ENDIF
        CALL ZBM_UHAU(OP,HL%EVEC)
        CALL ZBM_SIMPLIFY(OP,HL%SEC_J%DIM,HL%SEC_J%ID)
      ENDIF
      RETURN
      END SUBROUTINE CALC_LOCAL_OPR
!
!****************************************************************************
      SUBROUTINE CALC_RED_EVECNJ(HL,PJ,ECUT)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      REAL(gq) ECUT
! LOCAL
      INTEGER JBASE,I,J,NCOL,NCOL_SUM,ISECN
!
      ! 'CALC_RED_EVECNJ'
      JBASE=1; ISECN=1; NCOL_SUM=0
      DO I=1,PJ%RHO%NBK 
      DO J=JBASE,JBASE+PJ%RHO%BK(I)%NROW-1
        IF(PJ%RHO_EV(J)<RLBOUND)EXIT
        IF(-LOG(PJ%RHO_EV(J))<=ECUT)CYCLE
        IF(J>JBASE)THEN
         IF(-LOG(PJ%RHO_EV(J))+LOG(PJ%RHO_EV(J-1))<1.E-6_gq)CYCLE
        ENDIF
        EXIT
      ENDDO
      NCOL=J-JBASE
      PJ%RHO%BK(I)%NCOL=NCOL; NCOL_SUM=NCOL_SUM+NCOL
      HL%SEC_J%ID(I+1)=HL%SEC_J%ID(I)+NCOL
      JBASE=JBASE+PJ%RHO%BK(I)%NROW
      IF(JBASE>=HL%SEC_N%ID(ISECN+1))THEN
        HL%SEC_N%ID(ISECN+1)=HL%SEC_N%ID(ISECN)+NCOL_SUM
        ISECN=ISECN+1; NCOL_SUM=0
      ENDIF
      ENDDO
      ! 'CALC_RED_EVECNJ'
      RETURN
!
      END SUBROUTINE CALC_RED_EVECNJ
!
!****************************************************************************
! Calculate mean-field configuration (Fock states) probabilities
!****************************************************************************
      SUBROUTINE CALC_P0_FS(CO,FS,HL,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER NBASE,IFS,IA,BS
!
      NBASE=FS%DIM-HL%DIM
      IF(.NOT.ASSOCIATED(PJ%RHO0))THEN
        ALLOCATE(PJ%RHO0(HL%DIM))
      ENDIF
      PJ%RHO0=1._gq
      DO IFS=1,HL%DIM; BS=FS%BS(NBASE+IFS)
      DO IA=1,CO%DIM2
        IF(BTEST(BS,IA-1))THEN
          PJ%RHO0(IFS)=PJ%RHO0(IFS)*CO%NC_PHY(IA,IA)
        ELSE
          PJ%RHO0(IFS)=PJ%RHO0(IFS)*(1-CO%NC_PHY(IA,IA))
        ENDIF
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE CALC_P0_FS
!
!****************************************************************************
      SUBROUTINE SET_NVAR_KKP(CO,FS,HL,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER NR,I,N1
      TYPE(ZCSR_MATRIX) :: NCSR
      !
!
      ! "SET_NVAR_KKP"
      !
      NR=CO%SYMH%NR
      ALLOCATE(PJ%NVAR_KKP(NR))
      DO I=1,NR
        N1=CO%SYMH%DIMP(I)
        CALL CALC_NIJ_SYMGP(NCSR,FS,FS%DIM-HL%DIM+1,FS%DIM,CO%SYMH%IDP(:,1:N1,I),N1)
        CALL SET_NKKP(PJ%NVAR_KKP(I),NCSR,HL,PJ,-1,N1)
        CALL DEALLOC_ZCSR(NCSR)
      ENDDO
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
      TYPE(CORR_ORB)    CO
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX):: UKKH
      TYPE(ZBD_MATRIX) :: U
      !
!
      ! "SET_UKKH_FOCK"
      !
      CALL ALLOC_ZCSR_KKH_FOCK(UKKH,U,HL%LDIAPJ,PJ%N_PHIK)
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
      INTEGER I,J,K,KL,KP,N,INZ,IR,JR,NFK,JPHIK
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
! MODE = -1: DIAGONAL N_VAR    ; 1: DIAGONAL N_PHY
!        -2: OFF-DIAGONAL N_VAR; 2: OFF-DIAGONAL N_PHY
!****************************************************************************
      SUBROUTINE SET_NKKP(NKKP,NCSR,HL,PJ,MODE,NN)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(ZCSR_MATRIX) :: NKKP
      TYPE(ZCSR_MATRIX) :: NCSR
      INTEGER MODE,NN
! LOCAL
      !
!
      ! "SET_NKKP"
      !
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)THEN
        CALL ALLOC_ZCSR_KKP(NKKP,PJ%ID_PHIK_J,HL%SEC_J%DIM,NN,MODE)
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
      REAL(gq) RES
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
      INTEGER NR,I,IL,IA1,IA2,MODE
      !
!
      ! "SET_MKKP"
      !
      NR=CO%SYMG%NR
      ALLOCATE(PJ%MKKP(NR))
      DO I=1,NR
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)THEN
        MODE=CO%SYMG%DIMP(I)
      ELSE
        MODE=0
      ENDIF
      CALL ALLOC_ZBM(PJ%MKKP(I),HL%SEC_N%DIM,PJ%ID_PHIK_N,-1,MODE)
      IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)THEN
        CALL CALC_MKKP_FOCK(I,CO,HL,PJ,FS)
      ELSE
        DO IL=1,CO%SYMG%DIMP(I)
          IA1=CO%SYMG%IDP(1,IL,I); IA2=CO%SYMG%IDP(2,IL,I)
          CALL CALC_MKKP(PJ%MKKP(I),IA1,IA2,HL,PJ,FS)
        ENDDO
      ENDIF
      IF(HL%LGPRJ==4.OR.HL%LGPRJ==12)THEN
        CALL ZBMTOBKZCSR(PJ%MKKP(I))
      ENDIF
      ENDDO ! I
      !
      ! "SET_MKKP"
      RETURN
!
      END SUBROUTINE SET_MKKP
!
!****************************************************************************
! M_aA =
! Trace( S(k)^\dagger f^\dagger_A  SP(k') f_a )
!       |G2><G1| f+_A |G3><G4| f_a |G2><G2|
!****************************************************************************
      SUBROUTINE CALC_MKKP_FOCK(IR,CO,HL,PJ,FS)
      TYPE(CORR_ORB)    CO
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      TYPE(FOCK_STATE)  FS
      INTEGER IR
! LOCAL
      INTEGER N,K,KP,KO,INZ,IL,IA,BS1,BS2,BS3,BS4,SGN1,SGN2,NBASE,IFK3,IFK4
      REAL(gq) RES
      !
!
      ! "CALC_MKKP_FOCK"
      !
      NBASE=FS%DIM-HL%DIM
      DO N=2,HL%SEC_N%DIM
      PJ%MKKP(IR)%BK(N)%ACSR%I(1)=1; INZ=0
      DO K=PJ%ID_PHIK_N(N),PJ%ID_PHIK_N(N+1)-1
      KO=K-PJ%ID_PHIK_N(N)+1
      BS1=FS%BS(PJ%SKIJ(K)%I(1)+NBASE)
      BS2=FS%BS(PJ%SKIJ(K)%J(1)+NBASE)
      DO IL=1,CO%SYMG%DIMP(IR)
        IA=CO%SYMG%IDP(2,IL,IR) ! A
        BS3=BS1; SGN1=1
        CALL A_STAT(BS3,IA-1,.FALSE.,SGN1)
        IF(SGN1==0)CYCLE
        IFK3=FS%IBS(BS3+1)-NBASE-HL%SEC_N%ID(N-1)+1
        IF(IFK3<=0)CYCLE
        BS4=BS2; SGN2=1
        IA=CO%SYMG%IDP(1,IL,IR) ! a
        CALL A_STAT(BS4,IA-1,.FALSE.,SGN2)
        IF(SGN2==0)CYCLE
        IFK4=FS%IBS(BS4+1)-NBASE-HL%SEC_N%ID(N-1)+1
        IF(IFK4<=0)CYCLE
        KP=PJ%ID_SKIJ%BK(N-1)%A(IFK3,IFK4)-PJ%ID_PHIK_N(N-1)+1
        IF(KP<=0)CYCLE
        RES=PJ%SKIJ(K)%A(1)*SGN1*SGN2*PJ%SKIJ(KP+PJ%ID_PHIK_N(N-1)-1)%A(1)
        INZ=INZ+1
        PJ%MKKP(IR)%BK(N)%ACSR%J(INZ)=KP
        PJ%MKKP(IR)%BK(N)%ACSR%A(INZ)=RES
      ENDDO ! IL
      PJ%MKKP(IR)%BK(N)%ACSR%I(KO+1)=INZ+1
      ENDDO ! K
      CALL ZCSR_COMPACT(PJ%MKKP(IR)%BK(N)%ACSR)
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
! Contribution of N,M matrix
!****************************************************************************
      SUBROUTINE ACT_ZBMC(N,V1,V2,CO,HL,PJ)
      TYPE(CORR_ORB)    CO
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      INTEGER N
      COMPLEX(gq) V1(N),V2(N)
! LOCAL
      INTEGER I,IBK,NROW,NCOL,I1,I2,J1,J2
      REAL(gq) RES
      COMPLEX(gq) ZES
      !
!
      !
      DO I=1,CO%SYMH%NR
      ZES=CO%LA2R(I)
      IF(ABS(ZES)<1.E-10_gq)CYCLE
      CALL ZCSR_GESAMUX_SK('N',ZES,PJ%NVAR_KKP(I),V1,V2)
      IF(CO%SYMH%IDP(1,1,I)/=CO%SYMH%IDP(2,1,I))THEN ! H.C.
        ZES=CONJG(ZES)
        CALL ZCSR_GESAMUX_SK('C',ZES,PJ%NVAR_KKP(I),V1,V2)
      ENDIF
      ENDDO
!
      DO I=1,CO%SYMG%NR
      DO IBK=2,HL%SEC_N%DIM
        IF(PJ%MKKP(I)%BK(IBK)%LZERO)CYCLE
        I1=PJ%ID_PHIK_N(IBK); J1=PJ%ID_PHIK_N(IBK-1)
        NROW=PJ%ID_PHIK_N(IBK+1)-I1; NCOL=I1-J1
        I2=I1+NROW-1; J2=J1+NCOL-1
        ZES=CO%DR(I)
        IF(ABS(ZES)<1.E-10_gq)CYCLE
        IF(PJ%MKKP(I)%BK(IBK)%LSPARSE)THEN
          CALL ZCSR_GESAMUX_SK('N',ZES,PJ%MKKP(I)%BK(IBK)%ACSR,V1(J1:J2),V2(I1:I2))
          ZES=CONJG(ZES)
          CALL ZCSR_GESAMUX_SK('C',ZES,PJ%MKKP(I)%BK(IBK)%ACSR,V1(I1:I2),V2(J1:J2))
        ELSE
          CALL ZGEMV('N',NROW,NCOL,ZES,PJ%MKKP(I)%BK(IBK)%A,NROW,V1(J1:J2),1,Z1,V2(I1:I2),1) ! M
          ZES=CONJG(ZES)
          CALL ZGEMV('C',NROW,NCOL,ZES,PJ%MKKP(I)%BK(IBK)%A,NROW,V1(I1:I2),1,Z1,V2(J1:J2),1) ! M^H
        ENDIF
      ENDDO; ENDDO
      !
      RETURN
!
      END SUBROUTINE ACT_ZBMC
!
!****************************************************************************
      SUBROUTINE CALC_OCC_CONFIG_WT(HL,PJ)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER IBK_J,IBK_N,I
      REAL(gq),POINTER :: CW(:)
!
      ALLOCATE(PJ%CW_N(HL%SEC_N%DIM)); CW=>PJ%CW_N; CW=0
      IBK_N=1
      DO IBK_J=1,HL%SEC_J%DIM
      IF(HL%SEC_J%ID(IBK_J)>=HL%SEC_N%ID(IBK_N+1))THEN
        DO I=IBK_N+1,HL%SEC_N%DIM
          IF(HL%SEC_J%ID(IBK_J)<HL%SEC_N%ID(I+1))THEN
            EXIT
          ENDIF
        ENDDO
        IBK_N=I
      ENDIF
      DO I=1,PJ%RHO%BK(IBK_J)%NROW
        CW(IBK_N)=CW(IBK_N)+REAL(PJ%RHO%BK(IBK_J)%A(I,I),gq)
      ENDDO; ENDDO
      NULLIFY(CW)
      RETURN
!
      END SUBROUTINE CALC_OCC_CONFIG_WT
!
!****************************************************************************
      SUBROUTINE CALC_ENS(HL,PJ,ENS)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      REAL(gq) ENS
! LOCAL
      INTEGER IBK,NDIM,I
      COMPLEX(gq),ALLOCATABLE::LOGRHO(:,:)
!
      ENS=0
      DO IBK=1,HL%SEC_J%DIM
        NDIM=PJ%RHO%BK(IBK)%NROW
        ALLOCATE(LOGRHO(NDIM,NDIM)); LOGRHO=0
        CALL ZATOFA(PJ%RHO%BK(IBK)%A,LOGRHO,NDIM,2,1._gq,.TRUE.)
        DO I=1,NDIM
          ENS=ENS-SUM(PJ%RHO%BK(IBK)%A(I,:)*LOGRHO(:,I))
        ENDDO
        DEALLOCATE(LOGRHO)
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_ENS
!
!****************************************************************************
! Calculate projection entropy
!****************************************************************************
      SUBROUTINE CALC_PROJ_ENS(HL,PJ)
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
! LOCAL
      INTEGER IBK,NDIM,I,IND,NROW,NCOL,N1,N2
      COMPLEX(gq),ALLOCATABLE::LOGRHO(:,:),RHO0(:,:),V(:,:)
!
      !"CALC_PROJ_ENS"
      PJ%PJ_ENS=0; IND=0
      DO IBK=1,HL%SEC_J%DIM
        IF(HL%SEC_J%ID(IBK)>=HL%SEC_N%ID(IND+1))THEN
          DO I=IND+1,HL%SEC_N%DIM
            IF(HL%SEC_J%ID(IBK)<HL%SEC_N%ID(I+1))THEN
              EXIT
            ENDIF
          ENDDO
          IF(IND>0)THEN
            DEALLOCATE(RHO0)
          ENDIF
          IND=I
          NROW=HL%SEF_N%ID(IND+1)-HL%SEF_N%ID(IND)
          NCOL=HL%SEC_N%ID(IND+1)-HL%SEC_N%ID(IND)
          ALLOCATE(RHO0(NCOL,NCOL),V(NCOL,NROW)); RHO0=0
          IF(HL%LGPRJ==11.OR.HL%LGPRJ==14)THEN
            V=0; DO I=1,NROW; V(I,I)=1._gq; ENDDO
          ELSE
            V=CONJG(TRANSPOSE(HL%EVEC%BK(IND)%A))
          ENDIF
          N1=HL%SEF_N%ID(IND)
          CALL ZATOFA1(RHO0,PJ%RHO0(N1:N1+NROW-1),V,NROW,NCOL,-1,D1) ! 1/(P0*Tr[Rho])
          DEALLOCATE(V)
        ENDIF
        NDIM=PJ%RHO%BK(IBK)%NROW
        N1=HL%SEC_J%ID(IBK)-HL%SEC_N%ID(IND)+1; N2=N1+NDIM-1
        ALLOCATE(LOGRHO(NDIM,NDIM)); LOGRHO=0
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,Z1,PJ%RHO%BK(IBK)%A,NDIM,RHO0(N1:N2,N1:N2),NDIM,Z0,LOGRHO,NDIM)
        CALL ZATOFA(LOGRHO,LOGRHO,NDIM,2,1._gq,.FALSE.)
        DO I=1,NDIM
          PJ%PJ_ENS=PJ%PJ_ENS-SUM(PJ%RHO%BK(IBK)%A(I,:)*LOGRHO(:,I))
        ENDDO
        DEALLOCATE(LOGRHO)
      ENDDO
      !"CALC_PROJ_ENS"
      RETURN
!
      END SUBROUTINE CALC_PROJ_ENS
!
!****************************************************************************
      SUBROUTINE DUMP_PJ_C(PJ,NI)
      TYPE(PROJ)        PJ
      INTEGER NI
!
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,15))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
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
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,15))),STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
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
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,1501))),STATUS='REPLACE')
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
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(NI,0,1501))),STATUS='OLD')
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
!****************************************************************************
      SUBROUTINE CALC_SLJ_1(CO,HL,PJ)
      TYPE(FOCK_STATE)  FS
      TYPE(CORR_ORB)    CO
      TYPE(LOCAL_HAMIL) HL
      TYPE(PROJ)        PJ
      INTEGER IO
! LOCAL
      INTEGER LGPRJ,IBK,I
      COMPLEX(gq) RES
      LOGICAL LFILE
!
      LGPRJ=HL%LGPRJ
      IF(LGPRJ==1.OR.LGPRJ==2.OR.LGPRJ==4.OR.LGPRJ==21)THEN
        RES=0
        DO IBK=1,PJ%RHO%NBK; DO I=1,PJ%RHO%BK(IBK)%NROW
          RES=RES+REAL(PJ%RHO%BK(IBK)%A(I,I),gq)*HL%SEC_J%VAL(IBK)
        ENDDO; ENDDO
        IF(LGPRJ==1.OR.LGPRJ==2.OR.LGPRJ==21)THEN
          CO%J2=RES
        ELSEIF(LGPRJ==4)THEN
          CO%SZ=RES
        ENDIF
      ENDIF
      IF(LGPRJ==1.OR.LGPRJ==2.OR.LGPRJ==21)THEN
        CO%S2=-1._gq
        INQUIRE(FILE=TRIM(ADJUSTL(FILE_NAME(HL%NI,0,11))),EXIST=LFILE)
        IF(LFILE)THEN
          CALL ZBM_LOAD(HL%S2,HL%NI,0,11)
          CALL ZBM_UHAU(HL%S2,HL%EVEC)
          CALL ZBM_SIMPLIFY(HL%S2,HL%SEC_J%DIM,HL%SEC_J%ID)
          CALL ZBM_TR_RHOA(PJ%RHO,HL%S2,RES)
          CO%S2=RES
          CALL DEALLOC_ZBM(HL%S2)
        ENDIF
        CO%L2=-1._gq
        INQUIRE(FILE=TRIM(ADJUSTL(FILE_NAME(HL%NI,0,12))),EXIST=LFILE)
        IF(LFILE)THEN
          CALL ZBM_LOAD(HL%L2,HL%NI,0,12)
          CALL ZBM_UHAU(HL%L2,HL%EVEC)
          CALL ZBM_SIMPLIFY(HL%L2,HL%SEC_J%DIM,HL%SEC_J%ID)
          CALL ZBM_TR_RHOA(PJ%RHO,HL%L2,RES)
          CO%L2=RES
          CALL DEALLOC_ZBM(HL%L2)
        ENDIF
      ENDIF
      RETURN
!
      END SUBROUTINE CALC_SLJ_1
!
!****************************************************************************
      SUBROUTINE CALC_CZCSRC(C,A,ZES,MODE)
      TYPE(ZCSR_MATRIX):: A
      COMPLEX(gq) C(A%NROW),ZES
      INTEGER MODE
! LOCAL
      COMPLEX(gq),ALLOCATABLE::V1(:)
      COMPLEX(gq),EXTERNAL::ZDOTC
!
      ALLOCATE(V1(A%NROW)); V1=0
      IF(MODE>=0)THEN
        CALL ZCSR_SYAMUX_SK('L',A,C,V1)
      ELSE
        CALL ZCSR_GEAMUX_SK('N',A,C,V1)
      ENDIF
      ZES=ZDOTC(A%NROW,C,1,V1,1)
      DEALLOCATE(V1)
      RETURN
!
      END SUBROUTINE CALC_CZCSRC
!
!****************************************************************************
      SUBROUTINE CALC_CZBMC(C,A,ZES,LHM)
      TYPE(ZBD_MATRIX):: A
      COMPLEX(gq) C(A%DIM),ZES
      LOGICAL LHM ! Hermitian matrix
! LOCAL
      COMPLEX(gq),ALLOCATABLE::V1(:)
      COMPLEX(gq),EXTERNAL::ZDOTC
!
      ALLOCATE(V1(A%DIM)); V1=C
      CALL ZBM_AMUX(A,V1,LHM)
      ZES=ZDOTC(A%DIM,C,1,V1,1)
      DEALLOCATE(V1)
      RETURN
!
      END SUBROUTINE CALC_CZBMC
!
!
      END MODULE GPROJECTOR

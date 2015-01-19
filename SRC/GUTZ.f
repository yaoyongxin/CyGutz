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
      MODULE GUTZ
      USE gprec; USE gconstant; USE SPARSE; USE FOCKSTATE
      USE CORRORB; USE BANDSTRU; USE GPROJECTOR
      USE WAREHOUSE; USE GUTIL; USE GREENFUN; USE GMPI
      IMPLICIT NONE
!
      TYPE GL_YS
        TYPE (FOCK_STATE)  FS
        TYPE (CORR_ORB)    CO
        TYPE(LOCAL_HAMIL)  HL
        TYPE (PROJ)        PJ
      END TYPE GL_YS
!
      TYPE GL_YZ
        INTEGER NIMAP,NT,NI0,NIL
        TYPE (FOCK_STATE)  FS
        TYPE (CORR_ORB)    CO
        TYPE(LOCAL_HAMIL)  HL
        TYPE (PROJ)        PJ
      END TYPE GL_YZ
!
      TYPE GL_INFO
        INTEGER IO,IU
        INTEGER NI,ITER1,ITER2,ITER_RHO,ITER_GF,ITER_LA1,ITER_ETA
        INTEGER NASOMAX,NA2MAX,NAT_TOT,NEV_F,NCV_F
        INTEGER LDC,LSCF,LGREEN,LGPRJ,LDIAPJ,LHUB,LUNIT,LPSICOMP,LEBANDCOMP,LEL0,LBSCODE
! LBSCODE=0: eV and real Harmonics (VASP); 1: Ryd and complex Harmonics (Wien2K)
        INTEGER LENTANGLES,LSOLVER,LEIGV,LPJC,LBNDU
        INTEGER LMUSHIFT,LCHKLOC,LV2AO,LCLUSTER
        INTEGER LMODEL,LENSEMBLE,LPLTGF,LRROT,LIADBOCC,LCURRENT
        INTEGER NMAX_ITER,NMAX_MIX,LETA,LMCFLY,LXGUESS
        REAL(gq) DCMIX_A,RHO_CUT,RMIX_A,LA1SHIFT
      END TYPE GL_INFO
!
! VARIABLES
      TYPE (GL_INFO) ,SAVE :: GL
      TYPE (GL_YS),ALLOCATABLE,SAVE :: GYS(:)
      TYPE (GL_YZ),ALLOCATABLE,SAVE :: GYZ(:)
!
      CONTAINS
!
!********************************************************************************
! CALLED BY OUTSIDE DRIVER ROUTINES
!********************************************************************************
      SUBROUTINE GUTZ1_INI(NAT,IU,IO,LBSCODE)
      INTEGER IU,IO,NAT
      INTEGER,OPTIONAL::LBSCODE
! LOCAL
      INTEGER NIONS,NTYP,NAT_TOT
      INTEGER NT,NI,NA,NA2,NA2MAX,NA2TOT,NI0,NIMAP
      CHARACTER*256 STR
!
      !" GUTZ1_INI"
      GL%IO=IO; GL%IU=IU; GL%ITER_RHO=0; GMEM_SIZE=0
      GL%NAT_TOT=NAT
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" CyGutz VMPI.3.1.6.3 BUILT ON",A12," AT ",A8)')"Jan 12 2015","14:27:05"
      GL%LBSCODE=0 ! VASP
      GL%LSCF=6; GL%LRROT=0; GL%LBNDU=0
      GL%LENSEMBLE=0 ! Canonical ensemble
      GL%LEL0=0; BND%LFERWE=0; GL%RHO_CUT=0.99_gq; GL%RMIX_A=0.1_gq; GL%DCMIX_A=0.1_gq
      BND%ISPIN=1; GL%LGPRJ=1; GL%LHUB=1; GL%LDC=12; GL%NMAX_ITER=300; GL%NMAX_MIX=3
      GL%LPSICOMP=0; GL%LEBANDCOMP=0; GL%LENTANGLES=1; GL%LDIAPJ=0
      GL%LSOLVER=1; GL%NEV_F=10; GL%NCV_F=0; GL%LEIGV=0; GL%LPJC=0
      GL%LMUSHIFT=0; SYM%MODE=1; GL%LCHKLOC=0
      GL%LV2AO=-1 ! 1: complex Harmonics; 0: real Harmonics
      GL%LGREEN=0; GL%LMODEL=0 ! 0: Lattice model; Other: impurity model
      GF%WTYP=0 ! Complex frequency
      GF%DS=1._gq; GF%TSP=0.5_gq; GF%WMIN=0._gq; GF%WMAX=0._gq 
      GF%NW=1001; GF%TEV=1.E-8_gq; GF%T=1._gq ! Kelvin
      GF%ETA=1.E-12_gq; GL%LIADBOCC=0; GL%LCURRENT=0
      GL%LPLTGF=0; GL%LETA=0; GL%LMCFLY=.FALSE.
      GL%LA1SHIFT=0; GL%LCLUSTER=0; GL%LXGUESS=1
      GF%MU=0; BND%ISO=1
!
      IF(PRESENT(LBSCODE))GL%LBSCODE=LBSCODE
      IF(GL%LBSCODE==0)THEN
        GL%LUNIT=0 ! eV
      ELSE ! Wien2K
        GL%LUNIT=1 ! Ryd.
      ENDIF
      OPEN(IU,FILE='GL.INP',STATUS='OLD')
      DO 
        READ(IU,'(A256)',END=100)STR
        IF(INDEX(STR,'LSCF')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LSCF
        ELSEIF(INDEX(STR,'RHO_CUT')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%RHO_CUT
        ELSEIF(INDEX(STR,'LDIAPJ')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LDIAPJ
        ELSEIF(INDEX(STR,'LSOLVER')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LSOLVER
        ELSEIF(INDEX(STR,'LEIGV')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LEIGV
        ELSEIF(INDEX(STR,'ISPIN')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)BND%ISPIN
        ELSEIF(INDEX(STR,'LGPRJ')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LGPRJ
        ELSEIF(INDEX(STR,'LRROT')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LRROT
        ELSEIF(INDEX(STR,'LBNDU')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LBNDU
        ELSEIF(INDEX(STR,'LHUB')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LHUB
        ELSEIF(INDEX(STR,'BNDISO')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)BND%ISO
        ELSEIF(INDEX(STR,'LGREEN')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LGREEN
        ELSEIF(INDEX(STR,'LMCFLY')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LMCFLY
        ELSEIF(INDEX(STR,'LXGUESS')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LXGUESS
        ELSEIF(INDEX(STR,'LA1SHIFT')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LA1SHIFT
        ELSEIF(INDEX(STR,'LIADBOCC')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LIADBOCC
        ELSEIF(INDEX(STR,'LENSEMBLE')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LENSEMBLE
        ELSEIF(INDEX(STR,'LCLUSTER')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LCLUSTER
        ELSEIF(INDEX(STR,'GFWTYP')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GF%WTYP
        ELSEIF(INDEX(STR,'LCURRENT')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LCURRENT
        ELSEIF(INDEX(STR,'VBIAS')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GF%MU(1)
          GF%MU(2)=-GF%MU(1)/2; GF%MU(1)=GF%MU(1)/2
        ELSEIF(INDEX(STR,'GFNW')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GF%NW
          IF(MOD(GF%NW,2)==0)THEN
            GF%NW=GF%NW+1
          ENDIF
        ELSEIF(INDEX(STR,'GFWMIN')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GF%WMIN
        ELSEIF(INDEX(STR,'GFWMAX')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GF%WMAX
        ELSEIF(INDEX(STR,'GFDS')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GF%DS
        ELSEIF(INDEX(STR,'GFTSP')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GF%TSP
        ELSEIF(INDEX(STR,'LPLTGF')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LPLTGF
        ELSEIF(INDEX(STR,'GFTEMP')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GF%T
        ELSEIF(INDEX(STR,'GFTEV')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GF%TEV
        ELSEIF(INDEX(STR,'GFETA')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GF%ETA
        ELSEIF(INDEX(STR,'LMODEL')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LMODEL
        ELSEIF(INDEX(STR,'LETA')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LETA
        ELSEIF(INDEX(STR,'LPJC')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LPJC
        ELSEIF(INDEX(STR,'LEL0')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LEL0
        ELSEIF(INDEX(STR,'LV2AO')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LV2AO
        ELSEIF(INDEX(STR,'SYM_MODE')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)SYM%MODE
        ELSEIF(INDEX(STR,'LCHKLOC')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LCHKLOC
        ELSEIF(INDEX(STR,'LMUSHIFT')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LMUSHIFT
        ELSEIF(INDEX(STR,'LDC')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LDC
        ELSEIF(INDEX(STR,'RMIX_A')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%RMIX_A
        ELSEIF(INDEX(STR,'DCMIX_A')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%DCMIX_A
        ELSEIF(INDEX(STR,'NEV_F')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%NEV_F
        ELSEIF(INDEX(STR,'NCV_F')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%NCV_F
        ELSEIF(INDEX(STR,'NMAX_ITER')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%NMAX_ITER
        ELSEIF(INDEX(STR,'NMAX_MIX')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%NMAX_MIX
        ELSEIF(INDEX(STR,'LPSICOMP')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LPSICOMP
        ELSEIF(INDEX(STR,'LEBANDCOMP')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LEBANDCOMP
        ELSEIF(INDEX(STR,'LENTANGLES')>0)THEN
          READ(STR(INDEX(STR,'=')+1:),*)GL%LENTANGLES
        ELSEIF(INDEX(STR,'ATOM TYPE INFO')>0)THEN
          READ(IU,*)NTYP
          WH%NTYP=NTYP
          ALLOCATE(GYS(NTYP))
          DO NT=1,NTYP
            READ(IU,*)
            GYS(NT)%CO%U=0; GYS(NT)%CO%J=0
            GYS(NT)%CO%F0=0; GYS(NT)%CO%F2=0
            GYS(NT)%CO%F4=0; GYS(NT)%CO%F6=0
            SELECT CASE(GL%LHUB)
            CASE(1,2); READ(IU,*) GYS(NT)%CO%U,GYS(NT)%CO%J
            CASE(3  ); READ(IU,*) GYS(NT)%CO%F0,GYS(NT)%CO%F2,GYS(NT)%CO%F4,GYS(NT)%CO%F6
            CASE DEFAULT; READ(IU,*)
            END SELECT
            IF(GL%LUNIT==1)THEN ! Rydberg
              GYS(NT)%CO%U =GYS(NT)%CO%U /RYTOEV; GYS(NT)%CO%J =GYS(NT)%CO%J /RYTOEV
              GYS(NT)%CO%F0=GYS(NT)%CO%F0/RYTOEV; GYS(NT)%CO%F2=GYS(NT)%CO%F2/RYTOEV
              GYS(NT)%CO%F4=GYS(NT)%CO%F4/RYTOEV; GYS(NT)%CO%F6=GYS(NT)%CO%F6/RYTOEV
            ENDIF
            READ(IU,*) GYS(NT)%FS%NA
            GYS(NT)%CO%NELF1=-1._gq; GYS(NT)%CO%NELF2=-1._gq
            SELECT CASE(GL%LDC)
            CASE(2,12)
              READ(IU,*) GYS(NT)%HL%OCC,GYS(NT)%CO%NELF1
            CASE(3)
              READ(IU,*) GYS(NT)%HL%OCC,GYS(NT)%CO%NELF1,GYS(NT)%CO%NELF2
            CASE DEFAULT
              READ(IU,*) GYS(NT)%HL%OCC
            END SELECT
          ENDDO
        ELSEIF(INDEX(STR,'ATOM INFO')>0)THEN
          READ(IU,*)NIONS
          ALLOCATE(GYZ(NIONS))
          DO NI=1,NIONS
            READ(IU,*)
            READ(IU,*)NT,NI0,NIMAP ! Correlated atom type and index in the total atoms
            IF(NI0>GL%NAT_TOT)THEN
              WRITE(0,'(" TOTAL INEQUIVALENT ATOMS in STRUCT FILE=",I4," NI0=",I4)')GL%NAT_TOT,NI0
              STOP ' ERROR: NI0 IS THE INDEX OF INEQUIVALENT ATOM IN STRUCT FILE!'
            ENDIF
            GYZ(NI)%NT=NT; GYZ(NI)%NI0=NI0; GYZ(NI)%NIMAP=NIMAP
            NA=GYS(NT)%FS%NA; NA2=NA*2
            GYZ(NI)%FS%NA2 =NA2
            ALLOCATE(GYZ(NI)%CO%SYMG%ID(NA2,NA2))
            READ(IU,*)GYZ(NI)%CO%SYMG%ID
            GYZ(NI)%CO%SYMG%ID=TRANSPOSE(GYZ(NI)%CO%SYMG%ID) ! reading convention
          ENDDO
        ENDIF
      ENDDO
100   CLOSE(IU)
!
      IF(GL%LSCF==2)THEN
        GL%LEL0=1
      ENDIF
      IF(GL%LMODEL/=0.AND.GL%LGREEN==0)THEN
        GL%LGREEN=1
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" WARNING: GL%LMODEL/=0 -> SET GL%LGREEN=0")')
      ENDIF
      IF(GL%LV2AO==-1)THEN
        IF(GL%LHUB==2.OR.GL%LBSCODE==0)THEN ! Kanamori
          GL%LV2AO=0
        ELSE
          GL%LV2AO=1
        ENDIF
      ENDIF
      DO NT=1,NTYP
        NA=GYS(NT)%FS%NA; NA2=NA*2
        GYS(NT)%FS%NA2 =NA2
        GYS(NT)%CO%DIM =NA
        GYS(NT)%CO%DIM2=NA2
        GYS(NT)%FS%OCC(1)=MAX(GYS(NT)%HL%OCC(1)-1,0) ! No use C+C+CC (-2), instead C+C (-1)
        GYS(NT)%FS%OCC(2)=MIN(GYS(NT)%HL%OCC(2),NA2)
        CALL INI_FOCKSTATE(GYS(NT)%FS)
        CALL SET_U_MAT(GYS(NT)%CO,GL%LHUB,NT,GL%IU)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NT=",I3," FS%SEC_N%DIM=",I5," WITH ID:")')NT,GYS(NT)%FS%SEC_N%DIM
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(4X,10I8)')GYS(NT)%FS%SEC_N%ID
        IF(GL%LBSCODE==0)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" MEAN U/J:",2F8.4," eV.")')GYS(NT)%CO%UB,GYS(NT)%CO%JB
        ELSE
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" MEAN U/J:",2F8.4," Ryd.",2F8.4," eV.")')GYS(NT)%CO%UB,GYS(NT)%CO%JB,GYS(NT)%CO%UB*RYTOEV,GYS(NT)%CO%JB*RYTOEV
        ENDIF
        IF(GL%LDC.GE.2)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NEL_FIX_VDC1=",F10.2)')GYS(NT)%CO%NELF1
        ENDIF
        IF(GL%LDC.EQ.3)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NEL_FIX_EDC2=",F10.2)')GYS(NT)%CO%NELF2
        ENDIF
      ENDDO
      CALL DUMP_V2AO_ALL()
!
      NA2MAX=MAXVAL(GYS(:)%FS%NA2); GL%NA2MAX=NA2MAX
      NA2TOT=SUM(GYZ(:)%FS%NA2)
!
      CALL ALLOC_WAREHOUSE(NA2MAX,NIONS,NA2TOT)
      IF(GL%LBSCODE==0)THEN
        WH%B2N=>WH%R2N ! real Harmonics
      ELSE
        WH%B2N=>WH%C2N ! complex Harmonics
      ENDIF
      WH%DIMXG=0; WH%DIMXH=0; WH%DIMXHO=0
      DO NI=1,NIONS
        NT=GYZ(NI)%NT
        GYZ(NI)%CO%U =GYS(NT)%CO%U;  GYZ(NI)%CO%J =GYS(NT)%CO%J
        GYZ(NI)%CO%UB=GYS(NT)%CO%UB; GYZ(NI)%CO%JB=GYS(NT)%CO%JB
        NA2=GYS(NT)%FS%NA2
        GYZ(NI)%FS%NA  =GYS(NT)%FS%NA; GYZ(NI)%FS%NA2 =NA2
        GYZ(NI)%CO%DIM =GYS(NT)%FS%NA; GYZ(NI)%CO%DIM2=NA2
        GYZ(NI)%FS%DIM =GYS(NT)%FS%DIM
        GYZ(NI)%FS%OCC =GYS(NT)%FS%OCC; GYZ(NI)%HL%OCC =GYS(NT)%HL%OCC
        GYZ(NI)%CO%NELF1=GYS(NT)%CO%NELF1
        GYZ(NI)%CO%NELF2=GYS(NT)%CO%NELF2
        GYZ(NI)%FS%SEC_N=GYS(NT)%FS%SEC_N
        GYZ(NI)%FS%C    =>GYS(NT)%FS%C 
        GYZ(NI)%FS%BS   =>GYS(NT)%FS%BS
        GYZ(NI)%FS%BS_SZ=>GYS(NT)%FS%BS_SZ
        GYZ(NI)%FS%IBS  =>GYS(NT)%FS%IBS
        GYZ(NI)%FS%ID_FSC_N=>GYS(NT)%FS%ID_FSC_N
        GYZ(NI)%CO%NKS     =>WH%NKS   (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%ISIMIX  =>WH%ISIMIX(1:NA2,1:NA2,NI)
        GYZ(NI)%CO%NC_VAR  =>WH%NC_VAR(1:NA2,1:NA2,NI)
        GYZ(NI)%CO%NC_PHY  =>WH%NC_PHY(1:NA2,1:NA2,NI)
        GYZ(NI)%CO%R       =>WH%R     (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%R0      =>WH%R0    (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%C2N     =>WH%C2N   (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%R2N     =>WH%R2N   (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%B2N     =>WH%B2N   (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%N2N     =>WH%N2N   (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%Z       =>WH%Z     (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%D0      =>WH%D0    (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%D       =>WH%D     (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%LA1     =>WH%LA1   (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%ETA     =>WH%ETA   (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%LA2     =>WH%LA2   (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%EL0     =>WH%EL0   (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%VDC2    =>WH%VDC2  (1:NA2,1:NA2,NI)
        GYZ(NI)%CO%NPHY_FIX=>WH%NPHY_FIX(1:NA2,1:NA2,NI)
        ALLOCATE(GYZ(NI)%CO%V2H(1:NA2,1:NA2,1:NA2,1:NA2))
        GYZ(NI)%CO%V2H=GYS(NT)%CO%V2H
        GYZ(NI)%CO%SYMG%N=NA2
        CALL GROUP_SYM_ID(GYZ(NI)%CO%SYMG)
        CALL GET_ISYM_BK_MATRIX(GYZ(NI)%CO%SYMG%IDL,GYZ(NI)%CO%SYMG_BK,NA2)
        CALL OUT_CO_SYM(GYZ(NI)%CO%SYMG,GL%IO,NI,'SYMG')
        CALL REDUCE_SYM(GYZ(NI)%CO%SYMG,GYZ(NI)%CO%SYMH,1)
        CALL GROUP_SYM_ID(GYZ(NI)%CO%SYMH)
        CALL SYMH_TO_MATRIX_BASIS(GYZ(NI)%CO%SYMH)
        CALL OUT_CO_SYM(GYZ(NI)%CO%SYMH,GL%IO,NI,'SYMH')
        CALL REDUCE_SYM(GYZ(NI)%CO%SYMH,GYZ(NI)%CO%SYMHO,0)
        CALL OUT_CO_SYM(GYZ(NI)%CO%SYMHO,GL%IO,NI,'SYMHO')
        CALL INI_SYM_LOC_ARRAY(GYZ(NI)%CO)
        WH%DIMXG=MAX(WH%DIMXG,MAXVAL(GYZ(NI)%CO%SYMG%ID))
        WH%DIMXH=MAX(WH%DIMXH,MAXVAL(GYZ(NI)%CO%SYMH%ID))
        WH%DIMXHO=MAX(WH%DIMXHO,MAXVAL(GYZ(NI)%CO%SYMHO%ID))
        GYZ(NI)%PJ%LENTANGLES=GL%LENTANGLES
        GYZ(NI)%PJ%RHO%NBK=0
        GYZ(NI)%PJ%NEV_F=GL%NEV_F; GYZ(NI)%PJ%LEIGV=GL%LEIGV
        GYZ(NI)%PJ%NCV_F=GL%NCV_F
      ENDDO
!
      CALL READ_NELF1()
      CALL SET_NIL()
      CALL ORDER_COSYM_IDSUM()
      WH%DIMX1=WH%DIMXG; WH%DIMX2=WH%DIMX1+WH%DIMXH; WH%DIMX3=WH%DIMX2+WH%DIMXH
      WH%DIMXT=WH%DIMX3+WH%DIMXHO
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" DIM(X_Gutz)=",I3)')WH%DIMXT
      ALLOCATE(WH%X(WH%DIMXT),WH%XNC(WH%DIMXH),WH%XN0(WH%DIMXH),WH%LXREAL(WH%DIMXT))
      CALL GF_INI(GL%LGREEN,GL%LUNIT,GL%LMODEL,GL%IU,GL%IO)
      CALL OUT_PAR_GL(NIONS,NTYP,IO)
      CALL SET_WH_LXREAL()
      IF(SYM%MODE==10)THEN
        CALL READ_PAWSYM(GL%IU)
        CALL LOCMAP_SYM()
      ENDIF
      !"GUTZ1_INI"
      RETURN
!
      END SUBROUTINE GUTZ1_INI
!
!****************************************************************************
      SUBROUTINE GUTZ2_SET_NISO(ISO_IN,ISPIN_IN,NBMAX,NKPT)
      INTEGER ISO_IN,ISPIN_IN,NBMAX,NKPT
!
      !" GUTZ2_SET_NISO"
      CALL SET_BND_NISO(ISO_IN,ISPIN_IN,NBMAX,NKPT,GL%IO)
      GYS(:)%CO%DIMSO=GYS(:)%CO%DIM*BND%ISO
      GYZ(:)%CO%DIMSO=GYZ(:)%CO%DIM*BND%ISO
      GL%NASOMAX=MAXVAL(GYS(:)%CO%DIMSO)
      WH%NASOMAX=GL%NASOMAX
      BND%NASOTOT=SUM(GYZ(:)%CO%DIMSO)
      WH%NASOTOT =BND%NASOTOT
      !" GUTZ2_SET_NISO"
      RETURN
!
      END SUBROUTINE GUTZ2_SET_NISO
!
!******************************************************************************
      SUBROUTINE GUTZ4_SET_C2N_UH(B2N,NASOMAX,NIONS)
      INTEGER NASOMAX,NIONS
      COMPLEX(gq) B2N(NASOMAX,NASOMAX,NIONS)
! LOCAL
      INTEGER NI,NA2,NA1,NASO,I,L
      REAL(gq) MAXERR
      INTEGER,PARAMETER::LMAX=3
      COMPLEX(gq),POINTER :: C2N_U(:,:)
      COMPLEX(gq),ALLOCATABLE :: YLM_CR(:,:,:)
!
      !"GUTZ4_SET_C2N_UH"
! Catch: assume single l (e.g., d of f)
      ALLOCATE(YLM_CR(2*LMAX+1,2*LMAX+1,0:LMAX))
      CALL GET_YLM_CR_ALL(YLM_CR,LMAX)
      CALL SET_WH_N2N(GL%IU,GL%IO)
      WH%R2N=0; WH%C2N=0
      DO NI=1,WH%NIONS
        NA2=GYZ(NI)%CO%DIM2
        NA1=GYZ(NI)%CO%DIM; L=(NA1-1)/2
        NASO=NA1*BND%ISO_IN
        GYZ(NI)%CO%C2N=0
        IF(GL%LBSCODE==0)THEN ! REAL harmonics as basis -> complex
          GYZ(NI)%CO%R2N(1:NASO,1:NASO)=B2N(1:NASO,1:NASO,NI)
          IF(GL%LCLUSTER==0)THEN
            GYZ(NI)%CO%C2N(1:NA1,1:NA1)=MATMUL(YLM_CR(1:NA1,1:NA1,L),B2N(1:NA1,1:NA1,NI))
            IF(NASO==NA2)THEN
              GYZ(NI)%CO%C2N(1+NA1:NASO,1+NA1:NASO)=MATMUL(YLM_CR(1:NA1,1:NA1,L),B2N(1+NA1:NASO,1+NA1:NASO,NI))
            ENDIF
          ELSE
            GYZ(NI)%CO%C2N(1:NASO,1:NASO)=B2N(1:NASO,1:NASO,NI)
          ENDIF
        ELSE
          GYZ(NI)%CO%C2N(1:NASO,1:NASO)=B2N(1:NASO,1:NASO,NI)
          IF(GL%LCLUSTER==0)THEN
            GYZ(NI)%CO%R2N(1:NA1,1:NA1)=MATMUL(TRANSPOSE(CONJG(YLM_CR(1:NA1,1:NA1,L))),B2N(1:NA1,1:NA1,NI))
            IF(NASO==NA2)THEN
              GYZ(NI)%CO%R2N(1+NA1:NASO,1+NA1:NASO)=MATMUL(TRANSPOSE(CONJG(YLM_CR(1:NA1,1:NA1,L))),B2N(1+NA1:NASO,1+NA1:NASO,NI))
            ENDIF
          ELSE
            GYZ(NI)%CO%R2N(1:NASO,1:NASO)=B2N(1:NASO,1:NASO,NI)
          ENDIF
        ENDIF
        IF(NASO<NA2)THEN
          GYZ(NI)%CO%C2N(1+NASO:,1+NASO:)=GYZ(NI)%CO%C2N(1:NASO,1:NASO)
          GYZ(NI)%CO%R2N(1+NASO:,1+NASO:)=GYZ(NI)%CO%R2N(1:NASO,1:NASO)
        ENDIF
        CALL ZSPIN_BLK2U_TRANS(GYZ(NI)%CO%N2N,NA2,NA2,.TRUE.,BND%ISO,LURIGHT=.TRUE.) ! Additional transformation.
        GYZ(NI)%CO%C2N=MATMUL(GYZ(NI)%CO%C2N,GYZ(NI)%CO%N2N)
        GYZ(NI)%CO%R2N=MATMUL(GYZ(NI)%CO%R2N,GYZ(NI)%CO%N2N)
        IF(GL%LV2AO==0)THEN ! REAL HARMONICS
          C2N_U=>GYZ(NI)%CO%R2N
        ELSE
          C2N_U=>GYZ(NI)%CO%C2N
        ENDIF        
        CALL A4_TRANS_C(GYZ(NI)%CO%V2H,NA2,C2N_U)
        NULLIFY(C2N_U)
      ENDDO
      CALL ZOUT_MAT('U_N2N',WH%N2N,GL%IO)
      CALL ZOUT_MAT('U_C2N',WH%C2N,GL%IO)
      CALL ZOUT_MAT('U_B2N',WH%B2N,GL%IO)
      !"GUTZ4_SET_C2N_UH"
      RETURN
!
      END SUBROUTINE GUTZ4_SET_C2N_UH
!
!********************************************************************************
      SUBROUTINE GUTZ5_INI_CORR_UVEK(NKP,NKPL,WT,NELET,NE,KX,KY,KZ,KNAME)
      INTEGER NKP,NKPL
      REAL(gq) WT(NKP),NELET
      INTEGER,OPTIONAL :: NE(3,NKP)
      REAL(gq),OPTIONAL :: KX(NKP),KY(NKP),KZ(NKP)
      CHARACTER*10,OPTIONAL::KNAME(NKP)
! LOCAL
      INTEGER NI,NBASE,NASO
!
      !"GUTZ5_INI_CORR_UVEK"
      IF(GL%LMODEL==0)THEN
        ALLOCATE(BND%NE(3,NKP))
        IF(PRESENT(NE))THEN
          BND%NE=NE
        ELSE
          BND%NE(1,:)=BND%NMAX; BND%NE(2,:)=1; BND%NE(3,:)=BND%NMAX
        ENDIF
        BND%NE(1,:)= BND%NE(1,:)*BND%ISO/BND%ISO_IN
        BND%NE(3,:)= BND%NE(3,:)*BND%ISO/BND%ISO_IN
        BND%NE(2,:)=(BND%NE(2,:)-1)*BND%ISO/BND%ISO_IN+1
        BND%NMAXIN=MAXVAL(BND%NE(3,:)-BND%NE(2,:)+1)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" MIN/MAX(BND%NE(1,:))=",2I8)')MINVAL(BND%NE(1,:)),MAXVAL(BND%NE(1,:))
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" MIN/MAX(BND%NE(2,:))=",2I8)')MINVAL(BND%NE(2,:)),MAXVAL(BND%NE(2,:))
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" MIN/MAX(BND%NE(3,:))=",2I8)')MINVAL(BND%NE(3,:)),MAXVAL(BND%NE(3,:))
      ELSE
        BND%NMAXIN=BND%NASOTOT
      ENDIF
!
      IF(PRESENT(KX))THEN
        CALL SET_BND_UVEK(NKP,NKPL,WT,WH%NIONS,GL%LBNDU,KX=KX,KY=KY,KZ=KZ,KNAME=KNAME)
      ELSE
        CALL SET_BND_UVEK(NKP,NKPL,WT,WH%NIONS,GL%LBNDU)
      ENDIF
!
      IF(GL%LMODEL==0)THEN
        BND%NELEC=SUM((BND%NE(2,:)-1)*KPT%WT)*BND%RSPO
        BND%NELEL=NELET-BND%NELEC
      ELSE
        BND%NELEC=0; BND%NELEL=NELET
      ENDIF
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" VALENCE ELECTRONS: TOTAL=",F8.1)')NELET
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("                    CORRELATED BLOCK=",F8.1)')BND%NELEL
!
      GF%NBMAX=BND%NMAXIN; GF%NSPIN=BND%NSPIN
      IF(GL%LGREEN>0)THEN 
        ALLOCATE(GF%G(BND%NMAXIN ,BND%NMAXIN ,GF%NW,GF%NSPIN)); GF%G=0
        ALLOCATE(GF%M(BND%NASOTOT,BND%NASOTOT,GF%NW,GF%NSPIN)); GF%M=0
      ENDIF
      IF(GL%LMODEL==100.OR.GL%LMODEL==102)THEN
        CALL GF_HYBRD_INI() ! 1-band 
      ENDIF
!
      NBASE=1
      DO NI=1,WH%NIONS
        NASO=GYZ(NI)%CO%DIMSO
        IF(GL%LMODEL==0)THEN
          GYZ(NI)%CO%UK =>BND%UK (:,NBASE:NBASE+NASO-1,:,:)
          GYZ(NI)%CO%VK =>BND%VK (NBASE:NBASE+NASO-1,:,:,:,:)
          GYZ(NI)%CO%HK0=>BND%HK0(NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,:,:)
        ENDIF
        GYZ(NI)%CO%RB =>BND%R  (NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,:)
        GYZ(NI)%CO%DB =>BND%D0 (NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,:)
        GYZ(NI)%CO%LB1=>BND%LA1(NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,:)
        GYZ(NI)%CO%ETB=>BND%ETA(NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,:)
        IF(GL%LGREEN>0)THEN
          GYZ(NI)%CO%GF%T =GF%T; GYZ(NI)%CO%GF%NW=GF%NW
          GYZ(NI)%CO%GF%WMIN=GF%WMIN
          GYZ(NI)%CO%GF%WMAX=GF%WMAX
          GYZ(NI)%CO%GF%WTYP=GF%WTYP
          GYZ(NI)%CO%GF%MU  =GF%MU
          GYZ(NI)%CO%GF%W =>GF%W
          GYZ(NI)%CO%GF%G =>GF%G(NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,:,:)
          GYZ(NI)%CO%GF%M =>GF%M(NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,:,:)
          GYZ(NI)%CO%GF%NBMAX=NASO; GYZ(NI)%CO%GF%NSPIN=GF%NSPIN
          GYZ(NI)%CO%GF%NLR=GF%NLR
          IF(GL%LMODEL==100.OR.GL%LMODEL==102)THEN
            GYZ(NI)%CO%GF%D=>GF%D(NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,:,:,:)
            IF(GF%WTYP==1)THEN
              GYZ(NI)%CO%GF%GM=>GF%GM(NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,:,:,:)
              GYZ(NI)%CO%GF%WT=>GF%WT
            ENDIF
          ENDIF
        ENDIF
        NBASE=NBASE+NASO
      ENDDO
      BND%NELET=NELET
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SYM%IDI/IDF=",2I3)')SYM%IDI,SYM%IDF
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" BND%NMAXIN=",I5)')BND%NMAXIN
      !"GUTZ5_INI_CORR_UVEK"
      RETURN
!
      END SUBROUTINE GUTZ5_INI_CORR_UVEK
!
!****************************************************************************
      SUBROUTINE GUTZ6_SOLVE()
      INTEGER NI
      LOGICAL::LGLPRJ
      REAL :: TA1,TA2,TB1,TB2; INTEGER TIB1,TIB2,TIRATE
!
      !"GUTZ6_SOLVE"
      IF(GL%LGREEN>0)THEN
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" APPROXIMATED FERMI FUNCTION EWIN=",F8.3)')GF%EWIN
      ENDIF
      IF(GL%LMODEL==0)THEN
        CALL PRE_ANALYSIS()
        CALL ORTH_LOC_PROJ()
        CALL SETUP_TB(GL%LBNDU)
        CALL GET_EL0()
        CALL RM_TB_EL0() ! Remove local 1p part.
        CALL CALC_NKS(0)
        CALL CALC_DA0()
        CALL CALC_BAND_WIDGAP(BND%EK0,BND%NMAX,KPT%DIM,BND%EFLDA,GL%IO)
        CALL CALC_CORR_EBWIDTH(GL%IO)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" CORRELATED BAND WIDTH=",F8.3)')BND%EBWIDTH
        IF(GL%LGREEN>0.AND.BND%EBWIDTH>GF%EWIN)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" WARNING: MATSUBARA FREQUENCIES NOT SUFFICIENT!")')
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(0    ,'(" WARNING: MATSUBARA FREQUENCIES NOT SUFFICIENT!")')
        ENDIF
      ELSE
        CALL GET_EL0()
      ENDIF
      INQUIRE(FILE='GLMKKP_1_1.VSP',EXIST=LGLPRJ)
      IF(GL%LSCF==1.OR.GL%LSCF==3.OR.GL%LSCF==4.OR.GL%LSCF==5.OR.GL%LSCF==-11.OR.(.NOT.LGLPRJ.AND.GL%LSCF==6)) THEN
        CALL CPU_TIME(TA1); CALL SYSTEM_CLOCK(TIB1,TIRATE); TB1=REAL(TIB1,4)/REAL(TIRATE,4)
        CALL GUTZ_SET_PROJ()
        CALL CPU_TIME(TA2); CALL SYSTEM_CLOCK(TIB2,TIRATE); TB2=REAL(TIB2,4)/REAL(TIRATE,4)
        CALL OUT_TIME_USE('GUTZ_SET_PRJ',TA2-TA1,TB2-TB1,GL%IO)
      ELSEIF(GL%LSCF==6.OR.GL%LSCF==105.OR.GL%LSCF==106.OR.GL%LSCF==107)THEN
        CALL GUTZ_READ_PROJ()
      ELSEIF(GL%LSCF==2)THEN
      ELSE
        WRITE(0,'(" ERROR: ILLEGAL LSCF=",I2)')GL%LSCF; STOP
      ENDIF
      CALL OUT_SEC_ALL()
      CALL GMPI_BARRIER()
      IF(GL%LSCF==4.OR.GL%LSCF==5)STOP
      CALL CPU_TIME(TA1); CALL SYSTEM_CLOCK(TIB1,TIRATE); TB1=REAL(TIB1,4)/REAL(TIRATE,4)
      GL%ITER1=0; GL%ITER2=0
      IF(GL%LSOLVER==3.OR.GL%LSOLVER==101)THEN
        CALL  GUTZ_SOLVE_MIX()
      ELSE
        CALL GUTZ_SOLVE_HYBRD() ! HYBRD1
      ENDIF
      CALL CPU_TIME(TA2); CALL SYSTEM_CLOCK(TIB2,TIRATE); TB2=REAL(TIB2,4)/REAL(TIRATE,4)
      CALL OUT_TIME_USE('  GUTZ_SOLVE',TA2-TA1,TB2-TB1,GL%IO)
      CALL CHK_W_ETA_ALL()
      IF(GL%LPLTGF==1)THEN
        CALL PLOT_GF(GL%IU)
      ENDIF
      WH%EF=BND%EF
      CALL WRT_WH_RLNEF(GL%IU)
      CALL CALC_RNRL()
      IF(GL%LSCF==3)GOTO 100
      CALL WRT_PROJ_C()
      CALL CALC_PJRHO() ! LOCAL CONFIGURATION DENSITY MATRIX
      CALL CALC_NCPHY()
      CALL CALC_ENG()
      CALL CALC_RCW()
      CALL CALC_SLJ()
      IF(GL%LENTANGLES>=0)THEN
        CALL CALC_FSP0()  ! Fock states probabilities (diagonal)
        CALL CALC_PROJENS()
      ENDIF
      IF(GL%LEBANDCOMP>0.AND.GL%LMODEL==0)THEN
        CALL CALC_EBAND_COMP()
      ENDIF
      IF(GL%LENTANGLES>0)THEN
        CALL CALC_ENTANGLES()
      ENDIF
      CALL CALC_RED_VSP()
      IF(GL%LSCF==107)THEN
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SEC_N/J AFTER TRUNCATION BASED ON RHO:")')
        CALL OUT_SEC_ALL()
      ENDIF
100   CONTINUE
      CALL UPDATE_NELF1()
      IF(GL%LSOLVER==2)CALL UPDATE_MIX_R()
      WH%EF=BND%EF
      CALL WRT_WH_X(GL%IU)
      GL%ITER_RHO=GL%ITER_RHO+1
      CALL CLR_PROJ() ! RELEASE SOME MEMORY
      !"GUTZ6_SOLVE"
      RETURN
!
      END SUBROUTINE GUTZ6_SOLVE
!
!****************************************************************************
      SUBROUTINE GUTZ_SOLVE_HYBRD()
      INTEGER INFO,NX,NX1,NX2
      REAL(gq),ALLOCATABLE :: X(:),FVEC(:)
      REAL(gq),PARAMETER :: RTOL=1.E-10_gq,EPSFCN=1.E-10_gq
      EXTERNAL :: GUTZ_FCN1
!
      GL%ITER1=0
      IF(GL%ITER2<=1)THEN
        CALL INI_WH_X(GL%IU)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" INITIAL CHEMICAL POTENTIAL=",F10.4)')BND%EF
      ENDIF
      IF(GL%LSOLVER==2)THEN
        NX=SUM(WH%LXREAL(WH%DIMX1+1:WH%DIMX2))
        ALLOCATE(X(NX),FVEC(NX))
        CALL SET_XTOWHX(X,WH%X(WH%DIMX1+1:WH%DIMX2),WH%LXREAL(WH%DIMX1+1:WH%DIMX2),NX,WH%DIMXH,.TRUE.)
      ELSEIF(GL%LSOLVER==30)THEN
        NX1=SUM(WH%LXREAL(1:WH%DIMX1)); NX2=SUM(WH%LXREAL(1+WH%DIMX2:WH%DIMX3))
        NX=NX1+NX2
        ALLOCATE(X(NX),FVEC(NX))
        CALL SET_XTOWHX(X(1:NX1),WH%X(1:WH%DIMX1),WH%LXREAL(1:WH%DIMX1),NX1,WH%DIMXG,.TRUE.)
        CALL SET_XTOWHX(X(1+NX1:NX),WH%X(1+WH%DIMX2:WH%DIMX3),WH%LXREAL(1+WH%DIMX2:WH%DIMX3),NX2,WH%DIMXH,.TRUE.)
      ELSE
        NX=SUM(WH%LXREAL(1:WH%DIMX2))
        ALLOCATE(X(NX),FVEC(NX))
        CALL SET_XTOWHX(X,WH%X(1:WH%DIMX2),WH%LXREAL(1:WH%DIMX2),NX,WH%DIMX2,.TRUE.)
      ENDIF
      CALL HYBRD1(GUTZ_FCN1,NX,X,FVEC,RTOL,EPSFCN,INFO)
      DEALLOCATE(X,FVEC)
      SELECT CASE(INFO)
      CASE(0)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE1-HYBRD1: ERROR! Improper input parameters.")')
        WRITE(0     ,'(" SOLVE1-HYBRD1: ERROR! Improper input parameters.")')
      CASE(1)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE1-HYBRD1: Success.")')
      CASE(2)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE1-HYBRD1: WARNING! Number of calls to FCN has reached or exceeded 200*(N+1).")')
        WRITE(0     ,'(" SOLVE1-HYBRD1: WARNING! Number of calls to FCN has reached or exceeded 200*(N+1).")')
      CASE(3)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE1-HYBRD1: WARNING! TOL is too small. No further improvement in the approximate solution X is possible.")')
        WRITE(0     ,'(" SOLVE1-HYBRD1: WARNING! TOL is too small. No further improvement in the approximate solution X is possible.")')
      CASE(4)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE1-HYBRD1: WARNING! Iteration is not making good progress.")')
        WRITE(0     ,'(" SOLVE1-HYBRD1: WARNING! Iteration is not making good progress.")')
      END SELECT
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,  '(" SOLVE1 FINISHED WITH ",I7," ITERATIONS.")')GL%ITER1
      RETURN
!
      END SUBROUTINE GUTZ_SOLVE_HYBRD
!
!****************************************************************************
      SUBROUTINE GUTZ_SOLVE_LA1()
      INTEGER INFO,NX
      REAL(gq),ALLOCATABLE :: X(:),FVEC(:)
      REAL(gq),PARAMETER :: RTOL=1.E-10_gq,EPSFCN=1.E-10_gq
      EXTERNAL :: GUTZ_FCN_NK
!
      GL%ITER_LA1=0
      NX=SUM(WH%LXREAL(WH%DIMX1+1:WH%DIMX2))
      ALLOCATE(X(NX),FVEC(NX))
      CALL SET_XTOWHX(X,WH%X(WH%DIMX1+1:WH%DIMX2),WH%LXREAL(WH%DIMX1+1:WH%DIMX2),NX,WH%DIMXH,.TRUE.)
      CALL HYBRD1(GUTZ_FCN_NK,NX,X,FVEC,RTOL,EPSFCN,INFO)
      DEALLOCATE(X,FVEC)
      SELECT CASE(INFO)
      CASE(0)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE-LA1-HYBRD1: ERROR! Improper input parameters.")')
        WRITE(0     ,'(" SOLVE-LA1-HYBRD1: ERROR! Improper input parameters.")')
      CASE(1)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE-LA1-HYBRD1: Success.")')
      CASE(2)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE-LA1-HYBRD1: WARNING! Number of calls to FCN has reached or exceeded 200*(N+1).")')
        WRITE(0     ,'(" SOLVE-LA1-HYBRD1: WARNING! Number of calls to FCN has reached or exceeded 200*(N+1).")')
      CASE(3)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE-LA1-HYBRD1: WARNING! TOL is too small. No further improvement in the approximate solution X is possible.")')
        WRITE(0     ,'(" SOLVE-LA1-HYBRD1: WARNING! TOL is too small. No further improvement in the approximate solution X is possible.")')
      CASE(4)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE-LA1-HYBRD1: WARNING! Iteration is not making good progress.")')
        WRITE(0     ,'(" SOLVE-LA1-HYBRD1: WARNING! Iteration is not making good progress.")')
      END SELECT
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,  '(" SOLVE-LA1 FINISHED WITH ",I7," ITERATIONS.")')GL%ITER_LA1
      RETURN
!
      END SUBROUTINE GUTZ_SOLVE_LA1
!
!****************************************************************************
      SUBROUTINE GUTZ_SOLVE_ETA()
      INTEGER INFO,NX
      REAL(gq),ALLOCATABLE :: X(:),FVEC(:)
      REAL(gq),PARAMETER :: RTOL=1.E-10_gq,EPSFCN=1.E-10_gq
      EXTERNAL :: GUTZ_FCN_NKND
!
      GL%ITER_ETA=0
      NX=SUM(WH%LXREAL(WH%DIMX3+1:WH%DIMXT))
      ALLOCATE(X(NX),FVEC(NX)); X=0 ! ALWAYS START FROM 0.
      CALL HYBRD1(GUTZ_FCN_NKND,NX,X,FVEC,RTOL,EPSFCN,INFO)
      DEALLOCATE(X,FVEC)
      SELECT CASE(INFO)
      CASE(0)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE-ETA-HYBRD1: ERROR! Improper input parameters.")')
        WRITE(0     ,'(" SOLVE-ETA-HYBRD1: ERROR! Improper input parameters.")')
      CASE(1)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE-ETA-HYBRD1: Success.")')
      CASE(2)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE-ETA-HYBRD1: WARNING! Number of calls to FCN has reached or exceeded 200*(N+1).")')
        WRITE(0     ,'(" SOLVE-ETA-HYBRD1: WARNING! Number of calls to FCN has reached or exceeded 200*(N+1).")')
      CASE(3)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE-ETA-HYBRD1: WARNING! TOL is too small. No further improvement in the approximate solution X is possible.")')
        WRITE(0     ,'(" SOLVE-ETA-HYBRD1: WARNING! TOL is too small. No further improvement in the approximate solution X is possible.")')
      CASE(4)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE-ETA-HYBRD1: WARNING! Iteration is not making good progress.")')
        WRITE(0     ,'(" SOLVE-ETA-HYBRD1: WARNING! Iteration is not making good progress.")')
      END SELECT
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,  '(" SOLVE-ETA FINISHED WITH ",I7," ITERATIONS.")')GL%ITER_ETA
      RETURN
!
      END SUBROUTINE GUTZ_SOLVE_ETA
!
!****************************************************************************
      SUBROUTINE GUTZ_SOLVE_MIX()
      INTEGER INFO,NX,NX1,NX2
      REAL(gq),ALLOCATABLE :: X(:),FVEC(:)
      EXTERNAL :: GUTZ_FCN1
!
      GL%ITER1=0
      IF(GL%ITER2<=1)THEN
        CALL INI_WH_X(GL%IU)
      ENDIF
      NX1=SUM(WH%LXREAL(1:WH%DIMX1))
      IF(GL%LSOLVER==3)THEN
        NX2=SUM(WH%LXREAL(1+WH%DIMX2:WH%DIMX3))
      ELSEIF(GL%LSOLVER==101)THEN
        NX2=SUM(WH%LXREAL(1+WH%DIMX1:WH%DIMX2))
      ENDIF
      NX=NX1+NX2
      ALLOCATE(X(NX),FVEC(NX))
      CALL SET_XTOWHX(X(1:NX1),WH%X(1:WH%DIMX1),WH%LXREAL(1:WH%DIMX1),NX1,WH%DIMXG,.TRUE.)
      IF(GL%LSOLVER==3)THEN
        CALL SET_XTOWHX(X(1+NX1:NX),WH%X(1+WH%DIMX2:WH%DIMX3),WH%LXREAL(1+WH%DIMX2:WH%DIMX3),NX2,WH%DIMXH,.TRUE.)
      ELSEIF(GL%LSOLVER==101)THEN
        CALL SET_XTOWHX(X(1+NX1:NX),WH%X(1+WH%DIMX1:WH%DIMX2),WH%LXREAL(1+WH%DIMX1:WH%DIMX2),NX2,WH%DIMXH,.TRUE.)
      ENDIF
      CALL SIMPLE_MIX(GUTZ_FCN1,NX,X,FVEC,GL%RMIX_A,GL%NMAX_MIX,INFO)
      CALL SET_XTOWHX(X(1:NX1),WH%X(1:WH%DIMX1),WH%LXREAL(1:WH%DIMX1),NX1,WH%DIMXG,.FALSE.)
      IF(GL%LSOLVER==3)THEN
        CALL SET_XTOWHX(X(1+NX1:NX),WH%X(1+WH%DIMX2:WH%DIMX3),WH%LXREAL(1+WH%DIMX2:WH%DIMX3),NX2,WH%DIMXH,.FALSE.)
      ELSEIF(GL%LSOLVER==101)THEN
        CALL SET_XTOWHX(X(1+NX1:NX),WH%X(1+WH%DIMX1:WH%DIMX2),WH%LXREAL(1+WH%DIMX1:WH%DIMX2),NX2,WH%DIMXH,.FALSE.)
      ENDIF
      SELECT CASE(INFO)
      CASE(0)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE1-MIX: Sucess.")')
      CASE(1)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SOLVE1-MIX: Not yet converge.")')
        WRITE(0     ,'(" SOLVE1-MIX: Not yet converge.")')
      END SELECT
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,  '(" SOLVE1-MIX FINISHED WITH ",I7," ITERATIONS.")')GL%ITER1
      RETURN
!
      END SUBROUTINE GUTZ_SOLVE_MIX
!
!****************************************************************************
! Naive way to get Ef, not efficient at all.
!****************************************************************************
      SUBROUTINE GREEN_FERMI_NV()
      INTEGER INFO
      REAL(gq) X(1),FVEC(1)
      REAL(gq),PARAMETER :: RTOL=1.E-10_gq,EPSFCN=1.E-10_gq
      EXTERNAL :: GFERMI_FCN
!
      GL%ITER_GF=0
      X(1)=BND%EF
      CALL HYBRD1(GFERMI_FCN,1,X,FVEC,RTOL,EPSFCN,INFO)
      SELECT CASE(INFO)
      CASE(0)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GREEN_FERMI_NV: ERROR! Improper input parameters.")')
        WRITE(0     ,'(" GREEN_FERMI_NV: ERROR! Improper input parameters.")')
      CASE(1)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GREEN_FERMI_NV: Success.")')
      CASE(2)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GREEN_FERMI_NV: WARNING! Number of calls to FCN has reached or exceeded 200*(N+1).")')
        WRITE(0     ,'(" GREEN_FERMI_NV: WARNING! Number of calls to FCN has reached or exceeded 200*(N+1).")')
      CASE(3)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GREEN_FERMI_NV: WARNING! TOL is too small. No further improvement in the approximate solution X is possible.")')
        WRITE(0     ,'(" GREEN_FERMI_NV: WARNING! TOL is too small. No further improvement in the approximate solution X is possible.")')
      CASE(4)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GREEN_FERMI_NV: WARNING! Iteration is not making good progress.")')
        WRITE(0     ,'(" GREEN_FERMI_NV: WARNING! Iteration is not making good progress.")')
      END SELECT
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,  '(" GREEN_FERMI_NV FINISHED WITH ",I7," ITERATIONS.")')GL%ITER_GF
      BND%EF=X(1)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GUTZ_FERMI=",F12.5)')BND%EF
      RETURN
!
      END SUBROUTINE GREEN_FERMI_NV
!
!****************************************************************************
! INTERNAL ROUTINES
!****************************************************************************
      SUBROUTINE CALC_VDC()
      INTEGER NI
      REAL(gq) UB,JB,NT
!
      IF(GL%LDC.EQ.0)RETURN
      IF(GL%LDC<0)THEN
        IF(GL%LDC==-1.OR.ABS(WH%NPHY_FIX(1,1,1))>1.1_gq)THEN
          WH%NPHY_FIX=WH%NKS
        ENDIF
      ENDIF
      DO NI=1,WH%NIONS
        CALL CALC_CO_VDC(GYZ(NI)%CO,GL%LDC)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3," VDC=",F10.4)')NI,GYZ(NI)%CO%VDC
      ENDDO
      IF(GL%LDC<0)THEN
        CALL ZOUT_MAT('VDC2',WH%VDC2,GL%IO)
      ENDIF
      RETURN
!
      END SUBROUTINE CALC_VDC
!
!****************************************************************************
      SUBROUTINE CALC_LC()
      INTEGER NI,NIMAP
      REAL :: TA1,TA2,TB1,TB2; INTEGER TIB1,TIB2,TIRATE
!
      !"CALC_LC"
      CALL CPU_TIME(TA1); CALL SYSTEM_CLOCK(TIB1,TIRATE); TB1=REAL(TIB1,4)/REAL(TIRATE,4)
      WH%NC_VAR=0
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.EQ.NI)THEN
          GL%NI=NI; GYZ(NI)%HL%NI=NI
          CALL CALC_PJ_LC(GYZ(NI)%CO,GYZ(NI)%HL,GYZ(NI)%PJ,GL%IO,GL%LSCF)
          WRITE(GL%IO,'(" NI=",I3," PHI-step gap=",E10.2)')NI,GYZ(NI)%PJ%EVL(2)-GYZ(NI)%PJ%EVL(1)
          IF(GYZ(NI)%PJ%ITER==-1)THEN
            GL%ITER1=-1
          ENDIF
        ELSE
          GYZ(NI)%CO%NC_VAR=GYZ(NIMAP)%CO%NC_VAR
        ENDIF
      ENDDO
      CALL ZSUM_ALL_MPI(WH%NC_VAR,WH%NA2TOT)
      CALL ZOUT_MAT('NC_VAR',WH%NC_VAR,GL%IO,1)
      CALL CPU_TIME(TA2); CALL SYSTEM_CLOCK(TIB2,TIRATE); TB2=REAL(TIB2,4)/REAL(TIRATE,4)
      CALL OUT_TIME_USE('CALC_LC',TA2-TA1,TB2-TB1,GL%IO)
      !"CALC_LC"
      RETURN
!
      END SUBROUTINE CALC_LC
!
!****************************************************************************
      SUBROUTINE REGULARIZE_NKS(MODE)
      INTEGER MODE
! LOCAL
      INTEGER NI,NIMAP
!
      DO NI=1,WH%NIONS
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.EQ.NI)THEN
          IF(MODE==0)THEN
            CALL ZH_REGULARIZE(GYZ(NI)%CO%NKS   ,GYZ(NI)%CO%DIM2,1.E-10_gq)
          ELSE
            CALL ZH_REGULARIZE(GYZ(NI)%CO%NC_VAR,GYZ(NI)%CO%DIM2,1.E-10_gq)
          ENDIF
        ELSE
          IF(MODE==0)THEN
            GYZ(NI)%CO%NKS   =GYZ(NIMAP)%CO%NKS
          ELSE
            GYZ(NI)%CO%NC_VAR=GYZ(NIMAP)%CO%NC_VAR
          ENDIF
        ENDIF
      ENDDO
      IF(MODE==0)THEN
        CALL ZOUT_MAT('NKS-REGU',WH%NKS,GL%IO,1)
      ELSE
        CALL ZOUT_MAT('NC_VAR-REG',WH%NKS,GL%IO,1)
      ENDIF
      RETURN
!
      END SUBROUTINE REGULARIZE_NKS
!
!****************************************************************************
      SUBROUTINE CALC_ISIMIX(MODE)
      INTEGER MODE
! LOCAL
      INTEGER NI,NIMAP
!
      DO NI=1,WH%NIONS
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.EQ.NI)THEN
          IF(MODE==0)THEN
            CALL CALC_CO_ISIMIX(GYZ(NI)%CO,GYZ(NI)%CO%NKS   )
          ELSE
            CALL CALC_CO_ISIMIX(GYZ(NI)%CO,GYZ(NI)%CO%NC_VAR)
          ENDIF
        ELSE
          GYZ(NI)%CO%ISIMIX=GYZ(NIMAP)%CO%ISIMIX
        ENDIF
      ENDDO
      CALL ZOUT_MAT('ISIMIX-VAR',WH%ISIMIX,GL%IO,1)
      RETURN
!
      END SUBROUTINE CALC_ISIMIX
!
!****************************************************************************
      SUBROUTINE UPDATE_MIX_R()
      REAL(gq) :: R_OLD(WH%NA2MAX,WH%NA2MAX,WH%NIONS),R_DIF(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
!
      R_OLD=WH%R
      CALL UPDATE_R_ALL()
      R_DIF=WH%R-R_OLD
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" MAX R DIFF=",F14.8)')MAXVAL(ABS(R_DIF))
      WH%R=R_OLD+GL%RMIX_A*R_DIF
      CALL SET_WH_X123(.TRUE.)
      RETURN
!
      END SUBROUTINE UPDATE_MIX_R
!
!****************************************************************************
      SUBROUTINE UPDATE_R_ALL()
      INTEGER NI,NIMAP
      REAL :: TA1,TA2,TB1,TB2; INTEGER TIB1,TIB2,TIRATE
!
      !"UPDATE_R_ALL"
      CALL CPU_TIME(TA1); CALL SYSTEM_CLOCK(TIB1,TIRATE); TB1=REAL(TIB1,4)/REAL(TIRATE,4)
      WH%R=0; WH%R0=0
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.EQ.NI)THEN
          CALL CALC_R(GYZ(NI)%CO,GYZ(NI)%PJ)
        ELSE
          GYZ(NI)%CO%R =GYZ(NIMAP)%CO%R
          GYZ(NI)%CO%R0=GYZ(NIMAP)%CO%R0
        ENDIF
      ENDDO
      CALL ZSUM_ALL_MPI(WH%R ,WH%NA2MAX*WH%NA2MAX*WH%NIONS)
      CALL ZSUM_ALL_MPI(WH%R0,WH%NA2MAX*WH%NA2MAX*WH%NIONS)
      CALL ZOUT_MAT('R0-OUT',WH%R0,GL%IO)
      CALL ZOUT_MAT('R-OUT',WH%R,GL%IO)
      DO NI=1,WH%NIONS
        GYZ(NI)%CO%Z=MATMUL(GYZ(NI)%CO%R,CONJG(TRANSPOSE(GYZ(NI)%CO%R)))
      ENDDO
      CALL ZOUT_MAT('Z_NATURAL',WH%Z,GL%IO)
      DO NI=1,WH%NIONS
        GYZ(NI)%CO%Z=MATMUL(CONJG(TRANSPOSE(GYZ(NI)%CO%R)),GYZ(NI)%CO%R)
      ENDDO
      CALL ZOUT_MAT('Z_ORIG',WH%Z,GL%IO)
      IF(GL%LCHKLOC>0)THEN
        CALL CHK_EIGS_ANN(WH%Z,'Z_ORG','V')
      ENDIF
      CALL CPU_TIME(TA2); CALL SYSTEM_CLOCK(TIB2,TIRATE); TB2=REAL(TIB2,4)/REAL(TIRATE,4)
      CALL OUT_TIME_USE('UPDATE_R_ALL',TA2-TA1,TB2-TB1,GL%IO)
      IF(GL%LSCF.EQ.-11)STOP
      !"UPDATE_R_ALL"
      RETURN
!
      END SUBROUTINE UPDATE_R_ALL
!
!****************************************************************************
      SUBROUTINE CALC_LA(MODE)
      INTEGER NI,MODE
!
      DO NI=1,WH%NIONS
        CALL CALC_CO_LA(GYZ(NI)%CO,GL%LDC,MODE)
      ENDDO
      IF(MODE==2)THEN
        CALL SYM_AN(WH%LA2,1)
        CALL ZOUT_MAT('LA2-SYM',WH%LA2,GL%IO)
      ELSE
        CALL SYM_AN(WH%LA1,1)
        CALL ZOUT_MAT('LA1-SYM',WH%LA1,GL%IO)
      ENDIF
      RETURN
!
      END SUBROUTINE CALC_LA
!
!****************************************************************************
      SUBROUTINE CALC_ENG()
!
      CALL CALC_ONSITE()
      CALL CALC_EDBL()
      IF(GL%LMODEL==0)THEN
        CALL CALC_EBND()
      ELSE
        ENG%BAND=0
      ENDIF
      CALL CALC_E_HYBRD()
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%HYBRD:",F12.5)')ENG%HYBRD
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%BAND:",F12.5)')ENG%BAND
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%TS2 :",F12.5)')ENG%TS2
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%DL(DEEP LEVEL):",F12.5)')ENG%DL
      ENG%TB=ENG%BAND+ENG%DBLC+ENG%GAMM
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%TB(ENG%BAND+ENG%DBLC+ENG%GAMM)=",F12.5)')ENG%TB
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%TB-DC2=",F12.5)')ENG%TB-ENG%DC2
      RETURN
!
      END SUBROUTINE CALC_ENG
!
!****************************************************************************
      SUBROUTINE WRT_PROJ_C()
      INTEGER NI,NIMAP
!
      ! "WRT_PROJ_C"
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP; IF(NIMAP.NE.NI)CYCLE
        CALL DUMP_PJ_C(GYZ(NI)%PJ,NI)
        IF(GL%LPJC==1)THEN
          CALL DUMP_PJ_C01(GYZ(NI)%PJ,NI)
        ENDIF
      ENDDO
      ! "WRT_PROJ_C"
!
      END SUBROUTINE WRT_PROJ_C
!
!****************************************************************************
! Calculate reduced local configuration density matrix 
!****************************************************************************
      SUBROUTINE CALC_PJRHO()
      INTEGER NI,NIMAP
!
      ! "CALC_PJRHO"
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP; IF(NIMAP.NE.NI)CYCLE
        GYZ(NI)%HL%NI=NI
        CALL CALC_PJ_RHO(GYZ(NI)%HL,GYZ(NI)%PJ)
        IF(GL%LSCF==106)THEN
          CALL ZBM_DUMP(GYZ(NI)%PJ%RHO,NI,0,9)
        ENDIF
        CALL OUT_LOCAL_DBOCC1(GYZ(NI)%PJ,GYZ(NI)%FS,GYZ(NI)%HL,GL%IO,GL%LIADBOCC)
        CALL OUT_LOCAL_CURRENT1(GYZ(NI)%PJ,GYZ(NI)%FS,GYZ(NI)%HL,GL%IO,GL%LCURRENT)
      ENDDO
      ! "CALC_PJRHO"
      RETURN
!
      END SUBROUTINE CALC_PJRHO
!
!****************************************************************************
      SUBROUTINE OUT_SEC_ALL()
      INTEGER NI,NIMAP,NT
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      NT=0
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP; IF(NIMAP.NE.NI)CYCLE
        IF(NT==GYZ(NI)%NT)CYCLE
        NT=GYZ(NI)%NT
        CALL OUT_SEC_A(GL%IO,GYZ(NI)%HL%SEC_N,'SEC_N ',1)
        SELECT CASE (GL%LGPRJ)
        CASE(1,2,21)
          CALL OUT_SEC_A(GL%IO,GYZ(NI)%HL%SEC_J,'SEC_J ',2)
        CASE(3)
          CALL OUT_SEC_A(GL%IO,GYZ(NI)%HL%SEC_J,'SEC_S ',2)
        CASE(4)
          CALL OUT_SEC_A(GL%IO,GYZ(NI)%HL%SEC_J,'SEC_J ',1)
        END SELECT
      ENDDO
      RETURN
!
      END SUBROUTINE OUT_SEC_ALL
!
!******************************************************************************
      SUBROUTINE DUMP_V2AO_ALL()
      INTEGER NT
!
      OPEN(GL%IU,FILE='V2AO.OUT',STATUS='REPLACE')
      DO NT=1,WH%NTYP
        CALL WRITE_V2AO(GYS(NT)%CO,NT,GL%IU)
        DEALLOCATE(GYS(NT)%CO%V2AO)
      ENDDO
      CLOSE(GL%IU)
      RETURN
!
      END SUBROUTINE DUMP_V2AO_ALL
!
!****************************************************************************
      SUBROUTINE CALC_RED_VSP()
      INTEGER NI,NIMAP
!
      IF(GL%LSCF/=107)RETURN
      ! "CALC_RED_VSP"
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP; IF(NIMAP.NE.NI)CYCLE
        ALLOCATE(GYZ(NI)%PJ%RHO_EV(GYZ(NI)%PJ%RHO%DIM))
        CALL ZBM_EIGENSYS(GYZ(NI)%PJ%RHO,GYZ(NI)%PJ%RHO_EV,-1)
        CALL CALC_RED_EVECNJ(GYZ(NI)%HL,GYZ(NI)%PJ,GL%RHO_CUT)
        CALL SEC1_SIMPLIFY(GYZ(NI)%HL%SEC_J)
        IF(GL%LGPRJ==11.OR.GL%LGPRJ==14)THEN
          CALL ALLOC_ZBM(GYZ(NI)%HL%EVEC,GYZ(NI)%HL%SEC_N%DIM,GYZ(NI)%HL%SEC_N%ID,0,0)
          CALL SET_IDEN_ZBM(GYZ(NI)%HL%EVEC)
        ENDIF
        CALL ZBM_AXB(GYZ(NI)%HL%EVEC,GYZ(NI)%PJ%RHO)
        CALL ZBM_DUMP(GYZ(NI)%HL%EVEC,NI,0,16)
        CALL DUMP_SEC1(GYZ(NI)%HL%SEC_N,NI,17)
        CALL DUMP_SEC1(GYZ(NI)%HL%SEC_J,NI,18)
      ENDDO
      ! "CALC_RED_VSP"
      RETURN
!
      END SUBROUTINE CALC_RED_VSP
!
!****************************************************************************
      SUBROUTINE CALC_FSP0()
      INTEGER NI,NIMAP
!
      ! "CALC_FSP0"
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP; IF(NIMAP.NE.NI)CYCLE
        CALL CALC_P0_FS(GYZ(NI)%CO,GYZ(NI)%FS,GYZ(NI)%HL,GYZ(NI)%PJ)
      ENDDO
      ! "CALC_FSP0"
      RETURN
!
      END SUBROUTINE CALC_FSP0
!
!****************************************************************************
      SUBROUTINE CLR_PROJ()
      INTEGER NI,NIMAP
!
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP; IF(NIMAP.NE.NI)CYCLE
        CALL CLR_PJ(GYZ(NI)%CO,GYZ(NI)%PJ)
      ENDDO
      RETURN
!
      END SUBROUTINE CLR_PROJ
!
!****************************************************************************
! Calculate reduced configuration weight
!****************************************************************************
      SUBROUTINE CALC_RCW()
      INTEGER NI,NIMAP,I,NMAX
      INTEGER,POINTER::NDIM(:),NF(:,:)
      REAL(gq),POINTER::CWN(:,:)
!
      ! "CALC_RCW"
      NMAX=0
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP; IF(NIMAP.NE.NI)CYCLE
        NMAX=MAX(NMAX,GYZ(NI)%HL%SEC_N%DIM)
      ENDDO
      CALL IMAX1_ALL_MPI(NMAX)
      ALLOCATE(NDIM(WH%NIONS),NF(NMAX,WH%NIONS),CWN(NMAX,WH%NIONS))
      NDIM=0; NF=0; CWN=0
!
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP; IF(NIMAP.NE.NI)CYCLE
        NDIM(NI)=GYZ(NI)%HL%SEC_N%DIM
        CALL CALC_OCC_CONFIG_WT(GYZ(NI)%HL,GYZ(NI)%PJ)
        NF(1:NDIM(NI),NI)=NINT(GYZ(NI)%HL%SEC_N%VAL)
        CWN(1:NDIM(NI),NI)=GYZ(NI)%PJ%CW_N
      ENDDO
!
      CALL ISUM_MASTER_MPI(NF,WH%NIONS*NMAX)
      CALL ISUM_MASTER_MPI(NDIM,WH%NIONS)
      CALL DSUM_MASTER_MPI(CWN,NMAX*WH%NIONS)
      DO NI=1,WH%NIONS
        NIMAP=GYZ(NI)%NIMAP; IF(NIMAP.NE.NI)CYCLE
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3," OCC_CONFIG_WEIGHT:")')NI
        DO I=1,NDIM(NI)
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(8X,"OCC=",I2," WEIGHT=",F10.5)')NF(I,NI),CWN(I,NI)
        ENDDO
      ENDDO
      DEALLOCATE(NDIM,NF,CWN)
      ! "CALC_RCW"
      RETURN
!
      END SUBROUTINE CALC_RCW
!
!****************************************************************************
! Calculate entanglement entropy
!****************************************************************************
      SUBROUTINE CALC_ENTANGLES()
      INTEGER NI,NIMAP
      REAL(gq),ALLOCATABLE::ENS(:)
!
      ! "CALC_ENTANGLES"
      ALLOCATE(ENS(WH%NIONS)); ENS=0
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.NE.NI)CYCLE
        CALL CALC_ENS(GYZ(NI)%HL,GYZ(NI)%PJ,ENS(NI))
      ENDDO
      CALL DSUM_MASTER_MPI(ENS,WH%NIONS)
      DO NI=1,WH%NIONS
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.NE.NI)CYCLE
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3," ENTANGLEMENT ENTROPY=",F12.4)')NI,ENS(NI)
      ENDDO
      DEALLOCATE(ENS)
      ! "CALC_ENTANGLES"
      RETURN
!
      END SUBROUTINE CALC_ENTANGLES
!
!****************************************************************************
      SUBROUTINE CALC_PROJENS()
      INTEGER NI,NIMAP
      REAL(gq) S_MIX
      REAL(gq),ALLOCATABLE::ENS(:)
!
      ! "CALC_PROJENS"
      ALLOCATE(ENS(WH%NIONS)); ENS=0
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.NE.NI)CYCLE
        CALL CALC_PROJ_ENS(GYZ(NI)%HL,GYZ(NI)%PJ)
        ENS(NI)=GYZ(NI)%PJ%PJ_ENS
      ENDDO
      CALL DSUM_MASTER_MPI(ENS,WH%NIONS)
      DO NI=1,WH%NIONS
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP==NI)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3," PROJECTION ENTROPY=",F12.4)')NI,ENS(NI)
        ELSE
          ENS(NI)=ENS(NIMAP)
        ENDIF
      ENDDO
      IF(KPT%ISMEAR/=-1)THEN
        S_MIX=0
      ELSE
        S_MIX=-ENG%TS2/KPT%DELTA
      ENDIF
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" TOTAL PROJECTION ENTROPY=",F12.4," MERMIN ENTROPY=",F12.4)')SUM(ENS),S_MIX
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" TOTAL ELECTRONIC ENTROPY=",F12.4)')SUM(ENS)+S_MIX
      DEALLOCATE(ENS)
      ! "CALC_PROJENS"
      RETURN
!
      END SUBROUTINE CALC_PROJENS
!
!****************************************************************************
! Calculate <S^2>, <L^2> or <J^2>
!****************************************************************************
      SUBROUTINE CALC_SLJ()
      INTEGER NI
      REAL :: TA1,TA2,TB1,TB2; INTEGER TIB1,TIB2,TIRATE
!
      !"CALC_SLJ"
      CALL CPU_TIME(TA1); CALL SYSTEM_CLOCK(TIB1,TIRATE); TB1=REAL(TIB1,4)/REAL(TIRATE,4)
      GYZ(:)%CO%S2=0; GYZ(:)%CO%L2=0; GYZ(:)%CO%J2=0; GYZ(:)%CO%SZ=0
      IF(GL%LGPRJ==11.OR.GL%LGPRJ==14)RETURN
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        IF(GYZ(NI)%NIMAP.NE.NI)CYCLE
        CALL CALC_SLJ_1(GYZ(NI)%CO,GYZ(NI)%HL,GYZ(NI)%PJ)
      ENDDO
      IF(GL%LGPRJ==1.OR.GL%LGPRJ==2.OR.GL%LGPRJ==21)THEN
        CALL DSUM_MASTER_MPI(GYZ(:)%CO%J2,WH%NIONS)
      ENDIF
      IF(GL%LGPRJ==1.OR.GL%LGPRJ==2.OR.GL%LGPRJ==3.OR.GL%LGPRJ==21)THEN
        CALL DSUM_MASTER_MPI(GYZ(:)%CO%S2,WH%NIONS)
      ENDIF
      IF(GL%LGPRJ==1.OR.GL%LGPRJ==2.OR.GL%LGPRJ==21)THEN
        CALL DSUM_MASTER_MPI(GYZ(:)%CO%L2,WH%NIONS)
      ENDIF
      IF(GL%LGPRJ==4)THEN
        CALL DSUM_MASTER_MPI(GYZ(:)%CO%SZ,WH%NIONS)
      ENDIF
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIMAP.NE.NI)CYCLE
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3,\)')NI
        IF(GL%LGPRJ==1.OR.GL%LGPRJ==2.OR.GL%LGPRJ==21)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" J2=",F10.4,\)')GYZ(NI)%CO%J2
        ENDIF
        IF(GL%LGPRJ==1.OR.GL%LGPRJ==2.OR.GL%LGPRJ==3.OR.GL%LGPRJ==21)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" S2=",F10.4,\)')GYZ(NI)%CO%S2
        ENDIF
        IF(GL%LGPRJ==1.OR.GL%LGPRJ==2.OR.GL%LGPRJ==21)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" L2=",F10.4,\)')GYZ(NI)%CO%L2
        ENDIF
        IF(GL%LGPRJ==4)THEN
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" SZ=",F10.4,\)')GYZ(NI)%CO%SZ
        ENDIF
      ENDDO
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ")')
      CALL CPU_TIME(TA2); CALL SYSTEM_CLOCK(TIB2,TIRATE); TB2=REAL(TIB2,4)/REAL(TIRATE,4)
      CALL OUT_TIME_USE('    CALC_SLJ',TA2-TA1,TB2-TB1,GL%IO)
      !"CALC_SLJ"
      RETURN
!
      END SUBROUTINE CALC_SLJ
!
!****************************************************************************
      SUBROUTINE CALC_EBAND_COMP()
      INTEGER ISYM,IVEC,IKS,IKP,IKPL,NKP,NBANDS,ISP
      INTEGER NI,NEMIN,NEMAX,NASO,NASP,IA,NIP,NSYM,NBASE,NBASEP
      REAL(gq)WTK,WTK0
      COMPLEX(gq),POINTER :: VK(:,:),HK0(:,:)
      REAL(gq),POINTER    :: FERWE(:)
!
      !"CALC_EBAND_COMP"
      ALLOCATE(WH%TFC(BND%NASOTOT,BND%NSPIN),WH%TFF(WH%NASOTOT,WH%NASOTOT,BND%NSPIN))
      WH%TFC=0; WH%TFF=0; WH%TCC=0
      NSYM=SYM%IDF-SYM%IDI+1
!
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      NEMIN=BND%NE(2,IKP); NEMAX=BND%NE(3,IKP)
      NBANDS=NEMAX-NEMIN+1
      DO ISP=1,BND%NSPIN
      FERWE=>BND%FERWE(NEMIN:NEMAX,IKP,ISP)
      WTK=1._gq/NSYM/BND%RSPO; WTK0=KPT%WT(IKP)
      DO ISYM=SYM%IDI,SYM%IDF
      HK0=>BND%HK0(1:NBANDS,1:NBANDS,ISYM,IKPL)
      VK =>BND%VK (1:NBANDS,1:NBANDS,ISYM,IKPL,ISP)
      CALL CALC_ECOMP_1K(HK0,VK,FERWE,NBANDS,WTK,WTK0,ISP)
      ENDDO; ENDDO ! ISYM,ISP
      ENDDO; ENDDO
      NULLIFY(HK0,VK,FERWE)
      CALL DSUM_ALL_MPI(WH%TFF,WH%NASOTOT*WH%NASOTOT*BND%NSPIN)
      CALL DSUM_ALL_MPI(WH%TFC,WH%NASOTOT*BND%NSPIN)
      CALL DSUM1_ALL_MPI(WH%TCC)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" BAND ENERGY COMPONENT ANALYSIS (ORBITAL INDEX FAST)")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%NONL_BAND=",F12.5)')ENG%BAND+ENG%DC1
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%FCTOT=",F12.5)')SUM(WH%TFC)
      NBASE=1
      DO NI=1,WH%NIONS
        NASO=GYZ(NI)%CO%DIMSO
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3," ENG%FC =",14F12.5)')NI,WH%TFC(NBASE:NBASE+NASO-1,:)
        NBASE=NBASE+NASO
      ENDDO
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%FFTOT=",F12.5)')SUM(WH%TFF)
      NBASE=1
      DO NI=1,WH%NIONS
        NASO=GYZ(NI)%CO%DIMSO
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3," ENG%FF1=",14F12.5)')NI,((SUM(WH%TFF(IA,:,ISP)),IA=NBASE,NBASE+NASO-1),ISP=1,BND%NSPIN)
        NBASE=NBASE+NASO
      ENDDO
      NBASE=1
      DO NI=1,WH%NIONS; NASO=GYZ(NI)%CO%DIMSO
      DO IA=1,NASO
      NBASEP=NBASE
      DO NIP=NI,WH%NIONS
        NASP=GYZ(NIP)%CO%DIMSO
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3," IA=",I3," NIP=",I3," ENG%FF2=",14F12.5)') &
              &NI,IA,NIP,WH%TFF(NBASE+IA-1,NBASEP:NBASEP+NASP-1,:)
      NBASEP=NBASEP+NASP
      ENDDO; ENDDO 
      NBASE=NBASE+NASO
      ENDDO ! NI
!
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%CCTOT=",F12.5)')WH%TCC
      !"CALC_EBAND_COMP"
      RETURN
!
      END SUBROUTINE CALC_EBAND_COMP
!
!****************************************************************************
      SUBROUTINE CALC_ECOMP_1K(HK0,VK,FERWE,NBANDS,WTK,WTK0,ISP)
      INTEGER NBANDS,ISP
      REAL(gq) WTK,WTK0,FERWE(NBANDS)
      COMPLEX(gq) HK0(NBANDS,NBANDS),VK(NBANDS,NBANDS)
! LOCAL
      INTEGER IA
      INTEGER NASOT
      COMPLEX(gq),ALLOCATABLE :: NABR(:,:)
!
      ALLOCATE(NABR(NBANDS,NBANDS)); NABR=0
      CALL CALC_NABR_1K(NABR,VK,FERWE,NBANDS,WTK,WTK0,ISP)
      NABR=NABR*HK0
      NASOT=BND%NASOTOT
      DO IA=1,NASOT; WH%TFC(IA,ISP)=WH%TFC(IA,ISP)+SUM(NABR(IA,NASOT+1:))*2; ENDDO
      WH%TFF(:,:,ISP)=WH%TFF(:,:,ISP)+NABR(1:NASOT,1:NASOT)
      WH%TCC=WH%TCC+SUM(NABR(1+NASOT:,1+NASOT:))
      DEALLOCATE(NABR)
      RETURN
!
      END SUBROUTINE CALC_ECOMP_1K
!
!****************************************************************************
! For energy component analysis
! f_n <Psi_n|a> R_{a,A} R^+_{B,b} <b|Psi_n> 
!=R^+_{B,b} <b|Psi_n> f_n <Psi_n|a> R_{a,A}
!****************************************************************************
      SUBROUTINE CALC_NABR_1K(NABR,VK,FERWE,NBANDS,WTK,WTK0,ISP)
      INTEGER NBANDS,ISP
      REAL(gq) WTK,WTK0,FERWE(NBANDS)
      COMPLEX(gq) VK(NBANDS,NBANDS),NABR(NBANDS,NBANDS)
! LOCAL
      INTEGER IA,NASOT
      COMPLEX(gq) R(BND%NASOTOT,BND%NASOTOT)
!
      CALL CALC_NABR1_1K(NABR,VK,FERWE,NBANDS,WTK,ISP) ! R^+ <b|psi>f<psi|a>
      NASOT=BND%NASOTOT
      R=BND%R(:,:,ISP)
      CALL ZANMXBMM('N',NABR(:,1:NASOT),R,NBANDS,NASOT) ! R^+ <b|psi>f<psi|a> R  (B,A)
      NABR=TRANSPOSE(NABR) ! (A,B)
      NABR(1:NASOT,1:NASOT)=NABR(1:NASOT,1:NASOT)+(-BND%NRL(:,:,ISP)+BND%NC_PHY(:,:,ISP))*WTK0
      RETURN
!
      END SUBROUTINE CALC_NABR_1K
!
!****************************************************************************
! R^+ applied to right side only, for D
! R^+ <b|psi>f<psi|a>
!****************************************************************************
      SUBROUTINE CALC_NABR1_1K(NABR,VK,FERWE,NBANDS,WTK,ISP)
      INTEGER NBANDS,ISP
      REAL(gq) WTK,FERWE(NBANDS)
      COMPLEX(gq) VK(NBANDS,NBANDS),NABR(NBANDS,NBANDS)
! LOCAL
      INTEGER NASOT
      COMPLEX(gq) ZR(BND%NASOTOT,BND%NASOTOT)
!
      CALL CALC_NAB_1K(NABR,VK,FERWE,NBANDS,WTK)
      NASOT=BND%NASOTOT
      ZR=BND%R(:,:,ISP)
      CALL ZANNXBNM('C',ZR,NABR(1:NASOT,:),NASOT,NBANDS) ! R^+ <b|psi>f<psi|a>
      RETURN
!
      END SUBROUTINE CALC_NABR1_1K
!
!****************************************************************************
      SUBROUTINE CALC_NAB_1K(NAB,VK,FERWE,NBANDS,WTK)
      INTEGER NBANDS
      REAL(gq) WTK,FERWE(NBANDS)
      COMPLEX(gq) VK(NBANDS,NBANDS),NAB(NBANDS,NBANDS)
! LOCAL
      INTEGER IB
      COMPLEX(gq),ALLOCATABLE::VF(:,:)
!
      NAB=0
      ALLOCATE(VF(NBANDS,NBANDS))
      DO IB=1,NBANDS; VF(:,IB)=VK(:,IB)*FERWE(IB)*WTK; ENDDO
      CALL ZGEMM('N','C',NBANDS,NBANDS,NBANDS,Z1,VF,NBANDS,VK,NBANDS,Z0,NAB,NBANDS) ! <b|psi>f<psi|a>
      DEALLOCATE(VF)
      RETURN
!
      END SUBROUTINE CALC_NAB_1K
!
!****************************************************************************
! Calc converged onsite terms
!****************************************************************************
      SUBROUTINE CALC_ONSITE()
      INTEGER NI,NIMAP
!
      GYZ(:)%PJ%EGAMM=0; GYZ(:)%PJ%EPOT2=0
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.EQ.NI)THEN
          CALL CALC_EGAMM(GYZ(NI)%PJ)
          CALL DEALLOC_ZCSR(GYZ(NI)%PJ%U_KKH)
          IF(GYZ(NI)%PJ%LMCFLY)THEN
            CALL ZBM_LOAD(GYZ(NI)%HL%U,NI,0,1)
          ENDIF
          CALL CALC_UPOT (GYZ(NI)%HL,GYZ(NI)%PJ)
          IF(GYZ(NI)%PJ%LMCFLY)THEN
            CALL DEALLOC_ZBM(GYZ(NI)%HL%U)
          ENDIF
        ELSE
          GYZ(NI)%PJ%EGAMM=GYZ(NIMAP)%PJ%EGAMM
          GYZ(NI)%PJ%EPOT2=GYZ(NIMAP)%PJ%EPOT2
        ENDIF
      ENDDO
      CALL DSUM_ALL_MPI(GYZ(:)%PJ%EGAMM,WH%NIONS)
      CALL DSUM_ALL_MPI(GYZ(:)%PJ%EPOT2,WH%NIONS)
      ENG%GAMM=SUM(GYZ(:)%PJ%EGAMM); ENG%POT2=SUM(GYZ(:)%PJ%EPOT2)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GYZ%PJ%EGAMM:")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10F12.5)')GYZ(:)%PJ%EGAMM
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%GAMM:",F12.5)')ENG%GAMM
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%POT2:",F12.5)')ENG%POT2
      RETURN
!
      END SUBROUTINE CALC_ONSITE
!
!****************************************************************************
! Calc physical occupations
!****************************************************************************
      SUBROUTINE CALC_NCPHY()
      INTEGER NI,NIMAP,I
      REAL(gq) SUM_KS,SUM_PHY,SUM_VAR
!
      ! "CALC_NCPHY"
      WH%NC_PHY=0
      DO NI=1,WH%NIONS
        IF(GYZ(NI)%NIL.EQ.0)CYCLE
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.EQ.NI)THEN
          IF((GYZ(NI)%HL%EVEC%NBK<=0).AND.(GL%LGPRJ/=11).AND.(GL%LGPRJ/=14))THEN
            CALL ZBM_LOAD(GYZ(NI)%HL%EVEC,NI,0,2)
          ENDIF
          CALL CALC_NCPHY_COSYM(GYZ(NI)%FS,GYZ(NI)%HL,GYZ(NI)%CO,GYZ(NI)%PJ)
          SUM_KS=0; SUM_PHY=0; SUM_VAR=0
          DO I=1,GYZ(NI)%CO%DIM2
            SUM_KS =SUM_KS +GYZ(NI)%CO%NKS   (I,I)
            SUM_PHY=SUM_PHY+GYZ(NI)%CO%NC_PHY(I,I)
            SUM_VAR=SUM_VAR+GYZ(NI)%CO%NC_VAR(I,I)
          ENDDO
          WRITE(GL%IO,'(" NI=",I3," SUM_PHY-SUM_VAR=",F16.8                          )')NI,SUM_PHY-SUM_VAR
          WRITE(GL%IO,'(7X       ," SUM_PHY-SUM_KS =",F16.8," WOULD BE RENORMALIZED!")')   SUM_PHY-SUM_KS
          GYZ(NI)%CO%NC_PHY=GYZ(NI)%CO%NC_PHY/SUM_PHY*SUM_KS ! Renormalized
        ELSE
          GYZ(NI)%CO%NC_PHY=GYZ(NIMAP)%CO%NC_PHY
        ENDIF
      ENDDO
      CALL ZSUM_ALL_MPI(WH%NC_PHY,WH%NA2TOT)
      CALL MAP_WH_BND_R(WH%NC_PHY,BND%NC_PHY,.FALSE.)
      CALL ZOUT_MAT('NC_PHY_REN',WH%NC_PHY,GL%IO,1)
      ! "CALC_NCPHY"
      RETURN
!
      END SUBROUTINE CALC_NCPHY
!
!****************************************************************************
! Calc double counting energy
!****************************************************************************
      SUBROUTINE CALC_EDBL()
      INTEGER NI
      REAL(gq) UB,JB,NT,EDC2V
!
      DO NI=1,WH%NIONS
        CALL CALC_CO_EDCLA1(GYZ(NI)%CO)
      ENDDO
      IF(GL%LDC.NE.0)THEN
        DO NI=1,WH%NIONS
          CALL CALC_CO_EDCUJ(GYZ(NI)%CO,GL%LDC)
        ENDDO
      ENDIF
!
      GYZ(:)%CO%EDC=GYZ(:)%CO%EDCLA1+GYZ(:)%CO%EDCUJ
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" EDCLA1(NI) ARRAY:")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10F10.3)') GYZ(:)%CO%EDCLA1
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" EDCUJ(NI) ARRAY:")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10F10.3)') GYZ(:)%CO%EDCUJ
      IF(GL%LDC.EQ.2.OR.GL%LDC.EQ.12)THEN
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" EDCUJV(NI) ARRAY:")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10F10.3)') GYZ(:)%CO%EDCUJV
      ENDIF
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" EDC(NI) ARRAY:")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10F10.3)') GYZ(:)%CO%EDC
      ENG%DBLC=SUM(GYZ(:)%CO%EDC)
      ENG%DC1 =SUM(GYZ(:)%CO%EDCLA1)
      ENG%DC2 =SUM(GYZ(:)%CO%EDCUJ)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%DBLC:",F12.5)')ENG%DBLC
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%DC1 :",F12.5)')ENG%DC1
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%DC2 :",F12.5)')ENG%DC2
      IF(GL%LDC.EQ.2.OR.GL%LDC.EQ.12)THEN
      EDC2V=SUM(GYZ(:)%CO%EDCUJV)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ENG%DC2V:",F12.5)')EDC2V
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" DC_E_STD_TO_E_V:",F12.5)')EDC2V-ENG%DC2
      ENDIF
      RETURN
!
      END SUBROUTINE CALC_EDBL
!
!****************************************************************************
      SUBROUTINE CALC_DA0()
      INTEGER I
!
      BND%R=0._gq
      DO I=1,BND%NASOTOT; BND%R(I,I,:)=1._gq; ENDDO
      CALL CALC_DA_BND()
      CALL MAP_WH_BND_R(WH%D0,BND%D0,.TRUE.)
      CALL ZOUT_MAT('D0-UNSYM',WH%D0,GL%IO)
      IF(SYM%MODE==10)THEN
        CALL PAWSYM_ANN(WH%D0)
        CALL ZOUT_MAT('D0-PAWSYM',WH%D0,GL%IO,1)
      ENDIF
      IF(GL%LCHKLOC>0)THEN
        CALL CHK_EIGS_ANN(WH%D0,'D','N')
      ENDIF
      CALL SYM_AN(WH%D0,0)
      CALL ZOUT_MAT('D0-SYM',WH%D0,GL%IO)
      RETURN
!
      END SUBROUTINE CALC_DA0
!
!****************************************************************************
      SUBROUTINE CALC_DA()
      INTEGER NI
      COMPLEX(gq),ALLOCATABLE::MU(:,:,:)
!
      IF(GL%LGREEN==0)THEN
        CALL CALC_DA_BND()
      ELSE
        CALL CALC_DA_GREEN(GL%LGREEN)
      ENDIF
      CALL MAP_WH_BND_R(WH%D0,BND%D0,.TRUE.)
      CALL ZOUT_MAT('D0-UNSYM',WH%D0,GL%IO)
      CALL SYM_AN(WH%D0,0)
      CALL ZOUT_MAT('D0-SYM',WH%D0,GL%IO)
      CALL D0_TO_D()
      CALL ZOUT_MAT('D-SYM',WH%D,GL%IO)
      DO NI=1,WH%NIONS
        CALL SET_SYM_LOC_ARRAY(GYZ(NI)%CO%D,GYZ(NI)%CO%DR,GYZ(NI)%CO,.FALSE.,0)
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_DA
!
!****************************************************************************
      SUBROUTINE CALC_DA_GREEN(MODE)
      INTEGER MODE
! LOCAL
      INTEGER I,NI,NSPIN,NASO
      COMPLEX(gq),ALLOCATABLE::MU(:,:,:)
!
      NSPIN=BND%NSPIN
      DO NI=1,WH%NIONS
        NASO=GYZ(NI)%CO%DIMSO
        IF(GF%WTYP==0)THEN
          IF(MODE==1)THEN
            ALLOCATE(MU(NASO,NASO,NSPIN)); MU=GYZ(NI)%CO%LB1+GYZ(NI)%CO%ETB
            DO I=1,NASO; MU(I,I,:)=MU(I,I,:)-BND%EF; ENDDO
            CALL GREEN_DA(GYZ(NI)%CO%GF,MU,GYZ(NI)%CO%RB,GYZ(NI)%CO%DB,BND%EF)
            DEALLOCATE(MU)
          ELSE
            CALL GREEN_DA2(GYZ(NI)%CO%GF,GYZ(NI)%CO%RB,GYZ(NI)%CO%DB,BND%EF)
          ENDIF
        ELSEIF(GF%WTYP==1)THEN
          CALL GREEN_DA_R1(GYZ(NI)%CO%GF,GYZ(NI)%CO%RB,GYZ(NI)%CO%DB,1)
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_DA_GREEN
!
!****************************************************************************
      SUBROUTINE CALC_SELF_ENERGY_FULL()
      INTEGER NI,NSPIN,NASO
      COMPLEX(gq),ALLOCATABLE::MU(:,:,:)
!
      NSPIN=BND%NSPIN
      DO NI=1,WH%NIONS
        NASO=GYZ(NI)%CO%DIMSO
        ALLOCATE(MU(NASO,NASO,NSPIN)); MU=GYZ(NI)%CO%LB1+GYZ(NI)%CO%ETB
        CALL CALC_SELF_ENERGY(GYZ(NI)%CO%GF,GYZ(NI)%CO%RB,MU,BND%EF)
        DEALLOCATE(MU)
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_SELF_ENERGY_FULL
!
!****************************************************************************
      SUBROUTINE CALC_DA_BND()
      INTEGER ISYM,IVEC,IKS,IKP,IKPL,ISP,NKP,NBANDS,IA
      INTEGER NEMIN,NEMAX,NSYM
      REAL(gq)WTK
      COMPLEX(gq),POINTER :: VK(:,:),HK0(:,:)
      REAL(gq),POINTER    :: FERWE(:)
!
      NSYM=SYM%IDF-SYM%IDI+1
      BND%D0=0
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      NEMIN=BND%NE(2,IKP); NEMAX=BND%NE(3,IKP)
      NBANDS=NEMAX-NEMIN+1
      WTK=1._gq/NSYM/BND%RSPO
      DO ISP=1,BND%NSPIN
      FERWE=>BND%FERWE(NEMIN:NEMAX,IKP,ISP)
      DO ISYM=SYM%IDI,SYM%IDF
      HK0=>BND%HK0(1:NBANDS,1:NBANDS,ISYM,IKPL)
      VK =>BND%VK (1:NBANDS,1:NBANDS,ISYM,IKPL,ISP)
      CALL CALC_DA_1K(HK0,VK,FERWE,NBANDS,WTK,ISP)
      ENDDO; ENDDO ! ISYM,ISP
      ENDDO; ENDDO
      CALL ZSUM_ALL_MPI(BND%D0,WH%NASOTOT*WH%NASOTOT*BND%NSPIN)
      DO ISP=1,BND%NSPIN; BND%D0(:,:,ISP)=TRANSPOSE(BND%D0(:,:,ISP)); ENDDO
      RETURN
!
      END SUBROUTINE CALC_DA_BND
!
!****************************************************************************
! f_n <Psi_n|a> H_{A,B} R^+_{B,b} <b|Psi_n>
!=H_{A,B} R^+_{B,b} <b|Psi_n> f_n <Psi_n|a>
!****************************************************************************
      SUBROUTINE CALC_DA_1K(HK0,VK,FERWE,NBANDS,WTK,ISP)
      INTEGER NBANDS,ISP
      REAL(gq) WTK,FERWE(NBANDS)
      COMPLEX(gq) HK0(NBANDS,NBANDS),VK(NBANDS,NBANDS)
! LOCAL
      INTEGER IA
      INTEGER NASOT
      COMPLEX(gq),ALLOCATABLE :: NABR(:,:),ZD(:,:)
!
      NASOT=BND%NASOTOT
      ALLOCATE(NABR(NBANDS,NBANDS),ZD(NASOT,NASOT))
      CALL CALC_NABR1_1K(NABR,VK,FERWE,NBANDS,WTK,ISP) ! R^+ <b|psi>f<psi|a>
      ZD=0
      CALL ZGEMM('N','N',NASOT,NASOT,NBANDS,Z1,HK0(1:NASOT,:),NASOT,NABR(:,1:NASOT),NBANDS,Z0,ZD,NASOT)
      BND%D0(:,:,ISP)=BND%D0(:,:,ISP)+ZD
      DEALLOCATE(ZD)
      RETURN
!
      END SUBROUTINE CALC_DA_1K
!
!****************************************************************************
! LGREEN=1: VK --H_eff^G, calc GF%G--quasiparticle G
! LGREEN=2: VK0--H_bare,  calc GF%G--Coherent part G
!****************************************************************************
      SUBROUTINE SUM_GFK_FULL(MU)
      REAL(gq) MU
! LOCAL
      INTEGER ISYM,IVEC,IKS,IKP,IKPL,ISP,NKP,NBANDS,IA
      INTEGER NEMIN,NEMAX,NSYM,NBMAX,IW
      REAL(gq)WTK
      COMPLEX(gq) OMG
      COMPLEX(gq),POINTER :: HK(:,:)
      COMPLEX(gq),ALLOCATABLE :: GK(:,:)
!
      ! 'SUM_GFK_FULL'
      NSYM=SYM%IDF-SYM%IDI+1
      NBMAX=BND%NMAXIN
      ALLOCATE(GK(NBMAX,NBMAX)); GK=0
      GF%G=0
!
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      NEMIN=BND%NE(2,IKP); NEMAX=BND%NE(3,IKP)
      NBANDS=NEMAX-NEMIN+1
      WTK=KPT%WT(IKP)/NSYM
      DO ISP=1,BND%NSPIN
      DO ISYM=SYM%IDI,SYM%IDF
      IF(GL%LGREEN==1)THEN
        HK=>BND%VK(1:NBANDS,1:NBANDS,ISYM,IKPL,ISP)
      ELSEIF(GL%LGREEN==2)THEN
        HK=>BND%HK0(1:NBANDS,1:NBANDS,ISYM,IKPL)
      ENDIF
      DO IW=1,GF%NW
        OMG=GF%W(IW)
        IF(GL%LGREEN==1)THEN
          CALL CALC_GF_1K(HK,GK,NBANDS,NBMAX,MU,OMG,BND%EBMAX,BND%NASOTOT)
        ELSEIF(GL%LGREEN==2)THEN
          CALL CALC_GF_1K(HK,GK,NBANDS,NBMAX,MU,OMG,BND%EBMAX,BND%NASOTOT,GF%M(:,:,IW,ISP))
        ENDIF
        GF%G(:,:,IW,ISP)=GF%G(:,:,IW,ISP)+GK*WTK
      ENDDO
      ENDDO; ENDDO ! ISYM,ISP
      ENDDO; ENDDO
      CALL ZSUM_ALL_MPI(GF%G,NBMAX*NBMAX*GF%NW*BND%NSPIN)
      NULLIFY(HK); DEALLOCATE(GK)
      ! 'SUM_GFK_FULL'
      RETURN
!
      END SUBROUTINE SUM_GFK_FULL
!
!******************************************************************************
! Quasi-particle impurity Green function
!******************************************************************************
      SUBROUTINE CALC_GF_IMP_ALL()
      INTEGER I,IB
!
      IF(GF%WTYP==1)THEN
        DO I=1,GF%NLR
          IF(GF%NLR==1.OR.I==2)THEN
            IB=GF%NBMAX
          ELSE
            IB=1
          ENDIF
          CALL GET_RENORM_GAMMA1(GF%GM(:,:,:,:,I),BND%R,GF%NBMAX,IB)
        ENDDO
      ENDIF
      DO I=1,GF%NLR
        IF(GF%NLR==1.OR.I==2)THEN
          IB=GF%NBMAX
        ELSE
          IB=1
        ENDIF
        CALL GET_RENORM_HYBRIZ1(GF%D1,BND%R,GF%D(:,:,:,:,I),GF%NBMAX,IB)
      ENDDO
      IF(GL%LENSEMBLE==1)THEN ! Grand canonical ensemble
        CALL CALC_GF_IMP(0._gq)
      ELSE
        CALL GREEN_FERMI_NV()
      ENDIF
      RETURN
!
      END SUBROUTINE CALC_GF_IMP_ALL
!
!****************************************************************************
      SUBROUTINE CALC_GF_IMP(MU)
      REAL(gq) MU
! LOCAL
      INTEGER ISP,IW,ILR
      COMPLEX(gq) H(GF%NBMAX,GF%NBMAX)
!
      DO ISP=1,GF%NSPIN; DO IW=1,GF%NW
        H=BND%LA1(:,:,ISP)+BND%ETA(:,:,ISP)
        DO ILR=1,GF%NLR
          H=H+GF%D(:,:,IW,ISP,ILR)
        ENDDO
        CALL CALC_GF_1K(H,GF%G(:,:,IW,ISP),GF%NBMAX,GF%NBMAX,MU,GF%W(IW),100._gq,1) ! No \mu in G_imp
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE CALC_GF_IMP
!
!****************************************************************************
      SUBROUTINE CALC_NKS(MODE)
      INTEGER MODE
! LOCAL
      INTEGER IA,NI,NA2
!
      CALL CALC_NKS_ORTH(MODE)
      CALL ZSUM_ALL_MPI(WH%NKS,WH%NIONS*WH%NA2MAX*WH%NA2MAX)
      CALL ZOUT_MAT('NKS-UNSYM',WH%NKS,GL%IO,1)
      IF(SYM%MODE==10)THEN
        CALL PAWSYM_ANN(WH%NKS)
        CALL ZOUT_MAT('NKS-PAWSYM',WH%NKS,GL%IO,1)
      ENDIF
      IF(GL%LCHKLOC>0)THEN
        CALL CHK_EIGS_ANN(WH%NKS,'NKS','V')
      ENDIF
      CALL SYM_AN(WH%NKS,1)
      CALL ZOUT_MAT('NKS-SYM',WH%NKS,GL%IO,1)
      DO NI=1,WH%NIONS
        NA2=GYZ(NI)%CO%DIM2
        GYZ(NI)%CO%NET=0
        DO IA=1,NA2
          GYZ(NI)%CO%NET=GYZ(NI)%CO%NET+ WH%NKS(IA,IA,NI)
        ENDDO
        CALL SET_SYM_LOC_N0R(GYZ(NI)%CO)
      ENDDO
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NELE_LOC TOTAL:")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(4X,5F)')GYZ(:)%CO%NET
      RETURN
!
      END SUBROUTINE CALC_NKS
!
!****************************************************************************
      SUBROUTINE CALC_NKS_ORTH(MODE)
      INTEGER MODE
!LOCAL
      INTEGER NI
!
      WH%NKS=0
      DO NI=1,WH%NIONS
        CALL CALC_CO_NKS(GYZ(NI)%CO,MODE)
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_NKS_ORTH
!
!****************************************************************************
      SUBROUTINE D0_TO_D()
      INTEGER NI
!
      DO NI=1,WH%NIONS
        CALL CO_D0_TO_D(GYZ(NI)%CO)
      ENDDO
      RETURN
!
      END SUBROUTINE D0_TO_D
!
!****************************************************************************
      SUBROUTINE CHK_W_ETA_ALL()
      INTEGER NI
      REAL(gq) NTOT(WH%NASOMAX,BND%NSPIN)
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      IF(GF%WTYP/=1)THEN
        RETURN
      ENDIF
      WRITE(GL%IO,'(" GF%W NORMALIZATION CHECK:")')
      DO NI=1,WH%NIONS
        CALL GREEN_OCC_R1(0,1,GYZ(NI)%CO%GF,GYZ(NI)%CO%RB,NELE=NTOT(1:GYZ(NI)%CO%DIMSO,:))
        WRITE(GL%IO,'(" NI=",I3," NORM=",14F10.6)')NI,NTOT(1:GYZ(NI)%CO%DIMSO,:)
        IF(MAXVAL(ABS(NTOT(1:GYZ(NI)%CO%DIMSO,:)-1))>.001_gq)THEN
          WRITE(GL%IO,'(" WARNING! PLEASE INCREASE GF%NW!")')
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE CHK_W_ETA_ALL
!
!****************************************************************************
      SUBROUTINE CALC_PSI_COMP()
      INTEGER IVEC,IKS,IKP,IKPL,ISP,NKP,IB,IB_,NI1,NI2,I,J,ISYM,NATM,NATMP
      INTEGER NBL,NBU,NBMAX,NASO,IU,NBANDS,NSYM,NBASE,NASOT,IR,IRL,ISUM
      REAL(gq) RES,RSUM,RESG(14)
      REAL(gq),ALLOCATABLE::PSIA(:,:,:,:)
      COMPLEX(gq),POINTER::VK(:)
!
      IU=GL%IU
      IF(GP%MYRANK.EQ.GP%MASTER)OPEN(IU,FILE='GLQTL.OUT',STATUS='REPLACE')
      NBL=MAXVAL(BND%NE(2,:)); NBU=MINVAL(BND%NE(3,:)); NBMAX=MAXVAL(BND%NE(1,:))
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" CALC_PSI_COMP:")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("   NBL,NBU,NBMAX=",3I4)')NBL,NBU,NBMAX
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("   MIN/MAX(E)=",2F8.2)')MINVAL(BND%EK(NBL,:,:)),MAXVAL(BND%EK(NBU,:,:))
      ALLOCATE(PSIA(GL%NASOMAX,WH%NIONS,KPT%DIM,BND%NSPIN))
      NASOT=BND%NASOTOT
      NSYM=SYM%IDF-SYM%IDI+1
!
      DO IB=NBL,NBU
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(IU,'(" BAND:",I4)')IB-NBL+1
      PSIA=0
!
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      IB_=IB-BND%NE(2,IKP)+1
      DO ISP=1,BND%NSPIN; DO ISYM=SYM%IDI,SYM%IDF
      VK=>BND%VK(1:NASOT,IB_,ISYM,IKPL,ISP)
      NBASE=1
      DO NI1=1,WH%NIONS
        NASO=GYZ(NI1)%CO%DIMSO
        IF(NI1==GYZ(NI1)%NIMAP)THEN
          PSIA(1:NASO,NI1,IKP,ISP)=PSIA(1:NASO,NI1,IKP,ISP)+VK(NBASE:NBASE+NASO-1)*CONJG(VK(NBASE:NBASE+NASO-1))/NSYM
        ENDIF
        NBASE=NBASE+NASO
      ENDDO; ENDDO; ENDDO
      ENDDO; ENDDO
!
      CALL DSUM_MASTER_MPI(PSIA,WH%NASOMAX*WH%NIONS*KPT%DIM*BND%NSPIN)
!
      IF(GP%MYRANK.EQ.GP%MASTER)THEN
      DO ISP=1,BND%NSPIN
      DO IKP=1,KPT%DIM
      RSUM=SUM(PSIA(:,:,IKP,ISP))
      DO NI1=1,GL%NAT_TOT
      DO NI2=1,WH%NIONS
        IF(NI2/=GYZ(NI2)%NIMAP)THEN
          CYCLE
        ENDIF
        IF(GYZ(NI2)%NI0==NI1)THEN
          GOTO 100
        ENDIF
      ENDDO ! NI2
      WRITE(IU,1050)BND%EK(IB,IKP,ISP),NI1,0._gq
      CYCLE
100   CONTINUE
      RES=SUM(PSIA(:,NI2,IKP,ISP))
      RESG=0; ISUM=0
      DO IR=1,GYZ(NI2)%CO%SYMH%NR; DO IRL=1,GYZ(NI2)%CO%SYMH%DIMP(IR)
        I=GYZ(NI2)%CO%SYMH%IDP(1,IRL,IR); J=GYZ(NI2)%CO%SYMH%IDP(2,IRL,IR)
        IF(I/=J)CYCLE
        RESG(IR)=RESG(IR)+PSIA(I,NI2,IKP,ISP)
      ENDDO; ENDDO
!
      WRITE(IU,1050)BND%EK(IB,IKP,ISP),NI1,RES,RES,RESG(1:GYZ(NI2)%CO%SYMH%NR)
      ENDDO ! NI1
      WRITE(IU,1050)BND%EK(IB,IKP,ISP),NI1,1-RSUM
      ENDDO; ENDDO ! IKP,ISP
      ENDIF
!
      ENDDO ! IB
      CLOSE(IU)
!
      DEALLOCATE(PSIA); NULLIFY(VK)
      RETURN
1050  FORMAT(F10.5,I3,F8.5,3X,39F8.5)
!
      END SUBROUTINE CALC_PSI_COMP
!
!****************************************************************************
      SUBROUTINE OUT_BANDS()
      INTEGER IU,I,NBL,NBU,NBMAX,ISP
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      NBL=MINVAL(BND%NE(2,:)); NBU=MAXVAL(BND%NE(3,:)); NBMAX=MAXVAL(BND%NE(1,:))
      NBU=MIN(NBU,NBMAX)
      IU=GL%IU
      OPEN(IU,FILE='GLBANDS.OUT',STATUS='REPLACE')
      DO ISP=1,BND%NSPIN
      DO I=1,KPT%DIM
        WRITE(IU,6000)KPT%X(I),KPT%Y(I),KPT%Z(I),KPT%NAME(I)
        WRITE(IU,*)
        WRITE(IU,'("     EIGENVALUES ARE:")')
        WRITE(IU,6001)BND%EK(NBL:NBU,I,ISP)
        WRITE(IU,6002)
      ENDDO; ENDDO
      CLOSE(IU)
      RETURN
6000  FORMAT(5X,'K=',3F10.5,3X,A10)
6001  FORMAT(8(/2X,5F13.7))
6002  FORMAT(13X,' EIGENVALUES BELOW THE ENERGY ')
!
      END SUBROUTINE OUT_BANDS
!
!****************************************************************************
      SUBROUTINE INI_WH_X(IU)
      INTEGER IU
! LOCAL
      INTEGER I
      LOGICAL LEXIST
      COMPLEX(gq) LA1(WH%DIMX1)
!
      IF(GL%ITER_RHO.GT.0)RETURN
      WH%X=0._gq; LA1=0; BND%EF=0
      INQUIRE(FILE='GLX.INP',EXIST=LEXIST)
      IF(LEXIST)THEN
        OPEN(IU,FILE='GLX.INP',STATUS='OLD')
        READ(IU,*)WH%X(1:WH%DIMX2)
        READ(IU,*,ERR=100,END=100)WH%X(1+WH%DIMX2:)
        READ(IU,*,ERR=100,END=100)BND%EF
100     CLOSE(IU)
      ENDIF
      IF(BND%ISPIN==2)THEN
        IF(GL%LMUSHIFT==1)THEN
          LA1=WH%X(WH%DIMX1+1:WH%DIMX1*2); WH%X=0._gq
        ELSEIF(.NOT.LEXIST)THEN
          STOP "ERROR: GLX.INP MUST BE PROVIDED FOR MAGNEIC CALCULATIONS."
        ENDIF
      ENDIF
      IF(.NOT.LEXIST.OR.GL%LMUSHIFT==1)THEN
        CALL SET_WH_X_DEFAULT()
      ENDIF
      IF(GL%LMUSHIFT==1)THEN
        WH%X(WH%DIMX1+1:WH%DIMX2)=WH%X(WH%DIMX1+1:WH%DIMX2)+LA1
      ENDIF
      RETURN
!
      END SUBROUTINE INI_WH_X
!
!****************************************************************************
      SUBROUTINE SET_WH_X_DEFAULT()
      INTEGER I,NI
      LOGICAL LEXIST
!
      WH%LA1=WH%EL0; WH%ETA=0
      CALL SYM_AN(WH%LA1,1)
      WH%LA1=WH%LA1+GL%LA1SHIFT  ! For sake of hybrd1, somehow it helps a bit.
      WH%R=0
      DO I=1,WH%NA2MAX
        IF(GL%LXGUESS==0)THEN
          WH%R(I,I,:)=1._gq
        ELSE
          WH%R(I,I,:)=0.999_gq
        ENDIF
      ENDDO
      INQUIRE(FILE='WH_RLNEF.INP',EXIST=LEXIST)
      IF(LEXIST)THEN
        CALL READ_WH_RLNEF(GL%IU); BND%EF=WH%EF
      ENDIF
      CALL SYM_AN(WH%R,0); CALL SYM_AN(WH%LA1,1); CALL SYM_AN(WH%NKS,1); CALL SYM_AN(WH%ETA,2)
      CALL SET_WH_X123(.TRUE.)
      RETURN
!
      END SUBROUTINE SET_WH_X_DEFAULT
!
!****************************************************************************
      SUBROUTINE SET_WH_X123(LRED)
      LOGICAL LRED
!
      CALL SET_WH_X1(WH%X(1         :WH%DIMX1),WH%R,GL%NA2MAX,LRED,0)
      CALL SET_WH_X1(WH%X(1+WH%DIMX1:WH%DIMX2),WH%LA1,GL%NA2MAX,LRED,1)
      CALL SET_WH_X1(WH%X(1+WH%DIMX2:WH%DIMX3),WH%NKS,GL%NA2MAX,LRED,1)
      IF(WH%DIMXHO>0)THEN
        CALL SET_WH_X1(WH%X(1+WH%DIMX3:WH%DIMXT),WH%ETA,GL%NA2MAX,LRED,2)
      ELSE
        IF(.NOT.LRED)THEN
          WH%ETA=0._gq
        ENDIF
      ENDIF
      RETURN
!
      END SUBROUTINE SET_WH_X123
!
!****************************************************************************
! MODE = 0: -->SYMG
!        1: -->SYMH
!        2: -->SYMHO
!****************************************************************************
      SUBROUTINE SET_WH_X1(X,Y,NY,LRED,MODE)
      INTEGER NY,MODE
      COMPLEX(gq) X(*),Y(NY,NY,WH%NIONS)
      LOGICAL LRED
! LOCAL
      INTEGER NI,I1,I2,DIM2,ID,NX,ISUM(WH%DIMXG)
      INTEGER,POINTER::ID_(:,:)
!
      SELECT CASE(MODE)
      CASE(0); NX=WH%DIMXG
      CASE(1); NX=WH%DIMXH
      CASE(2); NX=WH%DIMXHO
      END SELECT
      IF(LRED)THEN
        X(1:NX)=0; ISUM=0
      ELSE
        Y=0
      ENDIF
!
      DO NI=1,WH%NIONS
      SELECT CASE(MODE)
      CASE(0); ID_=>GYZ(NI)%CO%SYMG%ID
      CASE(1); ID_=>GYZ(NI)%CO%SYMH%ID
      CASE(2); ID_=>GYZ(NI)%CO%SYMHO%ID
      END SELECT
      DIM2=GYZ(NI)%CO%DIM2
      DO I1=1,DIM2; DO I2=1,DIM2
        ID=ID_(I1,I2)
        IF(ID.EQ.0)CYCLE
        IF(LRED)THEN
          X(ID)=X(ID)+Y(I1,I2,NI)
          ISUM(ID)=ISUM(ID)+1
        ELSE
          Y(I1,I2,NI)=X(ID)
          IF(MODE==0.OR.I1==I2)CYCLE
          Y(I2,I1,NI)=CONJG(X(ID))
        ENDIF
      ENDDO; ENDDO; ENDDO
      IF(LRED)THEN
        DO I1=1,NX
          IF(ISUM(I1).EQ.0)CYCLE
          X(I1)=X(I1)/ISUM(I1)
        ENDDO
      ENDIF
      RETURN
!
      END SUBROUTINE SET_WH_X1
!
!****************************************************************************
      SUBROUTINE SET_XTOWHX(X,WHX,LXREAL,NX,N,LBACK)
      INTEGER NX,N
      LOGICAL LBACK
      REAL(gq) X(NX)
      INTEGER LXREAL(N)
      COMPLEX(gq) WHX(N)
! LOCAL
      INTEGER I,J,ISUM
!
      ISUM=0
      IF(LBACK)THEN
        DO I=1,N
          ISUM=ISUM+1
          X(ISUM)=REAL(WHX(I),gq)
          IF(LXREAL(I)==2)THEN
            ISUM=ISUM+1
            X(ISUM)=AIMAG(WHX(I))
          ENDIF
        ENDDO
      ELSE
        DO I=1,N
          ISUM=ISUM+1
          WHX(I)=X(ISUM)
          IF(LXREAL(I)==2)THEN
            ISUM=ISUM+1
            WHX(I)=WHX(I)+DCMPLX(0._gq,X(ISUM))
          ENDIF
        ENDDO
      ENDIF
      RETURN
!
      END SUBROUTINE SET_XTOWHX
!
!****************************************************************************
      SUBROUTINE SET_WH_LXREAL()
      INTEGER NI,I1,I2,DIM2,ID,J
!
      WH%LXREAL=0
      DO NI=1,WH%NIONS; DIM2=GYZ(NI)%CO%DIM2
      DO I1=1,DIM2; DO I2=1,DIM2
        ID=GYZ(NI)%CO%SYMG%ID(I1,I2)
        IF(ID>0)THEN
          WH%LXREAL(ID)=1 *2  ! R
        ENDIF
        ID=GYZ(NI)%CO%SYMH%ID(I1,I2)
        IF(ID>0)THEN
          IF(I1==I2)THEN
            WH%LXREAL(ID+WH%DIMX1)=1 ! LA1
            WH%LXREAL(ID+WH%DIMX2)=1 ! NKS
          ELSE
            WH%LXREAL(ID+WH%DIMX1)=1 *2
            WH%LXREAL(ID+WH%DIMX2)=1 *2
          ENDIF
        ENDIF
        ID=GYZ(NI)%CO%SYMHO%ID(I1,I2)
        IF(ID>0)THEN
          WH%LXREAL(ID+WH%DIMX3)=1 *2
        ENDIF
      ENDDO; ENDDO; ENDDO
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" WH%LXREAL: R_X")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10I2)')WH%LXREAL(1:WH%DIMX1)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" WH%LXREAL: LA1_X")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10I2)')WH%LXREAL(1+WH%DIMX1:WH%DIMX2)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" WH%LXREAL: NA0_X")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10I2)')WH%LXREAL(1+WH%DIMX2:WH%DIMX3)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" WH%LXREAL: ETA_X")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10I2)')(WH%LXREAL(J),J=1+WH%DIMX3,WH%DIMXT)
      RETURN
!
      END SUBROUTINE SET_WH_LXREAL
!
!****************************************************************************
      SUBROUTINE SYM_AN(A,MODE)
      COMPLEX(gq) A(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
      INTEGER MODE
! LOCAL
      COMPLEX(gq) X(WH%DIMXG)
!
      CALL SET_WH_X1(X,A,WH%NA2MAX,.TRUE.,MODE)
      CALL SET_WH_X1(X,A,WH%NA2MAX,.FALSE.,MODE)
      RETURN
!
      END SUBROUTINE SYM_AN
!
!****************************************************************************
      SUBROUTINE GUTZ_SET_PROJ()
      INTEGER NI
!
      !"GUTZ_SET_PROJ"
      DO NI=1,WH%NIONS
        GYZ(NI)%HL%LGPRJ=GL%LGPRJ; GYZ(NI)%HL%LDIAPJ=GL%LDIAPJ; GYZ(NI)%HL%NI=NI
        GYZ(NI)%PJ%LMCFLY=GL%LMCFLY
        IF(GYZ(NI)%NIMAP.NE.NI)CYCLE; IF(GYZ(NI)%NIL.EQ.0)CYCLE
        CALL SET_PROJ(GYZ(NI)%CO,GYZ(NI)%FS,GYZ(NI)%HL,GYZ(NI)%PJ,GL%ITER_RHO)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3," DIM_PHIK=",I10," HL%DIM=",I10)')NI,GYZ(NI)%PJ%N_PHIK,GYZ(NI)%HL%DIM
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(6X," ID_PHIK_N=")')
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10I10)')GYZ(NI)%PJ%ID_PHIK_N
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(6X," ID_PHIK_J=")')
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10I10)')GYZ(NI)%PJ%ID_PHIK_J
      ENDDO
      !"GUTZ_SET_PROJ"
      RETURN
!
      END SUBROUTINE GUTZ_SET_PROJ
!
!****************************************************************************
      SUBROUTINE GUTZ_READ_PROJ()
      INTEGER NI
      REAL :: TA1,TA2,TB1,TB2; INTEGER TIB1,TIB2,TIRATE
!
      !"GUTZ_READ_PROJ"
      CALL CPU_TIME(TA1); CALL SYSTEM_CLOCK(TIB1,TIRATE); TB1=REAL(TIB1,4)/REAL(TIRATE,4)
      DO NI=1,WH%NIONS
        GYZ(NI)%HL%LGPRJ=GL%LGPRJ; GYZ(NI)%HL%LDIAPJ=GL%LDIAPJ; GYZ(NI)%HL%NI=NI
        GYZ(NI)%PJ%LMCFLY=GL%LMCFLY
        IF(GYZ(NI)%NIMAP.NE.NI)CYCLE; IF(GYZ(NI)%NIL.EQ.0)CYCLE
        CALL READ_PROJ(GL%IO,GYZ(NI)%CO,GYZ(NI)%FS,GYZ(NI)%HL,GYZ(NI)%PJ,GL%ITER_RHO,GL%LSCF)
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3," DIM_PHIK=",I10," HL%DIM=",I10)')NI,GYZ(NI)%PJ%N_PHIK,GYZ(NI)%HL%DIM
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(6X," ID_PHIK_N=")')
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10I10)')GYZ(NI)%PJ%ID_PHIK_N
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(6X," ID_PHIK_J=")')
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(10I10)')GYZ(NI)%PJ%ID_PHIK_J
      ENDDO
      CALL CPU_TIME(TA2); CALL SYSTEM_CLOCK(TIB2,TIRATE); TB2=REAL(TIB2,4)/REAL(TIRATE,4)
      CALL OUT_TIME_USE('GUTZ_RED_PRJ',TA2-TA1,TB2-TB1,GL%IO)
      !"GUTZ_READ_PROJ"
      RETURN
!
      END SUBROUTINE GUTZ_READ_PROJ
!
!******************************************************************************
      SUBROUTINE OUT_PAR_GL(NIONS,NTYP,IO)
      INTEGER NIONS,NTYP,IO
! LOCAL
      INTEGER NI,NT
      CHARACTER*10 LCALC(6)
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      IF(GL%LMODEL==0)THEN
        WRITE(GL%IO,'(" Lattice model.")')
      ELSE
        WRITE(GL%IO,'(" Impurity model.")')
      ENDIF      
      LCALC=''
      SELECT CASE(GL%LSCF)
      CASE(  1); LCALC(1)="GA"
      CASE(  2); LCALC(1)="GA-BAND" ! Plot bands
      CASE(  3); LCALC(1)="1Iter-GA"
      CASE(  4); LCALC(1)="PRT_PROJ" ! Print variational setup stuff.
      CASE(  5); LCALC(1)="PRT_PRJ_VE" 
      CASE(  6); LCALC(1)="GAPJ-SAVE" ! Run with fixed variational setup (generate if necessary)
      CASE( -1); LCALC(1)="HF"
      CASE(-11); LCALC(1)="GA:B-IN" ! Check R with B read in
      CASE(105); LCALC(1)="GA FIX-C"
      CASE(106); LCALC(1)="MEE FIT"
      CASE(107); LCALC(1)="GEN R-VSP" ! Reduced variational setup
      END SELECT
      SELECT CASE(GL%LDIAPJ)
      CASE( 1); LCALC(2)="DIAGONAL "
      CASE DEFAULT; LCALC(2)="FULL "
      END SELECT
      SELECT CASE(GL%LGPRJ)
      CASE( 1); LCALC(3)=" -NJ-PROJ"
      CASE( 2); LCALC(3)=" -NJC-PRJ" ! Crystal field as perturbation
      CASE( 3); LCALC(3)=" -NS-PROJ"
      CASE( 4); LCALC(3)="  -NSZ-PJ"
      CASE(11); LCALC(3)="    -N-PJ"
      CASE(12); LCALC(3)=" -N-PJ(DG)"
      CASE(14); LCALC(3)="-NSZ-PJ-II"
      CASE(21); LCALC(3)="  -NJR-PJ"
      END SELECT
      SELECT CASE(GL%LHUB)
      CASE( 0); LCALC(4)="INP-HUB"
      CASE( 1); LCALC(4)="SLATER"
      CASE( 2); LCALC(4)="KANAMORI"
      END SELECT
      SELECT CASE(GL%LDC)
      CASE(-1); LCALC(5)="HF DC"
      CASE( 0); LCALC(5)="NO DC"
      CASE( 1); LCALC(5)="STD DC"
      CASE( 2); LCALC(5)="FIXDC1"
      CASE( 3); LCALC(5)="FIXDC2"
      CASE(12); LCALC(5)="VAR-FIXDC"
      END SELECT
      SELECT CASE(GL%LUNIT)
      CASE( 0); LCALC(6)="   eV"
      CASE DEFAULT; LCALC(6)="   Ryd."
      END SELECT
      WRITE(GL%IO,'(" CALC SET: ",6A10)')LCALC
      IF(GL%LGREEN>0)THEN
        WRITE(GL%IO,'(" WITH GREEN FUNCTION APPROACH.")')
      ENDIF
!
      DO NT=1,NTYP
        WRITE(GL%IO,'(" NT=",I2)')NT
        CALL OUT_PAR(GYS(NT)%CO,GYS(NT)%FS,GYS(NT)%HL,IO)
      ENDDO
      RETURN
!
      END SUBROUTINE OUT_PAR_GL
!
!******************************************************************************
      SUBROUTINE OUT_PAR(CO,FS,HL,IO)
      INTEGER IO
      TYPE (CORR_ORB)    CO
      TYPE (FOCK_STATE)  FS
      TYPE(LOCAL_HAMIL)  HL
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      WRITE(IO,'(" LOCAL INTERACTION PART   :")')
      IF(GL%LBSCODE==0)THEN
        WRITE(IO,'(" U              :",F9.4," eV.")')CO%U
        WRITE(IO,'(" J              :",F9.4," eV.")')CO%J
      ELSE
        WRITE(IO,'(" U              :",F9.4," Ryd.",F9.4," eV.")')CO%U,CO%U*RYTOEV
        WRITE(IO,'(" J              :",F9.4," Ryd.",F9.4," eV.")')CO%J,CO%J*RYTOEV
      ENDIF
      WRITE(IO,'(" FOCK SPACE PART:")')
      WRITE(IO,'(" NA             :",I3)')FS%NA
      WRITE(IO,'(" FS%OCC(1) L.B. :",I3)')FS%OCC(1)
      WRITE(IO,'(" FS%OCC(2) U.B. :",I3)')FS%OCC(2)
      WRITE(IO,'(" FS%DIM         :",I10)')FS%DIM
      WRITE(IO,'(" HL%OCC(1) L.B. :",I3)')HL%OCC(1)
      WRITE(IO,'(" HL%OCC(2) U.B. :",I3)')HL%OCC(2)
      RETURN
!
      END SUBROUTINE OUT_PAR
!
!*************************************************************************************
      SUBROUTINE GET_EL0()
      INTEGER NI,IA,NA2
      LOGICAL LREL0
      CHARACTER(20) FMT
!
      !"GET_EL0"
      IF(GL%LEL0==1)THEN
        CALL READ_WHEL0(GL%IU)
        WRITE(GL%IO,'(" WH%EL0 READ IN.")')
      ELSEIF(GL%LMODEL==0)THEN
        WH%EL0=0
        DO NI=1,WH%NIONS
!          CALL CALC_CO_EL0_UV(GYZ(NI)%CO)
          CALL CALC_CO_EL0(GYZ(NI)%CO)
        ENDDO
        CALL ZSUM_ALL_MPI(WH%EL0,WH%NIONS*WH%NA2MAX*WH%NA2MAX)
      ENDIF
      CALL ZOUT_MAT('EL-UNSYM',WH%EL0,GL%IO,-1)
      IF(SYM%MODE==10)THEN
        CALL PAWSYM_ANN(WH%EL0)
        CALL ZOUT_MAT('EL-PAWSYM',WH%EL0,GL%IO,-1)
      ENDIF
      IF(GL%LCHKLOC>0)THEN
        CALL CHK_EIGS_ANN(WH%EL0,'EL','N')
      ENDIF
      CALL SYM_AN(WH%EL0,1)
      DO NI=1,WH%NIONS
        NA2=GYZ(NI)%CO%DIM2
        WRITE(FMT,'("(4X",I0,"F8.4)")')NA2 *2
        IF(GL%LSCF==5)THEN
          WRITE(*,'(" NI=",I3," CURRENT EL0 =")')NI
          WRITE(*,FMT)GYZ(NI)%CO%EL0
          WRITE(*,'("PLEASE ENTER NEW EL0:")')
          READ(*,*)GYZ(NI)%CO%EL0
        ENDIF
        GYZ(NI)%CO%ELC=0
        DO IA=1,NA2; GYZ(NI)%CO%ELC=GYZ(NI)%CO%ELC+GYZ(NI)%CO%EL0(IA,IA); ENDDO
        GYZ(NI)%CO%ELC=GYZ(NI)%CO%ELC/NA2 ! Level center
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" NI=",I3," LOCAL ORBITAL LEVEL CENTER(ELC)=",F8.3)')NI,GYZ(NI)%CO%ELC
      ENDDO
      CALL ZOUT_MAT('EL-SYM',WH%EL0,GL%IO,-1)
      IF(GL%LEL0==0)THEN
        CALL WRT_WHEL0(GL%IU)
      ENDIF
      !"GET_EL0"
      RETURN
!
      END SUBROUTINE GET_EL0
!
!*************************************************************************************
      SUBROUTINE PAWSYM_ANN(ANN)
      COMPLEX(gq) ANN(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
! LOCAL
      INTEGER NI,NA2,NA1,L,ISP
      COMPLEX(gq) R2N(WH%NA2MAX,WH%NA2MAX)
!
      IF(BND%ISO==2)THEN
        STOP ' ERROR: PAWSYM_ANN NOT READY FOR ISO=2!'
      ENDIF
!
      DO NI=1,WH%NIONS
        NA2=GYZ(NI)%CO%DIM2
        R2N=0
        R2N(1:NA2,1:NA2)=GYZ(NI)%CO%R2N
        IF(MAXVAL(ABS(R2N(1:NA2,1:NA2)-GYZ(NI)%CO%R2N))>1.E-10_gq)THEN ! REAL/CMPLX
          STOP ' ERROR IN PAWSYM_ANN: GYZ(NI)%CO%R2N IS NOT COMPATIBLE!'
        ENDIF
        CALL ZA2_TRANS(ANN(1:NA2,1:NA2,NI),NA2,R2N,-1) ! Goto real spherical harmonics basis
      ENDDO
      NA1=GYZ(1)%CO%DIM; NA2=GYZ(1)%CO%DIM2  ! Roughly
      L=(NA1-1)/2
      DO ISP=1,BND%NSPIN
        CALL ANNSYM1(ANN((ISP-1)*NA1+1:ISP*NA1,(ISP-1)*NA1+1:ISP*NA1,:),NA1,WH%NIONS,L)
      ENDDO
      IF(BND%NSPIN==1)THEN
        ANN(NA1+1:NA2,NA1+1:NA2,:)=ANN(1:NA1,1:NA1,:)
      ENDIF
      DO NI=1,WH%NIONS
        NA2=GYZ(NI)%CO%DIM2
        R2N=0
        R2N(1:NA2,1:NA2)=GYZ(NI)%CO%R2N
        CALL ZA2_TRANS(ANN(1:NA2,1:NA2,NI),NA2,R2N,1) ! Goto "Natural" basis
      ENDDO
      RETURN
!
      END SUBROUTINE PAWSYM_ANN
!
!****************************************************************************
      SUBROUTINE LOCMAP_SYM()
      INTEGER NI
      INTEGER ROTMAP(SYM%NIOND,SYM%NROTK,SYM%NPCELL)
!
      ROTMAP=SYM%ROTMAP
      DEALLOCATE(SYM%ROTMAP)
      ALLOCATE(SYM%ROTMAP(WH%NIONS,SYM%NROTK,SYM%NPCELL)); SYM%ROTMAP=0
      DO NI=1,WH%NIONS
        SYM%ROTMAP(NI,:,:)=ROTMAP(GYZ(NI)%NI0,:,:)
      ENDDO
      RETURN
!
      END SUBROUTINE LOCMAP_SYM
!
!*************************************************************************************
      SUBROUTINE CHK_EIGS_ANN(ANN,STR,JOBZ)
      COMPLEX(gq) ANN(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
      CHARACTER(*) STR
      CHARACTER JOBZ
! LOCAL
      INTEGER NI,NA2,NASO
      REAL(gq) W(WH%NA2MAX)
      CHARACTER(20) FMT
      COMPLEX(gq) TMP(WH%NA2MAX,WH%NA2MAX)
      COMPLEX(gq) ZTMP(WH%NASOMAX,WH%NASOMAX)
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      OPEN(GL%IU,FILE='GUTZ4.OUT',STATUS='REPLACE')
      WRITE(GL%IU,*)WH%NASOMAX,WH%NIONS
      DO NI=1,WH%NIONS
        NA2=GYZ(NI)%CO%DIM2; NASO=GYZ(NI)%CO%DIMSO
        WRITE(FMT,'("(",I0,"F18.14)")')NASO
        TMP(1:NA2,1:NA2)=-ANN(1:NA2,1:NA2,NI)
        CALL ZSPIN_BLK2_TRANS(TMP,NA2,.FALSE.,BND%ISO) ! Orbital-fast <- Spin-fast
        CALL ZHEEV_(JOBZ,'L',TMP(1:NASO,1:NASO),W(1:NASO),NASO)
        W=-W
        WRITE(GL%IO,'(" NI=",I3," EIGEN VALUES OF ",A7,":")')NI,STR
        WRITE(GL%IO,'(14F10.4)')W(1:NASO)
        ZTMP=0; ZTMP(1:NASO,1:NASO)=TMP(1:NASO,1:NASO)
        WRITE(GL%IU,*)ZTMP
      ENDDO
      CLOSE(GL%IU)
      RETURN
!
      END SUBROUTINE CHK_EIGS_ANN
!
!*************************************************************************************
      SUBROUTINE RM_TB_EL0()
      INTEGER NI,ISYM,IKP
      INTEGER NASO,NA2,NBASE,NASOT
      COMPLEX(gq) EL0(BND%NASOTOT,BND%NASOTOT)
!
      EL0=0; NBASE=1
      DO NI=1,WH%NIONS
        NA2=GYZ(NI)%CO%DIM2; NASO=GYZ(NI)%CO%DIMSO
        CALL ZSPIN_BLK2_TRANS(GYZ(NI)%CO%EL0,NA2,.FALSE.,BND%ISO)
        EL0(NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1)=GYZ(NI)%CO%EL0(1:NASO,1:NASO)
        CALL ZSPIN_BLK2_TRANS(GYZ(NI)%CO%EL0,NA2,.TRUE.,BND%ISO)
        NBASE=NBASE+NASO
      ENDDO
      NASOT=BND%NASOTOT
      DO ISYM=SYM%IDI,SYM%IDF; DO IKP=1,KPT%DIML
        BND%HK0(1:NASOT,1:NASOT,ISYM,IKP)=BND%HK0(1:NASOT,1:NASOT,ISYM,IKP)-EL0
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE RM_TB_EL0
!
!*************************************************************************************
      SUBROUTINE OUT_VEC(STR,VEC,IO,MODE)
      INTEGER IO
      CHARACTER(*) STR
      REAL(gq) VEC(WH%NA2TOT)
      INTEGER,OPTIONAL::MODE
! LOCAL
      INTEGER NI,DIM2,NBASE
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      WRITE(GL%IO,'("************",2X,A10,2X,"************")')STR
      NBASE=1
      DO NI=1,WH%NIONS; DIM2=GYZ(NI)%CO%DIM2
        WRITE(GL%IO,'(" NI=",I3," NI_GLOBAL=",I3)')NI,GYZ(NI)%NI0
        CALL WRT_AN(VEC(NBASE:NBASE+DIM2-1),DIM2,IO)
        IF(PRESENT(MODE))THEN
          IF(MODE==1)THEN
            WRITE(GL%IO,'("   SUB_TOT=",F10.6)')SUM(VEC(NBASE:NBASE+DIM2-1))
          ENDIF
        ENDIF
        NBASE=NBASE+DIM2
      ENDDO
      RETURN
!
      END SUBROUTINE OUT_VEC
!
!*************************************************************************************
      SUBROUTINE OUT_VEC_BK(STR,VEC,IO,MODE)
      INTEGER IO
      CHARACTER(10) STR
      REAL(gq) VEC(BND%NASOTOT)
      INTEGER,OPTIONAL::MODE
! LOCAL
      INTEGER NI,DIMSO,NBASE
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      WRITE(GL%IO,'("************",2X,A10,2X,"************")')STR
      NBASE=1
      DO NI=1,WH%NIONS; DIMSO=GYZ(NI)%CO%DIMSO
        WRITE(GL%IO,'(" NI=",I3," NI_GLOBAL=",I3)')NI,GYZ(NI)%NI0
        CALL WRT_AN(VEC(NBASE:NBASE+DIMSO-1),DIMSO,IO)
        IF(PRESENT(MODE))THEN
          IF(MODE==1)THEN
            WRITE(GL%IO,'("   SUB_TOT=",F10.6)')SUM(VEC(NBASE:NBASE+DIMSO-1))
          ENDIF
        ENDIF
        NBASE=NBASE+DIMSO
      ENDDO
      RETURN
!
      END SUBROUTINE OUT_VEC_BK
!
!*************************************************************************************
      SUBROUTINE DOUT_MAT(STR,AM,IO,MODE)
      INTEGER IO
      CHARACTER(*) STR
      REAL(gq) AM(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
      INTEGER,OPTIONAL::MODE
! LOCAL
      INTEGER NI,DIM2,I
      REAL(gq) RES
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      WRITE(GL%IO,'("************",2X,A10,2X,"************")')STR
      DO NI=1,WH%NIONS; DIM2=GYZ(NI)%CO%DIM2
        WRITE(GL%IO,'(" NI=",I3," NI_GLOBAL=",I3)')NI,GYZ(NI)%NI0
        CALL DWRT_ANN(AM(1:DIM2,1:DIM2,NI),DIM2,IO)
        IF(PRESENT(MODE))THEN
          IF(MODE==1)THEN
            RES=0; DO I=1,DIM2; RES=RES+AM(I,I,NI); ENDDO
            WRITE(GL%IO,'("   SUB_TOT=",2F10.6)')RES
          ELSEIF(MODE==-1)THEN
            WRITE(GL%IO,'("   ORBITAL SPLITTING:")')
            WRITE(GL%IO,'(14F10.5)')(REAL(AM(I,I,NI)-AM(1,1,NI),gq),I=1,DIM2)
          ENDIF
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE DOUT_MAT
!
!*************************************************************************************
      SUBROUTINE ZOUT_MAT(STR,AM,IO,MODE)
      INTEGER IO
      CHARACTER(*) STR
      COMPLEX(gq) AM(WH%NA2MAX,WH%NA2MAX,WH%NIONS)
      INTEGER,OPTIONAL::MODE
! LOCAL
      INTEGER NI,DIM2,I
      COMPLEX(gq) RES
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      WRITE(GL%IO,'("************",2X,A10,2X,"************")')STR
      DO NI=1,WH%NIONS; DIM2=GYZ(NI)%CO%DIM2
        WRITE(GL%IO,'(" NI=",I3," NI_GLOBAL=",I3)')NI,GYZ(NI)%NI0
        CALL ZWRT_ANN(AM(1:DIM2,1:DIM2,NI),DIM2,IO)
        IF(PRESENT(MODE))THEN
          IF(MODE==1)THEN
            RES=0; DO I=1,DIM2; RES=RES+AM(I,I,NI); ENDDO
            WRITE(GL%IO,'("   SUB_TOT=",2F10.6)')RES
          ELSEIF(MODE==-1)THEN
            WRITE(GL%IO,'("   ORBITAL SPLITTING:")')
            WRITE(GL%IO,'(14F10.5)')(REAL(AM(I,I,NI)-AM(1,1,NI),gq),I=1,DIM2)
          ENDIF
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE ZOUT_MAT
!
!*************************************************************************************
      SUBROUTINE DOUT_MAT_BK(STR,AM,IO,MODE)
      INTEGER IO
      CHARACTER(*) STR
      REAL(gq) AM(BND%NASOTOT,BND%NASOTOT)
      INTEGER,OPTIONAL::MODE
! LOCAL
      INTEGER NI,DIMSO,NBASE,I
      REAL(gq) RES
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      WRITE(GL%IO,'("************",2X,A10,2X,"************")')STR
      NBASE=0
      DO NI=1,WH%NIONS; DIMSO=GYZ(NI)%CO%DIMSO
        WRITE(GL%IO,'(" NI=",I3," NI_GLOBAL=",I3)')NI,GYZ(NI)%NI0
        CALL DWRT_ANN(AM(NBASE+1:NBASE+DIMSO,NBASE+1:NBASE+DIMSO),DIMSO,IO)
        IF(PRESENT(MODE))THEN
          IF(MODE==1)THEN
            RES=0; DO I=1,DIMSO; RES=RES+AM(NBASE+I,NBASE+I); ENDDO
            WRITE(GL%IO,'("   SUB_TOT=",2F10.6)')RES
          ELSEIF(MODE==-1)THEN
            WRITE(GL%IO,'("   ORBITAL SPLITTING:")')
            WRITE(GL%IO,'(14F10.5)')(REAL(AM(NBASE+I,NBASE+I)-AM(NBASE+1,NBASE+1),gq),I=1,DIMSO)
          ENDIF
        ENDIF
        NBASE=NBASE+DIMSO
      ENDDO
      RETURN
!
      END SUBROUTINE DOUT_MAT_BK
!
!*************************************************************************************
      SUBROUTINE ZOUT_MAT_BK(STR,AM,IO,MODE)
      INTEGER IO
      CHARACTER(*) STR
      COMPLEX(gq) AM(BND%NASOTOT,BND%NASOTOT)
      INTEGER,OPTIONAL::MODE
! LOCAL
      INTEGER NI,DIMSO,NBASE,I
      COMPLEX(gq) RES
!
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      WRITE(GL%IO,'("************",2X,A10,2X,"************")')STR
      NBASE=0
      DO NI=1,WH%NIONS; DIMSO=GYZ(NI)%CO%DIMSO
        WRITE(GL%IO,'(" NI=",I3," NI_GLOBAL=",I3)')NI,GYZ(NI)%NI0
        CALL ZWRT_ANN(AM(NBASE+1:NBASE+DIMSO,NBASE+1:NBASE+DIMSO),DIMSO,IO)
        IF(PRESENT(MODE))THEN
          IF(MODE==1)THEN
            RES=0; DO I=1,DIMSO; RES=RES+AM(NBASE+I,NBASE+I); ENDDO
            WRITE(GL%IO,'("   SUB_TOT=",2F10.6)')RES
          ELSEIF(MODE==-1)THEN
            WRITE(GL%IO,'("   ORBITAL SPLITTING:")')
            WRITE(GL%IO,'(14F10.5)')(REAL(AM(NBASE+I,NBASE+I)-AM(NBASE+1,NBASE+1),gq),I=1,DIMSO)
          ENDIF
        ENDIF
        NBASE=NBASE+DIMSO
      ENDDO
      RETURN
!
      END SUBROUTINE ZOUT_MAT_BK
!
!*************************************************************************************
      SUBROUTINE CALC_BAND_ALL(LDIAG)
      LOGICAL LDIAG
      INTEGER NKPT,ISYMI,ISYMF,ISP
      INTEGER ISYM,IVEC,IKS,IKP,IKPL,NKP,NBANDS
      COMPLEX(gq),POINTER :: HK(:,:),HK0(:,:),HK1(:,:)
      REAL(gq),POINTER    :: EK(:)
      COMPLEX(gq) :: R(BND%NASOTOT,BND%NASOTOT),LA1(BND%NASOTOT,BND%NASOTOT)
      REAL :: TA1,TA2,TB1,TB2; INTEGER TIB1,TIB2,TIRATE
!
      CALL CPU_TIME(TA1); CALL SYSTEM_CLOCK(TIB1,TIRATE); TB1=REAL(TIB1,4)/REAL(TIRATE,4)
      NKPT=KPT%DIM; ISYMI=SYM%IDI; ISYMF=SYM%IDF
      BND%EK=0
!
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      DO ISYM=ISYMI,ISYMF
      NBANDS=BND%NE(3,IKP)-BND%NE(2,IKP)+1
      HK0=>BND%HK0(1:NBANDS,1:NBANDS,ISYM,IKPL)
      IF(GL%LBNDU==1)THEN
        HK1=>BND%HK1(1:NBANDS,1:NBANDS,ISYM,IKPL)
      ENDIF
      IF(ISYM/=SYM%IE)THEN
        ALLOCATE(EK(NBANDS)); EK=0
      ENDIF
      DO ISP=1,BND%NSPIN
      HK =>BND%VK(1:NBANDS,1:NBANDS,ISYM,IKPL,ISP)
      IF(ISYM==SYM%IE)THEN
        BND%EK(:,IKP,ISP)=BND%EK0(:,IKP)
        EK=>BND%EK(BND%NE(2,IKP):BND%NE(3,IKP),IKP,ISP)
      ENDIF
      R  =BND%R  (1:BND%NASOTOT,1:BND%NASOTOT,ISP)
      LA1=BND%LA1(1:BND%NASOTOT,1:BND%NASOTOT,ISP)+BND%ETA(1:BND%NASOTOT,1:BND%NASOTOT,ISP)
      IF(GL%LBNDU==0)THEN
        CALL CALC_BAND_1K(HK0,EK,HK,R,LA1,NBANDS,BND%NASOTOT,LDIAG)
      ELSE
        CALL CALC_BAND_1K(HK0,EK,HK,R,LA1,NBANDS,BND%NASOTOT,LDIAG,HK1=HK1)
      ENDIF
      ENDDO ! ISP
      IF(ISYM/=SYM%IE)DEALLOCATE(EK)
      ENDDO ! ISYM
      ENDDO; ENDDO
!
      NULLIFY(HK,HK0,EK)
      CALL DSUM_ALL_MPI(BND%EK,BND%NMAX*KPT%DIM*BND%NSPIN)
      CALL CPU_TIME(TA2); CALL SYSTEM_CLOCK(TIB2,TIRATE); TB2=REAL(TIB2,4)/REAL(TIRATE,4)
      CALL OUT_TIME_USE('CALC_BAND_AL',TA2-TA1,TB2-TB1,GL%IO)
      RETURN
!
      END SUBROUTINE CALC_BAND_ALL
!
!*************************************************************************************
      SUBROUTINE CALC_BAND_1K(HK0,EK,HK,R,LA1,NBANDS,NASOT,LDIAG,HK1)
      LOGICAL LDIAG
      INTEGER NBANDS,NASOT
      REAL(gq) EK(NBANDS)
      COMPLEX(gq) R(NASOT,NASOT),LA1(NASOT,NASOT)
      COMPLEX(gq) HK(NBANDS,NBANDS),HK0(NBANDS,NBANDS)
      COMPLEX(gq),OPTIONAL :: HK1(NBANDS,NBANDS)
!
      CALL CALC_HAMIL_1K(HK0,HK,R,LA1,NBANDS,NASOT)
      IF(PRESENT(HK1)) HK=HK+HK1
      EK=0
      IF(LDIAG) CALL ZHEEV_('V','U',HK,EK,NBANDS)
      RETURN
!
      END SUBROUTINE CALC_BAND_1K
!
!*************************************************************************************
      SUBROUTINE CALC_HAMIL_1K(HK0,HK,R,LA1,NBANDS,NASOT)
      INTEGER NBANDS,NASOT
      COMPLEX(gq) HK(NBANDS,NBANDS),HK0(NBANDS,NBANDS)
      COMPLEX(gq) R(NASOT,NASOT),LA1(NASOT,NASOT)
! Local
      INTEGER IA
!
      HK=HK0
      CALL ZANNXBNM('N',R,HK(1:NASOT,:),NASOT,NBANDS)
      CALL ZANMXBMM('C',HK(:,1:NASOT),R,NBANDS,NASOT) ! RHR^+
      HK(1:NASOT,1:NASOT)=HK(1:NASOT,1:NASOT)+LA1
      RETURN
!
      END SUBROUTINE CALC_HAMIL_1K
!
!=============================================================================
! Possibly non-orthogonal projector analysis
!=============================================================================
      SUBROUTINE PRE_ANALYSIS()
!
      CALL CALC_LOC_OVLP()
      IF(GL%LBNDU==1)THEN
        CALL CALC_LOC_HR2()
      ELSE
        CALL CALC_LOC_HR()
      ENDIF
      CALL GUTZ_FERMI(GL%IO)
      CALL CALC_NKS_NO()
      RETURN
!
      END SUBROUTINE PRE_ANALYSIS
!
!=============================================================================
      SUBROUTINE CALC_LOC_OVLP()
      INTEGER IVEC,IKS,IKP,IKPL,NKP,ISYM,NBANDS,NSYM,I
      REAL(gq) WTK
      COMPLEX(gq) SR0(BND%NASOTOT,BND%NASOTOT)
!
      SR0=0
      NSYM=SYM%IDF-SYM%IDI+1
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      NBANDS=BND%NE(3,IKP)-BND%NE(2,IKP)+1
      WTK=KPT%WT(IKP)/NSYM
      DO ISYM=SYM%IDI,SYM%IDF
        CALL ZGEMM('C','N',BND%NASOTOT,BND%NASOTOT,NBANDS,Z1, &
                  &BND%UK(1:NBANDS,:,ISYM,IKPL),NBANDS, &
                  &BND%UK(1:NBANDS,:,ISYM,IKPL),NBANDS, &
                  &Z0,BND%SK(:,:,ISYM,IKPL),BND%NASOTOT)
        SR0=SR0+BND%SK(:,:,ISYM,IKPL)*WTK
      ENDDO ! ISYM
      ENDDO; ENDDO
      CALL ZSUM_ALL_MPI(SR0,BND%NASOTOT*BND%NASOTOT)
      CALL ZOUT_MAT_BK('SR0',SR0,GL%IO)
      RETURN
!
      END SUBROUTINE CALC_LOC_OVLP
!
!=============================================================================
      SUBROUTINE ROTATE_BNDU(TRANSB,MODE)
      INTEGER MODE
      CHARACTER*1 TRANSB
! LOCAL
      INTEGER IVEC,IKS,IKP,IKPL,NKP,ISYM,NBANDS,I,NI,NBASE,NASO,NA2
!
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      NBANDS=BND%NE(3,IKP)-BND%NE(2,IKP)+1
      DO ISYM=SYM%IDI,SYM%IDF
      NBASE=1
      DO NI=1,WH%NIONS
        NASO=GYZ(NI)%CO%DIMSO; NA2=GYZ(NI)%CO%DIM2
        IF(MODE==0)THEN
          CALL ZANMXBMM(TRANSB,BND%UK(1:NBANDS,NBASE:NBASE-1+NASO,ISYM,IKPL),GYZ(NI)%CO%B2N(1:NASO,1:NASO),NBANDS,NASO)
        ELSE
          CALL ZSPIN_BLK2U_TRANS(GYZ(NI)%CO%N2N,NA2,NA2,.FALSE.,BND%ISO,LURIGHT=.TRUE.) ! orbital faster
          CALL ZANMXBMM(TRANSB,BND%UK(1:NBANDS,NBASE:NBASE-1+NASO,ISYM,IKPL),GYZ(NI)%CO%N2N(1:NASO,1:NASO),NBANDS,NASO)
          CALL ZSPIN_BLK2U_TRANS(GYZ(NI)%CO%N2N,NA2,NA2,.TRUE.,BND%ISO,LURIGHT=.TRUE.) ! spin faster
        ENDIF
      ENDDO; ENDDO
      ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE ROTATE_BNDU
!
!=============================================================================
! Calcuate local one-body using eigen-value and eigen0vecotr of the Hamilt.
!=============================================================================
      SUBROUTINE CALC_LOC_HR()
      INTEGER IVEC,IKS,IKP,IKPL,NKP,ISYM,NBANDS,NSYM,I
      REAL(gq) WTK
      COMPLEX(gq) :: HR0(BND%NASOTOT,BND%NASOTOT)
      COMPLEX(gq) :: EV(BND%NMAXIN,BND%NASOTOT),HKL(BND%NASOTOT,BND%NASOTOT)
!
      HR0=0
      NSYM=SYM%IDF-SYM%IDI+1
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      NBANDS=BND%NE(3,IKP)-BND%NE(2,IKP)+1
      WTK=KPT%WT(IKP)/NSYM
      DO ISYM=SYM%IDI,SYM%IDF
        DO I=1,NBANDS; EV(I,:)=BND%EK0(BND%NE(2,IKP)-1+I,IKP)*BND%UK(I,:,ISYM,IKPL); ENDDO
        CALL ZGEMM('C','N',BND%NASOTOT,BND%NASOTOT,NBANDS,Z1,&
                  &BND%UK(1:NBANDS,:,ISYM,IKPL),NBANDS, &
                  &EV(1:NBANDS,:),NBANDS,Z0,HKL,BND%NASOTOT)
        HR0=HR0+HKL*WTK
      ENDDO ! ISYM
      ENDDO; ENDDO
      CALL ZSUM_ALL_MPI(HR0,BND%NASOTOT*BND%NASOTOT)
      CALL ZOUT_MAT_BK('HR0',HR0,GL%IO,-1)
      RETURN
!
      END SUBROUTINE CALC_LOC_HR
!
!=============================================================================
! Calcuate local one-body using Hamiltonian matrix
!=============================================================================
      SUBROUTINE CALC_LOC_HR2()
      INTEGER IVEC,IKS,IKP,IKPL,NKP,ISYM,NBANDS,NSYM,I
      REAL(gq) WTK
      COMPLEX(gq) :: HR0(BND%NASOTOT,BND%NASOTOT)
      COMPLEX(gq) :: HKL(BND%NASOTOT,BND%NASOTOT)
!
      HR0=0
      NSYM=SYM%IDF-SYM%IDI+1
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      NBANDS=BND%NE(3,IKP)-BND%NE(2,IKP)+1
      WTK=KPT%WT(IKP)/NSYM
      DO ISYM=SYM%IDI,SYM%IDF
        CALL ZUHAU(BND%HK0(1:NBANDS,1:NBANDS,ISYM,IKPL),BND%UK(1:NBANDS,:,ISYM,IKPL),NBANDS,BND%NASOTOT,HKL)
        HR0=HR0+HKL*WTK
      ENDDO ! ISYM
      ENDDO; ENDDO
      CALL ZSUM_ALL_MPI(HR0,BND%NASOTOT*BND%NASOTOT)
      CALL ZOUT_MAT_BK('HR0',HR0,GL%IO,-1)
      RETURN
!
      END SUBROUTINE CALC_LOC_HR2
!
!=============================================================================
      SUBROUTINE CALC_NKS_NO()
      INTEGER IVEC,IKS,IKP,IKPL,NKP,ISYM,NBANDS,NSYM,I
      REAL(gq) WTK
      COMPLEX(gq) :: NKS0(BND%NASOTOT,BND%NASOTOT)
      COMPLEX(gq) :: NKL(BND%NASOTOT,BND%NASOTOT)
!
      NKS0=0
      NSYM=SYM%IDF-SYM%IDI+1
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      NBANDS=BND%NE(3,IKP)-BND%NE(2,IKP)+1
      WTK=1._gq/NSYM/BND%RSPO
      DO ISYM=SYM%IDI,SYM%IDF
        CALL CALC_NKS_NO1K(BND%FERWE(BND%NE(2,IKP):BND%NE(3,IKP),IKP,1), &
                           &BND%UK(1:NBANDS,:,ISYM,IKPL),BND%SK(:,:,ISYM,IKPL), &
                           &NKL,NBANDS,BND%NASOTOT)
        NKS0=NKS0+NKL*WTK
      ENDDO ! ISYM
      ENDDO; ENDDO
      CALL ZSUM_ALL_MPI(NKS0,BND%NASOTOT*BND%NASOTOT)
      CALL ZOUT_MAT_BK('NK0',NKS0,GL%IO,1)
      RETURN
!
      END SUBROUTINE CALC_NKS_NO
!
!=============================================================================
      SUBROUTINE CALC_E_HYBRD()
      INTEGER NI
!
      ENG%HYBRD=SUM(WH%R*(WH%D0))
      RETURN
!
      END SUBROUTINE CALC_E_HYBRD
!
!=============================================================================
!   f_n <psi_n|a'> [S^(-1)]_{a',a} S^(-1)]_{b,b'} <b'|psi_n>
! = S^(-1)]_{b,b'} <b'|psi_n> f_n <psi_n|a'> [S^(-1)]_{a',a}
!=============================================================================
      SUBROUTINE CALC_NKS_NO1K(FERWE,UK,SK,NKL,NBANDS,NAST)
      INTEGER NBANDS,NAST
      REAL(gq) FERWE(NBANDS)
      COMPLEX(gq) UK(NBANDS,NAST),SK(NAST,NAST),NKL(NAST,NAST)
! LOCAL
      INTEGER I
      COMPLEX(gq) :: EV(NBANDS,NAST),SKI(NAST,NAST),ZBUF(NAST,NAST)
!
      DO I=1,NBANDS; EV(I,:)=FERWE(I)*UK(I,:); ENDDO
      CALL ZGEMM('C','N',NAST,NAST,NBANDS,Z1,UK,NBANDS,EV,NBANDS,Z0,NKL,NAST) ! <b'|Psi_n> f_n <psi_n|a'>
      SKI=SK
      CALL ZINV_(SKI,NAST) ! SK^(-1)
      CALL ZGEMM('N','N',NAST,NAST,NAST,Z1,SKI,NAST,NKL,NAST,Z0,ZBUF,NAST) 
      CALL ZGEMM('N','N',NAST,NAST,NAST,Z1,ZBUF,NAST,SKI,NAST,Z0,NKL,NAST)
      NKL=TRANSPOSE(NKL)
      RETURN
!
      END SUBROUTINE CALC_NKS_NO1K
!
!=============================================================================
! ORTHOGONALIZE LOCAL PROJECTOR
!=============================================================================
      SUBROUTINE ORTH_LOC_PROJ()
      INTEGER IVEC,IKS,IKP,IKPL,NKP,ISYM,NBANDS
!
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
      NBANDS=BND%NE(3,IKP)-BND%NE(2,IKP)+1
      DO ISYM=SYM%IDI,SYM%IDF
        CALL ORTH_LOC_PROJ1(BND%UK(1:NBANDS,:,ISYM,IKPL),BND%SK(:,:,ISYM,IKPL),NBANDS,BND%NASOTOT)
      ENDDO ! ISYM
      ENDDO; ENDDO
      DEALLOCATE(BND%SK)
      RETURN
!
      END SUBROUTINE ORTH_LOC_PROJ
!
!=============================================================================
      SUBROUTINE ORTH_LOC_PROJ1(U,SK,NBANDS,NAST)
      INTEGER NBANDS,NAST
      COMPLEX(gq) U(NBANDS,NAST),SK(NAST,NAST)
! LOCAL
      INTEGER J1,J2
      COMPLEX(gq) U_(NBANDS,NAST)
      COMPLEX(gq),ALLOCATABLE :: SNH(:,:)
!
      ALLOCATE(SNH(NAST,NAST)); SNH=0
      CALL ZATOFA(SK,SNH,NAST,-12,D1,.TRUE.)
      U_=U; U=0
      DO J1=1,NAST; DO J2=1,NAST
        U(:,J1)=U(:,J1)+U_(:,J2)*SNH(J2,J1)
      ENDDO; ENDDO
      DEALLOCATE(SNH)
      RETURN
!
      END SUBROUTINE ORTH_LOC_PROJ1
!
!=============================================================================
      SUBROUTINE READ_NELF1()
      INTEGER NI
      COMPLEX(gq) NPHY(WH%NA2MAX,WH%NA2MAX)
!
      IF(GL%LDC/=2.AND.GL%LDC/=12.AND.GL%LDC/=-12)RETURN
      OPEN(GL%IU,FILE="GL_NELF1.INP",STATUS='OLD',ERR=100)
      IF(GL%LDC==2.OR.GL%LDC==12)THEN
        READ(GL%IU,*) GYZ(:)%CO%NELF1
      ELSE
        DO NI=1,WH%NIONS
          READ(GL%IU,*)NPHY
          WH%NPHY_FIX(:,:,NI)=NPHY
        ENDDO
        CALL ZOUT_MAT('NPHY_FIX',WH%NPHY_FIX,GL%IO,1)
      ENDIF
      CLOSE(GL%IU)
100   CONTINUE
      IF(GL%LDC==2.OR.GL%LDC==12)THEN
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ATOM DEPENDENT NELF1:")')
        IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(5F12.4)')GYZ(:)%CO%NELF1
      ENDIF
      RETURN
!
      END SUBROUTINE READ_NELF1
!
!=============================================================================
      SUBROUTINE UPDATE_NELF1()
      INTEGER NI
      REAL(gq) NDIFF(WH%NIONS)
      COMPLEX(gq) NDIFF2(WH%NA2MAX*WH%NA2MAX*WH%NIONS),DIFF(1)
      COMPLEX(gq) NPHY(WH%NA2MAX,WH%NA2MAX)
      CHARACTER(40) FMT
!
      WRITE(FMT,'(A,I3,A)')"(",WH%NA2MAX,'("(",F18.12,",",F18.12,")"))'
      IF(GL%LDC/=2.AND.GL%LDC/=12.AND.GL%LDC/=-12)RETURN
      IF(GL%LDC==2.OR.GL%LDC==12)THEN
        NDIFF=GYZ(:)%CO%NET-GYZ(:)%CO%NELF1
        DIFF=NDIFF(MAXLOC(ABS(NDIFF)))
      ELSE
        NDIFF2=RESHAPE(WH%NC_PHY-WH%NPHY_FIX,(/WH%NA2MAX*WH%NA2MAX*WH%NIONS/))
        DIFF=NDIFF2(MAXLOC(ABS(NDIFF2)))
      ENDIF
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" MAX NELF1 DIFF=",2F14.8)')DIFF
      IF(GL%LDC/=12.AND.GL%LDC/=-12)RETURN
      IF(GL%LDC==12)THEN
        GYZ(:)%CO%NELF1=GYZ(:)%CO%NELF1+GL%DCMIX_A*NDIFF
      ELSE
        WH%NPHY_FIX=WH%NPHY_FIX+GL%DCMIX_A*RESHAPE(NDIFF2,(/WH%NA2MAX,WH%NA2MAX,WH%NIONS/))
      ENDIF
      IF(GP%MYRANK.NE.GP%MASTER)RETURN
      OPEN(GL%IU,FILE="GL_NELF1.OUT",STATUS='REPLACE')
      IF(GL%LDC==12)THEN
        WRITE(GL%IU,*) GYZ(:)%CO%NELF1
      ELSE
        DO NI=1,WH%NIONS
          NPHY=WH%NPHY_FIX(:,:,NI)
          WRITE(GL%IU,FMT)NPHY
        ENDDO
      ENDIF
      CLOSE(GL%IU)
      RETURN
!
      END SUBROUTINE UPDATE_NELF1
!
!=============================================================================
      SUBROUTINE SET_DMFTU_BNDU(U,NBANDS,NASOMAX,NIONS,ISYM,IKP)
      INTEGER NBANDS,NASOMAX,NIONS,ISYM,IKP
      COMPLEX(gq) U(NBANDS,NASOMAX,NIONS)
! LOCAL
      INTEGER NI,NASO,ISYM_
!
      ISYM_=ISYM-SYM%IDI+1
      DO NI=1,WH%NIONS
        NASO=GYZ(NI)%CO%DIMSO
        GYZ(NI)%CO%UK(1:NBANDS,:,ISYM_,IKP)=U(:,1:NASO,NI)
      ENDDO
      RETURN
!
      END SUBROUTINE SET_DMFTU_BNDU
!
!****************************************************************************
! Calculate bare/renormalized occupation matrix of the original KS orbitals 
!****************************************************************************
      SUBROUTINE CALC_KSWT(MODE)
      INTEGER MODE
! LOCAL
      INTEGER IVEC,IKS,IKP,IKPL,NKP,NBANDS,I,ISYM,ISP
      INTEGER NEMIN,NEMAX
      REAL(gq) SUMWT(2),SUM1,SUM2,WTK,WTK0,MAXOFFDIAG
      COMPLEX(gq),POINTER :: VK(:,:),KSWT(:,:),UK(:,:)
      REAL(gq),POINTER    :: FERWE(:)
!
      ! "CALC_KSWT"
      IF(MODE==0) ALLOCATE(BND%CWT(BND%NMAXIN,BND%NMAXIN,KPT%DIML))
      BND%CWT=0; SUMWT=0; MAXOFFDIAG=0
!
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
        NEMIN=BND%NE(2,IKP); NEMAX=BND%NE(3,IKP)
        NBANDS=NEMAX-NEMIN+1
        KSWT =>BND%CWT  (1:NBANDS,1:NBANDS,IKPL)
        WTK=1._gq/BND%RSPO; WTK0=KPT%WT(IKP)
        DO ISP=1,BND%NSPIN
        FERWE=>BND%FERWE(NEMIN:NEMAX,IKP,ISP)
        ISYM=SYM%IE
        VK   =>BND%VK   (1:NBANDS,1:NBANDS,ISYM,IKPL,ISP)
        UK   =>BND%UK0  (1:NBANDS,1:NBANDS,ISYM,IKPL)
        CALL CALC_KSWT_1K(KSWT,VK,UK,FERWE,NBANDS,WTK,WTK0,ISP,MODE)
        SUM1=SUM(FERWE)
        SUMWT(1)=SUMWT(1)+SUM1
        ENDDO ! ISP
        KSWT=KSWT*BND%RSPO
        SUM2=0
        DO I=1,NBANDS; SUM2=SUM2+KSWT(I,I); ENDDO
        SUMWT(2)=SUMWT(2)+SUM2
        DO I=1,NBANDS-1; MAXOFFDIAG=MAX(MAXOFFDIAG,MAXVAL(ABS(KSWT(I+1:,I)))); ENDDO
        IF(BND%ISO>BND%ISO_IN)THEN
          CALL ZSPIN_BLK2_TRANS(BND%CWT(1:NBANDS,1:NBANDS,IKPL),NBANDS,.FALSE.,BND%ISO_IN)
          BND%CWT(1:NBANDS,1:NBANDS,IKPL)=BND%CWT(1:NBANDS,1:NBANDS,IKPL)*BND%ISO/BND%ISO_IN
        ENDIF
      ENDDO; ENDDO
      CALL DSUM_ALL_MPI(SUMWT,2)
      CALL DMAX1_ALL_MPI(MAXOFFDIAG)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("MODE=",I2," CORRELATED SUBSPACE, SUM_FERWT=",F15.7," SUM_KSWT=",F15.7)')MODE,SUMWT
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("      MAX OFF DIAGONAL KS.WT=",F12.6)')MAXOFFDIAG
      IF(ABS(SUMWT(1)-SUMWT(2)).GT.1.E-3_gq)STOP ' ERROR: SUM_FERWT \= SUM_KSWT!'
      NULLIFY(VK,FERWE,KSWT)
      ! "CALC_KSWT"
      RETURN
!
      END SUBROUTINE CALC_KSWT
!
!****************************************************************************
! Density matrix <a|psi>f_psi<psi|b>
!****************************************************************************
      SUBROUTINE CALC_KSWT_1K(KSWT,VK,UK,FERWE,NBANDS,WTK,WTK0,ISP,MODE)
      INTEGER NBANDS,ISP,MODE
      REAL(gq) FERWE(NBANDS),WTK,WTK0
      COMPLEX(gq) VK(NBANDS,NBANDS),KSWT(NBANDS,NBANDS),UK(NBANDS,NBANDS)
! LOCAL
      COMPLEX(gq),ALLOCATABLE :: NABR(:,:)
!
      ALLOCATE(NABR(NBANDS,NBANDS)); NABR=0
      IF(MODE==0)THEN
        CALL CALC_NABR_1K(NABR,VK,FERWE,NBANDS,WTK,WTK0,ISP) ! <psi|a><b|psi>
        NABR=TRANSPOSE(NABR) ! <a|psi><psi|b>
      ELSE
        CALL CALC_NAB_1K(NABR,VK,FERWE,NBANDS,WTK)
      ENDIF
      CALL ZUHAU(NABR,UK,NBANDS,NBANDS,TRUL='N',TRUR='C')
      KSWT=KSWT+NABR ! Density matrix <a|psi>f_psi<psi|b>
      DEALLOCATE(NABR)
      RETURN
!
      END SUBROUTINE CALC_KSWT_1K
!
!****************************************************************************
! FOR WIEN_DMFT2
!****************************************************************************
      SUBROUTINE WRT_KSWT()
      REAL(8) VNORM1,FAK
      INTEGER IVEC,IKP,IKPL,IKS,NKP,NBANDS,I,IU,ISYM,ID
      INTEGER NEMIN,NEMAX
      COMPLEX(8),POINTER :: KSWT(:,:)
!
      VNORM1=1.d0; IU=GL%IU
      IF(BND%ISO_IN.EQ.2.AND.BND%ISPIN_IN.EQ.1)VNORM1=.5d0
!
      IKPL=0
      DO IVEC=1,GP%NVEC
      IF(GP%LKPVEC)THEN; ID=GP%KVEC(IVEC,1); ELSE; ID=GP%MYRANK; ENDIF
      OPEN(IU,FILE=TRIM(ADJUSTL(FILE_NAME(ID,0,-2))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      WRITE(IU)BND%EF
      WRITE(IU)ENG%DBLC+ENG%GAMM
      WRITE(IU)ENG%BAND
      IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF
      DO IKS=1,NKP
        IF(GP%LKPVEC)THEN
          IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1
        ELSE
          IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML
        ENDIF
        IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
        NEMIN=(BND%NE(2,IKP)-1)*BND%ISO_IN/BND%ISO+1
        NEMAX=BND%NE(3,IKP)*BND%ISO_IN/BND%ISO
        NBANDS =NEMAX-NEMIN+1
        KSWT=>BND%CWT(1:NBANDS,1:NBANDS,IKPL)
        FAK=1/VNORM1
        KSWT=KSWT*FAK ! To be consistent with dmft2
        WRITE(IU)NEMIN,NEMAX
        WRITE(IU)KSWT
      ENDDO; ENDDO
      CLOSE(IU)
      NULLIFY(KSWT)
      RETURN
!
      END SUBROUTINE WRT_KSWT
!
!*************************************************************
! FOR VASP
!*************************************************************
      SUBROUTINE DIAG_KSWT()
      INTEGER IVEC,IKP,IKPL,IKS,NKP,NBANDS,I,IU,ISYM,ID
      INTEGER NEMIN,NEMAX
      COMPLEX(8),POINTER :: KSWT(:,:)
!
      BND%FERWER=0
      IKPL=0; DO IVEC=1,GP%NVEC; IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF; DO IKS=1,NKP;   IF(GP%LKPVEC)THEN; IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1; ELSE; IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF;   IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
        NEMIN=(BND%NE(2,IKP)-1)*BND%ISO_IN/BND%ISO+1
        NEMAX=BND%NE(3,IKP)*BND%ISO_IN/BND%ISO
        NBANDS=NEMAX-NEMIN+1
        BND%FERWER(1:NEMIN-1,IKP,1)=1._gq
        KSWT =>BND%CWT(1:NBANDS,1:NBANDS,IKPL)
        KSWT=-KSWT ! Eigen-value correct order
        CALL ZHEEV_('V','L',KSWT,BND%FERWER(NEMIN:NEMAX,IKP,1),NBANDS)
        BND%FERWER(NEMIN:NEMAX,IKP,1)=-BND%FERWER(NEMIN:NEMAX,IKP,1)/BND%RSPO/KPT%WT(IKP)
      ENDDO; ENDDO
      NULLIFY(KSWT)
      RETURN
!
      END SUBROUTINE DIAG_KSWT
!
!****************************************************************************
! FOR VASP
!****************************************************************************
      SUBROUTINE WRT0_KSWT(MODE)
      INTEGER MODE
! LOCAL
      INTEGER IVEC,IKP,IKPL,IKS,NKP,NBANDS,I,IU,ISYM,ID
      INTEGER NEMIN,NEMAX,NBTOT
      COMPLEX(8),POINTER :: KSWT(:,:)
!
      IKPL=0
      DO IVEC=1,GP%NVEC
      IF(GP%LKPVEC)THEN; ID=GP%KVEC(IVEC,1); ELSE; ID=GP%MYRANK; ENDIF
      OPEN(IU,FILE=TRIM(ADJUSTL(FILE_NAME(ID,0,MODE))),STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      WRITE(IU)ENG%TB
      IF(GP%LKPVEC)THEN; NKP=GP%KVEC(IVEC,2); ELSE; NKP=KPT%DIM; ENDIF
      DO IKS=1,NKP
        IF(GP%LKPVEC)THEN
          IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1
        ELSE
          IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML
        ENDIF
        IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
        NEMIN=(BND%NE(2,IKP)-1)*BND%ISO_IN/BND%ISO+1
        NEMAX=BND%NE(3,IKP)*BND%ISO_IN/BND%ISO
        NBTOT=BND%NE(1,IKP)*BND%ISO_IN/BND%ISO
        NBANDS=NEMAX-NEMIN+1
        KSWT=>BND%CWT(1:NBANDS,1:NBANDS,IKPL)
        WRITE(IU)NBTOT,NEMIN,NEMAX
        WRITE(IU)BND%FERWER(1:NBTOT,IKP,1)
        WRITE(IU)KSWT
      ENDDO; ENDDO
      CLOSE(IU)
      NULLIFY(KSWT)
      RETURN
!
      END SUBROUTINE WRT0_KSWT
!
!*************************************************************
      SUBROUTINE SET_NIL()
      INTEGER NI,NIMAP,NIS,NT
!
      NIS=0; NT=0; GYZ(:)%NIL=0
      DO NI=1,WH%NIONS
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.EQ.NI)THEN
          IF(GP%MYRANK.EQ.MOD(NIS,GP%NPROCS))GYZ(NI)%NIL=NI
          NIS=NIS+1
          IF(NT/=GYZ(NI)%NT)THEN
            NT=GYZ(NI)%NT
          ENDIF
        ELSEIF(NIMAP.LT.NI)THEN
          IF(GYZ(NIMAP)%NIL.GT.0)GYZ(NI)%NIL=NI
        ELSE
          WRITE(0,'(" ERROR: NIMAP>NI IN GL.INP!")'); STOP
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE SET_NIL
!
!*************************************************************
      SUBROUTINE ROT_RLA1()
      INTEGER NI,NIMAP
!
      DO NI=1,WH%NIONS
        NIMAP=GYZ(NI)%NIMAP
        IF(NIMAP.EQ.NI)THEN
          CALL ROT_CO_RLA1(GYZ(NI)%CO)
        ELSE
          GYZ(NI)%CO%R  =GYZ(NIMAP)%CO%R
          GYZ(NI)%CO%LA1=GYZ(NIMAP)%CO%LA1
        ENDIF
      ENDDO
      RETURN
!
      END SUBROUTINE ROT_RLA1
!
!*************************************************************
      SUBROUTINE ORDER_COSYM_IDSUM()
      INTEGER NI,ISUM,NIP
      INTEGER ID(WH%NIONS)
!
      ID=0; ISUM=1; ID(1)=ISUM
      DO NI=2,WH%NIONS
        CALL LOCATE_MI(GYZ(1:NI-1)%CO%SYMG%IDSUM,NI-1,GYZ(NI)%CO%SYMG%IDSUM,NIP)
        IF(NIP>0)THEN
          ID(NI)=ID(NIP)
        ELSE
          ISUM=ISUM+1
          ID(NI)=ISUM
        ENDIF
      ENDDO
      GYZ(:)%CO%SYMG%IDSUM=ID
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GYZ%SYM%IDSUM:")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(20I3)')ID
      RETURN
!
      END SUBROUTINE ORDER_COSYM_IDSUM
!
!*************************************************************
      SUBROUTINE MAP_WH_BND_R(X,Y,LBACK)
      COMPLEX(gq) X(WH%NA2MAX,WH%NA2MAX,WH%NIONS),Y(WH%NASOTOT,WH%NASOTOT,BND%NSPIN)
      LOGICAL LBACK
! LOCAL
      INTEGER NI,ISP
      INTEGER NA2,NASO,NBASE,NSPIN
      COMPLEX(gq) BUF(GL%NA2MAX,GL%NA2MAX)
!
      NBASE=1; NSPIN=BND%NSPIN
      DO NI=1,WH%NIONS
        NA2=GYZ(NI)%CO%DIM2; NASO=GYZ(NI)%CO%DIMSO; BUF=0
        IF(.NOT.LBACK)THEN
          BUF(1:NA2,1:NA2)=X(1:NA2,1:NA2,NI)
          CALL ZSPIN_BLK2_TRANS(BUF(1:NA2,1:NA2),NA2,.FALSE.,BND%ISO)
          DO ISP=1,NSPIN
            Y(NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,ISP)=BUF((ISP-1)*NASO+1:ISP*NASO,(ISP-1)*NASO+1:ISP*NASO)
          ENDDO
        ELSE
          DO ISP=1,NSPIN
            BUF((ISP-1)*NASO+1:ISP*NASO,(ISP-1)*NASO+1:ISP*NASO)=Y(NBASE:NBASE+NASO-1,NBASE:NBASE+NASO-1,ISP)
          ENDDO
          IF(NSPIN==1)BUF(1+NASO:NA2,1+NASO:NA2)=BUF(1:NASO,1:NASO)
          CALL ZSPIN_BLK2_TRANS(BUF(1:NA2,1:NA2),NA2,.TRUE.,BND%ISO)
          X(1:NA2,1:NA2,NI)=BUF(1:NA2,1:NA2)
        ENDIF
        NBASE=NBASE+NASO
      ENDDO
      RETURN
!
      END SUBROUTINE MAP_WH_BND_R
!
!
      END MODULE GUTZ
!
!*************************************************************
      SUBROUTINE GUTZ_FCN1(N,X,FVEC,IFLAG)
      USE gprec; USE GUTZ
      IMPLICIT NONE
      INTEGER N,IFLAG
      REAL(gq) X(N),FVEC(N)
! LOCAL
      REAL(gq) MAXERR
      INTEGER I,NX1,NX2
      COMPLEX(gq) WHX(WH%DIMX2)
      REAL :: TA1,TA2,TB1,TB2; INTEGER TIB1,TIB2,TIRATE
      REAL(gq),PARAMETER::RCUT=1.E-8_gq
!
      ! 'GUTZ_FCN1'
      CALL CPU_TIME(TA1); CALL SYSTEM_CLOCK(TIB1,TIRATE); TB1=REAL(TIB1,4)/REAL(TIRATE,4)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("************ WH%X ************")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(5F14.8)')X
      GL%ITER1=GL%ITER1+1
      NX1=SUM(WH%LXREAL(1:WH%DIMX1)); NX2=SUM(WH%LXREAL(1+WH%DIMX2:WH%DIMX3))
      IF(GL%LSOLVER==2)THEN
        CALL SET_XTOWHX(X,WH%X(WH%DIMX1+1:WH%DIMX2),WH%LXREAL(WH%DIMX1+1:WH%DIMX2),N,WH%DIMXH,.FALSE.)
      ELSEIF(GL%LSOLVER==3.OR.GL%LSOLVER==30)THEN
        CALL SET_XTOWHX(X(1:NX1),WH%X(1:WH%DIMX1),WH%LXREAL(1:WH%DIMX1),NX1,WH%DIMXG,.FALSE.)
        CALL SET_XTOWHX(X(1+NX1:N),WH%X(1+WH%DIMX2:WH%DIMX3),WH%LXREAL(1+WH%DIMX2:WH%DIMX3),NX2,WH%DIMXH,.FALSE.)
        WHX(1:WH%DIMX1)=WH%X(1:WH%DIMX1)
      ELSE  ! 1,101,102,103
        CALL SET_XTOWHX(X,WH%X(1:WH%DIMX2),WH%LXREAL(1:WH%DIMX2),N,WH%DIMX2,.FALSE.)
        WHX=WH%X(1:WH%DIMX2)
      ENDIF
      CALL SET_WH_X123(.FALSE.) ! UPDATE R,LA,etc.
      CALL MAP_WH_BND_R(WH%R,BND%R,.FALSE.)
      CALL ZOUT_MAT('R-IN',WH%R,GL%IO)
!
      IF(GL%LETA==0)WH%ETA=0
      IF(GL%LSOLVER==3.OR.GL%LSOLVER==30)THEN
        CALL MAP_WH_BND_R(WH%ETA,BND%ETA,.FALSE.)
        CALL ZOUT_MAT('ETA',WH%ETA,GL%IO)
        CALL GUTZ_SOLVE_LA1()
        IF(WH%DIMXHO>0.AND.GL%LETA==1)THEN
          CALL GUTZ_SOLVE_ETA()
        ENDIF
      ELSE
        CALL MAP_WH_BND_R(WH%LA1,BND%LA1,.FALSE.)
        CALL ZOUT_MAT('LA1',WH%LA1,GL%IO)
        IF((GL%LSCF/=2).AND.WH%DIMXHO>0.AND.GL%LETA==1)THEN
          CALL GUTZ_SOLVE_ETA()
        ELSE
          CALL MAP_WH_BND_R(WH%ETA,BND%ETA,.FALSE.)
          CALL ZOUT_MAT('ETA',WH%ETA,GL%IO)
          IF(GL%LMODEL==0)THEN
            CALL CALC_BAND_ALL(GL%LGREEN==0)
            IF(GL%LPSICOMP.GT.0)THEN
              CALL CALC_PSI_COMP()
              CALL OUT_BANDS()
            ENDIF
            IF(GL%LSCF==2)THEN
              STOP ' BAND STRUCTURE DONE!'
            ENDIF
            IF(GL%LGREEN==0)THEN
              CALL GUTZ_FERMI(GL%IO)
            ELSE
              CALL GREEN_FERMI_NV()
              IF(GL%LGREEN==2)THEN
                CALL CALC_SELF_ENERGY_FULL()
              ENDIF
              CALL SUM_GFK_FULL(BND%EF)
            ENDIF
            CALL CALC_MUP_DN()
            IF(BND%NSPIN==2)THEN
              IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" TOTAL MAGNETIC MOMENT: ",F8.3)')BND%MUP_DN
            ENDIF
          ELSE
            CALL CALC_GF_IMP_ALL()
          ENDIF
          IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GA FERMI LEVEL=",F16.8)')BND%EF
          CALL CALC_NKS(GL%LGREEN)
          CALL SET_WH_X1(WH%X(WH%DIMX2+1:WH%DIMX3),WH%NKS,WH%NA2MAX,.TRUE.,1)
        ENDIF
      ENDIF
!
      CALL REGULARIZE_NKS(0)
      CALL CALC_ISIMIX(0)
      CALL ZOUT_MAT('R-IN',WH%R,GL%IO)
!
      CALL MAP_WH_BND_R(WH%NKS,BND%NKS,.FALSE.)
      CALL CALC_DA()
!
      CALL CALC_VDC() ! Calculate V_DC
      CALL CALC_LA(2) ! Solve Lambda^{c} (-V_DC)
      CALL CALC_LC()  ! Solve {c}
      IF(GL%LSOLVER==102.OR.GL%LSOLVER==103)THEN
        CALL REGULARIZE_NKS(1)
        CALL CALC_ISIMIX(1)
      ENDIF
      IF(GL%LSOLVER==103)THEN
        CALL CALC_LA(1)
      ENDIF
      IF(GL%LSOLVER==1.OR.GL%LSOLVER==101.OR.GL%LSOLVER==3.OR.GL%LSOLVER==30.OR.GL%LSOLVER==102.OR.GL%LSOLVER==103)THEN
        CALL UPDATE_R_ALL()
      ENDIF
      CALL MAP_WH_BND_R(WH%NC_VAR,BND%NC_VAR,.FALSE.)
      CALL SET_WH_X123(.TRUE.) ! Update X
      CALL SET_WH_X1(WH%XN0,WH%NKS,WH%NA2MAX,.TRUE.,1)
      CALL SET_WH_X1(WH%XNC,WH%NC_VAR,WH%NA2MAX,.TRUE.,1)
!
      IF(GL%LSOLVER==2)THEN
        CALL SET_XTOWHX(FVEC,WH%XNC-WH%XN0(1:WH%DIMX1),WH%LXREAL(WH%DIMX1+1:WH%DIMX2),N,WH%DIMXH,.TRUE.)
      ELSE
        CALL SET_XTOWHX(FVEC(1:NX1),WH%X(1:WH%DIMX1)-WHX(1:WH%DIMX1),WH%LXREAL(1:WH%DIMX1),NX1,WH%DIMXG,.TRUE.)
        IF(GL%LSOLVER==103)THEN
          CALL SET_XTOWHX(FVEC(1+NX1:N),WH%X(WH%DIMX1+1:WH%DIMX2)-WHX(WH%DIMX1+1:WH%DIMX2),WH%LXREAL(WH%DIMX1+1:WH%DIMX2),NX2,WH%DIMXH,.TRUE.)
        ELSE
          CALL SET_XTOWHX(FVEC(1+NX1:N),WH%XNC-WH%XN0,WH%LXREAL(WH%DIMX1+1:WH%DIMX2),NX2,WH%DIMXH,.TRUE.)
        ENDIF
      ENDIF
      CALL CPU_TIME(TA2); CALL SYSTEM_CLOCK(TIB2,TIRATE); TB2=REAL(TIB2,4)/REAL(TIRATE,4)
      CALL OUT_TIME_USE('GUTZ_FCN1',TA2-TA1,TB2-TB1,GL%IO)
! Set back the input values
      IF(GL%LSOLVER==2)THEN
        CALL SET_XTOWHX(X,WH%X(WH%DIMX1+1:WH%DIMX2),WH%LXREAL(WH%DIMX1+1:WH%DIMX2),N,WH%DIMXH,.FALSE.)
      ELSEIF(GL%LSOLVER==3.OR.GL%LSOLVER==30)THEN
        CALL SET_XTOWHX(X(1:NX1),WH%X(1:WH%DIMX1),WH%LXREAL(1:WH%DIMX1),NX1,WH%DIMXG,.FALSE.)
        CALL SET_XTOWHX(X(1+NX1:N),WH%X(1+WH%DIMX2:WH%DIMX3),WH%LXREAL(1+WH%DIMX2:WH%DIMX3),NX2,WH%DIMXH,.FALSE.)
      ELSE
        CALL SET_XTOWHX(X,WH%X(1:WH%DIMX2),WH%LXREAL(1:WH%DIMX2),N,WH%DIMX2,.FALSE.)
      ENDIF
      CALL SET_WH_X123(.FALSE.) ! UPDATE R,LA,etc.
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("************ DIF_X ************")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(5F14.8)')FVEC
      MAXERR=MAXVAL(ABS(FVEC))
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(*,'(" ITER1=",I3," MAXERR=",F14.8)')GL%ITER1,MAXERR
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ITER1=",I3," MAXERR=",F14.8)')GL%ITER1,MAXERR
      IF(MAXERR.LT.RCUT)IFLAG=-1
      IF(GL%ITER1.GE.GL%NMAX_ITER.OR.GL%ITER1==-1)IFLAG=-1
      IF(GL%LSCF==3)IFLAG=-1
      ! 'GUTZ_FCN1'
      RETURN
!
      END SUBROUTINE GUTZ_FCN1
!
!*************************************************************
      SUBROUTINE GFERMI_FCN(N,X,FVEC,IFLAG)
      USE gprec; USE GUTZ
      IMPLICIT NONE
      INTEGER N
      INTEGER IFLAG
      REAL(gq) X(N),FVEC(N)
! LOCAL
      REAL(gq) DIFF,NELE
!
      GL%ITER_GF=GL%ITER_GF+1
      IF(GL%LMODEL==0)THEN
        CALL SUM_GFK_FULL(X(1))
        CALL GREEN_OCC(0,GF,NELE=NELE)
      ELSE
        CALL CALC_GF_IMP(X(1))
        CALL GREEN_OCC_R1(1,1,GF,BND%R,NELET=NELE)
      ENDIF
      NELE=NELE*BND%RSPO
      FVEC(1)=NELE-BND%NELEL
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(*    ,'(" ITER_GF=",I3," MU=",F12.6," NELE=",F14.8," BND%NELEL=",F14.8)')GL%ITER_GF,X(1),NELE,BND%NELEL
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ITER_GF=",I3," MU=",F12.6," NELE=",F14.8," BND%NELEL=",F14.8)')GL%ITER_GF,X(1),NELE,BND%NELEL
      RETURN
!
      END SUBROUTINE GFERMI_FCN
!
!====================================
! V2=A*V1
!====================================
      SUBROUTINE PHIK_AV(N,V1,V2)
      USE GUTZ; USE gprec
      IMPLICIT NONE
      INTEGER N
      COMPLEX(gq) V1(N),V2(N)
! LOCAL
      INTEGER NI
      !
!
      !
      NI=GL%NI
      GYZ(NI)%PJ%ITER=GYZ(NI)%PJ%ITER+1
      IF(.NOT.GYZ(NI)%PJ%LMCFLY)THEN
        CALL ZCSR_SYAMUX_SK('L',GYZ(NI)%PJ%FHL,V1,V2)
      ELSE
        CALL ZCSR_SYAMUX_SK('L',GYZ(NI)%PJ%U_KKH,V1,V2)
        CALL ACT_ZBMC(N,V1,V2,GYZ(NI)%CO,GYZ(NI)%HL,GYZ(NI)%PJ)
      ENDIF
      !
      RETURN
!
      END SUBROUTINE PHIK_AV
!
!*************************************************************
      SUBROUTINE GUTZ_FCN_PJS(V,N,FUN)
      USE gprec; USE GUTZ
      IMPLICIT NONE
      INTEGER N
      COMPLEX(gq) V(N)
      REAL(gq) FUN
! LOCAL
      INTEGER NI
      REAL(gq) TR
      COMPLEX(gq) VP(N)
      COMPLEX(gq),EXTERNAL::ZDOTC
!
      NI=GL%NI
      TR=ZDOTC(N,V,1,V,1)
      GYZ(NI)%PJ%C=V/SQRT(TR)
      CALL CALC_PJ_RHO(GYZ(NI)%HL,GYZ(NI)%PJ)
      IF((GYZ(NI)%HL%EVEC%NBK<=0).AND.(GL%LGPRJ/=11).AND.(GL%LGPRJ/=14))THEN
        CALL ZBM_LOAD(GYZ(NI)%HL%EVEC,NI,0,2)
      ENDIF
      CALL CALC_NCPHY_COSYM(GYZ(NI)%FS,GYZ(NI)%HL,GYZ(NI)%CO,GYZ(NI)%PJ)
      CALL CALC_P0_FS(GYZ(NI)%CO,GYZ(NI)%FS,GYZ(NI)%HL,GYZ(NI)%PJ)
      CALL CALC_PROJ_ENS(GYZ(NI)%HL,GYZ(NI)%PJ)
      CALL PHIK_AV(N,GYZ(NI)%PJ%C,VP)
      FUN=ZDOTC(N,GYZ(NI)%PJ%C,1,VP,1)-GYZ(NI)%PJ%PJ_ENS*KPT%DELTA
      RETURN
!
      END SUBROUTINE GUTZ_FCN_PJS
!
!*************************************************************
      SUBROUTINE GUTZ_FCN_NK(N,X,FVEC,IFLAG)
      USE gprec; USE GUTZ
      IMPLICIT NONE
      INTEGER N,IFLAG
      REAL(gq) X(N),FVEC(N)
! LOCAL
      INTEGER I
      REAL(gq),PARAMETER::RCUT=1.E-8_gq
      REAL(gq) MAXERR
!
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("************ WH%LA1 ************")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(5F14.8)')X
      GL%ITER_LA1=GL%ITER_LA1+1
      CALL SET_XTOWHX(X,WH%X(WH%DIMX1+1:WH%DIMX2),WH%LXREAL(WH%DIMX1+1:WH%DIMX2),N,WH%DIMXH,.FALSE.)
      CALL SET_WH_X1(WH%X(1+WH%DIMX1:WH%DIMX2),WH%LA1,GL%NA2MAX,.FALSE.,1)
      CALL MAP_WH_BND_R(WH%LA1,BND%LA1,.FALSE.)
      CALL ZOUT_MAT('LA1',WH%LA1,GL%IO)
      IF(GL%LMODEL==0)THEN
        CALL CALC_BAND_ALL(GL%LGREEN==0)
        IF(GL%LGREEN==0)THEN
          CALL GUTZ_FERMI(GL%IO)
        ELSE
          CALL GREEN_FERMI_NV()
          IF(GL%LGREEN==2)THEN
            CALL CALC_SELF_ENERGY_FULL()
          ENDIF
          CALL SUM_GFK_FULL(BND%EF)
        ENDIF
      ELSE
        CALL CALC_GF_IMP_ALL()
      ENDIF
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GA FERMI LEVEL=",F16.8)')BND%EF
      CALL CALC_NKS(GL%LGREEN)
      CALL SET_WH_X1(WH%XN0,WH%NKS,WH%NA2MAX,.TRUE.,1)
      CALL SET_XTOWHX(FVEC,WH%XN0-WH%X(1+WH%DIMX2:WH%DIMX3),WH%LXREAL(1+WH%DIMX2:WH%DIMX3),N,WH%DIMXH,.TRUE.)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("************ DIF_NK ************")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(5F14.8)')FVEC
      MAXERR=MAXVAL(ABS(FVEC))
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(*    ,'(" ITER_LA1=",I3," MAXERR=",F14.8)')GL%ITER_LA1,MAXERR
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ITER_LA1=",I3," MAXERR=",F14.8)')GL%ITER_LA1,MAXERR
      IF(MAXERR.LT.RCUT)IFLAG=-1
      IF(GL%ITER_LA1.GE.GL%NMAX_ITER)IFLAG=-1
      RETURN
!
      END SUBROUTINE GUTZ_FCN_NK
!
!*************************************************************
! ETA, Non-diaginal elements
!*************************************************************
      SUBROUTINE GUTZ_FCN_NKND(N,X,FVEC,IFLAG)
      USE gprec; USE GUTZ
      IMPLICIT NONE
      INTEGER N,IFLAG
      REAL(gq) X(N),FVEC(N)
! LOCAL
      INTEGER I
      COMPLEX(gq) :: XN1(WH%DIMXHO)
      REAL(gq),PARAMETER::RCUT=1.E-8_gq
      REAL(gq) MAXERR
!
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("************ WH%ETA ************")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(5F14.8)')X
      GL%ITER_ETA=GL%ITER_ETA+1
      CALL SET_XTOWHX(X,WH%X(WH%DIMX3+1:WH%DIMXT),WH%LXREAL(WH%DIMX3+1:WH%DIMXT),N,WH%DIMXHO,.FALSE.)
      CALL SET_WH_X1(WH%X(1+WH%DIMX3:WH%DIMXT),WH%ETA,GL%NA2MAX,.FALSE.,2)
      CALL MAP_WH_BND_R(WH%ETA,BND%ETA,.FALSE.)
      CALL ZOUT_MAT('ETA',WH%ETA,GL%IO)
      IF(GL%LMODEL==0)THEN
        CALL CALC_BAND_ALL(GL%LGREEN==0)
        IF(GL%LGREEN==0)THEN
          CALL GUTZ_FERMI(GL%IO)
        ELSE
          CALL GREEN_FERMI_NV()
          IF(GL%LGREEN==2)THEN
            CALL CALC_SELF_ENERGY_FULL()
          ENDIF
          CALL SUM_GFK_FULL(BND%EF)
        ENDIF
      ELSE
        CALL CALC_GF_IMP_ALL()
      ENDIF
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" GA FERMI LEVEL=",F16.8)')BND%EF
      CALL CALC_NKS(GL%LGREEN)
      CALL SET_WH_X1(XN1,WH%NKS,WH%NA2MAX,.TRUE.,2)
      CALL SET_XTOWHX(FVEC,XN1,WH%LXREAL(1+WH%DIMX3:WH%DIMXT),N,WH%DIMXHO,.TRUE.)
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'("************ DIF_NK ************")')
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(5F14.8)')FVEC
      MAXERR=MAXVAL(ABS(FVEC))
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(*    ,'(" ITER_ETA=",I3," MAXERR=",F14.8)')GL%ITER_ETA,MAXERR
      IF(GP%MYRANK.EQ.GP%MASTER)WRITE(GL%IO,'(" ITER_ETA=",I3," MAXERR=",F14.8)')GL%ITER_ETA,MAXERR
      IF(MAXERR.LT.RCUT)IFLAG=-1
      IF(GL%ITER_ETA.GE.GL%NMAX_ITER)IFLAG=-1
      RETURN
!
      END SUBROUTINE GUTZ_FCN_NKND

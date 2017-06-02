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
      MODULE BANDSTRU
      USE gprec; USE GMPI; USE GUTIL
      USE GCONSTANT, ONLY: D0, D1, Z0, Z1, ZI
      IMPLICIT NONE
! DEFINE TYPE
      TYPE BAND_STRU
        INTEGER ISPIN,ISPIN_IN,NSPIN,ISO,ISO_IN,ISPO,RSPO,LFERWE    ! ISO=2 => SOC; LFERWE=1:Pre-store FERWE
        INTEGER NASOTOT    ! Maximal number of bands for the local projector
        INTEGER N_FROZEN
        INTEGER,POINTER :: NE(:,:)          ! 1: TOTAL NUMBER OF BANDS, 2-3: CORRELATED BANDS INTERVAL
        REAL(gq),POINTER :: EK0(:,:),EK(:,:,:) ! Kohn-Sham / Gutz Eigen values 
        COMPLEX(gq),POINTER :: R(:,:,:),LA1(:,:,:),NC_PHY(:,:,:),NKS(:,:,:),D0(:,:,:)
        COMPLEX(gq),POINTER :: ETA(:,:,:)
        COMPLEX(gq),POINTER :: NRL(:,:,:) ! For onsite local 1PDM
        COMPLEX(gq),POINTER :: UK (:,:,:,:) ! (MAX_IN_BANDS,NASOTOT,NSYMOP,NKP), <Psi_k|Phi_loc>
        COMPLEX(gq),POINTER :: SK (:,:,:,:) ! <Phi_loc|Phi_loc>
        COMPLEX(gq),POINTER :: VK (:,:,:,:,:) ! Gutz Eigen-vector
        COMPLEX(gq),POINTER :: UK0(:,:,:,:) ! BANDS unitary transformation
        COMPLEX(gq),POINTER :: HK0(:,:,:,:) ! LDA DISPERSION
        COMPLEX(gq),POINTER :: HK1(:,:,:,:) ! Additional bare Hamiltonian term which is not to be normalized. (CMR-LDA)
        INTEGER NMAX,NMAXIN,NVMAX                   ! NMAX: MIXIMAL NUMBER OF BANDS; NVMAX: Highest occupied band
        REAL(gq),POINTER :: FERWE(:,:,:),FERWER(:,:,:)       ! fermi-weight, CONVENTION INCLUDING KPT%WT AND SPIN DEGENERACY
        REAL(gq) EFLDA,NELET,NELET_FROZEN,NELEL,NELEC,MUP_DN
        REAL(gq) :: EF=0._gq
        REAL(gq) EBMIN,EBMAX,EBWIDTH
        COMPLEX(gq),POINTER::CWT(:,:,:)   ! Occupation matrix for \Rho
      END TYPE BAND_STRU
!
      TYPE K_POINTS
        INTEGER DIM,DIML,ISMEAR,ICOR
        REAL(gq),POINTER :: WT(:)=>NULL(),X(:)=>NULL(),Y(:)=>NULL(),Z(:)=>NULL()
        CHARACTER*10,POINTER :: NAME(:)
        REAL(gq) TWT,DELTA,CORDEG
        CHARACTER::FILE_NAME*128
      END TYPE K_POINTS
!
      TYPE ENERGY
        REAL(gq) HYBRD,BAND,DBLC,GAMM,TOT,TB,DIFF,TS2,DL,POT2,DC1,DC2 ! TS2: BZ INTE CORRECTION. DL: DEEP LEVEL
      END TYPE ENERGY
!
      TYPE SYM_INFO
        INTEGER IE,NOP,IDI,IDF,MODE 
! VASP_PAWSYM; MODE=10, REAL HARMONICS BASIS
        INTEGER LMAX,MMAX,NROTK,NIOND,NPCELL
        INTEGER,ALLOCATABLE :: ROTMAP(:,:,:)
        REAL(gq),ALLOCATABLE :: SL(:,:,:,:) ! (L,M) <- (L,MP)
      END TYPE SYM_INFO
!
! VARIABLES
      TYPE (BAND_STRU),SAVE :: BND
      TYPE (ENERGY)   ,SAVE :: ENG
      TYPE (SYM_INFO) ,SAVE :: SYM
      TYPE (K_POINTS) ,SAVE :: KPT
!
! SUBROUTINE
      CONTAINS
!****************************************************************************
      SUBROUTINE SET_BND_NISO(ISO_IN,ISPIN_IN,NBMAX,NKPT,IO)
      INTEGER ISO_IN,ISPIN_IN,NBMAX,NKPT,IO
!
      BND%ISO_IN=ISO_IN
      IF(BND%ISO<ISO_IN)THEN
        IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" WARNING: BND%ISO<ISO_IN! SET BND%ISO=ISO_IN!")')
        BND%ISO=ISO_IN
      ENDIF
      BND%ISPIN_IN=ISPIN_IN; KPT%DIM=NKPT; BND%ISPO=MAX(BND%ISO,BND%ISPIN)
      BND%NSPIN=MAX(1,BND%ISPIN/BND%ISO)
      BND%RSPO=3-BND%ISPO
      BND%NMAX=NBMAX*BND%ISO/ISO_IN
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" ISO_IN=",I2," ISPIN_IN=",I2," ISO=",I2," ISPIN=",I2," ISPO=",I2)')BND%ISO_IN,BND%ISPIN_IN,BND%ISO,BND%ISPIN,BND%ISPO
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" BND%NMAX=",I4)')BND%NMAX
      RETURN
!
      END SUBROUTINE SET_BND_NISO
!
!********************************************************************************
      SUBROUTINE SET_BND_UVEK(NKP,NKPL,WT,MODE,KX,KY,KZ,KNAME)
      INTEGER NKP,NKPL,MODE
      REAL(gq) WT(NKP)
      REAL(gq),OPTIONAL::KX(NKP),KY(NKP),KZ(NKP)
      CHARACTER*10,OPTIONAL::KNAME(NKP)
!
      KPT%DIML=NKPL
      IF(NKPL>0)THEN
        ALLOCATE(KPT%WT(NKP)); KPT%TWT=SUM(WT); KPT%WT=WT/KPT%TWT
        IF(BND%LFERWE==0)THEN
          ALLOCATE(BND%FERWE (BND%NMAX,NKP,BND%NSPIN)); BND%FERWE=0
          ALLOCATE(BND%FERWER(BND%NMAX,NKP,BND%NSPIN)); BND%FERWER=0
        ENDIF
        ALLOCATE(BND%EK(BND%NMAX,NKP,BND%NSPIN),BND%EK0(BND%NMAX,NKP)); BND%EK=0; BND%EK0=0
        ALLOCATE(BND%UK (BND%NMAXIN,BND%NASOTOT,SYM%IDI:SYM%IDF,NKPL)); BND%UK=0
        ALLOCATE(BND%SK (BND%NASOTOT,BND%NASOTOT,SYM%IDI:SYM%IDF,NKPL)); BND%SK=0
        ALLOCATE(BND%UK0(BND%NMAXIN,BND%NMAXIN,SYM%IDI:SYM%IDF,NKPL)); BND%UK0=0
        ALLOCATE(BND%HK0(BND%NMAXIN,BND%NMAXIN,SYM%IDI:SYM%IDF,NKPL)); BND%HK0=0
        IF(MODE==1)THEN
          ALLOCATE(BND%HK1(BND%NMAXIN,BND%NMAXIN,SYM%IDI:SYM%IDF,NKPL)); BND%HK1=0
        ENDIF
        ALLOCATE(BND%VK (BND%NMAXIN,BND%NMAXIN,SYM%IDI:SYM%IDF,NKPL, &
            &BND%NSPIN)); BND%VK=0
      ENDIF
      ALLOCATE(BND%R  (BND%NASOTOT,BND%NASOTOT,BND%NSPIN)); BND%R  =0
      ALLOCATE(BND%D0 (BND%NASOTOT,BND%NASOTOT,BND%NSPIN)); BND%D0 =0
      ALLOCATE(BND%LA1(BND%NASOTOT,BND%NASOTOT,BND%NSPIN)); BND%LA1=0
      ALLOCATE(BND%ETA(BND%NASOTOT,BND%NASOTOT,BND%NSPIN)); BND%ETA=0
      ALLOCATE(BND%NC_PHY(BND%NASOTOT,BND%NASOTOT,BND%NSPIN)); BND%NC_PHY=0
      ALLOCATE(BND%NKS(BND%NASOTOT,BND%NASOTOT,BND%NSPIN)); BND%NKS=0
      ALLOCATE(BND%NRL(BND%NASOTOT,BND%NASOTOT,BND%NSPIN)); BND%NRL=0
!
      IF(PRESENT(KX))THEN
        ALLOCATE(KPT%X(NKP),KPT%Y(NKP),KPT%Z(NKP),KPT%NAME(NKP))
        KPT%X=KX; KPT%Y=KY; KPT%Z=KZ; KPT%NAME=KNAME
      ENDIF
!
      RETURN
!
      END SUBROUTINE SET_BND_UVEK
!
!********************************************************************************
      SUBROUTINE READ_BNDU(IU,MODE)
      INTEGER IU,MODE
! LOCAL
      INTEGER IVEC,IFILE,NKP,IUP
      CHARACTER FNAME*32
      LOGICAL LEXIST
!
      ! 'READ_BNDU'
      DO IVEC=1,GP%NVEC
        IF(GP%LKPVEC)THEN
          IFILE=GP%KVEC(IVEC,1); NKP=GP%KVEC(IVEC,2)
        ELSE
          IFILE=GP%MYRANK; NKP=KPT%DIM
        ENDIF
        FNAME = FILE_NAME(IFILE,0,0,-1)
        IUP=IU+IVEC
        INQUIRE(FILE=FNAME, EXIST=LEXIST)
        IF (LEXIST) THEN
          CALL READ_BNDU_UNFORMATTED(FNAME, IVEC, NKP, IUP, MODE)
          CYCLE
        ENDIF
        IF(MODE > 0)THEN
          STOP " ERROR IN READ_BNDU: MODE > 0 NOT AVAILABLE!"
        ENDIF
        FNAME = FILE_NAME(IFILE,0,0,-11)
        INQUIRE(FILE=FNAME, EXIST=LEXIST)
        IF (LEXIST) THEN 
          CALL READ_BNDU_FORMATTED(FNAME, IVEC, NKP, IUP)
          CYCLE 
        ENDIF
        STOP " ERROR: NO BNDU FILES!"
      ENDDO
      ! 'READ_BNDU'
      RETURN
!
      END SUBROUTINE READ_BNDU
!
!********************************************************************************
      SUBROUTINE READ_BNDU_UNFORMATTED(FNAME, IVEC, NKP, IU, MODE)
      INTEGER,INTENT(IN)::IVEC, NKP, IU, MODE
      CHARACTER,INTENT(IN)::FNAME*32
! LOCAL
      INTEGER I,IKP,IKPL,IKS,ISYM,NBANDS,NBANDS2,NBTOT,NBTOT2,NASOTOT,NASOTOT2
      COMPLEX(gq),ALLOCATABLE::UK(:,:)
!
      NASOTOT2=BND%NASOTOT; NASOTOT=NASOTOT2*BND%ISO_IN/BND%ISO
      ALLOCATE(UK(BND%NMAXIN,NASOTOT))
      OPEN(IU, FILE = FNAME, STATUS = 'OLD', FORM = "UNFORMATTED", ACCESS = "SEQUENTIAL")
      IKPL=0
      DO I=1,IVEC-1; IKPL=IKPL+GP%KVEC(I,2); ENDDO
      DO IKS=1,NKP
        IF(GP%LKPVEC)THEN
          IKP=GP%KVEC(IVEC,3)+IKS
          IF(GP%LOMP)THEN; IKPL=IKP; ELSE; IKPL=IKPL+1; ENDIF
        ELSE
          IKP=IKS
          IF(GP%LOMP)THEN; IKPL=IKP; ELSE; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF
        ENDIF
        IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
        NBTOT2=BND%NE(1,IKP)
        NBTOT =NBTOT2*BND%ISO_IN/BND%ISO
        READ(IU)BND%EK0(1:NBTOT,IKP)
        NBANDS2=BND%NE(3,IKP)-BND%NE(2,IKP)+1
        NBANDS =NBANDS2*BND%ISO_IN/BND%ISO
        DO ISYM=1,SYM%NOP
          IF(ISYM<SYM%IDI.OR.ISYM>SYM%IDF)THEN
            DO I=1,NASOTOT
              READ(IU)UK(1:NBANDS,I)
            ENDDO
          ELSE
            IF(MODE>0)THEN
              READ(IU)BND%HK0(1:NBTOT,1:NBTOT,ISYM,IKPL)
              READ(IU)BND%HK1(1:NBTOT,1:NBTOT,ISYM,IKPL)  ! In diagonal basis of the whole H
            ENDIF
            DO I=1,NASOTOT
              READ(IU)BND%UK(1:NBANDS,I,ISYM,IKPL)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(UK)
      RETURN
!
      END SUBROUTINE READ_BNDU_UNFORMATTED
!
!********************************************************************************
      SUBROUTINE READ_BNDU_FORMATTED(FNAME, IVEC, NKP, IU)
      INTEGER,INTENT(IN)::IVEC, NKP, IU
      CHARACTER,INTENT(IN)::FNAME*32
! LOCAL
      INTEGER I,IKP,IKPL,IKS,ISO,ISYM,NBANDS,NBANDS2,NBTOT,NBTOT2,NASOTOT,NASOTOT2
      REAL(gq),ALLOCATABLE::RBUF(:),IBUF(:)
!
      NASOTOT2=BND%NASOTOT; NASOTOT=NASOTOT2*BND%ISO_IN/BND%ISO
      ALLOCATE(RBUF(BND%NMAX), IBUF(BND%NMAX))
      OPEN(IU, FILE = FNAME, STATUS = 'OLD')
      IKPL=0
      DO I=1,IVEC-1; IKPL=IKPL+GP%KVEC(I,2); ENDDO
      DO IKS=1,NKP
        IF(GP%LKPVEC)THEN
          IKP=GP%KVEC(IVEC,3)+IKS
          IF(GP%LOMP)THEN; IKPL=IKP; ELSE; IKPL=IKPL+1; ENDIF
        ELSE
          IKP=IKS
          IF(GP%LOMP)THEN; IKPL=IKP; ELSE; IKPL=IKP-GP%MYRANK*KPT%DIML; ENDIF
        ENDIF
        IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
        NBTOT2=BND%NE(1,IKP)
        NBTOT =NBTOT2*BND%ISO_IN/BND%ISO
        READ(IU,*)BND%EK0(1:NBTOT,IKP)
        NBANDS2=BND%NE(3,IKP)-BND%NE(2,IKP)+1
        NBANDS =NBANDS2*BND%ISO_IN/BND%ISO
        DO ISYM=1,SYM%NOP; DO ISO = 1, NASOTOT
          READ(IU,*)RBUF(1:NBANDS); READ(IU,*)IBUF(1:NBANDS)
          IF(ISYM<SYM%IDI.OR.ISYM>SYM%IDF)CYCLE
          BND%UK(1:NBANDS,ISO,ISYM,IKPL)=RBUF(1:NBANDS)+ZI*IBUF(1:NBANDS)
        ENDDO; ENDDO
      ENDDO
      DEALLOCATE(RBUF, IBUF)
      RETURN
!
      END SUBROUTINE READ_BNDU_FORMATTED
!
!********************************************************************************
      SUBROUTINE PROC_BNDU(MODE)
      INTEGER MODE
! LOCAL
      INTEGER IVEC,IKP,IKPL,NKP,IKS,ISYM,NBANDS,NBANDS2,ISP,NBTOT,NBTOT2,NASOTOT,NASOTOT2
!
      NASOTOT2=BND%NASOTOT; NASOTOT=NASOTOT2*BND%ISO_IN/BND%ISO
      IKPL=0
      DO IVEC=1,GP%NVEC
        IF(GP%LKPVEC)THEN
          NKP=GP%KVEC(IVEC,2)
        ELSE
          NKP=KPT%DIM
        ENDIF
        DO IKS=1,NKP
          IF(GP%LKPVEC)THEN
            IKP=GP%KVEC(IVEC,3)+IKS; IKPL=IKPL+1
          ELSE
            IKP=IKS; IKPL=IKP-GP%MYRANK*KPT%DIML
          ENDIF
          IF(IKPL.LE.0)CYCLE; IF(IKPL.GT.KPT%DIML)EXIT
          NBTOT2=BND%NE(1,IKP)
          NBTOT =NBTOT2*BND%ISO_IN/BND%ISO
          IF(BND%ISO>BND%ISO_IN)THEN
            BND%EK0(1+NBTOT:NBTOT2,IKP)=BND%EK0(1:NBTOT,IKP)
            CALL ORBITAL_SPIN_TRANS(BND%EK0(1:NBTOT2,IKP),NBTOT2,.TRUE.,BND%ISO_IN)
          ENDIF
          NBANDS2=BND%NE(3,IKP)-BND%NE(2,IKP)+1
          NBANDS =NBANDS2*BND%ISO_IN/BND%ISO
          DO ISYM=1,SYM%NOP
            IF(ISYM<SYM%IDI.OR.ISYM>SYM%IDF)CYCLE
            IF(MODE>0)THEN
              IF(BND%ISO>BND%ISO_IN)THEN
                BND%HK0(1+NBTOT:NBTOT2,1+NBTOT:NBTOT2,ISYM,IKPL)=BND%HK0(1:NBTOT,1:NBTOT,ISYM,IKPL)
                BND%HK1(1+NBTOT:NBTOT2,1+NBTOT:NBTOT2,ISYM,IKPL)=BND%HK1(1:NBTOT,1:NBTOT,ISYM,IKPL)
                CALL ORBITAL_SPIN_TRANS(BND%HK0(:,:,ISYM,IKPL),NBTOT2,NBTOT2,.TRUE.,BND%ISO_IN,LURIGHT=.FALSE.)
                CALL ORBITAL_SPIN_TRANS(BND%HK1(:,:,ISYM,IKPL),NBTOT2,NBTOT2,.TRUE.,BND%ISO_IN,LURIGHT=.FALSE.)
              ENDIF
            ENDIF
            IF(BND%ISO>BND%ISO_IN)THEN
              BND%UK(1+NBANDS:NBANDS2,1+NASOTOT:NASOTOT2,ISYM,IKPL)=BND%UK(1:NBANDS,1:NASOTOT,ISYM,IKPL)
              CALL ORBITAL_SPIN_TRANS(BND%UK(1:NBANDS2,1:NASOTOT2,ISYM,IKPL),NBANDS2,NASOTOT2,.TRUE.,BND%ISO_IN,LURIGHT=.FALSE.)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CALL DSUM_ALL_MPI(BND%EK0,BND%NMAX*KPT%DIM)
      DO ISP=1,BND%NSPIN; BND%EK(:,:,ISP)=BND%EK0(:,:); ENDDO
      RETURN
!
      END SUBROUTINE PROC_BNDU
!
!******************************************
! Locate the index of identity operator
!******************************************
      SUBROUTINE GUTZ3_LOCATE_SYM_IE(IORD,IZ,TAU)
      INTEGER :: IORD
      INTEGER,OPTIONAL :: IZ(3,3,IORD)
      REAL(gq),OPTIONAL :: TAU(3,IORD)
! LOCAL
      INTEGER ISYM,I,J
!
      IF(.NOT.PRESENT(IZ))THEN
        SYM%IE=1; SYM%NOP=1; SYM%IDI=1; SYM%IDF=1
        RETURN
      ENDIF
!
      DO ISYM=1,IORD
      DO I=1,3; DO J=1,3
      IF(I.EQ.J)THEN
        IF(IZ(I,J,ISYM).NE.1)GOTO 100
      ELSE
        IF(IZ(I,J,ISYM).NE.0)GOTO 100
      ENDIF
      ENDDO; ENDDO
      DO I=1,3; IF(ABS(TAU(I,ISYM))>1.E-16_gq)GOTO 100; ENDDO
      SYM%IE=ISYM; SYM%NOP=IORD
      IF(SYM%MODE==0)THEN
        SYM%IDI=SYM%IE; SYM%IDF=SYM%IE
      ELSE
        SYM%IDI=1; SYM%IDF=IORD
      ENDIF
      RETURN
100   CONTINUE
      ENDDO ! ISYM
      STOP ' ERROR: FAILED TO LOCATE IDENTITY OPERATION!'
!
      END SUBROUTINE GUTZ3_LOCATE_SYM_IE
!
!****************************************************************************
! Calc band energy
!****************************************************************************
      SUBROUTINE CALC_EBND()
      INTEGER IK
!
      ENG%BAND=ENG%TS2; ENG%DL=0
      DO IK=1,KPT%DIM
        ENG%DL  =ENG%DL  +SUM(BND%EK(1:BND%NE(2,IK)-1,IK,:)*BND%FERWE(1:BND%NE(2,IK)-1,IK,:))
        ENG%BAND=ENG%BAND+SUM(BND%EK(:,IK,:)*BND%FERWE(:,IK,:))
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_EBND
!
!****************************************************************************
! f_n <Psi_n|a> R_{a,A} R^+_{B,b} <b|Psi_n>
!=R^+_{B,b} <b|Psi_n> f_n <Psi_n|a> R_{a,A}
!****************************************************************************
      SUBROUTINE CALC_RNRL()
      INTEGER ISP,N
!
      N=BND%NASOTOT
      DO ISP=1,BND%NSPIN
        CALL ZGEMM('T','N',N,N,N,Z1,BND%NKS(:,:,ISP),N,BND%R(:,:,ISP),N,Z0,BND%NRL(:,:,ISP),N) ! <b|Psi_n> f_n <Psi_n|a> R_{a,A}
        CALL ANNXB('C',BND%R(:,:,ISP),BND%NRL(:,:,ISP),N,N) ! R^+_{B,b} <b|Psi_n> f_n <Psi_n|a> R_{a,A}
        BND%NRL(:,:,ISP)=TRANSPOSE(BND%NRL(:,:,ISP))
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_RNRL
!
!****************************************************************************
! Calc total magnetic moment
!****************************************************************************
      SUBROUTINE CALC_MUP_DN()
      INTEGER IK
!
      BND%MUP_DN=0
      IF(BND%NSPIN==1)RETURN
      DO IK=1,KPT%DIM
        BND%MUP_DN=BND%MUP_DN+SUM(BND%FERWE(:,IK,1)-BND%FERWE(:,IK,2))
      ENDDO
      RETURN
!
      END SUBROUTINE CALC_MUP_DN
!
!****************************************************************************
! Band width, gaps
!****************************************************************************
      SUBROUTINE CALC_BAND_WIDGAP(EK,NMAX,NKP,EF,IO)
      INTEGER ::NMAX,NKP,IO
      REAL(gq)::EK(NMAX,NKP),EF
! LOCAL
      INTEGER IB,IBLK
      INTEGER ::NBLK(2,NMAX)
      REAL(gq)::EBOT1,EBOT2,ETOP1,ETOP2,EBLK(2,NMAX)
!
      IBLK=1; EBOT1=MINVAL(EK(1,:)); ETOP1=MAXVAL(EK(1,:))
      EBLK(1,1)=EBOT1;  EBLK(2,1)=ETOP1
      NBLK(1,1)=1    ;  NBLK(2,1)=1
      DO IB=2,NMAX
        EBOT2=MINVAL(EK(IB,:)); ETOP2=MAXVAL(EK(IB,:))
        IF(EBOT2.GE.ETOP1)THEN
          EBLK(2,IBLK)=ETOP1; NBLK(2,IBLK)=IB-1
          IBLK=IBLK+1
          EBLK(1,IBLK)=EBOT2; EBLK(2,IBLK)=ETOP2
          NBLK(1,IBLK)=IB   ; NBLK(2,IBLK)=IB
        ENDIF
        EBOT1=EBOT2; ETOP1=ETOP2
      ENDDO
      IF(NBLK(2,IBLK).NE.NMAX)THEN
        EBLK(2,IBLK)=ETOP1; NBLK(2,IBLK)=NMAX
      ENDIF
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" There are ",I2," blockes of bands, with band energy windows(wrt E_Fermi):")')IBLK
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,100)EBLK(:,1:IBLK)-EF
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" Band window indices:")')
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,102)NBLK(:,1:IBLK)
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" Band width")')
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,101)EBLK(2,1:IBLK)-EBLK(1,1:IBLK)
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" Band gap")')
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,101)(EBLK(1,IB+1)-EBLK(2,IB),IB=1,IBLK-1)
      RETURN
100   FORMAT(5(2F8.3," | "))
101   FORMAT(5(4X,F8.3,"   |  "))
102   FORMAT(5(3X,2I7," | "))
!
      END SUBROUTINE CALC_BAND_WIDGAP
!
!****************************************************************************
      SUBROUTINE CALC_CORR_EBWIDTH(IO)
      INTEGER IO
! LOCAL
      INTEGER IK
      REAL(gq) EMIN,EMAX
!
      EMIN=100._gq; EMAX=-100._gq
      DO IK=1,KPT%DIM
        EMIN=MIN(EMIN,BND%EK0(BND%NE(2,IK),IK))
        EMAX=MAX(EMAX,BND%EK0(BND%NE(3,IK),IK))
      ENDDO
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" CORRELATED BLOCK: EMIN/EMIX=",2F10.4)')EMIN,EMAX
      BND%EBMIN=EMIN; BND%EBMAX=EMAX
      BND%EBWIDTH=EMAX-EMIN
      RETURN
!
      END SUBROUTINE CALC_CORR_EBWIDTH
!
!****************************************************************************
      SUBROUTINE SETUP_TB(MODE)
      INTEGER MODE
! LOCAL
      INTEGER ISYM,IVEC,IKP,IKPL,IKS,IB
      INTEGER NEMIN,NEMAX,NKP,NBMAX,NBANDS
      COMPLEX(gq),POINTER :: UK(:,:),VK(:,:),HK(:,:)
!
      !"SETUP_TB"
      NBMAX=BND%NMAXIN; IKPL=0
!
      !$OMP PARALLEL DO FIRSTPRIVATE(IKPL) PRIVATE(IVEC,ISYM,NKP,IB,IKS,IKP,NEMIN,NEMAX,NBANDS,UK,VK,HK) SCHEDULE(STATIC,1)
      DO IVEC=1,GP%NVEC;   IF(GP%LKPVEC)THEN;     NKP=GP%KVEC(IVEC,2);   ELSE;     NKP=KPT%DIM;   ENDIF;   DO IKS=1,NKP;   IF(GP%LKPVEC)THEN;     IKP=GP%KVEC(IVEC,3)+IKS;     IF(GP%LOMP)THEN;       IKPL=IKP;     ELSE;       IKPL=IKPL+1;    ENDIF;   ELSE;     IKP=IKS;     IF(GP%LOMP)THEN;       IKPL=IKP;     ELSE;       IKPL=IKP-GP%MYRANK*KPT%DIML;     ENDIF;   ENDIF;   IF(IKPL.LE.0) CYCLE;   IF(IKPL.GT.KPT%DIML) EXIT
      NEMIN=BND%NE(2,IKP); NEMAX=BND%NE(3,IKP)
      NBANDS=NEMAX-NEMIN+1
      DO ISYM=SYM%IDI,SYM%IDF
        UK=>BND%UK (1:NBANDS,:,ISYM,IKPL)  ! <psi_0 | loc>
        VK=>BND%UK0(1:NBANDS,1:NBANDS,ISYM,IKPL) ! <psi_0|loc + non loc>
        HK=>BND%HK0(1:NBANDS,1:NBANDS,ISYM,IKPL)
        CALL SETUP_TB_1K(UK,VK,NBANDS,BND%NASOTOT)
        IF(MODE==0)THEN
          DO IB=1,NBANDS; HK(IB,:)=BND%EK0(NEMIN+IB-1,IKP)*VK(IB,:); ENDDO
          CALL ANNXB('C',VK,HK,NBANDS,NBANDS)
        ELSE
          CALL UHAU(HK,VK,NBANDS,NBANDS)
          CALL UHAU(BND%HK1(1:NBANDS,1:NBANDS,ISYM,IKPL),VK,NBANDS,NBANDS)
        ENDIF
        BND%VK(1:NBANDS,1:NBANDS,ISYM,IKPL,1)=TRANSPOSE(CONJG(VK)) ! <loc. orb. |psi0>
      ENDDO ! ISYM
      ENDDO; ENDDO
      !$OMP END PARALLEL DO
      NULLIFY(UK,VK,HK)
      !"SETUP_TB"
      RETURN
!
      END SUBROUTINE SETUP_TB
!
!****************************************************************************
      SUBROUTINE SETUP_TB_1K(UK,VK,NBANDS,NASOT)
      USE GUTIL; USE GCONSTANT
      INTEGER NBANDS,NASOT
      COMPLEX(gq) UK(NBANDS,NASOT),VK(NBANDS,NBANDS)
! LOCAL
      LOGICAL LOK
      REAL(gq),ALLOCATABLE :: W(:)
!
      LOK=.TRUE.
      CALL ZGEMM('N','C',NBANDS,NBANDS,NASOT,-Z1,UK,NBANDS,UK,NBANDS,Z0,VK,NBANDS) ! <psi0|a><a|psi0>
      ALLOCATE(W(NBANDS)); W=0
      CALL HERMEV('V','U',VK,W,NBANDS)
      IF(ABS(W(NASOT)+1)>1.E-6_gq) LOK=.FALSE.
      IF(NASOT<NBANDS)THEN
        IF(ABS(W(NASOT+1))>1.E-6_gq) LOK=.FALSE.
      ENDIF
      IF(.NOT.LOK)THEN
        WRITE(0,'(" FETAL ERROR: NBANDS=",I5," NASOT=",I5," WITH W OF LOCAL PROJECTOR:")')NBANDS,NASOT
        WRITE(0,*)W; STOP
      ENDIF
      VK(:,1:NASOT)=UK  ! The lowest are local orbitals; <psi0|loc. orb.>
      DEALLOCATE(W)
      RETURN
!
      END SUBROUTINE SETUP_TB_1K
!
!****************************************************************************
      SUBROUTINE READ_PAWSYM(IU)
      INTEGER IU
!
      OPEN(IU,FILE='PAWSYM.INP',STATUS='OLD')
      READ(IU,*)SYM%LMAX,SYM%NROTK
      SYM%MMAX=2*SYM%LMAX+1
      ALLOCATE(SYM%SL(SYM%MMAX,SYM%MMAX,0:SYM%LMAX,SYM%NROTK))
      READ(IU,*)SYM%SL(1:,1:,0:,1:)
      READ(IU,*)SYM%NIOND,SYM%NPCELL
      ALLOCATE(SYM%ROTMAP(SYM%NIOND,SYM%NROTK,SYM%NPCELL))
      READ(IU,*)SYM%ROTMAP
      CLOSE(IU)
      RETURN
!
      END SUBROUTINE READ_PAWSYM
!
!*******************************************************************************
      SUBROUTINE ANNSYM1(ANN,N,NIONS,L)
      INTEGER N,NIONS,L
      COMPLEX(gq) ANN(N,N,NIONS)
! LOCAL
      INTEGER IROT,M,MP,NI,NIP,ITRANS
      COMPLEX(gq) BUF1(N,N,NIONS),BUF2(N,N,NIONS),TMP(N,N,NIONS)
!
      TMP=0
      DO IROT=1,SYM%NROTK
      BUF1=0
      DO NI=1,NIONS; DO M=1,2*L+1; DO MP=1,2*L+1
        BUF1(M,:,NI)=BUF1(M,:,NI)+SYM%SL(M,MP,L,IROT)*ANN(MP,:,NI) ! Left trans
      ENDDO; ENDDO; ENDDO
      BUF2=0
      DO NI=1,NIONS; DO M=1,2*L+1; DO MP=1,2*L+1
        BUF2(:,M,NI)=BUF2(:,M,NI)+BUF1(:,MP,NI)*SYM%SL(M,MP,L,IROT) ! Right trans
      ENDDO; ENDDO; ENDDO
      DO NI=1,NIONS; DO ITRANS=1,SYM%NPCELL
        NIP=SYM%ROTMAP(NI,IROT,ITRANS)
        TMP(:,:,NI)=TMP(:,:,NI)+BUF2(:,:,NIP)
      ENDDO; ENDDO
      ENDDO ! IROT
      ANN=TMP/SYM%NROTK/SYM%NPCELL
      RETURN
!
      END SUBROUTINE ANNSYM1
!
!*******************************************************************************
      SUBROUTINE SET_KPT_ICOR(IO)
      INTEGER IO
!
      IF(KPT%ISMEAR/=-5)RETURN
      KPT%ICOR=1
      IF(ABS(KPT%DELTA)>100._gq)KPT%ICOR=0
      KPT%CORDEG=-1.D-6
      IF(KPT%DELTA>100._gq)THEN
        KPT%CORDEG=-KPT%DELTA+100._gq
      ELSEIF(KPT%DELTA>0._gq)THEN
        KPT%CORDEG=-KPT%DELTA
      ELSEIF(KPT%DELTA<-100._gq)THEN
        KPT%CORDEG= KPT%DELTA-100._gq
      ELSEIF(KPT%DELTA<0._gq)THEN
        KPT%CORDEG= KPT%DELTA
      ENDIF
      IF(ABS(KPT%CORDEG)>.01_gq)KPT%CORDEG=-1.D-6
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" BZ INTEGRATION WITH TETRA METHOD, ICOR=",I2)')KPT%ICOR
      IF(KPT%CORDEG<0._gq)THEN
        IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" EQUAL OCCUPANCY OF DEGENERATE STATES, TOL=",E8.1)')-KPT%CORDEG
      ENDIF
      RETURN
!
      END SUBROUTINE SET_KPT_ICOR
!
!*******************************************************************************
      SUBROUTINE GUTZ_FERMI(IO)
      INTEGER IO
!
      IF(KPT%ISMEAR==-5)THEN
        CALL GUTZ_FERMI_TETRA_W2K()
      ELSEIF(KPT%ISMEAR==0.OR.KPT%ISMEAR==-1)THEN
        CALL GET_EF_FUN()
      ELSE
        WRITE(0,'(" ISMEAR=",I2)')KPT%ISMEAR; STOP ' ERROR: UNSUPPORTED ISMEAR!'
      ENDIF
      IF((GP%MYRANK.EQ.GP%MASTER).AND.(OMP_GET_THREAD_NUM()==0)) WRITE(IO,'(" GUTZ FERMI LEVEL=",F16.8)')BND%EF
      RETURN
!
      END SUBROUTINE GUTZ_FERMI
!
!*******************************************************************************
      SUBROUTINE GUTZ_FERMI_TETRA_W2K()
      INTEGER,PARAMETER::NW=250000
      INTEGER :: IW(NW)
      REAL(gq),ALLOCATABLE::EB(:,:,:),E_(:,:),WEIGHT_(:)
      INTEGER,ALLOCATABLE::NEHELP(:,:)
      REAL(gq)ELECN,EF
      INTEGER NKPT,NEMAX,NSPIN,JSPIN,ISO,NBMAX
      INTEGER IK,ISP,N1,N2
!
      NKPT=KPT%DIM; NEMAX=BND%NMAX; NSPIN=BND%NSPIN; ISO=BND%ISO
      ALLOCATE(NEHELP(NKPT,2))
      NEHELP(:,1)=BND%NE(1,:); NEHELP(:,2)=NEHELP(:,1)
      ALLOCATE(EB(NEMAX,NKPT,2))
      EB=3._gq
      DO IK=1,NKPT; EB(1:BND%NE(1,IK),IK,1:NSPIN)=BND%EK(1:BND%NE(1,IK),IK,1:NSPIN); ENDDO
      IF(NSPIN==1)EB(:,:,2)=EB(:,:,1)
      ELECN=BND%NELET-1.D-10-BND%NELET_FROZEN
      JSPIN=BND%ISPIN
      ALLOCATE(E_(NEMAX*JSPIN,NKPT)); E_=0
      E_(1:NEMAX,:) = EB(:,:,1)
      IF(JSPIN==2)THEN
        E_(1+NEMAX:,:) = EB(:,:,2)
      ENDIF
      ALLOCATE(WEIGHT_(NEMAX*JSPIN*NKPT)); WEIGHT_=0
      CALL DOS(NEMAX*JSPIN,NKPT,E_,WEIGHT_,ELECN/2.d0*ISO*JSPIN,EF,IW,NW)
      CALL EWEIGH(EF,WEIGHT_,NEMAX,NKPT,JSPIN,NEHELP,EB,NBMAX)
      DO ISP=1,NSPIN; DO IK=1,NKPT
        N1=(IK-1)*NEMAX*JSPIN+(ISP-1)*NEMAX+1
        N2=(IK-1)*NEMAX*JSPIN+ ISP   *NEMAX
        BND%FERWE(:,IK,ISP)=WEIGHT_(N1:N2)*BND%RSPO  ! Convention
      ENDDO; ENDDO
      ENG%TS2=0; BND%EF=EF
      DEALLOCATE(EB,E_,WEIGHT_,NEHELP)
      RETURN
!
      END SUBROUTINE GUTZ_FERMI_TETRA_W2K
!
!*******************************************************************************
      SUBROUTINE BND_MODIFY_FROZEN()
      INTEGER IK,I,ISP,NTOP,RSPO
!
      IF(BND%N_FROZEN==0)RETURN
      RSPO=3-MAX(BND%ISO,BND%NSPIN)
      DO IK=1,KPT%DIM
        NTOP=BND%NE(3,IK)
        DO I=NTOP-BND%N_FROZEN*BND%ISO/2+1,NTOP; DO ISP=1,BND%NSPIN
          IF(ABS(BND%EK(I,IK,ISP)-30._gq)>1.E-6_gq)THEN
            WRITE(0,'(" FETAL ERROR IN BND_MODIFY_FROZEN (!=30): FIXED BAND LEVEL = ",2F12.4)')BND%EK(I,IK,ISP)
            STOP
          ENDIF
          BND%FERWE(I,IK,ISP)=REAL(BND%NELET_FROZEN,gq)/BND%N_FROZEN*RSPO*KPT%WT(IK)
        ENDDO; ENDDO
      ENDDO
      RETURN
!
      END SUBROUTINE BND_MODIFY_FROZEN
!
!--------------------------------------------------------------------
      SUBROUTINE SET_FERMI_WEIGHT(MU)
      IMPLICIT NONE
      REAL(gq) MU
! LOCAL
      INTEGER ISP,IKP,IB,RSPO
      REAL(gq) DT
!
      BND%FERWE=0
      RSPO=3-MAX(BND%ISO,BND%NSPIN)
      DO ISP=1,BND%NSPIN; DO IKP=1,KPT%DIM; DO IB=1,BND%NE(1,IKP)
        DT=(BND%EK(IB,IKP,ISP)-MU)/KPT%DELTA
        SELECT CASE(KPT%ISMEAR)
        CASE(-1)
          BND%FERWE(IB,IKP,ISP)=FERMI_FUN(DT)
        CASE(0)
          BND%FERWE(IB,IKP,ISP)=GAUSS_FUN(DT)
        CASE DEFAULT
          STOP ' ERROR IN NELECT_FUN: KPT%ISMEAR NOT SUPPORTED!'
        END SELECT
        BND%FERWE(IB,IKP,ISP)=BND%FERWE(IB,IKP,ISP)*RSPO*KPT%WT(IKP)
      ENDDO; ENDDO; ENDDO
      RETURN
!
      END SUBROUTINE SET_FERMI_WEIGHT
!
!*******************************************************************************
      SUBROUTINE GUTZ_FERMI_TETRA_VSP()
!
!
      END SUBROUTINE GUTZ_FERMI_TETRA_VSP
!
!
      END MODULE BANDSTRU
!
!***********************************************************************
      SUBROUTINE GET_EF_FUN()
      USE GPREC; USE BANDSTRU
      IMPLICIT NONE
      REAL(gq) EMIN,EMAX
      REAL(gq),PARAMETER::EPS=1.E-10_gq
      REAL(gq),EXTERNAL::RTBIS,DIF_NELE_FUN
!
      BND%FERWE=0
      IF(ABS(BND%NELET-BND%NELET_FROZEN)<1.E-6_gq)THEN
        BND%EF=0; RETURN
      ENDIF
      EMIN=MINVAL(BND%EK); EMAX=MAXVAL(BND%EK)
      BND%EF=RTBIS(DIF_NELE_FUN,EMIN,EMAX,EPS)
!
      SELECT CASE(KPT%ISMEAR)
      CASE(-1)
        CALL CALC_ENTROPY_FERMI()
      CASE(0)
        CALL CALC_CORRECTION_GAUSS()
      END SELECT
      RETURN
!
      END SUBROUTINE GET_EF_FUN
!
!***********************************************************************
      SUBROUTINE CALC_CORRECTION_GAUSS()
      USE GPREC; USE GCONSTANT; USE BANDSTRU
      IMPLICIT NONE
      INTEGER IKP,ISP,IB,RSPO
      REAL(gq) ETA,DE
!
      ETA=0
      RSPO=3-MAX(BND%ISO,BND%NSPIN)
      DO ISP=1,BND%NSPIN; DO IKP=1,KPT%DIM; DO IB=1,BND%NE(1,IKP)
        DE=(BND%EK(IB,IKP,ISP)-BND%EF)/KPT%DELTA
        DE=DE*DE
        IF(DE<15._gq)ETA=ETA+0.5*KPT%DELTA*EXP(-DE)*KPT%WT(IKP)*RSPO
      ENDDO; ENDDO; ENDDO
      ETA=-ETA*2._gq/SQRT(PI)
      ENG%TS2=ETA/2._gq
      RETURN
!
      END SUBROUTINE CALC_CORRECTION_GAUSS
!
!***********************************************************************
      SUBROUTINE CALC_ENTROPY_FERMI()
      USE GPREC; USE GCONSTANT; USE BANDSTRU
      IMPLICIT NONE
      INTEGER IKP,ISP,IB,RSPO
      REAL(gq) ENTR,FOCC,F1,F2,EINT
!
      ENTR=0
      RSPO=3-MAX(BND%ISO,BND%NSPIN)
      DO ISP=1,BND%NSPIN; DO IKP=1,KPT%DIM
      DO IB=1,BND%NE(1,IKP)
        FOCC=BND%FERWE(IB,IKP,ISP)/RSPO/KPT%WT(IKP)
        F1=FOCC; F2=1._gq-FOCC
        IF(F1>0._gq.AND.F2>0._gq)THEN
          EINT=F1*LOG(F1)+F2*LOG(F2)
          ENTR=ENTR+EINT*KPT%WT(IKP)*RSPO
        ENDIF
      ENDDO; ENDDO; ENDDO
      ENG%TS2=KPT%DELTA*ENTR
      RETURN
!
      END SUBROUTINE CALC_ENTROPY_FERMI
!
!***********************************************************************
      FUNCTION DIF_NELE_FUN(MU)
      USE GPREC; USE BANDSTRU; USE GUTIL
      IMPLICIT NONE
      REAL(gq) MU,DIF_NELE_FUN
! LOCAL
      REAL(gq) NELET
!
      CALL SET_FERMI_WEIGHT(MU)
      NELET=SUM(BND%FERWE)
      DIF_NELE_FUN=NELET-(BND%NELET-BND%NELET_FROZEN)
      RETURN
!
      END FUNCTION DIF_NELE_FUN
!
!--------------------------------------------------------------------
      FUNCTION RTBIS(FUNC,X1,X2,TOL)
      USE GPREC
      IMPLICIT NONE
      REAL(gq) X1,X2,TOL,DX,XMID,FMID,F,RTBIS
      INTEGER J
      INTEGER,PARAMETER::JMAX=500
      REAL(gq),EXTERNAL::FUNC
!
      FMID=FUNC(X2)
      F=FUNC(X1)
      IF(F*FMID.GE.0._gq) GOTO 900
      IF(F.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5_gq
        XMID=RTBIS+DX
        FMID=FUNC(XMID)
        IF(FMID.LE.0._gq)RTBIS=XMID
        IF(ABS(FMID)<TOL)THEN
          RTBIS=XMID
          RETURN
        ENDIF
 11   CONTINUE
      GOTO 910
 900  if(abs(FMID).lt.TOL)then
        RTBIS=X2
        return
      else if(abs(F).lt.TOL)then
        RTBIS=X1
        return
      endif
      WRITE(0,'(" FMID AND F=",2F16.8)')FMID,F
      STOP ' ERROR IN RTBIS: Root must be bracketed for bisection.'
 910  STOP ' ERROR IN RTBIS: too many bisections in rtbis'
!
      END FUNCTION RTBIS

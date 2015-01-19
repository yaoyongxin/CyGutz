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
      MODULE LOCALHAMIL
      USE gprec; USE SPARSE; USE FOCKSTATE; USE CORRORB; USE GUTIL; USE gconstant; USE GMPI
      IMPLICIT NONE
!
      TYPE LOCAL_HAMIL
        INTEGER DIM,LGPRJ,LDIAPJ,NI
        TYPE(ZBD_MATRIX) :: H,U ! Onsite Coulomb interaction
        TYPE(ZBD_MATRIX) :: EVEC ! local Hamiltonian eigen vector
        REAL(gq),POINTER :: EVAL(:) ! LOCAL MANY-BODY EIGEN VALUE
        TYPE(SECTOR1)    :: SEF_N,SEC_N,SEC_S,SEC_L,SEC_J,SEC_SZ,SEC_JZ
        INTEGER,POINTER :: ID_FSC_N(:) ! N-Block index for C
        TYPE(ZBD_MATRIX) :: S2,L2,J2,SZ,SP,LZ,LP,JZ,JP
        INTEGER          :: OCC(2)
      END TYPE LOCAL_HAMIL
!
      CONTAINS
! SUBROUTINES
!******************************************************************************
      SUBROUTINE MAP_HL_FS(FS,HL)
      TYPE(FOCK_STATE) FS
      TYPE(LOCAL_HAMIL) HL
! LOCAL
      INTEGER N1,NST
!
      ! "MAP_HL_FS"
      N1=HL%OCC(1)-FS%OCC(1)+1
      HL%ID_FSC_N=>FS%ID_FSC_N(N1:)
      HL%DIM=FS%DIM-(FS%SEC_N%ID(N1)-FS%SEC_N%ID(1))
      CALL CP_HLFS_SEC(HL%SEF_N,FS%SEC_N,N1)
      IF(HL%LGPRJ==21)THEN
        CALL LOAD_SEC1(HL%SEC_N,HL%NI,17)
      ELSE
        HL%SEC_N=HL%SEF_N
      ENDIF
      ! "MAP_HL_FS"
      RETURN
!
      END SUBROUTINE MAP_HL_FS
!
!******************************************************************************
      SUBROUTINE CP_HLFS_SEC(SEC_A,SEC_B,N1)
      INTEGER N1
      TYPE(SECTOR1)::SEC_A,SEC_B
!
      SEC_A%DIM=SEC_B%DIM-N1+1
      ALLOCATE(SEC_A%ID(SEC_A%DIM+1), SEC_A%VAL(SEC_A%DIM))
      SEC_A%ID=SEC_B%ID(N1:)-SEC_B%ID(N1)+1; SEC_A%VAL=SEC_B%VAL(N1:)
      RETURN
!
      END SUBROUTINE CP_HLFS_SEC
!
!******************************************************************************
! LGPRJ: 1-> J; 2-> S,L; 3->S
!******************************************************************************
      SUBROUTINE SET_HL_SLJ2(CO,FS,HL)
      TYPE(CORR_ORB) CO
      TYPE (FOCK_STATE) FS
      TYPE(LOCAL_HAMIL) HL
! LOCAL
      INTEGER I1,I2,IDX,N1,N2,NDIF
      INTEGER NA2,NA1,LL,MM,LGPRJ
      REAL(gq),PARAMETER::RCUT=1.E-12_gq
      COMPLEX(gq) C2MS(CO%DIM2,CO%DIM2),ZFAC
      TYPE(ZCSR_MATRIX),ALLOCATABLE :: F(:),CC(:) ! Up: m=-l...l; Dn: m=-l...l
      TYPE(ZCSR_MATRIX)::FUFD,FUFU,FDFD,SP,SZ,LP,LZ,JP,JZ,J2
!
      LGPRJ=HL%LGPRJ
      IF(LGPRJ==0)RETURN
      ! "SET_HL_SLJ2"
      NA2=CO%DIM2; NA1=CO%DIM
      NDIF=HL%OCC(1)-FS%OCC(1)
      N1=FS%SEC_N%ID(NDIF+1); N2=FS%SEC_N%ID(FS%SEC_N%DIM+1)-1
      ALLOCATE(F(NA2),CC(NA2))
      C2MS=TRANSPOSE(CONJG(CO%C2N)) ! Generally can be complex
      DO I1=1,NA2; CALL SM_DCOOZCSR(FS%C(I1),CC(I1)); ENDDO
      DO I1=1,NA2; IDX=0
      !(0,'(" CALC F_MS, IORB=",I3)')I1
      DO I2=1,NA2
        ZFAC=C2MS(I2,I1)
        IF(ABS(ZFAC)<RCUT)CYCLE
        IF(IDX==0)THEN
          CALL ZCSR_COPY(CC(I2),F(I1))
          F(I1)%A=F(I1)%A*ZFAC
        ELSE
          CALL ZCSR_APLSB_SK('N',F(I1),ZFAC,CC(I2))
        ENDIF
        IDX=IDX+1
      ENDDO; ENDDO ! I2,I1
      DO I1=1,NA2; CALL DEALLOC_ZCSR(CC(I1)); ENDDO
      DEALLOCATE(CC)
      LL=(NA1-1)/2
! S_z=0.5*(c+_u c_u - c+_d c_d); S^+=c+_u c_d
      DO I1=1,NA1
        !(0,'(" CALC S, IORB=",I3)')I1
        MM=I1-LL-1
        CALL ZCSR_AMUB_SK('C',F(I1),F(I1+NA1),C=FUFD)
        CALL ZCSR_AMUB_SK('C',F(I1),F(I1),C=FUFU)
        CALL ZCSR_AMUB_SK('C',F(I1+NA1),F(I1+NA1),C=FDFD)
        IF(N1>1)THEN
          CALL ZCSR_SUBMAT(FUFD,N1,N2,N1,N2)
          CALL ZCSR_SUBMAT(FUFU,N1,N2,N1,N2)
          CALL ZCSR_SUBMAT(FDFD,N1,N2,N1,N2)
        ENDIF
        IF(I1==1)THEN
          ZFAC=DCMPLX(.5_gq)
          CALL ZCSR_COPY(FUFU,SZ); SZ%A=SZ%A*ZFAC
          CALL ZCSR_COPY(FUFD,SP)
          CALL ZCSR_COPY(FUFU,LZ); LZ%A=LZ%A*MM
        ELSE
          ZFAC=DCMPLX(.5_gq)
          CALL ZCSR_APLSB_SK('N',SZ,ZFAC,FUFU)
          ZFAC=DCMPLX(1._gq)
          CALL ZCSR_APLSB_SK('N',SP,ZFAC,FUFD)
          ZFAC=DCMPLX(MM)
          CALL ZCSR_APLSB_SK('N',LZ,ZFAC,FUFU)
        ENDIF
        ZFAC=DCMPLX(-.5_gq)
        CALL ZCSR_APLSB_SK('N',SZ,ZFAC,FDFD)
        ZFAC=DCMPLX(MM)
        CALL ZCSR_APLSB_SK('N',LZ,ZFAC,FDFD)
        CALL DEALLOC_ZCSR(FUFD); CALL DEALLOC_ZCSR(FUFU); CALL DEALLOC_ZCSR(FDFD)
      ENDDO ! I1
!
! Many Body Angular Momentum operators, Schirmer, J. and Cederbaum, L.S., Phys.Rev.A (16)1575
! http://en.wikipedia.org/wiki/Ladder_operator
      DO I1=1,NA1-1
        !(0,'(" CALC L, IORB=",I3)')I1
        MM=I1-LL-1
        ZFAC=SQRT(DBLE((LL-MM)*(LL+MM+1)))
        CALL ZCSR_AMUB_SK('C',F(I1+1),F(I1),C=FUFU)
        CALL ZCSR_AMUB_SK('C',F(I1+1+NA1),F(I1+NA1),C=FDFD)
        IF(N1>1)THEN
          CALL ZCSR_SUBMAT(FUFU,N1,N2,N1,N2)
          CALL ZCSR_SUBMAT(FDFD,N1,N2,N1,N2)
        ENDIF
        IF(I1==1)THEN
          CALL ZCSR_COPY(FUFU,LP); LP%A=LP%A*ZFAC
        ELSE
          CALL ZCSR_APLSB_SK('N',LP,ZFAC,FUFU)
        ENDIF
        CALL ZCSR_APLSB_SK('N',LP,ZFAC,FDFD)
        CALL DEALLOC_ZCSR(FUFU); CALL DEALLOC_ZCSR(FDFD)
      ENDDO ! I1
!
! Output Sz,Sp; Lz,Lp
!      CALL ZCSRTOZBM(SZ,HL%SZ,HL%SEF_N%DIM,HL%SEF_N%ID,0)
!      CALL ZCSRTOZBM(SP,HL%SP,HL%SEF_N%DIM,HL%SEF_N%ID,0)
!      CALL ZCSRTOZBM(LZ,HL%LZ,HL%SEF_N%DIM,HL%SEF_N%ID,0)
!      CALL ZCSRTOZBM(LP,HL%LP,HL%SEF_N%DIM,HL%SEF_N%ID,0)
!      WRITE(6,'("SZ")')
!      CALL ZBM_WRT2(HL%SZ,6)
!      WRITE(6,'("S+")')
!      CALL ZBM_WRT2(HL%SP,6)
!      WRITE(6,'("LZ")')
!      CALL ZBM_WRT2(HL%LZ,6)
!      WRITE(6,'("L+")')
!      CALL ZBM_WRT2(HL%LP,6)
!      STOP
!
      IF(LGPRJ==3)THEN
        CALL ZCSRTOZBM(SZ,HL%SZ,HL%SEF_N%DIM,HL%SEF_N%ID,0)
        CALL ZCSRTOZBM(SP,HL%SP,HL%SEF_N%DIM,HL%SEF_N%ID,0)
      ENDIF
! S 
      IF(LGPRJ==1.OR.LGPRJ==2.OR.LGPRJ==3.OR.LGPRJ==21)THEN
        CALL JZJP_TO_J2(SZ,SP,J2)
        CALL ZCSRTOZBM(J2,HL%S2,HL%SEF_N%DIM,HL%SEF_N%ID,0)
        CALL DEALLOC_ZCSR(J2)
        IF(LGPRJ==1.OR.LGPRJ==2.OR.LGPRJ==21)THEN
          CALL ZBM_DUMP(HL%S2,HL%NI,0,11)
          CALL DEALLOC_ZBM(HL%S2)
        ENDIF
      ENDIF
! Sz
      IF(LGPRJ==4)THEN
        CALL ZCSRTOZBM(SZ,HL%SZ,HL%SEF_N%DIM,HL%SEF_N%ID,0)
      ENDIF
! L
      IF(LGPRJ==1.OR.LGPRJ==2.OR.LGPRJ==21)THEN
        CALL JZJP_TO_J2(LZ,LP,J2)
        CALL ZCSRTOZBM(J2,HL%L2,HL%SEF_N%DIM,HL%SEF_N%ID,0)
        CALL DEALLOC_ZCSR(J2)
        CALL ZBM_DUMP(HL%L2,HL%NI,0,12)
        CALL DEALLOC_ZBM(HL%L2)
      ENDIF
! J
      IF(LGPRJ==1.OR.LGPRJ==2.OR.LGPRJ==21)THEN
        ZFAC=DCMPLX(1._gq)
        CALL ZCSR_APLSB_SK('N',SZ,ZFAC,LZ,C=JZ)
        CALL ZCSR_APLSB_SK('N',SP,ZFAC,LP,C=JP)
        CALL JZJP_TO_J2(JZ,JP,J2)
        CALL ZCSRTOZBM(J2,HL%J2,HL%SEF_N%DIM,HL%SEF_N%ID,0)
        CALL ZCSRTOZBM(JZ,HL%JZ,HL%SEF_N%DIM,HL%SEF_N%ID,0)
        CALL ZCSRTOZBM(JP,HL%JP,HL%SEF_N%DIM,HL%SEF_N%ID,0)
        CALL DEALLOC_ZCSR(J2); CALL DEALLOC_ZCSR(JZ); CALL DEALLOC_ZCSR(JP)
        CALL ZBM_DUMP(HL%J2,HL%NI,0,13)
      ENDIF
      CALL DEALLOC_ZCSR(SZ); CALL DEALLOC_ZCSR(SP)
      CALL DEALLOC_ZCSR(LZ)
      IF(NA1>1)CALL DEALLOC_ZCSR(LP)
!
      ! "SET_HL_SLJ2"
      RETURN
!
      END SUBROUTINE SET_HL_SLJ2
!
!******************************************************************************
      SUBROUTINE JZJP_TO_J2(SZ,SP,S2)
      TYPE(ZCSR_MATRIX)::SZ,SP,S2
! LOCAL
      COMPLEX(gq) ZFAC
      TYPE(ZCSR_MATRIX)::SX,SY,ZBUF
!
      ZFAC=DCMPLX( 1._gq)
      CALL ZCSR_APLSB_SK('C',SP,ZFAC,SP,C=SX); SX%A=SX%A/2
      ZFAC=DCMPLX(-1._gq)
      CALL ZCSR_APLSB_SK('C',SP,ZFAC,SP,C=SY); SY%A=SY%A/(-2*DCMPLX(0._gq,1._gq))
      CALL ZCSR_AMUB_SK('N',SX,SX,C=S2)
      CALL DEALLOC_ZCSR(SX)
      CALL ZCSR_AMUB_SK('N',SY,SY,C=ZBUF)
      CALL DEALLOC_ZCSR(SY)
      ZFAC=DCMPLX( 1._gq)
      CALL ZCSR_APLSB_SK('N',S2,ZFAC,ZBUF)
      CALL DEALLOC_ZCSR(ZBUF)
      CALL ZCSR_AMUB_SK('N',SZ,SZ,C=ZBUF)
      CALL ZCSR_APLSB_SK('N',S2,ZFAC,ZBUF)
      CALL DEALLOC_ZCSR(ZBUF)
      RETURN
!
      END SUBROUTINE JZJP_TO_J2
!
!******************************************************************************
      SUBROUTINE SET_HLOC(CO,FS,HL,LSKIPU)
      TYPE(CORR_ORB) CO
      TYPE(FOCK_STATE) FS
      TYPE(LOCAL_HAMIL) HL
      LOGICAL LSKIPU
! LOCAL
      INTEGER DIM2,N1,NBASE,NBK
      INTEGER ISEC,IFS,BS1,SGN1,BS2,SGN2,BS3,SGN3,BS4,SGN4,IBS
      INTEGER IAS1,IAS2,IAS3,IAS4
      COMPLEX(gq) FAC
      INTEGER,POINTER::ID(:)
      COMPLEX(gq),POINTER::H(:,:),U(:,:)
!
      ! "SET_HLOC"
      DIM2=CO%DIM2
      NBK=HL%SEF_N%DIM; ID=>HL%SEF_N%ID
      CALL ALLOC_ZBM(HL%H,NBK,ID,0,0)
      IF(.NOT.LSKIPU) CALL ALLOC_ZBM(HL%U,NBK,ID,0,0)
!
      DO ISEC=1,NBK
      NBASE=ID(ISEC)-1+FS%DIM-HL%DIM
      N1=ID(ISEC+1)-ID(ISEC)
      H=>HL%H%BK(ISEC)%A
      IF(.NOT.LSKIPU) U=>HL%U%BK(ISEC)%A
! ONE-BODY TERM: C2+C1
      DO IFS=1,N1
      DO IAS1=1,DIM2
      BS1=FS%BS(NBASE+IFS); SGN1=1
      CALL A_STAT(BS1,IAS1-1,.FALSE.,SGN1)
      IF(SGN1.EQ.0)CYCLE
      DO IAS2=IAS1,DIM2
      FAC=CO%EL0(IAS2,IAS1)
      IF(ABS(FAC).LT.SMALL)CYCLE
      BS2=BS1; SGN2=SGN1
      CALL A_STAT(BS2,IAS2-1,.TRUE.,SGN2)
      IF(SGN2.EQ.0)CYCLE
      IBS=FS%IBS(BS2+1)-NBASE
      IF(IBS.LE.0)CYCLE
      H(IBS,IFS)=H(IBS,IFS)+FAC*SGN2
      IF(IAS1.EQ.IAS2)CYCLE
      H(IFS,IBS)=H(IFS,IBS)+CONJG(FAC*SGN2)
      ENDDO; ENDDO; ENDDO ! IAS2,IAS1,IFS
! TWO-BODY TERM: C4+C3+C2C1
      DO IFS=1,N1
      DO IAS1=1,DIM2
      BS1=FS%BS(NBASE+IFS); SGN1=1
      CALL A_STAT(BS1,IAS1-1,.FALSE.,SGN1)
      IF(SGN1.EQ.0)CYCLE
      DO IAS2=IAS1+1,DIM2  ! Take care of 1/2
      BS2=BS1; SGN2=SGN1
      CALL A_STAT(BS2,IAS2-1,.FALSE.,SGN2)
      IF(SGN2.EQ.0)CYCLE
      DO IAS3=1,DIM2
      BS3=BS2; SGN3=SGN2
      CALL A_STAT(BS3,IAS3-1,.TRUE.,SGN3)
      IF(SGN3.EQ.0)CYCLE
      DO IAS4=1,DIM2
      FAC=CO%V2H(IAS4,IAS3,IAS1,IAS2)
      IF(ABS(FAC).LT.SMALL)CYCLE
      BS4=BS3; SGN4=SGN3
      CALL A_STAT(BS4,IAS4-1,.TRUE.,SGN4)
      IF(SGN4.EQ.0)CYCLE
      IBS=FS%IBS(BS4+1)-NBASE
      IF(IBS.LE.0)CYCLE
      H(IBS,IFS)=H(IBS,IFS)+FAC*SGN4
      IF(.NOT.LSKIPU) U(IBS,IFS)=U(IBS,IFS)+FAC*SGN4
      ENDDO; ENDDO; ENDDO; ENDDO; ENDDO
      ENDDO ! ISEC
!
      NULLIFY(H,ID)
      IF(.NOT.LSKIPU) NULLIFY(U)
      ! "SET_HLOC"
      RETURN
!
      END SUBROUTINE SET_HLOC
!
!****************************************************************************
! LGPRJ: 1-> J; 2-> S,L; 3->S; 4->Sz
!****************************************************************************
      SUBROUTINE DIAG_HL_SLJ2(HL)
      TYPE(LOCAL_HAMIL) HL
! LOCAL
      INTEGER I1,I2,IDX,LGPRJ
      INTEGER NA2,NA1,LL,MM
      REAL(gq),PARAMETER :: RCUT=1.E-6_gq
!
      LGPRJ=HL%LGPRJ
      IF(LGPRJ==0)RETURN
      ! "DIAG_HL_SLJ2"
      ! OUT_SEC_A(0,HL%SEC_N,'SEC_N ',1)
      IF(LGPRJ==1.OR.LGPRJ==2)THEN
      CALL DIAG_J2(HL%SEC_N,HL%SEC_J,HL%J2,HL%EVEC,-1)
      CALL DEALLOC_ZBM(HL%J2)
      ! OUT_SEC_A(0,HL%SEC_J,'SEC_J ',2)
      ENDIF
      IF(LGPRJ==3)THEN
      CALL DIAG_J2(HL%SEC_N,HL%SEC_S,HL%S2,HL%EVEC,-1)
      CALL DEALLOC_ZBM(HL%S2)
      ! OUT_SEC_A(0,HL%SEC_S,'SEC_S ',2)
      ENDIF
      IF(LGPRJ==4)THEN
      CALL DIAG_J2(HL%SEC_N,HL%SEC_SZ,HL%SZ,HL%EVEC,-1)
      CALL DEALLOC_ZBM(HL%SZ)
      ! OUT_SEC_A(0,HL%SEC_SZ,'SEC_SZ',1)
      ENDIF
      ! "DIAG_HL_SLJ2"
      RETURN
!
      END SUBROUTINE DIAG_HL_SLJ2
!
!****************************************************************************
      SUBROUTINE SET_VNJ_JZ(JZ,JP,SEC_J,SEC_JZ,HL)
      TYPE(SECTOR1) SEC_J,SEC_JZ
      TYPE(ZBD_MATRIX) JZ,JP
      TYPE(LOCAL_HAMIL) HL
! LOCAL
      INTEGER HRDIM
!
      ! "SET_VNJ_JZ"
      CALL ZBM_UHAU(JZ,HL%EVEC); CALL ZBM_UHAU(JP,HL%EVEC)
      HRDIM=SEC_J%ID(SEC_J%DIM+1)-SEC_J%ID(1)
      CALL DIAG_JZJP_ES(SEC_J,SEC_JZ,JZ,JP,HL%EVEC,HRDIM)
      CALL DEALLOC_ZBM(JZ); CALL DEALLOC_ZBM(JP)
      ! "SET_VNJ_JZ"
      RETURN
!
      END SUBROUTINE SET_VNJ_JZ
!
!****************************************************************************
      SUBROUTINE OUT_SEC_A(IO,SEC_A,CHSEC,MODE)
      INTEGER IO,MODE
      CHARACTER*6 CHSEC
      TYPE(SECTOR1)::SEC_A
! LOCAL
      INTEGER I,IMAX,J,NBK
      INTEGER NMAX_BK
!
      NBK=SEC_A%DIM
      WRITE(IO,'(1X,A6,"%DIM=",I4)')CHSEC,NBK
      DO I=1,NBK,10; IMAX=MIN(NBK,I+9)
        WRITE(IO,'(1X,A6,"%ID      =",11I10)')CHSEC,SEC_A%ID(I:IMAX+1)
        WRITE(IO,'(1X,A6,"%VAL     =",\)')CHSEC
        IF(MODE==1)THEN
          WRITE(IO,'(5X,10F10.4)')SEC_A%VAL(I:IMAX)
        ELSE
          WRITE(IO,'(5X,10F10.4)')SQRT(SEC_A%VAL(I:IMAX)+0.25_gq)-0.5_gq
          WRITE(IO,'(1X,A5,"%DEG/MULT=",\)')CHSEC
          WRITE(IO,'(5X,10F10.4)')((SEC_A%ID(J+1)-SEC_A%ID(J))/(1+2*(SQRT(SEC_A%VAL(J)+0.25_gq)-0.5_gq)),J=I,IMAX)
        ENDIF
      ENDDO
      NMAX_BK=0
      DO I=1,SEC_A%DIM-1; NMAX_BK=MAX(SEC_A%ID(I+1)-SEC_A%ID(I),NMAX_BK); ENDDO
      WRITE(IO,'(" NMAX_BK=",I6)')NMAX_BK
      RETURN
!
      END SUBROUTINE OUT_SEC_A
!
!****************************************************************************
      SUBROUTINE OUT_SEC(HL)
      TYPE(LOCAL_HAMIL) HL
!
      IF(HL%LGPRJ==1.OR.HL%LGPRJ==21)THEN
      OPEN(GIU,FILE=TRIM(ADJUSTL(FILE_NAME(HL%NI,0,14))),STATUS='REPLACE')
      WRITE(GIU,'(" NOCC_STR,NOCC_END:")')
      WRITE(GIU,'(2I7)')HL%OCC
      WRITE(GIU,'(" SEC_N%ID:")')
      WRITE(GIU,'(10I7)')HL%SEC_N%ID
      WRITE(GIU,'(" SEC_J%DIM:")')
      WRITE(GIU,'(I7)')HL%SEC_J%DIM
      WRITE(GIU,'(" SEC_J%ID:")')
      WRITE(GIU,'(10I7)')HL%SEC_J%ID
      WRITE(GIU,'(" SEC_J%VAL:")')
      WRITE(GIU,'(10F7.2)')SQRT(HL%SEC_J%VAL+0.25_gq)-0.5_gq
      CLOSE(GIU)
      ENDIF
      RETURN
!
      END SUBROUTINE OUT_SEC
!
!****************************************************************************
      SUBROUTINE DIAG_J2(SEC_A,SEC_B,J2,EVEC,MODE)
      TYPE(SECTOR1)::SEC_A,SEC_B
      TYPE(ZBD_MATRIX)::J2,EVEC
      INTEGER MODE
! LOCAL
      INTEGER NBK,HLDIM
      REAL(gq),PARAMETER::RCUT=1.E-6_gq
      REAL(gq),ALLOCATABLE::EVAL(:)
!
      NBK=SEC_A%DIM; HLDIM=SEC_A%ID(NBK+1)-1
      ALLOCATE(EVAL(HLDIM)); EVAL=0
      CALL BK_DIAG(NBK,SEC_A%ID,J2,EVEC,EVAL,HLDIM,1,MODE)
      CALL SET_EIGEN_SPACE(SEC_A,SEC_B,EVAL,HLDIM,RCUT)
      DEALLOCATE(EVAL)
      RETURN
!
      END SUBROUTINE DIAG_J2
!
!******************************************************************************
      SUBROUTINE SET_EIGEN_SPACE(SEC_A,SEC_B,EVAL,NDIM,RCUT)
      INTEGER NDIM
      REAL(gq) RCUT,EVAL(NDIM)
      TYPE(SECTOR1)::SEC_A,SEC_B
! LOCAL
      INTEGER I,ISUM,IBK,N1,N2
      REAL(gq) RSUM
      INTEGER,ALLOCATABLE::IDX(:)
!
      !"SET_EIGEN_SPACE"
      ALLOCATE(IDX(NDIM+1)); IDX(1)=SEC_A%ID(1); ISUM=1
      DO IBK=1,SEC_A%DIM
      N1=SEC_A%ID(IBK); N2=SEC_A%ID(IBK+1)-1
      DO I=N1,N2
        IF(I.EQ.IDX(1))CYCLE
        IF(I.NE.N1.AND.ABS(EVAL(I)-EVAL(I-1)).LT.RCUT.AND.RCUT.GT.0)CYCLE
        ISUM=ISUM+1
        IDX(ISUM)=I
      ENDDO; ENDDO
      IDX(ISUM+1)=N2+1
      SEC_B%DIM=ISUM
      ALLOCATE(SEC_B%ID(ISUM+1))
      SEC_B%ID=IDX(1:ISUM+1)
      ALLOCATE(SEC_B%VAL(ISUM))
      CALL AVG_EIGEN_SPACE(ISUM,SEC_B%ID,SEC_B%VAL,EVAL)
      DEALLOCATE(IDX)
      !"SET_EIGEN_SPACE"
      RETURN
!
      END SUBROUTINE SET_EIGEN_SPACE
!
!******************************************************************************
      SUBROUTINE AVG_EIGEN_SPACE(NDIM,ID,EL,EVAL)
      INTEGER NDIM,ID(NDIM+1)
      REAL(gq) EL(NDIM),EVAL(ID(NDIM+1)-1)
! LOCAL
      INTEGER I
      REAL(gq) RES
!
      !"AVG_EIGEN_SPACE"
      DO I=1,NDIM
        RES=SUM(EVAL(ID(I):ID(I+1)-1))/(ID(I+1)-ID(I))
        EVAL(ID(I):ID(I+1)-1)=RES
        EL(I)=RES
      ENDDO
      !"AVG_EIGEN_SPACE"
      RETURN
!
      END SUBROUTINE AVG_EIGEN_SPACE
!
!******************************************************************************
! Start from |N,J,-J> and apply Jp to get |Gamma> (State mixing avoided)
!******************************************************************************
      SUBROUTINE DIAG_JZJP_ES(SEC_A,SEC_B,JZ,JP,EVEC,HRDIM)
      TYPE(SECTOR1)::SEC_A,SEC_B
      TYPE(ZBD_MATRIX)::JZ,JP,EVEC
      INTEGER HRDIM
! LOCAL
      REAL(gq),PARAMETER::RCUT=1.E-6_gq
      REAL(gq),ALLOCATABLE::JZVAL(:)
!
      ALLOCATE(JZVAL(HRDIM)); JZVAL=0
      CALL BK_DIAG_JZJP(SEC_A%DIM,SEC_A%ID,SEC_A%VAL,JZ,JP,EVEC,JZVAL,HRDIM)
      CALL SET_EIGEN_SPACE(SEC_A,SEC_B,JZVAL,HRDIM,RCUT)
      DEALLOCATE(JZVAL)
      RETURN
!
      END SUBROUTINE DIAG_JZJP_ES
!
!******************************************************************************
      SUBROUTINE BK_DIAG(NBK,ID,H,EVEC,EVAL,HLDIM,LGPRJ,MODE)
      INTEGER NBK,ID(NBK+1),HLDIM,LGPRJ,MODE
      REAL(gq) EVAL(HLDIM)
      TYPE(ZBD_MATRIX) H,EVEC
! LOCAL
      INTEGER IBK,N1,N2,NDIM,I1
      COMPLEX(gq),ALLOCATABLE::AM(:,:)
      REAL(gq),ALLOCATABLE::WM(:)
      !
!
      DO IBK=1,NBK
        N1=ID(IBK); N2=ID(IBK+1)-1
        NDIM=N2-N1+1
        ALLOCATE(AM(NDIM,NDIM)); AM=0
        CALL ZBM_DIA_SUBMAT(H,AM,N1,N2,N1,N2,1)
        ALLOCATE(WM(NDIM)); WM=0
        IF(LGPRJ==0)THEN ! FOCK STATE
          DO I1=1,NDIM; AM(1:I1-1,I1)=0; AM(I1+1:NDIM,I1)=0; ENDDO
        ENDIF
        !
        CALL ZHEEV_('V','L',AM,WM,NDIM)
        !
        EVAL(N1:N2)=WM
        CALL ZBM_DIA_SUBMAT(EVEC,AM,N1,N2,N1,N2,MODE)
        DEALLOCATE(AM,WM)
      ENDDO
      RETURN
!
      END SUBROUTINE BK_DIAG
!
!******************************************************************************
      SUBROUTINE BK_DIAG_JZJP(NBK,ID,JVAL,JZ,JP,EVEC,JZVAL,HRDIM)
      INTEGER NBK,ID(NBK+1),HRDIM
      REAL(gq) JVAL(NBK),JZVAL(HRDIM)
      TYPE(ZBD_MATRIX) JZ,JP,EVEC
! LOCAL
      INTEGER IBK,N1,N2,NDIM,I1
      COMPLEX(gq),ALLOCATABLE::AM(:,:),AP(:,:)
      !
!
      DO IBK=1,NBK
        N1=ID(IBK); N2=ID(IBK+1)-1
        NDIM=N2-N1+1
        ALLOCATE(AM(NDIM,NDIM),AP(NDIM,NDIM)); AM=0; AP=0
        CALL ZBM_DIA_SUBMAT(JZ,AM,N1,N2,N1,N2,1)
        CALL ZBM_DIA_SUBMAT(JP,AP,N1,N2,N1,N2,1)
        !
        CALL DIAG_JZJP('L',JVAL(IBK),AM,AP,JZVAL(N1:N2),NDIM)
        !
        CALL ZBM_DIA_SUBMAT(EVEC,AM,N1,N2,N1,N2,-2) ! Update HL%EVEC
        DEALLOCATE(AM,AP)
      ENDDO
      RETURN
!
      END SUBROUTINE BK_DIAG_JZJP
!
!******************************************************************************
      SUBROUTINE CALC_DUMP_HL_NAS(CO,HL,FS,NI)
      TYPE(CORR_ORB) CO
      TYPE(LOCAL_HAMIL) HL
      TYPE(FOCK_STATE) FS
      INTEGER NI
! LOCAL
      INTEGER I,J,NBK,IBK,NROW,NCOL,IBASE,N1,N2
      TYPE(ZBD_MATRIX) :: NAS
      INTEGER,POINTER::ID(:)
      TYPE(ZCSR_MATRIX) :: NCSR
      COMPLEX(gq),ALLOCATABLE::VU(:,:)
!
      N1=HL%OCC(1)-FS%OCC(1)
      NBK=HL%SEC_N%DIM; ID=>HL%SEC_N%ID
      DO I=1,CO%SYMH%NR
        N2=CO%SYMH%DIMP(I)
        CALL ALLOC_ZBM(NAS,NBK,ID,0,0)
        IBASE=0
        DO IBK=1,NBK
          CALL CALC_NIJ_SYMGP(NCSR,FS,FS%SEC_N%ID(IBK+N1),FS%SEC_N%ID(IBK+N1+1)-1,CO%SYMH%IDP(:,1:N2,I),N2)
          NROW=HL%EVEC%BK(IBK)%NROW; NCOL=HL%EVEC%BK(IBK)%NCOL
          ALLOCATE(VU(NROW,NCOL))
          CALL ZCSRMUDEN_SK(NCSR,HL%EVEC%BK(IBK)%A,VU,NCOL)
          CALL DEALLOC_ZCSR(NCSR)
          CALL ZGEMM('C','N',NCOL,NCOL,NROW,Z1,HL%EVEC%BK(IBK)%A,NROW,VU,NROW,Z0,NAS%BK(IBK)%A,NCOL)
          DEALLOCATE(VU)
          IBASE=IBASE+NROW
        ENDDO
        CALL ZBM_SIMPLIFY(NAS,HL%SEC_J%DIM,HL%SEC_J%ID)
        CALL ZBM_DUMP(NAS,NI,I,10)
        CALL DEALLOC_ZBM(NAS)
      ENDDO
      NULLIFY(ID)
      RETURN
!
      END SUBROUTINE CALC_DUMP_HL_NAS
!
!
      END MODULE LOCALHAMIL

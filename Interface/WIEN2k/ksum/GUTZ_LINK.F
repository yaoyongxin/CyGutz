!*********************************************************
! Set BND_NE: Correlated band interval
!*********************************************************
      SUBROUTINE GUTZ_SET_EF_LDA(NUME,NKPT,PROJECTOR)
      USE COM, ONLY: emin,emax,ef,iso,iso2,BND_NE,BND_EK
      USE STRUCT, ONLY: NAT
      USE com_mpi,ONLY: nvector,vector_para,fvectors,vectors,DSUM_ALL_MPI,ISUM_ALL_MPI,myrank,master
      IMPLICIT NONE
      INTEGER NUME,NKPT,PROJECTOR
!LOCAL
      INTEGER K,IKP,N,NEN,I,II,IOS,NUM,NEMIN,NEMAX
      INTEGER ivector,itape,is
      REAL(8) SS,T,Z,E1
      LOGICAL LINEF
      CHARACTER*10 KNAME
      CHARACTER*200   FNAME
      REAL(8),PARAMETER :: RYTOEV=13.60569193_8
      INTEGER,PARAMETER :: IU=90
      
      ALLOCATE(BND_EK(NUME,NKPT)); BND_EK=0
      ALLOCATE(BND_NE(3,NKPT)); BND_NE=0

      DO ivector=1,nvector
         if (vector_para) then
            DO is=1,iso    !------ over up/dn ---------------------!
               FNAME = fvectors(ivector, is+2)
               itape=61-is
               open(itape,FILE=FNAME,STATUS='old',FORM='formatted')
            ENDDO
         endif
         do is=1,iso
            itape = 61-is
            DO I=1,NAT
               READ(itape,*); READ(itape,*)
            ENDDO
            ios=0; K=0
            DO WHILE (ios == 0)
               READ(itape,'(3e19.12,a10,2i6)',IOSTAT=ios) SS,T,Z,KNAME,N,NEn
               IF (ios /= 0) CYCLE
               k=k+1
               if (vector_para) then
                 ikp = vectors(ivector,3)+k
               else
                 ikp=k
               endif
               BND_NE(1,IKP)=NEN ! Global
               DO ii=1,nen
                  READ(itape,*) NUM,E1
                  BND_EK(II,IKP)=E1
               ENDDO
            ENDDO
         end do
         DO is=1,iso    !------ over up/dn ---------------------!
            itape=61-is
            if (vector_para) then
               close(itape)
            else
               REWIND(itape)
            endif
         END DO
      END DO
#ifdef MPI
      CALL DSUM_ALL_MPI(BND_EK,NUME*NKPT)
      CALL ISUM_ALL_MPI(BND_NE,3*NKPT)
#endif
      INQUIRE(FILE='EFLDA.INP',EXIST=LINEF)
      IF(LINEF)THEN
        CALL GUTZ_READ_FERMI(EF,IU)
        if (myrank.EQ.master) write(6,'(" Fermi level read in from EFLDA.INP")')
      ELSE
        CALL GUTZ_FERMI(ISO,ISO2); EF=EF-0.5D0 ! UNDO SHIFT
        CALL GUTZ_WRT_FERMI(EF,IU)
      ENDIF
      if (myrank.EQ.master) WRITE(6,'(" FERMI-LDA=",F8.3," Ryd,",F8.3," eV.")')EF,EF*RYTOEV
      
      EMIN=EMIN+EF; EMAX=EMAX+EF
      if (myrank.EQ.master)THEN
        WRITE(6,'(" ABSOLUTE EMIN=",F8.3," Ryd.",F8.3," eV.")')EMIN,EMIN*RYTOEV
        WRITE(6,'(" ABSOLUTE EMAX=",F8.3," Ryd.",F8.3," eV.")')EMAX,EMAX*RYTOEV
      endif     
 
      DO K=1,NKPT
      NEMIN=1; NEMAX=0
      DO II=1,BND_NE(1,K)
        IF(BND_EK(II,K).LT.EMIN) NEMIN=NEMIN+1
        IF(BND_EK(II,K).LT.EMAX) NEMAX=NEMAX+1
      ENDDO ! II
      BND_NE(2,K)=NEMIN; BND_NE(3,K)=NEMAX
      ENDDO ! K
      IF(PROJECTOR>0)THEN  ! Keep same number of bands over k
        NEMIN=MINVAL(BND_NE(2,:)); BND_NE(2,:)=NEMIN
        NEMAX=MAXVAL(BND_NE(3,:))
        DO K=1,NKPT
          BND_NE(3,K)=MIN(NEMAX,BND_NE(1,K))
        ENDDO
      ELSE
        PROJECTOR=-PROJECTOR
      ENDIF
      PROJECTOR=MOD(PROJECTOR,10)
      DEALLOCATE(BND_EK)
      RETURN
      
      END SUBROUTINE GUTZ_SET_EF_LDA

!*******************************************************************************
      SUBROUTINE GUTZ_FERMI(ISO,JSPIN)
! NELEC IS THE NUMBER OF ELECTRONS IN THIS COMPOUND                
! EF IS THE FERMI-ENERGY TO BE CALCULATED                          
      USE param,   ONLY: NUME,NKPT
      USE com_mpi, ONLY: myrank,master
      USE com,only : ts2,efmod,weigh,ef,delta,BND_NE,BND_EK
      USE bandm
      IMPLICIT NONE
      INTEGER ISO,JSPIN
! LOCAL
      INTEGER IK,IB,ISP,NSPIN,NVMAX
      REAL(8) E1
    
      CALL init_bandm(nume)

      NSPIN=MAX(1,JSPIN/ISO)
      DO ISP=1,NSPIN; DO IK=1,NKPT; DO IB=1,BND_NE(1,IK); E1=BND_EK(IB,IK)
        if(e1.gt.ebmax(IB)) ebmax(IB)=e1
        if(e1.lt.ebmin(IB)) ebmin(IB)=e1
      ENDDO; ENDDO; ENDDO

      if(efmod.eq.'GAUSS'.OR.efmod.eq.'TEMPS') then ! Just to get the energy window
!         ef=ef-int(ef)
!         if(ef.gt.0.1d0) ef=0.002d0
         ef=delta
         if (myrank.EQ.master) write(6,'(" GAUSS-SMEARING WITH",f10.5,4h Ry )') ef
         call GUTZ_fermig(ISO,JSPIN)
         ef=ef+0.5d0
      elseif(efmod.eq.'TETRA') then	
         !     tetraeder method of bloechl (for ef.lt.100  corr is used!!)
         call GUTZ_fermi_tetra(NVMAX,ISO,JSPIN)
         ts2=0.d0
      ELSE
         WRITE(0,'(" efmod=",A5)')efmod
         STOP 'Not implimented yet!'
      end if

      RETURN
    
      END SUBROUTINE GUTZ_FERMI

!************************************************************************************
      SUBROUTINE GUTZ_fermi_tetra(nbmax,ISO,JSPIN)
! NELEC IS THE NUMBER OF ELECTRONS IN THIS COMPOUND                
! EF IS THE FERMI-ENERGY TO BE CALCULATED                          
      USE com_mpi, ONLY: myrank,master
      USE param,ONLY: NKPT,NUME
      USE com,only : weigh,ef,elecn,BND_NE,BND_EK
    
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (nw=250000)                                            
      INTEGER          :: iw(nw)
      REAL*8           :: eb(nume,nkpt,2)
      REAL(8),ALLOCATABLE :: E_(:,:),WEIGHT_(:)
      INTEGER          NEHELP(NKPT,2)
      common /correct/ cordeg,icor
      !--------------------------------------------------------------------- 
      !Clas0
      !     icor switches non-linear correction  on/off
      icor=1
      if(abs(ef).ge.100.d0) icor=0
      !     cordeg is a switch to correct occupancy of degenerate states
      cordeg=-1.d-6
      NSPIN=MAX(1,JSPIN/ISO)
      if(ef.gt.0.d0)  cordeg=-ef   
      if(ef-100.d0.gt.0.d0) cordeg=-ef+100.d0
      if(ef.lt.0.d0)        cordeg=ef
      if(ef.lt.-100.d0)     cordeg=ef+100.d0
      if(abs(cordeg).gt.0.01d0)  cordeg=-1.d-6
      !
      if (myrank.EQ.master) write(6,'(" BZ-integration with TETRA-program.   icor=:",I2  )') icor
      if(cordeg.lt.0.d0 .and. (myrank.EQ.master) ) write(6,'(" Equal occupancy of degenerate states, tol=:",E8.1)') -cordeg
      
      NEMAX=MAXVAL(BND_NE(1,:))
      NEHELP(:,1)=BND_NE(1,:); NEHELP(:,2)=NEHELP(:,1)
      EB=3.d0
      DO IK=1,NKPT; EB(1:BND_NE(1,IK),IK,1)=BND_EK(1:BND_NE(1,IK),IK); ENDDO
      IF(NSPIN==1)EB(:,:,2)=EB(:,:,1)
      !.....MORE DIMENS.REPRESENTATION,IF DELTA E LESS THEN TEST1            
      ELECN=ELECN-1.D-10                                               
      !para begin
      ! ensure to mimick the proper vector file
      ! if running on multiple processors
      if (myrank.EQ.master) write(6,*)'call eord...'
      ALLOCATE(E_(nemax*jspin,NKPT)); E_=0
      call eord(E_,eb,nemax,NKPT,jspin)
      if (myrank.EQ.master) write(6,*)'call dos...' 
      ALLOCATE(WEIGHT_(nemax*jspin*NKPT)); WEIGHT_=0
      call dos(nemax*jspin,NKPT,E_,weight_,elecn/2.d0*ISO*JSPIN,d,ef,iw,nw) 
      if (myrank.EQ.master) write(6,*)'call eweigh...'
      call eweigh(ef,weight_,nemax,NKPT,jspin,nehelp,eb,nbmax)
      ELECN=ELECN+1.D-10

      ISUM=0
      DO IK=1,NKPT; DO ISP=1,JSPIN; DO IB=1,NEMAX
        ISUM=ISUM+1
        WEIGH(IK+(ISP-1)*NKPT,IB)=WEIGHT_(ISUM)*2/JSPIN/ISO
      ENDDO; ENDDO; ENDDO

      DEALLOCATE(E_,WEIGHT_)
      RETURN            
                                 
      END SUBROUTINE GUTZ_fermi_tetra

!***********************************************************************
      SUBROUTINE GUTZ_FERMIG(ISO,JSPIN)
      USE com_mpi, ONLY: myrank,master
      USE param,ONLY: NKPT,NUME
      USE kpts,ONLY: weight=>MWEIGHT
      USE com,only : weigh,ef,ts2,ELECN,BND_EK,BND_NE
    
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ztot,del1
      integer i, j, k,NSPIN, NBTOP
      INTEGER,POINTER :: NE(:)
      REAL(8),ALLOCATABLE :: E(:,:,:),WE(:,:,:)
   
      ztot=ELECN; NSPIN=MAX(1,JSPIN/ISO)
      DEL1=EF; K=NKPT
      NE=>BND_NE(1,:)
      ALLOCATE(E(NUME,K,NSPIN)); E=4.d0
      DO I=1,K; E(1:NE(I),I,1)=BND_EK(1:NE(I),I); ENDDO
      ALLOCATE(WE(NUME,K,NSPIN)); WE=0
      CALL FERMI_GS(ZTOT,DEL1,ISO,JSPIN,NSPIN,E,WE,NUME,WEIGHT,NKPT,3.9D0,1.D-4,EF,TS2,NBTOP)
      if (myrank.eq.master) write(6,154) TS2
      DO J=1,NSPIN; DO I=1,K; WEIGH((J-1)*NKPT+I,1:NE(I))=WE(1:NE(I),I,J); ENDDO; ENDDO
      DEALLOCATE(E,WE); NULLIFY(NE)
      return
154   format(/,9x,'BANDENERGY CORRECTION:',f10.6)

      end SUBROUTINE GUTZ_FERMIG

!***********************************************************************
      SUBROUTINE FERMI_GS(ZTOT,DELTA,ISO,ISPIN,NSPIN,E,WE,NBMAX,WTK,NKPT,&
                         &ECUT,ESTEP,EF,TS2,NBTOP)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::NBMAX,NKPT,ISO,ISPIN,NSPIN
      REAL(8),INTENT(IN)::ZTOT,DELTA,ECUT,ESTEP,E(NBMAX,NKPT,NSPIN),WTK(NKPT)
      REAL(8),INTENT(OUT)::WE(NBMAX,NKPT,NSPIN),EF,TS2
      INTEGER,INTENT(OUT)::NBTOP
! LOCAL
      REAL(8) EFERMI,EMIN,EMAX,NEL,NEL0,DNNT,EF0,EF1,DE,WT,FAC,E1,ETA
      INTEGER IB,IKP,ISP,NNT,NC,RSO,RSPO
      REAL(8),PARAMETER::TOL=1.E-6_8,PI =3.141592653589793238_8
      INTEGER,PARAMETER::NCMAX=1000

      RSO=3-ISO; RSPO=3-MAX(ISO,ISPIN); NBTOP=0
      NNT=ZTOT/RSO; DNNT=ZTOT/RSO-NNT
      IF(NNT==0)THEN
        EF0=E(1,1,1)
      ELSE
        EF0=E(NNT,1,1)+DNNT*(E(NNT+1,1,1)-E(NNT,1,1))
      endif
      EFERMI=EF0

      DO NC=1,NCMAX
      WE=0; EMAX=EF0-1000._8; EMIN=EF0+1000_8
      DO ISP=1,NSPIN; DO IKP=1,NKPT; DO IB=1,NBMAX
        IF(E(IB,IKP,ISP)>ECUT)CYCLE
        E1=E(IB,IKP,ISP)
        EMIN=MIN(EMIN,E1)
        DE=(E1-EFERMI)/DELTA
        IF(DE<-3._8)THEN
          WT=2._8
        ELSEIF(DE<0._8)THEN
          WT=2-ERFC(-DE)
        ELSEIF(DE<3._8)THEN
          WT=ERFC(DE)
        ELSE
          WT=0._8
        ENDIF
        WT=WT/2*RSPO
        WE(IB,IKP,ISP)=WT*WTK(IKP)
        IF(WT>1.E-5_8)THEN
          EMAX=MAX(E1,EMAX); NBTOP=MAX(NBTOP,IB)
        ENDIF
      ENDDO; ENDDO; ENDDO
      NEL=SUM(WE)
! special treatment of the first iteration step
      IF(NC==1)THEN
        EFERMI=EFERMI+ESTEP
        NEL0=NEL
        CYCLE
      ENDIF
! all other steps
      IF(ABS(ZTOT-NEL)<TOL)GOTO 100 ! LOCATED EF
      EF1=EFERMI
      FAC=(NEL0-NEL)/(ZTOT-NEL)
      IF(ABS(FAC)>.1_8)THEN
        EFERMI=EFERMI+(EF0-EFERMI)/FAC
        EF0=EF1; NEL0=NEL
      ELSE
        EF0=EF1; NEL0=NEL
        IF((ZTOT-NEL)<0)THEN
          EFERMI=EFERMI-ESTEP
        ELSE
          EFERMI=EFERMI+ESTEP
        ENDIF
      ENDIF
      ENDDO ! NC

      WRITE(0,*) 'ERROR IN FERMI_GS: FERMILEVEL NOT CONVERGED!'; STOP
100   CONTINUE
      EF=EFERMI
      EMAX=EMAX+0.0001_8
      EMIN=EMIN-0.0001_8
! Caluculate energycorrection caused by gaussian smearing
! Start loop over k-points (for spinpol.systems divide by 2
      ETA=0
      DO IKP=1,NKPT; DO ISP=1,NSPIN; DO IB=1,NBMAX
        IF(E(IB,IKP,ISP)>ECUT)CYCLE
        DE=(E(IB,IKP,ISP)-EFERMI)/DELTA
        DE=DE*DE
        IF(DE<15._8)ETA=ETA+0.5*DELTA*EXP(-DE)*WTK(IKP)
      ENDDO; ENDDO; ENDDO
      ETA=-ETA*2._8/SQRT(PI)/ISPIN
      TS2=ETA/2._8
      RETURN

      END SUBROUTINE FERMI_GS

!****************************************************************************
      SUBROUTINE GUTZ1_WRT(NAT)
      IMPLICIT NONE
      INTEGER NAT
! LOCAL
      INTEGER,PARAMETER :: IU=90

      OPEN(IU,FILE='GUTZ1.INP',STATUS='REPLACE')
      WRITE(IU,*)NAT,1
      WRITE(IU,'(" NAT,LUNIT!")')
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ1_WRT

!****************************************************************************
      SUBROUTINE GUTZ2_WRT(ISO,ISO2,NUME,NKPT)
      IMPLICIT NONE
      INTEGER ISO,ISO2,NUME,NKPT
! LOCAL
      INTEGER,PARAMETER :: IU=90

      OPEN(IU,FILE='GUTZ2.INP',STATUS='REPLACE')
      WRITE(IU,*)ISO,ISO2,NUME,NKPT
      WRITE(IU,'(" ISO,ISO2,NUME,NKPT")')
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ2_WRT

!****************************************************************************
      SUBROUTINE GUTZ3_WRT(IORD,IZ,TAU)
      IMPLICIT NONE
      INTEGER IORD
      INTEGER IZ(3,3,IORD)
      REAL(8) TAU(3,IORD)
! LOCAL
      INTEGER I
      INTEGER,PARAMETER :: IU=90

      OPEN(IU,FILE='GUTZ3.INP',STATUS='REPLACE')
      WRITE(IU,*)IORD
      WRITE(IU,'(3(3I6,/),3F6.2,/)')(IZ(:,:,I),TAU(:,I),I=1,IORD)
      WRITE(IU,*)' IORD, (IZ(:,:,I),TAU(:,I),I=1,IORD)'
      RETURN

      END SUBROUTINE GUTZ3_WRT

!*****************************************************************
! Spherical Harmonics->Projector
!*****************************************************************
      SUBROUTINE GUTZ4_WRT(CFX,MAXDIM2,NORBITALS,NINDO)
      IMPLICIT NONE
      INTEGER MAXDIM2,NORBITALS,NINDO(NORBITALS)
      COMPLEX(8) CFX(MAXDIM2,MAXDIM2,NORBITALS,NORBITALS)
! LOCAL
      INTEGER I1,I2,NINDO_
      COMPLEX(8) C2N(MAXDIM2,MAXDIM2)
      CHARACTER(40) FMT
      INTEGER,PARAMETER :: IU=90

      WRITE(FMT,'(A,I3,A)')"(",MAXDIM2,"F20.16)"
      DO I1=1,NORBITALS; DO I2=1,NORBITALS; IF(I1.EQ.I2)CYCLE
        IF(SUM(ABS(CFX(:,:,I1,I2))).LT.1.E-6_8)CYCLE
        WRITE(*,'(" NOT IMPLIMENTED: CFX HAS NON-ZERO CROSS TERMS OVER NORBITALS!")'); STOP
      ENDDO; ENDDO

      OPEN(IU,FILE='GUTZ4.INP',STATUS='REPLACE')
      WRITE(IU,*)MAXDIM2,NORBITALS
      DO I1=1,NORBITALS; NINDO_=NINDO(I1)
        C2N=0
        C2N(1:NINDO_,1:NINDO_)=TRANSPOSE(CONJG(CFX(1:NINDO_,1:NINDO_,I1,I1)))
        WRITE(IU,FMT)REAL(C2N); WRITE(IU,FMT)AIMAG(C2N); WRITE(IU,*)
      ENDDO
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ4_WRT

!****************************************************************************
      SUBROUTINE GUTZ5_WRT(NKP,WT,NELET,NE,KX,KY,KZ,KNAME,EFMOD,DELTA)
      IMPLICIT NONE
      INTEGER NKP
      REAL(8) WT(NKP),NELET,DELTA
      INTEGER :: NE(3,NKP)
      REAL(8) :: KX(NKP),KY(NKP),KZ(NKP)
      CHARACTER*10::KNAME(NKP)
      CHARACTER*5::EFMOD
! LOCAL
      INTEGER I,ISMEAR
      CHARACTER::FILE_NAME*512
      INTEGER,PARAMETER :: IU=90

      OPEN(IU,FILE='GUTZ5.INP',STATUS='REPLACE')
      WRITE(IU,*)NKP
      WRITE(IU,*)WT
      IF(EFMOD=='GAUSS')THEN
        ISMEAR=0
      ELSEIF(EFMOD=='TETRA')THEN
        ISMEAR=-5
      ELSEIF(EFMOD=='TEMPS')THEN
        ISMEAR=-1
      ELSE
        WRITE(0,'(" EFMOD=",A6)')EFMOD; STOP ' ERROR: UNSUPPORTED EFMOD!'
      ENDIF
      WRITE(IU,*)ISMEAR,DELTA
      WRITE(IU,*)NELET
      WRITE(IU,*)NE
      INQUIRE(UNIT=14,NAME=FILE_NAME)
      FILE_NAME=FILE_NAME(INDEX(FILE_NAME,'/',.TRUE.)+1:)
      WRITE(IU,*)TRIM(FILE_NAME)
      WRITE(IU,*)NKP
      DO I=1,NKP
        WRITE(IU,'(3F16.12,A10)') KX(I),KY(I),KZ(I),KNAME(I)
      ENDDO
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ5_WRT

!****************************************************************************
      SUBROUTINE GMPI_WRT()
      USE com_mpi
      IMPLICIT NONE
! LOCAL
      INTEGER I
      INTEGER,PARAMETER :: IU=90
      CHARACTER*7 STR1

      WRITE(STR1,'(I7)')MYRANK
      OPEN(IU,FILE='GMPI_'//TRIM(ADJUSTL(STR1))//'.INP',STATUS='REPLACE')
      WRITE(IU,'(3(I3,2X),L1,I5)')MYRANK,NPROCS,MASTER,VECTOR_PARA,NVECTOR
      IF(VECTOR_PARA)THEN
        DO I=1,NVECTOR
          WRITE(IU,'(3I10)')VECTORS(I,1:3)
        ENDDO
      ENDIF
      WRITE(IU,'(" MYRANK,NPROCS,MASTER,VECTOR_PARA,NVECTOR,VECTORS")')
      CLOSE(IU)
      RETURN

      END SUBROUTINE GMPI_WRT

!****************************************************************************
      SUBROUTINE GUTZ_READ_FERMI(EF,IU)
      INTEGER IU
      REAL(8) EF

      OPEN(IU,FILE='EFLDA.INP',STATUS='OLD')
      READ(IU,*)EF
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ_READ_FERMI

!****************************************************************************
      SUBROUTINE GUTZ_WRT_FERMI(EF,IU)
      INTEGER IU
      REAL(8) EF

      OPEN(IU,FILE='EFLDA.OUT',STATUS='REPLACE')
      WRITE(IU,*)EF,' ! EF_LDA'
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ_WRT_FERMI

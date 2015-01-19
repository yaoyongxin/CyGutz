#include "symbol.inc"
!**********************************************************************
! LDAG_READER!
! THIS SUBROUTINE READS THE REQUIRED  PARAMETERS FROM THE INCAR FILE
! LDAG     = .TRUE. OR .FALSE.
! LDAGL   L-QUANTUM NUMBER ON WHICH U ACTS (ONE VALUE FOR EACH SPECIES)
!**********************************************************************
      SUBROUTINE LDAG_READER(NTYP,IU5,IU0)
      USE VASPXML; USE TB_MODULE; USE prec
      IMPLICIT NONE

      INTEGER IU5,IU0
      INTEGER ITYP,NTYP
      INTEGER IDUM,N,IERR,LL
      REAL(q) RDUM
      COMPLEX(q) CDUM
      LOGICAL LOPEN,LDUM
      CHARACTER*1 :: CHARAC

      LOPEN=.FALSE.
      OPEN(UNIT=IU5,FILE='INCAR',STATUS='OLD')
      TB%IU=81
      ALLOCATE(TYS(NTYP))

      TB%KEY=0
      CALL RDATAB(LOPEN,'INCAR',IU5,'LCMR','=','#',';','I', &
     &   TB%KEY,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) THEN
            WRITE(IU0,*)'ERROR READING ITEM ''LCMR'' FROM FILE INCAR.'
         ENDIF
      ENDIF
      CLOSE(IU5)
      RETURN

      END SUBROUTINE LDAG_READER

!=======================================================================
! INITIALIZE LDAG: QOCAR
!=======================================================================
      SUBROUTINE INITIALIZE_LDAG(T_INFO,INFO,P,IU0,IU6,WDES,KPOINTS,LMDIM,NEDOS)
      USE PSEUDO; USE POSCAR; USE WAVE; USE PREC; USE mkpoints; USE base; USE CONSTANT
      USE TB_MODULE
      IMPLICIT NONE
      TYPE (TYPE_INFO)   T_INFO
      TYPE (info_struct) INFO
      TYPE (POTCAR)      P(T_INFO%NTYP)
      TYPE (WAVEDES) :: WDES
      TYPE (wavespin):: W 
      TYPE (kpoints_struct) :: KPOINTS
      INTEGER IU0,IU6,LMDIM,NEDOS
! LOCAL
      INTEGER MAXDIM2

      IF(TB%KEY==0)RETURN
      TB%Y00=1._q/2/SQRT(PI) 
      TB%NTYP = T_INFO%NTYP; TB%NIONS= T_INFO%NIONS
      TB%ITYP =>T_INFO%ITYP; TB%NITYP=>T_INFO%NITYP
      TB%IU0=IU0; TB%IU6=IU6; TB%LMDIM=LMDIM
      TB%RSPIN =WDES%RSPIN ; TB%NKPTS =WDES%NKPTS; TB%ISPIN =WDES%ISPIN
      TB%NRPLWV=WDES%NRPLWV; TB%NPROD =WDES%NPROD; TB%NBANDS=WDES%NBANDS
      TB%NPSNL =NPSNL
      TB%NELECT=INFO%NELECT; TB%NUP_DOWN=INFO%NUP_DOWN
      TB%NEDOS=NEDOS
      CALL INITIALIZE_LDAG_ELEMENT(P)
      CALL GUTZ1_WRT(TB%NIONS)
      CALL GUTZ2_WRT(WDES%NRSPINORS,WDES%ISPIN,WDES%NBANDS,TB%NKPTS)
      CALL GUTZ3_WRT()
      MAXDIM2=MAXVAL(TYS(:)%LMDIM*WDES%ISPIN)
      CALL GUTZ5_WRT(TB%NBANDS,MAXDIM2,WDES%NKPTS,WDES%WTKPT,TB%NELECT,KPOINTS%ISMEAR,KPOINTS%SIGMA)
      CALL INITIALIZE_LDAG_TQB()
      CALL INITIALIZE_LDAG_ATOM(P)
      MAXDIM2=MAXDIM2/WDES%ISPIN*2
      CALL GUTZ4_WRT(MAXDIM2,TB%NIONS)
      RETURN

      END SUBROUTINE INITIALIZE_LDAG

!=======================================================================
      SUBROUTINE CALC_LOCAL_PROJ_BND(LMDIM,INFO,CQIJ,W,MODE)
      USE prec; USE base; USE wave; USE TB_MODULE
      IMPLICIT NONE
      INTEGER            LMDIM
      TYPE (info_struct) INFO          ! INFO structure of VASP
      TYPE (wavespin):: W
      OVERLAP CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      INTEGER MODE
! Local
      TYPE (wavedes1)    WDES1
      INTEGER ISP,NK,IQ,NI,IQL,NB
      INTEGER NPL
      GDEF,ALLOCATABLE       :: CPROW(:)  ! stores the projected coefficients of wave functions Q

      IF(TB%KEY<MODE)RETURN
      ALLOCATE(CPROW(TB%NPROD)); CPROW=0

      KPOINTS: DO NK=1,TB%NKPTS
      CALL SETWDES(W%WDES,WDES1,NK)
      WDES1%NBANDS=1    ! USED THIS ONLY HERE NOT QUITE CLEAN
      NPL=WDES1%NGVECTOR
! CALC <PSI|AO>
      DO IQ=1,QB%DIM
      CPROW=0; NI=QB%NI(IQ)
      CALL OVERL(WDES1,INFO%LOVERL,LMDIM,CQIJ(1,1,1,1),QB%CPROJ(1,IQ,NK),CPROW(1))
      DO ISP=1,TB%ISPIN; IQL=QB%IQL(IQ,ISP)
      DO NB=1,TB%NBANDS
        CALL WSWP(QB%CPTWFP(:,IQ,NK),W%CPTWFP(:,NB,NK,ISP), &
                & CPROW             ,W%CPROJ (:,NB,NK,ISP), &
                & NPL,W%WDES%NPROD,TYZ(NI)%PSIA(NB,IQL,NK)) !<PSI|AO>
      ENDDO; ENDDO ! ISP,NB
      ENDDO ! IQ
      ENDDO KPOINTS
      DEALLOCATE(CPROW)
      RETURN

      END SUBROUTINE CALC_LOCAL_PROJ_BND

!=======================================================================
      SUBROUTINE QB_U_TO_W(WL,WR,MODE)
      USE prec; USE TB_MODULE; USE wave
      IMPLICIT NONE
      TYPE (wavespin):: WL,WR
      INTEGER MODE
! LOCAL
      INTEGER NK,NPL,NBANDS,NQ,NPROD
      COMPLEX(q),ALLOCATABLE::CTMP(:,:)
      COMPLEX(q),PARAMETER :: Z1=(1._q,0._q),Z0=(0._q,0._q)

      IF(TB%KEY/=MODE)RETURN
      NPL=WL%WDES%NRPLWV; NBANDS=TB%NBANDS; NQ=QB%DIM
      ALLOCATE(CTMP(NPL,NQ))
      DO NK=1,WL%WDES%NKPTS
        CALL ZGEMM('N','N',NPL,NQ,NBANDS,Z1,WR%CPTWFP(:,:,NK,1),NPL,QB%PSIA(:,:,NK),NBANDS,Z0,CTMP,NPL)
        WL%CPTWFP(:,1:NQ,NK,1)=CTMP
      ENDDO
      NPROD=WL%WDES%NPROD
      DEALLOCATE(CTMP); ALLOCATE(CTMP(NPROD,NQ))
      DO NK=1,WL%WDES%NKPTS
        CALL ZGEMM('N','N',NPROD,NQ,NBANDS,Z1,WR%CPROJ(:,:,NK,1),NPROD,QB%PSIA(:,:,NK),NBANDS,Z0,CTMP,NPROD)
        WL%CPROJ (:,1:NQ,NK,1)=CTMP
      ENDDO
      DEALLOCATE(CTMP)
      RETURN

      END SUBROUTINE QB_U_TO_W

!=======================================================================
      SUBROUTINE DUMP_BNDU(W)
      USE prec; USE TB_MODULE; USE wave; USE WAVEBASIS
      IMPLICIT NONE
      TYPE (wavespin):: W
! LOCAL
      INTEGER ISP,NK,NI
      INTEGER,PARAMETER :: IU=90

      IF(TB%KEY==0)RETURN
      OPEN(IU,FILE='BNDU_0.INP',STATUS='REPLACE',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      DO ISP=1,TB%ISPIN; DO NK=1,TB%NKPTS
        WRITE(IU)REAL(W%CELEN(:,NK,ISP),q)
        IF(TB%KEY==2.OR.TB%KEY==12)THEN
          WRITE(IU)W_BASIS%H_CMR (:,:,NK,ISP)
          WRITE(IU)W_BASIS%H_FOCK(:,:,NK,ISP)
        ENDIF
        WRITE(IU)QB%PSIA(:,:,NK)
      ENDDO; ENDDO
      CLOSE(IU)
      RETURN

      END SUBROUTINE DUMP_BNDU

!=======================================================================
! Assume serial and ISPIN=1
!=======================================================================
      SUBROUTINE UPDATE_BANDS(W,ETB,FNAME)
      USE prec; USE TB_MODULE; USE wave
      IMPLICIT NONE
      TYPE (wavespin):: W
      REAL(q) ETB
      CHARACTER*4 FNAME
! LOCAL
      INTEGER NK,NBTOT,NB1,NB2,NBAND,NPL
      REAL(q) NELE
      COMPLEX(q),ALLOCATABLE::U(:,:)
      INTEGER,PARAMETER :: IU=90

      IF(TB%KEY==0)RETURN
      IF(TB%ISPIN/=1)THEN
        STOP' ERROR: ISPIN/=1!'
      ENDIF
      OPEN(IU,FILE=FNAME//'_0.DAT',STATUS='OLD',FORM="UNFORMATTED",ACCESS="SEQUENTIAL")
      READ(IU)ETB
      NELE=0; W%FERWE=0
      DO NK=1,TB%NKPTS
        NPL=W%WDES%NGVECTOR(NK)
        READ(IU)NBTOT,NB1,NB2
        NBAND=NB2-NB1+1
        READ(IU)W%FERWE(1:NBTOT,NK,1)
        NELE=NELE+SUM(W%FERWE(:,NK,1))*TB%RSPIN*W%WDES%WTKPT(NK)
        ALLOCATE(U(NBAND,NBAND))
        READ(IU)U
        CALL ZANMXBMM('N',W%CPTWFP(1:NPL     ,NB1:NB2,NK,1),U,NPL     ,NBAND)
        CALL ZANMXBMM('N',W%CPROJ (1:TB%NPROD,NB1:NB2,NK,1),U,TB%NPROD,NBAND)
        DEALLOCATE(U)
      ENDDO
      CLOSE(IU)
      IF(ABS(NELE-TB%NELECT)>1.E-6_q)THEN
        WRITE(0,'(" ERROR: TB%NELECT vs NELE=",2F12.5)')TB%NELECT,NELE
      ENDIF
      RETURN

      END SUBROUTINE UPDATE_BANDS

!************************************************************
      SUBROUTINE ZANMXBMM(TRANSB,A,B,N,M)
      USE prec
      IMPLICIT NONE
      CHARACTER*1 TRANSB
      INTEGER N,M
      COMPLEX(q) A(N,M),B(M,M)
! LOCAL
      COMPLEX(q),ALLOCATABLE::A_(:,:)
      COMPLEX(q),PARAMETER :: Z1=(1._q,0._q),Z0=(0._q,0._q)

      ALLOCATE(A_(N,M)); A_=0
      CALL ZGEMM('N',TRANSB,N,M,M,Z1,A,N,B,M,Z0,A_,N)
      A=A_
      DEALLOCATE(A_)
      RETURN

      END SUBROUTINE ZANMXBMM

!=======================================================================
      SUBROUTINE INITIALIZE_LDAG_ELEMENT(P)
      USE prec; USE pseudo; USE TB_MODULE
      IMPLICIT NONE
      TYPE (POTCAR)      P(TB%NTYP)
! LOCAL
      INTEGER NT,IR,IL,LL,IADD
      INTEGER LDIM,NMAX,LDIMR
      REAL(q) ALPHA(20)
      LOGICAL LIN(20)
      CHARACTER CSEL*1,LABEL*2,UNITS*1

      OPEN(TB%IU,FILE='QOCAR',STATUS='OLD')
      DO NT=1,TB%NTYP
        TYS(NT)%PSMAXN=P(NT)%PSMAXN
        READ(TB%IU,*)CSEL,LABEL,LDIM,NMAX,UNITS
        IF(UNITS.NE.'a')THEN
          WRITE(0,*)'Illegal unit:',UNITS; STOP
        ENDIF
        IF(INDEX(LABEL,P(NT)%ELEMENT,.TRUE.).EQ.0)THEN
          WRITE(*,'(A9,2(2x,A2))')' ELEMENT=',LABEL,P(NT)%ELEMENT
          WRITE(*,*)'QOCAR and POTCAR mis-match!'
          STOP
        ENDIF
        IF(NMAX.LT.P(NT)%R%NMAX)THEN
          WRITE(0,'(" NMAXT=",I4," NMAXP=",I4)')NMAX,P(NT)%R%NMAX
          WRITE(0,'(" TOT NUMBER OF RADIAL GRIDS < # IN AUG. SHERE!")')
          STOP
        ENDIF
        TYS(NT)%R%NMAX=NMAX
        READ(TB%IU,*)CSEL,LIN(1:LDIM)
        READ(TB%IU,*)CSEL,ALPHA(1:LDIM)
        LDIMR=LDIM
        DO IL=1,LDIM; IF(LIN(IL)) CYCLE; LDIMR=LDIMR-1; ENDDO
        TYS(NT)%LDIM=LDIMR
        ALLOCATE(TYS(NT)%L(LDIMR),TYS(NT)%ALPHA(LDIMR))
        IADD=0
        DO IL=1,LDIM
          IF(.NOT.LIN(IL)) CYCLE
          IADD=IADD+1
          TYS(NT)%ALPHA(IADD)=ALPHA(IL)
        ENDDO
        ALLOCATE(TYS(NT)%R%R(NMAX),TYS(NT)%PW(NMAX,LDIMR),TYS(NT)%QW(0:NPSNL,LDIMR))
        TYS(NT)%R%R=0; TYS(NT)%PW=0; TYS(NT)%QW=0
        IADD=0
        DO IL=1,LDIM
          READ(TB%IU,*)CSEL,LL
          IF(LIN(IL))THEN
            IADD=IADD+1
            TYS(NT)%L(IADD)=LL
            DO IR=1,NMAX; READ(TB%IU,*)TYS(NT)%R%R(IR),TYS(NT)%PW(IR,IADD); ENDDO
          ELSE
            DO IR=1,NMAX; READ(TB%IU,*); ENDDO
          ENDIF
        ENDDO ! IL
        READ(TB%IU,*)
        TYS(NT)%LMDIM=SUM(2*TYS(NT)%L+1)
        ALLOCATE(TYS(NT)%LMCHAR(TYS(NT)%LMDIM))
        IADD=0
        DO IL=1,LDIMR
          LL=TYS(NT)%L(IL)
          SELECT CASE (LL)
           CASE (0)
            TYS(NT)%LMCHAR(IADD+1)='  s'
           CASE (1)
            TYS(NT)%LMCHAR(IADD+1)=' py'
            TYS(NT)%LMCHAR(IADD+2)=' pz'
            TYS(NT)%LMCHAR(IADD+3)=' px'
           CASE (2)
            TYS(NT)%LMCHAR(IADD+1)='dxy'
            TYS(NT)%LMCHAR(IADD+2)='dyz'
            TYS(NT)%LMCHAR(IADD+3)='dz2'
            TYS(NT)%LMCHAR(IADD+4)='dxz'
            TYS(NT)%LMCHAR(IADD+5)='dx2'
           CASE (3)
            TYS(NT)%LMCHAR(IADD+1)='f-3'
            TYS(NT)%LMCHAR(IADD+2)='f-2'
            TYS(NT)%LMCHAR(IADD+3)='f-1'
            TYS(NT)%LMCHAR(IADD+4)='f-0'
            TYS(NT)%LMCHAR(IADD+5)='f+1'
            TYS(NT)%LMCHAR(IADD+6)='f+2'
            TYS(NT)%LMCHAR(IADD+7)='f+3'
           CASE DEFAULT
            STOP 'LMCHAR: LM OUT OF RANGE!'
          END SELECT
          IADD=IADD+2*LL+1
        ENDDO
        CALL CHK_AO_NORM(NT,P(NT))
        CALL QOGREP(P(NT),NT)
      ENDDO ! NT
      CLOSE(TB%IU)
      QB%DIM=SUM(TYS(:)%LMDIM*TB%NITYP(:))
      WRITE(TB%IU6,'(" QB%DIM=",I4)')QB%DIM
      RETURN

      END SUBROUTINE INITIALIZE_LDAG_ELEMENT

!=======================================================================
      SUBROUTINE INITIALIZE_LDAG_TQB()
      USE prec; USE TB_MODULE
      IMPLICIT NONE
! LOCAL
      INTEGER NQ,NRPLWV,NPROD,NKPTS

      NQ=QB%DIM; NRPLWV=TB%NRPLWV; NPROD=TB%NPROD
      NKPTS=TB%NKPTS

      ALLOCATE(QB%NI(NQ),QB%IQL(NQ,TB%ISPIN),QB%CPTWFP(NRPLWV,NQ,NKPTS), &
              &QB%CPROJ(NPROD,NQ,NKPTS), &
              &QB%PSIA(TB%NBANDS,NQ*TB%ISPIN,NKPTS))
      QB%NI=0; QB%CPTWFP=0; QB%CPROJ=0; QB%IQL=0; QB%PSIA=0
      RETURN

      END SUBROUTINE INITIALIZE_LDAG_TQB

!=======================================================================
      SUBROUTINE CHK_AO_NORM(NT,P)
      USE prec; USE pseudo; USE radial; USE TB_MODULE
      IMPLICIT NONE
      TYPE (POTCAR) P
      INTEGER NT
! LOCAL
      INTEGER NMAX,LDIM,IL,LL
      REAL(q) NORM
      REAL(q),EXTERNAL::GET_RNORM

      NMAX=TYS(NT)%R%NMAX
      TYS(NT)%R%RMAX  = TYS(NT)%R%R(NMAX)
      TYS(NT)%R%REND  = TYS(NT)%R%R(NMAX)
      TYS(NT)%R%D     =(TYS(NT)%R%R(NMAX  )-TYS(NT)%R%R(NMAX-1))/ &
                       &(TYS(NT)%R%R(NMAX-1)-TYS(NT)%R%R(NMAX-2))
      TYS(NT)%R%H     =LOG(TYS(NT)%R%D)
      TYS(NT)%R%RSTART=TYS(NT)%R%R(2)/(EXP(TYS(NT)%R%H)-1)
      IF(ABS((TYS(NT)%R%H-P%R%H)/P%R%H).GT.1.E-6)THEN
        WRITE(0,*)TYS(NT)%R%R(1:2),'v.s.'
        WRITE(0,*)P%R%R(1:2)
        STOP 'MISMATCH RADIAL GRID!'
      ENDIF
      CALL SET_SIMP(TYS(NT)%R)
      LDIM=TYS(NT)%LDIM
      WRITE(TB%IU6,'(" NT=",I2)')NT
! Check original norm
      DO IL=1,LDIM
        LL=TYS(NT)%L(IL)
        NORM=GET_RNORM(LL,IL,P,NT,1)
        WRITE(TB%IU6,'(2(a3,I2),2x,a10,f6.3)')'IL=',IL,' L=',LL,'NORM_ORIG=',NORM
        TYS(NT)%PW(:,IL)=TYS(NT)%PW(:,IL)*EXP(-TYS(NT)%ALPHA(IL)*TYS(NT)%R%R)
        NORM=GET_RNORM(LL,IL,P,NT,1)
        WRITE(TB%IU6,'(5X,a3,I2,2x,a10,f6.3)')' L=',LL,'NORM_SCAL=',NORM
        TYS(NT)%PW(:,IL)=TYS(NT)%PW(:,IL)*SQRT(ABS(1/NORM))
! Target orbital might be negative due to output.
        NORM=SUM(TYS(NT)%R%SI(:)*TYS(NT)%PW(:,IL))
        IF(NORM.LT.0) TYS(NT)%PW(:,IL)=-TYS(NT)%PW(:,IL)
      ENDDO ! IL
      RETURN

      END SUBROUTINE CHK_AO_NORM

!=======================================================================
! Get <AO~|S|AO~> for radial part. L=0,1,2...
!=======================================================================
      FUNCTION GET_RNORM(L,IL,P,NT,MODE)
      USE pseudo; USE constant; USE TB_MODULE; USE prec
      IMPLICIT NONE
      INTEGER L,IL,NT,MODE
      TYPE (potcar)      P
      REAL(q) GET_RNORM
! Local
      INTEGER CH,IND,I,IS,J,CHS(2)
      INTEGER,PARAMETER::IO=77
      REAL(q) ARGSC,DIST,REM,VPS,SNORM(2)

      IS=0; SNORM=0; CHS=1
      ARGSC=NPSRNL/P%PSRMAX
      DO CH=1,P%LMAX
      IF(P%LPS(CH).NE.L)CYCLE
      IS=IS+1; CHS(IS)=CH
      IF(MODE.EQ.1.AND.IS.EQ.1)THEN
      WRITE(IO,'("# WPS NT=",I2," L=",I2)')NT,L
      DO IND=1,P%R%NMAX
        WRITE(IO,*)P%R%R(IND),P%WPS(IND,CH)
      ENDDO
      WRITE(IO,*)
      WRITE(IO,'("# WAE NT=",I2," L=",I2)')NT,L
      DO IND=1,P%R%NMAX
        WRITE(IO,*)P%R%R(IND),P%WAE(IND,CH)
      ENDDO
      WRITE(IO,*)
      ENDIF ! MODE.EQ.1.AND.IS.EQ.1
      IF(IS.GT.2)THEN
        STOP " ERROR IN GET_RNORM: IS>2!"
      ENDIF
      DO IND=1,TYS(NT)%R%NMAX
      DIST=TYS(NT)%R%R(IND)
      IF(DIST.GT.P%PSRMAX)EXIT
      I=MIN(INT(DIST*ARGSC)+1,NPSRNL-1)
      REM=DIST-P%PSPRNL(I,1,CH)
      VPS=P%PSPRNL(I,2,CH)+REM*(P%PSPRNL(I,3,CH)+REM*(P%PSPRNL(I,4,CH)+REM*P%PSPRNL(I,5,CH)))
! Ahh: PSPRNL is only the R_l(r) instead of rR_l(r)
      SNORM(IS)=SNORM(IS)+TYS(NT)%R%SI(IND)*VPS*TYS(NT)%PW(IND,IL)*DIST
      ENDDO
      ENDDO ! CH

      IF(MODE.EQ.1)THEN
      WRITE(IO,'("# UNSCALED PW NT=",I2," L=",I2)')NT,L
      DO IND=1,TYS(NT)%R%NMAX
        WRITE(IO,*)TYS(NT)%R%R(IND),TYS(NT)%PW(IND,IL)
      ENDDO
      WRITE(IO,*)
      ENDIF ! MODE.EQ.1

      GET_RNORM=SUM(TYS(NT)%R%SI(:)*TYS(NT)%PW(:,IL)**2)
      DO I=1,2; DO J=1,2
      GET_RNORM=GET_RNORM+SNORM(I)*SNORM(J)*P%QION(CHS(I),CHS(J))
      ENDDO; ENDDO

      END FUNCTION GET_RNORM

!=======================================================================
! RADIAL FOURIER TRANSFORMATION OF RW TO GET QW. GMAX=P(NT)%PSMAXN
!=======================================================================
      SUBROUTINE QOGREP(P,NT)
      USE PSEUDO; USE prec; USE TB_MODULE
      IMPLICIT NONE
      INTEGER NT
      TYPE(POTCAR) P
! LOCAL
      INTEGER NMAX

      NMAX=TYS(NT)%R%NMAX
      CALL RAD_FOURIER(TYS(NT)%PW,TYS(NT)%L,NMAX,TYS(NT)%R%R,TYS(NT)%R%SI, &
                      &P%PSMAXN,NPSNL,TYS(NT)%QW)

      END SUBROUTINE QOGREP

!=======================================================================
! General subroutine to perform radial fourier transformation
!=======================================================================
      SUBROUTINE RAD_FOURIER(RW,L,NMAX,R,SI,PSMAXN,NPSNL,QW)
      USE prec
      IMPLICIT NONE
      INTEGER NMAX,NPSNL,L
      REAL(q) PSMAXN
      REAL(q) RW(NMAX),R(NMAX),QW(0:NPSNL),SI(NMAX)
! Local
      INTEGER IG,IR
      REAL(q) GLEN,RLEN,GSTEP,BJ
      REAL(q),POINTER :: AUX(:)

! FOLLOW THE RECIPE FOR PSPRNL.
      ALLOCATE(AUX(NMAX)); AUX=0
      GSTEP=PSMAXN/NPSNL

      DO IG=1,NPSNL
      GLEN=(IG-1)*GSTEP
      DO IR=1,NMAX
        RLEN=R(IR)
        CALL SBESSEL(GLEN*RLEN,BJ,L)
        AUX(IR)=RLEN*RW(IR)*BJ
      ENDDO ! IR
      QW(IG)=SUM(AUX*SI)
      ENDDO ! IG
      QW(0)=QW(2)*(-1)**L
      DEALLOCATE(AUX)

      END SUBROUTINE RAD_FOURIER

!**********************************************************************
! INITIALIZE HDFT FOR ATOM INFO
!**********************************************************************
      SUBROUTINE INITIALIZE_LDAG_ATOM()
      USE prec; USE pseudo; USE TB_MODULE
      IMPLICIT NONE
! LOCAL
      INTEGER NIONS,NI,NT,NH,NHADD,NBASE,NH1,NH2,NHS,ISPIN,ISP,IQ

      NIONS=TB%NIONS; ISPIN=TB%ISPIN
      ALLOCATE(TYZ(NIONS))

      NHADD=0; NBASE=0
      DO NI=1,NIONS
      NT=TB%ITYP(NI)
      TYZ(NI)%LDIM   =TYS(NT)%LDIM 
      TYZ(NI)%LMDIM   =TYS(NT)%LMDIM
      NH=TYS(NT)%LMDIM
      IF(NH.EQ.0)CYCLE
      NH1=NHADD+1; NH2=NHADD+NH; NHS=NH*ISPIN
      QB%NI(NH1:NH2)=NI
      DO ISP=1,ISPIN; DO IQ=1,NH
        QB%IQL(NHADD+IQ,ISP)=IQ+(ISP-1)*NH
      ENDDO; ENDDO
      TYZ(NI)%CPTWFP=>QB%CPTWFP(:,NH1:NH2,:)
      TYZ(NI)%CPROJ =>QB%CPROJ (:,NH1:NH2,:)
      TYZ(NI)%PSIA  =>QB%PSIA(:,NBASE+1:NBASE+NHS,:)
      TYZ(NI)%PW    =>TYS(NT)%PW
      TYZ(NI)%QW    =>TYS(NT)%QW
      TYZ(NI)%L     =>TYS(NT)%L
      TYZ(NI)%LMCHAR=>TYS(NT)%LMCHAR
      NHADD=NHADD+NH
      NBASE=NBASE+NHS
      ENDDO ! NI
      RETURN

      END SUBROUTINE INITIALIZE_LDAG_ATOM

!=======================================================================
! GENERATE LCAO for one k point and calc their CPROJ
!=======================================================================
      SUBROUTINE GEN_LCAO_1K(LATT_CUR,WDES,NONL_S,NK)
      USE lattice; USE wave; USE nonl; USE constant; USE asa; USE prec
      USE TB_MODULE
      IMPLICIT NONE

      TYPE (latt)        LATT_CUR
      TYPE (nonl_struct) NONL_S        ! descriptor for non local part of PP (reciprocal space)
      TYPE (wavedes)     WDES          ! descriptor for wavefunction
      INTEGER NK
! Local
      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W1
      REAL(q) KG(3),KGX,KGY,KGZ,QW1,FAKT
      COMPLEX(q)         CFKT
      INTEGER NIS,IQ,LM,NPL,MM,L,NITYP,NIONS,IL
      INTEGER IPL,NT,NI,IPRO,LMBASE
      REAL(q),ALLOCATABLE:: YLM(:,:),KGNORM(:)
      COMPLEX(q),ALLOCATABLE:: KGREXP(:,:) ! EXP(i(K+G)R) for atoms with same type.

!=======================================================================
! number of G-Vectors on the k-Point
!=======================================================================
      FAKT= 1/SQRT(LATT_CUR%OMEGA)
      NPL= WDES%NGVECTOR(NK)

      ALLOCATE(YLM(NPL,16),KGNORM(NPL)); YLM=0; KGNORM=0 ! maxL=3

! Set up phase and YLM
      NIONS=NONL_S%NIONS
      ALLOCATE(KGREXP(NPL,NIONS)); KGREXP=0

      DO IPL=1,NPL
      KGX=WDES%IGX(IPL,NK)+WDES%VKPT(1,NK); KG(1)=KGX
      KGY=WDES%IGY(IPL,NK)+WDES%VKPT(2,NK); KG(2)=KGY
      KGZ=WDES%IGZ(IPL,NK)+WDES%VKPT(3,NK); KG(3)=KGZ
      DO NI=1,NIONS
        KGREXP(IPL,NI)=EXP(CITPI*SUM(KG(:)*NONL_S%POSION(:,NI)))
      ENDDO ! NI
      KG(:)=(KGX*LATT_CUR%B(:,1)+KGY*LATT_CUR%B(:,2)+KGZ*LATT_CUR%B(:,3))*TPI
      KGNORM(IPL)=MAX(SQRT(SUM(KG(:)**2)),1E-10_q)
      KG=KG/KGNORM(IPL)
      CALL SETYLM(3,1,YLM(IPL:IPL,:),KG(1:1),KG(2:2),KG(3:3))
      ENDDO ! IPL
! Begin generating LCAO
      NIS=0
      DO NT=1,NONL_S%NTYP
      NITYP=NONL_S%NITYP(NT)
      LMBASE=0
      DO IL=1,TYS(NT)%LDIM
        L=TYS(NT)%L(IL); CFKT=4*PI*(0._q,1._q)**L*FAKT
        DO IPL=1,NPL
        CALL INTERPOLATE4(TYS(NT)%PSMAXN,TB%NPSNL,TYS(NT)%QW(:,IL),KGNORM(IPL),QW1)
        DO MM=1,2*L+1
        LM=L**2+MM
        DO NI=1,NITYP
          TYZ(NIS+NI)%CPTWFP(IPL,MM+LMBASE,NK)=CONJG(QW1*YLM(IPL,LM)*KGREXP(IPL,NIS+NI)*CFKT)  ! |LCAO>
        ENDDO; ENDDO; ENDDO ! NI,NM,IPL
        LMBASE=LMBASE+2*L+1
      ENDDO ! IL
      ENDDO ! NT
! Calc CPROJ
      CALL PHASE(WDES,NONL_S,NK)
! Set phase for NK. ONly exp(igr) since exp(ikr) would be cancelled considering |Pi>Qij<Pj|
      CALL SETWDES(WDES,WDES1,NK)
      DO NI=1,NIONS
      DO IQ=1,TYZ(NI)%LMDIM
      W1%CPTWFP=>TYZ(NI)%CPTWFP(:,IQ,NK)
      W1%CPROJ =>TYZ(NI)%CPROJ (:,IQ,NK)
      CALL PROJ1(NONL_S,WDES1,W1)  ! Calc <P~|LCAO>
      ENDDO ! IQ
      ENDDO ! NI
      DEALLOCATE(YLM,KGNORM,KGREXP)
      RETURN

      END SUBROUTINE GEN_LCAO_1K

!++++++++++++++++++++++++++++++++++++++++++++++++++
! Micmic the storage of CRHODE for symm.
!++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE WQ2CRHODE(WQQW,LMDIM,P,T_INFO,LBACK,ISPIN)
      USE poscar; USE pseudo; USE prec
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      INTEGER LMDIM,ISPIN
      REAL(q) WQQW(LMDIM,LMDIM,T_INFO%NIONS,ISPIN)
      LOGICAL LBACK
! Local
      INTEGER NT,NI,NIS,CH,LL,LL0,LM,QLM,NM,ISP
      OVERLAP WQBUF(LMDIM,LMDIM)

      DO ISP=1,ISPIN
      NIS=0
      DO NT=1,T_INFO%NTYP; DO NI=1,T_INFO%NITYP(NT)
      NIS=NIS+1; LL0=-10; LM=0; WQBUF=WQQW(:,:,NIS,ISP); WQQW(:,:,NIS,ISP)=0
      DO CH=1,P(NT)%LMAX
      LL=P(NT)%LPS(CH); NM=2*LL+1; LM=LM+NM
      IF(LL.EQ.LL0)CYCLE
      QLM=(LL+1)**2
      IF(LBACK)THEN
        WQQW(QLM-NM+1:QLM,QLM-NM+1:QLM,NIS,ISP)=WQBUF(LM-NM+1:LM,LM-NM+1:LM)
      ELSE
        WQQW(LM-NM+1:LM,LM-NM+1:LM,NIS,ISP)=WQBUF(QLM-NM+1:QLM,QLM-NM+1:QLM)
      ENDIF
      LL0=LL
      ENDDO; ENDDO; ENDDO
      ENDDO ! ISP

      END SUBROUTINE WQ2CRHODE

!++++++++++++++++++++++++++++++++++++++++++++++++++
! Given F(NMAX) and the maximum x, xmax, four-point
! interpolation of f @ x. used in VASP
!++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE INTERPOLATE4(XMAX,NMAX,F,X,FX)
      USE prec
      IMPLICIT NONE

      INTEGER NMAX
      REAL(q) XMAX,F(0:NMAX),X,FX
! Local
      REAL(q) ARGSC,ARG,V1,V2,V3,V4,T0,T1,T2,T3,REM
      INTEGER NADDR

      ARGSC=NMAX/XMAX
      ARG=X*ARGSC+1
      NADDR=INT(ARG)
      IF(NADDR<NMAX-2)THEN
        REM=MOD(ARG,1.0_q)
        V1=F(NADDR-1)
        V2=F(NADDR  )
        V3=F(NADDR+1)
        V4=F(NADDR+2)
        T0=V2
        T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
        T2=(V1+V3-(2*V2))/2._q
        T3=(V4-V1+(3*(V2-V3)))/6._q
        FX= T0+REM*(T1+REM*(T2+REM*T3))
      ELSE
        FX=0
      ENDIF

      END SUBROUTINE INTERPOLATE4

!=======================================================================
! INITIALIZE LCAO
!=======================================================================
      SUBROUTINE INI_LCAO(LATT_CUR,WDES,NONL_S)
      USE LATTICE; USE NONL; USE WAVE; USE prec; USE TB_MODULE
      IMPLICIT NONE
      TYPE (LATT)        LATT_CUR
      TYPE (NONL_STRUCT) NONL_S        ! DESCRIPTOR FOR NON LOCAL PART OF PP (RECIPROCAL SPACE)
      TYPE (WAVEDES)     WDES          ! DESCRIPTOR FOR WAVEFUNCTION
! LOCAL 
      INTEGER NK

      IF(TB%KEY==0)RETURN
      DO NK =1,TB%NKPTS
        CALL GEN_LCAO_1K(LATT_CUR,WDES,NONL_S,NK)
      ENDDO ! NK
      RETURN
      
      END SUBROUTINE INI_LCAO

!++++++++++++++++++++++++++++++++++++++++++++++++++
! Calc <WFP|S|WF>, S=1+Q
!++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE WSWP(WF,WFP,CPROW,CPROJ,NPL,NPROD,RES)
      USE prec
      IMPLICIT NONE

      INTEGER    NPL,NPROD
      COMPLEX(q) WF(*),WFP(*),CPROJ(*),CPROW(*),RES

      RES=DOT_PRODUCT(WFP(1:NPL),WF(1:NPL))   ! <WF'|WF>
      RES=RES+DOT_PRODUCT(CPROJ(1:NPROD),CPROW(1:NPROD))   ! Q included in CPROW

      END SUBROUTINE WSWP

!****************************************************************************
      SUBROUTINE ORTH_LOC_PROJ_BND(MODE)
      USE prec; USE TB_MODULE
      IMPLICIT NONE
      INTEGER MODE
!LOCAL
      INTEGER NK
      COMPLEX(q)::SK(QB%DIM,QB%DIM)
      COMPLEX(q),PARAMETER :: Z1=(1._q,0._q),Z0=(0._q,0._q)

      IF(TB%KEY/=MODE)RETURN
      IF(TB%ISPIN/=1)THEN
        STOP ' ERROR IN ORTH_BNDU: TB%ISPIN/=1!'
      ENDIF
      DO NK=1,TB%NKPTS
        CALL ZGEMM('C','N',QB%DIM,QB%DIM,TB%NBANDS,Z1, &
                  &QB%PSIA(:,:,NK),TB%NBANDS, &
                  &QB%PSIA(:,:,NK),TB%NBANDS, &
                  &Z0,SK,QB%DIM)
        CALL ORTH_LOC_PROJ_BND1(QB%PSIA(:,:,NK),SK,TB%NBANDS,QB%DIM)
      ENDDO
      RETURN

      END SUBROUTINE ORTH_LOC_PROJ_BND

!=============================================================================
      SUBROUTINE ORTH_LOC_PROJ_BND1(U,SK,NBANDS,NAST)
      USE prec
      IMPLICIT NONE
      INTEGER NBANDS,NAST
      COMPLEX(q) U(NBANDS,NAST),SK(NAST,NAST)
! LOCAL
      INTEGER J1,J2
      COMPLEX(q) U_(NBANDS,NAST)
      COMPLEX(q),ALLOCATABLE :: SNH(:,:)

      ALLOCATE(SNH(NAST,NAST)); SNH=0
      CALL CALC_S2SNH(SK,NAST,SNH)
      U_=U; U=0
      DO J1=1,NAST; DO J2=1,NAST
        U(:,J1)=U(:,J1)+U_(:,J2)*SNH(J2,J1)
      ENDDO; ENDDO
      DEALLOCATE(SNH)
      RETURN

      END SUBROUTINE ORTH_LOC_PROJ_BND1

!=======================================================================
      SUBROUTINE ORTH_LOCAL_PROJ(LMDIM,INFO,CQIJ,W)
      USE prec; USE base; USE wave; USE TB_MODULE
      IMPLICIT NONE
      INTEGER            LMDIM
      TYPE (info_struct) INFO          ! INFO structure of VASP
      TYPE (wavespin):: W
      OVERLAP CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
! Local
      TYPE (wavedes1)    WDES1
      INTEGER ISP,NK,IQ,NB
      INTEGER NPL
      COMPLEX(q)::SK(QB%DIM,QB%DIM),SNH(QB%DIM,QB%DIM)
      COMPLEX(q)::CPTWFP(TB%NRPLWV,QB%DIM),CPROJ(TB%NPROD,QB%DIM)
      GDEF,ALLOCATABLE :: CPROW(:)  ! stores the projected coefficients of wave functions Q
      COMPLEX(q),PARAMETER :: Z1=(1._q,0._q),Z0=(0._q,0._q)

      ALLOCATE(CPROW(TB%NPROD)); CPROW=0
      KPOINTS: DO NK=1,TB%NKPTS
      CALL SETWDES(W%WDES,WDES1,NK)
      WDES1%NBANDS=1    ! USED THIS ONLY HERE NOT QUITE CLEAN
      NPL=WDES1%NGVECTOR
      DO IQ=1,QB%DIM
      CPROW=0
      CALL OVERL(WDES1,INFO%LOVERL,LMDIM,CQIJ(1,1,1,1),QB%CPROJ(1,IQ,NK),CPROW(1))
      DO NB=1,IQ
        CALL WSWP(QB%CPTWFP(:,IQ,NK),QB%CPTWFP(:,NB,NK), &
                & CPROW             ,QB%CPROJ (:,NB,NK), &
                & NPL,TB%NPROD,SK(NB,IQ)) !<NB|IQ>
        IF(NB/=IQ)SK(IQ,NB)=CONJG(SK(NB,IQ))
      ENDDO; ENDDO ! IQ
      CALL CALC_S2SNH(SK,QB%DIM,SNH)
      CPTWFP=QB%CPTWFP(:,:,NK); CPROJ=QB%CPROJ(:,:,NK)
      CALL ZGEMM('N','N',TB%NRPLWV,QB%DIM,QB%DIM,Z1,CPTWFP,TB%NRPLWV,SNH,QB%DIM,Z0,QB%CPTWFP(:,:,NK),TB%NRPLWV)
      CALL ZGEMM('N','N',TB%NPROD ,QB%DIM,QB%DIM,Z1,CPROJ ,TB%NPROD ,SNH,QB%DIM,Z0,QB%CPROJ (:,:,NK),TB%NPROD )
      ENDDO KPOINTS
      DEALLOCATE(CPROW)
      RETURN

      END SUBROUTINE ORTH_LOCAL_PROJ

!=======================================================================
      SUBROUTINE MERG_LCAO_BND(W,INFO)
      USE prec; USE TB_MODULE; USE wave; USE base
      IMPLICIT NONE
      TYPE (wavespin):: W
      TYPE (info_struct) INFO  
!LOCAL
      INTEGER ISP,NK,NBAND,I,NINS,NREP
      REAL(q) VAL(QB%DIM)
      COMPLEX(q)::SK(QB%DIM,QB%DIM)
      REAL(q),PARAMETER :: TOL=1.E-6_q
      COMPLEX(q),PARAMETER :: Z1=(1._q,0._q),Z0=(0._q,0._q)

      DO ISP=1,TB%ISPIN; DO NK=1,TB%NKPTS
      NREP=QB%DIM
100   CONTINUE
      NBAND=TB%NBANDS-NREP
      CALL ZGEMM('C','N',QB%DIM,QB%DIM,NBAND,Z1, &
                &QB%PSIA(1:NBAND,:,NK),NBAND, &
                &QB%PSIA(1:NBAND,:,NK),NBAND, &
                &Z0,SK,QB%DIM)
      CALL ZHEEV_('V','U',SK,VAL,QB%DIM)
      NINS=0
      DO I=1,QB%DIM
        IF(ABS(1._q-VAL(I))>TOL)THEN
          NINS=NINS+1
        ELSE
          EXIT
        ENDIF
      ENDDO ! I

      IF(NINS==0)CYCLE ! All well expanded
      IF(NINS<NREP)THEN
        NREP=NINS; GOTO 100
      ELSEIF(NINS>NREP)THEN
        WRITE(0,'(" N_INSERT=",I4," VS N_REPLACE=",I4)')NINS,NREP
        STOP ' ERROR IN MERG_LCAO_BND!'
      ENDIF

      CALL ZGEMM('N','N',TB%NRPLWV,NREP,QB%DIM,Z1,QB%CPTWFP(:,:,NK),TB%NRPLWV, &
                &SK(:,1:NREP),QB%DIM,Z0,W%CPTWFP(:,NBAND+1:TB%NBANDS,NK,ISP),TB%NRPLWV)
      CALL ZGEMM('N','N',TB%NPROD ,NREP,QB%DIM,Z1,QB%CPROJ (:,:,NK),TB%NPROD , &
                &SK(:,1:NREP),QB%DIM,Z0,W%CPROJ (:,NBAND+1:TB%NBANDS,NK,ISP),TB%NPROD )
      ENDDO; ENDDO ! NK,ISP
      RETURN

      END SUBROUTINE MERG_LCAO_BND

!=======================================================================
      SUBROUTINE CALC_S2SNH(S,N,SNH)
      USE prec
      IMPLICIT NONE
      INTEGER N
      COMPLEX(q) S(N,N),SNH(N,N)
!LOCAL
      INTEGER I1,I2,I3
      REAL(q) W(N)

      CALL ZHEEV_('V','U',S,W,N)
      DO I1=1,N
        IF(W(I1).LT.1.E-6_q)THEN
          WRITE(*,'(/)')
          WRITE(*,'(" FETAL ERROR: EIGEN VALUES OF S IN CALC_S2SNH TOO SMALL!")')
          WRITE(*,'(14F8.3)')W
          STOP 'CHECK CHECK ZHEEV-II IN CALC_S2SNH!'
        ENDIF
      ENDDO
! SET UP S^(-1/2)
      SNH=0
      DO I1=1,N; DO I2=1,N; DO I3=1,N
        SNH(I2,I1)=SNH(I2,I1)+S(I2,I3)*CONJG(S(I1,I3))/SQRT(W(I3))
      ENDDO; ENDDO; ENDDO
      RETURN

      END SUBROUTINE CALC_S2SNH

!****************************************************************************
      SUBROUTINE QB_SET_AS_MLWF()
      USE TB_MODULE; USE mlwf; USE full_kpoints
      IMPLICIT NONE
!LOCAL
      INTEGER num_wann,num_bands,num_bands_tot,num_kpts

      MLWF_num_wann=QB%DIM
      num_wann=QB%DIM
      num_kpts=KPOINTS_FULL%NKPX*KPOINTS_FULL%NKPY*KPOINTS_FULL%NKPZ
      num_bands_tot=TB%NBANDS; num_bands=TB%NBANDS
      ALLOCATE(U_matrix(num_wann,num_wann,num_kpts,TB%ISPIN), &
          &   U_matrix_opt(num_bands,num_wann,num_kpts,TB%ISPIN), &
          &   lwindow(num_bands,num_kpts,TB%ISPIN), &
          &   lexclude_band(num_bands_tot))
      lwindow = .TRUE.
      lexclude_band=.FALSE.
      RETURN

      END SUBROUTINE QB_SET_AS_MLWF

!****************************************************************************
      SUBROUTINE OUT_PAWSYM(P,NTYP,A,B,ROTMAP,NIOND)
      USE pseudo; USE prec; USE TB_MODULE
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER NTYP,NIOND,NPCELL
      REAL(q) A(3,3),B(3,3)
      TYPE (potcar) P(NTYP)
      INTEGER ROTMAP(NIOND,NROTK,NPCELL)
! LOCAL
      INTEGER LMAX,MMAX
      REAL(q),ALLOCATABLE::SL(:,:,:,:) ! (L,M) <- (L,MP)
      INTEGER,EXTERNAL :: MAXL
      INTEGER,PARAMETER :: IU=90

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      IF(TB%KEY==0)RETURN
      LMAX=MAXL(NTYP,P)         ! maximum l quantum number
      MMAX=(2*LMAX+1)             ! maximum m quantum number
      ALLOCATE(SL(MMAX,MMAX,0:LMAX,NROTK)); SL=0
      CALL CALC_PAWSYM_SL(SL,LMAX,MMAX,NROTK,ISYMOP,A,B)
      OPEN(IU,FILE='PAWSYM.OUT',STATUS='REPLACE')
      WRITE(IU,'(2I5)')LMAX,NROTK
      WRITE(IU,*)SL(1:,1:,0:,1:)
      WRITE(IU,'(3I5)')NIOND,NPCELL
      WRITE(IU,'(10I5)')ROTMAP
      CLOSE(IU)
      CONTINUE

      END SUBROUTINE OUT_PAWSYM
!
!****************************************************************************
      SUBROUTINE DUMP_KPT_TET(KPOINTS)
      USE prec; USE TB_MODULE; USE mkpoints
      IMPLICIT NONE
      TYPE (kpoints_struct) :: KPOINTS
! LOCAL
      INTEGER,PARAMETER :: IU=90
     
      OPEN(IU,FILE='KPTET.OUT',STATUS='REPLACE')
      WRITE(IU,*)KPOINTS%NTET,KPOINTS%VOLWGT
      WRITE(IU,*)KPOINTS%IDTET
      CLOSE(IU)
      RETURN

      END SUBROUTINE DUMP_KPT_TET

!****************************************************************************
      SUBROUTINE CALC_PAWSYM_SL(SL,LMAX,MMAX,NR,S,A,B)
      USE TB_MODULE; USE pawsym; USE prec
      IMPLICIT NONE
      INTEGER LMAX,MMAX,NR
      REAL(q) SL(MMAX,MMAX,0:LMAX,NR)
      INTEGER S(3,3,48)
      REAL(q) A(3,3),B(3,3)
!LOCAL
      INTEGER IROT

      DO IROT=1,NR
        ! setup rotation matrices for L=0,...,LMAX
        CALL SETUP_SYM_LL(MMAX,LMAX,S(1,1,IROT),SL(1,1,0,IROT),A,B)
      ENDDO
      RETURN

      END SUBROUTINE CALC_PAWSYM_SL

!****************************************************************************
      SUBROUTINE GMPI_WRT(MYRANK,NPROCS,MASTER,VECTOR_PARA,NVECTOR)
      USE TB_MODULE
      IMPLICIT NONE
      INTEGER MYRANK,NPROCS,MASTER,NVECTOR
      LOGICAL VECTOR_PARA
! LOCAL
      INTEGER I
      INTEGER,PARAMETER :: IU=90
      CHARACTER*7 STR1

      IF(TB%KEY<=0)RETURN
      WRITE(STR1,'(I7)')MYRANK
      OPEN(IU,FILE='GMPI_'//TRIM(ADJUSTL(STR1))//'.INP',STATUS='REPLACE')
      WRITE(IU,*)MYRANK,NPROCS,MASTER,VECTOR_PARA,NVECTOR
      WRITE(IU,*)' MYRANK,NPROCS,MASTER,VECTOR_PARA,NVECTOR,VECTORS'
      CLOSE(IU)
      RETURN

      END SUBROUTINE GMPI_WRT

!****************************************************************************
      SUBROUTINE GUTZ1_WRT(NAT)
      USE TB_MODULE
      IMPLICIT NONE
      INTEGER NAT
! LOCAL
      INTEGER,PARAMETER :: IU=90

      IF(TB%KEY<=0)RETURN
      OPEN(IU,FILE='GUTZ1.INP',STATUS='REPLACE')
      WRITE(IU,*)NAT,0
      WRITE(IU,'(" NAT,LUNIT!")')
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ1_WRT

!****************************************************************************
      SUBROUTINE GUTZ2_WRT(ISO,ISO2,NUME,NKPT)
      USE TB_MODULE
      USE prec
      IMPLICIT NONE
      INTEGER ISO,ISO2,NUME,NKPT
! LOCAL
      INTEGER,PARAMETER :: IU=90

      IF(TB%KEY<=0)RETURN
      OPEN(IU,FILE='GUTZ2.INP',STATUS='REPLACE')
      WRITE(IU,*)ISO,ISO2,NUME,NKPT
      WRITE(IU,'(" ISO,ISO2,NUME,NKPT")')
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ2_WRT

!****************************************************************************
      SUBROUTINE GUTZ3_WRT()
      USE TB_MODULE
      IMPLICIT NONE
! LOCAL
      INTEGER,PARAMETER :: IU=90

      IF(TB%KEY<=0)RETURN
      OPEN(IU,FILE='GUTZ3.INP',STATUS='REPLACE')
      WRITE(IU,*)0,0
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ3_WRT

!*****************************************************************
      SUBROUTINE GUTZ4_WRT(MAXDIM2,NIONS)
      USE TB_MODULE
      IMPLICIT NONE
      INTEGER MAXDIM2,NIONS
! LOCAL
      COMPLEX(q) C2N(MAXDIM2,MAXDIM2,NIONS)
      INTEGER NI,I
      INTEGER,PARAMETER :: IU=90

      IF(TB%KEY<=0)RETURN
      C2N=0
      DO NI=1,NIONS; DO I=1,TYZ(NI)%LMDIM*2
        C2N(I,I,NI)=1._q
      ENDDO; ENDDO
      OPEN(IU,FILE='GUTZ4.INP',STATUS='NEW',ERR=100)
      WRITE(IU,*)MAXDIM2,NIONS
      WRITE(IU,*)C2N
      WRITE(IU,*)' MAXDIM2,NORBITALS,C2N'
      CLOSE(IU)
100   RETURN

      END SUBROUTINE GUTZ4_WRT

!****************************************************************************
      SUBROUTINE GUTZ5_WRT(NBANDS,MAXDIM2,NKP,WT,NELET,ISMEAR,DELTA)
      USE TB_MODULE
      IMPLICIT NONE
      INTEGER NBANDS,MAXDIM2,NKP,ISMEAR
      REAL(8) WT(NKP),NELET,DELTA
! LOCAL
      INTEGER I
      INTEGER,PARAMETER :: IU=90

      IF(TB%KEY<=0)RETURN
      OPEN(IU,FILE='GUTZ5.INP',STATUS='REPLACE')
      WRITE(IU,*)MAXDIM2,NKP
      WRITE(IU,*)WT
      WRITE(IU,*)ISMEAR,DELTA
      WRITE(IU,*)NELET
      DO I=1,NKP; WRITE(IU,*)NBANDS,1,NBANDS; ENDDO ! NE
      WRITE(IU,*)'KPOINTS'
      WRITE(IU,*)0
      CLOSE(IU)
      RETURN

      END SUBROUTINE GUTZ5_WRT

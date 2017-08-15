SUBROUTINE LoMAIN(nemin,nemax,lfirst,latom,n,jatom,isym,LC,iso,icase)
  USE param
  USE struct
  USE sym2
  USE abc
  USE case
  IMPLICIT NONE
  INTEGER, intent(in) :: nemin, nemax, lfirst, latom, n, jatom, isym, LC, iso, icase
  !IMPLICIT REAL*8 (A-H,O-Z)
  COMPLEX*16 :: YL((LMAX2+1)*(LMAX2+1))
  COMPLEX*16 :: PHSHEL,PH_SPIN(2)
  COMPLEX*16,ALLOCATABLE   :: phs(:)
  COMMON /GENER/   BR1(3,3),BR2(3,3)
  REAL*8        :: BR1, BR2
  COMMON /XA/      R(NRAD),BK(3),BKROT(3),BKROT2(3),BKROT3(3),BKRLOC(3)
  REAL*8        :: R, BK, BKROT, BKROT2, BKROT3, BKRLOC
  common /loabc/  alo(0:lomax,2,nloat,nrf)
  REAL*8        :: alo
  common /lolog/  nlo,nlov,nlon,loor(0:lomax),ilo(0:lomax),lapw(0:lmax2)
  INTEGER       :: nlo, nlov, nlon, ilo
  logical       :: loor, lapw  
  DATA           IMAG/(0.0D0,1.0D0)/         
  COMPLEX*16 :: IMAG
  REAL*8     :: PI, TWOPI, ARG1, ARG2, ARG3, ARGT, ARGT2
  INTEGER    :: I, L, jlo, jneq, M1, is, NUM, M, ind_yl, irf
!-----------------------------------------------------------------------------     
  PI=ACOS(-1.0D0)
  TWOPI=2.D0*PI
  ALLOCATE(phs(nume))
  i=n-(nlo+nlon)
  DO L=0,LoMAX                ! 10
     do jlo=1,ilo(l)          ! 10
        do jneq=1,mult(jatom) ! 20
           DO M1=-l,+l        ! 25
              i=i+1  
              if (.not.(L.eq.LC)) CYCLE !goto 25
              BK(1)=BKX(I) ! BK = K+k==G
              BK(2)=BKY(I)
              BK(3)=BKZ(I)
              ! BKROT = R_a.(k+K) transforms to the reducible k-point
              CALL ROTATE (BK,TMAT(1,1,isym),BKROT)
              ! BKROT2 = R_n.R_a.(k+K), transformation from the first atom to an equivalent atom               
              CALL ROTATE (BKROT, rotij(1,1,latom), BKROT2)
              !---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
              ! BKROT3 = R_n.R_a.(k+K), but in cartesian coordinate system
              BKROT3(1)=BKROT2(1)*BR1(1,1)+BKROT2(2)*BR1(1,2)+BKROT2(3)*BR1(1,3)
              BKROT3(2)=BKROT2(1)*BR1(2,1)+BKROT2(2)*BR1(2,2)+BKROT2(3)*BR1(2,3)
              BKROT3(3)=BKROT2(1)*BR1(3,1)+BKROT2(2)*BR1(3,2)+BKROT2(3)*BR1(3,3)
              !---- BKRLOC = crotloc.R_n.R_a.(k+K),  rotates according to the user specified local coordinate system.
              CALL ROTATE (BKROT3,crotloc(1,1,icase),BKRLOC) ! BKRLOC = Rotloc.R_g.(k+K)
              !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
              CALL YLM (BKRLOC,LMAX2,YL)
              ! (R_n.R_a.(k+K)) *  R(first) * 2pi
              ARG1=BKROT2(1)*(POS(1,lfirst))*TWOPI          ! ARG1 + ARG2 + ARG3 = (R_g.(k+K)) *  R(latom) * 2pi
              ARG2=BKROT2(2)*(POS(2,lfirst))*TWOPI
              ARG3=BKROT2(3)*(POS(3,lfirst))*TWOPI
              ! ARGT = (k+K)*tau(isym) * 2pi
              ARGT=(BKX(I)*TAU(1,isym)+BKY(I)*TAU(2,isym)+ BKZ(I)*TAU(3,isym))*TWOPI
              ! ARGT2 = (R_a.(k+K)).tau_n * 2pi
              ARGT2=(BKROT(1)*tauij(1,latom)+BKROT(2)*tauij(2,latom)+BKROT(3)*tauij(3,latom))*TWOPI 
              ! PHSEHL = e^{I*2pi*( (R_a.(k+K))*tau_n + (K+k)*tau(isym) + (R_n.R_a.(k+K)*R(first)))}
              PHSHEL=EXP(IMAG*(ARG1+ARG2+ARG3+ARGT+ARGT2))
              DO is=1,iso
                 PH_SPIN(IS)=EXP((2*is-3)*IMAG*PHASE(ISYM)/2)
                 DO NUM=NEMIN,NEMAX
                    PHS(NUM)=PHSHEL*A(I,NUM,is)*PH_SPIN(IS)
                 ENDDO
                 DO M=1,2*L+1
                    ind_yl=M+L*L                                                
                    DO NUM=NEMIN,NEMAX
                       DO irf=1,nrf
                          ALM(m,num,irf,is)=ALM(m,num,irf,is)+ALo(l,is,jlo,irf)*dconjg(YL(IND_YL))*PHS(NUM)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        enddo
     enddo
  ENDDO
  return                   
END SUBROUTINE LoMAIN

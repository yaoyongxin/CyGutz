SUBROUTINE lomain(XROTLOC,is,nemin,nemax,lfirst,latom,n,jatom,alm,blm,clm)
  USE param,ONLY: LoMAX, LMAX2, nume, nloat
  USE defs, ONLY: TWO, PI, IMAG
  USE struk,ONLY: mult, rotij, tauij, BR1, POS, ROTLOC
  USE lo, ONLY: nlo, nlov, nlon, ilo, Alo, Blo, Clo
  USE xa, ONLY: BK, BKROT, BKRLOC, PHS
  USE xa3,ONLY: BKXlo, BKYlo, BKZlo, As_lo
  USE com,ONLY: 
  IMPLICIT NONE
  REAL*8,  INTENT(in) :: XROTLOC(3,3)
  INTEGER, INTENT(in) :: is, nemin,nemax,lfirst,latom,n,jatom
  COMPLEX*16, INTENT(inout) :: ALM((LMAX2+1)*(LMAX2+1),nume),BLM((LMAX2+1)*(LMAX2+1),nume),cLM((LMAX2+1)*(LMAX2+1),nume,nloat)  
  COMPLEX*16 :: YL((LMAX2+1)*(LMAX2+1))
  COMPLEX*16 :: PHSHEL                     
  INTEGER    :: i,l,m,m1,jlo,jneq,num,index
  real*8     :: k(3),krot(3),twopi,arg1,arg2,arg3,argt, BKROT3(3)
  !.initialise a,b,c of lo                                      
  TWOPI=TWO*PI
  !print *, 'nlov=', nlov
  i=nlov
  DO L=0,LoMAX
     do jlo=1,ilo(l)
        do jneq=1,mult(jatom)
           DO M1=-l,+l 
              i=i+1
              BK(1)=BKXlo(I)
              BK(2)=BKYlo(I)
              BK(3)=BKZlo(I)
              ! BKROT2 = R_n.R_a.(k+K), transformation from the first atom to an equivalent atom
              CALL ROTATE (BK,ROTIJ(1,1,LATOM),BKROT)
              !---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
              ! BKROT3 = R_n.R_a.(k+K), but in cartesian coordinate system
              BKROT3(1)=BKROT(1)*BR1(1,1)+BKROT(2)*BR1(1,2)+BKROT(3)*BR1(1,3)   
              BKROT3(2)=BKROT(1)*BR1(2,1)+BKROT(2)*BR1(2,2)+BKROT(3)*BR1(2,3)   
              BKROT3(3)=BKROT(1)*BR1(3,1)+BKROT(2)*BR1(3,2)+BKROT(3)*BR1(3,3)   
              !---- BKRLOC = crotloc.R_n.R_a.(k+K),  rotates according to the user specified local coordinate system.
              CALL ROTATE (BKROT3,XROTLOC,BKRLOC)
              !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
              CALL YLM (BKRLOC,LoMAX,YL)
              ! (R_n.R_a.(k+K)) *  R(first) * 2pi
              ARG1=BKROT(1)*(POS(1,LFIRST))*TWOPI
              ARG2=BKROT(2)*(POS(2,LFIRST))*TWOPI
              ARG3=BKROT(3)*(POS(3,LFIRST))*TWOPI
              ! ARGT = (R_a.(k+K)).tau_n * 2pi
              ARGT=(BKXlo(I)*TAUIJ(1,LATOM)+BKYlo(I)*TAUIJ(2,LATOM)+BKZlo(I)*TAUIJ(3,LATOM))*TWOPI
              ! PHSEHL = e^{I*2pi*( (R_a.(k+K))*tau_n + (K+k)*tau(isym) + (R_n.R_a.(k+K)*R(first)))}
              PHSHEL=EXP( IMAG*(ARG1+ARG2+ARG3+ARGT) )

              !WRITE(*,'(A,1x,I4,1x,3f10.5)') 'BK0=', i+n-(nlo+nlov+nlon), BK
              !print *, 'pos=', POS(:,lfirst)
              !print *, 'BK=', BKROT

              !print *, 'p=', ARG1, ARG2, ARG3, ARGT
              DO NUM=NEMIN,NEMAX
                 PHS(NUM) = PHSHEL*As_lo(I,NUM,is)   
              ENDDO
              DO M=-l,+l
                 index=l*(l+1)+m+1
                 DO NUM=NEMIN,NEMAX
                    ALM(index,num)=ALM(index,num)+Alo(l,jlo)*conjg(YL(index))*PHS(NUM)
                    BLM(index,num)=BLM(index,num)+Blo(l,jlo)*conjg(YL(index))*PHS(NUM)
                    CLM(index,num,jlo)=CLM(index,num,jlo)+Clo(l,jlo)*CONJG(YL(index))*PHS(NUM)
             
                    !if (abs(Alo(l,jlo)).gt.1e-10) WRITE(*,'(A,1x,I4,1x,I2,1x,I2,1x,f20.15,1x,2f20.15,1x,2f20.15)') 'Q', I, m+l+1, num, Alo(l,jlo), YL(index), PHS(NUM)
                    !if (abs(Blo(l,jlo)).gt.1e-10) WRITE(*,'(A,1x,I4,1x,I2,1x,I2,1x,f20.15,1x,2f20.15,1x,2f20.15)') 'Q', I, m+l+1, num, Blo(l,jlo), YL(index), PHS(NUM)
                    !if (abs(Clo(l,jlo)).gt.1e-10) WRITE(*,'(A,1x,I4,1x,I2,1x,I2,1x,f20.15,1x,2f20.15,1x,2f20.15)') 'Q', I, m+l+1, num, Clo(l,jlo), YL(index), PHS(NUM)

                    !if (num.eq.52) WRITE(*,'(A,1x,I4,1x,I2,1x,I2,1x,2f20.15,1x,2f20.15,1x,2f20.15,1x,2f20.15)') 'Q', I, m+l+1, num, Alm(index,num), Blm(index,num), Clm(index,num,1), Clm(index,num,2)

                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  return
END SUBROUTINE lomain

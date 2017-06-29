USE param
subroutine xSrm3(jspin,jri,llmax,dx,rx)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /RADFU/   RR1(NRAD,0:LMAX2,2),RR2(NRAD,0:LMAX2,2), &
                       RE1(NRAD,0:LMAX2,2),RE2(NRAD,0:LMAX2,2),  &
                       Rlo1(NRAD,0:LOMAX,2),Rlo2(NRAD,0:LMAX2,2)          
      COMMON /RINTEG/  RAA(0:lmax2,2,2),RBB(0:lmax2,2,2), &
                       RCC(0:lomax,2,2),RAB(0:lmax2,2,2), &
                       RAC(0:lomax,2,2),RBC(0:lomax,2,2)
      DIMENSION  AMEAAM(NRAD),AMEABM(NRAD),AMEBBM(NRAD)
      DIMENSION  AMECCM(NRAD),AMEACM(NRAD),AMEBCM(NRAD)
      dimension rx(nrad)
      logical          loor(0:lomax)     
      common /lolog/   nlo,nlov,nlon,loor  
       do 91 isi=1,jspin
       do 91 isj=isi,jspin
!.........................................................................
!.....CALCULATE THE MATRIX ELEMENT
!
      DO 100 L=0,LLMAX
      LK=L+1
! perform <X> over the large component only
      DO 106 IR=1,JRI
       vde=1./rx(ir)**3
      AMEAAM(IR)=RR1(IR,L,isi)*RR1(IR,L,isj)*vde
      AMEBBM(IR)=RE1(IR,L,isi)*RE1(IR,L,isj)*vde
      AMEABM(IR)=RR1(IR,L,isi)*RE1(IR,L,isj)*vde
  106 CONTINUE
!
      CALL CALI(AMEAAM,TOTAA,RX(1),DX,JRI)
      CALL CALI(AMEBBM,TOTBB,RX(1),DX,JRI)
      CALL CALI(AMEABM,TOTAB,RX(1),DX,JRI)
!
      RAA(L,isi,isj)=TOTAA
      RBB(L,isi,isj)=TOTBB
      RAB(L,isi,isj)=TOTAB
  100 CONTINUE
!
!.........................................................................
!.....CALCULATE THE MATRIX ELEMENT of LO
!
      DO 200 L=0,llmax
      if ((l.le.lomax).and.loor(l)) then
      LK=L+1
      DO 206 IR=1,JRI
       vde=1./rx(ir)**3
      AMECCM(IR)=Rlo1(IR,L,isi)*Rlo1(IR,L,isj)*vde
      AMEACM(IR)=RR1(IR,L,isi)*Rlo1(IR,L,isj)*vde
      AMEBCM(IR)=RE1(IR,L,isi)*Rlo1(IR,L,isj)*vde
  206 CONTINUE
      CALL CALI(AMECCM,TOTCC,RX(1),DX,JRI)
      CALL CALI(AMEACM,TOTAC,RX(1),DX,JRI)
      CALL CALI(AMEBCM,TOTBC,RX(1),DX,JRI)
!
      RCC(L,isi,isj)=TOTCC
      RAC(L,isi,isj)=TOTAC
      RBC(L,isi,isj)=TOTBC
      else
      RCC(L,isi,isj)=0.0d0
      RAC(L,isi,isj)=0.0d0
      RBC(L,isi,isj)=0.0d0
      endif
  200 CONTINUE
91    continue
!
      RETURN
      END

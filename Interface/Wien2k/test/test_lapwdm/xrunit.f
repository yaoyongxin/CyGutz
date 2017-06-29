USE param
subroutine xrunit(jspin,jri,llmax,dx,rx)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /RADFU/   RR1(NRAD,0:LMAX2,2),RR2(NRAD,0:LMAX2,2), &
                       RE1(NRAD,0:LMAX2,2),RE2(NRAD,0:LMAX2,2),  &
                       Rlo1(NRAD,0:LOMAX,2),Rlo2(NRAD,0:LMAX2,2)          
      COMMON /RINTEG/  RAA(0:lmax2,2,2),RBB(0:lmax2,2,2), &
                       RCC(0:lomax,2,2),RAB(0:lmax2,2,2), &
                       RAC(0:lomax,2,2),RBC(0:lomax,2,2)
      logical          loor(0:lomax)     
      common /lolog/   nlo,nlov,nlon,loor  
      dimension totaa(2),totbb(2),totab(2),totac(2),totbc(2), &
                totcc(2)
      dimension rx(nrad)
      DIMENSION  AMEAAM(NRAD,2),AMEABM(NRAD,2),AMEBBM(NRAD,2)
      DIMENSION  AMECCM(NRAD,2),AMEACM(NRAD,2),AMEBCM(NRAD,2)
      COE=1.D0
      CIN=1.0D0/137.0359895d0**2
       do 91 isi=1,jspin
       do 91 isj=1,jspin
!.........................................................................
!.....CALCULATE THE MATRIX ELEMENT
!
      DO 100 L=0,LLMAX
      DO 104 IC=1,2
      DO 106 IR=1,JRI
!      if(isi.eq.isj)then
       vde=1.d0            
!     else
!      vde=0.d0                           
!      endif
      AMEAAM(IR,IC)=RR1(IR,L,isi)*RR1(IR,L,isj)*vde
      AMEBBM(IR,IC)=RE1(IR,L,isi)*RE1(IR,L,isj)*vde
      AMEABM(IR,IC)=RR1(IR,L,isi)*RE1(IR,L,isj)*vde
  106 CONTINUE
!
      CALL CALI(AMEAAM(1,IC),TOTAA(IC),RX(1),DX,JRI)
      CALL CALI(AMEBBM(1,IC),TOTBB(IC),RX(1),DX,JRI)
      CALL CALI(AMEABM(1,IC),TOTAB(IC),RX(1),DX,JRI)
  104 CONTINUE
!
      RAA(L,isi,isj)=COE*(TOTAA(1)+CIN*TOTAA(2))
      RBB(L,isi,isj)=COE*(TOTBB(1)+CIN*TOTBB(2))
      RAB(L,isi,isj)=COE*(TOTAB(1)+CIN*TOTAB(2))
  100 CONTINUE
!
!.........................................................................
!.....CALCULATE THE MATRIX ELEMENT of LO
!
      DO 200 L=0,llmax
      if ((l.le.lomax).and.loor(l)) then
      DO 204 IC=1,2
      DO 206 IR=1,JRI
       if(isi.eq.isj)then
       vde=1.d0            
       else
       vde=0.d0                           
       endif
      AMECCM(IR,IC)=Rlo1(IR,L,isi)*Rlo1(IR,L,isj)*vde
      AMEACM(IR,IC)=RR1(IR,L,isi)*Rlo1(IR,L,isj)*vde
      AMEBCM(IR,IC)=RE1(IR,L,isi)*Rlo1(IR,L,isj)*vde
  206 CONTINUE
      CALL CALI(AMECCM(1,IC),TOTCC(IC),RX(1),DX,JRI)
      CALL CALI(AMEACM(1,IC),TOTAC(IC),RX(1),DX,JRI)
      CALL CALI(AMEBCM(1,IC),TOTBC(IC),RX(1),DX,JRI)
  204 CONTINUE
!
      RCC(L,isi,isj)=COE*(TOTCC(1)+CIN*TOTCC(2))
      RAC(L,isi,isj)=COE*(TOTAC(1)+CIN*TOTAC(2))
      RBC(L,isi,isj)=COE*(TOTBC(1)+CIN*TOTBC(2))
      else
      RCC(L,isi,isj)=0.0d0
      RAC(L,isi,isj)=0.0d0
      RBC(L,isi,isj)=0.0d0
      endif
  200 CONTINUE
91    continue
      RETURN
      END

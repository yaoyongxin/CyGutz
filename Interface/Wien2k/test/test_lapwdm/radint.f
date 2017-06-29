SUBROUTINE RADINT(JATOM,ISPIN,vr,e)
USE param
USE struct
        IMPLICIT REAL*8 (A-H,O-Z)
      real*8           :: VR(NRAD)
      real*8  e(0:lmax2)
      COMMON /RADFU/   RF1(NRAD,0:LMAX2,2,nrf),RF2(NRAD,0:LMAX2,2,nrf)
      COMMON /RINTEG/  RI_MAT(0:lmax2,nrf,nrf,2,2)
      DIMENSION  rx(nrad),a11(nrad),a12(nrad),a21(nrad),a22(nrad)
      common /aver/ krad,kls,cx(-20:20,20),iprx
!.....set up radial mesh.
      l=1
      DO  I=1,JRJ(JATOM)
      RX(I)=R0(JATOM)*EXP(DX(JATOM)*(i-1))
      enddo
    
      ri_mat(0:lmax2,1:nrf,1:nrf,1:2,1:2)=0.d0

      DO 19 L=0,LMAX2
       E(L)=E(L)/2.d0
      DO 20 IS1=1,ISPIN
      DO 20 IS2=1,ISPIN
      DO 20 IF1=1,NRF
      DO 20 IF2=1,IF1
          DO  I=1,JRJ(JATOM)
          if((krad.eq.0).or.(krad.eq.1))then
          A11(i)=rf1(i,l,is1,if1)
          A21(i)=rf2(i,l,is1,if1)
          A12(i)=rf1(i,l,is2,if2)
          A22(i)=rf2(i,l,is2,if2)
          else if((krad.eq.2).or.(krad.eq.3))then
            dnom=1.d0+(E(L)-vr(i))/(2.d0*clight*clight)
            A11(i)=rf1(i,l,is1,if1)/(dnom*rx(i)**1.5)
            A21(i)=0.d0
            A12(i)=rf1(i,l,is2,if2)/(dnom*rx(i)**1.5)
            A22(i)=0.d0
          else if(krad.eq.4)then  ! WIEN2k_5 approx.
           A11(i)=rf1(i,l,is1,if1)/rx(i)**1.5
           A21(i)=0.d0
           A12(i)=rf1(i,l,is2,if2)/rx(i)**1.5
           A22(i)=0.d0
          else if(krad.gt.10)then  ! <r**krad>
           rkrad=dfloat(krad-10)/2.d0
           A11(i)=rf1(i,l,is1,if1)*rx(i)**rkrad
           A21(i)=rf2(i,l,is1,if1)*rx(i)**rkrad
           A12(i)=rf1(i,l,is2,if2)*rx(i)**rkrad
           A22(i)=rf2(i,l,is2,if2)*rx(i)**rkrad
          else if(krad.lt.-10)then  ! <1/r**krad>
           rkrad=-dfloat(10+krad)/2.d0
           A11(i)=rf1(i,l,is1,if1)/rx(i)**rkrad
           A21(i)=rf2(i,l,is1,if1)/rx(i)**rkrad
           A12(i)=rf1(i,l,is2,if2)/rx(i)**rkrad
           A22(i)=rf2(i,l,is2,if2)/rx(i)**rkrad
          endif
          enddo
      CALL RINT13(A11,A21,A12,A22, &
           ri_mat(l,if1,if2,is1,is2),JATOM)
      ri_mat(l,if2,if1,is2,is1)=ri_mat(l,if1,if2,is1,is2)
 20   CONTINUE
      E(L)=2.d0*E(L)
 19   CONTINUE
      END


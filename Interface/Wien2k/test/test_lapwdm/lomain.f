SUBROUTINE LoMAIN(nemin,nemax,lfirst,latom,n,jatom,is)
  USE param
  USE struct
  USE xxa
  USE xa
  USE xa3
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16       YL((LMAX2+1)*(LMAX2+1))
      COMPLEX*16       PHSHEL,IMAG                     
!                                                                       
      COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
      logical         loor(0:lomax),lapw(0:lomax)                                        
      common /loabc/  alo(0:lomax,2,nloat,nrf)
      common /lolog/  nlo,nlov,nlon,loor,ilo(0:lomax),lapw  
     
      DATA            CZERO/(0.0D0,0.0D0)/,IMAG/(0.0D0,1.0D0)/         
                                                                       
!------------------------------------------------------------------     
      PI=ACOS(-1.0D0)                                                   
      TWOPI=2.D0*PI
	mu=latom-lfirst+1
      i=n-(nlo+nlon)                                                      
      DO 10 L=0,LoMAX
       do 10 jlo=1,ilo(l)
        do 20 jneq=1,mult(jatom)
          DO 25 M1=-l,+l                                                    
          i=i+1    
        BK(1)=BKX(I)*BR1(1,1)+BKY(I)*BR1(1,2)+BKZ(I)*BR1(1,3)
        BK(2)=BKX(I)*BR1(2,1)+BKY(I)*BR1(2,2)+BKZ(I)*BR1(2,3)
        BK(3)=BKX(I)*BR1(3,1)+BKY(I)*BR1(3,2)+BKZ(I)*BR1(3,3)
          CALL YLM (BK,LOMAX,YL)   
        ARG1=BKX(I)*POS(1,LATOM)*TWOPI
        ARG2=BKY(I)*POS(2,LATOM)*TWOPI
        ARG3=BKZ(I)*POS(3,LATOM)*TWOPI
!         PHSHEL=EXP(IMAG*(ARG1+ARG2+ARG3))
          PHSHEL=DCMPLX(DCOS(ARG1+ARG2+ARG3),DSIN(ARG1+ARG2+ARG3))
          DO 50 NUM=NEMIN,NEMAX                                           
            PHS(NUM)=PHSHEL*A(I,NUM)         
  50      CONTINUE                                                        
          DO 30 M=-l,+l                                                    
            index=l*(l+1)+m+1 
            DO 40 NUM=NEMIN,NEMAX     
              DO 40 irf=1,nrf  
      ALM(index,num,mu,irf,is)=ALM(INDEX,num,mu,irf,is)+ &
                      ALo(l,is,jlo,irf)* &
                        dconjg(YL(INDEX))*PHS(NUM)
  40        CONTINUE                              
  30      CONTINUE
  25    CONTINUE                                                          
  20    CONTINUE                                                          
  10  CONTINUE   
      return                   
      END        

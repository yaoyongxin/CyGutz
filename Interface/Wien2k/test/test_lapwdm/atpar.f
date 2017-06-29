SUBROUTINE ATPAR (JATOM,LATOM,itape,jtape,is,ISPIN,notcalc)
USE param
USE struct
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8           :: VR(NRAD)
      LOGICAL          loor(0:lomax),lapw(0:lmax2),rlo(1:nloat, &
                       0:lomax),notcalc                                       
      DIMENSION     emist(0:lomax,nloat),E(0:LMAX2),elo(0:LOMAX,nloat), &
                       pei(0:lmax2),ilo(0:lomax) 
                                                                       
      COMMON /UHELP/   A(NRAD),B(NRAD),AE(NRAD),BE(NRAD)
      COMMON /ATSPDT/  P(0:LMAX2,2,nrf),DP(0:LMAX2,2,nrf)
      COMMON /RADFU/   RF1(NRAD,0:LMAX2,2,nrf),RF2(NRAD,0:LMAX2,2,nrf)
      COMMON /RINTEG/  RI_MAT(0:lmax2,nrf,nrf,2,2)
      common /loabc/   alo(0:lomax,2,nloat,nrf)
      common /lolog/   nlo,nlov,nlon,loor,ilo,lapw  
!---------------------------------------------------------------------  
 2022 FORMAT(3X,4E19.12) 
                                                                       
!.....READ TOTAL SPHERICAL POTENTIAL V(0,0) OF TAPEjtape=VSP               
!     NORM OF V(0,0)=V00(R)*R/SQRT(4.D0*PI)     
                        
      READ(jtape,1980)                                                     
      READ(jtape,2000)IDUMMY                                              
      READ(jtape,2031)                                                     
      READ(jtape,2022)(VR(J),J=1,JRJ(JATOM))                           
      READ(jtape,2031)                                                     
      READ(jtape,2030)
  
      DO 115 J=1,JRJ(JATOM)                                             
 115  VR(J)=VR(J)/2.0D0       
       if(.not.notcalc)then
       if(is.eq.1)then
	write(6,*)'ATPAR for spin up'  
       else
	write(6,*)'ATPAR for spin down'  
       endif
       endif
      
      nlo=0
      nlov=0
      nlon=0
      do i=0,lomax
      ilo(i)=0
      end do
      DO 15 I=1,JATOM                                                   
        READ(itape)E 
!       write(6,*)E
        READ(itape)elo
        do jlo=1,nloat
!       write(6,*)(elo(l,jlo),l=0,lomax)
        end do
        IF(i.EQ.jatom) THEN
           DO l=0,lmax2
              lapw(l)=.TRUE.
              IF(e(l).GT.150.) THEN
                 e(l)=e(l)-200.d+0
                 lapw(l)=.FALSE.
              ENDIF
           ENDDO
        ENDIF
        do 15 l = 0,lomax
              loor(l)=.FALSE.
           DO k=1,nloat
              rlo(k,l)=.FALSE.
              IF (i.EQ.jatom) THEN
                 IF (elo(l,k).LT.(995.d+0)) THEN
                    ilo(l)=ilo(l)+1
                    nlo=nlo+((2*l+1))*mult(i)
                    IF(.NOT.lapw(l).AND.k.EQ.1) GOTO 666
                    IF(k.EQ.nloat) THEN
                          rlo(ilo(l),l)=.TRUE.
                          GOTO 666
                    ENDIF
                    loor(l)=.TRUE.
 666                CONTINUE
                 ENDIF
              ELSE
                 IF (elo(l,k).LT.(995.d+0)) nlov=nlov+((2*l+1))*mult(i)
              ENDIF
           ENDDO
  15    continue
      IF(JATOM.EQ.NAT) GOTO 30                                          
      DO 20 I=JATOM+1,NAT  
        READ(itape) EMIST                                                    
        READ(itape) EMIST                                                    
       do 20 l=0,lomax
	do 20 k=1,nloat
          if (emist(l,k).lt.(995.0D+0))  &
              nlon=nlon+((2*l+1))*mult(i)
  20  continue
  30  CONTINUE     
      if(notcalc)return  
      WRITE(6,7) ANAME(JATOM)
      WRITE(6,5) E
      WRITE(6,14)

!                                                                       
      DO 70 l=0,LMAX2                                                  
      DELE=2.0D-3                                                       
      DELEI=0.25D0/DELE                                                 
      FL=L                                                              
      EI=E(l)/2.0d0        
!     CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE                  
!     DELE IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN HARTREES          
!                                                                       
      E1=EI-DELE  
      CALL OUTWIN(REL,VR,R0(JATOM),DX(JATOM),JRJ(JATOM),E1,            &
      FL,UVB,DUVB,NODEL,ZZ)                                                
      CALL RINT13(A,B,A,B,OVLP,JATOM)                               
      TRX=1.0D0/SQRT(OVLP)                                              
      IMAX=JRJ(JATOM)                                                   
      DO 45 M=1,IMAX               
      AE(M)=TRX*A(M)                                                    
      BE(M)=TRX*B(M)                                                    
   45 CONTINUE                                                          
      UVB=TRX*UVB                                                       
      DUVB=TRX*DUVB                                                     
      E1=EI+DELE                                                        
      CALL OUTWIN(REL,VR,R0(JATOM),DX(JATOM),JRJ(JATOM),E1,            &
      FL,UVE,DUVE,NODE,ZZ)                                                 
      CALL RINT13(A,B,A,B,OVLP,JATOM)                               
      TRX=1.0d0/SQRT(OVLP)                                                
      UVE=DELEI*(TRX*UVE-UVB)                                           
      DUVE=DELEI*(TRX*DUVE-DUVB)                                        
      IMAX=JRJ(JATOM)                                                   
      DO 50 M=1,IMAX                                                    
      AE(M)=DELEI*(TRX*A(M)-AE(M))                                      
      BE(M)=DELEI*(TRX*B(M)-BE(M))                                      
   50 CONTINUE                                                          
!                                                                       
!     CALCULATE FUNCTION AT EI                                          
!                                                                       
      CALL OUTWIN(REL,VR(1),R0(JATOM),DX(JATOM),JRJ(JATOM),EI,         &
      FL,UV,DUV,NODES,ZZ)                                                  
      CALL RINT13(A,B,A,B,OVLP,JATOM)                               
      TRX=1.0d0/SQRT(OVLP)                                                
      P(l,is,1)=TRX*UV                                                       
      DP(l,is,1)=TRX*DUV                                                     
      IMAX=JRJ(JATOM)                                                   
      DO 60 M=1,IMAX                                                    
      A(M)=TRX*A(M)                                                     
   60 B(M)=TRX*B(M)                                                     
!                                                                       
!     INSURE ORTHOGONALIZATION                                          
!                                                                       
      CALL RINT13(A,B,AE,BE,CROSS,JATOM)                            
      TRY=-CROSS                                                        
      IMAX=JRJ(JATOM)                                                   
      DO 55 M=1,IMAX                                                    
      AE(M)=(AE(M)+TRY*A(M))                                            
   55 BE(M)=(BE(M)+TRY*B(M))                                            
      IMAX=JRJ(JATOM)                                                   
      DO 80 I=1,IMAX                                                    
          RF1(I,l,is,1)=A(I)                                                   
          RF2(I,l,is,1)=B(I)                                                   
          RF1(I,l,is,2)=AE(I)                                                  
          RF2(I,l,is,2)=BE(I)                                                  
  80    CONTINUE                                                          
        P(l,is,2)=UVE+TRY*P(l,is,1)                                                
        DP(l,is,2)=DUVE+TRY*DP(l,is,1) 
        CALL RINT13(AE,BE,AE,BE,PEI(l),JATOM)                            
        WRITE(6,8) L,P(l,is,1),DP(l,is,1),P(l,is,2),DP(l,is,2)
  70    continue                                               
!                         
! nun fur lo
!
      DO 170 l=0,lomax
        irf=2
	do 170 jlo=1,ilo(l)
        if ((.not.lapw(l)).and.(jlo.eq.1)) goto 180
	irf=irf+1
        DELE=2.0D-3                                                       
        DELEI=0.25D0/DELE                                                 
        FL=L                                                              
        EI=elo(l,jlo)/2.d0                                                         
!                                                                       
!     CALCULATE FUNCTION AT EI 

         IF(rlo(jlo,l)) THEN
            ei=elo(l,nloat)/2.d0
            kappa=l
            CALL diracout(rel,vr(1),r0(jatom),dx(jatom),jrj(jatom),    &
                     ei,fl,kappa,uv,duv,nodes,zz)
            CALL dergl(a,b,r0(jatom),dx(jatom),jrj(jatom))
            DO m = 1, jrj(jatom)
            r_m = r0(jatom)*exp(dx(jatom)*(m-1))
            b(m) = b(m)*r_m/(2.d0*clight+(elo(l,jlo)- &
                   2.d0*vr(m)/r_m)/(2.d0*clight))
            b(m)=b(m)*clight
            ENDDO
         ELSE
            CALL outwin(rel,vr(1),r0(jatom),dx(jatom),jrj(jatom),   &
                    ei,fl,uv,duv,nodes,zz)
         ENDIF                                         
!                                                                       
        CALL RINT13(A,B,A,B,OVLP,JATOM)                               
        TRX=1.0d0/SQRT(OVLP)   
        P(l,is,irf)=TRX*UV  
        DP(l,is,irf)=TRX*DUV 
        IMAX=JRJ(JATOM)                                                   
        DO 160 M=1,IMAX                                                    
        rf1(M,l,is,irf)=TRX*A(M)                                                     
  160   rf2(M,l,is,irf)=TRX*B(M) 
                                                    
        CALL RINT13(rf1(1,l,is,1),rf2(1,l,is,1), &
                    rf1(1,l,is,irf),rf2(1,l,is,irf),pi12lo,JATOM)
        CALL RINT13(rf1(1,l,is,2),rf2(1,l,is,2), &
                    rf1(1,l,is,irf),rf2(1,l,is,irf),pe12lo,JATOM)
  180 continue         
      call abc (l,jatom,pei(l),pi12lo,pe12lo,is,jlo,lapw(l))
      WRITE(6,*)
 170  continue
!... CALCULATION OF RADIAL INTEGRALS
      IF (ISPIN.EQ.2.AND.IS.EQ.1) GOTO 469
      CALL RADINT(JATOM,ISPIN,vr,e)
 469  RETURN                                                            
                                                                       
    2 FORMAT(5E14.7)                                                    
    3 FORMAT(16I5)                                                      
    4 FORMAT(I4,E16.7)                                                  
    5 FORMAT(10X,' ENERGY PARAMETERS ARE',7F7.2)                        
    6 FORMAT(10X,'E(',I2,2H)=,F10.4)                                    
    7 FORMAT(/10X,'ATOMIC PARAMETERS FOR ',A10/)                        
   14 FORMAT(/11X,1HL,5X,4HU(R),10X,                                     &
             5HU'(R),9X,5HDU/DE,8X,6HDU'/DE,6X,7HNORM-U')               
    8 FORMAT(10X,I2,5E14.6,5X,3I2)                                      
    9 FORMAT(10X,'FOR L=',I2,' CORRECTION=',E14.5,' OVERLAP=',E14.5)    
   11 FORMAT(7F10.5)                                                    
   13 FORMAT(////,':POS',i3.3,':',1x,'AT.NR.',I3,2X,'POSITION =', &
             3F8.5,2X,'MULTIPLICITY =',I3)
 1040 FORMAT(8I2)                                                       
 1050 FORMAT(A10,I5,5X,2F10.9,I5,F5.2)                                  
 1060 FORMAT(//,3X,'NOT EQUIV ATOM ',A10,'  LOCAL ROTATION MATRIX')     
 1070 FORMAT(30X,3F10.5)                                                
 1080 FORMAT(13X,'EQUIV ATOM ',I3,3X,'POSITION: ',3F8.3)                
 1090 FORMAT(30X,3F10.5,5X,F10.5)                                       
 1980 FORMAT(3X)                                                        
 2000 FORMAT(15X,I3//)                                                  
 2030 FORMAT(///)                                                       
 2031 FORMAT(/)                                                         
      END                                                               

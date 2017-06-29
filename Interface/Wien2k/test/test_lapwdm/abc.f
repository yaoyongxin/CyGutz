      subroutine abc(l,jatom,pei,pi12lo,pe12lo,is,jlo,lapw)                           
        USE param
        USE struct
!     abc calculates the cofficients a,b,c of the lo    
      IMPLICIT REAL*8 (A-H,O-Z)
      common /loabc/   alo(0:lomax,2,nloat,nrf)
      COMMON /ATSPDT/  P(0:LMAX2,2,nrf),DP(0:LMAX2,2,nrf)
      logical lapw     
!---------------------------------------------------------------------  
!      
      do irf=1,nrf
	alo(l,is,jlo,irf)=0.d0
      end do                             
      if (lapw) then
      irf=2+jlo                                    
      xac=p(l,is,irf)*dp(l,is,2)-dp(l,is,irf)*p(l,is,2)
      xac=xac*rmt(jatom)*rmt(jatom)
      xbc=p(l,is,irf)*dp(l,is,1)-dp(l,is,irf)*p(l,is,1)
      xbc=-xbc*rmt(jatom)*rmt(jatom)
      clo=xac*(xac+2.0D0*pi12lo)+xbc* &
                (xbc*pei+2.0D0*pe12lo)+1.0D0  
      clo=1.0D0/sqrt(clo)
      write(6,*)clo
      if (clo.gt.2.0D2) clo=2.0d2
      alo(l,is,jlo,1)=clo*xac
      alo(l,is,jlo,2)=clo*xbc
      alo(l,is,jlo,irf)=clo
      write (6,10) l,alo(l,is,jlo,1),alo(l,is,jlo,2), &
                  alo(l,is,jlo,irf)
      else
      if (jlo.eq.1) then
	alonorm=sqrt(1.d0+(P(l,is,1)/P(l,is,2))**2 &
                 *PEI)
      alo(l,is,jlo,1)=1.d0/alonorm
      alo(l,is,jlo,2)=-P(l,is,1)/P(l,is,2)/alonorm
      else
      xbc=-P(l,is,1)/P(l,is,1+jlo)
      xac=sqrt(1+xbc**2+2*xbc*PI12LO)
      alo(l,is,jlo,1)=1.d0/xac
      alo(l,is,jlo,1+jlo)=xbc/xac
      end if
      write (6,10)l,(alo(l,is,jlo,jrf),jrf=1,nrf)
      end if	
 10   FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,6F12.5)
      return
      end



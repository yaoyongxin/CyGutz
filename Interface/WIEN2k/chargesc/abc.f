SUBROUTINE abc(l,jlo,jatom)
  !                                                                       
  ! abc calculates the cofficients a,b,c of the lo    
  ! Ouput: alo(l,jlo), blo(l,jlo), clo(l,jlo)
  !                                                                       
  USE param, ONLY: nrad, fastFilesystem
  USE defs, ONLY: 
  USE struk, ONLY: rmt
  USE lo, ONLY: lapw, plo, dplo, alo, blo, clo, pi12lo, pe12lo
  USE atspdt, ONLY: e=>el,p,dp,pe,dpe,pei
  USE com_mpi
  !USE parallel
  IMPLICIT NONE
  INTEGER, intent(in) :: l,jatom,jlo
  ! locals
  REAL*8              :: xac,xbc,alonorm

  !---------------------------------------------------------------------  
  !                                                                       
  if(lapw(l)) then
     ! DPE are indexed from 0 to LMAX2 in lapw2 usw
     ! DPE are indexed from 1 to LMAX2+1 in lapw1 usw
     xac=plo(jlo,l)*dpe(l)-dplo(jlo,l)*pe(l)
     xac=xac*rmt(jatom)*rmt(jatom)
     xbc=plo(jlo,l)*dp(l)-dplo(jlo,l)*p(l)
     xbc=-xbc*rmt(jatom)*rmt(jatom)
     clo(l,jlo)=xac*(xac+2.0D0*pi12lo(jlo,l))+xbc*(xbc*pei(l)+2.0D0*pe12lo(jlo,l))+1.0D0  
     clo(l,jlo)=1.0D0/sqrt(clo(l,jlo))
     if (clo(l,jlo).gt.2.0D2) clo(l,jlo)=2.0d2
     alo(l,jlo)=clo(l,jlo)*xac
     blo(l,jlo)=clo(l,jlo)*xbc
  else
     if(jlo.eq.1) then
        alonorm=sqrt(1.d0+(P(l)/PE(l))**2*PEI(l))
        ALO(l,jlo) = 1.d0 /alonorm 
        BLO(l,jlo) = -P(l)/PE(l)/alonorm
        CLO(l,jlo) = 0.d0
     else 
        xbc=-P(l)/PLO(jlo,l)
        xac=sqrt(1+xbc**2+2*xbc*PI12LO(jlo,l))
        !         xac=1.d0
        ALO(l,jlo) = 1.d0/xac 
        BLO(l,jlo) = 0.d0
        CLO(l,jlo) = xbc/xac
     endif
  endif
  if (myrank.EQ.master .OR. fastFilesystem/=0) WRITE(6,10) l,alo(l,jlo),blo(l,jlo),clo(l,jlo)
10 FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,3F12.5)
  return
END SUBROUTINE abc

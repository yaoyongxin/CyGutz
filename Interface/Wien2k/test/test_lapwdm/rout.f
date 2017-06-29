SUBROUTINE ROUT(LL)
USE param
	IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /RINTEG/  RI_MAT(0:lmax2,nrf,nrf,2,2)

	write(6,*)'LL=',ll
 3      FORMAT(8F9.5)
	WRITE(6,*)'RADIAL INTEGRALS:'
	DO 10 IRF1=1,NRF
        DO 10 IS1=1,2
	DO 10 IS2=1,2
	DO 10 IRF2=1,NRF
        if (abs(RI_MAT(LL,IRF1,IRF2,IS1,IS2)).gt.1d-4) &
      	WRITE(6,*)irf1,irf2,'     ',is1,is2,RI_MAT(LL,IRF1,IRF2,IS1,IS2)
 10     continue 
	end

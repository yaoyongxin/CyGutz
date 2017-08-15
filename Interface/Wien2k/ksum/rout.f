SUBROUTINE ROUT(L1,L2,iso)
  USE param
  IMPLICIT REAL*8 (A-H,O-Z)

  COMMON /RINTEG/  RI_MAT(0:lmax2,0:lmax2,nrf,nrf,2,2)

  write(6,*)'L1=',L1,' L2=',L2 
  WRITE(6,*)'RADIAL INTEGRALS:'
  DO IRF1=1,NRF      ! 10
     DO IS1=1,iso    ! 10
	DO IS2=1,iso ! 10
          DO IRF2=1,NRF ! 10
             if (abs(RI_MAT(L1,L2,IRF1,IRF2,IS1,IS2)).gt.1d-4) WRITE(6,*)irf1,irf2,'     ',is1,is2,RI_MAT(L1,L2,IRF1,IRF2,IS1,IS2)
          ENDDO !10     continue 
       ENDDO
    ENDDO
 ENDDO
END SUBROUTINE ROUT

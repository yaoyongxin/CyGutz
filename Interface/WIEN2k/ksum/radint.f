SUBROUTINE RADINT(JATOM,ISPIN,QoffDiag)
  USE param
  USE com_mpi
  USE struct
  USE case
  IMPLICIT REAL*8 (A-H,O-Z)
  logical  QoffDiag  ! off-diagonal l needed
  COMMON /RADFU/   RF1(NRAD,0:LMAX2,2,nrf),RF2(NRAD,0:LMAX2,2,nrf)
  COMMON /RINTEG/  RI_MAT(0:lmax2,0:lmax2,nrf,nrf,2,2)
  DIMENSION  rx(nrad)
  ! rf1 and rf2 are two parts of the solution of the Dirac equation : large and small component
  ! rf1 has the following indexes:
  !    rf1(r,l,s,irf)  where r is radial distance, s is spin,
  !                    and irf is one of the radial functions 
  !    for LAPW      we have u(E1),dot{u}(E1)
  !    for LAPW+LO   we have u(E1),dot{u}(E1), u(E2)
  !    for APW+lo    we have u(E1),dot{u}(E1)
  !    for APW+lo+LO we have u(E1),dot{u}(E1), u(E2)
  ri_mat(0:lmax2,0:lmax2,1:nrf,1:nrf,1:ispin,1:ispin)=0.d0
  if (QoffDiag) then
     do  L1=0,LMAX2
        do  L2=0,LMAX2
           do  IS1=1,ISPIN
              do  IS2=1,ISPIN
                 do  IF1=1,NRF
                    do  IF2=1,NRF
                       ! It integrates 
                       ! ri_mat(l1,l2,irf1,irf2,is1,is2) = integrate( rf1(1,l1,is1,irf1)*rf1(1,l2,is2,irf2) + c*rf2(1,l1,is1,irf1)*rf2(1,l2,is2,irf2) )
                       CALL RINT13(rf1(1,l1,is1,if1),rf2(1,l1,is1,if1), rf1(1,l2,is2,if2),rf2(1,l2,is2,if2), ri_mat(l1,l2,if1,if2,is1,is2), JATOM)
                    enddo
                 enddo
              enddo
           enddo
           if (myrank.EQ.master .OR. fastFilesystem) call rout(l1,l2,ispin)
        enddo
     enddo
  else
     do l=0,3
        do  IS1=1,ISPIN
           do  IS2=1,ISPIN
              do  IF1=1,NRF
                 do  IF2=1,NRF
                    CALL RINT13(rf1(1,l,is1,if1),rf2(1,l,is1,if1), rf1(1,l,is2,if2),rf2(1,l,is2,if2), ri_mat(l,l,if1,if2,is1,is2),JATOM)
                 enddo
              enddo
           enddo
        enddo
        if (myrank.EQ.master .OR. fastFilesystem) call rout(l,l,ispin)
     enddo
  endif
END SUBROUTINE RADINT


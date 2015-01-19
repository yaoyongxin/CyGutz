SUBROUTINE TestDoubleTransformation(DMFTrans, nbands, nind, maxdim2)
  ! It tests the following identity:
  ! sum_{ij} DMFTrans(i,j,L1,L2) * DMFTrans(j,i,L2p,L1p) = delta(L1,L1p)*delta(L2,L2p)
  ! which makes shure that a quantity G_{LL'} when transformed to G_{ij} and back to G_{LL'} takes the same form
  IMPLICIT NONE
  COMPLEX*16, intent(in) :: DMFTrans(nbands,nbands,maxdim2,maxdim2)
  INTEGER, intent(in)    :: nbands, nind, maxdim2
  ! local variables
  COMPLEX*16 :: csum
  REAL*8     :: small
  INTEGER    :: ind1, ind2, ind1p, ind2p, i, j
  small = 1e-3
  do ind1=1,nind
     do ind2=1,nind
        do ind1p=1,nind
           do ind2p=1,nind
              csum=0
              do i=1,nbands
                 do j=1,nbands
                    csum = csum + DMFTrans(i,j,ind1,ind2) * DMFTrans(j,i,ind2p,ind1p)
                 enddo
              enddo
              if (abs(csum)>small) WRITE(*,'(A1,2x,4I3,2x,2f10.4)') 'w', ind1, ind1p, ind2, ind2p, csum
           enddo
        enddo
     enddo
  enddo
END SUBROUTINE TestDoubleTransformation

SUBROUTINE TestOrtgogonality(DMFTrans, nbands, nind, maxdim2)
  ! It tests the following idenity 
  ! \sum_i DMFTrans(i,i,L1,L2) = delta(L1,L2)
  IMPLICIT NONE
  COMPLEX*16, intent(in) :: DMFTrans(nbands,nbands,maxdim2,maxdim2)
  INTEGER, intent(in)    :: nbands, nind, maxdim2
  ! local variables
  COMPLEX*16 :: csum
  REAL*8     :: small
  INTEGER    :: ind1, ind2, i
  small = 1e-3
  do ind1=1,nind
     do ind2=1,nind
        csum=0
        do i=1,nbands
           csum = csum + DMFTrans(i,i,ind1,ind2)
        enddo
        if (abs(csum)>1e-3) WRITE(*,'(A,2x,2I4,2x,2f10.5)') 't', ind1, ind2, csum
     enddo
  enddo
END SUBROUTINE TestOrtgogonality

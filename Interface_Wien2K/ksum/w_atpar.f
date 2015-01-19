module w_atpar
  IMPLICIT NONE
  !*******************************************
  ! atpar changes the following common blocks:
  REAL*8, allocatable :: w_alo(:,:,:,:,:)
  INTEGER, allocatable:: w_nlo(:), w_nlov(:), w_nlon(:), w_ilo(:,:)
  LOGICAL, allocatable:: w_loor(:,:), w_lapw(:,:)
  REAL*8, allocatable :: w_ri_mat(:,:,:,:,:,:,:)
  !REAL*8, allocatable :: w_rf1(:,:,:,:,:), w_rf2(:,:,:,:,:)
  REAL*8, allocatable :: w_P(:,:,:,:), w_DP(:,:,:,:)
  !INTEGER,allocatable :: w_lfirst(:)
  INTEGER, allocatable:: w_jatom(:)
  REAL*8, allocatable :: w_FJ(:,:,:)
  REAL*8, allocatable :: w_DFJ(:,:,:)
  CONTAINS

  SUBROUTINE w_allocate(maxucase)
    use param
    use case
    IMPLICIT NONE
    integer :: maxucase
    allocate( w_alo(0:lomax,2,nloat,nrf,maxucase) )
    allocate( w_nlo(maxucase), w_nlov(maxucase), w_nlon(maxucase), w_ilo(0:lomax,maxucase) )
    allocate( w_loor(0:lomax,maxucase), w_lapw(0:lmax2,maxucase) )
    allocate( w_ri_mat(0:lmax2,0:lmax2,nrf,nrf,2,2,maxucase) )
    !allocate( w_rf1(NRAD,0:LMAX2,2,nrf,maxucase), w_rf2(NRAD,0:LMAX2,2,nrf,maxucase) ) ! It seems to me that this one does not need to be saved!
    allocate( w_P(0:LMAX2,2,nrf,maxucase), w_DP(0:LMAX2,2,nrf,maxucase) )              ! Needs to be saved!
    !allocate( w_lfirst(maxucase) )
    allocate( w_jatom(maxucase))    
    allocate( w_FJ(0:lmax2,nmat,maxucase) )
    allocate( w_DFJ(0:lmax2,nmat,maxucase) )
    w_alo=0; w_nlo=0; w_nlov=0; w_nlon=0; w_ilo=0; w_loor=0; w_lapw=0; w_ri_mat=0; w_P=0; w_DP=0; 
    w_jatom=0; w_FJ=0; w_DFJ=0
  END SUBROUTINE w_allocate
  
  SUBROUTINE w_deallocate
    IMPLICIT NONE
    deallocate( w_alo )
    deallocate( w_nlo, w_nlov, w_nlon, w_ilo ) 
    deallocate( w_loor, w_lapw )
    deallocate( w_ri_mat )
    !deallocate( w_rf1, w_rf2 )
    deallocate( w_P, w_DP )
    !deallocate( w_lfirst )
    deallocate( w_jatom )
    deallocate( w_FJ )
    deallocate( w_DFJ )
  END SUBROUTINE w_deallocate
  
end module w_atpar

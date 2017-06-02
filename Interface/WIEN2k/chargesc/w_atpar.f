module w_atpar
  IMPLICIT NONE
  INTEGER, ALLOCATABLE :: w_lfirst(:)
  INTEGER, ALLOCATABLE :: w_nlo(:)
  INTEGER, ALLOCATABLE :: w_nlov(:)
  INTEGER, ALLOCATABLE :: w_nlon(:)
  INTEGER, ALLOCATABLE :: w_ilo(:,:)
  LOGICAL, ALLOCATABLE :: w_loor(:,:,:)
  LOGICAL, ALLOCATABLE :: w_lapw(:,:)
  REAL*8 , ALLOCATABLE :: w_alo(:,:,:)
  REAL*8 , ALLOCATABLE :: w_blo(:,:,:)
  REAL*8 , ALLOCATABLE :: w_clo(:,:,:)

  REAL*8 , ALLOCATABLE :: w_P(:,:)
  REAL*8 , ALLOCATABLE :: w_DP(:,:)
  REAL*8 , ALLOCATABLE :: w_PE(:,:)
  REAL*8 , ALLOCATABLE :: w_DPE(:,:)
  REAL*8 , ALLOCATABLE :: w_PEI(:,:)
  REAL*8 , ALLOCATABLE :: w_pi12lo(:,:,:)
  REAL*8 , ALLOCATABLE :: w_pe12lo(:,:,:)
  REAL*8 , ALLOCATABLE :: w_pr12lo(:,:,:,:)

  INTEGER, ALLOCATABLE :: w_lm(:,:,:)
  INTEGER, ALLOCATABLE :: w_lmmax(:)
  REAL*8 , ALLOCATABLE :: w_R(:,:)

  REAL*8 , ALLOCATABLE :: w_xwt1(:,:)
  REAL*8 , ALLOCATABLE :: w_xwt1l(:,:)
  REAL*8 , ALLOCATABLE :: w_xwt1h(:,:)
  REAL*8 , ALLOCATABLE :: w_xwteh(:,:)
  REAL*8 , ALLOCATABLE :: w_xwtel(:,:)

  REAL*8 , ALLOCATABLE :: w_a1lo(:,:,:,:)
  REAL*8 , ALLOCATABLE :: w_b1lo(:,:,:,:)
  REAL*8 , ALLOCATABLE :: w_RRAD1(:,:,:)
  REAL*8 , ALLOCATABLE :: w_RRAD2(:,:,:)
  REAL*8 , ALLOCATABLE :: w_RADE1(:,:,:)
  REAL*8 , ALLOCATABLE :: w_RADE2(:,:,:)

  REAL*8 , ALLOCATABLE :: w_vrholm(:,:,:)
  REAL*8 , ALLOCATABLE :: w_rholm(:,:,:)
  
  CONTAINS

  SUBROUTINE w_allocate0(nat)
    USE param
    IMPLICIT NONE
    INTEGER, intent(in) :: nat
    ! 
    ALLOCATE( w_lm(1:2,1:NCOM,1:nat) )
    ALLOCATE( w_lmmax(1:nat) )
  END SUBROUTINE w_allocate0
  
  SUBROUTINE w_allocate(LM_MAX, nat)
    USE param
    IMPLICIT NONE
    INTEGER, intent(in) :: LM_MAX, nat
    
    ALLOCATE( w_lfirst(1:nat) )
    ALLOCATE( w_nlo(1:nat) )
    ALLOCATE( w_nlov(1:nat) )
    ALLOCATE( w_nlon(1:nat) )
    ALLOCATE( w_ilo(0:lmax2,1:nat) )
    ALLOCATE( w_loor(1:nloat,0:lomax,1:nat) )
    ALLOCATE( w_lapw(0:lmax2,1:nat) )
    ALLOCATE( w_alo(0:lomax,1:nloat,1:nat) )
    ALLOCATE( w_blo(0:lomax,1:nloat,1:nat) )
    ALLOCATE( w_clo(0:lomax,1:nloat,1:nat) )
    
    ALLOCATE( w_P(0:LMAX2,1:nat) )
    ALLOCATE( w_DP(0:LMAX2,1:nat) )
    ALLOCATE( w_PE(0:LMAX2,1:nat) )
    ALLOCATE( w_DPE(0:LMAX2,1:nat) )
    ALLOCATE( w_PEI(0:LMAX2,1:nat) )
    ALLOCATE( w_pi12lo(1:nloat,0:lomax,1:nat) )
    ALLOCATE( w_pe12lo(1:nloat,0:lomax,1:nat) )
    ALLOCATE( w_pr12lo(1:nloat,1:nloat,0:lomax,1:nat) )
    
    ALLOCATE( w_R(1:nrad,1:nat) )
    
    ALLOCATE( w_xwt1(0:21,1:nat) )
    ALLOCATE( w_xwt1l(0:3,1:nat) )
    ALLOCATE( w_xwt1h(0:3,1:nat) )
    ALLOCATE( w_xwteh(0:3,1:nat) )
    ALLOCATE( w_xwtel(0:3,1:nat) )
    
    ALLOCATE( w_a1lo(1:nrad,1:nloat,0:lomax,1:nat) )
    ALLOCATE( w_b1lo(1:nrad,1:nloat,0:lomax,1:nat) )
    ALLOCATE( w_RRAD1(1:nrad,0:lmax2,1:nat) )
    ALLOCATE( w_RRAD2(1:nrad,0:lmax2,1:nat) )
    ALLOCATE( w_RADE1(1:nrad,0:lmax2,1:nat) )
    ALLOCATE( w_RADE2(1:nrad,0:lmax2,1:nat) )
    
    ALLOCATE( w_vrholm(1:NRAD,1:LM_MAX,1:nat) )
    ALLOCATE( w_rholm(1:NRAD,1:LM_MAX,1:nat) )
    
  END SUBROUTINE w_allocate
  
  SUBROUTINE w_deallocate
    IMPLICIT NONE
    DEALLOCATE( w_lfirst )
    DEALLOCATE( w_nlo )
    DEALLOCATE( w_nlov )
    DEALLOCATE( w_nlon )
    DEALLOCATE( w_ilo )
    DEALLOCATE( w_loor )
    DEALLOCATE( w_lapw )
    DEALLOCATE( w_alo )
    DEALLOCATE( w_blo )
    DEALLOCATE( w_clo )

    DEALLOCATE( w_P )
    DEALLOCATE( w_DP )
    DEALLOCATE( w_PE )
    DEALLOCATE( w_DPE )
    DEALLOCATE( w_PEI )
    DEALLOCATE( w_pi12lo )
    DEALLOCATE( w_pe12lo )
    DEALLOCATE( w_pr12lo )

    DEALLOCATE( w_lm )
    DEALLOCATE( w_lmmax )
    DEALLOCATE( w_R )

    DEALLOCATE( w_xwt1 )
    DEALLOCATE( w_xwt1l )
    DEALLOCATE( w_xwt1h )
    DEALLOCATE( w_xwteh )
    DEALLOCATE( w_xwtel )

    DEALLOCATE( w_a1lo )
    DEALLOCATE( w_b1lo )
    DEALLOCATE( w_RRAD1 )
    DEALLOCATE( w_RRAD2 )
    DEALLOCATE( w_RADE1 )
    DEALLOCATE( w_RADE2 )

    DEALLOCATE( w_vrholm )
    DEALLOCATE( w_rholm )
  END SUBROUTINE w_deallocate
end module w_atpar

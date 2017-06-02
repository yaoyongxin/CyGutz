      MODULE GHDF5
      USE GPREC; USE GHDF5_BASE; use sparse; use hdf5
      IMPLICIT NONE
      private
      public :: gh5_write_compound
!
      INTERFACE gh5_write_compound
        MODULE PROCEDURE gh5_write_dcsr_matrix, gh5_write_zcsr_matrix, &
            & gh5_write_zbd_csr_matrix, gh5_write_dcoo_zcsr_matrix, &
            & gh5_write_dbd_csr_matrix
      END INTERFACE
!
      contains
!
      subroutine gh5_write_dcsr_matrix(dcsr, path, f_id)
      type(dcsr_matrix), target, intent(in) :: dcsr
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer nnz
!
      call gh5_write(dcsr%nrow, path//".nrow", f_id)
      call gh5_write(dcsr%ncol, path//".ncol", f_id)
      call gh5_write(1, path//".base", f_id)
      call gh5_write(dcsr%i, dcsr%nrow + 1, path//".indptr", f_id)
      nnz = dcsr%I(dcsr%nrow + 1) -1
      call gh5_write(dcsr%j, nnz, path//".indices", f_id)
      call gh5_write(dcsr%a, nnz, path//".data", f_id)
      RETURN
!
      end subroutine gh5_write_dcsr_matrix
!
      subroutine gh5_write_zcsr_matrix(zcsr, path, f_id)
      type(zcsr_matrix), target, intent(in) :: zcsr
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer nnz
!
      call gh5_write(zcsr%nrow, path//".nrow", f_id)
      call gh5_write(zcsr%ncol, path//".ncol", f_id)
      call gh5_write(1, path//".base", f_id)
      call gh5_write(zcsr%i, zcsr%nrow + 1, path//".indptr", f_id)
      nnz = zcsr%I(zcsr%nrow + 1) -1
      call gh5_write(zcsr%j, nnz, path//".indices", f_id)
      call gh5_write(zcsr%a, nnz, path//".data", f_id)
      RETURN
!
      end subroutine gh5_write_zcsr_matrix
!
      subroutine gh5_write_zbd_csr_matrix(zbd, path, f_id)
      type(zbd_matrix) zbd
      character(len = *) path
      INTEGER(HID_T), intent(in) :: f_id
! local
      type(zcsr_matrix) zcsr
!
      CALL ZBMTOZCSR(zbd, zcsr)
      CALL gh5_write_zcsr_matrix(zcsr, path, f_id)
      CALL DEALLOC_ZCSR(zcsr)
      return
!
      end subroutine gh5_write_zbd_csr_matrix
!
      subroutine gh5_write_dbd_csr_matrix(dbd, path, f_id)
      type(dbd_matrix) dbd
      character(len = *) path
      INTEGER(HID_T), intent(in) :: f_id
! local
      type(dcsr_matrix) dcsr
!
      CALL DBMTODCSR(dbd, dcsr)
      CALL gh5_write_dcsr_matrix(dcsr, path, f_id)
      CALL DEALLOC_DCSR(dcsr)
      return
!
      end subroutine gh5_write_dbd_csr_matrix
!
      subroutine gh5_write_dcoo_zcsr_matrix(dcoo, path, f_id)
      type(dcoo_matrix) dcoo
      character(len = *) path
      INTEGER(HID_T), intent(in) :: f_id
! local
      type(zcsr_matrix) zcsr
!
      CALL SM_DCOOZCSR(dcoo, zcsr)
      CALL gh5_write_zcsr_matrix(zcsr, path, f_id)
      CALL DEALLOC_ZCSR(zcsr)
      return
!
      end subroutine gh5_write_dcoo_zcsr_matrix
!
!
      END MODULE GHDF5

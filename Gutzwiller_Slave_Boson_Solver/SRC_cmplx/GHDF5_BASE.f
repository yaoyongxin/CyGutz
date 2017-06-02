      MODULE GHDF5_BASE
      use hdf5; USE GPREC
      USE GUTIL, ONLY: int_to_str
      IMPLICIT NONE
      private 
      public :: gh5_init, gh5_end, gh5_write, gh5_create_groups, &
          & gh5_open, gh5_close
      INTEGER(HID_T),public :: log_file_id, embedH_file_id
      INTEGER(HID_T),private :: group_id, dset_id, dspace_id
      INTEGER(HID_T),private :: dcomplex_id, string_id
      INTEGER(HSIZE_T),private :: gh5_dims(7), gh5_size, gh5_offset
      TYPE(C_PTR),private :: f_ptr
      INTEGER,private :: gh5_err
!
      INTERFACE gh5_write
        MODULE PROCEDURE gh5_write_iscalar, gh5_write_dscalar, gh5_write_sarray, &
        & gh5_write_1d_iarray, gh5_write_2d_iarray, &
        & gh5_write_1d_darray, gh5_write_2d_darray, gh5_write_4d_darray, &
        & gh5_write_1d_zarray, gh5_write_2d_zarray, gh5_write_4d_zarray, &
        & gh5_write_3d_darray, gh5_write_3d_zarray, gh5_write_5d_darray, &
        & gh5_write_5d_zarray, gh5_write_1d_sarray
      END INTERFACE
      contains
!
      subroutine gh5_init(gh5_file)
      character(*) gh5_file
!
      ! Initialize FORTRAN interface.
      call h5open_f(gh5_err)
      ! Create a new file using default properties.
      CALL h5fcreate_f(gh5_file, H5F_ACC_TRUNC_F, log_file_id, gh5_err)
      CALL set_dcomplex_id()
      CALL H5Tcopy_f(H5T_FORTRAN_S1, string_id, gh5_err)
      return
!
      end subroutine gh5_init
!     
      subroutine gh5_open(gh5_file, f_id)
      character(*),intent(in) :: gh5_file
      INTEGER(HID_T),intent(out) :: f_id
!
      ! Create a new file using default properties.
      CALL h5fcreate_f(gh5_file, H5F_ACC_TRUNC_F, f_id, gh5_err)
      return
!      
      end subroutine gh5_open
!     
      subroutine gh5_close(f_id)
      INTEGER(HID_T),intent(in) :: f_id
!
      CALL h5fclose_f(f_id, gh5_err)
      return
!      
      end subroutine gh5_close
!     
      subroutine gh5_end()
!
      CALL H5Tclose_f(dcomplex_id, gh5_err)
      CALL H5Tclose_f(string_id, gh5_err)
      ! Close the file.
      call h5fclose_f(log_file_id, gh5_err)
      ! Close FORTRAN interface.
      call h5close_f(gh5_err)
      return
!
      end subroutine gh5_end
!     
      subroutine set_dcomplex_id()
!
      gh5_size = gq*2; gh5_offset = 0
      CALL H5Tcreate_f(H5T_COMPOUND_F, gh5_size, dcomplex_id, gh5_err)
      CALL H5Tinsert_f(dcomplex_id, "r", gh5_offset, h5kind_to_type(gq, H5_REAL_KIND), gh5_err)
      gh5_offset = gh5_offset + gq
      CALL H5Tinsert_f(dcomplex_id, "i", gh5_offset, h5kind_to_type(gq, H5_REAL_KIND), gh5_err)
      RETURN
!
      end subroutine set_dcomplex_id
!
      subroutine gh5_write_iscalar(i, path, f_id)
      integer, intent(in) :: i
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! LOCAL
      integer ibuf(1), N(1)
      
      ibuf = i; N = 1
      call gh5_write_iarray(ibuf, 1, N, path, f_id)
      RETURN
!
      end subroutine gh5_write_iscalar
!
      subroutine gh5_write_dscalar(d, path, f_id)
      real(gq), intent(in) :: d
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! LOCAL
      real(gq) dbuf(1)
      integer N(1)
      
      dbuf = d; N = 1
      call gh5_write_darray(dbuf, 1, N, path, f_id)
      RETURN
!
      end subroutine gh5_write_dscalar
!
      subroutine gh5_write_1d_iarray(iarray, N1, path, f_id)
      integer, intent(in) :: N1, iarray(N1)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(1)
!
      N = N1
      call gh5_write_iarray(iarray, 1, N, path, f_id)
      return
!
      end subroutine gh5_write_1d_iarray
!
      subroutine gh5_write_2d_iarray(iarray, N1, N2, path, f_id)
      integer, intent(in) :: N1, N2, iarray(N1, N2)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(2)
!
      N(1) = N1; N(2) = N2
      call gh5_write_iarray(iarray, 2, N, path, f_id)
      return
!
      end subroutine gh5_write_2d_iarray
!
      subroutine gh5_write_1d_darray(darray, N1, path, f_id)
      integer, intent(in) :: N1
      real(gq), intent(in) :: darray(N1)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(1)
!
      N = N1
      call gh5_write_darray(darray, 1, N, path, f_id)
      return
!
      end subroutine gh5_write_1d_darray
!
      subroutine gh5_write_2d_darray(darray, N1, N2, path, f_id)
      integer, intent(in) :: N1, N2
      real(gq), intent(in) :: darray(N1, N2)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(2)
!
      N(1) = N1; N(2) = N2
      call gh5_write_darray(darray, 2, N, path, f_id)
      return
!
      end subroutine gh5_write_2d_darray
!
      subroutine gh5_write_3d_darray(darray, N1, N2, N3, path, f_id)
      integer, intent(in) :: N1, N2, N3
      real(gq), intent(in) :: darray(N1, N2, N3)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(3)
!
      N(1) = N1; N(2) = N2; N(3) = N3
      call gh5_write_darray(darray, 3, N, path, f_id)
      return
!
      end subroutine gh5_write_3d_darray
!
      subroutine gh5_write_4d_darray(darray, N1, N2, N3, N4, path, f_id)
      integer, intent(in) :: N1, N2, N3, N4
      real(gq), intent(in) :: darray(N1, N2, N3, N4)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(4)
!
      N(1) = N1; N(2) = N2; N(3) = N3; N(4) = N4
      call gh5_write_darray(darray, 4, N, path, f_id)
      return
!
      end subroutine gh5_write_4d_darray
!
      subroutine gh5_write_5d_darray(darray, N1, N2, N3, N4, N5, path, f_id)
      integer, intent(in) :: N1, N2, N3, N4, N5
      real(gq), intent(in) :: darray(N1, N2, N3, N4, N5)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(5)
!
      N(1) = N1; N(2) = N2; N(3) = N3; N(4) = N4; N(5) = N5
      call gh5_write_darray(darray, 5, N, path, f_id)
      return
!
      end subroutine gh5_write_5d_darray
!
      subroutine gh5_write_1d_zarray(zarray, N1, path, f_id)
      integer, intent(in) :: N1
      complex(gq), intent(in) :: zarray(N1)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(1)
!
      N = N1
      call gh5_write_zarray(zarray, 1, N, path, f_id)
      return
!
      end subroutine gh5_write_1d_zarray
!
      subroutine gh5_write_1d_sarray(sarray, NS, N1, path, f_id)
      integer, intent(in) :: NS,N1
      character(NS), intent(in) :: sarray(N1)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(1)
!
      N = N1
      call gh5_write_sarray(sarray, NS, 1, N, path, f_id)
      return
!
      end subroutine gh5_write_1d_sarray
!
      subroutine gh5_write_2d_zarray(zarray, N1, N2, path, f_id)
      integer, intent(in) :: N1, N2
      complex(gq), intent(in) :: zarray(N1, N2)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(2)
!
      N(1) = N1; N(2) = N2
      call gh5_write_zarray(zarray, 2, N, path, f_id)
      return
!
      end subroutine gh5_write_2d_zarray
!
      subroutine gh5_write_3d_zarray(zarray, N1, N2, N3, path, f_id)
      integer, intent(in) :: N1, N2, N3
      complex(gq), intent(in) :: zarray(N1, N2, N3)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(3)
!
      N(1) = N1; N(2) = N2; N(3) = N3
      call gh5_write_zarray(zarray, 3, N, path, f_id)
      return
!
      end subroutine gh5_write_3d_zarray
!
      subroutine gh5_write_4d_zarray(zarray, N1, N2, N3, N4, path, f_id)
      integer, intent(in) :: N1, N2, N3, N4
      complex(gq), intent(in) :: zarray(N1, N2, N3, N4)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(4)
!
      N(1) = N1; N(2) = N2; N(3) = N3; N(4) = N4
      call gh5_write_zarray(zarray, 4, N, path, f_id)
      return
!
      end subroutine gh5_write_4d_zarray
!
      subroutine gh5_write_5d_zarray(zarray, N1, N2, N3, N4, N5, path, f_id)
      integer, intent(in) :: N1, N2, N3, N4, N5
      complex(gq), intent(in) :: zarray(N1, N2, N3, N4, N5)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer N(5)
!
      N(1) = N1; N(2) = N2; N(3) = N3; N(4) = N4; N(5) = N5
      call gh5_write_zarray(zarray, 5, N, path, f_id)
      return
!
      end subroutine gh5_write_5d_zarray
!
      subroutine gh5_write_iarray(iarray, ND, N, path, f_id)
      integer, intent(in) :: ND, N(ND), iarray(*)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
!
      gh5_dims(1:ND) = N
      ! Create the dataspace.
      CALL h5screate_simple_f(ND, gh5_dims(1:ND), dspace_id, gh5_err)
      ! Create the dataset with default properties.
      CALL h5dcreate_f(f_id, path, H5T_NATIVE_INTEGER, dspace_id, dset_id, gh5_err)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, iarray, gh5_dims(1:ND), gh5_err)
      ! End access to the dataset and release resources used by it.
      CALL h5dclose_f(dset_id, gh5_err)
      ! Terminate access to the data space.
      CALL h5sclose_f(dspace_id, gh5_err)
      RETURN
!
      end subroutine gh5_write_iarray
!
      subroutine gh5_write_darray(darray, ND, N, path, f_id)
      integer, intent(in) :: ND, N(ND)
      real(gq), intent(in) :: darray(*)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
!
      gh5_dims(1:ND) = N
      ! Create the dataspace.
      CALL h5screate_simple_f(ND, gh5_dims(1:ND), dspace_id, gh5_err)
      ! Create the dataset with default properties.
      CALL h5dcreate_f(f_id, path, H5T_NATIVE_DOUBLE, dspace_id, dset_id, gh5_err)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, darray, gh5_dims(1:ND), gh5_err)
      ! End access to the dataset and release resources used by it.
      CALL h5dclose_f(dset_id, gh5_err)
      ! Terminate access to the data space.
      CALL h5sclose_f(dspace_id, gh5_err)
      RETURN
!
      end subroutine gh5_write_darray
!
      subroutine gh5_write_zarray(zarray, ND, N, path, f_id)
      integer, intent(in) :: ND, N(ND)
      complex(gq), target, intent(in) :: zarray(*)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
!
      gh5_dims(1:ND) = N
      ! Create the dataspace.
      CALL h5screate_simple_f(ND, gh5_dims(1:ND), dspace_id, gh5_err)
      ! Create the dataset with default properties.
      CALL h5dcreate_f(f_id, path, dcomplex_id, dspace_id, dset_id, gh5_err)
      f_ptr = C_LOC(zarray)
      CALL h5dwrite_f(dset_id, dcomplex_id, f_ptr, gh5_err)
      ! End access to the dataset and release resources used by it.
      CALL h5dclose_f(dset_id, gh5_err)
      ! Terminate access to the data space.
      CALL h5sclose_f(dspace_id, gh5_err)
      RETURN
!
      end subroutine gh5_write_zarray
!
      subroutine gh5_write_sarray(sarray, NS, ND, N, path, f_id)
      integer, intent(in) :: NS, ND, N(ND)
      character(len=NS), target, intent(in) :: sarray(*)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
!
      gh5_size = NS
      CALL H5Tset_size_f(string_id, gh5_size, gh5_err)
      gh5_dims(1:ND) = N
      ! Create the dataspace.
      CALL h5screate_simple_f(ND, gh5_dims(1:ND), dspace_id, gh5_err)
      ! Create the dataset with default properties.
      CALL h5dcreate_f(f_id, path, string_id, dspace_id, dset_id, gh5_err)
      f_ptr = C_LOC(sarray(1)(1:1))
      CALL h5dwrite_f(dset_id, string_id, f_ptr, gh5_err)
      ! End access to the dataset and release resources used by it.
      CALL h5dclose_f(dset_id, gh5_err)
      ! Terminate access to the data space.
      CALL h5sclose_f(dspace_id, gh5_err)
      RETURN
!
      end subroutine gh5_write_sarray
!
      subroutine gh5_create_groups(iarray, N, path, f_id)
      integer, intent(in) :: N, iarray(N)
      character(len=*), intent(in) :: path
      INTEGER(HID_T), intent(in) :: f_id
! local
      integer i
!
      do i = 1, N
        CALL h5gcreate_f(f_id, path//trim(int_to_str(iarray(i))), group_id, gh5_err)
        CALL h5gclose_f(group_id, gh5_err)
      enddo
      RETURN
!
      end subroutine gh5_create_groups
!
      END MODULE GHDF5_BASE

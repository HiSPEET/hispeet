!> summary:  Write an SEM variable to HDF5
!> author:   Joerg Stiller
!> date:     2024/06/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(SEM__Variable__3D) MP_WriteHDF5
  use, intrinsic ::  ISO_C_Binding
  use HDF5_Binding
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Write SEM variable to HDF5 file

  module subroutine WriteHDF5_F(this, file) result(this)
    class(SEM_Variable_3D), intent(in) :: this
    character(len=*),       intent(in) :: file !< name of HDF5 file

    integer(hid_t) :: file_id, group_id
    integer :: err

    call Init_HDF5_Binding()

    call H5Fcreate_f(trim(file), H5F_ACC_TRUNC_F, file_id, err)
    call H5Gcreate_f(file_id, 'sev', group_id, err)

    call WriteHDF5_G(this, group_id)

    call H5Gclose_f(group_id, err)
    call H5Fclose_f(file_id, err)

  end subroutine WriteHDF5_F

  !-----------------------------------------------------------------------------
  !> Write SEM variable into given HDF5 group

  module subroutine WriteHDF5_G(this, group_id)
    class(SEM_Variable_3D), intent(in) :: this     !< SEM variable
    integer(HID_T),         intent(in) :: group_id !< ID of related HDF5 group

    ! HDF5 datatype, dataspace, dimensions and buffer adress pointer
    integer(hid_t)   :: space_id
    integer(hid_t)   :: data_id
    integer(hsize_t) :: dims(5), maxdims(5)
    type(C_Ptr)      :: buf

    integer :: err

    call Init_HDF5_Binding()

    if (.not. associated(this%val)) then
      call Error('WriteHDF5_G','val not associated','SEM__Variable__3D')
    end if

    dims = shape(this%val)
    buf  = C_Loc(this%val(0,0,0,1,1))

    call H5Screate_simple_f(1, dims, space_id, err)
    call H5Dcreate_f(group_id, 'val', H5T_REAL_RNP, space_id, data_id, err)
    call H5Dwrite_f(data_id, H5T_REAL_RNP, buf, err)
    call H5Sclose_f(space_id, err)
    call H5Dclose_f(data_id, err)

  end subroutine WriteHDF5_G

  !=============================================================================

end submodule MP_WriteHDF5

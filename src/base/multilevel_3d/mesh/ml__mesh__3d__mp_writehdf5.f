!> summary:  Writing a multilevel mesh partition to HDF5
!> author:   Joerg Stiller
!> date:     2024/10/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(ML__Mesh__3D) MP_WriteHDF5
  use, intrinsic :: ISO_C_Binding
  use HDF5_Binding
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Write multilevel mesh partition into HDF5 file

  module subroutine WriteHDF5(this, file)
    class(ML_Mesh_3D), intent(in) :: this !< multilevel mesh partition
    character(len=*),  intent(in) :: file !< name of HDF5 file

    character(len=:), allocatable :: file_pr
    character(len=80) :: tag
    integer(HID_T)    :: file_id, group_id
    integer, target   :: l_top
    integer           :: err, l

    ! preliminaries ............................................................

    ! safeguard
    call Init_HDF5_Binding()

    ! append process rank to file name
    write(tag,'(I0)') this % mesh(1) % proc
    file_pr = trim(file)//'_'//trim(tag)//'.h5'

    l_top = size(this % mesh)

    ! HDF5 file and base group .................................................

    call H5Fcreate_f(file_pr, H5F_ACC_TRUNC_F, file_id, err)
    call H5Gcreate_f(file_id, '/ml_mesh', group_id, err)

    ! attributes ...............................................................

    block
      integer(hid_t)   :: data_id, space_id
      integer(hsize_t) :: dims(1)

      call H5Gcreate_f(file_id, '/ml_mesh/attrib', group_id, err)

      ! l_top
      dims = 1
      call H5Screate_simple_f(size(dims), dims, space_id, err)
      call H5Dcreate_f(group_id, 'l_top', H5T_INTEGER, space_id, data_id, err)
      call H5Dwrite_f(data_id, H5T_INTEGER, C_Loc(l_top), err)
      call H5Sclose_f(space_id, err)
      call H5Dclose_f(data_id, err)

    end block

    ! mesh levels ..............................................................

    do l = 1, l_top
      write(tag,'(I0)') l
      call H5Gcreate_f(file_id, '/ml_mesh/level_'//trim(tag), group_id, err)
      call this % mesh(l) % WriteHDF5(group_id)
      call H5Gclose_f(group_id, err)
    end do

    ! finalize .................................................................

    call H5Fclose_f(file_id, err)

  end subroutine WriteHDF5

  !=============================================================================

end submodule MP_WriteHDF5

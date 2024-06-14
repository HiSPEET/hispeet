!> summary:  Creating an SEM variable from HDF5 data
!> author:   Joerg Stiller, Erik Pfister
!> date:     2024/06/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_New_HDF5
  use, intrinsic ::  ISO_C_Binding
  use HDF5_Binding
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> New SEM variable from spectral-element mesh and data given in HDF5 file
  !>
  !> An empty variable with zero bounds is created if `sem` has no elements.

  module function New_HDF5_F(sem, group)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    character(len=*), intent(in) :: file !< name of HDF5 file
    type(SEM_Variable_3D) :: this

    ! TBD

  end subroutine New_HDF5_F

  !-----------------------------------------------------------------------------
  !> New SEM variable from spectral-element mesh and data given in HDF5 group
  !>
  !> Creates a new spectral-element variable associated with `sem` with data
  !> provided by the HDF5 `group`.

  module function New_HDF5_G(sem, group)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    integer(hid_t), intent(in) :: group !< ID of HDF5 group
    type(SEM_Variable_3D) :: this

    ! HDF5 datatype, dataspace, dimensions and buffer adress pointer
    integer(hid_t)   :: type_id
    integer(hid_t)   :: space_id
    integer(hid_t)   :: data_id
    integer(hsize_t) :: dims(5), maxdims(5)
    type(C_Ptr)      :: buf

    integer :: err
    integer :: nc

    call Init_HDF5_Binding()

    if (group == H5I_INVALID_HID_F) then
      call Error('New_HDF5_G','invalid group given', 'SEM__Variable__3D')
    end if

    if (allocated(this % mem)) then
      deallocate(this % mem)
    end if

    call H5Dopen_f(group, 'mem', data_id, err)
    call H5Dget_space_f(data_id, space_id, err)
    call H5Sget_simple_extent_dims_f(space_id, dims, maxdims, err)

    associate(po => sem%std_op%po, ne => sem%mesh%n_elem)

      if (dims(1) == po+1 .and. dims(4) == ne) then

        nc = dims(5)
        allocate(this % mem(0:po, 0:po, 0:po, ne, nc))

        buf = C_Loc(this % mem)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Tclose_f(type_id, err)
        call H5Dclose_f(data_id, err)

        this % sem => sem
        this % val(0:,0:,0:,1:,1:) => this % mem

      else

        call Error('New_HDF5_G','HDF5 data not matching','SEM__Variable__3D')

      end if

    end associate

  end subroutine New_HDF5_G

  !=============================================================================

end submodule MP_New_HDF5

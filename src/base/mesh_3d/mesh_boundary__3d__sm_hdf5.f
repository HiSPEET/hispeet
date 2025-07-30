!> summary:  ...
!> author:   Joerg Stiller, Erik Pfister, Moritz Kreuseler
!> date:     2024/06/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_Boundary__3D) SM_HDF5
  use, intrinsic ::  ISO_C_Binding
  use HDF5_Binding
  implicit none

  !-----------------------------------------------------------------------------
  !> HDF5 datatype for mesh boundary attributes

  integer(HID_T) :: H5T_BoundaryAttributes = -1

contains

  !-----------------------------------------------------------------------------
  !> Get HDF5 datatype for mesh boundary attributes

  module subroutine Get_H5T_MeshBoundaryAttributes_3D &
      ( H5T_MeshBoundaryAttributes_3D )

    integer(HID_T), intent(out) :: H5T_MeshBoundaryAttributes_3D

    !$omp master

    if (H5T_BoundaryAttributes < 0) then
      call Init_H5T_BoundaryAttributes()
    end if

    H5T_MeshBoundaryAttributes_3D = H5T_BoundaryAttributes

    !$omp end master

  end subroutine Get_H5T_MeshBoundaryAttributes_3D

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_BoundaryAttributes

  subroutine Init_H5T_BoundaryAttributes()

    type(MeshBoundaryAttributes_3D), target :: attributes(2)
    integer(SIZE_T)  :: offset
    integer(HID_T)   :: type_id
    integer :: err

    ! initialize datatype
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_BoundaryAttributes, err)

    ! insert name
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%name(1:1)))
    call H5Tcopy_f(H5T_CHARACTER, type_id, err)
    call H5Tset_size_f(type_id, int(len(attributes(1)%name), SIZE_T), err)
    call H5Tinsert_f(H5T_BoundaryAttributes, 'name', offset, type_id, err)
    call H5Tclose_f(type_id, err)

    ! insert id
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%id))
    call H5Tinsert_f(H5T_BoundaryAttributes, 'id', offset, H5T_INTEGER, err)

    ! insert coupled
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%coupled))
    call H5Tinsert_f(H5T_BoundaryAttributes, 'coupled', offset, H5T_INTEGER, err)

    ! insert polarity
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%polarity))
    call H5Tinsert_f(H5T_BoundaryAttributes, 'polarity', offset, H5T_INTEGER, err)

    ! insert map
    offset = H5offsetof(C_Loc(attributes(1)), C_Loc(attributes(1)%map(1,1)))
    call H5Tarray_create_f(H5T_REAL_RNP, 2, int([4,4], HSIZE_T), type_id, err)
    call H5Tinsert_f(H5T_BoundaryAttributes, 'map', offset, type_id, err)
    call H5Tclose_f(type_id, err)

  end subroutine Init_H5T_BoundaryAttributes

  !=============================================================================

end submodule SM_HDF5

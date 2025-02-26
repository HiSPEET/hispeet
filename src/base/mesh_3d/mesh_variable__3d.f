!> summary:  3D mesh variable
!> author:   Joerg Stiller
!> date:     2024/09/27
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_Variable__3D
  use Kind_Parameters
  use Mesh__3D
  implicit none
  private

  public :: MeshVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D mesh variable

  type MeshVariable_3D

    integer :: po = 0 !< polynomial order
    integer :: nc = 0 !< number of components

    class(Mesh_3D), pointer :: mesh => null()
    real(RNP), contiguous, pointer :: val(:,:,:,:,:) => null() !< value access
    real(RNP), allocatable :: mem(:,:,:,:,:) !< value storage

  contains

    procedure :: Init => Init_MeshVariable_3D
    procedure :: GetSlice

  end type MeshVariable_3D

contains

  !-----------------------------------------------------------------------------
  !> 3D mesh variable initialization

  subroutine Init_MeshVariable_3D(this, mesh, po, nc)
    class(MeshVariable_3D), target, intent(inout) :: this
    class(Mesh_3D), target, intent(in) :: mesh
    integer, intent(in) :: po
    integer, intent(in) :: nc

    allocate(this % mem(0:po,0:po,0:po,mesh%n_elem,nc))

    this % po = po
    this % nc = nc
    this % mesh => mesh
    this % val => this % mem

  end subroutine Init_MeshVariable_3D

  !-----------------------------------------------------------------------------
  !> Create a new mesh variable as a slice of the given one
  !>
  !> The values of the new variable refer to `this % val` if `copy` is false
  !> or absent. Otherwise they are stored in fresh memory, i.e. `slice % val`.

  subroutine GetSlice(this, slice, first, last, copy)
    class(MeshVariable_3D), target, intent(in) :: this
    class(MeshVariable_3D), target, intent(inout) :: slice
    integer,           intent(in) :: first !< first component of slice
    integer,           intent(in) :: last  !< last component of slice
    logical, optional, intent(in) :: copy  !< copy into fresh memory [F]

    integer :: nc
    logical :: copy_

    associate(po => this%po, ne => this%mesh%n_elem)

      nc = 1 + last - first

      slice % po = po
      slice % nc = nc
      slice % mesh => this % mesh

      if (present(copy)) then
        copy_ = copy
      else
        copy_ = .false.
      end if

      if (copy_) then
        allocate( slice % mem(0:po,0:po,0:po,ne,nc)       &
                , source = this % val(:,:,:,:,first:last) )
        slice % val => slice % mem
      else
        slice % val(0:,0:,0:,1:,1:) => this % val(:,:,:,:,first:last)
      end if

    end associate

  end subroutine GetSlice

  !=============================================================================

end module Mesh_Variable__3D

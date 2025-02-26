!> summary:  3D single-level spacetime mesh variable
!> author:   Joerg Stiller
!> date:     2024/11/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spacetime_Variable__3D
  use Mesh__3D
  use Mesh_Variable__3D
  implicit none
  private

  public :: SpacetimeVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D single-level spacetime mesh variable
  !>
  !> The values of a spacetime variable `u` are accessed via
  !>
  !>       u % var(m,n) % val(i,j,k,e,c)
  !>
  !>  where
  !>
  !>    - `0 ≤   m   ≤ po_t`:  temporal collocation point ID
  !>    - `1 ≤   n   ≤ ne_t`:  temporal element ID in given time slice
  !>    - `0 ≤ i,j,k ≤ po_x`:  spatial collocation point triple index
  !>    - `1 ≤   e   ≤ ne_x`:  spatial element ID

  type SpacetimeVariable_3D
    integer :: po_t = 0 !< polynomial order in time
    integer :: ne_t = 0 !< number of time elements/steps
    integer :: nc       !< number of components
    type(MeshVariable_3D), allocatable :: var(:,:)
  contains
    procedure :: Init => Init_SpacetimeVariable_3D
    procedure :: GetSlice
  end type SpacetimeVariable_3D

contains

  !-----------------------------------------------------------------------------
  !> Initialization of 3D spacetime variable

  subroutine Init_SpacetimeVariable_3D(this, mesh, po_x, po_t, ne_t, nc)
    class(SpacetimeVariable_3D), intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh !< mesh partition
    integer,        intent(in) :: po_x !< polynomial order in space
    integer,        intent(in) :: po_t !< polynomial order in time
    integer,        intent(in) :: ne_t !< number of time elements/steps
    integer,        intent(in) :: nc   !< number of components

    integer :: m, n

    allocate(this % var(0:po_t,1:ne_t))

    this % po_t = po_t
    this % ne_t = ne_t
    this % nc   = nc

    do n = 1, ne_t
    do m = 0, po_t
      call this % var(m,n) % Init(mesh, po_x, nc)
    end do
    end do

  end subroutine Init_SpacetimeVariable_3D

  !-----------------------------------------------------------------------------
  !> Create a new spacetime variable as a slice of the given one
  !>
  !> The values of the new variable refer to `this % var(m,n) % val` if `copy`
  !> is false or absent. Otherwise they are stored in fresh memory.

  subroutine GetSlice(this, slice, first, last, copy)
    class(SpacetimeVariable_3D), intent(in)    :: this
    class(SpacetimeVariable_3D), intent(inout) :: slice
    integer,           intent(in) :: first !< first component of slice
    integer,           intent(in) :: last  !< last component of slice
    logical, optional, intent(in) :: copy  !< copy into fresh memory [F]

    integer :: m, n

    if (allocated(slice % var)) then
      deallocate(slice % var)
    end if

    allocate(slice % var(0:this%po_t,1:this%ne_t))

    slice % po_t = this % po_t
    slice % ne_t = this % ne_t
    slice % nc   = last - first + 1

    do n = 1, this % ne_t
    do m = 0, this % po_t
      call this % var(m,n) % GetSlice(slice%var(m,n), first, last, copy)
    end do
    end do

  end subroutine GetSlice

  !=============================================================================

end module Spacetime_Variable__3D

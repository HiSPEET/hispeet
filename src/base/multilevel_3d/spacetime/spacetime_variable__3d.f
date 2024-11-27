!> summary:  3D single-level spacetime mesh variable
!> author:   Joerg Stiller
!> date:     2024/11/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spacetime_Variable__3D
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

    type(MeshVariable_3D), allocatable :: var(:,:)

  end type SpacetimeVariable_3D

contains

  !=============================================================================

end module Spacetime_Variable__3D

!> summary:  3D multilevel spacetime mesh variable
!> author:   Joerg Stiller
!> date:     2024/11/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Spacetime_Variable__3D
  use Spacetime_Variable__3D
  implicit none
  private

  public :: ML_SpacetimeVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel spacetime mesh variable
  !>
  !> The values of a multilevel spacetime variable `u` are accessed via
  !>
  !>       u % level(l) % var(m,n) % val(i,j,k,e,c)
  !>
  !>  where
  !>
  !>    - `1 ≤   l   ≤ l_top`:  level
  !>    - `0 ≤   m   ≤ pt(l)`:  temporal collocation point ID
  !>    - `1 ≤   n   ≤ nt(l)`:  temporal element ID in given time slice
  !>    - `0 ≤ i,j,k ≤ po(l)`:  spatial collocation point triple index
  !>    - `1 ≤   e   ≤ ne(l)`:  spatial element ID

  type ML_SpacetimeVariable_3D
    type(SpacetimeVariable_3D), allocatable :: level(:) !< variable per level
  end type ML_SpacetimeVariable_3D

contains

  !-----------------------------------------------------------------------------
  !> Constructor of 3D multilevel spacetime mesh variable



  !-----------------------------------------------------------------------------
  !> Initialization of 3D multilevel spacetime mesh variable


  !=============================================================================

end module ML__Spacetime_Variable__3D

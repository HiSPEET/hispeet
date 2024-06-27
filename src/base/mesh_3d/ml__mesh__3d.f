!> summary:  3D multilevel mesh
!> author:   Joerg Stiller
!> date:     2024/06/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Mesh__3D
  use Mesh__3D
  implicit none
  private

  public :: ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel mesh
  !>
  !> The `mesh` component contains a sequence of meshes the from level `1` up to
  !> `l_top = size(mesh)`. The relation between mesh levels is described by the
  !> `refinement` pattern. It ranges from `1` to `l_top-1` and can take the
  !> following values
  !>
  !>   | `refinement` | child level generation | application  |
  !>   |:------------:| ---------------------- | ------------ |
  !>   |     'c'      | global cloning         | p-refinement |
  !>   |     'z'      | zonal cloning          | p-adaptivity |
  !>   |     'r'      | regular refinement     | h-refinement |
  !>   |     'l'      | local refinement       | h-adaptivity |
  !>

  type ML_Mesh_3D
    type(Mesh_3D), allocatable :: mesh(:)       !< mesh partitions
    character,     allocatable :: refinement(:) !< refinement pattern
  end type ML_Mesh_3D

contains

  !-----------------------------------------------------------------------------
  !> New multilevel mesh by global refinement of given mesh



  !=============================================================================

end module ML__Mesh__3D

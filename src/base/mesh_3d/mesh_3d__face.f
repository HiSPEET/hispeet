!> summary:  3D mesh face type
!> author:   Joerg Stiller
!> date:     2020/11/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_3d__Face
  implicit none
  private

  public :: Mesh3d_Face

  !-----------------------------------------------------------------------------
  !> Adjacent mesh element data

  type AdjacentElement
    integer :: id   = 0 !< element id
    integer :: face = 0 !< corresponding element face {1:6}
  end type AdjacentElement

  !-----------------------------------------------------------------------------
  !> 3D mesh face
  !>
  !> The main purpose is to access the adjacent elements. The latter are sorted
  !> according to the normal direction of the face. In case of a structured mesh
  !> the latter coincides with the corresponding coordinate direction.

  type Mesh3d_Face
    type(AdjacentElement) :: element(2) !< adjacent elements, including ghosts
  end type Mesh3d_Face

end module Mesh_3d__Face

!> summary:  Cartesian mesh face
!> author:   Joerg Stiller
!> date:     2017/04/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Cartesian mesh face
!===============================================================================

module CART__Mesh_Face
  implicit none
  private

  public :: MeshFace

  !-----------------------------------------------------------------------------
  !> Cartesian mesh face
  !>
  !> The elements are sorted following to the coordinate direction. For example,
  !> if the normal vector is aligned with the x-direction, element(1) precedes
  !> element(2) on the x-axis. Since the orientation of the mesh faces is known
  !> (see MeshPartition), this convention allows to identify the corresponding
  !> element face: In case of an x-face, element(1) abuts with face 2 and
  !> element(2) with face 1. In case of y- and z-faces these are faces 4 and 3,
  !> 6 and 5, respectively.

  type MeshFace
    integer :: element(2) = 0 !< adjacent element ID, including ghost, 0 if none
  end type MeshFace

!===============================================================================

end module CART__Mesh_Face

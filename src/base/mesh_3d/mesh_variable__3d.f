!> summary:  3D mesh variable
!> author:   Joerg Stiller
!> date:     2024/06/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_Variable__3D
  use Kind_Parameters
  implicit none
  private

  public :: MeshVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D mesh variable

  type MeshVariable_3D
    real(RNP), allocatable :: var(:,:,:,:,:) ! [0:po,0:po,0:po,1:ne,1:nc]
  end type MeshVariable_3D

  !=============================================================================

end module Mesh_Variable__3D

!> summary:  3D multilevel spectral element mesh
!> author:   Joerg Stiller
!> date:     2024/06/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Spectral_Element_Mesh__3D
  use Spectral_Element_Mesh__3D
  use Coarse_To_Fine_Interpolation__1D
  use Fine_To_Coarse_Projection__1D
  implicit none
  private

  public :: ML_SpectralElementMesh_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel spectral element mesh
  !>
  !> `sem(1:l_top)`    sequence of spectral element meshes
  !> `ifc(1:l_top-1)`  fine-to-coarse interpolation operators
  !> `icf(2:l_top)`    coarse-to-fine interpolation operators
  !> `pcf(2:l_top)`    coarse-to-fine L2-projection operators

  type ML_SpectralElementMesh_3D
    type(SpectralElementMesh_3D),       allocatable :: sem(:)
    type(CoarseToFineInterpolation_1D), allocatable :: icf(:)
    type(FineToCoarseProjection_1D),    allocatable :: ifc(:)
    type(FineToCoarseProjection_1D),    allocatable :: pfc(:)
  end type ML_SpectralElementMesh_3D

! contains

  !=============================================================================

end module ML__Spectral_Element_Mesh__3D

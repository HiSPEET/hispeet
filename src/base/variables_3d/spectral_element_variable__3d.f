!> summary:  3D spectral element variable
!> author:   Joerg Stiller
!> date:     2021/7/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spectral_Element_Variable__3D
  use Kind_Parameters, only: RNP
  use Spectral_Element_Mesh__3D
  implicit none
  private

  public :: SpectralElementVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D spectral element variable

  type SpectralElementVariable_3D
    class(SpectralElementMesh_3D), pointer :: se_mesh
    real(RNP), contiguous,  pointer :: value(:,:,:,:,:)  ! (0:po,0:po,0:po,n_elem,n_comp)
    real(RNP), allocatable, private :: storage(:,:,:,:,:)
  end type SpectralElementVariable_3D

! Remarks + issues:
! - original if storage is allocated
! - boundary conditions optional
! - extract slice (subset) as new variable
! – easy access to values, e.g: v(0:,0:,0:,1:,1:) => var % val(:,:,:,:,1:3)
! - extended types: scalar, vector
! - whether/how to include BC ?

contains

  !-----------------------------------------------------------------------------
  !> 3D spectral element variable initialization

  subroutine Init_SpectralElementVariable_3D(this, se_mesh)
    class(SpectralElementVariable_3D), intent(inout) :: this

  end subroutine Init_SpectralElementVariable_3D

  !=============================================================================

end module Spectral_Element_Variable__3D

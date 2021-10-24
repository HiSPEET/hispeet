!> summary:  3D spectral element boundary variable
!> author:   Joerg Stiller
!> date:     2021/??/??
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spectral_Element_Boundary_Variable__3D
  use Kind_Parameters  , only: RNP
  implicit none
  private

  public :: SpectralElementBoundaryVariable_3D

  !-----------------------------------------------------------------------------
  !> Type for keeping the spectral element data associated with one boundary

  type SpectralElementBoundaryData_3D
    real(RNP), contiguous,  pointer :: val(:,:,:,:) !< accessable values
    real(RNP), allocatable, private :: mem(:,:,:,:) !< memory allocated to val
  end type SpectralElementBoundaryData_3D

  !-----------------------------------------------------------------------------
  !> 3D spectral element boundary variable
  !>
  !> Base type for keeping a set of variables on the boundaries of a spectral
  !> element mesh. The mesh, standard operators and metric coefficients are
  !> provided by the pointer `sem`, which can be shared with other entities.
  !> The variables are accessed through
  !>
  !>     bnd(1:nb) % val(0:po,0:po,1:nf,1:nc),
  !>
  !> where
  !>
  !>   - `nb` is the number of boundaries,
  !>   - `po` the polynomial order,
  !>   - `nf` the number of faces, and
  !>   - `nc` the number of components.
  !>
  !> The first three of these dimensions correspond to the spectral element
  !> mesh `sem`, which implies that `nb` is generally different for each
  !> boundary, whereas `po` is constant.
  !> In the current implementation, `nc` is identical for all boundaries.
  !>
  !> Depending on the creation of the boundary variable, the values `val` in
  !> `bnd` can be stored in its own `mem` component or refer to the `mem`
  !> component of another instance.

  type SpectralElementBoundaryVariable_3D
    class(SpectralElementMesh_3D), pointer :: sem
    class(SpectralElementBoundaryData_3D), allocatable :: bnd(:)
  end type SpectralElementBoundaryVariable_3D

  !=============================================================================

end module Spectral_Element_Boundary_Variable__3D

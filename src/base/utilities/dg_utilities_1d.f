!> summary:  Discontinuous 1D spectral element utilities
!> author:   Johannes Stojanow, Joerg Stiller
!> date:     2019/05/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Discontinuous 1D spectral element utilities
!>
!> Provides common routines for Lobatto-based nodal DG-SEM, including
!>
!>   *  mesh generation (`GetMeshPoints`)
!>
!> Some routines require conditions for the left and right boundary points,
!> which are passed in the character array `bc(1:2)`. The following boundary
!> types are supported:
!>
!>   *  periodic  (`'P'`)
!>   *  Dirichlet (`'D'`)
!>   *  Neumann   (`'N'`)
!>
!> While the latter two can be combined as appropriate, periodic conditions
!> must be always specified on both sides, i.e. `bc(1:2) ='P'`.
!>
!===============================================================================

module DG_Utilities_1D
  use Kind_Parameters,   only: RNP
  use Constants,         only: HALF
  use Standard_Operators__1D
  implicit none
  private

  public :: GetMeshPoints

contains

  !-----------------------------------------------------------------------------
  !> Computes the mesh points for a given interval [a,b]
  !>
  !> Requires the initialized standard operators `sop` providing the collocation
  !> points `sop%x`. The polynomial order `sop%po` must match the upper bound of
  !> the first dimension in `x`. The second dimension of `x` determines the
  !> number of elements that are generated.
  !>
  !> Works with all nodal bases.

  subroutine GetMeshPoints(sop, a, b, dx, x)
    class(StandardOperators_1D), intent(in)  :: sop     !< standard operators
    real(RNP),                  intent(in)  :: a       !< left border
    real(RNP),                  intent(in)  :: b       !< right border
    real(RNP),                  intent(out) :: dx      !< element length
    real(RNP), contiguous,      intent(out) :: x(0:,:) !< mesh points x(0:po,1:ne)

    integer   :: k, ne
    real(RNP) :: xe

    ne = size(x,2)
    dx = (b - a) / ne

    associate(xi => sop%x)
      !$omp do
      do k = 1, ne
         xe = a + (k - HALF) * dx      ! element midpoint
         x(:,k) = xe + HALF * dx * xi  ! transformed Lobatto points
      end do
    end associate

  end subroutine GetMeshPoints

  !=============================================================================

end module DG_Utilities_1D

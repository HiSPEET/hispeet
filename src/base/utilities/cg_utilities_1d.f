!> summary:  Continuous 1D spectral element utilities
!> author:   Joerg Stiller
!> date:     2018/09/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Continuous 1D spectral element utilities
!>
!> Provides common routines for GLL-based nodal CG-SEM, including
!>
!>   *  mesh generation (`GetMeshPoints`)
!>   *  computation of the global mass matrix (`GetMassMatrix`)
!>   *  averaging of discontinuous data (`MakeContinuous`)
!>   *  assembly of element integrals (`Assembly`)
!>   *  provision of point weights based on valency (`GetPointWeights`)
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

module CG_Utilities_1D
  use Kind_Parameters,   only: RNP
  use Constants,         only: ONE, HALF
  use Standard_Operators__1D
  implicit none
  private

  public :: GetMeshPoints
  public :: GetMassMatrix
  public :: MakeContinuous
  public :: Assembly
  public :: GetPointWeights

contains

!-------------------------------------------------------------------------------
!> Computes the mesh points for a given interval [a,b]
!>
!> Requires the initialized standard operators `sop` providing the collocation
!> points `sop%x`. The polynomial order `sop%po` must match the upper bound of
!> the first dimension in `x`. The second dimension of `x` determines the number
!> of elements that are generated.
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
       x(:,k) = xe + HALF * dx * xi  ! transformed GLL points
    end do
  end associate

end subroutine GetMeshPoints

!-------------------------------------------------------------------------------
!> Returns the diagonal global mass matrix distributed to element points
!>
!> Current version restricted to GLL bases.

subroutine GetMassMatrix(sop, dx, M, periodic)
  class(StandardOperators_1D), intent(in)  :: sop      !< standard operators
  real(RNP),                  intent(in)  :: dx       !< element length
  real(RNP), contiguous,      intent(out) :: M(0:,:)  !< mass matrix
  logical,     optional,      intent(in)  :: periodic !< assume periodicity [F]

  integer :: e

  !$omp do
  do e = 1, size(M,2)
    M(:,e) = dx/2 * sop % w
  end do
  call Assembly(M, periodic)

end subroutine GetMassMatrix

!-------------------------------------------------------------------------------
!> Average variable over element boundaries
!>
!> Current version restricted to GLL bases.

subroutine MakeContinuous(u, periodic)
  real(RNP), contiguous, intent(inout) :: u(0:,:)  !< element contributions
  logical,     optional, intent(in)    :: periodic !< assume periodic domain [F]

  integer   :: l, po, ne
  real(RNP) :: ua

  po = ubound(u,1)
  ne = ubound(u,2)

  !$omp do
  do l = 1, ne-1
    ua = HALF * (u(po,l) + u(0,l+1))
    u(po,l  ) = ua
    u(0 ,l+1) = ua
  end do

  !$omp single
  if (present(periodic)) then
    if (periodic) then
      ua = HALF * (u(po,ne) + u(0,1))
      u(po,ne) = ua
      u(0 , 1) = ua
    end if
  end if
  !$omp end single

end subroutine MakeContinuous

!-------------------------------------------------------------------------------
!> Assembly of element contributions

subroutine Assembly(v, periodic)
  real(RNP), contiguous, intent(inout) :: v(0:,:) !< element contributions
  logical,     optional, intent(in)    :: periodic !< assume periodic domain [F]

  integer   :: l, po, ne
  real(RNP) :: va

  po = ubound(v,1)
  ne = ubound(v,2)

  !$omp do
  !$acc parallel loop present(v)
  do l = 1, ne-1
    va = v(po,l) + v(0,l+1)
    v(po, l  ) = va
    v( 0, l+1) = va
  end do
  !$acc end parallel

  !$omp single
  if (present(periodic)) then
    if (periodic) then
      va = v(po,ne) + v(0,1)
      v(po, ne) = va
      v( 0,  1) = va
    end if
  end if
  !$omp end single

end subroutine Assembly

!-------------------------------------------------------------------------------
!> Node weights based on inverse valency
!>
!> The provided node weights `w` allow to evaluate scalar products of global
!> coefficient vectors using the local (element based) coefficients.
!> Periodic boundary points are treated as interior ones.

subroutine GetPointWeights(w, periodic)
  real(RNP), contiguous, intent(out) :: w(0:,:)  !< node weights
  logical,     optional, intent(in)  :: periodic !< assume periodic domain [F]

  integer :: j, k, po, ne

  po = ubound(w,1)
  ne = ubound(w,2)

  ! initialize all weights to w = 1 ............................................

  !$omp do
  do k = 1, ne
    do j = 0, po
      w(j,k) = 1
    end do
  end do

  ! set interior element boundaries to w = 1/2 .................................

  !$omp do
  do k = 1, ne-1
    w(po, k  ) = HALF
    w( 0, k+1) = HALF
  end do

  ! periodic points ............................................................

  !$omp single
  if (present(periodic)) then
    if (periodic) then
      w(0 , 1 ) = HALF
      w(po, ne) = HALF
    end if
  end if
  !$omp end single

end subroutine GetPointWeights

!===============================================================================

end module CG_Utilities_1D

!> summary:  Continuous 1D spectral element utilities
!> author:   Joerg Stiller
!> date:     2018/09/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Continuous 1D spectral element utilities
!>
!> @todo
!>   *  OpenMP parallelization
!> @endtodo
!===============================================================================

module CG_Spectral_Element_Utils_1D
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, HALF
  use Standard_Operators_1D
  implicit none
  private

  public :: GetMeshPoints
  public :: GetMassMatrix
  public :: MakeContinuous
  public :: Assembly
  public :: GetPointWeights

contains

!-------------------------------------------------------------------------------
!> Computes the mesh points

subroutine GetMeshPoints(sop, dx, x)
  class(StandardOperators1D), intent(in)  :: sop     !< standard operators
  real(RNP),                  intent(in)  :: dx      !< element length
  real(RNP), contiguous,      intent(out) :: x(0:,:) !< mesh points x(0:po,1:ne)

  integer   :: l
  real(RNP) :: xe

  associate(xi => sop%x)
    do l = 1, ubound(x,2)
       xe = (l - HALF) * dx - ONE    ! element midpoint
       x(:,l) = xe + HALF * dx * xi  ! transformed GLL points
    end do
  end associate

end subroutine GetMeshPoints

!-------------------------------------------------------------------------------
!> Returns the diagonal global mass matrix distributed to element points

subroutine GetMassMatrix(sop, dx, bc, M)
  class(StandardOperators1D), intent(in)  :: sop      !< standard operators
  real(RNP),                  intent(in)  :: dx       !< element length
  character,                  intent(in)  :: bc(:)    !< left/right BC
  real(RNP), contiguous,      intent(out) :: M(0:,:)  !< mass matrix

  integer :: e

  do e = 1, size(M,2)
    M(:,e) = dx/2 * sop % w
  end do
  call Assembly(bc, M)

end subroutine GetMassMatrix

!-------------------------------------------------------------------------------
!> Average variable over element boundaries

subroutine MakeContinuous(bc, u)
  character,             intent(in)    :: bc(:)   !< left/right BC
  real(RNP), contiguous, intent(inout) :: u(0:,:) !< element contributions

  integer   :: l, po, ne
  real(RNP) :: ua

  po = ubound(u,1)
  ne = ubound(u,2)

  do l = 1, ne-1
    ua = HALF * (u(po,l) + u(0,l+1))
    u(po,l  ) = ua
    u(0 ,l+1) = ua
  end do

  if (all(bc == 'P')) then
    ua = HALF * (u(po,ne) + u(0,1))
    u(po,ne) = ua
    u(0 , 1) = ua
  end if

end subroutine MakeContinuous

!-------------------------------------------------------------------------------
!> Assembly of element contributions

subroutine Assembly(bc, v)
  character,             intent(in)    :: bc(:)   !< left/right BC
  real(RNP), contiguous, intent(inout) :: v(0:,:) !< element contributions

  integer   :: l, po, ne
  real(RNP) :: va

  po = ubound(v,1)
  ne = ubound(v,2)

  do l = 1, ne-1
    va = v(po,l) + v(0,l+1)
    v(po, l  ) = va
    v( 0, l+1) = va
  end do

  if (all(bc == 'P')) then
    va = v(po,ne) + v(0,1)
    v(po, ne) = va
    v( 0,  1) = va
  end if

end subroutine Assembly

!------------------------------------------------------------------------------
!> Node weights based on inverse multiplicity

subroutine GetPointWeights(bc, w)
  character,             intent(in)  :: bc(:)   !< left/right BC
  real(RNP), contiguous, intent(out) :: w(0:,:) !< node weights

  integer :: po, ne

  po = ubound(w,1)
  ne = ubound(w,2)

  w = 1

  ! element interfaces
  w(po, 1:ne-1) = HALF
  w( 0, 2:ne  ) = HALF

  ! left boundary
  select case(bc(1))
  case('D')
    w(0, 1) = 0
  case('P')
    w(0, 1) = HALF
  end select

  ! right boundary
  select case(bc(2))
  case('D')
    w(po, ne) = 0
  case('P')
    w(po, ne) = HALF
  end select

end subroutine GetPointWeights

!===============================================================================

end module CG_Spectral_Element_Utils_1D

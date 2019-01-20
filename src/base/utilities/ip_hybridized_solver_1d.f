!> summary:  Direct 1D elliptic solver for hybrid IP/DG-SEM
!> author:   Joerg Stiller
!> date:     2019/01/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>###  Direct 1D elliptic solver for hybrid IP/DG-SEM
!===============================================================================

module IP_Hybridized_Solver_1D
  use Kind_Parameters,  only: RNP
  use Linear_Equations, only: TridiagonalSolver, CyclicTridiagonalSolver
  use IP_Element_Operators_1D
  implicit none
  private

  public :: HybridEllipticSolver

contains

!-------------------------------------------------------------------------------
!> Direct elliptic solver based on hybridization

subroutine HybridEllipticSolver(eop, dx, c, nu, bc, ub, f, u, standby)
  class(IP_ElementOperators1D), intent(in) :: eop !< element operators
  real(RNP), intent(in)  :: dx       !< element width
  real(RNP), intent(in)  :: c        !< coefficient of linear term
  real(RNP), intent(in)  :: nu       !< diffusivity
  character, intent(in)  :: bc(2)    !< boundary conditions {'D','N','P'}
  real(RNP), intent(in)  :: ub(2)    !< Dirichlet boundary values
  real(RNP), intent(in)  :: f(0:,:)  !< source including Neumann BC
  real(RNP), intent(out) :: u(0:,:)  !< solution

  !> optionally keep suboperators for repeated application [F]
  logical, optional, intent(in) :: standby

  ! variables ..................................................................

  ! suboperators
  real(RNP), allocatable, save :: Aib(:,:)     ! interior-boundary op, Â(0:P,1:2)
  real(RNP), allocatable, save :: Aii_inv(:,:) ! inv interior op, Ã⁻¹(0:P,0:P)

  ! flux system
  real(RNP), allocatable :: Af(:,:) ! flux system matrix
  real(RNP), allocatable :: ff(:)   ! flux RHS / solution

  real(RNP) :: dxe(-1:1)

  ! preprocessing ............................................................

  if (allocated(Aib)) then
    if (ubound(Aib,1) /= eop%po) deallocate(Aib, Aii_inv)
  end if

  if (.not. allocated(Aib)) then
    allocate(Aib(0:eop%po, 2), Aii_inv(0:eop%po, 0:eop%po))
    dxe = dx
    call eop % GetEllipticSuboperators(dxe, bc, c, nu, Aib, Aii_inv)
  end if

  ! solution ...................................................................

  call BuildFluxSystem(Aib, Aii_inv, bc, ub, f, Af, ff)

  ! clean-up ...................................................................

  if (present(standby)) then
    if (standby) return
  end if

  deallocate(Aib, Aii_inv)

end subroutine HybridEllipticSolver

!-------------------------------------------------------------------------------
!> Build the flux system

subroutine BuildFluxSystem(Aib, Aii_inv, bc, ub, f, Af, ff)
  real(RNP), intent(in) :: Aib(:,:)     !< interior-boundary operator
  real(RNP), intent(in) :: Aii_inv(:,:) !< inverse interior operator
  character, intent(in) :: bc(2)        !< boundary conditions
  real(RNP), intent(in) :: ub(2)        !< Dirichlet boundary values
  real(RNP), intent(in) :: f(0:,:)      !< source including Neumann BC
  real(RNP), allocatable, intent(out) :: Af(:,:) !< flux system matrix
  real(RNP), allocatable, intent(out) :: ff(:)   !< flux system RHS


end subroutine BuildFluxSystem

!===============================================================================

end module IP_Hybridized_Solver_1D

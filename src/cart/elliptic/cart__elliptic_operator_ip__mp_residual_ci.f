!> summary:  Residual of the IP/DG elliptic operator: const isotropic
!> author:   Joerg Stiller
!> date:     2018/11/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Residual of the IP/DG elliptic operator: const isotropic diffusivity
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_Residual_CI
  use Constants, only: ONE
  use Array_Assignments
  implicit none

contains

!--------------------------------------------------------------------------
!> Computes the residual to given approximation: const isotropic

module subroutine Residual_CI(this, lambda, nu, bc, u, f, r)
  class(EllipticOperator3D_IP), intent(in) :: this
  real(RNP), intent(in)  :: lambda         !< Helmholtz parameter
  real(RNP), intent(in)  :: nu             !< diffusivity
  character, intent(in)  :: bc(:)          !< BC types {P,D,N}
  real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)  :: f(0:,0:,0:,:)  !< right hand side
  real(RNP), intent(out) :: r(0:,0:,0:,:)  !< result

  call this % Apply(lambda, nu, bc, u, r)
  call MergeArrays(-ONE, r, ONE, f)

end subroutine Residual_CI

!===============================================================================

end submodule MP_Residual_CI

!> summary:  Residual of the IP/DG elliptic operator
!> author:   Joerg Stiller
!> date:     2018/11/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Residual of the IP/DG elliptic operator: const isotropic diffusivity
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_Residual
  implicit none

contains

!--------------------------------------------------------------------------
!> Computes the residual to given approximation: const isotropic

module subroutine Residual(this, u, f, r)
  class(EllipticOperator3D_IP), intent(in) :: this
  real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)  :: f(0:,0:,0:,:)  !< right hand side
  real(RNP), intent(out) :: r(0:,0:,0:,:)  !< result

  call this % Apply(u, r)
  call MergeArrays(-ONE, r, ONE, f)

end subroutine Residual

!===============================================================================

end submodule MP_Residual

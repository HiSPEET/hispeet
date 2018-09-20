!> summary:  ISP flow: pressure step
!> author:   Joerg Stiller
!> date:     2018/04/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### ISP flow: pressure step
!===============================================================================

module CART__ISP_Flow__Pressure

  use Kind_Parameters, only: RNP
  use ISP_Flow_Problem
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Pressure__DG_EO

  implicit none
  private

  public :: PressureSolver

contains

!-------------------------------------------------------------------------------
!>  Pressure solver

subroutine PressureSolver(problem, flow_op, tau, f, p, w, consistent, i_max)
  class(FlowProblem),   intent(in)    :: problem      !< flow problem
  class(FlowOperators), intent(inout) :: flow_op      !< flow operators
  real(RNP),            intent(in)    :: tau          !< dt or t
  real(RNP),            intent(in)    :: f(:,:,:,:,:) !< v* or ∂v/∂t*
  real(RNP),            intent(inout) :: p(:,:,:,:)   !< pressure
  real(RNP),            intent(out)   :: w(:,:,:,:,:) !< workspace
  logical,    optional, intent(in)    :: consistent   !< switch to consistent BC
  integer,    optional, intent(in)    :: i_max        !< num MG/CG cycles


  if (flow_op % po_p == flow_op % po_u) then

    call PressureSolver_DG_EO(problem, flow_op, tau, f, p, w, consistent, i_max)

  else
  end if

end subroutine PressureSolver

!===============================================================================

end module CART__ISP_Flow__Pressure

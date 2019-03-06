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
  use CART__ISP_Flow__Pressure__EO
  use CART__ISP_Flow__Pressure__MO

  implicit none
  private

  public :: PressureSolver

contains

!-------------------------------------------------------------------------------
!>  Pressure solver

subroutine PressureSolver(problem, flow_op, dt, v_i, p, w, i_max)
  class(FlowProblem),   intent(in)    :: problem        !< flow problem
  class(FlowOperators), intent(inout) :: flow_op        !< flow operators
  real(RNP),            intent(in)    :: dt             !< time-step size
  real(RNP),            intent(in)    :: v_i(:,:,:,:,:) !< v*
  real(RNP),            intent(inout) :: p(:,:,:,:)     !< pressure
  real(RNP),            intent(out)   :: w(:,:,:,:,:)   !< workspace
  integer,    optional, intent(in)    :: i_max          !< num MG/CG cycles

  if (flow_op % po_p == flow_op % po_u) then
    call PressureSolver_EO(problem, flow_op, dt, v_i, p, w, i_max)
  else
    call PressureSolver_MO(problem, flow_op, dt, v_i, p, w, i_max)
  end if

end subroutine PressureSolver

!===============================================================================

end module CART__ISP_Flow__Pressure

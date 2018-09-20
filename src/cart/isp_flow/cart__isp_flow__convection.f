!> summary:  ISP flow: convective fluxes
!> author:   Joerg Stiller
!> date:     2018/04/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### ISP flow: convective fluxes
!===============================================================================

module CART__ISP_Flow__Convection

  use Kind_Parameters, only: RNP
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Convection__EQ
  use CART__ISP_Flow__Convection__IQ

  implicit none
  private

  public :: WeakConvectiveFlux


contains

!-------------------------------------------------------------------------------
!> Weak divergence of convective fluxes for incompressible flow

subroutine WeakConvectiveFlux(flow_op, u, div_f)

  ! arguments ..................................................................

  class(FlowOperators), intent(in)  :: flow_op    !< flow operators
  real(RNP),            intent(in)  :: u          !< flow variables
  real(RNP),            intent(out) :: div_f      !< flux divergence

  dimension :: u     (:,:,:,:,:)
  dimension :: div_f (:,:,:,:,:)

  ! computation ................................................................

  if (flow_op % po_q == flow_op % po_u) then
      call WeakConvectiveFlux_EQ(flow_op%mesh, flow_op%eop_u, u, div_f)
  else
      call WeakConvectiveFlux_IQ(flow_op, u, div_f)
  end if

end subroutine WeakConvectiveFlux

!===============================================================================

end module CART__ISP_Flow__Convection

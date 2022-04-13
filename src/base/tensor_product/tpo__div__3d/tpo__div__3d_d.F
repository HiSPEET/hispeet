!> summary:  3D Cartesian divergence operator
!> author:   Jerome Michel, Erik Pfister
!> date:     2021/10/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Div__3D_D
  use Kind_Parameters, only: RDP
  use TPO__Div__3D_D__Hand
  use TPO__Div__3D_D__XSMM
  implicit none
  private

  public :: TPO_Div

  interface TPO_Div
    module procedure TPO_Div_D_RDP
  end interface

contains

  subroutine TPO_Div_D_RDP(Ds, Ji, u, v)
    real(RDP), intent(in)  :: Ds(:,:)         !< 1D standard stiffness matrix
    real(RDP), intent(in)  :: Ji(:,:,:,:,:,:) !< inverse Jacobi matrix
    real(RDP), intent(in)  :: u(:,:,:,:,:)    !< 3D vector field
    real(RDP), intent(out) :: v(:,:,:,:)      !< element-wise divergence of u

    integer :: np, ne

    np = size(u,1)
    ne = size(u,4)

    select case(np)

#ifdef __LIBXSMM__
    case(3:)
      call TPO_Div_D_XSMM(np, ne, Ds, Ji, u, v)
#endif

    case default
      call TPO_Div_D_Hand(np, ne, Ds, Ji, u, v)
    end select

  end subroutine TPO_Div_D_RDP

  !=============================================================================

end module TPO__Div__3D_D

!> summary:  3D Cartesian divergence operator
!> author:   Erik Pfister
!> date:     2020/03/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Div__3D_R
  use Kind_Parameters, only: RDP
  use TPO__Div__3D_R__Hand
  use TPO__Div__3D_R__XSMM
  implicit none
  private

  public :: TPO_Div

  interface TPO_Div
    module procedure TPO_Div_R_RDP
  end interface

contains

  subroutine TPO_Div_R_RDP(Ds, dx, u, v)
    real(RDP), intent(in)  :: Ds(:,:)          !< 1D standard stiffness matrix
    real(RDP), intent(in)  :: dx(3)            !< element extensions
    real(RDP), intent(in)  :: u(:,:,:,:,:)     !< 3D vector field
    real(RDP), intent(out) :: v(:,:,:,:)       !< element-wise divergence of u

    integer :: np, ne

    np = size(u,1)
    ne = size(u,4)

    select case(np)

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Div_R_XSMM(np, ne, Ds, dx, u, v)
#endif

    case default
      call TPO_Div_R_Hand(np, ne, Ds, dx, u, v)

    end select

  end subroutine TPO_Div_R_RDP

  !=============================================================================

end module TPO__Div__3D_R

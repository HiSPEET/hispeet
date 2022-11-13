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

  public :: TPO_Div_R

  interface TPO_Div_R
    module procedure TPO_Div_R_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D divergence element operator (R)

  subroutine TPO_Div_R_RDP(Ms, Ds, dx, u, up, v)
    real(RDP), intent(in)  :: Ms(:)         !< 1D standard mass matrix
    real(RDP), intent(in)  :: Ds(:,:)       !< 1D standard stiffness matrix
    real(RDP), intent(in)  :: dx(3)         !< element extensions
    real(RDP), intent(in)  :: u(:,:,:,:,:)  !< 3D vector field
    real(RDP), intent(in)  :: up(:,:,:,:,:) !< exterior traces u⁺ at faces
    real(RDP), intent(out) :: v(:,:,:,:)    !< divergence of u

    optional :: Ms, up

    integer :: np, ne

    np = size(u,1)
    ne = size(u,4)

    select case(np)

#ifdef __LIBXSMM__
    case(4:)
      call TPO_Div_R_XSMM(np, ne, Ms, Ds, dx, u, up, v)
#endif

    case default
      call TPO_Div_R_Hand(np, ne, Ms, Ds, dx, u, up, v)

    end select

  end subroutine TPO_Div_R_RDP

  !=============================================================================

end module TPO__Div__3D_R

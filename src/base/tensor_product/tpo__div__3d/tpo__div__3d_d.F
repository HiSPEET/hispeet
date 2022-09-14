!> summary:  3D Cartesian divergence operator
!> author:   Jerome Michel, Joerg Stiller
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

  !-----------------------------------------------------------------------------
  !> 3D divergence element operator based on LIBXSMM (R)

  subroutine TPO_Div_D_RDP(Ms, Ds, Jd, Ji, a, n, u, up, v)
    real(RDP), intent(in)  :: Ms(:)           !< 1D standard mass matrix
    real(RDP), intent(in)  :: Ds(:,:)         !< 1D standard diff matrix
    real(RDP), intent(in)  :: Jd(:,:,:,:)     !< element Jacobian determinant
    real(RDP), intent(in)  :: Ji(:,:,:,:,:,:) !< inverse Jacobi matrix
    real(RDP), intent(in)  :: a(:,:,:,:)      !< face area coefficients
    real(RDP), intent(in)  :: n(:,:,:,:,:)    !< face unit normal vectors
    real(RDP), intent(in)  :: u(:,:,:,:,:)    !< 3D vector field
    real(RDP), intent(in)  :: up(:,:,:,:,:)   !< exterior traces u⁺ at faces
    real(RDP), intent(out) :: v(:,:,:,:)      !< divergence of u

    optional :: Ms, Jd, a, n, up

    integer :: np, ne

    np = size(u,1)
    ne = size(u,4)

    select case(np)

#ifdef __LIBXSMM__
    case(3:)
      call TPO_Div_D_XSMM(np, ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
#endif

    case default
      call TPO_Div_D_Hand(np, ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    end select

  end subroutine TPO_Div_D_RDP

  !=============================================================================

end module TPO__Div__3D_D

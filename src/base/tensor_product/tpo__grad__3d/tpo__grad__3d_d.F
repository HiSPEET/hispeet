!> summary:   3D Cartesian gradient operator (D)
!> author:    Jerome Michel, Erik Pfister
!> date:      2021/10/12
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Grad__3D_D
  use Kind_Parameters, only: RDP
  use TPO__Grad__3D_D__Hand
  use TPO__Grad__3D_D__XSMM
  implicit none
  private

  public :: TPO_Grad

  interface TPO_Grad
    module procedure TPO_Grad_D_RDP
  end interface

contains

  subroutine TPO_Grad_D_RDP(Ds, Ji, u, v)
    real(RDP), intent(in)  :: Ds(:,:)         !< 1D standard stiffness matrix
    real(RDP), intent(in)  :: Ji(:,:,:,:,:,:) !< inverse Jacobi matrix
    real(RDP), intent(in)  :: u(:,:,:,:)      !< 3D Scalar Field
    real(RDP), intent(out) :: v(:,:,:,:,:)    !< element-wise gradient of u

    integer :: np, ne

    np = size(u,1)
    ne = size(u,4)

    select case(np)

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Grad_D_XSMM(np, ne, Ds, Ji, u, v)
#endif

    case default
      call TPO_Grad_D_Hand(np, ne, Ds, Ji, u, v)
    end select

  end subroutine TPO_Grad_D_RDP

  !=============================================================================

end module TPO__Grad__3D_D

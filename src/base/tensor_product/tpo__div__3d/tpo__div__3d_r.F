!> summary:  3D Cartesian divergence operator
!> author:   Erik Pfister
!> date:     2020/03/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Div__3D_R
  use Kind_Parameters, only: RNP
  use TPO__Div__3D_R__Hand
  use TPO__Div__3D_R__XSMM
  implicit none
  private

  public :: TPO_Div_R

contains

  subroutine TPO_Div_R(Ds, dx, u, v)
    use Kind_Parameters,  only: RNP
    real(RNP), intent(in)  :: Ds(:,:)          !< 1D standard stiffness matrix
    real(RNP), intent(in)  :: dx(3)            !< element extensions
    real(RNP), intent(in)  :: u(:,:,:,:,:)     !< 3D vector field
    real(RNP), intent(out) :: v(:,:,:,:)       !< element-wise divergence of u

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

  end subroutine TPO_Div_R

  !=============================================================================

end module TPO__Div__3D_R

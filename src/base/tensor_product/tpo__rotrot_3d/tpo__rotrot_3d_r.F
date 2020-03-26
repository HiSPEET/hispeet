!> summary:  3D Cartesian rotrot operator
!> author:   Erik Pfister
!> date:     2020/03/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__RotRot_3d_R
  use Kind_Parameters, only: RNP
  use TPO__RotRot_3d_R__Hand
  use TPO__RotRot_3d_R__XSMM
  implicit none
  private

  public :: TPO_RotRot_R

contains

  subroutine TPO_RotRot_R(Ds, dx, u, v)
    use Kind_Parameters,  only: RNP
    real(RNP), intent(in)  :: Ds(:,:)          !< 1D standard stiffness matrix
    real(RNP), intent(in)  :: dx(3)            !< element extensions
    real(RNP), intent(in)  :: u(:,:,:,:,:)     !< operand
    real(RNP), intent(out) :: v(:,:,:,:,:)     !< result

    integer :: np, ne

    np = size(u,1)
    ne = size(u,4)

    select case(np)

#ifdef __LIBXSMM__
    case(8:)
      call TPO_RotRot_R_XSMM(np, ne, Ds, dx, u, v)
#endif

    case default
      call TPO_RotRot_R_Hand(np, ne, Ds, dx, u, v)

    end select

  end subroutine TPO_RotRot_R

  !=============================================================================

end module TPO__RotRot_3d_R

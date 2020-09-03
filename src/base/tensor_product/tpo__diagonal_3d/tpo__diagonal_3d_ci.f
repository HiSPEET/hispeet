!> summary:  3D scaled constant isotropic diagonal operator
!> author:   Joerg Stiller
!> date:     2020/06/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diagonal_3d_CI
  use Kind_Parameters, only: RNP
  use TPO__Diagonal_3d_CI__Gen
  implicit none
  private

  public :: TPO_Diagonal_CI

contains

  !-----------------------------------------------------------------------------
  !> 3D scaled constant isotropic diagonal operator, v = s DxDxD u

  subroutine TPO_Diagonal_CI(s, D, u, v)
    use Kind_Parameters,  only: RNP
    real(RNP), intent(in)  :: s          !< scaling factor
    real(RNP), intent(in)  :: D(:)       !< diagonal 1D operator
    real(RNP), intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), intent(out) :: v(:,:,:,:) !< result

    integer :: np, ne

    np = size(D)
    ne = size(u,4)

    call TPO_Diagonal_CI_Gen(np, ne, s, D, u, v)

  end subroutine TPO_Diagonal_CI

  !=============================================================================

end module TPO__Diagonal_3d_CI

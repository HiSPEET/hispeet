!> summary:  3D scaled constant isotropic diagonal operator
!> author:   Joerg Stiller
!> date:     2020/06/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diagonal__3D_CI
  use Kind_Parameters, only: RDP
  use TPO__Diagonal__3D_CI__Gen
  implicit none
  private

  public :: TPO_Diagonal_CI_RDP

contains

  !-----------------------------------------------------------------------------
  !> 3D scaled constant isotropic diagonal operator, v = s DxDxD u

  subroutine TPO_Diagonal_CI_RDP(s, D, u, v)
    real(RDP), intent(in)  :: s          !< scaling factor
    real(RDP), intent(in)  :: D(:)       !< diagonal 1D operator
    real(RDP), intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), intent(out) :: v(:,:,:,:) !< result

    integer :: np, ne

    np = size(D)
    ne = size(u,4)

    call TPO_Diagonal_CI_Gen(np, ne, s, D, u, v)

  end subroutine TPO_Diagonal_CI_RDP

  !=============================================================================

end module TPO__Diagonal__3D_CI

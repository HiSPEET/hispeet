!> summary:  3D generic isotropic averaging operator
!> author:   Joerg Stiller
!> date:     2022/10/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Average__3D_I
  use Kind_Parameters, only: RDP
  use TPO__Average__3D_I__Gen
  implicit none
  private

  public :: TPO_Average_I_RDP

contains

  !-----------------------------------------------------------------------------
  !> 3D scaled isotropic averaging operator

  subroutine TPO_Average_I_RDP(s, w, u, v)
    real(RDP), intent(in)  :: s          !< scaling factor
    real(RDP), intent(in)  :: w(:)       !< 1D averaging operator
    real(RDP), intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), intent(out) :: v(:)       !< result

    integer :: np, ne

    np = size(w)
    ne = size(u,4)

    call TPO_Average_I_Gen(np, ne, s, w, u, v)

  end subroutine TPO_Average_I_RDP

  !=============================================================================

end module TPO__Average__3D_I

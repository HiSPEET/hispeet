!> summary:  Evaluation of  v = (AxAxA) u   with 3D mesh variables
!> author:   Joerg Stiller
!> date:     2020/05/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__AAA__3D
  use Kind_Parameters, only: RDP
  use TPO__AAA__3D__Hand
  use TPO__AAA__3D__XSMM
  implicit none
  private

  public :: TPO_AAA

  interface TPO_AAA
    module procedure TPO_AAA_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D AxAxA (interpolation) operator

  subroutine TPO_AAA_RDP(A, u, v)
    real(RDP), intent(in)  :: A(:,:)     !< 1D operator
    real(RDP), intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), intent(out) :: v(:,:,:,:) !< result

    integer :: np_max

    np_max = max( size(u,1), size(v,1) )

    select case(np_max)

#ifdef __LIBXSMM__
    case(8:)
      call TPO_AAA_XSMM(A, u, v)
#endif

    case default
      call TPO_AAA_Hand(A, u, v)

    end select

  end subroutine TPO_AAA_RDP

  !=============================================================================

end module TPO__AAA__3D

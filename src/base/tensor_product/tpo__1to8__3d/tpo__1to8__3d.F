!> summary:  3d interpolation operator for regular 1:8 h-refinement
!> author:   Joerg Stiller
!> date:     2024/07/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__1to8__3D
  use Kind_Parameters, only: RDP
  use TPO__1to8__3D__Hand
  use TPO__1to8__3D__XSMM
  implicit none
  private

  public :: TPO_1to8

  interface TPO_1to8
    module procedure TPO_1to8_RDP, TPO_1to8_RDP_Array
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D 1:8 interpolation operator for RDP

  subroutine TPO_1to8_RDP(A, u, v)
    real(RDP), intent(in)  :: A(:,:,:)   !< 1D operator
    real(RDP), intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), intent(out) :: v(:,:,:,:) !< result

    integer :: np_max

    np_max = max( size(u,1), size(v,1) )

    select case(np_max)

#ifdef __LIBXSMM__
    case(8:)
      call TPO_1to8_XSMM(A, u, v)
#endif

    case default
      call TPO_1to8_Hand(A, u, v)

    end select

  end subroutine TPO_1to8_RDP

  !-----------------------------------------------------------------------------
  !> 3D 1:8 interpolation operator for RDP array variables

  subroutine TPO_1to8_RDP_Array(A, u, v)
    real(RDP), intent(in)  :: A(:,:,:)     !< 1D operator
    real(RDP), intent(in)  :: u(:,:,:,:,:) !< operand
    real(RDP), intent(out) :: v(:,:,:,:,:) !< result

    integer :: nc, ne, ni, np

    ni = size(A,1)
    np = size(A,2)
    ne = size(u,4)
    nc = size(u,5)

    call ArrayOperator(ni, np, ne, nc, A, u, v)

  contains

    subroutine ArrayOperator(ni, np, ne, nc, A, u, v)
      integer,   intent(in)  :: ni, np, ne, nc
      real(RDP), intent(in)  :: A(ni,np,2)          !< 1D operator
      real(RDP), intent(in)  :: u(np,np,np,ne*nc)   !< operand
      real(RDP), intent(out) :: v(ni,ni,ni,ne*nc*8) !< result

      call TPO_1to8_RDP(A, u, v)

    end subroutine ArrayOperator

  end subroutine TPO_1to8_RDP_Array

  !=============================================================================

end module TPO__1to8__3D

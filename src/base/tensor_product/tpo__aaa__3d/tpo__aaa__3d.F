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
    module procedure TPO_AAA_RDP, TPO_AAA_RDP_Array
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D AxAxA (interpolation) operator for RDP

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

  !-----------------------------------------------------------------------------
  !> 3D AxAxA (interpolation) operator for RDP array variables

  subroutine TPO_AAA_RDP_Array(A, u, v)
    real(RDP), intent(in)  :: A(:,:)       !< 1D operator
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
      integer, intent(in) :: ni, np, ne, nc
      real(RDP), intent(in)  :: A(ni,np)          !< 1D operator
      real(RDP), intent(in)  :: u(np,np,np,ne*nc) !< operand
      real(RDP), intent(out) :: v(ni,ni,ni,ne*nc) !< result

      call TPO_AAA_RDP(A, u, v)

    end subroutine ArrayOperator

  end subroutine TPO_AAA_RDP_Array

  !=============================================================================

end module TPO__AAA__3D

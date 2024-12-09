!> summary:  3d projection operator for regular 8:1 h- or hp-coarsening
!> author:   Joerg Stiller
!> date:     2024/08/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__8To1__3D
  use Kind_Parameters, only: RDP
  use TPO__8To1__3D__Gen
  use TPO__8To1__3D__XSMM
  implicit none
  private

  public :: TPO_8To1

  interface TPO_8To1
    module procedure TPO_8To1_RDP, TPO_8To1_RDP_Array
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D 8:1 projection operator for RDP scalar variables

  subroutine TPO_8To1_RDP(A, u, v)
    real(RDP), contiguous, intent(in)  :: A(:,:,:)   !< 1D projection operator
    real(RDP), contiguous, intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), contiguous, intent(out) :: v(:,:,:,:) !< result

    integer :: np_max

    np_max = max( size(A,1), size(A,2) )

    select case(np_max)

#ifdef __LIBXSMM__
    case(8:)
      call TPO_8To1_XSMM(A, u, v)
#endif

    case default
      call TPO_8To1_Gen(A, u, v)

    end select

  end subroutine TPO_8To1_RDP

  !-----------------------------------------------------------------------------
  !> 3D 1:8 interpolation operator for RDP array variables

  subroutine TPO_8To1_RDP_Array(A, u, v)
    real(RDP), contiguous, intent(in)  :: A(:,:,:)     !< 1D projection operator
    real(RDP), contiguous, intent(in)  :: u(:,:,:,:,:) !< operand
    real(RDP), contiguous, intent(out) :: v(:,:,:,:,:) !< result

    integer :: nc, ne, np, nq

    np = size(u,1)
    nq = size(v,1)
    ne = size(v,4)
    nc = size(v,5)

    call ArrayOperator(np, nq, ne, nc, A, u, v)

  contains

    subroutine ArrayOperator(np, nq, ne, nc, A, u, v)
      integer,   intent(in)  :: np, nq, ne, nc
      real(RDP), intent(in)  :: A(nq,np,2)          !< 1D projection operator
      real(RDP), intent(in)  :: u(np,np,np,ne*nc*8) !< operand
      real(RDP), intent(out) :: v(nq,nq,nq,ne*nc)   !< result

      call TPO_8To1_RDP(A, u, v)

    end subroutine ArrayOperator

  end subroutine TPO_8To1_RDP_Array

  !=============================================================================

end module TPO__8To1__3D

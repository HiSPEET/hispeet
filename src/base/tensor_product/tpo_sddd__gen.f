!> summary:  Generic s(DxDxD) operator with loops
!> author:   Joerg Stiller
!> date:     2017/05/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Generic s(DxDxD) operator with loops
!===============================================================================

subroutine TPO_sDDD__gen(np, ne, s, D, u, v)

  !-----------------------------------------------------------------------------
  ! used modules

  use Kind_Parameters, only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: np             !< number of points per direction
  integer,   intent(in)  :: ne             !< number of elements
  real(RNP), intent(in)  :: s              !< scaling factor
  real(RNP), intent(in)  :: D(np)          !< diagonal 1D operator
  real(RNP), intent(in)  :: u(np,np,np,ne) !< collapsed 3D operand
  real(RNP), intent(out) :: v(np,np,np,ne) !< collapsed 3D result

  !-----------------------------------------------------------------------------
  ! local variables

  integer :: e, i, j, k

  !$acc data present(u,v) copyin(D) async
  !$acc parallel loop collapse(4) async
  !$omp do
  do e = 1, ne
    do k = 1, np
    do j = 1, np
    do i = 1, np
      v(i,j,k,e) = s * D(k) * D(j) * D(i) * u(i,j,k,e)
    end do
    end do
    end do
  end do
  !$omp end do
  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine TPO_sDDD__gen

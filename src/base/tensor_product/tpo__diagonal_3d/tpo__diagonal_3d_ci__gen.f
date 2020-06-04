!> summary:  3D generic scaled constant isotropic diagonal operator
!> author:   Joerg Stiller
!> date:     2020/06/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diagonal_3d_CI__Gen
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_Diagonal_CI_Gen

contains

  !-----------------------------------------------------------------------------
  !> 3D generic scaled constant isotropic diagonal operator, v = s DxDxD u

  subroutine TPO_Diagonal_CI_Gen(np, ne, s, D, u, v)
    integer,   intent(in)  :: np             !< number of points per direction
    integer,   intent(in)  :: ne             !< number of elements
    real(RNP), intent(in)  :: s              !< scaling factor
    real(RNP), intent(in)  :: D(np)          !< diagonal 1D operator
    real(RNP), intent(in)  :: u(np,np,np,ne) !< collapsed 3D operand
    real(RNP), intent(out) :: v(np,np,np,ne) !< collapsed 3D result

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

  end subroutine TPO_Diagonal_CI_Gen

  !=============================================================================

end module TPO__Diagonal_3d_CI__Gen


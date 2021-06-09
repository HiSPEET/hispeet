!> summary:  3D generic scaled constant isotropic diagonal operator
!> author:   Joerg Stiller
!> date:     2020/06/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

!-----------------------------------------------------------------------------
!> 3D generic scaled constant isotropic diagonal operator, v = s DxDxD u

subroutine TPO_Diagonal_CI_Gen_RWP(np, ne, s, D, u, v)
  integer,   intent(in)  :: np             !< number of points per direction
  integer,   intent(in)  :: ne             !< number of elements
  real(RWP), intent(in)  :: s              !< scaling factor
  real(RWP), intent(in)  :: D(np)          !< diagonal 1D operator
  real(RWP), intent(in)  :: u(np,np,np,ne) !< collapsed 3D operand
  real(RWP), intent(out) :: v(np,np,np,ne) !< collapsed 3D result

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

end subroutine TPO_Diagonal_CI_Gen_RWP

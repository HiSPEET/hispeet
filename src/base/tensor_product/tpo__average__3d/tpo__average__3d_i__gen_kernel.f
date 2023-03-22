!> summary:  3D generic isotropic averaging operator
!> author:   Joerg Stiller
!> date:     2022/10/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

!-------------------------------------------------------------------------------
!> 3D generic isotropic averaging operator

subroutine TPO_Average_I_Gen_RWP(np, ne, w, u, v)
  integer,   intent(in)  :: np             !< number of points per direction
  integer,   intent(in)  :: ne             !< number of elements
  real(RWP), intent(in)  :: w(np)          !< 1D averaging operator
  real(RWP), intent(in)  :: u(np,np,np,ne) !< operand
  real(RWP), intent(out) :: v(ne)          !< elementwise average

  real(RWP) :: s
  integer   :: e, i, j, k

  s = 1 / sum(w)**3

  !$omp do
  do e = 1, ne
    v(e) = 0
    do k = 1, np
    do j = 1, np
    do i = 1, np
      v(e) = v(e) + w(k) * w(j) * w(i) * u(i,j,k,e)
    end do
    end do
    end do
    v(e) = s * v(e)
  end do

end subroutine TPO_Average_I_Gen_RWP

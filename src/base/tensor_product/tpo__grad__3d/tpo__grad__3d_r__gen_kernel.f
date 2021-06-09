!> summary:  Generic gradient of a vector field: 3D Cartesian equidistant
!> author:   Joerg Stiller, Erik Pfister
!> date:     2017/01/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Grad_R_Gen_RWP(np, ne, Ds, dx, u, v)
  integer,   intent(in)  :: np               !< number of points per direction
  integer,   intent(in)  :: ne               !< number of elements
  real(RWP), intent(in)  :: Ds(np,np)        !< 1D standard diff matrix
  real(RWP), intent(in)  :: dx(3)            !< element extensions
  real(RWP), intent(in)  :: u(np,np,np,ne)   !< 3D scalar field
  real(RWP), intent(out) :: v(np,np,np,ne,3) !< element-wise gradient of u

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP) :: A(np,np)
  real(RWP) :: g(3), tmp1, tmp2, tmp3

  integer :: e, i, j, k, p
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (np < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  ! Initialization
  A = transpose(Ds)

  ! metric coefficients
  g = 2 / dx

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(A,g) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker

  !$omp do private(e)
  do e = 1, ne

    ! v1 = du/dx1, v2 = du/dx2 .................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do p = 1, np
        tmp1 = tmp1 + A(p,i) * u(p,j,k,e)
        tmp2 = tmp2 + A(p,j) * u(i,p,k,e)
      end do
      v(i,j,k,e,1) = g(1) * tmp1
      v(i,j,k,e,2) = g(2) * tmp2

    end do
    end do
    end do

    ! v3 = du/dx3 ..............................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp3 = 0
      do p = 1, np
        tmp3 = tmp3 + A(p,k) * u(i,j,p,e)
      end do
      v(i,j,k,e,3) = g(3) * tmp3

    end do
    end do
    end do

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine TPO_Grad_R_Gen_RWP

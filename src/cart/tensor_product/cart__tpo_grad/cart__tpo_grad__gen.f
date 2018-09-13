!> summary:  Generic gradient of a scalar field: 3D Cartesian equidistant
!> author:   Joerg Stiller
!> date:     2017/02/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Generic gradient of a scalar field: 3D Cartesian equidistant
!===============================================================================

subroutine CART__TPO_Grad__gen(np, ne, Ds, dx, u, v)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters, only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: np               !< number of points per direction
  integer,   intent(in)  :: ne               !< number of elements
  real(RNP), intent(in)  :: Ds(np,np)        !< 1D standard diff matrix
  real(RNP), intent(in)  :: dx(3)            !< element extensions
  real(RNP), intent(in)  :: u(np,np,np,ne)   !< 3D scalar field
  real(RNP), intent(out) :: v(np,np,np,ne,3) !< element-wise gradient of u

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), allocatable :: Dt(:,:)
  real(RNP) :: g(3), tmp1, tmp2, tmp3

  integer :: e, i, j, k, m
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (np < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  ! workspace
  allocate(Dt(np,np))

  ! transposed 1D standard diff matrix
  do j = 1, np
  do i = 1, np
    Dt(j,i) = Ds(i,j)
  end do
  end do

  ! metric coefficients
  g = 2 / dx

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(Dt,g) async
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
      do m = 1, np
        tmp1 = tmp1 + Dt(m,i) * u(m,j,k,e)
        tmp2 = tmp2 + Dt(m,j) * u(i,m,k,e)
      end do
      v(i,j,k,e,1) = g(1) * tmp1
      v(i,j,k,e,2) = g(2) * tmp2

    end do
    end do
    end do

    ! v += du3/dx3 .............................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp3 = 0
      do m = 1, np
        tmp3 = tmp3 + Dt(m,k) * u(i,j,m,e)
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

end subroutine CART__TPO_Grad__gen

!> summary:  Generic rotrot of a vector field: 3D Cartesian equidistant
!> author:   Joerg Stiller, Erik Pfister
!> date:     2017/01/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Generic rotrot of a vector field: 3D Cartesian equidistant
!===============================================================================


module TPO__RotRot_3d_R__Gen
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_RotRot_R_Gen

contains


subroutine TPO_RotRot_R_Gen(np, ne, Ds, dx, u, v)

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
  real(RNP), intent(in)  :: u(np,np,np,ne,3) !< 3D vector field
  real(RNP), intent(out) :: v(np,np,np,ne,3) !< element-wise divergence of u

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP) :: A(np,np), DA(np,np)
  real(RNP) :: div_u(np,np,np)
  real(RNP) :: g(3), tmp, tmp1, tmp2

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
  A  = transpose(Ds)

  do j = 1, np
  do i = 1, np
    DA(j,i) = sum(Ds(i,:) * Ds(:,j))
  end do
  end do

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

    ! div_u = du1/dx1 ..........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,i) * u(p,j,k,e,1)
      end do
      div_u(i,j,k) = g(1) * tmp

    end do
    end do
    end do

    ! div_u += du2/dx2 .........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,j) * u(i,p,k,e,2)
      end do
      div_u(i,j,k) = div_u(i,j,k) + g(2) * tmp

    end do
    end do
    end do

    ! div_u += du3/dx3 .........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,k) * u(i,j,p,e,3)
      end do
      div_u(i,j,k) = div_u(i,j,k) + g(3) * tmp

    end do
    end do
    end do

    ! v1 = d(div_u)/dx1 ........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,i) * div_u(p,j,k)
      end do
      v(i,j,k,e,1) = g(1) * tmp

    end do
    end do
    end do

    ! v1 -= d( du1/dx1 )/dx1 + d( du1/dx2 )/dx2 ................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do p = 1, np
        tmp1 = tmp1 + DA(p,i) * u(p,j,k,e,1)
        tmp2 = tmp2 + DA(p,j) * u(i,p,k,e,1)
      end do
      v(i,j,k,e,1) = v(i,j,k,e,1) - (g(1)*g(1)*tmp1 + g(2)*g(2)*tmp2)

    end do
    end do
    end do

    ! v1 -= d( du1/dx3 )/dx3 ...................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + DA(p,k) * u(i,j,p,e,1)
      end do
      v(i,j,k,e,1) = v(i,j,k,e,1) - g(3)*g(3)*tmp

    end do
    end do
    end do

    ! v2 = d(div_u)/dx2 ........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,j) * div_u(i,p,k)
      end do
      v(i,j,k,e,2) = g(2) * tmp

    end do
    end do
    end do

    ! v2 -= d( du2/dx1 )/dx1 + d( du2/dx2 )/dx2 ................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do p = 1, np
        tmp1 = tmp1 + DA(p,i) * u(p,j,k,e,2)
        tmp2 = tmp2 + DA(p,j) * u(i,p,k,e,2)
      end do
      v(i,j,k,e,2) = v(i,j,k,e,2) - (g(1)*g(1)*tmp1 + g(2)*g(2)*tmp2)

    end do
    end do
    end do

    ! v2 -= d( du2/dx3 )/dx3 ...................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + DA(p,k) * u(i,j,p,e,2)
      end do
      v(i,j,k,e,2) = v(i,j,k,e,2) - g(3)*g(3)*tmp

    end do
    end do
    end do

    ! v3 = d(div_u)/dx3 ........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,k) * div_u(i,j,p)
      end do
      v(i,j,k,e,3) = g(3) * tmp

    end do
    end do
    end do

    ! v3 -= d( du3/dx1 )/dx1 + d( du3/dx2 )/dx2 ................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do p = 1, np
        tmp1 = tmp1 + DA(p,i) * u(p,j,k,e,3)
        tmp2 = tmp2 + DA(p,j) * u(i,p,k,e,3)
      end do
      v(i,j,k,e,3) = v(i,j,k,e,3) - (g(1)*g(1)*tmp1 + g(2)*g(2)*tmp2)

    end do
    end do
    end do

    ! v3 -= d( du3/dx3 )/dx3 ...................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + DA(p,k) * u(i,j,p,e,3)
      end do
      v(i,j,k,e,3) = v(i,j,k,e,3) - g(3)*g(3)*tmp

    end do
    end do
    end do

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine TPO_RotRot_R_Gen

end module TPO__RotRot_3d_R__Gen

!> summary:  Generic divergence of a vector field: 3D Cartesian equidistant
!> author:   Joerg Stiller, Erik Pfister
!> date:     2017/01/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Generic divergence of a vector field: 3D Cartesian equidistant
!===============================================================================


module TPO__Div_3d_R__Gen
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_Div_R_Gen

contains


subroutine TPO_Div_R_Gen(np, ne, Ds, dx, u, v)

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
  real(RNP), intent(out) :: v(np,np,np,ne)   !< element-wise divergence of u
!  logical, optional, intent(in) :: transp    !< if present, apply Ds^T

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP) :: A(np,np)
  real(RNP) :: g(3), tmp

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

    ! v = du1/dx1 ..............................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do m = 1, np
        tmp = tmp + A(m,i) * u(m,j,k,e,1)
      end do
      v(i,j,k,e) = g(1) * tmp

    end do
    end do
    end do

    ! v += du2/dx2 .............................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do m = 1, np
        tmp = tmp + A(m,j) * u(i,m,k,e,2)
      end do
      v(i,j,k,e) = v(i,j,k,e) + g(2) * tmp

    end do
    end do
    end do

    ! v += du3/dx3 .............................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do m = 1, np
        tmp = tmp + A(m,k) * u(i,j,m,e,3)
      end do
      v(i,j,k,e) = v(i,j,k,e) + g(3) * tmp

    end do
    end do
    end do

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine TPO_Div_R_Gen

end module TPO__Div_3d_R__Gen

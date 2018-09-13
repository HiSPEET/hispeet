!> summary:  Generic element Diffusion operator: 3D Cartesian equidistant
!> author:   Immo Huismann, Joerg Stiller
!> date:     2015/04/07, revised 2017/01/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>###   Generic element Diffusion operator: 3D Cartesian equidistant
!>
!> @note
!> When compiled with OpenACC and run on a GPU, any call to this
!> procedure must be followed by an "acc wait" before accessing u
!> or v from CPU or any other than the default accelerator queue.
!> @endnote
!===============================================================================

subroutine CART__TPO_Diffusion__gen(np, ne, Ms, Ls, lambda, nu, dx, u, v)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters,  only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: np             !< number of points per direction
  integer,   intent(in)  :: ne             !< number of elements
  real(RNP), intent(in)  :: Ms(np)         !< 1D standard mass matrix
  real(RNP), intent(in)  :: Ls(np,np)      !< 1D standard stiffness matrix
  real(RNP), intent(in)  :: lambda         !< Helmholtz parameter
  real(RNP), intent(in)  :: nu             !< diffusivity
  real(RNP), intent(in)  :: dx(3)          !< element extensions
  real(RNP), intent(in)  :: u(np,np,np,ne) !< operand
  real(RNP), intent(out) :: v(np,np,np,ne) !< result

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), allocatable :: M(:,:,:), M_u(:,:,:), Lm(:,:)
  real(RNP) :: g(3), tmp

  integer :: e, i, j, k, p
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if(np < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  ! workspace
  allocate(M(np,np,np), M_u(np,np,np), Lm(np,np))

  ! element mass matrix
  tmp = product(dx) / 8
  do k = 1, np
  do j = 1, np
  do i = 1, np
    M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  ! mass-weighted stiffness matrix: Lm = Ms^-1 Ls = (Ls Ms^-1)^T
  do j = 1, np
  do i = 1, np
    Lm(i,j) = Ls(i,j) / Ms(i)
  end do
  end do

  ! coefficients
  g = 4 * nu / dx**2

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(g,M,Lm) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker private(M_u)

  !$omp do private(e)
  do e = 1, ne

    ! M u and lambda M u .......................................................

    !$acc loop collapse(3) independent vector
    do k = 1, np
    do j = 1, np
    do i = 1, np
      M_u(i,j,k) = M(i,j,k) * u(i,j,k,e)
      v(i,j,k,e) = lambda * M_u(i,j,k)
    end do
    end do
    end do

    ! direction 1 ..............................................................

    !$acc loop collapse(3) independent vector
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Lm(p,i) * M_u(p,j,k)
      end do
      v(i,j,k,e) = v(i,j,k,e) + g(1) * tmp
    end do
    end do
    end do

    ! direction 2 ..............................................................

    !$acc loop collapse(3) independent vector
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Lm(p,j) * M_u(i,p,k)
      end do
      v(i,j,k,e) = v(i,j,k,e) + g(2) * tmp
    end do
    end do
    end do

    ! direction 3 ..............................................................

    !$acc loop collapse(3) independent vector
    do k = 1, np
    do j = 1, np
    do i = 1, np
      tmp = 0
      do p = 1, np
        tmp = tmp + Lm(p,k) * M_u(i,j,p)
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

end subroutine CART__TPO_Diffusion__gen

!> summary:  Isotropic Schwarz operator: bf-marching, separate, generic
!> author:   Joerg Stiller
!> date:     2017/10/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Isotropic Schwarz operator: bf-marching, separate, generic
!===============================================================================

subroutine CART__TPO_Schwarz_Old__gen( ns, nc, ne, S, V1, V2, V3, W, &
                                       ec, dx, lambda, nu, f, u      )

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters, only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: ns              !< dimension of S
  integer,   intent(in)  :: nc              !< number of configurations
  integer,   intent(in)  :: ne              !< number of elements
  real(RNP), intent(in)  :: S(ns,ns,nc)     !< 1D eigenvectors/subdomains
  real(RNP), intent(in)  :: V1(ns,nc)       !< 1D eigenvalues in direction 1
  real(RNP), intent(in)  :: V2(ns,nc)       !< 1D eigenvalues in direction 2
  real(RNP), intent(in)  :: V3(ns,nc)       !< 1D eigenvalues in direction 3
  real(RNP), intent(in)  :: W(ns,nc)        !< 1D weights
  integer,   intent(in)  :: ec(3,ne)        !< element-boundary configurations
  real(RNP), intent(in)  :: dx(3)           !< element extensions
  real(RNP), intent(in)  :: lambda          !< Helmholtz parameter
  real(RNP), intent(in)  :: nu(ne)          !< subdomain diffusivities
  real(RNP), intent(in)  :: f(ns,ns,ns,ne)  !< RHS
  real(RNP), intent(out) :: u(ns,ns,ns,ne)  !< solution

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), allocatable :: WS_t(:,:,:), z(:,:,:)
  real(RNP) :: g0, g1, g2, g3
  real(RNP) :: tmp

  integer :: e, i, j, k, m
  integer :: c, c1, c2, c3
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (ns < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  allocate(WS_t(ns,ns,nc), z(ns,ns,ns))

  ! transposed 1D operator
  do c = 1, nc
  do j = 1, ns
  do i = 1, ns
    WS_t(j,i,c) = W(i,c) * S(i,j,c)
  end do
  end do
  end do

  ! metric coefficients
  g0 = dx(1) * dx(2) * dx(3)
  g1 = dx(2) * dx(3) / dx(1)
  g2 = dx(3) * dx(1) / dx(2)
  g3 = dx(1) * dx(2) / dx(3)

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(ec,nu,f,u) copyin(V1,V2,V3,S,WS_t) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker private(z)

  !$omp do private(e)
  elements: do e = 1, ne

    ! element-boundary configurations
    c1 = ec(1,e)
    c2 = ec(2,e)
    c3 = ec(3,e)

    ! z = S^t x I x I f ........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + S(m,k,c3) * f(i,j,m,e)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    ! u = I x S^t x I z ........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + S(m,j,c2) * z(i,m,k)
      end do
      u(i,j,k,e) = tmp
    end do
    end do
    end do

    ! z = (d₀IxIxI + d₁IxIxV₁ + d₂IxIxV₂ + d₃IxIxV₃)⁻¹ (I x I x S^t) u .........

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + S(m,i,c1) * u(m,j,k,e)
      end do
      z(i,j,k) = tmp / ( lambda *   g0             &
                       + nu(e)  * ( g1 * V1(i,c1)  &
                                  + g2 * V2(j,c2)  &
                                  + g3 * V3(k,c3)  &
                                  )                &
                       )
    end do
    end do
    end do

    ! u = I x I x WS z .........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + WS_t(m,i,c1) * z(m,j,k)
      end do
      u(i,j,k,e) = tmp
    end do
    end do
    end do

    ! z = I x WS x I u .........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + WS_t(m,j,c2) * u(i,m,k,e)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    ! u = WS x I x I z .........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + WS_t(m,k,c3) * z(i,j,m)
      end do
      u(i,j,k,e) = tmp
    end do
    end do
    end do

  end do elements
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine CART__TPO_Schwarz_Old__gen

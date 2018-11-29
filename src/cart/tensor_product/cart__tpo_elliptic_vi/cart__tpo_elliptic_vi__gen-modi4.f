!> summary:  Elliptic element operator with variable diffusivity: generic
!> author:   Karl Schoppmann, Joerg Stiller
!> date:     2018/11/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Elliptic element operator with variable diffusivity: generic
!>
!>
!> @note
!> When compiled with OpenACC and run on a GPU, any call to this
!> procedure must be followed by an "acc wait" before accessing u
!> or v from CPU or any other than the default accelerator queue.
!> @endnote
!===============================================================================

subroutine CART__TPO_Elliptic_VI__gen(np, ne, Ms, Ds, lambda, nu, dx, u, v)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters,  only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: np              !< number of points per direction
  integer,   intent(in)  :: ne              !< number of elements
  real(RNP), intent(in)  :: Ms(np)          !< 1D standard mass matrix
  real(RNP), intent(in)  :: Ds(np,np)       !< 1D standard diff matrix
  real(RNP), intent(in)  :: lambda          !< Helmholtz parameter
  real(RNP), intent(in)  :: nu(np,np,np,ne) !< diffusivity
  real(RNP), intent(in)  :: dx(3)           !< element extensions
  real(RNP), intent(in)  :: u(np,np,np,ne)  !< operand
  real(RNP), intent(out) :: v(np,np,np,ne)  !< result


!######################### REWRITE STARTING FROM HERE ##########################

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), allocatable :: M(:,:,:), M_u(:,:,:), Ms_Ds(:,:), Msi_Dst(:,:), z(:)

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
  allocate(M(np,np,np), M_u(np,np,np), M_ujik(np,np,np), &
           M_ukij(np,np,np), Ms_Ds(np,np), Msi_Dst(np,np), z(np))

  ! element mass matrix
  tmp = product(dx) / 8
  do k = 1, np
  do j = 1, np
  do i = 1, np
    M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  ! modified differentiation operators
  !! to optimize cache use, the operators are used in transposed form
  !! to save loads, both operators are computed in one loop-nest
  do j = 1, np
  do i = 1, np
    Ms_Ds   (i,j) = Ms(i) * Ds(i,j)  ! = Ms Ds     =  [ (Ms Ds)ᵀ ]ᵀ
    Msi_Dst (i,j) = Ds(j,i) / Ms(i)  ! = Ms⁻¹ Dsᵀ  =  [  Ds Ms⁻¹ ]ᵀ
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

    !$acc loop collapse(2) independent vector
    do k = 1, np
    do j = 1, np

      do i = 1, np
        z(i) = 0
        do p = 1, np
          z(i) = z(i) + Msi_Dst(p,i) * M_u(p,j,k)
        end do
        z(i) = nu(i,j,k,e) * z(i)
      end do

      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ms_Ds(p,i) * z(p)
        end do
        v(i,j,k,e) = v(i,j,k,e) + g(1) * tmp
        M_ujik(i,j,k) = M_u(j,i,k)
      end do

    end do
    end do

    ! direction 2 ..............................................................
    
    do k = 1, np
    do j = 1, np

      do i = 1, np
        z(i) = 0
        do p = 1, np
          z(i) = z(i) + Msi_Dst(p,i) * M_ujik(p,j,k)
        end do
        z(i) = nu(j,i,k,e) * z(i)
      end do

      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ms_Ds(p,i) * z(p)
        end do
        v(j,i,k,e) = v(j,i,k,e) + g(2) * tmp
        M_ukij(i,j,k) = M_u(k,i,j)
      end do

    end do
    end do

    ! direction 3 ..............................................................
    
    do k = 1, np
    do j = 1, np

      do i = 1, np
        z(i) = 0
        do p = 1, np
          z(i) = z(i) + Msi_Dst(p,i) * M_ukij(p,j,k)
        end do
        z(i) = nu(j,k,i,e) * z(i)
      end do

      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ms_Ds(p,i) * z(p)
        end do
        v(j,k,i,e) = v(j,k,i,e) + g(3) * tmp
      end do

    end do
    end do

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

!######################### REWRITE STARTING UNTIL HERE #########################

!===============================================================================

end subroutine CART__TPO_Elliptic_VI__gen

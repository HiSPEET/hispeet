!> summary:  3d generic isotropic Schwarz operator
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Schwarz__3D_CI__Gen
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_Schwarz_CI_Gen

contains

  !-----------------------------------------------------------------------------
  !> 3d generic isotropic Schwarz operator

  subroutine TPO_Schwarz_CI_Gen(S, W, cfg, D_inv, f, u)
    real(RNP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors   (np,np,nc)
    real(RNP), intent(in)  :: W(:,:)         !< 1D weights        (np,nc)
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configs (3,nd)
    real(RNP), intent(in)  :: D_inv(:,:,:,:) !< inv. 3D eigenvals (np,np,np,nd)
    real(RNP), intent(in)  :: f(:,:,:,:)     !< RHS               (np,np,np,nd)
    real(RNP), intent(out) :: u(:,:,:,:)     !< solution          (np,np,np,nd)

    !-----------------------------------------------------------------------------
    ! local variables

    real(RNP) :: WS_t ( size(S,1), size(S,2), size(S,3) )
    real(RNP) :: y    ( size(f,1), size(f,2), size(f,3) )
    real(RNP) :: z    ( size(f,1), size(f,2), size(f,3) )
    real(RNP) :: tmp

    integer :: nc, nd, np
    integer :: i, j, k, l, p
    integer :: c, c1, c2, c3

    !-----------------------------------------------------------------------------
    ! initialization

    np = size(S,1)
    nc = size(S,3)
    nd = size(f,4)

    ! WS_t(:,:,c) = (WS)ᵀ(:,:,c)
    do c = 1, nc
    do j = 1, np
    do i = 1, np
      WS_t(j,i,c) = W(i,c) * S(i,j,c)
    end do
    end do
    end do

    !---------------------------------------------------------------------------
    ! evaluation

    !$acc data async present(cfg, D_inv, f, u) copyin(S, WS_t)
    !$acc parallel
    !$acc loop gang worker private(y,z)

    !$omp do private(l)
    subdomains: do l = 1, nd

      ! element-boundary configurations
      c1 = cfg(1,l)
      c2 = cfg(2,l)
      c3 = cfg(3,l)

      ! y = Sᵀ x I x I f .......................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + S(p,k,c3) * f(i,j,p,l)
          end do
          y(i,j,k) = tmp
        end do
      end do
      end do

      ! z = I x Sᵀ x I y .......................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + S(p,j,c2) * y(i,p,k)
          end do
          z(i,j,k) = tmp
        end do
      end do
      end do

      ! y = D⁻¹ (I x I x Sᵀ) z .................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + S(p,i,c1) * z(p,j,k)
          end do
          y(i,j,k) = tmp * D_inv(i,j,k,l)
        end do
      end do
      end do

      ! z = I x I x WS y .......................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + WS_t(p,i,c1) * y(p,j,k)
          end do
          z(i,j,k) = tmp
        end do
      end do
      end do

      ! y = I x WS x I z .......................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + WS_t(p,j,c2) * z(i,p,k)
          end do
          y(i,j,k) = tmp
        end do
      end do
      end do

      ! u = WS x I x I y .......................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + WS_t(p,k,c3) * y(i,j,p)
          end do
          u(i,j,k,l) = tmp
        end do
      end do
      end do

    end do subdomains

    !$acc end parallel
    !$acc end data

  end subroutine TPO_Schwarz_CI_Gen

  !=============================================================================

end module TPO__Schwarz__3D_CI__Gen

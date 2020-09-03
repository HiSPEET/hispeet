!> summary:  3d generic isotropic spectral operator
!> author:   Joerg Stiller
!> date:     2020/05/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral_3d_CI__Gen
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_Spectral_CI_Gen

contains

  !-----------------------------------------------------------------------------
  !> Generic procedure for evaluating  v = [(SxSxS) Λ (SxSxS)ᵀ] u

  subroutine TPO_Spectral_CI_Gen(S, Lambda, u, v)
    real(RNP), intent(in)  :: S(:,:)        !< 1D eigenvectors
    real(RNP), intent(in)  :: Lambda(:,:,:) !< 3D eigenvalues
    real(RNP), intent(in)  :: u(:,:,:,:)    !< operand
    real(RNP), intent(out) :: v(:,:,:,:)    !< result

    !---------------------------------------------------------------------------
    ! local variables

    real(RNP) :: St( size(S,1), size(S,1) )
    real(RNP) :: z ( size(S,1), size(S,1), size(S,1) )
    real(RNP) :: tmp
    integer   :: np, ne
    integer   :: e, i, j, k, p

    !---------------------------------------------------------------------------
    ! initialization

    np = size(S,1)
    ne = size(u,4)

    St = transpose(S)

    !---------------------------------------------------------------------------
    ! evaluation

    !$acc data present(u,v) copyin(S,St,Lambda)
    !$acc parallel
    !$acc loop gang worker private(z)

    !$omp do private(e)
    elements: do e = 1, ne

      ! z = Sᵀ x I x I u .......................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + S(p,k) * u(i,j,p,e)
          end do
          z(i,j,k) = tmp
        end do
      end do
      end do

      ! v = I x Sᵀ x I z .......................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + S(p,j) * z(i,p,k)
          end do
          v(i,j,k,e) = tmp
        end do
      end do
      end do

      ! z = Lambda I x I x Sᵀ v ................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + S(p,i) * v(p,j,k,e)
          end do
          z(i,j,k) = tmp * Lambda(i,j,k)
        end do
      end do
      end do

      ! v = I x I x S z ........................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + St(p,i) * z(p,j,k)
          end do
          v(i,j,k,e) = tmp
        end do
      end do
      end do

      ! z = I x S x I v ........................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + St(p,j) * v(i,p,k,e)
          end do
          z(i,j,k) = tmp
        end do
      end do
      end do

      ! v = S x I x I z ........................................................

      !$acc loop collapse(3) vector
      do k = 1, np
      do j = 1, np
        !DIR$ SIMD
        do i = 1, np
          tmp = 0
          do p = 1, np
            tmp = tmp + St(p,k) * z(i,j,p)
          end do
          v(i,j,k,e) = tmp
        end do
      end do
      end do

    end do elements
    !$omp end do

    !$acc end parallel
    !$acc end data

  end subroutine TPO_Spectral_CI_Gen

  !=============================================================================

end module TPO__Spectral_3d_CI__Gen


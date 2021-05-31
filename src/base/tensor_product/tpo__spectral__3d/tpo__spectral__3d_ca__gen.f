!> summary:  3d generic anisotropic spectral operator
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral__3D_CA__Gen
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_Spectral_CA_Gen

contains

  !-----------------------------------------------------------------------------
  !> Generic procedure for evaluating  v = [(S₃ x S₂ x S₁) Λ (S₃ x S₂ x S₁)ᵀ] u

  subroutine TPO_Spectral_CA_Gen(S1, S2, S3, Lambda, u, v)
    real(RNP), intent(in)  :: S1(:,:)       !< eigenvectors for dir 1: S₁(n1,n1)
    real(RNP), intent(in)  :: S2(:,:)       !< eigenvectors for dir 2: S₂(n2,n2)
    real(RNP), intent(in)  :: S3(:,:)       !< eigenvectors for dir 3: S₃(n3,n3)
    real(RNP), intent(in)  :: Lambda(:,:,:) !< 3D eigenvalues: Λ(n1,n2,n3)
    real(RNP), intent(in)  :: u(:,:,:,:)    !< operand: u(n1.n2,n3,ne)
    real(RNP), intent(out) :: v(:,:,:,:)    !< result:  v(n1,n2,n3,ne)

    !-----------------------------------------------------------------------------
    ! local variables

    real(RNP) :: z( size(S1,1), size(S2,1), size(S3,1) )
    real(RNP) :: tmp
    integer   :: n1, n2, n3, ne
    integer   :: e, i, j, k, p

    !---------------------------------------------------------------------------
    ! initialization

    n1 = size(S1,1)
    n2 = size(S2,1)
    n3 = size(S3,1)
    ne = size(u ,4)

    !---------------------------------------------------------------------------
    ! evaluation

    !$acc data present(u,v) copyin(S1,S2,S3,Lambda)
    !$acc parallel
    !$acc loop gang worker private(z1,z2,z3)

    !$omp do private(e)
    elements: do e = 1, ne

      ! z = S₃ᵀ x I x I u ......................................................

      !$acc loop collapse(3) vector
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n3
            tmp = tmp + S3(p,k) * u(i,j,p,e)
          end do
          z(i,j,k) = tmp
        end do
      end do
      end do

      ! vᵉ = I x S₂ᵀ x I z .....................................................

      !$acc loop vector collapse(3)
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n2
            tmp = tmp + S2(p,j) * z(i,p,k)
          end do
          v(i,j,k,e) = tmp
        end do
      end do
      end do

      ! z = Λ (I x I x S₁ᵀ) vᵉ .................................................

      !$acc loop vector collapse(3)
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n1
            tmp = tmp + S1(p,i) * v(p,j,k,e)
          end do
          z(i,j,k) = tmp * Lambda(i,j,k)
        end do
      end do
      end do

      ! vᵉ = I x I x S₁ z ......................................................

      !$acc loop vector collapse(3)
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n1
            tmp = tmp + S1(i,p) * z(p,j,k)
          end do
          v(i,j,k,e) = tmp
        end do
      end do
      end do

      ! z = I x S₂ x I vᵉ ......................................................

      !$acc loop vector collapse(3)
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n2
            tmp = tmp + S2(j,p) * v(i,p,k,e)
          end do
          z(i,j,k) = tmp
        end do
      end do
      end do

      ! vᵉ = S₃ x I x I z ......................................................

      !$acc loop collapse(3) vector
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n3
            tmp = tmp + S3(k,p) * z(i,j,p)
          end do
          v(i,j,k,e) = tmp
        end do
      end do
      end do

    end do elements

    !$acc end parallel
    !$acc end data

  end subroutine TPO_Spectral_CA_Gen

  !=============================================================================

end module TPO__Spectral__3D_CA__Gen

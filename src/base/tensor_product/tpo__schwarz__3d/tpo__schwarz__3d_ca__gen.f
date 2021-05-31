!> summary:  3d generic anisotropic Schwarz operator
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Schwarz__3D_CA__Gen
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_Schwarz_CA_Gen

contains

  !-----------------------------------------------------------------------------
  !> 3d generic anisotropic Schwarz operator

  subroutine TPO_Schwarz_CA_Gen(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, u)
    real(RNP), intent(in)  :: S1(:,:,:)      !< eigenvecs, dir 1   (n1,n1,nc)
    real(RNP), intent(in)  :: S2(:,:,:)      !< eigenvecs, dir 2   (n2,n2,nc)
    real(RNP), intent(in)  :: S3(:,:,:)      !< eigenvecs, dir 3   (n3,n3,nc)
    real(RNP), intent(in)  :: W1(:,:)        !< weights,   dir 1   (n1,nc)
    real(RNP), intent(in)  :: W2(:,:)        !< weights,   dir 2   (n2,nc)
    real(RNP), intent(in)  :: W3(:,:)        !< weights,   dir 3   (n3,nc)
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configs  (3,nd)
    real(RNP), intent(in)  :: D_inv(:,:,:,:) !< inv 3D eigenvals   (n1,n2,n3,nd)
    real(RNP), intent(in)  :: f(:,:,:,:)     !< RHS                (n1,n2,n3,nd)
    real(RNP), intent(out) :: u(:,:,:,:)     !< solution           (n1,n2,n3,nd)

    !---------------------------------------------------------------------------
    ! local variables

    real(RNP) :: WS1_t ( size(S1,1), size(S1,2), size(S1,3) )
    real(RNP) :: WS2_t ( size(S2,1), size(S2,2), size(S2,3) )
    real(RNP) :: WS3_t ( size(S3,1), size(S3,2), size(S3,3) )
    real(RNP) :: y     ( size(f ,1), size(f ,2), size(f ,3) )
    real(RNP) :: z     ( size(f ,1), size(f ,2), size(f ,3) )
    real(RNP) :: tmp

    integer :: n1, n2, n3, nc, nd
    integer :: i, j, k, l, p
    integer :: c, c1, c2, c3

    !-----------------------------------------------------------------------------
    ! initialization

    n1 = size(S1, 1)
    n2 = size(S2, 1)
    n3 = size(S3, 1)
    nc = size(S3, 3)
    nd = size(f , 4)

    ! WS1_t(:,:,c) = (W₁S₁)ᵀ(:,:,c)
    do c = 1, nc
    do j = 1, n1
    do i = 1, n1
      WS1_t(j,i,c) = W1(i,c) * S1(i,j,c)
    end do
    end do
    end do

    ! WS2_t(:,:,c) = (W₂S₂)ᵀ(:,:,c)
    do c = 1, nc
    do j = 1, n2
    do i = 1, n2
      WS2_t(j,i,c) = W2(i,c) * S2(i,j,c)
    end do
    end do
    end do

    ! WS3_t(:,:,c) = (W₃S₃)ᵀ(:,:,c)
    do c = 1, nc
    do j = 1, n3
    do i = 1, n3
      WS3_t(j,i,c) = W3(i,c) * S3(i,j,c)
    end do
    end do
    end do

    !---------------------------------------------------------------------------
    ! evaluation

    !$acc data async &
    !$acc & present(cfg, D_inv, f, u) &
    !$acc & copyin(S1, S2, S3, WS1_t, WS2_t, WS3_t)

    !$acc parallel async

    !$acc loop gang worker private(y,z)

    !$omp do private(l)
    subdomains: do l = 1, nd

      ! element-boundary configurations
      c1 = cfg(1,l)
      c2 = cfg(2,l)
      c3 = cfg(3,l)

      ! y = S₃ᵀ x I x I f ......................................................

      !$acc loop collapse(3) vector
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n3
            tmp = tmp + S3(p,k,c3) * f(i,j,p,l)
          end do
          y(i,j,k) = tmp
        end do
      end do
      end do

      ! z = I x S₂ᵀ x I y ......................................................

      !$acc loop vector collapse(3)
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n2
            tmp = tmp + S2(p,j,c2) * y(i,p,k)
          end do
          z(i,j,k) = tmp
        end do
      end do
      end do

      ! y = D⁻¹ (I x I x S₁ᵀ) z ................................................

      !$acc loop vector collapse(3)
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n1
            tmp = tmp + S1(p,i,c1) * z(p,j,k)
          end do
          y(i,j,k) = tmp * D_inv(i,j,k,l)
        end do
      end do
      end do

      ! z = I x I x W₁S₁ y .....................................................

      !$acc loop vector collapse(3)
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n1
            tmp = tmp + WS1_t(p,i,c1) * y(p,j,k)
          end do
          z(i,j,k) = tmp
        end do
      end do
      end do

      ! y = I x W₂S₂ x I z .....................................................

      !$acc loop vector collapse(3)
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n2
            tmp = tmp + WS2_t(p,j,c2) * z(i,p,k)
          end do
          y(i,j,k) = tmp
        end do
      end do
      end do

      ! u = W₃S₃ x I x I y .....................................................

      !$acc loop collapse(3) vector
      do k = 1, n3
      do j = 1, n2
        !DIR$ SIMD
        do i = 1, n1
          tmp = 0
          do p = 1, n3
            tmp = tmp + WS3_t(p,k,c3) * y(i,j,p)
          end do
          u(i,j,k,l) = tmp
        end do
      end do
      end do

    end do subdomains
    !$omp end do

    !$acc end parallel
    !$acc end data

  end subroutine TPO_Schwarz_CA_Gen

  !=============================================================================

end module TPO__Schwarz__3D_CA__Gen

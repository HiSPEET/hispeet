!> summary:   3d generic AxAxA operator
!> author:    Jörg Stiller
!> date:      2020/05/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__AAA__3D__Gen
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_AAA_Gen

contains

  !-----------------------------------------------------------------------------
  !> Generic AxAxA operator

  subroutine TPO_AAA_Gen(A, u, v)
    real(RNP), intent(in)  :: A(:,:)     !< 1D operator
    real(RNP), intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), intent(out) :: v(:,:,:,:) !< result

    !---------------------------------------------------------------------------
    ! local variables

    real(RNP) :: At( size(A,2), size(A,1) )
    real(RNP) :: z2( size(A,2), size(A,1), size(A,1))
    real(RNP) :: z3( size(A,2), size(A,2), size(A,1))
    real(RNP) :: tmp
    integer   :: na1, na2, ne
    integer   :: e, i, j, k, p

    !---------------------------------------------------------------------------
    ! initialization

    na1 = size(A,1)
    na2 = size(A,2)
    ne  = size(u,4)

    At = transpose(A)

    !---------------------------------------------------------------------------
    ! evaluation

    !$acc data present(u,v) copyin(At)
    !$acc parallel
    !$acc loop gang worker private(z2,z3)

    !$omp do private(e)
    do e = 1, ne

      ! z3 = AxIxI u^e .........................................................

      !$acc loop collapse(3) vector
      do k = 1, na1
      do j = 1, na2
      !DIR$ SIMD
      do i = 1, na2
        tmp = 0
        do p = 1, na2
          tmp = tmp + At(p,k) * u(i,j,p,e)
        end do
        z3(i,j,k) = tmp
      end do
      end do
      end do

      ! z2 = IxAxI z3 ..........................................................

      !$acc loop collapse(3) vector
      do k = 1, na1
      do j = 1, na1
      !DIR$ SIMD
      do i = 1, na2
        tmp = 0
        do p = 1, na2
          tmp = tmp + At(p,j) * z3(i,p,k)
        end do
        z2(i,j,k) = tmp
      end do
      end do
      end do

      ! v^e = IxIxA z2 .........................................................

      !$acc loop collapse(3) vector
      do k = 1, na1
      do j = 1, na1
      !DIR$ SIMD
      do i = 1, na1
        tmp = 0
        do p = 1, na2
          tmp = tmp + At(p,i) * z2(p,j,k)
        end do
        v(i,j,k,e) = tmp
      end do
      end do
      end do

    end do
    !$omp end do

    !$acc end parallel
    !$acc end data

  end subroutine TPO_AAA_Gen

  !=============================================================================

end module TPO__AAA__3D__Gen

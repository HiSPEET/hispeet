!> summary:   3d regular 8:1 projection operator using LIBXSMM for real(RDP)
!> author:    Joerg Stiller
!> date:      2024/08/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__8To1__3D__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_8To1_XSMM

  interface TPO_8To1_XSMM
    module procedure TPO_8To1_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> XSMM AxAxA operator

  subroutine TPO_8To1_XSMM_RDP(A, B, smooth, u, v)
    real(RDP), contiguous, intent(in)  :: A(:,:,:)   !< 1D projection operator
    real(RDP), contiguous, intent(in)  :: B(:,:)     !< 1D smoothing operator
    integer,               intent(in)  :: smooth     !< discontinuity handling
    real(RDP), contiguous, intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), contiguous, intent(out) :: v(:,:,:,:) !< result


    real(RDP), parameter :: alpha = 1
    real(RDP), parameter :: beta  = 1
    real(RDP), parameter :: HALF  = 0.5

    !---------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP) :: At   ( size(A,2), size(A,1), 2 )
    real(RDP) :: uc   ( size(u,1), size(u,1), size(u,1), 2, 2, 2 )

    real(RDP) :: z3   ( size(u,1), size(u,1), size(v,1), 2, 2 )
    real(RDP) :: jmp3 ( size(u,1), size(u,1))

    real(RDP) :: z2   ( size(u,1), size(v,1), size(v,1), 2 )
    real(RDP) :: y2   ( size(u,1), size(A,1), 2 )
    real(RDP) :: jmp2 ( size(u,1))

    real(RDP) :: y1   ( size(A,1), size(v,1), size(v,1), 2 )
    real(RDP) :: jmp1

    integer   :: na1, na2, npu, npv
    integer   :: e, i, j, k, l, m, o, p, q1, q2

    !---------------------------------------------------------------------------
    ! initialization

    na1 = size(A,1)
    na2 = size(A,2)

    npu = size(u,1) ! = na2
    npv = size(v,1)

    At(:,:,1) = transpose(A(:,:,1))
    At(:,:,2) = transpose(A(:,:,2))

    ! point offsetsfor contribution of first and second fine element
    q1 = 0
    q2 = npv - na1

    !-----------------------------------------------------------------------------
    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_3, npu**2, na1   , na2, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, npu   , na1   , na2, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_1, na1   , npv**2, na2, alpha=alpha, beta=beta)
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, size(u,4)

      ! z3 = AxIxI u ...........................................................

      ! offset of fine element index
      o = 8 * (e-1)

      ! extract 2x2x2 cluster of fine element variable
      uc(:,:,:,1,1,1) = u(:,:,:,o+1)
      uc(:,:,:,2,1,1) = u(:,:,:,o+2)
      uc(:,:,:,1,2,1) = u(:,:,:,o+3)
      uc(:,:,:,2,2,1) = u(:,:,:,o+4)
      uc(:,:,:,1,1,2) = u(:,:,:,o+5)
      uc(:,:,:,2,1,2) = u(:,:,:,o+6)
      uc(:,:,:,1,2,2) = u(:,:,:,o+7)
      uc(:,:,:,2,2,2) = u(:,:,:,o+8)

      ! discontinuity treatment in direction 3
      select case(smooth)
      case(1)
        ! linear compensation
        do m = 1, 2
        do l = 1, 2
          jmp3 = uc(:,:,npu,l,m,1) - uc(:,:,1,l,m,2)
          do k = 1, npu
            uc(:,:,k,l,m,1) = uc(:,:,k,l,m,1) + jmp3 * B(k,1)
            uc(:,:,k,l,m,2) = uc(:,:,k,l,m,2) + jmp3 * B(k,2)
          end do
        end do
        end do
      case(2)
        ! averaging interface coefficients
        do m = 1, 2
        do l = 1, 2
          jmp3 = uc(:,:,npu,l,m,1) - uc(:,:,1,l,m,2)
          uc(:,:,npu,l,m,1) = uc(:,:,npu,l,m,1) - HALF * jmp3
          uc(:,:,  1,l,m,2) = uc(:,:,  1,l,m,2) + HALF * jmp3
        end do
        end do
      end select

      ! projection in direction 3
      do m = 1, 2
      do l = 1, 2
        z3(:,:,:,l,m) = 0
        call LIBXSMM_DMMCall(xmm_3, uc(:,:,:,l,m,1), At(:,:,1), z3(:,:,q1+1:,l,m))
        call LIBXSMM_DMMCall(xmm_3, uc(:,:,:,l,m,2), At(:,:,2), z3(:,:,q2+1:,l,m))
      end do
      end do

      ! z2 = IxAxI z3 ..........................................................

      ! discontinuity treatment in direction 2
      select case(smooth)
      case(1)
        ! linear compensation
        do k = 1, npv
        do l = 1, 2
          jmp2 = z3(:,npu,k,l,1) - z3(:,1,k,l,2)
          do j = 1, npu
            z3(:,j,k,l,1) = z3(:,j,k,l,1) + jmp2 * B(j,1)
            z3(:,j,k,l,2) = z3(:,j,k,l,2) + jmp2 * B(j,2)
          end do
        end do
        end do
      case(2)
        ! averaging interface coefficients
        do k = 1, npv
        do l = 1, 2
          jmp2 = z3(:,npu,k,l,1) - z3(:,1,k,l,2)
          z3(:,npu,k,l,1) = z3(:,npu,k,l,1) - HALF * jmp2
          z3(:,  1,k,l,2) = z3(:,  1,k,l,2) + HALF * jmp2
        end do
        end do
      end select

      ! projection in direction 2
      do l = 1, 2
      do k = 1, npv
        ! use y2 for temporary storage as LIBXSMM expects dense memory operands
        call LIBXSMM_DMMCall(xmm_2, z3(:,:,k,l,1), At(:,:,1), y2(:,:,1))
        call LIBXSMM_DMMCall(xmm_2, z3(:,:,k,l,2), At(:,:,2), y2(:,:,2))
        z2(:,:,k,l) = 0
        do p = 1, na1
          z2(:,q1+p,k,l) = z2(:,q1+p,k,l) + y2(:,p,1)
          z2(:,q2+p,k,l) = z2(:,q2+p,k,l) + y2(:,p,2)
        end do
      end do
      end do

      ! v^e = IxIxA z2 ............................................................

      ! discontinuity treatment in direction 1
      select case(smooth)
      case(1)
        ! linear compensation
        do k = 1, npv
        do j = 1, npv
          jmp1 = z2(npu,j,k,1) - z2(1,j,k,2)
          do i = 1, npu
            z2(q1+i,j,k,1) = z2(q1+i,j,k,1) + jmp1 * B(i,1)
            z2(q2+i,j,k,2) = z2(q2+i,j,k,2) + jmp1 * B(i,2)
          end do
        end do
        end do
      case(2)
        ! averaging interface coefficients
        do k = 1, npv
        do j = 1, npv
          jmp1 = z2(npu,j,k,1) - z2(1,j,k,2)
          z2(npu,j,k,1) = z2(npu,j,k,1) - HALF * jmp1
          z2(  1,j,k,2) = z2(  1,j,k,2) + HALF * jmp1
        end do
        end do
      end select

      ! projection in direction 1
      call LIBXSMM_DMMCall(xmm_1, A(:,:,1), z2(:,:,:,1), y1(:,:,:,1))
      call LIBXSMM_DMMCall(xmm_1, A(:,:,2), z2(:,:,:,2), y1(:,:,:,2))
      do k = 1, npv
      do j = 1, npv
        v(:,j,k,e) = 0
        do p = 1, na1
          v(q1+p,j,k,e) = v(q1+p,j,k,e) + y1(p,j,k,1)
          v(q2+p,j,k,e) = v(q2+p,j,k,e) + y1(p,j,k,2)
        end do
      end do
      end do

    end do

  end subroutine TPO_8To1_XSMM_RDP

  !=============================================================================

#endif

end module TPO__8To1__3D__XSMM_RDP

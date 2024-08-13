!> summary:   3d generic projection operator for regular 8:1 h-coarsening
!> author:    Jörg Stiller
!> date:      2024/08/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> Discontinuities between the fine elements are treated according to the
!> `smooth` parameter:
!>   - `0`  no discontinuity handling
!>   - `1`  jump removal using antisymmetric linear correction
!>   - `2`  jump removal by averaging interface coefficients
!>
!> This variant runs through the directions in the order 3,2,1
!===============================================================================

!-------------------------------------------------------------------------------
!> Generic 3d projection operator for regular 8:1 h-coarsening

subroutine TPO_8To1_Gen_RWP(A, B, smooth, u, v)
  real(RWP), contiguous, intent(in)  :: A(:,:,:)   !< 1D projection operator
  real(RWP), contiguous, intent(in)  :: B(:,:)     !< 1D smoothing operator
  integer,               intent(in)  :: smooth     !< discontinuity handling
  real(RWP), contiguous, intent(in)  :: u(:,:,:,:) !< operand
  real(RWP), contiguous, intent(out) :: v(:,:,:,:) !< result

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP), parameter :: HALF = 0.5

  real(RWP) :: At   ( size(A,2), size(A,1), 2 )
  real(RWP) :: w    ( size(u,1), size(u,1), size(u,1), 2, 2, 2 )
  real(RWP) :: z3   ( size(u,1), size(u,1), size(v,1), 2, 2 )
  real(RWP) :: z2   ( size(u,1), size(v,1), size(v,1), 2 )
  real(RWP) :: jmp3 ( size(u,1), size(u,1))
  real(RWP) :: jmp2 ( size(u,1))
  real(RWP) :: jmp1

  integer   :: na1, na2, npu, npv
  integer   :: e, i, j, k, l, m, o, p, q1, q2

  !-----------------------------------------------------------------------------
  ! initialization

  na1 = size(A,1)
  na2 = size(A,2)

  npu = size(u,1) ! = na2
  npv = size(v,1)

  At(:,:,1) = transpose(A(:,:,1))
  At(:,:,2) = transpose(A(:,:,2))

  ! point offsets for contribution of first and second fine element
  q1 = 0
  q2 = npv - na1

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e,o)
  do e = 1, size(v,4)

    ! z3 = AxIxI u .............................................................

    ! offset of fine element index
    o = 8 * (e-1)

    ! extract fine element variable
    w(:,:,:,1,1,1) = u(:,:,:,o+1)
    w(:,:,:,2,1,1) = u(:,:,:,o+2)
    w(:,:,:,1,2,1) = u(:,:,:,o+3)
    w(:,:,:,2,2,1) = u(:,:,:,o+4)
    w(:,:,:,1,1,2) = u(:,:,:,o+5)
    w(:,:,:,2,1,2) = u(:,:,:,o+6)
    w(:,:,:,1,2,2) = u(:,:,:,o+7)
    w(:,:,:,2,2,2) = u(:,:,:,o+8)

    ! discontinuity treatment in direction 3
    select case(smooth)
    case(1)
      ! linear compensation
      do m = 1, 2
      do l = 1, 2
        jmp3 = w(:,:,npu,l,m,1) - w(:,:,1,l,m,2)
        do k = 1, npu
          w(:,:,k,l,m,1) = w(:,:,k,l,m,1) + jmp3 * B(k,1)
          w(:,:,k,l,m,2) = w(:,:,k,l,m,2) + jmp3 * B(k,2)
        end do
      end do
      end do
    case(2)
      ! averaging interface coefficients
      do m = 1, 2
      do l = 1, 2
        jmp3 = w(:,:,npu,l,m,1) - w(:,:,1,l,m,2)
        w(:,:,npu,l,m,1) = w(:,:,npu,l,m,1) - HALF * jmp3
        w(:,:,  1,l,m,2) = w(:,:,  1,l,m,2) + HALF * jmp3
      end do
      end do
    end select

    ! projection in direction 3
    z3 = 0
    do k = 1, na1
    do j = 1, npu
    do i = 1, npu
      do p = 1, na2
        do m = 1, 2
        do l = 1, 2
          z3(i,j,q1+k,l,m) = z3(i,j,q1+k,l,m) + w(i,j,p,l,m,1) * At(p,k,1)
          z3(i,j,q2+k,l,m) = z3(i,j,q2+k,l,m) + w(i,j,p,l,m,2) * At(p,k,2)
        end do
        end do
      end do
    end do
    end do
    end do

    ! z2 = IxAxI z3 ............................................................

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
    z2 = 0
    do k = 1, npv
    do j = 1, na1
    do i = 1, npu
      do p = 1, na2
        do l = 1, 2
          z2(i,q1+j,k,l) = z2(i,q1+j,k,l) + z3(i,p,k,l,1) * At(p,j,1)
          z2(i,q2+j,k,l) = z2(i,q2+j,k,l) + z3(i,p,k,l,2) * At(p,j,2)
        end do
      end do
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
    v(:,:,:,e) = 0
    do k = 1, npv
    do j = 1, npv
    do i = 1, na1
      do p = 1, na2
        v(q1+i,j,k,e) = v(q1+i,j,k,e) + z2(p,j,k,1) * At(p,i,1)
        v(q2+i,j,k,e) = v(q2+i,j,k,e) + z2(p,j,k,2) * At(p,i,2)
      end do
    end do
    end do
    end do

  end do
  !$omp end do

end subroutine TPO_8To1_Gen_RWP

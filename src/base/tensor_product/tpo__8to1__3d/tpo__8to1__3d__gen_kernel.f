!> summary:   3d generic projection operator for regular 8:1 h-coarsening
!> author:    Jörg Stiller
!> date:      2024/08/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

!-------------------------------------------------------------------------------
!> Generic 3d projection operator for regular 8:1 h-coarsening

subroutine TPO_8To1_Gen_RWP(A, u, v)
  real(RWP), contiguous, intent(in)  :: A(:,:,:)   !< 1D projection operator
  real(RWP), contiguous, intent(in)  :: u(:,:,:,:) !< operand
  real(RWP), contiguous, intent(out) :: v(:,:,:,:) !< result

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP), parameter :: HALF = 0.5

  real(RWP) :: At ( size(A,2), size(A,1), 2 )
  real(RWP) :: w  ( size(u,1), size(u,1), size(u,1), 2, 2, 2 )
  real(RWP) :: z3 ( size(u,1), size(u,1), size(v,1), 2, 2 )
  real(RWP) :: z2 ( size(u,1), size(v,1), size(v,1), 2 )

  integer   :: np, nq
  integer   :: e, i, j, k, l, m, o, p

  !-----------------------------------------------------------------------------
  ! initialization

  np = size(u,1) ! = size(A,2)
  nq = size(v,1) ! = size(A,1)

  At(:,:,1) = transpose(A(:,:,1))
  At(:,:,2) = transpose(A(:,:,2))

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

    ! projection in direction 3
    z3 = 0
    do k = 1, nq
    do j = 1, np
    do i = 1, np
      do p = 1, np
        do m = 1, 2
        do l = 1, 2
          z3(i,j,k,l,m) = z3(i,j,k,l,m) + w(i,j,p,l,m,1) * At(p,k,1) &
                                        + w(i,j,p,l,m,2) * At(p,k,2)
        end do
        end do
      end do
    end do
    end do
    end do

    ! z2 = IxAxI z3 ............................................................

    ! projection in direction 2
    z2 = 0
    do k = 1, nq
    do j = 1, nq
    do i = 1, np
      do p = 1, np
        do l = 1, 2
          z2(i,j,k,l) = z2(i,j,k,l) + z3(i,p,k,l,1) * At(p,j,1) &
                                    + z3(i,p,k,l,2) * At(p,j,2)
        end do
      end do
    end do
    end do
    end do

    ! v^e = IxIxA z2 ............................................................

    ! projection in direction 1
    v(:,:,:,e) = 0
    do k = 1, nq
    do j = 1, nq
    do i = 1, nq
      do p = 1, np
        v(i,j,k,e) = v(i,j,k,e) + z2(p,j,k,1) * At(p,i,1) &
                                + z2(p,j,k,2) * At(p,i,2)
      end do
    end do
    end do
    end do

  end do
  !$omp end do

end subroutine TPO_8To1_Gen_RWP

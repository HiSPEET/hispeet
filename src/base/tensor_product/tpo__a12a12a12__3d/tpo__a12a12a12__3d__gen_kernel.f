!> summary:   3d generic 1:2 interpolation operator
!> author:    Jörg Stiller
!> date:      2024/07/24
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

!-----------------------------------------------------------------------------
!> Generic 3d 1:2 interpolation operator

subroutine TPO_A12A12A12_Gen_RWP(A12, u, v)
  real(RWP), intent(in)  :: A12(:,:,:) !< 1D operator (na1,na2,2)
  real(RWP), intent(in)  :: u(:,:,:,:) !< operand
  real(RWP), intent(out) :: v(:,:,:,:) !< result

  !---------------------------------------------------------------------------
  ! local variables

  real(RWP) :: z2( size(A,2), size(A,1), size(A,1), 2, 2 )
  real(RWP) :: z3( size(A,2), size(A,2), size(A,1), 2 )
  integer   :: na1, na2, nec
  integer   :: e, i, j, k, l, p

  !---------------------------------------------------------------------------
  ! initialization

  na1 = size(A,1)
  na2 = size(A,2)
  nec = size(u,4) ! number of coarse elements / fine element clusters

  !---------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e,l)
  do e = 1, nec

    ! z3 = A21xIxI u^e .........................................................

    do k = 1, na1
    do j = 1, na2
    do i = 1, na2
      z3(i,j,k,1) = 0
      z3(i,j,k,2) = 0
      do p = 1, na2
        z3(i,j,k,1) = z3(i,j,k,1) + A12(k,p,1) * u(i,j,p,e)
        z3(i,j,k,2) = z3(i,j,k,2) + A12(k,p,2) * u(i,j,p,e)
      end do
    end do
    end do
    end do

    ! z2 = IxA12xI z3 ..........................................................

    do k = 1, na1
    do j = 1, na1
    do i = 1, na2
      z2(i,j,k,1,1) = 0
      z2(i,j,k,2,1) = 0
      z2(i,j,k,1,2) = 0
      z2(i,j,k,2,2) = 0
      do p = 1, na2
        z2(i,j,k,1,1) = z2(i,j,k,1,1) + A12(j,p,1) * z3(i,p,k,1)
        z2(i,j,k,2,1) = z2(i,j,k,2,1) + A12(j,p,2) * z3(i,p,k,1)
        z2(i,j,k,1,2) = z2(i,j,k,1,2) + A12(j,p,1) * z3(i,p,k,2)
        z2(i,j,k,2,2) = z2(i,j,k,2,2) + A12(j,p,2) * z3(i,p,k,2)
      end do
    end do
    end do
    end do

    ! v^e = IxIxA12 z2 .........................................................

    ! fine element offset
    l = 8 * (e-1)

    do k = 1, na1
    do j = 1, na1
    do i = 1, na1
      v(i,j,k,l+1) = 0
      v(i,j,k,l+2) = 0
      v(i,j,k,l+3) = 0
      v(i,j,k,l+4) = 0
      v(i,j,k,l+5) = 0
      v(i,j,k,l+6) = 0
      v(i,j,k,l+7) = 0
      v(i,j,k,l+8) = 0
      do p = 1, na2                                            ! position
        v(i,j,k,l+1) = v(i,j,k,l+1) + A(i,p,1) * z2(p,j,k,1,1) !  1,1,1
        v(i,j,k,l+2) = v(i,j,k,l+2) + A(i,p,2) * z2(p,j,k,1,1) !  2,1,1
        v(i,j,k,l+3) = v(i,j,k,l+3) + A(i,p,1) * z2(p,j,k,2,1) !  1,2,1
        v(i,j,k,l+4) = v(i,j,k,l+4) + A(i,p,2) * z2(p,j,k,2,1) !  2,2,1
        v(i,j,k,l+5) = v(i,j,k,l+5) + A(i,p,1) * z2(p,j,k,1,2) !  1,1,2
        v(i,j,k,l+6) = v(i,j,k,l+6) + A(i,p,2) * z2(p,j,k,1,2) !  2,1,2
        v(i,j,k,l+7) = v(i,j,k,l+7) + A(i,p,1) * z2(p,j,k,2,2) !  1,2,2
        v(i,j,k,l+8) = v(i,j,k,l+8) + A(i,p,2) * z2(p,j,k,2,2) !  2,2,2
      end do
    end do
    end do
    end do

  end do
  !$omp end do

end subroutine TPO_A12A12A12_Gen_RWP

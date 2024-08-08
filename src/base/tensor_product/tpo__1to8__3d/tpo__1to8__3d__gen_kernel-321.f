!> summary:   3d generic interpolation operator for regular 1:8 h-refinement
!> author:    Jörg Stiller
!> date:      2024/07/24
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> This variant runs through the directions in the order 3,2,1
!===============================================================================

!-----------------------------------------------------------------------------
!> Generic 3d interpolation operator for regular 1:8 h-refinement

subroutine TPO_1To8_Gen_RWP(A, u, v)
  real(RWP), contiguous, intent(in)  :: A(:,:,:)   !< 1D operator (na1,na2,2)
  real(RWP), contiguous, intent(in)  :: u(:,:,:,:) !< operand
  real(RWP), contiguous, intent(out) :: v(:,:,:,:) !< result

  !---------------------------------------------------------------------------
  ! local variables

  real(RWP) :: y ( size(A,2), size(A,1), size(A,1), 2, 2 )
  real(RWP) :: z ( size(A,2), size(A,2), size(A,1), 2 )
  real(RWP) :: At( size(A,2), size(A,1), 2 )
  integer   :: na1, na2
  integer   :: e, i, j, k, p, o

  !---------------------------------------------------------------------------
  ! initialization

  na1 = size(A,1)
  na2 = size(A,2)

  At(:,:,1) = transpose(A(:,:,1))
  At(:,:,2) = transpose(A(:,:,2))

  !---------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e,o)
  do e = 1, size(u,4)

    ! z = AxIxI u^e .........................................................

    ! z([:,:],[:,:]) = u([:,:],[:],e) • At([:],[:,:])
    do k = 1, na1
    do j = 1, na2
    do i = 1, na2
      z(i,j,k,1:2) = 0
      do p = 1, na2
        z(i,j,k,1) = z(i,j,k,1) + u(i,j,p,e) * At(p,k,1)
        z(i,j,k,2) = z(i,j,k,2) + u(i,j,p,e) * At(p,k,2)
      end do
    end do
    end do
    end do

    ! y = IxAxI z ..............................................................

    ! y([:],[:],k,m,n) = z([:],[:],k,n) • At([:],[:],m)
    do k = 1, na1
    do j = 1, na1
    do i = 1, na2
      y(i,j,k,1:2,1:2) = 0
      do p = 1, na2
        y(i,j,k,1,1) = y(i,j,k,1,1) + z(i,p,k,1) * At(p,j,1)
        y(i,j,k,2,1) = y(i,j,k,2,1) + z(i,p,k,1) * At(p,j,2)
        y(i,j,k,1,2) = y(i,j,k,1,2) + z(i,p,k,2) * At(p,j,1)
        y(i,j,k,2,2) = y(i,j,k,2,2) + z(i,p,k,2) * At(p,j,2)
      end do
    end do
    end do
    end do

    ! v^e = IxIxA y ............................................................

    ! offset of fine element index
    o = 8 * (e-1)

    ! v([:],[:,:],o+1) = A([:],[:],1) • y([:],[:,:],1,1)
    !                  ⋮
    ! v([:],[:,:],o+8) = A([:],[:],2) • y([:],[:,:],2,2)
    do k = 1, na1
    do j = 1, na1
    do i = 1, na1
      v(i,j,k,o+1:o+8) = 0
      do p = 1, na2                                           ! position
        v(i,j,k,o+1) = v(i,j,k,o+1) + A(i,p,1) * y(p,j,k,1,1) !  1,1,1
        v(i,j,k,o+2) = v(i,j,k,o+2) + A(i,p,2) * y(p,j,k,1,1) !  2,1,1
        v(i,j,k,o+3) = v(i,j,k,o+3) + A(i,p,1) * y(p,j,k,2,1) !  1,2,1
        v(i,j,k,o+4) = v(i,j,k,o+4) + A(i,p,2) * y(p,j,k,2,1) !  2,2,1
        v(i,j,k,o+5) = v(i,j,k,o+5) + A(i,p,1) * y(p,j,k,1,2) !  1,1,2
        v(i,j,k,o+6) = v(i,j,k,o+6) + A(i,p,2) * y(p,j,k,1,2) !  2,1,2
        v(i,j,k,o+7) = v(i,j,k,o+7) + A(i,p,1) * y(p,j,k,2,2) !  1,2,2
        v(i,j,k,o+8) = v(i,j,k,o+8) + A(i,p,2) * y(p,j,k,2,2) !  2,2,2
      end do
    end do
    end do
    end do

  end do
  !$omp end do

end subroutine TPO_1To8_Gen_RWP

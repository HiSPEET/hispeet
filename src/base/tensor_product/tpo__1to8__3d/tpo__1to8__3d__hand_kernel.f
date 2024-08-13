!> summary:   3d regular 1:8 interpolation using hand-crafted suboperators
!> author:    Jörg Stiller
!> date:      2024/07/30
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> This variant runs through the directions in the order 3,2,1
!===============================================================================

!-----------------------------------------------------------------------------
!> 3d regular 1:8 interpolation using hand-crafted suboperators

subroutine PROC(TPO_1To8_Hand__,_NA1_,_NA2_)(ne, A, u, v)
  integer,   intent(in)  :: ne                        !< num coarse elements
  real(RWP), intent(in)  :: A(_NA1_,_NA2_,2)          !< 1D operator
  real(RWP), intent(in)  :: u(_NA2_,_NA2_,_NA2_,ne)   !< operand
  real(RWP), intent(out) :: v(_NA1_,_NA1_,_NA1_,ne*8) !< result

  !---------------------------------------------------------------------------
  ! local variables

  real(RWP), parameter :: alpha = 1
  real(RWP), parameter :: beta  = 0
  real(RWP) :: y (_NA2_,_NA1_,_NA1_,2,2)
  real(RWP) :: z (_NA2_,_NA2_,_NA1_,2)
  real(RWP) :: At(_NA2_,_NA1_,2)

  integer, parameter :: NA2_NA2 = _NA2_ * _NA2_
  integer :: e, o

  !---------------------------------------------------------------------------
  ! initialization

  At(:,:,1) = transpose(A(:,:,1))
  At(:,:,2) = transpose(A(:,:,2))

  ! make sure that aux array does not contain NaNs
  y = 0
  z = 0

  !---------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e,o)
  do e = 1, ne

    ! offset of fine element index
    o = 8 * (e-1)

    ! avoid trouble with NaNs
    v(:,:,:,o+1:o+8) = 0

    ! z = AxIxI uᵉ
    call PROC(CtxIxI__,_NA2_,_NA1_)(NA2_NA2, At(:,:,1), alpha, beta, u(:,:,:,e), z(:,:,:,1))
    call PROC(CtxIxI__,_NA2_,_NA1_)(NA2_NA2, At(:,:,2), alpha, beta, u(:,:,:,e), z(:,:,:,2))

    ! y = IxAxI z
    call PROC(IxBtxI__,_NA2_,_NA1_)(_NA2_, _NA1_, At(:,:,1), alpha, beta, z(:,:,:,1), y(:,:,:,1,1))
    call PROC(IxBtxI__,_NA2_,_NA1_)(_NA2_, _NA1_, At(:,:,2), alpha, beta, z(:,:,:,1), y(:,:,:,2,1))
    call PROC(IxBtxI__,_NA2_,_NA1_)(_NA2_, _NA1_, At(:,:,1), alpha, beta, z(:,:,:,2), y(:,:,:,1,2))
    call PROC(IxBtxI__,_NA2_,_NA1_)(_NA2_, _NA1_, At(:,:,2), alpha, beta, z(:,:,:,2), y(:,:,:,2,2))

    ! v = IxIxA y
    call PROC(IxIxAt__,_NA2_,_NA1_)(_NA1_, _NA1_, At(:,:,1), alpha, beta, y(:,:,:,1,1), v(:,:,:,o+1))
    call PROC(IxIxAt__,_NA2_,_NA1_)(_NA1_, _NA1_, At(:,:,2), alpha, beta, y(:,:,:,1,1), v(:,:,:,o+2))
    call PROC(IxIxAt__,_NA2_,_NA1_)(_NA1_, _NA1_, At(:,:,1), alpha, beta, y(:,:,:,2,1), v(:,:,:,o+3))
    call PROC(IxIxAt__,_NA2_,_NA1_)(_NA1_, _NA1_, At(:,:,2), alpha, beta, y(:,:,:,2,1), v(:,:,:,o+4))
    call PROC(IxIxAt__,_NA2_,_NA1_)(_NA1_, _NA1_, At(:,:,1), alpha, beta, y(:,:,:,1,2), v(:,:,:,o+5))
    call PROC(IxIxAt__,_NA2_,_NA1_)(_NA1_, _NA1_, At(:,:,2), alpha, beta, y(:,:,:,1,2), v(:,:,:,o+6))
    call PROC(IxIxAt__,_NA2_,_NA1_)(_NA1_, _NA1_, At(:,:,1), alpha, beta, y(:,:,:,2,2), v(:,:,:,o+7))
    call PROC(IxIxAt__,_NA2_,_NA1_)(_NA1_, _NA1_, At(:,:,2), alpha, beta, y(:,:,:,2,2), v(:,:,:,o+8))

  end do
  !$omp end do

end subroutine PROC(TPO_1To8_Hand__,_NA1_,_NA2_)

!> summary:  3D divergence element operator using LIBXSSM (R)
!> author:   Erik Pfister
!> date:     2019/11/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Div_3d_R__XSMM

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RNP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Div_R_XSMM

  type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

contains

  !-----------------------------------------------------------------------------
  !> 3D divergence element operator based on LIBXSMM (R)

  subroutine TPO_Div_R_XSMM(np, ne, Ds, dx, u, v)
    integer,   intent(in)  :: np               !< num points per direction
    integer,   intent(in)  :: ne               !< num elements
    real(RNP), intent(in)  :: Ds(np, np)       !< standard diff matrix
    real(RNP), intent(in)  :: dx(3)            !< element extensions
    real(RNP), intent(in)  :: u(np,np,np,ne,3) !< operand
    real(RNP), intent(out) :: v(np,np,np,ne)   !< result

    real(RNP) :: A(np,np), A_1T(np,np)
    real(RNP) :: A_2(np,np), A_3(np, np)
    real(RNP), parameter :: alpha  = 1
    real(RNP), parameter :: beta   = 1

    real(RNP) :: g(3)
    integer   :: e, i, j, k

    !-----------------------------------------------------------------------------
    ! initialization

    A = transpose(Ds)

    ! metric coefficients
    g = 2 / dx

    A_1T  = g(1) * transpose(A)
    A_2   = g(2) * A
    A_3   = g(3) * A

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_Elliptic_RLCI_XSMM: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, ne
    
      v(:,:,:,e)=0     

      ! apply operator in direction 1
      call LIBXSMM_DMMCall(xmm_1, A_1T, u(:,:,1,e,1) , v(:,:,:,e))

      ! direction 2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,2), A_2, v(:,:,k,e))
      end do

      ! direction 3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,1,e,3), A_3, v(:,:,:,e))

    end do

  end subroutine TPO_Div_R_XSMM

  !=============================================================================

#endif

end module TPO__Div_3d_R__XSMM

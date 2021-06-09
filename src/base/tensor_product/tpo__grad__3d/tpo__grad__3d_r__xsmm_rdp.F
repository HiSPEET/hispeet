!> summary:  3D gradient element operator using LIBXSSM (R)
!> author:   Erik Pfister
!> date:     2019/11/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Grad__3D_R__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Grad_R_XSMM

  interface TPO_Grad_R_XSMM
    module procedure TPO_Grad_R_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D gradient element operator based on LIBXSMM (R)

  subroutine TPO_Grad_R_XSMM_RDP(np, ne, Ds, dx, u, v)
    integer,   intent(in)  :: np               !< num points per direction
    integer,   intent(in)  :: ne               !< num elements
    real(RDP), intent(in)  :: Ds(np, np)       !< standard diff matrix
    real(RDP), intent(in)  :: dx(3)            !< element extensions
    real(RDP), intent(in)  :: u(np,np,np,ne)   !< 3D scalar field
    real(RDP), intent(out) :: v(np,np,np,ne,3) !< element-wise gradient of u

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP) :: Ds_1(np,np)
    real(RDP) :: Ds_2t(np,np), Ds_3t(np, np)
    real(RDP), parameter :: alpha  = 1
    real(RDP), parameter :: beta   = 0

    real(RDP) :: g(3)
    integer   :: e, k

    !-----------------------------------------------------------------------------

    ! metric coefficients
    g = 2 / dx

    Ds_1  = g(1) * Ds
    Ds_2t = g(2) * transpose(Ds)
    Ds_3t = g(3) * transpose(Ds)

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_Grad_R_XSMM: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, ne

      ! v1 = du/dx1
      call LIBXSMM_DMMCall(xmm_1, Ds_1, u(:,:,:,e) , v(:,:,:,e,1))

      ! v2 = du/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e), Ds_2t, v(:,:,k,e,2))
      end do

      ! v3 = du/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e), Ds_3t, v(:,:,:,e,3))

    end do

  end subroutine TPO_Grad_R_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Grad__3D_R__XSMM_RDP

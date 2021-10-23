!> summary:  3D divergence element operator using LIBXSMM for real(RDP) (D)
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/08/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Div__3D_D__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Div_D_XSMM

  interface TPO_Div_D_XSMM
    module procedure TPO_Div_D_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D divergence element operator based on LIBXSMM (R)

  subroutine TPO_Div_D_XSMM_RDP(np, ne, Ds, Ji, u, v)
    integer,   intent(in)  :: np               !< num points per direction
    integer,   intent(in)  :: ne               !< num elements
    real(RDP), intent(in)  :: Ds(np, np)       !< standard diff matrix
    real(RDP), intent(in)  :: Ji(:,:,:,:,:,:)  !< inverse Jacobi matrix
    real(RDP), intent(in)  :: u(np,np,np,ne,3) !< 3D vector field
    real(RDP), intent(out) :: v(np,np,np,ne)   !< element-wise divergence of u

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP), parameter :: alpha  = 1
    real(RDP), parameter :: beta   = 0

    real(RDP) :: Ds_t(np,np)

    real(RDP), dimension(np,np,np) :: r, s, t, w

    integer :: e, k

    !-----------------------------------------------------------------------------
    ! Initialization

    Ds_t = transpose(Ds)

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! Evaluation

    !$omp do
    do e = 1, ne

      ! direction 1 ............................................................

      ! r = ∂u₁/∂ξ
      call LIBXSMM_DMMCall(xmm_1, Ds, u(:,:,:,e,1), r)

      ! s = ∂u₁/∂η
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,1), Ds_t, s(:,:,k))
      end do

      ! t = ∂u₁/∂ζ
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,1), Ds_t, t)

      w = r * Ji(:,:,:,e,1,1) &
        + s * Ji(:,:,:,e,2,1) &
        + t * Ji(:,:,:,e,3,1)

      ! direction 2 ............................................................

      ! r = ∂u₂/∂ξ
      call LIBXSMM_DMMCall(xmm_1, Ds, u(:,:,:,e,2), r)

      ! s = ∂u₂/∂η
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,2), Ds_t, s(:,:,k))
      end do

      ! t = ∂u₂/∂ζ
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,2), Ds_t, t)

      w = w + r * Ji(:,:,:,e,1,2) &
            + s * Ji(:,:,:,e,2,2) &
            + t * Ji(:,:,:,e,3,2)

      ! direction 3 ............................................................

      ! r = ∂u₃/∂ξ
      call LIBXSMM_DMMCall(xmm_1, Ds, u(:,:,:,e,3), r)

      ! s = ∂u₂/∂η
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,3), Ds_t, s(:,:,k))
      end do

      ! t = ∂u₂/∂ζ
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,3), Ds_t, t)

      w = w + r * Ji(:,:,:,e,1,3) &
            + s * Ji(:,:,:,e,2,3) &
            + t * Ji(:,:,:,e,3,3)

      v(:,:,:,e) = w

    end do

  end subroutine TPO_Div_D_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Div__3D_D__XSMM_RDP

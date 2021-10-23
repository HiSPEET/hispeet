!> summary:  3D gradient element operator using LIBXSMM (D)
!> author:   Jerome Michel, Erik Pfister
!> date:     2021/08/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Grad__3D_D__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  ! like this? if use directly, maybe not needed
  public :: TPO_Grad_D_XSMM

  interface TPO_Grad_D_XSMM
    module procedure TPO_Grad_D_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D gradient element operator based on LIBXSMM (R)

  subroutine TPO_Grad_D_XSMM_RDP(np, ne, Ds, Ji, u, v)
    integer,   intent(in)  :: np               !< num points per direction
    integer,   intent(in)  :: ne               !< num elements
    real(RDP), intent(in)  :: Ds(np, np)       !< standard diff matrix
    real(RDP), intent(in)  :: Ji(:,:,:,:,:,:)  !< inverse Jacobi matrix
    real(RDP), intent(in)  :: u(np,np,np,ne)   !< 3D scalar field
    real(RDP), intent(out) :: v(np,np,np,ne,3) !< element-wise gradient of u

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP), parameter :: alpha  = 1
    real(RDP), parameter :: beta   = 0

    real(RDP) :: Ds_t(np,np)


    real(RDP), dimension(np,np,np) :: r,s,t

    integer   :: e, k

    !-----------------------------------------------------------------------------
    ! intialization

    Ds_t  = transpose(Ds)

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, ne


      ! v1 = du/dx1
      call LIBXSMM_DMMCall(xmm_1, Ds, u(:,:,:,e) , r)

      ! v2 = du/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e), Ds_t, s(:,:,k))
      end do

      ! v3 = du/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e), Ds_t, t)


      ! is there a faster way?
      v(:,:,:,e,1) = r(:,:,:) * Ji(:,:,:,e,1,1) + &
                     s(:,:,:) * Ji(:,:,:,e,2,1) + &
                     t(:,:,:) * Ji(:,:,:,e,3,1)

      v(:,:,:,e,2) = r(:,:,:) * Ji(:,:,:,e,1,2) + &
                     s(:,:,:) * Ji(:,:,:,e,2,2) + &
                     t(:,:,:) * Ji(:,:,:,e,3,2)

      v(:,:,:,e,3) = r(:,:,:) * Ji(:,:,:,e,1,3) + &
                     s(:,:,:) * Ji(:,:,:,e,2,3) + &
                     t(:,:,:) * Ji(:,:,:,e,3,3)

    end do

  end subroutine TPO_Grad_D_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Grad__3D_D__XSMM_RDP

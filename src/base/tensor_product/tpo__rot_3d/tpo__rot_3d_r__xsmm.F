!> summary:  3D rotation element operator using LIBXSSM (R)
!> author:   Erik Pfister
!> date:     2019/11/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Rot_3d_R__XSMM

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RNP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Rot_R_XSMM

  type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

contains

  !-----------------------------------------------------------------------------
  !> 3D rotation element operator based on LIBXSMM (R)

  subroutine TPO_Rot_R_XSMM(np, ne, Ds, dx, u, v)
    integer,   intent(in)  :: np               !< num points per direction
    integer,   intent(in)  :: ne               !< num elements
    real(RNP), intent(in)  :: Ds(np, np)       !< standard diff matrix
    real(RNP), intent(in)  :: dx(3)            !< element extensions
    real(RNP), intent(in)  :: u(np,np,np,ne,3) !< 3D vector field
    real(RNP), intent(out) :: v(np,np,np,ne,3) !< element-wise rotation of u

    real(RNP) :: Ds_1(np,np)
    real(RNP) :: Ds_2T(np,np), Ds_3T(np, np)
    real(RNP), parameter :: alpha  = 1
    real(RNP), parameter :: beta   = 1

    real(RNP) :: g(3)
    integer   :: e, k

    !-----------------------------------------------------------------------------
    
    ! metric coefficients
    g = 2 / dx

    Ds_1  = g(1) * Ds
    Ds_2T = g(2) * transpose(Ds)
    Ds_3T = g(3) * transpose(Ds)

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
    
      v(:,:,:,e,1) = 0
      v(:,:,:,e,2) = 0
      v(:,:,:,e,3) = 0     

      ! v2 = du1/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,1,e,1), Ds_3T,  v(:,:,:,e,2))

      ! v3 =-du1/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,1), -Ds_2T, v(:,:,k,e,3))
      end do

      ! v3+= du2/dx1
      call LIBXSMM_DMMCall(xmm_1, Ds_1,  u(:,:,1,e,2), v(:,:,:,e,3))

      ! v1 =-du1/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,1,e,2), -Ds_3T,  v(:,:,:,e,1))

      ! v1+= du1/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,3), Ds_2T, v(:,:,k,e,1))
      end do

      ! v2+=-du3/dx1
      call LIBXSMM_DMMCall(xmm_1, -Ds_1,  u(:,:,1,e,3), v(:,:,:,e,2))

    end do

  end subroutine TPO_Rot_R_XSMM

  !=============================================================================

#endif

end module TPO__Rot_3d_R__XSMM

!> summary:  3D divergence element operator using LIBXSMM for real(RDP) (D)
!> author:   Jerome Michel,Erik Pfister
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

    real(RDP) :: Ds_1(np,np)
    real(RDP) :: Ds_2t(np,np), Ds_3t(np, np)

    real(RDP), dimension(np,np,np) :: r, s, t

    integer   :: e, k

    !-----------------------------------------------------------------------------
    ! metric coefficients

    Ds_1  = Ds
    Ds_2t = transpose(Ds)
    Ds_3t = transpose(Ds)

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_Div_D_XSMM: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, ne

      v(:,:,:,e) = 0


      ! direction 1

      r = 0;
      s = 0;
      t = 0;

      ! r = du1/dx1
      call LIBXSMM_DMMCall(xmm_1, Ds_1, u(:,:,:,e,1) , r(:,:,:))

      ! s= du1/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,1), Ds_2t, s(:,:,k))
      end do

      ! t= du1/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,1), Ds_3t, t(:,:,:))


      ! consider SIMD vectorization using !$omp simd or LIBXSMM

      v(:,:,:,e) =  r(:,:,:) * Ji(:,:,:,e,1,1) &
                  + s(:,:,:) * Ji(:,:,:,e,2,1) &
                  + t(:,:,:) * Ji(:,:,:,e,3,1)


      ! direction 2

      r = 0;
      s = 0;
      t = 0;

      ! r = du2/dx1
      call LIBXSMM_DMMCall(xmm_1, Ds_1, u(:,:,:,e,2) , r(:,:,:))

      ! s= du2/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,2), Ds_2t, s(:,:,k))
      end do

      ! t= du2/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,2), Ds_3t, t(:,:,:))


      ! consider SIMD vectorization using !$omp simd or LIBXSMM

      v(:,:,:,e) = v(:,:,:,e) + r(:,:,:) * Ji(:,:,:,e,1,2) &
                              + s(:,:,:) * Ji(:,:,:,e,2,2) &
                              + t(:,:,:) * Ji(:,:,:,e,3,2)


      ! direction 3

      r = 0;
      s = 0;
      t = 0;

      ! r= du3/dx1
      call LIBXSMM_DMMCall(xmm_1, Ds_1, u(:,:,:,e,3) , r(:,:,:))

      ! s= du3/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,3), Ds_2t, s(:,:,k))
      end do

      ! t= du3/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,3), Ds_3t, t(:,:,:))


      ! consider SIMD vectorization using !$omp simd or LIBXSMM

      v(:,:,:,e) = v(:,:,:,e) + r(:,:,:) * Ji(:,:,:,e,1,3) &
                              + s(:,:,:) * Ji(:,:,:,e,2,3) &
                              + t(:,:,:) * Ji(:,:,:,e,3,3)


    end do

  end subroutine TPO_Div_D_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Div__3D_D__XSMM_RDP

!> summary:  3D rotrot element operator using LIBXSSM (R)
!> author:   Erik Pfister
!> date:     2020/03/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__RotRot_3d_R__XSMM

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RNP
  use LIBXSMM
  implicit none
  private

  public :: TPO_RotRot_R_XSMM

contains

  !-----------------------------------------------------------------------------
  !> 3D rotrot element operator based on LIBXSMM (R)

  subroutine TPO_RotRot_R_XSMM(np, ne, Ds, dx, u, v)
    integer,   intent(in)  :: np               !< num points per direction
    integer,   intent(in)  :: ne               !< num elements
    real(RNP), intent(in)  :: Ds(np, np)       !< standard diff matrix
    real(RNP), intent(in)  :: dx(3)            !< element extensions
    real(RNP), intent(in)  :: u(np,np,np,ne,3) !< 3D vector field
    real(RNP), intent(out) :: v(np,np,np,ne,3) !< v = rot(rot(u))

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RNP) :: Ds_1(np,np), DDs_1(np,np)
    real(RNP) :: Ds_2t(np,np), Ds_3t(np,np)
    real(RNP) :: DDs_2t(np,np), DDs_3t(np,np)
    real(RNP) :: DDs_t(np,np)
    real(RNP) :: div_u(np,np,np)
    real(RNP), parameter :: alpha  = 1
    real(RNP), parameter :: beta   = 1

    real(RNP) :: g(3)
    integer   :: e, i, j, k

    !-----------------------------------------------------------------------------

    ! metric coefficients
    g = 2 / dx

    Ds_1  = g(1) * Ds
    Ds_2t = g(2) * transpose(Ds)
    Ds_3t = g(3) * transpose(Ds)

    do j = 1, np
    do i = 1, np
      DDs_t(j,i) = sum(Ds(i,:) * Ds(:,j))
    end do
    end do

    DDs_1  = -(g(1) * g(1)) * transpose(DDs_t)
    DDs_2t = -(g(2) * g(2)) * DDs_t
    DDs_3t = -(g(3) * g(3)) * DDs_t

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_RotRot_R_XSMM: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, ne

      div_u(:,:,:) = 0
      v(:,:,:,e,1) = 0
      v(:,:,:,e,2) = 0
      v(:,:,:,e,3) = 0

      ! div(u) = du1/dx1
      call LIBXSMM_DMMCall(xmm_1, Ds_1, u(:,:,:,e,1),  div_u(:,:,:))

      ! div(u)+= du2/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,2), Ds_2t, div_u(:,:,k))
      end do

      ! div(u)+= du3/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,3), Ds_3t, div_u(:,:,:))


      ! direction 1
      ! v1 = d(div_u)/dx1
      call LIBXSMM_DMMCall(xmm_1, Ds_1, div_u(:,:,:),  v(:,:,:,e,1))

      ! v1+=id( du1/dx1 )/dx1
      call LIBXSMM_DMMCall(xmm_1, DDs_1, u(:,:,:,e,1), v(:,:,:,e,1))

      ! v1+=id( du1/dx2 )/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,1), DDs_2t, v(:,:,k,e,1))
      end do

      ! v1+=id( du1/dx3 )/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,1), DDs_3t, v(:,:,:,e,1))


      !direction 2
      ! v2 = d(div_u)/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, div_u(:,:,k), Ds_2t, v(:,:,k,e,2))
      end do

      ! v2+=-d( du2/dx1 )/dx1
      call LIBXSMM_DMMCall(xmm_1, DDs_1, u(:,:,:,e,2), v(:,:,:,e,2))

      ! v2+=-d( du2/dx2 )/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,2), DDs_2t, v(:,:,k,e,2))
      end do

      ! v2+=-d( du2/dxx )/dxx
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,2), DDs_3t, v(:,:,:,e,2))


      !direction 3
      ! v3 = d(div_u)/dx3
      call LIBXSMM_DMMCall(xmm_3, div_u(:,:,:), Ds_3t, v(:,:,:,e,3))

      ! v3+=-d( du3/dx1 )/dx1
      call LIBXSMM_DMMCall(xmm_1, DDs_1, u(:,:,:,e,3), v(:,:,:,e,3))

      ! v3+=-d( du3/dx2 )/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,3), DDs_2t, v(:,:,k,e,3))
      end do

      ! v3+=-d( du3/dx3 )/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,3), DDs_3t, v(:,:,:,e,3))

    end do

  end subroutine TPO_RotRot_R_XSMM

  !=============================================================================

#endif

end module TPO__RotRot_3d_R__XSMM

!> summary:  3D element diffusion operator using LIBXSMM (DLCI) for real(RDP)
!> author:   Jerome Michel, Erik Pfister, Joerg Stiller
!> date:     2021/08/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diffusion__3D_DLCI__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Diffusion_DLCI_XSMM

  interface TPO_Diffusion_DLCI_XSMM
    module procedure TPO_Diffusion_DLCI_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element diffusion operator based on LIBXSMM (RLCI)

  subroutine TPO_Diffusion_DLCI_XSMM_RDP(Ms, Ds, Jd, G, lambda, nu, u, v)
    real(RDP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
    real(RDP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
    real(RDP), intent(in)  :: Jd(:,:,:,:) !< Jacobian determinant  (np,np,np,ne)
    real(RDP), intent(in)  :: G(:,:,:,:,:)!< Laplacian metrics   (np,np,np,ne,6)
    real(RDP), intent(in)  :: lambda      !< Helmholtz parameter
    real(RDP), intent(in)  :: nu          !< diffusivity
    real(RDP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
    real(RDP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP), parameter :: alpha = 1
    real(RDP), parameter :: beta  = 1

    real(RDP), dimension(size(Ms), size(Ms), size(Ms)) :: M, r, s, t, z

    real(RDP), dimension(size(Ms), size(Ms)) :: Ds_1, Ds_2t, Ds_3t
    real(RDP), dimension(size(Ms), size(Ms)) :: Ds_1t, Ds_2, Ds_3

    integer   :: np, ne
    integer   :: e, i, j, k

    !---------------------------------------------------------------------------
    ! initialization

    np = size(Ms)
    ne = size(u,4)

    ! element mass matrix
    do k = 1, np
    do j = 1, np
    do i = 1, np
      M(i,j,k) =  Ms(k) * Ms(j) * Ms(i)
    end do
    end do
    end do


    ! initialize operators

    ! first derivatives
    Ds_1 = Ds
    Ds_2t = transpose(Ds)
    Ds_3t = transpose(Ds)

    ! second derivatives
    Ds_1t = transpose(Ds)
    Ds_2 = Ds
    Ds_3 = Ds


    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_Diffusion_DLCI_XSMM: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, ne

      ! projection to element mass matrix and computation of the Helmholtz term
      do k = 1, np
      do j = 1, np
      do i = 1, np
        v(i,j,k,e) = lambda * M(i,j,k) * Jd(i,j,k,e) * u(i,j,k,e)
      end do
      end do
      end do

      r(:,:,:) = 0
      s(:,:,:) = 0
      t(:,:,:) = 0

      ! first derivatives

      call LIBXSMM_DMMCall(xmm_1, Ds_1, u(:,:,:,e), r(:,:,:))

      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e), Ds_2t, s(:,:,k))
      end do

      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e), Ds_3t, t(:,:,:))


      ! second derivative direction 1 and application of metric terms

      z(:,:,:) = 0

      do k = 1, np
      do j = 1, np
      do i = 1, np
        z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,1) * r(i,j,k) &
                                   + G(i,j,k,e,2) * s(i,j,k) &
                                   + G(i,j,k,e,3) * t(i,j,k) )
      end do
      end do
      end do

      call LIBXSMM_DMMCall(xmm_1, Ds_1t, z(:,:,:), v(:,:,:,e))


      ! second derivative direction 2 and application of metric terms

      z(:,:,:) = 0

      do k = 1, np
      do j = 1, np
      do i = 1, np
        z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,2) * r(i,j,k) &
                                   + G(i,j,k,e,4) * s(i,j,k) &
                                   + G(i,j,k,e,5) * t(i,j,k) )
      end do
      end do
      end do

      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, z(:,:,k), Ds_2, v(:,:,k,e))
      end do


      ! second derivative direction 3 and application of metric terms

      z(:,:,:) = 0

      do k = 1, np
      do j = 1, np
      do i = 1, np
        z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,3) * r(i,j,k) &
                                   + G(i,j,k,e,5) * s(i,j,k) &
                                   + G(i,j,k,e,6) * t(i,j,k) )
      end do
      end do
      end do

      call LIBXSMM_DMMCall(xmm_3, z(:,:,:), Ds_3, v(:,:,:,e))


    end do

  end subroutine TPO_Diffusion_DLCI_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Diffusion__3D_DLCI__XSMM_RDP

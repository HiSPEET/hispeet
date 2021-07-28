!> summary:  3d generic element diffusion operator (RLVI) using LIBXSMM
!> author:   Erik Pfister, Joerg Stiller
!> date:     2020/06/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diffusion__3D_RLVI__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Diffusion_RLVI_XSMM

  interface TPO_Diffusion_RLVI_XSMM
    module procedure TPO_Diffusion_RLVI_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Generic element diffusion operator (RLVI)

  subroutine TPO_Diffusion_RLVI_XSMM_RDP(Ms, Ds, dx, lambda, nu, u, v)
    real(RDP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
    real(RDP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
    real(RDP), intent(in)  :: lambda      !< Helmholtz parameter
    real(RDP), intent(in)  :: dx(3)       !< element extensions
    real(RDP), intent(in)  :: nu(:,:,:,:) !< diffusivity           (np,np,np,ne)
    real(RDP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
    real(RDP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP), parameter :: alpha = 1
    real(RDP), parameter :: beta  = 1

    !---------------------------------------------------------------------------
    ! local variables

    real(RDP), dimension(size(Ms), size(Ms), size(Ms)) :: M, M_u, z
    real(RDP), dimension(size(Ms), size(Ms))           :: Ms_Ds, Msi_Dst
    real(RDP), dimension(size(Ms), size(Ms))           :: Ms_Ds_2, Ms_Ds_3
    real(RDP), dimension(size(Ms), size(Ms))           :: Ms_Dst_1, Msi_Ds

    real(RDP) :: g(3), tmp
    integer   :: np, ne
    integer   :: e, i, j, k, p

    !---------------------------------------------------------------------------
    ! initialization

    np = size(Ms)
    ne = size(u,4)

    ! element mass matrix
    tmp = product(dx) / 8
    do k = 1, np
    do j = 1, np
    do i = 1, np
      M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
    end do
    end do
    end do

    ! modified differentiation operators
    !! to optimize cache use, the operators are used in transposed form
    !! to save loads, both operators are computed in one loop-nest
    do j = 1, np
    do i = 1, np
      Ms_Ds   (i,j) = Ms(i) * Ds(i,j)  ! = Ms Ds     =  [ (Ms Ds)ᵀ ]ᵀ
      Msi_Dst (i,j) = Ds(j,i) / Ms(i)  ! = Ms⁻¹ Dsᵀ  =  [  Ds Ms⁻¹ ]ᵀ
    end do
    end do

    ! coefficients
    g = 4 / dx**2

    Ms_Ds_2   = g(2) * Ms_Ds
    Ms_Ds_3   = g(3) * Ms_Ds
    Ms_Dst_1  = g(1) * transpose(Ms_Ds)
    Msi_Ds    =        transpose(Msi_Dst)

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_Diffusion_RLCI_XSMM: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do private(e)
    do e = 1, ne

      ! M_u = Mᵉuᵉ, vᵉ = λ M uᵉ ................................................

      do k = 1, np
      do j = 1, np
      do i = 1, np
        M_u(i,j,k) = M(i,j,k) * u(i,j,k,e)
        v(i,j,k,e) = lambda * M_u(i,j,k)
      end do
      end do
      end do

      ! direction 1 ............................................................
      z(:,:,:) = 0
      ! vᵉ += g₁ [I x I x (MˢDˢ)ᵀ] [νᵉ] [I x I x (Mˢ)⁻¹(Dˢ)ᵀ)] Mᵉuᵉ
      call LIBXSMM_DMMCall(xmm_1, Msi_Ds(:,:), M_u(:,:,:), z(:,:,:))
      z = nu(:,:,:,e) * z
      call LIBXSMM_DMMCall(xmm_1, Ms_Dst_1(:,:) , z(:,:,:), v(:,:,:,e))

      ! direction 2 ............................................................
      z(:,:,:) = 0
      ! vᵉ += g₂ [I x (MˢDˢ)ᵀ x I] [νᵉ] [I x (Mˢ)⁻¹(Dˢ)ᵀ x I] Mᵉuᵉ
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, M_u(:,:,k), Msi_Dst(:,:), z(:,:,k))
      end do
      z = nu(:,:,:,e) * z
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, z(:,:,k), Ms_Ds_2(:,:), v(:,:,k,e))
      end do

      ! direction 3 ............................................................
      z(:,:,:) = 0
      ! vᵉ += g₃ [(MˢDˢ)ᵀ x I x I] [νᵉ] [(Mˢ)⁻¹(Dˢ)ᵀ x I x I] Mᵉuᵉ
      ! z = [νᵉ] [(Mˢ)⁻¹(Dˢ)ᵀ x I x I] Mᵉuᵉ
      call LIBXSMM_DMMCall(xmm_3, M_u(:,:,:), Msi_Dst(:,:), z(:,:,:))
      z = nu(:,:,:,e) * z
      call LIBXSMM_DMMCall(xmm_3, z(:,:,:), Ms_Ds_3(:,:), v(:,:,:,e))

    end do
    !$omp end do

  end subroutine TPO_Diffusion_RLVI_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Diffusion__3D_RLVI__XSMM_RDP

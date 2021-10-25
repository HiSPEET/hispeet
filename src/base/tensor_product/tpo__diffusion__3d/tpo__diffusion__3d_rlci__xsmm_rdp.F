!> summary:  3D element diffusion operator using LIBXSSM (RLCI) for real(RDP)
!> author:   Erik Pfister, Joerg Stiller
!> date:     2019/11/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diffusion__3D_RLCI__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Diffusion_RLCI_XSMM

  interface TPO_Diffusion_RLCI_XSMM
    module procedure TPO_Diffusion_RLCI_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element diffusion operator based on LIBXSMM (RLCI)

  subroutine TPO_Diffusion_RLCI_XSMM_RDP(Ms, Ls, dx, lambda, nu, u, v, Ds, ub, qb)
    real(RDP), intent(in)  :: Ms(:)      !< 1D standard mass matrix         (np)
    real(RDP), intent(in)  :: Ls(:,:)    !< 1D standard stiffness matrix (np,np)
    real(RDP), intent(in)  :: dx(3)      !< element extensions
    real(RDP), intent(in)  :: lambda     !< Helmholtz parameter
    real(RDP), intent(in)  :: nu         !< diffusivity
    real(RDP), intent(in)  :: u(:,:,:,:) !< operand                (np,np,np,ne)
    real(RDP), intent(out) :: v(:,:,:,:) !< result                 (np,np,np,ne)

    real(RDP), optional, intent(in)    :: Ds(:,:)     !< 1D standard diff matrix
    real(RDP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
    real(RDP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

    !---------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP), parameter :: alpha = 1
    real(RDP), parameter :: beta  = 1

    real(RDP), dimension(size(Ms), size(Ms), size(Ms)) :: M, Mu_e
    real(RDP), dimension(size(Ms), size(Ms))           :: Lm, Lm_1t, Lm_2, Lm_3
    real(RDP), dimension(size(Ms))                     :: Bs1, Bs2

    real(RDP) :: c(3), tmp, tmp1, tmp2
    integer   :: np, ne
    integer   :: e, i, j, k, p
    logical   :: get_traces

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

    ! mass-weighted stiffness matrix: Lm = Ms^-1 Ls = (Ls Ms^-1)^T
    do j = 1, np
    do i = 1, np
      Lm(i,j) = Ls(i,j) / Ms(i)
    end do
    end do

    ! coefficients
    c = 4 * nu / dx**2

    ! initialize operators
    Lm_1t = c(1) * transpose(Lm)
    Lm_2  = c(2) * Lm
    Lm_3  = c(3) * Lm

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

    ! handling of boundary traces
    get_traces = present(Ds) .and. present(ub) .and. present(qb)

    if (get_traces) then
      Bs1 = -nu * Ds( 1,:)  ! n⋅νDˢ @ "left" boundary
      Bs2 =  nu * Ds(np,:)  ! n⋅νDˢ @ "right" boundary
    end if

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, ne

      ! projection to element mass matrix and computation of the Helmholtz term
      do k = 1, np
      do j = 1, np
      do i = 1, np
        Mu_e(i,j,k) = M(i,j,k) * u(i,j,k,e)
        v(i,j,k,e) = lambda * Mu_e(i,j,k)
      end do
      end do
      end do

      ! apply operator in direction 1
      call LIBXSMM_DMMCall(xmm_1, Lm_1t, Mu_e, v(:,:,:,e))

      ! direction 2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, Mu_e(:,:,k), Lm_2, v(:,:,k,e))
      end do

      ! direction 3
      call LIBXSMM_DMMCall(xmm_3, Mu_e, Lm_3, v(:,:,:,e))

      ! extract/compute traces .................................................

      if (get_traces) then

        ! ub = u, n⋅ν∇u  @ Γ₁ + Γ₂
        do k = 1, np
        do j = 1, np
          tmp1 = 0
          tmp2 = 0
          do p = 1, np
            tmp1 = tmp1 + Bs1(p) * u(p,j,k,e)
            tmp2 = tmp2 + Bs2(p) * u(p,j,k,e)
          end do
          ub(j,k,1,e) = u( 1,j,k,e)
          ub(j,k,2,e) = u(np,j,k,e)
          qb(j,k,1,e) = 2/dx(1) * tmp1
          qb(j,k,2,e) = 2/dx(1) * tmp2
        end do
        end do

        ! ub = u, n⋅ν∇u  @ Γ₃ + Γ₄
        do k = 1, np
        do i = 1, np
          tmp1 = 0
          tmp2 = 0
          do p = 1, np
            tmp1 = tmp1 + Bs1(p) * u(i,p,k,e)
            tmp2 = tmp2 + Bs2(p) * u(i,p,k,e)
          end do
          ub(i,k,3,e) = u(i, 1,k,e)
          ub(i,k,4,e) = u(i,np,k,e)
          qb(i,k,3,e) = 2/dx(2) * tmp1
          qb(i,k,4,e) = 2/dx(2) * tmp2
        end do
        end do

        ! ub = u, n⋅ν∇u  @ Γ₅ + Γ₆
        do j = 1, np
        do i = 1, np
          tmp1 = 0
          tmp2 = 0
          do p = 1, np
            tmp1 = tmp1 + Bs1(p) * u(i,j,p,e)
            tmp2 = tmp2 + Bs2(p) * u(i,j,p,e)
          end do
          ub(i,j,5,e) = u(i,j, 1,e)
          ub(i,j,6,e) = u(i,j,np,e)
          qb(i,j,5,e) = 2/dx(3) * tmp1
          qb(i,j,6,e) = 2/dx(3) * tmp2
        end do
        end do

      end if

    end do

  end subroutine TPO_Diffusion_RLCI_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Diffusion__3D_RLCI__XSMM_RDP

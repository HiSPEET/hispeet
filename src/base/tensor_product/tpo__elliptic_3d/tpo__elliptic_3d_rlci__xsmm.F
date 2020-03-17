!> summary:  3D elliptic element operator using LIBXSSM (RCLI)
!> author:   Erik Pfister
!> date:     2019/11/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Elliptic_3d_RLCI__XSMM

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RNP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Elliptic_RLCI_XSMM

  type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

contains

  !-----------------------------------------------------------------------------
  !> 3D elliptic element operator based on LIBXSMM (RCLI)

  subroutine TPO_Elliptic_RLCI_XSMM(np, ne, Ms, Ls, lambda, nu, dx, u, v)
    integer,   intent(in)  :: np             !< num points per direction
    integer,   intent(in)  :: ne             !< num elements
    real(RNP), intent(in)  :: Ms(np)         !< standard mass matrix
    real(RNP), intent(in)  :: Ls(np,np)      !< standard stiffness matrix
    real(RNP), intent(in)  :: lambda         !< Helmholtz parameter λ
    real(RNP), intent(in)  :: nu             !< diffusivity nu
    real(RNP), intent(in)  :: dx(3)          !< element extensions
    real(RNP), intent(in)  :: u(np,np,np,ne) !< operand
    real(RNP), intent(out) :: v(np,np,np,ne) !< result

    real(RNP), parameter :: alpha = 1
    real(RNP), parameter :: beta  = 1

    real(RNP) :: M(np,np,np), Lm(np,np), Mu_e(np,np,np)
    real(RNP) :: Lm_1T(np,np), Lm_2(np,np), Lm_3(np,np)
    real(RNP) :: c(3), tmp
    integer   :: e, i, j, k

    !---------------------------------------------------------------------------
    ! initialization

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
    Lm_1T = c(1) * transpose(Lm)
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

      stop "TPO_Elliptic_RLCI_XSMM: LIBXSMM_Dispatch failed"

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
        Mu_e(i,j,k) = M(i,j,k) * u(i,j,k,e)
        v(i,j,k,e) = lambda * Mu_e(i,j,k)
      end do
      end do
      end do

      ! apply operator in direction 1
      call LIBXSMM_DMMCall(xmm_1, Lm_1T, Mu_e, v(:,:,:,e))

      ! direction 2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, Mu_e(:,:,k), Lm_2, v(:,:,k,e))
      end do

      ! direction 3
      call LIBXSMM_DMMCall(xmm_3, Mu_e, Lm_3, v(:,:,:,e))

    end do

  end subroutine TPO_Elliptic_RLCI_XSMM

  !=============================================================================

#endif

end module TPO__Elliptic_3d_RLCI__XSMM

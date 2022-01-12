!> summary:  3d curvilinear element stress operator (DLCI)
!> author:   Joerg Stiller
!> date:     2022/01/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Stress__3D_DLCI
  use Kind_Parameters, only: RDP
! use TPO__Stress__3D_DLCI__Hand
! use TPO__Stress__3D_DLCI__XSMM
  use TPO__Stress__3D_DLCI__Gen ! workaround until the upper two become available
  implicit none
  private

  public :: TPO_Stress

  interface TPO_Stress
    module procedure TPO_Stress_DLCI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element stress operator

  subroutine TPO_Stress_DLCI_RDP(Ms, Ds, Jd, Ji, n, eta, zeta, v, fv, v_b, tau_b)

    real(RDP), contiguous, intent(in) :: Ms(:)
    !< 1D diagonal standard mass matrix (np)

    real(RDP), contiguous, intent(in)  :: Ds(:,:)
    !< 1D standard diff matrix (np,np)

    real(RDP), contiguous, intent(in)  :: Jd(:,:,:,:)
    !< Jacobian determinant (np,np,np,ne)

    real(RDP), contiguous, intent(in)  :: Ji(:,:,:,:,:,:)
    !< inverse Jacobian matrix  (np,np,np,ne,3,3)

    real(RDP), contiguous, intent(in)  :: n(:,:,:,:,:)
    !< unit normal vector at element boundary faces (np,np,6,ne,3)

    real(RDP), intent(in)  :: eta
    !< shear viscosity or modulus η

    real(RDP), intent(in)  :: zeta
    !< bulk viscosity or modulus ζ

    real(RDP), contiguous, intent(in)  :: v(:,:,:,:,:)
    !< velocity or displacement (np,np,np,ne,3)

    real(RDP), contiguous, intent(out) :: fv(:,:,:,:,:)
    !< element volume contribution to stress operator (np,np,np,ne,3)

    real(RDP), contiguous, intent(out) :: v_b(:,:,:,:)
    !< element boundary values of velocity or displacement (np,np,6,ne,3)

    real(RDP), contiguous, intent(out) :: tau_b(:,:,:,:)
    !< stress normal to element boundary fluxes, tau_b = n⁻⋅τ⁻

    ! workaround
    call TPO_Stress_DLCI_Gen_RWP(Ms, Ds, Jd, Ji, n, eta, zeta, v, &
                                 fv, v_b, tau_b)

  end subroutine TPO_Stress_DLCI_RDP

  !=============================================================================

end module TPO__Stress__3D_DLCI

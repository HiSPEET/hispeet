!> summary:  3d curvilinear element stress operator (DLCI)
!> author:   Joerg Stiller
!> date:     2022/01/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Stress__3D_DLCI
  use Kind_Parameters, only: RDP
  use TPO__Stress__3D_DLCI__Gen
  use TPO__Stress__3D_DLCI__XSMM

  implicit none
  private

  public :: TPO_Stress

  interface TPO_Stress
    module procedure TPO_Stress_DLCI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element stress operator -- constant viscosity

  subroutine TPO_Stress_DLCI_RDP(Ms, Ds, Jd, Ji, n, zeta, eta, v, fv, vb, sb)

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

    real(RDP), intent(in)  :: zeta
    !< bulk viscosity or modulus ζ

    real(RDP), intent(in)  :: eta
    !< shear viscosity or modulus η

    real(RDP), contiguous, intent(in)  :: v(:,:,:,:,:)
    !< velocity or displacement (np,np,np,ne,3)

    real(RDP), contiguous, intent(out) :: fv(:,:,:,:,:)
    !< element volume contribution to stress operator (np,np,np,ne,3)

    real(RDP), contiguous, optional, intent(out) :: vb(:,:,:,:,:)
    !< element boundary values of velocity or displacement (np,np,6,ne+ng,3),

    real(RDP), contiguous, optional, intent(out) :: sb(:,:,:,:,:)
    !< stress normal to element boundary fluxes, sb = n⁻⋅τ⁻ (np,np,6,ne+ng,3)

    select case(size(Ms))

#ifdef __LIBXSMM__
    case(3:)
      call TPO_Stress_DLCI_XSMM(Ms, Ds, Jd, Ji, n, zeta, eta, v, fv, vb, sb)
#endif

    case default
      call TPO_Stress_DLCI_Gen(Ms, Ds, Jd, Ji, n, zeta, eta, v, fv, vb, sb)
    end select

  end subroutine TPO_Stress_DLCI_RDP

  !=============================================================================

end module TPO__Stress__3D_DLCI

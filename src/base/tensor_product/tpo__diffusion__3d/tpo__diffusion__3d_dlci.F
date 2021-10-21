!> summary: 3D curvilinear diffusion operator with constant isotropic diffusivity
!> author:  Jerome Michel, Joerg Stiller
!> date:    2021/10/12
!> license: Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diffusion__3D_DLCI
  use Kind_Parameters, only: RDP
  use TPO__Diffusion__3D_DLCI__Hand
  use TPO__Diffusion__3D_DLCI__XSMM
  implicit none
  private

  public :: TPO_Diffusion

  interface TPO_Diffusion
    module procedure TPO_Diffusion_DLCI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element diffusion operator for scalars
  !>
  !> The procedure evaluates the diffusion operator with constant lambda and
  !> nu for deformed 3D elements. Given the scalar
  !> variable u(1:np,1:np,1:np,1:ne), the result is computed as
  !>   .......
  !>   .......
  !> for each point (i,j,k) in every element e.
  !>
  !> Ms and Ds are the 1D standard mass and differention matrices for polynomial
  !> order P = np-1, respectively

  subroutine TPO_Diffusion_DLCI_RDP(Ms, Ds, Jd, G, lambda, nu, u, v)
    real(RDP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
    real(RDP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
    real(RDP), intent(in)  :: Jd(:,:,:,:) !< Jacobian determinant  (np,np,np,ne)
    real(RDP), intent(in)  :: G(:,:,:,:,:)!< Laplacian metrics
    real(RDP), intent(in)  :: lambda      !< Helmholtz parameter
    real(RDP), intent(in)  :: nu          !< diffusivity
    real(RDP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
    real(RDP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

    select case(size(Ms))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Diffusion_DLCI_XSMM(Ms, Ds, Jd, G, lambda, nu, u, v)
#endif

    case default
      call TPO_Diffusion_DLCI_Hand(Ms, Ds, Jd, G, lambda, nu, u, v)
    end select

  end subroutine TPO_Diffusion_DLCI_RDP

  !=============================================================================

end module TPO__Diffusion__3D_DLCI

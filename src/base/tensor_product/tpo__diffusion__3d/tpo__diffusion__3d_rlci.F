!> summary:  3D Cartesian diffusion operator with constant isotropic diffusivity
!> author:   Joerg Stiller
!> date:     2020/06/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diffusion__3D_RLCI
  use Kind_Parameters, only: RDP
  use TPO__Diffusion__3D_RLCI__Hand
  use TPO__Diffusion__3D_RLCI__XSMM
  implicit none
  private

  public :: TPO_Diffusion

  interface TPO_Diffusion
    module procedure TPO_Diffusion_RLCI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element diffusion operator for scalars
  !>
  !> The procedure evaluates the diffusion operator with constant coefficients
  !> lambda and nu for equidistant 3D Cartesian elements. Given the scalar
  !> variable u(1:np,1:np,1:np,1:ne), the result is computed as
  !>
  !>     v(i,j,k,e)
  !>
  !>       = lambda * dx(1)*dx(2)*dx(3)/8 * Ms(i) * Ms(j) * Ms(k) * u(i,j,k,e)
  !>
  !>       + nu * dx(2)*dx(3)/(2*dx(1)) * Ms(j) * Ms(k)
  !>                                    * sum( Ms(:) * Ds(:,i) * u(:,j,k,e) )
  !>
  !>       + nu * dx(1)*dx(3)/(2*dx(2)) * Ms(i) * Ms(k)
  !>                                    * sum( Ms(:) * Ds(:,j) * u(i,:,k,e) )
  !>
  !>       + nu * dx(2)*dx(2)/(2*dx(3)) * Ms(i) * Ms(j)
  !>                                    * sum( Ms(:) * Ds(:,k) * u(i,j,:,e) )
  !>
  !> for each point (i,j,k) in every element e.
  !> Ms and Ls are the 1D standard mass and stiffness matrices for polynomial
  !> order P = np-1, respectively, and dx(1:3) the element dimensions.

  subroutine TPO_Diffusion_RLCI_RDP(Ms, Ls, dx, lambda, nu, u, v, Ds, ub, qb)
    real(RDP), intent(in)  :: Ms(:)      !< 1D standard mass matrix
    real(RDP), intent(in)  :: Ls(:,:)    !< 1D standard stiffness matrix
    real(RDP), intent(in)  :: dx(3)      !< element extensions
    real(RDP), intent(in)  :: lambda     !< Helmholtz parameter
    real(RDP), intent(in)  :: nu         !< diffusivity
    real(RDP), intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), intent(out) :: v(:,:,:,:) !< result

    real(RDP), optional, intent(in)    :: Ds(:,:)     !< 1D standard diff matrix
    real(RDP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
    real(RDP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

    select case(size(Ms))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Diffusion_RLCI_XSMM(Ms, Ls, dx, lambda, nu, u, v, Ds, ub, qb)
#endif

    case default
      call TPO_Diffusion_RLCI_Hand(Ms, Ls, dx, lambda, nu, u, v, Ds, ub, qb)

    end select

  end subroutine TPO_Diffusion_RLCI_RDP

  !=============================================================================

end module TPO__Diffusion__3D_RLCI

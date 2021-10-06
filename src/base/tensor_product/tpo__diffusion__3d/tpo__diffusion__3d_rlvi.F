!> summary:  3D Cartesian diffusion operator with variable isotropic diffusivity
!> author:   Joerg Stiller
!> date:     2020/06/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diffusion__3D_RLVI
  use Kind_Parameters, only: RDP
  use TPO__Diffusion__3D_RLVI__Hand
  use TPO__Diffusion__3D_RLVI__XSMM
  implicit none
  private

  public :: TPO_Diffusion

  interface TPO_Diffusion
    module procedure TPO_Diffusion_RLVI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element diffusion operator for scalars
  !>
  !> The procedure evaluates the diffusion operator with constant lambda and
  !> variable nu for equidistant 3D Cartesian elements. Given the scalar
  !> variable u(1:np,1:np,1:np,1:ne), the result is computed as
  !>
  !>     v(i,j,k,e)
  !>
  !>       = g0 * lambda * Ms(i) * Ms(j) * Ms(k) * u(i,j,k,e)
  !>
  !>       + g1 * Ms(j) * Ms(k)
  !>            * sum( nu(:,j,k) * Ms(:) * Ds(:,i) * matmul(Ds, u(:,j,k,e)) )
  !>
  !>       + g2 * Ms(k) * Ms(i)
  !>            * sum( nu(i,:,k) * Ms(:) * Ds(:,j) * matmul(Ds, u(i,:,k,e)) )
  !>
  !>       + g3 * Ms(i) * Ms(j)
  !>            * sum( nu(i,j,:) * Ms(:) * Ds(:,k) * matmul(Ds, u(i,j,:,e)) )
  !>
  !> for each point (i,j,k) in every element e.
  !>
  !> Ms and Ds are the 1D standard mass and differention matrices for polynomial
  !> order P = np-1, respectively, dx(1:3) the element dimensions, and
  !>
  !>     g0 = dx(1) * dx(2) * dx(3) / 8
  !>     g1 = dx(2) * dx(3) / (2 * dx(1))
  !>     g2 = dx(3) * dx(1) / (2 * dx(2))
  !>     g3 = dx(1) * dx(2) / (2 * dx(3))
  !>

  subroutine TPO_Diffusion_RLVI_RDP(Ms, Ds, dx, lambda, nu, u, v, ub, qb)
    real(RDP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
    real(RDP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
    real(RDP), intent(in)  :: dx(3)       !< element extensions
    real(RDP), intent(in)  :: lambda      !< Helmholtz parameter
    real(RDP), intent(in)  :: nu(:,:,:,:) !< diffusivity           (np,np,np,ne)
    real(RDP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
    real(RDP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

    real(RDP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
    real(RDP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

    select case(size(Ms))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Diffusion_RLVI_XSMM(Ms, Ds, dx, lambda, nu, u, v, ub, qb)
#endif

    case default
      call TPO_Diffusion_RLVI_Hand(Ms, Ds, dx, lambda, nu, u, v, ub, qb)

    end select

  end subroutine TPO_Diffusion_RLVI_RDP

  !=============================================================================

end module TPO__Diffusion__3D_RLVI

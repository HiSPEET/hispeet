!> summary:  3D Cartesian elliptic operator with variable isotropic diffusivity
!> author:   Joerg Stiller
!> date:     2020/06/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Elliptic__3D_RLVI
  use Kind_Parameters, only: RNP
  use TPO__Elliptic__3D_RLVI__Hand
  use TPO__Elliptic__3D_RLVI__XSMM
  implicit none
  private

  public :: TPO_Elliptic_RLVI

contains

  !-----------------------------------------------------------------------------
  !> 3D elliptic element operator for scalars
  !>
  !> The procedure evaluates the elliptic operator with constant lambda and
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

  subroutine TPO_Elliptic_RLVI(Ms, Ds, lambda, nu, dx, u, v)
    real(RNP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
    real(RNP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
    real(RNP), intent(in)  :: lambda      !< Helmholtz parameter
    real(RNP), intent(in)  :: nu(:,:,:,:) !< diffusivity           (np,np,np,ne)
    real(RNP), intent(in)  :: dx(3)       !< element extensions
    real(RNP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
    real(RNP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

    select case(size(Ms))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Elliptic_RLVI_XSMM(Ms, Ds, lambda, nu, dx, u, v)
#endif

    case default
      call TPO_Elliptic_RLVI_Hand(Ms, Ds, lambda, nu, dx, u, v)

    end select

  end subroutine TPO_Elliptic_RLVI

  !=============================================================================

end module TPO__Elliptic__3D_RLVI

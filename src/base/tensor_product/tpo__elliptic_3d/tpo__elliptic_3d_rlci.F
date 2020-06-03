!> summary:  3D Cartesian elliptic operator with constant isotropic diffusivity
!> author:   Joerg Stiller
!> date:     2020/06/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Elliptic_3d_RLCI
  use Kind_Parameters, only: RNP
  use TPO__Elliptic_3d_RLCI__Hand
  use TPO__Elliptic_3d_RLCI__XSMM
  implicit none
  private

  public :: TPO_Elliptic_RLCI

contains

  !-----------------------------------------------------------------------------
  !> 3D elliptic element operator for scalars
  !>
  !> The procedure evaluates the elliptic operator with constant coefficients
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

  subroutine TPO_Elliptic_RLCI(Ms, Ls, lambda, nu, dx, u, v)
    use Kind_Parameters,  only: RNP
    real(RNP), intent(in)  :: Ms(:)      !< 1D standard mass matrix
    real(RNP), intent(in)  :: Ls(:,:)    !< 1D standard stiffness matrix
    real(RNP), intent(in)  :: lambda     !< Helmholtz parameter
    real(RNP), intent(in)  :: nu         !< diffusivity
    real(RNP), intent(in)  :: dx(3)      !< element extensions
    real(RNP), intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), intent(out) :: v(:,:,:,:) !< result

    select case(size(Ms))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Elliptic_RLCI_XSMM(Ms, Ls, lambda, nu, dx, u, v)
#endif

    case default
      call TPO_Elliptic_RLCI_Hand(Ms, Ls, lambda, nu, dx, u, v)

    end select

  end subroutine TPO_Elliptic_RLCI

  !=============================================================================

end module TPO__Elliptic_3d_RLCI

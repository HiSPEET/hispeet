!> summary:  DG Schwarz operator for elliptic equations: base type
!> author:   Joerg Stiller
!> date:     2022/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module DG__Schwarz_Operator__3D__Base
  use Kind_Parameters
  use Mesh__3D
  use DG__Element_Operators__1D

  implicit none
  private

  public :: DG_SCHWARZ_BC_3D
  public :: DG_SchwarzOperator_3D

  !-----------------------------------------------------------------------------
  !> Supported boundary conditions
  !>
  !> Periodic boundaries are treated as interior: 'P' → ' '

  character, parameter :: DG_SCHWARZ_BC_3D(3) = [ ' ', 'D', 'N']

  !-----------------------------------------------------------------------------
  !> Schwarz operator
  !>
  !> In the Schwarz method we consider a rectangular subdomain surrounding an
  !> element located in its center. Curvilinear elements are approximated by
  !> their cuboids. The corresponding subdomains are constructed by adopting
  !> `no` layers of collocation points from the adjoining elements.
  !>
  !> The Schwarz operator is the inverse of the truncated diffusion operator,
  !> which is given in tensor-product form by
  !>
  !>     A  =  c0 M x M x M
  !>        +  c1 M x M x L
  !>        +  c2 M x L x M
  !>        +  c3 L x M x M
  !>
  !> where `M` is the 1D mass matrix and `L` the corresponding stiffness matrix.
  !> These operators are normalized to unit mesh spacing and, thus, depend only
  !> on the following parameters
  !>
  !>   -  polynomial order and, possibly, further discretization parameters
  !>   -  number of overlapped points `no`
  !>   -  Helmholtz and diffusion coefficients
  !>   -  boundary conditions.
  !>
  !> In the case with spectral vanishing viscosity (SVV), the stiffness matrices
  !> are defined as
  !>
  !>     L1 = A1 / (nu + nu_svv)
  !>
  !> where `A1` is the 1D element diffusion matrix for `dx=1` etc.
  !>
  !> The effect of element extensions `dx` is incorporated into the coefficients
  !>
  !>     c0 = dx(1) * dx(2) * dx(3) * lambda
  !>     c1 = dx(2) * dx(3) / dx(1) * (nu + nu_svv)
  !>     c2 = dx(3) * dx(1) / dx(2) * (nu + nu_svv)
  !>     c3 = dx(1) * dx(2) / dx(3) * (nu + nu_svv)
  !>
  !> where `lambda` represents the Helmholtz parameter λ, `nu` the physical
  !> diffusivity ν and `nu_svv` the spectral viscosity amplitude νˢ.
  !>
  !> The element-boundary configuration describes the conditions at the element
  !> faces:
  !>
  !>   *  In the standard configuration, the element is completely enclosed by
  !>      adjoining elements and, hence, every everywhere coated by the layer
  !>      of overlapped points.
  !>
  !>   *  In boundary configurations, one or more element faces coincide with
  !>      the boundary of the computational domain. At those faces, the
  !>      Helmholtz operator is modified according to the boundary conditions,
  !>      and no exterior points are adopted to the Schwarz subdomain.
  !>
  !> Considering interior (I), Dirichlet (D) and Neumann (N) faces, 9 different
  !> configurations have to be distinguished in each coordinate direction:
  !>
  !>    1.  I-I
  !>    2.  D-I
  !>    3.  N-I
  !>    4.  I-D
  !>    5.  D-D
  !>    6.  N-D
  !>    7.  I-N
  !>    8.  D-N
  !>    9.  N-N
  !>
  !> For computing the inverse operator, a generalized 1D eigenvalue problem is
  !> solved for every configuration in each coordinate direction, yielding the
  !> matrices of right eigenvectors `S1`, `S2`, `S3` and the diagonal matrices
  !> of eigenvalues `Λ1`, `Λ2`, `Λ3` such that
  !>
  !>     S1ᵀ L1 S1 = Λ1
  !>     S1ᵀ M1 S1 = I1
  !>
  !> where `I1` ist the matching unit matrix etc.
  !> The inverse Helmholtz operator can be expressed in the tensor-product form
  !>
  !>     A⁻¹  =  (S3 x S2 x S1 D⁻¹ (S3ᵀ x S2ᵀ x S1ᵀ)
  !>
  !> with the diagonal matrix
  !>
  !>     D  =  c0 I3 x I2 x I1
  !>        +  c1 I3 x I2 x Λ1
  !>        +  c2 I3 x Λ2 x I1
  !>        +  c3 Λ3 x I2 x I1
  !>
  !> Before assembling the global correction to an approximate solution, the
  !> subdomain correction is weighted according to
  !>
  !>     Δu = W (A⁻¹ r)
  !>
  !> The weights form a diagonal tensor-product matrix
  !>
  !>     W = W3 x W2 x W1
  !>
  !> with 1D distributions `W1`, `W2`, `W3` depending on the element-boundary
  !> configuration.

  type, abstract :: DG_SchwarzOperator_3D

    integer :: po = -1     !< polynomial order
    integer :: no = -1     !< number of overlapped layers
    integer :: nc = -1     !< number of 1D configurations
    logical :: restrictive !< T/F for ex/including neighbor results

    integer, allocatable :: cfg(:,:) !< subdomain configurations (nc,*)

  contains

    procedure, nopass :: ConfigurationID

    procedure(InitOperators), deferred :: InitOperators

    generic :: InitSubdomains => InitSubdomains_CI, &
                                 InitSubdomains_CS, &
                                 InitSubdomains_VI

    procedure(InitSubdomains_CI), deferred :: InitSubdomains_CI
    procedure(InitSubdomains_CS), deferred :: InitSubdomains_CS
    procedure(InitSubdomains_VI), deferred :: InitSubdomains_VI
    ! INTEL bug prevents declaring InitSubdomains_CI/CS/VI private :(

  end type DG_SchwarzOperator_3D

  !=============================================================================
  ! Interfaces to deferred procedures

  abstract interface

    !---------------------------------------------------------------------------
    !> Build the 1D eigenvalues, eigenvectors and weights

    subroutine InitOperators(this, eop, no, weighting, svv)
      import
      class(DG_SchwarzOperator_3D), intent(inout) :: this
      class(DG_ElementOperators_1D), intent(in) :: eop !< DG element operators
      integer,   intent(in) :: no        !< number of overlapped points
      integer,   intent(in) :: weighting !< weighting method {0,1,3,5,7,9}
      real(RNP), intent(in), optional :: svv !< ratio νˢ/(ν + νˢ)
    end subroutine InitOperators

    !---------------------------------------------------------------------------
    !> Set subdomain configurations and eigenvalues: const isotropic w/o SVV

    subroutine InitSubdomains_CI(this, mesh, lambda, nu, bc)
      import
      class(DG_SchwarzOperator_3D), intent(inout) :: this
      class(Mesh_3D),  intent(in) :: mesh   !< mesh partition
      real(RNP),       intent(in) :: lambda !< Helmholtz parameter
      real(RNP),       intent(in) :: nu     !< physical diffusivity
      character,       intent(in) :: bc(:)  !< BC {'D','N','P'}
    end subroutine InitSubdomains_CI

    !---------------------------------------------------------------------------
    !> Set subdomain configurations and eigenvalues: const isotropic with SVV

    subroutine InitSubdomains_CS(this, mesh, lambda, nu, nu_svv, bc)
      import
      class(DG_SchwarzOperator_3D), intent(inout) :: this
      class(Mesh_3D), intent(in) :: mesh   !< mesh partition
      real(RNP),      intent(in) :: lambda !< Helmholtz parameter
      real(RNP),      intent(in) :: nu     !< physical diffusivity
      real(RNP),      intent(in) :: nu_svv !< spectral diffusivity
      character,      intent(in) :: bc(:)  !< BC {'D','N','P'}
    end subroutine InitSubdomains_CS

    !---------------------------------------------------------------------------
    !> Set subdomain configurations and inverse 3D eigenvals: variable isotropic

    subroutine InitSubdomains_VI(this, eop, mesh, lambda, nu, bc)
      import
      class(DG_SchwarzOperator_3D), intent(inout) :: this
      class(DG_ElementOperators_1D), intent(in) :: eop !< DG element operators
      class(Mesh_3D), intent(in) :: mesh               !< mesh partition
      real(RNP),      intent(in) :: lambda             !< Helmholtz parameter
      real(RNP),      intent(in) :: nu(0:,0:,0:,:)     !< physical diffusivity
      character,      intent(in) :: bc(:)              !< BC {'D','N','P'}
    end subroutine InitSubdomains_VI

  end interface

contains

  !=============================================================================
  ! Utilities

  !-----------------------------------------------------------------------------
  !> Returns the 1D subdomain configuration ID corresponding to the given BCs

  pure integer function ConfigurationID(bc) result(cfg)
    character, intent(in) :: bc(2) !< left/right boundary types {' ','D','N','P'}

    integer :: i1, i2

    select case(bc(1))
    case('D')
      i1 = 2
    case('N')
      i1 = 3
    case default ! ' ' and 'P'
      i1 = 1
    end select

    select case(bc(2))
    case('D')
      i2 = 2
    case('N')
      i2 = 3
    case default ! ' ' and 'P'
      i2 = 1
    end select

    cfg = i1 + (i2 - 1) * size(DG_SCHWARZ_BC_3D)

  end function ConfigurationID

  !=============================================================================

end module DG__Schwarz_Operator__3D__Base

!> summary:  Schwarz operator for elliptic equations
!> author:   Joerg Stiller
!> date:     2018/12/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Schwarz operator for elliptic equations
!===============================================================================
module CART__Schwarz_Operator
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use Standard_Operators_1D
  use IP_Element_Operators_1D
  use XMPI
  use CART__Mesh_Partition

  implicit none
  private

  public :: SchwarzOperator3D
  public :: SchwarzOptions3D

  ! ConfigurationID is made public only to circumvent an error of GCC 8.2 which
  ! makes private entities invisible to submodules :(
  public :: ConfigurationID

  !-----------------------------------------------------------------------------
  !> Schwarz operator
  !>
  !> In the Schwarz method we consider a rectangular subdomain surrounding an
  !> element located in its center. The subdomain is constructed by adopting a
  !> layer of collocation points from the adjoining elements. The thickness of
  !> this layer is a directional property, which depends on two parameters:
  !> the relative thickness `delta` and the minimal number of overlapped points
  !> `no_min`. Typically, the thickness assumes a value between 0 and 1, though
  !> a negative value can be chosen for restricting the subdomain to the element
  !> alone.
  !>
  !> The Schwarz operator is the inverse of the truncated elliptic operator,
  !> which is given in tensor-product form by
  !>
  !>       A  =  c0 M3 x M2 x M1
  !>          +  c1 M3 x M2 x L1
  !>          +  c2 M3 x L2 x M1
  !>          +  c3 L3 x M2 x M1
  !>
  !> where `M1`, `M2`, `M3` are the 1D mass matrices and  `L1`, `L2`, `L3` the
  !> corresponding stiffness matrices. These operators are normalized to unit
  !> mesh spacing and, thus, depend only on the following parameters
  !>
  !>   *  polynomial order and, possibly, further discretization parameters
  !>   *  number of overlapped points `no`
  !>   *  Helmholtz and diffusion coefficients
  !>   *  boundary conditions.
  !>
  !> The effect of element extensions `dx` is incorporated into the coefficients
  !>
  !>     c0 = dx(1) * dx(2) * dx(3) * lambda
  !>     c1 = dx(2) * dx(3) / dx(1) * nu
  !>     c2 = dx(3) * dx(1) / dx(2) * nu
  !>     c3 = dx(1) * dx(2) / dx(3) * nu
  !>
  !> `lambda` represents the Helmholtz parameter and `nu` the diffusivity.
  !>
  !> The element-boundary configuration describes the conditions met at the
  !> element faces:
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
  !>     A⁻¹  =  (S3 x S2 x S1) D⁻¹ (S3ᵀ x S2ᵀ x S1ᵀ)
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
  !>
  !> It is worth noting, that eigenvectors, eigenvalues and weights coincide,
  !> if the number of overlapped points is identical in each direction.
  !> To benefit from possible optimizations, this quasi-isotropic case is
  !> indicated by setting component `isotropic` to true.

  !-----------------------------------------------------------------------------
  !> Options for the 3D Schwarz operator

  type SchwarzOptions3D
    real(RNP) :: delta(3)  =  0.125 !< relative overlap
    integer   :: no_min    = -1     !< min overlap in points
    integer   :: weighting =  5     !< weighting method
  contains
    procedure :: Bcast => SchwarzOptions3D_Bcast
  end type SchwarzOptions3D

  type SchwarzOperator3D

    ! 1D eigensystems ..........................................................

    integer :: no(3) = -1                    !< overlapped node layers
    logical :: isotropic                     !< switch to isotropic operator

    integer :: n1 = -1                       !< number of points in direction 1
    integer :: n2 = -1                       !< number of points in direction 2
    integer :: n3 = -1                       !< number of points in direction 3
    integer :: nc = -1                       !< number of 1D configurations

    real(RNP), allocatable :: S1(:,:,:)      !< eigenvectors for direction 1
    real(RNP), allocatable :: S2(:,:,:)      !< eigenvectors for direction 2
    real(RNP), allocatable :: S3(:,:,:)      !< eigenvectors for direction 3

    real(RNP), allocatable :: V1(:,:)        !< eigenvalues for direction 1 (Λ1)
    real(RNP), allocatable :: V2(:,:)        !< eigenvalues for direction 2 (Λ2)
    real(RNP), allocatable :: V3(:,:)        !< eigenvalues for direction 3 (Λ3)

    real(RNP), allocatable :: W1(:,:)        !< weights for direction 1
    real(RNP), allocatable :: W2(:,:)        !< weights for direction 2
    real(RNP), allocatable :: W3(:,:)        !< weights for direction 3

    ! subdomain configuration and inverse eigenvalues ..........................

    integer,   allocatable :: cfg(:,:)       !< subdomain configurations
    real(RNP), allocatable :: D_inv(:,:,:,:) !< subdomain inverse 3D eigenvalues

    ! variable for the SVV case ................................................

    class(SchwarzOptions3D), allocatable :: opt ! options for the Schwarz operator
    real(RNP) :: svv_ratio = ONE !< ν/(ν+νˢ), if this changes the Eigensystem
                                 !! has to be built again

  contains
    private

    generic, public :: Init_SchwarzOperator3D => Init_Base,   Init_CI,         &
                                                 Init_CI_svv, Init_VI
    procedure :: Init_Base
    procedure :: Init_CI
    procedure :: Init_CI_svv
    procedure :: Init_VI

    generic, public :: SetProblem => SetProblem_CI, SetProblem_CI_svv,         &
                                     SetProblem_VI
    procedure :: SetProblem_CI
    procedure :: SetProblem_CI_svv
    procedure :: SetProblem_VI

  end type SchwarzOperator3D

  ! constructor interface
  interface SchwarzOperator3D
    module procedure New_Base
    module procedure New_CI
    module procedure New_CI_svv
    module procedure New_VI
  end interface

  !=============================================================================
  !> Interfaces to procedures for generating or updating Schwarz operators

  interface

    !---------------------------------------------------------------------------
    !> Build 1D eigensystems for IP/DG-SEM

    module subroutine BuildEigensystems_IP(this, eop, delta, no_min, weighting)
      use IP_Element_Operators_1D

      class(SchwarzOperator3D),     intent(inout) :: this !< Schwarz operator
      class(IP_ElementOperators1D), intent(in)    :: eop  !< IP-DG SE operators

      real(RNP),         intent(in) :: delta(3)  !< relative overlap
      integer, optional, intent(in) :: no_min    !< min overlap in points [-1]
      integer, optional, intent(in) :: weighting !< weighting method      [ 5]

    end subroutine BuildEigensystems_IP

    !---------------------------------------------------------------------------
    !> Build subdomains for constant isotropic coefficients

    module subroutine BuildSubdomains_CI(this, mesh, lambda, nu, bc)
      class(SchwarzOperator3D), intent(inout) :: this    !< Schwarz operator
      class(MeshPartition),     intent(in) :: mesh       !< mesh partition
      real(RNP),                intent(in) :: lambda     !< Helmholtz parameter
      real(RNP),                intent(in) :: nu         !< diffusivity
      character,                intent(in) :: bc(:)      !< BC {'D','N'}
    end subroutine BuildSubdomains_CI

    !---------------------------------------------------------------------------
    !> Build subdomains for variable isotropic coefficients

    module subroutine BuildSubdomains_VI(this, eop, mesh, lambda, nu, bc)
      class(SchwarzOperator3D),   intent(inout) :: this     !< Schwarz operator
      class(StandardOperators1D), intent(in) :: eop         !< 1D SE operators
      class(MeshPartition),       intent(in) :: mesh        !< mesh partition
      real(RNP),                  intent(in) :: lambda      !< Helmholtz parameter
      real(RNP),                  intent(in) :: nu(:,:,:,:) !< diffusivity
      character,                  intent(in) :: bc(:)       !< BC {'D','N'}
    end subroutine BuildSubdomains_VI

  end interface

  !=============================================================================
  ! Private module variables

  character, parameter :: boundary_type(3)   = [ ' ', 'D', 'N']
  integer,   parameter :: num_boundary_types = size(boundary_type)

contains

!===============================================================================
! SchwarzOperator3D :: constructors

!-------------------------------------------------------------------------------
!> New Schwarz operator without problem data

function New_Base(opt, eop) result(this)
  class(SchwarzOptions3D),    intent(in) :: opt  !< Schwarz options
  class(StandardOperators1D), intent(in) :: eop  !< 1D standard SE ops

  type(SchwarzOperator3D) :: this

  call Init_Base(this, opt, eop)

end function New_Base

!-------------------------------------------------------------------------------
!> New Schwarz operator with constant isotropic diffusivity

function New_CI(opt, eop, mesh, lambda, nu, bc) result(this)
  class(SchwarzOptions3D),    intent(in) :: opt    !< Schwarz options
  class(StandardOperators1D), intent(in) :: eop    !< 1D standard SE ops
  class(MeshPartition),       intent(in) :: mesh   !< mesh partition
  real(RNP),                  intent(in) :: lambda !< Helmholtz parameter
  real(RNP),                  intent(in) :: nu     !< diffusivity
  character,                  intent(in) :: bc(:)  !< BC {'D','N'}

  type(SchwarzOperator3D) :: this

  call Init_CI(this, opt, eop, mesh, lambda, nu, bc)

end function New_CI

!-------------------------------------------------------------------------------
!> New Schwarz operator with a combination of constant isotropic diffusivity
!> and constant isotropic spectral diffusivity

function New_CI_svv(opt, eop, mesh, lambda, nu, nu_svv, bc) result(this)
  class(SchwarzOptions3D),    intent(in) :: opt    !< Schwarz options
  class(StandardOperators1D), intent(in) :: eop    !< 1D standard SE ops
  class(MeshPartition),       intent(in) :: mesh   !< mesh partition
  real(RNP),                  intent(in) :: lambda !< Helmholtz parameter
  real(RNP),                  intent(in) :: nu     !< diffusivity
  real(RNP),                  intent(in) :: nu_svv !< spectral diffusivity
  character,                  intent(in) :: bc(:)  !< BC {'D','N'}

  type(SchwarzOperator3D) :: this

  this % opt = opt
  call Init_CI_svv(this, opt, eop, mesh, lambda, nu, nu_svv, bc)

end function New_CI_svv

!-------------------------------------------------------------------------------
!> New Schwarz operator with variable isotropic diffusivity

function New_VI(opt, eop, mesh, lambda, nu, bc) result(this)
  class(SchwarzOptions3D),    intent(in) :: opt         !< Schwarz options
  class(StandardOperators1D), intent(in) :: eop         !< 1D standard SE ops
  class(MeshPartition),       intent(in) :: mesh        !< mesh partition
  real(RNP),                  intent(in) :: lambda      !< Helmholtz parameter
  real(RNP),                  intent(in) :: nu(:,:,:,:) !< diffusivity
  character,                  intent(in) :: bc(:)       !< BC {'D','N'}

  type(SchwarzOperator3D) :: this

  call Init_VI(this, opt, eop, mesh, lambda, nu, bc)

end function New_VI

!===============================================================================
! SchwarzOperator3D :: initialization

!-------------------------------------------------------------------------------
!> Initialize Schwarz operator without problem data

subroutine Init_Base(this, opt, eop)
  class(SchwarzOperator3D),   intent(inout) :: this !< Schwarz operator
  class(SchwarzOptions3D),    intent(in)    :: opt  !< Schwarz options
  class(StandardOperators1D), intent(in)    :: eop  !< 1D standard SE ops

  select type(eop)
  class is(IP_ElementOperators1D)
    call BuildEigensystems_IP(this, eop, opt%delta, opt%no_min, opt%weighting)
  end select

end subroutine Init_Base

!-------------------------------------------------------------------------------
!> Initialize Schwarz operator with constant isotropic diffusivity

subroutine Init_CI(this, opt, eop, mesh, lambda, nu, bc)

  class(SchwarzOperator3D),   intent(inout) :: this !< Schwarz operator
  class(SchwarzOptions3D),    intent(in)    :: opt  !< Schwarz options
  class(StandardOperators1D), intent(in)    :: eop  !< 1D standard SE ops
  class(MeshPartition),       intent(in)    :: mesh !< mesh partition

  real(RNP), intent(in) :: lambda    !< Helmholtz parameter
  real(RNP), intent(in) :: nu        !< diffusivity
  character, intent(in) :: bc(:)     !< BC {'D','N'}

  call Init_Base(this, opt, eop)
  call SetProblem_CI(this, mesh, lambda, nu, bc)

end subroutine Init_CI

!-------------------------------------------------------------------------------
!> Initialize Schwarz operator with a combination of a constant isotropic
!> diffusivity and a constant isotropic spectral diffusivity

subroutine Init_CI_svv(this, opt, eop, mesh, lambda, nu, nu_svv, bc)

  class(SchwarzOperator3D),   intent(inout) :: this !< Schwarz operator
  class(SchwarzOptions3D),    intent(in)    :: opt  !< Schwarz options
  class(StandardOperators1D), intent(in)    :: eop  !< 1D standard SE ops
  class(MeshPartition),       intent(in)    :: mesh !< mesh partition

  real(RNP), intent(in) :: lambda    !< Helmholtz parameter
  real(RNP), intent(in) :: nu        !< diffusivity
  real(RNP), intent(in) :: nu_svv    !< spectral diffusivity
  character, intent(in) :: bc(:)     !< BC {'D','N'}

  call Init_Base(this, opt, eop)
  call SetProblem_CI_svv(this, eop, mesh, lambda, nu, nu_svv, bc)

end subroutine Init_CI_svv

!-------------------------------------------------------------------------------
!> Initialize Schwarz operator with variable isotropic diffusivity

subroutine Init_VI(this, opt, eop, mesh, lambda, nu, bc)

  class(SchwarzOperator3D),   intent(inout) :: this !< Schwarz operator
  class(SchwarzOptions3D),    intent(in)    :: opt  !< Schwarz options
  class(StandardOperators1D), intent(in)    :: eop  !< 1D standard SE ops
  class(MeshPartition),       intent(in)    :: mesh !< mesh partition

  real(RNP), intent(in) :: lambda      !< Helmholtz parameter
  real(RNP), intent(in) :: nu(:,:,:,:) !< diffusivity
  character, intent(in) :: bc(:)       !< BC {'D','N'}

  call Init_Base(this, opt, eop)
  call SetProblem_VI(this, eop, mesh, lambda, nu, bc)

end subroutine Init_VI

!===============================================================================
! SchwarzOperator3D :: SetProblem

!-------------------------------------------------------------------------------
!> (Re)Set problem parameters for constant isotropic viscosity

subroutine SetProblem_CI(this, mesh, lambda, nu, bc)
  class(SchwarzOperator3D), intent(inout) :: this   !< Schwarz operator
  class(MeshPartition),     intent(in)    :: mesh   !< mesh partition
  real(RNP),                intent(in)    :: lambda !< Helmholtz parameter
  real(RNP),                intent(in)    :: nu     !< diffusivity
  character,                intent(in)    :: bc(:)  !< BC {'D','N'}

  call BuildSubdomains_CI(this, mesh, lambda, nu, bc)

end subroutine SetProblem_CI

!-------------------------------------------------------------------------------
!> (Re)Set problem parameters for constant isotropic viscosity

subroutine SetProblem_CI_svv(this, eop, mesh, lambda, nu, nu_svv, bc)
  class(SchwarzOperator3D),   intent(inout) :: this   !< Schwarz operator
  class(StandardOperators1D), intent(in)    :: eop    !< 1D standard SE ops
  class(MeshPartition),       intent(in)    :: mesh   !< mesh partition
  real(RNP),                  intent(in)    :: lambda !< Helmholtz parameter
  real(RNP),                  intent(in)    :: nu     !< diffusivity
  real(RNP),                  intent(in)    :: nu_svv !< diffusivity
  character,                  intent(in)    :: bc(:)  !< BC {'D','N'}

  real(RNP) :: svv_ratio
  
  svv_ratio = nu / (nu + nu_svv)

  ! if the ratio ν/(ν+νˢ) changed compared to the time where the Eigensystem was
  ! built, the Eigensystem has to be built again
  if (this%svv_ratio /= svv_ratio) then
    this%svv_ratio = svv_ratio

    select type(eop)
    class is(IP_ElementOperators1D)
      call BuildEigensystems_IP(this, eop, this%opt%delta, this%opt%no_min,    &
                                this%opt%weighting)
    end select
  end if

  ! use ν+νˢ as prefactor for the D_inv computation as this is the common
  ! denominator due to the SVV ratio ν/(ν+νˢ)
  call BuildSubdomains_CI(this, mesh, lambda, nu + nu_svv, bc)

end subroutine SetProblem_CI_svv

!-------------------------------------------------------------------------------
!> (Re)Set problem parameters for variable isotropic viscosity

subroutine SetProblem_VI(this, eop, mesh, lambda, nu, bc)
  class(SchwarzOperator3D),   intent(inout) :: this        !< Schwarz operator
  class(StandardOperators1D), intent(in)    :: eop         !< 1D standard SE ops
  class(MeshPartition),       intent(in)    :: mesh        !< mesh partition
  real(RNP),                  intent(in)    :: lambda      !< Helmholtz parameter
  real(RNP),                  intent(in)    :: nu(0:,0:,0:,:) !< diffusivity
  character,                  intent(in)    :: bc(:)       !< BC {'D','N'}

  call BuildSubdomains_VI(this, eop, mesh, lambda, nu, bc)

end subroutine SetProblem_VI

!===============================================================================
! SchwarzOptions3D :: Bcast

!-------------------------------------------------------------------------------
!> Extension of MPI_Bcast to objects of type SchwarzOptions3D

subroutine SchwarzOptions3D_Bcast(this, root, comm)
  class(SchwarzOptions3D), intent(inout) :: this
  integer,                 intent(in)    :: root !< rank of broadcast root
  type(MPI_Comm),          intent(in)    :: comm !< MPI communicator

  call XMPI_Bcast( this % delta     , root, comm )
  call XMPI_Bcast( this % no_min    , root, comm )
  call XMPI_Bcast( this % weighting , root, comm )

end subroutine SchwarzOptions3D_Bcast

!===============================================================================
! Utilities

!-------------------------------------------------------------------------------
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

  cfg = i1 + num_boundary_types * (i2 - 1)

end function ConfigurationID

!===============================================================================

end module CART__Schwarz_Operator

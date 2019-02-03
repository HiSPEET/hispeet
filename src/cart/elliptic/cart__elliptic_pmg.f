!> summary:  Polynomial multigrid for use with elliptic solvers
!> author:   Joerg Stiller
!> date:     2019/02/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Polynomial multigrid for use with elliptic solvers
!===============================================================================

module CART__Elliptic_PMG
  use Kind_Parameters, only: RNP
  use XMPI
  use IP_Element_Operators_1D
  use CART__DG_PMG_Transfer_Operators  !!! yet to be revised
  use CART__Mesh_Partition
  use CART__Schwarz_Operator
  use CART__Elliptic_PMG_Level

  implicit none
  private

  public :: PolynomialMultigrid
  public :: PMG_Options

  !-----------------------------------------------------------------------------
  !> Polynomial multigrid operators and related procedures

  type PolynomialMultigrid

    ! mesh, levels and operators
    type(MeshPartition),         pointer     :: mesh        !< mesh partition
    type(PMG_Level),             allocatable :: level(:)    !< levels
    type(PMG_TransferOperators), allocatable :: transfer(:) !< transfer ops

    ! p-MG/CG solver settings
    integer   :: i_max   !< max number iterations (cycles)
    real(RNP) :: r_red   !< min residual reduction
    real(RNP) :: r_max   !< max admissible residual
    real(RNP) :: dr_min  !< termination threshold for Δr

    ! coarse grid solver settings
    character :: solver  !< coarse grid solver, 'C': CG, 'S': Schwarz
    integer   :: i0_max  !< max num of iterations on coarse grid
    real(RNP) :: r0_red  !< min residual reduction on coarse grid

    ! control
    logical   :: monitor !< switch for monitoring

  contains

    generic :: Init_PolynomialMultigrid => Init__IP_CI
    procedure, private :: Init__IP_CI

  end type PolynomialMultigrid

  !-----------------------------------------------------------------------------
  !> Type bundling polynomial multigrid options

  type PMG_Options

    ! levels
    integer   :: po_top    = -1  !< polynomial order at top level
    integer   :: po_bot    =  1  !< polynomial order at bottom4 level
    real(RNP) :: cr        =  2  !< target coarsening ratio
    real(RNP) :: cr_max    =  2  !< max coarsening ratio

    ! p-MG/CG solver settings
    integer   :: i_max     =  1  !< max number iterations (cycles)
    real(RNP) :: r_red     = -1  !< min residual reduction
    real(RNP) :: r_max     = -1  !< max admissible residual
    real(RNP) :: dr_min    = -1  !< termination threshold for Δr

    ! smoothing settings
    integer   :: ns1       =  1  !< num pre-smoothing  steps on top level
    integer   :: ns2       =  1  !< num post-smoothing steps on top level
    integer   :: mvs       =  1  !< multiplier for variable smoothing

    ! coarse grid solver settings
    character :: solver    = 'C' !< coarse grid solver, 'C': CG, 'S': Schwarz
    integer   :: i0_max    =  1  !< max number coarse grid iterations
    real(RNP) :: r0_red    = -1  !< min coarse grid residual reduction

    ! control
    logical   :: monitor   = .false. !< switch for monitoring

  contains

    procedure :: Bcast => PMG_Options_Bcast

  end type PMG_Options

  !=============================================================================
  ! Separate procedures

  interface

    !---------------------------------------------------------------------------
    !> PolynomialMultigrid initialization: IP, constant isotropic diffusivity

    module subroutine Init__IP_CI( this, mesh, lambda, nu, bc,  &
                                   pmg_opt, ip_opt, schwarz_opt )

      class(PolynomialMultigrid), intent(inout) :: this
      class(MeshPartition), target, intent(in) :: mesh   !< mesh partition
      real(RNP), intent(in) :: lambda                    !< Helmholtz parameter
      real(RNP), intent(in) :: nu                        !< diffusivity
      character, intent(in) :: bc(:)                     !< boundary conditions
      class(PMG_Options), intent(in) :: pmg_opt          !< PMG options
      class(IP_ElementOptions1D), intent(in) :: ip_opt   !< IP/DG options
      class(SchwarzOptions3D), intent(in) :: schwarz_opt !< Schwarz options

    end subroutine Init__IP_CI

    !---------------------------------------------------------------------------
    !> PolynomialMultigrid initialization: IP, variable isotropic diffusivity

    module subroutine Init__IP_VI( this, mesh, lambda, nu, bc,  &
                                   pmg_opt, ip_opt, schwarz_opt )

      class(PolynomialMultigrid), intent(inout) :: this
      class(MeshPartition), target, intent(in) :: mesh   !< mesh partition
      real(RNP), intent(in) :: lambda                    !< Helmholtz parameter
      real(RNP), intent(in) :: nu(0:,0:,0:,:)            !< diffusivity
      character, intent(in) :: bc(:)                     !< boundary conditions
      class(PMG_Options), intent(in) :: pmg_opt          !< PMG options
      class(IP_ElementOptions1D), intent(in) :: ip_opt   !< IP/DG options
      class(SchwarzOptions3D), intent(in) :: schwarz_opt !< Schwarz options

    end subroutine Init__IP_VI

  end interface

contains

!===============================================================================
! PolynomialMultigrid: type-bound procedures

!===============================================================================
! PMG_Options: type-bound procedures

subroutine PMG_Options_Bcast(this, root, comm)
  class(PMG_Options), intent(inout) :: this
  integer,        intent(in) :: root !< rank of broadcast root
  type(MPI_Comm), intent(in) :: comm !< MPI communicator

  type(MPI_Request) :: request(15)
  type(MPI_Status)  :: stat(size(request))
  integer :: n

  n = 1
  call XMPI_Ibcast( this % po_top    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % po_bot    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % cr        , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % cr_max    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % i_max     , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % r_red     , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % r_max     , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % dr_min    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % ns1       , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % ns2       , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % mvs       , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % solver    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % i0_max    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % r0_red    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % monitor   , root, comm, request(n) )

  call MPI_Waitall(n, request, stat)

end subroutine PMG_Options_Bcast

!===============================================================================

end module CART__Elliptic_PMG

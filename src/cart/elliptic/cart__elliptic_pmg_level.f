!> summary:  Polynomial multigrid level for use with elliptic solvers
!> author:   Joerg Stiller
!> date:     2019/01/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Polynomial multigrid level for use with DG elliptic solvers
!===============================================================================

module CART__Elliptic_PMG_Level
  use Kind_Parameters, only: RNP
  use IP_Element_Operators_1D
  use CART__Mesh_Partition
  use CART__Schwarz_Operator
  use CART__Elliptic_Operator
  use CART__Elliptic_Operator_IP

  implicit none
  private

  public :: PMG_Level

  !-----------------------------------------------------------------------------
  !> Polynomial multigrid level

  type PMG_Level

    integer :: po = -1  !< polynomial order
    integer :: ne = -1  !< number of local elements
    logical :: bottom   !< switch for bottom level
    logical :: top      !< switch for top level
    integer :: ns1 = 1  !< number of pre-smoothing steps
    integer :: ns2 = 1  !< number of post-smoothing steps

    class(EllipticOperator3D), allocatable :: eop !< elliptic operator

    real(RNP), allocatable :: u(:,:,:,:) !< solution
    real(RNP), allocatable :: f(:,:,:,:) !< RHS
    real(RNP), allocatable :: v(:,:,:,:) !< residual / correction

  contains

    generic :: Init => Init__IP_CI, Init__IP_VI
    procedure, private :: Init__IP_CI
    procedure, private :: Init__IP_VI

  end type PMG_Level

contains

!===============================================================================
! Initialization of PMG_Level

!-------------------------------------------------------------------------------
!> Initialization of PMG_Level with IP/DG-SEM and constant isotropic diffusivity

subroutine Init__IP_CI( this, bottom, top, ns1, ns2, mesh, lambda, nu, bc, &
                        ip_opt, schwarz_opt )

  ! arguments ..................................................................

  class(PMG_Level), intent(inout) :: this

  logical,   intent(in) :: bottom           !< switch for bottom level
  logical,   intent(in) :: top              !< switch for top level
  integer,   intent(in) :: ns1              !< number of pre-smoothing steps
  integer,   intent(in) :: ns2              !< number of post-smoothing steps

  class(MeshPartition), intent(in) :: mesh  !< mesh partition

  real(RNP), intent(in) :: lambda           !< Helmholtz parameter
  real(RNP), intent(in) :: nu               !< diffusivity
  character, intent(in) :: bc(:)            !< boundary conditions

  class(IP_ElementOptions1D), intent(in) :: ip_opt
  !< options for the IP/DG method, including polynomial order `po` and `penalty`

  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt
  !< options for the Schwarz method

  ! initialization of components ...............................................

  this % po     = ip_opt % po
  this % ne     = mesh   % ne
  this % bottom = bottom
  this % top    = top
  this % ns1    = ns1
  this % ns2    = ns2

  this % eop = EllipticOperator3D_IP(mesh, lambda, nu, bc, ip_opt, schwarz_opt)

  call GetWorkspace(this)

end subroutine Init__IP_CI

!-------------------------------------------------------------------------------
!> Initialization of PMG_Level with IP/DG-SEM and variable isotropic diffusivity

subroutine Init__IP_VI( this, bottom, top, ns1, ns2, mesh, lambda, nu, bc, &
                        ip_opt, schwarz_opt )

  ! arguments ..................................................................

  class(PMG_Level), intent(inout) :: this

  logical,   intent(in) :: bottom          !< switch for bottom level
  logical,   intent(in) :: top             !< switch for top level
  integer,   intent(in) :: ns1             !< number of pre-smoothing steps
  integer,   intent(in) :: ns2             !< number of post-smoothing steps

  class(MeshPartition), intent(in) :: mesh !< mesh partition

  real(RNP), intent(in) :: lambda          !< Helmholtz parameter
  real(RNP), intent(in) :: nu(0:,0:,0:,:)  !< diffusivity
  character, intent(in) :: bc(:)           !< boundary conditions

  class(IP_ElementOptions1D), intent(in) :: ip_opt
  !< options for the IP/DG method, including polynomial order `po` and `penalty`

  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt
  !< options for the Schwarz method

  ! initialization of components ...............................................

  this % po     = ip_opt % po
  this % ne     = mesh   % ne
  this % bottom = bottom
  this % top    = top
  this % ns1    = ns1
  this % ns2    = ns2

  this % eop = EllipticOperator3D_IP(mesh, lambda, nu, bc, ip_opt, schwarz_opt)

  call GetWorkspace(this)

end subroutine Init__IP_VI

!===============================================================================
! Helpers

!-------------------------------------------------------------------------------
!> Allocate workspace

subroutine GetWorkspace(this)
  class(PMG_Level), intent(inout) :: this

  associate(po => this%po, ne => this%ne)

    if (allocated(this % u)) then
      if (any(shape(this % u) /= [po+1, po+1, po+1, ne])) then
        deallocate(this % u)
        deallocate(this % f)
        deallocate(this % v)
      end if
    end if

    if (.not. allocated(this % u)) then
      allocate(this % u(0:po, 0:po, 0:po, ne))
      allocate(this % f, mold = this % u)
      allocate(this % v, mold = this % u)
    end if

  end associate

end subroutine GetWorkspace

!===============================================================================

end module CART__Elliptic_PMG_Level

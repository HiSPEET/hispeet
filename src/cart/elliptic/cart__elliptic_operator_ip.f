!> summary:  3D Cartesian elliptic operator for interior penalty (IP) DG-SEM
!> author:   Joerg Stiller
!> date:     2018/11/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### 3D Cartesian elliptic operator for interior penalty (IP) DG-SEM
!===============================================================================

module CART__Elliptic_Operator_IP
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF
  use Array_Assignments
  use Array_Reductions
  use IP_Element_Operators_1D
  use XMPI, only: XMPI_Bcast
  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Schwarz_Operator
  use CART__Elliptic_Operator

  implicit none
  private

  public :: EllipticOperator3D_IP

  !-----------------------------------------------------------------------------
  !> Type accommodating 3D Cartesian elliptic operators for IP/DG-SEM

  type, extends(EllipticOperator3D) :: EllipticOperator3D_IP

    real(RNP), allocatable :: nu_hat(:,:,:) !< diffusivity on faces

  contains

    generic :: Init_EllipticOperator3D_IP => Init_Base, Init_CI, Init_VI
    procedure, private :: Init_Base
    procedure, private :: Init_CI
    procedure, private :: Init_VI

    ! specific procedures for generic SetProblem
    procedure :: SetProblem_CI  ! should be private
    procedure :: SetProblem_VI  ! but fails with ifort

    procedure :: Apply
    procedure :: BcToRHS
    procedure :: Residual
    procedure :: ConjugateGradients
    procedure :: SchwarzMethod

  end type EllipticOperator3D_IP

  ! constructor interface
  interface EllipticOperator3D_IP
    module procedure New_Base
    module procedure New_CI
    module procedure New_VI
  end interface

  !=============================================================================
  ! Separate procedures

  interface

    !---------------------------------------------------------------------------
    !> Initialize problem with constant isotropic diffusivity

    module subroutine SetProblem_CI(this, lambda, nu, bc)
      class(EllipticOperator3D_IP), intent(inout) :: this
      real(RNP), intent(in) :: lambda !< Helmholtz parameter
      real(RNP), intent(in) :: nu     !< diffusivity
      character, intent(in) :: bc(:)  !< boundary conditions
    end subroutine SetProblem_CI

    !--------------------------------------------------------------------------
    !> Initialize problem with variable isotropic diffusivity

    module subroutine SetProblem_VI(this, lambda, nu, bc)
      class(EllipticOperator3D_IP), intent(inout) :: this
      real(RNP), intent(in) :: lambda         !< Helmholtz parameter
      real(RNP), intent(in) :: nu(0:,0:,0:,:) !< diffusivity
      character, intent(in) :: bc(:)          !< boundary conditions
    end subroutine SetProblem_VI

    !---------------------------------------------------------------------------
    !> Application of the IP/DG elliptic operator

    module subroutine Apply(this, u, v)
      class(EllipticOperator3D_IP), intent(in) :: this
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result
    end subroutine Apply

    !---------------------------------------------------------------------------
    !> Adds the boundary contributions of the right hand side

    module subroutine BcToRHS(this, bv, c, f)
      class(EllipticOperator3D_IP), intent(in)    :: this
      type(BoundaryVariable),       intent(in)    :: bv(:)      !< BC
      integer,            optional, intent(in)    :: c          !< component [1]
      real(RNP),                    intent(inout) :: f(:,:,:,:) !< RHS
    end subroutine BcToRHS

    !--------------------------------------------------------------------------
    !> Computes the residual to given approximation

    module subroutine Residual(this, u, f, r)
      class(EllipticOperator3D_IP), intent(in) :: this
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(in)  :: f(0:,0:,0:,:)  !< right hand side
      real(RNP), intent(out) :: r(0:,0:,0:,:)  !< result
    end subroutine Residual

    !---------------------------------------------------------------------------
    !> Element-centered overlapping Schwarz method with constant coefficients

    module subroutine SchwarzMethod(this, u, f, i_max, r_red, r_max, ni)
      class(EllipticOperator3D_IP), intent(in) :: this
      real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
      integer,   intent(in)    :: i_max          !< max num iterations
      real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
      real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
      integer,   optional, intent(out) :: ni     !< exec num iterations
    end subroutine SchwarzMethod

  end interface

contains

!===============================================================================
! Constructors

!-------------------------------------------------------------------------------
!> New EllipticOperator3D_IP with no problem data

function New_Base(mesh, ip_opt, schwarz_opt) result(this)
  class(MeshPartition), target,      intent(in) :: mesh
  class(IP_ElementOptions1D),        intent(in) :: ip_opt
  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt

  type(EllipticOperator3D_IP) :: this

  call Init_Base(this, mesh, ip_opt, schwarz_opt)

end function New_Base

!-------------------------------------------------------------------------------
!> New EllipticOperator3D_IP with constant isotropic diffusivity

function New_CI(mesh, lambda, nu, bc, ip_opt, schwarz_opt) result(this)
  class(MeshPartition), target,      intent(in) :: mesh
  real(RNP),                         intent(in) :: lambda
  real(RNP),                         intent(in) :: nu
  character,                         intent(in) :: bc(:)
  class(IP_ElementOptions1D),        intent(in) :: ip_opt
  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt

  type(EllipticOperator3D_IP) :: this

  call Init_CI(this, mesh, lambda, nu, bc, ip_opt, schwarz_opt)

end function New_CI

!-------------------------------------------------------------------------------
!> New EllipticOperator3D_IP with variable isotropic diffusivity

function New_VI(mesh, lambda, nu, bc, ip_opt, schwarz_opt) result(this)
  class(MeshPartition), target,      intent(in) :: mesh
  real(RNP),                         intent(in) :: lambda
  real(RNP),                         intent(in) :: nu(0:,0:,0:,:)
  character,                         intent(in) :: bc(:)
  class(IP_ElementOptions1D),        intent(in) :: ip_opt
  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt

  type(EllipticOperator3D_IP) :: this

  call Init_VI(this, mesh, lambda, nu, bc, ip_opt, schwarz_opt)

end function New_VI

!===============================================================================
! Type-bound procedures

!-------------------------------------------------------------------------------
!> Initialize operator with no problem data

subroutine Init_Base(this, mesh, ip_opt, schwarz_opt)
  class(EllipticOperator3D_IP),      intent(inout) :: this
  class(MeshPartition), target,      intent(in)    :: mesh
  class(IP_ElementOptions1D),        intent(in)    :: ip_opt
  class(SchwarzOptions3D), optional, intent(in)    :: schwarz_opt

  this % mesh => mesh

  allocate(this % eop, source = IP_ElementOperators1D(ip_opt))

  if (present(schwarz_opt)) then
    this % schwarz = SchwarzOperator3D(schwarz_opt, this%eop)
  end if

end subroutine Init_Base

!-------------------------------------------------------------------------------
!> Initialize operator with constant isotropic diffusivity

subroutine Init_CI(this, mesh, lambda, nu, bc, ip_opt, schwarz_opt)
  class(EllipticOperator3D_IP),      intent(inout) :: this
  class(MeshPartition), target,      intent(in)    :: mesh
  real(RNP),                         intent(in)    :: lambda
  real(RNP),                         intent(in)    :: nu
  character,                         intent(in)    :: bc(:)
  class(IP_ElementOptions1D),        intent(in)    :: ip_opt
  class(SchwarzOptions3D), optional, intent(in)    :: schwarz_opt

  call Init_Base(this, mesh, ip_opt, schwarz_opt)
  call SetProblem_CI(this, lambda, nu, bc)

end subroutine Init_CI

!-------------------------------------------------------------------------------
!> Initialize operator with variable isotropic diffusivity

subroutine Init_VI(this, mesh, lambda, nu, bc, ip_opt, schwarz_opt)
  class(EllipticOperator3D_IP),      intent(inout) :: this
  class(MeshPartition), target,      intent(in)    :: mesh
  real(RNP),                         intent(in)    :: lambda
  real(RNP),                         intent(in)    :: nu(0:,0:,0:,:)
  character,                         intent(in)    :: bc(:)
  class(IP_ElementOptions1D),        intent(in)    :: ip_opt
  class(SchwarzOptions3D), optional, intent(in)    :: schwarz_opt

  call Init_Base(this, mesh, ip_opt, schwarz_opt)
  call SetProblem_VI(this, lambda, nu, bc)

end subroutine Init_VI

!-------------------------------------------------------------------------------
!> Conjugate gradient method

subroutine ConjugateGradients(this, u, f, i_max, r_red, r_max, ni)
  class(EllipticOperator3D_IP), intent(in) :: this
  real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
  integer,   intent(in)    :: i_max          !< max num iterations
  real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
  real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
  integer,   optional, intent(out) :: ni     !< exec num iterations

  ! local variables ............................................................

  real(RNP), dimension(:,:,:,:), allocatable, save :: g, r, p, q
  real(RNP), save :: rr_term
  logical  , save :: converged

  real(RNP) :: alpha, pq, rr, rr_old
  integer   :: i

  ! initialization .............................................................

  associate(mesh => this%mesh)

    ! work space
    !$omp single
    allocate(g, mold = u)
    allocate(r, mold = u)
    allocate(p, mold = u)
    allocate(q, mold = u)
    !$omp end single

    !$acc data create(g,r,p,q) present(u,f)

    ! RHS
    call SetArray(g, f)

    ! calibrate RHS of singular problem
    if (abs(this%lambda) < epsilon(ONE) .and. all(this%bc /= 'D')) then
      call CalibrateArray(g, mesh%comm)
    end if

    ! initial residual .........................................................

    call this % Residual(u, g, r)
    call SetArray(p, r)

    rr = ScalarProduct(r, r, mesh%comm)

    !$omp single
    if (present(r_red)) then
      rr_term  = max(ZERO, sqrt(rr) * r_red)**2
      if (present(r_max)) then
        rr_term = max(rr_term, max(ZERO, r_max)**2)
      end if
    else
      rr_term = 0
    end if
    !$omp end single

    rr_old = 0

    ! iteration ................................................................

    do i = 1, i_max

      ! termination check  . . . . . . . . . . . . . . . . . . . . . . . . . . .

      ! MPI master decides about termination
      !$omp master
      if (mesh%part == 0) then
        converged = rr <= rr_term
      end if
      rr_old = rr
      call XMPI_Bcast(converged, root=0, comm=mesh%comm)
      !$omp end master
      !$omp barrier

      if (converged) exit

      ! next iteration . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      call this % Apply(p, q)

      pq = ScalarProduct(p, q, mesh%comm)
      alpha = rr_old / pq

      call MergeArrays(ONE, u,  alpha, p)

      if (mod(i,50) == 0) then
        ! compute true residual to get rid of round-off errors
        call this % Residual(u, g, r)
      else
        call MergeArrays(ONE, r, -alpha, q)
      end if

      rr = ScalarProduct(r, r, mesh%comm)

      call  MergeArrays(rr/rr_old, p, ONE, r)

    end do

    if (present(ni)) ni = i - 1

    !$acc end data

    !$omp barrier
    !$omp master
    deallocate(g, r, p, q)
    !$omp end master

  end associate

end subroutine ConjugateGradients

!===============================================================================

end module CART__Elliptic_Operator_IP

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

    type(IP_ElementOperators1D) :: eop           !< 1D operators for IP-DG-SEM
    real(RNP), allocatable      :: nu_hat(:,:,:) !< diffusivity on faces

  contains

    generic :: New => New_CI, New_VI
    procedure, private :: New_CI
    procedure, private :: New_VI

    procedure :: Apply
    procedure :: BcToRHS
    procedure :: Residual
    procedure :: ConjugateGradients
    procedure :: SchwarzMethod

  end type EllipticOperator3D_IP

  !=============================================================================
  ! Separate procedures

  interface

    !---------------------------------------------------------------------------
    !> New operator with constant isotropic diffusivity

    module subroutine New_CI(this, mesh, lambda, nu, bc, po, penalty, &
                             schwarz_opt)

      class(EllipticOperator3D_IP), intent(inout) :: this

      class(MeshPartition), target, &
                            intent(in) :: mesh        !< mesh partition
      real(RNP),            intent(in) :: lambda      !< Helmholtz parameter
      real(RNP),            intent(in) :: nu          !< diffusivity
      character,            intent(in) :: bc(:)       !< boundary conditions
      integer,              intent(in) :: po          !< polynomial order
      real(RNP),  optional, intent(in) :: penalty     !< penalty parameter

      class(SchwarzOptions3D), &
                  optional, intent(in) :: schwarz_opt !< Schwarz options

    end subroutine New_CI

    !---------------------------------------------------------------------------
    !> New operator with variable isotropic diffusivity

    module subroutine New_VI(this, mesh, lambda, nu, bc, po, penalty, &
                             schwarz_opt)

      class(EllipticOperator3D_IP), intent(inout) :: this

      class(MeshPartition), target, &
                            intent(in) :: mesh           !< mesh partition
      real(RNP),            intent(in) :: lambda         !< Helmholtz parameter
      real(RNP),            intent(in) :: nu(0:,0:,0:,:) !< diffusivity
      character,            intent(in) :: bc(:)          !< boundary conditions
      integer,              intent(in) :: po             !< polynomial order
      real(RNP),  optional, intent(in) :: penalty        !< penalty parameter

      class(SchwarzOptions3D), &
                  optional, intent(in) :: schwarz_opt    !< Schwarz options

    end subroutine New_VI

    !---------------------------------------------------------------------------
    !> Application of the IP/DG elliptic operator

    module subroutine Apply(this, u, v)
      class(EllipticOperator3D_IP), intent(in) :: this
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result
    end subroutine Apply

    !---------------------------------------------------------------------------
    !> Adds the boundary contributions of the right hand side

    module subroutine BcToRHS(this, bv, f)
      class(EllipticOperator3D_IP), intent(in)    :: this
      type(BoundaryVariable),       intent(in)    :: bv(:)      !< BC
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

!===============================================================================

contains

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
    call AssignArray(g, f)

    ! calibrate RHS of singular problem
    if (abs(this%lambda) < epsilon(ONE) .and. all(this%bc /= 'D')) then
      call CalibrateArray(g, mesh%comm)
    end if

    ! initial residual .........................................................

    call this % Residual(u, g, r)
    call AssignArray(p, r)

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

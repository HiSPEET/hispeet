!> summary:  3D Cartesian elliptic operator for interior penalty (IP) DG-SEM
!> author:   Joerg Stiller, Gustav Tschirschnitz
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
  use DG__Element_Operators__1D
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

    generic :: Init_EllipticOperator3D_IP => Init_Base, Init_CI, Init_CI_svv,  &
                                             Init_VI


    procedure, private :: Init_Base
    procedure, private :: Init_CI
    procedure, private :: Init_CI_svv
    procedure, private :: Init_VI

    ! specific procedures for generic SetProblem
    procedure :: SetProblem_CI     ! should be private
    procedure :: SetProblem_CI_svv ! but fails with ifort
    procedure :: SetProblem_VI
    procedure :: Apply
    procedure :: BcToRHS
    procedure :: Residual
    procedure :: ConjugateGradients
    procedure :: SchwarzMethod
    procedure :: SchwarzPreConjugateGradients

  end type EllipticOperator3D_IP

  ! constructor interface
  interface EllipticOperator3D_IP
    module procedure New_Base
    module procedure New_CI
    module procedure New_CI_svv
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

    !---------------------------------------------------------------------------
    !> Initialize problem with a combination of constant isotropic diffusivity
    !> and constant isotropic spectral diffusivity

    module subroutine SetProblem_CI_svv(this, lambda, nu, nu_svv, bc)
      class(EllipticOperator3D_IP), intent(inout) :: this
      real(RNP), intent(in) :: lambda !< Helmholtz parameter
      real(RNP), intent(in) :: nu     !< diffusivity
      real(RNP), intent(in) :: nu_svv !< spectral diffusivity
      character, intent(in) :: bc(:)  !< boundary conditions
    end subroutine SetProblem_CI_svv

    !---------------------------------------------------------------------------
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

    !---------------------------------------------------------------------------
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
  class(DG_ElementOptions_1D),        intent(in) :: ip_opt
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
  class(DG_ElementOptions_1D),        intent(in) :: ip_opt
  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt

  type(EllipticOperator3D_IP) :: this

  call Init_CI(this, mesh, lambda, nu, bc, ip_opt, schwarz_opt)

end function New_CI

!-------------------------------------------------------------------------------
!> New EllipticOperator3D_IP with a combination of constant isotropic
!> diffusivity and constant isotropic spectral diffusivity

function New_CI_svv(mesh, lambda, nu, nu_svv, bc, ip_opt, schwarz_opt) &
    result(this)

  class(MeshPartition), target,      intent(in) :: mesh
  real(RNP),                         intent(in) :: lambda
  real(RNP),                         intent(in) :: nu
  real(RNP),                         intent(in) :: nu_svv
  character,                         intent(in) :: bc(:)
  class(DG_ElementOptions_1D),        intent(in) :: ip_opt
  class(SchwarzOptions3D), optional, intent(in) :: schwarz_opt

  type(EllipticOperator3D_IP) :: this

  call Init_CI_svv(this, mesh, lambda, nu, nu_svv, bc, ip_opt, schwarz_opt)

end function New_CI_svv

!-------------------------------------------------------------------------------
!> New EllipticOperator3D_IP with variable isotropic diffusivity

function New_VI(mesh, lambda, nu, bc, ip_opt, schwarz_opt) result(this)
  class(MeshPartition), target,      intent(in) :: mesh
  real(RNP),                         intent(in) :: lambda
  real(RNP),                         intent(in) :: nu(0:,0:,0:,:)
  character,                         intent(in) :: bc(:)
  class(DG_ElementOptions_1D),        intent(in) :: ip_opt
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
  class(DG_ElementOptions_1D),       intent(in)    :: ip_opt
  class(SchwarzOptions3D), optional, intent(in)    :: schwarz_opt

  this % mesh => mesh

  allocate(this % eop, source = DG_ElementOperators_1D(ip_opt))

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
  class(DG_ElementOptions_1D),        intent(in)    :: ip_opt
  class(SchwarzOptions3D), optional, intent(in)    :: schwarz_opt

  call Init_Base(this, mesh, ip_opt, schwarz_opt)
  call SetProblem_CI(this, lambda, nu, bc)

end subroutine Init_CI

!-------------------------------------------------------------------------------
!> Initialize operator with a combination of constant isotropic diffusivity and
!> constant isotropic spectral diffusivity

subroutine Init_CI_svv(this, mesh, lambda, nu, nu_svv, bc, ip_opt, schwarz_opt)
  class(EllipticOperator3D_IP),      intent(inout) :: this
  class(MeshPartition), target,      intent(in)    :: mesh
  real(RNP),                         intent(in)    :: lambda
  real(RNP),                         intent(in)    :: nu
  real(RNP),                         intent(in)    :: nu_svv
  character,                         intent(in)    :: bc(:)
  class(DG_ElementOptions_1D),        intent(in)    :: ip_opt
  class(SchwarzOptions3D), optional, intent(in)    :: schwarz_opt

  call Init_Base(this, mesh, ip_opt, schwarz_opt)
  call SetProblem_CI_svv(this, lambda, nu, nu_svv, bc)

end subroutine Init_CI_svv

!-------------------------------------------------------------------------------
!> Initialize operator with variable isotropic diffusivity

subroutine Init_VI(this, mesh, lambda, nu, bc, ip_opt, schwarz_opt)
  class(EllipticOperator3D_IP),      intent(inout) :: this
  class(MeshPartition), target,      intent(in)    :: mesh
  real(RNP),                         intent(in)    :: lambda
  real(RNP),                         intent(in)    :: nu(0:,0:,0:,:)
  character,                         intent(in)    :: bc(:)
  class(DG_ElementOptions_1D),        intent(in)    :: ip_opt
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

  real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
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
      call XMPI_Bcast(converged, root=0, comm=mesh%comm)
      !$omp end master
      !$omp barrier

      if (converged) exit

      rr_old = rr

      ! next iteration . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      call this % Apply(p, q)

      pq = ScalarProduct(p, q, mesh%comm)
!print '(A,G0)', 'rr_old = ', rr_old
!print '(A,G0)', 'pq1 = ', pq
      pq = sign(max(abs(pq),eps), pq)
!print '(A,G0)', 'pq2 = ', pq
      alpha = rr_old / pq
!print '(A,G0)', 'alpha = ', alpha
      call MergeArrays(ONE, u,  alpha, p)  !?! nowait in case mod(i,50) /= 0

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

    !$omp master
    deallocate(g, r, p, q)
    !$omp end master

  end associate

end subroutine ConjugateGradients

!-------------------------------------------------------------------------------
!> Schwarz-preconditioned conjugate gradient method

subroutine SchwarzPreConjugateGradients(this, u, f, i_max, r_red, r_max, ni)
  class(EllipticOperator3D_IP), intent(in) :: this
  real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
  integer,   intent(in)    :: i_max          !< max num iterations
  real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
  real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
  integer,   optional, intent(out) :: ni     !< exec num iterations

  ! local variables ............................................................

  real(RNP), dimension(:,:,:,:), allocatable, save :: g, r, p, q, s, z
  real(RNP), save :: rr_term
  logical  , save :: converged

  real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
  logical   :: check_convergence, singular
  real(RNP) :: alpha, beta, delta, rr
  integer   :: i, i_max_

  ! initialization .............................................................

  check_convergence = present(r_red) .or. present(r_max)
  singular = abs(this%lambda) < epsilon(ONE) .and. all(this%bc /= 'D')

  associate(mesh => this%mesh)

    ! work space
    !$omp single
    allocate(g, mold = u)
    allocate(p, mold = u)
    allocate(q, mold = u)
    allocate(r, mold = u)
    allocate(s, mold = u)
    allocate(z, mold = u)
    !$omp end single

    ! RHS
    call SetArray(g, f)

    ! calibrate RHS of singular problem
    if (singular) then
      call CalibrateArray(g, mesh%comm)
    end if

    ! initial residual .........................................................

    call this % Residual(u, g, r)
    call SetArray(p, r)

    ! termination conditions
    if (check_convergence) then
      rr = ScalarProduct(r, r, mesh%comm)
      !$omp master
      if (present(r_red)) then
        rr_term  = max(ZERO, sqrt(rr) * r_red)**2
      else
        rr_term = 0
      end if
      if (present(r_max)) then
        rr_term = max(rr_term, max(ZERO, r_max)**2)
      end if
      converged = rr <= rr_term
      call XMPI_Bcast(converged, root=0, comm=mesh%comm)
      !$omp end master
      !$omp barrier
    else
      !$omp master
      converged = .false.
      !$omp end master
      !$omp barrier
    end if

    if (converged) then
      i_max_ = 0
      i      = 0
    else
      i_max_ = i_max
    end if

    ! iteration ................................................................

    do i = 1, i_max_

      ! Schwarz preconditioner
      call SetArray(z, ZERO)
      call this % SchwarzMethod(z, r, 1)

      ! set/update search vector
      if (i == 1) then
        if (singular) then
          call CalibrateArray(z, mesh%comm)
        end if
        call SetArray(p, z)                               ! p = z
      else
        call SetArray(q, r)                               ! q = r
        call MergeArrays(ONE, q, -ONE, s)                 ! q = r - s
        beta = ScalarProduct(q, z, mesh%comm) / delta
        call MergeArrays(beta, p, ONE, z)                 ! p = beta p + z
      end if

      ! save old residual
      call SetArray(s, r)

      ! correction
      call this % Apply(p, q)
      delta = ScalarProduct(r, z, mesh%comm)
      alpha = delta / ScalarProduct(p, q, mesh%comm)
      call MergeArrays(ONE, u,  alpha, p)                 ! u = u + alpha p
      call MergeArrays(ONE, r, -alpha, q)                 ! r = r - alpha q

      if (check_convergence) then
        rr = ScalarProduct(r, r, mesh%comm)
        !$omp master
        converged = rr <= rr_term
        call XMPI_Bcast(converged, root=0, comm=mesh%comm)
        !$omp end master
        !$omp barrier
      end if

      if (converged .or. i == i_max_) exit

    end do

    ! finalization .............................................................

    if (present(ni)) ni = i

    !$omp master
    deallocate(g, p, q, r, s, z)
    !$omp end master

  end associate

end subroutine SchwarzPreConjugateGradients

!===============================================================================

end module CART__Elliptic_Operator_IP

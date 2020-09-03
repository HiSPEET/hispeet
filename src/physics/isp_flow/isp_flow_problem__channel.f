!> summary:  3D turbulent channel flow
!> author:   Joerg Stiller
!> date:     2020/08/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> This module provides an environment for turbulent plane channel flow.
!> ...
!>
!>
!>
!> References
!>
!>   *  J Kim, P Moin, RD Moser, Turbulence statistics in fully developed
!>      turbulent channel flow at low Reynolds number. JFM 177:133-166, 1987
!>
!>   *  KS Nelson, OB Fringer, Reducing spin-up time for simulations of
!>      turbulent channel flow. Phys Fluids 29:105101, 2017
!>
!===============================================================================

module ISP_Flow_Problem__Channel

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, TWO, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use ISP_Flow_Problem

  implicit none
  private

  public :: FlowProblem_Channel

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the problem

  type, extends(FlowProblem) :: FlowProblem_Channel

    real(RNP) :: re_t   !< friction Reynolds number, Re_τ = u_τ δ / ν
    real(RNP) :: delta  !< channel half height δ
    real(RNP) :: alpha  !< max relative perturbation of velocity components

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetBoundaryTimeDerivative
    procedure :: GetExternalSources

  end type FlowProblem_Channel

contains

!===============================================================================
! type bound procedures

!-------------------------------------------------------------------------------
!> Initialization
!>
!> The default initialization adopts the data from Kim et al. (1987), i.e.
!> Re_τ = u_τ δ / ν = 180, δ = 1, ν = 1/Re_τ and, hence,  u_τ = 1 for the
!> friction velocity.

subroutine SetProblem(problem, file, comm)
  class(FlowProblem_Channel), intent(inout) :: problem
  character(len=*), optional, intent(in) :: file  !< input file
  type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

  ! local variables ............................................................

  ! problem parameters according to Kim et al. (1987)
  real(RNP) :: re_t  = 180  ! friction Reynolds number, Re_τ = u_τ δ / ν
  real(RNP) :: alpha = 1    ! max relative perturbation of velocity components
  real(RNP) :: l     = 4*PI ! channel length
  real(RNP) :: h     = 2    ! channel height, h = 2δ
  real(RNP) :: w     = 2*PI ! channel width
  real(RNP) :: r_svv = 0    ! relative SVV amplitude, r_svv = ν_svv / ν
  character, allocatable :: bc(:,:)

  namelist /parameters/ re_t, alpha, l, h, w, r_svv

  logical :: exists
  integer :: prm, rank

  ! preliminaries ..............................................................

  if (present(comm)) then
    call MPI_Comm_rank(comm, rank)
  else
    rank = 0
  end if

  ! check for input file
  if (rank == 0 .and. present(file)) then

    exists = len_trim(file) > 0
    if (exists) then
      inquire(file=trim(file)//'.prm', exist=exists)
    end if

    if (exists) then
      open(newunit=prm, file=trim(file)//'.prm')
    else
      call Warning('SetProblem', 'Input file "'//trim(file)//'.prm" not found')
    end if

  else
    exists = .false.
  end if

  ! parameters .................................................................

  ! default BC (pressure is ignored)
  allocate(bc(6, problem % nc), source = 'P')
  bc(3:4,:) = 'D' ! bottom and top wall

  if (rank == 0 .and. exists) then
    read(prm, nml=parameters)
    close(prm)
  end if

  if (present(comm)) then
    call XMPI_Bcast(re_t , 0, comm)
    call XMPI_Bcast(alpha, 0, comm)
    call XMPI_Bcast(l    , 0, comm)
    call XMPI_Bcast(h    , 0, comm)
    call XMPI_Bcast(w    , 0, comm)
    call XMPI_Bcast(r_svv, 0, comm)
  end if

  allocate(problem % nu_ref(     problem%nc ), source = ZERO)
  allocate(problem % nu_svv_ref( problem%nc ), source = ZERO)

  problem % stokes          = .false.   ! Navier-Stokes flow
  problem % exact_solution  = .false.   ! exact solution is not provided
  problem % re_t            =  re_t
  problem % delta           =  h / 2
  problem % alpha           =  alpha
  problem % v_ref           =  13.26 / h * re_t**(ONE/7) ! bulk velo (Dean 1978)
  problem % nu_ref(1:3)     =  ONE   / re_t
  problem % nu_svv_ref(1:3) =  r_svv / re_t
  problem % x0              =  ZERO
  problem % x1              = [l, h, w]

  call move_alloc(bc, problem % bc)
  call problem % GeneratePressureBC()

end subroutine SetProblem

!-------------------------------------------------------------------------------
!> Provides the initial values u(x,0) for mesh points x

subroutine GetInitialValues(problem, x, u)
  class(FlowProblem_Channel), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

  integer :: m, n

  n = size(x(:,:,:,:,1))

  call GetInitialVelocity(problem, n, x(:,:,:,:,2), u(:,:,:,:,1:3))

  ! remaining variables get zero
  do m = 4, size(u,5)
    call SetArray(u(:,:,:,:,m), ZERO)
  end do

end subroutine GetInitialValues

!-------------------------------------------------------------------------------
!> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

subroutine GetBoundaryValues(problem, b, xb, t, ub)
  class(FlowProblem_Channel), intent(in) :: problem
  integer,   intent(in)  :: b             !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
  real(RNP), intent(in)  :: t             !< time
  real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

  call SetArray(ub, ZERO, multi=.true.)

  ! silence the compiler ;)
  if (b < 0 .or. size(xb) < 0 .or. t < 0) return

end subroutine GetBoundaryValues

!-------------------------------------------------------------------------------
!> Provides the values `dt_ub = ∂u/∂t(xb,t)` for points `xb` on boundary `b`

subroutine GetBoundaryTimeDerivative(problem, b, xb, t, dt_ub)
  class(FlowProblem_Channel), intent(in)  :: problem
  integer,   intent(in)  :: b              !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)    !< boundary points
  real(RNP), intent(in)  :: t              !< time
  real(RNP), intent(out) :: dt_ub(:,:,:,:) !< ∂u/∂t

  call SetArray(dt_ub, ZERO, multi=.true.)

  ! silence the compiler ;)
  if (b < 0 .or. size(xb) < 0 .or. t < 0) return

end subroutine GetBoundaryTimeDerivative

!-------------------------------------------------------------------------------
!> Provides the external sources for all variables at points x and time t

subroutine GetExternalSources(problem, x, t, f)
  class(FlowProblem_Channel), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(in)  :: t            !< time
  real(RNP), intent(out) :: f(:,:,:,:,:) !< external sources

  ! f₁ = (u_τ)² / (2δ)
  call SetArray(f(:,:,:,:,1 ), HALF / problem % delta ** 3)
  call SetArray(f(:,:,:,:,2:), ZERO, multi=.true.)

  ! silence the compiler ;)
  if (size(x) < 0 .or. t < 0) return

end subroutine GetExternalSources

!===============================================================================
! problem-specific procedures

!-------------------------------------------------------------------------------
!> Velocity

subroutine GetInitialVelocity(problem, n, y, v)
  class(FlowProblem_Channel), intent(in) :: problem
  integer,   intent(in)  :: n      !< number of points
  real(RNP), intent(in)  :: y(n)   !< y-coordinate of points
  real(RNP), intent(out) :: v(n,3) !< velocity at mesh points

  real(RNP) :: a, yh
  integer   :: i

  associate( u_m  => problem % v_ref )

    call SetArray(v, ZERO, multi = .true.)

    ! initialize v with random numbers in the range [0,1)
    call random_number(v)

    ! rescale random values to ±α u_m and add linear profile to v₁
    !$omp do
    do i = 1, n

      yh = y(i) / (2 * problem % delta)

      ! eliminate fluctuations on the wall
      if (min(yh, 1-yh) < epsilon(yh)) then
        a = 0
      else
        a = problem % alpha * u_m
      end if

      v(i,1) = a * (2*v(i,1) - 1)  +  6 * u_m * (1 - yh) * yh
      v(i,2) = a * (2*v(i,2) - 1)
      v(i,3) = a * (2*v(i,3) - 1)

    end do

  end associate

end subroutine GetInitialVelocity

!===============================================================================

end module ISP_Flow_Problem__Channel

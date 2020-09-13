!> summary:  Taylor-Green vortex adopted from Minion & Saye (2018)
!> author:   Joerg Stiller
!> date:     2018/06/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!> The solution is dimensionally inconsistent.
!> @endnote
!===============================================================================

module ISP_Flow_Problem__Vortex_TG

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use ISP_Flow_Problem

  implicit none
  private

  public :: FlowProblem_Vortex_TG

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the Minion-Saye vortex problem

  type, extends(FlowProblem) :: FlowProblem_Vortex_TG

    real(RNP) :: vt(3) = 0  !< translation velocity
    real(RNP) :: xt(3) = 0  !< initial displacement

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetBoundaryTimeDerivative
    procedure :: GetExternalSources
    procedure :: GetExactSolution
    procedure :: GetExactTimeDerivative

  end type FlowProblem_Vortex_TG

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, file, comm)
    class(FlowProblem_Vortex_TG), intent(inout) :: problem
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    logical   :: stokes = .false.
    real(RNP) :: nu     =  0.01         ! kinematic viscosity
    real(RNP) :: vt(3)  =  [1,1,0]      ! translation velocity
    real(RNP) :: xt(3)  =  [0,1,0] / 8. ! initial displacement
    real(RNP) :: x0(3)  = -HALF         ! bounding box: corner nearest to -∞
    real(RNP) :: x1(3)  =  HALF         ! bounding box: corner nearest to +∞
    character, allocatable :: bc(:,:)

    namelist /parameters/ stokes, nu, vt, xt, x0, x1, bc

    logical :: exists
    integer :: prm, rank
    integer :: b1, b2, d

    ! preliminaries ............................................................

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
        call Warning('SetProblem','Input file "'//trim(file)//'.prm" not found')
      end if

    else
      exists = .false.
    end if

    ! parameters ...............................................................

    ! default BC (pressure is ignored)
    allocate(bc(6, problem % nc), source = 'P')

    if (rank == 0 .and. exists) then
      read(prm, nml=parameters)
      close(prm)
    end if

    if (present(comm)) then
      call XMPI_Bcast(stokes, 0, comm)
      call XMPI_Bcast(nu    , 0, comm)
      call XMPI_Bcast(vt    , 0, comm)
      call XMPI_Bcast(xt    , 0, comm)
      call XMPI_Bcast(x0    , 0, comm)
      call XMPI_Bcast(x1    , 0, comm)
      call XMPI_Bcast(bc    , 0, comm)
    end if

    problem % stokes = stokes
    problem % exact_solution = .true.

    allocate(problem % nu_ref( problem%nc ), source = ZERO)
    problem % nu_ref(1:3) = nu

    problem % vt = vt
    problem % xt = xt
    problem % x0 = x0
    problem % x1 = x1

    ! pressure BC
    do d = 1, 3
      b1 = 2*d - 1
      b2 = b1  + 1
      if (all(bc(b1:b2, d) == 'P')) then
        bc(b1:b2, 4) = 'P'
      else
        bc(b1:b2, 4) = 'N'
      end if
    end do

    call move_alloc(bc, problem % bc)

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(FlowProblem_Vortex_TG), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref(1), vt => problem%vt, xt => problem%xt)
      call GetVelocity(n, x, ZERO, nu, vt, xt, u(:,:,:,:,1:3))
      call GetPressure(n, x, ZERO, nu, vt, xt, u(:,:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(FlowProblem_Vortex_TG), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary ID
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    integer :: m, n

    n = size(xb(:,:,:,1))

    associate(nu => problem%nu_ref(1), vt => problem%vt, xt => problem%xt)
      call GetVelocity(n, xb, t, nu, vt, xt, ub(:,:,:,1:3))
      call GetPressure(n, xb, t, nu, vt, xt, ub(:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(ub,4)
      call SetArray(ub(:,:,:,m), ZERO)
    end do

    ! silence the compiler ;)
    if (b < 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the values `dt_ub = ∂u/∂t(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryTimeDerivative(problem, b, xb, t, dt_ub)
    class(FlowProblem_Vortex_TG), intent(in)  :: problem
    integer,   intent(in)  :: b              !< boundary ID
    real(RNP), intent(in)  :: xb(:,:,:,:)    !< boundary points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: dt_ub(:,:,:,:) !< ∂u/∂t

    integer :: m, n

    n = size(xb(:,:,:,1))

    associate(nu => problem%nu_ref(1), vt => problem%vt, xt => problem%xt)
      call GetVelocityTimeDerivative(n, xb, t, nu, vt, xt, dt_ub(:,:,:,1:3))
      call GetPressureTimeDerivative(n, xb, t, nu, vt, xt, dt_ub(:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(dt_ub,4)
      call SetArray(dt_ub(:,:,:,m), ZERO)
    end do

    ! silence the compiler ;)
    if (b < 0) return

  end subroutine GetBoundaryTimeDerivative

  !-----------------------------------------------------------------------------
  !> Provides the external sources for all variables at points x and time t

  subroutine GetExternalSources(problem, x, t, f)
    class(FlowProblem_Vortex_TG), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(out) :: f(:,:,:,:,:) !< external sources

    integer :: m, n

    n = size(x(:,:,:,:,1))

    if (problem % stokes) then
    associate(nu => problem%nu_ref(1), vt => problem%vt, xt => problem%xt)
        call GetStokesSource(n, x, t, nu, vt, xt, f(:,:,:,:,1:3))
        do m = 4, size(f,5)
          call SetArray(f(:,:,:,:,m), ZERO)
        end do
      end associate
    else
      do m = 1, size(f,5)
        call SetArray(f(:,:,:,:,m), ZERO)
      end do
    end if

  end subroutine GetExternalSources

  !---------------------------------------------------------------------------
  !> Provides the exact solution u(x,t), or zero if not available

  subroutine GetExactSolution(problem, x, t, u)
    class(FlowProblem_Vortex_TG), intent(in) :: problem
    real(RNP),          intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP),          intent(in)  :: t            !< time
    real(RNP),          intent(out) :: u(:,:,:,:,:) !< solution

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref(1), vt => problem%vt, xt => problem%xt)
      call GetVelocity(n, x, t, nu, vt, xt, u(:,:,:,:,1:3))
      call GetPressure(n, x, t, nu, vt, xt, u(:,:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactSolution

  !---------------------------------------------------------------------------
  !> Provides the time derivative of the exact solution, ∂u/∂t(x,t)

  subroutine GetExactTimeDerivative(problem, x, t, dt_u)
    class(FlowProblem_Vortex_TG), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)    !< mesh points
    real(RNP), intent(in)  :: t               !< time
    real(RNP), intent(out) :: dt_u(:,:,:,:,:) !< time derivative

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref(1), vt => problem%vt, xt => problem%xt)
      call GetVelocityTimeDerivative(n, x, t, nu, vt, xt, dt_u(:,:,:,:,1:3))
      call GetPressureTimeDerivative(n, x, t, nu, vt, xt, dt_u(:,:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(dt_u,5)
      call SetArray(dt_u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactTimeDerivative

  !=============================================================================
  ! problem-specific procedures

  !-----------------------------------------------------------------------------
  !> Velocity

  subroutine GetVelocity(n, x, t, nu, vt, xt, v)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(in)  :: nu      !< kinematic viscosity
    real(RNP), intent(in)  :: vt(3)   !< translation velocity
    real(RNP), intent(in)  :: xt(3)   !< initial displacement
    real(RNP), intent(out) :: v(n,3)  !< velocity at mesh points

    real(RNP) :: a, c, phi1, phi2
    integer   :: i

    a = -8 * PI**2 * nu
    c =  exp(a * t)

    !$omp do
    do i = 1, n
      phi1 = 2 * PI * (x(i,1) - vt(1)*t - xt(1))
      phi2 = 2 * PI * (x(i,2) - vt(2)*t - xt(2))
      v(i,1) = vt(1) + c * sin(phi1) * cos(phi2)
      v(i,2) = vt(2) - c * cos(phi1) * sin(phi2)
      v(i,3) = vt(3)
    end do

  end subroutine GetVelocity

  !-----------------------------------------------------------------------------
  !> Velocity time derivative

  subroutine GetVelocityTimeDerivative(n, x, t, nu, vt, xt, dt_v)
    integer,   intent(in)  :: n          !< number of points
    real(RNP), intent(in)  :: x(n,3)     !< mesh points
    real(RNP), intent(in)  :: t          !< time
    real(RNP), intent(in)  :: nu         !< kinematic viscosity
    real(RNP), intent(in)  :: vt(3)      !< translation velocity
    real(RNP), intent(in)  :: xt(3)      !< initial displacement
    real(RNP), intent(out) :: dt_v(n,3)  !< velocity at mesh points

    real(RNP) :: a, c
    real(RNP) :: phi1, sin1, cos1
    real(RNP) :: phi2, sin2, cos2
    integer   :: i

    a = -8 * PI**2 * nu
    c =  exp(a * t)

    !$omp do
    do i = 1, n

      phi1 = 2 * PI * (x(i,1) - vt(1)*t - xt(1))
      sin1 = sin(phi1)
      cos1 = cos(phi1)

      phi2 = 2 * PI * (x(i,2) - vt(2)*t - xt(2))
      sin2 = sin(phi2)
      cos2 = cos(phi2)

      dt_v(i,1) =  c * (a*sin1*cos2 - 2*PI*(vt(1)*cos1*cos2 - vt(2)*sin1*sin2))
      dt_v(i,2) = -c * (a*cos1*sin2 + 2*PI*(vt(1)*sin1*sin2 - vt(2)*cos1*cos2))
      dt_v(i,3) =  ZERO

    end do

  end subroutine GetVelocityTimeDerivative

  !-----------------------------------------------------------------------------
  !> Pressure

  subroutine GetPressure(n, x, t, nu, vt, xt, p)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(in)  :: nu      !< kinematic viscosity
    real(RNP), intent(in)  :: vt(3)   !< translation velocity
    real(RNP), intent(in)  :: xt(3)   !< initial displacement
    real(RNP), intent(out) :: p(n)    !< pressure at mesh points

    real(RNP) :: a, c, phi1, phi2
    integer   :: i

    a = -16 * PI**2 * nu
    c =  exp(a * t) / 4

    !$omp do
    do i = 1, n
      phi1 = 4 * PI * (x(i,1) - vt(1)*t - xt(1))
      phi2 = 4 * PI * (x(i,2) - vt(2)*t - xt(2))
      p(i) = c * (cos(phi1) + cos(phi2))
    end do

  end subroutine GetPressure

  !-----------------------------------------------------------------------------
  !> Pressure time derivative

  subroutine GetPressureTimeDerivative(n, x, t, nu, vt, xt, dt_p)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(in)  :: nu      !< kinematic viscosity
    real(RNP), intent(in)  :: vt(3)   !< translation velocity
    real(RNP), intent(in)  :: xt(3)   !< initial displacement
    real(RNP), intent(out) :: dt_p(n) !< pressure time derivative

    real(RNP) :: a, c, phi1, phi2
    integer   :: i

    a = -16 * PI**2 * nu
    c = -exp(a * t)

    !$omp do
    do i = 1, n
      phi1 = 4 * PI * (x(i,1) - vt(1)*t - xt(1))
      phi2 = 4 * PI * (x(i,2) - vt(2)*t - xt(2))
      dt_p(i) = c * ( a * cos(phi1) + 4*PI*vt(1) * sin(phi1) &
                    + a * cos(phi2) + 4*PI*vt(2) * sin(phi2) )
    end do

  end subroutine GetPressureTimeDerivative

  !-----------------------------------------------------------------------------
  !> Stokes source

  subroutine GetStokesSource(n, x, t, nu, vt, xt, f)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(in)  :: nu      !< kinematic viscosity
    real(RNP), intent(in)  :: vt(3)   !< translation velocity
    real(RNP), intent(in)  :: xt(3)   !< initial displacement
    real(RNP), intent(out) :: f(n,3)  !< nonlinear term, -v · ∇v

    real(RNP) :: a, c
    real(RNP) :: phi1, sin1, cos1, v1
    real(RNP) :: phi2, sin2, cos2, v2
    integer   :: i

    a = -8 * PI**2 * nu
    c =  exp(a * t)

    !$omp do
    do i = 1, n

      phi1 = 2 * PI * (x(i,1) - vt(1)*t - xt(1))
      sin1 = sin(phi1)
      cos1 = cos(phi1)

      phi2 = 2 * PI * (x(i,2) - vt(2)*t - xt(2))
      sin2 = sin(phi2)
      cos2 = cos(phi2)

      v1 = vt(1) + c * sin1 * cos2
      v2 = vt(2) - c * cos1 * sin2

      f(i,1) = 2*PI*c * ( -v1 * cos1*cos2 + v2 * sin1*sin2 )
      f(i,2) = 2*PI*c * ( -v1 * sin1*sin2 + v2 * cos1*cos2 )
      f(i,3) = ZERO

    end do

  end subroutine GetStokesSource

  !=============================================================================

end module ISP_Flow_Problem__Vortex_TG

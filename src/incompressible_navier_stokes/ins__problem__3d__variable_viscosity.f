!> summary:  Method of Manufactured Solution for vortex field with different
!>           cases of variable diffusivity
!> author:   Karl Schoppmann, Jörg Stiller
!> date:     2019/05/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> Provides procedures for setting
!>
!>   *  the problem
!>
!> and provides procedures for getting
!>
!>   *  initial values
!>   *  boundary values
!>   *  the time derivative of boundary values
!>   *  external sources
!>   *  an exact solution
!>   *  the time derivative of the exact solution
!>   *  a variable viscosity
!>
!> Four cases of variable viscosity are considered:
!>
!>   1. Position dependent viscosity `nu(x)`,
!>      implementation follows
!>      Niemann, M.: Bouyancy effects in turbulent liquid metal flow, 2017
!>   2. Position and time dependent viscosity `nu(x,t)`,
!>      implementation based on case 1 with further development
!>   3. Nonlinear viscosity based on kinetic energy.
!>
!> For more details [see](http://arxiv.org/abs/2001.11902) and
!> Karl Schoppmann, Verifikation von Methoden hoher Ordnung für die numerische
!> Simulation inkompressibler Strömungen mit variabler Viskosität,
!> Projektarbeit zum Forschungspraktikum, ISM, TU Dresden, 2019.
!===============================================================================

module INS__Problem__3D__Variable_Viscosity

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_3D_VariableViscosity

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the vortex field problem

  type, extends(INS_Problem_3D) :: INS_Problem_3D_VariableViscosity

    integer   :: test_case = 1    !< variable viscosity case 1|2|3|4
    real(RNP) :: nu_0      = ONE  !< constant base viscosity
    real(RNP) :: nu_1      = ONE  !< coefficient of variable viscosity

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetExternalSources
    procedure :: GetViscosity

    procedure :: GetExactSolution
    procedure :: GetExactTimeDerivative
    procedure :: GetExactConvectiveTerm
    procedure :: GetExactPressureTerm
    procedure :: GetExactDiffusiveTerm

  end type INS_Problem_3D_VariableViscosity

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization
  !>
  !> The viscosity is split into a constant base part and a variable part

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_3D_VariableViscosity), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    logical   :: stokes    = .false.
    integer   :: test_case =  1         ! 1|2|3 (`nu(x)`|`nu(x,t)`|`nu(x,t,u)`)
    real(RNP) :: nu_0      =  ONE       ! constant base viscosity
    real(RNP) :: nu_1      =  ONE       ! coefficient of variable viscosity
    real(RNP) :: x0(3)     = -HALF      ! bounding box: corner nearest to -∞
    real(RNP) :: x1(3)     =  HALF      ! bounding box: corner nearest to +∞
    character, allocatable :: bc_v(:)

    namelist /parameters/ stokes, test_case, nu_0, nu_1, x0, x1, bc_v

    logical :: exists
    integer :: prm, rank

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
    allocate(bc_v(nb), source = 'D')

    if (rank == 0) then

      ! read parameters
      if (exists) then
        read(prm, nml=parameters)
        close(prm)
      end if

      ! normalize turbulent viscosity to nu_1
      select case(test_case)
      case(3)
        nu_1 = nu_1 / (48 * PI**2)
      case(4)
        nu_1 = nu_1 / 4
      end select

    end if

    if (present(comm)) then
      call XMPI_Bcast(stokes   , 0, comm)
      call XMPI_Bcast(test_case, 0, comm)
      call XMPI_Bcast(nu_0     , 0, comm)
      call XMPI_Bcast(nu_1     , 0, comm)
      call XMPI_Bcast(x0       , 0, comm)
      call XMPI_Bcast(x1       , 0, comm)
      call XMPI_Bcast(bc_v     , 0, comm)
    end if

    problem % stokes         = stokes
    problem % exact_solution = .true.
    problem % variable_props = .true.
    problem % test_case      = test_case
    problem % nu_ref         = nu_0 + HALF * nu_1
    problem % nu_0           = nu_0
    problem % nu_1           = nu_1
    problem % x0             = x0
    problem % x1             = x1

    call move_alloc(bc_v, problem % bc_v)

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values `u(x,0)` for mesh points `x`

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_3D_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetVelocity(n, x, ZERO, u(:,:,:,:,1:3))
    call GetPressure(n, x, ZERO, u(:,:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` in points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_3D_VariableViscosity), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary identifier
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    integer :: m, n

    n = size(xb(:,:,:,1))

    call GetVelocity(n, xb, t, ub(:,:,:,1:3))
    call GetPressure(n, xb, t, ub(:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(ub,4)
      call SetArray(ub(:,:,:,m), ZERO)
    end do

    ! ignore boundary identifier
    if (b < 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the external sources for all variables at points `x` and time `t`

  subroutine GetExternalSources(problem, x, t, F_s)
    class(INS_Problem_3D_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_s(:,:,:,:,:) !< external sources

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call Get_MMS_Source(problem, n, x, t, F_s(:,:,:,:,1:3))

    do m = 4, size(F_s,5)
      call SetArray(F_s(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExternalSources

  !-----------------------------------------------------------------------------
  !> Provides a variable viscosity depending on input parameter `test_case`
  !>
  !> 1. `nu(x)`
  !> 2. `nu(x,t)`
  !> 3. `nu(x,t,u)`, where `u` is the approximate solution

  subroutine GetViscosity(problem, x, t, u, nu)
    class(INS_Problem_3D_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x (:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(in)  :: u (:,:,:,:,:) !< flow variables
    real(RNP), intent(out) :: nu(:,:,:,:)   !< viscosity

    integer :: m, n

    n = size(nu)

    associate(nu_0 => problem % nu_0, nu_1 => problem % nu_1)
      select case(problem % test_case)
        case(1)
          call GetViscosity_Case_1(n, x, nu_0, nu_1, nu)     ! nu(x)
        case(2)
          call GetViscosity_Case_2(n, x, t, nu_0, nu_1, nu)  ! nu(x,t)
        case(3)
          call GetViscosity_Case_3(n, nu_0, nu_1, u, nu)     ! nu(x,t,u)
      end select
    end associate

  end subroutine GetViscosity

  !-----------------------------------------------------------------------------
  !> Provides the exact solution `u(x,t)`

  subroutine GetExactSolution(problem, x, t, u)
    class(INS_Problem_3D_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: u(:,:,:,:,:)   !< solution

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetVelocity(n, x, t, u(:,:,:,:,1:3))
    call GetPressure(n, x, t, u(:,:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactSolution

  !-----------------------------------------------------------------------------
  !> Provides the time derivative of the exact solution, `∂u/∂t(x,t)`

  subroutine GetExactTimeDerivative(problem, x, t, dt_u)
    class(INS_Problem_3D_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)    !< mesh points
    real(RNP), intent(in)  :: t               !< time
    real(RNP), intent(out) :: dt_u(:,:,:,:,:) !< time derivative

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetVelocityTimeDerivative(n, x, t, dt_u(:,:,:,:,1:3))
    call GetPressureTimeDerivative(n, x, t, dt_u(:,:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(dt_u,5)
      call SetArray(dt_u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactTimeDerivative

  !---------------------------------------------------------------------------
  !> Exact convection term, F_c = -[∇⋅vv,0]

  subroutine GetExactConvectiveTerm(problem, x, t, F_c)
    class(INS_Problem_3D_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_c(:,:,:,:,:) !< convection term

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetConvectiveTerm(n, x, t, F_c(:,:,:,:,1:3))
    do m = 4, size(F_c,5)
      call SetArray(F_c(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactConvectiveTerm

  !---------------------------------------------------------------------------
  !> Exact pressure term, F_p = [-∇p,0]

  subroutine GetExactPressureTerm(problem, x, t, F_p)
    class(INS_Problem_3D_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_p(:,:,:,:,:) !< diffusion term

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetPressureTerm(n, x, t, F_p(:,:,:,:,1:3))
    do m = 4, size(F_p,5)
      call SetArray(F_p(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactPressureTerm

  !---------------------------------------------------------------------------
  !> Exact diffusion term, F_d = [∇⋅τ,0]

  subroutine GetExactDiffusiveTerm(problem, x, t, F_d)
    class(INS_Problem_3D_VariableViscosity), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_d(:,:,:,:,:) !< diffusion term

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetDiffusiveTerm(problem, n, x, t, F_d(:,:,:,:,1:3))
    do m = 4, size(F_d,5)
      call SetArray(F_d(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactDiffusiveTerm

  !=============================================================================
  ! problem-specific procedures

  !-----------------------------------------------------------------------------
  !> Velocity

  subroutine GetVelocity(n, x, t, v)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: v(n,3)  !< velocity at mesh points

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)
      v(i,1) = ( sin(kxt) + cos(kyt) ) * sin(kzt)
      v(i,2) = ( cos(kxt) + sin(kyt) ) * sin(kzt)
      v(i,3) = ( cos(kxt) + cos(kyt) ) * cos(kzt)
    end do

  end subroutine GetVelocity

  !-----------------------------------------------------------------------------
  !> Velocity time derivative

  subroutine GetVelocityTimeDerivative(n, x, t, dt_v)
    integer,   intent(in)  :: n          !< number of points
    real(RNP), intent(in)  :: x(n,3)     !< mesh points
    real(RNP), intent(in)  :: t          !< time
    real(RNP), intent(out) :: dt_v(n,3)  !< time derivative of velocity

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)

      dt_v(i,1) =  k * ( ( cos(kxt) - sin(kyt)) * sin(kzt) &
                       + ( sin(kxt) + cos(kyt)) * cos(kzt) )
      dt_v(i,2) =  k * ( (-sin(kxt) + cos(kyt)) * sin(kzt) &
                       + ( cos(kxt) + sin(kyt)) * cos(kzt) )
      dt_v(i,3) = -k * ( ( sin(kxt) + sin(kyt)) * cos(kzt) &
                       + ( cos(kxt) + cos(kyt)) * sin(kzt) )
    end do

  end subroutine GetVelocityTimeDerivative

  !-----------------------------------------------------------------------------
  !> Pressure

  subroutine GetPressure(n, x, t, p)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: p(n)    !< pressure at mesh points

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)
      p(i) = sin(kxt) * sin(kyt) * sin(kzt)
    end do

  end subroutine GetPressure

  !-----------------------------------------------------------------------------
  !> Pressure time derivative

  subroutine GetPressureTimeDerivative(n, x, t, dt_p)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: dt_p(n) !< time derivative of pressure

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)
      dt_p(i) = k * ( cos(kxt) * sin(kyt) * sin(kzt) &
                    + sin(kxt) * cos(kyt) * sin(kzt) &
                    + sin(kxt) * sin(kyt) * cos(kzt) )
    end do

  end subroutine GetPressureTimeDerivative

  !-----------------------------------------------------------------------------
  !> Viscosity, case 1: position dependent viscosity `nu(x)`

  subroutine GetViscosity_Case_1(n, x, nu_0, nu_1, nu)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: nu_0    !< constant base viscosity
    real(RNP), intent(in)  :: nu_1    !< coefficient of variable viscosity
    real(RNP), intent(out) :: nu(n)   !< kinematic viscosity

    real(RNP) :: k, kx, ky, kz
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kx = k * x(i,1)
      ky = k * x(i,2)
      kz = k * x(i,3)

      nu(i) = nu_0 + nu_1 * sin(kx)**2 * sin(ky)**2 * sin(kz)**2
    end do

  end subroutine GetViscosity_Case_1

  !-----------------------------------------------------------------------------
  !> Viscosity, case 2: position and time dependent viscosity `nu(x,t)`

  subroutine GetViscosity_Case_2(n, x, t, nu_0, nu_1, nu)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(in)  :: nu_0    !< constant base viscosity
    real(RNP), intent(in)  :: nu_1    !< coefficient of variable viscosity
    real(RNP), intent(out) :: nu(n)   !< kinematic viscosity

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) - t)
      kyt = k * (x(i,2) - t)
      kzt = k * (x(i,3) - t)

      nu(i) = nu_0 + nu_1 * sin(kxt)**2 * sin(kyt)**2 * sin(kzt)**2
    end do

  end subroutine GetViscosity_Case_2

  !-----------------------------------------------------------------------------
  !> Viscosity, case 3: solution dependent viscosity

  subroutine GetViscosity_Case_3(n, nu_0, nu_1, v, nu)
    integer,   intent(in)  :: n        !< number of points
    real(RNP), intent(in)  :: nu_0     !< constant base viscosity
    real(RNP), intent(in)  :: nu_1     !< coefficient of variable viscosity
    real(RNP), intent(in)  :: v(n,3)   !< velocity
    real(RNP), intent(out) :: nu(n)    !< kinematic viscosity

    integer :: i

    !$omp do
    do i = 1, n
      nu(i) = nu_0 + nu_1 * (v(i,1)**2 + v(i,2)**2 + v(i,3)**2)
    end do

  end subroutine GetViscosity_Case_3

  !-----------------------------------------------------------------------------
  !> Viscosity gradient, case 1

  subroutine GetViscosityGradient_Case_1(n, x, nu_1, grad_nu)
    integer,   intent(in)  :: n            !< number of points
    real(RNP), intent(in)  :: x(n,3)       !< mesh points
    real(RNP), intent(in)  :: nu_1         !< coefficient of variable viscosity
    real(RNP), intent(out) :: grad_nu(n,3) !< viscosity gradient

    real(RNP) :: k, kx, ky, kz
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kx = k * x(i,1)
      ky = k * x(i,2)
      kz = k * x(i,3)

      grad_nu(i,1) = 2 * k * nu_1 * sin(kx) * cos(kx) * sin(ky)**2 * sin(kz)**2
      grad_nu(i,2) = 2 * k * nu_1 * sin(ky) * cos(ky) * sin(kx)**2 * sin(kz)**2
      grad_nu(i,3) = 2 * k * nu_1 * sin(kz) * cos(kz) * sin(kx)**2 * sin(ky)**2
    end do

  end subroutine GetViscosityGradient_Case_1

  !-----------------------------------------------------------------------------
  !> Viscosity gradient, case 2

  subroutine GetViscosityGradient_Case_2(n, x, t, nu_1, grad_nu)
    integer,   intent(in)  :: n            !< number of points
    real(RNP), intent(in)  :: x(n,3)       !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(in)  :: nu_1         !< coefficient of variable viscosity
    real(RNP), intent(out) :: grad_nu(n,3) !< viscosity gradient

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) - t)
      kyt = k * (x(i,2) - t)
      kzt = k * (x(i,3) - t)

      grad_nu(i,1) = 2 * k * nu_1 * sin(kxt) * cos(kxt) &
                     * sin(kyt)**2 * sin(kzt)**2
      grad_nu(i,2) = 2 * k * nu_1 * sin(kyt) * cos(kyt) &
                     * sin(kxt)**2 * sin(kzt)**2
      grad_nu(i,3) = 2 * k * nu_1 * sin(kzt) * cos(kzt) &
                     * sin(kxt)**2 * sin(kyt)**2
    end do

  end subroutine GetViscosityGradient_Case_2

  !------------------------------------------------------------------------------
  !> Viscosity gradient, case 3

  subroutine GetViscosityGradient_Case_3(n, x, t, nu_1, grad_nu)
    integer,   intent(in)  :: n            !< number of points
    real(RNP), intent(in)  :: x(n,3)       !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(in)  :: nu_1         !< coefficient of variable viscosity
    real(RNP), intent(out) :: grad_nu(n,3) !< viscosity gradient

    real(RNP), allocatable, save :: v(:,:), grad_v(:,:,:)
    integer :: i, l

    !$omp master
    allocate(v(n,3), grad_v(n,3,3))
    !$omp end master
    !$omp barrier

    call GetVelocity(n, x, t, v)
    call GetVelocityGradient(n, x, t, grad_v)

    do l = 1, 3
      !$omp do
      do i = 1, n
        grad_nu(i,l) = nu_1 * 2 * ( grad_v(i,l,1) * v(i,1) &
                                  + grad_v(i,l,2) * v(i,2) &
                                  + grad_v(i,l,3) * v(i,3) )
      end do
    end do

    !$omp master
    deallocate(v, grad_v)
    !$omp end master

  end subroutine GetViscosityGradient_Case_3

  !-----------------------------------------------------------------------------
  !> MMS source

  subroutine Get_MMS_Source(problem, n, x, t, F_s)
    class(INS_Problem_3D_VariableViscosity), intent(in)  :: problem
    integer,   intent(in)  :: n         !< number of points
    real(RNP), intent(in)  :: x(n,3)    !< mesh points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(out) :: F_s(n,3)  !< MMS source

    real(RNP), allocatable, save :: dt_v(:,:), F_c(:,:), F_p(:,:), F_d(:,:)
    integer :: i, l

    !$omp master
    allocate(dt_v(n,3), F_c(n,3), F_d(n,3), F_p(n,3))
    !$omp end master
    !$omp barrier

    call GetVelocityTimeDerivative(n, x, t, dt_v)
    call GetDiffusiveTerm(problem, n, x, t, F_d)

    if(problem % stokes) then

      do l = 1, 3
        !$omp do
        do i = 1, n
          F_s(i,l) = dt_v(i,l) - F_d(i,l)
        end do
      end do

    else ! Navier-Stokes

      call GetConvectiveTerm (n, x, t, F_c)
      call GetPressureTerm   (n, x, t, F_p)
      do l = 1, 3
        !$omp do
        do i = 1, n
          F_s(i,l) = dt_v(i,l) - (F_c(i,l) + F_d(i,l) + F_p(i,l))
        end do
      end do

    end if

    !$omp master
    deallocate(dt_v, F_c, F_p, F_d)
    !$omp end master

  end subroutine Get_MMS_Source

  !-----------------------------------------------------------------------------
  !> Pressure term, -∇p

  subroutine GetPressureTerm(n, x, t, F_p)
    integer,   intent(in)  :: n          !< number of points
    real(RNP), intent(in)  :: x(n,3)     !< mesh points
    real(RNP), intent(in)  :: t          !< time
    real(RNP), intent(out) :: F_p(n,3)   !< pressure gradient

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)
      F_p(i,1) = -k * cos(kxt) * sin(kyt) * sin(kzt)
      F_p(i,2) = -k * sin(kxt) * cos(kyt) * sin(kzt)
      F_p(i,3) = -k * sin(kxt) * sin(kyt) * cos(kzt)
    end do

  end subroutine GetPressureTerm

  !-----------------------------------------------------------------------------
  !> Velocity Gradient
  !>
  !>       grad_v(i,l,m) = ∂v_m/∂x_l
  !>
  !> for mesh point `x(i,:)` and time `t`

  subroutine GetVelocityGradient(n, x, t, grad_v)
    integer,   intent(in)  :: n               !< number of points
    real(RNP), intent(in)  :: x(n,3)          !< mesh points
    real(RNP), intent(in)  :: t               !< time
    real(RNP), intent(out) :: grad_v(n,3,3)   !< velocity gradient

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)
      grad_v(i,1,1) =  k * cos(kxt) * sin(kzt)                  ! ∂u/∂x
      grad_v(i,2,1) = -k * sin(kyt) * sin(kzt)                  ! ∂u/∂y
      grad_v(i,3,1) =  k * ( sin(kxt) + cos(kyt) ) * cos(kzt)   ! ∂u/∂z

      grad_v(i,1,2) = -k * sin(kxt) * sin(kzt)                  ! ∂v/∂x
      grad_v(i,2,2) =  k * cos(kyt) * sin(kzt)                  ! ∂v/∂y
      grad_v(i,3,2) =  k * ( cos(kxt) + sin(kyt) ) * cos(kzt)   ! ∂v/∂z

      grad_v(i,1,3) = -k * sin(kxt) * cos(kzt)                  ! ∂w/∂x
      grad_v(i,2,3) = -k * sin(kyt) * cos(kzt)                  ! ∂w/∂y
      grad_v(i,3,3) = -k * ( cos(kxt) + cos(kyt) ) * sin(kzt)   ! ∂w/∂z
    end do

  end subroutine GetVelocityGradient

  !-----------------------------------------------------------------------------
  !> Laplacian of velocity components
  !>
  !>       lapl_v(i,j) = ∇²v_j
  !>
  !> for mesh point `x(i,:)` and time `t`

  subroutine GetVelocityLaplacian(n, x, t, lapl_v)
    integer,   intent(in)  :: n           !< number of points
    real(RNP), intent(in)  :: x(n,3)      !< mesh points
    real(RNP), intent(in)  :: t           !< time
    real(RNP), intent(out) :: lapl_v(n,3) !< Laplacian of velocity

    real(RNP) :: k, kxt, kyt, kzt
    integer   :: i

    k = 2 * PI  ! wave number

    !$omp do
    do i = 1, n
      kxt = k * (x(i,1) + t)
      kyt = k * (x(i,2) + t)
      kzt = k * (x(i,3) + t)

      lapl_v(i,1) = -2 * k**2 * (sin(kxt)+cos(kyt)) * sin(kzt)  ! ∇²u(x,t)
      lapl_v(i,2) = -2 * k**2 * (cos(kxt)+sin(kyt)) * sin(kzt)  ! ∇²v(x,t)
      lapl_v(i,3) = -2 * k**2 * (cos(kxt)+cos(kyt)) * cos(kzt)  ! ∇²w(x,t)
    end do

  end subroutine GetVelocityLaplacian

!   !-----------------------------------------------------------------------------
!   !> Hessian of velocity components -- currently not used
!   !>
!   !>       hess_v(i,l,m,j) = ∂²v_j /(∂x_l ∂x_m)
!   !>
!   !> for mesh point `x(i,:)` and time `t`
!
!   subroutine GetVelocityHessian(n, x, t, hess_v)
!     integer,   intent(in)  :: n               !< number of points
!     real(RNP), intent(in)  :: x(n,3)          !< mesh points
!     real(RNP), intent(in)  :: t               !< time
!     real(RNP), intent(out) :: hess_v(n,3,3,3) !< hessian of velocity
!
!     real(RNP) :: k, kxt, kyt, kzt
!     integer   :: i
!
!     k = 2 * PI  ! wave number
!
!     !$omp do
!     do i = 1, n
!       kxt = k * (x(i,1) + t)
!       kyt = k * (x(i,2) + t)
!       kzt = k * (x(i,3) + t)
!
!       ! hessian for velocity component u
!       hess_v(i,1,1,1) = -k**2 * sin(kxt) * sin(kzt)
!       hess_v(i,2,2,1) = -k**2 * cos(kyt) * sin(kzt)
!       hess_v(i,3,3,1) = -k**2 * (sin(kxt) + cos(kyt)) * sin(kzt)
!       hess_v(i,2,3,1) = -k**2 * sin(kyt) * cos(kzt)
!       hess_v(i,1,3,1) =  k**2 * cos(kxt) * cos(kzt)
!       hess_v(i,1,2,1) = 0
!       hess_v(i,3,2,1) = hess_v(i,2,3,1)
!       hess_v(i,3,1,1) = hess_v(i,1,3,1)
!       hess_v(i,2,1,1) = hess_v(i,1,2,1)
!
!       ! hessian for velocity component v
!       hess_v(i,1,1,2) = -k**2 * cos(kxt) * sin(kzt)
!       hess_v(i,2,2,2) = -k**2 * sin(kyt) * sin(kzt)
!       hess_v(i,3,3,2) = -k**2 * (cos(kxt) + sin(kyt)) * sin(kzt)
!       hess_v(i,2,3,2) =  k**2 * cos(kyt) * cos(kzt)
!       hess_v(i,1,3,2) = -k**2 * sin(kxt) * cos(kzt)
!       hess_v(i,1,2,2) = 0
!       hess_v(i,3,2,2) = hess_v(i,2,3,2)
!       hess_v(i,3,1,2) = hess_v(i,1,3,2)
!       hess_v(i,2,1,2) = hess_v(i,1,2,2)
!
!       ! hessian for velocity component w
!       hess_v(i,1,1,3) = -k**2 * cos(kxt) * cos(kzt)
!       hess_v(i,2,2,3) = -k**2 * cos(kyt) * cos(kzt)
!       hess_v(i,3,3,3) = -k**2 * (cos(kxt) + cos(kyt)) * cos(kzt)
!       hess_v(i,2,3,3) =  k**2 * sin(kyt) * sin(kzt)
!       hess_v(i,1,3,3) =  k**2 * sin(kxt) * sin(kzt)
!       hess_v(i,1,2,3) = 0
!       hess_v(i,3,2,3) = hess_v(i,2,3,3)
!       hess_v(i,3,1,3) = hess_v(i,1,3,3)
!       hess_v(i,2,1,3) = hess_v(i,1,2,3)
!
!     end do
!
!   end subroutine GetVelocityHessian

  !-----------------------------------------------------------------------------
  !> Diffusive term, F_d = ∇⋅τ

  subroutine GetDiffusiveTerm(problem, n, x, t, F_d)
    class(INS_Problem_3D_VariableViscosity), intent(in) :: problem
    integer,   intent(in)  :: n            !< number of points
    real(RNP), intent(in)  :: x(n,3)       !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(out) :: F_d(n,3)     !< diffusive term

    real(RNP), allocatable, save :: v(:,:), lapl_v(:,:), strain_rate(:,:,:)
    real(RNP), allocatable, save :: nu(:,:), grad_nu(:,:)
    real(RNP) :: tmp
    integer   :: i, l, m

    !$omp master
    allocate(v(n,3), lapl_v(n,3), strain_rate(n,3,3))
    allocate(nu(n,3), grad_nu(n,3))
    !$omp end master
    !$omp barrier

    call GetVelocity          (n, x, t, v)
    call GetVelocityLaplacian (n, x, t, lapl_v)
    call GetStrainRate        (n, x, t, strain_rate)

    associate(nu_0 => problem % nu_0,  nu_1 => problem % nu_1)
      select case(problem % test_case)
        case(1)
          call GetViscosity_Case_1(n, x, nu_0, nu_1, nu)
          call GetViscosityGradient_Case_1(n, x, nu_1, grad_nu)
        case(2)
          call GetViscosity_Case_2(n, x, t, nu_0, nu_1, nu)
          call GetViscosityGradient_Case_2(n, x, t, nu_1, grad_nu)
        case(3)
          call GetViscosity_Case_3(n, nu_0, nu_1, v, nu)
          call GetViscosityGradient_Case_3(n, x, t, nu_1, grad_nu)
      end select
    end associate

    ! diffusive term: ∇(ν 2 S) = ∇ν · 2 S + ν ∇²v
    do l = 1, 3
      !$omp do
      do i = 1, n
        tmp = ZERO
        do m = 1, 3
          tmp = tmp + grad_nu(i,m) * 2 * strain_rate(i,m,l)
        end do
        F_d(i,l) = tmp + nu(i,1) * lapl_v(i,l)
      end do
    end do

    !$omp master
    deallocate(v, lapl_v, strain_rate, nu, grad_nu)
    !$omp end master

  end subroutine GetDiffusiveTerm

  !-----------------------------------------------------------------------------
  !> Convective term, F_c = -∇⋅vv

  subroutine GetConvectiveTerm(n, x, t, F_c)
    integer,   intent(in)  :: n          !< number of points
    real(RNP), intent(in)  :: x(n,3)     !< mesh points
    real(RNP), intent(in)  :: t          !< time
    real(RNP), intent(out) :: F_c(n,3)   !< convective term

    real(RNP), allocatable, save :: grad_v(:,:,:), v(:,:)
    real(RNP) :: tmp
    integer   :: i, l, m

    !$omp master
    allocate(grad_v(n,3,3), v(n,3))
    !$omp end master
    !$omp barrier

    call GetVelocity        (n, x, t, v      )
    call GetVelocityGradient(n, x, t, grad_v )

    do m = 1, 3
      !$omp do
      do i = 1, n
        tmp = ZERO
        do l = 1, 3
          tmp = tmp + v(i,l) * grad_v(i,l,m)
        end do
        F_c(i,m) = -tmp
      end do
    end do

    !$omp master
    deallocate(grad_v, v)
    !$omp end master

  end subroutine GetConvectiveTerm

  !-----------------------------------------------------------------------------
  !> Strain rate tensor

  subroutine GetStrainRate(n, x, t, strain_rate)
    integer,   intent(in)  :: n                    !< number of points
    real(RNP), intent(in)  :: x(n,3)               !< mesh points
    real(RNP), intent(in)  :: t                    !< time
    real(RNP), intent(out) :: strain_rate(n,3,3)   !< strain rate tensor

    real(RNP), allocatable, save :: grad_v(:,:,:)
    integer   :: i, l, m

    !$omp single
    allocate(grad_v(n,3,3))
    !$omp end single

    call GetVelocityGradient(n, x, t, grad_v)

    do m = 1, 3
    do l = 1, 3
      !$omp do
      do i = 1, n
        ! symmetric part of velocity gradient tensor: S = 1/2 * (∇v + (∇v)ᵀ)
        strain_rate(i,l,m) = HALF * (grad_v(i,l,m) + grad_v(i,m,l))
      end do
    end do
    end do

    !$omp barrier
    !$omp master
    deallocate(grad_v)
    !$omp end master

  end subroutine GetStrainRate

!===============================================================================

end module INS__Problem__3D__Variable_Viscosity


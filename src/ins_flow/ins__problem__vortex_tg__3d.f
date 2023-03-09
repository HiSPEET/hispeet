!> summary:  Taylor-Green vortex adopted from Minion & Saye (2018)
!> author:   Joerg Stiller
!> date:     2018/06/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!> The solution is dimensionally inconsistent.
!> @endnote
!===============================================================================

module INS__Problem__Vortex_TG__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_Vortex_TG_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the Minion-Saye vortex problem

  type, extends(INS_Problem_3D) :: INS_Problem_Vortex_TG_3D

    real(RNP) :: vt(3) = 0  !< translation velocity
    real(RNP) :: xt(3) = 0  !< initial displacement

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetExternalSources

    procedure :: GetExactSolution
    procedure :: GetExactTimeDerivative
    procedure :: GetExactConvectiveTerm
    procedure :: GetExactPressureTerm
    procedure :: GetExactDiffusiveTerm

  end type INS_Problem_Vortex_TG_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_Vortex_TG_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    logical   :: stokes = .false.
    real(RNP) :: nu     =  0.01         ! kinematic viscosity
    real(RNP) :: vt(3)  =  [1,1,0]      ! translation velocity
    real(RNP) :: xt(3)  =  [0,1,0] / 8. ! initial displacement
    real(RNP) :: x0(3)  = -HALF         ! bounding box: corner nearest to -∞
    real(RNP) :: x1(3)  =  HALF         ! bounding box: corner nearest to +∞
    character, allocatable :: bc_v(:)

    namelist /parameters/ stokes, nu, vt, xt, x0, x1, bc_v

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

    ! default BC
    allocate(bc_v(nb), source = 'D')

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
      call XMPI_Bcast(bc_v  , 0, comm)
    end if

    problem % stokes = stokes
    problem % exact_solution = .true.

    ! set v_ref to reasonable approximation of initial velocity magnitude
    problem % v_ref  = sqrt( (sqrt(1 + vt(1)**2 + vt(2)**2))**2 + vt(3)**2 )
    problem % nu_ref = nu

    problem % vt = vt
    problem % xt = xt
    problem % x0 = x0
    problem % x1 = x1

    call move_alloc(bc_v, problem % bc_v)

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_Vortex_TG_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, vt => problem%vt, xt => problem%xt)
      call GetVelocity(n, x, ZERO, nu, vt, xt, u(:,:,:,:,1:3))
      call GetPressure(n, x, ZERO, nu, vt, xt, u(:,:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` in points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_Vortex_TG_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary identifier
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    integer :: m, n

    n = size(xb(:,:,:,1))

    associate(nu => problem%nu_ref, vt => problem%vt, xt => problem%xt)
      call GetVelocity(n, xb, t, nu, vt, xt, ub(:,:,:,1:3))
      call GetPressure(n, xb, t, nu, vt, xt, ub(:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(ub,4)
      call SetArray(ub(:,:,:,m), ZERO)
    end do

    ! ignore boundary identifier
    if (b < 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the external sources for all variables at points x and time t

  subroutine GetExternalSources(problem, x, t, F_s)
    class(INS_Problem_Vortex_TG_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_s(:,:,:,:,:) !< external sources

    integer :: m, n

    n = size(x(:,:,:,:,1))

    if (problem % stokes) then
      associate(nu => problem%nu_ref, vt => problem%vt, xt => problem%xt)
        call GetConvectiveTerm(n, x, t, nu, vt, xt, F_s(:,:,:,:,1:3))
        do m = 4, size(F_s,5)
          call SetArray(F_s(:,:,:,:,m), ZERO)
        end do
      end associate
    else
      call SetArray(F_s, ZERO, multi = .true.)
    end if

  end subroutine GetExternalSources

  !---------------------------------------------------------------------------
  !> Provides the exact solution u(x,t)

  subroutine GetExactSolution(problem, x, t, u)
    class(INS_Problem_Vortex_TG_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(out) :: u(:,:,:,:,:) !< solution

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, vt => problem%vt, xt => problem%xt)
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
    class(INS_Problem_Vortex_TG_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)    !< mesh points
    real(RNP), intent(in)  :: t               !< time
    real(RNP), intent(out) :: dt_u(:,:,:,:,:) !< time derivative

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, vt => problem%vt, xt => problem%xt)
      call GetVelocityTimeDerivative(n, x, t, nu, vt, xt, dt_u(:,:,:,:,1:3))
      call GetPressureTimeDerivative(n, x, t, nu, vt, xt, dt_u(:,:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(dt_u,5)
      call SetArray(dt_u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactTimeDerivative

  !---------------------------------------------------------------------------
  !> Exact convection term, F_c = -[∇⋅vv,0]

  subroutine GetExactConvectiveTerm(problem, x, t, F_c)
    class(INS_Problem_Vortex_TG_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_c(:,:,:,:,:) !< convection term

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, vt => problem%vt, xt => problem%xt)
      call GetConvectiveTerm(n, x, t, nu, vt, xt, F_c(:,:,:,:,1:3))
      do m = 4, size(F_c,5)
        call SetArray(F_c(:,:,:,:,m), ZERO)
      end do
    end associate

  end subroutine GetExactConvectiveTerm

  !---------------------------------------------------------------------------
  !> Exact pressure term, F_p = [-∇p,0]

  subroutine GetExactPressureTerm(problem, x, t, F_p)
    class(INS_Problem_Vortex_TG_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_p(:,:,:,:,:) !< diffusion term

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, vt => problem%vt, xt => problem%xt)
      call GetPressureTerm(n, x, t, nu, vt, xt, F_p(:,:,:,:,1:3))
      do m = 4, size(F_p,5)
        call SetArray(F_p(:,:,:,:,m), ZERO)
      end do
    end associate

  end subroutine GetExactPressureTerm

  !---------------------------------------------------------------------------
  !> Exact diffusion term, F_d = [∇⋅τ,0]

  subroutine GetExactDiffusiveTerm(problem, x, t, F_d)
    class(INS_Problem_Vortex_TG_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_d(:,:,:,:,:) !< diffusion term

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, vt => problem%vt, xt => problem%xt)
      call GetDiffusiveTerm(n, x, t, nu, vt, xt, F_d(:,:,:,:,1:3))
      do m = 4, size(F_d,5)
        call SetArray(F_d(:,:,:,:,m), ZERO)
      end do
    end associate

  end subroutine GetExactDiffusiveTerm

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
  !> Convection term

  subroutine GetConvectiveTerm(n, x, t, nu, vt, xt, F_c)
    integer,   intent(in)  :: n         !< number of points
    real(RNP), intent(in)  :: x(n,3)    !< mesh points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(in)  :: nu        !< kinematic viscosity
    real(RNP), intent(in)  :: vt(3)     !< translation velocity
    real(RNP), intent(in)  :: xt(3)     !< initial displacement
    real(RNP), intent(out) :: F_c(n,3)  !< convection term, -v · ∇v

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

      F_c(i,1) = 2*PI*c * ( -v1 * cos1*cos2 + v2 * sin1*sin2 )
      F_c(i,2) = 2*PI*c * ( -v1 * sin1*sin2 + v2 * cos1*cos2 )
      F_c(i,3) = ZERO

    end do

  end subroutine GetConvectiveTerm

  !-----------------------------------------------------------------------------
  !> Pressure term

  subroutine GetPressureTerm(n, x, t, nu, vt, xt, F_p)
    integer,   intent(in)  :: n         !< number of points
    real(RNP), intent(in)  :: x(n,3)    !< mesh points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(in)  :: nu        !< kinematic viscosity
    real(RNP), intent(in)  :: vt(3)     !< translation velocity
    real(RNP), intent(in)  :: xt(3)     !< initial displacement
    real(RNP), intent(out) :: F_p(n,3)  !< pressure term, -∇p

    real(RNP) :: a, c, phi1, phi2
    integer   :: i

    a = -16 * PI**2 * nu
    c =  PI * exp(a * t)

    !$omp do
    do i = 1, n

      phi1 = 4 * PI * (x(i,1) - vt(1)*t - xt(1))
      phi2 = 4 * PI * (x(i,2) - vt(2)*t - xt(2))

      F_p(i,1) = c * sin(phi1)
      F_p(i,2) = c * sin(phi2)
      F_p(i,3) = ZERO

    end do

  end subroutine GetPressureTerm

  !-----------------------------------------------------------------------------
  !> Diffusion term, F_d = [ν∇²v,0]

  subroutine GetDiffusiveTerm(n, x, t, nu, vt, xt, F_d)
    integer,   intent(in)  :: n         !< number of points
    real(RNP), intent(in)  :: x(n,3)    !< mesh points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(in)  :: nu        !< kinematic viscosity
    real(RNP), intent(in)  :: vt(3)     !< translation velocity
    real(RNP), intent(in)  :: xt(3)     !< initial displacement
    real(RNP), intent(out) :: F_d(n,3)  !< diffusion term

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

      F_d(i,1) =  a * c * sin1*cos2
      F_d(i,2) = -a * c * cos1*sin2
      F_d(i,3) =  ZERO

    end do

  end subroutine GetDiffusiveTerm

  !=============================================================================

end module INS__Problem__Vortex_TG__3D

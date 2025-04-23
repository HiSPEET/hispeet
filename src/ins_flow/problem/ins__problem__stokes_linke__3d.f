!> summary:  Stokes test case of Linke (2014)
!> author:   Joerg Stiller
!> date:     2023/03/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> 2D Stokes problem that is suited to investigate pressure robustness. The
!> velocity is defined by the stream function ψ = g(x)g(y) with g(x) = x²(x-1)²
!> and the pressure is p = x⁵ + y⁵ - 1/3. The source is set to f = ∇p - ν∇²v,
!> so that the solution is steady.
!>
!> ### References
!>
!> 1. Linke A, Comput Meth Appl Mech Engrg 268:782-800, 2014
!===============================================================================

module INS__Problem__Stokes_Linke__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, THIRD
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_Stokes_Linke_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the Linke Stokes problem

  type, extends(INS_Problem_3D) :: INS_Problem_Stokes_Linke_3D
  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetExternalSources

    procedure :: GetExactSolution
    procedure :: GetExactPressureTerm
    procedure :: GetExactDiffusiveTerm

  end type INS_Problem_Stokes_Linke_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_Stokes_Linke_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    real(RNP) :: nu    = 1 ! kinematic viscosity
    real(RNP) :: x0(3) = 0 ! bounding box corner nearest to to -∞
    real(RNP) :: x1(3) = 1 ! bounding box corner nearest to to +∞

    namelist /parameters/ nu, x0, x1

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

    if (rank == 0 .and. exists) then
      read(prm, nml=parameters)
      close(prm)
    end if

    if (present(comm)) then
      call XMPI_Bcast(nu, 0, comm)
      call XMPI_Bcast(x0, 0, comm)
      call XMPI_Bcast(x1, 0, comm)
    end if

    problem % stokes         = .true.
    problem % exact_solution = .true.

    problem % nu_ref = nu
    problem % v_ref  = ONE / 72
    problem % x0     = x0
    problem % x1     = x1

    allocate(problem % bc_v(nb), source = 'D')

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_Stokes_Linke_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetVelocity(n, x, u(:,:,:,:,1:3))
    call GetPressure(n, x, u(:,:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` in points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_Stokes_Linke_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary identifier
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    integer :: m, n

    n = size(xb(:,:,:,1))

    call GetVelocity(n, xb, ub(:,:,:,1:3))
    call GetPressure(n, xb, ub(:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(ub,4)
      call SetArray(ub(:,:,:,m), ZERO)
    end do

    ! ignore boundary identifier and time
    if (b < 0 .or. t > 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the external sources for all variables at points x and time t

  subroutine GetExternalSources(problem, x, t, F_s)
    class(INS_Problem_Stokes_Linke_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_s(:,:,:,:,:) !< external sources

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref)
      call GetSources(n, x, nu, F_s(:,:,:,:,1:3))
    end associate

    do m = 5, size(F_s,5)
      call SetArray(F_s(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExternalSources

  !---------------------------------------------------------------------------
  !> Provides the exact solution u(x,t)

  subroutine GetExactSolution(problem, x, t, u)
    class(INS_Problem_Stokes_Linke_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(out) :: u(:,:,:,:,:) !< solution

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetVelocity(n, x, u(:,:,:,:,1:3))
    call GetPressure(n, x, u(:,:,:,:,4)  )

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

    if (t > 0) return

  end subroutine GetExactSolution

  !---------------------------------------------------------------------------
  !> Exact pressure term, F_p = [-∇p,0]

  subroutine GetExactPressureTerm(problem, x, t, F_p)
    class(INS_Problem_Stokes_Linke_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_p(:,:,:,:,:) !< diffusion term

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetPressureTerm(n, x, F_p(:,:,:,:,1:3))
    do m = 4, size(F_p,5)
      call SetArray(F_p(:,:,:,:,m), ZERO)
    end do

    if (t > 0) return

  end subroutine GetExactPressureTerm

  !---------------------------------------------------------------------------
  !> Exact diffusion term, F_d = [∇⋅τ,0]

  subroutine GetExactDiffusiveTerm(problem, x, t, F_d)
    class(INS_Problem_Stokes_Linke_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_d(:,:,:,:,:) !< diffusion term

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref)
      call GetDiffusiveTerm(n, x, nu, F_d(:,:,:,:,1:3))
      do m = 4, size(F_d,5)
        call SetArray(F_d(:,:,:,:,m), ZERO)
      end do
    end associate

    if (t > 0) return

  end subroutine GetExactDiffusiveTerm

  !=============================================================================
  ! problem-specific procedures

  !-----------------------------------------------------------------------------
  !> Velocity

  subroutine GetVelocity(n, x, v)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(out) :: v(n,3)  !< velocity at mesh points

    real(RNP) :: x1, x2
    integer   :: i

    !$omp do
    do i = 1, n
      x1 = x(i,1)
      x2 = x(i,2)
      v(i,1) =  g0(x1) * g1(x2)
      v(i,2) = -g1(x1) * g0(x2)
      v(i,3) = ZERO
    end do

  end subroutine GetVelocity

  !-----------------------------------------------------------------------------
  !> Pressure

  subroutine GetPressure(n, x, p)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(out) :: p(n)    !< pressure at mesh points

    integer :: i

    !$omp do
    do i = 1, n
      p(i) = x(i,1)**5 + x(i,2)**5 - THIRD
    end do

  end subroutine GetPressure

  !-----------------------------------------------------------------------------
  !> Pressure term

  subroutine GetPressureTerm(n, x, F_p)
    integer,   intent(in)  :: n         !< number of points
    real(RNP), intent(in)  :: x(n,3)    !< mesh points
    real(RNP), intent(out) :: F_p(n,3)  !< pressure term, -∇p

    integer :: i

    !$omp do
    do i = 1, n
      F_p(i,1) = -5 * x(i,1)**4
      F_p(i,2) = -5 * x(i,2)**4
      F_p(i,3) =  ZERO
    end do

  end subroutine GetPressureTerm

  !-----------------------------------------------------------------------------
  !> Diffusion term, F_d = [ν∇²v,0]

  subroutine GetDiffusiveTerm(n, x, nu, F_d)
    integer,   intent(in)  :: n         !< number of points
    real(RNP), intent(in)  :: x(n,3)    !< mesh points
    real(RNP), intent(in)  :: nu        !< kinematic viscosity
    real(RNP), intent(out) :: F_d(n,3)  !< diffusion term

    real(RNP) :: x1, x2
    integer   :: i

    !$omp do
    do i = 1, n
      x1 = x(i,1)
      x2 = x(i,2)
      F_d(i,1) =  nu * (g2(x1) * g1(x2) + g0(x1) * g3(x2))
      F_d(i,2) = -nu * (g3(x1) * g0(x2) + g1(x1) * g2(x2))
      F_d(i,3) =  ZERO
    end do

  end subroutine GetDiffusiveTerm

  !-----------------------------------------------------------------------------
  !> Sources

  subroutine GetSources(n, x, nu, f)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: nu      !< kinematic viscosity
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(out) :: f(n,3)  !< source at mesh points

    real(RNP) :: x1, x2
    integer   :: i

    !$omp do
    do i = 1, n
      x1 = x(i,1)
      x2 = x(i,2)
      f(i,1) = 5 * x1**4 - nu * (g2(x1) * g1(x2) + g0(x1) * g3(x2))
      f(i,2) = 5 * x2**4 + nu * (g3(x1) * g0(x2) + g1(x1) * g2(x2))
      f(i,3) = ZERO
    end do

  end subroutine GetSources

  !-----------------------------------------------------------------------------
  !> 1D function used to build the stream function

  pure real(RNP) function g0(x)
    real(RNP), intent(in) :: x
    g0 = (x * (x - 1))**2
  end function g0

  !-----------------------------------------------------------------------------
  !> First derivative, g1 = g'(x)

  pure real(RNP) function g1(x)
    real(RNP), intent(in) :: x
    g1 = 2 * x * (x - 1) * (2*x - 1)
  end function g1

  !-----------------------------------------------------------------------------
  !> Second derivative, g2 = g"(x)

  pure real(RNP) function g2(x)
    real(RNP), intent(in) :: x
    g2 = 12 * x * (x - 1) + 2
  end function g2

  !-----------------------------------------------------------------------------
  !> Third derivative, g3 = g'''(x)

  pure real(RNP) function g3(x)
    real(RNP), intent(in) :: x
    g3 = 24 * x  - 12
  end function g3

  !=============================================================================

end module INS__Problem__Stokes_Linke__3D

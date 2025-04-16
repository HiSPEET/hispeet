!> summary:  Stokes test case of Deville, Kleiser & Montigny-Rannou (1984)
!> author:   Joerg Stiller
!> date:     2018/06/08, revised 2023/03/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> This is a 2D Stokes problem based on a family of analytic solutions derived
!> in [1]. The present test case was formulated in [2], but giving a wrong
!> expression for the exact solution. The correct form is found, e.g., in [3].
!> As the solution is 2π-periodic in x, the original test case was stated in the
!> domain [0,2π]×[-1,1] with periodicity in x and Dirichlet conditions in y
!> direction. Later studies used the domain [-1,1]² with Dirichlet conditions
!> throughout [3-5].
!> In HiSPEET, the problem is embedded 3D by adding z as the third dimension.
!>
!> 1. Deville MO, Kleiser L & Montigny-Rannou F, Int J Meth Fluids 4:1149-1163,
!>    1984
!>
!> 2. Maday Y, Patera A & Rønquist EM, J Sci Comput 5(4):263-292, 1990
!>
!> 3. Shahbazi K, Fischer PF & Ethier CR, J Comput Phys 222:391-407, 2007
!>
!> 4. Ferrer E, Moxey D, Willden RHJ &  Sherwin SJ, Commun Comput Phys
!>    16:817-840, 2014
!>
!> 5. Krank B, Fehn N, Wall WA & Kronbichler M, J Comput Phys 348:634-659, 2017
!===============================================================================

module INS__Problem__Stokes_DKM__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_Stokes_DKM_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the DKM Stokes problem

  type, extends(INS_Problem_3D) :: INS_Problem_Stokes_DKM_3D

    real(RNP) :: a = 2.883356_RNP !< wave length

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues

    procedure :: GetExactSolution
    procedure :: GetExactTimeDerivative
    procedure :: GetExactPressureTerm
    procedure :: GetExactDiffusiveTerm

  end type INS_Problem_Stokes_DKM_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_Stokes_DKM_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    logical   :: stokes = .true.
    real(RNP) :: nu    = 1                   ! kinematic viscosity
    real(RNP) :: a     = 2.883356_RNP        ! wave length

    ! bounding box corner nearest to ∓∞
    real(RNP) :: x0(3) = [ -real(PI, RNP), -ONE, -ONE ]
    real(RNP) :: x1(3) = [  real(PI, RNP),  ONE,  ONE ]

    namelist /parameters/ stokes, nu, a, x0, x1

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
      call XMPI_Bcast(stokes, 0, comm)
      call XMPI_Bcast(nu, 0, comm)
      call XMPI_Bcast(a , 0, comm)
      call XMPI_Bcast(x0, 0, comm)
      call XMPI_Bcast(x1, 0, comm)
    end if

    problem % stokes = stokes
    problem % exact_solution = .true.

    problem % nu_ref = nu
    problem % v_ref  = 3.46307
    problem % a      = a
    problem % x0     = x0
    problem % x1     = x1

    allocate(problem % bc_v(nb), source = 'D')

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_Stokes_DKM_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, a => problem%a)
      call GetVelocity(n, x, ZERO, nu, a, u(:,:,:,:,1:3))
      call GetPressure(n, x, ZERO, nu, a, u(:,:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` in points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_Stokes_DKM_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary identifier
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    integer :: m, n

    n = size(xb(:,:,:,1))

    associate(nu => problem%nu_ref, a => problem%a)
      call GetVelocity(n, xb, t, nu, a, ub(:,:,:,1:3))
      call GetPressure(n, xb, t, nu, a, ub(:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(ub,4)
      call SetArray(ub(:,:,:,m), ZERO)
    end do

    ! ignore boundary identifier
    if (b < 0) return

  end subroutine GetBoundaryValues

  !---------------------------------------------------------------------------
  !> Provides the exact solution u(x,t)

  subroutine GetExactSolution(problem, x, t, u)
    class(INS_Problem_Stokes_DKM_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(out) :: u(:,:,:,:,:) !< solution

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, a => problem%a)
      call GetVelocity(n, x, t, nu, a, u(:,:,:,:,1:3))
      call GetPressure(n, x, t, nu, a, u(:,:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactSolution

  !---------------------------------------------------------------------------
  !> Provides the time derivative of the exact solution, ∂u/∂t(x,t)

  subroutine GetExactTimeDerivative(problem, x, t, dt_u)
    class(INS_Problem_Stokes_DKM_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)    !< mesh points
    real(RNP), intent(in)  :: t               !< time
    real(RNP), intent(out) :: dt_u(:,:,:,:,:) !< time derivative

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, a => problem%a)
      call GetVelocityTimeDerivative(n, x, t, nu, a, dt_u(:,:,:,:,1:3))
      call GetPressureTimeDerivative(n, x, t, nu, a, dt_u(:,:,:,:,4)  )
    end associate

    ! remaining variables get zero
    do m = 5, size(dt_u,5)
      call SetArray(dt_u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetExactTimeDerivative

  !---------------------------------------------------------------------------
  !> Exact pressure term, F_p = [-∇p,0]

  subroutine GetExactPressureTerm(problem, x, t, F_p)
    class(INS_Problem_Stokes_DKM_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_p(:,:,:,:,:) !< diffusion term

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, a => problem%a)
      call GetPressureTerm(n, x, t, nu, a, F_p(:,:,:,:,1:3))
      do m = 4, size(F_p,5)
        call SetArray(F_p(:,:,:,:,m), ZERO)
      end do
    end associate

  end subroutine GetExactPressureTerm

  !---------------------------------------------------------------------------
  !> Exact diffusion term, F_d = [∇⋅τ,0]

  subroutine GetExactDiffusiveTerm(problem, x, t, F_d)
    class(INS_Problem_Stokes_DKM_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_d(:,:,:,:,:) !< diffusion term

    integer :: m, n

    n = size(x(:,:,:,:,1))

    associate(nu => problem%nu_ref, a => problem%a)
      call GetDiffusiveTerm(n, x, t, nu, a, F_d(:,:,:,:,1:3))
      do m = 4, size(F_d,5)
        call SetArray(F_d(:,:,:,:,m), ZERO)
      end do
    end associate

  end subroutine GetExactDiffusiveTerm

  !=============================================================================
  ! problem-specific procedures

  !-----------------------------------------------------------------------------
  !> Velocity

  subroutine GetVelocity(n, x, t, nu, a, v)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(in)  :: nu      !< kinematic viscosity
    real(RNP), intent(in)  :: a       !< wave length
    real(RNP), intent(out) :: v(n,3)  !< velocity at mesh points

    real(RNP) :: lambda, b, c, x1, x2
    integer   :: i

    lambda = nu * (1 + a*a)
    b = cos(a)
    c = exp(-lambda * t)

    !$omp do
    do i = 1, n
      x1 = x(i,1)
      x2 = x(i,2)
      v(i,1) = c * sin(x1) * (a * sin(a*x2) - b * sinh(x2))
      v(i,2) = c * cos(x1) * (    cos(a*x2) + b * cosh(x2))
      v(i,3) = ZERO
    end do

  end subroutine GetVelocity

  !-----------------------------------------------------------------------------
  !> Velocity time derivative

  subroutine GetVelocityTimeDerivative(n, x, t, nu, a, dt_v)
    integer,   intent(in)  :: n          !< number of points
    real(RNP), intent(in)  :: x(n,3)     !< mesh points
    real(RNP), intent(in)  :: t          !< time
    real(RNP), intent(in)  :: nu         !< kinematic viscosity
    real(RNP), intent(in)  :: a          !< wave length
    real(RNP), intent(out) :: dt_v(n,3)  !< velocity time derivative

    real(RNP) :: lambda, b, c, x1, x2
    integer   :: i

    lambda = nu * (1 + a*a)
    b = cos(a)
    c = -lambda * exp(-lambda * t)

    !$omp do
    do i = 1, n
      x1 = x(i,1)
      x2 = x(i,2)
      dt_v(i,1) = c * sin(x1) * (a * sin(a*x2) - b * sinh(x2))
      dt_v(i,2) = c * cos(x1) * (    cos(a*x2) + b * cosh(x2))
      dt_v(i,3) = ZERO
    end do

  end subroutine GetVelocityTimeDerivative

  !-----------------------------------------------------------------------------
  !> Pressure

  subroutine GetPressure(n, x, t, nu, a, p)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(in)  :: nu      !< kinematic viscosity
    real(RNP), intent(in)  :: a       !< wave length
    real(RNP), intent(out) :: p(n)    !< pressure at mesh points

    real(RNP) :: lambda, c
    integer   :: i

    lambda = nu * (1 + a*a)
    c = lambda * cos(a) * exp(-lambda * t)

    !$omp do
    do i = 1, n
      p(i) = c * cos(x(i,1)) * sinh(x(i,2))
    end do

  end subroutine GetPressure

  !-----------------------------------------------------------------------------
  !> Pressure time derivative

  subroutine GetPressureTimeDerivative(n, x, t, nu, a, dt_p)
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(in)  :: nu      !< kinematic viscosity
    real(RNP), intent(in)  :: a       !< wave length
    real(RNP), intent(out) :: dt_p(n) !< pressure time derivative

    real(RNP) :: lambda, c
    integer   :: i

    lambda = nu * (1 + a*a)
    c = -lambda**2 * cos(a) * exp(-lambda * t)

    !$omp do
    do i = 1, n
      dt_p(i) = c * cos(x(i,1)) * sinh(x(i,2))
    end do

  end subroutine GetPressureTimeDerivative

  !-----------------------------------------------------------------------------
  !> Pressure term

  subroutine GetPressureTerm(n, x, t, nu, a, F_p)
    integer,   intent(in)  :: n         !< number of points
    real(RNP), intent(in)  :: x(n,3)    !< mesh points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(in)  :: nu        !< kinematic viscosity
    real(RNP), intent(in)  :: a         !< wave length
    real(RNP), intent(out) :: F_p(n,3)  !< pressure term, -∇p

    real(RNP) :: lambda, c
    integer   :: i

    lambda = nu * (1 + a*a)
    c = lambda * cos(a) * exp(-lambda * t)

    !$omp do
    do i = 1, n
      F_p(i,1) =  c * sin(x(i,1)) * sinh(x(i,2))
      F_p(i,2) = -c * cos(x(i,1)) * cosh(x(i,2))
      F_p(i,3) =  ZERO
    end do

  end subroutine GetPressureTerm

  !-----------------------------------------------------------------------------
  !> Diffusion term, F_d = [ν∇²v,0]

  subroutine GetDiffusiveTerm(n, x, t, nu, a, F_d)
    integer,   intent(in)  :: n         !< number of points
    real(RNP), intent(in)  :: x(n,3)    !< mesh points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(in)  :: nu        !< kinematic viscosity
    real(RNP), intent(in)  :: a         !< wave length
    real(RNP), intent(out) :: F_d(n,3)  !< diffusion term

    real(RNP) :: lambda, a2, a3, c
    integer   :: i

    a2 = a * a
    a3 = a * a2
    lambda = nu * (1 + a2)
    c = -nu * exp(-lambda * t)

    !$omp do
    do i = 1, n
      F_d(i,1) = c * (a + a3) * sin(x(i,1)) * sin(a*x(i,2))
      F_d(i,2) = c * (1 + a2) * cos(x(i,1)) * cos(a*x(i,2))
      F_d(i,3) = ZERO
    end do

  end subroutine GetDiffusiveTerm

  !=============================================================================

end module INS__Problem__Stokes_DKM__3D

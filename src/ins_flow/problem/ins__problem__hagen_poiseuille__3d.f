!> summary:  Hagen Poiseuille flow
!> author:   Matthias Frey, Joerg Stiller
!> date:     2023/08/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> For Hagen-Poiseuille flow the Navier-Stokes equations simplify to
!>
!>       ρν (d²v/dr² + 1/r dv/dr) = -f
!>
!> with BC
!>
!>       v(r=R) = 0
!>
!> where v = v₃ and r² = x² + y². Here we use the nondimensional form
!>
!>       1/Re (d²v/dr² + 1/r dv/dr) = -f
!>
!> based on the bulk Reynolds number `Re = v_m R / ν`. For simplicity we
!> further impose `ρ = 1`, `R = 1` and `v_m = 1` so that `ν = 1/Re`.
!> With these assumptions the axial volume force becomes
!>
!>       f = 8 / Re
!>
!> and the exact solution takes the form
!>
!>       v(r) = 2(1 - r²)
!>
!===============================================================================

module INS__Problem__Hagen_Poiseuille__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, TWO, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_HagenPoiseuille_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the problem

  type, extends(INS_Problem_3D) :: INS_Problem_HagenPoiseuille_3D

    real(RNP) :: re     !< bulk Reynolds number, Re = 1/ν
    real(RNP) :: alpha  !< max relative perturbation

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetExternalSources
    procedure :: GetExactSolution

  end type INS_Problem_HagenPoiseuille_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_HagenPoiseuille_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    logical   :: stokes = .false. ! Navier-Stokes is default
    real(RNP) :: re     = 1       ! bulk flow Reynolds number, Re = 1/ν
    real(RNP) :: alpha  = 0       ! max relative perturbation
    character, allocatable :: bc_v(:)

    namelist /parameters/ stokes, re, alpha, bc_v

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

    if (rank == 0 .and. exists) then
      read(prm, nml=parameters)
      close(prm)
    end if

    if (present(comm)) then
      call XMPI_Bcast(stokes, 0, comm)
      call XMPI_Bcast(re    , 0, comm)
      call XMPI_Bcast(alpha , 0, comm)
      call XMPI_Bcast(bc_v  , 0, comm)
    end if

    problem % stokes          =  stokes
    problem % exact_solution  =  .true.
    problem % re              =  re
    problem % alpha           =  alpha
    problem % v_ref           =  ONE
    problem % nu_ref          =  ONE / re
    problem % x0              = -ONE
    problem % x1              =  ONE

    call move_alloc(bc_v, problem % bc_v)

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_HagenPoiseuille_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetVelocity( n, alpha = problem % alpha &
                    , x = x(:,:,:,:,1)           &
                    , y = x(:,:,:,:,2)           &
                    , v = u(:,:,:,:,1:3)         )

    ! remaining variables get zero
    do m = 4, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_HagenPoiseuille_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary identifier
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    integer :: m, n

    n = size(xb(:,:,:,1))

    call GetVelocity( n, alpha = ZERO   &
                    , x = xb(:,:,:,1)   &
                    , y = xb(:,:,:,2)   &
                    , v = ub(:,:,:,1:3) )

    ! remaining variables get zero
    do m = 4, size(ub,4)
      call SetArray(ub(:,:,:,m), ZERO)
    end do

    ! silence the compiler ;)
    if (b < 0 .or. t < 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the external sources for all variables at points x and time t

  subroutine GetExternalSources(problem, x, t, F_s)
    class(INS_Problem_HagenPoiseuille_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_s(:,:,:,:,:) !< external sources

    call SetArray(F_s(:,:,:,:,1), ZERO)
    call SetArray(F_s(:,:,:,:,2), ZERO)
    call SetArray(F_s(:,:,:,:,3), 8 / problem % re)

    ! silence the compiler ;)
    if (size(x) < 0 .or. t < 0) return

  end subroutine GetExternalSources

  !---------------------------------------------------------------------------
  !> Provides the exact solution u(x,t)

  subroutine GetExactSolution(problem, x, t, u)
    class(INS_Problem_HagenPoiseuille_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(out) :: u(:,:,:,:,:) !< solution

    integer :: m, n

    n = size(x(:,:,:,:,1))

    call GetVelocity (n, alpha = ZERO    &
                    , x = x(:,:,:,:,1)   &
                    , y = x(:,:,:,:,2)   &
                    , v = u(:,:,:,:,1:3) )

    ! remaining variables get zero
    do m = 4, size(u,5)
      call SetArray(u(:,:,:,:,m), ZERO)
    end do

    if (t < 0) return

  end subroutine GetExactSolution

  !=============================================================================
  ! Problem-specific procedures

  !-----------------------------------------------------------------------------
  !> Velocity

  subroutine GetVelocity(n, alpha, x, y, v)
    integer,   intent(in)  :: n      !< number of points
    real(RNP), intent(in)  :: alpha  !< max relative perturbation
    real(RNP), intent(in)  :: x(n)   !< x-coordinate of points
    real(RNP), intent(in)  :: y(n)   !< y-coordinate of points
    real(RNP), intent(out) :: v(n,3) !< velocity at mesh points

    real(RNP) :: a, rr, tol
    integer   :: i

    if (alpha > ZERO) then
      call random_number(v)  ! random numbers in the range [0,1)
    else
      call SetArray(v, ZERO, multi = .true.)
    end if

    tol = ONE - epsilon(ONE)

    !$omp do
    do i = 1, n

      rr = x(i)**2 + y(i)**2

      ! eliminate fluctuations on the wall
      if (rr > tol) then
        a = 0
      else
        a = alpha
      end if

      v(i,1) = a * (2*v(i,1) - 1)
      v(i,2) = a * (2*v(i,2) - 1)
      v(i,3) = a * (2*v(i,3) - 1) + 2 * (1 - rr)

    end do

  end subroutine GetVelocity

  !=============================================================================

end module INS__Problem__Hagen_Poiseuille__3D

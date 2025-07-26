!> summary:  Cylinder flow
!> author:   Joerg Stiller
!> date:     2025/07/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Problem__Linear_Cascade__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_LinearCascade_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the problem

  type, extends(INS_Problem_3D) :: INS_Problem_LinearCascade_3D

    real(RNP) :: beta !< inflow angle in DEG
    integer   :: init !< initialize with zero (0) or inflow velocity (1)

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues

  end type INS_Problem_LinearCascade_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_LinearCascade_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    character, allocatable :: bc_v(:)
    real(RNP) :: beta = 60.75_RNP
    real(RNP) :: v    = 1.000_RNP
    real(RNP) :: nu   = 0.001_RNP
    integer   :: init = 0

    namelist /parameters/ beta, v, nu, init, bc_v

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

    bc_v(2) = 'O'

    if (rank == 0 .and. exists) then
      read(prm, nml=parameters)
      close(prm)
    end if

    if (present(comm)) then
      call XMPI_Bcast(beta, 0, comm)
      call XMPI_Bcast(v   , 0, comm)
      call XMPI_Bcast(nu  , 0, comm)
      call XMPI_Bcast(init, 0, comm)
      call XMPI_Bcast(bc_v, 0, comm)
    end if

    problem % stokes          =  .false.
    problem % exact_solution  =  .false.
    problem % v_ref           =  v
    problem % nu_ref          =  nu
    problem % beta            =  beta

    call move_alloc(bc_v, problem % bc_v)
    call problem % SetPressureBC()

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_LinearCascade_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    real(RNP) :: v_x, v_y

    v_x = problem%v_ref * cos(problem%beta * PI/180)
    v_y = problem%v_ref * sin(problem%beta * PI/180)

    select case(problem % init)
    case(0)
      call SetArray(u, ZERO, multi = .true.)
    case(1)
      call SetArray(u(:,:,:,:,1), v_x)
      call SetArray(u(:,:,:,:,2), v_y)
      call SetArray(u(:,:,:,:,3:), ZERO, multi = .true.)
    end select

    if (size(x) == 0) return

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_LinearCascade_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary identifier
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    real(RNP) :: v_x, v_y

    if (b == 1) then
      v_x = problem%v_ref * cos(problem%beta * PI/180)
      v_y = problem%v_ref * sin(problem%beta * PI/180)
      call SetArray(ub(:,:,:,1), v_x)
      call SetArray(ub(:,:,:,2), v_y)
      call SetArray(ub(:,:,:,3:), ZERO, multi = .true.)
    else
      call SetArray(ub, ZERO, multi = .true.)
    end if

    if (size(xb) == 0 .or. t > 0) return

  end subroutine GetBoundaryValues

  !=============================================================================

end module INS__Problem__Linear_Cascade__3D

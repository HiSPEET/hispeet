!> summary:  Cylinder flow
!> author:   Matthias Frey, Joerg Stiller
!> date:     2023/10/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Problem__Cylinder_2D__3D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Problem_Cylinder2D_3D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the problem

  type, extends(INS_Problem_3D) :: INS_Problem_Cylinder2D_3D

    integer :: test_case !< 1 for steady, and 2 for unsteady flow

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues

  end type INS_Problem_Cylinder2D_3D

contains

  !=============================================================================
  ! type bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine SetProblem(problem, nb, file, comm)
    class(INS_Problem_Cylinder2D_3D), intent(inout) :: problem
    integer,                    intent(in) :: nb    !< number of boundaries
    character(len=*), optional, intent(in) :: file  !< input file
    type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

    ! local variables ..........................................................

    integer :: test_case = 1
    character, allocatable :: bc_v(:)

    namelist /parameters/ test_case, bc_v

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
      call XMPI_Bcast(test_case, 0, comm)
      call XMPI_Bcast(bc_v     , 0, comm)
    end if

    problem % stokes          =  .false.
    problem % exact_solution  =  .false.
    problem % nu_ref          =  1E-3_RNP
    problem % test_case       =  test_case

    call move_alloc(bc_v, problem % bc_v)

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  subroutine GetInitialValues(problem, x, u)
    class(INS_Problem_Cylinder2D_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

    call SetArray(u, ZERO, multi = .true.)

    if (problem % nc == 4 .or. size(x) == 0) return

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

  subroutine GetBoundaryValues(problem, b, xb, t, ub)
    class(INS_Problem_Cylinder2D_3D), intent(in) :: problem
    integer,   intent(in)  :: b             !< boundary identifier
    real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

    real(RNP), parameter :: hc = 0.41_RNP
    real(RNP), parameter :: vm = 0.30_RNP

    real(RNP) :: c0, c1, c2, y
    integer   :: np, nf, nc
    integer   :: i, j, k

    np = size(ub, 1)
    nf = size(ub, 3)
    nc = size(ub, 4)

    call SetArray(ub, ZERO, multi = .true.)

    if (b /= 1) return

    select case(problem % test_case)
    case(1)

      ! steady case ..........................................................

      c0 = 4 * vm / hc ** 2
      c1 = 0.20_RNP
      c2 = 0.21_RNP

      !$omp do collapse(3)
      do k = 1, nf
      do j = 1, np
      do i = 1, np
        y = xb(i,j,k,2)
        ub(i,j,k,1) = c0 * (y + c1) * (hc - (y + c2))
      end do
      end do
      end do

    case(2)

      ! unsteady case ........................................................

      c0 = real(6 * sin(PI * t/8) / 0.41_RNP ** 2, RNP)
      c1 = 0.20_RNP
      c2 = 0.21_RNP

      !$omp do collapse(3)
      do k = 1, nf
      do j = 1, np
      do i = 1, np
        y = xb(i,j,k,2)
        ub(i,j,k,1) = c0 * (y + c1) * (c2 - y)
      end do
      end do
      end do

    end select

  end subroutine GetBoundaryValues

  !=============================================================================

end module INS__Problem__Cylinder_2D__3D

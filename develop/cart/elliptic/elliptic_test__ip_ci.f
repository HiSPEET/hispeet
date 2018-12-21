!> summary:  Test elliptic IP/DG operator with constant isotropic diffusivity
!> author:   Joerg Stiller
!> date:     2018/12/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Test elliptic solvers
!===============================================================================

program Elliptic_Test__IP_CI

  use Kind_Parameters, only: RDP, RNP, IXL
  use Constants, only: PI, ZERO, ONE
  use Array_Assignments
  use TPO_sDDD
  use Export_Volume_Data_To_VTK

  use XMPI

  use CART__Mesh_Partition
  use CART__Generate_Structured_Mesh
  use CART__Boundary_Variable
  use CART__Schwarz_Operator
  use CART__Elliptic_Operator_IP

  use Elliptic_Test_Case

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! problem parameters .........................................................

  real(RNP) :: lambda  = 0       ! Helmholtz parameter
  real(RNP) :: nu      = 1       ! diffusivity
  integer   :: kappa   = 1       ! wave number
  real(RNP) :: xo(3)   = 0       ! corner closest to -infinity
  real(RNP) :: lx(3)   = 2*PI    ! domain extensions
  character :: bc(6)   = 'P'     ! boundary conditions {'P'|'D'|'N'}

  namelist /problem/ lambda, nu, kappa, xo, lx, bc

  ! discretization parameters ..................................................

  integer   :: np(3)   = 1       ! number of partitions in directions 1:3
  integer   :: ep(3)   = 2       ! elements per partition and direction

  integer   :: po      = 2       ! polynomial order
  real(RNP) :: penalty = 2       ! penalty parameter (> 1)

  namelist /discretization/ np, ep, po, penalty

  ! solution parameters ........................................................

  integer   :: method  = 1       ! CG,, Schwarz, p-MG or p-MG/CG {1|2|3|4}

  integer   :: i_max   = huge(1) ! max number of iterations/cycles
  real(RNP) :: r_red   = 1E-6    ! min residual reduction

  type(SchwarzOptions3D) :: schwarz_opt

  namelist /solver/ method
  namelist /solver_cg/ i_max, r_red
  namelist /solver_schwarz/ i_max, r_red, schwarz_opt

  ! MPI ........................................................................

  type(MPI_Comm) :: comm         ! communicator
  integer        :: comm_size    ! number of processes
  integer        :: rank         ! local rank

  ! mesh .......................................................................

  real(RNP)              :: dx(3)          ! spacing in directions 1:3
  logical                :: periodic(3)    ! periodic directions set true
  type(MeshPartition)    :: mesh           ! mesh partition
  real(RNP), allocatable :: x(:,:,:,:,:)   ! mesh points

  ! variables ..................................................................

  real(RNP), allocatable, target :: scalars(:)       ! storage for scalar fields
  character(len=80), allocatable :: scalar_names(:)  ! names of scalars

  real(RNP), dimension(:,:,:,:), pointer, contiguous :: s  ! exact solution
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: u  ! numeric solution
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: f  ! source
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: r  ! residual
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: e  ! error

  ! operators ..................................................................

  type(EllipticOperator3D_IP) :: elliptic_op

  ! input / output .............................................................

  character(len=80) :: parameter_file = 'elliptic_test__ip_ci'
  character(len=80) :: plot_file      = ''
  integer :: prm
  logical :: exists

  namelist /control/ plot_file

  ! auxiliary ..................................................................

  type(BoundaryVariable) :: bv(6)
  real(RNP), allocatable :: grad_u(:,:,:,:,:)
  real(RNP), allocatable :: laplace_u(:,:,:,:)
  real(RNP) :: r_max, r_max_loc
  real(RNP) :: e_min, e_min_loc
  real(RNP) :: e_max, e_max_loc
  real(RNP) :: c0
  real(RDP) :: time, time0

  integer :: nt = 10
  integer :: b, n
  integer :: i, ni
  integer(IXL) :: dof

  !-----------------------------------------------------------------------------
  ! Initialization

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_size(comm, comm_size)
  call MPI_Comm_rank(comm, rank)

  if (rank == 0) then

    parameter_file = trim(parameter_file) // '.prm'
    inquire(file=parameter_file, exist=exists)
    if (exists) then
      open(newunit=prm, file=parameter_file)
      read(prm, nml=problem)
      read(prm, nml=discretization)
      read(prm, nml=solver)
      select case(method)
      case(1)
        read(prm, nml=solver_cg)
      case(2)
        read(prm, nml=solver_schwarz)
      !case(3:)
      !  read(prm, nml=solver_pmg)
      end select
      read(prm, nml=control)
      close(prm)
    end if

  end if

  ! problem parameters
  call XMPI_Bcast(lambda   , 0, comm)
  call XMPI_Bcast(nu       , 0, comm)
  call XMPI_Bcast(kappa    , 0, comm)
  call XMPI_Bcast(xo       , 0, comm)
  call XMPI_Bcast(lx       , 0, comm)
  call XMPI_Bcast(bc       , 0, comm)

  ! discretization parameters
  call XMPI_Bcast(np       , 0, comm)
  call XMPI_Bcast(ep       , 0, comm)
  call XMPI_Bcast(po       , 0, comm)
  call XMPI_Bcast(penalty  , 0, comm)

  ! solver
  call XMPI_Bcast(method   , 0, comm)

  ! CG/Schwarz options
  call XMPI_Bcast(i_max    , 0, comm)
  call XMPI_Bcast(r_red    , 0, comm)

  ! Schwarz options
  call schwarz_opt % Bcast(0, comm)

  ! control parameters
  call XMPI_Bcast(plot_file, 0, comm)

  ! mesh and variables .........................................................

  ! spacing
  dx = lx/ (np * ep)

  ! periodicity
  periodic(1) = all(bc(1:2) == 'P')
  periodic(2) = all(bc(3:4) == 'P')
  periodic(3) = all(bc(5:6) == 'P')

  ! mesh partition and points
  call GenerateStructuredMesh(mesh, np, ep, xo, dx, periodic, comm)
  call mesh % GetPoints(po, 'GLL', x)

  ! mesh variables
  call InitializeMeshVariables()

  ! auxiliary variables
  allocate(grad_u(0:po, 0:po, 0:po, mesh%ne, 3))
  allocate(laplace_u(0:po, 0:po, 0:po, mesh%ne))

  ! operators ..................................................................

  call elliptic_op % New(mesh, lambda, nu, bc, po, penalty, schwarz_opt)

  !-----------------------------------------------------------------------------
  ! Tests

  ! exact solution and RHS  ....................................................

  n = size(u)

  if (n > 0) then

    ! exact solution, gradient and Laplacian
    call GetExactSolution(  kappa, n, x, u           )
    call GetExactGradient(  kappa, n, x, grad_u      )
    call GetExactLaplacian( kappa, n, x, laplace_u=r )

    ! s = u
    call AssignArray(s, u)

    ! r = -nu laplace u + lambda u
    call MergeArrays(-nu, r, lambda, u)

    ! project source:  f = M r
    c0 = product(dx) / 8
    call TPO_sDDD_Eval(po+1, mesh%ne, c0, elliptic_op%eop%w, r, f)

    ! boundary values
    do b = 1, size(bc)
      select case(bc(b))
      case('D')
        call bv(b) % Extract(mesh, u, b, bc(b))
      case('N')
        call bv(b) % ExtractNormalComponent(mesh, grad_u, b, bc(b))
      case default
        call bv(b) % New(mesh, po, b, bc(b))
      end select
    end do

    call elliptic_op % BcToRHS(bv, f)

  end if

  ! operator ...................................................................

  if (rank == 0) then
    write(*,'(/,A)') repeat('-',80)
    write(*,'(A,/)') 'IP/DG EllipticOperator: Apply'
  end if

  !$omp parallel
  !$acc data copyin(u) copyout(r)

  ! setup call
  call elliptic_op % Apply(u, r)
  !$acc wait

  if (rank == 0) then
    time0 = MPI_Wtime()
  end if

  do i = 1, nt
    call elliptic_op % Apply(u, r)
    !$acc wait
  end do

  !$acc end data
  !$omp end parallel

  if (rank == 0) then
    time = MPI_Wtime()
    time = (time - time0) / nt
    dof  = product(np) * product(ep) * (po + 1)**3
    write(*,'(A,1X,I0,A,ES10.3,A)') 'dof      =', dof, ' (', real(dof), ' )'
    write(*,'(A,ES10.3)')           'time/dof =', time / dof
    write(*,'(A,ES10.3)')           'dof/time =', dof / time
  end if

  ! residual ...................................................................

  if (rank == 0) then
    write(*,'(/,A)') repeat('-',80)
    write(*,'(A,/)') 'IP/DG EllipticOperator: Residual'
  end if

  !$omp parallel
  !$acc data copyin(u,f) copyout(r)

  call elliptic_op % Residual(u, f, r)

  !$acc end data
  !$omp end parallel

  if (mesh%part >= 0) then
    r_max_loc = maxval(abs(r))
  else
    r_max_loc = 0
  end if
  call XMPI_Reduce(r_max_loc, r_max, MPI_MAX, 0, mesh%comm)

  if (rank == 0) then
    write(*,'(A,ES10.3)')  'consistency:  r_max =', r_max
  end if

  ! solution ...................................................................

  if (rank == 0) then
    write(*,'(/,A)') repeat('-',80)
    select case(method)
    case(1)
      write(*,'(A,/)') 'IP/DG EllipticOperator: Conjugate Gradients'
    case(2)
      write(*,'(A,/)') 'IP/DG EllipticOperator: Schwarz Method'
    case(3)
      write(*,'(A,/)') 'IP/DG EllipticOperator: p-Multigrid'
    case(4)
      write(*,'(A,/)') 'IP/DG EllipticOperator: p-MG/CG'
    end select
  end if

  !$omp parallel
  !$acc data copyin(f) copyout(u) create(r)

  !call AssignScalar(u, ZERO)
  call random_number(u)
  u = 2*u - 1

  call elliptic_op % Residual(u, f, r)
  if (mesh%part >= 0) then
    r_max_loc = maxval(abs(r))
  else
    r_max_loc =  0
  end if
  call XMPI_Reduce(r_max_loc, r_max, MPI_MAX, 0, mesh%comm)
  if (rank == 0) then
    write(*,'(A,ES10.3)') 'initial residual:  r_0 =', r_max
    write(*,*)
  end if

  !$acc end data
  !$omp end parallel

  select case(method)

  case(1) ! conjugate gradients

    call elliptic_op % ConjugateGradients(u, f, i_max, r_red, ni=ni)

  end select

  call elliptic_op % Residual(u, f, r)

  if (mesh%part >= 0) then
    r_max_loc = maxval(abs(r))
    e = u - s
    e_min_loc = minval(e)
    e_max_loc = maxval(e)
  else
    r_max_loc =  0
    e_min_loc = -huge(ONE)
    e_max_loc =  huge(ONE)
  end if
  call XMPI_Reduce(r_max_loc, r_max, MPI_MAX, 0, mesh%comm)
  call XMPI_Reduce(e_min_loc, e_min, MPI_MIN, 0, mesh%comm)
  call XMPI_Reduce(e_max_loc, e_max, MPI_MAX, 0, mesh%comm)

  if (rank == 0) then
    write(*,'(A,1X,I0)')  'iterations:   ni    =', ni
    write(*,'(A,ES10.3)') 'residual:     r_max =', r_max
    write(*,'(A,ES10.3)') 'error:        e_max =', (e_max - e_min)/2
    write(*,*)
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

!-------------------------------------------------------------------------------
!> Initialization of mesh variables

subroutine InitializeMeshVariables()

  integer :: ns = 5  ! number of scalar fields
  integer :: ls      ! length of one scalar field
  integer :: i, j, k

  ! provide memory .............................................................

  ls = (po + 1)**3 * mesh%ne

  allocate( scalar_names(ns) )
  allocate( scalars(ls * ns) )

  ! assign scalars .............................................................

  i = 1
  j = 1
  k = ls

  scalar_names(i) = 's'
  s(0:po, 0:po, 0:po, 1:mesh%ne) => scalars(j:k)

  i = i + 1
  j = j + ls
  k = k + ls

  scalar_names(i) = 'u'
  u(0:po, 0:po, 0:po, 1:mesh%ne) => scalars(j:k)

  i = i + 1
  j = j + ls
  k = k + ls

  scalar_names(i) = 'f'
  f(0:po, 0:po, 0:po, 1:mesh%ne) => scalars(j:k)

  i = i + 1
  j = j + ls
  k = k + ls

  scalar_names(i) = 'r'
  r(0:po, 0:po, 0:po, 1:mesh%ne) => scalars(j:k)

  i = i + 1
  j = j + ls
  k = k + ls

  scalar_names(i) = 'e'
  e(0:po, 0:po, 0:po, 1:mesh%ne) => scalars(j:k)

end subroutine InitializeMeshVariables

!===============================================================================

end program Elliptic_Test__IP_CI

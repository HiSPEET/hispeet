!> summary:  Test elliptic IP/DG operator with constant isotropic diffusivity
!> author:   Joerg Stiller, Gustav Tschirschnitz
!> date:     2018/12/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Test elliptic solvers
!===============================================================================

program Elliptic_Test__IP_CI

  use Kind_Parameters, only: RDP, RNP, IXL
  use Constants, only: PI, ZERO, ONE
  use Array_Assignments
  use Array_Reductions
  use TPO__Diagonal__3D
  use OpenMP_Binding
  use XMPI
  use DG__Element_Operators__1D
  use Export_Volume_Data_To_VTK

  use CART__Mesh_Partition
  use CART__Generate_Structured_Mesh
  use CART__Boundary_Variable
  use CART__Schwarz_Operator
  use CART__Elliptic_Operator_IP
  use CART__Elliptic_PMG

  use Elliptic_Problem
  use Elliptic_Problem__Simple_1D
  use Elliptic_Problem__Simple_2D
  use Elliptic_Problem__Simple_3D
  use Elliptic_Problem__Knotty

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! problem ....................................................................

  class(EllipticProblem), allocatable :: problem

  integer   :: test      =  1     ! case {1,2,3,4} = {simple 1d/2d/3d, knotty}

  real(RNP) :: lambda    =  0      ! Helmholtz parameter
  real(RNP) :: nu        =  1      ! diffusivity
  real(RNP) :: nu_svv    = -1      ! spectral diffusivity amplitude
  integer   :: k_u       =  1      ! solution wave number

  real(RNP) :: xo(3)     =  0      ! corner closest to -infinity
  real(RNP) :: lx(3)     =  2*PI   ! domain extensions
  character :: bc(6)     = 'P'     ! boundary conditions {'P'|'D'|'N'}

  namelist /test_case/ test
  namelist /test_case/ lambda, nu, nu_svv, k_u
  namelist /test_case/ xo, lx, bc

  ! discretization .............................................................

  integer   :: np(3)     = 1       ! number of partitions in directions 1:3
  integer   :: ep(3)     = 2       ! elements per partition and direction
  integer   :: po        = 2       ! polynomial order
  real(RNP) :: penalty   = 2       ! penalty parameter > 1
  logical   :: adjust_dx = .false. ! adjust mesh spacing: Δx₃ = max(Δx₁,Δx₂)
  logical   :: svv       = .false. ! switch for the SVV model

  namelist /discretization/ np, ep, po, penalty, adjust_dx, svv

  ! solution ...................................................................

  integer :: method = 1
  ! 0  none
  ! 1  CG
  ! 2  Schwarz
  ! 3  Schwarz-PCG
  ! 4  p-MG
  ! 5  p-MG/CG

  integer   :: i_max     = huge(1) ! max number of iterations/cycles
  real(RNP) :: r_red     = 1E-6    ! min residual reduction

  type(SchwarzOptions3D) :: schwarz_opt
  type(PMG_Options3D)    :: pmg_opt

  namelist /solver/             method
  namelist /solver_cg/          i_max, r_red
  namelist /solver_schwarz/     i_max, r_red, schwarz_opt
  namelist /solver_schwarz_pcg/ i_max, r_red, schwarz_opt
  namelist /solver_pmg/         pmg_opt

  ! MPI and OpenMP .............................................................

  type(MPI_Comm) :: comm         ! MPI communicator
  integer        :: rank         ! local MPI rank
  integer        :: n_proc       ! number of MPI processes
  integer        :: n_thread     ! number of OpenMP threads

  ! mesh .......................................................................

  real(RNP)              :: dx(3)          ! spacing in directions 1:3
  logical                :: periodic(3)    ! periodic directions set true
  type(MeshPartition)    :: mesh           ! mesh partition
  real(RNP), allocatable :: x(:,:,:,:,:)   ! mesh points

  ! variables ..................................................................

  real(RNP), allocatable, target :: scalars(:)       ! storage for scalar fields
  character(len=80), allocatable :: scalar_names(:)  ! names of scalars

  real(RNP), pointer, contiguous :: s(:,:,:,:)  ! exact solution
  real(RNP), pointer, contiguous :: u(:,:,:,:)  ! numeric solution
  real(RNP), pointer, contiguous :: f(:,:,:,:)  ! source
  real(RNP), pointer, contiguous :: r(:,:,:,:)  ! residual
  real(RNP), pointer, contiguous :: e(:,:,:,:)  ! error

  ! operators / methods ........................................................

  type(EllipticOperator3D_IP) :: elliptic_op

  type(PMG_Method3D) :: pmg

  ! input / output .............................................................

  character(len=80) :: parameter_file = 'elliptic_test__ip_ci'
  character(len=80) :: plot_file      = ''
  integer :: nt = 10
  integer :: prm
  logical :: exists

  namelist /control/ nt, plot_file

  ! auxiliary ..................................................................

  type(DG_ElementOptions_1D) :: ip_opt
  type(BoundaryVariable) :: bv(6)
  real(RNP), allocatable :: grad_u(:,:,:,:,:)
  real(RNP), allocatable :: laplace_u(:,:,:,:)
  real(RNP) :: r_max, r_max_loc, r_l2, r_l2_0
  real(RNP) :: e_min, e_min_loc
  real(RNP) :: e_max, e_max_loc
  real(RNP) :: c0
  real(RDP) :: time, time0

  integer :: b, n
  integer :: i, ni
  integer(IXL) :: dof

  !-----------------------------------------------------------------------------
  ! Initialization

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  !$omp parallel
  n_thread = OMP_Num_Threads()
  !$omp end parallel

  if (rank == 0) then

    parameter_file = trim(parameter_file) // '.prm'
    inquire(file=parameter_file, exist=exists)
    if (exists) then
      open(newunit=prm, file=parameter_file)
      read(prm, nml=test_case)
      read(prm, nml=discretization)
      read(prm, nml=solver)
      select case(method)
      case(1)
        read(prm, nml=solver_cg)
      case(2)
        read(prm, nml=solver_schwarz)
      case(3)
        read(prm, nml=solver_schwarz_pcg)
      case(4,5)
        read(prm, nml=solver_pmg)
        pmg_opt % po_top = po
      end select
      read(prm, nml=control)
      close(prm)
    end if

  end if

  ! problem
  call XMPI_Bcast(test   , 0, comm)
  call XMPI_Bcast(lambda , 0, comm)
  call XMPI_Bcast(nu     , 0, comm)
  call XMPI_Bcast(nu_svv , 0, comm)
  call XMPI_Bcast(k_u    , 0, comm)
  call XMPI_Bcast(xo     , 0, comm)
  call XMPI_Bcast(lx     , 0, comm)
  call XMPI_Bcast(bc     , 0, comm)

  ! discretization
  call XMPI_Bcast(np        , 0, comm)
  call XMPI_Bcast(ep        , 0, comm)
  call XMPI_Bcast(po        , 0, comm)
  call XMPI_Bcast(penalty   , 0, comm)
  call XMPI_Bcast(adjust_dx , 0, comm)

  ! solver
  call XMPI_Bcast(method , 0, comm)
  call XMPI_Bcast(svv    , 0, comm)

  ! CG/Schwarz options
  call XMPI_Bcast(i_max , 0, comm)
  call XMPI_Bcast(r_red , 0, comm)

  ! Schwarz, PCG and PMG options
  call schwarz_opt % Bcast(0, comm)
  call pmg_opt     % Bcast(0, comm)

  ! control parameters
  call XMPI_Bcast(nt        , 0, comm)
  call XMPI_Bcast(plot_file , 0, comm)

  ! problem ....................................................................

  select case(test)
  case(1)
    allocate(EllipticProblem_Simple1D :: problem)
  case(2)
    allocate(EllipticProblem_Simple2D :: problem)
  case(3)
    allocate(EllipticProblem_Simple3D :: problem)
  case default
    allocate(EllipticProblem_Knotty   :: problem)
  end select

  call problem % SetProblem(lambda, nu_0 = nu, k_u = k_u)

  ! mesh and variables .........................................................

  ! spacing
  dx = lx / (np * ep)

  if (adjust_dx) then
    dx(3) = maxval(dx(1:2))
  end if

  ! periodicity
  periodic(1) = all(bc(1:2) == 'P')
  periodic(2) = all(bc(3:4) == 'P')
  periodic(3) = all(bc(5:6) == 'P')

  ! mesh partition and points
  call GenerateStructuredMesh(mesh, np, ep, xo, dx, periodic, comm)
  call mesh % GetPoints(po, 'L', x)

  ! mesh variables
  call InitializeMeshVariables()
  n = size(u)

  ! auxiliary variables
  allocate(grad_u(0:po, 0:po, 0:po, mesh%ne, 3))
  allocate(laplace_u(0:po, 0:po, 0:po, mesh%ne))

  ! set spectral diffusivity to default value if not given as parameter
  if (nu_svv == -1) nu_svv = ONE / real(po,RNP)

  ! operators ..................................................................

  ip_opt = DG_ElementOptions_1D(po = po, penalty = penalty, svv = svv)

  if (svv) then
    elliptic_op = EllipticOperator3D_IP(mesh, lambda, nu, nu_svv, bc, ip_opt,  &
                                        schwarz_opt)
  else
    elliptic_op = EllipticOperator3D_IP(mesh, lambda, nu,         bc, ip_opt,  &
                                        schwarz_opt)
  end if

  !-----------------------------------------------------------------------------
  ! Tests

  ! exact solution and RHS  ....................................................

  if (n > 0) then

    ! exact solution, gradient and Laplacian
    call problem % GetExactSolution  (x, u           )
    call problem % GetExactGradient  (x, grad_u      )
    call problem % GetExactLaplacian (x, laplace_u=r )

    ! s = u
    call SetArray(s, u)

    ! r = -ν ∇²u +  λ u
    call MergeArrays(-nu, r, lambda, u)

    ! project source:  f = M r
    c0 = product(dx) / 8
    call TPO_Diagonal(c0, elliptic_op%eop%w, r, f)

    ! boundary values
    do b = 1, size(bc)
      select case(bc(b))
      case('D')
        call bv(b) % Extract(mesh, u, b, bc(b))
      case('N')
        call bv(b) % ExtractNormalComponent(mesh, grad_u, b, bc(b))
      case default
        bv(b) = BoundaryVariable(mesh, po, b, bc(b))
      end select
    end do

    call elliptic_op % BcToRHS(bv, 1, f)

  end if

  ! methods
  select case(method)
  case(4,5)
    pmg = PMG_Method3D(mesh, ip_opt, pmg_opt)
    if (svv) then
      call pmg % SetProblem(lambda, nu, nu_svv, bc)
    else
      call pmg % SetProblem(lambda, nu,         bc)
    end if
  end select

  if (rank == 0) then
    dof = product(np) * product(ep) * (po + 1)**3
    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'IP/DG Elliptic Solver with constant diffusivity'
    write(*,*)
    write(*,'(2X,A,1X,I0)')        'n_proc   = ', n_proc
    write(*,'(2X,A,1X,I0)')        'n_thread = ', n_thread
    write(*,*)
    write(*,'(2X,A,1X,I0)')        'P   = ', po
    write(*,'(2X,A,1X,I0)')        'ne  = ', product(np) * product(ep)
    write(*,'(2X,A,1X,I0)')        'DOF = ', dof
    write(*,'(2X,A,3(ES12.5,1X))') 'dx  = ', dx
    write(*,*)
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
    write(*,'(2X,A,ES10.3)') 'time/DOF =', time / dof
    write(*,'(2X,A,ES10.3)') 'DOF/time =', dof / time
  end if

  ! residual ...................................................................

  if (rank == 0) then
    write(*,'(/,A)') repeat('-',80)
    write(*,'(A,/)') 'IP/DG EllipticOperator: Consistency'
  end if

  !$omp parallel
  !$acc data copyin(u,f) copyout(r)

  call elliptic_op % Residual(u, f, r)
  r_l2 = ScalarProduct(r, r, mesh%comm)
  r_l2 = sqrt(r_l2)

  !$acc end data
  !$omp end parallel

  if (mesh%part >= 0) then
    r_max_loc = maxval(abs(r))
  else
    r_max_loc = 0
  end if
  call XMPI_Reduce(r_max_loc, r_max, MPI_MAX, 0, mesh%comm)

  if (rank == 0) then
    write(*,'(2X,A,ES10.3)') 'r_L2  =', r_l2
    write(*,'(2X,A,ES10.3)') 'r_max =', r_max
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
      write(*,'(A,/)') 'IP/DG EllipticOperator: Schwarz-PCG Method'
    case(4)
      write(*,'(A,/)') 'IP/DG EllipticOperator: p-Multigrid'
    case(5)
      write(*,'(A,/)') 'IP/DG EllipticOperator: p-MG/CG'
    case default
      write(*,'(A,/)') 'skipping solver test'
    end select
  end if

  if (method == 0) then
    call MPI_Finalize()
    stop
  end if

  !$omp parallel
  !$acc data copyin(f) copyout(u) create(r)

  call SetArray(u, ZERO)
  !call random_number(u)
  !u = 2*u - 1

  call elliptic_op % Residual(u, f, r)
  r_l2_0 = ScalarProduct(r, r, mesh%comm)
  r_l2_0 = sqrt(r_l2_0)

  !$omp master
  if (mesh%part >= 0) then
    r_max_loc = maxval(abs(r))
  else
    r_max_loc =  0
  end if
  call XMPI_Reduce(r_max_loc, r_max, MPI_MAX, 0, mesh%comm)
  if (rank == 0) then
    write(*,'(A)') 'initial residual:'
    write(*,'(2X,A,ES10.3)') 'r_L2  =', r_l2_0
    write(*,'(2X,A,ES10.3)') 'r_max =', r_max
    select case(method)
    case(4:5)
      if (pmg_opt % monitor) write(*,*)
    end select
  end if

  if (rank == 0) then
    time0 = MPI_Wtime()
  end if
  !$omp end master

  select case(method)
  case(1) ! conjugate gradients
    call elliptic_op % ConjugateGradients(u, f, i_max, r_red, ni=ni)
  case(2) ! Schwarz method
    call elliptic_op % SchwarzMethod(u, f, i_max, r_red, ni=ni)
  case(3) ! Schwarz-PCG method
    call elliptic_op % SchwarzPreConjugateGradients(u, f, i_max, r_red, ni=ni)
  case(4) ! p-MG method
    call pmg % MG_Solver(u, f, ni=ni)
  case(5) ! p-MG/CG method
    call pmg % MG_CG_Solver(u, f, ni=ni)
  end select

  !$omp master
  if (rank == 0) then
    time = MPI_Wtime()
    time = time - time0
  end if
  !$omp end master

  call elliptic_op % Residual(u, f, r)
  r_l2 = ScalarProduct(r, r, mesh%comm)
  r_l2 = sqrt(r_l2)

  !$acc end data
  !$omp end parallel

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
    write(*,'(/,A)')      'solution:'
    write(*,'(A,1X,I0)')  '   ni    =', ni
    write(*,'(A,ES10.3)') '   r_L2  =', r_l2
    write(*,'(A,ES10.3)') '   r_max =', r_max
    write(*,'(A,ES10.3)') '   e_max =', (e_max - e_min)/2
    if (ni > 0) then
      write(*,'(A,ES10.3)') '   -lg ρ =', log10(r_l2_0 / r_l2) / ni
    end if
    write(*,'(/,A)')      'performance:'
    write(*,'(A,ES10.3)') '   time     =', time
    write(*,'(A,ES10.3)') '   time/DOF =', time / dof
    write(*,'(A,ES10.3)') '   DOF/time =', dof / time
    write(*,*)
  end if

  !  plot file .................................................................

  if (len_trim(plot_file) > 0 .and. mesh%part >= 0) then
    call ExportVolumeDataToVTK( po, mesh%ne,              &
                                size(scalar_names),       &
                                0,                        &
                                x,                        &
                                scalars,                  &
                                scalar_names,             &
                                file   = trim(plot_file), &
                                part   = mesh%part,       &
                                n_part = mesh%n_part      )
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

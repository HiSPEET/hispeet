!> summary:  Test elliptic IP/DG operator with variable isotropic diffusivity
!> author:   Joerg Stiller
!> date:     2019/01/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Test elliptic solvers
!===============================================================================

program Elliptic_Test__IP_VI

  use Kind_Parameters, only: RDP, RNP, IXL
  use Constants, only: PI, ZERO, ONE
  use Array_Assignments
  use Array_Reductions
  use TPO_sDDD
  use XMPI
  use IP_Element_Operators_1D
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

  integer   :: test      = 1       ! case {1,2,3,4} = {simple 1d/2d/3d, knotty}

  real(RNP) :: lambda    = 0       ! Helmholtz parameter
  real(RNP) :: nu_0      = 1       ! diffusivity mean
  real(RNP) :: nu_1      = 0.1     ! diffusivity fluctuation
  integer   :: k_nu      = 1       ! diffusivity wave number
  real(RNP) :: d_nu      = 0       ! diffusivity phase shift
  integer   :: k_u       = 1       ! solution wave number

  ! domain and boundary conditions
  real(RNP) :: xo(3)     = 0       ! corner closest to -infinity
  real(RNP) :: lx(3)     = 2*PI    ! domain extensions
  character :: bc(6)     = 'P'     ! boundary conditions {'P'|'D'|'N'}

  namelist /test_case/ test
  namelist /test_case/ lambda, nu_0, nu_1, d_nu, k_nu, k_u
  namelist /test_case/ xo, lx, bc

  ! discretization .............................................................

  integer   :: np(3)     = 1       ! number of partitions in directions 1:3
  integer   :: ep(3)     = 2       ! elements per partition and direction
  integer   :: po        = 2       ! polynomial order
  integer   :: penalty   = 2       ! penalty parameter > 1
  logical   :: adjust_dx = .false. ! adjust mesh spacing: Δx₃ = max(Δx₁,Δx₂)

  namelist /discretization/ np, ep, po, penalty, adjust_dx

  ! solution ...................................................................

  integer   :: method  = 1       ! CG,, Schwarz, p-MG or p-MG/CG {1|2|3|4}

  integer   :: i_max   = huge(1) ! max number of iterations/cycles
  real(RNP) :: r_red   = 1E-6    ! min residual reduction

  type(SchwarzOptions3D) :: schwarz_opt
  type(PMG_Options3D)    :: pmg_opt

  namelist /solver/         method
  namelist /solver_cg/      i_max, r_red
  namelist /solver_schwarz/ i_max, r_red, schwarz_opt
  namelist /solver_pmg/     pmg_opt

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
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: nu ! diffusivity
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: f  ! source
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: r  ! residual
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: e  ! error

  ! operators ..................................................................

  type(EllipticOperator3D_IP) :: elliptic_op

  type(PMG_Method3D) :: pmg

  ! input / output .............................................................

  character(len=80) :: parameter_file = 'elliptic_test__ip_vi'
  character(len=80) :: plot_file      = ''
  integer :: prm
  logical :: exists

  namelist /control/ plot_file

  ! auxiliary ..................................................................

  type(IP_ElementOptions1D) :: ip_opt
  type(BoundaryVariable) :: bv(6)
  real(RNP), allocatable :: grad_u(:,:,:,:,:)
  real(RNP) :: r_max, r_max_loc, r_l2, r_l2_0
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
      read(prm, nml=test_case)
      read(prm, nml=discretization)
      read(prm, nml=solver)
      select case(method)
      case(1)
        read(prm, nml=solver_cg)
      case(2)
        read(prm, nml=solver_schwarz)
      case(3,4)
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
  call XMPI_Bcast(nu_0   , 0, comm)
  call XMPI_Bcast(nu_1   , 0, comm)
  call XMPI_Bcast(d_nu   , 0, comm)
  call XMPI_Bcast(k_nu   , 0, comm)
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
  call XMPI_Bcast(method, 0, comm)

  ! CG/Schwarz options
  call XMPI_Bcast(i_max, 0, comm)
  call XMPI_Bcast(r_red, 0, comm)

  ! Schwarz and PMG options
  call schwarz_opt % Bcast(0, comm)
  call pmg_opt     % Bcast(0, comm)

  ! control parameters
  call XMPI_Bcast(plot_file, 0, comm)

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

  call problem % SetProblem(lambda, nu_0, nu_1, d_nu, k_nu, k_u)

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
  call mesh % GetPoints(po, 'GLL', x)

  ! mesh variables
  call InitializeMeshVariables()
  n = size(u)

  ! diffusivity
  call problem % GetDiffusivity(x, nu)

  ! auxiliary variables
  allocate(grad_u(0:po, 0:po, 0:po, mesh%ne, 3))

  ! operators ..................................................................

  ip_opt = IP_ElementOptions1D(po = po, penalty = penalty)
  elliptic_op = EllipticOperator3D_IP(mesh, lambda, nu, bc, ip_opt, schwarz_opt)

  !-----------------------------------------------------------------------------
  ! Tests

  ! exact solution and RHS  ....................................................

  if (n > 0) then

    ! exact solution and gradient
    call problem % GetExactSolution(x, u)
    call problem % GetExactGradient(x, grad_u)
    call SetArray(s, u)

    ! r = λ u - ∇·(ν ∇u)
    call problem % GetSource(x, r)

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
        bv(b) = BoundaryVariable(mesh, po, b, bc(b))
      end select
    end do

    call elliptic_op % BcToRHS(bv, 1, f)

  end if

  ! methods
  select case(method)
  case(3,4)
    pmg = PMG_Method3D(mesh, ip_opt, pmg_opt)
    call pmg % SetProblem(lambda, nu, bc)
  end select

  if (rank == 0) then
    dof = product(np) * product(ep) * (po + 1)**3
    write(*,'(/,A)') repeat('=',80)
    write(*,'(A,/)') 'IP/DG Elliptic Solver with variable diffusivity'
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
      write(*,'(A,/)') 'IP/DG EllipticOperator: p-Multigrid'
    case(4)
      write(*,'(A,/)') 'IP/DG EllipticOperator: p-MG/CG'
    end select
  end if

  !$omp parallel
  !$acc data copyin(f) copyout(u) create(r)

  !call SetArray(u, ZERO)
  call random_number(u)
  u = 2*u - 1

  call elliptic_op % Residual(u, f, r)
  r_l2_0 = ScalarProduct(r, r, mesh%comm)
  r_l2_0 = sqrt(r_l2_0)
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
    case(3:4)
      if (pmg_opt % monitor) write(*,*)
    end select
  end if

  !$acc end data
  !$omp end parallel

  if (rank == 0) then
    time0 = MPI_Wtime()
  end if

  select case(method)
  case(1) ! conjugate gradients
    call elliptic_op % ConjugateGradients(u, f, i_max, r_red, ni=ni)
  case(2) ! Schwarz method
    call elliptic_op % SchwarzMethod(u, f, i_max, r_red, ni=ni)
  case(3) ! p-MG method
    call pmg % MG_Solver(u, f, ni=ni)
  case(4) ! p-MG/CG method
    call pmg % MG_CG_Solver(u, f, ni=ni)
  end select

  if (rank == 0) then
    time = MPI_Wtime()
    time = (time - time0) / nt
  end if

  call elliptic_op % Residual(u, f, r)
  r_l2 = ScalarProduct(r, r, mesh%comm)
  r_l2 = sqrt(r_l2)

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

    !### CHECK
    block
      integer :: l, pl
      real(RNP), allocatable :: xl(:,:,:,:,:), vl(:,:,:,:,:)
      character(len=2) :: vl_names(4) = [ 'u ', 'f ', 'v ', 'nu' ]
      select case(method)
      case(3,4)
        do l = 0, ubound(pmg % level,1)
          pl = pmg % level(l) % po
          write(plot_file, '(A,I0)') 'elliptic_test__ip_vi_', l
          print '(2(A,I0))', 'level = ',l,': P_l = ', pl
          allocate(xl(0:pl,0:pl,0:pl,mesh%ne,4))
          allocate(vl(0:pl,0:pl,0:pl,mesh%ne,4))
          call mesh % GetPoints(pl, 'GLL', xl)
          vl(:,:,:,:,1) = pmg % level(l) % u
          vl(:,:,:,:,2) = pmg % level(l) % f
          vl(:,:,:,:,3) = pmg % level(l) % v
          if (allocated(pmg % level(l) % elliptic_op % nu_vi)) then
            vl(:,:,:,:,4) = pmg % level(l) % elliptic_op % nu_vi
          else
            vl(:,:,:,:,4) = 0
          end if
          call ExportVolumeDataToVTK( pl, mesh%ne, 4, 0, xl, vl, vl_names, &
                                      file   = trim(plot_file),            &
                                      part   = mesh%part,                  &
                                      n_part = mesh%n_part                 )
          deallocate(xl, vl)
        end do
      end select
    end block
    !### CHECK END

  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

!-------------------------------------------------------------------------------
!> Initialization of mesh variables

subroutine InitializeMeshVariables()

  integer :: ns = 6  ! number of scalar fields
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

  scalar_names(i) = 'nu'
  nu(0:po, 0:po, 0:po, 1:mesh%ne) => scalars(j:k)

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

end program Elliptic_Test__IP_VI

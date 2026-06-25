!> summary:  Testing elliptic multilevel solvers with adaptive refinement
!> author:   Joerg Stiller
!> date:     2025/01/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program ML_Elliptic_Test_Adaptive

  use, intrinsic :: ieee_arithmetic

  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Logging_Levels
  use Array_Assignments
  use Array_Reductions

  use TPO__Grad__3D
  use Elliptic_Problem__3D
  use Elliptic_Problem__Sphere__3D

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Mesh_Element__3D
  use Verify_Mesh__3D
  use Data_Exchange__3D
  use Export_VTK_Mesh_SFC__3D

  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  use ML__Boundary_Variable__3D
  use ML__Array_Reductions__3D
  use ML__DG__Elliptic_Solver__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP .............................................................

  type(MPI_Comm) :: comm      ! MPI communicator
  integer        :: rank      ! local MPI rank
  integer        :: n_proc    ! number of MPI processes
  integer        :: n_thread  ! number of OpenMP threads

  ! control .....................................................................

  ! NOTE
  ! The case name is defined by the first command line argument.
  ! If no argument is given, the default case is assumed.

  character(len=*), parameter :: default_case = 'ml_elliptic_test_adaptive'
  character(len=80) :: case_name ! case name
  character(len=80) :: case_file ! case input file: trim(case_name).prm

  logical :: use_solver = .true.  ! apply MG solver, else inject exact solution
  logical :: export_vtk = .false. ! switch for VTK export
  integer :: vtk_mode   =  3      ! 1/2/3: all/active/leaf elements

  namelist/control_prm/ use_solver, export_vtk, vtk_mode

  namelist/control_prm/ log_level
  namelist/control_prm/ log_level_inner_iteration
  namelist/control_prm/ log_level_outer_iteration
  namelist/control_prm/ log_level_multigrid_cycle

  ! domain .....................................................................

  integer :: test_domain = 1 ! computational domain
                             !   0  import from GMSH
                             !   1  cuboidal with Cartesian mesh
                             !   2  cuboidal with unstructured "diamond" mesh
                             !   3  cylindrical domain

  character(len=80) :: gmsh_file = '../../gmsh/cylinder_2d'

  namelist/domain_prm/ test_domain, gmsh_file

  ! mesh .......................................................................

  type(GenericMesh_3D)     , save :: generic_mesh
  type(Mesh_3D)            , save :: base_mesh
  type(ML_Mesh_Options_3D) , save :: ml_mesh_opt
  type(ML_Mesh_3D)         , save :: ml_mesh

  ! problem ....................................................................

  integer :: test_problem = 1 ! 1: sphere
  integer :: start_values = 0 ! 0/1/2: zero, exact, random

  namelist/problem_prm/ test_problem, start_values

  real(RNP) :: lambda =  0         ! Helmholtz parameter
  real(RNP) :: alpha  =  200       ! radial scaling factor
  real(RNP) :: x_c(3) = -0.05      ! front center
  real(RNP) :: r_f(2) = [0.7, 0.0] ! front radii
  integer   :: n_f

  namelist/problem_prm/ lambda, alpha, x_c, r_f

  ! boundary conditions: Dirichlet, if not periodic
  character, allocatable :: bc(:)

  class(EllipticProblem_3D), allocatable, save :: problem
  character(len=80) :: problem_name = ''

  ! solver parameters ..........................................................

  integer, save :: n_cycle = 1 ! number of solution cycles
  integer, allocatable, save :: po(:) ! sequence of polynomial orders

  type(ML_MeshOperators_3D),      save :: ml_op
  type(ML_DG_EllipticSolver_3D),  save :: ml_elliptic
  type(ML_DG_EllipticOptions_3D), save :: ml_elliptic_opt

  namelist/solver_prm/ n_cycle, po, ml_elliptic_opt

  ! adaptation parameters ......................................................

  integer   :: adapt_criterion = 4
  real(RNP) :: adapt_remove = 0.0
  real(RNP) :: adapt_refine = 0.7

  namelist/adaptation_prm/ adapt_criterion, adapt_remove, adapt_refine

  type(DataExchangePlan_3D), allocatable, save :: x_plan(:)

  ! variables ..................................................................

  character(len=:), allocatable, save :: name_var(:) ! names of variables
  type(ML_MeshVariable_3D), save :: var ! container of variables

  type(ML_MeshVariable_3D), save :: mm     ! diagonal mass matrix
  type(ML_MeshVariable_3D), save :: f      ! sources
  type(ML_MeshVariable_3D), save :: u      ! numerical solution
  type(ML_MeshVariable_3D), save :: s      ! exact solution
  type(ML_MeshVariable_3D), save :: e      ! H0/H1 errors per element
  type(ML_MeshVariable_3D), save :: r      ! residual or just workspace
  type(ML_MeshVariable_3D), save :: qi     ! quantity of interest
  type(ML_MeshVariable_3D), save :: grad_u ! gradient of numerical solution
  type(ML_MeshVariable_3D), save :: grad_s ! gradient of exact solution

  type(ML_BoundaryVariable_3D), save :: bv ! boundary values

  ! auxiliaries ................................................................

  character(:), allocatable, save :: domain_name
  real(RNP)   , allocatable, save :: t_adapt(:), t_setup(:), t_solve(:)

  character(80), save :: message
  real(RNP)    , save :: t_start
  logical      , save :: exists, passed, all_passed
  integer      , save :: l_top, l_max, l_adapt
  integer      , save :: io, stat

  integer :: i, l, m

  !-----------------------------------------------------------------------------
  ! Initialization

  ! MPI and OpenMP .............................................................

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  !$omp parallel
  n_thread = OMP_Num_Threads()
  !$omp end parallel

  ! control and domain parameters ..............................................

  if (rank == 0) then

    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Validation of elliptic multilevel solvers'
    write(*,*)
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:', n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'  , n_thread
    write(*,*)

    call get_command_argument(1, case_name, status=stat)
    if (stat /= 0 .or. len_trim(case_name) == 0) then
      case_name = default_case
    end if
    case_file = trim(case_name) // '.prm'

    inquire(file=case_file, exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(case_file)
      open(newunit = io, file = case_file)
      read(io, nml = control_prm)
      read(io, nml = domain_prm)
      ml_mesh_opt = ML_Mesh_Options_3D(io, n_proc)
      close(io)
    else
      call Error('ML_Elliptic_Test_Adaptive', &
                 'file "'// trim(case_file) //'" not found')
    end if

  end if

  ! globalize multilevel mesh options
  call ml_mesh_opt % Bcast(0, comm)

  ! globalize logging levels
  call XMPI_Bcast_LoggingLevels(0, comm)

  ! globalize remaining parameters
  call XMPI_Bcast(case_name   , 0, comm)
  call XMPI_Bcast(use_solver  , 0, comm)
  call XMPI_Bcast(export_vtk  , 0, comm)
  call XMPI_Bcast(vtk_mode    , 0, comm)
  call XMPI_Bcast(test_domain , 0, comm)
  call XMPI_Bcast(gmsh_file   , 0, comm)

  ! base mesh ..................................................................

  if (rank == 0) then
    write(*,'(/,A)') 'Base mesh'
  end if

  select case(test_domain)
  case(0)
    if (rank == 0) then
      call ImportGMSH_3D(gmsh_file, generic_mesh)
    end if
    call base_mesh % ImportGenericMesh(generic_mesh, comm)
    domain_name = gmsh_file
  case(1)
    call CreateCuboidCartesian(comm, case_file, base_mesh)
    domain_name = 'Cuboidal domain with Cartesian mesh'
  case(2)
    call CreateCuboidDiamonds(comm, case_file, base_mesh)
    domain_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCylinder(comm, case_file, base_mesh)
    domain_name = 'Cylindrical domain with unstructured mesh'
  case(4)
    call CreateAnnulus(comm, case_file, base_mesh)
    domain_name = 'Annular domain with unstructured mesh'
  end select

  if (rank == 0) then
    write(*,'(2X,2A)') 'created mesh for ', trim(domain_name)
  end if

  call VerifyMesh_3D(base_mesh, passed)
  call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, comm)

  if (rank == 0) then
    write(*,'(2X,A,G0)') 'verification: passed = ', all_passed
  end if

  call MPI_Barrier(comm)

  ! initial multilevel mesh ....................................................

  if (rank == 0) then
    write(*,'(/,A)') 'Initial multilevel mesh'
  end if

  ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)

  l_top   = size(ml_mesh % mesh)
  l_max   = ml_mesh_opt % l_max
  l_adapt = ml_mesh_opt % l_adapt

  associate(mesh => ml_mesh%mesh)
    ! verification
    do l = 1, l_top
      call VerifyMesh_3D(mesh(l), passed)
      call XMPI_Allreduce(passed, all_passed, MPI_LAND, comm)
      if (.not. all_passed) exit
    end do
    call MPI_Barrier(comm)
    if (rank == 0) then
      if (all_passed) then
        write(*,'(2X,A,/)') 'verification: all levels passed'
      else
        write(*,'(2X,A,I0,A,/)') 'verification of level ',l,' failed'
      end if
    end if
  end associate

  ! problem ....................................................................

  allocate(bc(base_mesh % n_bound), source = 'D')

  if (rank == 0) then
    write(*,'(/,A)') 'Elliptic problem'
    open(newunit = io, file = case_file)
    read(io, nml = problem_prm)
    close(io)
  end if

  ! globalize problem parameters
  call XMPI_Bcast( test_problem, 0, comm )
  call XMPI_Bcast( start_values, 0, comm )
  call XMPI_Bcast( lambda      , 0, comm )
  call XMPI_Bcast( alpha       , 0, comm )
  call XMPI_Bcast( x_c         , 0, comm )
  call XMPI_Bcast( r_f         , 0, comm )
  call XMPI_Bcast( bc          , 0, comm )
  n_f = count(r_f > 0)

  select case(test_problem)
  case default
    problem_name = 'Sphere'
    problem = EllipticProblem_Sphere_3D(lambda, alpha, x_c, r_f(1:n_f))
  end select

  ! enforce periodicity at coupled boundaries
  where(base_mesh % boundary % coupled > 0) bc = 'P'

  if (rank == 0) then
    write(*,'(T3,A,T26,9(G0,X))') 'problem name:'        , trim(problem_name)
    write(*,'(T3,A,T26,9(G0,X))') 'boundary conditions:' , bc
  end if

  ! solver and adaptation parameters ...........................................

  allocate(po(l_max), source = 1)

  if (rank == 0) then
    open(newunit = io, file = case_file)
    read(io, nml = solver_prm)
    read(io, nml = adaptation_prm)
    close(io)
  end if

  call XMPI_Bcast(n_cycle, 0, comm)
  call XMPI_Bcast(po     , 0, comm)
  call ml_elliptic_opt % Bcast(0, comm)

  call XMPI_Bcast(adapt_criterion, 0, comm)
  call XMPI_Bcast(adapt_remove   , 0, comm)
  call XMPI_Bcast(adapt_refine   , 0, comm)

  allocate(t_adapt(n_cycle), source = ZERO)
  allocate(t_setup, t_solve, source = t_adapt)

  !-----------------------------------------------------------------------------
  ! first iteration

  ! preliminaries ..............................................................

  if (rank == 0) then
    t_start = MPI_Wtime()
  end if

  ml_op = ML_MeshOperators_3D(ml_mesh, po)
  ml_elliptic = ML_DG_EllipticSolver_3D(ml_op, ml_elliptic_opt)

  call mm % Init(ml_op, 1)
  call f  % Init(ml_op, 1)
  call u  % Init(ml_op, 1)
  call s  % Init(ml_op, 1)
  call e  % Init(ml_op, 1)
  call r  % Init(ml_op, 1)
  call bv % Init(ml_op, 1)

  call qi % Init(ml_mesh, po = 0, nc = 1)

  call grad_u % Init(ml_op, 3)
  call grad_s % Init(ml_op, 3)


  ! RHS, BC and start values ...................................................

  do l = 1, l_top
    associate( x_l  => ml_op % sem(l) % metrics % x   &
             , mm_l => mm % level(l) % val(:,:,:,:,1) &
             , f_l  => f  % level(l) % val(:,:,:,:,1) &
             , s_l  => s  % level(l) % val(:,:,:,:,1) &
             , u_l  => u  % level(l) % val(:,:,:,:,1) )

      call ml_op % sem(l) % Get_DG_DiagonalMassMatrix(mm_l)

      call problem % GetExactSolution(x_l, s_l)
      call problem % GetExactGradient(x_l, grad_s % level(l) % val)
      call problem % GetSource(x_l, f_l)
      f_l = mm_l * f_l

      ! boundary conditions
      do i = 1, size(bc)
        select case(bc(i))
        case('D')
          call bv % level(l) % var(i) % Extract(s_l)
        end select
      end do

      ! start values
      select case(start_values)
      case(0)
        call SetArray(u_l, ZERO)
      case(1)
        call SetArray(u_l, s_l)
      case default
        call random_number(u_l)
        u_l = 2*u_l - 1
      end select

    end associate
  end do

  if (rank == 0) then
    t_setup(1) = MPI_Wtime() - t_start
  end if

  ! solution ...................................................................

  if (rank == 0) then
    write(*,'(/,A,I0)') 'Cycle ', 1
    t_start = MPI_Wtime()
  end if

  call ml_elliptic % FAS_MG_Solver(bc, lambda, problem%nu_0, bv, f, u)

  if (rank == 0) then
    t_solve(1) = MPI_Wtime() - t_start
  end if

  call Evaluation

  !------------------------------------------------------------------------------
  ! adaptation/solution cycles

  do m = 2, n_cycle

    if (rank == 0) then
       write(*,'(/,A,I0)') 'Cycle ', m
    end if

    if (l_top < l_max) then
      ml_mesh % mesh(l_top) % refinement = ml_mesh_opt % refinement(l_top)
    end if
    call SetAdaptationMarks

    if (rank == 0) then
      t_start = MPI_Wtime()
    end if

    call ml_mesh % Adapt(ml_mesh_opt % partition, x_plan)

    if (rank == 0) then
      t_adapt(m) = MPI_Wtime() - t_start
    end if

    ! verification
    associate(mesh => ml_mesh%mesh)
      do l = 1, size(mesh)
        call VerifyMesh_3D(mesh(l), passed)
        call XMPI_Allreduce(passed, all_passed, MPI_LAND, comm)
        if (all_passed) cycle
        write(message,'(A,I0,A)') 'verification of level ',l,' failed'
        call Error('ML_Elliptic_Test_Adaptive', message)
      end do
    end associate

    l_top = size(ml_mesh % mesh)
    ml_op = ML_MeshOperators_3D(ml_mesh, po)

    if (rank == 0) then
      t_start = MPI_Wtime()
    end if

    ! rebuild elliptic operators -- can be optimized
    ml_elliptic = ML_DG_EllipticSolver_3D(ml_op, ml_elliptic_opt)

    if (rank == 0) then
      t_setup(m) = MPI_Wtime() - t_start
      t_start = MPI_Wtime()
    end if

    ! interpolate/redistribute solution
    call u % FitAdapt(ml_op, x_plan)

    if (rank == 0) then
      t_adapt(m) = t_adapt(m) + MPI_Wtime() - t_start
    end if

    do l = 1, l_top
      associate(u_l => u  % level(l) % val(:,:,:,:,1))
        if (any(ieee_is_nan(u_l))) then
          print '(99G0)', '#1 [',ml_mesh%mesh(l)%proc,'] u_',l,' has NaN'
        end if
      end associate
    end do

    if (rank == 0) then
      t_start = MPI_Wtime()
    end if

    ! adjust remaining variables
    call mm % Init(ml_op, 1)
    call f  % Init(ml_op, 1)
    call s  % Init(ml_op, 1)
    call e  % Init(ml_op, 1)
    call r  % Init(ml_op, 1)
    call bv % Init(ml_op, 1)

    call qi % Init(ml_mesh, po = 0, nc = 1)

    call grad_u % Init(ml_op, 3)
    call grad_s % Init(ml_op, 3)

    ! RHS and BC
    do l = 1, l_top
      associate( x_l  => ml_op % sem(l) % metrics % x   &
               , mm_l => mm % level(l) % val(:,:,:,:,1) &
               , f_l  => f  % level(l) % val(:,:,:,:,1) &
               , s_l  => s  % level(l) % val(:,:,:,:,1) )

        ! right hand side
        call ml_op % sem(l) % Get_DG_DiagonalMassMatrix(mm_l)
        call problem % GetExactSolution(x_l, s_l)
        call problem % GetExactGradient(x_l, grad_s % level(l) % val)
        call problem % GetSource(x_l, f_l)
        f_l = mm_l * f_l

        ! boundary conditions
        do i = 1, size(bc)
          select case(bc(i))
          case('D')
            call bv % level(l) % var(i) % Extract(s_l)
          end select
        end do

      end associate
    end do

    if (rank == 0) then
      t_setup(m) = t_setup(m) + MPI_Wtime() - t_start
      t_start = MPI_Wtime()
    end if

    if (use_solver) then
      call ml_elliptic % FAS_MG_Solver(bc, lambda, problem%nu_0, bv, f, u)
    else
      ! inject exact solution
      do l = 1, l_top
        call problem % GetExactSolution( x = ml_op % sem(l) % metrics % x  &
                                       , u = u % level(l) % val(:,:,:,:,1) )
      end do
    end if

    if (rank == 0) then
      t_solve(m) = MPI_Wtime() - t_start
    end if

    do l = 1, l_top
      associate(u_l => u  % level(l) % val(:,:,:,:,1))
        if (any(ieee_is_nan(u_l))) then
          print '(99G0)', '#2 [',ml_mesh%mesh(l)%proc,'] u_',l,' has NaN'
        end if
      end associate
    end do

    call Evaluation

  end do

  ! print times
  if (rank == 0) then
    write(*,'(/,2X,A5,3(4X,A7,X))') 'cycle', 't_adapt', 't_setup', 't_solve'
    do m = 1, n_cycle
      write(*,'(I5,2X,3(2X,ES10.3))') m, t_adapt(m), t_setup(m), t_solve(m)
    end do
    write(*,'(A5,2X,3(2X,ES10.3))') &
        'sum', sum(t_adapt), sum(t_setup), sum(t_solve)
  end if

  !-----------------------------------------------------------------------------
  ! VTK export

  name_var = [ 'f ', 'u ', 's ', 'r ' ]
  call var % Init(ml_op, size(name_var), name_var)

  do l = 1, l_top
    call SetArray(var%level(l)%val(:,:,:,:,1), f%level(l)%val(:,:,:,:,1))
    call SetArray(var%level(l)%val(:,:,:,:,2), u%level(l)%val(:,:,:,:,1))
    call SetArray(var%level(l)%val(:,:,:,:,3), s%level(l)%val(:,:,:,:,1))
    call SetArray(var%level(l)%val(:,:,:,:,4), r%level(l)%val(:,:,:,:,1))
  end do

  !-----------------------------------------------------------------------------
  ! VTK export

  if (export_vtk) then
    vtk_mode = max(1, min(3, vtk_mode))
    call var % ExportVTK(ml_op, trim(case_name), vtk_mode)
    if (vtk_mode < 3) then
      associate(mesh => ml_mesh%mesh)
        block
          character(len=9) :: tag
          do l = 1, size(mesh)
            if (mesh(l) % has_sfc) then
              write(tag,'(A,I0,A)') '_sfc_l', l
              call ExportVTK_MeshSFC(mesh(l), file = trim(case_name)//trim(tag))
            end if
          end do
        end block
      end associate
    end if
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !> Evaluation

  subroutine Evaluation

    real(RNP), save :: e_0, e_1, r_2

    ! mesh metrics .............................................................

    call ml_op % Print_MeshCharacteristics()

    ! E2 composite leaf residual ...............................................

    call ml_elliptic % FAS_MG_Residual(bc, lambda, problem%nu_0, bv, f, u, r)
    r_2 = sqrt(ML_ScalarProduct_3D(r, r, leaf = .true.))

    ! H0 composite error ε₀ ....................................................

    do l = 1, l_top
      associate( e_l => e % level(l) % val(:,:,:,:,1) &
               , s_l => s % level(l) % val(:,:,:,:,1) &
               , u_l => u % level(l) % val(:,:,:,:,1) )

        e_l = u_l - s_l

      end associate
    end do

    ! composite H0 norm
    e_0 = sqrt(ML_WeightedScalarProduct_3D(mm, e, e, leaf = .true.))

    ! quantity of interest based on elemental H0 norm
    if (1 <= adapt_criterion .and. adapt_criterion <= 2) then
      do l = 1, l_top
        associate( mesh_l => ml_mesh % mesh(l)              &
                 , mm_l   => mm % level(l) % val(:,:,:,:,1) &
                 , qi_l   => qi % level(l) % val(0,0,0,:,1) &
                 , e_l    => e  % level(l) % val(:,:,:,:,1) )

          do i = 1, mesh_l % n_elem_active
            select case(adapt_criterion)
            case(1)
              qi_l(i) = sqrt( sum(mm_l(:,:,:,i) * e_l(:,:,:,i)**2) )
            case(2)
              qi_l(i) = sqrt( sum(mm_l(:,:,:,i) * e_l(:,:,:,i)**2) &
                            / sum(mm_l(:,:,:,i)) )
            end select
          end do
        end associate
      end do
    end if

    ! H1 error ε₁ ..............................................................

    do l = 1, l_top
      associate( ell_l => ml_elliptic % elliptic_op(l)  &
               , e_l   => e % level(l) % val(:,:,:,:,1) &
               , r_l   => r % level(l) % val(:,:,:,:,1) )

        ! r = Aε
        call ml_elliptic % elliptic_op(l) &
                 % Apply(bc, lambda, problem%nu_0, u = e_l, r = r_l)

      end associate
    end do

    ! ε₁ = √ εAε
    e_1 = sqrt(ML_ScalarProduct_3D(e, r, leaf = .true.))

    ! quantity of interest based on elemental H1 seminorm
    if (3 <= adapt_criterion .and. adapt_criterion <= 4) then
      do l = 1, l_top
        associate( mesh_l => ml_mesh % mesh(l)              &
                 , mm_l   => mm % level(l) % val(:,:,:,:,1) &
                 , qi_l   => qi % level(l) % val(0,0,0,:,1) &
                 , e_l    => e  % level(l) % val(:,:,:,:,1) &
                 , r_l    => r  % level(l) % val(:,:,:,:,1) )

          do i = 1, mesh_l % n_elem_active
            select case(adapt_criterion)
            case(3)
              qi_l(i) = sqrt( max(sum(e_l(:,:,:,i) * r_l(:,:,:,i)), ZERO) )
            case(4)
              qi_l(i) = sqrt( max(sum(e_l(:,:,:,i) * r_l(:,:,:,i)), ZERO) &
                            / sum(mm_l(:,:,:,i)) )
            end select
          end do
        end associate
      end do
    end if

    ! print error metrics ......................................................

    if (rank == 0) then
      write(*,*)
      write(*,'(2X,A)') 'error metrics'
      write(*,'(T5,A,T10,ES12.5)') 'e_0 =', e_0
      write(*,'(T5,A,T10,ES12.5)') 'e_1 =', e_1
      write(*,'(T5,A,T10,ES12.5)') 'r_2 =', r_2
    end if

  end subroutine Evaluation

  !-----------------------------------------------------------------------------
  !> Set the adaptation marks

  subroutine SetAdaptationMarks

    integer, allocatable, save :: n_refine_loc(:), n_refine(:)
    integer, allocatable, save :: n_remove_loc(:), n_remove(:)

    real(RNP) :: max_e = 0, max_e_loc = 0
    real(RNP) :: qi_refine, qi_remove

    associate(mesh => ml_mesh%mesh)

      allocate(n_refine_loc(l_top), n_refine(l_top), source = 0)
      allocate(n_remove_loc(l_top), n_remove(l_top), source = 0)

      ! mark elements in globally refined levels ...............................

      do l = 1, min(l_top,l_adapt-2)
        if (mesh(l)%n_elem > 0) then
          call mesh(l) % element % MarkForRefinement()
        end if
      end do

      ! criteria for refinement/removal ........................................

      max_e_loc = 0

      ! maximum error over leaf elements
      do l = l_adapt-1, min(l_top,l_max-1)
        associate(qi_l => qi % level(l) % val(0,0,0,:,1))
          do i = 1, mesh(l) % n_elem_active
            if (mesh(l) % element(i) % IsLeaf()) then
              max_e_loc = max(max_e_loc, qi_l(i))
            end if
          end do
        end associate
      end do

      call XMPI_Allreduce(max_e_loc, max_e, MPI_MAX, comm)

      qi_refine = adapt_refine * max_e
      qi_remove = adapt_remove * max_e

      ! mark for adaptation ....................................................

      ! l < l_max: mark elements for removal or refinement
      do l = l_adapt-1, min(l_top,l_max-1)
      do i = 1, mesh(l) % n_elem
        associate( element => mesh(l) % element(i)        &
                 , qi_e => qi % level(l) % val(0,0,0,i,1) )

          ! elements that are not leaves are marked for removal
          if (.not. element % IsLeaf()) then
            call element % MarkForRemoval()

          ! unmark elements close to front center
          else if (AtCenter(element, x_c)) then
            call element % Unmark()

          ! refine elements cut by the front or above refinement threshold
          else if (OnFront(element, x_c, r_f) .or. qi_e > qi_refine) then
            call element % MarkForRefinement()
            n_refine_loc(l) = n_refine_loc(l) + 1

          ! remove elements below removal threshold
          else if (qi_e < adapt_remove) then
            call element % MarkForRemoval()
            n_remove_loc(l) = n_remove_loc(l) + 1

          else
            call element % Unmark()
          end if

        end associate
      end do
      end do

      ! l = l_max: mark elements for removal
      if (l_top == l_max) then
        do i = 1, mesh(l_top) % n_elem
          associate( element => mesh(l_top) % element(i)        &
                   , qi_e => qi % level(l_top) % val(0,0,0,i,1) )

            if (qi_e < adapt_remove) then
              call element % MarkForRemoval()
              n_remove_loc(l_top) = n_remove_loc(l_top) + 1
            else
              call element % Unmark()
            end if

          end associate
        end do
      end if

      ! number of leaf elements marked for refinement
      call XMPI_Reduce(n_refine_loc, n_refine, MPI_SUM, 0, comm)
      call XMPI_Reduce(n_remove_loc, n_remove, MPI_SUM, 0, comm)

      ! control output .........................................................

      if (rank == 0) then
        write(*,'(/,2X,A,99(X,G0))') 'elements marked for refinement: ',n_refine
        write(*,'(  2X,A,99(X,G0))') 'elements marked for removal:    ',n_remove
      end if

      deallocate(n_refine_loc, n_refine)
      deallocate(n_remove_loc, n_remove)

    end associate

  end subroutine SetAdaptationMarks

  !-----------------------------------------------------------------------------
  !> Check if element is located close to front center

  pure logical function AtCenter(element, x_c)
    class(MeshElement_3D), intent(in) :: element
    real(RNP), intent(in) :: x_c(3) !< front center

    real(RNP), parameter :: tol = 0.6
    real(RNP) :: delta, dx_e(3)

    ! displacement between front center and element center
    delta = sqrt(sum( (element % geometry % x_c(0,1:3) - x_c)**2 ))

    ! element dimensions
    call element % GetCuboidDimensions(dx_e)

    AtCenter = delta < tol * maxval(dx_e)

  end function AtCenter

  !-----------------------------------------------------------------------------
  !> Check if element is cut by front

  pure logical function OnFront(element, x_c, r_f)
    class(MeshElement_3D), intent(in) :: element
    real(RNP), intent(in) :: x_c(3) !< front center
    real(RNP), intent(in) :: r_f(:) !< front radii

    real(RNP), parameter :: tol = 0.7
    real(RNP) :: dx_c(3), dx_e(3), r_e

    ! radius at element center
    dx_c = element % geometry % x_c(0,1:3) - x_c
    r_e  = sqrt(sum( dx_c**2 ))

    ! element dimensions
    call element % GetCuboidDimensions(dx_e)

    OnFront = any(r_f > 0 .and. abs(r_f - r_e) < tol * maxval(dx_e))

  end function OnFront

  !=============================================================================

end program ML_Elliptic_Test_Adaptive

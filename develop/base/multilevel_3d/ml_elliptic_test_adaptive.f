!> summary:  Testing elliptic multilevel solvers with adaptive refinement
!> author:   Joerg Stiller
!> date:     2025/01/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program ML_Elliptic_Test_Adaptive
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Logging_Levels
  use Array_Assignments
  use Array_Reductions

  use Elliptic_Problem__3D
  use Elliptic_Problem__Sphere__3D

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Verify_Mesh__3D
  use Data_Exchange__3D

  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  use ML__Boundary_Variable__3D
  use ML__Array_Reductions__3D
  use ML__DG__Elliptic_Solver__3D

!### CHECK
use, intrinsic :: ieee_arithmetic
!### CHECK END
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

  logical :: export_vtk = .false.  ! switch for VTK export

  namelist/control_prm/ export_vtk

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

  character(len=80) :: gmsh_file = '../gmsh_3d/cylinder_2d'

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

  real(RNP) :: lambda =  0    ! Helmholtz parameter
  real(RNP) :: x_c(3) = -0.05 ! sphere center
  real(RNP) :: r_0    =  0.7  ! sphere radius
  real(RNP) :: alpha  =  200  ! radial scaling factor

  namelist/problem_prm/ lambda, x_c, r_0, alpha

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

  integer   :: adapt_criterion = 1
  real(RNP) :: adapt_remove = 0.0
  real(RNP) :: adapt_refine = 0.7

  namelist/adaptation_prm/ adapt_criterion, adapt_remove, adapt_refine

  type(DataExchangePlan_3D), allocatable, save :: x_plan(:)

  ! variables ..................................................................

  character(len=:), allocatable, save :: name_var(:) ! names of variables
  type(ML_MeshVariable_3D), save :: var ! container of variables

  type(ML_MeshVariable_3D), save :: mm  ! diagonal mass matrix
  type(ML_MeshVariable_3D), save :: f   ! sources
  type(ML_MeshVariable_3D), save :: u   ! numerical solution
  type(ML_MeshVariable_3D), save :: s   ! exact solution
  type(ML_MeshVariable_3D), save :: e   ! error
  type(ML_MeshVariable_3D), save :: r   ! residual or just workspace

  type(ML_BoundaryVariable_3D), save :: bv ! boundary values

  ! auxiliaries ................................................................

  character(:), allocatable, save :: domain_name

  character(len=100), save :: message
  logical, save :: exists, passed, all_passed
  integer, save :: l_top, l_max, l_adapt
  integer, save :: io, stat

  logical :: singular
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
  call XMPI_Bcast(case_name       , 0, comm)
  call XMPI_Bcast(export_vtk      , 0, comm)
  call XMPI_Bcast(test_domain     , 0, comm)
  call XMPI_Bcast(gmsh_file       , 0, comm)

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

  l_top   = ml_mesh_opt % l_top
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
  call XMPI_Bcast( test_problem   , 0, comm )
  call XMPI_Bcast( start_values   , 0, comm )
  call XMPI_Bcast( lambda         , 0, comm )
  call XMPI_Bcast( x_c            , 0, comm )
  call XMPI_Bcast( r_0            , 0, comm )
  call XMPI_Bcast( alpha          , 0, comm )
  call XMPI_Bcast( bc             , 0, comm )

  select case(test_problem)
  case default
    problem_name = 'Sphere'
    problem = EllipticProblem_Sphere_3D(lambda, x_c, r_0, alpha)
  end select

  ! enforce periodicity at coupled boundaries
  where(base_mesh % boundary % coupled > 0) bc = 'P'

  singular = lambda == ZERO .and. all(bc == 'P' .or. bc == 'N')

  if (rank == 0) then
    write(*,'(T3,A,T26,9(G0,X))') 'problem name:'        , trim(problem_name)
    write(*,'(T3,A,T26,9(G0,X))') 'boundary conditions:' , bc
  end if

  ! solver and adaptation parameters ...........................................

  allocate(po(l_max), source = 1)

  if (rank == 0) then
    open(newunit = io, file = case_file)
    read(io, nml = solver_prm, iostat = stat)
    read(io, nml = adaptation_prm)
    close(io)
  end if

  call XMPI_Bcast(n_cycle, 0, comm)
  call XMPI_Bcast(po     , 0, comm)
  call ml_elliptic_opt % Bcast(0, comm)

  call XMPI_Bcast(adapt_criterion, 0, comm)
  call XMPI_Bcast(adapt_remove   , 0, comm)
  call XMPI_Bcast(adapt_refine   , 0, comm)

  !-----------------------------------------------------------------------------
  ! first iteration

  ! preliminaries ..............................................................

  ml_op = ML_MeshOperators_3D(ml_mesh, po)
  ml_elliptic = ML_DG_EllipticSolver_3D(ml_op, ml_elliptic_opt)

  call mm % Init(ml_op, 1)
  call f  % Init(ml_op, 1)
  call u  % Init(ml_op, 1)
  call s  % Init(ml_op, 1)
  call e  % Init(ml_op, 1)
  call r  % Init(ml_op, 1)
  call bv % Init(ml_op, 1)

  ! RHS, BC and start values ...................................................

  do l = 1, l_top
    associate( x_l   => ml_op % sem(l) % metrics % x   &
             , mm_l  => mm % level(l) % val(:,:,:,:,1) &
             , f_l   => f  % level(l) % val(:,:,:,:,1) &
             , s_l   => s  % level(l) % val(:,:,:,:,1) &
             , u_l   => u  % level(l) % val(:,:,:,:,1) )

      call ml_op % sem(l) % Get_DG_DiagonalMassMatrix(mm_l)

      call problem % GetExactSolution(x_l, s_l)
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

  ! solution ...................................................................

  if (rank == 0) then
     write(*,'(/,A,I0)') 'Cycle ', 1
  end if

  call ml_elliptic % FAS_MG_Solver(bc, lambda, problem%nu_0, u, f, bv)

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

    call ml_mesh % Adapt(ml_mesh_opt % partition, x_plan)

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

    ! rebuild elliptic operators -- can be optimized
    ml_elliptic = ML_DG_EllipticSolver_3D(ml_op, ml_elliptic_opt)

    ! interpolate/redistribute solution
    call u % FitAdapt(ml_op, x_plan)
!### CHECK
do l = 1, l_top
  associate(u_l => u  % level(l) % val(:,:,:,:,1))
    if (any(ieee_is_nan(u_l))) then
      print '(99G0)', '#1 [',ml_mesh%mesh(l)%proc,'] u_',l,' has NaN'
    end if
  end associate
end do
!### CHECK END

    ! adjust remaining variables
    call mm % Init(ml_op, 1)
    call f  % Init(ml_op, 1)
    call s  % Init(ml_op, 1)
    call e  % Init(ml_op, 1)
    call r  % Init(ml_op, 1)
    call bv % Init(ml_op, 1)

    ! RHS and BC
    do l = 1, l_top
      associate( x_l  => ml_op % sem(l) % metrics % x   &
               , mm_l => mm % level(l) % val(:,:,:,:,1) &
               , f_l  => f  % level(l) % val(:,:,:,:,1) &
               , s_l  => s  % level(l) % val(:,:,:,:,1) )

        ! right hand side
        call ml_op % sem(l) % Get_DG_DiagonalMassMatrix(mm_l)
        call problem % GetExactSolution(x_l, s_l)
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

    call ml_elliptic % FAS_MG_Solver(bc, lambda, problem%nu_0, u, f, bv)
!### CHECK
do l = 1, l_top
  associate(u_l => u  % level(l) % val(:,:,:,:,1))
    if (any(ieee_is_nan(u_l))) then
      print '(99G0)', '#2 [',ml_mesh%mesh(l)%proc,'] u_',l,' has NaN'
    end if
  end associate
end do
!### CHECK END

    call Evaluation

  end do

  !-----------------------------------------------------------------------------
  ! VTK export

  name_var = [ 'f ', 'u ', 's ', 'e ', 'r ' ]
  call var % Init(ml_op, size(name_var), name_var)

  do l = 1, l_top
    call SetArray(var%level(l)%val(:,:,:,:,1), f%level(l)%val(:,:,:,:,1))
    call SetArray(var%level(l)%val(:,:,:,:,2), u%level(l)%val(:,:,:,:,1))
    call SetArray(var%level(l)%val(:,:,:,:,3), s%level(l)%val(:,:,:,:,1))
    call SetArray(var%level(l)%val(:,:,:,:,4), e%level(l)%val(:,:,:,:,1))
    call SetArray(var%level(l)%val(:,:,:,:,5), r%level(l)%val(:,:,:,:,1))
  end do

  if (export_vtk) then
    call var % ExportVTK(ml_op, trim(case_name)//'_full', mode=1)
    call var % ExportVTK(ml_op, trim(case_name)//'_leaf', mode=3)
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !> Evaluation

  subroutine Evaluation

    real(RNP), parameter :: eps = epsilon(ONE) / 1000

    real(RNP), save :: int_1, int_1_loc
    real(RNP), save :: int_e, int_e_loc
    real(RNP), save :: e_l2, r_l2

    integer, save :: ne_max, ne_min, ne_tot, ne_leaf, ne_leaf_loc
    integer, save :: ne_tot_sum, ne_leaf_sum
    integer, save :: np_tot_sum, np_leaf_sum

    real(RNP) :: e_avg
    integer   :: i, l

    associate(mesh => ml_mesh%mesh)

      ! residual .................................................................

      call ml_elliptic % FAS_MG_Residual(bc, lambda, problem%nu_0, f, bv, u, r)
      r_l2 = sqrt(ML_WeightedScalarProduct_3D(mm, r, r, leaf = .true.))

      ! error ....................................................................

      int_1_loc = 0
      int_e_loc = 0

      do l = 1, l_top
        associate( mm_l   => mm % level(l) % val(:,:,:,:,1) &
                 , s_l    => s  % level(l) % val(:,:,:,:,1) &
                 , u_l    => u  % level(l) % val(:,:,:,:,1) &
                 , e_l    => e  % level(l) % val(:,:,:,:,1) )

          do i = 1, mesh(l) % n_elem
            e_l(:,:,:,i) = u_l(:,:,:,i) - s_l(:,:,:,i)
            if (mesh(l)%element(i)%IsLeaf()) then
              int_1_loc = int_1_loc + sum(mm_l(:,:,:,i))
              int_e_loc = int_e_loc + sum(mm_l(:,:,:,i) * e_l(:,:,:,i))
            end if
          end do

        end associate
      end do

      ! mean error for calibration in singular case
      if (singular) then
        call XMPI_Allreduce(int_1_loc, int_1, MPI_SUM, comm)
        call XMPI_Allreduce(int_e_loc, int_e, MPI_SUM, comm)
        e_avg = int_e / int_1
      else
        e_avg = 0
      end if

      ! calibration and maximum norm
      do l = 1, l_top
        associate(e_l => e % level(l) % val(:,:,:,:,1))
          do i = 1, mesh(l) % n_elem
            e_l(:,:,:,i) = e_l(:,:,:,i) - e_avg
          end do
        end associate
      end do

      ! L2 norm
      e_l2 = sqrt(ML_WeightedScalarProduct_3D(mm, e, e, leaf = .true.))

      ! mesh metrics .............................................................

      if (rank == 0) then
        write(*,'(/,A7,5A9,A5)') '  level'   &
                               , '  n_parts' &
                               , '   ne_min' &
                               , '   ne_max' &
                               , '   ne_tot' &
                               , '  ne_leaf' &
                               , '  po'
      end if

      ne_tot_sum  = 0
      np_tot_sum  = 0
      ne_leaf_sum = 0
      np_leaf_sum = 0

      do l = 1, l_top

        ne_leaf_loc = 0

        if (mesh(l)%part >= 0) then
          do i = 1, mesh(l) % n_elem_active
            if (mesh(l) % element(i) % IsLeaf()) then
              ne_leaf_loc = ne_leaf_loc + 1
            end if
          end do
          call XMPI_Reduce(mesh(l)%n_elem, ne_min , MPI_MIN, 0, mesh(l)%comm_parts)
          call XMPI_Reduce(mesh(l)%n_elem, ne_max , MPI_MAX, 0, mesh(l)%comm_parts)
          call XMPI_Reduce(mesh(l)%n_elem, ne_tot , MPI_SUM, 0, mesh(l)%comm_parts)
          call XMPI_Reduce(ne_leaf_loc   , ne_leaf, MPI_SUM, 0, mesh(l)%comm_parts)
        end if

        call XMPI_Bcast(ne_min , mesh(l)%proc_part(0), mesh(l)%comm_world)
        call XMPI_Bcast(ne_max , mesh(l)%proc_part(0), mesh(l)%comm_world)
        call XMPI_Bcast(ne_tot , mesh(l)%proc_part(0), mesh(l)%comm_world)
        call XMPI_Bcast(ne_leaf, mesh(l)%proc_part(0), mesh(l)%comm_world)

        if (rank == 0) then
          write(*,'(I7,5I9,I5)') l, mesh(l)%n_parts,             &
                                 ne_min, ne_max, ne_tot,ne_leaf, &
                                 po(l)

          ne_tot_sum  = ne_tot_sum  + ne_tot
          np_tot_sum  = np_tot_sum  + ne_tot  * (po(l) + 1)**3
          ne_leaf_sum = ne_leaf_sum + ne_leaf
          np_leaf_sum = np_leaf_sum + ne_leaf * (po(l) + 1)**3
        end if
      end do

      if (rank == 0) then
        write(*,*)
        write(*,'(2X,A)') 'overall metrics'
        write(*,'(T5,A,T16,I0)')     'ne_tot  =', ne_tot_sum
        write(*,'(T5,A,T16,I0)')     'np_tot  =', np_tot_sum
        write(*,'(T5,A,T16,I0)')     'ne_leaf =', ne_leaf_sum
        write(*,'(T5,A,T16,I0)')     'np_leaf =', np_leaf_sum
        write(*,'(T5,A,T15,ES12.5)') 'r_l2    =', r_l2
        write(*,'(T5,A,T15,ES12.5)') 'e_l2    =', e_l2
      end if

    end associate

  end subroutine Evaluation

  !-----------------------------------------------------------------------------
  !> Set the adaptation marks

  subroutine SetAdaptationMarks

    type(ML_MeshVariable_3D), save :: qi
    real(RNP), save :: qi_max, qi_max_loc
    integer,   save :: n_ref_loc, n_ref

    real(RNP) :: qi_remove, qi_refine
    integer   :: i, l

    associate(mesh => ml_mesh%mesh)

      n_ref_loc = 0

      ! mark for global refinement ..............................................

      do l = 1, min(l_top,l_adapt-2)
        if (mesh(l)%n_elem > 0) then
          call mesh(l) % element % MarkForRefinement()
          n_ref_loc = n_ref_loc + 1
        end if
      end do

      ! quantity of interest for adaptation .....................................

      call qi % Init(ml_mesh, po = 0, nc = 1)
      qi_max_loc = 0

      ! evaluation
      do l = l_adapt-1, l_top
        associate( mesh_l => ml_op % sem(l) % mesh          &
                 , mm_l   => mm % level(l) % val(:,:,:,:,1) &
                 , qi_l   => qi % level(l) % val(0,0,0,:,1) &
                 , e_l    => e  % level(l) % val(:,:,:,:,1) )

          call SetArray(qi_l, ZERO)

          do i = 1, mesh_l % n_elem_active
            select case(adapt_criterion)
            case(1)
              qi_l(i) = sum(mm_l(:,:,:,i) * e_l(:,:,:,i)**2)
            case default
              qi_l(i) = sum(mm_l(:,:,:,i) * e_l(:,:,:,i)**2) &
                      / sum(mm_l(:,:,:,i))
            end select
            qi_max_loc = max(qi_max_loc, qi_l(i))
         end do

        end associate
      end do

      ! global maximum
      call XMPI_Reduce(qi_max_loc, qi_max, MPI_MAX, 0, comm)

      ! mark for adaptation ....................................................

      select case(adapt_criterion)
      case(1:2)
        qi_remove = adapt_remove * qi_max
        qi_refine = adapt_refine * qi_max
      case default ! 3
        qi_remove = adapt_remove
        qi_refine = adapt_refine
      end select

      ! l < l_max: mark elements for removal or refinement
      do l = l_adapt-1, min(l_top,l_max-1)
        do i = 1, mesh(l) % n_elem
          associate( element => mesh(l) % element(i)      &
                   , qi => qi % level(l) % val(0,0,0,i,1) )
            if (qi < adapt_remove) then
              call element % MarkForRemoval()
            else if (qi > qi_refine) then
              call element % MarkForRefinement()
              n_ref_loc = n_ref_loc + 1
            else
              call element % Unmark()
            end if
          end associate
        end do
      end do

      ! l = l_max: mark elements for removal
      if (l_top == l_max) then
        do i = 1, mesh(l_top) % n_elem
          associate( element => mesh(l_top) % element(i)      &
                   , qi => qi % level(l_top) % val(0,0,0,i,1) )
            if (qi < adapt_remove) then
              call element % MarkForRemoval()
            else
              call element % Unmark()
            end if
          end associate
        end do
      end if

      ! number of leaf elements marked for refinement
      call XMPI_Reduce(n_ref_loc, n_ref, MPI_SUM, 0, comm)

      ! control output .........................................................

      if (rank == 0) then
        write(*,'(/,2X,A,I0)') 'elements marked for refinement: ', n_ref
      end if

    end associate

  end subroutine SetAdaptationMarks

  !=============================================================================

end program ML_Elliptic_Test_Adaptive

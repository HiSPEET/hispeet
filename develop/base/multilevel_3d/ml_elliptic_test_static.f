!> summary:  Testing elliptic multilevel solvers with static refinement
!> author:   Joerg Stiller
!> date:     2024/09/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program ML_Elliptic_Test_Static
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Logging_Levels
  use Array_Assignments
  use Array_Reductions

  use Elliptic_Problem__3D
  use Elliptic_Problem__Knotty__3D
  use Elliptic_Problem__Simple__3D
  use Elliptic_Problem__Sphere__3D
  use Elliptic_Problem__TGV_Pressure__3D

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Verify_Mesh__3D
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

  character(len=*), parameter :: default_case = 'ml_elliptic_test_static'
  character(len=80) :: case_name ! case name
  character(len=80) :: case_file ! case input file: trim(case_name).prm

  integer :: solution_method = 21
  ! 11:  CS-MG
  ! 12:  CS-MGCG
  ! 21:  FAS-MG

  logical :: check_hdf5 = .false.  ! write and re-read ML mesh before solving
  logical :: export_vtk = .false.  ! switch for VTK export

  namelist/control_prm/ solution_method, check_hdf5, export_vtk

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

  integer :: test_problem = 3 ! 1...6: simple_{1/2/3}d/knotty/TGV/sphere
  integer :: start_values = 0 ! 0/1/2: zero, exact, random

  namelist/problem_prm/ test_problem, start_values

  ! common
  real(RNP) :: lambda =  0    ! Helmholtz parameter
  real(RNP) :: nu_0   =  1    ! constant or mean diffusivity ν₀
  real(RNP) :: nu_1   =  0    ! diffusivity fluctuation amplitude ν₁
  real(RNP) :: d_nu   =  0    ! diffusivity fluctuation phase shift
  integer   :: k_nu   =  1    ! diffusivity fluctuation wave number
  integer   :: k_u    =  1    ! solution wave number

  namelist/problem_prm/ lambda, nu_0, nu_1, d_nu, k_nu, k_u

  ! sphere
  real(RNP) :: x_c(3) = -0.05 ! sphere center
  real(RNP) :: r_0    =  0.7  ! sphere radius
  real(RNP) :: alpha  =  200  ! radial scaling factor

  namelist/problem_prm/ x_c, r_0, alpha

  ! boundary conditions {'D','N','P'} ['D']
  character, allocatable :: bc(:)

  namelist/problem_prm/  bc

  class(EllipticProblem_3D), allocatable, save :: problem
  character(len=80) :: problem_name = ''

  ! solvers .....................................................................

  integer, allocatable, save :: po(:) ! sequence of polynomial orders
  integer :: fc_smooth = 0 ! fine-to-coarse discontinuity handling {0,1,2}
  integer :: fc_filter = 0 ! fine-to-coarse filter order, 0: no filtering

  type(ML_MeshOperators_3D),      save :: ml_op
  type(ML_DG_EllipticSolver_3D),  save :: ml_elliptic
  type(ML_DG_EllipticOptions_3D), save :: ml_elliptic_opt

  namelist/solver_prm/ po, fc_smooth, fc_filter, ml_elliptic_opt

  ! variables ...................................................................

  character(len=:), allocatable, save :: name_var(:) ! names of variables

  type(ML_MeshVariable_3D), save :: var ! container of variables
  type(ML_MeshVariable_3D), save :: mm  ! diagonal mass matrix
  type(ML_MeshVariable_3D), save :: nu  ! diffusivity
  type(ML_MeshVariable_3D), save :: f   ! sources
  type(ML_MeshVariable_3D), save :: s   ! exact solution
  type(ML_MeshVariable_3D), save :: u   ! numerical solution
  type(ML_MeshVariable_3D), save :: e   ! error
  type(ML_MeshVariable_3D), save :: r   ! residual or just workspace

  type(ML_BoundaryVariable_3D), save :: bv ! boundary values

  ! auxiliaries ................................................................

  character(:), allocatable, save :: domain_name
  real(RNP), dimension(:), allocatable, save :: e0_mx, r0_e2, r0_mx
  real(RNP), dimension(:), allocatable, save :: en_mx, rn_e2, rn_mx
  real(RNP), save :: e0_l2, en_l2, r0_l2, rn_l2
  real(RNP), save :: t_start, t_solve
  integer,   save :: n_i

  logical :: exists, passed, all_passed, singular
  integer :: io, stat
  integer :: dim
  integer :: l_top, ne_max, ne_min, ne_tot

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
      call Error('ML_Elliptic_Test_Static', &
                 'file "'// trim(case_file) //'" not found')
    end if

  end if

  ! globalize multilevel mesh options
  call ml_mesh_opt % Bcast(0, comm)

  ! globalize logging levels
  call XMPI_Bcast_LoggingLevels(0, comm)

  ! globalize remaining parameters
  call XMPI_Bcast(case_name       , 0, comm)
  call XMPI_Bcast(check_hdf5      , 0, comm)
  call XMPI_Bcast(export_vtk      , 0, comm)
  call XMPI_Bcast(solution_method , 0, comm)
  call XMPI_Bcast(test_domain     , 0, comm)
  call XMPI_Bcast(gmsh_file       , 0, comm)

  ! base mesh ..................................................................

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
    write(*,'(/,2A)') 'verifying base mesh for ', trim(domain_name)
  end if

  call VerifyMesh_3D(base_mesh, passed)
  call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, comm)

  if (rank == 0) then
    write(*,'(2X,A,G0)') 'passed = ', all_passed
  end if

  call MPI_Barrier(comm)

  ! multilevel mesh ............................................................

  ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)
  l_top   = size(ml_mesh%mesh)

  if (check_hdf5) then
    call ml_mesh % WriteHDF5(case_name)
    deallocate(ml_mesh % mesh)
    call ml_mesh % ReadHDF5(case_name, comm)
  end if

  associate(mesh => ml_mesh%mesh)
    block
      integer :: l

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

    end block
  end associate

  ! problem ....................................................................

  allocate(bc(base_mesh % n_bound), source = 'D')

  if (rank == 0) then
    write(*,'(/,A)') 'initializing elliptic problem'
    open(newunit = io, file = case_file)
    read(io, nml = problem_prm)
    close(io)
  end if

  ! globalize problem parameters
  call XMPI_Bcast( test_problem   , 0, comm )
  call XMPI_Bcast( start_values   , 0, comm )
  call XMPI_Bcast( lambda         , 0, comm )
  call XMPI_Bcast( nu_0           , 0, comm )
  call XMPI_Bcast( nu_1           , 0, comm )
  call XMPI_Bcast( d_nu           , 0, comm )
  call XMPI_Bcast( k_nu           , 0, comm )
  call XMPI_Bcast( k_u            , 0, comm )
  call XMPI_Bcast( x_c            , 0, comm )
  call XMPI_Bcast( r_0            , 0, comm )
  call XMPI_Bcast( alpha          , 0, comm )
  call XMPI_Bcast( bc             , 0, comm )

  select case(test_problem)
  case(1:3)
    dim = test_problem
    write(problem_name,'(A,I0,A)') 'Simple_', dim, 'D'
    problem = EllipticProblem_Simple_3D(lambda, nu_0, nu_1, d_nu, k_nu, k_u, dim)
  case(4)
    problem_name = 'Knotty'
    problem = EllipticProblem_Knotty_3D(lambda, nu_0, nu_1, d_nu, k_nu, k_u)
  case(5)
    problem_name = 'TGV_Pressure'
    problem = EllipticProblem_TGV_Pressure_3D(lambda, nu_0, nu_1, d_nu, k_nu)
  case(6)
    problem_name = 'Sphere'
    problem = EllipticProblem_Sphere_3D(lambda, x_c, r_0, alpha)
  end select

  ! enforce periodicity at coupled boundaries
  where(base_mesh % boundary % coupled > 0) bc = 'P'

  singular = lambda == ZERO .and. all(bc == 'P' .or. bc == 'N')

  if (rank == 0) then
    write(*,'(T3,A,T26,9(G0,X))') 'problem name:'        , trim(problem_name)
    write(*,'(T3,A,T26,9(G0,X))') 'variable diffusivity:', problem % nu_1 > 0
    write(*,'(T3,A,T26,9(G0,X))') 'boundary conditions:' , bc
  end if

  ! solver ......................................................................

  allocate(po(l_top), source = -1)

  if (rank == 0) then
    open(newunit = io, file = case_file)
    read(io, nml = solver_prm, iostat = stat)
    close(io)
    write(*,'(/,A)') 'initializing multilevel operators'
  end if

  call XMPI_Bcast(po       , 0, comm)
  call XMPI_Bcast(fc_smooth, 0, comm)
  call XMPI_Bcast(fc_filter, 0, comm)
  ml_op = ML_MeshOperators_3D(ml_mesh, po, smooth=fc_smooth, filter=fc_filter)

  if (rank == 0) then
    write(*,'(/,A)') 'initializing multilevel solver'
  end if

  call ml_elliptic_opt % Bcast(0, comm)
  ml_elliptic = ML_DG_EllipticSolver_3D(ml_op, ml_elliptic_opt)

  associate(mesh => ml_mesh%mesh)
    block
      integer :: l

      if (rank == 0) then
        write(*,'(4X,A,4(2X,A7),2(3X,A2),8X,A2)')  &
            'l','n_parts','min(ne)','max(ne)','sum(ne)','po','no','ns'
      end if

      do l = 1, l_top
        if (mesh(l)%part >= 0) then
          call XMPI_Reduce(mesh(l)%n_elem, ne_min, MPI_MIN, 0, mesh(l)%comm_parts)
          call XMPI_Reduce(mesh(l)%n_elem, ne_max, MPI_MAX, 0, mesh(l)%comm_parts)
          call XMPI_Reduce(mesh(l)%n_elem, ne_tot, MPI_SUM, 0, mesh(l)%comm_parts)
        end if
        if (mesh(l)%part == 0) then
          write(*,'(I5,4I9,2I5,I6)',advance='NO') &
              l, mesh(l)%n_parts, ne_min, ne_max, ne_tot, &
              po(l), ml_elliptic % elliptic_op(l) % schwarz % no
          if (l == 1) then
            write(*,'(I10)') ml_elliptic % i_crs
          else
            write(*,'(2X,2(I3,X,A))') &
                ml_elliptic % NumSmoothingSteps(l, stage = 1), '+', &
                ml_elliptic % NumSmoothingSteps(l, stage = 2)
          end if
        end if
      end do

    end block
  end associate

  ! variables ..................................................................

  if (rank == 0) then
    write(*,'(/,A)') 'initializing multilevel variables'
  end if

  name_var = [ 'nu', 'f ', 's ', 'u ', 'e ', 'r ' ]

  call var % Init(ml_op, size(name_var), name_var)

  ! handles for accessing individual variables
  call var % GetSlice(nu, first = 1, last = 1)
  call var % GetSlice(f , first = 2, last = 2)
  call var % GetSlice(s , first = 3, last = 3)
  call var % GetSlice(u , first = 4, last = 4)
  call var % GetSlice(e , first = 5, last = 5)
  call var % GetSlice(r , first = 6, last = 6)

  ! structure for keeping diagonal mass matrices
  call mm % Init(ml_op, 1)

  ! structure for keeping boundary values
  call bv % Init(ml_op, nc = 1)

  allocate(e0_mx(l_top), source = ZERO)
  allocate(en_mx, r0_e2, r0_mx, rn_e2, rn_mx, source = e0_mx)

  block
    real(RNP), allocatable, save :: q(:,:,:,:,:)
    integer :: i, l

    do l = 1, l_top
      associate( sem_l => ml_op % sem(l)                 &
               , x_l   => ml_op % sem(l) % metrics % x   &
               , mm_l  => mm % level(l) % val(:,:,:,:,1) &
               , nu_l  => nu % level(l) % val(:,:,:,:,1) &
               , f_l   => f  % level(l) % val(:,:,:,:,1) &
               , r_l   => r  % level(l) % val(:,:,:,:,1) &
               , s_l   => s  % level(l) % val(:,:,:,:,1) &
               , u_l   => u  % level(l) % val(:,:,:,:,1) &
               , bu_l  => bv % level(l) % var            )

        call sem_l % Get_DG_DiagonalMassMatrix(mm_l)

        call problem % GetExactSolution (x_l, s_l)
        call problem % GetDiffusivity   (x_l, nu_l)
        call problem % GetSource        (x_l, r_l)     ! r = λ u - ∇·(ν ∇u)
        f_l = mm_l * r_l                               ! f = M r

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

        if (any(bc == 'N')) then
          allocate(q, mold = x_l)
          call problem % GetExactGradient(x_l, q)
          q(:,:,:,:,1) = nu_l * q(:,:,:,:,1)
          q(:,:,:,:,2) = nu_l * q(:,:,:,:,2)
          q(:,:,:,:,3) = nu_l * q(:,:,:,:,3)
        end if

        ! extract boundary conditions
        do i = 1, size(bc)
          select case(bc(i))
          case('D')
            call bu_l(i) % Extract(s_l)
          case('N')
            call bu_l(i) % ExtractNormalComponent(sem_l, q)
          end select
        end do

        if (allocated(q)) deallocate(q)

      end associate
    end do

  end block

  !-----------------------------------------------------------------------------
  ! Initial error and residual norms

  call Evaluation(r0_e2, r0_mx, r0_l2, e0_mx, e0_l2)

  !-----------------------------------------------------------------------------
  ! Solution

  if (rank == 0) then
    write(*,'(/,A)') 'executing multilevel solver'
    t_start = MPI_Wtime()
  end if

  if (problem % nu_1 > 0) then
    select case(solution_method)
    case(11)
      call ml_elliptic % CS_MG_Solver(bc, lambda, nu, bv, f, u, ni=n_i)
    case(12)
      call ml_elliptic % CS_MGCG_Solver(bc, lambda, nu, bv, f, u, ni=n_i)
    case(21)
      call ml_elliptic % FAS_MG_Solver(bc, lambda, nu, bv, f, u, ni=n_i)
    end select
  else
    select case(solution_method)
    case(11)
      call ml_elliptic % CS_MG_Solver(bc, lambda, nu_0, bv, f, u, ni=n_i)
    case(12)
      call ml_elliptic % CS_MGCG_Solver(bc, lambda, nu_0, bv, f, u, ni=n_i)
    case(21)
      call ml_elliptic % FAS_MG_Solver(bc, lambda, nu_0, bv, f, u, ni=n_i)
    end select
  end if

  if (rank == 0) then
    t_solve = MPI_Wtime() - t_start
  end if

  !-----------------------------------------------------------------------------
  ! Evaluation

  if (rank == 0) then
    write(*,'(/,A)') 'evaluation'
  end if

  ! final error and residual norms
  call Evaluation(rn_e2, rn_mx, rn_l2, en_mx, en_l2)

  if (rank == 0) then
    write(*,*)
    write(*,'(A,I0)')     'num iterations:  n_i = ', n_i
    write(*,'(A,ES10.3)') 'solution time:   t_s =' , t_solve

    write(*,'(/,A,/)') 'L2 residual and error over leaf elements'

    write(*,'(2X,A,ES10.3)') 'r0_l2  = ', r0_l2
    write(*,'(2X,A,ES10.3)') 'rn_l2  = ', rn_l2
    write(*,'(2X,A,ES10.3)') 'e0_l2  = ', e0_l2
    write(*,'(2X,A,ES10.3)') 'en_l2  = ', en_l2

    write(*,'(/,A,/)') 'residuals and errors over active elements per level'
    block
      integer :: l, l_min

      select case(solution_method)
      case(11,12)
        l_min = l_top
      case default
        l_min = 1
      end select

      write(*,'(4X,A,3X,7(A7,5X))') 'l', &
                                    '  r0_e2', '  rn_e2', &
                                    '  r0_mx', '  rn_mx', &
                                    '  e0_mx', '  en_mx', &
                                    '-lg rho'
      do l = l_min, l_top
        write(*,'(I5,7(2X,ES10.3))') l, &
                                     r0_e2(l), rn_e2(l), &
                                     r0_mx(l), rn_mx(l), &
                                     e0_mx(l), en_mx(l), &
                                     log10(r0_e2(l) / rn_e2(l)) / n_i
      end do
    end block

  end if

  !-----------------------------------------------------------------------------
  ! VTK export

  if (export_vtk) then

    ! mesh and variables
    call var % ExportVTK(ml_op, trim(case_name)//'_full', mode=1)
    call var % ExportVTK(ml_op, trim(case_name)//'_leaf', mode=3)

    ! space filling curve
    associate(mesh => ml_mesh%mesh)
      block
        character(len=9) :: tag
        integer :: l
        do l = 1, size(mesh)
          if (mesh(l) % has_sfc) then
            write(tag,'(A,I0,A)') '_sfc_l', l
            call ExportVTK_MeshSFC(mesh(l), file = trim(case_name)//trim(tag))
          end if
        end do
      end block
    end associate

  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !> Evaluation

  subroutine Evaluation(r_e2, r_mx, r_l2, e_mx, e_l2)
    real(RNP), contiguous, intent(out) :: r_e2(:)
    real(RNP), contiguous, intent(out) :: r_mx(:)
    real(RNP),             intent(out) :: r_l2
    real(RNP), contiguous, intent(out) :: e_mx(:)
    real(RNP),             intent(out) :: e_l2

    real(RNP), allocatable, save :: e_mx_loc(:), r_mx_loc(:), r_e2_loc(:)
    real(RNP), save :: int_1, int_1_loc
    real(RNP), save :: int_e, int_e_loc

    real(RNP) :: e_avg
    integer   :: i, l

    allocate(e_mx_loc(l_top), source = ZERO)
    allocate(r_mx_loc, r_e2_loc, source = e_mx_loc)

    ! residual .................................................................

    if (problem % nu_1 > 0) then
      call ml_elliptic % FAS_MG_Residual(bc, lambda, nu, bv, f, u, r)
    else
      call ml_elliptic % FAS_MG_Residual(bc, lambda, nu_0, bv, f, u, r)
    end if

    ! maximum norm
    do l = 1, l_top
      associate( mesh_l => ml_op % sem(l) % mesh         &
               , r_l    => r % level(l) % val(:,:,:,:,1) )

        do i = 1, mesh_l % n_elem_active
          r_mx_loc(l) = max(r_mx_loc(l), maxval(abs(r_l(:,:,:,i))))
          r_e2_loc(l) = r_e2_loc(l) + sum(r_l(:,:,:,i)**2)
        end do

      end associate
    end do
    call XMPI_Reduce(r_mx_loc, r_mx, MPI_MAX, 0, comm)
    call XMPI_Reduce(r_e2_loc, r_e2, MPI_SUM, 0, comm)
    r_e2 = sqrt(r_e2)

    ! L2 norms
    r_l2 = sqrt(ML_WeightedScalarProduct_3D(mm, r, r, leaf = .true.))

    ! error ....................................................................

    int_1_loc = 0
    int_e_loc = 0

    do l = 1, l_top
      associate( mesh_l => ml_op % sem(l) % mesh          &
               , mm_l   => mm % level(l) % val(:,:,:,:,1) &
               , s_l    => s  % level(l) % val(:,:,:,:,1) &
               , u_l    => u  % level(l) % val(:,:,:,:,1) &
               , e_l    => e  % level(l) % val(:,:,:,:,1) )

        do i = 1, size(e_l,4)
          e_l(:,:,:,i) = u_l(:,:,:,i) - s_l(:,:,:,i)
          if (mesh_l%element(i)%IsLeaf()) then
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
      associate( mesh_l => ml_op % sem(l) % mesh         &
               , e_l    => e % level(l) % val(:,:,:,:,1) )

        do i = 1, mesh_l % n_elem_active
          e_l(:,:,:,i) = e_l(:,:,:,i) - e_avg
          e_mx_loc(l)  = max(e_mx_loc(l), maxval(abs(e_l(:,:,:,i))))
        end do

      end associate
    end do
    call XMPI_Reduce(e_mx_loc, e_mx, MPI_MAX, 0, comm)

    ! L2 norm
    e_l2 = sqrt(ML_WeightedScalarProduct_3D(mm, e, e, leaf = .true.))

    deallocate(r_mx_loc, r_e2_loc, e_mx_loc)

  end subroutine Evaluation

  !=============================================================================

end program ML_Elliptic_Test_Static

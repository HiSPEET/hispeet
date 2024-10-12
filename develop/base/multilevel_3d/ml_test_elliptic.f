!> summary:  Program for testing the elliptic multilevel solvers
!> author:   Joerg Stiller
!> date:     2024/09/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program ML_Test_Elliptic
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Logging_Levels
  use Array_Assignments
  use Array_Reductions

  use Elliptic_Problem__3D
  use Elliptic_Problem__Simple__3D
  use Elliptic_Problem__Knotty__3D
  use Elliptic_Problem__TGV_Pressure__3D

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Verify_Mesh__3D

  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  use ML__Boundary_Variable__3D
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

  character(len=*), parameter :: default_case = 'ml_test_elliptic'
  character(len=80) :: case_name ! case name
  character(len=80) :: case_file ! case input file: trim(test_case).prm

  namelist/control_prm/ log_level
  namelist/control_prm/ log_level_inner_iteration
  namelist/control_prm/ log_level_outer_iteration
  namelist/control_prm/ log_level_multigrid_cycle

  logical :: export_vtk = .false.  ! switch for VTK export

  namelist/control_prm/ export_vtk

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

  integer :: test_problem    =  3      ! 1/2/3/4/5: simple_{1/2/3}d/knotty/TGV
  logical :: has_variable_nu = .false. ! T/F: variable/constant ν

  namelist/problem_prm/ test_problem, has_variable_nu

  ! NOTE
  ! The fluctuation amplitude ν₁ is ignored in case of constant ν

  real(RNP) :: lambda = 0         ! Helmholtz parameter
  real(RNP) :: nu_0   = 1         ! diffusivity mean value ν₀
  real(RNP) :: nu_1   = 0         ! diffusivity fluctuation amplitude ν₁
  real(RNP) :: d_nu   = 0         ! diffusivity fluctuation phase shift
  integer   :: k_nu   = 1         ! diffusivity fluctuation wave number
  integer   :: k_u    = 1         ! solution wave number
  character, allocatable :: bc(:) ! boundary conditions {'D','N','P'} ['D']

  namelist/problem_prm/ lambda, nu_0, nu_1, d_nu, k_nu, k_u, bc

  class(EllipticProblem_3D), allocatable, save :: problem
  character(len=:),          allocatable, save :: problem_name

  ! solvers .....................................................................

  integer, allocatable, save :: po(:) ! sequence of polynomial orders

  type(ML_MeshOperators_3D),      save :: ml_op
  type(ML_DG_EllipticSolver_3D),  save :: ml_elliptic
  type(ML_DG_EllipticOptions_3D), save :: ml_elliptic_opt

  namelist/solver_prm/ po, ml_elliptic_opt

  ! variables ...................................................................

  character(len=:), allocatable, save :: name_var(:) ! names of variables

  type(ML_MeshVariable_3D), save :: ml_var ! container of variables
  type(ML_MeshVariable_3D), save :: ml_nu  ! diffusivity
  type(ML_MeshVariable_3D), save :: ml_f   ! sources
  type(ML_MeshVariable_3D), save :: ml_s   ! exact solution
  type(ML_MeshVariable_3D), save :: ml_u   ! numerical solution
  type(ML_MeshVariable_3D), save :: ml_e   ! error
  type(ML_MeshVariable_3D), save :: ml_r   ! residual or just workspace

  type(ML_BoundaryVariable_3D), save :: ml_bv ! boundary values

  ! auxiliaries ................................................................

  character(:), allocatable, save :: domain_name
  real(RNP), allocatable, save :: e_max(:), e_min(:)
  real(RNP), allocatable, save :: r_max(:), r_rms(:)
  real(RNP), allocatable, save :: v_loc(:)

  real(RNP) :: rr
  logical   :: exists, passed, all_passed
  integer   :: io, stat
  integer   :: dim
  integer   :: l_top, ne_max, ne_min, ne_tot, n_i
  integer   :: i, l

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
      call Error('ML_Test_Elliptic', 'file "'// trim(case_file) //'" not found')
    end if

  end if

  ! globalize multilevel mesh options
  call ml_mesh_opt % Bcast(0, comm)

  ! globalize logging levels
  call XMPI_Bcast_LoggingLevels(0, comm)

  ! globalize remaining parameters
  call XMPI_Bcast(export_vtk , 0, comm)
  call XMPI_Bcast(test_domain, 0, comm)
  call XMPI_Bcast(gmsh_file  , 0, comm)

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
    write(*,'(/,A)') 'initializing elliptic problem'
    open(newunit = io, file = case_file)
    read(io, nml = problem_prm)
    close(io)
    if (.not. has_variable_nu) then
      nu_1 = 0
    end if
  end if

  ! globalize problem parameters
  call XMPI_Bcast( test_problem   , 0, comm )
  call XMPI_Bcast( has_variable_nu, 0, comm )
  call XMPI_Bcast( lambda         , 0, comm )
  call XMPI_Bcast( nu_0           , 0, comm )
  call XMPI_Bcast( nu_1           , 0, comm )
  call XMPI_Bcast( d_nu           , 0, comm )
  call XMPI_Bcast( k_nu           , 0, comm )
  call XMPI_Bcast( k_u            , 0, comm )
  call XMPI_Bcast( bc             , 0, comm )

  select case(test_problem)
  case(1:3)
    dim = test_problem
    problem_name = 'Simple xD'
    write(problem_name(8:8),'(I1)') dim
    problem = EllipticProblem_Simple_3D &
                  (lambda, nu_0, nu_1, d_nu, k_nu, k_u, dim)
  case(4)
    problem_name = 'Knotty'
    problem = EllipticProblem_Knotty_3D &
                  (lambda, nu_0, nu_1, d_nu, k_nu, k_u)
  case default
    problem_name = 'TGV_Pressure'
    problem = EllipticProblem_TGV_Pressure_3D &
                  (lambda, nu_0, nu_1, d_nu, k_nu, k_u)
  end select

  ! enforce periodicity at coupled boundaries
  where(base_mesh % boundary % coupled > 0) bc = 'P'

  if (rank == 0) then
    write(*,'(T3,A,T26,9(G0,X))') 'problem name:'        , trim(problem_name)
    write(*,'(T3,A,T26,9(G0,X))') 'variable diffusivity:', has_variable_nu
    write(*,'(T3,A,T26,9(G0,X))') 'boundary conditions:' , bc
  end if

  ! solver ......................................................................

  allocate(po(l_top), source = 1)

  if (rank == 0) then
    open(newunit = io, file = case_file)
    read(io, nml = solver_prm)
    close(io)
    write(*,'(/,A)') 'initializing multilevel operators'
  end if

  call XMPI_Bcast(po, 0, comm)
  ml_op = ML_MeshOperators_3D(ml_mesh, po)

  associate(mesh => ml_mesh%mesh)
    do l = 1, l_top
      if (mesh(l)%part >= 0) then
        call XMPI_Reduce(mesh(l)%n_elem, ne_min, MPI_MIN, 0, mesh(l)%comm_parts)
        call XMPI_Reduce(mesh(l)%n_elem, ne_max, MPI_MAX, 0, mesh(l)%comm_parts)
        call XMPI_Reduce(mesh(l)%n_elem, ne_tot, MPI_SUM, 0, mesh(l)%comm_parts)
      end if
      if (mesh(l)%part == 0) then
        write(*,'(2X,A,I4,A,I5,A,3(X,A,X,I6),X,A,X,I3)') &
          'level ',l,': n_parts =',mesh(l)%n_parts,',  ', &
          'min/max/sum(n_elem)/po =',ne_min,'/',ne_max,'/',ne_tot,'/',po(l)
      end if
    end do
  end associate

  if (rank == 0) then
    write(*,'(/,A)') 'initializing multilevel solver'
  end if

  call ml_elliptic_opt % Bcast(0, comm)
  ml_elliptic = ML_DG_EllipticSolver_3D(ml_op, ml_elliptic_opt, bc)

  ! variables ..................................................................

  if (rank == 0) then
    write(*,'(/,A)') 'initializing multilevel variables'
  end if

  name_var = [ 'nu', 'f ', 's ', 'u ', 'e ', 'r ' ]

  ml_var = ML_MeshVariable_3D(ml_op, size(name_var), name_var)

  ! handles for accessing individual variables
  call ml_var % GetSlice(ml_nu, first = 1, last = 1)
  call ml_var % GetSlice(ml_f , first = 2, last = 2)
  call ml_var % GetSlice(ml_s , first = 3, last = 3)
  call ml_var % GetSlice(ml_u , first = 4, last = 4)
  call ml_var % GetSlice(ml_e , first = 5, last = 5)
  call ml_var % GetSlice(ml_r , first = 6, last = 6)

  ml_bv = ML_BoundaryVariable_3D(ml_op, nc = 1)

  allocate(v_loc(l_top), source = ZERO)
  allocate(e_min, e_max, r_max, r_rms, source = v_loc)

  block
    real(RNP), allocatable, save :: q(:,:,:,:,:)

    do l = 1, l_top
      associate( sem  => ml_op % sem(l)                    &
               , x    => ml_op % sem(l)   % metrics % x    &
               , nu   => ml_nu % level(l) % val(:,:,:,:,1) &
               , f    => ml_f  % level(l) % val(:,:,:,:,1) &
               , r    => ml_r  % level(l) % val(:,:,:,:,1) &
               , s    => ml_s  % level(l) % val(:,:,:,:,1) &
               , u    => ml_u  % level(l) % val(:,:,:,:,1) &
               , mm   => ml_e  % level(l) % val(:,:,:,:,1) &
               , bv_u => ml_bv % level(l) % var            )

        call problem % GetExactSolution (x, s)
        call problem % GetDiffusivity   (x, nu)
        call problem % GetSource        (x, r)
        call SetArray(u, ZERO)
!!         call SetArray(u, s)
!! if (l == l_top) call SetArray(u, s)

        ! r = λ u - ∇·(ν ∇u)
        call sem % Get_DG_DiagonalMassMatrix(mm)

        ! project source: f = M r
        f = mm * r

        if (any(bc == 'N')) then
          allocate(q, mold = x)
          call problem % GetExactGradient(x, q)
          q(:,:,:,:,1) = nu * q(:,:,:,:,1)
          q(:,:,:,:,2) = nu * q(:,:,:,:,2)
          q(:,:,:,:,3) = nu * q(:,:,:,:,3)
        end if

        ! extract and apply boundary conditions
        do i = 1, size(bc)
          select case(bc(i))
          case('D')
            call bv_u(i) % Extract(s)
          case('N')
            call bv_u(i) % ExtractNormalComponent(sem, q)
          end select
        end do

        if (allocated(q)) deallocate(q)

      end associate
    end do

  end block

  !-----------------------------------------------------------------------------
  ! Solution

  if (rank == 0) then
    write(*,'(/,A)') 'executing multilevel solver'
  end if

  if (has_variable_nu) then
    call ml_elliptic % MG_Solver(lambda, ml_nu, ml_u, ml_f, ml_bv, n_i)
  else
    call ml_elliptic % MG_Solver(lambda, nu_0, ml_u, ml_f, ml_bv, n_i)
  end if

  !-----------------------------------------------------------------------------
  ! Evaluation

  do l = 1, l_top
    associate( sem  => ml_op % sem(l)                    &
             , nu   => ml_nu % level(l) % val(:,:,:,:,1) &
             , f    => ml_f  % level(l) % val(:,:,:,:,1) &
             , s    => ml_s  % level(l) % val(:,:,:,:,1) &
             , u    => ml_u  % level(l) % val(:,:,:,:,1) &
             , e    => ml_e  % level(l) % val(:,:,:,:,1) &
             , r    => ml_r  % level(l) % val(:,:,:,:,1) &
             , bv_u => ml_bv % level(l) % var            )

      if (has_variable_nu) then
        call ml_elliptic % elliptic_op(l) % Residual(lambda, nu, f, bv_u, u, r)
      else
        call ml_elliptic % elliptic_op(l) % Residual(lambda, nu_0, f, bv_u, u, r)
      end if

      if (sem%mesh%part >= 0) then
        rr = ScalarProduct(r, r, sem%mesh%comm_parts)
      else
        rr = 0
      end if

      r_rms(l) = sqrt(rr)
      r_max(l) = maxval(abs(r))

      do i = 1, size(u,4)
        e(:,:,:,i) = u(:,:,:,i) - s(:,:,:,i)
        e_min(l) = min(e_min(l), minval(e(:,:,:,i)))
        e_max(l) = max(e_max(l), maxval(e(:,:,:,i)))
      end do

    end associate
  end do

  v_loc = r_rms; call XMPI_Reduce(v_loc, r_rms, MPI_MAX, 0, comm)
  v_loc = r_max; call XMPI_Reduce(v_loc, r_max, MPI_MAX, 0, comm)
  v_loc = e_min; call XMPI_Reduce(v_loc, e_min, MPI_MIN, 0, comm)
  v_loc = e_max; call XMPI_Reduce(v_loc, e_max, MPI_MAX, 0, comm)

  if (rank == 0) then
    write(*,'(/,A)') 'number of iterations'
    write(*,'(2X,A,I0)') 'n_i  = ', n_i

    write(*,'(/,A,/)') 'error metrics'
    write(*,'(4X,A,5(4X,A5,3X))') 'l', 'r_rms', 'r_max', 'e_min', 'e_max'
    do l = 1, l_top
      write(*,'(I5,5(2X,ES10.3))') l, r_rms(l), r_max(l), e_min(l), e_max(l)
    end do
  end if

  !-----------------------------------------------------------------------------
  ! VTK export

  if (export_vtk) then
    call ml_var % ExportVTK(ml_op, trim(case_file)//'_full', mode=1)
    call ml_var % ExportVTK(ml_op, trim(case_file)//'_leaf', mode=3)
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

  !=============================================================================

end program ML_Test_Elliptic

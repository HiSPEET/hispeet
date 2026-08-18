!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Validation of incompressible Navier-Stokes time integrators
!> author:   Joerg Stiller
!> date:     2022/10/19
!>
!> If present, the first argument of the invoking command will be interpreted
!> as the base name of the control file. If omitted, the program looks for
!> `ins_timeintegrator_3d_test.prm`.
!===============================================================================

program INS_Integrator_3D_Test
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Array_Assignments
  use Array_Reductions
  use Logging_Levels

  use Mesh__3D
  use Boundary_Variable__3D
  use TPO__Div__3D
  use Trace_Operators__3D
  use Volume_Integrals__3D
  use Surface_Integrals__3D
  use Smooth_Mesh_Data__3D
  use Export_VTK_Volume_Data__3D
  use Spectral_Element_Mesh__3D

  use INS__Problem__3D
  use INS__Problem__Test_Suite__3D

  use INS__Operator__3D
  use INS__Time_Scales__3D
  use INS__Integrator__3D
  use INS__Integrator__Euler__3D
  use INS__Integrator__BDF__3D
  use INS__Integrator__Runge_Kutta__3D
  use INS__Flow_Characteristics__3D

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  use Generic_Mesh__3D
  use Verify_Mesh__3D
  use Import_GMSH__3D

  use Root_Mesh_Partitioning__3D

  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__DG__Elliptic_Solver__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm      ! MPI communicator
  integer        :: rank      ! local MPI rank
  integer        :: n_proc    ! number of MPI processes
  integer        :: n_thread  ! number of OpenMP threads

  ! control parameters .........................................................

  ! control of logging levels
  namelist/control_prm/ log_level
  namelist/control_prm/ log_level_inner_iteration
  namelist/control_prm/ log_level_outer_iteration
  namelist/control_prm/ log_level_multigrid_cycle

  character(len=*), parameter :: default_case = 'ins_timeintegrator_3d_test'
  character(len=80) :: flow_case ! flow case name
  character(len=80) :: case_file ! flow case input file: trim(flow_case).prm

  character(len=80) :: flow_problem = 'Vortex_TG' ! problem name
  character(len=80) :: problem_file = 'vortex_tg' ! problem parameters file

  integer :: flow_domain = 1
    ! computational flow domain (u/s = un/structured, r = regular, d = deformed)
    !   1  cuboidal domain with Cartesian mesh                             (s+r)
    !   2  cuboidal domain with unstructured "diamond" mesh                (u+d)
    !   3  cylindrical domain                                              (u+d)
    !   4  annular domain                                                  (u+d)
    !  10  import from GMSH

  character(len=80) :: raw_mesh_file = ''

  namelist/control_prm/ flow_problem, problem_file, flow_domain, raw_mesh_file

  logical :: mesh_stat  = .true.   ! show mesh statistics
  integer :: char_freq  = 1        ! characteristics output frequency
  integer :: avg_rate   = 0        ! sampling rate for averaging, 0 if none
  logical :: export_vtk = .false.  ! generate VTK files
  logical :: subdiv_vtk = .true.   ! use quadratic subdivision for VTK export

  namelist/control_prm/ mesh_stat, char_freq, avg_rate, export_vtk, subdiv_vtk

  ! restart options
  character(len=80) :: restart_tag_in  = ''  ! tag for restart input files
  character(len=80) :: restart_tag_out = ''  ! tag for restart output files
  namelist/control_prm/ restart_tag_in, restart_tag_out
    !
    ! restart input is read from
    !   - trim(flow_case)_trim(restart_tag_in)_mesh_<rank>.h5  for the mesh
    !   - trim(flow_case)_trim(restart_tag_in)_data_<rank>.h5  for flow data
    ! where <rank> is the process rank in mesh%comm_world
    !
    ! output written to
    !   - trim(flow_case)_trim(restart_tag_out)_mesh_<rank>.h5  for the mesh
    !   - trim(flow_case)_trim(restart_tag_out)_data_<rank>.h5  for flow data
    !
    ! no restart data is read or written if the corresponding tag is empty

  ! spatial parameters .........................................................

  type(ML_Mesh_Options_3D), save :: ml_mesh_opt
    ! multilevel mesh options, defining top level, partitioning and refinement,
    ! stacic components defined in namelist `ml_mesh_options_3d__static` and
    ! dynamic components in `ml_mesh_options_3d__dynamic`

  integer, save :: po_u = 8
    ! polynomial order of solution u \ p on top level
  integer, allocatable, save :: po_p(:)
    ! polynomial orders for pressure on levels 1:l_top,
    ! matching condition: po_p(l_top) = ins_op % eop_p % po

  type(ML_DG_EllipticOptions_3D), save :: ml_solver_p_opt
  type(INS_OperatorOptions_3D),   save :: ins_op_opt

  namelist/spatial_prm/ po_u, po_p, ml_solver_p_opt, ins_op_opt

  ! temporal parameters ........................................................

  real(RNP) :: t_end       = 1  ! final time
  real(RNP) :: dt          = 1  ! time step size
  integer   :: nt_max      = 0  ! max num time steps, < 0 if no limit
  integer   :: time_method = 1  ! 1/2/3: Euler/BDF/Runge-Kutta

  logical   :: smooth_initial_data = .false.
  integer   :: smooth_filter = 3  ! 1/2/3: cut-off/erfc-log/exponential
  integer   :: smooth_order  = 0  ! filter order (0: auto)

  type(INS_Integrator_Euler_Options_3D)      :: ins_ti_euler_opt
  type(INS_Integrator_BDF_Options_3D)        :: ins_ti_bdf_opt
  type(INS_Integrator_RungeKutta_Options_3D) :: ins_ti_rk_opt

  namelist/temporal_prm/ t_end, dt, nt_max
  namelist/temporal_prm/ smooth_initial_data, smooth_filter, smooth_order
  namelist/temporal_prm/ time_method
  namelist/temporal_prm/ ins_ti_euler_opt
  namelist/temporal_prm/ ins_ti_bdf_opt
  namelist/temporal_prm/ ins_ti_rk_opt

  ! mesh .......................................................................

  type(GenericMesh_3D),          save :: generic_mesh
  type(Mesh_3D), allocatable,    save :: base_mesh
  type(ML_Mesh_3D),              save :: ml_mesh
  type(ML_MeshOperators_3D),     save :: ml_op_p
  type(ML_DG_EllipticSolver_3D), save :: ml_solver_p

  ! problem, operators and solvers .............................................

  class(INS_Problem_3D), allocatable, save :: problem

  type(SpectralElementMesh_3D), save :: sem_u
  type(INS_Operator_3D),        save :: ins_op

  class(INS_Integrator_3D), allocatable, save :: ins_ti

  type(INS_TimeScales_3D)          :: time_scales
  type(INS_FlowCharacteristics_3D) :: flow_char

  ! variables ..................................................................

  real(RNP) :: t = 0 ! problem time

  real(RNP), allocatable, target, save :: var(:,:,:,:,:)
  character(len=20), allocatable, save :: var_name(:)

  real(RNP), pointer, contiguous, save :: u(:,:,:,:,:)     ! u = [v, p]
  real(RNP), pointer, contiguous, save :: v(:,:,:,:,:)     ! velocity
  real(RNP), pointer, contiguous, save :: p(:,:,:,:)       ! pressure
  real(RNP), pointer, contiguous, save :: mu(:,:,:,:)      ! bulk viscosity μ
  real(RNP), pointer, contiguous, save :: nu(:,:,:,:)      ! shear viscosity ν
  real(RNP), pointer, contiguous, save :: div_v(:,:,:,:)   ! ∇⋅v

  real(RNP), pointer, contiguous, save :: u_ex(:,:,:,:,:)  ! u_ex = [v_ex, p_ex]
  real(RNP), pointer, contiguous, save :: v_ex(:,:,:,:,:)  ! exact velocity
  real(RNP), pointer, contiguous, save :: p_ex(:,:,:,:)    ! exact pressure

  real(RNP), pointer, contiguous, save :: err_u(:,:,:,:,:) ! error, u - u_ex
  real(RNP), pointer, contiguous, save :: err_v(:,:,:,:,:) ! velocity error
  real(RNP), pointer, contiguous, save :: err_p(:,:,:,:)   ! pressure error

  real(RNP), pointer, contiguous, save :: q_avg(:,:,:,:,:) ! averaged quantities

  real(RNP), allocatable, save :: vp(:,:,:,:,:)  ! velocity trace
  real(RNP), allocatable, save :: w(:,:,:,:,:)   ! workspace

  type(BoundaryVariable_3D), allocatable, save :: bv_vn(:) ! n⋅v on Γ=∂Ω
  real(RNP), allocatable, save :: int_vn(:,:) ! ∫n⋅v dΓ

  ! auxiliaries ................................................................

  character(:), allocatable :: domain_name
  character(:), allocatable :: mesh_file
  character(:), allocatable :: data_file

  real(RNP) :: domain_volume
  logical   :: exists, last, restart_in, restart_out
  integer   :: io, stat
  integer   :: n_avg, n_bound, n_elem, n_elem_tot, n_ghost, n_point, n_var
  integer   :: l_top
  integer   :: i, nt
  real(RDP) :: time

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

  ! parameters .................................................................

  ! read control parameters
  if (rank == 0) then

    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Validation of incompressible Navier-Stokes time integrators'
    write(*,*)
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:', n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'  , n_thread
    write(*,*)

    call get_command_argument(1, flow_case, status=stat)
    if (stat /= 0 .or. len_trim(flow_case) == 0) then
      flow_case = default_case
    end if
    case_file = trim(flow_case) // '.prm'

    inquire(file=case_file, exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(case_file)
      open(newunit = io, file = case_file)
      read(io, nml = control_prm)
      close(io)
      if (char_freq < 1) then
        ! disable intermediate control output
        char_freq = huge(1)
      end if
    else
       call Error( 'INS_Integrator_3D_Test', &
                   'input file "' // trim(case_file) // '" not found' )
    end if

  end if

  ! globalize logging levels
  call XMPI_Bcast_LoggingLevels(0, comm)

  ! globalize control parameters
  call XMPI_Bcast(flow_case      , 0, comm)
  call XMPI_Bcast(case_file      , 0, comm)
  call XMPI_Bcast(flow_problem   , 0, comm)
  call XMPI_Bcast(problem_file   , 0, comm)
  call XMPI_Bcast(flow_domain    , 0, comm)
  call XMPI_Bcast(raw_mesh_file  , 0, comm)
  call XMPI_Bcast(mesh_stat      , 0, comm)
  call XMPI_Bcast(char_freq      , 0, comm)
  call XMPI_Bcast(avg_rate       , 0, comm)
  call XMPI_Bcast(export_vtk     , 0, comm)
  call XMPI_Bcast(subdiv_vtk     , 0, comm)
  call XMPI_Bcast(restart_tag_in , 0, comm)
  call XMPI_Bcast(restart_tag_out, 0, comm)

  ! restart switches
  restart_in  = len_trim(restart_tag_in)  > 0
  restart_out = len_trim(restart_tag_out) > 0

  ! domain name ................................................................

  select case(flow_domain)
  case(1)
    domain_name = 'Cuboidal domain with Cartesian mesh'
  case(2)
    domain_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    domain_name = 'Cylindrical domain with unstructured mesh'
  case(4)
    domain_name = 'Annular domain with unstructured mesh'
  case(10)
    domain_name = raw_mesh_file
  end select

  if (restart_in) then

    ! read multilevel mesh .....................................................

    mesh_file = trim(flow_case) // '_' // trim(restart_tag_in) // '_mesh'
    call ml_mesh % ReadHDF5(mesh_file, comm)

  else

    ! create base mesh .........................................................

    allocate(base_mesh)

    select case(flow_domain)
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
    case(10)
      if (rank == 0) then
        call ImportGMSH_3D(raw_mesh_file, generic_mesh)
      end if
      call base_mesh % ImportGenericMesh(generic_mesh, comm)
      domain_name = trim(raw_mesh_file)
    end select

    ! create multilevel mesh ...................................................

    if (rank == 0) then
      write(*,'(2X,A)') 'creating multilevel mesh'
      ! read multilevel mesh options
      open(newunit = io, file = case_file)
      ml_mesh_opt = ML_Mesh_Options_3D(io, n_proc)
      close(io)
    end if

    call ml_mesh_opt % Bcast(0, comm)
    ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)
    deallocate(base_mesh)

  end if

  l_top   = size(ml_mesh % mesh)
  n_elem  = ml_mesh % mesh(l_top) % n_elem
  n_ghost = ml_mesh % mesh(l_top) % n_ghost
  n_bound = ml_mesh % mesh(l_top) % n_bound

  ! problem ....................................................................

  call Set_INS_TestProblem_3D(problem, flow_problem, problem_file, n_bound, comm)

  ! check & fix boundary conditions
  do i = 1, n_bound
    if (ml_mesh % mesh(l_top) % boundary(i) % coupled > 0) then
      if (problem % bc_v(i) /= 'P') then
        if (rank == 0) then
          write(*,'(A,I0)') '  *** enforcing periodic BC on coupled boundary ',i
        end if
        problem % bc_v(i) = 'P'
        problem % bc_p(i) = 'P'
      end if
    end if
  end do

  ! spatial ....................................................................

  allocate(po_p(l_top), source = 1)

  ! read
  if (rank == 0) then
    write(*,'(2X,A)') 'creating spatial operators'
    open(newunit = io, file = case_file)
    read(io, nml = spatial_prm)
    close(io)
  end if

  ! globalize parameters
  call XMPI_Bcast(po_u, 0, comm)
  call XMPI_Bcast(po_p, 0, comm)

  ! globalize options
  call ml_solver_p_opt % Bcast(0, comm)
  call ins_op_opt      % Bcast(0, comm)

  ! multilevel pressure solver
  ml_op_p     = ML_MeshOperators_3D(ml_mesh, po_p)
  ml_solver_p = ML_DG_EllipticSolver_3D(ml_op_p, ml_solver_p_opt)

  ! Navier-Stokes operator
  sem_u = SpectralElementMesh_3D(ml_mesh%mesh(l_top), po_u)
  call ins_op % Init( opt         = ins_op_opt           &
                    , problem     = problem              &
                    , sem_u       = sem_u                &
                    , sem_p       = ml_op_p % sem(l_top) &
                    , ml_solver_p = ml_solver_p          &
                    , level       = l_top                )

  ! temporal ...................................................................

  ! read
  if (rank == 0) then
    write(*,'(2X,A)') 'creating temporal operators'
    open(newunit = io, file = case_file)
    read(io, nml = temporal_prm)
    close(io)
  end if

  ! globalize parameters
  call XMPI_Bcast(t_end              , 0, comm)
  call XMPI_Bcast(dt                 , 0, comm)
  call XMPI_Bcast(nt_max             , 0, comm)
  call XMPI_Bcast(time_method        , 0, comm)
  call XMPI_Bcast(smooth_initial_data, 0, comm)


  ! globalize options
  call ins_ti_euler_opt % Bcast(0, comm)
  call ins_ti_bdf_opt   % Bcast(0, comm)
  call ins_ti_rk_opt % Bcast(0, comm)

  ! time integrator
  select case(time_method)
  case(1)
    ins_ti = INS_Integrator_Euler_3D(problem, ins_op, ins_ti_euler_opt)
  case(2)
    ins_ti = INS_Integrator_BDF_3D(problem, ins_op, ins_ti_bdf_opt)
  case(3)
    ins_ti = INS_Integrator_RungeKutta_3D(problem, ins_op, ins_ti_rk_opt)
  end select

  ! variables ..................................................................

  n_var = 7

  if (problem % HasExactSolution()) then
    n_var = n_var + 8
  end if

  if (avg_rate > 0) then
    n_var = n_var + 10
  end if

  allocate(var(0:po_u,0:po_u,0:po_u,1:n_elem,1:n_var), source = ZERO)
  allocate(var_name(1:n_var))

  u(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,1:4)
  v(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,1:3)
  p(0:,0:,0:,1:)     =>  var(:,:,:,:,4)
  mu(0:,0:,0:,1:)    =>  var(:,:,:,:,5)
  nu(0:,0:,0:,1:)    =>  var(:,:,:,:,6)
  div_v(0:,0:,0:,1:) =>  var(:,:,:,:,7)

  var_name(1:7) = [ 'v_x  ', 'v_y  ', 'v_z  ', 'p    ' &
                  , 'mu   ', 'nu   ', 'div_v'          ]

  i = 7

  if (problem % HasExactSolution()) then

    u_ex(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,i+1:i+4)
    v_ex(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,i+1:i+3)
    p_ex(0:,0:,0:,1:)     =>  var(:,:,:,:,i+1)

    var_name(i+1:i+4) = [ 'v_x__exact', 'v_y__exact', 'v_z__exact', 'p__exact  ']

    i = i + 4

    err_u(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,i+1:i+4)
    err_v(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,i+1:i+3)
    err_p(0:,0:,0:,1:)     =>  var(:,:,:,:,i+1)

    var_name(i+1:i+4) = [ 'error(v_x)', 'error(v_y)', 'error(v_z)', 'error(p)  ']

    i = i + 4

  end if

  if (avg_rate > 0) then
    q_avg(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,i+1:i+10)
    var_name(i+ 1) = '⟨v_x⟩'
    var_name(i+ 2) = '⟨v_y⟩'
    var_name(i+ 3) = '⟨v_z⟩'
    var_name(i+ 4) = '⟨p⟩'
    var_name(i+ 5) = '⟨v_x v_x⟩'
    var_name(i+ 6) = '⟨v_x v_y⟩'
    var_name(i+ 7) = '⟨v_x v_z⟩'
    var_name(i+ 8) = '⟨v_y v_y⟩'
    var_name(i+ 9) = '⟨v_y v_z⟩'
    var_name(i+10) = '⟨v_z v_z⟩'
  else
    q_avg => null()
  end if

  allocate(vp(0:po_u,0:po_u,1:6,1:n_elem,1:3))
  allocate(w (0:po_u,0:po_u,0:po_u,1:n_elem,1:4) )

  ! initial conditions .........................................................

  if (restart_in) then
    data_file = trim(flow_case) // '_' // trim(restart_tag_in) // '_data'
    call ReadRestartData(data_file, rank, t, u, n_avg, q_avg)
    call XMPI_Bcast(t, 0, comm)
  else
    t = 0
    n_avg = 0
    call SetArray(mu, ins_op % mu_0)
    call SetArray(nu, ins_op % nu_0)
    call problem % GetInitialValues(ins_op % sem_u % metrics % x, u)
    if (smooth_initial_data) then
      call SmoothMeshData_3D( ins_op%mesh, ins_op%eop_u, u &
                            , smooth_filter, smooth_order  )
    end if
  end if

  ! time scales
  call time_scales % Evaluate(ins_op, u, comm)

  ! info .......................................................................

  call ins_op % sem_u % Get_Volume(domain_volume)

  call XMPI_Reduce(n_elem, n_elem_tot, MPI_SUM, 0, comm)
  n_point = n_elem_tot * (po_u+1)**3

  if (rank == 0) then
    write(*,'(/,A)') 'problem and discretization parameters'
    write(*,'(T3,A,T30,9(G0,X))') 'flow problem:', trim(flow_problem)
    write(*,'(T3,A,T30,9(G0,X))') 'domain:',  domain_name
    write(*,'(T3,A,T30,9(G0,X))') 'domain volume:', domain_volume
    write(*,'(T3,A,T30,9(G0,X))') 'boundary conditions for v:', problem % bc_v
    write(*,'(T3,A,T30,9(G0,X))') 'boundary conditions for p:', problem % bc_p
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of v:', ins_op % eop_u % po
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of p:', ins_op % eop_p % po
    write(*,'(T3,A,T30,9(G0,X))') 'conv quadrature order:', ins_op % sop_q % po
    write(*,'(T3,A,T30,9(G0,X))') 'number of mesh points:', n_point
    write(*,'(T3,A,T30,9(G0,X))') 'time integrator:'      , trim(ins_ti % name)
    write(*,'(T3,A,T29,ES18.11)') 'time step size:'       , dt
    write(*,*)
    write(*,'(T3,A)') 'convective and diffusive CFL numbers'
    write(*,'(T5,A,T16,ES12.5)') 'C(v_0  ) =' , dt / time_scales % tau_conv_v
    write(*,'(T5,A,T16,ES12.5)') 'C(v_ref) =' , dt / time_scales % tau_conv_r
    write(*,'(T5,A,T17,ES12.5)') 'D(ν_ref) =' , dt / time_scales % tau_diff_r
    write(*,*)
  end if

  if (mesh_stat) then
    call PrintStatistics()
  end if

  !-----------------------------------------------------------------------------
  ! Time integration

  call flow_char % Evaluate(ins_op, t, u, dt, domain_volume)
  call flow_char % PrintHeader()
  call flow_char % PrintValues('#init#')

  if (rank == 0) then
    !$omp master
    time = MPI_Wtime()
    !$omp end master
  end if

  do nt = 1, nt_max
    last = t + dt >= t_end .or. nt == nt_max
    call ins_ti % TimeStep(t, dt, mu, nu, u, standby = .not. last)
    if (avg_rate > 0 .and. mod(nt, max(avg_rate,1)) == 0) then
      call TemporalAveraging(u, q_avg, n_avg)
    end if
    if (last) then
      call flow_char % Evaluate(ins_op, t, u, dt, domain_volume)
      call flow_char % PrintValues('#last#')
      exit
    else if (mod(nt, char_freq) == 0) then
      call flow_char % Evaluate(ins_op, t, u, dt, domain_volume)
      call flow_char % PrintValues()
    end if
  end do

  if (rank == 0) then
    !$omp master
    time = MPI_Wtime() - time
    if (nt_max > 0) then
      nt = max(min(nt, nt_max),1)
      write(*,'(/,A)') 'performance'
      write(*,'(A,ES10.3)') ' time               =', time
      write(*,'(A,ES10.3)') ' time       / step  =', time / nt
      write(*,'(A,ES10.3)') ' throughput / step  =', nt / time * size(u) * n_proc
    end if
    !$omp end master
  end if

  !-----------------------------------------------------------------------------
  ! evaluation

  ! errors .....................................................................

  if (problem % HasExactSolution() .and. ins_op % mesh % part >= 0) then

    call problem % GetExactSolution(ins_op % sem_u % metrics % x, t, u_ex)

    call SetArray(err_u, u, multi=.true.)
    call MergeArrays(ONE, err_u, -ONE, u_ex, multi=.true.)
    call CalibrateArray(err_p, comm = ins_op % mesh % comm_parts)

  end if

  ! boundary fluxes ............................................................

  if (ins_op % mesh % part >= 0) then

    !$omp master
    allocate(bv_vn(n_bound), int_vn(1,n_bound))
    !$omp end master
    !$omp barrier

    do i = 1, n_bound
      call bv_vn(i) % Init(ins_op % mesh % boundary(i), po_u, nc = 1)
      call bv_vn(i) % ExtractNormalComponent(ins_op % sem_u, v)
    end do

    call GetSurfaceIntegrals(ins_op % sem_u, bv_vn, int_vn)

  end if

  if (ins_op % mesh % part == 0) then
    !$omp master
    write(*,'(/,A)') 'boundary fluxes'
    do i = 1, n_bound
       write(*,'(T3,A,I3,A,T29,ES12.5)') &
           'boundary',i,',   int(vn)  =', int_vn(1,i)
    end do
    write(*,'(T3,A,T29,ES12.5)') 'total:     sum(int(vn)) =', sum(int_vn)
    write(*,*)
    !$omp end master
  end if

  !-----------------------------------------------------------------------------
  ! Write plot files

  if (export_vtk .and. ins_op % mesh%part >= 0) then

    ! compute divergence
    call GetOuterVectorTraces_3D(ins_op%mesh, v, vp)   ! vp = v⁺
    call TPO_Div(ins_op%eop_u, ins_op%sem_u, v, vp, div_v)

    call ExportVTK_VolumeData( x       = ins_op % sem_u % metrics % x &
                             , s       = var                          &
                             , sname   = var_name                     &
                             , file    = flow_case                    &
                             , part    = ins_op % mesh % part         &
                             , n_parts = ins_op % mesh % n_parts      &
                             , subdiv  = subdiv_vtk                   )
  end if

  !-----------------------------------------------------------------------------
  ! Write restart data

  if (restart_out) then
    mesh_file = trim(flow_case) // '_' // trim(restart_tag_out) // '_mesh'
    data_file = trim(flow_case) // '_' // trim(restart_tag_out) // '_data'
    ! mesh
    call ml_mesh % WriteHDF5(mesh_file)
    ! flow data
    call WriteRestartData(data_file, rank, t, u, n_avg, q_avg)
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !> Write restart data

  subroutine WriteRestartData(file, rank, t, u, n_avg, q_avg)
    use, intrinsic ::  ISO_C_Binding
    use Kind_Parameters
    use HDF5_Binding
    implicit none

    character(len=*),  intent(in) :: file                !< file base name
    integer,           intent(in) :: rank                !< process rank
    real(RNP), target, intent(in) :: t                   !< problem time
    real(RNP), target, intent(in) :: u(0:,0:,0:,:,:)     !< variables
    integer,   target, intent(in) :: n_avg               !< averaging counter
    real(RNP), target, intent(in) :: q_avg(0:,0:,0:,:,:) !< averaged quantities

    optional :: q_avg

    integer(HID_T)    :: data_id, file_id, group_id, space_id
    integer(HSIZE_T)  :: dims_t(1), dims_n(1), dims_var(5)
    integer           :: err
    character(len=80) :: tag
    character(len=:), allocatable :: file_pr

    call Init_HDF5_Binding()

    ! append process rank to file name
    write(tag,'(I0)') rank
    file_pr = trim(file)//'_'//trim(tag)//'.h5'

    ! create HDF5 file and group
    call H5Fcreate_f(file_pr, H5F_ACC_TRUNC_F, file_id, err)
    call H5Gcreate_f(file_id, 'data', group_id, err)

    ! time
    dims_t = 1
    call H5Screate_simple_f(size(dims_t), dims_t, space_id, err)
    call H5Dcreate_f(group_id, 't', H5T_REAL_RNP, space_id, data_id, err)
    call H5Dwrite_f(data_id, H5T_REAL_RNP, C_Loc(t), err)
    call H5Sclose_f(space_id, err)
    call H5Dclose_f(data_id, err)

    ! u
    dims_var = shape(u)
    call H5Screate_simple_f(size(dims_var), dims_var, space_id, err)
    call H5Dcreate_f(group_id, 'u', H5T_REAL_RNP, space_id, data_id, err)
    call H5Dwrite_f(data_id, H5T_REAL_RNP, C_Loc(u), err)
    call H5Sclose_f(space_id, err)
    call H5Dclose_f(data_id, err)

    ! n_avg
    dims_n = 1
    call H5Screate_simple_f(size(dims_n), dims_n, space_id, err)
    call H5Dcreate_f(group_id, 'n_avg', H5T_INTEGER, space_id, data_id, err)
    call H5Dwrite_f(data_id, H5T_INTEGER, C_Loc(n_avg), err)
    call H5Sclose_f(space_id, err)
    call H5Dclose_f(data_id, err)

    ! q_avg
    if (present(q_avg)) then
      dims_var = shape(q_avg)
      call H5Screate_simple_f(size(dims_var), dims_var, space_id, err)
      call H5Dcreate_f(group_id, 'q_avg', H5T_REAL_RNP, space_id, data_id, err)
      call H5Dwrite_f(data_id, H5T_REAL_RNP, C_Loc(q_avg), err)
      call H5Sclose_f(space_id, err)
      call H5Dclose_f(data_id, err)
    end if

    ! release resources
    call H5Gclose_f(group_id, err)
    call H5Fclose_f(file_id, err)

  end subroutine WriteRestartData

  !-----------------------------------------------------------------------------
  !> Read restart data

  subroutine ReadRestartData(file, rank, t, u, n_avg, q_avg)
    use, intrinsic ::  ISO_C_Binding
    use Kind_Parameters
    use Execution_Control
    use HDF5_Binding
    implicit none

    character(len=*),  intent(in)    :: file                !< file base name
    integer,           intent(in)    :: rank                !< process rank
    real(RNP), target, intent(inout) :: t                   !< problem time
    real(RNP), target, intent(inout) :: u(0:,0:,0:,:,:)     !< variables
    integer,   target, intent(inout) :: n_avg               !< averaging counter
    real(RNP), target, intent(inout) :: q_avg(0:,0:,0:,:,:) !< averaged quantities

    optional :: q_avg

    integer(HID_T)    :: data_id, file_id, group_id, space_id, type_id
    integer(HSIZE_T)  :: dims_var(5), maxdims_var(5)
    integer           :: err
    logical           :: exists
    type(C_Ptr)       :: buf
    character(len=80) :: tag
    character(len=:), allocatable :: file_pr

    call Init_HDF5_Binding()

    ! append process rank to file name and check if it exists
    write(tag,'(I0)') rank
    file_pr = trim(file)//'_'//trim(tag)//'.h5'
    inquire(file=file_pr, exist=exists)

    if (.not. exists) return

    ! open HDF5 file and group for reading
    call H5Fopen_f(file_pr, H5F_ACC_RDWR_F, file_id, err)
    call H5Gopen_f(file_id, 'data', group_id, err)

    ! time
    buf = C_Loc(t)
    call H5Dopen_f(group_id, 't', data_id, err)
    call H5Dget_type_f(data_id, type_id, err)
    call H5Dread_f(data_id, type_id, buf, err)
    call H5Tclose_f(type_id, err)
    call H5Dclose_f(data_id, err)

    ! u
    buf = C_Loc(u)
    call H5Dopen_f(group_id, 'u', data_id, err)
    call H5Dget_space_f(data_id, space_id, err)
    call H5Sget_simple_extent_dims_f(space_id, dims_var, maxdims_var, err)
    if (any(dims_var /= shape(u))) then
      call Error('ReadRestartData','shape(u) not matching')
    end if
    call H5Dget_type_f(data_id, type_id, err)
    call H5Dread_f(data_id, type_id, buf, err)
    call H5Tclose_f(type_id, err)
    call H5Dclose_f(data_id, err)

    ! n_avg
    buf = C_Loc(n_avg)
    call H5Dopen_f(group_id, 'n_avg', data_id, err)
    call H5Dget_type_f(data_id, type_id, err)
    call H5Dread_f(data_id, type_id, buf, err)
    call H5Tclose_f(type_id, err)
    call H5Dclose_f(data_id, err)

    ! q_avg
    if (n_avg > 0 .and. present(q_avg)) then
      buf = C_Loc(q_avg)
      call H5Dopen_f(group_id, 'q_avg', data_id, err)
      call H5Dget_space_f(data_id, space_id, err)
      call H5Sget_simple_extent_dims_f(space_id, dims_var, maxdims_var, err)
      if (any(dims_var /= shape(q_avg))) then
        call Error('ReadRestartData','shape(q_avg) not matching')
      end if
      call H5Dget_type_f(data_id, type_id, err)
      call H5Dread_f(data_id, type_id, buf, err)
      call H5Tclose_f(type_id, err)
      call H5Dclose_f(data_id, err)
    end if

    ! release resources
    call H5Gclose_f(group_id, err)
    call H5Fclose_f(file_id, err)

  end subroutine ReadRestartData

  !-----------------------------------------------------------------------------
  !> Temporal averaging

  subroutine TemporalAveraging(u, q_avg, n_avg)
    real(RNP), contiguous, intent(in)    :: u(:,:,:,:,:)     !< sample
    real(RNP), contiguous, intent(inout) :: q_avg(:,:,:,:,:) !< avg quantities
    integer,               intent(inout) :: n_avg            !< counter

    real(RNP) :: wa, ws
    integer   :: e, i, j, k, ne, np

    n_avg = n_avg + 1
    if (n_avg > 1) then
      wa = real(n_avg - 1, RNP) / n_avg   ! weight of current average
      ws = ONE - wa                       ! weight of sample
    else
      wa = 0
      ws = 1
    end if

    np = size(u,1)
    ne = size(u,4)

    !$omp do
    do e = 1, ne
      do k = 1, np
      do j = 1, np
      do i = 1, np
        q_avg(i,j,k,e, 1) = wa * q_avg(i,j,k,e, 1) + ws * u(i,j,k,e,1)
        q_avg(i,j,k,e, 2) = wa * q_avg(i,j,k,e, 2) + ws * u(i,j,k,e,2)
        q_avg(i,j,k,e, 3) = wa * q_avg(i,j,k,e, 3) + ws * u(i,j,k,e,3)
        q_avg(i,j,k,e, 4) = wa * q_avg(i,j,k,e, 4) + ws * u(i,j,k,e,4)
        q_avg(i,j,k,e, 5) = wa * q_avg(i,j,k,e, 5) + ws * u(i,j,k,e,1) * u(i,j,k,e,1)
        q_avg(i,j,k,e, 6) = wa * q_avg(i,j,k,e, 6) + ws * u(i,j,k,e,1) * u(i,j,k,e,2)
        q_avg(i,j,k,e, 7) = wa * q_avg(i,j,k,e, 7) + ws * u(i,j,k,e,1) * u(i,j,k,e,3)
        q_avg(i,j,k,e, 8) = wa * q_avg(i,j,k,e, 8) + ws * u(i,j,k,e,2) * u(i,j,k,e,2)
        q_avg(i,j,k,e, 9) = wa * q_avg(i,j,k,e, 9) + ws * u(i,j,k,e,2) * u(i,j,k,e,3)
        q_avg(i,j,k,e,10) = wa * q_avg(i,j,k,e,10) + ws * u(i,j,k,e,3) * u(i,j,k,e,3)
      end do
      end do
      end do
    end do

  end subroutine TemporalAveraging

  !-----------------------------------------------------------------------------
  !> Mesh statistics

  subroutine PrintStatistics()

    integer, allocatable, save :: n_elem(:,:), n_leaf(:,:)

    integer(IXL) :: dof_p_leaf, dof_p_tot, dof_v_leaf
    integer      :: l

    call ml_op_p % Print_MeshCharacteristics()
    call ml_op_p % Get_MeshCharacteristics(n_elem, n_leaf = n_leaf)

    if (rank == 0) then

      dof_p_tot  = 0

      do l = 1, l_top
        dof_p_tot = dof_p_tot + n_elem(l,4) * (po_p(l) + 1)**3
      end do

      dof_p_leaf = int(n_leaf(l_top,4), IXL) * (po_p(l_top) + 1)**3
      dof_v_leaf = int(n_leaf(l_top,4), IXL) * (po_u        + 1)**3 * 3

      write(*,*)
      write(*,'(4X,A)') 'flow solver'
      write(*,'(T7,A,T20,I0)') 'dof_p_tot  =', dof_p_tot
      write(*,'(T7,A,T20,I0)') 'dof_p_leaf =', dof_p_leaf
      write(*,'(T7,A,T20,I0)') 'dof_v_leaf =', dof_v_leaf
      write(*,*)

    end if

  end subroutine PrintStatistics

  !=============================================================================

end program INS_Integrator_3D_Test

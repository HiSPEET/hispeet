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

!> summary:  3D multilevel solver for incompressible Navier-Stokes problems
!> author:   Joerg Stiller
!> date:     2025/03/28
!>
!> If present, the first argument of the invoking command will be interpreted
!> as the base name of the control file. If omitted, the program looks for
!> `ml_ins_solver_3d.prm`.
!===============================================================================

program ML_INS_Solver_3D
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Array_Assignments
  use Execution_Control
  use Logging_Levels

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  use Mesh__3D
  use Generic_Mesh__3D
  use Verify_Mesh__3D
  use Import_GMSH__3D
  use Data_Exchange__3D

  use INS__Problem__3D
  use INS__Problem__Test_Suite__3D

  use ML__Mesh__3D
  use ML__Mesh_Variable__3D
  use ML__INS__Operator__3D
  use ML__INS__Integrator__BDF__3D
  use ML__INS__Flow_Characteristics__3D
  use ML__INS__Time_Scales__3D
  use ML__INS__Time_Averaging__3D

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

  character(len=*), parameter :: default_case = 'ml_ins_solver_3d'
  character(len=80) :: flow_case ! flow case name
  character(len=80) :: case_file ! flow case input file: trim(flow_case).prm

  character(len=80) :: flow_problem = 'VariableViscosity'  ! problem name
  character(len=80) :: problem_file = 'variable_viscosity' ! parameters file

  integer :: flow_domain = 1
    ! computational flow domain (u/s = un/structured, r = regular, d = deformed)
    !   1  cuboidal domain with Cartesian mesh                             (s+r)
    !   2  cuboidal domain with unstructured "diamond" mesh                (u+d)
    !   3  cylindrical domain                                              (u+d)
    !   4  annular domain                                                  (u+d)
    !  10  import from GMSH

  character(len=80) :: raw_mesh_file = ''

  namelist/control_prm/ flow_problem, problem_file, flow_domain, raw_mesh_file

  integer :: char_freq  = 1 ! characteristics output frequency
  integer :: avg_rate   = 0 ! sampling rate for averaging, 0 if none
  integer :: vtk_mode   = 0 ! VTK export mode, 0/1/2/3: none/all/active/leafs

  namelist/control_prm/ char_freq, avg_rate, vtk_mode

  ! restart options
  character(len=80) :: restart_tag_in  = ''  ! tag for restart input files
  character(len=80) :: restart_tag_out = ''  ! tag for restart output files
  namelist/control_prm/ restart_tag_in, restart_tag_out
    !
    ! restart input is read from
    !   - trim(flow_case)_trim(restart_tag_in)_restart.prm     for parameters
    !   - trim(flow_case)_trim(restart_tag_in)_mesh_<rank>.h5  for the mesh
    !   - trim(flow_case)_trim(restart_tag_in)_data_<rank>.h5  for flow data
    ! where <rank> is the process rank in mesh%comm_world
    !
    ! output written to
    !   - trim(flow_case)_trim(restart_tag_out)_restart.prm     for parameters
    !   - trim(flow_case)_trim(restart_tag_out)_mesh_<rank>.h5  for the mesh
    !   - trim(flow_case)_trim(restart_tag_out)_data_<rank>.h5  for flow data
    !
    ! no restart data is read or written if the corresponding tag is empty

  ! mesh handling
  logical :: restart_adjust = .false. ! adjust mesh after restart
  namelist/control_prm/ restart_adjust

  ! mesh .......................................................................

  type(GenericMesh_3D),       save :: generic_mesh
  type(Mesh_3D), allocatable, save :: base_mesh
  type(ML_Mesh_3D),           save :: ml_mesh

  type(ML_Mesh_Options_3D), save :: ml_mesh_opt
    ! multilevel mesh options, defining top level, partitioning and refinement,
    ! stacic components defined in namelist `ml_mesh_options_3d__static` and
    ! dynamic components in `ml_mesh_options_3d__dynamic`

  ! problem ....................................................................

  class(INS_Problem_3D), allocatable, save :: problem

  ! spatial operators ..........................................................

  integer, allocatable, save :: po(:)
  type(ML_INS_OperatorOptions_3D), save :: ml_ins_opt
  type(ML_INS_Operator_3D), save :: ml_ins

  namelist/spatial_prm/ po, ml_ins_opt

  ! time integration ...........................................................

  type(ML_INS_Integrator_BDF_Options_3D), save :: ml_bdf_opt
  type(ML_INS_Integrator_BDF_3D), save :: ml_bdf

  namelist/temporal_prm/ ml_bdf_opt

  real(RNP) :: t_end  = 0.25
  real(RNP) :: dt     = 1E-3
  integer   :: nt_max = 1

  namelist/temporal_prm/ t_end, dt, nt_max

  ! variables ..................................................................

  real(RNP) :: t ! problem time

  type(ML_MeshVariable_3D), save :: var   ! solution and averaged variables
  type(ML_MeshVariable_3D), save :: u     ! handle for solution
  type(ML_MeshVariable_3D), save :: mu    ! handle for bulk viscosity
  type(ML_MeshVariable_3D), save :: nu    ! handle for shear viscosity
  type(ML_MeshVariable_3D), save :: q_avg ! handle for averaged quantities
  type(ML_MeshVariable_3D), save :: vtk   ! variables exported to VTK

  ! re/start parameters ........................................................

  ! variables stored in the restart parameter file
  real(RNP) :: t_0 = 0      ! physical time at start
  integer   :: n_sample = 0 ! averaging weight

  namelist/restart_prm/ t_0, n_sample

  ! auxiliaries ................................................................

  type(DataExchangePlan_3D), allocatable :: x_plan(:)

  type(ML_INS_TimeScales_3D)          :: ml_time_scales
  type(ML_INS_FlowCharacteristics_3D) :: ml_flow_char

  character(:), allocatable :: domain_name
  character(:), allocatable :: restart_file ! restart parameter file
  character(:), allocatable :: mesh_file    ! multilevel mesh file
  character(:), allocatable :: data_file    ! multilevel data file

  character(len=20), allocatable :: var_name(:)
  character(len=20), allocatable :: vtk_name(:)

  real(RNP) :: domain_volume
  logical   :: exists, restart_in, restart_out
  logical   :: first, last
  logical   :: perform_average
  logical   :: print_flow_char
  integer   :: io, stat
  integer   :: l_top, l_max, n_bound
  integer   :: n_comp  ! number of solution components
  integer   :: n_aux   ! number of auxiliary variables
  integer   :: n_avg   ! number of averaged quantities
  integer   :: n_var   ! number of solution and averaged variables
  integer   :: n_vtk   ! number of VTK quantities
  integer   :: i, k, l, nt

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

  if (rank == 0) then
    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Multilevel Navier-Stokes solver for incompressible flow'
    write(*,*)
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:', n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'  , n_thread
    write(*,*)
  end if

  ! control parameters .........................................................

  if (rank == 0) then

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
      if (char_freq < 1) then
        ! disable intermediate control output
        char_freq = huge(1)
      end if
    else
       call Error( 'ML_INS_Solver_3D', &
                   'input file "' // trim(case_file) // '" not found' )
    end if

  end if

  ! globalize logging levels
  call XMPI_Bcast_LoggingLevels(0, comm)

  ! globalize control parameters
  call XMPI_Bcast( flow_case           , 0, comm)
  call XMPI_Bcast( case_file           , 0, comm)
  call XMPI_Bcast( flow_problem        , 0, comm)
  call XMPI_Bcast( problem_file        , 0, comm)
  call XMPI_Bcast( flow_domain         , 0, comm)
  call XMPI_Bcast( raw_mesh_file       , 0, comm)
  call XMPI_Bcast( char_freq           , 0, comm)
  call XMPI_Bcast( avg_rate            , 0, comm)
  call XMPI_Bcast( vtk_mode            , 0, comm)
  call XMPI_Bcast( restart_tag_in      , 0, comm)
  call XMPI_Bcast( restart_tag_out     , 0, comm)
  call XMPI_Bcast( restart_adjust      , 0, comm)

  ! restart switches
  restart_in  = len_trim(restart_tag_in)  > 0
  restart_out = len_trim(restart_tag_out) > 0

  ! multilevel mesh options ....................................................

  if (rank == 0) then
    ml_mesh_opt = ML_Mesh_Options_3D(io, n_proc)
  end if

  call ml_mesh_opt % Bcast(0, comm)

  l_max = max(ml_mesh_opt%l_top, ml_mesh_opt%l_max)

  ! spatial ....................................................................

  allocate(po(l_max), source = -1)

  if (rank == 0) then
    read(io, nml = spatial_prm)
  end if

  call XMPI_Bcast(po, 0, comm)
  call ml_ins_opt % Bcast(0, comm)

  ! temporal ...................................................................

  if (rank == 0) then
    read(io, nml = temporal_prm)
    close(io)
  end if

  call ml_bdf_opt % Bcast(0, comm)

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

    ! read restart parameters ..................................................

    if (rank == 0) then
      restart_file = trim(flow_case)//'_'//trim(restart_tag_in)//'_restart.prm'
      inquire(file=restart_file, exist=exists)
      if (exists) then
        write(*,'(2X,A)') 'reading restart parameters from '//trim(restart_file)
        open(newunit = io, file = restart_file)
        read(io, nml = restart_prm)
        close(io)
      else
         call Error( 'ML_INS_Solver_3D' &
                   , 'file "' // trim(restart_file) // '" not found' )
      end if
    end if

    call XMPI_Bcast(t_0     , 0, comm)
    call XMPI_Bcast(n_sample, 0, comm)

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
    end if

    ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)
    deallocate(base_mesh)

  end if

  l_top = size(ml_mesh % mesh)
  n_bound = ml_mesh % mesh(1) % n_bound

  ! problem ....................................................................

  call Set_INS_TestProblem_3D(problem, flow_problem, problem_file, n_bound, comm)

  ! check & fix boundary conditions
  do i = 1, n_bound
    if (ml_mesh % mesh(1) % boundary(i) % coupled > 0) then
      if (problem % bc_v(i) /= 'P') then
        if (rank == 0) then
          write(*,'(2X,A,X,I0)') '*** enforcing periodic BC on boundary', i
        end if
        problem % bc_v(i) = 'P'
        problem % bc_p(i) = 'P'
      end if
    end if
  end do

  ! multilevel Navier-Stokes operator ..........................................

  ml_ins = ML_INS_Operator_3D(ml_mesh, po, problem, ml_ins_opt)

  ! variables ..................................................................

  ! number of variables
  n_comp = problem % nc
  n_aux  = 2
  if (avg_rate > 0) then
    n_avg = 2 * n_comp  & ! ⟨uᵢ⟩, ⟨uᵢuᵢ⟩
          + 3             ! ⟨u₁u₂⟩, ⟨u₁u₃⟩, ⟨u₂u₃⟩
  else
    n_avg = 0
  end if
  n_var = n_comp + n_aux + n_avg

  allocate(character(len=20) :: var_name(n_var))

  ! names of solution components
  var_name(1:4) = [ 'v_x', 'v_y', 'v_z', 'p  ']
  do i = 5, n_comp
    write(var_name(i),'(A,I0)') 'u_',i
  end do

  ! auxiliary variables
  k = n_comp
  var_name(k+1:k+2) = [ 'mu', 'nu' ]

  ! names of averaged quantities
  if (avg_rate > 0) then
    k = n_comp + n_aux
    var_name(k+1) = '⟨v_x⟩'
    var_name(k+2) = '⟨v_y⟩'
    var_name(k+3) = '⟨v_z⟩'
    var_name(k+4) = '⟨p⟩'
    do i = 5, n_comp
      write(var_name(k+i),'(A,I0,A)') '⟨u_',i,'⟩'
    end do
    k = 2 * n_comp + n_aux
    k = k + 1;  var_name(k) = '⟨v_x v_x⟩'
    k = k + 1;  var_name(k) = '⟨v_x v_y⟩'
    k = k + 1;  var_name(k) = '⟨v_x v_z⟩'
    k = k + 1;  var_name(k) = '⟨v_y v_y⟩'
    k = k + 1;  var_name(k) = '⟨v_y v_z⟩'
    k = k + 1;  var_name(k) = '⟨v_z v_z⟩'
    k = k + 1;  var_name(k) = '⟨p p⟩'
    do i = 5, n_comp
      k = k + 1
      write(var_name(k),'(3(A,I0))') '⟨u_',i,', u_',i,'⟩'
    end do
  end if

  ! generate multilevel variables and handles
  call var % Init(ml_ins%ml_op_u, nc = n_var, name = var_name)
  k = n_comp
  call var % GetSlice(u , first =     1, last = k)
  call var % GetSlice(mu, first = k + 1, last = k + 1)
  call var % GetSlice(nu, first = k + 2, last = k + 2)
  k = n_comp + n_aux
  call var % GetSlice(q_avg, first = k + 1, last = k + n_avg)

  !-----------------------------------------------------------------------------
  ! Initial conditions and operators

  if (rank == 0) then
    write(*,'(/,A)') 'read/generate initial conditions and (re)create operators'
  end if

  if (restart_in) then

    ! read restart data
    data_file = trim(flow_case) // '_' // trim(restart_tag_in) // '_data'
    call var % ReadHDF5(data_file)

    ! adjust mesh, Navier-Stokes operator and variables
    if (restart_adjust) then
      call ml_mesh % MarkByOptions(ml_mesh_opt)
      call ml_mesh % Adapt(ml_mesh_opt%partition, x_plan)
      ml_ins = ML_INS_Operator_3D(ml_mesh, po, problem, ml_ins_opt)
      call var % FitAdapt(ml_ins % ml_op_u, x_plan)
      k = n_comp
      call var % GetSlice(u , first =     1, last = k)
      call var % GetSlice(mu, first = k + 1, last = k + 1)
      call var % GetSlice(nu, first = k + 2, last = k + 2)
      k = n_comp + n_aux
      call var % GetSlice(q_avg, first = k + 1, last = k + n_avg)
    end if
    l_top = size(ml_mesh % mesh)

  else

    ! set initial conditions
    do l = 1, l_top
      associate(ins_l => ml_ins % ins_op(l), u_l => u % level(l) % val)
        call problem % GetInitialValues(ins_l % sem_u % metrics % x, u_l)
      end associate
    end do

  end if

  ! temporal operators .........................................................

  ml_bdf = ML_INS_Integrator_BDF_3D(problem, ml_ins, ml_bdf_opt)

  !-----------------------------------------------------------------------------
  ! Print info

  ! mesh
  if (rank == 0) then
    write(*,'(/,A)') 'multilevel mesh characteristics'
  end if
  call ml_ins % ml_op_u % Print_MeshCharacteristics()

  ! time scales
  call ml_time_scales % Evaluate(ml_ins, u)
  if (rank == 0) then
    write(*,*)
    write(*,'(T1,A)') 'convective and diffusive CFL numbers'
    write(*,'(T3,A,T14,ES12.5)') 'C(v_0  ) =' , &
        dt / minval(ml_time_scales % level % tau_conv_v)
    write(*,'(T3,A,T14,ES12.5)') 'C(v_ref) =' , &
        dt / minval(ml_time_scales % level % tau_conv_r)
    write(*,'(T3,A,T15,ES12.5)') 'D(ν_ref) =' , &
        dt / minval(ml_time_scales % level % tau_diff_r)
    write(*,*)
  end if

  !-----------------------------------------------------------------------------
  ! Time integration

  t = t_0

  call XMPI_Bcast(dt    , 0, comm)
  call XMPI_Bcast(t_end , 0, comm)
  call XMPI_Bcast(nt_max, 0, comm)

  !  domain volume
  call ml_ins % ml_op_u % Get_Volume(domain_volume)

  ! initial flow characteristics
  call ml_flow_char % Evaluate(ml_ins, t, u, dt, domain_volume, leaf = .true.)
  call ml_flow_char % PrintHeader()
  call ml_flow_char % PrintValues('#init#')

  do nt = 1, nt_max
    first = nt == 1
    last  = t + dt >= t_end .or. nt == nt_max

    call ml_bdf % TimeStep(t, dt, u, first, last)

    perform_average = avg_rate > 0 .and. mod(nt, max(avg_rate,1)) == 0
    print_flow_char = mod(nt, char_freq) == 0 .or. last

    if (perform_average .or. print_flow_char) then
      call ml_ins % CalibratePressure(u)
    end if

    if (perform_average) then
      call ML_INS_TimeAveraging_3D(u, q_avg, n_avg)
    end if

    if (print_flow_char) then
      call ml_flow_char % Evaluate(ml_ins, t, u, dt, domain_volume, leaf=.true.)
      if (last) then
        call ml_flow_char % PrintValues('#last#')
      else
        call ml_flow_char % PrintValues()
      end if
    end if

    if (last) exit
  end do

  if (rank == 0) then
    write(*,*)
  end if

  !-----------------------------------------------------------------------------
  ! Write restart files

  if (restart_out) then

    restart_file = trim(flow_case)//'_'//trim(restart_tag_out) // '_restart.prm'
    mesh_file    = trim(flow_case)//'_'//trim(restart_tag_out) // '_mesh'
    data_file    = trim(flow_case)//'_'//trim(restart_tag_out) // '_data'

    if (rank == 0) then
      write(*,'(A)') 'writing restart data'
      open(newunit = io, file = restart_file)
      t_0 = t
      write(io, nml = restart_prm)
      close(io)
    end if

    call ml_mesh % WriteHDF5(mesh_file)
    call var % WriteHDF5(data_file)

  end if

  !-----------------------------------------------------------------------------
  ! Write plot files

  EXPORT_VTK: if (vtk_mode > 0) then

    if (rank == 0) then
      write(*,'(A)') 'writing VTK file'
    end if

    if (problem % HasExactSolution()) then

      ! export data including exact solution and error .........................

      n_vtk = n_var + 2 * n_comp

      allocate(vtk_name(n_vtk))
      vtk_name(1:n_var) = var_name
      k = n_var
      do i = 1, n_comp
        write(vtk_name(k + i         ),'(2A)') trim(var_name(i)), '__exact'
        write(vtk_name(k + i + n_comp),'(2A)') trim(var_name(i)), '__error'
      end do

      call vtk % Init(ml_ins%ml_op_u, n_vtk, vtk_name)

      do l = 1, l_top
        call SetArray(vtk%level(l)%val(:,:,:,:,1:n_var), var%level(l)%val)
        associate( x_l => ml_ins % ml_op_u % sem(l) % metrics % x           &
                 , u_l => vtk % level(l) % val(:,:,:,:,1+0*n_comp:1*n_comp) &
                 , s_l => vtk % level(l) % val(:,:,:,:,1+1*n_comp:2*n_comp) &
                 , e_l => vtk % level(l) % val(:,:,:,:,1+2*n_comp:3*n_comp) )

          call problem % GetExactSolution(x_l, t, s_l)
          e_l = u_l - s_l
        end associate
      end do

      call vtk % ExportVTK(ml_ins%ml_op_u, file = flow_case, mode = vtk_mode)

    else

      ! export solution and averaged variables .................................

      call var % ExportVTK(ml_ins%ml_op_u, file = flow_case, mode = vtk_mode)

    end if

  end if EXPORT_VTK

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

  !=============================================================================

end program ML_INS_Solver_3D

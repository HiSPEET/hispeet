!> summary:  Validation of incompressible Navier-Stokes time integrators
!> author:   Joerg Stiller
!> date:     2022/10/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> If present, the first argument of the invoking command will be interpreted
!> as the base name of the control file. If omitted, the program looks for
!> `ins_timeintegrator_3d_test.prm`.
!===============================================================================

program INS_TimeIntegrator_3D_Test
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Array_Assignments
  use Array_Reductions

  use Export_VTK_Volume_Data__3D
  use Mesh__3D

  use INS__Problem__3D
  use INS__Problem__Vortex_TG__3D
  use INS__Problem__Variable_Viscosity__3D
  use INS__Operator__3D
  use INS__Time_Integrator__3D
  use INS__Time_Integrator__Euler__3D

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm      ! MPI communicator
  integer        :: rank      ! local MPI rank
  integer        :: n_proc    ! number of MPI processes
  integer        :: n_thread  ! number of OpenMP threads

  ! control parameters .........................................................

  character(len=*), parameter :: default_case = 'ins_timeintegrator_3d_test'
  character(len=80) :: flow_case
  ! input file (*.prm)

  character(len=80) :: flow_problem = 'Vortex_TG'
  ! flow problem:
  !   'Vortex_TG'          2D Taylor-Green vortex
  !   'VariableViscosity'  3D Vortex array with variable viscosity

  character(len=80) :: problem_file = 'vortex_tg'
  ! file containing the problem parameters (*.prm)

  integer :: flow_domain = 1
  ! computational flow domain (u/s = un/structured, r = regular, d = deformed)
  !   1  cuboidal domain with Cartesian mesh                               (s+r)
  !   2  cuboidal domain with unstructured "diamond" mesh                  (u+d)
  !   3  cylindrical domain                                                (u+d)
  !   4  annular domain                                                    (u+d)

  namelist/control_prm/ flow_problem, problem_file, flow_domain

  type(INS_OperatorOptions_3D) :: ins_op_opts
  ! options for the incompressible Navier-Stokes operator

  type(INS_TimeIntegrator_Euler_Options_3D) :: ins_ti_euler_opts
  ! options for the Euler time integrator

  namelist/control_prm/ ins_op_opts, ins_ti_euler_opts

  real(RNP) :: t_end      = 1        ! final time
  real(RNP) :: dt         = 1        ! time step size
  integer   :: nt_max     = 0        ! max num time steps
  logical   :: export_vtk = .false.  ! generate VTK files

  namelist/control_prm/ t_end, dt, nt_max, export_vtk

  ! operators and variables ....................................................

  type(Mesh_3D), save :: mesh

  class(INS_Problem_3D), allocatable, save :: problem
  ! flow problem

  type(INS_Operator_3D), save :: ins_op
  ! incompressible Navier-Stokes operator

  class(INS_TimeIntegrator_3D), allocatable, save :: ins_ti
  ! incompressible Navier-Stokes time integrator

  real(RNP) :: t = 0 ! problem time

  real(RNP), allocatable, target, save :: var(:,:,:,:,:)
  character(len=20), allocatable, save :: var_name(:)

  real(RNP), pointer, contiguous, save :: u(:,:,:,:,:)    ! u = [v, p]
  real(RNP), pointer, contiguous, save :: v(:,:,:,:,:)    ! velocity
  real(RNP), pointer, contiguous, save :: p(:,:,:,:)      ! pressure

  real(RNP), pointer, contiguous, save :: u_ex(:,:,:,:,:) ! u_ex = [v_ex, p_ex]
  real(RNP), pointer, contiguous, save :: v_ex(:,:,:,:,:) ! exact velocity
  real(RNP), pointer, contiguous, save :: p_ex(:,:,:,:)   ! exact pressure

   real(RNP), allocatable, save :: w(:,:,:,:,:)   ! workspace

  ! auxiliaries ................................................................

  character(len=80) :: domain_name = ''
! real(RDP) :: time, time0
  real(RNP) :: e_p, e_v
  logical   :: exists
  integer   :: io, stat
  integer   :: n_bound, n_elem, n_elem_tot, n_ghost, n_point, n_var, po
  integer   :: i, nt

  !-----------------------------------------------------------------------------
  ! Initialization

  ! MPI and OpenMP .............................................................

  call XMPI_Init()
!### CHECK
print *, '$ 00'
!### CHECK END

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
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:'  , n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'    , n_thread
    write(*,*)

    call get_command_argument(1, flow_case, status=stat)
    if (stat /= 0 .or. len_trim(flow_case) == 0) then
      flow_case = default_case
    end if
    flow_case = trim(flow_case) // '.prm'

    inquire(file=trim(flow_case), exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(flow_case)
      open(newunit = io, file = flow_case)
      read(io, nml = control_prm)
      close(io)
    else
       call Error( 'INS_Operator_3D_Test', &
                   'input file "' // trim(flow_case) // '" not found' )
    end if

  end if

  ! globalize control parameters
  call XMPI_Bcast(flow_case   , 0, comm)
  call XMPI_Bcast(flow_problem, 0, comm)
  call XMPI_Bcast(flow_domain , 0, comm)
  call XMPI_Bcast(problem_file, 0, comm)
  call XMPI_Bcast(t_end       , 0, comm)
  call XMPI_Bcast(dt          , 0, comm)
  call XMPI_Bcast(nt_max      , 0, comm)
  call XMPI_Bcast(export_vtk  , 0, comm)

  ! globalize options
  call ins_op_opts       % Bcast(0, comm)
  call ins_ti_euler_opts % Bcast(0, comm)

  ! mesh .......................................................................

  select case(flow_domain)
  case(2)
    call CreateCuboidDiamonds(comm, flow_case, mesh)
    domain_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCylinder(comm, flow_case, mesh)
    domain_name = 'Cylindrical domain with unstructured mesh'
  case(4)
    call CreateAnnulus(comm, flow_case, mesh)
    domain_name = 'Annular domain with unstructured mesh'
  case default
    call CreateCuboidCartesian(comm, flow_case, mesh)
    domain_name = 'Cuboidal domain with Cartesian mesh'
  end select

  n_elem  = mesh % n_elem
  n_ghost = mesh % n_ghost
  n_bound = mesh % n_bound
!### CHECK
print *, '$ 01'
!### CHECK END

  ! problem ....................................................................

  select case(flow_problem)
  case('Vortex_TG')
    allocate(INS_Problem_Vortex_TG_3D         :: problem)
  case default
    allocate(INS_Problem_VariableViscosity_3D :: problem)
  end select

  call problem % SetProblem(n_bound, problem_file, comm)

  ! check & fix boundary conditions
  do i = 1, n_bound
    if (mesh % boundary(i) % coupled > 0) then
      if (problem % bc_v(i) /= 'P') then
        if (rank == 0) then
          write(*,'(A,I0)') '  *** enforcing periodic BC on coupled boundary ',i
        end if
        problem % bc_v(i) = 'P'
      end if
    end if
  end do

  ! operators ..................................................................

  ! Navier-Stokes operator, mesh will be copied
  ins_op = INS_Operator_3D(ins_op_opts, problem, mesh)
!### CHECK
print *, '$ 02'
print *, '$ 02, allocated(mesh % boundary)                  =', allocated(mesh % boundary)
print *, '$ 02, allocated(ins_op % mesh % boundary)         =', allocated(ins_op % mesh % boundary)
print *, '$ 02, allocated(ins_op % sem_v % mesh % boundary) =', allocated(ins_op % sem_v % mesh % boundary)
print *, '$ 02, allocated(ins_op % sem_p % mesh % boundary) =', allocated(ins_op % sem_p % mesh % boundary)
!### CHECK END

  ! time integrator: up to now only Euler
  ins_ti = INS_TimeIntegrator_Euler_3D(problem, ins_op, ins_ti_euler_opts)
!### CHECK
print *, '$ 03'
print *, '$ 03, shape(mesh % boundary)                   =', shape(mesh % boundary)
print *, '$ 03, shape(ins_op % mesh % boundary)          =', shape(ins_op % mesh % boundary)
print *, '$ 03, shape(ins_op % sem_v % mesh % boundary)  =', shape(ins_op % sem_v % mesh % boundary)
print *, '$ 03, shape(ins_ti % ins_op % mesh % boundary) =', shape(ins_ti % ins_op % mesh % boundary)
print *, '$ 03, shape(ins_ti % ins_op % sem_v % mesh % boundary) =', shape(ins_ti % ins_op % sem_p % mesh % boundary)
!### CHECK END

  ! variables ..................................................................

  po = ins_op % eop_v % po
  n_var = 8

  allocate(var(0:po,0:po,0:po,1:n_elem,1:n_var), source = ZERO)
  allocate(var_name(1:n_var))

  u(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,1:4)
  v(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,1:3)
  p(0:,0:,0:,1:)     =>  var(:,:,:,:,4)

  var_name(1:4) = [ 'v_x', 'v_y', 'v_z', 'p  ']

  u_ex(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,5:8)
  v_ex(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,5:7)
  p_ex(0:,0:,0:,1:)     =>  var(:,:,:,:,8)

  var_name(5:8) = [ 'v_x__exact', 'v_y__exact', 'v_z__exact', 'p__exact  ']

  allocate(w  (0:po,0:po,0:po,1:n_elem,1:4) )

  ! info .......................................................................

  call XMPI_Reduce(n_elem, n_elem_tot, MPI_SUM, 0, comm)
  n_point = n_elem_tot * (po+1)**3

  if (rank == 0) then
    write(*,'(/,A)') 'problem and discretization parameters'
    write(*,'(T3,A,T30,9(G0,X))') 'flow problem:', trim(flow_problem)
    write(*,'(T3,A,T30,9(G0,X))') 'domain:',  trim(domain_name)
    write(*,'(T3,A,T30,9(G0,X))') 'boundary conditions:'  , problem % bc_v
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of v:', ins_op % eop_v % po
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of p:', ins_op % eop_p % po
    write(*,'(T3,A,T30,9(G0,X))') 'conv quadrature order:', ins_op % sop_q % po
    write(*,'(T3,A,T30,9(G0,X))') 'conv quadrature type:' , ins_op % sop_q % basis
    write(*,'(T3,A,T30,9(G0,X))') 'number of mesh points:', n_point
    write(*,'(T3,A,T30,9(G0,X))') 'time step size:'       , dt
  end if
!### CHECK
print *, '$ 04'
!### CHECK END

  !-----------------------------------------------------------------------------
  ! Time integration

  ! initial conditions

  call problem % GetExactSolution(ins_op % sem_v % metrics % x, t, u)
!### CHECK
print *, '$ 05'
!### CHECK END

  do nt = 1, nt_max
    call ins_ti % TimeStep(t, dt, u)
    if (t >= t_end) exit
  end do
  nt = min(nt, nt_max)
!### CHECK
print *, '$ 06'
!### CHECK END

  !-----------------------------------------------------------------------------
  ! evaluation

  if (mesh % part >= 0) then

    call problem % GetExactSolution(ins_op % sem_v % metrics % x, t, u_ex)

    ! calibrate pressure to zero mean value
    call CalibrateArray(p   , comm=mesh%comm_parts)
    call CalibrateArray(p_ex, comm=mesh%comm_parts)

    ! w = u - u_ex,
    call SetArray(w, u, multi=.true.)
    call MergeArrays(ONE, w, -ONE, u_ex, multi=.true.)

    e_v = ScalarProduct(w(:,:,:,:,1:3), w(:,:,:,:,1:3), mesh%comm_parts)
    e_v = sqrt(e_v / n_point)

    e_p = ScalarProduct(w(:,:,:,:,4), w(:,:,:,:,4), mesh%comm_parts)
    e_p = sqrt(e_p / n_point)

  end if

  if (mesh % part == 0) then
    write(*,'(/,A)') 'time integration'
    write(*,'(T3,A,T29,ES12.5)') 'final time          t    =', t
    write(*,'(T3,A,T30,ES12.5)') 'velocity error      ε_v  =', e_v
    write(*,'(T3,A,T30,ES12.5)') 'diffusion           ε_p  =', e_p
    write(*,*)
  end if

  !-----------------------------------------------------------------------------
  ! Write plot files

  if (export_vtk .and. mesh%part >= 0) then
    call ExportVTK_VolumeData( x       = ins_op % sem_v % metrics % x &
                             , s       = var                          &
                             , sname   = var_name                     &
                             , file    = flow_case                    &
                             , part    = mesh % part                  &
                             , n_parts = mesh % n_parts               )
  end if

  !-----------------------------------------------------------------------------
  ! Finalization
!### CHECK
print *, '$ XX'
!### CHECK END

  call MPI_Finalize()

  !=============================================================================

end program INS_TimeIntegrator_3D_Test

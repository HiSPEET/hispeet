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

  use TPO__AAA__3D
  use TPO__Div__3D

  use Mesh__3D
  use Trace_Operators__3D
  use Volume_Integrals__3D
  use Export_VTK_Volume_Data__3D

  use INS__Problem__3D
  use INS__Problem__Vortex_TG__3D
  use INS__Problem__Variable_Viscosity__3D
  use INS__Operator__3D
  use INS__Time_Integrator__3D
  use INS__Time_Integrator__Euler__3D
  use INS__Time_Integrator__BDF2__3D

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
  character(len=80) :: case_file
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

  integer :: time_method = 1
  ! 1  Euler
  ! 2  BDF2

  namelist/control_prm/ time_method

  type(INS_OperatorOptions_3D) :: ins_op_opts
  type(INS_TimeIntegrator_Euler_Options_3D) :: ins_ti_euler_opts
  type(INS_TimeIntegrator_BDF2_Options_3D)  :: ins_ti_bdf2_opts

  namelist/control_prm/ ins_op_opts, ins_ti_euler_opts, ins_ti_bdf2_opts

  real(RNP) :: t_end      = 1        ! final time
  real(RNP) :: dt         = 1        ! time step size
  integer   :: nt_max     = 0        ! max num time steps
  logical   :: export_vtk = .false.  ! generate VTK files

  namelist/control_prm/ t_end, dt, nt_max, export_vtk

  ! operators and variables ....................................................

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

  real(RNP), allocatable, save :: w (:,:,:,:,:)   ! workspace

  ! auxiliaries ................................................................

  character(len=80) :: domain_name = ''
! real(RDP) :: time, time0
  real(RNP) :: domain_volume
  real(RNP) :: e_p, e_v, e_div_vv, e_div_vp, e_div_vq
  logical   :: exists, last
  integer   :: io, stat
  integer   :: n_bound, n_elem, n_elem_tot, n_ghost, n_point, n_var, po
  integer   :: i, nt

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
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:'  , n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'    , n_thread
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
    else
       call Error( 'INS_Operator_3D_Test', &
                   'input file "' // trim(case_file) // '" not found' )
    end if

  end if

  ! globalize control parameters
  call XMPI_Bcast(flow_case   , 0, comm)
  call XMPI_Bcast(case_file   , 0, comm)
  call XMPI_Bcast(flow_problem, 0, comm)
  call XMPI_Bcast(problem_file, 0, comm)
  call XMPI_Bcast(flow_domain , 0, comm)
  call XMPI_Bcast(time_method , 0, comm)
  call XMPI_Bcast(t_end       , 0, comm)
  call XMPI_Bcast(dt          , 0, comm)
  call XMPI_Bcast(nt_max      , 0, comm)
  call XMPI_Bcast(export_vtk  , 0, comm)

  ! globalize options
  call ins_op_opts       % Bcast(0, comm)
  call ins_ti_euler_opts % Bcast(0, comm)
  call ins_ti_bdf2_opts  % Bcast(0, comm)

  ! mesh .......................................................................

  select case(flow_domain)
  case(2)
    call CreateCuboidDiamonds(comm, case_file, ins_op % mesh)
    domain_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCylinder(comm, case_file, ins_op % mesh)
    domain_name = 'Cylindrical domain with unstructured mesh'
  case(4)
    call CreateAnnulus(comm, case_file, ins_op % mesh)
    domain_name = 'Annular domain with unstructured mesh'
  case default
    call CreateCuboidCartesian(comm, case_file, ins_op % mesh)
    domain_name = 'Cuboidal domain with Cartesian mesh'
  end select

  n_elem  = ins_op % mesh % n_elem
  n_ghost = ins_op % mesh % n_ghost
  n_bound = ins_op % mesh % n_bound

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
    if (ins_op % mesh % boundary(i) % coupled > 0) then
      if (problem % bc_v(i) /= 'P') then
        if (rank == 0) then
          write(*,'(A,I0)') '  *** enforcing periodic BC on coupled boundary ',i
        end if
        problem % bc_v(i) = 'P'
      end if
    end if
  end do

  ! operators ..................................................................

  ! Navier-Stokes operator
  call ins_op % Init(ins_op_opts, problem)

  ! time integrator: up to now only Euler
  select case(time_method)
  case(1)
    ins_ti = INS_TimeIntegrator_Euler_3D(problem, ins_op, ins_ti_euler_opts)
  case(2)
    ins_ti = INS_TimeIntegrator_BDF2_3D(problem, ins_op, ins_ti_bdf2_opts)
  end select

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

  call ins_op % sem_v % Get_Volume(domain_volume)

  call XMPI_Reduce(n_elem, n_elem_tot, MPI_SUM, 0, comm)
  n_point = n_elem_tot * (po+1)**3

  if (rank == 0) then
    write(*,'(/,A)') 'problem and discretization parameters'
    write(*,'(T3,A,T30,9(G0,X))') 'flow problem:', trim(flow_problem)
    write(*,'(T3,A,T30,9(G0,X))') 'domain:',  trim(domain_name)
    write(*,'(T3,A,T30,9(G0,X))') 'domain volume:', domain_volume
    write(*,'(T3,A,T30,9(G0,X))') 'boundary conditions:'  , problem % bc_v
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of v:', ins_op % eop_v % po
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of p:', ins_op % eop_p % po
    write(*,'(T3,A,T30,9(G0,X))') 'conv quadrature order:', ins_op % sop_q % po
    write(*,'(T3,A,T30,9(G0,X))') 'conv quadrature type:' , ins_op % sop_q % basis
    write(*,'(T3,A,T30,9(G0,X))') 'number of mesh points:', n_point
    write(*,'(T3,A,T30,9(G0,X))') 'time integrator:'      , trim(ins_ti % name)
    write(*,'(T3,A,T30,9(G0,X))') 'time step size:'       , dt
  end if

  !-----------------------------------------------------------------------------
  ! Time integration

  ! initial conditions
  call problem % GetExactSolution(ins_op % sem_v % metrics % x, t, u)

  do nt = 1, nt_max
    last = t + dt >= t_end .or. nt == nt_max
    call ins_ti % TimeStep(t, dt, u, standby = .not. last)
    if (last) exit
  end do

  !-----------------------------------------------------------------------------
  ! evaluation

  if (ins_op % mesh % part >= 0) then

    call problem % GetExactSolution(ins_op % sem_v % metrics % x, t, u_ex)

    ! calibrate pressure to zero mean value
    call CalibrateArray(p   , comm=ins_op % mesh % comm_parts)
    call CalibrateArray(p_ex, comm=ins_op % mesh % comm_parts)

    ! w = u - u_ex
    call SetArray(w, u, multi=.true.)
    call MergeArrays(ONE, w, -ONE, u_ex, multi=.true.)

    !$omp do
    do i = 1, ins_op % mesh % n_elem
      w(:,:,:,i,1) = w(:,:,:,i,1) ** 2 + w(:,:,:,i,2) ** 2 + w(:,:,:,i,3) ** 2
    end do
    call GetVolumeIntegral(ins_op%sem_v, w(:,:,:,:,1), e_v)
    e_v = sqrt(e_v)

    !$omp do
    do i = 1, ins_op % mesh % n_elem
      w(:,:,:,i,4) = w(:,:,:,i,4) ** 2
    end do
    call GetVolumeIntegral(ins_op%sem_v, w(:,:,:,:,4), e_p)
    e_p = sqrt(e_p)

    call EvalDivError(ins_op, v, e_div_vv, e_div_vp, e_div_vq)

    e_v      = e_v      / domain_volume
    e_p      = e_p      / domain_volume
    e_div_vv = e_div_vv / domain_volume
    e_div_vp = e_div_vp / domain_volume
    e_div_vq = e_div_vq / domain_volume

  end if

  if (ins_op % mesh % part == 0) then
    write(*,'(/,A)') 'time integration'
    write(*,'(T3,A,T29,ES12.5)') 'final time          t    =', t
    write(*,'(T3,A,T30,ES12.5)') 'velocity error      ε_v  =', e_v
    write(*,'(T3,A,T30,ES12.5)') 'pressure error      ε_p  =', e_p
    write(*,'(T3,A,T30,ES12.5)') 'div errors     ε_div_vv  =', e_div_vv
    write(*,'(T3,A,T30,ES12.5)') 'div errors     ε_div_vp  =', e_div_vp
    write(*,'(T3,A,T30,ES12.5)') 'div errors     ε_div_vq  =', e_div_vq
    write(*,*)
  end if

  !-----------------------------------------------------------------------------
  ! Write plot files

  if (export_vtk .and. ins_op % mesh%part >= 0) then
    call ExportVTK_VolumeData( x       = ins_op % sem_v % metrics % x &
                             , s       = var                          &
                             , sname   = var_name                     &
                             , file    = flow_case                    &
                             , part    = ins_op % mesh % part         &
                             , n_parts = ins_op % mesh % n_parts      )
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !> Evaluation of the divergence error

  subroutine EvalDivError(ins_op, v, e_div_vv, e_div_vp, e_div_vq)
    class(INS_Operator_3D), intent(in) :: ins_op
    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    real(RNP), intent(out) :: e_div_vv
    real(RNP), intent(out) :: e_div_vp
    real(RNP), intent(out) :: e_div_vq

    real(RNP), allocatable, save :: vp(:,:,:,:,:), div_v(:,:,:,:), q(:,:,:,:)

    integer :: ne, np, nq
    integer :: e

    ! L2 divergence in velocity space ..........................................

    ne = ins_op % mesh % n_elem
    np = ins_op % eop_v % po + 1

    !$omp master
    allocate(vp(np,np,6,ne,3), div_v(np,np,np,ne), q(np,np,np,ne))
    !$omp end master

    call GetOuterTraces_3D(ins_op%mesh, v, vp)
    call TPO_Div(ins_op % eop_v, ins_op % sem_v, v, vp, div_v)
    !$omp do
    do e = 1, ne
      q(:,:,:,e) = div_v(:,:,:,e) ** 2
    end do
    call GetVolumeIntegral(ins_op%sem_v, q, e_div_vv)
    e_div_vv = sqrt(e_div_vv)

    ! L2 divergence in pressure space ..........................................

    nq = ins_op % eop_p % po + 1

    !$omp master
    deallocate(q)
    allocate(q(nq,nq,nq,ne))
    !$omp end master

    call TPO_AAA(ins_op % iop_vp % A, div_v, q)
    !$omp do
    do e = 1, ne
      q(:,:,:,e) = q(:,:,:,e) ** 2
    end do
    call GetVolumeIntegral(ins_op%sem_p, q, e_div_vp)
    e_div_vp = sqrt(e_div_vp)

    ! L2 divergence in quadrature space ........................................

    nq = ins_op % sop_q % po + 1

    !$omp master
    deallocate(q)
    allocate(q(nq,nq,nq,ne))
    !$omp end master

    call TPO_AAA(ins_op % iop_vq % A, div_v, q)
    !$omp do
    do e = 1, ne
      q(:,:,:,e) = q(:,:,:,e) ** 2
    end do
    call GetVolumeIntegral(ins_op%sem_q, q, e_div_vq)
    e_div_vq = sqrt(e_div_vq)

    !$omp master
    deallocate(vp, div_v, q)
    !$omp end master

  end subroutine EvalDivError

  !=============================================================================

end program INS_TimeIntegrator_3D_Test

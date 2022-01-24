!> summary:  Validation of incompressible Navier-Stokes operators
!> author:   Joerg Stiller
!> date:     2022/01/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> If present, the first argument of the invoking command will be interpreted
!> as the base name of the control file. If omitted, the program looks for
!> `ins_operator_test.prm`.
!>
!> @todo
!>   - Extension to variable viscosity
!>   - Extension to free-slip boundary conditions
!===============================================================================
program INS_Operator_3D_Test
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Array_Assignments
  use Array_Reductions

  use Convert_Traces__3D
  use Spectral_Element_Variable__3D
  use Spectral_Element_Boundary_Variable__3D
  use Export_VTK_Volume_Data__3D
  use DG__Diffusion_Operator__3D

  use TPO__INS_Convection__3D_D__Gen

  use INS__Problem__3D
  use INS__Problem__3D__Vortex_TG
  use INS__Problem__3D__Variable_Viscosity
  use INS__Operator__3D

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

  character(len=*), parameter :: input_default = 'ins_operator_3d_test'
  character(len=80) :: input_file
  ! input file (*.prm)

  character(len=80) :: flow_problem = 'Vortex_TG'
  ! flow problem:
  !   'Vortex_TG'          2D Taylor-Green vortex
  !   'VariableViscosity'  3D Vortex array with variable viscosity

  character(len=80) :: problem_file = 'vortex_tg'
  ! file containing the problem parameters (*.prm)

  integer :: comp_domain = 1
  ! computational domain (u/s = un/structured, r = regular, d = deformed)
  !   1  cuboidal domain with Cartesian mesh                               (s+r)
  !   2  cuboidal domain with unstructured "diamond" mesh                  (u+d)
  !   3  cylindrical domain                                                (u+d)
  !   4  annular domain                                                    (u+d)
  ! configuration > 1 currently available only with one MPI process

  type(INS_Options_3D) :: ins_options
  ! options for the incompressible Navier-Stokes operator

  namelist/control_prm/ flow_problem, problem_file, comp_domain, ins_options

  integer :: n_test = 1            ! repetitions of performance test
  logical :: export_vtk = .false.  ! generate VTK files

  namelist/control_prm/ n_test, export_vtk

  ! operators and variables ....................................................

  class(INS_Problem_3D), allocatable :: problem
  ! flow problem

  type(INS_Operator_3D) :: ins_op
  ! incompressible Navier-Stokes operator

  type(DG_DiffusionOperator_3D) :: diffusion_op

  real(RNP) :: t = 0 ! problem time

  real(RNP), allocatable, target :: var(:,:,:,:,:)
  character(len=20), allocatable :: var_name(:)

  real(RNP), pointer, contiguous :: u(:,:,:,:,:)    ! u = [v,p]
  real(RNP), pointer, contiguous :: v(:,:,:,:,:)    ! velocity
  real(RNP), pointer, contiguous :: p(:,:,:,:)      ! pressure
  real(RNP), pointer, contiguous :: F_ce(:,:,:,:,:) ! exact  F_c = -∇⋅vv
  real(RNP), pointer, contiguous :: F_de(:,:,:,:,:) ! exact  F_d =  ∇⋅τ
  real(RNP), pointer, contiguous :: F_pe(:,:,:,:,:) ! exact  F_p = -∇p
  real(RNP), pointer, contiguous :: F_ch(:,:,:,:,:) ! approx F_c = -∇⋅vv
  real(RNP), pointer, contiguous :: F_dh(:,:,:,:,:) ! approx F_d =  ∇⋅τ
  real(RNP), pointer, contiguous :: F_ph(:,:,:,:,:) ! approx F_p = -∇p

  real(RNP), allocatable :: mm(:,:,:,:)    ! diagonal mass matrix
  real(RNP), allocatable :: um(:,:,:,:,:)  ! interior traces u⁻
  real(RNP), allocatable :: up(:,:,:,:,:)  ! exterior traces u⁺
  real(RNP), allocatable :: w(:,:,:,:,:)   ! workspace

  ! auxiliaries ................................................................

  type(SpectralElementVariable_3D) :: sev_u
  type(SpectralElementBoundaryVariable_3D), allocatable :: sebv_u(:)
  type(SpectralElementBoundaryVariable_3D), allocatable :: sebv_vi(:)

  character(len=80) :: domain_name = ''
! real(RDP) :: time, time0
  logical   :: exists
  integer   :: io, stat
  integer   :: n_bound, n_elem, n_elem_tot, n_ghost, n_point, n_var, po
  integer   :: i

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
    write(*,'(A)') 'Validation of incompressible Navier-Stokes operators'
    write(*,*)

    call get_command_argument(1, input_file, status=stat)
    if (stat /= 0 .or. len_trim(input_file) == 0) then
      input_file = input_default
    end if
    input_file = trim(input_file) // '.prm'

    inquire(file=trim(input_file), exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(input_file)
      open(newunit = io, file = input_file)
      read(io, nml = control_prm)
      close(io)
    else
       call Error( 'INS_Operator_3D_Test', &
                   'input file "' // trim(input_file) // '" not found' )
    end if

  end if

  ! globalize control parameters
  call XMPI_Bcast(flow_problem, 0, comm)
  call XMPI_Bcast(problem_file, 0, comm)
  call XMPI_Bcast(comp_domain , 0, comm)
  call ins_options % Bcast(0, comm)

  ! mesh .......................................................................

  select case(comp_domain)
  case(2)
    call CreateCuboidDiamonds(comm, input_file, ins_op % mesh)
    domain_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCylinder( comm, input_file, ins_op % mesh)
    domain_name = 'Cylindrical domain with unstructured mesh'
  case(4)
    call CreateAnnulus( comm, input_file, ins_op % mesh)
    domain_name = 'Annular domain with unstructured mesh'
  case default
    call CreateCuboidCartesian(comm, input_file, ins_op % mesh)
    domain_name = 'Cuboidal domain with Cartesian mesh'
  end select

  n_elem  = ins_op % mesh % n_elem
  n_ghost = ins_op % mesh % n_ghost
  n_bound = ins_op % mesh % n_bound

  ! problem ....................................................................

  select case(flow_problem)
  case('Vortex_TG')
    allocate(INS_Problem_3D_Vortex_TG         :: problem)
  case default
    allocate(INS_Problem_3D_VariableViscosity :: problem)
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

  call ins_op % Init(ins_options)

  ! variables ..................................................................

  po = ins_op % eop_v % po
  n_var = 22

  allocate(var(0:po,0:po,0:po,1:n_elem,1:n_var), source = ZERO)
  allocate(var_name(1:n_var))

  u(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,1:4)
  v(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,1:3)
  p(0:,0:,0:,1:)     =>  var(:,:,:,:,4)

  var_name(1:4) = [ 'v_x', 'v_y', 'v_z', 'p  ']

  F_ce(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,  5 :  7)
  F_de(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,  8 : 10)
  F_pe(0:,0:,0:,1:,1:)  =>  var(:,:,:,:, 11 : 13)

  var_name(5:13) = [ 'F_ce_x', 'F_ce_y', 'F_ce_z' &
                   , 'F_de_x', 'F_de_y', 'F_de_z' &
                   , 'F_pe_x', 'F_pe_y', 'F_pe_z' ]

  F_ch(0:,0:,0:,1:,1:)  =>  var(:,:,:,:, 14 : 16)
  F_dh(0:,0:,0:,1:,1:)  =>  var(:,:,:,:, 17 : 19)
  F_ph(0:,0:,0:,1:,1:)  =>  var(:,:,:,:, 20 : 22)

  var_name(14:22) = [ 'F_ch_x', 'F_ch_y', 'F_ch_z' &
                    , 'F_dh_x', 'F_dh_y', 'F_dh_z' &
                    , 'F_ph_x', 'F_ph_y', 'F_ph_z' ]

  allocate(mm (0:po,0:po,0:po,1:n_elem)     )
  allocate(w  (0:po,0:po,0:po,1:n_elem,1:4) )

  allocate(um (0:po,0:po,1:6,1:n_elem+n_ghost,1:4), source = ZERO )
  allocate(up (0:po,0:po,1:6,1:n_elem        ,1:4), source = ZERO )

  call ins_op % sem_v % Get_DG_DiagonalMassMatrix(mm)

  ! exact solution and terms ...................................................

  associate(x => ins_op % sem_v % metrics % x)

    call problem % GetExactSolution       (x, t, u)
    call problem % GetExactConvectiveTerm (x, t, F_ce)
    call problem % GetExactDiffusiveTerm  (x, t, F_de)
    call problem % GetExactPressureTerm   (x, t, F_pe)

  end associate

  ! traces and boundary values .................................................

  associate(sem => ins_op % sem_v, mesh => ins_op % mesh)

    ! wrap solution in spectral element variable
    sev_u  = SpectralElementVariable_3D(sem, u)

    ! traces
    call sev_u % GetTraces(um)
    call ConvertTraces(mesh, um, up)

    ! extract boundary values
    sebv_u = SpectralElementBoundaryVariable_3D(sem, mesh % boundary, nc = 4)
    call sebv_u % Extract(sev_u, mesh % boundary)

  end associate

  ! convection term: F_c = -∇·vv ...............................................

  associate(metrics => ins_op % sem_q % metrics)

    call TPO_INS_Convection_D_Gen( nv   = ins_op % eop_v % po + 1  &
                                 , nq   = ins_op % sop_q % po + 1  &
                                 , ne   = n_elem                   &
                                 , D_v  = ins_op % eop_v  % D      &
                                 , I_vq = ins_op % iop_vq % A      &
                                 , w_q  = ins_op % sop_q  % w      &
                                 , Jd_q = metrics % Jd             &
                                 , Ji_q = metrics % Ji             &
                                 , a_q  = metrics % a              &
                                 , n_q  = metrics % n              &
                                 , v    = v                        &
                                 , vp   = up(:,:,:,:,1:3)          &
                                 , F_c  = w (:,:,:,:,1:3)          )
  end associate


  do i = 1,3
    F_ch(:,:,:,:,i) = w(:,:,:,:,i) / mm
  end do

  ! pressure term: F_p = -∇p ...................................................

  ! viscous term, diffusion part: F_d1 = ∇·(ν ∇v) ..............................

  ! so far ν is constant and boundaries are periodic of have Dirichlet BC

  diffusion_op = DG_DiffusionOperator_3D( sem    = ins_op % sem_v      &
                                        , dg_opt = ins_options % opt_v &
                                        , lambda = ZERO                &
                                        , nu_p   = problem % nu_ref    &
                                        , bc     = problem % bc_v      )

  allocate(sebv_vi(n_bound))

  do i = 1, 3
    associate(r => w(:,:,:,:,i), f => w(:,:,:,:,4))

      ! store boundary contribution in f and compute residual r
      f = 0
      call sebv_u % GetSlice(sebv_vi, first=i, last=i)
      call diffusion_op % AddBC(sebv_vi, f)
      call diffusion_op % Apply(v(:,:,:,:,i), r, f)  ! r = -M ∇·(ν ∇vᵢ)

      ! compute nodal values
      F_dh(:,:,:,:,i) = -r / mm

    end associate
  end do

  ! viscous term, complete: F_d = ∇·τ ..........................................


  ! info .......................................................................

  call XMPI_Reduce(n_elem, n_elem_tot, MPI_SUM, 0, comm)
  n_point = n_elem_tot * (po+1)**3

  if (rank == 0) then
    write(*,*)
    write(*,'(T3,A,T30,9(G0,X))') 'problem:', trim(flow_problem)
    write(*,'(T3,A,T30,9(G0,X))') 'domain:',  trim(domain_name)
    write(*,'(T3,A,T30,9(G0,X))') 'boundary conditions:'  , problem % bc_v
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:'  , n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'    , n_thread
    write(*,'(T3,A,T30,9(G0,X))') 'number of elements:'   , n_elem_tot
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of v:', ins_op % eop_v % po
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of p:', ins_op % eop_p % po
    write(*,'(T3,A,T30,9(G0,X))') 'conv quadrature order:', ins_op % sop_q % po
    write(*,'(T3,A,T30,9(G0,X))') 'conv quadrature type:' , ins_op % sop_q % basis
    write(*,'(T3,A,T30,9(G0,X))') 'number of mesh points:', n_point
  end if

!?!  !-----------------------------------------------------------------------------
!?!  ! Consistency test
!?!
!?!  if (rank == 0) then
!?!    write(*,'(/,A,/)') 'Consistency'
!?!  end if
!?!
!?!  !$omp parallel
!?!
!?!  ! setup call
!?!  call diffusion_op % Apply(u, r)
!?!
!?!  !$omp master
!?!  if (rank == 0) time0 = MPI_Wtime()
!?!  !$omp end master
!?!
!?!  do i = 1, n_test
!?!    call diffusion_op % Apply(u, r, f)  ! r = A u - f
!?!!   call diffusion_op % Apply(u, r)     ! r = A u
!?!!   call MergeArrays(ONE, r, -ONE, f)   ! r = r - f
!?!  end do
!?!
!?!  !$omp master
!?!  if (rank == 0) time = MPI_Wtime()
!?!  !$omp end master
!?!
!?!  r_l2 = ScalarProduct(r, r, comm)
!?!  r_l2 = sqrt(r_l2)
!?!
!?!  !$omp end parallel
!?!
!?!  if (sem % mesh % part >= 0) then
!?!    r_max_loc = maxval(abs(r))
!?!  else
!?!    r_max_loc = 0
!?!  end if
!?!  call XMPI_Reduce(r_max_loc, r_max, MPI_MAX, 0, comm)
!?!
!?!  if (rank == 0) then
!?!    write(*,'(T3,A,T11,ES10.3)') 'r_L2  =', r_l2
!?!    write(*,'(T3,A,T11,ES10.3)') 'r_max =', r_max
!?!    if (n_test > 0) then
!?!      time = (time - time0) / n_test
!?!      write(*,'(T3,A,T11,ES10.3)') 't/DOF =', time / dof
!?!      write(*,'(T3,A,T11,ES10.3)') 'DOF/t =', dof / time
!?!    end if
!?!  end if

  !-----------------------------------------------------------------------------
  ! Write plot files

  if (export_vtk) then
    call ExportVTK_VolumeData( x      = ins_op % sem_v % metrics % x &
                             , s      = var                          &
                             , sname  = var_name                     &
                             , file   = 'ins_operator_3d_test'       &
                             , part   = ins_op % mesh % part         &
                             , n_part = ins_op % mesh % n_part       )
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !>
  !=============================================================================

end program INS_Operator_3D_Test

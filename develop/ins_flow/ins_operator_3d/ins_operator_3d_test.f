!> summary:  Validation of incompressible Navier-Stokes operators
!> author:   Joerg Stiller
!> date:     2022/01/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> If present, the first argument of the invoking command will be interpreted
!> as the base name of the control file. If omitted, the program looks for
!> `ins_operator_test.prm`.
!===============================================================================

program INS_Operator_3D_Test
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Array_Assignments
  use Array_Reductions

  use Mesh__3D
  use Trace_Operators__3D
  use Spectral_Element_Mesh__3D
  use Boundary_Variable__3D
  use Export_VTK_Volume_Data__3D
  use DG__Elliptic_Operator__3D
  use DG__Schwarz_Operator__3D

  use TPO__Div__3D
  use TPO__INS_Convection__3D_D__Gen

  use INS__Problem__3D
  use INS__Problem__Vortex_TG__3D
  use INS__Problem__Variable_Viscosity__3D
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

  integer :: flow_domain = 1
  ! computational domain (u/s = un/structured, r = regular, d = deformed)
  !   1  cuboidal domain with Cartesian mesh                               (s+r)
  !   2  cuboidal domain with unstructured "diamond" mesh                  (u+d)
  !   3  cylindrical domain                                                (u+d)
  !   4  annular domain                                                    (u+d)
  ! configuration > 1 currently available only with one MPI process

  namelist/control_prm/ flow_problem, problem_file, flow_domain

  integer :: po_u = 8 ! polynomial order of solution u \ p
  integer :: po_p = 7 ! polynomial order of pressure p
  type(INS_OperatorOptions_3D) :: ins_op_opt

  namelist/control_prm/ po_u, po_p, ins_op_opt

  integer   :: n_test = 1            ! repetitions of performance test
  integer   :: i_max  = 10
  real(RNP) :: r_red  = -1
  real(RNP) :: r_max  = -1
  logical   :: export_vtk = .false.  ! generate VTK files

  namelist/control_prm/ n_test, i_max, r_red, r_max, export_vtk

  ! operators and variables ....................................................

  class(INS_Problem_3D), allocatable :: problem

  type(Mesh_3D) :: mesh
  type(SpectralElementMesh_3D) :: sem_u
  type(SpectralElementMesh_3D) :: sem_p
  type(INS_Operator_3D) :: ins_op

  type(DG_EllipticOperator_3D) :: elliptic_op

  real(RNP) :: t = 0 ! problem time

  real(RNP), allocatable, target :: var(:,:,:,:,:)
  character(len=20), allocatable :: var_name(:)

  real(RNP), pointer, contiguous :: u(:,:,:,:,:)    ! u = [v,p]
  real(RNP), pointer, contiguous :: v(:,:,:,:,:)    ! velocity
  real(RNP), pointer, contiguous :: p(:,:,:,:)      ! pressure
  real(RNP), pointer, contiguous :: nu(:,:,:,:)     ! shear viscosity ν
  real(RNP), pointer, contiguous :: F_ce(:,:,:,:,:) ! exact  F_c = -∇⋅vv
  real(RNP), pointer, contiguous :: F_de(:,:,:,:,:) ! exact  F_d =  ∇⋅τ
  real(RNP), pointer, contiguous :: F_pe(:,:,:,:,:) ! exact  F_p = -∇p
  real(RNP), pointer, contiguous :: F_ch(:,:,:,:,:) ! approx F_c = -∇⋅vv
  real(RNP), pointer, contiguous :: F_dh(:,:,:,:,:) ! approx F_d =  ∇⋅τ
  real(RNP), pointer, contiguous :: F_ph(:,:,:,:,:) ! approx F_p = -∇p
  real(RNP), pointer, contiguous :: p_h(:,:,:,:)    ! approx p

  real(RNP), allocatable :: mm(:,:,:,:)    ! diagonal mass matrix
  real(RNP), allocatable :: mu(:,:,:,:)    ! bulk viscosity μ
  real(RNP), allocatable :: up(:,:,:,:,:)  ! exterior traces u⁺
  real(RNP), allocatable :: sp(:,:,:,:,:)  ! exterior traces s⁺ = n⋅τ⁺
  real(RNP), allocatable :: w(:,:,:,:,:)   ! workspace

  ! auxiliaries ................................................................

  type(BoundaryVariable_3D), allocatable :: bv_u(:)
  type(BoundaryVariable_3D), allocatable :: bv_vi(:)

  character(len=80) :: domain_name = ''
! real(RDP) :: time, time0
  real(RNP) :: e_c, e_d, e_d1, d_d1
  logical   :: exists
  integer   :: io, stat
  integer   :: n_bound, n_elem, n_elem_tot, n_ghost, n_point, n_var
  integer   :: b, i

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
  call XMPI_Bcast(flow_domain , 0, comm)

  call XMPI_Bcast(po_u, 0, comm)
  call XMPI_Bcast(po_p, 0, comm)

  call ins_op_opt % Bcast(0, comm)

  ! mesh .......................................................................

  select case(flow_domain)
  case(2)
    call CreateCuboidDiamonds(comm, input_file, mesh)
    domain_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCylinder(comm, input_file, mesh)
    domain_name = 'Cylindrical domain with unstructured mesh'
  case(4)
    call CreateAnnulus(comm, input_file, mesh)
    domain_name = 'Annular domain with unstructured mesh'
  case default
    call CreateCuboidCartesian(comm, input_file, mesh)
    domain_name = 'Cuboidal domain with Cartesian mesh'
  end select

  n_elem  = mesh % n_elem
  n_ghost = mesh % n_ghost
  n_bound = mesh % n_bound

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
        problem % bc_p(i) = 'P'
      end if
    end if
  end do

  ! operators ..................................................................

  sem_u = SpectralElementMesh_3D(mesh, po_u)
  sem_p = SpectralElementMesh_3D(mesh, po_p)
  call ins_op % Init(ins_op_opt, problem, sem_u, sem_p)

  ! variables ..................................................................

  n_var = 24

  allocate(var(0:po_u,0:po_u,0:po_u,1:n_elem,1:n_var), source = ZERO)
  allocate(var_name(1:n_var))

  u (0:,0:,0:,1:,1:)  =>  var(:,:,:,:,1:4)
  v (0:,0:,0:,1:,1:)  =>  var(:,:,:,:,1:3)
  p (0:,0:,0:,1:)     =>  var(:,:,:,:,4)
  nu(0:,0:,0:,1:)     =>  var(:,:,:,:,5)

  var_name(1:5) = [ 'v_x', 'v_y', 'v_z', 'p  ', 'nu ']

  F_ce(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,  6 :  8)
  F_de(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,  9 : 11)
  F_pe(0:,0:,0:,1:,1:)  =>  var(:,:,:,:, 12 : 14)

  var_name(6:14) = [ 'F_ce_x', 'F_ce_y', 'F_ce_z' &
                   , 'F_de_x', 'F_de_y', 'F_de_z' &
                   , 'F_pe_x', 'F_pe_y', 'F_pe_z' ]

  F_ch(0:,0:,0:,1:,1:)  =>  var(:,:,:,:, 15 : 17)
  F_dh(0:,0:,0:,1:,1:)  =>  var(:,:,:,:, 18 : 20)
  F_ph(0:,0:,0:,1:,1:)  =>  var(:,:,:,:, 21 : 23)

  var_name(15:23) = [ 'F_ch_x', 'F_ch_y', 'F_ch_z' &
                    , 'F_dh_x', 'F_dh_y', 'F_dh_z' &
                    , 'F_ph_x', 'F_ph_y', 'F_ph_z' ]

  p_h(0:,0:,0:,1:)  =>  var(:,:,:,:,24)

  var_name(24) = 'p_h'

  allocate(mm (0:po_u,0:po_u,0:po_u,1:n_elem)     )
  allocate(mu (0:po_u,0:po_u,0:po_u,1:n_elem)     )
  allocate(w  (0:po_u,0:po_u,0:po_u,1:n_elem,1:4) )

  allocate(up (0:po_u,0:po_u,1:6,1:n_elem,1:4), source = ZERO )
  allocate(sp (0:po_u,0:po_u,1:6,1:n_elem,1:3), source = ZERO )

  call ins_op % sem_u % Get_DG_DiagonalMassMatrix(mm)

  ! info .......................................................................

  call XMPI_Reduce(n_elem, n_elem_tot, MPI_SUM, 0, comm)
  n_point = n_elem_tot * (po_u+1)**3

  if (rank == 0) then
    write(*,'(/,A)') 'initialization of operators'
    write(*,'(T3,A,T30,9(G0,X))') 'problem:', trim(flow_problem)
    write(*,'(T3,A,T30,9(G0,X))') 'domain:',  trim(domain_name)
    write(*,'(T3,A,T30,9(G0,X))') 'boundary conditions for v:', problem % bc_v
    write(*,'(T3,A,T30,9(G0,X))') 'boundary conditions for p:', problem % bc_p
    write(*,'(T3,A,T30,9(G0,X))') 'variable properties:'  , &
                                               problem % HasVariableProperties()
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:'  , n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'    , n_thread
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of v:', po_u
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of p:', po_p
    write(*,'(T3,A,T30,9(G0,X))') 'conv quadrature order:', ins_op % sop_q % po
    write(*,'(T3,A,T30,9(G0,X))') 'number of mesh points:', n_point
  end if

  ! exact solution and terms ...................................................

  associate(x => ins_op % sem_u % metrics % x)

    call problem % GetExactSolution       (x, t, u)
    call problem % GetViscosity           (x, t, u, nu)
    call problem % GetExactConvectiveTerm (x, t, F_ce)
    call problem % GetExactDiffusiveTerm  (x, t, F_de)
    call problem % GetExactPressureTerm   (x, t, F_pe)

  end associate

  call SetArray(mu, ins_op % mu_0)

  ! traces and boundary values .................................................

  associate(sem => ins_op % sem_u)

    ! outer traces u⁺
    call GetOuterTraces_3D(mesh, u, up)

    ! boundary values
    allocate(bv_u(mesh % n_bound))
    do i = 1, mesh % n_bound
      call bv_u(i) % Init(mesh % boundary(i), po_u, nc = 4)
      call bv_u(i) % Extract(u)
    end do

  end associate

  ! convection term: F_c = -∇·vv ...............................................

  associate(metrics => ins_op % sem_q % metrics)

    call TPO_INS_Convection_D_Gen( nv   = ins_op % eop_u % po + 1  &
                                 , nq   = ins_op % sop_q % po + 1  &
                                 , ne   = n_elem                   &
                                 , D_v  = ins_op % eop_u  % D      &
                                 , I_vq = ins_op % iop_uq % A      &
                                 , w_q  = ins_op % sop_q  % w      &
                                 , Jd_q = metrics % Jd             &
                                 , Ji_q = metrics % Ji             &
                                 , a_q  = metrics % a              &
                                 , n_q  = metrics % n              &
                                 , v    = v                        &
                                 , vp   = up(:,:,:,:,1:3)          &
                                 , F_c  = w (:,:,:,:,1:3)          )
  end associate

  do i = 1, 3
    F_ch(:,:,:,:,i) = w(:,:,:,:,i) / mm
  end do

  ! error
  w(:,:,:,:,1:3) = F_ch - F_ce
  e_c = ScalarProduct(w(:,:,:,:,1:3), w(:,:,:,:,1:3), comm)
  e_c = sqrt(e_c / n_point)

  ! viscous term: F_d = ∇·τ ....................................................

  if (problem % HasVariableProperties()) then
    call ins_op % GetDiffusionTerm_V(mu, nu, v, up, sp, w, bv_u)
  else
    call ins_op % GetDiffusionTerm_C(v, up, sp, w, bv_u)
  end if

  do i = 1, 3
    F_dh(:,:,:,:,i) = w(:,:,:,:,i) / mm
  end do

  ! error
  w(:,:,:,:,1:3) = F_dh - F_de
  e_d = ScalarProduct(w(:,:,:,:,1:3), w(:,:,:,:,1:3), comm)
  e_d = sqrt(e_d / n_point)

  if (.not. problem % HasVariableProperties()) then

    ! diffusion part: w = ∇·(ν ∇v) .............................................
    ! using F_ph as workspace for F_d1h

    elliptic_op = DG_EllipticOperator_3D(ins_op%sem_u, DG_SchwarzOptions_3D())

    allocate(bv_vi(n_bound))

    do i = 1, 3
      associate(r => w(:,:,:,:,1), f => w(:,:,:,:,4), v_i => v(:,:,:,:,i))

        f = 0
        do b = 1, n_bound
          call bv_u(b) % GetSlice(bv_vi(b), first=i, last=i)
        end do
        ! r = M ∇·(ν ∇vᵢ)
        if (problem % HasVariableProperties()) then
          call elliptic_op % Residual( ins_op % problem % bc_v, ZERO &
                                     , nu, f, bv_vi, v_i, r          )
        else
          call elliptic_op % Residual( ins_op % problem % bc_v, ZERO      &
                                     , problem % nu_ref, f, bv_vi, v_i, r )
        end if
        ! compute nodal values
        F_ph(:,:,:,:,i) = r / mm

      end associate
    end do

    ! error (if ν is constant)
    w(:,:,:,:,1:3) = F_ph - F_de
    e_d1 = ScalarProduct(w(:,:,:,:,1:3), w(:,:,:,:,1:3), comm)
    e_d1 = sqrt(e_d1 / n_point)

    ! deviation between F_dh and F_d1h
    w(:,:,:,:,1:3) = F_ph - F_dh
    d_d1 = ScalarProduct(w(:,:,:,:,1:3), w(:,:,:,:,1:3), comm)
    d_d1 = sqrt(d_d1 / n_point)

  end if

  !-----------------------------------------------------------------------------
  ! Result info

  if (rank == 0) then
    write(*,'(/,A)') 'operator evaluation'
    write(*,'(T3,A,T30,ES12.5)') 'convection          ε_c  =', e_c
    write(*,'(T3,A,T30,ES12.5)') 'diffusion           ε_d  =', e_d
    if (.not. problem % HasVariableProperties()) then
      write(*,'(T3,A,T30,ES12.5)') '                    ε_d1 =', e_d1
      write(*,'(T3,A,T30,ES12.5)') '                    δ_d1 =', d_d1
    end if
    write(*,*)
  end if

  !-----------------------------------------------------------------------------
  ! Write plot files

  if (export_vtk) then
    call ExportVTK_VolumeData( x       = sem_u % metrics % x    &
                             , s       = var                    &
                             , sname   = var_name               &
                             , file    = 'ins_operator_3d_test' &
                             , part    = mesh % part            &
                             , n_parts = mesh % n_parts         )
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

end program INS_Operator_3D_Test

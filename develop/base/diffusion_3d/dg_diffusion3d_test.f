program DG_Diffusion3D_Test
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Array_Assignments
  use Array_Reductions

  use TPO__Diagonal__3D
  use Mesh__3D
  use Spectral_Element_Mesh__3D
  use Export_VTK_Volume_Data__3D
  use DG__Diffusion_Operator__3D

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  use Elliptic_Problem
  use Elliptic_Problem__Simple_1D
  use Elliptic_Problem__Simple_2D
  use Elliptic_Problem__Simple_3D
  use Elliptic_Problem__Knotty

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! control parameters .........................................................

  character(len=*), parameter :: input = 'dg_diffusion3d_test.prm'
  ! input file

  integer :: config = 1
  ! configuration (u/s = un/structured, r = regular, d = deformed)
  !   1  cuboidal domain with Cartesian mesh                               (s+r)
  !   2  cuboidal domain with unstructured "diamond" mesh                  (u+d)
  !   3  cylindrical domain                                                (u+d)
  !   4  annular domain                                                    (u+d)
  ! configuration > 1 currently available only with one MPI process

  integer :: n_test = 1            ! repetitions of consistency test
  logical :: export_vtk = .false.  ! generate VTK files

  namelist/control_prm/ config, n_test, export_vtk

  ! problem parameters .........................................................

  integer :: test_case  = 3            ! set 1/2/3/4 for simple_1/2/3d or knotty
  logical :: has_variable_nu = .false. ! set T/F for variable/constant ν
  logical :: has_spectral_nu = .false. ! set T/F to use/discard spectral ν

  namelist/problem_prm/ test_case, has_variable_nu, has_spectral_nu

  ! NOTE:
  !   -  spectral diffusivity  (nu_s) available in case of constant ν only
  !   -  fluctuation amplitude (nu_1) ignored   in case of constant ν

  real(RNP) :: lambda = 0     ! Helmholtz parameter
  real(RNP) :: nu_0   = 1     ! diffusivity mean value ν₀
  real(RNP) :: nu_1   = 0     ! diffusivity fluctuation amplitude ν₁
  real(RNP) :: nu_s   = 0.1   ! spectral diffusivity amplitude
  real(RNP) :: d_nu   = 0     ! diffusivity fluctuation phase shift
  integer   :: k_nu   = 1     ! diffusivity fluctuation wave number
  integer   :: k_u    = 1     ! solution wave number

  namelist/problem_prm/ lambda, nu_0, nu_1, nu_s, d_nu, k_nu, k_u

  ! boundary conditions {'D','N','P'} -- initialized with configuration
  character, allocatable :: bc(:)

  ! discretization parameters ..................................................

  ! NOTE: domain and mesh parameters are set via test case

  integer   :: pg      = 1    ! polynomial degree of geometry description
  integer   :: po      = 7    ! polynomial degree of spectral elements
  real(RNP) :: penalty = 2    ! penalty parameter > 1

  namelist/dicretization_prm/ pg, po, penalty

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm      ! MPI communicator
  integer        :: rank      ! local MPI rank
  integer        :: n_proc    ! number of MPI processes
  integer        :: n_thread  ! number of OpenMP threads

  ! problem variables ..........................................................

  type(Mesh_3D)                 :: mesh
  type(SpectralElementMesh_3D)  :: sem
  type(DG_ElementOptions_1D)    :: dg_opt
  type(DG_DiffusionOperator_3D) :: diffusion_op

  class(EllipticProblem), allocatable :: problem

  real(RNP), allocatable, target :: var(:,:,:,:,:)
  real(RNP), pointer, contiguous :: s (:,:,:,:)  ! exact solution
  real(RNP), pointer, contiguous :: u (:,:,:,:)  ! approximate solution
  real(RNP), pointer, contiguous :: nu(:,:,:,:)  ! diffusivity
  real(RNP), pointer, contiguous :: f (:,:,:,:)  ! RHS
  real(RNP), pointer, contiguous :: r (:,:,:,:)  ! residual
  real(RNP), pointer, contiguous :: e (:,:,:,:)  ! error

  real(RNP), allocatable :: mm(:,:,:,:) ! diagonal mass matrix

  ! auxiliary variables ........................................................

  character(len=80) :: config_name = '', test_case_name = ''
  real(RDP) :: time, time0
  real(RNP) :: r_l2, r_max, r_max_loc
  integer   :: io, n_bound, n_elem, n_var
  integer   :: n_elem_tot, dof
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
    open(newunit = io, file = input)
    read(io, nml = control_prm)
    read(io, nml = problem_prm)
    read(io, nml = dicretization_prm)
    close(io)

    if (has_variable_nu) then
      has_spectral_nu = .false.
      nu_s = 0
    else
      has_spectral_nu = has_spectral_nu .and. nu_s > 0
      nu_1 = 0
    end if

  end if

  ! globalize control parameters
  call XMPI_Bcast( config          , 0, comm )
  call XMPI_Bcast( n_test          , 0, comm )
  call XMPI_Bcast( export_vtk      , 0, comm )

  ! globalize problem parameters
  call XMPI_Bcast( test_case       , 0, comm )
  call XMPI_Bcast( has_variable_nu , 0, comm )
  call XMPI_Bcast( has_spectral_nu , 0, comm )
  call XMPI_Bcast( lambda          , 0, comm )
  call XMPI_Bcast( nu_0            , 0, comm )
  call XMPI_Bcast( nu_1            , 0, comm )
  call XMPI_Bcast( nu_s            , 0, comm )
  call XMPI_Bcast( d_nu            , 0, comm )
  call XMPI_Bcast( k_nu            , 0, comm )
  call XMPI_Bcast( k_u             , 0, comm )

  ! globalize discretization parameters
  call XMPI_Bcast( pg              , 0, comm )
  call XMPI_Bcast( po              , 0, comm )
  call XMPI_Bcast( penalty         , 0, comm )

  ! mesh .......................................................................

  select case(config)
  case(2)
    call CreateCuboidDiamonds( comm, input, pg, mesh, bc )
    config_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCylinder( comm, input, pg, mesh, bc )
    config_name = 'Cylindrical domain with unstructured mesh'
  case(4)
    call CreateAnnulus( comm, input, pg, mesh, bc )
    config_name = 'Annular domain with unstructured mesh'
  case default
    call CreateCuboidCartesian( comm, input, pg, mesh, bc )
    config_name = 'Cuboidal domain with Cartesian mesh'
  end select

  sem = SpectralElementMesh_3D(mesh, po)

  ! problem ....................................................................

  select case(test_case)
  case(1)
    allocate(EllipticProblem_Simple1D :: problem)
    test_case_name = 'Simple 1D'
  case(2)
    allocate(EllipticProblem_Simple2D :: problem)
    test_case_name = 'Simple 2D'
  case(3)
    allocate(EllipticProblem_Simple3D :: problem)
    test_case_name = 'Simple 3D'
  case default
    allocate(EllipticProblem_Knotty   :: problem)
    test_case_name = 'Knotty'
  end select

  call problem % SetProblem(lambda, nu_0, nu_1, d_nu, k_nu, k_u)

  ! variables ..................................................................

  n_elem  = sem % mesh % n_elem
  n_bound = sem % mesh % n_bound
  n_var   = 6

  allocate(var(0:po,0:po,0:po,1:n_elem,1:n_var), source = ZERO)

  s (0:,0:,0:,1:) => var(:,:,:,:,1)
  u (0:,0:,0:,1:) => var(:,:,:,:,2)
  nu(0:,0:,0:,1:) => var(:,:,:,:,3)
  f (0:,0:,0:,1:) => var(:,:,:,:,4)
  r (0:,0:,0:,1:) => var(:,:,:,:,5)
  e (0:,0:,0:,1:) => var(:,:,:,:,6)

  allocate(mm(0:po,0:po,0:po,1:n_elem))
  call sem % Get_DG_DiagonalMassMatrix(mm)

  ! solution and RHS ...........................................................

  associate(x => sem % metrics % x)

    ! exact solution and diffusivity
    call problem % GetExactSolution(x, s)
    call problem % GetDiffusivity(x, nu)
    u = s

    ! r = λ u - ∇·(ν ∇u)
    call problem % GetSource(x, r)

    ! project source:  f = M r
    f = mm * r

    ! apply boundary values:  f = f + f_bc
    !***TBD***!

  end associate

  ! operators ..................................................................

  dg_opt = DG_ElementOptions_1D(po, svv = has_spectral_nu, penalty = penalty)

  if (has_variable_nu) then
    diffusion_op = DG_DiffusionOperator_3D(sem, dg_opt, lambda, nu, bc)
  else if (has_spectral_nu) then
    diffusion_op = DG_DiffusionOperator_3D(sem, dg_opt, lambda, nu_0, nu_s, bc)
  else
    diffusion_op = DG_DiffusionOperator_3D(sem, dg_opt, lambda, nu_0, bc)
  end if

  ! info .......................................................................

  call XMPI_Reduce(n_elem, n_elem_tot, MPI_SUM, 0, comm)
  dof = n_elem_tot * (po+1)**3

  if (rank == 0) then
    write(*,'(/,A,/)') 'DG Diffusion Test'
    write(*,'(T3,A,T25,9(G0,X))') 'configuration:',config,' ',trim(config_name)
    write(*,'(T3,A,T25,9(G0,X))') 'test case:',test_case,' ',trim(test_case_name)
    write(*,'(T3,A,T25,9(G0,X))') 'spectral diffusivity:', has_spectral_nu
    write(*,'(T3,A,T25,9(G0,X))') 'variable diffusivity:', has_variable_nu
    write(*,'(T3,A,T25,9(G0,X))') 'boundary conditions:' , bc
    write(*,'(T3,A,T25,9(G0,X))') 'number of processes:' , n_proc
    write(*,'(T3,A,T25,9(G0,X))') 'number of threads:'   , n_thread
    write(*,'(T3,A,T25,9(G0,X))') 'number of elements:'  , n_elem_tot
    write(*,'(T3,A,T25,9(G0,X))') 'polynomial order:'    , po
    write(*,'(T3,A,T25,9(G0,X))') 'degrees of freedom:'  , dof
  end if

  !-----------------------------------------------------------------------------
  ! Consistency test

  if (rank == 0) then
    write(*,'(/,A,/)') 'Consistency'
  end if

  !$omp parallel

  ! setup call
  call diffusion_op % Apply(u, r)

  !$omp master
  if (rank == 0) time0 = MPI_Wtime()
  !$omp end master

  do i = 1, n_test
    call diffusion_op % Apply(u, r, f)  ! r = A u - f
!   call diffusion_op % Apply(u, r)     ! r = A u
!   call MergeArrays(ONE, r, -ONE, f)   ! r = r - f
  end do

  !$omp master
  if (rank == 0) time = MPI_Wtime()
  !$omp end master

  r_l2 = ScalarProduct(r, r, comm)
  r_l2 = sqrt(r_l2)

  !$omp end parallel

  if (sem % mesh % part >= 0) then
    r_max_loc = maxval(abs(r))
  else
    r_max_loc = 0
  end if
  call XMPI_Reduce(r_max_loc, r_max, MPI_MAX, 0, comm)

  if (rank == 0) then
    write(*,'(T3,A,T11,ES10.3)') 'r_L2  =', r_l2
    write(*,'(T3,A,T11,ES10.3)') 'r_max =', r_max
    if (n_test > 0) then
      time = (time - time0) / n_test
      write(*,'(T3,A,T11,ES10.3)') 't/DOF =', time / dof
      write(*,'(T3,A,T11,ES10.3)') 'DOF/t =', dof / time
    end if
  end if

  !-----------------------------------------------------------------------------
  ! ...

!?!    allocate(area(n_bound))
!?!    call se_mesh % GetVolume(vol)
!?!    call se_mesh % GetSurfaceAreas(area)
!?!    write(*,'(A)') 'Spectral element mesh'
!?!    write(*,'(2X,A,G0)') 'volume  = ', vol
!?!    do i = 1, mesh % n_bound
!?!      write(*,'(2X,A,I0,A,G0)') 'area(',i,') = ', area(i)
!?!    end do

!?!    if (export_vtk) then
!?!      call ExportVTK_VolumeData( x, var                  &
!?!                               , sname  = ['r','e']      &
!?!                               , file   = 'element_mesh' &
!?!                               , part   = mesh % part    &
!?!                               , n_part = mesh % n_part  )
!?!
!?!    end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !>
  !=============================================================================

end program DG_Diffusion3D_Test

program DG_Diffusion3D_Test
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Spectral_Element_Mesh__3D
  use Export_VTK_Volume_Data__3D

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

  logical :: export_vtk = .false. ! generate VTK files

  namelist/control_prm/ config, export_vtk

  ! problem parameters .........................................................

  integer :: test_case   = 3       ! set 1/2/3 for simple 1/2/3d or 4 for knotty
  logical :: variable_nu = .false. ! set T/F for variable/constant ν

  namelist/problem_prm/ test_case, variable_nu

  ! NOTE:
  !   -  spectral diffusivity  (nu_s) available in case of constant ν only
  !   -  fluctuation amplitude (nu_1) ignored   in case of constant ν

  real(RNP) :: lambda = 0     ! Helmholtz parameter
  real(RNP) :: nu_0   = 1     ! diffusivity mean value ν₀
  real(RNP) :: nu_1   = 0     ! diffusivity fluctuation amplitude ν₁
  real(RNP) :: nu_s   = 0     ! spectral  diffusivity amplitude
  real(RNP) :: d_nu   = 0     ! diffusivity fluctuation phase shift
  integer   :: k_nu   = 1     ! diffusivity fluctuation wave number
  integer   :: k_u    = 1     ! solution wave number

  namelist/problem_prm/ lambda, nu_0, nu_1, nu_s, d_nu, k_nu, k_u

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

  type(SpectralElementMesh_3D) :: sem

  class(EllipticProblem), allocatable :: problem

  real(RNP), allocatable, target :: var(:,:,:,:,:)
  real(RNP), pointer :: s (:,:,:,:)  ! exact solution
  real(RNP), pointer :: u (:,:,:,:)  ! approximate solution
  real(RNP), pointer :: nu(:,:,:,:)  ! diffusivity
  real(RNP), pointer :: f (:,:,:,:)  ! RHS
  real(RNP), pointer :: r (:,:,:,:)  ! residual
  real(RNP), pointer :: e (:,:,:,:)  ! error


  ! auxiliary variables ........................................................

  integer :: io
  integer :: n_bound, n_elem, n_var

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

    if (n_proc > 1) config = 1  ! so far
  end if

  ! globalize control parameters
  call XMPI_Bcast( config      , 0, comm )
  call XMPI_Bcast( export_vtk  , 0, comm )

  ! globalize problem parameters
  call XMPI_Bcast( test_case   , 0, comm )
  call XMPI_Bcast( variable_nu , 0, comm )
  call XMPI_Bcast( lambda      , 0, comm )
  call XMPI_Bcast( nu_0        , 0, comm )
  call XMPI_Bcast( nu_1        , 0, comm )
  call XMPI_Bcast( nu_s        , 0, comm )
  call XMPI_Bcast( d_nu        , 0, comm )
  call XMPI_Bcast( k_nu        , 0, comm )
  call XMPI_Bcast( k_u         , 0, comm )

  ! globalize discretization parameters
  call XMPI_Bcast( pg          , 0, comm )
  call XMPI_Bcast( po          , 0, comm )
  call XMPI_Bcast( penalty     , 0, comm )

  ! mesh .......................................................................

  select case(config)
  case(2)
    call CreateCuboidDiamonds( comm, input, pg, po, sem )
  case(3)
    call CreateCylinder( comm, input, pg, po, sem )
  case(4)
    call CreateAnnulus( comm, input, pg, po, sem )
  case default
    call CreateCuboidCartesian( comm, input, pg, po, sem )
  end select

  ! problem ....................................................................

  select case(test_case)
  case(1)
    allocate(EllipticProblem_Simple1D :: problem)
  case(2)
    allocate(EllipticProblem_Simple2D :: problem)
  case(3)
    allocate(EllipticProblem_Simple3D :: problem)
  case default
    allocate(EllipticProblem_Knotty   :: problem)
  end select

  call problem % SetProblem(lambda, nu_0, nu_1, d_nu, k_nu, k_u)

  ! variables ..................................................................

  n_elem  = sem % mesh % n_elem
  n_bound = sem % mesh % n_bound
  n_var   = 6

  allocate(var(0:po,0:po,0:po,1:n_elem,1:n_var))

  s (0:,0:,0:,1:) => var(:,:,:,:,1)
  u (0:,0:,0:,1:) => var(:,:,:,:,2)
  nu(0:,0:,0:,1:) => var(:,:,:,:,3)
  f (0:,0:,0:,1:) => var(:,:,:,:,4)
  r (0:,0:,0:,1:) => var(:,:,:,:,5)
  e (0:,0:,0:,1:) => var(:,:,:,:,6)

  ! solution and RHS ...........................................................

  associate(x => sem % metrics % x)

    ! exact solution
    call problem % GetExactSolution(x, s)

    ! r = λ u - ∇·(ν ∇u)
    call problem % GetSource(x, r)

    ! f = M r
    ! TBD: DG_Projection

  end associate

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

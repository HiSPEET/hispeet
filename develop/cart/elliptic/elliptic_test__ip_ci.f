!> summary:  Test elliptic IP/DG operator with constant isotropic diffusivity
!> author:   Joerg Stiller
!> date:     2018/12/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Test elliptic solvers
!===============================================================================

program Elliptic_Test__IP_CI

  use Kind_Parameters, only: RDP, RNP, IXL
  use Constants, only: PI, ZERO, ONE
  use Array_Assignments
  use TPO_sDDD
  use Export_Volume_Data_To_VTK

  use XMPI

  use CART__Mesh_Partition
  use CART__Generate_Structured_Mesh
  use CART__Elliptic_Operator_IP

  use Elliptic_Test_Case

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! problem parameters .........................................................

  real(RNP) :: lambda  = 0       ! Helmholtz parameter
  real(RNP) :: nu      = 1       ! diffusivity
  integer   :: kappa   = 1       ! wave number
  real(RNP) :: xo(3)   = 0       ! corner closest to -infinity
  real(RNP) :: lx(3)   = 2*PI    ! domain extensions
  character :: bc(6)   = 'P'     ! boundary conditions {'P'|'D'|'N'}

  namelist /problem/ lambda, nu, kappa, xo, lx, bc

  ! discretization parameters ..................................................

  integer   :: np(3)   = 1       ! number of partitions in directions 1:3
  integer   :: ep(3)   = 2       ! elements per partition and direction

  integer   :: po      = 2       ! polynomial order
  real(RNP) :: penalty = 2       ! penalty parameter (> 1)

  namelist /discretization/ np, ep, po, penalty

  ! solution parameters ........................................................

  integer   :: method  = 1       ! CG,, Schwarz, p-MG or p-MG/CG {1|2|3|4}

  integer   :: i_max   = huge(1) ! max number of iterations/cycles
  real(RNP) :: r_red   = 1E-6    ! min residual reduction

  type(SchwarzOptions3D) :: schwarz_opt

  namelist /solver/ method
  namelist /solver_cg/ i_max, r_red
  namelist /solver_schwarz/ i_max, r_red, schwarz_opt

  ! MPI ........................................................................

  type(MPI_Comm) :: comm         ! communicator
  integer        :: comm_size    ! number of processes
  integer        :: rank         ! local rank

  ! mesh .......................................................................

  real(RNP)              :: dx(3)          ! spacing in directions 1:3
  logical                :: periodic(3)    ! periodic directions set true
  type(MeshPartition)    :: mesh           ! mesh partition
  real(RNP), allocatable :: x(:,:,:,:,:)   ! mesh points

  ! variables ..................................................................

  real(RNP), allocatable, target :: scalars(:)       ! storage for scalar fields
  character(len=80), allocatable :: scalar_names(:)  ! names of scalars

  real(RNP), dimension(:,:,:,:), pointer, contiguous :: s  ! exact solution
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: u  ! numeric solution
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: f  ! source
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: r  ! residual
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: e  ! error

  ! operators ..................................................................

  type(EllipticOperator3D_IP) :: elliptic_op

  ! input / output .............................................................

  character(len=80) :: parameter_file = 'elliptic_test__ip_ci'
  character(len=80) :: plot_file      = ''
  integer :: prm
  logical :: exists

  namelist /control/ plot_file

  !-----------------------------------------------------------------------------
  ! Initialization

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_size(comm, comm_size)
  call MPI_Comm_rank(comm, rank)

  if (rank == 0) then

    parameter_file = trim(parameter_file) // '.prm'
    inquire(file=parameter_file, exist=exists)
    if (exists) then
      open(newunit=prm, file=parameter_file)
      read(prm, nml=problem)
      read(prm, nml=discretization)
      read(prm, nml=solver)
      select case(method)
      case(1)
        read(prm, nml=solver_cg)
      case(2)
        read(prm, nml=solver_schwarz)
      !case(3:)
      !  read(prm, nml=solver_pmg)
      end select
      read(prm, nml=control)
      close(prm)
    end if

  end if

  ! problem parameters
  call XMPI_Bcast(lambda   , 0, comm)
  call XMPI_Bcast(nu       , 0, comm)
  call XMPI_Bcast(kappa    , 0, comm)
  call XMPI_Bcast(xo       , 0, comm)
  call XMPI_Bcast(lx       , 0, comm)
  call XMPI_Bcast(bc       , 0, comm)

  ! discretization parameters
  call XMPI_Bcast(np       , 0, comm)
  call XMPI_Bcast(ep       , 0, comm)
  call XMPI_Bcast(po       , 0, comm)
  call XMPI_Bcast(penalty  , 0, comm)

  ! solver
  call XMPI_Bcast(method   , 0, comm)

  ! CG/Schwarz options
  call XMPI_Bcast(i_max    , 0, comm)
  call XMPI_Bcast(r_red    , 0, comm)

  ! Schwarz options
  call schwarz_opt % Bcast(0, comm)

!===============================================================================

end program Elliptic_Test__IP_CI

!> summary:  Creation of a cuboid domain with regular Cartesian mesh
!> author:   Joerg Stiller
!> date:     2021/09/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Create_Cuboid_Cartesian
  use Kind_Parameters
  use Constants
  use XMPI
  use Mesh__3D
  use Generate_Regular_Mesh__3D
  implicit none
  private

  public :: CreateCuboidCartesian

contains

  !-----------------------------------------------------------------------------
  !> Creation of a cuboid domain with regular Cartesian mesh

  subroutine CreateCuboidCartesian(comm, input, pg, mesh, bc)

    ! arguments ................................................................

    type(MPI_Comm),         intent(in)  :: comm  !< MPI communicator
    character(len=*),       intent(in)  :: input !< input file
    integer,                intent(in)  :: pg    !< polynomial order of geometry
    type(Mesh_3D),          intent(out) :: mesh  !< spectral-element mesh
    character, allocatable, intent(out) :: bc(:) !< boundary conditions

    ! input parameters .........................................................

    real(RNP) :: xo(3) = 0     ! corner closest to -infinity
    real(RNP) :: lx(3) = 2*PI  ! domain extensions
    integer   :: np(3) = 1     ! number of partitions per direction
    integer   :: ep(3) = 2     ! elements per partition and direction

    logical   :: regular    = .true. ! define mesh to be (not) regular
    logical   :: structured = .true. ! define mesh to be (un)structured

    namelist/cuboid_cartesian_prm/ xo, lx, np, ep, bc, regular, structured

    ! auxiliary variables ......................................................

    integer   :: rank
    integer   :: io
    integer   :: nb = 6        ! number of boundaries
    logical   :: periodic(3)   ! F/T for non/periodic directions
    real(RNP) :: dx(3)         ! element spacing in directions 1:3

    ! initialization ...........................................................

    call MPI_Comm_rank(comm, rank)

    allocate(bc(nb), source = 'P')
    if (rank == 0) then
      open(newunit = io, file = input)
      read(io, nml = cuboid_cartesian_prm)
      close(io)
    end if

    call XMPI_Bcast(xo, 0, comm)
    call XMPI_Bcast(lx, 0, comm)
    call XMPI_Bcast(np, 0, comm)
    call XMPI_Bcast(ep, 0, comm)
    call XMPI_Bcast(bc, 0, comm)

    call XMPI_Bcast(regular   , 0, comm)
    call XMPI_Bcast(structured, 0, comm)

    periodic(1) = all(bc(1:2) == 'P')
    periodic(2) = all(bc(3:4) == 'P')
    periodic(3) = all(bc(5:6) == 'P')

    dx = lx / (np * ep)

    ! create mesh ..............................................................

    call GenerateRegularMesh(mesh, np, ep, xo, dx, periodic, comm, pg)

    mesh % regular = regular
    mesh % structured = structured

  end subroutine CreateCuboidCartesian

  !=============================================================================

end module Create_Cuboid_Cartesian

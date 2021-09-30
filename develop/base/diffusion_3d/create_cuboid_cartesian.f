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
  use Spectral_Element_Mesh__3D
  implicit none
  private

  public :: CreateCuboidCartesian

contains

  !-----------------------------------------------------------------------------
  !> Creation of a cuboid domain with regular Cartesian mesh

  subroutine CreateCuboidCartesian(comm, input, pg, po, sem)
    type(MPI_Comm),   intent(in)  :: comm  !< MPI communicator
    character(len=*), intent(in)  :: input !< input file
    integer,          intent(in)  :: pg    !< polynomial order of geometry
    integer,          intent(in)  :: po    !< polynomial order of elements
    type(SpectralElementMesh_3D), intent(out) :: sem !< spectral-element mesh

    ! input parameters .........................................................

    real(RNP) :: xo(3) = 0             ! corner closest to -infinity
    real(RNP) :: lx(3) = 2*PI          ! domain extensions
    integer   :: np(3) = 1             ! number of partitions per direction
    integer   :: ep(3) = 2             ! elements per partition and direction
    logical   :: periodic(3) = .false. !  F/T for non/periodic directions

    namelist/cuboid_cartesian_prm/ xo, lx, np, ep, periodic

    ! auxiliary variables ......................................................

    integer   :: rank
    integer   :: io
    real(RNP) :: dx(3)                 ! element spacing in directions 1:3

    type(Mesh_3D) :: mesh              ! mesh partition

    ! initialization ...........................................................

    call MPI_Comm_rank(comm, rank)

    if (rank == 0) then
      open(newunit = io, file = input)
      read(io, nml = cuboid_cartesian_prm)
      close(io)
    end if

    call XMPI_Bcast(xo      , 0, comm)
    call XMPI_Bcast(lx      , 0, comm)
    call XMPI_Bcast(np      , 0, comm)
    call XMPI_Bcast(ep      , 0, comm)
    call XMPI_Bcast(periodic, 0, comm)

    dx = lx / (np * ep)

    ! create mesh ..............................................................

    call GenerateRegularMesh(mesh, np, ep, xo, dx, periodic, comm, pg)
    sem = SpectralElementMesh_3D(mesh, po)

  end subroutine CreateCuboidCartesian

  !=============================================================================

end module Create_Cuboid_Cartesian

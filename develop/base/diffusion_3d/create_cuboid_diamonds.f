!> summary:  Creation of a cuboid with unstructured "diamond" mesh
!> author:   Joerg Stiller
!> date:     2021/09/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Create_Cuboid_Diamonds
  use Kind_Parameters
  use Constants
  use XMPI
  use Mesh__3D
  use Generic_Mesh__3D
  use Spectral_Element_Mesh__3D
  implicit none
  private

  public :: CreateCuboidDiamonds

contains

  !-----------------------------------------------------------------------------
  !> Creation of a cuboid domain with unstructured "diamond" mesh

  subroutine CreateCuboidDiamonds(comm, input, pg, po, sem)
    type(MPI_Comm),   intent(in)  :: comm  !< MPI communicator
    character(len=*), intent(in)  :: input !< input file
    integer,          intent(in)  :: pg    !< polynomial order of geometry
    integer,          intent(in)  :: po    !< polynomial order of elements
    type(SpectralElementMesh_3D), intent(out) :: sem !< spectral-element mesh

    ! input parameters .........................................................

    real(RNP) :: lx = 2*PI             ! domain length in x-direction
    real(RNP) :: ly = 2*PI             ! domain length in y-direction
    real(RNP) :: lz = 2*PI             ! domain length in z-direction
    integer   :: nx = 2                ! number of cells in x-direction
    integer   :: ny = 2                ! number of cells in y-direction
    integer   :: nz = 2                ! number of elements in z-direction
    logical   :: periodic(3) = .false. ! F/T for non/periodic directions

    namelist/cuboid_diamonds_prm/ lx, ly, lz, nx, ny, nz, periodic

    ! auxiliary variables ......................................................

    integer   :: rank
    integer   :: io

    type(GenericMesh_3D) :: generic_mesh
    type(Mesh_3D)        :: mesh

    ! initialization ...........................................................

    call MPI_Comm_rank(comm, rank)

    if (rank == 0) then
      open(newunit = io, file = input)
      read(io, nml = cuboid_diamonds_prm)
      close(io)
    end if

    call XMPI_Bcast(lx      , 0, comm)
    call XMPI_Bcast(ly      , 0, comm)
    call XMPI_Bcast(lz      , 0, comm)
    call XMPI_Bcast(nx      , 0, comm)
    call XMPI_Bcast(ny      , 0, comm)
    call XMPI_Bcast(nz      , 0, comm)
    call XMPI_Bcast(periodic, 0, comm)

    ! create mesh ..............................................................

    call generic_mesh % CreateDiamonds(lx, ly, lz, nx, ny, nz, pg, periodic)
    call mesh % ImportGenericMesh(generic_mesh, comm = comm)

    sem = SpectralElementMesh_3D(mesh, po)

  end subroutine CreateCuboidDiamonds

  !=============================================================================

end module Create_Cuboid_Diamonds

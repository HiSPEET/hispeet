!> summary:  Creation of a 3x3x3 cuboid with a rotated center element
!> author:   Joerg Stiller
!> date:     2022/03/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Create_Cuboid_OneRotated
  use Kind_Parameters
  use Constants
  use XMPI
  use Mesh__3D
  use Generic_Mesh__3D
  implicit none
  private

  public :: CreateCuboidOneRotated

contains

  !-----------------------------------------------------------------------------
  !> Creation of a cylindrical domain with an unstructured mesh

  subroutine CreateCuboidOneRotated(comm, input, mesh)

    ! arguments ................................................................

    type(MPI_Comm),   intent(in)  :: comm  !< MPI communicator
    character(len=*), intent(in)  :: input !< input file
    type(Mesh_3D),    intent(out) :: mesh  !< spectral-element mesh

    ! input parameters .........................................................

    real(RNP) :: xo(3) = 0             ! corner closest to -infinity
    real(RNP) :: lx(3) = 2*PI          ! domain extensions
    integer   :: pg = 3                ! polynomial order of geometry
    integer   :: rotation(3) = 0       ! x/y/z-rotation of center element
    logical   :: periodic(3) = .false. ! F/T for non/periodic directions

    namelist/cuboid_onerotated_prm/ xo, lx, pg, rotation, periodic

    ! auxiliary variables ......................................................

    integer :: rank
    integer :: io

    type(GenericMesh_3D) :: generic_mesh

    ! initialization ...........................................................

    call MPI_Comm_rank(comm, rank)

    if (rank == 0) then
      open(newunit = io, file = input)
      read(io, nml = cuboid_onerotated_prm)
      close(io)
    end if

    call XMPI_Bcast(xo      , 0, comm)
    call XMPI_Bcast(lx      , 0, comm)
    call XMPI_Bcast(pg      , 0, comm)
    call XMPI_Bcast(rotation, 0, comm)
    call XMPI_Bcast(periodic, 0, comm)

    ! create mesh ..............................................................

    call generic_mesh % CreateOneRotated(xo, lx, pg, rotation, periodic)
    call mesh % ImportGenericMesh(generic_mesh, comm = comm)

  end subroutine CreateCuboidOneRotated

  !=============================================================================

end module Create_Cuboid_OneRotated

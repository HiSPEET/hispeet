!> summary:  Creation of a annular domain with an unstructured mesh
!> author:   Joerg Stiller
!> date:     2021/09/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Create_Annulus
  use Kind_Parameters
  use Constants
  use XMPI
  use Mesh__3D
  use Generic_Mesh__3D
  implicit none
  private

  public :: CreateAnnulus

contains

  !-----------------------------------------------------------------------------
  !> Creation of a cylindrical domain with an unstructured mesh

  subroutine CreateAnnulus(comm, input, mesh)

    ! arguments ................................................................

    type(MPI_Comm),   intent(in)  :: comm  !< MPI communicator
    character(len=*), intent(in)  :: input !< input file
    type(Mesh_3D),    intent(out) :: mesh  !< spectral-element mesh

    ! input parameters .........................................................

    real(RNP) :: r0 = 0.8*PI        ! inner radius
    real(RNP) :: r1 = 1.0*PI        ! outer radius
    real(RNP) :: h  = 2.0*PI        ! height = axial extension
    integer   :: nr = 2             ! num elements in radial direction
    integer   :: np = 4             ! num elements in azimuthal direction ≥ 3
    integer   :: nz = 4             ! num elements in axial  direction
    integer   :: pg = 3             ! polynomial order of geometry
    logical   :: periodic = .false. ! switch for periodicity in axial direction

    namelist/annulus_prm/ r0, r1, h, nr, np, nz, pg, periodic

    ! auxiliary variables ......................................................

    integer   :: rank
    integer   :: io

    type(GenericMesh_3D) :: generic_mesh

    ! initialization ...........................................................

    call MPI_Comm_rank(comm, rank)

    if (rank == 0) then
      open(newunit = io, file = input)
      read(io, nml = annulus_prm)
      close(io)
    end if

    call XMPI_Bcast(r0, 0, comm)
    call XMPI_Bcast(r1, 0, comm)
    call XMPI_Bcast(h , 0, comm)
    call XMPI_Bcast(nr, 0, comm)
    call XMPI_Bcast(np, 0, comm)
    call XMPI_Bcast(nz, 0, comm)
    call XMPI_Bcast(pg, 0, comm)

    call XMPI_Bcast(periodic, 0, comm)

    ! create mesh ..............................................................

    call generic_mesh % CreateAnnulus(r0, r1, h, nr, np, nz, pg, periodic)
    call mesh % ImportGenericMesh(generic_mesh, comm = comm)

  end subroutine CreateAnnulus

  !=============================================================================

end module Create_Annulus

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
  use Spectral_Element_Mesh__3D
  implicit none
  private

  public :: CreateAnnulus

contains

  !-----------------------------------------------------------------------------
  !> Creation of a cylindrical domain with an unstructured mesh

  subroutine CreateAnnulus(comm, input, pg, po, sem)
    type(MPI_Comm),   intent(in)  :: comm  !< MPI communicator
    character(len=*), intent(in)  :: input !< input file
    integer,          intent(in)  :: pg    !< polynomial order of geometry
    integer,          intent(in)  :: po    !< polynomial order of elements
    type(SpectralElementMesh_3D), intent(out) :: sem !< spectral-element mesh

    ! input parameters .........................................................

    real(RNP) :: r0   = 0.8*PI      ! inner radius
    real(RNP) :: r1   = 1.0*PI      ! outer radius
    real(RNP) :: h    = 2.0*PI      ! height = axial extension
    integer   :: nr   = 2           ! num elements in radial direction
    integer   :: np   = 4           ! num elements in azimuthal direction ≥ 3
    integer   :: nz   = 4           ! num elements in axial  direction
    logical   :: periodic = .false. ! F/T for non/periodicity in axial direction

    namelist/annulus/ r0, r1, h, nr, np, nz, periodic

    ! auxiliary variables ......................................................

    integer   :: mpi_rank
    integer   :: io

    type(GenericMesh_3D) :: generic_mesh
    type(Mesh_3D)        :: mesh

    ! initialization ...........................................................

    call MPI_Comm_rank(comm, mpi_rank)

    if (mpi_rank == 0) then
      open(newunit = io, file = input)
      read(io, nml = annulus)
      close(io)
    end if

    call XMPI_Bcast(r0      , 0, comm)
    call XMPI_Bcast(r1      , 0, comm)
    call XMPI_Bcast(h       , 0, comm)
    call XMPI_Bcast(nr      , 0, comm)
    call XMPI_Bcast(np      , 0, comm)
    call XMPI_Bcast(nz      , 0, comm)
    call XMPI_Bcast(periodic, 0, comm)

    ! create mesh ..............................................................

    call generic_mesh % CreateAnnulus(r0, r1, h, nr, np, nz, pg, periodic)
    call mesh % ImportGenericMesh(generic_mesh, comm = comm)

    sem = SpectralElementMesh_3D(mesh, po)

  end subroutine CreateAnnulus

  !=============================================================================

end module Create_Annulus

!> summary:  Creation of a cylindrical domain with an unstructured mesh
!> author:   Joerg Stiller
!> date:     2021/09/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Create_Cylinder
  use Kind_Parameters
  use Constants
  use XMPI
  use Mesh__3D
  use Generic_Mesh__3D
  implicit none
  private

  public :: CreateCylinder

contains

  !-----------------------------------------------------------------------------
  !> Creation of a cylindrical domain with an unstructured mesh

  subroutine CreateCylinder(comm, input, pg, mesh, bc)

    ! arguments ................................................................

    type(MPI_Comm),         intent(in)  :: comm  !< MPI communicator
    character(len=*),       intent(in)  :: input !< input file
    integer,                intent(in)  :: pg    !< polynomial order of geometry
    type(Mesh_3D),          intent(out) :: mesh  !< spectral-element mesh
    character, allocatable, intent(out) :: bc(:) !< boundary conditions

    ! input parameters .........................................................

    real(RNP) :: r  = 1*PI        ! radius
    real(RNP) :: h  = 2*PI        ! height = axial extension
    integer   :: nr = 2           ! num elements in radial direction
    integer   :: nz = 2           ! num elements in axial  direction

    namelist/cylinder_prm/ r, h, nr, nz, bc

    ! auxiliary variables ......................................................

    integer   :: rank
    integer   :: io
    integer   :: nb = 3    ! num boundaries {shell, bottom, top}
    logical   :: periodic  ! F/T for non/periodicity in axial direction

    type(GenericMesh_3D) :: generic_mesh

    ! initialization ...........................................................

    call MPI_Comm_rank(comm, rank)

    allocate(bc(nb), source = 'P')
    if (rank == 0) then
      open(newunit = io, file = input)
      read(io, nml = cylinder_prm)
      close(io)
    end if

    call XMPI_Bcast(r , 0, comm)
    call XMPI_Bcast(h , 0, comm)
    call XMPI_Bcast(nr, 0, comm)
    call XMPI_Bcast(nz, 0, comm)
    call XMPI_Bcast(bc, 0, comm)

    periodic = all(bc(2:3) == 'P')

    ! create mesh ..............................................................

    call generic_mesh % CreateCylinder(r, h, nr, nz, pg, periodic)
    call mesh % ImportGenericMesh(generic_mesh, comm = comm)

  end subroutine CreateCylinder

  !=============================================================================

end module Create_Cylinder

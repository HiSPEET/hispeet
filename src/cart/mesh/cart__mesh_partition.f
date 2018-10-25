!> summary:  Cartesian mesh partition
!> author:   Joerg Stiller
!> date:     2017/04/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Cartesian mesh partition
!===============================================================================

module CART__Mesh_Partition

  use Kind_Parameters, only: RNP
  use XMPI
  use CART__Mesh_Element
  use CART__Mesh_Face
  use CART__Mesh_Boundary
  use CART__Mesh_Link

  implicit none
  private

  public :: MeshPartition
  public :: EmptyMeshPartition

  !-----------------------------------------------------------------------------
  !> Definition a Cartesian mesh partition
  !>
  !> ...
  !> In the structured case, all entities are numbered as defined in the
  !> Structured_Mesh module. To foster efficient flux computation, the
  !> mesh faces are always partitioned into in three groups according to
  !> their normal direction and ordered accordingly.

  type MeshPartition

    ! global attributes
    integer :: level      = 0        !< mesh level
    integer :: n_boundary = 0        !< number of domain boundaries
    integer :: n_part     = 0        !< number of non-empty partitions

    ! local attributes
    integer :: part       = -1       !< partition ID
    logical :: structured = .false.  !< mapping to structured mesh exists

    ! spacing
    real(RNP) :: dx(3)               !< mesh spacing in directions 1:3

    ! dimensions
    integer :: nf = 0  !< number of mesh faces
    integer :: ne = 0  !< number of mesh elements
    integer :: ng = 0  !< number of ghost elements

    ! number of mesh faces per orientation
    integer :: nf1 = 0  !< number of faces normal to direction 1
    integer :: nf2 = 0  !< number of faces normal to direction 2
    integer :: nf3 = 0  !< number of faces normal to direction 3

    ! structured mesh dimensions
    integer :: ne1 = 0  !< number of elements in direction 1
    integer :: ne2 = 0  !< number of elements in direction 2
    integer :: ne3 = 0  !< number of elements in direction 3

    ! mesh components
    type(MeshElement),  allocatable :: element(:)  !< mesh elements
    type(MeshFace),     allocatable :: face(:)     !< mesh faces
    type(MeshBoundary), allocatable :: boundary(:) !< boundaries

    ! MPI
    type(MPI_Comm) :: comm  !< communicator

    ! mesh links
    type(MeshLink), allocatable :: link(:)

  contains

    procedure :: GetPoints
    procedure :: GetPointValency

  end type MeshPartition

  !-----------------------------------------------------------------------------
  !> Generates Gauss-Lobatto or Gauss points to all elements of a mesh partition

  interface
    module subroutine GetPoints(mesh, po, basis, x)
      use Kind_Parameters, only: RNP

      class(MeshPartition),   intent(in)  :: mesh         !< mesh parition
      integer,                intent(in)  :: po           !< polynomial order
      character(len=*),       intent(in)  :: basis        !< 'GL' or 'GLL'
      real(RNP), allocatable, intent(out) :: x(:,:,:,:,:) !< mesh points

    end subroutine GetPoints
  end interface

  !-----------------------------------------------------------------------------
  !> Computes the valency of mesh points

  interface
    module subroutine GetPointValency(mesh, basis, v)
      class(MeshPartition), intent(in)  :: mesh          !< mesh parition
      character(len=*),     intent(in)  :: basis         !< 'GL' or 'GLL'
      integer,              intent(out) :: v(0:,0:,0:,:) !< point valency
    end subroutine GetPointValency
  end interface

contains

!-------------------------------------------------------------------------------
!> Generate an empty mesh partition

function EmptyMeshPartition(level, n_boundary, n_part, dx, comm) result(mesh)
  integer,        intent(in) :: level       !< mesh level
  integer,        intent(in) :: n_boundary  !< number of external boundaries
  integer,        intent(in) :: n_part      !< number of non-empty partitions
  real(RNP),      intent(in) :: dx(3)       !< mesh spacing in directions 1:3
  type(MPI_Comm), intent(in) :: comm        !< MPI communicator
  optional :: comm

  type(MeshPartition) :: mesh

  mesh%level      = level
  mesh%n_boundary = n_boundary
  mesh%n_part     = n_part
  mesh%dx         = dx

  allocate(mesh%element(1:0))
  allocate(mesh%face(1:0))
  allocate(mesh%boundary(1:n_boundary))

  if (present(comm)) then
    mesh%comm = comm
  else
    mesh%comm = MPI_COMM_WORLD
  end if

end function EmptyMeshPartition

!===============================================================================

end module CART__Mesh_Partition

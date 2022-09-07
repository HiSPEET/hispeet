!> summary:  3D mesh boundary
!> author:   Joerg Stiller
!> date:     2020/11/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_Boundary__3D
  use XMPI
  implicit none
  private

  public :: MeshBoundary_3D
  public :: MeshBoundaryAttributes_3D

  !-----------------------------------------------------------------------------
  !> 3d mesh boundary face

  type MeshBoundaryFace_3D
    integer :: element_id   = 0 !< local element ID
    integer :: element_face = 0 !< adjacent face of the element {1:6}
  end type MeshBoundaryFace_3D

  !-----------------------------------------------------------------------------
  !> Mesh boundary
  !>
  !> This type provides the properties of a boundary and serves for accessing
  !> adjacent mesh faces and elements of the local partition.
  !> The components `id`, `name`, `coupled` and `polarity` represent global
  !> properties and, therefore, are same in all mesh partitions even if they
  !> do not share any part of the boundary.
  !> The numeric identifier `id` equals the index in the `boundary` component
  !> of the related `Mesh_3D` object, whereas `name` represents a textual label.
  !>
  !> Periodicity is supported via the `coupled` and `polarity` components.
  !> While the former specifies the ID of the coupled boundary, the latter
  !> allows to identify their relative location respective to the pertinent
  !> periodic direction. More precisely, `polarity` can assume the values
  !> 0, ±1, ±2 or ±3. For example, a polarity of -1 indicates a boundary
  !> limiting the domain toward -∞ in periodic direction 1. Similarly,
  !> +1 refers to a position on the opposite side, i.e. toward +∞. Thus, then
  !> `polarity` component enables the identification of all boundaries forming
  !> an extremity or "pole" of one periodic direction.
  !>
  !> The `face` component provides the required information for accessing the
  !> local mesh face and the element face adjoining a given boundary face.
  !> Once initialized, the number of local boundary faces is stored in `n_face`.

  type MeshBoundary_3D

    character(len=80) :: name     = ''  !< name
    integer           :: id       =  0  !< boundary identifier
    integer           :: coupled  =  0  !< ID of coupled boundary, 0 if none
    integer           :: polarity =  0  !< position WRT to periodic direction
    integer           :: n_face   =  0  !< number of faces

    type(MeshBoundaryFace_3D), allocatable :: face(:) !< boundary faces

  end type MeshBoundary_3D

  ! constructor
  interface MeshBoundary_3D
    module procedure New_MeshBoundary_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for collecting and transmitting the mesh boundary attributes

  type MeshBoundaryAttributes_3D
    character(len=80) :: name     = ''  !< name
    integer           :: id       =  0  !< boundary identifier
    integer           :: coupled  =  0  !< ID of coupled boundary, 0 if none
    integer           :: polarity =  0  !< position WRT to periodic direction
  contains
    procedure :: Bcast => Bcast_MeshBoundaryAttributes_3D
  end type MeshBoundaryAttributes_3D

  ! constructor
  interface MeshBoundaryAttributes_3D
    module procedure ExtractBoundaryAttributes
  end interface

contains

  !=============================================================================
  ! MeshBoundary_3D constructor and type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for MeshBoundary_3D

  function New_MeshBoundary_3D(attrib, n_face) result(this)
    class(MeshBoundaryAttributes_3D), intent(in) :: attrib !< attributes
    type(MeshBoundary_3D) :: this !< 3d mesh boundary object
    integer, optional, intent(in) :: n_face !< number of mesh faces [-1]

    call Init_MeshBoundary_3D(this, attrib, n_face)

  end function New_MeshBoundary_3D

  !-----------------------------------------------------------------------------
  !> Initialize a new 3d mesh boundary

  subroutine Init_MeshBoundary_3D(this, attrib, n_face)
    class(MeshBoundary_3D), intent(inout) :: this !< 3d mesh boundary object
    class(MeshBoundaryAttributes_3D), intent(in) :: attrib !< attributes
    integer, optional, intent(in) :: n_face !< number of mesh faces [-1]

    if (allocated(this%face)) deallocate(this%face)

    this % name      =  attrib % name
    this % id        =  attrib % id
    this % coupled   =  attrib % coupled
    this % polarity  =  attrib % polarity

    if (present(n_face)) then
      allocate(this % face( n_face ))
    end if

  end subroutine Init_MeshBoundary_3D

  !=============================================================================
  ! MeshBoundaryAttributes_3D constructor and type-bound procedures

  elemental function ExtractBoundaryAttributes(boundary) result(this)
    class(MeshBoundary_3D), intent(in) :: boundary
    type(MeshBoundaryAttributes_3D) :: this

    this % name      =  boundary % name
    this % id        =  boundary % id
    this % coupled   =  boundary % coupled
    this % polarity  =  boundary % polarity

  end function ExtractBoundaryAttributes

  !-----------------------------------------------------------------------------
  !> Broadcast mesh boundary attributes

  subroutine Bcast_MeshBoundaryAttributes_3D(this, root, comm)
    class(MeshBoundaryAttributes_3D), intent(inout) :: this
    integer       , intent(in) :: root !< MPI root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    integer :: attrib_int(3)

    call XMPI_Bcast(this % name, root, comm)

    attrib_int(1) = this % id
    attrib_int(2) = this % coupled
    attrib_int(3) = this % polarity

    call XMPI_Bcast(attrib_int, root, comm)

    this % id       = attrib_int(1)
    this % coupled  = attrib_int(2)
    this % polarity = attrib_int(3)

  end subroutine Bcast_MeshBoundaryAttributes_3D

  !=============================================================================

end module Mesh_Boundary__3D

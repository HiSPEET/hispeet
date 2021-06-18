!> summary:  3D mesh boundary
!> author:   Joerg Stiller
!> date:     2020/11/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_Boundary__3D
  implicit none
  private

  public :: MeshBoundary_3D

  !-----------------------------------------------------------------------------
  !> Adjacent mesh face data

  type AdjacentFace
    integer :: id   = 0 !< mesh face ID
    integer :: side = 0 !< corresponding side of the face {1,2}
  end type AdjacentFace

  !-----------------------------------------------------------------------------
  !> Adjacent mesh element data

  type AdjacentElement
    integer :: id   = 0 !< element ID
    integer :: face = 0 !< corresponding face of the element {1:6}
  end type AdjacentElement

  !-----------------------------------------------------------------------------
  !> 3d mesh boundary face

  type MeshBoundaryFace_3D
    type(AdjacentFace)    :: mesh_face    !< adjacent mesh face
    type(AdjacentElement) :: mesh_element !< adjacent mesh element
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
  !> of the related `Mesh3D_Partition` object, whereas `name` represents a
  !> textual label.
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

    integer           :: id       =  0  !< boundary identifier
    character(len=80) :: name     = ''  !< name
    integer           :: coupled  =  0  !< ID of coupled boundary, 0 if none
    integer           :: polarity =  0  !< position WRT to periodic direction
    integer           :: n_face   = -1  !< number of faces

    type(MeshBoundaryFace_3D), allocatable :: face(:) !< boundary faces

  end type MeshBoundary_3D

  ! constructor
  interface MeshBoundary_3D
    module procedure New_MeshBoundary_3D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor for MeshBoundary_3D

  function New_MeshBoundary_3D(id, name, coupled, polarity, n_face) result(this)
    integer,           intent(in) :: id       !< identifier
    character(len=*),  intent(in) :: name     !< name
    integer, optional, intent(in) :: coupled  !< ID coupled boundary, 0 if none
    integer, optional, intent(in) :: polarity !< pos WRT to periodic direction
    integer, optional, intent(in) :: n_face   !< number of faces
    type(MeshBoundary_3D)         :: this     !< 3d mesh boundary object

    call Init_MeshBoundary_3D(this, id, name, coupled, polarity, n_face)

  end function New_MeshBoundary_3D

  !-----------------------------------------------------------------------------
  !> Initialize a new 3d mesh boundary

  subroutine Init_MeshBoundary_3D(this, id, name, coupled, polarity, n_face)
    class(MeshBoundary_3D), intent(inout) :: this !< 3d mesh boundary object
    integer,           intent(in) :: id       !< identifier
    character(len=*),  intent(in) :: name     !< name
    integer, optional, intent(in) :: coupled  !< ID coupled boundary, 0 if none
    integer, optional, intent(in) :: polarity !< pos WRT to periodic direction
    integer, optional, intent(in) :: n_face   !< number of faces

    if (allocated(this%face)) deallocate(this%face)

    this % id      = id
    this % name    = name

    if (present(coupled )) this % coupled  = coupled
    if (present(polarity)) this % polarity = polarity

    if (present(n_face)) then
      this % n_face  = n_face
      allocate(this%face(n_face))
    end if

  end subroutine Init_MeshBoundary_3D

  !=============================================================================

end module Mesh_Boundary__3D

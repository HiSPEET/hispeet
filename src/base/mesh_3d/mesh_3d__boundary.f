!> summary:  3D mesh boundary
!> author:   Joerg Stiller
!> date:     2020/11/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_3d__Boundary
  implicit none
  private

  public :: Mesh3d_Boundary

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

  type Mesh_3d_BoundaryFace
    type(AdjacentFace)    :: mesh_face    !< adjacent mesh face
    type(AdjacentElement) :: mesh_element !< adjacent mesh element
  end type Mesh_3d_BoundaryFace

  !-----------------------------------------------------------------------------
  !> Cartesian mesh boundary

  type Mesh3d_Boundary

    integer           :: id      =  0  !< boundary identifier
    character(len=80) :: name    = ''  !< name
    integer           :: n_face  = -1  !< number of faces
    integer           :: coupled =  0  !< ID of coupled boundary, 0 if none

    type(Mesh_3d_BoundaryFace), allocatable :: face(:) !< boundary faces

  end type Mesh3d_Boundary

  ! constructor
  interface Mesh3d_Boundary
    module procedure New_Mesh3d_Boundary
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor for Mesh3d_Boundary

  function New_Mesh3d_Boundary(id, name, n_face, coupled) result(this)
    integer,          intent(in) :: id       !< identifier
    character(len=*), intent(in) :: name     !< name
    integer,          intent(in) :: n_face   !< number of faces
    integer,          intent(in) :: coupled  !< ID coupled boundary, 0 if none
    type(Mesh3d_Boundary)        :: this     !< 3d mesh boundary object

    call Init_Mesh3d_Boundary(this, id, name, n_face, coupled)

  end function New_Mesh3d_Boundary

  !-----------------------------------------------------------------------------
  !> Initialize a new 3d mesh boundary

  subroutine Init_Mesh3d_Boundary(this, id, name, n_face, coupled)
    class(Mesh3d_Boundary), intent(inout) :: this !< 3d mesh boundary object
    integer,          intent(in) :: id       !< identifier
    character(len=*), intent(in) :: name     !< name
    integer,          intent(in) :: n_face   !< number of faces
    integer,          intent(in) :: coupled  !< ID coupled boundary, 0 if none

    if (allocated(this%face)) deallocate(this%face)

    this % id      = id
    this % name    = name
    this % n_face  = n_face
    this % coupled = coupled

    allocate(this%face(n_face))

  end subroutine Init_Mesh3d_Boundary

  !=============================================================================

end module Mesh_3d__Boundary

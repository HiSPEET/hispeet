!> summary:  Cartesian mesh boundary
!> author:   Joerg Stiller
!> date:     2017/04/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Cartesian mesh boundary
!>
!> @todo
!> Accelerate
!> @endtodo
!===============================================================================

module CART__Mesh_Boundary

  implicit none
  private

  public :: MeshBoundary

  !-----------------------------------------------------------------------------
  !> Mesh face data

  type MeshFaceData
    integer :: id = -1 !< face id
  end type MeshFaceData

  !-----------------------------------------------------------------------------
  !> Mesh element face data

  type MeshElementData
    integer :: id   = -1 !< element id
    integer :: face = -1 !< corresponding element face
  end type MeshElementData

  !-----------------------------------------------------------------------------
  !> Boundary face
  !>
  !> The orientation specifies the coordinate direction {1|2|3} and the sense
  !> {+|-} of the exterior normal vector. E.g., orientation = 1 indicates that
  !> the normal coincides with the direction vector e1.

  type BoundaryFace
    integer               :: orientation =  0  !< normal orientation
    type(MeshFaceData)    :: mesh_face         !< corresponding mesh face
    type(MeshElementData) :: mesh_element      !< adjacent mesh element
  end type BoundaryFace

  !-----------------------------------------------------------------------------
  !> Cartesian mesh boundary

  type MeshBoundary

    integer           :: id   =  0      !< boundary identifier
    character(len=80) :: name =  ''     !< name
    integer           :: nf   = -1      !< number of faces

    logical :: is_periodic = .false.    !< indicates periodic boundary
    logical :: is_interior = .false.    !< indicates interior boundary

    type(BoundaryFace), allocatable :: face(:) !< boundary faces

  contains

    generic :: New => New_MeshBoundary
    procedure, private :: New_MeshBoundary

  end type MeshBoundary

contains

!-------------------------------------------------------------------------------
!> Create a new mesh boundary structure

subroutine New_MeshBoundary(this, id, name, nf, periodic)
  class(MeshBoundary), intent(inout) :: this     !< mesh boundary object
  integer,             intent(in)    :: id       !< identifier
  character(len=*),    intent(in)    :: name     !< name
  integer,             intent(in)    :: nf       !< number of faces
  logical,             intent(in)    :: periodic !< switch for periodicity

  if (allocated(this%face)) deallocate(this%face)

  this % id          = id
  this % name        = name
  this % nf          = nf
  this % is_periodic = periodic

  allocate(this%face(nf))

end subroutine New_MeshBoundary

!===============================================================================

end module CART__Mesh_Boundary

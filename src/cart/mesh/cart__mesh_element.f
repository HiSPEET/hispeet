!> summary:  Cartesian mesh element
!> author:   Joerg Stiller
!> date:     2017/04/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Cartesian mesh boundary
!===============================================================================

module CART__Mesh_Element
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: MeshElement

  !-----------------------------------------------------------------------------
  !> Structure for keeping element face data

  type ElementFace
    integer :: id       = 0  !< mesh face ID
    integer :: neighbor = 0  !< neighbor element, including virtual one
    integer :: boundary = 0  !< adjacent boundary, if any
  end type ElementFace

  !-----------------------------------------------------------------------------
  !> Structure for keeping element edge data

  type ElementEdge
    integer :: neighbor = 0  !< neighbor element, including virtual one
  end type ElementEdge

  !-----------------------------------------------------------------------------
  !> Structure for keeping element vertex data

  type ElementVertex
    real(RNP) :: x(3)          !< coordinates
    integer   :: neighbor = 0  !< neighbor element, including virtual one
  end type ElementVertex

  !-----------------------------------------------------------------------------
  !> Cartesian mesh element

  type MeshElement
    integer             :: global_id = 0 !< global ID
    type(ElementFace)   :: face(6)       !< face ID, neighbor and boundary info
    type(ElementEdge)   :: edge(12)      !< edge neighbor
    type(ElementVertex) :: vertex(8)     !< vertex coordinates and neighbor
  end type MeshElement

!===============================================================================

end module CART__Mesh_Element
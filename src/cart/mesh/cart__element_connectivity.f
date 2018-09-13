!> summary:  Structures for temporary storage of element connectivity
!> author:   Joerg Stiller
!> date:     2017/04/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Structures for temporary storage of element connectivity
!===============================================================================

module CART__Element_Connectivity
  private

  public :: ElementConnectivity

  !-----------------------------------------------------------------------------
  !> Structure for identifying a neighbor element

  type NeighborElement
    integer :: part = -1  !< host partition
    integer :: id   = -1  !< element ID within host partition
  end type NeighborElement

  !-----------------------------------------------------------------------------
  !> Structure for storing the connectivity of a Cartesian mesh element

  type ElementConnectivity
    type(NeighborElement) :: face(6)    !< face neighbors
    type(NeighborElement) :: edge(12)   !< edge neighbors
    type(NeighborElement) :: vertex(8)  !< vertex neighbors
  end type ElementConnectivity

!===============================================================================

end module CART__Element_Connectivity

!> summary:  Cartesian mesh link
!> author:   Joerg Stiller
!> date:     2017/04/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Cartesian mesh link
!===============================================================================

module CART__Mesh_Link
  use Kind_Parameters, only: IXS
  implicit none
  private

  public :: MeshLink
  public :: ElementLink

  !-----------------------------------------------------------------------------
  !> Structure for keeping element linking information

  type ElementLink
    integer      :: id = 0        !< local (virtual) element ID
    integer(IXS) :: face(6)   = 0 !< set to 1(0) for (un)linked faces
    integer(IXS) :: edge(12)  = 0 !< set to 1(0) for (un)linked edges
    integer(IXS) :: vertex(8) = 0 !< set to 1(0) for (un)linked vertices
  end type ElementLink

  !-----------------------------------------------------------------------------
  !> Structure for keeping mesh coupling information
  !>
  !> We maintain a list of master elements and a list of ghost elements.
  !> The former belong to the local partition, whereas the latter represent
  !> virtual entities for keeping data that refers to remote elements.

  type MeshLink
    integer :: part = -1                        !< remote partition ID
    integer :: nf = 0                           !< number of linked faces
    integer :: nm = 0                           !< number of master elements
    integer :: ng = 0                           !< number of ghost elements
    integer, allocatable :: face(:)             !< list of linked faces
    type(ElementLink), allocatable :: master(:) !< master elements
    type(ElementLink), allocatable :: ghost(:)  !< ghost elements
  end type MeshLink

!===============================================================================

end module CART__Mesh_Link

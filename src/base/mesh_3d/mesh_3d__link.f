!> summary:  3D mesh link
!> author:   Joerg Stiller
!> date:     2020/11/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_3d__Link
  use Kind_Parameters, only: IXS
  implicit none
  private

  public :: Mesh3d_Link
  public :: Mesh3d_ElementLink

  !-----------------------------------------------------------------------------
  !> Element linking information
  !>
  !> The structure may refer to a local element or a virtual "ghost" element.
  !> The components `vertex`, `edge` and `face` indicate element components,
  !> for which data needs to be transferred. This specification allows for
  !> implementing partial transfers moving only a selected subset of the given
  !> element data.

  type Mesh3d_ElementLink
    integer      :: id = 0        !< local/ghost element ID
    integer(IXS) :: vertex(8) = 0 !< set to 1(0) for (un)linked vertices
    integer(IXS) :: edge(12)  = 0 !< set to 1(0) for (un)linked edges
    integer(IXS) :: face(6)   = 0 !< set to 1(0) for (un)linked faces
  end type Mesh3d_ElementLink

  !-----------------------------------------------------------------------------
  !> Structure for keeping information for coupling mesh partitions
  !>
  !> Provides the following lists:
  !>
  !>   - `face`: the IDs of faces linked with partition `part`
  !>   - `coupled_face`: the IDs of local faces linked to those in `face`
  !>   - `master`: link to local elements with a ghost in partition `part`
  !>   - `ghost` : link to ghost elements with their master in partition `part`
  !>
  !> Note that `coupled_face` is empty unless `part` refers to the local
  !> partition.

  type Mesh3d_Link
    integer :: part     = -1                 !< remote partition ID
    integer :: n_face   =  0                 !< number of linked faces
    integer :: n_master =  0                 !< number of master elements
    integer :: n_ghost  =  0                 !< number of ghost elements
    integer, allocatable :: face(:)          !< list of linked faces
    integer, allocatable :: coupled_face(:)  !< list of coupled local faces
    type(Mesh3d_ElementLink), allocatable :: master(:) !< master elements
    type(Mesh3d_ElementLink), allocatable :: ghost(:)  !< ghost elements
  end type Mesh3d_Link

  !=============================================================================

end module Mesh_3d__Link

!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Numbering and identification of element components
!> author:   Joerg Stiller
!> date:     2017/04/30; revised 2020/11/21
!>
!> Hexahedral elements can be divided into different components
!>
!>   -  1  volume
!>   -  6  faces
!>   -  12 edges
!>   –  8  vertices
!>
!> Each component can be identified by the
!>
!>   -  index within its group
!>         * face:   1 ≤ f ≤ 6
!>         * edge:   1 ≤ e ≤ 12
!>         * vertex: 1 ≤ v ≤ 8
!>
!>   –  element component index l, 0 ≤ l ≤ 26
!>         * center: l = 0
!>         * face:   l = f
!>         * edge:   l = e + 6
!>         * face:   l = v + 18
!>
!>   –  triple position index (i,j,k), with -1 ≤ i,j,k ≤ 1
!>
!> Identifying (i,j,k) as a point in the standard element [-1,1]³ allows to
!> associate the triple positition index with the element component containing
!> this point. For example, (0,0,1) selects the face located at ζ=1.
!> More generally, the triplet corresponds to
!>
!>   –  the volume, if i=j=k=0,
!>   –  a face,     if containing 2 zeros,
!>   –  an edge,    if containing 1 zero, and
!>   –  a vertex,   if i,j,k ≠ 0
!>
!> These indices are used to establish mappings between the element components.
!> For example, V_FACE(1:4,f) lists the vertices v on face f. Within these
!> mappings, the subordinate components are ordered according to their linear
!> group index.
!>
!===============================================================================

module Mesh_Element_Indexing__3D
  use Kind_Parameters, only: IXS
  implicit none
  private

  !-----------------------------------------------------------------------------
  ! Public procedures

  public :: ElementComponentIndex
  public :: ElementComponentType
  public :: ElementFaceID
  public :: ElementEdgeID
  public :: ElementVertexID

  !-----------------------------------------------------------------------------
  !> `F_EDGE(:,e)` -- indices of element faces sharing edge `e`

  integer, parameter, public ::                               &
      F_EDGE(2,12) = reshape( [ 3,5, 4,5, 3,6, 4,6,           &
                                1,5, 2,5, 1,6, 2,6,           &
                                1,3, 2,3, 1,4, 2,4 ],  [2,12] )

  !-----------------------------------------------------------------------------
  !> `V_EDGE(i,e)` -- index `v` of vertex `i` of edge `e`

  integer, parameter, public ::                               &
      V_EDGE(2,12) = reshape( [ 1,2, 3,4, 5,6, 7,8,           &
                                1,3, 2,4, 5,7, 6,8,           &
                                1,5, 2,6, 3,7, 4,8 ],  [2,12] )

  !-----------------------------------------------------------------------------
  !> `V_FACE(i,f)` -- index `v` of vertex `i` of face `f`

  integer, parameter, public ::                               &
      V_FACE(4,6)  = reshape( [ 1, 3, 5, 7,                   &
                                2, 4, 6, 8,                   &
                                1, 2, 5, 6,                   &
                                3, 4, 7, 8,                   &
                                1, 2, 3, 4,                   &
                                5, 6, 7, 8 ],  [4,6] )

  !-----------------------------------------------------------------------------
  !> `E_FACE(i,f)` -- index `e` of edge `i` of face `f`

  integer, parameter, public ::                               &
      E_FACE(4,6)  = reshape( [ 5,  7,  9, 11,                &
                                6,  8, 10, 12,                &
                                1,  3,  9, 10,                &
                                2,  4, 11, 12,                &
                                1,  2,  5,  6,                &
                                3,  4,  7,  8 ],  [4,6] )

  !-----------------------------------------------------------------------------
  !> `F_VERT(:,v)` -- indices of element faces sharing vertex `v`

  integer, parameter, public ::                               &
      F_VERT(3,8)  = reshape( [ 1,3,5,  2,3,5,                &
                                1,4,5,  2,4,5,                &
                                1,3,6,  2,3,6,                &
                                1,4,6,  2,4,6 ],  [3,8] )

  !-----------------------------------------------------------------------------
  !> `E_VERT(:,v)` -- indices of element edges sharing vertex `v`

  integer, parameter, public ::                               &
      E_VERT(3,8)  = reshape( [ 1, 5,  9,  1, 6, 10,          &
                                2, 5, 11,  2, 6, 12,          &
                                3, 7,  9,  3, 8, 10,          &
                                4, 7, 11,  4, 8, 12 ],  [3,8] )

  !-----------------------------------------------------------------------------
  !> Element component index `l` associated with linear position index `p`

  integer, parameter, public ::         &
      L_POS(-1:1,-1:1,-1:1) = reshape( [ 19,     & ! -1, -1, -1 : vertex WSB
                                          7,     & !  0, -1, -1 : x1-edge SB
                                         20,     & !  1, -1, -1 : vertex ESB
                                         11,     & ! -1,  0, -1 : x2-edge WB
                                          5,     & !  0,  0, -1 : face B
                                         12,     & !  1,  0, -1 : x2-edge EB
                                         21,     & ! -1,  1, -1 : vertex WNB
                                          8,     & !  0,  1, -1 : x1-edge NB
                                         22,     & !  1,  1, -1 : vertex ENB
                                         15,     & ! -1, -1,  0 : x3-edge WS
                                          3,     & !  0, -1,  0 : face S
                                         16,     & !  1, -1,  0 : x3-edge ES
                                          1,     & ! -1,  0,  0 : face W
                                          0,     & !  0,  0,  0 : center
                                          2,     & !  1,  0,  0 : face E
                                         17,     & ! -1,  1,  0 : x3-edge WN
                                          4,     & !  0,  1,  0 : face N
                                         18,     & !  1,  1,  0 : x3-edge EN
                                         23,     & ! -1, -1,  1 : vertex WST
                                          9,     & !  0, -1,  1 : x1-edge ST
                                         24,     & !  1, -1,  1 : vertex EST
                                         13,     & ! -1,  0,  1 : x2-edge WT
                                          6,     & !  0,  0,  1 : face T
                                         14,     & !  1,  0,  1 : x2 edge ET
                                         25,     & ! -1,  1,  1 : vertex WNT
                                         10,     & !  0,  1,  1 : x1-edge NT
                                         26 ],   & !  1,  1,  1 : vertex ENT
                                         [3,3,3] )

  !-----------------------------------------------------------------------------
  ! Component identifiers

  integer, parameter, public :: IS_VOLUME = 0  !< volume/cell identifier
  integer, parameter, public :: IS_FACE   = 1  !< face identifier
  integer, parameter, public :: IS_EDGE   = 2  !< edge identifier
  integer, parameter, public :: IS_VERTEX = 3  !< vertex identifier

  !-----------------------------------------------------------------------------
  ! Generic interfaces

  interface ElementComponentType
    module procedure :: ElementComponentType_L_IDK
    module procedure :: ElementComponentType_L_IXS
    module procedure :: ElementComponentType_T
  end interface

  interface ElementFaceID
    module procedure :: ElementFaceID_L_IDK
    module procedure :: ElementFaceID_L_IXS
    module procedure :: ElementFaceID_T
  end interface

  interface ElementEdgeID
    module procedure :: ElementEdgeID_L_IDK
    module procedure :: ElementEdgeID_L_IXS
    module procedure :: ElementEdgeID_T
  end interface

  interface ElementVertexID
    module procedure :: ElementVertexID_L_IDK
    module procedure :: ElementVertexID_L_IXS
    module procedure :: ElementVertexID_T
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Converts triple position index (i,j,k) into linear element component index

  pure integer function ElementComponentIndex(i,j,k) result(l)
    integer, intent(in) :: i   !< position in direction 1,  -1 <= i <= 1
    integer, intent(in) :: j   !< position in direction 2,  -1 <= j <= 1
    integer, intent(in) :: k   !< position in direction 3,  -1 <= k <= 1

    if (abs(i) <= 1 .and. abs(j) <= 1 .and. abs(k) <= 1) then
      l = L_POS(i,j,k)
    else
      l = -1
    end if

  end function ElementComponentIndex

  !-----------------------------------------------------------------------------
  !> Returns the component type for given component index

  pure integer function ElementComponentType_L_IDK(l) result(typ)
    integer, intent(in) :: l !< component index

    select case(l)
    case(0)
      typ = IS_VOLUME
    case(1:6)
      typ = IS_FACE
    case(7:18)
      typ = IS_EDGE
    case default ! 19:26
      typ = IS_VERTEX
    end select

  end function ElementComponentType_L_IDK

  !-----------------------------------------------------------------------------
  !> Returns the component type for given component index

  pure integer function ElementComponentType_L_IXS(l) result(typ)
    integer(IXS), intent(in) :: l !< component index

    select case(l)
    case(0)
      typ = IS_VOLUME
    case(1:6)
      typ = IS_FACE
    case(7:18)
      typ = IS_EDGE
    case default ! 19:26
      typ = IS_VERTEX
    end select

  end function ElementComponentType_L_IXS

  !-----------------------------------------------------------------------------
  !> Returns the component type for given triple position index

  pure integer function ElementComponentType_T(i,j,k) result(typ)
    integer, intent(in) :: i !< position in direction 1,  -1 <= i <= 1
    integer, intent(in) :: j !< position in direction 2,  -1 <= j <= 1
    integer, intent(in) :: k !< position in direction 3,  -1 <= k <= 1

    typ = ElementComponentType(L_POS(i,j,k))

  end function ElementComponentType_T

  !-----------------------------------------------------------------------------
  !> Returns the element face `f` associated with component index `l`

  pure integer function ElementFaceID_L_IDK(l) result(f)
    integer, intent(in) :: l !< component index

    f = l

  end function ElementFaceID_L_IDK

  !-----------------------------------------------------------------------------
  !> Returns the element face `f` associated with component index `l`

  pure integer function ElementFaceID_L_IXS(l) result(f)
    integer(IXS), intent(in) :: l !< component index

    f = l

  end function ElementFaceID_L_IXS

  !-----------------------------------------------------------------------------
  !> Returns the element face `f` associated with position `(i,j,k)`

  pure integer function ElementFaceID_T(i,j,k) result(f)
    integer, intent(in) :: i !< position in direction 1,  -1 <= i <= 1
    integer, intent(in) :: j !< position in direction 2,  -1 <= j <= 1
    integer, intent(in) :: k !< position in direction 3,  -1 <= k <= 1

    f = L_POS(i,j,k)

  end function ElementFaceID_T

  !-----------------------------------------------------------------------------
  !> Returns the element edge `e` associated with component index `l`

  pure integer function ElementEdgeID_L_IDK(l) result(e)
    integer, intent(in) :: l !< component index

    e = l - 6

  end function ElementEdgeID_L_IDK

  !-----------------------------------------------------------------------------
  !> Returns the element edge `e` associated with component index `l`

  pure integer function ElementEdgeID_L_IXS(l) result(e)
    integer(IXS), intent(in) :: l !< component index

    e = l - 6

  end function ElementEdgeID_L_IXS

  !-----------------------------------------------------------------------------
  !> Returns the element edge `e` associated with position `(i,j,k)`

  pure integer function ElementEdgeID_T(i,j,k) result(e)
    integer, intent(in) :: i !< position in direction 1,  -1 <= i <= 1
    integer, intent(in) :: j !< position in direction 2,  -1 <= j <= 1
    integer, intent(in) :: k !< position in direction 3,  -1 <= k <= 1

    e = L_POS(i,j,k) - 6

  end function ElementEdgeID_T

  !-----------------------------------------------------------------------------
  !> Returns the element vertex `v` associated with component index `l`

  pure integer function ElementVertexID_L_IDK(l) result(v)
    integer, intent(in) :: l !< component index

    v = l - 18

  end function ElementVertexID_L_IDK

  !-----------------------------------------------------------------------------
  !> Returns the element vertex `v` associated with component index `l`

  pure integer function ElementVertexID_L_IXS(l) result(v)
    integer(IXS), intent(in) :: l !< component index

    v = l - 18

  end function ElementVertexID_L_IXS

  !-----------------------------------------------------------------------------
  !> Returns the element vertex `v` associated with position `(i,j,k)`

  pure integer function ElementVertexID_T(i,j,k) result(v)
    integer, intent(in) :: i !< position in direction 1,  -1 <= i <= 1
    integer, intent(in) :: j !< position in direction 2,  -1 <= j <= 1
    integer, intent(in) :: k !< position in direction 3,  -1 <= k <= 1

    v = L_POS(i,j,k) - 18

  end function ElementVertexID_T

  !=============================================================================

end module Mesh_Element_Indexing__3D

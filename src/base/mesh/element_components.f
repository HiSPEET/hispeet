!> summary:  Numbering and identification of element components
!> author:   Joerg Stiller
!> date:     2017/04/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Element_Components
  implicit none
  private

  public :: ElementComponentIndex
  public :: ElementComponentType
  public :: ElementComponentID

  integer, parameter, public :: IS_VOLUME = 0
  integer, parameter, public :: IS_FACE   = 1
  integer, parameter, public :: IS_EDGE   = 2
  integer, parameter, public :: IS_VERTEX = 3

  interface ElementComponentType
    module procedure ElementComponentType_L
    module procedure ElementComponentType_T
  end interface

  interface ElementComponentID
    module procedure ElementComponentID_L
    module procedure ElementComponentID_T
  end interface

contains

!-------------------------------------------------------------------------------
!> Converts triple position index (i,j,k) into linear element component index l
!>
!>     l        |   component(s)
!>    --------- | -----------------------------------------------------
!>     0        |   element center
!>     1:6      |   faces: W, E, S, N, B, T
!>     7:10     |   x1-edges: SB, NB, ST, NT
!>     11:14    |   x2-edges: WB, EB, WT, ET
!>     15:18    |   x3-edges: WS, ES, WN, EN
!>     19:26    |   vertices: WSB, ESB, WNB, ENB, WST, EST, WNT, ENT

pure integer function ElementComponentIndex(i,j,k) result(l)
  integer, intent(in) :: i   !< position in direction 1,  -1 <= i <= 1
  integer, intent(in) :: j   !< position in direction 2,  -1 <= j <= 1
  integer, intent(in) :: k   !< position in direction 3,  -1 <= k <= 1

  integer, parameter :: pos(-1:1, -1:1, -1:1)   & !  i,  j,  k : component
                        = reshape( [ 19,        & ! -1, -1, -1 : vertex WSB
                                      7,        & !  0, -1, -1 : x1-edge SB
                                     20,        & !  1, -1, -1 : vertex ESB
                                     11,        & ! -1,  0, -1 : x2-edge WB
                                      5,        & !  0,  0, -1 : face B
                                     12,        & !  1,  0, -1 : x2-edge EB
                                     21,        & ! -1,  1, -1 : vertex WNB
                                      8,        & !  0,  1, -1 : x1-edge NB
                                     22,        & !  1,  1, -1 : vertex ENB
                                     15,        & ! -1, -1,  0 : x3-edge WS
                                      3,        & !  0, -1,  0 : face S
                                     16,        & !  1, -1,  0 : x3-edge ES
                                      1,        & ! -1,  0,  0 : face W
                                      0,        & !  0,  0,  0 : center
                                      2,        & !  1,  0,  0 : face E
                                     17,        & ! -1,  1,  0 : x3-edge WN
                                      4,        & !  0,  1,  0 : face N
                                     18,        & !  1,  1,  0 : x3-edge EN
                                     23,        & ! -1, -1,  1 : vertex WST
                                      9,        & !  0, -1,  1 : x1-edge ST
                                     24,        & !  1, -1,  1 : vertex EST
                                     13,        & ! -1,  0,  1 : x2-edge WT
                                      6,        & !  0,  0,  1 : face T
                                     14,        & !  1,  0,  1 : x2 edge ET
                                     25,        & ! -1,  1,  1 : vertex WNT
                                     10,        & !  0,  1,  1 : x1-edge NT
                                     26         & !  1,  1,  1 : vertex ENT
                                   ], [ 3,3,3 ] )

  if (abs(i) <= 1 .and. abs(j) <= 1 .and. abs(k) <= 1) then
    l = pos(i,j,k)
  else
    l = -1
  end if

end function ElementComponentIndex

!-------------------------------------------------------------------------------
!> Returns the component type for given linear component index

pure integer function ElementComponentType_L(l) result(typ)
  integer, intent(in) :: l !< linear element component index

  select case(l)
  case(0)
    typ = IS_VOLUME
  case(1:6)
    typ = IS_FACE
  case(7:18)
    typ = IS_EDGE
  case(19:26)
    typ = IS_VERTEX
  case default
    typ = -1
  end select

end function ElementComponentType_L

!-------------------------------------------------------------------------------
!> Returns the component type for given triple position index

pure integer function ElementComponentType_T(i,j,k) result(typ)
  integer, intent(in) :: i   !< position in direction 1,  -1 <= i <= 1
  integer, intent(in) :: j   !< position in direction 2,  -1 <= j <= 1
  integer, intent(in) :: k   !< position in direction 3,  -1 <= k <= 1

  typ = ElementComponentType_L( ElementComponentIndex(i,j,k) )

end function ElementComponentType_T

!-------------------------------------------------------------------------------
!> Returns the component ID for given linear component index

pure integer function ElementComponentID_L(l) result(id)
  integer, intent(in) :: l !< linear element component index

  select case(l)
  case(0)
    id = 1
  case(1:6)
    id = l
  case(7:18)
    id = l - 6
  case(19:26)
    id = l - 18
  case default
    id = -1
  end select

end function ElementComponentID_L

!-------------------------------------------------------------------------------
!> Returns the component ID for given triple position index

pure integer function ElementComponentID_T(i,j,k) result(id)
  integer, intent(in) :: i   !< position in direction 1,  -1 <= i <= 1
  integer, intent(in) :: j   !< position in direction 2,  -1 <= j <= 1
  integer, intent(in) :: k   !< position in direction 3,  -1 <= k <= 1

  id = ElementComponentID_L( ElementComponentIndex(i,j,k) )

end function ElementComponentID_T

!===============================================================================

end module Element_Components

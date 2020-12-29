!> summary:  Procedures for lexical ordering
!> author:   Joerg Stiller
!> date:     2013/05/24; revised 2014/11/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Structured_Mesh
  implicit none
  private

  public :: NumberOfVertices
  public :: NumberOfEdges
  public :: NumberOfFaces
  public :: NumberOfElements

  public :: LexicalIndex
  public :: LexicalVertexIndex
  public :: LexicalEdgeIndex
  public :: LexicalFaceIndex
  public :: LexicalElementIndex

  public :: TripleIndex
  public :: TripleVertexIndex
  public :: TripleEdgeIndex
  public :: TripleFaceIndex
  public :: TripleElementIndex

contains

!-------------------------------------------------------------------------------
!> Number of Vertices.

pure integer function NumberOfVertices(n1, n2, n3, nb) result(n)
  integer, intent(in) :: n1  !< number of intervals in direction 1
  integer, intent(in) :: n2  !< number of intervals in direction 2
  integer, intent(in) :: n3  !< number of intervals in direction 3
  integer, intent(in) :: nb  !< number of buffer layers

  integer :: n1b, n2b, n3b

  n1b = n1 + 2*nb
  n2b = n2 + 2*nb
  n3b = n3 + 2*nb

  n = (n1b + 1) * (n2b + 1) * (n3b + 1)

end function NumberOfVertices

!-------------------------------------------------------------------------------
!> Calculates the number of edges in a structured mesh.
!>
!> This function returns the number of edges in a structured mesh depending on
!> the number of intervals used in each direction and the number of buffer
!> layers.

pure integer function NumberOfEdges(n1, n2, n3, nb) result(n)
  integer, intent(in) :: n1  !< number of intervals in direction 1
  integer, intent(in) :: n2  !< number of intervals in direction 2
  integer, intent(in) :: n3  !< number of intervals in direction 3
  integer, intent(in) :: nb  !< number of buffer layers

  integer :: n1b, n2b, n3b

  ! number of elements in each direction
  n1b = n1 + 2*nb
  n2b = n2 + 2*nb
  n3b = n3 + 2*nb

  ! The number of edges parallel to the first direction is equal to the number
  ! of vertices in a plane times the number of intervals in the first direction,
  ! the same applies to the other directions.
  n =  n1b      * (n2b + 1) * (n3b + 1)  &
    + (n1b + 1) *  n2b      * (n3b + 1)  &
    + (n1b + 1) * (n2b + 1) *  n3b

end function NumberOfEdges

!-------------------------------------------------------------------------------
!> Calculates the number of faces in the mesh.
!>
!> This function calculates the number of faces of the structured mesh depending
!> on the number of intervals in each direction and the number of buffer layers.
!>
!> As the mesh is structured, the number of faces in one direction equals the
!> number of vertices in that direction times the number of faces in one plane.
!> Summing up over all three directions gives the total number of faces.

pure integer function NumberOfFaces(n1, n2, n3, nb) result(n)
  integer, intent(in) :: n1  !< number of intervals in direction 1
  integer, intent(in) :: n2  !< number of intervals in direction 2
  integer, intent(in) :: n3  !< number of intervals in direction 3
  integer, intent(in) :: nb  !< number of buffer layers

  integer :: n1b, n2b, n3b

  ! number of elements in each direction
  n1b = n1 + 2*nb
  n2b = n2 + 2*nb
  n3b = n3 + 2*nb

  ! number of faces
  n = (n1b + 1) *  n2b      *  n3b       & ! x1-faces
    +  n1b      * (n2b + 1) *  n3b       & ! x2-faces
    +  n1b      *  n2b      * (n3b + 1)    ! x3-faces

end function NumberOfFaces

!-------------------------------------------------------------------------------
!> Calculates the number of elements in a structured mesh.
!>
!> This function returns the number of elements in a structured mesh depending
!> on the number of intervals used in each direction and the number of buffer
!> layers.

pure integer function NumberOfElements(n1, n2, n3, nb) result(n)
  integer, intent(in) :: n1  !< number of intervals in direction 1
  integer, intent(in) :: n2  !< number of intervals in direction 2
  integer, intent(in) :: n3  !< number of intervals in direction 3
  integer, intent(in) :: nb  !< number of buffer layers

  n = (n1 + 2*nb) * (n2 + 2*nb) * (n3 + 2*nb)

end function NumberOfElements

!-------------------------------------------------------------------------------
!> Returns the lexical index of a mesh component.
!>
!> Given the triple index (i,j,k) the function returns the single-valued lexical
!> index l starting from One. The directional indices i, j and k start from i0,
!> j0 and k0, respectively. The number of items is n1 in the first direction, n2
!> in the second, and unspecified in the third direction.

pure integer function LexicalIndex(i, j, k, n1, n2, i0, j0, k0) result(l)
  integer, intent(in) :: i  !< index in direction 1
  integer, intent(in) :: j  !< index in direction 2
  integer, intent(in) :: k  !< index in direction 3
  integer, intent(in) :: n1 !< number of items in direction 1
  integer, intent(in) :: n2 !< number of items in direction 2
  integer, intent(in) :: i0 !< start index in direction 1
  integer, intent(in) :: j0 !< start index in direction 2
  integer, intent(in) :: k0 !< start index in direction 3

  l = 1 + (i-i0) + n1*(j-j0) + n1*n2*(k-k0)

end function LexicalIndex

!-------------------------------------------------------------------------------
!> Recovers the triple index of a from the given lexical index.
!>
!> Converts the lexical index l into the corresponding triple index i, j, k for
!> given start indices i0, j0, k0 and item numbers n1, n2 in the first two
!> directions.

pure subroutine TripleIndex(i, j, k, n1, n2, i0, j0, k0, l)
  integer, intent(out) :: i   !< index in direction 1
  integer, intent(out) :: j   !< index in direction 2
  integer, intent(out) :: k   !< index in direction 3
  integer, intent(in)  :: n1  !< number of items in direction 1
  integer, intent(in)  :: n2  !< number of items in direction 2
  integer, intent(in)  :: i0  !< start index in direction 1
  integer, intent(in)  :: j0  !< start index in direction 2
  integer, intent(in)  :: k0  !< start index in direction 3
  integer, intent(in)  :: l   !< lexical index

  k = k0 + (l - 1) / (n1*n2)
  j = j0 + (l - 1 - n1*n2 * (k - k0)) / n1
  i = i0 +  l - 1 - n1 * (j - j0) - n1*n2 * (k - k0)

end subroutine TripleIndex

!-------------------------------------------------------------------------------
!> Lexical index of a vertex.

pure integer function LexicalVertexIndex(i, j, k, n1, n2, nb) result(l)
  integer, intent(in) :: i   !< vertex index in direction 1
  integer, intent(in) :: j   !< vertex index in direction 2
  integer, intent(in) :: k   !< vertex index in direction 3
  integer, intent(in) :: n1  !< number of intervals in direction 1
  integer, intent(in) :: n2  !< number of intervals in direction 2
  integer, intent(in) :: nb  !< number of buffer layers

  l = LexicalIndex(i, j, k, n1+2*nb+1 , n2+2*nb+1, -nb, -nb, -nb)

end function LexicalVertexIndex

!-------------------------------------------------------------------------------
!> Triple index of a vertex.

pure subroutine TripleVertexIndex(i, j, k, n1, n2, nb, l)
  integer, intent(out) :: i   !< vertex index in direction 1
  integer, intent(out) :: j   !< vertex index in direction 2
  integer, intent(out) :: k   !< vertex index in direction 3
  integer, intent(in)  :: n1  !< number of intervals in direction 1
  integer, intent(in)  :: n2  !< number of intervals in direction 2
  integer, intent(in)  :: nb  !< number of buffer layers
  integer, intent(in)  :: l   !< lexical index

  call TripleIndex(i, j, k, n1+2*nb+1 , n2+2*nb+1, -nb, -nb, -nb, l)

end subroutine TripleVertexIndex

!-------------------------------------------------------------------------------
!> Lexical index of an edge.

pure integer function LexicalEdgeIndex(i, j, k, d, n1, n2, n3, nb) result(l)
  integer, intent(in) :: i   !< vertex index in direction 1
  integer, intent(in) :: j   !< vertex index in direction 2
  integer, intent(in) :: k   !< vertex index in direction 3
  integer, intent(in) :: d   !< direction (1,2, or 3)
  integer, intent(in) :: n1  !< number of intervals in direction 1
  integer, intent(in) :: n2  !< number of intervals in direction 2
  integer, intent(in) :: n3  !< number of intervals in direction 3
  integer, intent(in) :: nb  !< number of buffer layers

  integer :: n1b, n2b, n3b

  n1b = n1 + 2*nb
  n2b = n2 + 2*nb
  n3b = n3 + 2*nb

  select case(d)
  case(1)
    l = LexicalIndex(i, j, k, n1b, n2b+1, -nb, -nb, -nb)
  case(2)
    l = LexicalIndex(i, j, k, n1b+1, n2b, -nb, -nb, -nb)  &
      + n1b * (n2b + 1) * (n3b + 1)
  case(3)
    l = LexicalIndex(i, j, k, n1b+1, n2b+1, -nb, -nb, -nb)  &
      + n1b * (n2b + 1) * (n3b + 1)  &
      + n2b * (n3b + 1) * (n1b + 1)
  case default
    l = -1
  end select

end function LexicalEdgeIndex

!-------------------------------------------------------------------------------
!> Triple index of an edge.

pure subroutine TripleEdgeIndex(i, j, k, d, n1, n2, n3, nb, l)
  integer, intent(out) :: i   !< vertex index in direction 1
  integer, intent(out) :: j   !< vertex index in direction 2
  integer, intent(out) :: k   !< vertex index in direction 3
  integer, intent(out) :: d   !< orientation (1,2, or 3)
  integer, intent(in)  :: n1  !< number of intervals in direction 1
  integer, intent(in)  :: n2  !< number of intervals in direction 2
  integer, intent(in)  :: n3  !< number of intervals in direction 3
  integer, intent(in)  :: nb  !< number of buffer layers
  integer, intent(in)  :: l   !< lexical index

  integer :: n1b, n2b, n3b, l1, l2

  n1b = n1 + 2*nb
  n2b = n2 + 2*nb
  n3b = n3 + 2*nb

  l1 = n1b * (n2b + 1) * (n3b + 1)        ! number of x1-edges
  l2 = n2b * (n3b + 1) * (n1b + 1) + l1   ! number of x1- and x2-edges

  if (l > l2) then ! x3-edge
    d = 3
    call TripleIndex(i, j, k, n1b+1, n2b+1, -nb, -nb, -nb, l-l2)
  else if (l > l1) then ! x2-edge
    d = 2
    call TripleIndex(i, j, k, n1b+1, n2b, -nb, -nb, -nb, l-l1)
  else ! x1-edge
    d = 1
    call TripleIndex(i, j, k, n1b, n2b+1, -nb, -nb, -nb, l)
  end if

end subroutine TripleEdgeIndex

!-------------------------------------------------------------------------------
!> Lexical index of a face.
!>
!> The face are grouped according to their normal direction. Lexical numbering
!> is applied within each group. Numbering starts at 1 with faces normal to
!> direction 1, continues with direction 2 and, finally, direction 3.
!> Thus, a unique ID is assigned to every mesh face, while retaining the natural
!> ordering and easy accessibility of individual directions.

pure integer function LexicalFaceIndex(i, j, k, d, n1, n2, n3, nb) result(l)
  integer, intent(in) :: i   !< vertex index in direction 1
  integer, intent(in) :: j   !< vertex index in direction 2
  integer, intent(in) :: k   !< vertex index in direction 3
  integer, intent(in) :: d   !< normal orientation (1,2, or 3)
  integer, intent(in) :: n1  !< number of intervals in direction 1
  integer, intent(in) :: n2  !< number of intervals in direction 2
  integer, intent(in) :: n3  !< number of intervals in direction 3
  integer, intent(in) :: nb  !< number of buffer layers

  integer :: n1b, n2b, n3b

  n1b = n1 + 2*nb
  n2b = n2 + 2*nb
  n3b = n3 + 2*nb

  select case(d)
  case(1)
    l = LexicalIndex(i, j, k, n1b+1, n2b, -nb, -nb, -nb)
  case(2)
    l = LexicalIndex(i, j, k, n1b, n2b+1, -nb, -nb, -nb)  &
      + (n1b + 1) * n2b * n3b
  case(3)
    l = LexicalIndex(i, j, k, n1b, n2b, -nb, -nb, -nb)  &
      + (n1b + 1) * n2b * n3b  &
      + (n2b + 1) * n3b * n1b
  case default
    l = -1
  end select

end function LexicalFaceIndex

!-------------------------------------------------------------------------------
!> Triple index of a face.

pure subroutine TripleFaceIndex(i, j, k, d, n1, n2, n3, nb, l)
  integer, intent(out) :: i   !< vertex index in direction 1
  integer, intent(out) :: j   !< vertex index in direction 2
  integer, intent(out) :: k   !< vertex index in direction 3
  integer, intent(out) :: d   !< normal orientation (1,2, or 3)
  integer, intent(in)  :: n1  !< number of intervals in direction 1
  integer, intent(in)  :: n2  !< number of intervals in direction 2
  integer, intent(in)  :: n3  !< number of intervals in direction 3
  integer, intent(in)  :: nb  !< number of buffer layers
  integer, intent(in)  :: l   !< lexical index

  integer :: n1b, n2b, n3b, l1, l2

  n1b = n1 + 2*nb
  n2b = n2 + 2*nb
  n3b = n3 + 2*nb

  l1 = (n1b + 1) * n2b * n3b       ! number of x1-faces
  l2 = (n2b + 1) * n3b * n1b + l1  ! number of x1- and x2-faces

  if (l > l2) then ! x3-face
    d = 3
    call TripleIndex(i, j, k, n1b, n2b, -nb, -nb, -nb, l-l2)
  else if (l > l1) then ! x2-face
    d = 2
    call TripleIndex(i, j, k, n1b, n2b+1, -nb, -nb, -nb, l-l1)
  else ! x1-face
    d = 1
    call TripleIndex(i, j, k, n1b+1, n2b, -nb, -nb, -nb, l)
  end if

end subroutine TripleFaceIndex

!-------------------------------------------------------------------------------
!> Lexical index of an element.

pure integer function LexicalElementIndex(i, j, k, n1, n2, nb) result(l)
  integer, intent(in) :: i   !< element index in direction 1
  integer, intent(in) :: j   !< element index in direction 2
  integer, intent(in) :: k   !< element index in direction 3
  integer, intent(in) :: n1  !< number of intervals in direction 1
  integer, intent(in) :: n2  !< number of intervals in direction 2
  integer, intent(in) :: nb  !< number of buffer layers

  l = LexicalIndex(i, j, k, n1+2*nb, n2+2*nb, i0=1-nb, j0=1-nb, k0=1-nb)

end function LexicalElementIndex

!-------------------------------------------------------------------------------
!> Triple index of an element.

pure subroutine TripleElementIndex(i, j, k, n1, n2, nb, l)
  integer, intent(out) :: i   !< element index in direction 1
  integer, intent(out) :: j   !< element index in direction 2
  integer, intent(out) :: k   !< element index in direction 3
  integer, intent(in)  :: n1  !< number of intervals in direction 1
  integer, intent(in)  :: n2  !< number of intervals in direction 2
  integer, intent(in)  :: nb  !< number of buffer layers
  integer, intent(in)  :: l   !< lexical index

  call TripleIndex(i, j, k, n1+2*nb, n2+2*nb, i0=1-nb, j0=1-nb, k0=1-nb, l=l)

end subroutine TripleElementIndex

!===============================================================================

end module Structured_Mesh

!> summary:  Lexical ordering applied to 3d structured meshes
!> author:   Joerg Stiller
!> date:     2013/05/24; revised 2020/11/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> This modules provides procedures for lexical numbering of structured meshes.
!> Every mesh component is identified by a triple index, which is supplemented
!> its orientation in case of edges and faces. For each category of components
!> the triple index is mapped to a unique linear index starting at 1.
!>
!> The range of the triple index `(i,j,k)` varies from component to component:
!>
!>   - for elements
!>       *  `1 ≤ i ≤ n1`
!>       *  `1 ≤ j ≤ n2`
!>       *  `1 ≤ k ≤ n3`
!>
!>   - for faces with their normal pointing in direction 1
!>       *  `0 ≤ i ≤ n1`
!>       *  `1 ≤ j ≤ n2`
!>       *  `1 ≤ k ≤ n3`
!>
!>   - for faces with their normal pointing in direction 2
!>       *  `1 ≤ i ≤ n1`
!>       *  `0 ≤ j ≤ n2`
!>       *  `1 ≤ k ≤ n3`
!>
!>   - for faces with their normal pointing in direction 3
!>       *  `1 ≤ i ≤ n1`
!>       *  `1 ≤ j ≤ n2`
!>       *  `0 ≤ k ≤ n3`
!>
!>   - for edges aligned with direction 1
!>       *  `1 ≤ i ≤ n1`
!>       *  `0 ≤ j ≤ n2`
!>       *  `0 ≤ k ≤ n3`
!>
!>   - for edges aligned with direction 2
!>       *  `0 ≤ i ≤ n1`
!>       *  `1 ≤ j ≤ n2`
!>       *  `0 ≤ k ≤ n3`
!>
!>   - for edges aligned with direction 3
!>       *  `0 ≤ i ≤ n1`
!>       *  `0 ≤ j ≤ n2`
!>       *  `1 ≤ k ≤ n3`
!>
!>   - for vertices
!>       *  `0 ≤ i ≤ n1`
!>       *  `0 ≤ j ≤ n2`
!>       *  `0 ≤ k ≤ n3`
!>
!> The linear indices of faces and edges depend on their orientation. Numbering
!> starts with components aligned with direction 1, followed by direction 2 and,
!> finally, direction 3.
!>
!===============================================================================

module Mesh_3d__Structured_Indexing
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

  !-----------------------------------------------------------------------------
  !> Number of Vertices.

  pure integer function NumberOfVertices(n1, n2, n3) result(n)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3

    n = (n1 + 1) * (n2 + 1) * (n3 + 1)

  end function NumberOfVertices

  !-----------------------------------------------------------------------------
  !> Calculates the number of edges in a structured mesh.
  !>
  !> This function returns the number of edges in a structured mesh depending on
  !> the number of intervals used in each direction.

  pure integer function NumberOfEdges(n1, n2, n3) result(n)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3

    ! The number of edges parallel to the first direction is equal to the number
    ! of vertices in a plane times the number of intervals in the first direction,
    ! the same applies to the other directions.
    n =  n1      * (n2 + 1) * (n3 + 1)  &
      + (n1 + 1) *  n2      * (n3 + 1)  &
      + (n1 + 1) * (n2 + 1) *  n3

  end function NumberOfEdges

  !-----------------------------------------------------------------------------
  !> Calculates the number of faces in the mesh.
  !>
  !> This function calculates the number of faces of the structured mesh
  !> depending on the number of intervals in each direction.
  !>
  !> As the mesh is structured, the number of faces in one direction equals
  !> the number of vertices in that direction times the number of faces in
  !> one plane. Summing up over all three directions gives the total number
  !> of faces.

  pure integer function NumberOfFaces(n1, n2, n3) result(n)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3

    n = (n1 + 1) *  n2      *  n3       & ! x1-faces
      +  n1      * (n2 + 1) *  n3       & ! x2-faces
      +  n1      *  n2      * (n3 + 1)    ! x3-faces

  end function NumberOfFaces

  !-----------------------------------------------------------------------------
  !> Calculates the number of elements in a structured mesh.
  !>
  !> This function returns the number of elements in a structured mesh depending
  !> on the number of intervals used in each direction.

  pure integer function NumberOfElements(n1, n2, n3) result(n)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3

    n = n1 * n2 * n3

  end function NumberOfElements

  !-----------------------------------------------------------------------------
  !> Returns the lexical index of a mesh component.
  !>
  !> Given the triple index (i,j,k) the function returns the single-valued
  !> lexical index l starting from 1. The directional indices i, j and k
  !> start from i0, j0 and k0, respectively. The number of items is n1 in the
  !> first direction, n2 in the second, and unspecified in the third direction.

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

  !-----------------------------------------------------------------------------
  !> Recovers the triple index of a from the given lexical index.
  !>
  !> Converts the lexical index l into the corresponding triple index i, j, k
  !> for given start indices i0, j0, k0 and item numbers n1, n2 in the first two
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

  !-----------------------------------------------------------------------------
  !> Lexical index of a vertex.

  pure integer function LexicalVertexIndex(i, j, k, n1, n2) result(l)
    integer, intent(in) :: i   !< vertex index in direction 1
    integer, intent(in) :: j   !< vertex index in direction 2
    integer, intent(in) :: k   !< vertex index in direction 3
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2

    l = LexicalIndex(i, j, k, n1+1 , n2+1, 0, 0, 0)

  end function LexicalVertexIndex

  !-----------------------------------------------------------------------------
  !> Triple index of a vertex.

  pure subroutine TripleVertexIndex(i, j, k, n1, n2, l)
    integer, intent(out) :: i   !< vertex index in direction 1
    integer, intent(out) :: j   !< vertex index in direction 2
    integer, intent(out) :: k   !< vertex index in direction 3
    integer, intent(in)  :: n1  !< number of intervals in direction 1
    integer, intent(in)  :: n2  !< number of intervals in direction 2
    integer, intent(in)  :: l   !< lexical index

    call TripleIndex(i, j, k, n1+1 , n2+1, 0, 0, 0, l)

  end subroutine TripleVertexIndex

  !-----------------------------------------------------------------------------
  !> Lexical index of an edge.

  pure integer function LexicalEdgeIndex(i, j, k, d, n1, n2, n3) result(l)
    integer, intent(in) :: i   !< vertex index in direction 1
    integer, intent(in) :: j   !< vertex index in direction 2
    integer, intent(in) :: k   !< vertex index in direction 3
    integer, intent(in) :: d   !< direction (1, 2, or 3)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3

    select case(d)
    case(1)
      l = LexicalIndex(i, j, k, n1, n2+1, 1, 0, 0)
    case(2)
      l = LexicalIndex(i, j, k, n1+1, n2, 0, 1, 0)  &
        + n1 * (n2 + 1) * (n3 + 1)
    case(3)
      l = LexicalIndex(i, j, k, n1+1, n2+1, 0, 0, 1)  &
        + n1 * (n2 + 1) * (n3 + 1)  &
        + n2 * (n3 + 1) * (n1 + 1)
    case default
      l = -1
    end select

  end function LexicalEdgeIndex

  !-----------------------------------------------------------------------------
  !> Triple index of an edge.

  pure subroutine TripleEdgeIndex(i, j, k, d, n1, n2, n3, l)
    integer, intent(out) :: i   !< vertex index in direction 1
    integer, intent(out) :: j   !< vertex index in direction 2
    integer, intent(out) :: k   !< vertex index in direction 3
    integer, intent(out) :: d   !< orientation (1, 2, or 3)
    integer, intent(in)  :: n1  !< number of intervals in direction 1
    integer, intent(in)  :: n2  !< number of intervals in direction 2
    integer, intent(in)  :: n3  !< number of intervals in direction 3
    integer, intent(in)  :: l   !< lexical index

    integer :: l1, l2

    l1 = n1 * (n2 + 1) * (n3 + 1)        ! number of x1-edges
    l2 = n2 * (n3 + 1) * (n1 + 1) + l1   ! number of x1- and x2-edges

    if (l > l2) then ! x3-edge
      d = 3
      call TripleIndex(i, j, k, n1+1, n2+1, 0, 0, 1, l-l2)
    else if (l > l1) then ! x2-edge
      d = 2
      call TripleIndex(i, j, k, n1+1, n2  , 0, 1, 0, l-l1)
    else ! x1-edge
      d = 1
      call TripleIndex(i, j, k, n1  , n2+1, 1, 0, 0, l)
    end if

  end subroutine TripleEdgeIndex

  !-----------------------------------------------------------------------------
  !> Lexical index of a face.
  !>
  !> The face are grouped according to their normal direction. Lexical numbering
  !> is applied within each group. Numbering starts at 1 with faces normal to
  !> direction 1, continues with direction 2 and, finally, direction 3.
  !> Thus, a unique ID is assigned to every mesh face, while retaining the
  !> natural ordering and easy accessibility of individual directions.

  pure integer function LexicalFaceIndex(i, j, k, d, n1, n2, n3) result(l)
    integer, intent(in) :: i   !< vertex index in direction 1
    integer, intent(in) :: j   !< vertex index in direction 2
    integer, intent(in) :: k   !< vertex index in direction 3
    integer, intent(in) :: d   !< normal orientation (1,2, or 3)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3

    select case(d)
    case(1)
      l = LexicalIndex(i, j, k, n1+1, n2  , 0, 1, 1)
    case(2)
      l = LexicalIndex(i, j, k, n1  , n2+1, 1, 0, 1)  &
        + (n1 + 1) * n2 * n3
    case(3)
      l = LexicalIndex(i, j, k, n1  , n2  , 1, 1, 0)  &
        + (n1 + 1) * n2 * n3  &
        + (n2 + 1) * n3 * n1
    case default
      l = -1
    end select

  end function LexicalFaceIndex

  !-----------------------------------------------------------------------------
  !> Triple index of a face.

  pure subroutine TripleFaceIndex(i, j, k, d, n1, n2, n3, l)
    integer, intent(out) :: i   !< vertex index in direction 1
    integer, intent(out) :: j   !< vertex index in direction 2
    integer, intent(out) :: k   !< vertex index in direction 3
    integer, intent(out) :: d   !< normal orientation (1,2, or 3)
    integer, intent(in)  :: n1  !< number of intervals in direction 1
    integer, intent(in)  :: n2  !< number of intervals in direction 2
    integer, intent(in)  :: n3  !< number of intervals in direction 3
    integer, intent(in)  :: l   !< lexical index

    integer :: l1, l2

    l1 = (n1 + 1) * n2 * n3       ! numer of x1-faces
    l2 = (n2 + 1) * n3 * n1 + l1  ! numer of x1- and x2-faces

    if (l > l2) then ! x3-face
      d = 3
      call TripleIndex(i, j, k, n1  , n2  , 1, 1, 0, l-l2)
    else if (l > l1) then ! x2-face
      d = 2
      call TripleIndex(i, j, k, n1  , n2+1, 1, 0, 1, l-l1)
    else ! x1-face
      d = 1
      call TripleIndex(i, j, k, n1+1, n2  , 0, 1, 1, l)
    end if

  end subroutine TripleFaceIndex

  !-----------------------------------------------------------------------------
  !> Lexical index of an element.

  pure integer function LexicalElementIndex(i, j, k, n1, n2) result(l)
    integer, intent(in) :: i   !< element index in direction 1
    integer, intent(in) :: j   !< element index in direction 2
    integer, intent(in) :: k   !< element index in direction 3
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2

    l = LexicalIndex(i, j, k, n1, n2, i0=1, j0=1, k0=1)

  end function LexicalElementIndex

  !-----------------------------------------------------------------------------
  !> Triple index of an element.

  pure subroutine TripleElementIndex(i, j, k, n1, n2, l)
    integer, intent(out) :: i   !< element index in direction 1
    integer, intent(out) :: j   !< element index in direction 2
    integer, intent(out) :: k   !< element index in direction 3
    integer, intent(in)  :: n1  !< number of intervals in direction 1
    integer, intent(in)  :: n2  !< number of intervals in direction 2
    integer, intent(in)  :: l   !< lexical index

    call TripleIndex(i, j, k, n1, n2, i0=1, j0=1, k0=1, l=l)

  end subroutine TripleElementIndex

  !=============================================================================

end module Mesh_3d__Structured_Indexing

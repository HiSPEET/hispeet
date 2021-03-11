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
  use Kind_Parameters, only: IXL
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

  interface LexicalElementIndex
    module procedure :: LexicalElementIndex_IDK
    module procedure :: LexicalElementIndex_IXL
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Number of Vertices.

  pure integer function NumberOfVertices(n1, n2, n3, periodic) result(n)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3
    logical, intent(in) :: periodic(3) !< T (F) for (non)periodic directions

    integer :: m(3)

    where (periodic)
      m = 0
    elsewhere
      m = 1
    end where

    n = (n1 + m(1)) * (n2 + m(2)) * (n3 + m(3))

  end function NumberOfVertices

  !-----------------------------------------------------------------------------
  !> Calculates the number of edges in a structured mesh.
  !>
  !> This function returns the number of edges in a structured mesh depending on
  !> the number of intervals used in each direction.

  pure integer function NumberOfEdges(n1, n2, n3, periodic) result(n)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3
    logical, intent(in) :: periodic(3) !< T (F) for (non)periodic directions

    integer :: m(3)

    where (periodic)
      m = 0
    elsewhere
      m = 1
    end where

    ! The number of edges parallel to the first direction is equal to the number
    ! of vertices in a plane times the number of intervals in the first direction,
    ! the same applies to the other directions.
    n =  n1         * (n2 + m(2)) * (n3 + m(3))  &
      + (n1 + m(1)) *  n2         * (n3 + m(3))  &
      + (n1 + m(1)) * (n2 + m(2)) *  n3

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

  pure integer function NumberOfFaces(n1, n2, n3, periodic) result(n)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3
    logical, intent(in) :: periodic(3) !< T (F) for (non)periodic directions

    integer :: m(3)

    where (periodic)
      m = 0
    elsewhere
      m = 1
    end where

    n = (n1 + m(1)) *  n2         *  n3          & ! x1-faces
      +  n1         * (n2 + m(2)) *  n3          & ! x2-faces
      +  n1         *  n2         * (n3 + m(3))    ! x3-faces

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
  !> Unified handling of indices for periodic and non-periodic cases

  pure integer function NormalizedIndex(m, n, periodic) result(i)
    integer, intent(in) :: m        !< index ranging from `0` to `n-1`
    integer, intent(in) :: n        !< number of items
    logical, intent(in) :: periodic !< T (F) for (non)periodic case

    if (periodic) then
      i = mod(m, max(n,1))
    else
      i = m
    end if

  end function NormalizedIndex

  !-----------------------------------------------------------------------------
  !> Returns the lexical index of a mesh component.
  !>
  !> Given the triple index (i,j,k) the function returns the single-valued
  !> lexical index l starting from 1. The directional indices i, j and k
  !> start from i0, j0 and k0, respectively. The number of items is n1 in the
  !> first direction, n2 in the second, and n3 in the third direction.

  pure integer function LexicalIndex(i, j, k, i0, j0, k0, n1, n2, n3, periodic) &
      result(l)

    integer          , intent(in) :: i  !< index in direction 1
    integer          , intent(in) :: j  !< index in direction 2
    integer          , intent(in) :: k  !< index in direction 3
    integer          , intent(in) :: i0 !< start index in direction 1
    integer          , intent(in) :: j0 !< start index in direction 2
    integer          , intent(in) :: k0 !< start index in direction 3
    integer          , intent(in) :: n1 !< number of items in direction 1
    integer          , intent(in) :: n2 !< number of items in direction 2
    integer, optional, intent(in) :: n3 !< number of items in direction 3
    logical, optional, intent(in) :: periodic(3) !< T for periodic directions [F]

    integer :: ii, jj, kk, m(2)

    if (present(n3) .and. present(periodic)) then
      ii = NormalizedIndex(i-i0, n1-1, periodic(1))
      jj = NormalizedIndex(j-j0, n2-1, periodic(2))
      kk = NormalizedIndex(k-k0, n3-1, periodic(3))
      where(periodic(1:2))
        m = -1
      elsewhere
        m  = 0
      end where
    else
      ii = i - i0
      jj = j - j0
      kk = k - k0
      m  = 0
    end if

    l = 1 + ii + (n1 + m(1)) * (jj + (n2 + m(2))*kk)

  end function LexicalIndex

  !-----------------------------------------------------------------------------
  !> Lexical index of a vertex.

  pure integer function LexicalVertexIndex(i, j, k, n1, n2, n3, periodic) &
      result(l)

    integer          , intent(in) :: i   !< vertex index in direction 1
    integer          , intent(in) :: j   !< vertex index in direction 2
    integer          , intent(in) :: k   !< vertex index in direction 3
    integer          , intent(in) :: n1  !< number of intervals in direction 1
    integer          , intent(in) :: n2  !< number of intervals in direction 2
    integer, optional, intent(in) :: n3  !< number of intervals in direction 3
    logical, optional, intent(in) :: periodic(3) !< T periodic directions [F]

    if (present(n3) .and. present(periodic)) then
      l = LexicalIndex(i, j, k, 0, 0, 0, n1+1, n2+1, n3+1, periodic)
    else
      l = LexicalIndex(i, j, k, 0, 0, 0, n1+1, n2+1)
    end if

  end function LexicalVertexIndex

  !-----------------------------------------------------------------------------
  !> Lexical index of an edge.

  pure integer function LexicalEdgeIndex(i, j, k, d, n1, n2, n3, periodic) &
      result(l)

    integer, intent(in) :: i   !< vertex index in direction 1
    integer, intent(in) :: j   !< vertex index in direction 2
    integer, intent(in) :: k   !< vertex index in direction 3
    integer, intent(in) :: d   !< direction (1, 2, or 3)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3
    logical, optional, intent(in) :: periodic(3) !< T periodic directions [F]

    logical :: p(3)
    integer :: m(3)

    if (present(periodic)) then
      p = periodic
    else
      p = .false.
    end if

    where (p)
      m = 0
    elsewhere
      m = 1
    end where

    select case(d)
    case(1)
      l = LexicalIndex(i,j,k, 1,0,0, n1, n2+1, n3+1, [ .false., p(2), p(3) ])
    case(2)
      l = LexicalIndex(i,j,k, 0,1,0, n1+1, n2, n3+1, [ p(1), .false., p(3) ]) &
        + n1 * (n2 + m(2)) * (n3 + m(3))
    case(3)
      l = LexicalIndex(i,j,k, 0,0,1, n1+1, n2+1, n3, [ p(1), p(2), .false. ]) &
        + n1 * (n2 + m(2)) * (n3 + m(3))  &
        + n2 * (n3 + m(3)) * (n1 + m(1))
    case default
      l = -1
    end select

  end function LexicalEdgeIndex

  !-----------------------------------------------------------------------------
  !> Lexical index of a face.
  !>
  !> The face are grouped according to their normal direction. Lexical numbering
  !> is applied within each group. Numbering starts at 1 with faces normal to
  !> direction 1, continues with direction 2 and, finally, direction 3.
  !> Thus, a unique ID is assigned to every mesh face, while retaining the
  !> natural ordering and easy accessibility of individual directions.

  pure integer function LexicalFaceIndex(i, j, k, d, n1, n2, n3, periodic) &
      result(l)

    integer, intent(in) :: i   !< vertex index in direction 1
    integer, intent(in) :: j   !< vertex index in direction 2
    integer, intent(in) :: k   !< vertex index in direction 3
    integer, intent(in) :: d   !< normal orientation (1,2, or 3)
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2
    integer, intent(in) :: n3  !< number of intervals in direction 3
    logical, optional, intent(in) :: periodic(3) !< T periodic directions [F]

    logical :: p(3)
    integer :: m(2)

    if (present(periodic)) then
      p = periodic
    else
      p = .false.
    end if

    where (p(1:2))
      m = 0
    elsewhere
      m = 1
    end where

    select case(d)
    case(1)
      l = LexicalIndex(i,j,k, 0,1,1, n1+1, n2, n3, [ p(1), .false., .false. ])
    case(2)
      l = LexicalIndex(i,j,k, 1,0,1, n1, n2+1, n3, [ .false., p(2), .false. ]) &
        + (n1 + m(1)) * n2 * n3
    case(3)
      l = LexicalIndex(i,j,k, 1,1,0, n1, n2, n3+1, [ .false., .false., p(3) ]) &
        + (n1 + m(1)) * n2 * n3  &
        + (n2 + m(2)) * n3 * n1
    case default
      l = -1
    end select

  end function LexicalFaceIndex

  !-----------------------------------------------------------------------------
  !> Lexical index of an element -- Integer Default Kind

  pure integer function LexicalElementIndex_IDK(i, j, k, n1, n2) result(l)

    integer, intent(in) :: i   !< element index in direction 1
    integer, intent(in) :: j   !< element index in direction 2
    integer, intent(in) :: k   !< element index in direction 3
    integer, intent(in) :: n1  !< number of intervals in direction 1
    integer, intent(in) :: n2  !< number of intervals in direction 2

    l = i + n1*(j-1) + n1*n2*(k-1)

  end function LexicalElementIndex_IDK

  !-----------------------------------------------------------------------------
  !> Lexical index of an element -- Integer eXtra Large

  pure integer(IXL) function LexicalElementIndex_IXL(i, j, k, n1, n2) result(l)

    integer(IXL), intent(in) :: i   !< element index in direction 1
    integer(IXL), intent(in) :: j   !< element index in direction 2
    integer(IXL), intent(in) :: k   !< element index in direction 3
    integer(IXL), intent(in) :: n1  !< number of intervals in direction 1
    integer(IXL), intent(in) :: n2  !< number of intervals in direction 2

    l = i + n1*(j-1) + n1*n2*(k-1)

  end function LexicalElementIndex_IXL

  !-----------------------------------------------------------------------------
  !> Recovers the triple index of a from the given lexical index.
  !>
  !> Converts the lexical index l into the corresponding triple index i, j, k
  !> for given start indices i0, j0, k0 and item numbers n1, n2 in the first two
  !> directions.
  !>
  !> This routine cannot handle periodicity!

  pure subroutine TripleIndex(i, j, k, i0, j0, k0, n1, n2, l)
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

  !=============================================================================

end module Mesh_3d__Structured_Indexing

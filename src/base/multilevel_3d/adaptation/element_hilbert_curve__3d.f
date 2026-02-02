!> summary:  Element Hilbert curve utilities
!> author:   Joerg Stiller
!> date:     2026/01/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> On each level, the elements are ordered according to their position on a
!> Hilbert curve that is defined in a bounding box enclosing the computational
!> domain. Inside a given element, the curve runs through the midpoints of all
!> children which (would) result from a regular subdivision of the element. The
!> path of this element Hilbert curve is defined by its start and end points or,
!> equivalently, the closest vertices. Generally, the path of the curve running
!> from vertex `m` to vertex `n`is encoded a single integer as `10*m + n`.
!> Every path is associated with a specific sequence of the child elements.
!> The subdivision also yields the next finer level of the Hilbert curve, which
!> is defined by the paths through the child elements.
!===============================================================================

module Element_Hilbert_Curve__3D
  implicit none
  private

  public :: EHC_ChildSequence
  public :: EHC_ChildCurvePaths

  !-----------------------------------------------------------------------------
  ! Child sequence and paths reference configurations

  !> child sequence along element Hilbert curve running from vertex 1 to 2
  integer, parameter :: CHILD_SEQ_12(8) = [ 1, 3, 7, 5, 6, 8, 4, 2 ]

  !> child curve paths for element Hilbert curve running from vertex 1 to 2
  integer, parameter :: CHILD_PATH_12(8) = [ 13, 42, 15, 62, 78, 73, 15, 12 ]


  !> triple child index `[ic(l), jc(l), kc(l)]` for given linear index
  integer, parameter :: ic(8) = [ 1, 2, 1, 2, 1, 2, 1, 2 ], &
                        jc(8) = [ 1, 1, 2, 2, 1, 1, 2, 2 ], &
                        kc(8) = [ 1, 1, 1, 1, 2, 2, 2, 2 ]

contains

  !-----------------------------------------------------------------------------
  !> Returns the sequence of children along given Hilbert curve path

  pure function EHC_ChildSequence(path) result(seq)
    integer, intent(in) :: path !< pass 10*m+n for start/end at vertex m/n
    integer :: seq(8)

    integer :: c, i, j, k, l

    do l = 1, 8
      c = CHILD_SEQ_12(l)
      i = ic(c)
      j = jc(c)
      k = kc(c)
      call TransformTripleIndex(path, i, j, k)
      seq(l) = i + 2*(j - 1) + 4*(k - 1)
    end do

  end function EHC_ChildSequence

  !-----------------------------------------------------------------------------
  !> Returns the Hilbert curve paths in children for given parent path

  pure function EHC_ChildCurvePaths(path_p) result(path_c)
    integer, intent(in) :: path_p !< pass 10*m+n for start/end at vertex m/n
    integer :: path_c(8)

    integer :: i, j, k, l, m, n, seq(8)

    seq = EHC_ChildSequence(path_p)

    ! processing child elements along parent SFC
    do l = 1, 8

      ! transform start point from path 12
      m = CHILD_PATH_12( CHILD_SEQ_12(l) ) / 10
      i = ic(m)
      j = jc(m)
      k = kc(m)
      call TransformTripleIndex(path_p, i, j, k)
      m = i + 2*(j - 1) + 4*(k - 1)

      ! transform end point from path 12
      n = mod(CHILD_PATH_12( CHILD_SEQ_12(l) ), 10)
      i = ic(n)
      j = jc(n)
      k = kc(n)
      call TransformTripleIndex(path_p, i, j, k)
      n = i + 2*(j - 1) + 4*(k - 1)

      ! insert transformed child path
      path_c(seq(l)) = 10*m + n

    end do

  end function EHC_ChildCurvePaths

  !-----------------------------------------------------------------------------
  !> Transformation of child triple index to given path
  !>
  !> Transforms a given triple index 1 ≤ i,j,k ≤ 2 to a coordinate system in
  !> which the start and end points of the given path correspond to vertex 1
  !> and 2, respectively.

  pure subroutine TransformTripleIndex(path, i, j, k)
    integer, intent(in)    :: path !< pass 10*m+n for start/end at vertex m/n
    integer, intent(inout) :: i    !< index 1 in original/transformed CS
    integer, intent(inout) :: j    !< index 2 in original/transformed CS
    integer, intent(inout) :: k    !< index 3 in original/transformed CS

    integer :: i0, j0, k0

    i0 = i
    j0 = j
    k0 = k

    select case(path)
    case(12);  i =     i0;  j =     j0;  k =     k0
    case(13);  i =     k0;  j =     i0;  k =     j0
    case(15);  i =     j0;  j =     k0;  k =     i0

    case(21);  i = 3 - i0;  j =     k0;  k =     j0
    case(24);  i = 3 - j0;  j =     i0;  k =     k0
    case(26);  i = 3 - k0;  j =     j0;  k =     i0

    case(34);  i =     i0;  j = 3 - k0;  k =     j0
    case(31);  i =     j0;  j = 3 - i0;  k =     k0
    case(37);  i =     k0;  j = 3 - j0;  k =     i0

    case(43);  i = 3 - i0;  j = 3 - k0;  k =     j0
    case(42);  i = 3 - k0;  j = 3 - i0;  k =     j0
    case(48);  i = 3 - j0;  j = 3 - k0;  k =     i0

    case(56);  i =     i0;  j =     k0;  k = 3 - j0
    case(57);  i =     j0;  j =     i0;  k = 3 - k0
    case(51);  i =     k0;  j =     j0;  k = 3 - i0

    case(65);  i = 3 - i0;  j =     j0;  k = 3 - k0
    case(68);  i = 3 - k0;  j =     i0;  k = 3 - j0
    case(62);  i = 3 - j0;  j =     k0;  k = 3 - i0

    case(78);  i =     i0;  j = 3 - j0;  k = 3 - k0
    case(75);  i =     k0;  j = 3 - i0;  k = 3 - j0
    case(73);  i =     j0;  j = 3 - k0;  k = 3 - i0

    case(87);  i = 3 - i0;  j = 3 - k0;  k = 3 - j0
    case(86);  i = 3 - j0;  j = 3 - i0;  k = 3 - k0
    case(84);  i = 3 - k0;  j = 3 - j0;  k = 3 - i0

    case default
      error stop 'TransformTripleIndex: invalid path'
    end select

  end subroutine TransformTripleIndex

  !=============================================================================

end module Element_Hilbert_Curve__3D

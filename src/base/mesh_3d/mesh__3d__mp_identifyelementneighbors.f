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

!> summary:  Identification of neighbor elements
!> author:   Joerg Stiller
!> date:     2021/09/14
!>
!> @note
!>   - works with unpartitioned mesh only
!>   - edge neighbor sorting seems to be obsolete
!>
!===============================================================================

submodule(Mesh__3D:MP_ImportGenericMesh) MP_IdentifyElementNeighbors
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Identification of neighbor elements

  module subroutine IdentifyElementNeighbors(mesh)
    class(Mesh_3D), intent(inout) :: mesh !< mesh partition

    integer, allocatable :: ne_face(:)  ! num elements per face
    integer, allocatable :: ne_edge(:)  ! num elements per edge
    integer, allocatable :: ne_vert(:)  ! num elements per vertex

    integer, allocatable :: nl_face(:)  ! num local elements per face
    integer, allocatable :: nl_edge(:)  ! num local elements per edge
    integer, allocatable :: nl_vert(:)  ! num local elements per vertex

    type NeighborElement
      integer      :: id        = 0 !< ID
      integer(IXS) :: component = 0 !< coupled neighbor component
    end type NeighborElement

    ! adjoining elements
    type(NeighborElement), allocatable :: elem_face(:,:) ! per face
    type(NeighborElement), allocatable :: elem_edge(:,:) ! per edge
    type(NeighborElement), allocatable :: elem_vert(:,:) ! per vertex

    type(MeshElementNeighbor_3D), allocatable :: neighbor(:)

    ! opposite element faces, edges and vertices
    integer, parameter :: opp_face( 6) = [ 2, 1, 4, 3, 6, 5 ]
    integer, parameter :: opp_edge(12) = [ 4, 3, 2, 1, 8, 7, 6, 5, 12, 11, 10, 9 ]
    integer, parameter :: opp_vert( 8) = [ 8, 7, 6, 5, 4, 3, 2, 1 ]

    logical, allocatable :: mask(:)
    logical :: right
    integer :: f(6), e(12), v(8)
    integer :: i, i1, i2, j, k, l, m, n, ne, nf, nn, self, step

    !---------------------------------------------------------------------------
    ! count number of adjoining elements for each vertex, edge and face

    ! local elements ...........................................................

    allocate(nl_face(mesh % n_face), source = 0)
    allocate(nl_edge(mesh % n_edge), source = 0)
    allocate(nl_vert(mesh % n_vert), source = 0)

    do l = 1, mesh % n_elem

      f = mesh % element(l) % face   % id
      e = mesh % element(l) % edge   % id
      v = mesh % element(l) % vertex % id

      do i = 1, 6
        nl_face(f(i)) = nl_face(f(i)) + 1
      end do

      do i = 1, 12
        nl_edge(e(i)) = nl_edge(e(i)) + 1
      end do

      do i = 1, 8
        nl_vert(v(i)) = nl_vert(v(i)) + 1
      end do

    end do

    !---------------------------------------------------------------------------
    ! establish lists of adjoining elements

    allocate(elem_face(maxval(nl_face), mesh%n_face))
    allocate(elem_edge(maxval(nl_edge), mesh%n_edge))
    allocate(elem_vert(maxval(nl_vert), mesh%n_vert))

    ! workspace based on upper bound for number of neighbors
    nn = 8 * maxval(nl_vert)
    allocate(neighbor(nn), mask(nn))

    ! counters
    allocate(ne_face(mesh % n_face), source = 0)
    allocate(ne_edge(mesh % n_edge), source = 0)
    allocate(ne_vert(mesh % n_vert), source = 0)

    ! adjoining element and components IDs .....................................

    do l = 1, mesh % n_elem

      f = mesh % element(l) % face   % id
      e = mesh % element(l) % edge   % id
      v = mesh % element(l) % vertex % id

      do j = 1, 6
        ne_face(f(j)) = ne_face(f(j)) + 1
        elem_face(ne_face(f(j)), f(j)) = NeighborElement(l, int(j, IXS))
      end do

      do j = 1, 12
        ne_edge(e(j)) = ne_edge(e(j)) + 1
        elem_edge(ne_edge(e(j)), e(j)) = NeighborElement(l, int(j+6, IXS))
      end do

      do j = 1, 8
        ne_vert(v(j)) = ne_vert(v(j)) + 1
        elem_vert(ne_vert(v(j)), v(j)) = NeighborElement(l, int(j+18, IXS))
      end do

    end do

    !---------------------------------------------------------------------------
    ! identify element neighbors

    do l = 1, mesh % n_elem
      associate(element => mesh % element(l))

        f = element % face   % id
        e = element % edge   % id
        v = element % vertex % id

        i = 0  ! neighbor counter

        ! faces ................................................................

        nf = 0  ! number of face neighbors

        do j = 1, 6
          associate(face => element % face(j))
            face % i_neighbor = i + 1
            face % n_neighbor = 0
            self = 0
            m = f(j)
            do k = 1, nl_face(m)

              ! accept self-reference if opposite face is matching
              if (elem_face(k,m) % id == l) then
                self = self + 1
                if (self /= 2) cycle ! skip first occurence
                if (element % face(opp_face(j)) % id /= face % id) cycle
              end if

              i = i + 1
              neighbor(i) % id        = elem_face(k,m) % id
              neighbor(i) % component = elem_face(k,m) % component
              neighbor(i) % part      = 0

              face % n_neighbor = face % n_neighbor + 1

              nf = nf + 1
            end do
            if (face % n_neighbor == 0) face % i_neighbor = -1
          end associate
        end do

        ! edges ................................................................

        ne = 0  ! number of edge neighbors

        do j = 1, 12
          associate(edge => element % edge(j))
            edge % i_neighbor = i + 1
            edge % n_neighbor = 0
            self = 0
            m = e(j)

            EDGE_NEIGHBORS: do k = 1, nl_edge(m)

              ! accept self-reference if opposite edge is matching
              if (elem_edge(k,m) % id == l) then
                self = self + 1
                if (self /= 2) cycle ! skip first occurence
                if (element % edge(opp_edge(j)) % id /= edge % id) cycle
              end if

              ! skip face neighbors
              do n = 1, nf
                if (elem_edge(k,m) % id == neighbor(n) % id) then
                  cycle EDGE_NEIGHBORS
                end if
              end do

              i = i + 1
              neighbor(i) % id        = elem_edge(k,m) % id
              neighbor(i) % component = elem_edge(k,m) % component
              neighbor(i) % part      = 0

              edge % n_neighbor = edge % n_neighbor + 1

              ne = ne + 1
            end do EDGE_NEIGHBORS
            if (edge % n_neighbor == 0) edge % i_neighbor = -1
          end associate
        end do

        ! vertices
        do j = 1, 8
          associate(vertex => element % vertex(j))
            vertex % i_neighbor = i + 1
            vertex % n_neighbor = 0
            self = 0
            m = v(j)

            VERTEX_NEIGHBORS: do k = 1, nl_vert(m)

              ! accept self-reference if opposite vertex is matching
              if (elem_vert(k,m) % id == l) then
                self = self + 1
                if (self /= 2) cycle ! skip first occurence
                if (element % vertex(opp_vert(j)) % id /= vertex % id) cycle
              end if

              ! skip face and edge neighbors
              do n = 1, nf+ne
                if (elem_vert(k,m) % id == neighbor(n) % id) then
                  cycle VERTEX_NEIGHBORS
                end if
              end do

              i = i + 1
              neighbor(i) % id        = elem_vert(k,m) % id
              neighbor(i) % component = elem_vert(k,m) % component
              neighbor(i) % part      = 0

              vertex % n_neighbor = vertex % n_neighbor + 1

            end do VERTEX_NEIGHBORS
            if (vertex % n_neighbor == 0) vertex % i_neighbor = -1
          end associate
        end do

        element % neighbor = neighbor(1:i)

      end associate
    end do

    !---------------------------------------------------------------------------
    ! order element edge neighbors according to right-hand rule

    do l = 1, mesh % n_elem
      do k = 1, 12
        if (mesh % element(l) % edge(k) % n_neighbor < 2) cycle
        call InitEdgeNeighborSorting(mesh % element(l), k, n, right)
        if (n == 0) cycle ! skip degenerate case
        nn = mesh % element(l) % edge(k) % n_neighbor
        i1 = mesh % element(l) % edge(k) % i_neighbor
        i2 = i1 + nn - 1
        neighbor(1:nn) = mesh % element(l) % neighbor(i1:i2)
        mask(1:nn) = .true.
        if (right) then
          step =  1
        else ! reverse search
          step = -1
          i2   = i1
          i1   = i2 + nn - 1
        end if
        do i = i1, i2, step
          SEARCH: do j = 1, nn
            if (.not. mask(j)) cycle
            m = neighbor(j) % id
            if (IsFaceNeighbor(mesh % element(m), n)) then
              mesh % element(l) % neighbor(i) = neighbor(j)
              mask(j) = .false.
              n = m
              exit SEARCH
            end if
          end do SEARCH
        end do
        if (any(mask(1:nn))) then
          call Error('IdentifyElementNeighbors','Failed sorting egde neighbors')
        end if
      end do
    end do

    !---------------------------------------------------------------------------
    ! identification of neighbor element orientations

    associate(element => mesh % element)

      do l = 1, mesh % n_elem

        ! orientation of face neighbors ..........................................

        do j = 1, 6
          if (element(l) % face(j) % n_neighbor /= 1) cycle
          i = element(l) % face(j) % i_neighbor
          n = element(l) % neighbor(i) % id
          k = ElementFaceID(element(l) % neighbor(i) % component)

          element(l) % neighbor(i) % orientation =  &
              FaceNeighborOrientation(element(l), j, element(n), k)

        end do

        ! orientation of edge neighbors ..........................................

        do j = 1, 12
          if (element(l) % edge(j) % n_neighbor < 1) cycle
          i1 = element(l) % edge(j) % i_neighbor
          i2 = element(l) % edge(j) % n_neighbor + i1 - 1
          do i = i1, i2
            n = element(l) % neighbor(i) % id
            k = ElementEdgeID(element(l) % neighbor(i) % component)

            element(l) % neighbor(i) % orientation =  &
                EdgeNeighborOrientation(element(l), j, element(n), k)

          end do
        end do

        ! orientation of regular vertex neighbors ................................

        do j = 1, 8
          if (element(l) % vertex(j) % n_neighbor /= 1) cycle

          i = element(l) % vertex(j) % i_neighbor
          n = element(l) % neighbor(i) % id
          k = ElementVertexID(element(l) % neighbor(i) % component)

          element(l) % neighbor(i) % orientation =  &
              VertexNeighborOrientation(element(l), j, element(n), k)

        end do

      end do

    end associate

  end subroutine IdentifyElementNeighbors

  !=============================================================================
  ! Procedures for edge sorting

  !-----------------------------------------------------------------------------
  !> Provides the initial data for edge neighbor sorting
  !>
  !> `n > 0` is the ID of neighbor element adjacent to the starting face and
  !> `right` the search direction (T: positive, F: negative orientation).
  !> `n = 0` indicates the degenerate case with no neighbors at the adjoining
  !> faces.

  subroutine InitEdgeNeighborSorting(element, k, n, right)
    class(MeshElement_3D), intent(in)  :: element !< mesh element
    integer,               intent(in)  :: k       !< element edge ID
    integer,               intent(out) :: n       !< initial neighbor, 0 if none
    logical,               intent(out) :: right   !< search direction

    integer :: i

    associate(face => element % face)
      select case(k)
      case(1)
        right = face(3) % n_neighbor > 0
        if (right) then
          i = face(3) % i_neighbor
        else
          i = face(5) % i_neighbor
        end if
      case(2)
        right = face(5) % n_neighbor > 0
        if (right) then
          i = face(5) % i_neighbor
        else
          i = face(4) % i_neighbor
        end if
      case(3)
        right = face(6) % n_neighbor > 0
        if (right) then
          i = face(6) % i_neighbor
        else
          i = face(3) % i_neighbor
        end if
      case(4)
        right = face(4) % n_neighbor > 0
        if (right) then
          i = face(4) % i_neighbor
        else
          i = face(6) % i_neighbor
        end if
      case(5)
        right = face(5) % n_neighbor > 0
        if (right) then
          i = face(5) % i_neighbor
        else
          i = face(1) % i_neighbor
        end if
      case(6)
        right = face(2) % n_neighbor > 0
        if (right) then
          i = face(2) % i_neighbor
        else
          i = face(5) % i_neighbor
        end if
      case(7)
        right = face(1) % n_neighbor > 0
        if (right) then
          i = face(1) % i_neighbor
        else
          i = face(6) % i_neighbor
        end if
      case(8)
        right = face(6) % n_neighbor > 0
        if (right) then
          i = face(6) % i_neighbor
        else
          i = face(2) % i_neighbor
        end if
      case(9)
        right = face(1) % n_neighbor > 0
        if (right) then
          i = face(1) % i_neighbor
        else
          i = face(3) % i_neighbor
        end if
      case(10)
        right = face(3) % n_neighbor > 0
        if (right) then
          i = face(3) % i_neighbor
        else
          i = face(2) % i_neighbor
        end if
      case(11)
        right = face(4) % n_neighbor > 0
        if (right) then
          i = face(4) % i_neighbor
        else
          i = face(1) % i_neighbor
        end if
      case(12)
        right = face(2) % n_neighbor > 0
        if (right) then
          i = face(2) % i_neighbor
        else
          i = face(4) % i_neighbor
        end if
      end select
    end associate

    if (i > 0) then
      n = element % neighbor(i) % id
    else
      n = 0
    end if

  end subroutine InitEdgeNeighborSorting

  !-----------------------------------------------------------------------------
  !> Checks whether element `e` is flushed with the given element

  pure logical function IsFaceNeighbor(element, e) result(r)
    class(MeshElement_3D), intent(in)  :: element !< mesh element
    integer,               intent(in)  :: e       !< element ID

    integer :: i

    r = .false.
    do i = 1, 6
      if (element % face(i) % n_neighbor < 1) cycle
      r = element % neighbor(element % face(i) % i_neighbor) % id == e
      if (r) exit
    end do

  end function IsFaceNeighbor

  !=============================================================================
  ! Common procedures for identifying neighbor orientation

  pure integer(IXS) function Orientation(de, dn) result(o)

    integer, intent(in) :: de(3,3)  ! element directions,  `de(:,i) = di`
    integer, intent(in) :: dn(3,3)  ! neighbor directions, `dn(:,i) = di`

    integer :: i, j, k

    ! identify element direction parallel to ξ
    do k = 1, 3
      if (abs(de(1,k)) == 1) exit
    end do

    ! identify matching neighbor direction
    do i = 1, 3
      if (abs(dn(i,k)) == 1) exit
    end do

    if (de(1,k) * dn(i,k) > 0) then ! parallel
      o = 10 * i
    else ! antiparallel
      o = 10 * i + 30
    end if

    ! identify element direction parallel to η
    do k = 1, 3
      if (abs(de(2,k)) == 1) exit
    end do

    ! identify matching neighbor direction
    do j = 1, 3
      if (abs(dn(j,k)) == 1) exit
    end do

    if (de(2,k) * dn(j,k) > 0) then ! parallel
      o = o + j
    else ! antiparallel
      o = o + j + 3
    end if

  end function Orientation

  !=============================================================================
  ! Procedures for identifying face neighbor orientation

  integer(IXS) function FaceNeighborOrientation(element, f, neighbor, fc)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer,               intent(in) :: f        !< element face
    class(MeshElement_3D), intent(in) :: neighbor !< neighbor at face f
    integer,               intent(in) :: fc       !< coupled neighbor face

    integer :: de(3,3)  ! element face directions
    integer :: dn(3,3)  ! neighbor face directions

    call FaceDirections( f                              &
                       , element % face(f) % normal     &
                       , element % face(f) % rotation   &
                       , de(:,1), de(:,2), de(:,3)      )

    call FaceDirections( fc                             &
                       , neighbor % face(fc) % normal   &
                       , neighbor % face(fc) % rotation &
                       , dn(:,1), dn(:,2), dn(:,3)      )

    FaceNeighborOrientation = Orientation(de, dn)

  end function FaceNeighborOrientation

  !-----------------------------------------------------------------------------
  !> Returns the directions of an element face matching the adjoining mesh face

  pure subroutine FaceDirections(f, normal, rotation, d1, d2, d3)
    integer,      intent(in)  :: f        !< element face
    integer(IXS), intent(in)  :: normal   !< normal orientation against MF (±1)
    integer(IXS), intent(in)  :: rotation !< 1/4-rotation for aligning with MF
    integer,      intent(out) :: d1(3)    !< tangential direction 1
    integer,      intent(out) :: d2(3)    !< tangential direction 2
    integer,      intent(out) :: d3(3)    !< normal direction

    integer, dimension(3), parameter :: xi   = [1,0,0]
    integer, dimension(3), parameter :: eta  = [0,1,0]
    integer, dimension(3), parameter :: zeta = [0,0,1]

    if (normal == 1) then
      select case(f)
      case(1,2)
        select case(int( rotation ))
        case(0)
          d1 =  eta
          d2 =  zeta
          d3 =  xi
        case(1)
          d1 =  zeta
          d2 = -eta
          d3 =  xi
        case(2)
          d1 = -eta
          d2 = -zeta
          d3 =  xi
        case default ! 3
          d1 = -zeta
          d2 =  eta
          d3 =  xi
        end select
      case(3,4)
        select case(int( rotation ))
        case(0)
          d1 =  xi
          d2 =  zeta
          d3 = -eta
        case(1)
          d1 =  zeta
          d2 = -xi
          d3 = -eta
        case(2)
          d1 = -xi
          d2 = -zeta
          d3 = -eta
        case default ! 3
           d1 = -zeta
           d2 =  xi
           d3 = -eta
        end select
      case default ! 5,6
        select case(int( rotation ))
        case(0)
          d1 =  xi
          d2 =  eta
          d3 =  zeta
        case(1)
          d1 =  eta
          d2 = -xi
          d3 =  zeta
        case(2)
          d1 = -xi
          d2 = -eta
          d3 =  zeta
        case default ! 3
          d1 = -eta
          d2 =  xi
          d3 =  zeta
        end select
     end select
    else ! normal = -1
      select case(f)
      case(1,2)
        select case(int( rotation ))
        case(0)
          d1 =  eta
          d2 = -zeta
          d3 = -xi
        case(1)
          d1 = -zeta
          d2 = -eta
          d3 = -xi
        case(2)
          d1 = -eta
          d2 =  zeta
          d3 = -xi
        case default ! 3
          d1 =  zeta
          d2 =  eta
          d3 = -xi
        end select
      case(3,4)
        select case(int( rotation ))
        case(0)
          d1 =  xi
          d2 = -zeta
          d3 =  eta
        case(1)
          d1 = -zeta
          d2 = -xi
          d3 =  eta
        case(2)
          d1 = -xi
          d2 =  zeta
          d3 =  eta
        case default ! 3
          d1 =  zeta
          d2 =  xi
          d3 =  eta
        end select
      case default ! 5,6
        select case(int( rotation ))
        case(0)
          d1 =  xi
          d2 = -eta
          d3 = -zeta
        case(1)
          d1 = -eta
          d2 = -xi
          d3 = -zeta
        case(2)
          d1 = -xi
          d2 =  eta
          d3 = -zeta
        case default ! 3
          d1 =  eta
          d2 =  xi
          d3 = -zeta
        end select
      end select
    end if

  end subroutine FaceDirections

  !=============================================================================
  ! Procedures for identifying edge neighbor orientation

  !-----------------------------------------------------------------------------
  !> Returns the orientation of an edge neighbor

  integer(IXS) function EdgeNeighborOrientation(element, e, neighbor, ec)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer,               intent(in) :: e        !< element edge
    class(MeshElement_3D), intent(in) :: neighbor !< neighbor at edge e
    integer,               intent(in) :: ec       !< coupled neighbor edge

    integer :: de(3,3)  ! element edge directions
    integer :: dn(3,3)  ! neighbor edge directions

    ! given element edge directions
    call ParallelEdgeDirections(e, de(:,1), de(:,2), de(:,3))

    ! neighbor element edge directions
    if (element%edge(e)%orientation == neighbor%edge(ec)%orientation) then
      call ParallelEdgeDirections(ec, dn(:,1), dn(:,2), dn(:,3))
    else
      call AntiparallelEdgeDirections(ec, dn(:,1), dn(:,2), dn(:,3))
    end if

    ! rotate neighbor trihedron to match given element
    dn(:,2) = -dn(:,2)
    dn(:,3) = -dn(:,3)

    EdgeNeighborOrientation = Orientation(de, dn)

  end function EdgeNeighborOrientation

  !-----------------------------------------------------------------------------
  !> Returns the directions of the right-handed trihedron formed by element edge
  !> `e` and the edges emerging from its starting point

  pure subroutine ParallelEdgeDirections(e, d1, d2, d3)
    integer, intent(in)  :: e     !< element edge
    integer, intent(out) :: d1(3) !< direction 1: aligned with edge
    integer, intent(out) :: d2(3) !< direction 2
    integer, intent(out) :: d3(3) !< direction 3

    integer :: v0(3)

    select case(e)
    case(1)
      v0 = [ -1, -1, -1 ]         ! vertex 1
      d1 = [  1, -1, -1 ] - v0    ! vertex 1 → 2
      d2 = [ -1,  1, -1 ] - v0    ! vertex 1 → 3
      d3 = [ -1, -1,  1 ] - v0    ! vertex 1 → 5
    case(2)
      v0 = [ -1,  1, -1 ]         ! vertex 3
      d1 = [  1,  1, -1 ] - v0    ! vertex 3 → 4
      d2 = [ -1,  1,  1 ] - v0    ! vertex 3 → 7
      d3 = [ -1, -1, -1 ] - v0    ! vertex 3 → 1
    case(3)
      v0 = [ -1, -1,  1 ]         ! vertex 5
      d1 = [  1, -1,  1 ] - v0    ! vertex 5 → 6
      d2 = [ -1, -1, -1 ] - v0    ! vertex 5 → 1
      d3 = [ -1,  1,  1 ] - v0    ! vertex 5 → 7
    case(4)
      v0 = [ -1,  1,  1 ]         ! vertex 7
      d1 = [  1,  1,  1 ] - v0    ! vertex 7 → 8
      d2 = [ -1, -1,  1 ] - v0    ! vertex 7 → 5
      d3 = [ -1,  1, -1 ] - v0    ! vertex 7 → 3
    case(5)
      v0 = [ -1, -1, -1 ]         ! vertex 1
      d1 = [ -1,  1, -1 ] - v0    ! vertex 1 → 3
      d2 = [ -1, -1,  1 ] - v0    ! vertex 1 → 5
      d3 = [  1, -1, -1 ] - v0    ! vertex 1 → 2
    case(6)
      v0 = [  1, -1, -1 ]         ! vertex 2
      d1 = [  1,  1, -1 ] - v0    ! vertex 2 → 4
      d2 = [ -1, -1, -1 ] - v0    ! vertex 2 → 1
      d3 = [  1, -1,  1 ] - v0    ! vertex 2 → 6
    case(7)
      v0 = [ -1, -1,  1 ]         ! vertex 5
      d1 = [ -1,  1,  1 ] - v0    ! vertex 5 → 7
      d2 = [  1, -1,  1 ] - v0    ! vertex 5 → 6
      d3 = [ -1, -1, -1 ] - v0    ! vertex 5 → 1
    case(8)
      v0 = [  1, -1,  1 ]         ! vertex 6
      d1 = [  1,  1,  1 ] - v0    ! vertex 6 → 8
      d2 = [  1, -1, -1 ] - v0    ! vertex 6 → 2
      d3 = [ -1, -1,  1 ] - v0    ! vertex 6 → 5
    case(9)
      v0 = [ -1, -1, -1 ]         ! vertex 1
      d1 = [ -1, -1,  1 ] - v0    ! vertex 1 → 5
      d2 = [  1, -1, -1 ] - v0    ! vertex 1 → 2
      d3 = [ -1,  1, -1 ] - v0    ! vertex 1 → 3
    case(10)
      v0 = [  1, -1, -1 ]         ! vertex 2
      d1 = [  1, -1,  1 ] - v0    ! vertex 2 → 6
      d2 = [  1,  1, -1 ] - v0    ! vertex 2 → 4
      d3 = [ -1, -1, -1 ] - v0    ! vertex 2 → 1
    case(11)
      v0 = [ -1,  1, -1 ]         ! vertex 3
      d1 = [ -1,  1,  1 ] - v0    ! vertex 3 → 7
      d2 = [ -1, -1, -1 ] - v0    ! vertex 3 → 1
      d3 = [  1,  1, -1 ] - v0    ! vertex 3 → 4
    case default ! 12
      v0 = [  1,  1, -1 ]         ! vertex 4
      d1 = [  1,  1,  1 ] - v0    ! vertex 4 → 8
      d2 = [ -1,  1, -1 ] - v0    ! vertex 4 → 3
      d3 = [  1, -1, -1 ] - v0    ! vertex 4 → 2
    end select

    d1 = d1 / 2
    d2 = d2 / 2
    d3 = d3 / 2

  end subroutine ParallelEdgeDirections

  !-----------------------------------------------------------------------------
  !> Returns the directions of the right-handed trihedron formed by element edge
  !> `e` and the edges emerging from its end point

  pure subroutine AntiparallelEdgeDirections(e, d1, d2, d3)
    integer, intent(in)  :: e     !< element edge
    integer, intent(out) :: d1(3) !< direction 1: aligned with reversed edge
    integer, intent(out) :: d2(3) !< direction 2
    integer, intent(out) :: d3(3) !< direction 3

    integer :: v0(3)

    select case(e)
    case(1)
      v0 = [  1, -1, -1 ]         ! vertex 2
      d1 = [ -1, -1, -1 ] - v0    ! vertex 2 → 1
      d2 = [  1, -1,  1 ] - v0    ! vertex 2 → 6
      d3 = [  1,  1, -1 ] - v0    ! vertex 2 → 4
    case(2)
      v0 = [  1,  1, -1 ]         ! vertex 4
      d1 = [ -1,  1, -1 ] - v0    ! vertex 4 → 3
      d2 = [  1, -1, -1 ] - v0    ! vertex 4 → 2
      d3 = [  1,  1,  1 ] - v0    ! vertex 4 → 8
    case(3)
      v0 = [  1, -1,  1 ]         ! vertex 6
      d1 = [ -1, -1,  1 ] - v0    ! vertex 6 → 5
      d2 = [  1,  1,  1 ] - v0    ! vertex 6 → 8
      d3 = [  1, -1, -1 ] - v0    ! vertex 6 → 2
    case(4)
      v0 = [  1,  1,  1 ]         ! vertex 8
      d1 = [ -1,  1,  1 ] - v0    ! vertex 8 → 7
      d2 = [  1,  1, -1 ] - v0    ! vertex 8 → 4
      d3 = [  1, -1,  1 ] - v0    ! vertex 8 → 6
    case(5)
      v0 = [ -1,  1, -1 ]         ! vertex 3
      d1 = [ -1, -1, -1 ] - v0    ! vertex 3 → 1
      d2 = [  1,  1, -1 ] - v0    ! vertex 3 → 4
      d3 = [ -1,  1,  1 ] - v0    ! vertex 3 → 7
    case(6)
      v0 = [  1,  1, -1 ]         ! vertex 4
      d1 = [  1, -1, -1 ] - v0    ! vertex 4 → 2
      d2 = [  1,  1,  1 ] - v0    ! vertex 4 → 8
      d3 = [ -1,  1, -1 ] - v0    ! vertex 4 → 3
    case(7)
      v0 = [ -1,  1,  1 ]         ! vertex 7
      d1 = [ -1, -1,  1 ] - v0    ! vertex 7 → 5
      d2 = [ -1,  1, -1 ] - v0    ! vertex 7 → 3
      d3 = [  1,  1,  1 ] - v0    ! vertex 7 → 8
    case(8)
      v0 = [  1,  1,  1 ]         ! vertex 8
      d1 = [  1, -1,  1 ] - v0    ! vertex 8 → 6
      d2 = [ -1,  1,  1 ] - v0    ! vertex 8 → 7
      d3 = [  1,  1, -1 ] - v0    ! vertex 8 → 4
    case(9)
      v0 = [ -1, -1,  1 ]         ! vertex 5
      d1 = [ -1, -1, -1 ] - v0    ! vertex 5 → 1
      d2 = [ -1,  1,  1 ] - v0    ! vertex 5 → 7
      d3 = [  1, -1,  1 ] - v0    ! vertex 5 → 6
    case(10)
      v0 = [  1, -1,  1 ]         ! vertex 6
      d1 = [  1, -1, -1 ] - v0    ! vertex 6 → 2
      d2 = [ -1, -1,  1 ] - v0    ! vertex 6 → 5
      d3 = [  1,  1,  1 ] - v0    ! vertex 6 → 8
    case(11)
      v0 = [ -1,  1,  1 ]         ! vertex 7
      d1 = [ -1,  1, -1 ] - v0    ! vertex 7 → 3
      d2 = [  1,  1,  1 ] - v0    ! vertex 7 → 8
      d3 = [ -1, -1,  1 ] - v0    ! vertex 7 → 5
    case default ! 12
      v0 = [  1,  1,  1 ]         ! vertex 8
      d1 = [  1,  1, -1 ] - v0    ! vertex 8 → 4
      d2 = [  1, -1,  1 ] - v0    ! vertex 8 → 6
      d3 = [ -1,  1,  1 ] - v0    ! vertex 8 → 7
    end select

    d1 = d1 / 2
    d2 = d2 / 2
    d3 = d3 / 2

  end subroutine AntiparallelEdgeDirections

  !=============================================================================
  ! Procedures for identifying vertex neighbor orientation

  !-----------------------------------------------------------------------------
  !> Returns the orientation of a regular vertex neighbor
  !>
  !> The orientation is determined only in the regular case, where the given
  !> element has a neighbor at the three adjoining faces. It is encoded in a
  !> two-digit integer `ab`, where `a` and `b` indicate the directions of the
  !> neighbor coinciding with the ξ- and η-directions of the given element.
  !> Values 1,2,3 correspond to ξ,η,ζ and 4,5,6 to -ξ,-η,-ζ, respectively.
  !> For example, `ab = 53` means that ξ and η in element match -η and ζ in
  !> neighbor.
  !>
  !> A return value of 0 indicates an irregular configuration and -1 an error.

  integer(IXS) function VertexNeighborOrientation(element, v, neighbor, vc)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer,               intent(in) :: v        !< element vertex
    class(MeshElement_3D), intent(in) :: neighbor !< neighbor at vertex v
    integer,               intent(in) :: vc       !< coupled neighbor vertex

    integer :: d1, d2, n, s

    ! preliminaries ............................................................

    ! check for regularity
    if ( element % vertex(v) % n_neighbor /= 1 .or.           &
         any(element % face(F_VERT(:,v)) % n_neighbor /= 1) ) &
    then
      VertexNeighborOrientation = -1
      return
    end if

    ! neighbor direction matching ±ξ ...........................................

    ! face neighbor adjacent to the prolongated ξ-edge through the vertex
    select case(v)
    case(1,3,5,7) ! face 1
      s = -1      ! prolongated edge points in -ξ direction
      n = element % neighbor(element % face(1) % i_neighbor) % id
    case default  ! face 2
      s =  1      ! prolongated edge points in +ξ direction
      n = element % neighbor(element % face(2) % i_neighbor) % id
    end select

    ! matching neighbor direction
    select case(AdjoiningEdge(neighbor, n))

    case(1:4)         ! ξ-edge
      select case(vc) ! position of coupled vertex
      case(1,3,5,7)   ! on face 1
        d1 =  1 * s
      case default    ! on face 2
        d1 = -1 * s
      end select

    case(5:8)         ! η-edge
      select case(vc) ! position of coupled vertex
      case(1,2,5,6)   ! on face 3
        d1 =  2 * s
      case default    ! on face 4
        d1 = -2 * s
      end select

    case(9:12)        ! ζ-edge
      select case(vc) ! position of coupled vertex
      case(1,2,3,4)   ! on face 5
        d1 =  3 * s
      case default    ! on face 6
        d1 = -3 * s
      end select

    case default      ! search failed
      d1 = 0

    end select

    ! neighbor direction matching ±η ...........................................

    ! face neighbor adjacent to the prolongated η-edge through the vertex
    select case(v)
    case(1,2,5,6) ! face 3
      s = -1      ! edge prolongated in -η direction
      n = element % neighbor(element % face(3) % i_neighbor) % id
    case default  ! face 4
      s =  1      ! edge prolongated in +η direction
      n = element % neighbor(element % face(4) % i_neighbor) % id
    end select

    ! matching neighbor direction
    select case(AdjoiningEdge(neighbor, n))

    case(1:4)         ! ξ-edge
      select case(vc) ! position of coupled vertex
      case(1,3,5,7)   ! on face 1
        d2 =  1 * s
      case default    ! on face 2
        d2 = -1 * s
      end select

    case(5:8)         ! η-edge
      select case(vc) ! position of coupled vertex
      case(1,2,5,6)   ! on face 3
        d2 =  2 * s
      case default    ! on face 4
        d2 = -2 * s
      end select

    case(9:12)        ! ζ-edge
      select case(vc) ! position of coupled vertex
      case(1,2,3,4)   ! on face 5
        d2 =  3 * s
      case default    ! on face 6
        d2 = -3 * s
      end select

    case default      ! search failed
      d2 = 0

    end select

    ! orientation ..............................................................

    if (d1 == 0 .or. d2 == 0) then
      VertexNeighborOrientation = -1
    else
      if (d1 < 0) d1 = 3 - d1
      if (d2 < 0) d2 = 3 - d2
      VertexNeighborOrientation = 10 * d1 + d2
    end if

  end function VertexNeighborOrientation

  !-----------------------------------------------------------------------------
  !> Edge of the vertex neighbor adjoining the face neighbor element

  pure integer function AdjoiningEdge(element, n) result(e)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer,               intent(in) :: n        !< vertex neighbor element ID

    integer :: i, i1, i2, j

    EDGES: do j = 1, 12
      if (element % edge(j) % n_neighbor < 1) then
        cycle
      else
        i1 = element % edge(j) % i_neighbor
        i2 = element % edge(j) % n_neighbor + i1 - 1
        do i = i1, i2
          if (element % neighbor(i) % id == n) exit EDGES
        end do
      end if
    end do EDGES

    if (j <= 12) then
      e = j
    else
      e = -1
    end if

  end function AdjoiningEdge

  !=============================================================================

end submodule MP_IdentifyElementNeighbors

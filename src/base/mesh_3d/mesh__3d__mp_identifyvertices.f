!> summary:  Identification of mesh vertices
!> author:   Joerg Stiller
!> date:     2022/07/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_IdentifyVertices
  use Mesh_Element_Indexing__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Identification of mesh edges
  !>
  !> Requires
  !>   - `mesh % element % neighbor` for local elements
  !>
  !> Generates
  !>   - mesh % n_vert
  !>   - mesh % element % vertex % id

  module subroutine IdentifyVertices(mesh)
    class(Mesh_3D), intent(inout) :: mesh !< mesh partition

    ! local data ...............................................................

    integer :: vid(8,mesh%n_elem), vid_ef(2,2), vid_nf(2,2)
    integer :: e, i, j, k, l, m, v
    logical :: periodic(3)

    associate(element => mesh % element)

      ! initialization .........................................................

      vid = 0

      v = 0
      do e = 1, mesh % n_elem

        ! adopt vertex IDs from face neighbors .................................

        periodic = .false.

        do k = 1, 6
          if (element(e) % face(k) % n_neighbor == 1) then
            i = element(e) % face(k) % i_neighbor
            if (element(e) % neighbor(i) % part /= mesh % part) cycle ! ghosts
            if (element(e) % neighbor(i) % id > e) cycle ! untouched neighbors

            l = element(e) % neighbor(i) % id
            m = ElementFaceID(element(e) % neighbor(i) % component)

            if (l == e) then
              if (k == 2) periodic(1) = .true.
              if (k == 4) periodic(2) = .true.
              if (k == 6) periodic(3) = .true.
              cycle
            end if

            vid_nf(1,1) = vid(V_FACE(1,m), l)
            vid_nf(2,1) = vid(V_FACE(2,m), l)
            vid_nf(1,2) = vid(V_FACE(3,m), l)
            vid_nf(2,2) = vid(V_FACE(4,m), l)

            call element(l) % AlignFromNeighborFace(k, i, vid_nf, vid_ef)

            vid(V_FACE(1,k), e) = vid_ef(1,1)
            vid(V_FACE(2,k), e) = vid_ef(2,1)
            vid(V_FACE(3,k), e) = vid_ef(1,2)
            vid(V_FACE(4,k), e) = vid_ef(2,2)
          end if
        end do

        ! adopt vertex IDs from edge neighbors .................................

        do k = 1, 12
          j = element(e) % edge(k) % i_neighbor
          do i = j, j + element(e) % edge(k) % n_neighbor - 1
            if (element(e) % neighbor(i) % part /= mesh % part) cycle
            if (element(e) % neighbor(i) % id >= e) cycle

            l = element(e) % neighbor(i) % id
            m = ElementEdgeID(element(e) % neighbor(i) % component)

            if (element(e) % NeigborEdgeIsAligned(k, i)) then
              vid(V_EDGE(1,k), e) = vid(V_EDGE(1,m), l)
              vid(V_EDGE(2,k), e) = vid(V_EDGE(2,m), l)
            else
              vid(V_EDGE(1,k), e) = vid(V_EDGE(2,m), l)
              vid(V_EDGE(2,k), e) = vid(V_EDGE(1,m), l)
            end if

            exit ! skip remaining neighbors on this edge
          end do
        end do

        ! adopt vertex IDs from vertex neighbors ...............................

        do k = 1, 8
          j = element(e) % vertex(k) % i_neighbor
          do i = j, j + element(e) % vertex(k) % n_neighbor - 1
            if (element(e) % neighbor(i) % part /= mesh % part) cycle ! ghosts
            if (element(e) % neighbor(i) % id >= e) cycle ! untouched neighbors

            l = element(e) % neighbor(i) % id
            m = ElementVertexID(element(e) % neighbor(i) % component)

            vid(k,e) = vid(m,l)

            exit ! skip remaining neighbors on this vertex
          end do
        end do

        ! identify remaining vertices ..........................................

        k = 1
        if (vid(k,e) < 1) then
          v = v + 1
          vid(k,e) = v
        end if

        k = 2
        if (periodic(1)) then
          vid(k,e) = vid(1,e)
        else if (vid(k,e) < 1) then
          v = v + 1
          vid(k,e) = v
        end if

        k = 3
        if (periodic(2)) then
          vid(k,e) = vid(1,e)
        else if (vid(k,e) < 1) then
          v = v + 1
          vid(k,e) = v
        end if

        k = 4
        if (periodic(1)) then
          vid(k,e) = vid(3,e)
        else if (periodic(2)) then
          vid(k,e) = vid(2,e)
        else if (vid(k,e) < 1) then
          v = v + 1
          vid(k,e) = v
        end if

        if (periodic(3)) then

          do k = 5, 8
            vid(k,e) = vid(k-4,e)
          end do

        else

          k = 5
          if (vid(k,e) < 1) then
            v = v + 1
            vid(k,e) = v
          end if

          k = 6
          if (periodic(1)) then
            vid(k,e) = vid(5,e)
          else if (vid(k,e) < 1) then
            v = v + 1
            vid(k,e) = v
          end if

          k = 7
          if (periodic(2)) then
            vid(k,e) = vid(5,e)
          else if (vid(k,e) < 1) then
            v = v + 1
            vid(k,e) = v
          end if

          k = 8
          if (periodic(1)) then
            vid(k,e) = vid(7,e)
          else if (periodic(2)) then
            vid(k,e) = vid(6,e)
          else if (vid(k,e) < 1) then
            v = v + 1
            vid(k,e) = v
          end if

        end if

        ! set element vertex IDs ...............................................

        do k = 1, 8
          element(e) % vertex(k) % id = vid(k,e)
        end do

      end do

      mesh % n_vert = v

    end associate

  end subroutine IdentifyVertices

  !=============================================================================

end submodule MP_IdentifyVertices

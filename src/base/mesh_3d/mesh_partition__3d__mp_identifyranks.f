!> summary:  Identification of element component ranks
!> author:   Joerg Stiller
!> date:     2021/02/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_Partition__3D) MP_IdentifyRanks
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Identification of element component ranks
  !>
  !> Requires
  !>   - mesh % element % face   % {n_neighbor, i_neighbor}
  !>   - mesh % element % edge   % {n_neighbor, i_neighbor}
  !>   - mesh % element % vertex % {n_neighbor, i_neighbor}
  !>   - mesh % element % neighbor
  !>
  !> Generates
  !>   - mesh % max_vert_val
  !>   - mesh % max_edge_val
  !>   - mesh % element % face   % {rank, val}
  !>   - mesh % element % edge   % {rank, val}
  !>   - mesh % element % vertex % {rank, val}

  module subroutine IdentifyRanks(mesh)
    class(MeshPartition_3D), intent(inout) :: mesh !< mesh partition

    integer, allocatable :: m_vert(:), m_edge(:), m_face(:)
    integer :: max_vert_val, max_edge_val
    integer :: e, i, k

    allocate(m_vert(mesh % n_vert), source = 0)
    allocate(m_edge(mesh % n_edge), source = 0)
    allocate(m_face(mesh % n_face), source = 0)
    max_vert_val = 0
    max_edge_val = 0

    ! ranks ....................................................................

    do e = 1, mesh % n_elem
      associate(element => mesh % element(e))

        do i = 1, 8
          k = element % vertex(i) % id
          m_vert(k) = m_vert(k) + 1
          element % vertex(i) % rank = m_vert(k)
        end do

        do i = 1, 12
          k = element % edge(i) % id
          m_edge(k) = m_edge(k) + 1
          element % edge(i) % rank = m_edge(k)
        end do

        do i = 1, 6
          k = element % face(i) % id
          m_face(k) = m_face(k) + 1
          element % face(i) % rank = m_face(k)
        end do

      end associate
    end do

    do e = 1, mesh % n_ghost
      associate(ghost => mesh % ghost(e))

        do i = 1, 8
          k = ghost % vertex(i) % id
          if (k > 0) then
            m_vert(k) = m_vert(k) + 1
            ghost % vertex(i) % rank = m_vert(k)
          end if
        end do

        do i = 1, 12
          k = ghost % edge(i) % id
          if (k > 0) then
            m_edge(k) = m_edge(k) + 1
            ghost % edge(i) % rank = m_edge(k)
          end if
        end do

        do i = 1, 6
          k = ghost % face(i) % id
          if (k > 0) then
            m_face(k) = m_face(k) + 1
            ghost % face(i) % rank = m_face(k)
          end if
        end do

      end associate
    end do

    ! valency ..................................................................

    do e = 1, mesh % n_elem
      associate(element => mesh % element(e))

        do i = 1, 8
          k = element % vertex(i) % id
          element % vertex(i) % val = m_vert(k)
          max_vert_val = max(max_vert_val, m_vert(k))
        end do

        do i = 1, 12
          k = element % edge(i) % id
          element % edge(i) % val = m_edge(k)
          max_edge_val = max(max_edge_val, m_edge(k))
        end do

        do i = 1, 6
          element % face(i) % val = m_face(element % face(i) % id)
        end do

      end associate
    end do

    do e = 1, mesh % n_ghost
      associate(ghost => mesh % ghost(e))

        do i = 1, 8
          k = ghost % vertex(i) % id
          if (k > 0) then
            ghost % vertex(i) % val = m_vert(k)
            max_vert_val = max(max_vert_val, m_vert(k))
          end if
        end do

        do i = 1, 12
          k = ghost % edge(i) % id
          if (k > 0) then
            ghost % edge(i) % val = m_edge(k)
            max_edge_val = max(max_edge_val, m_edge(k))
          end if
        end do

        do i = 1, 6
          k = ghost % face(i) % id
          if (k > 0) then
            ghost % face(i) % val = m_face(k)
          end if
        end do

      end associate
    end do

    mesh % max_vert_val = max_vert_val
    mesh % max_edge_val = max_edge_val

  end subroutine IdentifyRanks

  !=============================================================================

end submodule MP_IdentifyRanks

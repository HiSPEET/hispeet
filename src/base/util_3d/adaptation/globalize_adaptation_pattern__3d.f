module Globalize_Adaptation_Pattern__3D
  use Mesh__3D
  use Mesh_Element_Indexing__3D
  use Element_Transfer_Buffer__3D

  implicit none
  private

  public :: GlobalizeAdaptationPattern_3D

contains

  !-----------------------------------------------------------------------------
  !> Make adaptation pattern globally consistent
  !>
  !> Upgrades the adaptation mark to `max(mark,0)` for all active elements which
  !> possess a neighbor with `mark > 0`

  subroutine GlobalizeAdaptationPattern_3D(mesh)
    class(Mesh_3D), intent(inout) :: mesh

    integer, allocatable, target :: mark(:)
    integer, contiguous, pointer :: mark_val(:,:,:,:)
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: mark_buf

    integer :: e, i, j, k, l, m, n, o

    allocate(mark(mesh%n_elem + mesh%n_ghost), source = -1)

    ! transfer adaptation marks to ghosts ......................................

    ! extract adaptation marks
    do e = 1, mesh % n_elem
      if (mesh % element(e) % frozen) cycle
      mark(e) = mesh % element(e) % adaptation % mark
    end do

    ! transfer marks to ghosts
    if (mesh % n_ghost > 0) then
      mark_val(1:1,1:1,1:1,1:mesh%n_elem+mesh%n_ghost) => mark
      mark_buf = ElementTransferBuffer_3D(mesh, mark_val)
      call mark_buf % Transfer(mesh, mark_val, tag = 1437)
      call mark_buf % Merge(mark_val)
    end if

    ! upgrade according to neighbor refinement .................................

    NEIGHBORS: do e = 1, mesh % n_elem
      associate( element  => mesh % element(e)            &
               , neighbor => mesh % element(e) % neighbor )

        if (element % frozen .or. element % adaptation % mark > 0) cycle

        ! check face neighbors
        do k = 1, 6
          i = element % face(k) % i_neighbor
          if (i <= 0) cycle
          m = mark(neighbor(i)%id)
          if (m <= 0) cycle
          element % adaptation % mark = max(element % adaptation % mark, 0)
          if (m == 0) cycle
          ! coupled neighbor face
          l = neighbor(i)%component
          ! bit-encoded neighbor child refinement
          o = mod(m,1000)
          if (any( BTest(o, V_FACE(:,l) - 1) )) then
            element % adaptation % mark = 1000
            cycle NEIGHBORS
          end if
        end do

        ! check edge neighbors
        do k = 1, 12
          i = element % edge(k) % i_neighbor
          n = element % edge(k) % n_neighbor
          if (n < 1) cycle
          do j = i, i + n - 1
            m = mark(neighbor(j)%id)
            if (m <= 0) cycle
            element % adaptation % mark = max(element % adaptation % mark, 0)
            if (m == 0) cycle
            ! coupled neighbor edge
            l = ElementEdgeID(neighbor(j)%component)
            ! bit-encoded neighbor child refinement
            o = mod(m,1000)
            if (any( BTest(o, V_EDGE(:,l) - 1) )) then
              element % adaptation % mark = 1000
              cycle NEIGHBORS
            end if
          end do
        end do

        ! check vertex neighbors
        do k = 1, 8
          i = element % vertex(k) % i_neighbor
          n = element % vertex(k) % n_neighbor
          if (n < 1) cycle
          do j = i, i + n - 1
            m = mark(neighbor(j)%id)
            if (m <= 0) cycle
            element % adaptation % mark = max(element % adaptation % mark, 0)
            if (m == 0) cycle
            ! coupled neighbor vertex
            l = ElementVertexID(neighbor(j)%component)
            ! bit-encoded neighbor child refinement
            o = mod(m,1000)
            if (BTest(o, l - 1)) then
              element % adaptation % mark = 1000
              cycle NEIGHBORS
            end if
          end do
        end do

      end associate
    end do NEIGHBORS

    ! upgrade according to own refinement .......................................

    do e = 1, mesh % n_elem
      associate(mark => mesh % element(e) % adaptation % mark)
        if (mark > 0) then
          mark = (mark / 1000 + min(1, mod(mark,1000))) * 1000
        end if
      end associate
    end do

  end subroutine GlobalizeAdaptationPattern_3D

  !=============================================================================

end module Globalize_Adaptation_Pattern__3D

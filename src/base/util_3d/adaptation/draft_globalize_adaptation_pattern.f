
  use Mesh__3D
  use Element_Transfer_Buffer__3D

  implicit none
  private

contains


  subroutine GlobalizeAdaptationPattern(mesh)
    class(Mesh_3D), intent(inout) :: mesh

    integer, allocatable, target :: mark(:)
    integer, contiguous, pointer :: mark_val(:,:,:,:)
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: tp_elem_buf

    integer :: e

    allocate(mark(old_mesh%n_elem + old_mesh%n_ghost), source = -1)

    ! extract adaptation marks
    do e = 1, mesh % n_elem
      if (mesh % element(e) % frozen) cycle
      mark(e) = mesh % element(e) % adaptation % mark
    end do

    ! transfer marks to ghosts
    if (mesh % n_ghost > 0) then
      mark_val(1:1,1:1,1:1,1:mesh%n_elem+mesh%n_ghost) => mark
      mark_buf = ElementTransferBuffer_3D(old_mesh, mark_val)
      call mark_buf % Transfer(mesh, mark_val, tag = 1437)
      call mark_buf % Merge(mark_val)
    end if

    ! upgrade marks
    do e = 1, mesh % n_elem
      if (mesh % element(e) % frozen) cycle
      associate(neighbor => mesh % element(e) % neighbor)
        Neighbors: do i = 1, size(neighbor)
          if (mark(neighbor(i)%id) > 0) then
            mark(e) = max(mark(e), 0)
            exit Neighbors
          end if
        end do Neighbors
      end
    end do

  end subroutine GlobalizeAdaptationPattern

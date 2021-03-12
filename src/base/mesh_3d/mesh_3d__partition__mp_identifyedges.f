!> summary:  Identification of mesh edges
!> author:   Joerg Stiller
!> date:     2020/12/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_3d__Partition) MP_IdentifyEdges
  use Quick_Sort
  use Mesh_3d__Element_Indexing, only: V_EDGE
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Identification of mesh edges
  !>
  !> Requires
  !>   - mesh % element % vertex % id
  !>
  !> Generates
  !>   - mesh % n_edge
  !>   - mesh % element % edge % id
  !>   - mesh % element % edge % orientation

  module subroutine IdentifyEdges(mesh)
    class(Mesh3d_Partition), intent(inout) :: mesh !< mesh partition

    ! local data ...............................................................

    integer :: element_edge(5, 12*mesh%n_elem)
    ! encoding:
    ! - element_edge(1,:) :  first vertex ID
    ! - element_edge(2,:) :  second vertex ID
    ! - element_edge(3,:) :  element ID
    ! - element_edge(4,:) :  corresponding element edge (1 .. 12)
    ! - element_edge(5,:) :  orientation (±1)

    integer :: e, i, j, k, l, nee
    integer(IXS) :: o

    ! build ordered list of element edges ......................................

    ! extract element edges, orientation = 1 (aligned) by definition
    j = 0
    do i = 1, size(mesh%element)
      associate(element => mesh % element(i))
        do k = 1, 12
          element_edge(1:5,j+k) = [ element % vertex(V_EDGE(1:2,k)) % id, i, k, 1 ]
        end do
      end associate
      j = j + 12
    end do
    nee = size(element_edge, 2)

    ! flip the edges to achieve an ascending order vertex IDs
    do i = 1, nee
      if (element_edge(1,i) > element_edge(2,i)) then
        call FlipElementEdge(element_edge(:,i))
      end if
    end do

    ! sort the element edges according to 1) first, 2) second vertex ID
    call SortPairs(element_edge)

    ! count and align mesh edges ...............................................

    i = 1
    j = 1
    k = 1
    COUNT_EDGES: do
      !
      ! identify element edges i:j-1 corresponding to mesh edge k
      do
        j = j + 1
        if (j > nee) exit
        if (element_edge(1,j) /= element_edge(1,i)) exit
        if (element_edge(2,j) /= element_edge(2,i)) exit
      end do
      !
      ! set mesh edge orientation to match as many element edges as possible
      call AlignEdge(element_edge(:,i:j-1))
      if (j > nee) exit COUNT_EDGES
      i = j
      k = k + 1
    end do COUNT_EDGES

    mesh % n_edge = k

    ! map element to mesh edges ................................................

    i = 1
    j = 1
    k = 1
    MAP_EDGES: do
      do
        l = element_edge(3,j)           ! corresponding element ID
        e = element_edge(4,j)           ! element edge
        o = int(element_edge(5,j), IXS) ! element edge orientation
        mesh % element(l) % edge(e) % id = k
        mesh % element(l) % edge(e) % orientation = o
        j = j + 1
        if (j > nee) exit
        if (element_edge(1,j) /= element_edge(1,i)) exit
        if (element_edge(2,j) /= element_edge(2,i)) exit
      end do
      if (j > nee) exit MAP_EDGES
      i = j
      k = k + 1
    end do MAP_EDGES

  end subroutine IdentifyEdges

  !---------------------------------------------------------------------------
  !> Reverse edge orientation

  pure subroutine FlipElementEdge(element_edge)
    integer, intent(inout) :: element_edge(5)
    element_edge(1:2) =  element_edge([2,1])  ! swap vertices
    element_edge(5)   = -element_edge(5)      ! switch orientation
  end subroutine FlipElementEdge

  !---------------------------------------------------------------------------
  !> Align edge with as many elements as possible

  pure subroutine AlignEdge(element_edge)
    integer, intent(inout) :: element_edge(:,:)
    integer :: i, n
    n = size(element_edge,2)
    if (2*count(element_edge(5,:) < 0) > n) then
      do i = 1, n
        call FlipElementEdge(element_edge(:,i))
      end do
    end if
  end subroutine AlignEdge

  !=============================================================================

end submodule MP_IdentifyEdges

!> summary:  Identification of mesh edges
!> author:   Joerg Stiller
!> date:     2020/12/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_IdentifyEdges
  use Quick_Sort
  use Mesh_Element_Indexing__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Identification of mesh edges
  !>
  !> Requires
  !>   - mesh % element % vertex % id
  !>   - mesh % element % face   % boundary
  !>
  !> Generates
  !>   - mesh % n_edge
  !>   - mesh % element % edge % id
  !>   - mesh % element % edge % orientation

  module subroutine IdentifyEdges(mesh)
    class(Mesh_3D), intent(inout) :: mesh !< mesh partition

    ! local data ...............................................................

    integer :: element_edge(6, 12*mesh%n_elem)
    ! encoding:
    ! - element_edge(1,:) :  first vertex ID
    ! - element_edge(2,:) :  second vertex ID
    ! - element_edge(3,:) :  polarity indicator
    ! - element_edge(4,:) :  element ID
    ! - element_edge(5,:) :  corresponding element edge (1 .. 12)
    ! - element_edge(6,:) :  orientation (±1)

    integer :: b, e, i, j, k, l, m, p, nee
    integer :: v1, v2, vp(3,8)
    integer(IXS) :: o

    ! build ordered list of element edges ......................................

    ! extract element edges, orientation = 1 (aligned) by definition
    k = 0
    do l = 1, mesh % n_elem
      associate(element => mesh % element(l))

        ! element vertex polarity
        vp = 0
        do m = 1, 6
          b = element % face(m) % boundary
          if (b <= 0) cycle
          p = mesh % boundary(b) % polarity
          select case(p)
          case(-1,1)
            vp(1, V_FACE(:,m)) = sign(100,p)
          case(-2,2)
            vp(2, V_FACE(:,m)) = sign(10,p)
          case(-3,3)
            vp(3, V_FACE(:,m)) = sign(1,p)
          end select
        end do

        ! element edge properties
        do m = 1, 12
          i  = k + m
          v1 = V_EDGE(1,m)
          v2 = V_EDGE(2,m)

          element_edge(1,i) = element%vertex(v1)%id
          element_edge(2,i) = element%vertex(v2)%id
          element_edge(3,i) = 0
          element_edge(4,i) = l
          element_edge(5,i) = m
          element_edge(6,i) = 1

          ! adjust polarity
          do j = 1, 3
            if (vp(j,v1) == vp(j,v2)) then
              ! vertices on same or no extremity
              cycle
            else if (vp(j,v1) == 0) then
              ! vertex 1 inside, vertex 2 on extremity
              element_edge(3,i) = element_edge(3,i) + vp(j,v2)
            else if (vp(j,v2) == 0) then
              ! vertex 1 on extremity, vertex 2 inside
              element_edge(3,i) = element_edge(3,i) + vp(j,v1)
            else if (vp(j,v1) < vp(j,v2)) then
              cycle
            else if (element_edge(1,i) == element_edge(2,i)) then
              ! edge spans over full length, flip to align with polarity
              element_edge(1:2,i) =  element_edge([2,1],i) ! swap vertices
              element_edge( 6 ,i) = -element_edge(  6  ,i) ! switch orientation
            end if
          end do

          ! adjust orientation
          if (element_edge(1,i) > element_edge(2,i)) then
            element_edge(1:2,i) =  element_edge([2,1],i) ! swap vertices
            element_edge( 6 ,i) = -element_edge(  6  ,i) ! switch orientation
          end if

        end do

      end associate
      k = k + 12
    end do
    nee = size(element_edge, 2)

    ! sort the element edges according to
    !   1. first vertex ID
    !   2. second vertex ID
    !   3. polarity indicator
    call SortTriplets(element_edge)

    ! count and align mesh edges ...............................................

    i = 1
    j = 1
    k = 1
    COUNT_EDGES: do
      ! identify element edges i:j-1 corresponding to mesh edge k
      do
        j = j + 1
        if (j > nee) exit
        if (element_edge(1,j) /= element_edge(1,i)) exit ! vertex 1 differs
        if (element_edge(2,j) /= element_edge(2,i)) exit ! vertex 2 differs
        if (element_edge(3,j) /= element_edge(3,i)) exit ! polarity differs
      end do
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
        l = element_edge(4,j)           ! corresponding element ID
        e = element_edge(5,j)           ! element edge
        o = int(element_edge(6,j), IXS) ! element edge orientation
        mesh % element(l) % edge(e) % id = k
        mesh % element(l) % edge(e) % orientation = o
        j = j + 1
        if (j > nee) exit
        if (element_edge(1,j) /= element_edge(1,i)) exit ! vertex 1 differs
        if (element_edge(2,j) /= element_edge(2,i)) exit ! vertex 2 differs
        if (element_edge(3,j) /= element_edge(3,i)) exit ! polarity differs
      end do
      if (j > nee) exit MAP_EDGES
      i = j
      k = k + 1
    end do MAP_EDGES

  end subroutine IdentifyEdges

  !=============================================================================

end submodule MP_IdentifyEdges

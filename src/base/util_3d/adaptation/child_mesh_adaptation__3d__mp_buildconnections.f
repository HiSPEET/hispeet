submodule(Child_Mesh_Adaptation__3D) MP_BuildConnections
  use Quick_Sort
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Build connections between parent, new child, old child and grandchildren
  !>
  !> Entities being built:
  !>
  !> - mapping from parent to new children
  !>     + `parent % n_child`
  !>     + `parent % map_child`
  !>
  !> - mapping from grandchildren to new children
  !>     + `grandchild % n_parent`
  !>     + `grandchild % map_parent`
  !>     + `grandchild % element % adaptation % parent_proc`
  !>     + `grandchild % element % adaptation % parent_id`
  !>
  !>  - mapping between old and new children
  !>     + `rd_send_map` for sending retained data
  !>     + `rd_recv_map` for receiving retained data

  module subroutine BuildConnections( parent, new_child_map, new_child &
                                    , old_child, grandchild            &
                                    , rd_send_map, rd_recv_map         )

    ! arguments ................................................................

    class(Mesh_3D),                 intent(inout) :: parent
    class(ChildDistributionMap_3D), intent(in)    :: new_child_map
    class(Mesh_3D),                 intent(in)    :: new_child
    class(Mesh_3D),       optional, intent(in)    :: old_child
    class(Mesh_3D),       optional, intent(inout) :: grandchild

    type(DataExchangeMap_3D), allocatable, optional, intent(out) :: &
        rd_send_map(:), rd_recv_map(:)

    !$omp master

    ! connect old child and grandchild with new child ..........................

    if (present(old_child) .and. present(rd_send_map)) then
      call ConnectOldChild( parent, new_child_map, new_child   &
                          , old_child, grandchild, rd_send_map )
    end if

    ! connect parent with new child ............................................

    call ConnectParent(parent, new_child_map, new_child)

    ! connect new child with old child .........................................

    if (present(rd_recv_map)) then
      call ConnectNewWithOldChild(new_child, rd_recv_map)
    end if

    !$omp end master
    !$omp barrier

  end subroutine BuildConnections

  !-----------------------------------------------------------------------------
  !> Build connections between new child, old child and grandchildren

  subroutine ConnectOldChild( parent, new_child_map, new_child   &
                            , old_child, grandchild, rd_send_map )

    ! arguments ................................................................

    class(Mesh_3D),                 intent(inout) :: parent
    class(ChildDistributionMap_3D), intent(in)    :: new_child_map
    class(Mesh_3D),                 intent(in)    :: new_child
    class(Mesh_3D),                 intent(in)    :: old_child
    class(Mesh_3D),       optional, intent(inout) :: grandchild

    type(DataExchangeMap_3D), allocatable, intent(out) :: rd_send_map(:)

    ! internal variables .......................................................

    type(DataExchangeMap_3D), allocatable :: send_map(:)
    type(DataExchangeMap_3D), allocatable :: recv_map(:)

    type(DataExchangeSendBuf_3D), allocatable :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable :: recv_buf(:)

    integer, allocatable :: send_attrib(:,:), recv_attrib(:,:)
    integer, allocatable :: cluster_rank(:)
    integer, allocatable :: e(:), m(:)
    integer :: n_recv, n_send, n_attrib
    integer :: a, c, i, j, k, l
    logical :: has_grandchild

    ! initialization ...........................................................

    n_send = parent % n_child
    allocate(send_map(n_send), send_buf(n_send))
    do i = 1, n_send
      send_map(i) = DataExchangeMap_3D(parent % map_child(i))
    end do

    n_recv = old_child % n_parent
    allocate(recv_map(n_recv), recv_buf(n_recv))
    do i = 1, n_recv
      recv_map(i) = DataExchangeMap_3D(old_child % map_parent(i))
    end do

    has_grandchild = present(grandchild)

    ! extract new child attributes for transfer ................................

    if (has_grandchild) then
      n_attrib = 9  ! transfer new child part and IDs
    else
      n_attrib = 1  ! transfer new child part only
    end if

    allocate(send_attrib(n_attrib, parent % n_elem, source = -1))
    do l = 1, parent % n_elem
      if (parent % element(l) % adaptation % refinement /= 100) cycle
      if (parent % element(l) % adaptation % mark       /= 100) cycle
      ! new child target partition
      send_attrib(1,l) = new_child_map % tp_child(l)
      if (has_grandchild) then
        a = 2
        do k = 1, 2
        do j = 1, 2
        do i = 1, 2
          ! new child IDs needed by grandchildren
          send_attrib(a,l) = new_child_map % id_child(i,j,k,l)
          a = a + 1
        end do
        end do
        end do
      end if
    end do

    ! pass new child attributes to old children ................................

    allocate(recv_attrib(n_attrib, old_child % n_cluster), source = -1)

    do i = 1, n_send
      call send_buf(i) % Extract_Data(send_map(i), send_attrib)
      call send_buf(i) % Send_Start()
    end do

    do i = 1, n_recv
      call recv_buf(i) % Init(recv_map(i), recv_attrib)
      call recv_buf(i) % Recv_Start()
    end do

    do i = 1, n_send
      call send_buf(i) % Send_Finish()
    end do

    do i = 1, n_recv
      call recv_buf(i) % Recv_Finish()
      call recv_buf(i) % Assign_Data(recv_attrib)
    end do

    deallocate(recv_buf, send_buf, send_attrib)

    ! create child data send map ...............................................

    ! number of retained element clusters per new partition
    allocate(m(0:new_child%n_parts-1), source = 0)
    associate(tp_child => recv_attrib(1,:))
      do i = 1, old_child % n_cluster
        p = tp_child(i)
        if (p >= 0) then
          m(p) = m(p) + 1
        end if
      end do
    end associate

    ! sort clusters according to 1) parent proc and 2) parent ID
    allocate(cluster_rank(old_child % n_cluster))
    call SortElementClusters(old_child, cluster_rank)

    ! prepare exchange maps
    allocate(rd_send_map(count(m > 0)))
    n = 0
    do p = 0, new_child%n_parts-1
      if (m(p) > 0) then
        n = n + 1
        rd_send_map(n) % comm = new_child % comm_world
        rd_send_map(n) % proc = new_child % proc_part(p)
        allocate(rd_send_map(n) % id_child(8 * m(p)))
        m(p) = n
      end if
    end do

    ! initialize element counter
    allocate(e(size(rd_send_map)), source = 0)

    ! build maps
    do i = 1, old_child % n_cluster
      c = cluster_rank(i)
      p = recv_attrib(c)
      if (p < 0) cycle
      n = m(p)        ! map index
      l = e(n)        ! map element ID offset
      k = 8 * (c - 1) ! mesh element ID offset
      do j = 1, 8
        rd_send_map(n) % id_elem(l + j) = k + j
      end do
      e(n) = e(n) + 8
    end do

    ! connect grandchild .......................................................

    if (has_grandchild) then

      ! convert received new child attributes for transfer to grandchild
      allocate(send_attrib(2, old_child%n_elem), source = -1)
      associate(proc_part => new_child % proc_part)
        do i = 1, old_child % n_elem
          c = old_child % element(i) % cluster_id
          k = old_child % element(i) % cluster_oct + 1
          p = recv_attrib(1,c)
          if (p < 0) cycle ! element not retained
          send_attrib(1,i) = proc_part(p)     ! new element process
          send_attrib(2,i) = recv_attrib(k,c) ! new element ID
        end do
      end associate

      ! received attributes are sent to grandchild
      call ConnectGrandchild(old_child, send_attrib, grandchild)

    end if

  end subroutine ConnectOldChild

  !-----------------------------------------------------------------------------
  !> Sort active clusters according parent process and parent ID
  !>
  !> The returned array `cluster_rank` provides a ranking of the element
  !> clusters. Its first part contains a list of the active clusters sorted
  !> according to 1) the parent element process and 2) the parent element ID.
  !> The second part lists the frozen clusters in the original order.
  !> It is exploited that the active clusters always precede the frozen ones.
  !> As the latter are not retained, they do not need to be sorted.

  subroutine SortElementClusters(mesh, cluster_rank)
    class(Mesh_3D), intent(in) :: mesh
    integer, intent(out) :: cluster_rank(mesh%n_cluster)

    integer, allocatable :: attrib(:,:)
    integer :: n_cluster_active
    integer :: c, e

    n_cluster_active = mesh % n_elem_active / 4

    allocate(attrib(3,n_cluster_active))

    do c = 1, n_cluster_active
      e = 1 + 8 * (c - 1)
      attrib(1,c) = mesh % element(e) % adaptation % parent_proc
      attrib(2,c) = mesh % element(e) % adaptation % parent_id
      attrib(3,c) = c
    end do

    ! sort active clusters according to parent process and parent ID
    call SortPairs(attrib)

    ! build list of active clusters with ascending rank
    do c = 1, n_cluster_active
      cluster_rank(c) = attrib(3,c)
    end do

    ! append unsorted frozen clusters
    do c = n_cluster_active+1, mesh % n_cluster
      cluster_rank(c) = c
    end do

  end subroutine SortElementClusters

  !-----------------------------------------------------------------------------
  !> Update parent info and maps in grandchild partitions

  subroutine ConnectGrandchild(old_child, send_attrib, grandchild)
    class(Mesh_3D), intent(in)    :: old_child
    integer,        intent(in)    :: send_attrib
    class(Mesh_3D), intent(inout) :: grandchild

    type(DataExchangeMap_3D), allocatable :: send_map(:)
    type(DataExchangeMap_3D), allocatable :: recv_map(:)

    type(DataExchangeSendBuf_3D), allocatable :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable :: recv_buf(:)

    integer, allocatable :: recv_attrib(:,:)
    integer :: n_recv, n_send, n_attrib
    integer :: i

    ! initialization ...........................................................

    n_send = old_child % n_child
    allocate(send_map(n_send), send_buf(n_send))
    do i = 1, n_send
      send_map(i) = DataExchangeMap_3D(old_child % map_child(i))
    end do

    n_recv = grandchild % n_parent
    allocate(recv_map(n_recv), recv_buf(n_recv))
    do i = 1, n_recv
      recv_map(i) = DataExchangeMap_3D(grandchild % map_parent(i))
    end do

    n_attrib = size(send_attrib, 1)
    allocate(recv_attrib(n_attrib, grandchild % n_cluster), source = -1)

    ! pass new child attributes to grandchildren ...............................

    do i = 1, n_send
      call send_buf(i) % Extract_Data(send_map(i), send_attrib)
      call send_buf(i) % Send_Start()
    end do

    do i = 1, n_recv
      call recv_buf(i) % Init(recv_map(i), recv_attrib)
      call recv_buf(i) % Recv_Start()
    end do

    do i = 1, n_send
      call send_buf(i) % Send_Finish()
    end do

    do i = 1, n_recv
      call recv_buf(i) % Recv_Finish()
      call recv_buf(i) % Assign_Data(recv_attrib)
    end do

    ! extract received attributes ..............................................

    do i = 1, grandchild % n_elem
      associate(element => grandchild % element(i))
        element % adaptation % parent_proc = recv_attrib(1, element%cluster_id)
        element % adaptation % parent_id   = recv_attrib(2, element%cluster_id)
      end associate
    end do

    ! rebuild map from grandchild to parent .........................................

    call grandchild % BuildMapToParent()

  end subroutine ConnectGrandchild

  !-----------------------------------------------------------------------------
  !> Connect parent with new child

  subroutine ConnectParent(parent, new_child_map, new_child)
    class(Mesh_3D),                 intent(inout) :: parent
    class(Mesh_3D),                 intent(in)    :: new_child
    class(ChildDistributionMap_3D), intent(in)    :: new_child_map

    integer :: i

    if (parent % part >= 0) then

      associate(proc_part => new_child % proc_part)
        do i = 1, parent % n_elem
          associate(adaptation => parent % element(i) % adaptation)
            adaptation % refinement = adaptation % mark
            if (adaptation % refinement > 0) then
              adaptation % child_proc = proc_part(new_child_map % tp_child(i))
            else
              adaptation % child_proc = -1
            end if
            adaptation % mark = 0
          end associate
        end do
      end associate

      call parent % BuildMapToChild()

    end if
  end subroutine ConnectParent

  !-----------------------------------------------------------------------------
  !> Build map indicating the old location of retained child elements
  !>
  !> Only active child elements can be retained. This implies that old and the
  !> new child both emerge from regular refinement and are not frozen.

  subroutine ConnectNewWithOldChild(new_child, rd_recv_map)
    class(Mesh_3D), intent(in) :: new_child
    type(DataExchangeMap_3D), allocatable, intent(out) :: rd_recv_map(:)

    integer, allocatable :: e(:), m(:)
    integer :: i, p, n, n_proc

    call MPI_Comm_size(new_child % comm_world)

    allocate(e(0:n_proc-1), source = 0)
    allocate(m(0:n_proc-1), source = 0)

    ! count retained elements per process
    do i = 1, new_child % n_elem_active
      ! process of retained child is stored in mark
      p = new_child % element(i) % adaptation % mark
      if (p >= 0) then
        m(p) = m(p) + 1
      end if
    end do

    ! prepare exchange maps
    allocate(rd_recv_map(count(m > 0)))
    n = 0
    do p = 0, n_proc-1
      if (m(p) > 0) then
        n = n + 1
        rd_recv_map(n) % comm = new_child % comm_world
        rd_recv_map(n) % proc = p
        allocate(rd_recv_map(n) % id_child(m(p)))
        m(p) = n
      end if
    end do

    ! build maps
    do i = 1, new_child % n_elem_active
      p = new_child % element(i) % adaptation % mark
      if (p < 0) cycle
      n    = m(p)
      e(n) = e(n) + 1
      rd_recv_map(n) % id_elem(e(n)) = i
    end do

  end subroutine ConnectNewWithOldChild

  !=============================================================================

end submodule MP_BuildConnections

!> summary:  Generation of 3d mesh links
!> author:   Joerg Stiller
!> date:     2017/04/24; revised 2020/11/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_BuildLinks

  use Quick_Sort
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of mesh links from global element neighbor information
  !>
  !> On entry, the mesh elements must be complete and `mesh%element%neighbor%id`
  !> set to the local ID of the neighbors within their home-partitions.
  !>
  !> Using this information
  !>
  !>   - `mesh % link` :
  !>      is built
  !>
  !>   – `mesh % element % neighbor % id` :
  !>      entries referring to remote neighbors are translated into ghost IDs
  !>
  !>   - `mesh % face % element` :
  !>     is completed by inserting the ghost ID of adjacent remote elements

  module subroutine BuildLinks(mesh)
    class(Mesh_3D), intent(inout) :: mesh !< local partition

    integer, allocatable :: nf(:), nm(:), ng(:), og(:)
    integer, allocatable :: link_face(:,:)
    integer, allocatable :: link_master(:,:,:)
    integer, allocatable :: link_ghost(:,:)
    integer, allocatable :: map(:), perm(:)
    integer :: c, e, f, i, j, j1, j2, k, l, m, n, p, q, r, s

    !---------------------------------------------------------------------------
    ! Initialization
!### CHECK
print '(9(G0,1X))', 'proc =',mesh%proc,' LINKS  0'
!### CHECK END

    if (mesh%part < 0) then
      allocate(mesh%link(0))
      return
    end if
!### CHECK
print '(9(G0,1X))', 'proc =',mesh%proc,' LINKS  1'
!### CHECK END

    p = mesh%n_parts - 1

    allocate(nf(0:p), source = 0)
    allocate(nm(0:p), source = 0)
    allocate(ng(0:p), source = 0)
    allocate(og(0:p), source = 0)

    !---------------------------------------------------------------------------
    ! Identification of linked entities
    !
    ! In periodic domains with only two element layers, two elements can be
    ! linked twice (at opposite sides). To match the face links to each other,
    ! the entries must ordered in an unique manner. This is achieved by sorting
    ! according to the adjoining element ID from the lowest rank partition

    ! faces ....................................................................

    ! link_face(1,:) : remote partition ID
    ! link_face(2,:) : adjacent element ID from lowest rank partition
    ! link_face(3,:) : local face ID

    allocate(link_face(3, mesh%n_face), source = -1)

    ! identify and count
    k = 0
    do e = 1, mesh % n_elem
      associate(element => mesh % element(e))

        do i = 1, 6

          j = element % face(i) % i_neighbor

          if (j < 1) cycle

          ! identify linked partition
          p = element % neighbor(j) % part

          ! skip adjacent local elements
          if (p == mesh % part) cycle

          ! increase counters
          k     = k     + 1  ! linked faces
          nf(p) = nf(p) + 1  ! faces linked to partition p

          ! store remote partition ID
          link_face(1,k) = p

          ! store adjacent element ID and face from lowest rank partition
          if (p <= mesh%part) then
            link_face(2,k) = element % neighbor(j) % id
          else
            link_face(2,k) = e
          end if

          ! store local face ID
          link_face(3,k) = element % face(i) % id

        end do

      end associate
    end do
!### CHECK
print '(9(G0,1X))', 'proc =',mesh%proc,' LINKS  2: nf =', nf
!### CHECK END

    ! sort linked faces according to
    !   1) remote partition ID
    !   2) element ID
    call SortPairs(link_face(:,1:k))

    ! identify master elements and number of ghosts ............................

    ! nm(p)              : number of elements (masters) linked to partition p
    ! link_master(:,1,e) : IDs of remote partitions to which element e is linked
    ! link_master(:,2,e) : linear index of the linked element component
    ! n                  : number of ghost elements

    ! max number of links per element
    n = 0
    do e = 1, mesh % n_elem
      n = max(n, size(mesh % element(e) % neighbor))
    end do

    allocate(perm(max(n, mesh % n_face)), source = 0)
    allocate(link_master(n, 2, mesh % n_elem), source = -1)

    n = 0
    do e = 1, mesh % n_elem
      associate(element => mesh % element(e))

        k = 0

        ! face neighbors
        do i = 1, 6
          j = element % face(i) % i_neighbor
          if (j < 1) cycle
          p = element % neighbor(j) % part
          if (p == mesh%part) cycle
          k = k + 1
          link_master(k,1,e) = p
          link_master(k,2,e) = i
        end do

        ! edge neighbors
        do i = 1, 12
          j1 = element % edge(i) % i_neighbor
          if (j1 < 1) cycle
          j2 = element % edge(i) % n_neighbor + j1 - 1
          do j = j1, j2
            p = element % neighbor(j) % part
            if (p == mesh%part) cycle
            k = k + 1
            link_master(k,1,e) = p
            link_master(k,2,e) = i + 6
          end do
        end do

        ! vertex neighbors
        do i = 1, 8
          j1 = element % vertex(i) % i_neighbor
          if (j1 < 1) cycle
          j2 = element % vertex(i) % n_neighbor + j1 - 1
          do j = j1, j2
            p = element % neighbor(j) % part
            if (p == mesh%part) cycle
            k = k + 1
            link_master(k,1,e) = p
            link_master(k,2,e) = i + 18
          end do
        end do

        ! sort according to partitions
        call SortIndex(link_master(1:k,1,e), perm(1:k))
        link_master(1:k,1,e) = link_master(perm(1:k),1,e)
        link_master(1:k,2,e) = link_master(perm(1:k),2,e)

        ! count
        q = -1
        do i = 1, k
          p = link_master(i,1,e)
          if (p /= q) then
            nm(p) = nm(p) + 1
            q = p
          end if
        end do

        n = n + k

      end associate
    end do
!### CHECK
print '(9(G0,1X))', 'proc =',mesh%proc,' LINKS  2: nm =', nm
!### CHECK END

    ! ghost elements ...........................................................

    ! link_ghost(1,i) : remote partition ID
    ! link_ghost(2,i) : remote element ID
    ! link_ghost(3,i) : linked local element ID
    ! link_ghost(4,i) : linked local element component index
    ! link_ghost(5,i) : linked local element neighbor index
    ! ng(p)           : number of elements (ghosts) linked from partition p

    allocate(link_ghost(5,n), source = -1)

    k = 0

    do e = 1, mesh % n_elem
      associate(element => mesh % element(e))

        ! via faces
        do i = 1, 6
          j = element % face(i) % i_neighbor
          if (j < 1) cycle
          p = element % neighbor(j) % part
          if (p == mesh%part) cycle
          k = k + 1
          link_ghost(1,k) = p
          link_ghost(2,k) = element % neighbor(j) % id
          link_ghost(3,k) = e
          link_ghost(4,k) = i
          link_ghost(5,k) = j
        end do

        ! via edges
        do i = 1, 12
          j1 = element % edge(i) % i_neighbor
          if (j1 < 1) cycle
          j2 = element % edge(i) % n_neighbor + j1 - 1
          do j = j1, j2
            p = element % neighbor(j) % part
            if (p == mesh%part) cycle
            k = k + 1
            link_ghost(1,k) = p
            link_ghost(2,k) = element % neighbor(j) % id
            link_ghost(3,k) = e
            link_ghost(4,k) = i + 6
            link_ghost(5,k) = j
          end do
        end do

        ! via vertices
        do i = 1, 8
          j1 = element % vertex(i) % i_neighbor
          if (j1 < 1) cycle
          j2 = element % vertex(i) % n_neighbor + j1 - 1
          do j = j1, j2
            p = element % neighbor(j) % part
            if (p == mesh%part) cycle
            k = k + 1
            link_ghost(1,k) = p
            link_ghost(2,k) = element % neighbor(j) % id
            link_ghost(3,k) = e
            link_ghost(4,k) = i + 18
            link_ghost(5,k) = j
          end do
        end do

      end associate
    end do

    ! sort according to 1) partition, 2) remote ID, 3) linked local element
    call SortTriplets(link_ghost)

    ! count
    q = -1
    s = -1
    do k = 1, n
      p = link_ghost(1,k)  ! remote partition ID
      r = link_ghost(2,k)  ! remote element ID
      if (p /= q .or. r /= s) then
        ng(p) = ng(p) + 1
        q = p
        s = r
      end if
    end do
!### CHECK
print '(9(G0,1X))', 'proc =',mesh%proc,' LINKS  2: ng =', ng
!### CHECK END

    !---------------------------------------------------------------------------
    ! generate links

    ! set up data structures ...................................................

    n = count(nf > 0 .or. nm > 0 .or. ng > 0)

    allocate(map(0:mesh%n_parts-1), source = -1)
    allocate(mesh%link(n))
    mesh % n_link = n

    k = 0
    do p = 0, mesh%n_parts - 1
      if (nf(p) > 0 .or. nm(p) > 0 .or. ng(p) > 0) then
        k = k + 1
        map(p) = k
        mesh % link(k) % part     = p
        mesh % link(k) % n_face   = nf(p)
        mesh % link(k) % n_master = nm(p)
        mesh % link(k) % n_ghost  = ng(p)
        allocate(mesh % link(k) % face   ( nf(p) ))
        allocate(mesh % link(k) % master ( nm(p) ))
        allocate(mesh % link(k) % ghost  ( ng(p) ))
      end if
    end do

    ! face links ...............................................................

    nf = 0

    do i = 1, mesh % n_face
      p = link_face(1,i)
      if (p < 0) exit
      nf(p) = nf(p) + 1
      mesh % link(map(p)) % face( nf(p) ) = link_face(3,i)
    end do

    ! master elements ..........................................................

    nm = 0
    do e = 1, mesh % n_elem
      q = -1
      do k = 1, size(link_master,1)

        p = link_master(k,1,e)
        j = link_master(k,2,e)

        if (p == -1) then
          exit
        else if (p /= q) then
          nm(p) = nm(p) + 1
          q = p
        end if

        associate(master => mesh % link( map(p) ) % master( nm(p) ))

          master % id = e

          select case(j)
          case(1:6)
            master%face(j) = 1
          case(7:18)
            master%edge(j-6) = 1
          case default
            master%vertex(j-18) = 1
          end select

        end associate

      end do
    end do

    ! ghost elements ...........................................................

    ! number of ghost elements
    mesh % n_ghost = sum(ng)

    ! offsets for numbering
    og(0) = mesh % n_elem
    do p = 1, mesh % n_parts - 1
      og(p) = og(p-1) + ng(p-1)
    end do

    ! reset counters
    ng = 0

    q = -1
    s = -1
    do k = 1, size(link_ghost,2)

      p = link_ghost(1,k)  ! remote partition ID
      r = link_ghost(2,k)  ! remote element ID
      e = link_ghost(3,k)  ! linked local element ID
      i = link_ghost(4,k)  ! linked local element component index
      j = link_ghost(5,k)  ! linked local element neighbor index

      if (p /= q .or. r /= s) then
        ng(p) = ng(p) + 1
        q = p
        s = r
      end if

      associate( ghost   => mesh % link(map(p)) % ghost(ng(p)) &
               , element => mesh % element(e)                  )

        l = og(p) + ng(p)

        ghost % id = l
        element % neighbor(j) % id = l

        if (i > 18) then ! vertex
          c = ElementVertexID(element % neighbor(j) % component)
          ghost % vertex(c) = 1

        else if (i > 6) then ! edge
          c = ElementEdgeID(element % neighbor(j) % component)
          ghost % edge(c) = 1

        else ! face
          c = ElementFaceID(element % neighbor(j) % component)
          ghost % face(c) = 1
          m = mod(i,2)
          f = element % face(i) % id
          if (       m == 1 .and. element % face(i) % normal > 0    &
               .or.  m == 0 .and. element % face(i) % normal < 0  ) &
          then
            mesh % face(f) % element(1) % id   = l
            mesh % face(f) % element(1) % face = c
          else
            mesh % face(f) % element(2) % id   = l
            mesh % face(f) % element(2) % face = c
          end if

        end if

      end associate

    end do
!### CHECK
print '(9(G0,1X))', 'proc =',mesh%proc,' LINKS  X'
!### CHECK END

  end subroutine BuildLinks

  !=============================================================================

end submodule MP_BuildLinks

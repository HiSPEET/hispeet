!> summary:  Generation of Cartesian mesh links
!> author:   Joerg Stiller
!> date:     2017/04/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Generation of Cartesian mesh links
!===============================================================================

module CART__Generate_Mesh_Links

  use Quick_Sort
  use CART__Element_Connectivity
  use CART__Mesh_Partition

  implicit none
  private

  public :: GenerateMeshLinks

contains

!-------------------------------------------------------------------------------
!> Generation of mesh links from global element neighbor information

subroutine GenerateMeshLinks(mesh, conn)
  type(MeshPartition),       intent(inout) :: mesh    !< mesh partition
  type(ElementConnectivity), intent(in)    :: conn(:) !< element connectivity

  integer, allocatable :: nf(:), nm(:), ng(:), og(:)
  integer, allocatable :: link_face(:,:)
  integer, allocatable :: link_master(:,:,:)
  integer, allocatable :: link_ghost(:,:)
  integer, allocatable :: map(:)
  integer :: e, i, j, k, l, n, p, q, r, s
  integer :: nef(6), perm(26)

  !-----------------------------------------------------------------------------
  ! initialization

  if (mesh%part < 0 .or. mesh%n_part == 1) then
    allocate(mesh%link(0))
    return
  end if

  p = mesh%n_part - 1

  allocate(nf(0:p), source = 0)
  allocate(nm(0:p), source = 0)
  allocate(ng(0:p), source = 0)
  allocate(og(0:p), source = 0)

  !-----------------------------------------------------------------------------
  ! identify and count linked entities

  ! faces ......................................................................

  ! In periodic domains with only two element layers, two elements can be
  ! linked twice (at opposite sides). To match the face links to each other,
  ! the entries must ordered in an unique manner. This is achieved by sorting
  ! according
  !            1) the adjoining element ID from the lowest rank partition
  !            2) the adjacent face of this element

  allocate(link_face(4, mesh%nf), source = -1)

  ! map to corresponding neigbor element face
  nef = [ 2, 1, 4, 3, 6, 5]

  ! identify and count
  k = 0
  do e = 1, mesh%ne
    do i = 1, 6
      p = conn(e)%face(i)%part
      if (p >= 0 .and. p /= mesh%part) then
        ! increase counters
        k = k + 1
        nf(p) = nf(p) + 1
        ! store partition of remote element
        link_face(1,k) = p
        ! store adjacent element ID and face from lowest rank partition
        if (p < mesh%part) then
          link_face(2,k) = conn(e)%face(i)%id
          link_face(3,k) = nef(i)
        else
          link_face(2,k) = e
          link_face(3,k) = i
        end if
        ! store local face ID
        link_face(4,k) = mesh%element(e)%face(i)%id
      end if
    end do
  end do

  ! sort linked faces according to 1) remote partition, 2) element ID
  call SortTriplets(link_face(:,1:k))

  ! partitions needing a copy of local elements  ...............................

  allocate(link_master(26, 2, mesh%ne), source = -1)

  n = 0
  do e = 1, mesh%ne

    k = 0

    ! face neighbors
    do i = 1, 6
      p = conn(e)%face(i)%part
      if (p >= 0 .and. p /= mesh%part) then
        k = k + 1
        link_master(k,1,e) = p
        link_master(k,2,e) = i
      end if
    end do

    ! edge neighbors
    do i = 1, 12
      p = conn(e)%edge(i)%part
      if (p >= 0 .and. p /= mesh%part) then
        k = k + 1
        link_master(k,1,e) = p
        link_master(k,2,e) = i + 6
      end if
    end do

    ! vertex neighbors
    do i = 1, 8
      p = conn(e)%vertex(i)%part
      if (p >= 0 .and. p /= mesh%part) then
        k = k + 1
        link_master(k,1,e) = p
        link_master(k,2,e) = i + 18
      end if
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

  end do

  ! ghost elements .............................................................

  allocate(link_ghost(4,n))

  k = 0

  do e = 1, mesh%ne

    ! via faces
    do i = 1, 6
      p = conn(e)%face(i)%part
      if (p >= 0 .and. p /= mesh%part) then
        k = k + 1
        link_ghost(1,k) = p
        link_ghost(2,k) = conn(e)%face(i)%id
        link_ghost(3,k) = e
        link_ghost(4,k) = i
      end if
    end do

    ! via edges
    do i = 1, 12
      p = conn(e)%edge(i)%part
      if (p >= 0 .and. p /= mesh%part) then
        k = k + 1
        link_ghost(1,k) = p
        link_ghost(2,k) = conn(e)%edge(i)%id
        link_ghost(3,k) = e
        link_ghost(4,k) = i + 6
      end if
    end do

    ! via vertices
    do i = 1, 8
      p = conn(e)%vertex(i)%part
      if (p >= 0 .and. p /= mesh%part) then
        k = k + 1
        link_ghost(1,k) = p
        link_ghost(2,k) = conn(e)%vertex(i)%id
        link_ghost(3,k) = e
        link_ghost(4,k) = i + 18
      end if
    end do

  end do

  ! sort according to 1) partition, 2) remote ID, 3) linked local element
  call SortTriplets(link_ghost)

  ! count
  q = -1
  s = -1
  do k = 1, n
    p = link_ghost(1,k)
    r = link_ghost(2,k)
    if (p /= q .or. r /= s) then
      ng(p) = ng(p) + 1
      q = p
      s = r
    end if
  end do

  !-----------------------------------------------------------------------------
  ! generate links

  ! set up data structures .....................................................

  n = count(nf > 0 .or. nm > 0 .or. ng > 0)

  allocate(map(0:mesh%n_part-1), source = -1)
  allocate(mesh%link(n))

  k = 0
  do p = 0, mesh%n_part - 1
    if (nf(p) > 0 .or. nm(p) > 0 .or. ng(p) > 0) then
      k = k + 1
      map(p) = k
      mesh % link(k) % part = p
      mesh % link(k) % nf = nf(p)
      mesh % link(k) % nm = nm(p)
      mesh % link(k) % ng = ng(p)
      allocate(mesh%link(k)%face ( nf(p) ))
      allocate(mesh%link(k)%master( nm(p) ))
      allocate(mesh%link(k)%ghost( ng(p) ))
    end if
  end do

  ! face links .................................................................

  nf = 0

  do i = 1, mesh%nf
    p = link_face(1,i)
    if (p < 0) exit
    nf(p) = nf(p) + 1
    mesh%link(map(p))%face ( nf(p) ) = link_face(4,i)
  end do

  ! master elements ............................................................

  nm = 0
  do e = 1, mesh%ne
    q = -1
    do k = 1, 26

      p = link_master(k,1,e)
      j = link_master(k,2,e)

      if (p == -1) then
        exit
      else if (p /= q) then
        nm(p) = nm(p) + 1
        q = p
      end if

      associate(master => mesh%link( map(p) )%master( nm(p) ))

        master%id = e

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

  ! ghost elements .............................................................

  ! number of ghost elements
  mesh%ng = sum(ng)

  ! offsets for numbering
  og(0) = mesh%ne
  do p = 1, mesh%n_part - 1
    og(p) = og(p-1) + ng(p-1)
  end do

  ! reset counters
  ng = 0

  q = -1
  s = -1
  do k = 1, size(link_ghost,2)

    p = link_ghost(1,k)
    r = link_ghost(2,k)
    e = link_ghost(3,k)
    j = link_ghost(4,k)

    if (p /= q .or. r /= s) then
      ng(p) = ng(p) + 1
      q = p
      s = r
    end if

    associate( ghost   => mesh%link(map(p))%ghost(ng(p)), &
               element => mesh%element(e)                 )

      l = og(p) + ng(p)

      ghost%id  = l

      select case(j)

      case(1) ! west face
        ghost   % face(2) = 1
        element % face(1) % neighbor = l
        mesh%face( element%face(1)%id ) % element(1) = l

      case(2) ! east face
        ghost   % face(1) = 1
        element % face(2) % neighbor = l
        mesh%face( element%face(2)%id ) % element(2) = l

      case(3) ! south face
        ghost   % face(4) = 1
        element % face(3) % neighbor = l
        mesh%face( element%face(3)%id ) % element(1) = l

      case(4) ! north face
        ghost   % face(3) = 1
        element % face(4) % neighbor = l
        mesh%face( element%face(4)%id ) % element(2) = l

      case(5) ! bottom face
        ghost   % face(6) = 1
        element % face(5) % neighbor = l
        mesh%face( element%face(5)%id ) % element(1) = l

      case(6) ! top face
        ghost   % face(5) = 1
        element % face(6) % neighbor = l
        mesh%face( element%face(6)%id ) % element(2) = l

      case(7:10) ! x1 edges
        i = j - 6
        ghost   % edge(5-i) = 1
        element % edge(  i) % neighbor = l

      case(11:14) ! x2 edges
        i = j - 6
        ghost   % edge(13-i) = 1
        element % edge(   i) % neighbor = l

      case(15:18) ! x3 edges
        i = j - 6
        ghost   % edge(21-i) = 1
        element % edge(   i) % neighbor = l

      case default ! vertices
        i = j - 18
        ghost   % vertex(9-i) = 1
        element % vertex(  i) % neighbor = l

      end select

    end associate

  end do

end subroutine GenerateMeshLinks

!===============================================================================

end module CART__Generate_Mesh_Links

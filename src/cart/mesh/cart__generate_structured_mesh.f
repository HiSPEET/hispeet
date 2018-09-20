!> summary:  Generate a distributed, structured Cartesian mesh
!> author:   Joerg Stiller
!> date:     2013/05/24; revised 2014/11/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Generate a distributed, structured Cartesian mesh
!===============================================================================

module CART__Generate_Structured_Mesh

  use Kind_Parameters,   only: RNP
  use Constants,         only: HALF
  use Execution_Control, only: Error
  use Element_Components
  use Structured_Mesh

  use XMPI
  use CART__Element_Connectivity
  use CART__Mesh_Boundary
  use CART__Mesh_Partition
  use CART__Generate_Mesh_Links

  implicit none
  private

  public :: GenerateStructuredMesh
  public :: GenerateStructuredMeshFaces
  public :: GenerateStructuredMeshBoundaries

contains

!-------------------------------------------------------------------------------
!> Creates structured mesh partitions

subroutine GenerateStructuredMesh(mesh, np, ep, xo, dx, periodic, comm)
  type(MeshPartition), intent(out) :: mesh  !< local partition
  integer,   intent(in) :: np(3)       !< number of partitions in directions 1:3
  integer,   intent(in) :: ep(3)       !< elements per partition and direction
  real(RNP), intent(in) :: xo(3)       !< corner closest to -infinity
  real(RNP), intent(in) :: dx(3)       !< mesh spacing in directions 1:3
  logical,   intent(in) :: periodic(3) !< periodic directions set true
  type(MPI_Comm), intent(in) :: comm   !< MPI communicator

  ! local variables ............................................................

  integer :: part       ! partition ID
  integer :: n_part     ! number of non-empty partitions
  integer :: n_boundary ! number of external boundaries
  integer :: ip, jp, kp ! partion triple index
  integer :: n1, n2, n3 ! local mesh dimensions
  integer :: i0, j0, k0 ! offset of local element indices WRT global numbering

  ! element corner closest to -infinity
  real(RNP) :: xe(3)

  ! mesh element connectivity
  type(ElementConnectivity), allocatable :: conn(:)

  ! MPI
  integer :: comm_size
  integer :: rank

  ! initialization .............................................................

  call MPI_Comm_size(comm, comm_size)
  call MPI_Comm_rank(comm, rank)

  ! number of partitions
  n_part = product(np)
  if (n_part > comm_size .and. rank == 0) then
    call Error('GenerateStructuredMesh', 'n_part > MPI communicator size')
  end if

  ! number of non-periodic domain boundaries
  n_boundary = 6

  ! partition ID and triple index
  if (rank < n_part) then
    part = rank
    call TripleIndex(ip, jp, kp, np(1), np(2), 1, 1, 1, l=part+1)
  else
    part = -1
    ip   = -1
    jp   = -1
    kp   = -1
  end if

  ! empty partitions
  if (part < 0) then
    mesh = EmptyMeshPartition(0, n_boundary, n_part, dx, comm=comm)
    call GenerateStructuredMeshBoundaries(mesh, periodic)
    return
  end if

  ! local mesh dimensions
  n1 = ep(1)
  n2 = ep(2)
  n3 = ep(3)

  ! offset of local element indices WRT global numbering
  i0 = ep(1) * (ip - 1)
  j0 = ep(2) * (jp - 1)
  k0 = ep(3) * (kp - 1)

  ! mesh attributes and properties .............................................

  mesh % n_boundary = n_boundary
  mesh % n_part     = n_part
  mesh % part       = part
  mesh % structured = .true.

  mesh % dx = dx

  mesh % ne = n1 * n2 * n3

  mesh % nf = (n1+1) * n2 * n3  &
            + (n2+1) * n3 * n1  &
            + (n3+1) * n1 * n2

  mesh % ne1 = n1
  mesh % ne2 = n2
  mesh % ne3 = n3

  mesh % comm = comm

  ! mesh components ............................................................

  call GenerateStructuredMeshElements
  call GenerateStructuredMeshFaces(mesh)
  call GenerateStructuredMeshBoundaries(mesh, periodic)
  call GenerateMeshLinks(mesh, conn)

contains

  !-----------------------------------------------------------------------------
  !> Builds the local mesh elements

  subroutine GenerateStructuredMeshElements

    integer :: ie, je, ke             ! local element indices
    integer :: ig, jg, kg, ng(3)      ! global elements indices and counts
    integer :: ies, jes, kes          ! shifted element indices
    integer :: ips, jps, kps          ! shifted partition indices
    integer :: i, e, l, p, r, s, t
    integer :: boundary_id(6) = [ 1,2,3,4,5,6 ]


    ! prerequisites ............................................................

    ! global mesh dimensions
    ng = np * ep

    ! generate elements ........................................................

    allocate(mesh%element( mesh%ne ), conn( mesh%ne ))

    do ke = 1, n3
    do je = 1, n2
    do ie = 1, n1

      l = LexicalElementIndex(ie, je, ke, n1, n2, 0)

      associate(element => mesh % element(l))

        ! global element ID ....................................................

        ig = i0 + ie
        jg = j0 + je
        kg = k0 + ke

        element % global_id = LexicalElementIndex(ig, jg, kg, ng(1), ng(2), 0)

        ! vertex coordinates ...................................................

        ! corner closest to -infinity
        xe(1) = xo(1) + (ig - 1) * dx(1)
        xe(2) = xo(2) + (jg - 1) * dx(2)
        xe(3) = xo(3) + (kg - 1) * dx(3)

        element % vertex(1) % x  =  xe  !  [ 0, 0, 0 ] * dx
        element % vertex(2) % x  =  xe  +  [ 1, 0, 0 ] * dx
        element % vertex(3) % x  =  xe  +  [ 0, 1, 0 ] * dx
        element % vertex(4) % x  =  xe  +  [ 1, 1, 0 ] * dx
        element % vertex(5) % x  =  xe  +  [ 0, 0, 1 ] * dx
        element % vertex(6) % x  =  xe  +  [ 1, 0, 1 ] * dx
        element % vertex(7) % x  =  xe  +  [ 0, 1, 1 ] * dx
        element % vertex(8) % x  =  xe  +  [ 1, 1, 1 ] * dx

        ! connectivity .........................................................

        do t = -1, 1
        do s = -1, 1
        do r = -1, 1

          if (r == 0 .and. s == 0 .and. t == 0) cycle

          ies = ie + r
          jes = je + s
          kes = ke + t

          ips = ip
          jps = jp
          kps = kp

          ! shift in direction 1
          if (ies < 1) then
            ies = n1
            ips = ip - 1
            if (ips < 1) then
              if (periodic(1)) then
                ips = np(1)
              else
                cycle
              end if
            end if
          else if (ies > n1) then
            ies = 1
            ips = ip + 1
            if (ips > np(1)) then
              if (periodic(1)) then
                ips = 1
              else
                cycle
              end if
            end if
          end if

          ! shift in direction 2
          if (jes < 1) then
            jes = n2
            jps = jp - 1
            if (jps < 1) then
              if (periodic(2)) then
                jps = np(2)
              else
                cycle
              end if
            end if
          else if (jes > n2) then
            jes = 1
            jps = jp + 1
            if (jps > np(2)) then
              if (periodic(2)) then
                jps = 1
              else
                cycle
              end if
            end if
          end if

          ! shift in direction 3
          if (kes < 1) then
            kes = n3
            kps = kp - 1
            if (kps < 1) then
              if (periodic(3)) then
                kps = np(3)
              else
                cycle
              end if
            end if
          else if (kes > n3) then
            kes = 1
            kps = kp + 1
            if (kps > np(3)) then
              if (periodic(3)) then
                kps = 1
              else
                cycle
              end if
            end if
          end if

          p = LexicalIndex(ips, jps, kps, np(1), np(2), 1, 1, 1)

          p = p - 1
          e = LexicalElementIndex(ies, jes, kes, n1, n2 ,0)
          i = ElementComponentID(r,s,t)

          ! store connectivity
          select case(ElementComponentType(r,s,t))
          case(IS_FACE)
            conn(l) % face(i) % part = p
            conn(l) % face(i) % id   = e
            if (p == mesh%part) then
              element % face(i) % neighbor = e
            end if
          case(IS_EDGE)
            conn(l) % edge(i) % part = p
            conn(l) % edge(i) % id   = e
            if (p == mesh%part) then
              element % edge(i) % neighbor = e
            end if
          case(IS_VERTEX)
            conn(l) % vertex(i) % part = p
            conn(l) % vertex(i) % id   = e
            if (p == mesh%part) then
              element % vertex(i) % neighbor = e
            end if
          end select

        end do
        end do
        end do

        ! face data ............................................................

        ! west
        if (ie == 1 .and. ip == 1) then
          element%face(1)%boundary = boundary_id(1)
        end if
        element%face(1)%id = LexicalFaceIndex(ie-1, je-1, ke-1, 1, n1,n2,n3, 0)

        ! east
        if (ie == n1 .and. ip == np(1)) then
          element%face(2)%boundary = boundary_id(2)
        end if
        element%face(2)%id = LexicalFaceIndex(ie, je-1, ke-1, 1, n1,n2,n3, 0)

        ! south
        if (je == 1 .and. jp == 1) then
          element%face(3)%boundary = boundary_id(3)
        end if
        element%face(3)%id = LexicalFaceIndex(ie-1, je-1, ke-1, 2, n1,n2,n3, 0)

        ! north
        if (je == n2 .and. jp == np(2)) then
          element%face(4)%boundary = boundary_id(4)
        end if
        element%face(4)%id = LexicalFaceIndex(ie-1, je, ke-1, 2, n1,n2,n3, 0)

        ! bottom
        if (ke == 1 .and. kp == 1) then
          element%face(5)%boundary = boundary_id(5)
        end if
        element%face(5)%id = LexicalFaceIndex(ie-1, je-1, ke-1, 3, n1,n2,n3, 0)

        ! top
        if (ke == n3 .and. kp == np(3)) then
          element%face(6)%boundary = boundary_id(6)
        end if
        element%face(6)%id = LexicalFaceIndex(ie-1, je-1, ke, 3, n1,n2,n3, 0)

      end associate

    end do
    end do
    end do

  end subroutine GenerateStructuredMeshElements

end subroutine GenerateStructuredMesh

!-------------------------------------------------------------------------------
!> Creates the faces of a structured mesh partitions

subroutine GenerateStructuredMeshFaces(mesh)
  type(MeshPartition), intent(inout) :: mesh

  integer :: n1, n2, n3
  integer :: i, j, k, l

  ! initialization .............................................................

  if (.not. mesh%structured) then
    call Error('GenerateStructuredMeshFaces', 'structured mesh required')
  end if

  n1 = mesh % ne1
  n2 = mesh % ne2
  n3 = mesh % ne3

  mesh % nf1 = (n1+1) * n2 * n3
  mesh % nf2 = (n2+1) * n3 * n1
  mesh % nf3 = (n3+1) * n1 * n2

  mesh % nf = mesh % nf1  &
            + mesh % nf2  &
            + mesh % nf3

  allocate(mesh%face( mesh%nf ))

  ! x1-faces ...................................................................

  do k = 0, n3-1
  do j = 0, n2-1
  do i = 0, n1

    l = LexicalFaceIndex(i, j, k, 1, n1, n2, n3, 0)

    associate(face => mesh%face(l))

      if (i > 0) then
        face % element(1) = LexicalElementIndex(i, j+1, k+1, n1, n2, 0)
      end if

      if (i < n1) then
        face % element(2) = LexicalElementIndex(i+1, j+1, k+1, n1, n2, 0)
      end if

    end associate

  end do
  end do
  end do

  ! x2-faces ...................................................................

  do k = 0, n3-1
  do j = 0, n2
  do i = 0, n1-1

    l = LexicalFaceIndex(i, j, k, 2, n1, n2, n3, 0)

    associate(face => mesh%face(l))

      if (j > 0) then
        face % element(1) = LexicalElementIndex(i+1, j, k+1, n1, n2, 0)
      end if

      if (j < n2) then
        face % element(2) = LexicalElementIndex(i+1, j+1, k+1, n1, n2, 0)
      end if

    end associate

  end do
  end do
  end do

  ! x3-faces ...................................................................

  do k = 0, n3
  do j = 0, n2-1
  do i = 0, n1-1

    l = LexicalFaceIndex(i, j, k, 3, n1, n2, n3, 0)

    associate(face => mesh%face(l))

      if (k > 0) then
        face % element(1) = LexicalElementIndex(i+1, j+1, k, n1, n2, 0)
      end if

      if (k < n3) then
        face % element(2) = LexicalElementIndex(i+1, j+1, k+1, n1, n2, 0)
      end if

    end associate

  end do
  end do
  end do

end subroutine GenerateStructuredMeshFaces

!-------------------------------------------------------------------------------
!> Creates the boundaries of a structured mesh partitions

subroutine GenerateStructuredMeshBoundaries(mesh, periodic)
  type(MeshPartition), intent(inout) :: mesh
  logical, intent(in) :: periodic(3) !< periodic directions set true

  integer :: b, e, e0, f, i, j, k

  ! initialization .............................................................

  if (allocated(mesh%boundary)) then
    deallocate(mesh%boundary)
  end if

  mesh%n_boundary = 6
  allocate(mesh%boundary(6))

  if (mesh%part < 0) then
    call mesh % boundary(1) % New(1, 'west',   nf = 0, periodic = periodic(1))
    call mesh % boundary(2) % New(2, 'east',   nf = 0, periodic = periodic(1))
    call mesh % boundary(3) % New(3, 'south',  nf = 0, periodic = periodic(2))
    call mesh % boundary(4) % New(4, 'north',  nf = 0, periodic = periodic(2))
    call mesh % boundary(5) % New(5, 'bottom', nf = 0, periodic = periodic(3))
    call mesh % boundary(6) % New(6, 'top',    nf = 0, periodic = periodic(3))
    return
  else if (.not. mesh%structured) then
    call Error('GenerateStructuredMeshBoundaries', 'mesh must be structured')
  end if

  associate( element => mesh%element, &
             n1 => mesh%ne1,          &
             n2 => mesh%ne2,          &
             n3 => mesh%ne3           )

    ! west .....................................................................

    b = 1
    e = 1
    associate(boundary => mesh%boundary(b))
      if (mesh%element(e)%face(b)%boundary == b) then
        call boundary % New(id=1, name='west', nf=n2*n3, periodic=periodic(1))
        f = 1
        do k = 0, n3-1
          e0 = 1 + n1*n2 * k
          do j = 0, n2-1
            e = e0 + n1 * j
            boundary % face(f) % orientation         = -1
            boundary % face(f) % mesh_face    % id   =  element(e)%face(b)%id
            boundary % face(f) % mesh_element % id   =  e
            boundary % face(f) % mesh_element % face =  b
            f = f + 1
          end do
        end do
      else
        call boundary % New(id=1, name='west', nf=0, periodic=periodic(1))
      end if
    end associate

    ! east .....................................................................

    b = 2
    e = n1
    associate(boundary => mesh%boundary(b))
      if (mesh%element(e)%face(b)%boundary == b) then
        call boundary % New(id=2, name='east', nf=n2*n3, periodic=periodic(1))
        f = 1
        do k = 0, n3-1
          e0 = n1 + n1*n2 * k
          do j = 0, n2-1
            e = e0 + n1 * j
            boundary % face(f) % orientation         =  1
            boundary % face(f) % mesh_face    % id   =  element(e)%face(b)%id
            boundary % face(f) % mesh_element % id   =  e
            boundary % face(f) % mesh_element % face =  b
            f = f + 1
          end do
        end do
      else
        call boundary % New(id=2, name='east', nf=0, periodic=periodic(1))
      end if
    end associate

    ! south ....................................................................

    b = 3
    e = 1
    associate(boundary => mesh%boundary(b))
      if (mesh%element(e)%face(b)%boundary == b) then
        call boundary % New(id=3, name='south', nf=n1*n3, periodic=periodic(2))
        f = 1
        do k = 0, n3-1
          e0 = 1 + n1*n2 * k
          do i = 0, n1-1
            e = e0 + i
            boundary % face(f) % orientation         = -2
            boundary % face(f) % mesh_face    % id   =  element(e)%face(b)%id
            boundary % face(f) % mesh_element % id   =  e
            boundary % face(f) % mesh_element % face =  b
            f = f + 1
          end do
        end do
      else
        call boundary % New(id=3, name='south', nf=0, periodic=periodic(2))
      end if
    end associate

    ! north ....................................................................

    b = 4
    e = 1 + n1*(n2-1)
    associate(boundary => mesh%boundary(b))
      if (mesh%element(e)%face(b)%boundary == b) then
        call boundary % New(id=4, name='north', nf=n1*n3, periodic=periodic(2))
        f = 1
        do k = 0, n3-1
          e0 = 1 + n1*(n2-1) + n1*n2 * k
          do i = 0, n1-1
            e = e0 + i
            boundary % face(f) % orientation         =  2
            boundary % face(f) % mesh_face    % id   =  element(e)%face(b)%id
            boundary % face(f) % mesh_element % id   =  e
            boundary % face(f) % mesh_element % face =  b
            f = f + 1
          end do
        end do
      else
        call boundary % New(id=4, name='north', nf=0, periodic=periodic(2))
      end if
    end associate

    ! bottom ...................................................................

    b = 5
    e = 1
    associate(boundary => mesh%boundary(b))
      if (mesh%element(e)%face(b)%boundary == b) then
        call boundary % New(id=5, name='bottom', nf=n1*n2, periodic=periodic(3))
        f = 1
        do j = 0, n2-1
          e0 = 1 + n1 * j
          do i = 0, n1-1
            e = e0 + i
            boundary % face(f) % orientation         = -3
            boundary % face(f) % mesh_face    % id   =  element(e)%face(b)%id
            boundary % face(f) % mesh_element % id   =  e
            boundary % face(f) % mesh_element % face =  b
            f = f + 1
          end do
        end do
      else
        call boundary % New(id=5, name='bottom', nf=0, periodic=periodic(3))
      end if
    end associate

    ! top ......................................................................

    b = 6
    e = 1 + n1*n2*(n3-1)
    associate(boundary => mesh%boundary(b))
      if (mesh%element(e)%face(b)%boundary == b) then
        call boundary % New(id=6, name='top', nf=n1*n2, periodic=periodic(3))
        f = 1
        do j = 0, n2-1
          e0 = 1 + n1 * j + n1*n2*(n3-1)
          do i = 0, n1-1
            e = e0 + i
            boundary % face(f) % orientation         =  3
            boundary % face(f) % mesh_face    % id   =  element(e)%face(b)%id
            boundary % face(f) % mesh_element % id   =  e
            boundary % face(f) % mesh_element % face =  b
            f = f + 1
          end do
        end do
      else
        call boundary % New(id=6, name='top', nf=0, periodic=periodic(3))
      end if
    end associate

  end associate

end subroutine GenerateStructuredMeshBoundaries

!===============================================================================

end module CART__Generate_Structured_Mesh

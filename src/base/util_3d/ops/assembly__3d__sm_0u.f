!> summary:  Assembly of unstructured 3D mesh variables
!> author:   Joerg Stiller
!> date:     2021/03/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Assembly__3D) SM_0U
  use Mesh_Element__3D
  use Mesh_Element_Indexing__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Unstructured assembly and, optionally, averaging of a scalar variable

  module subroutine Assembly_3D_0U(mesh, u, u_buf, avg)
    !> mesh partition
    class(Mesh_3D), intent(in)    :: mesh
    !> mesh variable, including ghost entries
    real(RNP), intent(inout) :: u(:,:,:,:)
    !> MPI transfer buffer
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: u_buf
    !> switch for averaging over element boundaries [F]
    logical, optional, intent(in) :: avg

    real(RNP), allocatable, save :: u_face(:,:,:,:), u_edge(:,:,:), u_vert(:,:)
    integer  , allocatable, save :: m_face(:)      , m_edge(:)    , m_vert(:)

    integer :: ni

    ! initialization ...........................................................

    !$omp master
    ni = size(u,1) - 2
    allocate(m_face(mesh%n_face), u_face(ni,ni,              2, mesh%n_face))
    allocate(m_edge(mesh%n_edge), u_edge(ni, mesh%max_edge_val, mesh%n_edge))
    allocate(m_vert(mesh%n_vert), u_vert(    mesh%max_vert_val, mesh%n_vert))
    !$omp end master
    !$omp barrier

    ! start transfer of ghost data .............................................

    call u_buf % Transfer(mesh, u, tag=10303)

    ! extract contributions of local elements ..................................

    call ExtractTraces( mesh % element, u, u_face, u_edge, u_vert &
                      , m_face , m_edge, m_vert                   )

    ! finish transfer and inject ghost data ....................................

    call u_buf % Merge(u)

    ! extract contributions data of ghost elements .............................

    call ExtractTraces(mesh % ghost, u, u_face, u_edge, u_vert)

    ! add and, optionally, average contributions ...............................

    call AssembleTraces(u_face, u_edge, u_vert, m_face, m_edge, m_vert, avg)

    ! distribute skeleton data to elements .....................................

    call DistributeSkeletonData(mesh % element, u, u_face, u_edge, u_vert)

    ! finalization .............................................................

    !$omp master
    deallocate(m_face, u_face)
    deallocate(m_edge, u_edge)
    deallocate(m_vert, u_vert)
    !$omp end master

  end subroutine Assembly_3D_0U

  !-----------------------------------------------------------------------------
  !> Extraction of scalar element boundary data

  subroutine ExtractTraces( element, u, u_face, u_edge, u_vert &
                          , m_face , m_edge, m_vert            )

    class(MeshElement_3D), intent(in)    :: element(:)
    real(RNP)            , intent(in)    :: u(0:,0:,0:,:)
    real(RNP)            , intent(inout) :: u_face(:,:,:,:)
    real(RNP)            , intent(inout) :: u_edge(:,:,:)
    real(RNP)            , intent(inout) :: u_vert(:,:)
    integer,     optional, intent(inout) :: m_face(:)
    integer,     optional, intent(inout) :: m_edge(:)
    integer,     optional, intent(inout) :: m_vert(:)

    integer :: i, k, l, p, q

    p = ubound(u,1)
    q = p - 1

    !$omp do schedule(static)
    do i = 1, size(element)
      associate( face   => element(i) % face   &
               , edge   => element(i) % edge   &
               , vertex => element(i) % vertex )

        l = element(i) % id

        ! vertex data
        if (vertex(1)%rank > 0) u_vert(vertex(1)%rank, vertex(1)%id) = u(0,0,0,l)
        if (vertex(2)%rank > 0) u_vert(vertex(2)%rank, vertex(2)%id) = u(p,0,0,l)
        if (vertex(3)%rank > 0) u_vert(vertex(3)%rank, vertex(3)%id) = u(0,p,0,l)
        if (vertex(4)%rank > 0) u_vert(vertex(4)%rank, vertex(4)%id) = u(p,p,0,l)
        if (vertex(5)%rank > 0) u_vert(vertex(5)%rank, vertex(5)%id) = u(0,0,p,l)
        if (vertex(6)%rank > 0) u_vert(vertex(6)%rank, vertex(6)%id) = u(p,0,p,l)
        if (vertex(7)%rank > 0) u_vert(vertex(7)%rank, vertex(7)%id) = u(0,p,p,l)
        if (vertex(8)%rank > 0) u_vert(vertex(8)%rank, vertex(8)%id) = u(p,p,p,l)

        ! vertex valencies
        if (present(m_vert)) then
          do k = 1, 8
            if (vertex(k) % rank == 1)  m_vert(vertex(k) % id) = vertex(k) % val
          end do
        end if

        if (q < 1) cycle

        ! edge data
        if (edge( 1)%rank > 0) call edge( 1) % AlignWithMesh(u(1:q,0,0,l), u_edge(:, edge( 1)%rank, edge( 1)%id))
        if (edge( 2)%rank > 0) call edge( 2) % AlignWithMesh(u(1:q,p,0,l), u_edge(:, edge( 2)%rank, edge( 2)%id))
        if (edge( 3)%rank > 0) call edge( 3) % AlignWithMesh(u(1:q,0,p,l), u_edge(:, edge( 3)%rank, edge( 3)%id))
        if (edge( 4)%rank > 0) call edge( 4) % AlignWithMesh(u(1:q,p,p,l), u_edge(:, edge( 4)%rank, edge( 4)%id))
        if (edge( 5)%rank > 0) call edge( 5) % AlignWithMesh(u(0,1:q,0,l), u_edge(:, edge( 5)%rank, edge( 5)%id))
        if (edge( 6)%rank > 0) call edge( 6) % AlignWithMesh(u(p,1:q,0,l), u_edge(:, edge( 6)%rank, edge( 6)%id))
        if (edge( 7)%rank > 0) call edge( 7) % AlignWithMesh(u(0,1:q,p,l), u_edge(:, edge( 7)%rank, edge( 7)%id))
        if (edge( 8)%rank > 0) call edge( 8) % AlignWithMesh(u(p,1:q,p,l), u_edge(:, edge( 8)%rank, edge( 8)%id))
        if (edge( 9)%rank > 0) call edge( 9) % AlignWithMesh(u(0,0,1:q,l), u_edge(:, edge( 9)%rank, edge( 9)%id))
        if (edge(10)%rank > 0) call edge(10) % AlignWithMesh(u(p,0,1:q,l), u_edge(:, edge(10)%rank, edge(10)%id))
        if (edge(11)%rank > 0) call edge(11) % AlignWithMesh(u(0,p,1:q,l), u_edge(:, edge(11)%rank, edge(11)%id))
        if (edge(12)%rank > 0) call edge(12) % AlignWithMesh(u(p,p,1:q,l), u_edge(:, edge(12)%rank, edge(12)%id))

        ! edge valencies
        if (present(m_edge)) then
          do k = 1, 12
            if (edge(k) % rank == 1)  m_edge(edge(k) % id) = edge(k) % val
          end do
        end if

        ! face data
        if (face(1)%rank > 0) call face(1) % AlignWithMesh(u(0,1:q,1:q,l), u_face(:,:, face(1)%rank, face(1)%id))
        if (face(2)%rank > 0) call face(2) % AlignWithMesh(u(p,1:q,1:q,l), u_face(:,:, face(2)%rank, face(2)%id))
        if (face(3)%rank > 0) call face(3) % AlignWithMesh(u(1:q,0,1:q,l), u_face(:,:, face(3)%rank, face(3)%id))
        if (face(4)%rank > 0) call face(4) % AlignWithMesh(u(1:q,p,1:q,l), u_face(:,:, face(4)%rank, face(4)%id))
        if (face(5)%rank > 0) call face(5) % AlignWithMesh(u(1:q,1:q,0,l), u_face(:,:, face(5)%rank, face(5)%id))
        if (face(6)%rank > 0) call face(6) % AlignWithMesh(u(1:q,1:q,p,l), u_face(:,:, face(6)%rank, face(6)%id))

        ! face valencies
        if (present(m_face)) then
          do k = 1, 6
            if (face(k) % rank == 1)  m_face(face(k) % id) = face(k) % val
          end do
        end if

      end associate
    end do

  end subroutine ExtractTraces

  !-----------------------------------------------------------------------------
  !> Summation and, optionally averaging of the traces

  subroutine AssembleTraces(u_face, u_edge, u_vert, m_face, m_edge, m_vert, avg)

    real(RNP)        , intent(inout) :: u_face(:,:,:,:)
    real(RNP)        , intent(inout) :: u_edge(:,:,:)
    real(RNP)        , intent(inout) :: u_vert(:,:)
    integer          , intent(in)    :: m_face(:)
    integer          , intent(in)    :: m_edge(:)
    integer          , intent(in)    :: m_vert(:)
    logical, optional, intent(in)    :: avg

    integer :: i, k
    logical :: avg_

    if (present(avg)) then
      avg_ = avg
    else
      avg_ = .false.
    end if

    ! vertices .................................................................

    !$omp do schedule(static)
    do k = 1, size(m_vert)
      do i = 2, m_vert(k)
        u_vert(1,k) = u_vert(1,k) + u_vert(i,k)
      end do
      if (avg_) u_vert(1,k) = ONE/m_vert(k) * u_vert(1,k)
    end do
    !$omp end do nowait

    if (size(u_edge,1) == 0) return

    ! edges ....................................................................

    !$omp do schedule(static)
    do k = 1, size(m_edge)
      do i = 2, m_edge(k)
        u_edge(:,1,k) = u_edge(:,1,k) + u_edge(:,i,k)
      end do
      if (avg_) u_edge(:,1,k) = ONE/m_edge(k) * u_edge(:,1,k)
    end do
    !$omp end do nowait

    ! faces ....................................................................

    !$omp do schedule(static)
    do k = 1, size(m_face)
      do i = 2, m_face(k)
        u_face(:,:,1,k) = u_face(:,:,1,k) + u_face(:,:,i,k)
      end do
      if (avg_) u_face(:,:,1,k) = ONE/m_face(k) * u_face(:,:,1,k)
    end do

  end subroutine AssembleTraces

  !-----------------------------------------------------------------------------
  !> Distribution of skeleton data to elements

  subroutine DistributeSkeletonData(element, u, u_face, u_edge, u_vert)

    class(MeshElement_3D), intent(in) :: element(:)
    real(RNP), intent(inout) :: u(0:,0:,0:,:)
    real(RNP), intent(in)    :: u_face(:,:,:,:)
    real(RNP), intent(in)    :: u_edge(:,:,:)
    real(RNP), intent(in)    :: u_vert(:,:)

    integer :: i, l, p, q

    p = ubound(u,1)
    q = p - 1

    !$omp do schedule(static)
    do i = 1, size(element)
      associate( face   => element(i) % face   &
               , edge   => element(i) % edge   &
               , vertex => element(i) % vertex )

        l = element(i) % id

        ! vertex data
        u(0,0,0,l) = u_vert(1, vertex(1)%id)
        u(p,0,0,l) = u_vert(1, vertex(2)%id)
        u(0,p,0,l) = u_vert(1, vertex(3)%id)
        u(p,p,0,l) = u_vert(1, vertex(4)%id)
        u(0,0,p,l) = u_vert(1, vertex(5)%id)
        u(p,0,p,l) = u_vert(1, vertex(6)%id)
        u(0,p,p,l) = u_vert(1, vertex(7)%id)
        u(p,p,p,l) = u_vert(1, vertex(8)%id)

        if (q < 1) cycle

        ! edge data
        call edge( 1) % AlignFromMesh(u_edge(:,1, edge( 1)%id), u(1:q,0,0,l))
        call edge( 2) % AlignFromMesh(u_edge(:,1, edge( 2)%id), u(1:q,p,0,l))
        call edge( 3) % AlignFromMesh(u_edge(:,1, edge( 3)%id), u(1:q,0,p,l))
        call edge( 4) % AlignFromMesh(u_edge(:,1, edge( 4)%id), u(1:q,p,p,l))
        call edge( 5) % AlignFromMesh(u_edge(:,1, edge( 5)%id), u(0,1:q,0,l))
        call edge( 6) % AlignFromMesh(u_edge(:,1, edge( 6)%id), u(p,1:q,0,l))
        call edge( 7) % AlignFromMesh(u_edge(:,1, edge( 7)%id), u(0,1:q,p,l))
        call edge( 8) % AlignFromMesh(u_edge(:,1, edge( 8)%id), u(p,1:q,p,l))
        call edge( 9) % AlignFromMesh(u_edge(:,1, edge( 9)%id), u(0,0,1:q,l))
        call edge(10) % AlignFromMesh(u_edge(:,1, edge(10)%id), u(p,0,1:q,l))
        call edge(11) % AlignFromMesh(u_edge(:,1, edge(11)%id), u(0,p,1:q,l))
        call edge(12) % AlignFromMesh(u_edge(:,1, edge(12)%id), u(p,p,1:q,l))

        ! face data
        call face(1) % AlignFromMesh(u_face(:,:,1, face(1)%id), u(0,1:q,1:q,l))
        call face(2) % AlignFromMesh(u_face(:,:,1, face(2)%id), u(p,1:q,1:q,l))
        call face(3) % AlignFromMesh(u_face(:,:,1, face(3)%id), u(1:q,0,1:q,l))
        call face(4) % AlignFromMesh(u_face(:,:,1, face(4)%id), u(1:q,p,1:q,l))
        call face(5) % AlignFromMesh(u_face(:,:,1, face(5)%id), u(1:q,1:q,0,l))
        call face(6) % AlignFromMesh(u_face(:,:,1, face(6)%id), u(1:q,1:q,p,l))

      end associate
    end do

  end subroutine DistributeSkeletonData

  !=============================================================================

end submodule SM_0U

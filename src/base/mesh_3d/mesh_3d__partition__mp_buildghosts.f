!> summary:  Generation of 3d mesh ghost elements
!> author:   Joerg Stiller
!> date:     2021/02/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_3d__Partition) MP_BuildGhosts
  use Element_Transfer_Buffer_3d
  use Mesh_3d__Element_Indexing
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of ghost elements
  !>
  !> On entry, the mesh elements and mesh links must be complete. Using this
  !> information, the ghosts are created in `mesh % ghost(1:n_ghost)` and
  !> initialized as follows:
  !>
  !>   - `ghost % global_id` :
  !>      is the global ID of the corresponding mesh element
  !>
  !>   - `ghost % local_id` :
  !>      is the virtual element ID in the local mesh partition. It holds
  !>      `ghost(i) % local_id = mesh % n_elem + i`
  !>

  module subroutine BuildGhosts(mesh)
    class(Mesh3d_Partition), intent(inout) :: mesh !< local partition

    type(ElementTransferBuffer3d), asynchronous, allocatable, save :: global_id_buf
    type(ElementTransferBuffer3d), asynchronous, allocatable, save :: orientation_buf
    integer(IXL), allocatable, save :: global_id(:,:,:,:)
    integer(IXS), allocatable, save :: orientation(:,:,:,:)

    integer :: lf_id(-1:1,-1:1), gf_id(-1:1,-1:1)
    integer :: e, f, g, i, j, k, l, n

    !---------------------------------------------------------------------------
    ! initialization

    ! set up data structures ...................................................

    !$omp master
    allocate(mesh % ghost(mesh % n_ghost))
    !$omp end master
    if (mesh % n_ghost == 0) return

    ! these could be two OpenMP tasks
    !$omp master
    allocate(global_id( 1, 1, 1, mesh%n_elem + mesh%n_ghost ))
    global_id_buf = ElementTransferBuffer3d(mesh, global_id)
    allocate(orientation( 24, 1, 1, mesh%n_elem + mesh%n_ghost ))
    orientation_buf = ElementTransferBuffer3d(mesh, orientation)
    !$omp end master
    !$omp barrier

    ! transfer global ID and orientation .......................................

    !$omp do schedule(static)
    do l = 1, mesh % n_elem

      global_id(1,1,1,l) = mesh % element(l) % global_id

      do k = 1, 6
        orientation(2*k-1,1,1,l) = mesh % element(l) % face(k) % normal
        orientation(2*k  ,1,1,l) = mesh % element(l) % face(k) % rotation
      end do

      do k = 1, 12
        orientation(k+12,1,1,l) = mesh % element(l) % edge(k) % orientation
      end do

    end do

    call global_id_buf   % Transfer(mesh, global_id  , tag=100)
    call orientation_buf % Transfer(mesh, orientation, tag=200)

    call global_id_buf   % Merge(global_id)
    call orientation_buf % Merge(orientation)

    ! assign received data to ghosts ...........................................

    !$omp do schedule(static)
    do g = 1, mesh % n_ghost
      associate( face => mesh % ghost(g) % face &
               , edge => mesh % ghost(g) % edge )

        l = mesh % n_elem + g

        mesh % ghost(g) % global_id = global_id(1,1,1,l)
        mesh % ghost(g) % local_id  = l

        ! faces
        do k = 1, 6
          face(k) % normal   = orientation(2*k-1,1,1,l)
          face(k) % rotation = orientation(2*k  ,1,1,l)
        end do

        ! edges
        do k = 1, 12
          edge(k) % orientation = orientation(k+12,1,1,l)
        end do

      end associate
    end do

    !---------------------------------------------------------------------------
    ! copy face, edge and vertex IDs to elements to adjoining ghosts

    !$omp do schedule(static)
    do l = 1, mesh % n_elem
      associate( face     => mesh % element(l) % face     &
               , edge     => mesh % element(l) % edge     &
               , vertex   => mesh % element(l) % vertex   &
               , neighbor => mesh % element(l) % neighbor )

        ! faces, including associated edges and vertices .......................

        do k = 1, 6
          if (face(k) % primary < 1) cycle

          n = face(k) % n_neighbor
          i = face(k) % i_neighbor
          do j = i, i+n-1

            ! check if neighbor(j) is a ghost
            g = neighbor(j) % id - mesh % n_elem
            if (g > 0) then

              ! identify corresponding ghost face
              f = neighbor(j) % cc

              ! get local element face mesh IDs at position ξ₁=i and ξ₂=j
              lf_id(-1,-1) = vertex( V_FACE(1,k) ) % id
              lf_id( 0,-1) = edge  ( E_FACE(1,k) ) % id
              lf_id( 1,-1) = vertex( V_FACE(2,k) ) % id
              lf_id(-1, 0) = edge  ( E_FACE(3,k) ) % id
              lf_id( 0, 0) = face  (          k  ) % id
              lf_id( 1, 0) = edge  ( E_FACE(4,k) ) % id
              lf_id(-1, 1) = vertex( V_FACE(3,k) ) % id
              lf_id( 0, 1) = edge  ( E_FACE(2,k) ) % id
              lf_id( 1, 1) = vertex( V_FACE(4,k) ) % id

              ! transform ID array to ghost face orientation
              if (mesh % structured) then
                gf_id = lf_id
              else if (face(k)%normal == 1 .and. face(k)%rotation == 0) then
                ! local element face is aligned with mesh face
                call mesh%ghost(g)%face(f)%AlignWithElementFace(lf_id, gf_id)
              else
                ! ghost face is aligned with mesh face
                call face(k) % AlignWithMeshFace(lf_id, gf_id)
              end if

              ! set ghost face mesh IDs
              mesh % ghost(g) % vertex( V_FACE(1,f) ) % id = gf_id(-1,-1)
              mesh % ghost(g) % edge  ( E_FACE(1,f) ) % id = gf_id( 0,-1)
              mesh % ghost(g) % vertex( V_FACE(2,f) ) % id = gf_id( 1,-1)
              mesh % ghost(g) % edge  ( E_FACE(3,f) ) % id = gf_id(-1, 0)
              mesh % ghost(g) % face  (          f  ) % id = gf_id( 0, 0)
              mesh % ghost(g) % edge  ( E_FACE(4,f) ) % id = gf_id( 1, 0)
              mesh % ghost(g) % vertex( V_FACE(3,f) ) % id = gf_id(-1, 1)
              mesh % ghost(g) % edge  ( E_FACE(2,f) ) % id = gf_id( 0, 1)
              mesh % ghost(g) % vertex( V_FACE(4,f) ) % id = gf_id( 1, 1)

            end if
          end do
        end do

        ! edges, including associated vertices .................................

        do k = 1, 12
          if (edge(k) % primary < 1) cycle
          n = edge(k) % n_neighbor
          i = edge(k) % i_neighbor
          do j = i, i+n-1

            ! check if neighbor(j) is a ghost
            g = neighbor(j) % id - mesh % n_elem
            if (g > 0) then
              associate( ghost_edge   => mesh % ghost(g) % edge   &
                       , ghost_vertex => mesh % ghost(g) % vertex )

                ! identify corresponding ghost edge
                e = neighbor(j) % cc

                ! copy edge and vertex IDs
                if ( mesh % structured    .or.                              &
                     edge(k) % orientation == ghost_edge(e) % orientation ) &
                then
                  ghost_edge  (          e  ) % id = edge  (          k  ) % id
                  ghost_vertex( V_EDGE(1,e) ) % id = vertex( V_EDGE(1,k) ) % id
                  ghost_vertex( V_EDGE(2,e) ) % id = vertex( V_EDGE(2,k) ) % id
                else
                  ghost_edge  (          e  ) % id = edge  (          k  ) % id
                  ghost_vertex( V_EDGE(1,e) ) % id = vertex( V_EDGE(2,k) ) % id
                  ghost_vertex( V_EDGE(2,e) ) % id = vertex( V_EDGE(1,k) ) % id
                end if

              end associate
            end if
          end do
        end do

        ! vertices .............................................................

        do k = 1, 8
          if (vertex(k) % primary < 1) cycle
          n = vertex(k) % n_neighbor
          i = vertex(k) % i_neighbor
          do j = i, i+n-1
            g = neighbor(j) % id - mesh % n_elem
            if (g > 0) then
              mesh % ghost(g) % vertex(neighbor(j) % cc) % id = vertex(k) % id
            end if
          end do
        end do

      end associate
    end do

    !---------------------------------------------------------------------------
    ! finalization

    !$omp master
    deallocate(orientation, orientation_buf)
    deallocate(global_id  , global_id_buf  )
    !$omp end master

  end subroutine BuildGhosts

  !=============================================================================

end submodule MP_BuildGhosts

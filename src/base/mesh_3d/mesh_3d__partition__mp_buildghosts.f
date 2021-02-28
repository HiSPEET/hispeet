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
  !>   - `ghost % face % {id, rank, val}` :
  !>      refer to
  !>        * the corresponding mesh face,
  !>        * the rank of the ghost element face and
  !>        * the valency of the face
  !>
  !>   - `ghost % edge % {id, rank, val}` :
  !>      refer to
  !>        * the corresponding mesh edge,
  !>        * the rank of the ghost element edge and
  !>        * the valency of the edge
  !>
  !>   - `ghost % vertex % {id, rank, val}` :
  !>      refer to
  !>        * the corresponding mesh vertex,
  !>        * the rank of the ghost element vertex and
  !>        * the valency of the vertex

  module subroutine BuildGhosts(mesh)
    class(Mesh3d_Partition), intent(inout) :: mesh !< local partition

    type(ElementTransferBuffer3d), asynchronous, allocatable, save :: global_id_buf
    type(ElementTransferBuffer3d), asynchronous, allocatable, save :: orientation_buf
    integer(IXL), allocatable, save :: global_id(:,:,:,:)
    integer(IXS), allocatable, save :: orientation(:,:,:,:)

    integer      :: lfi(-1:1,-1:1), gfi(-1:1,-1:1), lei(-1:1)
    integer(IXS) :: lfr(-1:1,-1:1), gfr(-1:1,-1:1), ler(-1:1)
    integer(IXS) :: lfv(-1:1,-1:1), gfv(-1:1,-1:1), lev(-1:1)
    integer      :: e, f, g, i, j, k, l, m, n, v

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
      associate(element => mesh % element(l))

        ! faces, including associated edges and vertices .......................

        do k = 1, 6
          if (element % face(k) % rank /= 1) cycle

          n = element % face(k) % n_neighbor
          i = element % face(k) % i_neighbor
          do j = i, i+n-1

            ! check if neighbor(j) is a ghost
            m = element % neighbor(j) % id
            g = m - mesh % n_elem
            if (g > 0) then
              associate(ghost => mesh % ghost(g))

                ! identify corresponding ghost face
                f = element % neighbor(j) % cc

                ! local element face IDs at position ξ₁=i and ξ₂=j
                lfi(-1,-1) = element % vertex( V_FACE(1,k) ) % id
                lfi( 0,-1) = element % edge  ( E_FACE(1,k) ) % id
                lfi( 1,-1) = element % vertex( V_FACE(2,k) ) % id
                lfi(-1, 0) = element % edge  ( E_FACE(3,k) ) % id
                lfi( 0, 0) = element % face  (          k  ) % id
                lfi( 1, 0) = element % edge  ( E_FACE(4,k) ) % id
                lfi(-1, 1) = element % vertex( V_FACE(3,k) ) % id
                lfi( 0, 1) = element % edge  ( E_FACE(2,k) ) % id
                lfi( 1, 1) = element % vertex( V_FACE(4,k) ) % id

                ! ranks of ghost components based on virtual element ID
                call element % DetermineVertexRank(V_FACE(1,k), m, lfr(-1,-1))
                call element % DetermineEdgeRank  (E_FACE(1,k), m, lfr( 0,-1))
                call element % DetermineVertexRank(V_FACE(2,k), m, lfr( 1,-1))
                call element % DetermineEdgeRank  (E_FACE(3,k), m, lfr(-1, 0))
                call element % DetermineFaceRank  (         k , m, lfr( 0, 0))
                call element % DetermineEdgeRank  (E_FACE(4,k), m, lfr( 1, 0))
                call element % DetermineVertexRank(V_FACE(3,k), m, lfr(-1, 1))
                call element % DetermineEdgeRank  (E_FACE(2,k), m, lfr( 0, 1))
                call element % DetermineVertexRank(V_FACE(4,k), m, lfr( 1, 1))

                ! local element face valencies at position ξ₁=i and ξ₂=j
                lfv(-1,-1) = element % vertex( V_FACE(1,k) ) % val
                lfv( 0,-1) = element % edge  ( E_FACE(1,k) ) % val
                lfv( 1,-1) = element % vertex( V_FACE(2,k) ) % val
                lfv(-1, 0) = element % edge  ( E_FACE(3,k) ) % val
                lfv( 0, 0) = element % face  (          k  ) % val
                lfv( 1, 0) = element % edge  ( E_FACE(4,k) ) % val
                lfv(-1, 1) = element % vertex( V_FACE(3,k) ) % val
                lfv( 0, 1) = element % edge  ( E_FACE(2,k) ) % val
                lfv( 1, 1) = element % vertex( V_FACE(4,k) ) % val

                ! transform ID and rank arrays to ghost face orientation
                if (mesh % structured) then
                  gfi = lfi
                  gfr = lfr
                  gfv = lfv
                else if ( element % face(k) % normal   == 1 .and. &
                          element % face(k) % rotation == 0 ) then
                  ! local element face is aligned with mesh face
                  call ghost % face(f) % AlignWithElementFace(lfi, gfi)
                  call ghost % face(f) % AlignWithElementFace(lfr, gfr)
                  call ghost % face(f) % AlignWithElementFace(lfv, gfv)
                else
                  ! ghost face is aligned with mesh face
                  call element % face(k) % AlignWithMeshFace(lfi, gfi)
                  call element % face(k) % AlignWithMeshFace(lfr, gfr)
                  call element % face(k) % AlignWithMeshFace(lfv, gfv)
                end if

                ! set ghost face IDs and ranks
                ghost % vertex( V_FACE(1,f) ) % id   = gfi(-1,-1)
                ghost % vertex( V_FACE(1,f) ) % rank = gfr(-1,-1)
                ghost % vertex( V_FACE(1,f) ) % val  = gfv(-1,-1)
                ghost % edge  ( E_FACE(1,f) ) % id   = gfi( 0,-1)
                ghost % edge  ( E_FACE(1,f) ) % rank = gfr( 0,-1)
                ghost % edge  ( E_FACE(1,f) ) % val  = gfv( 0,-1)
                ghost % vertex( V_FACE(2,f) ) % id   = gfi( 1,-1)
                ghost % vertex( V_FACE(2,f) ) % rank = gfr( 1,-1)
                ghost % vertex( V_FACE(2,f) ) % val  = gfv( 1,-1)
                ghost % edge  ( E_FACE(3,f) ) % id   = gfi(-1, 0)
                ghost % edge  ( E_FACE(3,f) ) % rank = gfr(-1, 0)
                ghost % edge  ( E_FACE(3,f) ) % val  = gfv(-1, 0)
                ghost % face  (          f  ) % id   = gfi( 0, 0)
                ghost % face  (          f  ) % rank = gfr( 0, 0)
                ghost % face  (          f  ) % val  = gfv( 0, 0)
                ghost % edge  ( E_FACE(4,f) ) % id   = gfi( 1, 0)
                ghost % edge  ( E_FACE(4,f) ) % rank = gfr( 1, 0)
                ghost % edge  ( E_FACE(4,f) ) % val  = gfv( 1, 0)
                ghost % vertex( V_FACE(3,f) ) % id   = gfi(-1, 1)
                ghost % vertex( V_FACE(3,f) ) % rank = gfr(-1, 1)
                ghost % vertex( V_FACE(3,f) ) % val  = gfv(-1, 1)
                ghost % edge  ( E_FACE(2,f) ) % id   = gfi( 0, 1)
                ghost % edge  ( E_FACE(2,f) ) % rank = gfr( 0, 1)
                ghost % edge  ( E_FACE(2,f) ) % val  = gfv( 0, 1)
                ghost % vertex( V_FACE(4,f) ) % id   = gfi( 1, 1)
                ghost % vertex( V_FACE(4,f) ) % rank = gfr( 1, 1)
                ghost % vertex( V_FACE(4,f) ) % val  = gfv( 1, 1)

              end associate
            end if
          end do
        end do

        ! edges, including associated vertices .................................

        do k = 1, 12
          if (element % edge(k) % rank /= 1) cycle
          n = element % edge(k) % n_neighbor
          i = element % edge(k) % i_neighbor
          do j = i, i+n-1

            ! check if neighbor(j) is a ghost
            m = element % neighbor(j) % id
            g = m - mesh % n_elem
            if (g > 0) then
              associate(ghost => mesh % ghost(g))

                ! identify corresponding ghost edge
                e = element % neighbor(j) % cc

                ! get element edge mesh IDs at position ξ=i
                lei(-1) = element % vertex( V_EDGE(1,k) ) % id
                lei( 0) = element % edge  (          k  ) % id
                lei( 1) = element % vertex( V_EDGE(2,k) ) % id

                ! ranks of the corresponding ghost components at position ξ=i
                call element % DetermineVertexRank( V_EDGE(1,k), m, ler(-1))
                call element % DetermineEdgeRank  (          k , m, ler( 0))
                call element % DetermineVertexRank( V_EDGE(2,k), m, ler( 1))

                ! get element edge mesh valencies at position ξ=i
                lev(-1) = element % vertex( V_EDGE(1,k) ) % val
                lev( 0) = element % edge  (          k  ) % val
                lev( 1) = element % vertex( V_EDGE(2,k) ) % val

                ! copy IDs and ranks
                if ( mesh % structured    .or.                                  &
                     element%edge(k)%orientation == ghost%edge(e)%orientation ) &
                then
                  ghost % vertex( V_EDGE(1,e) ) % id   = lei(-1)
                  ghost % vertex( V_EDGE(1,e) ) % rank = ler(-1)
                  ghost % vertex( V_EDGE(1,e) ) % val  = lev(-1)
                  ghost % edge  (          e  ) % id   = lei( 0)
                  ghost % edge  (          e  ) % rank = ler( 0)
                  ghost % edge  (          e  ) % val  = lev( 0)
                  ghost % vertex( V_EDGE(2,e) ) % id   = lei( 1)
                  ghost % vertex( V_EDGE(2,e) ) % rank = ler( 1)
                  ghost % vertex( V_EDGE(2,e) ) % val  = lev( 1)
                else
                  ghost % vertex( V_EDGE(1,e) ) % id   = lei( 1)
                  ghost % vertex( V_EDGE(1,e) ) % rank = ler( 1)
                  ghost % vertex( V_EDGE(1,e) ) % val  = lev( 1)
                  ghost % edge  (          e  ) % id   = lei( 0)
                  ghost % edge  (          e  ) % rank = ler( 0)
                  ghost % edge  (          e  ) % val  = lev( 0)
                  ghost % vertex( V_EDGE(2,e) ) % id   = lei(-1)
                  ghost % vertex( V_EDGE(2,e) ) % rank = ler(-1)
                  ghost % vertex( V_EDGE(2,e) ) % val  = lev(-1)
                end if

              end associate
            end if
          end do
        end do

        ! vertices .............................................................

        do k = 1, 8
          if (element % vertex(k) % rank /= 1) cycle
          n = element % vertex(k) % n_neighbor
          i = element % vertex(k) % i_neighbor
          do j = i, i+n-1
            m = element % neighbor(j) % id
            g = m - mesh % n_elem
            if (g > 0) then
              associate(ghost => mesh % ghost(g))
                v = element % neighbor(j) % cc
                ghost % vertex(v) % id  = element % vertex(k) % id
                ghost % vertex(v) % val = element % vertex(k) % val
                call element % DetermineVertexRank(k, m, ghost%vertex(v)%rank)
              end associate
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

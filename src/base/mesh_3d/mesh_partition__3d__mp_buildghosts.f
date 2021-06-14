!> summary:  Generation of 3d mesh ghost elements
!> author:   Joerg Stiller
!> date:     2021/02/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_Partition__3D) MP_BuildGhosts
  use Element_Transfer_Buffer_3d
  use Mesh_Element_Indexing__3D
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
  !>   - `ghost % face % id` :
  !>      is the corresponding mesh face,
  !>
  !>   - `ghost % edge % id` :
  !>      is the corresponding mesh edge,
  !>
  !>   - `ghost % vertex % id` :
  !>      is the corresponding mesh vertex,

  module subroutine BuildGhosts(mesh)
    class(MeshPartition_3D), intent(inout) :: mesh !< local partition

    type(ElementTransferBuffer3d), asynchronous, allocatable :: global_id_buf
    type(ElementTransferBuffer3d), asynchronous, allocatable :: orientation_buf
    integer(IXL), allocatable :: global_id(:,:,:,:)
    integer(IXS), allocatable :: orientation(:,:,:,:)

    integer :: lfi(-1:1,-1:1), gfi(-1:1,-1:1), lei(-1:1)
    integer :: e, f, g, i, j, k, l, m, n, v

    !---------------------------------------------------------------------------
    ! initialization

    ! set up data structures ...................................................

    allocate(mesh % ghost(mesh % n_ghost))
    if (mesh % n_ghost == 0) return

    ! these could be two OpenMP tasks
    allocate(global_id( 1, 1, 1, mesh%n_elem + mesh%n_ghost ))
    global_id_buf = ElementTransferBuffer3d(mesh, global_id)
    allocate(orientation( 24, 1, 1, mesh%n_elem + mesh%n_ghost ))
    orientation_buf = ElementTransferBuffer3d(mesh, orientation)

    ! transfer global ID and orientation .......................................

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

    do l = 1, mesh % n_elem
      associate(element => mesh % element(l))

        ! faces, including associated edges and vertices .......................

        do k = 1, 6
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

                ! transform ID arrays to ghost face orientation
                if (mesh % structured) then
                  gfi = lfi
                else if ( element % face(k) % normal   == 1 .and. &
                          element % face(k) % rotation == 0 ) then
                  ! local element face is aligned with mesh face
                  call ghost % face(f) % AlignWithELement(lfi, gfi)
                else
                  ! ghost face is aligned with mesh face
                  call element % face(k) % AlignWithMesh(lfi, gfi)
                end if

                ! set ghost face IDs and ranks
                ghost % vertex( V_FACE(1,f) ) % id   = gfi(-1,-1)
                ghost % edge  ( E_FACE(1,f) ) % id   = gfi( 0,-1)
                ghost % vertex( V_FACE(2,f) ) % id   = gfi( 1,-1)
                ghost % edge  ( E_FACE(3,f) ) % id   = gfi(-1, 0)
                ghost % face  (          f  ) % id   = gfi( 0, 0)
                ghost % edge  ( E_FACE(4,f) ) % id   = gfi( 1, 0)
                ghost % vertex( V_FACE(3,f) ) % id   = gfi(-1, 1)
                ghost % edge  ( E_FACE(2,f) ) % id   = gfi( 0, 1)
                ghost % vertex( V_FACE(4,f) ) % id   = gfi( 1, 1)

              end associate
            end if
          end do
        end do

        ! edges, including associated vertices .................................

        do k = 1, 12
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

                ! copy IDs and ranks
                if ( mesh % structured    .or.                                  &
                     element%edge(k)%orientation == ghost%edge(e)%orientation ) &
                then
                  ghost % vertex( V_EDGE(1,e) ) % id   = lei(-1)
                  ghost % edge  (          e  ) % id   = lei( 0)
                  ghost % vertex( V_EDGE(2,e) ) % id   = lei( 1)
                else
                  ghost % vertex( V_EDGE(1,e) ) % id   = lei( 1)
                  ghost % edge  (          e  ) % id   = lei( 0)
                  ghost % vertex( V_EDGE(2,e) ) % id   = lei(-1)
                end if

              end associate
            end if
          end do
        end do

        ! vertices .............................................................

        do k = 1, 8
          n = element % vertex(k) % n_neighbor
          i = element % vertex(k) % i_neighbor
          do j = i, i+n-1
            m = element % neighbor(j) % id
            g = m - mesh % n_elem
            if (g > 0) then
              associate(ghost => mesh % ghost(g))
                v = element % neighbor(j) % cc
                ghost % vertex(v) % id  = element % vertex(k) % id
              end associate
            end if
          end do
        end do

      end associate
    end do

  end subroutine BuildGhosts

  !=============================================================================

end submodule MP_BuildGhosts

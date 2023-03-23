!> summary:  Generation of 3d mesh ghost elements
!> author:   Joerg Stiller
!> date:     2021/02/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_BuildGhosts
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of ghost elements
  !>
  !> On entry, the mesh elements and mesh links must be complete. Using this
  !> information, the ghosts are created in `mesh % ghost(1:n_ghost)` and
  !> initialized as follows:
  !>
  !>   - `ghost % id` :
  !>      is the virtual element ID in the local mesh partition. It holds
  !>      `ghost(i) % id = mesh % n_elem + i`
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
    class(Mesh_3D), intent(inout) :: mesh !< local partition

    integer :: lfi(-1:1,-1:1), gfi(-1:1,-1:1), lei(-1:1)
    integer :: e, f, g, i, j, k, l, m, n, v

    !---------------------------------------------------------------------------
    ! initialization

    ! set up data structures ...................................................

    allocate(mesh % ghost(mesh % n_ghost))
    if (mesh % n_ghost == 0) return

    !---------------------------------------------------------------------------
    ! ghost IDs

    do l = 1, mesh % n_ghost
      mesh % ghost(l) % id = mesh % n_elem + l
    end do

    !---------------------------------------------------------------------------
    ! copy face, edge and vertex IDs from elements to adjoining ghosts

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
                f = ElementFaceID(element % neighbor(j) % component)

                ! local element-face component IDs at position ξ₁=i and ξ₂=j
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
                else
                  call element % AlignWithNeighborFace(k, j, lfi, gfi)
                end if

                ! set ghost face IDs and ranks
                ghost % vertex( V_FACE(1,f) ) % id = gfi(-1,-1)
                ghost % edge  ( E_FACE(1,f) ) % id = gfi( 0,-1)
                ghost % vertex( V_FACE(2,f) ) % id = gfi( 1,-1)
                ghost % edge  ( E_FACE(3,f) ) % id = gfi(-1, 0)
                ghost % face  (          f  ) % id = gfi( 0, 0)
                ghost % edge  ( E_FACE(4,f) ) % id = gfi( 1, 0)
                ghost % vertex( V_FACE(3,f) ) % id = gfi(-1, 1)
                ghost % edge  ( E_FACE(2,f) ) % id = gfi( 0, 1)
                ghost % vertex( V_FACE(4,f) ) % id = gfi( 1, 1)

                ! ghost face orientation, exploiting that the local element
                ! is aligned with the mesh face by construction
                call ghost % face(f) % SetOrientation(           &
                       efv = ghost   % vertex(V_FACE(:,f)) % id, &
                       mfv = element % vertex(V_FACE(:,k)) % id  )

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
                e = ElementEdgeID(element % neighbor(j) % component)

                ! get element edge mesh IDs at position ξ=i
                lei(-1) = element % vertex( V_EDGE(1,k) ) % id
                lei( 0) = element % edge  (          k  ) % id
                lei( 1) = element % vertex( V_EDGE(2,k) ) % id

                ! copy IDs and assign edge orientation
                if (mesh % structured .or. element % NeigborEdgeIsAligned(k,j)) &
                then
                  ghost % vertex( V_EDGE(1,e) ) % id = lei(-1)
                  ghost % vertex( V_EDGE(2,e) ) % id = lei( 1)
                  ghost % edge(e) % id          =  lei( 0)
                  ghost % edge(e) % orientation =  element % edge(k) % orientation
                else
                  ghost % vertex( V_EDGE(1,e) ) % id = lei( 1)
                  ghost % vertex( V_EDGE(2,e) ) % id = lei(-1)
                  ghost % edge(e) % id          =  lei( 0)
                  ghost % edge(e) % orientation = -element % edge(k) % orientation
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
                v = ElementVertexID(element % neighbor(j) % component)
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

!> summary:  Generation of 3d mesh ghost elements
!> author:   Joerg Stiller
!> date:     2021/02/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_3d__Partition) MP_BuildGhosts
  use Element_Transfer_Buffer_3d
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
    integer :: e, g, i, j, k, n

    ! initialization ...........................................................

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

    ! extract and adopt IDs ....................................................

    !$omp do schedule(static)
    do e = 1, mesh % n_elem
      associate( face     => mesh % element(e) % face     &
               , edge     => mesh % element(e) % edge     &
               , vertex   => mesh % element(e) % vertex   &
               , neighbor => mesh % element(e) % neighbor )

        ! global element ID
        global_id(1,1,1,e) = mesh % element(e) % global_id

        ! faces
        do k = 1, 6
          orientation(2*k-1,1,1,e) = face(k) % normal
          orientation(2*k  ,1,1,e) = face(k) % rotation
          if (face(k) % primary /= 1) cycle
          n = face(k) % n_neighbor
          i = face(k) % i_neighbor
          do j = i, i+n-1
            g = neighbor(j) % id - mesh % n_elem
            if (g > 0) then
              mesh % ghost(g) % face(neighbor(j) % cc) % id = face(k) % id
            end if
          end do
        end do

        ! edges
        do k = 1, 12
          orientation(k+12,1,1,e) = edge(k) % orientation
          if (edge(k) % primary /= 1) cycle
          n = edge(k) % n_neighbor
          i = edge(k) % i_neighbor
          do j = i, i+n-1
            g = neighbor(j) % id - mesh % n_elem
            if (g > 0) then
              mesh % ghost(g) % edge(neighbor(j) % cc) % id = edge(k) % id
            end if
          end do
        end do

        ! vertices
        do k = 1, 8
          if (vertex(k) % primary /= 1) cycle
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

    ! transfer & merge ghost data ..............................................

    call global_id_buf   % Transfer(mesh, global_id  , tag=100)
    call orientation_buf % Transfer(mesh, orientation, tag=200)

    call global_id_buf   % Merge(global_id)
    call orientation_buf % Merge(orientation)

    ! assign received data to ghosts ...........................................

    !$omp do schedule(static)
    do g = 1, mesh % n_ghost
      associate( face => mesh % ghost(g) % face &
               , edge => mesh % ghost(g) % edge )

        e = mesh % n_elem + g

        mesh % ghost(g) % global_id = global_id(1,1,1,e)

        ! faces
        do k = 1, 6
          face(k) % normal   = orientation(2*k-1,1,1,e)
          face(k) % rotation = orientation(2*k  ,1,1,e)
        end do

        ! edges
        do k = 1, 12
          edge(k) % orientation = orientation(k+12,1,1,e)
        end do

      end associate
    end do

    ! finalization .............................................................

    !$omp master
    deallocate(orientation, orientation_buf)
    deallocate(global_id  , global_id_buf  )
    !$omp end master

  end subroutine BuildGhosts

  !=============================================================================

end submodule MP_BuildGhosts

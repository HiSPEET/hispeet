!> summary:  Merging of periodic vertices
!> author:   Benedikt Wex, Joerg Stiller
!> date:     2025/03/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> The identification of periodic vertices is based on the following assumptions
!>   - coupled boundaries possess an identical subdivision in to surfaces
!>   - the surfaces are grouped in to master-slave pairs
!>   - surface nodes are set to related mesh vertex or to -1 if there is none
!>   - every master surface provides an affine transformation to the slave
!>
!> Periodic vertices are identified and merged as follows
!>   1. Select a master-slave surface pair
!>   2. Build a KD tree based on the coordinates of the slave vertices
!>   3. Use the affine transform to project a given master vertex to the
!>      slave surface
!>   4. Seach the KD tree to identify the closest slave vertex
!>   5. Merge the IDs (tags) of master and slave based on the following rules:
!>        - if the master has not been merged yet, adopt the slave ID
!>        - if only the master has been merged, adopt the master ID
!>        - if both have been merged with different IDs, master and all other
!>          nodes sharing its ID adopt the slave ID
!>        - if both have been merged with an identical ID, nothing is left to do
!>
!===============================================================================

submodule(Import_GMSH__3D) MP_MergePeriodicVertices
  use KD_Tree
  implicit none

  !-----------------------------------------------------------------------------
  !> number of neighbors to return from KD tree

  integer, parameter :: nn = 2

contains

  !-----------------------------------------------------------------------------
  !> Merges the IDs of periodic vertices

  module subroutine MergePeriodicVertices(mesh, surface)
    class(GenericMesh_3D), intent(inout) :: mesh       !< generic mesh
    type(MshSurface),      intent(in)    :: surface(:) !< GMSH surfaces

    type(KDTree), pointer :: tree        ! KD tree
    type(KDTree_Result)   :: results(nn) ! nearest neighbors found

    real(KDTree_RP), allocatable :: leaf_pos(:,:) ! KD tree leaf position
    integer,         allocatable :: leaf_id(:)    ! KD tree leaf ID
    logical,         allocatable :: merged(:)     ! indicator of merged node

    real(RNP) :: x_master(3)  ! node coordinates on master surface
    real(RNP) :: x_slave (3)  ! node coordinates on slave surface
    real(RNP) :: rel_dis, max_rel_dis = 0.1
    real(RNP) :: max_dis, mean_dis

    integer :: i, j, k, l, m, s
    integer :: id_master, id_slave
    integer :: n_leaf, n_match

    associate(vertex => mesh % vertex)

      ! initialization .........................................................

      allocate(merged(size(vertex)), source = .false.)

      n_match = 0
      max_dis = 0

      ! merge periodic node IDs ................................................

      do m = 1, size(surface)
        if (.not. surface(m)%master) cycle

        ! ID of coupled "slave" surface
        s = surface(m) % periodicID

        ! coords and tags on slave surface (kd-tree leaves)
        n_leaf = surface(s) % nNodes
        allocate(leaf_pos(3,n_leaf), leaf_id(n_leaf))
        do k = 1, n_leaf
          leaf_pos(:,k) = vertex(surface(s) % nodes(k)) % x
          leaf_id (  k) = vertex(surface(s) % nodes(k)) % id
        end do

        ! build kd-tree
        tree => KDTree_Create(leaf_pos, sort=.false., rearrange=.false.)

        ! find periodic node on slave surface s for node k of surface m
        do k = 1, n_leaf

          ! coordinates of node k on master surface
          x_master = vertex(surface(m) % nodes(k)) % x

          ! apply affinity transform to obtain coordinates on slave surface
          x_slave = matmul(surface(m)%A(1:3,:), [x_master, ONE])

          ! seach kd-tree for nn nearest neighbor points of projection of master
          call KDTree_N_Nearest(tree, qv = x_slave, nn = nn, results = results)
          call KDTree_SortResults(nn, results)

          ! check security margin of nearest neigbor and other points
          rel_dis = results(1)%dis / max(results(2)%dis, epsilon(ONE))
          if (rel_dis < max_rel_dis) then
            n_match  = n_match + 1
            mean_dis = mean_dis + rel_dis
            max_dis  = max(max_dis, rel_dis)
          else
            write(*,'(2X,A)') 'distance 1:', results(1)%dis
            write(*,'(2X,A)') 'distance 2:', results(2)%dis
            write(*,'(2X,A)') 'difference:', abs(results(1)%dis - results(2)%dis)
            write(*,'(2X,A)') 'ratio:     ', rel_dis
            call Error( 'MergePeriodicVertices'                                  &
                      , 'Failed to distinguish matching points on slave surface' &
                      , 'Import_GMSH__3D'                                        )
          end if

          id_master = vertex(surface(m) % nodes(k)) % id
          id_slave  = leaf_id(results(1) % idx)

          MERGING: if (merged(id_master)) then
            ! master has already been merged

            if (merged(id_slave)) then
              ! slave has already been merged
              if (id_master == id_slave) then
                ! a) master and slave merged with identical ID
                cycle
              else
                ! b) change master and related nodes to slave ID
                do l = 1, size(surface)
                do i = 1, surface(l)%nNodes
                  if (vertex(surface(l)%nodes(i)) % id == id_master) then
                    vertex(surface(l)%nodes(i)) % id = id_slave
                  end if
                end do
                end do
              end if
            else
              ! c) slave has not been merged yet and adopts master ID
              do l = 1, n_leaf
                if (vertex(surface(s)%nodes(l)) % id == id_slave) then
                  vertex(surface(s)%nodes(l)) % id = id_master
                end if
              end do
            end if

          else
            ! d) master not merged yet
            vertex(surface(m)%nodes(k)) % id = id_slave
            merged(id_master) = .true.
            merged(id_slave ) = .true.

          end if MERGING

        end do

        deallocate(leaf_pos, leaf_id)
        call KDTree_Destroy(tree)

      end do

      write(*,'(4X,A,G0)') 'number of matching pairs:  ', n_match
      write(*,'(4X,A,G0)') 'mean    relative distance: ', mean_dis / n_match
      write(*,'(4X,A,G0)') 'maximum relative distance: ', max_dis

    end associate

  end subroutine MergePeriodicVertices

  !=============================================================================

end submodule MP_MergePeriodicVertices

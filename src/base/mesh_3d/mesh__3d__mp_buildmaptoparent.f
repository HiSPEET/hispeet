!> summary:  Generation of map to parent elements
!> author:   Joerg Stiller
!> date:     2023/04/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_BuildMapToParent
  use Quick_Sort
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of map to parent elements
  !>
  !> The routine provides
  !>
  !>   - `mesh % n_parent`
  !>   - `mesh % map_parent(:)`
  !>
  !> and requires
  !>
  !>   -  `mesh % element % cluster_id`
  !>   -  `mesh % element % frozen`
  !>   -  `mesh % element % adaptation % parent_proc`
  !>   -  `mesh % element % adaptation % parent_id`
  !>
  !> The procedure is capable to cope with orphaned elements whose parents have
  !> been removed in course of an ongoing adaptation process.

  module subroutine BuildMapToParent(mesh)
    class(Mesh_3D), intent(inout) :: mesh  !< mesh parition

    integer, allocatable :: n_cluster(:)
    integer, allocatable :: n_active(:)
    integer, allocatable :: n_frozen(:)
    integer, allocatable :: parent_proc(:)

    type ParentRanking
      integer, allocatable :: id_active(:) ! active element parent ID
      integer, allocatable :: rk_active(:) ! active element parent rank
      integer, allocatable :: id_frozen(:) ! frozen element parent ID
      integer, allocatable :: rk_frozen(:) ! frozen element parent rank
    end type ParentRanking
    type(ParentRanking), allocatable :: parent(:)

    integer :: n_proc, cluster_id
    integer :: e, i, p, p_min, p_max

    !$omp master

    if (allocated(mesh % map_parent)) then
      deallocate(mesh % map_parent)
    end if

    if (mesh % is_root .or. mesh % n_elem < 1) then

      ! root or empty mesh .....................................................

      mesh % n_parent = 0
      allocate(mesh % map_parent(0))

    else

      ! initialization .........................................................

      call MPI_Comm_size(mesh % comm_world, n_proc)

      allocate( n_cluster   (0:n_proc-1))
      allocate( n_active    (0:n_proc-1), source =  0)
      allocate( n_frozen    (0:n_proc-1), source =  0)
      allocate( parent_proc (0:n_proc-1), source = -1)

      p_min = n_proc-1
      p_max = 0

      ! counts .................................................................

      cluster_id = 0
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))
          if (element % cluster_id == cluster_id)     cycle ! skip siblings
          if (element % adaptation % parent_proc < 0) cycle ! skip orphans
          cluster_id = element % cluster_id
          p = element % adaptation % parent_proc
          p_min = min(p_min, p)
          p_max = min(p_max, p)
          if (mesh % element(e) % frozen) then
            n_frozen(p) = n_frozen(p) + 1
          else
            n_active(p) = n_active(p) + 1
          end if
        end associate
      end do

      ! setup of maps ..........................................................

      i = 0
      do p = p_min, p_max
        n_cluster(p) = n_active(p) + n_frozen(p)
        if (n_cluster(p) > 0) then
          i = i + 1
          parent_proc(p) = i
        end if
      end do
      mesh % n_parent = i
      allocate(mesh % map_parent(i))

      i = 0
      do p = p_min, p_max
        if (n_cluster(p) > 0) then
          i = i + 1
          mesh % map_parent(i) % comm      = mesh % comm_world
          mesh % map_parent(i) % proc      = p
          mesh % map_parent(i) % n_cluster = n_cluster(p)
          mesh % map_parent(i) % n_active  = n_active(p)
          mesh % map_parent(i) % n_frozen  = n_frozen(p)
        end if
      end do

      do i = 1, mesh % n_parent
        allocate(mesh%map_parent(i)%id_cluster( mesh%map_parent(i)%n_cluster ))
      end do

      ! ranking of parent elements .............................................

      allocate(parent(mesh % n_parent))

      i = 0
      do i = 1, mesh % n_parent
        allocate(parent(i) % id_active( mesh % map_parent(i) % n_active ))
        allocate(parent(i) % rk_active( mesh % map_parent(i) % n_active ))
        allocate(parent(i) % id_frozen( mesh % map_parent(i) % n_frozen ))
        allocate(parent(i) % rk_frozen( mesh % map_parent(i) % n_frozen ))
      end do

      ! extract parent IDs
      cluster_id = 0
      n_active   = 0
      n_frozen   = 0
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))
          if (element % cluster_id == cluster_id)     cycle ! skip siblings
          if (element % adaptation % parent_proc < 0) cycle ! skip orphans
          cluster_id = element % cluster_id
          p = element % adaptation % parent_proc
          i = parent_proc(p)
          if (element % frozen) then
            n_frozen(p) = n_frozen(p) + 1
            parent(i) % id_frozen(n_frozen(p)) = element % adaptation % parent_id
          else
            n_active(p) = n_active(p) + 1
            parent(i) % id_active(n_active(p)) = element % adaptation % parent_id
          end if
        end associate
      end do

      ! parent ranking
      do i = 1, mesh % n_parent
        call SortIndex(parent(i) % id_active, parent(i) % rk_active)
        call SortIndex(parent(i) % id_frozen, parent(i) % rk_frozen)
      end do

      ! build element lists ....................................................

      cluster_id = 0
      n_active   = 0
      n_frozen   = 0

      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))
          if (element % cluster_id == cluster_id)     cycle ! skip siblings
          if (element % adaptation % parent_proc < 0) cycle ! skip orphans
          cluster_id = element % cluster_id
          p = element % adaptation % parent_proc
          i = parent_proc(p)
          associate(map => mesh % map_parent(i))
            if (element % frozen) then
              n_frozen(p) = n_frozen(p) + 1
              map % id_cluster(n_frozen(p) + map % n_active) = cluster_id
            else
              n_active(p) = n_active(p) + 1
              map % id_cluster(n_active(p)) = cluster_id
            end if
          end associate
        end associate
      end do

    end if

    !$omp end master
    !$omp barrier

  end subroutine BuildMapToParent

  !=============================================================================

end submodule MP_BuildMapToParent

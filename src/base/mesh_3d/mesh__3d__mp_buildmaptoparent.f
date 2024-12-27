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

    integer, allocatable :: n_cluster(:)   ! total num clusters  per proc
    integer, allocatable :: n_active(:)    ! num active clusters per proc
    integer, allocatable :: n_frozen(:)    ! num frozen clusters per proc
    integer, allocatable :: map_proc(:)    ! map entry per proc

    !---------------------------------------------------------------------------
    !> Structure to identify the rank of elements descending from given parent
    !>
    !> `active` comprises the list of active descendents where
    !>   - `active(1,:)` is the parent ID,
    !>   - `active(2,:)` is the parent process,
    !>   - `active(3,:)` is the ID of the related element cluster
    !>
    !> `frozen` is a corresponding list of frozen element clusters.

    type ParentRanking
      integer, allocatable :: active(:,:) ! active element parent ID
      integer, allocatable :: frozen(:,:) ! frozen element parent ID
    end type ParentRanking
    type(ParentRanking), allocatable :: parent(:)

    integer :: n_proc, cluster_id
    integer :: e, i, j, k, p, p_min, p_max

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

      allocate( n_cluster(0:n_proc-1), source = 0)
      allocate( n_active (0:n_proc-1), source = 0)
      allocate( n_frozen (0:n_proc-1), source = 0)
      allocate( map_proc (0:n_proc-1), source = 0)

      p_min = n_proc
      p_max = -1

      ! counts .................................................................

      cluster_id = 0
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))
          if (element % cluster_id == cluster_id)     cycle ! skip siblings
          if (element % adaptation % parent_proc < 0) cycle ! skip orphans
          cluster_id = element % cluster_id
          p = element % adaptation % parent_proc
          p_min = min(p_min, p)
          p_max = max(p_max, p)
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
          map_proc(p) = i
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

      do i = 1, mesh % n_parent
        allocate(parent(i) % active(3, mesh % map_parent(i) % n_active ))
        allocate(parent(i) % frozen(3, mesh % map_parent(i) % n_frozen ))
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
          i = map_proc(p)
          if (element % frozen) then
            n_frozen(p) = n_frozen(p) + 1
            parent(i) % frozen(1,n_frozen(p)) = element % adaptation % parent_id
            parent(i) % frozen(2,n_frozen(p)) = p
            parent(i) % frozen(3,n_frozen(p)) = cluster_id
          else
            n_active(p) = n_active(p) + 1
            parent(i) % active(1,n_active(p)) = element % adaptation % parent_id
            parent(i) % active(2,n_active(p)) = p
            parent(i) % active(3,n_active(p)) = cluster_id
          end if
        end associate
      end do

      ! parent ranking
      do i = 1, mesh % n_parent
        call SortPairs(parent(i) % active)
        call SortPairs(parent(i) % frozen)
      end do

      ! build cluster lists ....................................................

      do i = 1, size(parent)
        k = 0
        do j = 1, size(parent(i)%active, 2)
          k = k + 1
          mesh % map_parent(i) % id_cluster(k) = parent(i) % active(3,j)
        end do
        do j = 1, size(parent(i)%frozen, 2)
          k = k + 1
          mesh % map_parent(i) % id_cluster(k) = parent(i) % frozen(3,j)
        end do
      end do

    end if

    !$omp end master
    !$omp barrier

  end subroutine BuildMapToParent

  !=============================================================================

end submodule MP_BuildMapToParent

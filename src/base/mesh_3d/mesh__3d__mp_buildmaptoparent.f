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

    type ParentRanking
      integer, allocatable :: active_id(:)  ! active cluster ID
      integer, allocatable :: active_pid(:) ! active cluster parent element ID
      integer, allocatable :: active_rk(:)  ! active cluster rank
      integer, allocatable :: frozen_id(:)  ! frozen cluster ID
      integer, allocatable :: frozen_pid(:) ! frozen cluster parent element ID
      integer, allocatable :: frozen_rk(:)  ! frozen cluster rank
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
        allocate(parent(i) % active_id ( mesh % map_parent(i) % n_active ))
        allocate(parent(i) % active_pid( mesh % map_parent(i) % n_active ))
        allocate(parent(i) % active_rk ( mesh % map_parent(i) % n_active ))
        allocate(parent(i) % frozen_id ( mesh % map_parent(i) % n_frozen ))
        allocate(parent(i) % frozen_pid( mesh % map_parent(i) % n_frozen ))
        allocate(parent(i) % frozen_rk ( mesh % map_parent(i) % n_frozen ))
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
            k = n_frozen(p) + 1
            parent(i) % frozen_id(k) = cluster_id
            parent(i) % frozen_pid(k) = element % adaptation % parent_id
            n_frozen(p) = k
          else
            k = n_active(p) + 1
            parent(i) % active_id(k) = cluster_id
            parent(i) % active_pid(k) = element % adaptation % parent_id
            n_active(p) = k
          end if
        end associate
      end do

      ! parent ranking
      do i = 1, mesh % n_parent
        call SortIndex(parent(i) % active_pid, parent(i) % active_rk)
        call SortIndex(parent(i) % frozen_pid, parent(i) % frozen_rk)
      end do

      ! build cluster lists ....................................................

      do i = 1, size(parent)
        associate(id_cluster => mesh % map_parent(i) % id_cluster)
          k = 0
          do j = 1, size(parent(i)%active_rk)
            k = k + 1
            id_cluster(k) = parent(i) % active_id( parent(i) % active_rk(j) )
          end do
          do j = 1, size(parent(i)%frozen_rk)
            k = k + 1
            id_cluster(k) = parent(i) % frozen_id( parent(i) % frozen_rk(j) )
          end do
        end associate
      end do

    end if

    !$omp end master
    !$omp barrier

  end subroutine BuildMapToParent

  !=============================================================================

end submodule MP_BuildMapToParent

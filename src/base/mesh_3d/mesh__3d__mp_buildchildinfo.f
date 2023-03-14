!> summary:  Generation of child mesh information
!> author:   Joerg Stiller
!> date:     2023/03/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_BuildChildInfo
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of information on child mesh partitions
  !>
  !> The routine provides
  !>
  !>   - `mesh % n_child`
  !>   - `mesh % child(:)`
  !>
  !> and requires `mesh % element`
  !>
  !>   -  `mesh % element % adaptation % refinement`
  !>   -  `mesh % element % adaptation % child_proc`

  module subroutine BuildChildInfo(mesh)
    class(Mesh_3D), intent(inout) :: mesh  !< mesh parition

    integer, allocatable :: n_child_elem(:)
    integer, allocatable :: n_child_elem_active(:)
    integer, allocatable :: n_child_elem_frozen(:)
    integer :: n_proc
    integer :: c, e, p, p_min, p_max

    !$omp master

    if (mesh % n_elem > 0) then

      call MPI_Comm_size(mesh % comm_world, n_proc)

      allocate(n_child_elem(0:n_proc-1))
      allocate(n_child_elem_active (0:n_proc-1), source = 0)
      allocate(n_child_elem_frozen (0:n_proc-1), source = 0)

      p_min = n_proc-1
      p_max = 0

      ! count active and frozen child elements
      do e = 1, mesh % n_elem
        p = mesh % element(e) % adaptation % child_proc
        if (p >= 0) then
          p_min = min(p_min, p)
          p_max = min(p_max, p)
          select case(int( mesh % element(e) % adaptation % refinement ))
          case(100) ! regular active
            n_child_elem_active(p) = n_child_elem_active(p) + 8
          case(50)  ! regular frozen
            n_child_elem_frozen(p) = n_child_elem_frozen(p) + 8
          case(1:6) ! face 1:6 frozen
            n_child_elem_frozen(p) = n_child_elem_frozen(p) + 4
          case(7:18) ! edge 1:12 frozen
            n_child_elem_frozen(p) = n_child_elem_frozen(p) + 2
          case(19:26) ! vertex 1:8 frozen
            n_child_elem_frozen(p) = n_child_elem_frozen(p) + 1
          end select
        end if
      end do

      ! total number of child elements and number of child partitions
      mesh % n_child = 0
      do p = p_min, p_max
        n_child_elem(p) = n_child_elem_active(p) + n_child_elem_frozen(p)
        if (n_child_elem(p) > 0) then
          mesh % n_child = mesh % n_child + 1
        end if
      end do

      allocate(mesh % child(mesh % n_child))

      c = 0
      do p = p_min, p_max
        if (n_child_elem(p) > 0) then
          c = c + 1
          mesh % child(c) % proc          = p
          mesh % child(c) % n_elem        = n_child_elem(p)
          mesh % child(c) % n_elem_active = n_child_elem_active(p)
          mesh % child(c) % n_elem_frozen = n_child_elem_frozen(p)
        end if
      end do

    else

      mesh % n_child = 0
      allocate(mesh % child(0))

    end if

    !$omp end master
    !$omp barrier

  end subroutine BuildChildInfo

  !=============================================================================

end submodule MP_BuildChildInfo

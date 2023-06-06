!> summary:  Generation of map to child elements
!> author:   Joerg Stiller
!> date:     2023/04/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_BuildMapToChild
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of map to child elements
  !>
  !> The routine provides
  !>
  !>   - `mesh % n_child`
  !>   - `mesh % map_child(:)`
  !>
  !> and requires `mesh % element`
  !>
  !>   -  `mesh % element % adaptation % refinement`
  !>   -  `mesh % element % adaptation % child_proc`

  module subroutine BuildMapToChild(mesh)
    class(Mesh_3D), intent(inout) :: mesh  !< mesh parition

    integer, allocatable :: n_elem(:)
    integer, allocatable :: n_active(:)
    integer, allocatable :: n_frozen(:)
    integer, allocatable :: child_proc(:)
    integer :: n_proc
    integer :: c, e, p, p_min, p_max

    !$omp master

    if (allocated(mesh % map_child)) then
      deallocate(mesh % map_child)
    end if

    if (mesh % is_top .or. mesh % n_elem < 1) then

      ! top or empty mesh ......................................................

      mesh % n_child = 0
      allocate(mesh % map_child(0))

    else

      ! initialization .........................................................

      call MPI_Comm_size(mesh % comm_world, n_proc)

      allocate( n_elem     (0:n_proc-1))
      allocate( n_active   (0:n_proc-1), source = 0)
      allocate( n_frozen   (0:n_proc-1), source = 0)
      allocate( child_proc (0:n_proc-1), source = 0)

      p_min = n_proc-1
      p_max = 0

      ! counts .................................................................

      do e = 1, mesh % n_elem
        p = mesh % element(e) % adaptation % child_proc
        if (p >= 0) then
          p_min = min(p_min, p)
          p_max = min(p_max, p)
          select case(mesh % element(e) % adaptation % refinement)
          case(100)
            n_active(p) = n_active(p) + 1
          case(1:26, 50)
            n_frozen(p) = n_frozen(p) + 1
          end select
        end if
      end do

      ! setup of maps ..........................................................

      c = 0
      do p = p_min, p_max
        n_elem(p) = n_active(p) + n_frozen(p)
        if (n_elem(p) > 0) then
          c = c + 1
          child_proc(p) = c
        end if
      end do
      mesh % n_child = c
      allocate(mesh % map_child(c))

      c = 0
      do p = p_min, p_max
        if (n_elem(p) > 0) then
          c = c + 1
          mesh % map_child(c) % proc     = p
          mesh % map_child(c) % n_elem   = n_elem(p)
          mesh % map_child(c) % n_active = n_active(p)
          mesh % map_child(c) % n_frozen = n_frozen(p)
        end if
      end do

      do c = 1, mesh % n_child
        allocate(mesh % map_child(c) % id_elem( mesh % map_child(c) % n_elem ))
      end do

      ! build element lists ....................................................

      n_active = 0
      n_frozen = 0

      associate(map => mesh % map_child)
        do e = 1, mesh % n_elem
          p = mesh % element(e) % adaptation % child_proc
          if (p >= 0) then
            c = child_proc(p)
            select case(mesh % element(e) % adaptation % refinement)
            case(100)
              n_active(p) = n_active(p) + 1
              map(c) % id_elem(n_active(p)) = e
            case(1:26, 50)
              n_frozen(p) = n_frozen(p) + 1
              map(c) % id_elem(map(c) % n_active + n_frozen(p)) = e
            end select
          end if
        end do
      end associate

    end if

    !$omp end master
    !$omp barrier

  end subroutine BuildMapToChild

  !=============================================================================

end submodule MP_BuildMapToChild

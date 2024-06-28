!> summary:  3D multilevel mesh
!> author:   Joerg Stiller
!> date:     2024/06/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Mesh__3D
  use Mesh__3D
  use Execution_Control
  use Partitioner_Interface__3D
  implicit none
  private

  public :: ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel mesh
  !>
  !> The `mesh` component contains a sequence of meshes the from level `1` up to
  !> `l_top = size(mesh)`. The relation between mesh levels is described by the
  !> `refinement` pattern. It ranges from `1` to `l_top-1` and can take the
  !> following values
  !>
  !>   | `refinement` | child level generation | application  |
  !>   |:------------:| ---------------------- | ------------ |
  !>   |     'c'      | global cloning         | p-refinement |
  !>   |     'z'      | zonal cloning          | p-adaptivity |
  !>   |     'r'      | regular refinement     | h-refinement |
  !>   |     'l'      | local refinement       | h-adaptivity |
  !>

  type ML_Mesh_3D
    type(Mesh_3D), allocatable :: mesh(:)       !< mesh partitions
    character,     allocatable :: refinement(:) !< refinement pattern
  end type ML_Mesh_3D

contains

  !-----------------------------------------------------------------------------
  !> New multilevel mesh by global refinement of given mesh

  function New_ML_Mesh_3D(mesh, refinement, part_opt) result(this)
    type(Mesh_3D),                 intent(in) :: mesh
    integer,                       intent(in) :: refinement(:)
    class(PartitioningOptions_3D), intent(in) :: opt(:)
    type(ML_Mesh_3D) :: this

    integer :: l, l_top

    ! preliminaries ............................................................

    l_top = size(opt)

    if (size(refinement) /= l_top-1) then
      call Error( 'New_ML_Mesh_3D', 'incompatible arguments', 'ML__Mesh__3D')
    end if

    do l = 1, size(refinement)
      if (scan('cr',refinement(l)) == 0) then
        call Error('New_ML_Mesh_3D','invalid refinement pattern','ML__Mesh__3D')
      end if
    end do

    allocate(this % mesh(l_top))
    this % refinement = refinement

    ! partition root level .....................................................

    if (old_mesh%n_parts == opt%n_parts) then
      this % mesh(1) = mesh
    else
      call RootMeshPartitioning_3D( opt       = part_opt(1)  &
                                  , old_mesh  = mesh         &
                                  , new_mesh  = this%mesh(1) )
    end if

    ! create higher levels .....................................................

    do l = 1, l_top-1

       select case(refinement(l))
       case('c')


    end do

! partitioning?
    if (l_top == 1) then
      ! single level case
      allocate(this % mesh(1), source = mesh)
      allocate(this % refinement(0))
      return
    end if

    allocate(old_mesh(1), source = mesh)
    allocate(new_mesh(2))

  end function New_ML_Mesh_3D

  !=============================================================================

end module ML__Mesh__3D

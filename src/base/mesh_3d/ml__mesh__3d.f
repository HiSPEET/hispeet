!> summary:  3D multilevel mesh
!> author:   Joerg Stiller
!> date:     2024/06/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Mesh__3D
  use Kind_Parameters
  use Constants
  use Execution_Control
  use XMPI
  use Mesh__3D
  use Data_Exchange__3D
  use Child_Mesh_Adaptation__3D
  use Partitioner_Interface__3D
  use Process_Adaptation_Pattern__3D
  use Root_Mesh_Partitioning__3D
  implicit none
  private

  public :: ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel mesh

  type ML_Mesh_3D
    type(Mesh_3D), allocatable :: mesh(:) !< mesh partitions
  end type ML_Mesh_3D

  ! constructor
  interface ML_Mesh_3D
    module procedure New_ML_Mesh_3D
  end interface

  !-----------------------------------------------------------------------------
  !> 3D multilevel mesh generation options

  type ML_MeshOptions_3D
    integer :: l_top   = 1                     !< top level
    integer :: l_adapt = huge(1)               !< first level to be adapted ≥1
    character, allocatable :: refinement(:)    !< refinement type {'c','s'}
    integer,   allocatable :: adapt_bnd(:)     !< boundaries to be adapted
    real(RNP), allocatable :: adapt_box(:,:,:) !< boxes to be adapted [3,2,*]
    real(RNP), allocatable :: adapt_tol(:)     !< tolerance per level
    type(PartitioningOptions_3D), allocatable :: partition(:)
  end type ML_MeshOptions_3D

contains

  !-----------------------------------------------------------------------------
  !> New multilevel mesh by global refinement of given mesh

  function New_ML_Mesh_3D(mesh, opt) result(this)
    type(Mesh_3D),            intent(in) :: mesh
    class(ML_MeshOptions_3D), intent(in) :: opt
    type(ML_Mesh_3D) :: this

    if (opt%l_top > 1 .and. opt%l_adapt <= opt%l_top) then
      this = Adapted_ML_Mesh_3D(mesh, opt)
    else
      this = Global_ML_Mesh_3D(mesh, opt)
    end if

  end function New_ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> Multilevel mesh created by global refinement of given mesh

  function Global_ML_Mesh_3D(mesh, opt) result(this)
    type(Mesh_3D),            intent(in) :: mesh
    class(ML_MeshOptions_3D), intent(in) :: opt
    type(ML_Mesh_3D) :: this

    integer :: l

    ! preliminaries ............................................................

    if (mesh%part == 0) then
      write(*,'(/,2X,A)') 'creating global multilevel mesh'
    end if

    allocate(this % mesh(opt%l_top))

    ! partition root level .....................................................

    if (mesh%part == 0) then
      write(*,'(4X,A)') 'partitioning root mesh'
    end if

    if (mesh%n_parts == opt%partition(1)%n_parts) then
      this % mesh(1) = mesh
    else
      call RootMeshPartitioning_3D( opt       = opt%partition(1) &
                                  , old_mesh  = mesh             &
                                  , new_mesh  = this%mesh(1)     )
    end if

    ! create higher levels .....................................................

    do l = 1, opt%l_top-1
      if (mesh%part == 0) then
        write(*,'(4X,A,I0)') 'creating level ', l+1
      end if
      this % mesh(l) % refinement = opt % refinement(l)
      call this % mesh(l) % element % MarkForRefinement()
      call ProcessAdaptationPattern_3D(this%mesh(l))
      call ChildMeshAdaptation_3D( opt        = opt  % partition(l+1) &
                                 , parent     = this % mesh(l)        &
                                 , new_child  = this % mesh(l+1)      )
    end do

  end function Global_ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> Multilevel mesh created by adaptive refinement of given mesh

  function Adapted_ML_Mesh_3D(mesh, opt) result(this)
    type(Mesh_3D),            intent(in) :: mesh
    class(ML_MeshOptions_3D), intent(in) :: opt
    type(ML_Mesh_3D) :: this

    type(Mesh_3D),             allocatable, save :: old_mesh(:)
    type(DataExchangePlan_3D), allocatable, save :: exch_plan(:)
    logical, save :: finish, finish_loc

    integer :: n_box, n_bnd
    integer :: i, e, l, m

    ! preliminaries ............................................................

    if (mesh%part == 0) then
      write(*,'(/,2X,A)') 'creating adapted multilevel mesh'
    end if

    if (allocated(opt%adapt_bnd)) then
      n_bnd = size(opt%adapt_bnd, 1)
    else
      n_bnd = 0
    end if

    if (allocated(opt%adapt_box)) then
      n_box = size(opt%adapt_box, 3)
    else
      n_box = 0
    end if

    allocate(this % mesh(1), source = mesh)

    ADAPTATION: do m = 1, opt%l_top - 1

      if (mesh%part == 0) then
        write(*,'(4X,9G0)') 'adaptation cycle ',m,'/',opt%l_top-1
      end if

      ! set adaptation marks ...................................................

      do l = 1, m
        associate(parent => this % mesh(l))

          parent % refinement = opt % refinement(l)

          if (l < opt%l_adapt) then

            call parent % element % MarkForRefinement()

          else

            MARK_INIT: if (parent % is_root) then
              call parent % element % Unmark()
            else
              call parent % element % MarkForRemoval()
            end if MARK_INIT

            MARK_BND: if (n_bnd > 0) then
              do e = 1, parent % n_elem
                associate(element => parent%element(e), bnd => opt%adapt_bnd)

                  do i = 1, n_bnd
                    if (any(element%face%boundary == bnd(i))) then
                      call element % MarkForRefinement()
                      exit
                    end if
                  end do

                end associate
              end do
            end if MARK_BND

            MARK_BOX: if (n_box > 0) then
              do e = 1, parent % n_elem
                associate( element => parent % element(e)                  &
                         , x_cub   => parent % element(e) % geometry % x_c &
                         , x_box   => opt % adapt_box                      &
                         , tol     => opt % adapt_tol(l)                   )

                  if (ElementIntersectsBox(n_box, x_box, x_cub, tol)) then
                    call element % MarkForRefinement()
                  end if

                end associate
              end do
            end if MARK_BOX

          end if

        end associate
      end do

      ! check wether to continue adaptation ....................................

      finish_loc = .not. any(this%mesh(m)%element%adaptation%mark > 0)
      call XMPI_Allreduce(finish_loc, finish, MPI_LAND, mesh%comm_world)

      if (finish) then
        if (mesh%part == 0) then
          write(*,'(4X,9G0)') 'adaptation finished with ',m, 'levels'
        end if
        exit ADAPTATION
      end if

      ! make adaptation pattern consistent .....................................

      do l = m, 2, -1
        call GlobalizeAdaptationPattern_3D(this%mesh(l))
        call RestrictAdaptationPattern_3D(this%mesh(l), this%mesh(l-1))
      end do
      call GlobalizeAdaptationPattern_3D(this%mesh(1))

      ! adaptation cycle .......................................................

      ! data structures
      call move_alloc(this%mesh, old_mesh)  ! save current mesh
      allocate(this%mesh(m+1))              ! prepare new mesh
      allocate(exch_plan(m))                ! exchange plan for retained data

      ! root level
      call ProcessAdaptationPattern_3D(old_mesh(1))
      if (max(old_mesh(1)%n_parts, opt%partition(1)%n_parts) == 1) then
        this % mesh(1) = old_mesh(1)
      else if (old_mesh(1) % is_top) then
        call RootMeshPartitioning_3D( opt       = opt%partition(1) &
                                    , old_mesh  = old_mesh(1)      &
                                    , new_mesh  = this%mesh(1)     &
                                    , exch_plan = exch_plan(1)     )
      else
        call RootMeshPartitioning_3D( opt       = opt%partition(1) &
                                    , old_mesh  = old_mesh(1)      &
                                    , new_mesh  = this%mesh(1)     &
                                    , child     = old_mesh(2)      &
                                    , exch_plan = exch_plan(1)     )
      end if

      ! higher levels
      do l = 1, m

        if (l > 1) then
          call ProcessAdaptationPattern_3D(this%mesh(l))
        end if

        select case(m-l)
        case(0)
          call ChildMeshAdaptation_3D( opt        = opt%partition(l+1) &
                                     , parent     = this%mesh(l)       &
                                     , new_child  = this%mesh(l+1)     )
        case(1)
          call ChildMeshAdaptation_3D( opt        = opt%partition(l+1) &
                                     , parent     = this%mesh(l)       &
                                     , new_child  = this%mesh(l+1)     &
                                     , old_child  = old_mesh(l+1)      &
                                     , exch_plan  = exch_plan(l+1)     )
        case(2:)
          call ChildMeshAdaptation_3D( opt        = opt%partition(l+1) &
                                     , parent     = this%mesh(l)       &
                                     , new_child  = this%mesh(l+1)     &
                                     , old_child  = old_mesh(l+1)      &
                                     , grandchild = old_mesh(l+2)      &
                                     , exch_plan  = exch_plan(l+1)     )
        end select

      end do

      ! release workspace
      deallocate(old_mesh, exch_plan)

    end do ADAPTATION

  end function Adapted_ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> Check for intersection of element cuboid with any of given boxes

  logical function ElementIntersectsBox(n_box, x_box, x_cub, tol)
    integer,   intent(in) :: n_box            !< number of boxes
    real(RNP), intent(in) :: x_box(3,2,n_box) !< boxes
    real(RNP), intent(in) :: x_cub(0:3,3)     !< element cuboid
    real(RNP), intent(in) :: tol              !< tolerance

    real(RNP) :: x_cub_min(3), x_cub_max(3)
    real(RNP) :: x_box_min(3), x_box_max(3)
    real(RNP) :: dx(3), tol2
    integer   :: i

    ElementIntersectsBox = .false.

    tol2 = tol * tol

    ! cuboid center position ∓ approximate half width in directions 1:3
    x_cub_min = x_cub(0,:) - (x_cub(1,:) + x_cub(2,:) + x_cub(3,:))
    x_cub_max = x_cub(0,:) + (x_cub(1,:) + x_cub(2,:) + x_cub(3,:))

    do i = 1, n_box
      x_box_min = x_box(3,1,i)
      x_box_max = x_box(3,2,i)
      dx = max(x_box_min - x_cub_min, x_cub_max - x_box_min, ZERO)
      if (dx(1)**2 + dx(2)**2 + dx(3)**2 <= tol2) then
        ElementIntersectsBox = .true.
        exit
      end if
    end do

  end function ElementIntersectsBox

  !=============================================================================

end module ML__Mesh__3D

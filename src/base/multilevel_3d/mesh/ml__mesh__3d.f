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
  use Globalize_Adaptation_Pattern__3D
  use Partitioner_Interface__3D
  use Process_Adaptation_Pattern__3D
  use Restrict_Adaptation_Pattern__3D
  use Root_Mesh_Partitioning__3D
  implicit none
  private

  public :: ML_Mesh_3D
  public :: ML_Mesh_Options_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel mesh

  type ML_Mesh_3D
    type(Mesh_3D), allocatable :: mesh(:) !< mesh partitions
  contains
    procedure :: Init_ML_Mesh_3D
    procedure :: ReadHDF5
    procedure :: WriteHDF5
  end type ML_Mesh_3D

  ! constructor
  interface ML_Mesh_3D
    procedure New_ML_Mesh_3D
  end interface

  !-----------------------------------------------------------------------------
  !> 3D multilevel mesh generation options

  type ML_Mesh_Options_3D
    integer :: l_top   = 1                     !< top level
    integer :: l_adapt = huge(1)               !< first level to be adapted ≥1
    character, allocatable :: refinement(:)    !< refinement type {'c','s'}
    integer,   allocatable :: adapt_bnd(:)     !< boundaries to be adapted
    real(RNP), allocatable :: adapt_box(:,:,:) !< boxes to be adapted (3,2,*)
    real(RNP), allocatable :: adapt_tol(:)     !< tolerance per level
    type(PartitioningOptions_3D), allocatable :: partition(:)
  contains
    procedure :: Bcast => Bcast_ML_Mesh_Options_3D
  end type ML_Mesh_Options_3D

  ! constructor
  interface ML_Mesh_Options_3D
    procedure Read_ML_Mesh_Options_3D
  end interface

  !=============================================================================
  ! ML_Mesh_3D: separate type-bound procedures

  interface

    !---------------------------------------------------------------------------
    !> Read multilevel mesh partition from HDF5 file

    module subroutine ReadHDF5(this, file, comm)
      class(ML_Mesh_3D), intent(inout) :: this !< multilevel mesh partition
      character(len=*),  intent(in)    :: file !< name of HDF5 file
      type(MPI_Comm),    intent(in)    :: comm !< MPI "world" communicator
    end subroutine ReadHDF5

    !---------------------------------------------------------------------------
    !> Write multilevel mesh partition into HDF5 file

    module subroutine WriteHDF5(this, file)
      class(ML_Mesh_3D), intent(in) :: this !< multilevel mesh partition
      character(len=*),  intent(in) :: file !< name of HDF5 file
    end subroutine WriteHDF5

  end interface

contains

  !=============================================================================
  ! ML_Mesh_3D procedures

  !-----------------------------------------------------------------------------
  !> New multilevel mesh by global refinement of given mesh

  function New_ML_Mesh_3D(mesh, opt) result(this)
    type(Mesh_3D),             intent(in) :: mesh
    class(ML_Mesh_Options_3D), intent(in) :: opt
    type(ML_Mesh_3D) :: this

    call Init_ML_Mesh_3D(this, mesh, opt)

  end function New_ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> Create multilevel mesh by global refinement of given mesh

  subroutine Init_ML_Mesh_3D(this, mesh, opt)
    class(ML_Mesh_3D),         intent(inout) :: this
    type(Mesh_3D),             intent(in)    :: mesh
    class(ML_Mesh_Options_3D), intent(in)    :: opt

    integer :: l

    do l = 1, opt%l_top-1
      if (scan(opt%refinement(l),'cs') == 0) then
        call Error('Init_ML_Mesh_3D','invalid opt%refinement')
      end if
    end do

    if (opt%l_top > 1 .and. opt%l_adapt <= opt%l_top) then
      call Create_Adapted_ML_Mesh_3D(this, mesh, opt)
    else
      call Create_Global_ML_Mesh_3D(this, mesh, opt)
    end if

  end subroutine Init_ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> Create multilevel mesh by global refinement of given mesh

  subroutine Create_Global_ML_Mesh_3D(this, mesh, opt)
    class(ML_Mesh_3D),         intent(inout) :: this
    type(Mesh_3D),             intent(in)    :: mesh
    class(ML_Mesh_Options_3D), intent(in)    :: opt

    integer :: l

    ! preliminaries ............................................................

    if (mesh%part == 0) then
      write(*,'(/,A)') 'creating global multilevel mesh'
    end if

    allocate(this % mesh(opt%l_top))

    ! partition root level .....................................................

    if (mesh%part == 0) then
      write(*,'(2X,A)') 'partitioning root mesh'
    end if

    if (mesh%n_parts == opt%partition(1)%n_parts) then
      this % mesh(1) = mesh
    else
      call RootMeshPartitioning_3D( opt      = opt%partition(1) &
                                  , old_mesh = mesh             &
                                  , new_mesh = this%mesh(1)     )
    end if

    ! create higher levels .....................................................

    do l = 1, opt%l_top-1
      if (mesh%part == 0) then
        write(*,'(2X,A,I0)') 'creating level ', l+1
      end if
      this % mesh(l) % refinement = opt % refinement(l)
      if (this%mesh(l)%n_elem > 0) then
        call this % mesh(l) % element % MarkForRefinement()
      end if
      call ProcessAdaptationPattern_3D(this%mesh(l))
      call ChildMeshAdaptation_3D( opt       = opt  % partition(l+1) &
                                 , parent    = this % mesh(l)        &
                                 , new_child = this % mesh(l+1)      )
    end do

  end subroutine Create_Global_ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> Create multilevel mesh by adaptive refinement of given mesh

  subroutine Create_Adapted_ML_Mesh_3D(this, mesh, opt)
    class(ML_Mesh_3D),         intent(inout) :: this
    type(Mesh_3D),             intent(in)    :: mesh
    class(ML_Mesh_Options_3D), intent(in)    :: opt

    type(Mesh_3D),             allocatable, save :: old_mesh(:)
    type(DataExchangePlan_3D), allocatable, save :: exch_plan(:)
    logical, save :: finish, finish_loc

    integer :: n_box, n_bnd
    integer :: i, e, l, m

    ! preliminaries ............................................................

    if (mesh%part == 0) then
      write(*,'(/,A)') 'creating adapted multilevel mesh'
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
        write(*,'(2X,9G0)') 'adaptation cycle ',m,'/',opt%l_top-1
      end if

      ! set adaptation marks ...................................................

      do l = 1, m
        associate(parent => this % mesh(l))

          parent % refinement = opt % refinement(l)

          if (parent%part < 0) cycle

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

      finish_loc = .true.
      do e = 1, this%mesh(m)%n_elem
        if (this%mesh(m)%element(e)%adaptation%mark > 0) then
          finish_loc = .false.
          exit
        end if
      end do
      call XMPI_Allreduce(finish_loc, finish, MPI_LAND, mesh%comm_world)

      if (finish) then
        if (mesh%part == 0) then
          write(*,'(2X,9G0)') 'adaptation finished with ',m, ' levels'
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

        this%mesh(l)%refinement = opt%refinement(l)

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

  end subroutine Create_Adapted_ML_Mesh_3D

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
    dx = abs(x_cub(1,:)) + abs(x_cub(2,:)) + abs(x_cub(3,:))
    x_cub_min = x_cub(0,:) - dx
    x_cub_max = x_cub(0,:) + dx

    do i = 1, n_box
      x_box_min = x_box(:,1,i) - tol
      x_box_max = x_box(:,2,i) + tol
      if (all(x_box_min < x_cub_max .and. x_box_max > x_cub_min)) then
        ElementIntersectsBox = .true.
        exit
      end if
    end do

  end function ElementIntersectsBox

  !=============================================================================
  ! ML_Mesh_Options_3D procedures

  !-----------------------------------------------------------------------------
  !> Generate multilevel mesh options from namelist input

  function Read_ML_Mesh_Options_3D(unit, n_proc) result(this)
    integer,           intent(in) :: unit   !< unit connected for namelist input
    integer, optional, intent(in) :: n_proc !< number of processes [∞]
    type(ML_Mesh_Options_3D) :: this

    ! static input variables ...................................................

    ! adaptation
    integer :: l_top   = 1         ! top level
    integer :: l_adapt = huge(1)   ! first level to be adapted ≥1
    integer :: n_bnd   = 0         ! number of boundaries to be adapted
    integer :: n_box   = 0         ! number of boxes to be adapted

    ! partitioning
    integer :: n_parts_root   = 1  ! number of partitions at root level
    integer :: n_parts_growth = 1  ! partition number growth rate

    namelist/ml_mesh_options_3d__static/ l_top, l_adapt, n_bnd, n_box
    namelist/ml_mesh_options_3d__static/ n_parts_root, n_parts_growth

    ! dynamic input variables ..................................................

    character, allocatable :: refinement(:)    ! refinement type {'c','s'}
    integer,   allocatable :: adapt_bnd(:)     ! boundaries to be adapted
    real(RNP), allocatable :: adapt_box(:,:,:) ! boxes to be adapted (3,2,*)
    real(RNP), allocatable :: adapt_tol(:)     ! tolerance per level

    namelist /ml_mesh_options_3d__dynamic/ refinement
    namelist /ml_mesh_options_3d__dynamic/ adapt_bnd, adapt_box, adapt_tol

    ! auxiliary variables ......................................................

    integer :: l

    ! read static input variables ..............................................

    rewind(unit)
    read(unit, nml = ml_mesh_options_3d__static)

    ! read dynamic input variables .............................................

    allocate(refinement(l_top-1),   source = ' ')
    allocate(adapt_bnd (n_bnd),     source = 0 )
    allocate(adapt_box (3,2,n_box), source = huge(ONE))
    allocate(adapt_tol (l_top-1),   source = ZERO)

    read(unit, nml = ml_mesh_options_3d__dynamic)

    ! create multilevel mesh options ...........................................

    this % l_top   = l_top
    this % l_adapt = l_adapt

    call move_alloc( refinement, this % refinement )
    call move_alloc( adapt_bnd , this % adapt_bnd  )
    call move_alloc( adapt_box , this % adapt_box  )
    call move_alloc( adapt_tol , this % adapt_tol  )

    ! partitioning options
    allocate(this % partition(l_top))
    associate(partition => this % partition)
      partition(1) % mode    = 1
      partition(1) % n_parts = n_parts_root
      do l = 2, l_top
        partition(l) % mode    = 2
        partition(l) % n_parts = n_parts_growth * partition(l-1)%n_parts
      end do
      if (present(n_proc)) then
        partition % n_parts = min(partition%n_parts, n_proc)
      end if
    end associate

  end function Read_ML_Mesh_Options_3D

  !-----------------------------------------------------------------------------
  !> Broadcasting multilevel mesh options

  subroutine Bcast_ML_Mesh_Options_3D(this, root, comm)
    class(ML_Mesh_Options_3D), intent(inout) :: this !< options
    integer       , intent(in) :: root !< rank root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    integer :: l, rank, n_bnd, n_box

    associate(l_top => this%l_top, l_adapt => this%l_adapt)

      call MPI_Comm_rank(comm, rank)

      if (rank == root) then
        n_bnd = size(this % adapt_bnd,1)
        n_box = size(this % adapt_box,3)
      end if

      call XMPI_Bcast(l_top  , root, comm)
      call XMPI_Bcast(l_adapt, root, comm)
      call XMPI_Bcast(n_bnd  , root, comm)
      call XMPI_Bcast(n_box  , root, comm)

      if (rank /= root) then
        allocate( this % refinement(l_top-1)  )
        allocate( this % adapt_bnd(n_bnd)     )
        allocate( this % adapt_box(3,2,n_box) )
        allocate( this % adapt_tol(l_top-1)   )
        allocate( this % partition(l_top)     )
      end if

      call XMPI_Bcast(this % refinement, root, comm)
      call XMPI_Bcast(this % adapt_bnd , root, comm)
      call XMPI_Bcast(this % adapt_box , root, comm)
      call XMPI_Bcast(this % adapt_tol , root, comm)

      do l = 1, l_top
        call this % partition(l) % Bcast( 0, comm )
      end do

    end associate

  end subroutine Bcast_ML_Mesh_Options_3D

  !=============================================================================

end module ML__Mesh__3D

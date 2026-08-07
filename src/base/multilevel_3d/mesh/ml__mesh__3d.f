!> summary:  3D multilevel mesh
!> author:   Joerg Stiller
!> date:     2024/06/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Mesh__3D
  use, intrinsic :: ISO_Fortran_Env

  use Kind_Parameters
  use Constants
  use Logging_Levels
  use Execution_Control
  use XMPI
  use Mesh__3D
  use Data_Exchange__3D
  use Child_Mesh_Adaptation__3D
  use Globalize_Adaptation_Pattern__3D
  use Mesh_Partitioner__3D
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
    procedure :: MarkByOptions
    procedure :: Adapt
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
    integer :: l_top   =  1                    !< initial top level
    integer :: l_max   = -1                    !< ultimate top level
    integer :: l_adapt = huge(1)               !< first level to be adapted ≥1
    character, allocatable :: refinement(:)    !< refinement type {'c','s'}
    integer,   allocatable :: adapt_bnd(:)     !< boundaries to be adapted (*)
    real(RNP), allocatable :: adapt_box(:,:,:) !< boxes to be adapted (3,2,*)
    type(MeshPartitionerOptions_3D), allocatable :: partition(:)
  contains
    procedure :: SetUp => SetUp_ML_Mesh_Options_3D
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
    !> Adapt multilevel mesh

    module subroutine Adapt(this, part_opt, x_plan)
      class(ML_Mesh_3D),                      intent(inout) :: this
      class(MeshPartitionerOptions_3D),       intent(in)    :: part_opt(:)
      type(DataExchangePlan_3D), allocatable, intent(out)   :: x_plan(:)
    end subroutine Adapt

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
  !> Create multilevel mesh using prescribed refinement of given mesh

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
  !> Create multilevel mesh by local refinement of given mesh

  subroutine Create_Global_ML_Mesh_3D(this, mesh, opt)
    class(ML_Mesh_3D),         intent(inout) :: this
    type(Mesh_3D),             intent(in)    :: mesh
    class(ML_Mesh_Options_3D), intent(in)    :: opt

    type(MeshPartitionerOptions_3D) :: part_opt
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
      part_opt = opt%partition(1)
      part_opt % n_con_root = 1
      call RootMeshPartitioning_3D( opt      = part_opt     &
                                  , old_mesh = mesh         &
                                  , new_mesh = this%mesh(1) )
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
      part_opt = opt%partition(l+1)
      part_opt % n_con_child = 1
      call ProcessAdaptationPattern_3D(this%mesh(l))
      call ChildMeshAdaptation_3D( opt       = part_opt         &
                                 , parent    = this % mesh(l)   &
                                 , new_child = this % mesh(l+1) )
    end do

  end subroutine Create_Global_ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> Create multilevel mesh by adaptive refinement of given mesh

  subroutine Create_Adapted_ML_Mesh_3D(this, mesh, opt)
    class(ML_Mesh_3D),         intent(inout) :: this
    type(Mesh_3D),             intent(in)    :: mesh
    class(ML_Mesh_Options_3D), intent(in)    :: opt

    type(DataExchangePlan_3D), allocatable, save :: x_plan(:)

    integer :: l

    ! preliminaries ............................................................

    if (mesh%part == 0) then
      write(*,'(/,A)') 'creating adapted multilevel mesh'
    end if

    allocate(this % mesh(1), source = mesh)

    ADAPTATION: do l = 1, opt%l_top - 1

      if (mesh%part == 0) then
        write(*,'(2X,9G0)') 'adaptation cycle ',l,'/',opt%l_top-1
      end if

      call this % MarkByOptions(opt)
      call this % Adapt(opt%partition, x_plan)

    end do ADAPTATION

  end subroutine Create_Adapted_ML_Mesh_3D

  !-----------------------------------------------------------------------------
  !> Set adaptation marks according to given ML mesh options

  subroutine MarkByOptions(this, opt)
    class(ML_Mesh_3D),         intent(inout) :: this
    class(ML_Mesh_Options_3D), intent(in)    :: opt

    integer :: n_box, n_bnd
    integer :: i, e, l

    ! preliminaries ............................................................

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

    ! set marks ................................................................

    do l = 1, size(this % mesh)
      associate(mesh => this % mesh(l))

        mesh % refinement = opt % refinement(l)

        if (mesh%part < 0) then
          cycle

        else if (l > opt % l_max) then
          call mesh % element % MarkForRemoval()

        else if (l == opt % l_max) then
          call mesh % element % Unmark()

        else if (l < opt%l_adapt) then
          call mesh % element % MarkForRefinement()

        else

          MARK_INIT: if (mesh % is_root) then
            call mesh % element % Unmark()
          else
            call mesh % element % MarkForRemoval()
          end if MARK_INIT

          MARK_BND: if (n_bnd > 0) then
            do e = 1, mesh % n_elem
              associate(element => mesh%element(e), bnd => opt%adapt_bnd)
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
            do e = 1, mesh % n_elem
              associate( element => mesh % element(e)                  &
                       , x_cub   => mesh % element(e) % geometry % x_c &
                       , x_box   => opt % adapt_box                    )

                if (ElementIntersectsBox(n_box, x_box, x_cub)) then
                  call element % MarkForRefinement()
                end if

              end associate
            end do
          end if MARK_BOX

        end if

      end associate
    end do

  end subroutine MarkByOptions

  !-----------------------------------------------------------------------------
  !> Check for intersection of element cuboid with any of given boxes

  logical function ElementIntersectsBox(n_box, x_box, x_cub)
    integer,   intent(in) :: n_box            !< number of boxes
    real(RNP), intent(in) :: x_box(3,2,n_box) !< boxes
    real(RNP), intent(in) :: x_cub(0:3,3)     !< element cuboid

    real(RNP) :: dx(3), x_cub_min(3), x_cub_max(3)
    integer   :: i

    ! cuboid center position ∓ approximate half width in directions 1:3
    dx = abs(x_cub(1,:)) + abs(x_cub(2,:)) + abs(x_cub(3,:))
    x_cub_min = x_cub(0,:) - dx
    x_cub_max = x_cub(0,:) + dx

    ElementIntersectsBox = .false.
    do i = 1, n_box
      if (any(x_cub_max <= x_box(1:3,1,i))) cycle ! cuboid left  outside the box
      if (any(x_cub_min >= x_box(1:3,2,i))) cycle ! cuboid right outside the box
      ElementIntersectsBox = .true.
      exit
    end do

  end function ElementIntersectsBox

  !=============================================================================
  ! ML_Mesh_Options_3D procedures

  !-----------------------------------------------------------------------------
  !> Set up multilevel mesh options

  subroutine SetUp_ML_Mesh_Options_3D(this, l_top, l_max, l_adapt, n_bnd, n_box)
    class(ML_Mesh_Options_3D), intent(inout) :: this
    integer,           intent(in) :: l_top
    integer, optional, intent(in) :: l_max
    integer, optional, intent(in) :: l_adapt
    integer, optional, intent(in) :: n_bnd
    integer, optional, intent(in) :: n_box

    integer :: l_adapt_, l_max_, n_bnd_, n_box_

    ! process optional arguments ...............................................

    if (present(l_max)) then
       l_max_ = max(l_top, l_max)
    else
       l_max_ = l_top
    end if

    if (present(l_adapt)) then
       l_adapt_ = l_adapt
    else
       l_adapt_ = huge(1)
    end if

    if (present(n_bnd)) then
       n_bnd_ = n_bnd
    else
       n_bnd_ = 0
    end if

    if (present(n_box)) then
       n_box_ = n_box
    else
       n_box_ = 0
    end if

    ! deallocate dynamic components ............................................

    if (allocated( this % refinement )) deallocate( this % refinement )
    if (allocated( this % adapt_bnd  )) deallocate( this % adapt_bnd  )
    if (allocated( this % adapt_box  )) deallocate( this % adapt_box  )
    if (allocated( this % partition  )) deallocate( this % partition  )

    this % l_top   = l_top
    this % l_max   = l_max_
    this % l_adapt = l_adapt_

    allocate(this % refinement(l_max_ - 1)  , source = ' ' )
    allocate(this % adapt_bnd(n_bnd_))
    allocate(this % adapt_box(3, 2, n_box_))
    allocate(this % partition(l_max_))

  end subroutine SetUp_ML_Mesh_Options_3D

  !-----------------------------------------------------------------------------
  !> Generate multilevel mesh options from namelist input

  function Read_ML_Mesh_Options_3D(unit, n_proc) result(this)
    integer,           intent(in) :: unit   !< unit connected for namelist input
    integer, optional, intent(in) :: n_proc !< number of processes [∞]
    type(ML_Mesh_Options_3D) :: this

    ! static input variables ...................................................

    ! adaptation
    integer :: l_top   =  1         ! top level
    integer :: l_max   = -1         ! top level
    integer :: l_adapt =  huge(1)   ! first level to be adapted > 1
    integer :: n_bnd   =  0         ! number of boundaries to be adapted
    integer :: n_box   =  0         ! number of boxes to be adapted

    ! partitioning
    integer :: n_parts_root   = -1      ! number of partitions at root level
    integer :: n_parts_growth = -1      ! partition number growth rate
    integer :: partitioner    =  1      ! partitioning method, 1/2: SFC/graph
    integer :: c_active       = -1      ! cost of active child elements [preset]
    integer :: c_frozen       = -1      ! cost of frozen child elements [preset]
    logical :: split          = .false. ! use parent subdivision for child mesh
    integer :: n_con_root     =  1      ! num constraints for root     ≤ 3
    integer :: n_con_child    =  1      ! num constraints for children ≤ 2
    integer :: n_con_sub      = 10      ! max num sublevels to be weighted

    namelist/ml_mesh_options_3d__static/ l_top, l_max, l_adapt
    namelist/ml_mesh_options_3d__static/ n_bnd, n_box
    namelist/ml_mesh_options_3d__static/ n_parts_root, n_parts_growth
    namelist/ml_mesh_options_3d__static/ partitioner, c_active, c_frozen, split
    namelist/ml_mesh_options_3d__static/ n_con_root, n_con_child, n_con_sub

    ! dynamic input variables ..................................................

    character, allocatable :: refinement(:)    ! refinement type {c,s} (l_max-1)
    integer,   allocatable :: n_parts(:)       ! num partitions        (l_max)
    integer,   allocatable :: adapt_bnd(:)     ! boundaries to adapt   (n_bnd)
    real(RNP), allocatable :: adapt_box(:,:,:) ! boxes to be adapted (3,2,n_box)

    namelist /ml_mesh_options_3d__dynamic/ refinement
    namelist /ml_mesh_options_3d__dynamic/ n_parts
    namelist /ml_mesh_options_3d__dynamic/ adapt_bnd, adapt_box

    ! auxiliary variables ......................................................

    integer :: stat
    integer :: l

    ! read static input variables ..............................................

    rewind(unit)
    read(unit, nml = ml_mesh_options_3d__static, iostat = stat)

    if (stat == 0) then
      call this % SetUp(l_top, l_max, l_adapt, n_bnd, n_box)
    else if (stat == iostat_end) then
      call this % SetUp(l_top = 1)
      this % partition % n_parts = n_proc
      return
    else
      call Error( 'Read_ML_Mesh_Options_3D' &
                , 'failed reading static ML mesh options' &
                , 'ML__Mesh__3D')
    end if

    ! read dynamic input variables .............................................

    allocate(refinement(this%l_max-1), source = ' ')
    allocate(n_parts(this%l_max),      source = -1 )
    allocate(adapt_bnd(n_bnd),         source = -1 )
    allocate(adapt_box(3,2,n_box),     source = huge(ONE))

    read(unit, nml = ml_mesh_options_3d__dynamic)

    ! create multilevel mesh options ...........................................

    ! partitioning options
    associate(partition => this % partition)

      partition % method = partitioner

      if (c_active >= 0) then
        partition % c_active = c_active
      end if
      if (c_frozen >= 0) then
        partition % c_frozen = c_frozen
      end if

      partition(1 ) % child = .false.
      partition(2:) % child = .true.
      partition(2:) % split = split

      if (all(n_parts > 0)) then
        ! number of partitions given for all levels
        partition(1:this%l_max) % n_parts = n_parts
      else if (n_parts_growth > 0) then
        ! number of partitions grows geometrically
        partition(1) % n_parts = max(n_parts_root, 1)
        do l = 2, this%l_max
          partition(l) % n_parts = n_parts_growth * partition(l-1)%n_parts
        end do
      else
        ! constant number of partitions
        partition(1:this%l_max) % n_parts = max(n_parts_root, 1)
      end if

      ! enforce upper limit for number of partitions
      if (present(n_proc)) then
        partition % n_parts = min(partition%n_parts, n_proc)
      end if

      partition % n_con_root  = n_con_root
      partition % n_con_child = n_con_child
      partition % n_con_sub   = n_con_sub

    end associate

    this % refinement = refinement
    this % adapt_bnd  = adapt_bnd
    this % adapt_box  = adapt_box

  end function Read_ML_Mesh_Options_3D

  !-----------------------------------------------------------------------------
  !> Broadcasting multilevel mesh options

  subroutine Bcast_ML_Mesh_Options_3D(this, root, comm)
    class(ML_Mesh_Options_3D), intent(inout) :: this !< options
    integer       , intent(in) :: root !< rank root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    integer :: l, rank, n_bnd, n_box

    associate( l_top   => this%l_top   &
             , l_max   => this%l_max   &
             , l_adapt => this%l_adapt )

      call MPI_Comm_rank(comm, rank)

      if (rank == root) then
        n_bnd = size(this % adapt_bnd,1)
        n_box = size(this % adapt_box,3)
      end if

      call XMPI_Bcast(l_top  , root, comm)
      call XMPI_Bcast(l_max  , root, comm)
      call XMPI_Bcast(l_adapt, root, comm)
      call XMPI_Bcast(n_bnd  , root, comm)
      call XMPI_Bcast(n_box  , root, comm)

      if (rank /= root) then
        call this % SetUp(l_top, l_max, l_adapt, n_bnd, n_box)
      end if

      if (l_max > 1) then
        call XMPI_Bcast(this % refinement, root, comm)
      end if

      if (n_bnd > 0) then
        call XMPI_Bcast(this % adapt_bnd , root, comm)
      end if

      if (n_box > 0) then
        call XMPI_Bcast(this % adapt_box , root, comm)
      end if

      do l = 1, l_max
        call this % partition(l) % Bcast( 0, comm )
      end do

    end associate

  end subroutine Bcast_ML_Mesh_Options_3D

  !=============================================================================

end module ML__Mesh__3D

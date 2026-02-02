module Mesh_Partitioner__3D
  use Kind_Parameters
  use Logging_Levels
  use Execution_Control
  use XMPI
  use Mesh__3D

  implicit none
  private

  public :: MeshPartitionerOptions_3D
  public :: MeshPartitioner_3D

  !-----------------------------------------------------------------------------
  !> Partitioner options
  !>
  !> The computational mesh can be partitioned either using
  !>   - a space filling curve, by choosing `method = 1`, or
  !>   - an abstract graph, with `method = 2`.
  !>
  !> Use `n_parts` to specify the requested number of partitions.
  !>
  !> Set `child = T` to partition the child mesh and `child = F` to repartition
  !> the given `mesh`.
  !>
  !> Set `split` to obtain child partitions by subdividing parent partitions.
  !> This option is only available with graph partitioning and `child = T`.
  !>
  !> Use `c_active` and `c_frozen` to specify the weights of active and frozen
  !> elements, respectively. While for `c_active` a minimum of `1` is enforced,
  !> `c_frozen` can be zero.
  !>
  !> Additional options apply to graph partitioning. In this case, the mesh
  !> elements assume the role of the graph vertices, whereas the connections
  !> to neighbor elements are identified as the graph edges. Depending on the
  !> adjacency type, the edges are classified into three groups: element face,
  !> element edge and element vertex connections.
  !>
  !> The option `w_adj` determines the weights of the graph edges resulting
  !> from element adjacency at the faces, edges and vertices. While the faces
  !> are always taken into account, edges and vertex connections are included
  !> only if the corresponding weights do not vanish.
  !>
  !> Further, the graph partitioning is subject to `n_con` constraints which
  !> depend on the target mesh (given or child) and the number of sublevels
  !> considered: For any element, let be `c` the number of children, `g` the
  !> number of grand children and `q` further descendants in sublevels 3 up to
  !> `n_con_sub`. Then the constraints are defined as follows:
  !>
  !>   | n_con |  given mesh    | child mesh  |
  !>   | ----- | -------------- | ----------- |
  !>   |   1   |  `1`           |  `c+g+q`    |
  !>   |   2   |  `1, c+g+q`    |  `c, g+q`   |
  !>   |   3   |  `1, c, g+q`   |  `c, g, q`  |
  !>   !   4   !  `1, c, g, q`  |             |

  type MeshPartitionerOptions_3D

    ! general options
    integer :: method   = 1       !< 1/2 for SFC/graph-based partitioner
    integer :: n_parts  = 1       !< num partitions requested
    integer :: c_active = 20      !< cost of active child elements
    integer :: c_frozen = 1       !< cost of frozen child elements
    logical :: child    = .false. !< T/F for partitioning child/current mesh
    logical :: split    = .false. !< T for subdivision mode, ignored if child=F

    ! graph partitioning options
    integer :: w_adj(3)    = [1,0,0] !< face/edge/vertex adjacency weights
    integer :: n_con_root  = 1       !< max num constraints for mesh itself
    integer :: n_con_child = 1       !< max num constraints for child mesh
    integer :: n_con_sub   = 10      !< max num sublevels to be weighted

  contains
    procedure :: Bcast => Bcast_PartitionerOptions
  end type MeshPartitionerOptions_3D

  !=============================================================================
  ! Module procedures

  interface

    !----------------------------------------------------------------------------
    !> Partitioner based on space filling curves

    module subroutine SFC_Partitioner(opt, mesh, tp_elem, n_parts)
      class(MeshPartitionerOptions_3D), intent(in) :: opt
      class(Mesh_3D), intent(in) :: mesh !< current mesh partition
      integer, intent(out) :: tp_elem(:) !< new/child element target partitions
      integer, intent(out) :: n_parts    !< actual number of partitions
    end subroutine SFC_Partitioner

    !---------------------------------------------------------------------------
    !> Graph-based mesh partitioning using Metis or ParMetis

    module subroutine Graph_Partitioner(opt, mesh, tp_elem, n_parts)
      class(MeshPartitionerOptions_3D), intent(in) :: opt
      class(Mesh_3D), intent(in) :: mesh !< current mesh partition
      integer, intent(out) :: tp_elem(:) !< new/child element target partitions
      integer, intent(out) :: n_parts    !< actual number of partitions
    end subroutine Graph_Partitioner

  end interface

contains

  !=============================================================================
  ! MeshPartitionerOptions_3D: type-bound procedures

  subroutine Bcast_PartitionerOptions(opt, root, comm)
    class(MeshPartitionerOptions_3D), intent(inout) :: opt !< options
    integer       , intent(in) :: root !< rank root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(opt % method      , root, comm)
    call XMPI_Bcast(opt % n_parts     , root, comm)
    call XMPI_Bcast(opt % c_active    , root, comm)
    call XMPI_Bcast(opt % c_frozen    , root, comm)
    call XMPI_Bcast(opt % child       , root, comm)
    call XMPI_Bcast(opt % split       , root, comm)
    call XMPI_Bcast(opt % w_adj       , root, comm)
    call XMPI_Bcast(opt % n_con_root  , root, comm)
    call XMPI_Bcast(opt % n_con_child , root, comm)
    call XMPI_Bcast(opt % n_con_sub   , root, comm)

  end subroutine Bcast_PartitionerOptions

  !=============================================================================
  ! Mesh partitioner

  !-----------------------------------------------------------------------------
  !> Computes the new or child element target partitions for the given mesh

  subroutine MeshPartitioner_3D(opt, mesh, tp_elem, n_parts)
    class(MeshPartitionerOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh !< current mesh partition
    integer, intent(out) :: tp_elem(:) !< new/child element target partitions
    integer, intent(out) :: n_parts    !< actual number of partitions

    if (opt % method == 1 .and. mesh % has_sfc) then
      call SFC_Partitioner(opt, mesh, tp_elem, n_parts)
    else
      call Graph_Partitioner(opt, mesh, tp_elem, n_parts)
    end if

  end subroutine MeshPartitioner_3D

  !=============================================================================

end module Mesh_Partitioner__3D

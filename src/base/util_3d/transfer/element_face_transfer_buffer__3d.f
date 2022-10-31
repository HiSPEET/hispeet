!> summary:  Type and methods for transferring linked element face data
!> author:   Joerg Stiller
!> date:     2021/09/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Element_Face_Transfer_Buffer__3D

  use Kind_Parameters
  use Constants        , only: ZERO, ONE
  use Execution_Control, only: Error
  use XMPI
  use Mesh_Link__3D
  use Mesh__3D

  implicit none
  private

  public :: ElementFaceTransferBuffer_3D

  ! made public to circumvent an error with gfortran*9/10:
  public :: Init_ElementFaceTransferBuffer_Shared_3D

  !-----------------------------------------------------------------------------
  !> Auxiliary structure for keeping linked element-face data and metadata

  type ElementFaceTransferData
    integer,           allocatable :: element(:)   !< coupled element ID
    integer,           allocatable :: face(:)      !< coupled element face
    class(*),          allocatable :: buf(:,:,:,:) !< message buffer
    type(MPI_Request), allocatable :: request(:)   !< requests
  end type ElementFaceTransferData

  !-----------------------------------------------------------------------------
  !> Type for transferring face data from master to ghost elements
  !>
  !> The type supports two flavors of mesh variables
  !>
  !>    -  single variable of the shape `v(np,np,6,ne+ng)`
  !>    -  variable arrays of the shape `v(np,np,6,ne+ng,nc)`
  !>
  !> where `v` is of type
  !>
  !>    -  `real(RDP)`  (others may be added on demand)
  !>
  !> and
  !>
  !>    -  `np` is the number of element points per direction
  !>    -  `ne` is the number of local elements, i.e. `mesh % n_elem`
  !>    -  `ng` is the number of ghost elements, i.e. `mesh % n_ghost`
  !>    -  `nc` is the number of components/variables (auto-detected)
  !>
  !> The third dimension of `v` refers to the element faces.
  !> It is recommended to pass the variable arguments always in the shape
  !> given above, though the ghost entries are needed only for `Merge`.
  !>
  !> Use with single thread (no OpenMP):
  !>
  !>    -  the `asynchronous` attribute is required for non-blocking MPI send
  !>       and receive operations in `Transfer`
  !>
  !>         type(ElementFaceTransferBuffer_3D), asynchronous :: buf
  !>         ...
  !>         ! create and fill buffer, start transfer
  !>         buf = ElementFaceTransferBuffer_3D(mesh, v)
  !>         call buf % Transfer(mesh, v, tag)
  !>         ...
  !>         ! possibly perform some computations to hide communication costs
  !>         ...
  !>         ! merge received master data into ghosts and finish transfer
  !>         call buf % Merge(v)
  !>
  !> Use with OpenMP:
  !>
  !>    -  the buffer must be declared `save` to be shared among the threads
  !>    -  it must be `allocatable` and (de)allocated explicitly to ensure
  !>       correct finalization and, thus, release of component storage
  !>
  !>         type(ElementFaceTransferBuffer_3D), &
  !>             asynchronous, allocatable, save :: buf
  !>         ...
  !>         !$omp master
  !>         buf = ElementFaceTransferBuffer_3D(mesh, v)
  !>         !$omp end master
  !>         !$omp barrier
  !>         ...
  !>         call buf % Transfer(mesh, v, tag)
  !>         ...
  !>         call buf % Merge(v)
  !>         ...
  !>         !$omp master
  !>         deallocate(buf)
  !>         !$omp end master

  type ElementFaceTransferBuffer_3D
    integer :: np = 0                      !< number of points per direction
    integer :: nc = 0                      !< number of components
    integer :: ne = 0                      !< number of local (master) elements
    integer :: ng = 0                      !< number of ghost elements
    integer :: nf = 0                      !< number of transferred faces
    integer, allocatable :: start(:)       !< message start addresses
    integer, allocatable :: len(:)         !< message lengths
    type(ElementFaceTransferData) :: send  !< local element-face data
    type(ElementFaceTransferData) :: recv  !< remote element-face data

  contains

    generic, public :: Transfer => Transfer_RDP_S, Transfer_RDP_A
    generic, public :: Merge    => Merge_RDP_S, Merge_RDP_A

    procedure, private :: Transfer_RDP_S, Transfer_RDP_A
    procedure, private :: Merge_RDP_S, Merge_RDP_A

  end type ElementFaceTransferBuffer_3D

  !=============================================================================
  ! Constructor

  interface ElementFaceTransferBuffer_3D

    !---------------------------------------------------------------------------
    !> New element face transfer buffer for a single real variable

    module function New_TransferBuffer_RDP_S(mesh, v) result(this)
      type(ElementFaceTransferBuffer_3D) :: this
      class(Mesh_3D), intent(in) :: mesh  !< mesh partition
      real(RDP), intent(in) :: v(:,:,:,:) !< element-face variable (np,np,6,*)
    end function New_TransferBuffer_RDP_S

    !---------------------------------------------------------------------------
    !> New element face transfer buffer for an array of real variables

    module function New_TransferBuffer_RDP_A(mesh, v) result(this)
      type(ElementFaceTransferBuffer_3D) :: this
      class(Mesh_3D), intent(in) :: mesh    !< mesh partition
      real(RDP), intent(in) :: v(:,:,:,:,:) !< element-face var (np,np,6,*,nc)
    end function New_TransferBuffer_RDP_A

  end interface

  !=============================================================================
  ! Transfer procedures

  interface

    !--------------------------------------------------------------------------
    !> Send local data to ghosts and remote data -- real scalar

    module subroutine Transfer_RDP_S(this, mesh, v, tag)
      class(ElementFaceTransferBuffer_3D), asynchronous, intent(inout) :: this
      class(Mesh_3D), intent(in) :: mesh       !< mesh partition
      real(RDP),      intent(in) :: v(:,:,:,:) !< element-face variable
      integer,        intent(in) :: tag        !< message tag
    end subroutine Transfer_RDP_S

    !---------------------------------------------------------------------------
    !> Send local data to ghosts and remote data -- real array

    module subroutine Transfer_RDP_A(this, mesh, v, tag)
      class(ElementFaceTransferBuffer_3D), asynchronous, intent(inout) :: this
      class(Mesh_3D), intent(in) :: mesh         !< mesh partition
      real(RDP),      intent(in) :: v(:,:,:,:,:) !< element-face variable
      integer,        intent(in) :: tag          !< message tag
    end subroutine Transfer_RDP_A

  end interface

  !=============================================================================
  ! Merge procedures

  interface

    !---------------------------------------------------------------------------
    !> Complete receive and merge buffer into local data -- real(RDP) scalar

    module subroutine Merge_RDP_S(this, v)
      class(ElementFaceTransferBuffer_3D), intent(inout) :: this
      real(RDP), intent(inout) :: v(:,:,:,:) !< element-face variable
    end subroutine Merge_RDP_S

    !---------------------------------------------------------------------------
    !> Complete receive and merge buffer into local data -- real(RDP) array

    module subroutine Merge_RDP_A(this, v)
      class(ElementFaceTransferBuffer_3D), intent(inout) :: this
      real(RDP), intent(inout) :: v(:,:,:,:,:) !< element-face variable
    end subroutine Merge_RDP_A

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Common initialization of element face transfer buffers

  subroutine Init_ElementFaceTransferBuffer_Shared_3D(this, mesh)
    class(ElementFaceTransferBuffer_3D), intent(inout) :: this !< buffer
    class(Mesh_3D), intent(in) :: mesh !< mesh partition

    integer :: e, i, j, f, l, m

    if (mesh % part   < 0) return  ! inactive
    if (mesh % n_link < 1) return  ! no links

    this % ne = mesh % n_elem
    this % ng = mesh % n_ghost
    this % nf = sum( mesh % link % n_face )

    allocate(this % send % element( this % nf ))
    allocate(this % send % face   ( this % nf ))
    allocate(this % recv % element( this % nf ))
    allocate(this % recv % face   ( this % nf ))

    associate(send => this % send, recv => this % recv)

      ! set side to linked face IDs
      i = 1
      do l = 1, mesh % n_link
      do m = 1, mesh % link(l) % n_face
        associate(face => mesh % link(l) % face(m))
          e = face % element_id
          f = face % element_face
          associate(element => mesh % element(e))
            j = element % face(f) % i_neighbor ! > 0: must exist
            ! local element
            send % element(i) = e
            send % face(i)    = f
            ! remote element
            recv % element(i) = element % neighbor(j) % id
            recv % face(i)    = element % neighbor(j) % component
          end associate
        end associate
        i = i + 1
      end do
      end do

    end associate

    allocate(this % send % request(mesh % n_link))
    allocate(this % recv % request(mesh % n_link))

  end subroutine Init_ElementFaceTransferBuffer_Shared_3D

  !=============================================================================

end module Element_Face_Transfer_Buffer__3D
